#=
 = Perform a single or a series of analysis calibration and performance estimates,
 = with the series performed across a parameter space defined in terms of surface
 = current and collocation position (lat/lon; mainly because these are available
 = from buoy obs, although Lagrangian or Eulerian parameters such as MLD, SSH/SST/
 = SSS structure, bathymetry or nearness to the coast, and possibly date or season
 = could be employed).  Fixed-size subsets of available collocations are selected
 = based on geometrical closeness to the target parameters (by equating units and
 = ignoring pdf shape, for example) and employed to obtain each calibration and
 = performance estimate.  Target parameters are then varied with each subset,
 = yielding the actual mean parameters by a simple average (and where two averages
 = are close, these are combined).  The ungridded variations in analysis quality
 = are then gridded in order to permit a lookup table for maps of analysis quality
 = that depend on these parameters.  Simple polynomials (whose coefficients are
 = found by least squares) are employed to regrid - RD May 2016.
 =#

using My, Optim
const ODAT             = 1                              # identify indecies of the input data:
const OLAT             = 2                              # date/lat/lon on the collocation grid
const OLON             = 3
const OCUR             = 4                              # drifter current component as well as
const TOTB             = 11                             # the two total current estimates from
const TOTA             = 12                             # before and after

const FRAC             = 0.9                            # fractional update during iterations
const DELTA            = 0.001                          # generic convergence criterion
const SDTRIM           = 6.0                            # standard deviation trimming limit
const STORAGE          = true                           # retain the triple calibration metrics in a file
const ITERATE          = false                          # iterate on representativeness and calibration
const RECALIB          = false                          # perform an affine recalibration
const GLOBAL           = true                           # with STORAGE = ITERATE = true and RECALIB = false
const ANALYS           = 1                              # number of current analyses

const DIRS  = ["v2.0_global_025_deg_total_15m"]
if ITERATE
const UCUAC = [ -9999.00000000]
const UCUBC = [ -9999.00000000]
const VCUAC = [ -9999.00000000]
const VCUBC = [ -9999.00000000]

const UCUAV = [ -9999.00000000]
const UCUBV = [ -9999.00000000]
const VCUAV = [ -9999.00000000]
const VCUBV = [ -9999.00000000]
else
const UCUAC = [     1.00000000]
const UCUBC = [     1.00000000]
const VCUAC = [     1.00000000]
const VCUBC = [     1.00000000]

const UCUAV = [     1.00000000]
const UCUBV = [     1.00000000]
const VCUAV = [     1.00000000]
const VCUBV = [     1.00000000]
end
const UCURC = [     1.00000000]
const VCURC = [     1.00000000]
const UCURV = [     1.00000000]
const VCURV = [     1.00000000]

if size(ARGS) != (1,)
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb\n\n")
  exit(1)
end

#=
 = Function returning triple collocation cal/val measures for a group of analyses, following McColl
 = et al. (2014).  Inputs are an array of collocated values and an array of measures of short timescale
 = variability (rsqr).  Here we loop over many analyses (e.g., given by two types of extrapolations: one
 = analysis extrapolated from before and the other analyses extrapolated from after each date, and vice
 = versa, with the third dataset being in situ data that are valid and independent on the date).  Stats
 = are returned for each analysis, with possible consideration (using rsqr) of effective resolution after
 = a rescaling, given that the analysis of interest may be either low or intermediate resolution compared
 = to the other analysis and in situ.  Effective resolution is then monitored, as only at low resolution
 = does rsqr enter into the equations and require iteration [i.e., when the analysis of interest has the
 = lowest resolution, the covariance between the other analysis and in situ should be reduced by rsqr
 = (McColl et al. 2014); this impacts rescaling, and in turn, rsqr depends on rescaling (Vogelzang et al.
 = 2011, Eq. 10), so iteration is required as long as the analysis of interest remains lowest resolution;
 = the rsqr ranking is checked to ensure that the effective resolution ranking does not change; if the
 = analysis of interest is intermediate resolution, then its RMSE does not depend on rsqr and no iteration
 = is required (Vogelzang and Stoffelen 2012 Eq. 2.25)] - RD September, October 2015, February 2016.
 =#

function triple(curr::Array{Float64,3}, rsqr::Array{Float64,1})
  allalpa = Array(Float64, ANALYS, ANALYS)                                    # use in situ and (a) any one analysis
  allbeta = Array(Float64, ANALYS, ANALYS)                                    # in turn as references and store the
  allsiga = Array(Float64, ANALYS, ANALYS)                                    # third (b) error and correlation values
  allcora = Array(Float64, ANALYS, ANALYS)                                    # (do this for both extrapolations)
  allmasa = Array(Float64, ANALYS, ANALYS, 3)
  allalpb = Array(Float64, ANALYS, ANALYS)
  allbetb = Array(Float64, ANALYS, ANALYS)
  allsigb = Array(Float64, ANALYS, ANALYS)
  allcorb = Array(Float64, ANALYS, ANALYS)
  allmasb = Array(Float64, ANALYS, ANALYS, 3)
  allalpc = Array(Float64, ANALYS)
  allbetc = Array(Float64, ANALYS)
  allsigc = Array(Float64, ANALYS)
  allcorc = Array(Float64, ANALYS)
  allmass = Array(Float64, ANALYS, 3)

  for a = 1:ANALYS
    for b = 1:ANALYS                                                          # in addition to the "now" in situ obs,
      mask = masquextreme(curr[1,:,ANALYS+1], SDTRIM) &                       # use bef analysis "a" and aft analysis "b"
             masquextreme(curr[1,:,       a], SDTRIM) &                       # (having removed collocations that are beyond
             masquextreme(curr[2,:,       b], SDTRIM)                         # SDTRIM standard deviations from their mean)
      sampbuoy = curr[1,mask,ANALYS+1]                                        # and iterate if "b" is higher resolution
      sampsate = curr[1,mask,       a]                                        # then get the parametric center of mass of
      sampfore = curr[2,mask,       b]                                        # the resulting subset using its buoy values
      sampairt = mean(curr[1,mask,ANALYS+2])
      sampwspd = mean(curr[2,mask,ANALYS+1])
      sampsstt = mean(curr[2,mask,ANALYS+2])
      allmasa[a,b,:] = [sampairt sampwspd sampsstt]
      if GLOBAL  @show a, b  end

      deltasqr = rsqr[b] > rsqr[a] ? rsqr[b] - rsqr[a] : 0.0
      bet2 = bet3 = 1.0
      alp2 = alp3 = 0.0
#=
      flag = ITERATE
      while flag == true
        avg1 = mean(sampbuoy)
        avg2 = mean(sampsate)
        avg3 = mean(sampfore)
        cv11 = mean(sampbuoy.*sampbuoy) - avg1^2
        cv12 = mean(sampbuoy.*sampsate) - avg1 * avg2 - deltasqr              # write("cv12 = $cv12\n")
        cv13 = mean(sampbuoy.*sampfore) - avg1 * avg3
        cv22 = mean(sampsate.*sampsate) - avg2^2
        cv23 = mean(sampsate.*sampfore) - avg2 * avg3
        cv33 = mean(sampfore.*sampfore) - avg3^2

        alp2old = alp2
        alp3old = alp3
        bet2old = bet2
        bet3old = bet3
        subfrac = FRAC
        bet2 = subfrac * bet2old + (1.0 - subfrac) * (cv23 / cv13)            # write("bet2 = $bet2\n")
        bet3 = subfrac * bet3old + (1.0 - subfrac) * (cv23 / cv12)            # write("bet3 = $bet3\n")
        alp2 = subfrac * alp2old + (1.0 - subfrac) * (avg2 - bet2 * avg1)     # write("alp2 = $alp2\n")
        alp3 = subfrac * alp3old + (1.0 - subfrac) * (avg3 - bet3 * avg1)     # write("alp3 = $alp3\n")
        sampsate = (curr[1,mask,a] - alp2) / bet2                             # rescale analysis current
        sampfore = (curr[2,mask,b] - alp3) / bet3
        rsqrsate =      rsqr[a]         / bet2 / bet2
        rsqrfore =      rsqr[b]         / bet3 / bet3
#       print("cv23 $cv23 cv13 $cv13 cv23 / cv13 $(cv23 / cv13)\n")
#       print("loop-1 $a $b rsqr[a] $(rsqr[a]) / bet2 $bet2 = rsqrsate $rsqrsate   rsqr[b] $(rsqr[b]) / bet3 $bet3 = rsqrfore $rsqrfore\n")
#       print("$a $b rsqr[a] $(rsqr[a]) / bet2 $bet2 = rsqrsate $rsqrsate   rsqr[b] $(rsqr[b]) / bet3 $bet3 = rsqrfore $rsqrfore\n")
#       @printf("%d %d rsqr[a] %9.3f / bet2 %9.3f = rsqrsate %9.3f   rsqr[b] %9.3f / bet3 %9.3f = rsqrfore %9.3f\n",
#         a, b, rsqr[a], bet2, rsqrsate, rsqr[b], bet3, rsqrfore)

        deltaold = deltasqr
        deltasqr = rsqrfore > rsqrsate ? rsqrfore - rsqrsate : 0.0
        if abs(deltasqr - deltaold) < DELTA  flag = false  end
        if rsqr[a] > rsqr[b] || a == b  flag = false  end
      end
=#
      avg1 = mean(sampbuoy)
      avg2 = mean(sampsate)
      avg3 = mean(sampfore)
      cv11 = mean(sampbuoy.*sampbuoy) - avg1^2
      cv12 = mean(sampbuoy.*sampsate) - avg1 * avg2 - deltasqr
      cv13 = mean(sampbuoy.*sampfore) - avg1 * avg3
      cv22 = mean(sampsate.*sampsate) - avg2^2
      cv23 = mean(sampsate.*sampfore) - avg2 * avg3
      cv33 = mean(sampfore.*sampfore) - avg3^2

      bet2 = cv23 / cv13
      bet3 = cv23 / cv12
      alp2 = avg2 - bet2 * avg1
      alp3 = avg3 - bet3 * avg1

      tmpval = cv11 - cv12 * cv13 / cv23 ; sig1 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv22 - cv12 * cv23 / cv13 ; sig2 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv33 - cv13 * cv23 / cv12 ; sig3 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv12 * cv13 / cv11 / cv23 ; cor1 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv12 * cv23 / cv22 / cv13 ; cor2 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv13 * cv23 / cv33 / cv12 ; cor3 = tmpval > 0 ? sqrt(tmpval) : 0.0

      allalpa[a,b] = alp3
      allbeta[a,b] = bet3
      allsiga[a,b] = sig3
      allcora[a,b] = cor3
#=
      mask = masquextreme(curr[1,:,ANALYS+1], SDTRIM) &                       # use aft analysis "a" and bef analysis "b"
             masquextreme(curr[2,:,       a], SDTRIM) &                       # (having removed collocations that are beyond
             masquextreme(curr[1,:,       b], SDTRIM)                         # SDTRIM standard deviations from their mean)
      sampbuoy = curr[1,mask,ANALYS+1]                                        # and iterate if "b" is higher resolution
      sampsate = curr[2,mask,       a]                                        # then get the parametric center of mass of
      sampfore = curr[1,mask,       b]                                        # the resulting subset using its buoy values
      sampairt = mean(curr[1,mask,ANALYS+2])
      sampwspd = mean(curr[2,mask,ANALYS+1])
      sampsstt = mean(curr[2,mask,ANALYS+2])
      allmasb[a,b,:] = [sampairt sampwspd sampsstt]

      deltasqr = rsqr[b] > rsqr[a] ? rsqr[b] - rsqr[a] : 0.0
      bet2 = bet3 = 1.0
      alp2 = alp3 = 0.0

      flag = ITERATE
      while flag == true
        avg1 = mean(sampbuoy)
        avg2 = mean(sampsate)
        avg3 = mean(sampfore)
        cv11 = mean(sampbuoy.*sampbuoy) - avg1^2
        cv12 = mean(sampbuoy.*sampsate) - avg1 * avg2 - deltasqr              # write("cv12 = $cv12\n")
        cv13 = mean(sampbuoy.*sampfore) - avg1 * avg3
        cv22 = mean(sampsate.*sampsate) - avg2^2
        cv23 = mean(sampsate.*sampfore) - avg2 * avg3
        cv33 = mean(sampfore.*sampfore) - avg3^2

        alp2old = alp2
        alp3old = alp3
        bet2old = bet2
        bet3old = bet3
        subfrac = FRAC
        bet2 = subfrac * bet2old + (1.0 - subfrac) * (cv23 / cv13)            # write("bet2 = $bet2\n")
        bet3 = subfrac * bet3old + (1.0 - subfrac) * (cv23 / cv12)            # write("bet3 = $bet3\n")
        alp2 = subfrac * alp2old + (1.0 - subfrac) * (avg2 - bet2 * avg1)     # write("alp2 = $alp2\n")
        alp3 = subfrac * alp3old + (1.0 - subfrac) * (avg3 - bet3 * avg1)     # write("alp3 = $alp3\n")
        sampsate = (curr[2,mask,a] - alp2) / bet2                             # rescale analysis current
        sampfore = (curr[1,mask,b] - alp3) / bet3
        rsqrsate =      rsqr[a]         / bet2 / bet2
        rsqrfore =      rsqr[b]         / bet3 / bet3
#       print("cv23 $cv23 cv13 $cv13 cv23 / cv13 $(cv23 / cv13)\n")
#       print("loop-2 $a $b rsqr[a] $(rsqr[a]) / bet2 $bet2 = rsqrsate $rsqrsate   rsqr[b] $(rsqr[b]) / bet3 $bet3 = rsqrfore $rsqrfore\n")
#       print("$a $b rsqr[a] $(rsqr[a]) / bet2 $bet2 = rsqrsate $rsqrsate   rsqr[b] $(rsqr[b]) / bet3 $bet3 = rsqrfore $rsqrfore\n")
#       @printf("%d %d rsqr[a] %9.3f / bet2 %9.3f = rsqrsate %9.3f   rsqr[b] %9.3f / bet3 %9.3f = rsqrfore %9.3f\n",
#         a, b, rsqr[a], bet2, rsqrsate, rsqr[b], bet3, rsqrfore)

        deltaold = deltasqr
        deltasqr = rsqrfore > rsqrsate ? rsqrfore - rsqrsate : 0.0
        if abs(deltasqr - deltaold) < DELTA  flag = false  end
        if rsqr[a] > rsqr[b] || a == b  flag = false  end
      end

      avg1 = mean(sampbuoy)
      avg2 = mean(sampsate)
      avg3 = mean(sampfore)
      cv11 = mean(sampbuoy.*sampbuoy) - avg1^2
      cv12 = mean(sampbuoy.*sampsate) - avg1 * avg2 - deltasqr
      cv13 = mean(sampbuoy.*sampfore) - avg1 * avg3
      cv22 = mean(sampsate.*sampsate) - avg2^2
      cv23 = mean(sampsate.*sampfore) - avg2 * avg3
      cv33 = mean(sampfore.*sampfore) - avg3^2

      bet2 = cv23 / cv13
      bet3 = cv23 / cv12
      alp2 = avg2 - bet2 * avg1
      alp3 = avg3 - bet3 * avg1

      tmpval = cv11 - cv12 * cv13 / cv23 ; sig1 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv22 - cv12 * cv23 / cv13 ; sig2 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv33 - cv13 * cv23 / cv12 ; sig3 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv12 * cv13 / cv11 / cv23 ; cor1 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv12 * cv23 / cv22 / cv13 ; cor2 = tmpval > 0 ? sqrt(tmpval) : 0.0
      tmpval = cv13 * cv23 / cv33 / cv12 ; cor3 = tmpval > 0 ? sqrt(tmpval) : 0.0

      allalpb[a,b] = alp3
      allbetb[a,b] = bet3
      allsigb[a,b] = sig3
      allcorb[a,b] = cor3
=#
    end
  end

# for b = 1:ANALYS
#   tmpa = 0.0 ; for a = 1:ANALYS  if a != b  tmpa += allalpa[a,b]    + allalpb[a,b]    end  end ; allalpc[b]   = 0.5 * tmpa / (ANALYS - 1.0)
#   tmpa = 0.0 ; for a = 1:ANALYS  if a != b  tmpa += allbeta[a,b]    + allbetb[a,b]    end  end ; allbetc[b]   = 0.5 * tmpa / (ANALYS - 1.0)
#   tmpa = 0.0 ; for a = 1:ANALYS  if a != b  tmpa += allsiga[a,b]    + allsigb[a,b]    end  end ; allsigc[b]   = 0.5 * tmpa / (ANALYS - 1.0)
#   tmpa = 0.0 ; for a = 1:ANALYS  if a != b  tmpa += allcora[a,b]    + allcorb[a,b]    end  end ; allcorc[b]   = 0.5 * tmpa / (ANALYS - 1.0)
#   tmpa = 0.0 ; for a = 1:ANALYS  if a != b  tmpa += allmasa[a,b,:] .+ allmasb[a,b,:]  end  end ; allmass[b,:] = 0.5 * tmpa / (ANALYS - 1.0)
# end
  allalpc[1]   = allalpa[1,1]
  allbetc[1]   = allbeta[1,1]
  allsigc[1]   = allsiga[1,1]
  allcorc[1]   = allcora[1,1]
  allmass[1,:] = allmasa[1,1,:]

  return(allmass, allalpc, allbetc, allsigc, allcorc)                # then return the average stats
end

#=
 = main program
 =#

const MISS             = -9999.0                        # generic missing value
if GLOBAL
  const RANGA          =   0.0:10.0: 0.0                # target sampling range for AIRT dimension
  const RANGB          =   0.0:10.0: 0.0                # target sampling range for WSPD dimension
  const RANGC          =   0.0:10.0: 0.0                # target sampling range for SSTT dimension
  const CUTOFF         = 4500000000                     # number of collocations in a subset
else
  const RANGA          = -40.0:10.0:30.0                # target sampling range for AIRT dimension
  const RANGB          =   0.0:10.0:30.0                # target sampling range for WSPD dimension
  const RANGC          =   0.0:10.0:30.0                # target sampling range for SSTT dimension
  const CUTOFF         = 5000                           # number of collocations in a subset
end

const ALPH             = 1                              # error model x = ALPH + BETA * truth + error
const BETA             = 2                              # error model x = ALPH + BETA * truth + error
const SIGM             = 3                              # triple coll RMSE
const CORR             = 4                              # triple coll correlation coefficient
const MAIR             = 5                              # center-of-mass parameter
const MSPD             = 6                              # center-of-mass parameter
const MSST             = 7                              # center-of-mass parameter
const PARAMS           = 7                              # number of triple coll statistics

contains(ARGS[1], "calib.ucur") && (alph = UCUAV ; beta = UCUBV ; rsqr = UCURC ; varname = "UCU.C")
contains(ARGS[1], "calib.vcur") && (alph = VCUAV ; beta = VCUBV ; rsqr = VCURC ; varname = "VCU.C")
contains(ARGS[1], "valid.ucur") && (alph = UCUAC ; beta = UCUBC ; rsqr = UCURV ; varname = "UCU.V")
contains(ARGS[1], "valid.vcur") && (alph = VCUAC ; beta = VCUBC ; rsqr = VCURV ; varname = "VCU.V")

fpa = My.ouvre(ARGS[1], "r") ; lines = readlines(fpa) ; close(fpa)            # read and count all triple collocations
(linum,) = size(lines)                                                        # and allocate for the target parameters and
dist = zeros(linum)                                                           # distance to them, the resulting mean params,
chnk = linum < CUTOFF ? linum : CUTOFF                                        # and the triple collocation cal/val estimates
curr = zeros(2, chnk, ANALYS + 2)
target = [MISS for a = RANGA, b = RANGB, c = RANGC, d = 1:3]
calval = [MISS for a = RANGA, b = RANGB, c = RANGC, d = 1:ANALYS, e = 1:PARAMS]

for (a, rana) in enumerate(RANGA)                                             # loop through the target parameters and
  for (b, ranb) in enumerate(RANGB)                                           # isolate the nearest CUTOFF set of obs
    for (c, ranc) in enumerate(RANGC)
      for (d, line) in enumerate(lines)
        vals = float(split(line))
        dist[d] = (rana - vals[OCUR])^2.0 +
                  (ranb - vals[OLAT])^2.0 +
                  (ranc - cos(pi / 180.0 * vals[OLON]))^2.0
      end
      lims = sort(dist)[chnk]

      e = 1                                                                   # get cal/val parameters for this subset
      for (d, line) in enumerate(lines)                                       # (possibly after recalibrating)
        if dist[d] <= lims && e <= chnk
          vals = float(split(line))
          curr[1,e,:] = [vals[TOTB] vals[OCUR] vals[OCUR]]
          curr[2,e,:] = [vals[TOTA] vals[OLAT] cos(pi / 180.0 * vals[OLON])]
          if RECALIB
            for f = 1:ANALYS
              curr[1,e,f] = (curr[1,e,f] - alph[f]) / beta[f]
              curr[1,e,f] = (curr[1,e,f] - alph[f]) / beta[f]
            end
          end
          e += 1
        end
      end
      (allmas, allalp, allbet, allsig, allcor) = triple(curr, rsqr)

      target[a,b,c,:] = [rana ranb ranc]
      calval[a,b,c,:,ALPH] = allalp
      calval[a,b,c,:,BETA] = allbet
      calval[a,b,c,:,SIGM] = allsig
      calval[a,b,c,:,CORR] = allcor
      calval[a,b,c,:,MAIR] = allmas[:,1]
      calval[a,b,c,:,MSPD] = allmas[:,2]
      calval[a,b,c,:,MSST] = allmas[:,3]

      if STORAGE
        fpb = My.ouvre(ARGS[1] * ".cali", "w")
        form = @sprintf("const %sA%c = [%15.8lf,      0.0]\n", varname[1:3], varname[5], allalp[1])
        write(fpb, form)
        form = @sprintf("const %sB%c = [%15.8lf,      1.0]\n", varname[1:3], varname[5], allbet[1])
        write(fpb, form)
        form = @sprintf("\ntarget params   OCUR,OLAT,OLON are %6.2f %6.2f %6.2f\n",   rana, ranb, ranc)
        write(fpb, form)
        form = @sprintf("  mean params   OCUR,OLAT,OLON are %6.2f %6.2f %6.2f\n\n", mean(allmas[:,1]), mean(allmas[:,2]), mean(allmas[:,3]))
        write(fpb, form)
        form = @sprintf("%22s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
        write(fpb, form)
        for d = 1:ANALYS
          form = @sprintf("%22s %8.2f %8.2f %8.2f %8.2f\n", DIRS[d], allalp[d], allbet[d], allsig[d], allcor[d])
          write(fpb, form)
        end
        close(fpb)
      end

      @printf("cala = [%15.8lf]\n", allalp[1])
      @printf("calb = [%15.8lf]\n", allbet[1])
      @printf("\ntarget params   AIRT,WSPD,SSTT are %6.2f %6.2f %6.2f\n",   rana, ranb, ranc)
      @printf("  mean params   AIRT,WSPD,SSTT are %6.2f %6.2f %6.2f\n\n", mean(allmas[:,1]), mean(allmas[:,2]), mean(allmas[:,3]))
      @printf("%22s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
      for d = 1:ANALYS
        @printf("%22s %8.2f %8.2f %8.2f %8.2f\n", DIRS[d], allalp[d], allbet[d], allsig[d], allcor[d])
      end
    end
  end
end

if GLOBAL  exit(0)  end

varair = Array(Float64, 0)                                                    # the sqerror closure requires data
varspd = Array(Float64, 0)                                                    # arrays in global scope
varsst = Array(Float64, 0)
varcol = Array(Float64, 0)

function sqerror(coef::Array{Float64,1})                                      # define the least squares metric
  err = 0.0
  for i in 1:length(varair)
    res  = coef[1] * varair[i] * varair[i] + coef[2] * varspd[i] * varspd[i] + coef[3] * varsst[i] * varsst[i] +
           coef[4] * varair[i] * varspd[i] + coef[5] * varair[i] * varsst[i] + coef[6] * varspd[i] * varsst[i] +
           coef[7] * varair[i]             + coef[8] * varspd[i]             + coef[9] * varsst[i] + coef[10]
    err += (varcol[i] - res)^2
  end
  return err
end

for a = 1:ANALYS                                                              # for each analysis, store the dependence
  outfile = ARGS[1] * "." * DIRS[a]                                           # of the four triple-collocation statistics
  fpa = My.ouvre(outfile, "w")                                                # in terms of second-order poynomial coeffs
  for b = 1:4                                                                 # obtained using unconstrained optimization
    varair = Array(Float64, 0)                                                # (so the stats' "hypercubes" are reduced to
    varspd = Array(Float64, 0)                                                # variations on a multidimensional curve, but
    varsst = Array(Float64, 0)                                                # with values everywhere in parameter space)
    varcol = Array(Float64, 0)
    for (c, ranc) in enumerate(RANGA)
      for (d, rand) in enumerate(RANGB)
        for (e, rane) in enumerate(RANGC)
          push!(varair, calval[c,d,e,a,MAIR])
          push!(varspd, calval[c,d,e,a,MSPD])
          push!(varsst, calval[c,d,e,a,MSST])
          push!(varcol, calval[c,d,e,a,b])
        end
      end
    end
    res = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)

    line = @sprintf("%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n",
      res.minimum[1], res.minimum[2], res.minimum[3], res.minimum[4], res.minimum[5],
      res.minimum[6], res.minimum[7], res.minimum[8], res.minimum[9], res.minimum[10])
    write(fpa, line)
    print("$(DIRS[a]) $b $(show(res)) \n")
  end
  close(fpa)
end

exit(0)
