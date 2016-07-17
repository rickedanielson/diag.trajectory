#=
 = Perform a single or a series of analysis calibration and performance estimates,
 = with the single option employing all available collocations, or for a series of
 = estimates performed across a parameter space defined in terms of a set of
 = hypercube dimensions (e.g., surface current and lat/lon position, mainly because
 = these are readily available from in situ obs, inertial oscillation proxy, Lagrangian
 = or Eulerian parameters such as MLD, SSH/SST/SSS structure, bathymetry or nearness
 = to the coast, and possibly date or season).  For the hypercube option, fixed-size
 = subsets of available collocations are selected based on geometrical closeness to
 = the target parameters (by equating units and ignoring pdf shape, for example) and
 = employed to obtain each calibration and performance estimate.  Target parameters
 = are then varied with each subset, yielding the actual mean parameters by a simple
 = average (and where two averages are close, these should be combined).  The ungridded
 = variations in analysis quality are then gridded in order to permit a lookup table
 = for maps of analysis quality that depend on these parameters.  Simple polynomials
 = (whose coefficients are found by least squares) are employed to regrid - RD May 2016.
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
const RECALIB          = true                           # perform an affine recalibration
const GLOBAL           = false                          # employ all collocations or a targeted subsets
const ANALYS           = 1                              # number of current analyses
const HYPER            = 1                              # number of hypercube dimensions (independent variables)

const DIRS  = ["v2.0_global_025_deg_total_15m"]
const UCUAC = [     0.00312866]
const UCUBC = [     1.45594908]
const VCUAC = [     0.00055068]
const VCUBC = [     1.48561760]
const UCUAV = [     0.00387766]
const UCUBV = [     1.30775895]
const VCUAV = [     0.00135119]
const VCUBV = [     1.32680928]

const UCLAC = [     0.04417926     -0.10298563      0.02980934]
const UCLBC = [     4.01156863     -6.30790060      3.07021902]
const VCLAC = [    -0.05374986      0.06217939     -0.01455402]
const VCLBC = [     1.60881144     -3.01035929      2.12685222]
const UCLAV = [    -0.11382001      0.08965406      0.00104184]
const UCLBV = [     2.98988430     -4.72762832      2.52096027]
const VCLAV = [    -0.02881019      0.00653482      0.01137809]
const VCLBV = [     4.51764345     -6.78495552      3.07755129]

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb\n\n")
  exit(1)
end

#=
 = Function returning triple collocation cal/val measures for a group of analyses, following McColl
 = et al. (2014).  Inputs are an array of collocated values and stats are returned for each analysis,
 = where it is assumed that extrapolation from before and after is done using the same analysis, so
 = no consideration of relative effective resolution (and no iteration as in Vogelzang et al. 2011)
 = is necessary (i.e., in situ is highest resolution, but there is no representativeness error
 = associated with one analysis being intermediate resolution and another being even lower resolution).
 =#

function triple(curr::Array{Float64,3})
  allalp = Array(Float64, ANALYS)
  allbet = Array(Float64, ANALYS)
  allsig = Array(Float64, ANALYS)
  allcor = Array(Float64, ANALYS)
  allmas = Array(Float64, ANALYS, HYPER)

  for a = 1:ANALYS                                                            # get the parametric center of mass
    mask = masquextreme(curr[1,   :,ANALYS+1], SDTRIM) &                      # that defines each subset in terms of
           masquextreme(curr[1,   :,       a], SDTRIM) &                      # hypercube independent variable (after 
           masquextreme(curr[2,   :,       a], SDTRIM)                        # trimming extreme values first)
    sampsitu    =       curr[1,mask,ANALYS+1]
    samprefa    =       curr[1,mask,       a]
    samprefb    =       curr[2,mask,       a]
    allmas[a,1] =  mean(curr[2,mask,ANALYS+1])

    avg1 = mean(sampsitu)                                                     # and use a robust calculation of covariance
    avg2 = mean(samprefa)                                                     # (two-pass here, but more algorithms are at
    avg3 = mean(samprefb)                                                     # en.wikipedia.org/wiki/Algorithms_for_calculating_variance)
    cv11 = mean((sampsitu - avg1) .* (sampsitu - avg1))
    cv12 = mean((sampsitu - avg1) .* (samprefa - avg2))
    cv13 = mean((sampsitu - avg1) .* (samprefb - avg3))
    cv22 = mean((samprefa - avg2) .* (samprefa - avg2))
    cv23 = mean((samprefa - avg2) .* (samprefb - avg3))
    cv33 = mean((samprefb - avg3) .* (samprefb - avg3))

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

    allalp[a] = 0.5 * (alp2 + alp3)
    allbet[a] = 0.5 * (bet2 + bet3)
    allsig[a] = 0.5 * (sig2 + sig3)
    allcor[a] = 0.5 * (cor2 + cor3)
  end

  return(allmas, allalp, allbet, allsig, allcor)                              # then return the average stats
end

#=
 = main program
 =#

const MISS             = -9999.0                        # generic missing value
if GLOBAL || RECALIB
  const RANGE          = 0.0:10.0: 0.0                  # target sampling range for current speed
  const CUTOFF         = 4500000000                     # number of collocations in a subset
else
  const RANGE          = 0.1: 0.1: 1.3                  # target sampling range for current speed
  const CUTOFF         = 200                            # number of collocations in a subset
end

const ALPH             = 1                              # error model x = ALPH + BETA * truth + error
const BETA             = 2                              # error model x = ALPH + BETA * truth + error
const SIGM             = 3                              # triple coll RMSE
const CORR             = 4                              # triple coll correlation coefficient
const MSPD             = 5                              # center-of-mass parameter
const PARAMS           = 5                              # number of triple collocation parameters

contains(ARGS[1], "calib.ucur") && (alph = UCUAV ; beta = UCUBV ; alpl = UCLAV ; betl = UCLBV ; varname = "UC..C")
contains(ARGS[1], "calib.vcur") && (alph = VCUAV ; beta = VCUBV ; alpl = VCLAV ; betl = VCLBV ; varname = "VC..C")
contains(ARGS[1], "valid.ucur") && (alph = UCUAC ; beta = UCUBC ; alpl = UCLAC ; betl = UCLBC ; varname = "UC..V")
contains(ARGS[1], "valid.vcur") && (alph = VCUAC ; beta = VCUBC ; alpl = VCLAC ; betl = VCLBC ; varname = "VC..V")
contains(ARGS[1],      ".ucur") && (uvuv = replace(ARGS[1], "ucur", "vcur"))
contains(ARGS[1],      ".vcur") && (uvuv = replace(ARGS[1], "vcur", "ucur"))
contains(ARGS[1], "calib.ucur") && (alpu = UCUAV ; betu = UCUBV ; alpv = VCUAV ; betv = VCUBV)
contains(ARGS[1], "calib.vcur") && (alpu = VCUAV ; betu = VCUBV ; alpv = UCUAV ; betv = UCUBV)
contains(ARGS[1], "valid.ucur") && (alpu = UCUAC ; betu = UCUBC ; alpv = VCUAC ; betv = VCUBC)
contains(ARGS[1], "valid.vcur") && (alpu = VCUAC ; betu = VCUBC ; alpv = UCUAC ; betv = UCUBC)

fpa = My.ouvre(ARGS[1], "r") ; lines = readlines(fpa) ; close(fpa)            # read and count all triple collocations
fpa = My.ouvre(   uvuv, "r") ; linez = readlines(fpa) ; close(fpa)            # (including both current components) and
linum = length(lines) ; linuz = length(linez)                                 # allocate for the target parameters and
if linum != linuz
  print("\nERROR: number of lines in $(ARGS[1]) and $uvuv are $linum != $linuz\n\n")
  exit(-1)
end
dist = zeros(linum)                                                           # distance to them, the resulting mean params,
chnk = linum < CUTOFF ? linum : CUTOFF                                        # and the triple collocation cal/val estimates
curr = zeros(2, chnk, ANALYS + 1)
calval = [MISS for a = RANGE, b = 1:ANALYS, c = 1:PARAMS]

for a = 1:ANALYS                                                              # delete files so as to append below
  filename = ARGS[1] * "." * DIRS[a] * ".cali.locl"
  isfile(filename) && (print("rm $filename\n") ; rm(filename))
end

for (a, rana) in enumerate(RANGE)                                             # loop through the target parameters and
  for d = 1:linum                                                             # isolate the nearest CUTOFF set of obs
    vala = float(split(lines[d]))
    valb = float(split(linez[d]))
    cspd = (vala[OCUR]^2.0 + valb[OCUR]^2.0)^0.5
    dist[d] = abs(rana - cspd)
  end
  lims = sort(dist)[chnk]

  e = 1                                                                       # get cal/val parameters for this subset
  for d = 1:linum                                                             # (possibly after recalibrating)
    if dist[d] <= lims && e <= chnk
      vala = float(split(lines[d]))
      valb = float(split(linez[d]))
      cspd = (vala[OCUR]^2.0 + valb[OCUR]^2.0)^0.5
      curr[1,e,:] = [vala[TOTB] vala[OCUR]]
      curr[2,e,:] = [vala[TOTA]       cspd]
      if RECALIB
        if GLOBAL
          for f = 1:ANALYS
            curr[1,e,f] = (curr[1,e,f] - alph[f]) / beta[f]
            curr[2,e,f] = (curr[2,e,f] - alph[f]) / beta[f]
          end
        else
          for f = 1:ANALYS
cspd = (((vala[TOTB] - alpu[f]) / betu[f])^2.0 + ((valb[TOTB] - alpv[f]) / betv[f])^2.0)^0.5
            localph = alpl[1] * cspd^2 + alpl[2] * cspd + alpl[3]
            locbeta = betl[1] * cspd^2 + betl[2] * cspd + betl[3]
            curr[1,e,f] = (curr[1,e,f] - localph) / locbeta
cspd = (((vala[TOTA] - alpu[f]) / betu[f])^2.0 + ((valb[TOTA] - alpv[f]) / betv[f])^2.0)^0.5
            localph = alpl[1] * cspd^2 + alpl[2] * cspd + alpl[3]
            locbeta = betl[1] * cspd^2 + betl[2] * cspd + betl[3]
            curr[2,e,f] = (curr[2,e,f] - localph) / locbeta
          end
        end
      end
      e += 1
    end
  end
  (allmas, allalp, allbet, allsig, allcor) = triple(curr)

  calval[a,:,ALPH] = allalp
  calval[a,:,BETA] = allbet
  calval[a,:,SIGM] = allsig
  calval[a,:,CORR] = allcor
  calval[a,:,MSPD] = allmas[:,1]

  if GLOBAL || RECALIB
    fpb = My.ouvre(ARGS[1] * ".cali.glob", "w")
    form = @sprintf("const %sUA%c = [%15.8lf]\n", varname[1:2], varname[5], allalp[1])
    write(fpb, form)
    form = @sprintf("const %sUB%c = [%15.8lf]\n", varname[1:2], varname[5], allbet[1])
    write(fpb, form)
    form = @sprintf("\ntarget param   CSPD is %6.2f\n", rana)
    write(fpb, form)
    form = @sprintf("  mean param   CSPD is %6.2f\n\n", mean(allmas[:,1]))
    write(fpb, form)
    form = @sprintf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
    write(fpb, form)
    for d = 1:ANALYS
      form = @sprintf("%33s %8.3f %8.3f %8.3f %8.3f\n", DIRS[d], allalp[d], allbet[d], allsig[d], allcor[d])
      write(fpb, form)
    end
    close(fpb)
  else
    for d = 1:ANALYS
      fpb = My.ouvre(ARGS[1] * "." * DIRS[d] * ".cali.locl", "a")
      form = @sprintf("%33s %15s %15s %8s %8s %8s %8s\n",
        "analysis", "target value", "mean value", "allalp", "allbet", "allsig", "allcor")
      a == 1 && write(fpb, form)
      form = @sprintf("%33s %15.3f %15.3f %8.3f %8.3f %8.3f %8.3f\n",
        DIRS[d], rana, mean(allmas[d,1]), allalp[d], allbet[d], allsig[d], allcor[d])
      write(fpb, form)
      close(fpb)
    end
  end

  @printf("cala = [%15.8lf]\n", allalp[1])
  @printf("calb = [%15.8lf]\n", allbet[1])
  @printf("\ntarget param   CSPD is %6.2f\n", rana)
  @printf("  mean param   CSPD is %6.2f\n\n", mean(allmas[:,1]))
  @printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
  for d = 1:ANALYS
    @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", DIRS[d], allalp[d], allbet[d], allsig[d], allcor[d])
  end
end

if GLOBAL || RECALIB  exit(0)  end

varspd = Array(Float64, 0)                                                    # the sqerror closure requires data
varcol = Array(Float64, 0)                                                    # arrays in global scope

function sqerror(coef::Array{Float64,1})                                      # define the least squares metric
  err = 0.0
  for i in 1:length(varspd)
    res  = coef[1] * varspd[i]^2 + coef[2] * varspd[i] + coef[3]
    err += (varcol[i] - res)^2
  end
  return err
end

for a = 1:ANALYS                                                              # for each analysis, store the dependence
  fpa = My.ouvre(ARGS[1] * "." * DIRS[a] * ".cali.locl", "a")                 # of the four triple-collocation statistics
  for b = 1:4                                                                 # in terms of second-order poynomial coeffs
    varspd = Array(Float64, 0)                                                # obtained using unconstrained optimization
    varcol = Array(Float64, 0)                                                # (so the stats' "hypercubes" are reduced to
    for (c, ranc) in enumerate(RANGE)                                         # variations on a multidimensional curve, but
      push!(varspd, calval[c,a,MSPD])                                         # with values everywhere in parameter space)
      push!(varcol, calval[c,a,b])
    end
    res = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)

    b == 1 && (stra = "const" ; strb = "A")
    b == 2 && (stra = "const" ; strb = "B")
    b == 3 && (stra = "     " ; strb = "S")
    b == 4 && (stra = "     " ; strb = "C")
    line = @sprintf("%s %sL%s%c = [%15.8f %15.8f %15.8f]\n",
      stra, varname[1:2], strb, varname[5], res.minimum[1], res.minimum[2], res.minimum[3])
    write(fpa, line) ; print("$(DIRS[a]) $b $(show(res)) \n")
  end
  close(fpa)
end

exit(0)
