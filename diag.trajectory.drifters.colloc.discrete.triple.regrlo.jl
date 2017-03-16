#=
 = Perform a series of paired analysis validation (and performance) estimates by nonlinear (bin)
 = regression.  For each bin, calculate the bounding ordinary and reverse linear regression slope
 = parameters and all metrics used in the collocation technique (cf. McColl et al. 2014, Gruber et
 = al. 2016).
 = Each of the series involves a recalibration using the OLS slope and intercept from the opposing
 = collocations.  Fixed-size subsets of the available collocations are selected based on closeness
 = to a target value of the variable of interest (the actual mean value is obtained as a simple
 = average and if two averages are close, these might be combined).  The target value is then
 = varied over a reasonable range and calibration is obtained for each subset.  These variations
 = in calibration (as a function of the observed or analyzed variable of interest) are then used
 = to assess global performance - RD June, November 2016.
 =#

using My, Winston, RCall ; R"library(DetMCD)"
const ODAT             = 1                              # identify indecies of the input data (starting
const OLAT             = 2                              # with date/lat/lon/obs on the collocation grid)
const OLON             = 3
const OCUR             = 4
const MCDTRIM          = 0.95                           # Minimum Covariance Determinant trimming (nonoutlier percent)
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 2 && (argc = length(ARGS)) != 3
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb v2.0_global_025_deg_total_15m [.agulhas]\n\n")
  exit(1)
end
regionlab = "" ; length(ARGS) == 3 && (regionlab = ARGS[3])

shift = -1
ARGS[2] == "v2.0_global_025_deg_ekman_15m"   && (shift =  0)
ARGS[2] == "v2.0_global_025_deg_ekman_hs"    && (shift =  3)
ARGS[2] == "v2.0_global_025_deg_geostrophic" && (shift =  6)
ARGS[2] == "v2.0_global_025_deg_total_15m"   && (shift =  9)
ARGS[2] == "v2.0_global_025_deg_total_hs"    && (shift = 12)
if shift > -1
  const TOTB           = 5 + shift                      # identify the three analysis indecies
  const TOTN           = 6 + shift
  const TOTA           = 7 + shift
else
  print("\nERROR: $(ARGS[2]) is not recognized\n\n")
  exit(-1)
end
contains(ARGS[1], ".ucur") && (ARGS222 = replace(ARGS[1], "ucur", "vcur"))
contains(ARGS[1], ".vcur") && (ARGS222 = replace(ARGS[1], "vcur", "ucur"))

#=
 = Function returning three sets of calibration and performance metrics, following McColl et al. (2014) and Gruber
 = et al. (2016).  Two sets of metrics correspond to ordinary and reverse least squares regression and the third set
 = (somewhere between these two ordinary bounds) is defined by a parameter "confid" that is the relative confidence
 = in the uncalibrated error variance relative to the true variance.  Because the ordinary and reverse least squares
 = regressions bound all metrics, a bound on relative confidence is also checked.  First define the indecies of the
 = two-dimensional regression coefficient matrix that is returned.
 =#

const RCMR             = 3                              # number of regression coefficient matrix rows
const RCMC             = 7                              # number of regression coefficient matrix columns
const BOLS             = 1                              # ordinary  least squares solution bound
const BCON             = 2                              # confident least squares solution bound
const BRLS             = 3                              # reverse   least squares solution bound
const VART             = 1                              # variance in target/true parameter
const ALPH             = 2                              # additive       bias of uncalibrated parameter
const BETA             = 3                              # multiplicative bias of uncalibrated parameter
const SIGC             = 4                              # error variance in   calibrated parameter
const SIGU             = 5                              # error variance in uncalibrated parameter
const SNRC             = 6                              # signal to noise ratio of   calibrated parameter
const SNRU             = 7                              # signal to noise ratio of uncalibrated parameter
const DELT             = 1e-9                           # small number for comparison

function regress(refer::Array{Float64,1}, calib::Array{Float64,1}, uncal::Array{Float64,1}, confid::Float64)
  temp = [calib uncal]
  remp = rcopy(R"DetMCD($temp, alpha = $MCDTRIM)")
  mask = falses(length(temp[:,1])) ; for a in remp[:Hsubsets]  mask[a] = true  end

  avgr = mean(refer[mask])                                                    # avgr = mean(refer)  use all data (including outliers)
  avgc = remp[:center][1]                                                     # avgc = mean(calib)  and a robust calculation of covariance
  avgu = remp[:center][2]                                                     # avgu = mean(uncal)
  cvcc = remp[:cov][1,1]                                                      # cvcc = mean((calib - avgc) .* (calib - avgc))
  cvcu = remp[:cov][1,2]                                                      # cvcu = mean((calib - avgc) .* (uncal - avgr))
  cvuu = remp[:cov][2,2]                                                      # cvuu = mean((uncal - avgr) .* (uncal - avgr))

  confidlo = 0.0
  confidhi = cvuu / cvcc - cvcu^2 / cvcc^2
# if confid < confidlo  print("\nWARNING: resetting desired confidence $confid to $confidlo (the minimum) \n\n") ; confid = confidlo  end
# if confid > confidhi  print("\nWARNING: resetting desired confidence $confid to $confidhi (the maximum) \n\n") ; confid = confidhi  end

  rcm = Array(Float64, RCMR, RCMC)                                            # allocate the regression coefficient matrix
  beto =       cvcu / cvcc                                                    # and the slope parameter for all solutions
# betc = (cvuu + (cvuu^2 - 4.0 * confid * cvcu^2)^0.5) / 2.0 / cvcu           # including an intermediate solution based on
# betc = (cvuu - (cvuu^2 - 4.0 * confid * cvcu^2)^0.5) / 2.0 / cvcu
  betr =       cvuu / cvcu                                                    # variance ratio: uncalibrated-error / truth
#tmpval = cvuu^2 - 4.0 * confid * cvcu^2 ; betc = tmpval > 0 ? (cvuu + tmpval^0.5) / 2.0 / cvcu : cvuu / 2.0 / cvcu
#@show betc
#tmpval = cvuu^2 - 4.0 * confid * cvcu^2 ; betc = tmpval > 0 ? (cvuu - tmpval^0.5) / 2.0 / cvcu : cvuu / 2.0 / cvcu
#@show betc
  betc = 0.5 * cvuu / cvcu

  betc = (cvuu / cvcc)^0.5
  betc = betc < beto ? beto : betc
  betc = betc > betr ? betr : betc

  rcm[BOLS,VART] = cvcu / beto               # ***** NOTE variance here and SD below *****
  rcm[BCON,VART] = cvcu / betc
  rcm[BRLS,VART] = cvcu / betr
  rcm[BOLS,ALPH] = avgu - beto * avgc
  rcm[BCON,ALPH] = avgu - betc * avgc
  rcm[BRLS,ALPH] = avgu - betr * avgc

  rcm[BOLS,ALPH] = beto * (cvuu - beto * cvcu) / cvcu
  rcm[BCON,ALPH] = betc * (cvuu - betc * cvcu) / cvcu
  rcm[BRLS,ALPH] = betr * (cvuu - betr * cvcu) / cvcu

  rcm[BOLS,BETA] =        beto
  rcm[BCON,BETA] =        betc
  rcm[BRLS,BETA] =        betr
  tmpval = cvcc - cvcu / beto ; rcm[BOLS,SIGC] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvcc - cvcu / betc ; rcm[BCON,SIGC] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvcc - cvcu / betr ; rcm[BRLS,SIGC] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvuu - cvcu * beto ; rcm[BOLS,SIGU] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvuu - cvcu * betc ; rcm[BCON,SIGU] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvuu - cvcu * betr ; rcm[BRLS,SIGU] = tmpval > 0 ? sqrt(tmpval) : 0.0

  tmpval = cvcu / cvcc * beto ; rcm[BOLS,SNRC] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvcu / cvcc * betc ; rcm[BCON,SNRC] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvcu / cvcc * betr ; rcm[BRLS,SNRC] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvcu / cvuu * beto ; rcm[BOLS,SNRU] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvcu / cvuu * betc ; rcm[BCON,SNRU] = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cvcu / cvuu * betr ; rcm[BRLS,SNRU] = tmpval > 0 ? sqrt(tmpval) : 0.0

# tmpval = cvcu / cvcc * beto ; -DELT < tmpval < DELT && (tmpval = 0) ; -DELT < tmpval - 1 < DELT && (tmpval = 1) ; rcm[BOLS,SNRC] = tmpval > 0 ? sqrt(tmpval) : 0.0
# tmpval = cvcu / cvcc * betc ; -DELT < tmpval < DELT && (tmpval = 0) ; -DELT < tmpval - 1 < DELT && (tmpval = 1) ; rcm[BCON,SNRC] = tmpval > 0 ? sqrt(tmpval) : 0.0
# tmpval = cvcu / cvcc * betr ; -DELT < tmpval < DELT && (tmpval = 0) ; -DELT < tmpval - 1 < DELT && (tmpval = 1) ; rcm[BRLS,SNRC] = tmpval > 0 ? sqrt(tmpval) : 0.0
# tmpval = cvcu / cvuu * beto ; -DELT < tmpval < DELT && (tmpval = 0) ; -DELT < tmpval - 1 < DELT && (tmpval = 1) ; rcm[BOLS,SNRU] = tmpval > 0 ? sqrt(tmpval) : 0.0
# tmpval = cvcu / cvuu * betc ; -DELT < tmpval < DELT && (tmpval = 0) ; -DELT < tmpval - 1 < DELT && (tmpval = 1) ; rcm[BCON,SNRU] = tmpval > 0 ? sqrt(tmpval) : 0.0
# tmpval = cvcu / cvuu * betr ; -DELT < tmpval < DELT && (tmpval = 0) ; -DELT < tmpval - 1 < DELT && (tmpval = 1) ; rcm[BRLS,SNRU] = tmpval > 0 ? sqrt(tmpval) : 0.0
# tmpval = cvcu / (cvcc * beto - cvcu) ; rcm[BOLS,SNRC] = tmpval > 0 ? 10 * log(tmpval) : -9e9
# tmpval = cvcu / (cvcc * betc - cvcu) ; rcm[BCON,SNRC] = tmpval > 0 ? 10 * log(tmpval) : -9e9
# tmpval = cvcu / (cvcc * betr - cvcu) ; rcm[BRLS,SNRC] = tmpval > 0 ? 10 * log(tmpval) : -9e9
# tmpval = cvcu / (cvuu / beto - cvcu) ; rcm[BOLS,SNRU] = tmpval > 0 ? 10 * log(tmpval) : -9e9
# tmpval = cvcu / (cvuu / betc - cvcu) ; rcm[BCON,SNRU] = tmpval > 0 ? 10 * log(tmpval) : -9e9
# tmpval = cvcu / (cvuu / betr - cvcu) ; rcm[BRLS,SNRU] = tmpval > 0 ? 10 * log(tmpval) : -9e9
  return(avgr, rcm)
end

#=
 = main program
 =#

const RANGE            = 0.1:0.01:1.1                   # target sampling range for current speed
const CUTOFF           = 200                            # number of collocations in a subset

const MEMO             = 1                              # center-of-mass parameter
const MEMB             = 2                              # error model x = ALPH + BETA * truth + error
const MEMA             = 3                              # error model x = ALPH + BETA * truth + error
const MEMS             = 3                              # number of triple collocation members

const MASS             = 1                              # center-of-mass parameter (again...)
const ALPH             = 2                              # error model x = ALPH + BETA * truth + error
const BETA             = 3                              # error model x = ALPH + BETA * truth + error
const SIGM             = 4                              # triple coll RMSE
const CORR             = 5                              # triple coll correlation coefficient
const PARS             = 5                              # number of triple collocation parameters

ARGS333 = replace(ARGS[1], "calib", "valid")                                  # read both sets of collocations
fpa    = My.ouvre(ARGS[1], "r") ; tinea = readlines(fpa) ; close(fpa)
fpb    = My.ouvre(ARGS333, "r") ; tineb = readlines(fpb) ; close(fpb)
tinuma = length(tinea)
tinumb = length(tineb)

ARGS444 = replace(ARGS222, "calib", "valid")                                  # also read the other current component
fpa    = My.ouvre(ARGS222, "r") ; tinec = readlines(fpa) ; close(fpa)
fpb    = My.ouvre(ARGS444, "r") ; tined = readlines(fpb) ; close(fpb)
tinumc = length(tinec)
tinumd = length(tined)
if tinuma != tinumc || tinumb != tinumd
  print("\nERROR: number of lines in $(ARGS[1]) and $(ARGS222) are $tinuma != $tinumc\n\n")
  print("\nERROR: number of lines in $(ARGS333) and $(ARGS444) are $tinumb != $tinumd\n\n")
  exit(-1)
end

refa = Array(Float64, tinuma)                                                 # and calculate a pair of reference variables
refb = Array(Float64, tinumb)                                                 # (either from observations or from analyses)
for a = 1:tinuma
  vala = float(split(tinea[a]))
  valb = float(split(tinec[a]))
  refa[a] = (vala[OCUR]^2.0 + valb[OCUR]^2.0)^0.5
# refa[a] = (vala[TOTN]^2.0 + valb[TOTN]^2.0)^0.5
end
for a = 1:tinumb
  vala = float(split(tineb[a]))
  valb = float(split(tined[a]))
  refb[a] = (vala[OCUR]^2.0 + valb[OCUR]^2.0)^0.5
# refb[a] = (vala[TOTN]^2.0 + valb[TOTN]^2.0)^0.5
end

locmas = [MISS for z = RANGE]                                                 # allocate two sets of local metrics
lodmas = [MISS for z = RANGE]
locrcm = Array(Float64, length(RANGE), RCMR, RCMC)
lodrcm = Array(Float64, length(RANGE), RCMR, RCMC)
linuma = tinuma < CUTOFF ? tinuma : CUTOFF
linumb = tinumb < CUTOFF ? tinumb : CUTOFF
dista  = Array(Float64, tinuma)
distb  = Array(Float64, tinumb)
maska  = Array(Bool,    tinuma)
maskb  = Array(Bool,    tinumb)
tmpaa  = zeros(linuma) ; tmpab = zeros(linuma) ; tmpac = zeros(linuma)
tmpba  = zeros(linumb) ; tmpbb = zeros(linumb) ; tmpbc = zeros(linumb)

for (z, ranz) in enumerate(RANGE)                                             # then loop through the target parameter
  for a = 1:tinuma   dista[a] = abs(ranz - refa[a])  end                      # and isolate the nearest CUTOFF set of obs
  for a = 1:tinumb   distb[a] = abs(ranz - refb[a])  end
  lima = sort(dista)[linuma]
  limb = sort(distb)[linumb]
  b = 1 ; for a = 1:tinuma  if dista[a] <= lima && b <= linuma  maska[a] = true ; b += 1  else  maska[a] = false  end  end
  b = 1 ; for a = 1:tinumb  if distb[a] <= limb && b <= linumb  maskb[a] = true ; b += 1  else  maskb[a] = false  end  end
  linea = tinea[maska] ; lrefa = refa[maska]
  lineb = tineb[maskb] ; lrefb = refb[maskb]

  for a = 1:linuma                                                            # compute cal/val metrics for the first set
    vals = float(split(linea[a]))
    tmpaa[a] = lrefa[a]
    tmpab[a] = vals[OCUR]
    tmpac[a] = vals[TOTN]
  end
  (mass, rcm) = regress(tmpaa, tmpab, tmpac, 1.0)
  locmas[z] = mass
  locrcm[z,:,:] = rcm

  for a = 1:linumb                                                            # compute cal/val metrics for the second set
    vals = float(split(lineb[a]))
    tmpba[a] = lrefb[a]
    tmpbb[a] = vals[OCUR]
    tmpbc[a] = vals[TOTN]
  end
  (mass, rcm) = regress(tmpba, tmpbb, tmpbc, 1.0)
  lodmas[z] = mass
  lodrcm[z,:,:] = rcm
end

locfit = Array(Function, RCMR, RCMC)                                          # define an array of simple functions that fit
lodfit = Array(Function, RCMR, RCMC)                                          # the regression coefficient matrix values and
function funfit(a, b, c, min = -9e99, max = 9e99)                             # save their coefficients to a file
  x -> (y = a + b * exp(c * x) ; y >= min ? (y < max ? y : max) : min)
end

fpb = My.ouvre(ARGS[1] * regionlab * ".cali.regr.funfit", "w")
mask = trues(length(locrcm[:,BOLS,SNRU])) ; tmpa = locmas[mask]
for a = 1:RCMR, b = 1:RCMC
  tmpx = locrcm[1,a,b] * ones(locrcm[mask,a,b]) - DELT ; tmpy = tmpx + DELT * 2 ; tmpz = trues(locrcm[mask,a,b])
  if (tmpx .< locrcm[mask,a,b] .< tmpy) == tmpz  tmpc = tmpf = tmpg = locrcm[1,a,b] ; tmpd = tmpe = 0 ; else
  tmpb = locrcm[mask,a,b] ; tmpc, tmpd, tmpe = My.integralexp(tmpa, tmpb) ; tmpf = minimum(tmpb) ; tmpg = maximum(tmpb) ; end
  tmpp = @sprintf("%2d %2d %8.4f %8.4f %8.4f %8.4f %8.4f\n", a, b, tmpc, tmpd, tmpe, tmpf, tmpg) ; write(fpb, tmpp)
  locfit[a,b] = funfit(tmpc, tmpd, tmpe, tmpf, tmpg)
end
mask = trues(length(lodrcm[:,BOLS,SNRU])) ; tmpa = lodmas[mask]
for a = 1:RCMR, b = 1:RCMC
# if lodrcm[mask,a,b] == lodrcm[1,a,b] * ones(lodrcm[mask,a,b])  tmpc = tmpf = tmpg = lodrcm[1,a,b] ; tmpd = tmpe = 0 ; else
  tmpx = lodrcm[1,a,b] * ones(lodrcm[mask,a,b]) - DELT ; tmpy = tmpx + DELT * 2 ; tmpz = trues(lodrcm[mask,a,b])
  if (tmpx .< lodrcm[mask,a,b] .< tmpy) == tmpz  tmpc = tmpf = tmpg = lodrcm[1,a,b] ; tmpd = tmpe = 0 ; else
  tmpb = lodrcm[mask,a,b] ; tmpc, tmpd, tmpe = My.integralexp(tmpa, tmpb) ; tmpf = minimum(tmpb) ; tmpg = maximum(tmpb) ; end
  tmpp = @sprintf("%2d %2d %8.4f %8.4f %8.4f %8.4f %8.4f\n", a, b, tmpc, tmpd, tmpe, tmpf, tmpg) ; write(fpb, tmpp)
  lodfit[a,b] = funfit(tmpc, tmpd, tmpe, tmpf, tmpg)
end
close(fpb)

tars = collect(RANGE)                                                         # plot the binned sums 
tarv = [locmas locmas locmas]
alpv = zeros(length(tars), RCMR) ; alpf = zeros(length(tars), RCMR)
betv = zeros(length(tars), RCMR) ; betf = zeros(length(tars), RCMR)
sigv = zeros(length(tars), RCMR) ; sigf = zeros(length(tars), RCMR)
snrv = zeros(length(tars), RCMR) ; snrf = zeros(length(tars), RCMR)
for a = 1:RCMR
  alpv[:,a] = locrcm[:,a,ALPH]
  betv[:,a] = locrcm[:,a,BETA]
  sigv[:,a] = locrcm[:,a,SIGU]
  snrv[:,a] = locrcm[:,a,SNRU]
  for (b, ref) in enumerate(tars)
    alpf[b,a] = locfit[a,ALPH](ref)
    betf[b,a] = locfit[a,BETA](ref)
    sigf[b,a] = locfit[a,SIGU](ref)
    snrf[b,a] = locfit[a,SNRU](ref)
  end
end

zonmed = ""
contains(ARGS[1], ".ucur") && (zonmed = " for U")
contains(ARGS[1], ".vcur") && (zonmed = " for V")

ppp = Winston.Table(2,2) ; setattr(ppp, "cellpadding", -0.5)                  # and then create the plots
for z = 1:4
  z == 1 && (varname = "a) Additive Bias (ms^{-1})" * zonmed ; bound = tarv ; grid = alpv ; grie = alpf ; tpos = (1,1))
  z == 2 && (varname = "b) Multiplicative Bias"     * zonmed ; bound = tarv ; grid = betv ; grie = betf ; tpos = (1,2))
  z == 3 && (varname = "c) RMSE (ms^{-1})"          * zonmed ; bound = tarv ; grid = sigv ; grie = sigf ; tpos = (2,1))
  z == 4 && (varname = "d) Correlation"             * zonmed ; bound = tarv ; grid = snrv ; grie = snrf ; tpos = (2,2))

  z == 1 && (xmin = 0.05 ; xmax = 1.15 ; ymin = -0.1  ; ymax = 1.5 )
# z == 1 && (xmin = 0.05 ; xmax = 1.15 ; ymin = -0.1  ; ymax = 0.1 )          # and locate the plot limits
  z == 2 && (xmin = 0.05 ; xmax = 1.15 ; ymin =  0.0  ; ymax = 3.0 )
  z == 3 && (xmin = 0.05 ; xmax = 1.15 ; ymin = -0.05 ; ymax = 0.25)
  z == 4 && (xmin = 0.05 ; xmax = 1.15 ; ymin =  0.4  ; ymax = 1.1 )

  ump = Array(Any, 6)
  cols = [  "red", "green",  "blue",    "red",  "green",    "blue"]
  kynd = ["solid", "solid", "solid", "dashed", "dashed",  "dashed"]
  dirs = [  "OLS",  "CONF",  "ROLS", "OLSfit", "CONfit", "ROLSfit"]
# xmin = 0.0 ; xmax = 1.4 ; ymin = -0.5 ; ymax = 1.0

  tmp = Winston.FramedPlot(title="$varname", xrange = (xmin,xmax), yrange = (ymin,ymax))
  ppp[tpos...] = Winston.add(tmp)

  for a = 1:3
    ump[a]   = Winston.Curve(bound[:,a], grid[:,a], "color", parse(Winston.Colorant, cols[a]))
               style(ump[a], kind = kynd[a])
               setattr(ump[a], label = dirs[a])
               Winston.add(ppp[tpos...], ump[a])
    ump[a+3] = Winston.Curve(      tars, grie[:,a], "color", parse(Winston.Colorant, cols[a+3]))
               style(ump[a+3], kind = kynd[a+3])
               setattr(ump[a+3], label = dirs[a+3])
               Winston.add(ppp[tpos...], ump[a+3])
  end
  if z == 2 && argc == 2
    tmp = Winston.Legend(.45, .82, Any[ump[1], ump[2], ump[3]])
          Winston.add(ppp[tpos...], tmp)
#   tmp = Winston.Legend(.70, .82, Any[ump[5], ump[6], ump[7], ump[8]])
#         Winston.add(ppp[tpos...], tmp)
  end
  if z == 2 && argc == 3
    tmp = Winston.PlotLabel(.45, .82, regionlab[2:end], "texthalign", "left", "size", 3.0)
          Winston.add(ppp[tpos...], tmp)
  end
end

xyzzy = ARGS[1] * regionlab * ".cali.regr.png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
exit(0)

#=
 = Closure function that defines a set of optimal (in a least squares sense) polynomial coefficients.
 = Given xclose and yclose (e.g., reference current speed and corresponding triple collocation metric
 = values, respectively, for a set of collocation subsets), this closure returns optimal coefficients.

xclose = Array(Float64, 0)                                                    # (data arrays are in global scope)
yclose = Array(Float64, 0)

function sqerror(coef::Array{Float64,1})
  err = 0.0
  for i in 1:length(xclose)
    res  = coef[1] * xclose[i]^2 + coef[2] * xclose[i] + coef[3]
    err += (yclose[i] - res)^2
  end
  return err
end

function triple(curr::Array{Float64,3})
#=mask = masquextreme(curr[1,:,2], SDTRIM) &                                  # get the parametric center of mass
         masquextreme(curr[1,:,1], SDTRIM) &                                  # after trimming extreme values first
         masquextreme(curr[2,:,1], SDTRIM)
  mass     =     mean(curr[2,mask,2])
  sampsitu =          curr[1,mask,2]
  samprefa =          curr[1,mask,1]
  samprefb =          curr[2,mask,1]
# @show length(mask) length(mask[mask])
# mass     =     mean(curr[2,:,2])
# sampsitu =          curr[1,:,2] - mass
# samprefa =          curr[1,:,1] - mass
# samprefb =          curr[2,:,1] - mass

  avg1 = mean(sampsitu)                                                       # and use a robust calculation of covariance
  avg2 = mean(samprefa)                                                       # (two-pass here, but more algorithms are at
  avg3 = mean(samprefb)                                                       # en.wikipedia.org/wiki/Algorithms_for_calculating_variance)
  cv11 = mean((sampsitu - avg1) .* (sampsitu - avg1))
  cv12 = mean((sampsitu - avg1) .* (samprefa - avg2))
  cv13 = mean((sampsitu - avg1) .* (samprefb - avg3))
  cv22 = mean((samprefa - avg2) .* (samprefa - avg2))
  cv23 = mean((samprefa - avg2) .* (samprefb - avg3))
  cv33 = mean((samprefb - avg3) .* (samprefb - avg3))
=#
  mass = mean(curr[2,:,2])
# curr[1,:,2] -= mass
# curr[1,:,1] -= mass
# curr[2,:,1] -= mass
  temp = [curr[1,:,2]' curr[1,:,1]' curr[2,:,1]']
  remp = rcopy(R"DetMCD($temp, alpha = $MCDTRIM)")
  mask = falses(length(temp[:,1])) ; for a in remp[:Hsubsets]  mask[a] = true  end
# mass = mean(curr[2,mask,2])

  avg1 = remp[:center][1]
  avg2 = remp[:center][2]
  avg3 = remp[:center][3]
  cv11 = remp[:cov][1,1]
  cv12 = remp[:cov][1,2]
  cv13 = remp[:cov][1,3]
  cv22 = remp[:cov][2,2]
  cv23 = remp[:cov][2,3]
  cv33 = remp[:cov][3,3]

  bet2 = cv23 / cv13
  bet3 = cv23 / cv12
  alp2 = avg2 - bet2 * avg1 #+ mass * (1.0 - bet2)
  alp3 = avg3 - bet3 * avg1 #+ mass * (1.0 - bet3)

  tmpval = cv11 - cv12 * cv13 / cv23 ; sig1 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv22 - cv12 * cv23 / cv13 ; sig2 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv33 - cv13 * cv23 / cv12 ; sig3 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv12 * cv13 / cv11 / cv23 ; cor1 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv12 * cv23 / cv22 / cv13 ; cor2 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv13 * cv23 / cv33 / cv12 ; cor3 = tmpval > 0 ? sqrt(tmpval) : 0.0

  return(mass, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3)    # then return all statistics
end


statis = [MISS for a = 1:4, b = 1:MEMS, c = 1:PARS]                           # allocate a set of global cal/val arrays
glomas = [MISS for a = 1:4]
gloalp = [MISS for a = 1:4]
globet = [MISS for a = 1:4]
glosig = [MISS for a = 1:4]
glocor = [MISS for a = 1:4]
curga  = zeros(2, tinuma, 2)
curgb  = zeros(2, tinumb, 2)

for a = 1:tinuma                                                              # report cal/val metrics for the first set
  vals = float(split(tinea[a]))
  curga[1,a,:] = [vals[TOTB] vals[OCUR]]
  curga[2,a,:] = [vals[TOTA] refa[a]   ]
end
a = 1 ; (mass, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(curga)
statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        glomas[a] = mass
statis[a,MEMO,ALPH] =  0.0 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; gloalp[a] = 0.5 * (alp2 + alp3)
statis[a,MEMO,BETA] =  1.0 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; globet[a] = 0.5 * (bet2 + bet3)
statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; glosig[a] = 0.5 * (sig2 + sig3)
statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; glocor[a] = 0.5 * (cor2 + cor3)
#a = 1 ; (glomas[a], gloalp[a], globet[a], glosig[a], glocor[a]) = triple(curga)

@printf("\nnumb = %15d for %s\n", tinuma, ARGS[1])
@printf("cala = %15.8f mean(vals[TOTB]) = %15.8f\n",      gloalp[a],  mean(curga[1,:,1]))
@printf("calb = %15.8f mean(vals[TOTA]) = %15.8f\n",      globet[a],  mean(curga[2,:,1]))
@printf("mean = %15.8f mean(vals[OCUR]) = %15.8f\n", mean(glomas[a]), mean(curga[1,:,2]))
@printf("%33s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", gloalp[a], globet[a], glosig[a], glocor[a])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMO,ALPH], statis[a,MEMO,BETA], statis[a,MEMO,SIGM], statis[a,MEMO,CORR])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMB,ALPH], statis[a,MEMB,BETA], statis[a,MEMB,SIGM], statis[a,MEMB,CORR])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMA,ALPH], statis[a,MEMA,BETA], statis[a,MEMA,SIGM], statis[a,MEMA,CORR])

for a = 1:tinumb                                                              # report cal/val metrics for the second set
  vals = float(split(tineb[a]))
  curgb[1,a,:] = [vals[TOTB] vals[OCUR]]
  curgb[2,a,:] = [vals[TOTA] refb[a]   ]
end
a = 2 ; (mass, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(curgb)
statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        glomas[a] = mass
statis[a,MEMO,ALPH] =  0.0 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; gloalp[a] = 0.5 * (alp2 + alp3)
statis[a,MEMO,BETA] =  1.0 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; globet[a] = 0.5 * (bet2 + bet3)
statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; glosig[a] = 0.5 * (sig2 + sig3)
statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; glocor[a] = 0.5 * (cor2 + cor3)
#a = 2 ; (glomas[a], gloalp[a], globet[a], glosig[a], glocor[a]) = triple(curgb)

@printf("\nnumb = %15d for %s\n", tinumb, ARGS333)
@printf("cala = %15.8f mean(vals[TOTB]) = %15.8f\n",      gloalp[a],  mean(curgb[1,:,1]))
@printf("calb = %15.8f mean(vals[TOTA]) = %15.8f\n",      globet[a],  mean(curgb[2,:,1]))
@printf("mean = %15.8f mean(vals[OCUR]) = %15.8f\n", mean(glomas[a]), mean(curgb[1,:,2]))
@printf("%33s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", gloalp[a], globet[a], glosig[a], glocor[a])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMO,ALPH], statis[a,MEMO,BETA], statis[a,MEMO,SIGM], statis[a,MEMO,CORR])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMB,ALPH], statis[a,MEMB,BETA], statis[a,MEMB,SIGM], statis[a,MEMB,CORR])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMA,ALPH], statis[a,MEMA,BETA], statis[a,MEMA,SIGM], statis[a,MEMA,CORR])

fpb = My.ouvre(ARGS[1] * regionlab * ".cali.regr", "w")
form = @sprintf("  mean param   CSPD is %6.2f\n", mean(glomas[1]))
write(fpb, form)
form = @sprintf("  mean param   CSPD is %6.2f\n", mean(glomas[2]))
write(fpb, form)
form = @sprintf("%77s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", ARGS[1] * "." * ARGS[2], gloalp[1], globet[1], glosig[1], glocor[1])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", ARGS333 * "." * ARGS[2], gloalp[2], globet[2], glosig[2], glocor[2])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "obs", statis[1,MEMO,ALPH], statis[1,MEMO,BETA], statis[1,MEMO,SIGM], statis[1,MEMO,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "obs", statis[2,MEMO,ALPH], statis[2,MEMO,BETA], statis[2,MEMO,SIGM], statis[2,MEMO,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "bef", statis[1,MEMB,ALPH], statis[1,MEMB,BETA], statis[1,MEMB,SIGM], statis[1,MEMB,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "bef", statis[2,MEMB,ALPH], statis[2,MEMB,BETA], statis[2,MEMB,SIGM], statis[2,MEMB,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "aft", statis[1,MEMA,ALPH], statis[1,MEMA,BETA], statis[1,MEMA,SIGM], statis[1,MEMA,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "aft", statis[2,MEMA,ALPH], statis[2,MEMA,BETA], statis[2,MEMA,SIGM], statis[2,MEMA,CORR])
write(fpb, form)
close(fpb)


if POLY                                                                       # either solve polynomial coefficients or
  xclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(xclose, locmas[a])  end
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, localp[a])  end ; localpint = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, locbet[a])  end ; locbetint = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, locsig[a])  end ; locsigint = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, loccor[a])  end ; loccorint = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)
  xclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(xclose, lodmas[a])  end
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, lodalp[a])  end ; lodalpint = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, lodbet[a])  end ; lodbetint = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, lodsig[a])  end ; lodsigint = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, lodcor[a])  end ; lodcorint = optimize(sqerror, [0.0, 0.0, 0.0], iterations = 10000)
else
  localpint = interpolate((locmas, ), localp, Gridded(Linear()))              # create interpolation functions for the local
  locbetint = interpolate((locmas, ), locbet, Gridded(Linear()))              # dependence of alpha and beta on the _in situ_
  locsigint = interpolate((locmas, ), locsig, Gridded(Linear()))              # reference variable (note also that one can't
  loccorint = interpolate((locmas, ), loccor, Gridded(Linear()))              # usually assume to have in situ in practice)
  lodalpint = interpolate((lodmas, ), lodalp, Gridded(Linear()))
  lodbetint = interpolate((lodmas, ), lodbet, Gridded(Linear()))
  lodsigint = interpolate((lodmas, ), lodsig, Gridded(Linear()))              # (tinea = ucal  tineb = uval)
  lodcorint = interpolate((lodmas, ), lodcor, Gridded(Linear()))              # (tinec = vcal  tined = vval)
end


# localp[z] = locrcm[z,BOLS,ALPH]
# locbet[z] = locrcm[z,BOLS,BETA]
# locsig[z] = locrcm[z,BOLS,SIGU]
# loccor[z] = locrcm[z,BOLS,SNRU]

# @printf("\nnumb = %15.0f for subset of %s\n", linuma, ARGS[1])
# @printf("cala = %15.8f\n",                localp[z])
# @printf("calb = %15.8f\n",                locbet[z])
# @printf("mean = %15.8f target = %5.2f\n", locmas[z], ranz)
# @printf("%33s %8s %8s %8s %8s\n", " ", "localp", "locbet", "locsig", "loccor")
# @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", localp[z], locbet[z], locsig[z], loccor[z])

# lodalp[z] = lodrcm[z,BOLS,ALPH]
# lodbet[z] = lodrcm[z,BOLS,BETA]
# lodsig[z] = lodrcm[z,BOLS,SIGU]
# lodcor[z] = lodrcm[z,BOLS,SNRU]

# @printf("\nnumb = %15.0f for subset of %s\n", linumb, ARGS333)
# @printf("cala = %15.8f\n",                lodalp[z])
# @printf("calb = %15.8f\n",                lodbet[z])
# @printf("mean = %15.8f target = %5.2f\n", lodmas[z], ranz)
# @printf("%33s %8s %8s %8s %8s\n", " ", "lodalp", "lodbet", "lodsig", "lodcor")
# @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", lodalp[z], lodbet[z], lodsig[z], lodcor[z])


for a = 1:tinuma                                                              # recalibrate using the calibration parameters from
  vala = float(split(tinea[a]))                                               # the other set; first get a refbef from gloalp/bet
  valb = float(split(tinec[a]))                                               # (using that of u or v for both u and v) and then
# refbef = (  vala[TOTB]                          ^2.0 +   valb[TOTB]                          ^2.0)^0.5
# refaft = (  vala[TOTA]                          ^2.0 +   valb[TOTA]                          ^2.0)^0.5
# refbef = (((vala[TOTB] - gloalp[2]) / globet[2])^2.0 + ((valb[TOTB] - gloalp[2]) / globet[2])^2.0)^0.5
# refaft = (((vala[TOTA] - gloalp[2]) / globet[2])^2.0 + ((valb[TOTA] - gloalp[2]) / globet[2])^2.0)^0.5
  refnow = (  vala[TOTN]                          ^2.0 +   valb[TOTN]                          ^2.0)^0.5
# alpbef = lodalpint.minimum[1] * refbef^2 + lodalpint.minimum[2] * refbef + lodalpint.minimum[3]
# alpaft = lodalpint.minimum[1] * refaft^2 + lodalpint.minimum[2] * refaft + lodalpint.minimum[3]
# betbef = lodbetint.minimum[1] * refbef^2 + lodbetint.minimum[2] * refbef + lodbetint.minimum[3]
# betaft = lodbetint.minimum[1] * refaft^2 + lodbetint.minimum[2] * refaft + lodbetint.minimum[3]
##alpnow = lodalpint.minimum[1] * refnow^2 + lodalpint.minimum[2] * refnow + lodalpint.minimum[3]
##betnow = lodbetint.minimum[1] * refnow^2 + lodbetint.minimum[2] * refnow + lodbetint.minimum[3]
  alpnow = lodalpint(refnow)
  betnow = lodbetint(refnow)
# vala[TOTB] = (vala[TOTB] - alpbef) / betbef                                 # get alp/betbef from refbef, and similarly for aft
# vala[TOTA] = (vala[TOTA] - alpaft) / betaft
  vala[TOTB] = (vala[TOTB] - alpnow) / betnow
  vala[TOTA] = (vala[TOTA] - alpnow) / betnow
  curga[1,a,:] = [vala[TOTB] vala[OCUR]]
  curga[2,a,:] = [vala[TOTA] refa[a]   ]
end
a = 3 ; (mass, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(curga)
statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        glomas[a] = mass
statis[a,MEMO,ALPH] =  0.0 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; gloalp[a] = 0.5 * (alp2 + alp3)
statis[a,MEMO,BETA] =  1.0 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; globet[a] = 0.5 * (bet2 + bet3)
statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; glosig[a] = 0.5 * (sig2 + sig3)
statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; glocor[a] = 0.5 * (cor2 + cor3)
#a = 3 ; (glomas[a], gloalp[a], globet[a], glosig[a], glocor[a]) = triple(curga)

tmpstr = "after recalibration only (using alpha and beta from the other collocations)"

@printf("\nnumb = %15.0f for %s\n", tinuma, ARGS[1])
@printf("cala = %15.8f %s\n",      gloalp[a],  tmpstr)
@printf("calb = %15.8f %s\n",      globet[a],  tmpstr)
@printf("mean = %15.8f %s\n", mean(glomas[a]), tmpstr)
@printf("%33s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", gloalp[a], globet[a], glosig[a], glocor[a])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMO,ALPH], statis[a,MEMO,BETA], statis[a,MEMO,SIGM], statis[a,MEMO,CORR])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMB,ALPH], statis[a,MEMB,BETA], statis[a,MEMB,SIGM], statis[a,MEMB,CORR])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMA,ALPH], statis[a,MEMA,BETA], statis[a,MEMA,SIGM], statis[a,MEMA,CORR])

for a = 1:tinumb                                                              # recalibrate using the calibration parameters
  vala = float(split(tineb[a]))                                               # the other set; first get a refbef from gloalp/bet
  valb = float(split(tined[a]))                                               # (using that of u or v for both u and v) and then
# refbef = (  vala[TOTB]                          ^2.0 +   valb[TOTB]                          ^2.0)^0.5
# refaft = (  vala[TOTA]                          ^2.0 +   valb[TOTA]                          ^2.0)^0.5
# refbef = (((vala[TOTB] - gloalp[1]) / globet[1])^2.0 + ((valb[TOTB] - gloalp[1]) / globet[1])^2.0)^0.5
# refaft = (((vala[TOTA] - gloalp[1]) / globet[1])^2.0 + ((valb[TOTA] - gloalp[1]) / globet[1])^2.0)^0.5
  refnow = (  vala[TOTN]                          ^2.0 +   valb[TOTN]                          ^2.0)^0.5
# alpbef = localpint.minimum[1] * refbef^2 + localpint.minimum[2] * refbef + localpint.minimum[3]
# alpaft = localpint.minimum[1] * refaft^2 + localpint.minimum[2] * refaft + localpint.minimum[3]
# betbef = locbetint.minimum[1] * refbef^2 + locbetint.minimum[2] * refbef + locbetint.minimum[3]
# betaft = locbetint.minimum[1] * refaft^2 + locbetint.minimum[2] * refaft + locbetint.minimum[3]
##alpnow = localpint.minimum[1] * refnow^2 + localpint.minimum[2] * refnow + localpint.minimum[3]
##betnow = locbetint.minimum[1] * refnow^2 + locbetint.minimum[2] * refnow + locbetint.minimum[3]
  alpnow = localpint(refnow)
  betnow = locbetint(refnow)
# vala[TOTB] = (vala[TOTB] - alpbef) / betbef                                 # get alp/betbef from refbef, and similarly for aft
# vala[TOTA] = (vala[TOTA] - alpaft) / betaft
  vala[TOTB] = (vala[TOTB] - alpnow) / betnow
  vala[TOTA] = (vala[TOTA] - alpnow) / betnow
  curgb[1,a,:] = [vala[TOTB] vala[OCUR]]
  curgb[2,a,:] = [vala[TOTA] refb[a]   ]
end
a = 4 ; (mass, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(curgb)
statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        glomas[a] = mass
statis[a,MEMO,ALPH] =  0.0 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; gloalp[a] = 0.5 * (alp2 + alp3)
statis[a,MEMO,BETA] =  1.0 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; globet[a] = 0.5 * (bet2 + bet3)
statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; glosig[a] = 0.5 * (sig2 + sig3)
statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; glocor[a] = 0.5 * (cor2 + cor3)
#a = 4 ; (glomas[a], gloalp[a], globet[a], glosig[a], glocor[a]) = triple(curgb)

@printf("\nnumb = %15.0f for %s\n", tinumb, ARGS333)
@printf("cala = %15.8f %s\n",      gloalp[a],  tmpstr)
@printf("calb = %15.8f %s\n",      globet[a],  tmpstr)
@printf("mean = %15.8f %s\n", mean(glomas[a]), tmpstr)
@printf("%33s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", gloalp[a], globet[a], glosig[a], glocor[a])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMO,ALPH], statis[a,MEMO,BETA], statis[a,MEMO,SIGM], statis[a,MEMO,CORR])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMB,ALPH], statis[a,MEMB,BETA], statis[a,MEMB,SIGM], statis[a,MEMB,CORR])
@printf("%33s %8.4f %8.4f %8.4f %8.4f\n", " ", statis[a,MEMA,ALPH], statis[a,MEMA,BETA], statis[a,MEMA,SIGM], statis[a,MEMA,CORR])

fpb = My.ouvre(ARGS[1]  * regionlab * ".cali.regr", "a")
form = @sprintf("  mean param   CSPD is %6.2f %s\n", mean(glomas[3]), tmpstr)
write(fpb, form)
form = @sprintf("  mean param   CSPD is %6.2f %s\n", mean(glomas[4]), tmpstr)
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", ARGS[1] * "." * ARGS[2], gloalp[3], globet[3], glosig[3], glocor[3])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", ARGS222 * "." * ARGS[2], gloalp[4], globet[4], glosig[4], glocor[4])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "obs", statis[3,MEMO,ALPH], statis[3,MEMO,BETA], statis[3,MEMO,SIGM], statis[3,MEMO,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "obs", statis[4,MEMO,ALPH], statis[4,MEMO,BETA], statis[4,MEMO,SIGM], statis[4,MEMO,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "bef", statis[3,MEMB,ALPH], statis[3,MEMB,BETA], statis[3,MEMB,SIGM], statis[3,MEMB,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "bef", statis[4,MEMB,ALPH], statis[4,MEMB,BETA], statis[4,MEMB,SIGM], statis[4,MEMB,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "aft", statis[3,MEMA,ALPH], statis[3,MEMA,BETA], statis[3,MEMA,SIGM], statis[3,MEMA,CORR])
write(fpb, form)
form = @sprintf("%77s %8.4f %8.4f %8.4f %8.4f\n", "aft", statis[4,MEMA,ALPH], statis[4,MEMA,BETA], statis[4,MEMA,SIGM], statis[4,MEMA,CORR])
write(fpb, form)
close(fpb)

#localp = [MISS for b = RANGE] ; lodalp = [MISS for b = RANGE]
#locbet = [MISS for b = RANGE] ; lodbet = [MISS for b = RANGE]
#locsig = [MISS for b = RANGE] ; lodsig = [MISS for b = RANGE]
#loccor = [MISS for b = RANGE] ; lodcor = [MISS for b = RANGE]
mask = trues(length(loccor)) # (abs(locbet) .< 10) # trues(length(loccor)) # (0 .< loccor .< 1) # (abs(localp) .< 10000) & (abs(locbet) .< 500) & (0 .< loccor .< 1)

@show My.integralexp(locmas[msk], localp[msk])..., minimum(localp[msk]), maximum(localp[msk])
tmpa = locmas[msk] ; tmpb = locrcm[msk,BOLS,ALPH]
@show My.integralexp(tmpa, tmpb)..., minimum(tmpb), maximum(tmpb)
exit(0)

localpint =                                                                funb(My.integralexp(locmas[msk], localp[msk])..., minimum(localp[msk]), maximum(localp[msk]))
localpfrm = @sprintf("localpint a,b,c,min,max %8.4f %8.4f %8.4f %8.4f %8.4f\n", My.integralexp(locmas[msk], localp[msk])..., minimum(localp[msk]), maximum(localp[msk]))
locbetint =                                                                funb(My.integralexp(locmas[msk], locbet[msk])..., minimum(locbet[msk]), maximum(locbet[msk]))
locbetfrm = @sprintf("locbetint a,b,c,min,max %8.4f %8.4f %8.4f %8.4f %8.4f\n", My.integralexp(locmas[msk], locbet[msk])..., minimum(locbet[msk]), maximum(locbet[msk]))
locsigint =                                                                funb(My.integralexp(locmas[msk], locsig[msk])..., minimum(locsig[msk]), maximum(locsig[msk]))
locsigfrm = @sprintf("locsigint a,b,c,min,max %8.4f %8.4f %8.4f %8.4f %8.4f\n", My.integralexp(locmas[msk], locsig[msk])..., minimum(locsig[msk]), maximum(locsig[msk]))
loccorint =                                                                funb(My.integralexp(locmas[msk], loccor[msk])..., minimum(loccor[msk]), maximum(loccor[msk]))
loccorfrm = @sprintf("loccorint a,b,c,min,max %8.4f %8.4f %8.4f %8.4f %8.4f\n", My.integralexp(locmas[msk], loccor[msk])..., minimum(loccor[msk]), maximum(loccor[msk]))
      msk = trues(length(lodcor)) # (abs(lodbet) .< 10) # trues(length(lodcor)) # (0 .< lodcor .< 1) # (abs(lodalp) .< 10000) & (abs(lodbet) .< 500) & (0 .< lodcor .< 1)
lodalpint =                                                                funb(My.integralexp(lodmas[msk], lodalp[msk])..., minimum(lodalp[msk]), maximum(lodalp[msk]))
lodalpfrm = @sprintf("lodalpint a,b,c,min,max %8.4f %8.4f %8.4f %8.4f %8.4f\n", My.integralexp(lodmas[msk], lodalp[msk])..., minimum(lodalp[msk]), maximum(lodalp[msk]))
lodbetint =                                                                funb(My.integralexp(lodmas[msk], lodbet[msk])..., minimum(lodbet[msk]), maximum(lodbet[msk]))
lodbetfrm = @sprintf("lodbetint a,b,c,min,max %8.4f %8.4f %8.4f %8.4f %8.4f\n", My.integralexp(lodmas[msk], lodbet[msk])..., minimum(lodbet[msk]), maximum(lodbet[msk]))
lodsigint =                                                                funb(My.integralexp(lodmas[msk], lodsig[msk])..., minimum(lodsig[msk]), maximum(lodsig[msk]))
lodsigfrm = @sprintf("lodsigint a,b,c,min,max %8.4f %8.4f %8.4f %8.4f %8.4f\n", My.integralexp(lodmas[msk], lodsig[msk])..., minimum(lodsig[msk]), maximum(lodsig[msk]))
lodcorint =                                                                funb(My.integralexp(lodmas[msk], lodcor[msk])..., minimum(lodcor[msk]), maximum(lodcor[msk]))
lodcorfrm = @sprintf("lodcorint a,b,c,min,max %8.4f %8.4f %8.4f %8.4f %8.4f\n", My.integralexp(lodmas[msk], lodcor[msk])..., minimum(lodcor[msk]), maximum(lodcor[msk]))
fpb = My.ouvre(ARGS[1] * regionlab * ".cali.regr.funb", "w")
write(fpb, localpfrm)
write(fpb, locbetfrm)
write(fpb, locsigfrm)
write(fpb, loccorfrm)
write(fpb, lodalpfrm)
write(fpb, lodbetfrm)
write(fpb, lodsigfrm)
write(fpb, lodcorfrm)
close(fpb)
=#
