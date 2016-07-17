#=
 = Perform a series of paired analysis validation (and performance) estimates by triple collocation.
 = Each of the series involves a recalibration using the OLS slope and intercept from the opposing
 = collocations.  Fixed-size subsets of the available collocations are selected based on closeness
 = to a target value of the variable of interest (the actual mean value is obtained as a simple
 = average and if two averages are close, these might be combined).  The target value is then
 = varied over a reasonable range and calibration is obtained for each subset.  These variations
 = in calibration (as a function of the observed or analyzed variable of interest) are then used
 = to assess global performance - RD June 2016.
 =#

using My, Interpolations, Optim, Winston
const ODAT             = 1                              # identify indecies of the input data:
const OLAT             = 2                              # date/lat/lon on the collocation grid
const OLON             = 3
const OCUR             = 4                              # drifter current component as well as
const TOTB             = 14                             # the three total current estimates from
const TOTN             = 15                             # before, now, and after
const TOTA             = 16

const POLY             = true                           # polynomial interpolation (not linear)
const MISS             = -9999.0                        # generic missing value
const SDTRIM           = 6.0                            # standard deviation trimming limit

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb\n\n")
  exit(1)
end
contains(ARGS[1], ".ucur") && (ARGS222 = replace(ARGS[1], "ucur", "vcur"))
contains(ARGS[1], ".vcur") && (ARGS222 = replace(ARGS[1], "vcur", "ucur"))

#=
 = Closure function that defines a set of optimal (in a least squares sense) polynomial coefficients.
 = Given xclose and yclose (e.g., reference current speed and corresponding triple collocation metric
 = values, respectively, for a set of collocation subsets), this closure returns optimal coefficients.
 =#

xclose = Array(Float64, 0)                                                    # (data arrays are in global scope)
yclose = Array(Float64, 0)

function sqerror(coef::Array{Float64,1})
  err = 0.0
  for i in 1:length(xclose)
    res  = coef[1] * xclose[i]^4 + coef[2] * xclose[i]^3 + coef[3] * xclose[i]^2 + coef[4] * xclose[i] + coef[5]
    err += (yclose[i] - res)^2
  end
  return err
end

#=
 = Function returning triple collocation cal/val measures for a group of analyses, following McColl
 = et al. (2014).  Inputs are an array of collocated values and stats are returned for a collocation
 = set, where it is assumed that extrapolation from before and after is done using the same analysis,
 = so no consideration of relative effective resolution is necessary (cf. Vogelzang et al. 2011)
 =#

function triple(curr::Array{Float64,3})
  allalp = MISS
  allbet = MISS
  allsig = MISS
  allcor = MISS
  allmas = MISS

  mask = masquextreme(curr[1,   :,2], SDTRIM) &                               # get the parametric center of mass
         masquextreme(curr[1,   :,1], SDTRIM) &                               # after trimming extreme values first
         masquextreme(curr[2,   :,1], SDTRIM)
  sampsitu =          curr[1,mask,2]
  samprefa =          curr[1,mask,1]
  samprefb =          curr[2,mask,1]
  allmas   =     mean(curr[2,mask,2])

  avg1 = mean(sampsitu)                                                       # and use a robust calculation of covariance
  avg2 = mean(samprefa)                                                       # (two-pass here, but more algorithms are at
  avg3 = mean(samprefb)                                                       # en.wikipedia.org/wiki/Algorithms_for_calculating_variance)
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

  allalp = 0.5 * (alp2 + alp3)
  allbet = 0.5 * (bet2 + bet3)
  allsig = 0.5 * (sig2 + sig3)
  allcor = 0.5 * (cor2 + cor3)

  return(allmas, allalp, allbet, allsig, allcor)                              # then return the average stats
end

#=
 = main program
 =#

const RANGE            = 0.1:0.01:1.1                   # target sampling range for current speed
const CUTOFF           = 200                            # number of collocations in a subset

const ALPH             = 1                              # error model x = ALPH + BETA * truth + error
const BETA             = 2                              # error model x = ALPH + BETA * truth + error
const SIGM             = 3                              # triple coll RMSE
const CORR             = 4                              # triple coll correlation coefficient
const MSPD             = 5                              # center-of-mass parameter
const PARAMS           = 5                              # number of triple collocation parameters

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
# refa[a] = (vala[OCUR]^2.0 + valb[OCUR]^2.0)^0.5
  refa[a] = (vala[TOTN]^2.0 + valb[TOTN]^2.0)^0.5
end
for a = 1:tinumb
  vala = float(split(tineb[a]))
  valb = float(split(tined[a]))
# refb[a] = (vala[OCUR]^2.0 + valb[OCUR]^2.0)^0.5
  refb[a] = (vala[TOTN]^2.0 + valb[TOTN]^2.0)^0.5
end

glomas = [MISS for a = 1:4]                                                   # allocate a set of global cal/val arrays
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
a = 1 ; (glomas[a], gloalp[a], globet[a], glosig[a], glocor[a]) = triple(curga)

@printf("\nnumb = %15d for %s\n", tinuma, ARGS[1])
@printf("cala = %15.8f\n", gloalp[a])
@printf("calb = %15.8f\n", globet[a])
@printf("mean = %15.8f\n", glomas[a])
@printf("%33s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", gloalp[a], globet[a], glosig[a], glocor[a])

for a = 1:tinumb                                                              # report cal/val metrics for the second set
  vals = float(split(tineb[a]))
  curgb[1,a,:] = [vals[TOTB] vals[OCUR]]
  curgb[2,a,:] = [vals[TOTA] refb[a]   ]
end
a = 2 ; (glomas[a], gloalp[a], globet[a], glosig[a], glocor[a]) = triple(curgb)

@printf("\nnumb = %15d for %s\n", tinumb, ARGS333)
@printf("cala = %15.8f\n", gloalp[a])
@printf("calb = %15.8f\n", globet[a])
@printf("mean = %15.8f\n", glomas[a])
@printf("%33s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", gloalp[a], globet[a], glosig[a], glocor[a])

fpb = My.ouvre(ARGS[1] * ".cali.ploc", "w")
form = @sprintf("  mean param   CSPD is %6.2f\n", glomas[1])
write(fpb, form)
form = @sprintf("  mean param   CSPD is %6.2f\n", glomas[2])
write(fpb, form)
form = @sprintf("%77s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS[1], gloalp[1], globet[1], glosig[1], glocor[1])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS333, gloalp[2], globet[2], glosig[2], glocor[2])
write(fpb, form)
close(fpb)

locmas = [MISS for b = RANGE] ; lodmas = [MISS for b = RANGE]                 # allocate two sets of local metrics
localp = [MISS for b = RANGE] ; lodalp = [MISS for b = RANGE]
locbet = [MISS for b = RANGE] ; lodbet = [MISS for b = RANGE]
locsig = [MISS for b = RANGE] ; lodsig = [MISS for b = RANGE]
loccor = [MISS for b = RANGE] ; lodcor = [MISS for b = RANGE]
linuma = tinuma < CUTOFF ? tinuma : CUTOFF
linumb = tinumb < CUTOFF ? tinumb : CUTOFF
dista  = Array(Float64, tinuma)
distb  = Array(Float64, tinumb)
maska  = Array(Bool,    tinuma)
maskb  = Array(Bool,    tinumb)
curla  = zeros(2, linuma, 2)
curlb  = zeros(2, linumb, 2)

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
    curla[1,a,:] = [vals[TOTB] vals[OCUR]]
    curla[2,a,:] = [vals[TOTA] lrefa[a]  ]
  end
  (locmas[z], localp[z], locbet[z], locsig[z], loccor[z]) = triple(curla)

# @printf("\nnumb = %15.0f for subset of %s\n", linuma, ARGS[1])
# @printf("cala = %15.8f\n",                localp[z])
# @printf("calb = %15.8f\n",                locbet[z])
# @printf("mean = %15.8f target = %5.2f\n", locmas[z], ranz)
# @printf("%33s %8s %8s %8s %8s\n", " ", "localp", "locbet", "locsig", "loccor")
# @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", localp[z], locbet[z], locsig[z], loccor[z])

  for a = 1:linumb                                                            # compute cal/val metrics for the second set
    vals = float(split(lineb[a]))
    curlb[1,a,:] = [vals[TOTB] vals[OCUR]]
    curlb[2,a,:] = [vals[TOTA] lrefb[a]  ]
  end
  (lodmas[z], lodalp[z], lodbet[z], lodsig[z], lodcor[z]) = triple(curlb)

# @printf("\nnumb = %15.0f for subset of %s\n", linumb, ARGS333)
# @printf("cala = %15.8f\n",                lodalp[z])
# @printf("calb = %15.8f\n",                lodbet[z])
# @printf("mean = %15.8f target = %5.2f\n", lodmas[z], ranz)
# @printf("%33s %8s %8s %8s %8s\n", " ", "lodalp", "lodbet", "lodsig", "lodcor")
# @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", lodalp[z], lodbet[z], lodsig[z], lodcor[z])
end

if POLY                                                                       # either solve polynomial coefficients or
  xclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(xclose, locmas[a])  end
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, localp[a])  end ; localpint = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, locbet[a])  end ; locbetint = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, locsig[a])  end ; locsigint = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, loccor[a])  end ; loccorint = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)
  xclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(xclose, lodmas[a])  end
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, lodalp[a])  end ; lodalpint = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, lodbet[a])  end ; lodbetint = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, lodsig[a])  end ; lodsigint = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)
  yclose = Array(Float64, 0) ; for (a, rana) in enumerate(RANGE)  push!(yclose, lodcor[a])  end ; lodcorint = optimize(sqerror, [0.0, 0.0, 0.0, 0.0, 0.0], iterations = 10000)
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

for a = 1:tinuma                                                              # recalibrate using the calibration parameters from
  vala = float(split(tinea[a]))                                               # the other set; first get a refbef from gloalp/bet
  valb = float(split(tinec[a]))                                               # (using that of u or v for both u and v) and then
  refbef = (  vala[TOTB]                          ^2.0 +   valb[TOTB]                          ^2.0)^0.5
  refaft = (  vala[TOTA]                          ^2.0 +   valb[TOTA]                          ^2.0)^0.5
# refbef = (((vala[TOTB] - gloalp[2]) / globet[2])^2.0 + ((valb[TOTB] - gloalp[2]) / globet[2])^2.0)^0.5
# refaft = (((vala[TOTA] - gloalp[2]) / globet[2])^2.0 + ((valb[TOTA] - gloalp[2]) / globet[2])^2.0)^0.5
  for b = 1:100
    alpbef = lodalpint.minimum[1] * refbef^4 + lodalpint.minimum[2] * refbef^3 + lodalpint.minimum[3] * refbef^2 + lodalpint.minimum[4] * refbef + lodalpint.minimum[5]
    alpaft = lodalpint.minimum[1] * refaft^4 + lodalpint.minimum[2] * refaft^3 + lodalpint.minimum[3] * refaft^2 + lodalpint.minimum[4] * refaft + lodalpint.minimum[5]
    betbef = lodbetint.minimum[1] * refbef^4 + lodbetint.minimum[2] * refbef^3 + lodbetint.minimum[3] * refbef^2 + lodbetint.minimum[4] * refbef + lodbetint.minimum[5]
    betaft = lodbetint.minimum[1] * refaft^4 + lodbetint.minimum[2] * refaft^3 + lodbetint.minimum[3] * refaft^2 + lodbetint.minimum[4] * refaft + lodbetint.minimum[5]
    refbef = (((vala[TOTB] - alpbef) / betbef)^2.0 + ((valb[TOTB] - alpbef) / betbef)^2.0)^0.5
    refaft = (((vala[TOTA] - alpaft) / betaft)^2.0 + ((valb[TOTA] - alpaft) / betaft)^2.0)^0.5
  end
  if POLY
    alpbef = lodalpint.minimum[1] * refbef^4 + lodalpint.minimum[2] * refbef^3 + lodalpint.minimum[3] * refbef^2 + lodalpint.minimum[4] * refbef + lodalpint.minimum[5]
    alpaft = lodalpint.minimum[1] * refaft^4 + lodalpint.minimum[2] * refaft^3 + lodalpint.minimum[3] * refaft^2 + lodalpint.minimum[4] * refaft + lodalpint.minimum[5]
    betbef = lodbetint.minimum[1] * refbef^4 + lodbetint.minimum[2] * refbef^3 + lodbetint.minimum[3] * refbef^2 + lodbetint.minimum[4] * refbef + lodbetint.minimum[5]
    betaft = lodbetint.minimum[1] * refaft^4 + lodbetint.minimum[2] * refaft^3 + lodbetint.minimum[3] * refaft^2 + lodbetint.minimum[4] * refaft + lodbetint.minimum[5]
  else
    alpbef = lodalpint[refbef] ; alpaft = lodalpint[refaft]
    betbef = lodbetint[refbef] ; betaft = lodbetint[refaft]
  end
  vala[TOTB] = (vala[TOTB] - alpbef) / betbef                                 # get alp/betbef from refbef, and similarly for aft
  vala[TOTA] = (vala[TOTA] - alpaft) / betaft
  curga[1,a,:] = [vala[TOTB] vala[OCUR]]
  curga[2,a,:] = [vala[TOTA] refa[a]   ]
end
a = 3 ; (glomas[a], gloalp[a], globet[a], glosig[a], glocor[a]) = triple(curga)

tmpstr = "after recalibration only (using alpha and beta from the other collocations)"

@printf("\nnumb = %15.0f for %s\n", tinuma, ARGS[1])
@printf("cala = %15.8f %s\n", gloalp[a], tmpstr)
@printf("calb = %15.8f %s\n", globet[a], tmpstr)
@printf("mean = %15.8f %s\n", glomas[a], tmpstr)
@printf("%33s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", gloalp[a], globet[a], glosig[a], glocor[a])

for a = 1:tinumb                                                              # recalibrate using the calibration parameters
  vala = float(split(tineb[a]))                                               # the other set; first get a refbef from gloalp/bet
  valb = float(split(tined[a]))                                               # (using that of u or v for both u and v) and then
  refbef = (  vala[TOTB]                          ^2.0 +   valb[TOTB]                          ^2.0)^0.5
  refaft = (  vala[TOTA]                          ^2.0 +   valb[TOTA]                          ^2.0)^0.5
# refbef = (((vala[TOTB] - gloalp[1]) / globet[1])^2.0 + ((valb[TOTB] - gloalp[1]) / globet[1])^2.0)^0.5
# refaft = (((vala[TOTA] - gloalp[1]) / globet[1])^2.0 + ((valb[TOTA] - gloalp[1]) / globet[1])^2.0)^0.5
  for b = 1:100
    alpbef = localpint.minimum[1] * refbef^4 + localpint.minimum[2] * refbef^3 + localpint.minimum[3] * refbef^2 + localpint.minimum[4] * refbef + localpint.minimum[5]
    alpaft = localpint.minimum[1] * refaft^4 + localpint.minimum[2] * refaft^3 + localpint.minimum[3] * refaft^2 + localpint.minimum[4] * refaft + localpint.minimum[5]
    betbef = locbetint.minimum[1] * refbef^4 + locbetint.minimum[2] * refbef^3 + locbetint.minimum[3] * refbef^2 + locbetint.minimum[4] * refbef + locbetint.minimum[5]
    betaft = locbetint.minimum[1] * refaft^4 + locbetint.minimum[2] * refaft^3 + locbetint.minimum[3] * refaft^2 + locbetint.minimum[4] * refaft + locbetint.minimum[5]
    refbef = (((vala[TOTB] - alpbef) / betbef)^2.0 + ((valb[TOTB] - alpbef) / betbef)^2.0)^0.5
    refaft = (((vala[TOTA] - alpaft) / betaft)^2.0 + ((valb[TOTA] - alpaft) / betaft)^2.0)^0.5
  end
  if POLY
    alpbef = localpint.minimum[1] * refbef^4 + localpint.minimum[2] * refbef^3 + localpint.minimum[3] * refbef^2 + localpint.minimum[4] * refbef + localpint.minimum[5]
    alpaft = localpint.minimum[1] * refaft^4 + localpint.minimum[2] * refaft^3 + localpint.minimum[3] * refaft^2 + localpint.minimum[4] * refaft + localpint.minimum[5]
    betbef = locbetint.minimum[1] * refbef^4 + locbetint.minimum[2] * refbef^3 + locbetint.minimum[3] * refbef^2 + locbetint.minimum[4] * refbef + locbetint.minimum[5]
    betaft = locbetint.minimum[1] * refaft^4 + locbetint.minimum[2] * refaft^3 + locbetint.minimum[3] * refaft^2 + locbetint.minimum[4] * refaft + locbetint.minimum[5]
  else
    alpbef = localpint[refbef] ; alpaft = localpint[refaft]
    betbef = locbetint[refbef] ; betaft = locbetint[refaft]
  end
  vala[TOTB] = (vala[TOTB] - alpbef) / betbef                                 # get alp/betbef from refbef, and similarly for aft
  vala[TOTA] = (vala[TOTA] - alpaft) / betaft
  curgb[1,a,:] = [vala[TOTB] vala[OCUR]]
  curgb[2,a,:] = [vala[TOTA] refb[a]   ]
end
a = 4 ; (glomas[a], gloalp[a], globet[a], glosig[a], glocor[a]) = triple(curgb)

@printf("\nnumb = %15.0f for %s\n", tinumb, ARGS333)
@printf("cala = %15.8f %s\n", gloalp[a], tmpstr)
@printf("calb = %15.8f %s\n", globet[a], tmpstr)
@printf("mean = %15.8f %s\n", glomas[a], tmpstr)
@printf("%33s %8s %8s %8s %8s\n", " ", "gloalp", "globet", "glosig", "glocor")
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", gloalp[a], globet[a], glosig[a], glocor[a])

fpb = My.ouvre(ARGS[1] * ".cali.ploc", "a")
form = @sprintf("  mean param   CSPD is %6.2f %s\n", mean(glomas[3]), tmpstr)
write(fpb, form)
form = @sprintf("  mean param   CSPD is %6.2f %s\n", mean(glomas[4]), tmpstr)
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS[1], gloalp[3], globet[3], glosig[3], glocor[3])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS333, gloalp[4], globet[4], glosig[4], glocor[4])
write(fpb, form)
close(fpb)

# exit(0)

tars = collect(RANGE)                                                         # plot the binned sums 
tarn = zeros(length(tars), 2)
alpn = zeros(length(tars), 2) ; alpo = zeros(length(tars), 2)
betn = zeros(length(tars), 2) ; beto = zeros(length(tars), 2)
sign = zeros(length(tars), 2) ; sigo = zeros(length(tars), 2)
corn = zeros(length(tars), 2) ; coro = zeros(length(tars), 2)

tarn[:,1] = locmas ; tarn[:,2] = lodmas
alpn[:,1] = localp ; alpn[:,2] = lodalp
betn[:,1] = locbet ; betn[:,2] = lodbet
sign[:,1] = locsig ; sign[:,2] = lodsig
corn[:,1] = loccor ; corn[:,2] = lodcor
for (a, ref) in enumerate(tars)
  alpo[a,1] = localpint.minimum[1] * ref^4 + localpint.minimum[2] * ref^3 + localpint.minimum[3] * ref^2 + localpint.minimum[4] * ref + localpint.minimum[5]
  alpo[a,2] = lodalpint.minimum[1] * ref^4 + lodalpint.minimum[2] * ref^3 + lodalpint.minimum[3] * ref^2 + lodalpint.minimum[4] * ref + lodalpint.minimum[5]
  beto[a,1] = locbetint.minimum[1] * ref^4 + locbetint.minimum[2] * ref^3 + locbetint.minimum[3] * ref^2 + locbetint.minimum[4] * ref + locbetint.minimum[5]
  beto[a,2] = lodbetint.minimum[1] * ref^4 + lodbetint.minimum[2] * ref^3 + lodbetint.minimum[3] * ref^2 + lodbetint.minimum[4] * ref + lodbetint.minimum[5]
  sigo[a,1] = locsigint.minimum[1] * ref^4 + locsigint.minimum[2] * ref^3 + locsigint.minimum[3] * ref^2 + locsigint.minimum[4] * ref + locsigint.minimum[5]
  sigo[a,2] = lodsigint.minimum[1] * ref^4 + lodsigint.minimum[2] * ref^3 + lodsigint.minimum[3] * ref^2 + lodsigint.minimum[4] * ref + lodsigint.minimum[5]
  coro[a,1] = loccorint.minimum[1] * ref^4 + loccorint.minimum[2] * ref^3 + loccorint.minimum[3] * ref^2 + loccorint.minimum[4] * ref + loccorint.minimum[5]
  coro[a,2] = lodcorint.minimum[1] * ref^4 + lodcorint.minimum[2] * ref^3 + lodcorint.minimum[3] * ref^2 + lodcorint.minimum[4] * ref + lodcorint.minimum[5]
end

ppp = Winston.Table(2,2) ; setattr(ppp, "cellpadding", -0.5)                  # and then create the plots
for z = 1:4
  z == 1 && (varname = "a) Bias (ms^{-1})" ; bound = tarn ; grid = alpn ; tpos = (1,1) ; grie = alpo)
  z == 2 && (varname = "b) Slope"          ; bound = tarn ; grid = betn ; tpos = (1,2) ; grie = beto)
  z == 3 && (varname = "c) RMSE (ms^{-1})" ; bound = tarn ; grid = sign ; tpos = (2,1) ; grie = sigo)
  z == 4 && (varname = "d) Correlation"    ; bound = tarn ; grid = corn ; tpos = (2,2) ; grie = coro)

  z == 1 && (xmin = 0.05 ; xmax = 1.15 ; ymin = -0.1 ; ymax = 0.05)           # and locate the plot limits
  z == 2 && (xmin = 0.05 ; xmax = 1.15 ; ymin =  0.4 ; ymax = 3.5)
  z == 3 && (xmin = 0.05 ; xmax = 1.15 ; ymin =  0.0 ; ymax = 0.2)
  z == 4 && (xmin = 0.05 ; xmax = 1.15 ; ymin =  0.8 ; ymax = 1.0)

  ump = Array(Any, 4)
  cols = [  "red",  "blue",    "red",   "blue"]
  kynd = ["solid", "solid", "dashed", "dashed"]
  dirs = ["Grp-A", "Grp-B",  "Est-A",  "Est-B"]
# xmin = 0.0 ; xmax = 1.4 ; ymin = -0.5 ; ymax = 1.0

  tmp = Winston.FramedPlot(title="$varname", xrange = (xmin,xmax), yrange = (ymin,ymax))
  ppp[tpos...] = Winston.add(tmp)

  for a = 1:2
    ump[a]   = Winston.Curve(bound[:,a], grid[:,a], "color", parse(Winston.Colorant, cols[a]))
               style(ump[a], kind = kynd[a])
               setattr(ump[a], label = dirs[a])
               Winston.add(ppp[tpos...], ump[a])
    ump[a+2] = Winston.Curve(      tars, grie[:,a], "color", parse(Winston.Colorant, cols[a+2]))
               style(ump[a+2], kind = kynd[a+2])
               setattr(ump[a+2], label = dirs[a+2])
               Winston.add(ppp[tpos...], ump[a+2])
  end
  if z == 2
    tmp = Winston.Legend(.45, .82, Any[ump[1], ump[2], ump[3], ump[4]])
          Winston.add(ppp[tpos...], tmp)
#   tmp = Winston.Legend(.70, .82, Any[ump[5], ump[6], ump[7], ump[8]])
#         Winston.add(ppp[tpos...], tmp)
  end
end

xyzzy = ARGS[1] * ".cali.ploc.png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
exit(0)
