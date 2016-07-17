#=
 = Perform a series of paired analysis validation (and performance) estimates, each of which
 = involves a recalibration using the OLS slope and intercept from the opposing collocations,
 = followed by a global calibration and performance estimate.  Fixed-size subsets of the
 = available collocations are selected based on closeness to a target value of the variable
 = of interest (the actual mean value is obtained as a simple average and if two averages
 = are close, these might be combined).  The target value is then varied over a reasonable
 = range and calibration is obtained for each subset.  Smooth variations in calibration (as
 = a function of the variable of interest and expressed by a polynomial whose coefficients
 = are determined by least squares) is applied in the final global performance assessment
 = - RD June 2016.
 =#

using My, Morph, RCall                                  # , Winston
const ODAT             = 1                              # identify indecies of the input data:
const OLAT             = 2                              # date/lat/lon on the collocation grid
const OLON             = 3
const OCUR             = 4                              # drifter current component as well as
const TOTB             = 11                             # the two total current estimates from
const TOTA             = 12                             # before and after

const MISS             = -9999.0                        # generic missing value
const SDTRIM           = 6.0                            # standard deviation trimming limit
const MORPH            = false                          # perform Gaussian anamorphosis

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb\n\n")
  exit(1)
end
contains(ARGS[1], ".ucur") && (ARGS222 = replace(ARGS[1], "ucur", "vcur"))
contains(ARGS[1], ".vcur") && (ARGS222 = replace(ARGS[1], "vcur", "ucur"))
#=
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
=#
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

const RANGE            = 0.1:0.1:1.1                    # target sampling range for current speed
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

refa = Array(Float64, tinuma)                                                 # and calculate current speed from obs
refb = Array(Float64, tinumb)
for a = 1:tinuma
  vala = float(split(tinea[a]))
  valb = float(split(tinec[a]))
  refa[a] = (vala[OCUR]^2.0 + valb[OCUR]^2.0)^0.5
end
for a = 1:tinumb
  vala = float(split(tineb[a]))
  valb = float(split(tined[a]))
  refb[a] = (vala[OCUR]^2.0 + valb[OCUR]^2.0)^0.5
end

allmas = [MISS for a = 1:4, b = RANGE]                                        # allocate input/output arrays, including
allalp = [MISS for a = 1:4, b = RANGE]                                        # distance to obs, the resulting mean params,
allbet = [MISS for a = 1:4, b = RANGE]                                        # and the triple collocation cal/val estimates
allsig = [MISS for a = 1:4, b = RANGE]
allcor = [MISS for a = 1:4, b = RANGE]
linuma = tinuma < CUTOFF ? tinuma : CUTOFF
linumb = tinumb < CUTOFF ? tinumb : CUTOFF
dista  = Array(Float64, tinuma)
distb  = Array(Float64, tinumb)
maska  = Array(Bool,    tinuma)
maskb  = Array(Bool,    tinumb)
cura   = zeros(2, linuma, 2)
curb   = zeros(2, linumb, 2)

for (z, ranz) in enumerate(RANGE)                                             # loop through the target parameter and
  for a = 1:tinuma   dista[a] = abs(ranz - refa[a])  end                      # isolate the nearest CUTOFF set of obs
  for a = 1:tinumb   distb[a] = abs(ranz - refb[a])  end
  lima = sort(dista)[linuma]
  limb = sort(distb)[linumb]
  b = 1 ; for a = 1:tinuma  if dista[a] <= lima && b <= linuma  maska[a] = true ; b += 1  else  maska[a] = false  end  end
  b = 1 ; for a = 1:tinumb  if distb[a] <= limb && b <= linumb  maskb[a] = true ; b += 1  else  maskb[a] = false  end  end
  linea = tinea[maska]
  lineb = tineb[maskb]

  if MORPH                                                                    # define for both collocation sets a mapping
    tmpa = Array(Float64, linuma)                                             # (Gaussian anamorphosis) based on the in situ
    tmpb = Array(Float64, linumb)                                             # observations (which are taken to be unbiased)
    for a = 1:linuma  vals = float(split(linea[a])) ; tmpa[a] = vals[OCUR]  end     # and apply each map to the opposite set
    for a = 1:linumb  vals = float(split(lineb[a])) ; tmpb[a] = vals[OCUR]  end
    shfa = minimum(tmpa) < 0 ? -2.0 * minimum(tmpa) : 0.0
    shfb = minimum(tmpb) < 0 ? -2.0 * minimum(tmpb) : 0.0
    rawa = sort(tmpa) + shfa
    rawb = sort(tmpb) + shfb

    R"library(RGeostats)"
    R"dba  = db.create(aa = c($rawa))"
    R"fia  = anam.fit(dba, name='aa', type='emp')"
    R"ana <- anam.z2y(fia, dba, 'aa')"
    cdfa = unique(rcopy(R"db.extract(ana, names = c('Gaussian.aa'))"))
    R"rm(dba, fia, ana)"
    rawa = unique(rawa - shfa)
    if length(rawa) != length(cdfa)
      print("\nERROR: length(rawa) $(length(rawa)) != $(length(cdfa)) length(cdfa)\n\n")
      exit(1)
    end

    R"dbb  = db.create(bb = c($rawb))"
    R"fib  = anam.fit(dbb, name='bb', type='emp')"
    R"anb <- anam.z2y(fib, dbb, 'bb')"
    cdfb = unique(rcopy(R"db.extract(anb, names = c('Gaussian.bb'))"))
    R"rm(dbb, fib, anb)"
    rawb = unique(rawb - shfb)
    if length(rawb) != length(cdfb)
      print("\nERROR: length(rawb) $(length(rawb)) != $(length(cdfb)) length(cdfb)\n\n")
      exit(1)
    end                                                                       # finally extend and equate the tails

    rawa[  1] > rawb[  1] && (rawa[  1] = rawb[  1]) ; rawb[  1] > rawa[  1] && (rawb[  1] = rawa[  1])
    rawa[end] < rawb[end] && (rawa[end] = rawb[end]) ; rawb[end] < rawa[end] && (rawb[end] = rawa[end])
    cdfa[  1] > cdfb[  1] && (cdfa[  1] = cdfb[  1]) ; cdfb[  1] > cdfa[  1] && (cdfb[  1] = cdfa[  1])
    cdfa[end] < cdfb[end] && (cdfa[end] = cdfb[end]) ; cdfb[end] < cdfa[end] && (cdfb[end] = cdfa[end])
    cdfa[  1] >        -9 && (cdfa[  1] =        -9) ; cdfb[  1] >        -9 && (cdfb[  1] =        -9)
    cdfa[end] <         9 && (cdfa[end] =         9) ; cdfb[end] <         9 && (cdfb[end] =         9)

#   tmp      = Winston.FramedPlot(title="Empirical anamorphosis", xlabel="Gaussian Current", ylabel="Actual Current")
#   ppp      = Winston.add(tmp)
#   tmp      = Winston.Curve(cdfa, rawa, "color", parse(Winston.Colorant,   "black"))
#              Winston.add(ppp, tmp)
#   tmp      = Winston.Curve(cdfb, rawb, "color", parse(Winston.Colorant,   "red"))
#              Winston.add(ppp, tmp)
#   xyzzy = "gaussian_anamorphosis_test.png"
#   print("writing $xyzzy\n")
#   Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
  end

  for a = 1:linuma                                                            # report cal/val parameters for the first set
    vals = float(split(linea[a]))                                             # (in original units regardless of morphing)
    cura[1,a,:] = [vals[TOTB] vals[OCUR]]
    cura[2,a,:] = [vals[TOTA] vals[OCUR]]
  end

  a = 1 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(cura)
  @printf("\nnumb = %15d for %s\n", linuma, ARGS[1])
  @printf("cala = %15.8f\n", allalp[a])
  @printf("calb = %15.8f\n", allbet[a])
  @printf("mean = %15.8f\n", mean(allmas[a]))
  @printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
  @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", allalp[a], allbet[a], allsig[a], allcor[a])

  for a = 1:linumb                                                            # report cal/val parameters for the second set
    vals = float(split(lineb[a]))                                             # (in original units regardless of morphing)
    curb[1,a,:] = [vals[TOTB] vals[OCUR]]
    curb[2,a,:] = [vals[TOTA] vals[OCUR]]
  end

  a = 2 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(curb)
  @printf("\nnumb = %15d for %s\n", linumb, ARGS333)
  @printf("cala = %15.8f\n", allalp[a])
  @printf("calb = %15.8f\n", allbet[a])
  @printf("mean = %15.8f\n", mean(allmas[a]))
  @printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
  @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", allalp[a], allbet[a], allsig[a], allcor[a])

  if MORPH                                                                    # save linear mapping metrics
    rawavga11 = mean(cura[1,:,1]) ; rawavga21 = mean(cura[2,:,1]) ; rawavga22 = mean(cura[2,:,2])
    rawvara11 =  var(cura[1,:,1]) ; rawvara21 =  var(cura[2,:,1]) ; rawvara22 =  var(cura[2,:,2])
    rawavgb11 = mean(curb[1,:,1]) ; rawavgb21 = mean(curb[2,:,1]) ; rawavgb22 = mean(curb[2,:,2])
    rawvarb11 =  var(curb[1,:,1]) ; rawvarb21 =  var(curb[2,:,1]) ; rawvarb22 =  var(curb[2,:,2])
  end

  MORPH && (fpb = My.ouvre(ARGS[1] * ".cali.pair.morph", "w"))
  MORPH || (fpb = My.ouvre(ARGS[1] * ".cali.pair",       "w"))
  form = @sprintf("  mean param   CSPD is %6.2f\n", mean(allmas[1]))
  write(fpb, form)
  form = @sprintf("  mean param   CSPD is %6.2f\n", mean(allmas[2]))
  write(fpb, form)
  form = @sprintf("%77s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
  write(fpb, form)
  form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS[1], allalp[1], allbet[1], allsig[1], allcor[1])
  write(fpb, form)
  form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS333, allalp[2], allbet[2], allsig[2], allcor[2])
  write(fpb, form)
  close(fpb)

  for a = 1:linuma                                                            # compute cal/val parameters for the first set
    vals = float(split(linea[a]))                                             # (in Gaussian units if morphing)
    MORPH && (vals[TOTB] = rawtonorm(cdfb, rawb, vals[TOTB]))
    MORPH && (vals[TOTA] = rawtonorm(cdfb, rawb, vals[TOTA]))
    MORPH && (vals[OCUR] = rawtonorm(cdfb, rawb, vals[OCUR]))
    cura[1,a,:] = [vals[TOTB] vals[OCUR]]
    cura[2,a,:] = [vals[TOTA] vals[OCUR]]
  end
  a = 1 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(cura)

  for a = 1:linumb                                                            # compute cal/val parameters for the second set
    vals = float(split(lineb[a]))                                             # (in Gaussian units if morphing)
    MORPH && (vals[TOTB] = rawtonorm(cdfa, rawa, vals[TOTB]))
    MORPH && (vals[TOTA] = rawtonorm(cdfa, rawa, vals[TOTA]))
    MORPH && (vals[OCUR] = rawtonorm(cdfa, rawa, vals[OCUR]))
    curb[1,a,:] = [vals[TOTB] vals[OCUR]]
    curb[2,a,:] = [vals[TOTA] vals[OCUR]]
  end
  a = 2 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(curb)

  if MORPH                                                                    # save linear mapping metrics
    gauavga11 = mean(cura[1,:,1]) ; gauavga21 = mean(cura[2,:,1]) ; gauavga22 = mean(cura[2,:,2])
    gauvara11 =  var(cura[1,:,1]) ; gauvara21 =  var(cura[2,:,1]) ; gauvara22 =  var(cura[2,:,2])
    gauavgb11 = mean(curb[1,:,1]) ; gauavgb21 = mean(curb[2,:,1]) ; gauavgb22 = mean(curb[2,:,2])
    gauvarb11 =  var(curb[1,:,1]) ; gauvarb21 =  var(curb[2,:,1]) ; gauvarb22 =  var(curb[2,:,2])
  end

  for a = 1:linuma                                                            # apply recalibration either in Gaussian or original units
    vals = float(split(linea[a]))                                             # using calibration parameters (and anamorphosis) from the
    MORPH && (vals[TOTB] = rawtonorm(cdfb, rawb, vals[TOTB]))                 # other set and report new cal/val parameters in original
    MORPH && (vals[TOTA] = rawtonorm(cdfb, rawb, vals[TOTA]))                 # units
    MORPH && (vals[OCUR] = rawtonorm(cdfb, rawb, vals[OCUR]))
    vals[TOTB] = (vals[TOTB] - allalp[2]) / allbet[2]
    vals[TOTA] = (vals[TOTA] - allalp[2]) / allbet[2]
    MORPH && (vals[TOTB] = normtoraw(cdfb, rawb, vals[TOTB]))
    MORPH && (vals[TOTA] = normtoraw(cdfb, rawb, vals[TOTA]))
    MORPH && (vals[OCUR] = normtoraw(cdfb, rawb, vals[OCUR]))
    cura[1,a,:] = [vals[TOTB] vals[OCUR]]
    cura[2,a,:] = [vals[TOTA] vals[OCUR]]
  end
  a = 3 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(cura)

            tmpstr  = "after recalibration only (using alpha and beta from the other collocations)"
  MORPH && (tmpstr  = "after applying Gaussian anamorphosis, recalibration, and inverse Gaussian anamorphosis")

  @printf("\nnumb = %15.0f for %s\n", linuma, ARGS[1])
  @printf("cala = %15.8f %s\n", allalp[a], tmpstr)
  @printf("calb = %15.8f %s\n", allbet[a], tmpstr)
  @printf("mean = %15.8f %s\n", mean(allmas[a]), tmpstr)
  @printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
  @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", allalp[a], allbet[a], allsig[a], allcor[a])

  for a = 1:linumb                                                            # apply recalibration either in Gaussian or original units
    vals = float(split(lineb[a]))                                             # using calibration parameters (and anamorphosis) from the
    MORPH && (vals[TOTB] = rawtonorm(cdfa, rawa, vals[TOTB]))                 # other set and report new cal/val parameters in original
    MORPH && (vals[TOTA] = rawtonorm(cdfa, rawa, vals[TOTA]))                 # units
    MORPH && (vals[OCUR] = rawtonorm(cdfa, rawa, vals[OCUR]))
    vals[TOTB] = (vals[TOTB] - allalp[1]) / allbet[1]
    vals[TOTA] = (vals[TOTA] - allalp[1]) / allbet[1]
    MORPH && (vals[TOTB] = normtoraw(cdfa, rawa, vals[TOTB]))
    MORPH && (vals[TOTA] = normtoraw(cdfa, rawa, vals[TOTA]))
    MORPH && (vals[OCUR] = normtoraw(cdfa, rawa, vals[OCUR]))
    curb[1,a,:] = [vals[TOTB] vals[OCUR]]
    curb[2,a,:] = [vals[TOTA] vals[OCUR]]
  end
  a = 4 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(curb)

  @printf("\nnumb = %15.0f for %s\n", linumb, ARGS333)
  @printf("cala = %15.8f %s\n", allalp[a], tmpstr)
  @printf("calb = %15.8f %s\n", allbet[a], tmpstr)
  @printf("mean = %15.8f %s\n", mean(allmas[a]), tmpstr)
  @printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
  @printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", allalp[a], allbet[a], allsig[a], allcor[a])

  if MORPH                                                                    # get average linear mapping metrics (as
    linrawa = (abs(rawavga11 - normtoraw(cdfb, rawb, gauavga11)) / rawvara11 +  # mean difference normalized by variance)
               abs(rawavga21 - normtoraw(cdfb, rawb, gauavga21)) / rawvara21 +
               abs(rawavga22 - normtoraw(cdfb, rawb, gauavga22)) / rawvara22) * 100 / 3
    linrawb = (abs(rawavgb11 - normtoraw(cdfb, rawb, gauavgb11)) / rawvarb11 +
               abs(rawavgb21 - normtoraw(cdfb, rawb, gauavgb21)) / rawvarb21 +
               abs(rawavgb22 - normtoraw(cdfb, rawb, gauavgb22)) / rawvarb22) * 100 / 3
    lingaua = (abs(gauavga11 - rawtonorm(cdfb, rawb, rawavga11)) / gauvara11 +
               abs(gauavga21 - rawtonorm(cdfb, rawb, rawavga21)) / gauvara21 +
               abs(gauavga22 - rawtonorm(cdfb, rawb, rawavga22)) / gauvara22) * 100 / 3
    lingaub = (abs(gauavgb11 - rawtonorm(cdfb, rawb, rawavgb11)) / gauvarb11 +
               abs(gauavgb21 - rawtonorm(cdfb, rawb, rawavgb21)) / gauvarb21 +
               abs(gauavgb22 - rawtonorm(cdfb, rawb, rawavgb22)) / gauvarb22) * 100 / 3
    linstr  = @sprintf("linrawa = %5.1f linrawb = %5.1f lingaua = %5.1f lingaub = %5.1f", linrawa, linrawb, lingaua, lingaub)
    println("\n$linstr\n")
    tmpstr *= " $linstr"
  end

  MORPH && (fpb = My.ouvre(ARGS[1] * ".cali.pair.morph", "a"))
  MORPH || (fpb = My.ouvre(ARGS[1] * ".cali.pair",       "a"))
  form = @sprintf("  mean param   CSPD is %6.2f %s\n", mean(allmas[3]), tmpstr)
  write(fpb, form)
  form = @sprintf("  mean param   CSPD is %6.2f %s\n", mean(allmas[4]), tmpstr)
  write(fpb, form)
  form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS[1], allalp[3], allbet[3], allsig[3], allcor[3])
  write(fpb, form)
  form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS333, allalp[4], allbet[4], allsig[4], allcor[4])
  write(fpb, form)
  close(fpb)
end

exit(0)
