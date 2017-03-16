#=
 = Loop through the available collocations and grid the corresponding forward and backward
 = extrapolations relative to the actual (uninterpolated) values.  Note that BEF refers to
 = a forecast extrapolation using analysis data from before the target value and AFT refers
 = to a revcast extrapolation using data from afterward.  Relative to the values at the
 = extrapolation target time (TOTN), both local (binwise) and global regressions of the
 = two separate extrapolations (from before and after) are also saved.  The global
 = regression is then used to adjust the data and the same files are stored.  This defines
 = an extrapolation bias correction.  Extrapolation variance (relative to the uninterpolated
 = values) is also calculated and stored, both for all collocations and as a function of
 = current speed - RD May, November 2016, January, March 2017.
 =#

using My
const ODAT             = 1                              # identify indecies of the input data:
const OLAT             = 2                              # date/lat/lon/obs on the collocation grid
const OLON             = 3
const OCUR             = 4
const MISS             = -9999.0                        # generic missing value
const CUTOFF           = 200

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_extra.ucur.got2000_obs.comb v2.0_global_025_deg_total_15m\n\n")
  exit(1)
end

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

step = 0.02 ; bound = collect(-4.0:step:4.0)
gridbefnow = zeros(length(bound), length(bound)) ; conbefnow = zeros(length(bound))
gridaftnow = zeros(length(bound), length(bound)) ; conaftnow = zeros(length(bound))
gridbefobs = zeros(length(bound), length(bound)) ; conbefobs = zeros(length(bound)) ; conobsbef = zeros(length(bound))
gridnowobs = zeros(length(bound), length(bound)) ; connowobs = zeros(length(bound)) ; conobsnow = zeros(length(bound))
gridaftobs = zeros(length(bound), length(bound)) ; conaftobs = zeros(length(bound)) ; conobsaft = zeros(length(bound))
gridbefaft = zeros(length(bound), length(bound)) ; conbefaft = zeros(length(bound)) ; conaftbef = zeros(length(bound))
regressbef =      Array(Float64, 0)
regressnow =      Array(Float64, 0)
regressaft =      Array(Float64, 0)
regressobs =      Array(Float64, 0)
binbef     = fill(Array(Float64, 0), length(bound))
binnow     = fill(Array(Float64, 0), length(bound))
binaft     = fill(Array(Float64, 0), length(bound))
cutbef     =      Array(Float64,     length(bound), CUTOFF)
cutnow     =      Array(Float64,     length(bound), CUTOFF)
cutaft     =      Array(Float64,     length(bound), CUTOFF)

if contains(ARGS[1], "wcur")                                                  # substitute speed for current component
  fila = replace(ARGS[1], "wcur", "ucur")
  filb = replace(ARGS[1], "wcur", "vcur")
  fpa = My.ouvre(fila, "r") ; tinea = readlines(fpa) ; close(fpa) ; tinuma = length(tinea)
  fpb = My.ouvre(filb, "r") ; tineb = readlines(fpb) ; close(fpb) ; tinumb = length(tineb)
  tinuma == tinumb || (print("\nERROR: tinuma $tinuma != $tinumb tinumb\n\n") ; exit(-1))
  for a = 1:tinuma
    vala = float(split(tinea[a]))
    valb = float(split(tineb[a]))
    out = tinea[a][1:30]
    for b = 4:19
      valc = (vala[b]^2 + valb[b]^2)^0.5
      out *= @sprintf(" %9.3f", valc)
    end
    tinea[a] = out
  end
else
  fpa = My.ouvre(ARGS[1], "r") ; tinea = readlines(fpa) ; close(fpa) ; tinuma = length(tinea)
end

for a = 1:tinuma                                                              # grid the collocations and save values
  vals = float(split(tinea[a]))                                               # for a regression versus TOTN and OCUR
  if vals[TOTB] > -333 && vals[TOTB] < 333 &&
     vals[TOTN] > -333 && vals[TOTN] < 333 &&
     vals[TOTA] > -333 && vals[TOTA] < 333 &&
     vals[OCUR] > -333 && vals[OCUR] < 333
    delbef, indbef = findmin(abs(bound - vals[TOTB])) ; bound[indbef] > vals[TOTB] && indbef > 1 && (indbef -= 1)
    delnow, indnow = findmin(abs(bound - vals[TOTN])) ; bound[indnow] > vals[TOTN] && indnow > 1 && (indnow -= 1)
    delaft, indaft = findmin(abs(bound - vals[TOTA])) ; bound[indaft] > vals[TOTA] && indaft > 1 && (indaft -= 1)
    delobs, indobs = findmin(abs(bound - vals[OCUR])) ; bound[indobs] > vals[OCUR] && indobs > 1 && (indobs -= 1)
    gridbefnow[indbef,indnow] += 1 ; conbefnow[indnow] += vals[TOTB]
    gridaftnow[indaft,indnow] += 1 ; conaftnow[indnow] += vals[TOTA]
    gridbefobs[indbef,indobs] += 1 ; conbefobs[indobs] += vals[TOTB] ; conobsbef[indbef] += vals[OCUR]
    gridnowobs[indnow,indobs] += 1 ; connowobs[indobs] += vals[TOTN] ; conobsnow[indnow] += vals[OCUR]
    gridaftobs[indaft,indobs] += 1 ; conaftobs[indobs] += vals[TOTA] ; conobsaft[indaft] += vals[OCUR]
    gridbefaft[indbef,indaft] += 1 ; conbefaft[indaft] += vals[TOTB] ; conaftbef[indbef] += vals[TOTA]
    push!(regressbef, vals[TOTB])
    push!(regressnow, vals[TOTN])
    push!(regressaft, vals[TOTA])
    push!(regressobs, vals[OCUR])
    tmp = copy(binbef[indnow]) ; push!(tmp, vals[TOTB]) ; binbef[indnow] = tmp
    tmp = copy(binnow[indnow]) ; push!(tmp, vals[TOTN]) ; binnow[indnow] = tmp
    tmp = copy(binaft[indnow]) ; push!(tmp, vals[TOTA]) ; binaft[indnow] = tmp
  end
end

valsb  = Array(Float64, tinuma)                                               # also grid for the extrapolation variance
valsn  = Array(Float64, tinuma)                                               # (versus TOTN) using all values in each bin
valsa  = Array(Float64, tinuma)                                               # (above) and for the nearest CUTOFF values
dista  = Array(Float64, tinuma)                                               # (here, but only save nearest CUTOFFs below)
linuma = tinuma < CUTOFF ? tinuma : CUTOFF

for a = 1:tinuma  vals = float(split(tinea[a])) ; valsb[a] = vals[TOTB] ; valsn[a] = vals[TOTN] ; valsa[a] = vals[TOTA]  end
for (z, ranz) in enumerate(bound)
  for a = 1:tinuma  dista[a] = abs(ranz - valsn[a])  end ; lima = sort(dista)[linuma] ; b = 1
  for a = 1:tinuma
    if dista[a] <= lima && b <= linuma
      cutbef[z,b] = valsb[a]
      cutnow[z,b] = valsn[a]
      cutaft[z,b] = valsa[a]
      b += 1
    end
  end
end

fname = ARGS[1] * "." * ARGS[2] * ".extra.dat"
fpa = My.ouvre(fname, "w")                                                    # and save the grids
for (a, vala) in enumerate(bound)
  for (b, valb) in enumerate(bound)
    @printf(fpa, "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n", gridbefnow[b,a], gridaftnow[b,a], gridbefobs[b,a], gridaftobs[b,a], gridbefaft[b,a], gridnowobs[b,a])
  end
end
close(fpa)

sumbefnow = sum(gridbefnow, 1)                                                # as well as the corresponding count and sum
sumaftnow = sum(gridaftnow, 1)                                                # of extrapolation values in each TOTN interval
sumbefobs = sum(gridbefobs, 1) ; sumobsbef = sum(gridbefobs, 2)
sumnowobs = sum(gridnowobs, 1) ; sumobsnow = sum(gridnowobs, 2)
sumaftobs = sum(gridaftobs, 1) ; sumobsaft = sum(gridaftobs, 2)
sumbefaft = sum(gridbefaft, 1) ; sumaftbef = sum(gridbefaft, 2)
fname = ARGS[1] * "." * ARGS[2] * ".extra.sum"
fpa = My.ouvre(fname, "w")
for (a, vala) in enumerate(bound)
  @printf(fpa, "%22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %18.11e %18.11e %22d %8.4f\n",
    sumbefnow[a], conbefnow[a], sumaftnow[a], conaftnow[a],
    sumbefobs[a], conbefobs[a], sumaftobs[a], conaftobs[a],
    sumobsbef[a], conobsbef[a], sumobsaft[a], conobsaft[a],
    sumbefaft[a], conbefaft[a], sumaftbef[a], conaftbef[a],
    sumnowobs[a], connowobs[a], sumobsnow[a], conobsnow[a],
    cov(vec(cutbef[a,:]), vec(cutnow[a,:])), cov(vec(cutaft[a,:]), vec(cutnow[a,:])), length(binbef[a]), bound[a])
end
close(fpa)

fname = ARGS[1] * "." * ARGS[2] * ".extra.reg"                                # and finally save the regression coefficients
fpa = My.ouvre(fname, "w")
(intbefnow, slobefnow) = linreg(regressnow, regressbef)
(intaftnow, sloaftnow) = linreg(regressnow, regressaft)
(intbefobs, slobefobs) = linreg(regressobs, regressbef) ; (intobsbef, sloobsbef) = linreg(regressbef, regressobs)
(intnowobs, slonowobs) = linreg(regressobs, regressnow) ; (intobsnow, sloobsnow) = linreg(regressnow, regressobs)
(intaftobs, sloaftobs) = linreg(regressobs, regressaft) ; (intobsaft, sloobsaft) = linreg(regressaft, regressobs)
(intbefaft, slobefaft) = linreg(regressaft, regressbef) ; (intaftbef, sloaftbef) = linreg(regressbef, regressaft)
            covbefnow  =    cov(regressnow, regressbef)
            covaftnow  =    cov(regressnow, regressaft)
            covbefobs  =    cov(regressobs, regressbef) ;             covobsbef  =    cov(regressbef, regressobs)
            covnowobs  =    cov(regressobs, regressnow) ;             covobsnow  =    cov(regressnow, regressobs)
            covaftobs  =    cov(regressobs, regressaft) ;             covobsaft  =    cov(regressaft, regressobs)
            covbefaft  =    cov(regressaft, regressbef) ;             covaftbef  =    cov(regressbef, regressaft)
@printf(fpa, "%33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f\n",
  intbefnow, slobefnow, intaftnow, sloaftnow,
  intbefobs, slobefobs, intaftobs, sloaftobs,
  intobsbef, sloobsbef, intobsaft, sloobsaft,
  intbefaft, slobefaft, intaftbef, sloaftbef,
  intnowobs, slonowobs, intobsnow, sloobsnow,
             covbefnow,            covaftnow,
             covbefobs,            covaftobs,
             covobsbef,            covobsaft,
             covbefaft,            covaftbef,
             covnowobs,            covobsnow)
close(fpa)

gridbefnow = zeros(length(bound), length(bound)) ; conbefnow = zeros(length(bound))
gridaftnow = zeros(length(bound), length(bound)) ; conaftnow = zeros(length(bound))
gridbefobs = zeros(length(bound), length(bound)) ; conbefobs = zeros(length(bound)) ; conobsbef = zeros(length(bound))
gridnowobs = zeros(length(bound), length(bound)) ; connowobs = zeros(length(bound)) ; conobsnow = zeros(length(bound))
gridaftobs = zeros(length(bound), length(bound)) ; conaftobs = zeros(length(bound)) ; conobsaft = zeros(length(bound))
gridbefaft = zeros(length(bound), length(bound)) ; conbefaft = zeros(length(bound)) ; conaftbef = zeros(length(bound))
regressbef =      Array(Float64, 0)
regressnow =      Array(Float64, 0)
regressaft =      Array(Float64, 0)                                           # reinitialize the variables for calibration
regressobs =      Array(Float64, 0)
binbef     = fill(Array(Float64, 0), length(bound))
binnow     = fill(Array(Float64, 0), length(bound))
binaft     = fill(Array(Float64, 0), length(bound))
cutbef     =      Array(Float64,     length(bound), CUTOFF)
cutnow     =      Array(Float64,     length(bound), CUTOFF)
cutaft     =      Array(Float64,     length(bound), CUTOFF)

for a = 1:tinuma                                                              # now calibrate the bef and aft collocations
  vals = float(split(tinea[a]))                                               # using global regressions above and regrid
  if vals[TOTB] > -333 && vals[TOTB] < 333 &&                                 # (but omit calibration of now versus obs)
     vals[TOTN] > -333 && vals[TOTN] < 333 &&
     vals[TOTA] > -333 && vals[TOTA] < 333 &&
     vals[OCUR] > -333 && vals[OCUR] < 333
    tmpb = (vals[TOTB] - intbefnow) /slobefnow
    tmpn =  vals[TOTN]
    tmpa = (vals[TOTA] - intaftnow) /sloaftnow
    tmpo =  vals[OCUR]
    delbef, indbef = findmin(abs(bound - tmpb)) ; bound[indbef] > tmpb && indbef > 1 && (indbef -= 1)
    delnow, indnow = findmin(abs(bound - tmpn)) ; bound[indnow] > tmpn && indnow > 1 && (indnow -= 1)
    delaft, indaft = findmin(abs(bound - tmpa)) ; bound[indaft] > tmpa && indaft > 1 && (indaft -= 1)
    delobs, indobs = findmin(abs(bound - tmpo)) ; bound[indobs] > tmpo && indobs > 1 && (indobs -= 1)
    gridbefnow[indbef,indnow] += 1 ; conbefnow[indnow] += tmpb
    gridaftnow[indaft,indnow] += 1 ; conaftnow[indnow] += tmpa
    gridbefobs[indbef,indobs] += 1 ; conbefobs[indobs] += tmpb ; conobsbef[indbef] += tmpo
    gridnowobs[indnow,indobs] += 1 ; connowobs[indobs] += tmpn ; conobsnow[indnow] += tmpo
    gridaftobs[indaft,indobs] += 1 ; conaftobs[indobs] += tmpa ; conobsaft[indaft] += tmpo
    gridbefaft[indbef,indaft] += 1 ; conbefaft[indaft] += tmpb ; conaftbef[indbef] += tmpa
    push!(regressbef, tmpb)
    push!(regressnow, tmpn)
    push!(regressaft, tmpa)
    push!(regressobs, tmpo)
    tmp = copy(binbef[indnow]) ; push!(tmp, tmpb) ; binbef[indnow] = tmp
    tmp = copy(binnow[indnow]) ; push!(tmp, tmpn) ; binnow[indnow] = tmp
    tmp = copy(binaft[indnow]) ; push!(tmp, tmpa) ; binaft[indnow] = tmp
  end
end

valsb  = Array(Float64, tinuma)                                               # also grid for the extrapolation variance
valsn  = Array(Float64, tinuma)                                               # (versus TOTN) using all values in each bin
valsa  = Array(Float64, tinuma)                                               # (above) and for the nearest CUTOFF values
dista  = Array(Float64, tinuma)                                               # (here, but only save nearest CUTOFFs below)
linuma = tinuma < CUTOFF ? tinuma : CUTOFF

for a = 1:tinuma  vals = float(split(tinea[a])) ; valsb[a] = (vals[TOTB] - intbefnow) /slobefnow ; valsn[a] = vals[TOTN] ; valsa[a] = (vals[TOTA] - intaftnow) /sloaftnow  end
for (z, ranz) in enumerate(bound)
  for a = 1:tinuma  dista[a] = abs(ranz - valsn[a])  end ; lima = sort(dista)[linuma] ; b = 1
  for a = 1:tinuma
    if dista[a] <= lima && b <= linuma
      cutbef[z,b] = valsb[a]
      cutnow[z,b] = valsn[a]
      cutaft[z,b] = valsa[a]
      b += 1
    end
  end
end

fname = ARGS[1] * "." * ARGS[2] * ".extra.dau"
fpa = My.ouvre(fname, "w")                                                    # and save the grids
for (a, vala) in enumerate(bound)
  for (b, valb) in enumerate(bound)
    @printf(fpa, "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n", gridbefnow[b,a], gridaftnow[b,a], gridbefobs[b,a], gridaftobs[b,a], gridbefaft[b,a], gridnowobs[b,a])
  end
end
close(fpa)

sumbefnow = sum(gridbefnow, 1)                                                # as well as the corresponding count and sum
sumaftnow = sum(gridaftnow, 1)                                                # of extrapolation values in each TOTN interval
sumbefobs = sum(gridbefobs, 1) ; sumobsbef = sum(gridbefobs, 2)
sumnowobs = sum(gridnowobs, 1) ; sumobsnow = sum(gridnowobs, 2)
sumaftobs = sum(gridaftobs, 1) ; sumobsaft = sum(gridaftobs, 2)
sumbefaft = sum(gridbefaft, 1) ; sumaftbef = sum(gridbefaft, 2)
fname = ARGS[1] * "." * ARGS[2] * ".extra.sun"
fpa = My.ouvre(fname, "w")
for (a, vala) in enumerate(bound)
  @printf(fpa, "%22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %18.11e %18.11e %22d %8.4f\n",
    sumbefnow[a], conbefnow[a], sumaftnow[a], conaftnow[a],
    sumbefobs[a], conbefobs[a], sumaftobs[a], conaftobs[a],
    sumobsbef[a], conobsbef[a], sumobsaft[a], conobsaft[a],
    sumbefaft[a], conbefaft[a], sumaftbef[a], conaftbef[a],
    sumnowobs[a], connowobs[a], sumobsnow[a], conobsnow[a],
    cov(vec(cutbef[a,:]), vec(cutnow[a,:])), cov(vec(cutaft[a,:]), vec(cutnow[a,:])), length(binbef[a]), bound[a])
end
close(fpa)

fname = ARGS[1] * "." * ARGS[2] * ".extra.reh"                                # and finally save the regression coefficients
fpa = My.ouvre(fname, "w")
(intbefnow, slobefnow) = linreg(regressnow, regressbef)
(intaftnow, sloaftnow) = linreg(regressnow, regressaft)
(intbefobs, slobefobs) = linreg(regressobs, regressbef) ; (intobsbef, sloobsbef) = linreg(regressbef, regressobs)
(intnowobs, slonowobs) = linreg(regressobs, regressnow) ; (intobsnow, sloobsnow) = linreg(regressnow, regressobs)
(intaftobs, sloaftobs) = linreg(regressobs, regressaft) ; (intobsaft, sloobsaft) = linreg(regressaft, regressobs)
(intbefaft, slobefaft) = linreg(regressaft, regressbef) ; (intaftbef, sloaftbef) = linreg(regressbef, regressaft)
            covbefnow  =    cov(regressnow, regressbef)
            covaftnow  =    cov(regressnow, regressaft)
            covbefobs  =    cov(regressobs, regressbef) ;             covobsbef  =    cov(regressbef, regressobs)
            covnowobs  =    cov(regressobs, regressnow) ;             covobsnow  =    cov(regressnow, regressobs)
            covaftobs  =    cov(regressobs, regressaft) ;             covobsaft  =    cov(regressaft, regressobs)
            covbefaft  =    cov(regressaft, regressbef) ;             covaftbef  =    cov(regressbef, regressaft)
@printf(fpa, "%33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f\n",
  intbefnow, slobefnow, intaftnow, sloaftnow,
  intbefobs, slobefobs, intaftobs, sloaftobs,
  intobsbef, sloobsbef, intobsaft, sloobsaft,
  intbefaft, slobefaft, intaftbef, sloaftbef,
  intnowobs, slonowobs, intobsnow, sloobsnow,
             covbefnow,            covaftnow,
             covbefobs,            covaftobs,
             covobsbef,            covobsaft,
             covbefaft,            covaftbef,
             covnowobs,            covobsnow)
close(fpa)
exit(0)

#=
buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_extra.ucur.got2000_obs.comb
dat, lat, lon, obs, bef1, now1, aft1, bef2, now2, aft2, bef3, now3, aft3, bef4, now4, aft4, bef5, now5, aft5

count = 0
for a = 1:length(bound)
  b = length(binbef[a])
  c = b >= 10 ? cov(binbef[a], binnow[a]) : 0
  d = b >= 10 ? cov(binaft[a], binnow[a]) : 0
  f =           cov(vec(cutbef[a,:]), vec(cutnow[a,:]))
  g =           cov(vec(cutaft[a,:]), vec(cutnow[a,:]))
  @printf("%8.4f %8d %18.14f %18.14f %18.14e %18.14e\n", bound[a], b, c, d, f, g)
  count += b
end
exit(0)
=#
