#=
 = Split the drifter observations into calibration and validation groups, where the
 = calibration subset is taken to be at locations within each RESOL-degree box that
 = contains the largest number of available observations - RD March, April 2016.
 =#

using My
const RESOL            = 2.0                            # resolution at which to choose the calib subset
const CUTOFF           = 10                             # minimum number of obs at calib location
const LAT              = 1
const LON              = 2
const NUM              = 3
const PARAMS           = 3

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate\n\n")
  exit(1)
end

minlons = collect(-180.0:RESOL:180.0-RESOL)                                   # define a search grid for speed
midlons = minlons + RESOL / 2.0
maxlons = minlons + RESOL
minlats = collect( -90.0:RESOL: 90.0-RESOL)
midlats = minlats + RESOL / 2.0
maxlats = minlats + RESOL
numlons = length(minlons)
numlats = length(minlats)

fpa = My.ouvre(ARGS[1], "r")                                                  # read the available obs counts
lines = readlines(fpa) ; close(fpa)                                           # and convert to float
valn = length(lines)
vals = Array(Float64, valn, PARAMS)
for (a, line) in enumerate(lines)
  tmp = split(line)
  vals[a,LAT] = float(tmp[1])
  vals[a,LON] = float(tmp[2])
  vals[a,NUM] = float(tmp[3])
end

locs = Array(Tuple{Float64, Float64}, 0)                                      # find the calib locations
for a = 1:numlons                                                             # (largest number of drifter obs
  for b = 1:numlats                                                           #  in each gridbox if available)
#   @printf("%8.3f %8.3f\n", midlons[a], midlats[b])
    maxlat = maxlon = maxnum = -1.0
    for c = 1:valn
      if vals[c,NUM] > maxnum && minlons[a] <= vals[c,LON] < maxlons[a] &&
                                 minlats[b] <= vals[c,LAT] < maxlats[b]
        maxlat = vals[c,LAT]
        maxlon = vals[c,LON]
        maxnum = vals[c,NUM]
      end
    end
    if maxnum > CUTOFF
      push!(locs, (maxlat, maxlon))
    end
  end
end

tmp = @sprintf("%s_%.1f_calib", ARGS[1], RESOL) ; fpb = My.ouvre(tmp, "w")
tmp = @sprintf("%s_%.1f_valid", ARGS[1], RESOL) ; fpc = My.ouvre(tmp, "w")
for line in lines
  vals = split(line)
  vlat = float(vals[1])
  vlon = float(vals[2])
  if findfirst(locs, (vlat, vlon)) > 0  write(fpb, line)
  else                                  write(fpc, line)  end
end
close(fpb)
close(fpc)
exit(0)
