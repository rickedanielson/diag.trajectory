#=
 = Create an animation comparing an actual drifter with simulated drifters
 = that are launched at regular intervals from the actual drifter positions,
 = but whose drift is advected following a GlobCurrent combined near-surface
 = current.  The smallest display that captures all visible trajectories is
 = determined from the active tracks - RD July 2015
 =#

using My
const LEN              = 101
const LOTS             = 10000
const MISS             = -9999.0                        # generic missing value
const SKIP             = 3                              # number of 6-hourly simulated drifters to skip (3 = daily)

if (argc = length(ARGS)) != 1
  print("\nUsage: jj $(basename(@__FILE__)) 2010-10-04-1800.00070972\n\n")
  exit(1)
end
xyzzy = ARGS[1]*".gif"
if isfile(xyzzy)  write(STDERR, "ERROR : $xyzzy already exists\n") ; exit(-1)  end

count = SKIP + 1                                                              # identify all simulated drifter files in
tmpfiles = Array(UTF8String, 0)                                               # the current dir but skip every SKIP+1
driftfiles = filter(x -> (contains(x, ARGS[1]*".traj") &&                     # file in a driftfiles list
                         !contains(x,           ".nc") &&
                         !contains(x,         "plot.")), readdir("."))
for (a, fila) in enumerate(driftfiles)
  if count == SKIP + 1
    push!(tmpfiles,fila)
    count = 0
  end
  count += 1
end
driftfiles = tmpfiles

count = 0
xyzzy = ARGS[1]*".xyzzy" ; fpa = My.ouvre(xyzzy, "w")                         # for each unskipped drifter store the
for file in driftfiles                                                        # netcdf trajectory and list all files
  print("nc.template.track.drifters $file\n")                                 # in a temporary plot file
    run(`nc.template.track.drifters $file`)
  write(fpa, "$(file).nc\n")
  count += 1
end
close(fpa)
print("\ncreated $count track files to plot\n\n")

linlat = linlon =  361.0
laxlat = laxlon = -361.0
minlat = minlon =  361.0
maxlat = maxlon = -361.0
simlat = Array(Float64, count, LEN)                                           # read all the trajectory files and
simlon = Array(Float64, count, LEN)                                           # identify the limits of the actual
for (a, fila) in enumerate(driftfiles)                                        # drifter and all its simulations
  fpa = My.ouvre(fila, "r")                                                   # (set invalid positions to MISS)
  b = 0
  for (b, linb) in enumerate(eachline(fpa))
    tmp = split(linb)
    simlat[a,b] = float(tmp[5])
    simlon[a,b] = float(tmp[6])
    if simlon[a,b] > 180.0  simlon[a,b] -= 360.0  end
    if simlat[a,1] > -361 && simlat[a,1] < linlat  linlat = simlat[a,1]  end
    if                       simlat[a,1] > laxlat  laxlat = simlat[a,1]  end
    if simlon[a,1] > -361 && simlon[a,1] < linlon  linlon = simlon[a,1]  end
    if                       simlon[a,1] > laxlon  laxlon = simlon[a,1]  end
    if simlat[a,b] > -361 && simlat[a,b] < minlat  minlat = simlat[a,b]  end
    if                       simlat[a,b] > maxlat  maxlat = simlat[a,b]  end
    if simlon[a,b] > -361 && simlon[a,b] < minlon  minlon = simlon[a,b]  end
    if                       simlon[a,b] > maxlon  maxlon = simlon[a,b]  end
  end
  for c = b+1:LEN
    simlat[a,c] = simlon[a,c] = MISS
  end
  close(fpa)
end
dellat = 0.1 * (laxlat - linlat)                                              # define a region of interest somewhere
dellon = 0.1 * (laxlon - linlon)                                              # between the actual and simulated drifter
linlat -= dellat ; laxlat += dellat                                           # limits
linlon -= dellon ; laxlon += dellon
if minlat < linlat  linlat = (linlat + 2.0 * minlat) / 3.0  end
if maxlat > laxlat  laxlat = (laxlat + 2.0 * maxlat) / 3.0  end
if minlon < linlon  linlon = (linlon + 2.0 * minlon) / 3.0  end
if maxlon > laxlon  laxlon = (laxlon + 2.0 * maxlon) / 3.0  end

dellat = laxlat - linlat ; cenlat = (laxlat + linlat) / 2.0                   # and constrain the aspect ratio of small
dellon = laxlon - linlon ; cenlon = (laxlon + linlon) / 2.0                   # domains
asplat = dellon * 2.0 / 3.0 ; if asplat > 20  asplat = 20.0  end
asplon = dellat * 3.0 / 2.0 ; if asplon > 30  asplon = 30.0  end
if dellat < 20 && asplat > dellat  linlat = cenlat - asplat / 2.0 ; laxlat = cenlat + asplat / 2.0  end
if dellon < 30 && asplon > dellon  linlon = cenlon - asplon / 2.0 ; laxlon = cenlon + asplon / 2.0  end

minlat = Array(Float64, count)                                                # then set this as the fixed viewing domain
maxlat = Array(Float64, count)                                                # (except if a moving domain is required)
minlon = Array(Float64, count)
maxlon = Array(Float64, count)
for (a, file) in enumerate(driftfiles)
  minlat[a] = linlat
  maxlat[a] = laxlat
  minlon[a] = linlon
  maxlon[a] = laxlon
end

if laxlat - linlat > 20                                                       # but constrain large domains in longitude
  for (a, file) in enumerate(driftfiles)                                      # (i.e., as a moving domain)
    minlat[a] = simlat[a,1] - 10.0
    maxlat[a] = simlat[a,1] + 10.0
  end
  minlat = lisse(minlat, [1.0, 1, 1, 1, 1]) ; minlat = lisse(minlat, [1.0, 1, 1, 1, 1]) ; minlat = lisse(minlat, [1.0, 1, 1, 1, 1])
  maxlat = lisse(maxlat, [1.0, 1, 1, 1, 1]) ; maxlat = lisse(maxlat, [1.0, 1, 1, 1, 1]) ; maxlat = lisse(maxlat, [1.0, 1, 1, 1, 1])
  for (a, file) in enumerate(driftfiles)
    if minlat[a] < linlat  minlat[a] = linlat  end
    if maxlat[a] > laxlat  maxlat[a] = laxlat  end
  end
end

if laxlon - linlon > 30                                                       # as well as in latitude (less common)
  for (a, file) in enumerate(driftfiles)
    minlon[a] = simlon[a,1] - 15.0
    maxlon[a] = simlon[a,1] + 15.0
  end
  minlon = lisse(minlon, [1.0, 1, 1, 1, 1]) ; minlon = lisse(minlon, [1.0, 1, 1, 1, 1]) ; minlon = lisse(minlon, [1.0, 1, 1, 1, 1])
  maxlon = lisse(maxlon, [1.0, 1, 1, 1, 1]) ; maxlon = lisse(maxlon, [1.0, 1, 1, 1, 1]) ; maxlon = lisse(maxlon, [1.0, 1, 1, 1, 1])
  for (a, file) in enumerate(driftfiles)
    if minlon[a] < linlon  minlon[a] = linlon  end
    if maxlon[a] > laxlon  maxlon[a] = laxlon  end
  end
end

for (a, file) in enumerate(driftfiles)                                        # then create the snapshot images
  link  = "../links/$(file[31:34])$(file[36:37])"                             # and gif animation and clean up
  link *=          "$(file[39:40])$(file[42:43])"
  link *= "0000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v01.0-fv01.0.nc"
  lims  = @sprintf("%9.3f %9.3f %9.3f %9.3f", minlat[a], maxlat[a], minlon[a], maxlon[a])
  print("grads -blc \"plot.globcurrent.drifter $a $file $link $lims\"\n")
    run(`grads -blc  "plot.globcurrent.drifter $a $file $link $lims"`)
end

#print("ssh $(ENV["HOST"]) \"cd $(pwd()) ; convert -delay 0 -loop 0 plot.globcurrent.drift.$(ARGS[1])*png $(ARGS[1]).gif\"\n")
#  run(`ssh $(ENV["HOST"])  "cd $(pwd()) ; convert -delay 0 -loop 0 plot.globcurrent.drift.$(ARGS[1])*png $(ARGS[1]).gif"`)
print("convert -delay 0 -loop 0 plot.globcurrent.drift.$(ARGS[1])*png $(ARGS[1]).gif\n")
  run(`convert -delay 0 -loop 0 plot.globcurrent.drift.$(ARGS[1])*png $(ARGS[1]).gif`)
for file in My.trouve("plot.globcurrent.drift.$(ARGS[1])", ".")  rm(file)  end
rm(xyzzy)
exit(0)
