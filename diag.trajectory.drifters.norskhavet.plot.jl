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
const SKIP             = 0                              # number of 6-hourly simulated drifters to skip (3 = daily)

if size(ARGS) != (1,)
  print("\nUsage: jj $(basename(@__FILE__)) lofoten_basin\n\n")
  exit(1)
end
xyzzy = ARGS[1]*".gif"
if isfile(xyzzy)  write(STDERR, "ERROR : $xyzzy already exists\n") ; exit(-1)  end

count = SKIP + 1                                                              # identify all simulated drifter files in
tmpfiles = Array(ASCIIString, 0)                                              # the current dir but skip every SKIP+1
driftfiles = filter(x -> (contains(x, ".traj.") &&                            # file in a driftfiles list
                         !contains(x,    ".nc") && 
                          filesize(x) > 150), readdir("."))
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
print("\ncreated $count track files\n")

dount = 1
datefiles = Array(ASCIIString, 0)
refdate = "2010-01-01-12" ; push!(datefiles, refdate)
for a = 1:4500
  newdate = My.dateadd(refdate, 1, "dy")
  if newdate[1:4] == "2013"  break  end
  push!(datefiles, newdate) ; refdate = newdate
  dount += 1
end
print("created $dount dates to plot\n\n")

minlat =  65.0                                                                # then create the snapshot images
maxlat =  75.0                                                                # (in parallel) and gif animation
minlon = -10.0
maxlon =  20.0
para = "xyzzy_$(basename(@__FILE__))_$(ARGS[1])" ; fpa = My.ouvre(para, "w")
for (a, file) in enumerate(datefiles)
  link  = "../norskhavet/norskhavet_$(file[1:4])$(file[6:7])"
  link *=                         "$(file[9:10])$(file[12:13])"
  link *= "0000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v01.0-fv01.0.nc"
  lims  = @sprintf("%9.3f %9.3f %9.3f %9.3f", minlat, maxlat, minlon, maxlon)
  write(fpa, "grads -blc \"plot.globcurrent.drifter.norskhavet $xyzzy $count $file $link $lims\"\n")
#        run(`grads -blc  "plot.globcurrent.drifter.norskhavet $xyzzy $count $file $link $lims"`)
end
close(fpa)
print("parallel -j 7 -a $para >&/dev/null\n")
# run(`parallel -j 7 -a $para >&/dev/null`)

#print("ssh $(ENV["HOST"]) \"cd $(pwd()) ; convert -delay 0 -loop 0 plot.globcurrent.drift.$(ARGS[1])*png $(ARGS[1]).gif\"\n")
#  run(`ssh $(ENV["HOST"])  "cd $(pwd()) ; convert -delay 0 -loop 0 plot.globcurrent.drift.$(ARGS[1])*png $(ARGS[1]).gif"`)
print("convert -delay 0 -loop 0 plot.globcurrent.drift.*png $(ARGS[1]).gif\n")
# run(`convert -delay 0 -loop 0 plot.globcurrent.drift.*png $(ARGS[1]).gif`)
#for file in My.trouve("plot.globcurrent.drift", ".")  print("rm $file\n") ; rm(file)  end
print("rm $para $xyzzy\n")
#rm(para) ; rm(xyzzy)
exit(0)
