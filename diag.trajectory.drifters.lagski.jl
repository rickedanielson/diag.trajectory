#=
 = Plot the mean and standard deviation of Lagrangian separation distance
 = versus time from launch.  All available trajectory files in the current
 = dir are identified and if these number more than a few thousand, every
 = few files are skipped.  Possibly restrict the trajectory files to those
 = with data at all times (25 days), but in any case, omit stats that are
 = not meaningful (too few samples) - RD August 2015
 =#

using My
const LEN              = 101
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1
  print("\nUsage: jj $(basename(@__FILE__)) lagski\n\n")
  exit(1)
end

driftfiles = filter(x -> (contains(x, ".traj") &&                             # identify all simulated drifter files
                         !contains(x,   ".nc") &&                             # in the current dir
                         !contains(x, "plot.")), readdir("."))
(count,) = size(driftfiles)
print("\n  found $count trajectory files\n")

tmpfiles = Array(UTF8String, 0)
for (a, fila) in enumerate(driftfiles)                                        # first restrict the trajectory files to
  if filesize(fila) > 10700                                                   # those with data at all times (25 days)
    push!(tmpfiles,fila)
  end
end
driftfiles = tmpfiles
(count,) = size(driftfiles) ; skip = int(count / 2000)
print("  among $count files valid at all times, skipping every $skip\n")

count = 1
tmpfiles = Array(UTF8String, 0)                                               # then reduce the number of trajectory
for (a, fila) in enumerate(driftfiles)                                        # files used to compute statistics
  if count == skip + 1                                                        # (skip = 3 is daily)
    push!(tmpfiles,fila)
    count = 0
  end
  count += 1
end
driftfiles = tmpfiles
(count,) = size(driftfiles)

trkfiles = Set{UTF8String}()                                                  # and identify the number of drifters
for (a, fila) in enumerate(driftfiles)                                        # associated with these trajectories
  push!(trkfiles, fila[1:24])
end
dount = length(trkfiles) ; plural = dount == 1 ? "drifter" : "drifters"
print("reading $count trajectory files associated with $dount individual $plural\n")

hours  = [0.0:6:600] / 24.0                                                   # set up the arrays and loop through the
dissum = zeros(hours)                                                         # the files to get sum and sum of squares
dissqr = zeros(hours)                                                         # at each relative time
disnum = zeros(hours)
for (a, fila) in enumerate(driftfiles)
  fpa = My.ouvre(fila, "r", false)
  for (b, linb) in enumerate(eachline(fpa))
    tmp = split(linb)
    driind = int((int(tmp[4]) + 6) / 6)
    drichk =    float(tmp[6])
    dridis =    float(tmp[9])
    if dridis >= 0 && drichk > -361
      dissum[driind] += dridis
      dissqr[driind] += dridis^2.0
      disnum[driind] += 1.0
#     if driind == 45  print("$dridis\n")  end
    end
  end
  close(fpa)
end

disavg = zeros(hours)                                                         # get each mean and standard deviation
disstd = zeros(hours)                                                         # (but only if they are "meaningful")
for a = 1:LEN                                                                 # then plot the timeseries
  if disnum[a] > 10
    disavg[a] =                           dissum[a] / disnum[a]
    disstd[a] = ((dissqr[a] - dissum[a] * dissum[a] / disnum[a]) / (disnum[a] - 1.0))^0.5
  else
    disavg[a] = disstd[a] = NaN
  end
# print("$a $(disnum[a]) $(disavg[a]) $(disstd[a])\n")
end
ymin = minimum(disavg-disstd)
ymax = maximum(disavg+disstd)
disavg[1] = NaN
disstd[1] = NaN

xyzzy = split(pwd(),"/") ; curdir = replace(xyzzy[end], "_", "-") ; xyzzy = ARGS[1]*"-$curdir.png"
print("writing $xyzzy\n")
using Winston
tmp = Winston.FramedPlot(
        title="$curdir\n($dount $plural and $count trajectories)",
        xlabel="Days after launch", xrange = (0,25),
        ylabel="Cumulative skill",  yrange = (-0.1,1.0))
ppp = Winston.add(tmp)
#for a = 1:4:LEN
#  tmp = Winston.DataLabel(hours[a], 240, "$(int(disnum[a]))", "textangle", 90.0, "texthalign", "right", "size", 1.1)
#        Winston.add(ppp, tmp)
#end
tmp = Winston.Curve(hours, disavg, "color", parse(Winston.Colorant, "red"))
      Winston.add(ppp, tmp)
tmp = Winston.SymmetricErrorBarsY(hours[1:4:end], disavg[1:4:end], disstd[1:4:end], "color", parse(Winston.Colorant, "red"))
      Winston.add(ppp, tmp)
      Winston.savefig(ppp, xyzzy, "width", 700, "height", 700)
exit(0)


#       ylabel="Lagrangian separation distance (km)", yrange = (ymin,ymax))
#;      Winston.style(pppc, "linetype", "dotted")
#if isfile(xyzzy)  write(STDERR, "ERROR : $xyzzy already exists\n") ; exit(-1)  end
#pppc = Winston.Curve(hours, disavg+disstd, "color", parse(Winston.Colorant, "blue"))
#;      Winston.setattr(pppc, label="± Standard deviation")
#pppd = Winston.Curve(hours, disavg-disstd, "color", parse(Winston.Colorant, "blue"))
#pppe = Winston.Legend(.1, .7, {pppb,pppc})
