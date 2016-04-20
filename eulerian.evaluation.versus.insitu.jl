#=
 = Report on simple metrics of comparison (bias and RMSE) relative to drifters.
 = Rescaling GlobCurrent velocity prior to comparison permits a simple evaluation
 = of the impact of rescaling - RD April 2016.
 =#

using My
const UCUR             = 10                             # indecies of all data variables
const VCUR             = 11

const MISS             = -9999.0                        # generic missing value
const CALIB            = 0                              # flag determining whether to calibrate (1 = yes)

if size(ARGS) != (2,)
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs ucur\n\n")
  exit(1)
end

vind = 0
if     ARGS[2] == "ucur"  vind = UCUR
  cala = [ 2.2,          4.1,     3.4,           3.4,      3.3]
  calb = [0.96,         0.92,    0.82,          0.98,     0.91]
elseif ARGS[2] == "vcur"  vind = VCUR
  cala = [ 2.2,          4.1,     3.4,           3.4,      3.3]
  calb = [0.96,         0.92,    0.82,          0.98,     0.91]
end
dirs = ["v2.0_global_025_deg_ekman_15m", "v2.0_global_025_deg_ekman_hs", "v2.0_global_025_deg_geostrophic", "v2.0_global_025_deg_total_15m", "v2.0_global_025_deg_total_hs"]
(dirn,) = size(dirs)

fpa = My.ouvre(ARGS[1], "r")                                                  # open the insitu and analysis files
fpn = Array(IOStream, 0)
for (a, dira) in enumerate(dirs)
  fpz = My.ouvre(ARGS[1] * "." * dira, "r") ; push!(fpn, fpz)
end

data = Array(Float64, 0)                                                      # read sets of valid data
datb = Array(Float64, dirn + 1)
for line in eachline(fpa)
  flag = 1
  datb[dirn+1] = float(split(line)[vind])
  if datb[dirn+1] < -333.0 || datb[dirn+1] > 333.0  flag = 0  end
  for a = 1:dirn
    linz = readline(fpn[a])
    datb[a] = float(split(linz)[vind])
    if datb[a] < -333.0 || datb[a] > 333.0  flag = 0  end
  end

  if CALIB == 1                                                               # rescale valid analysis data as needed
    for a = 1:dirn  if datb[a] != MISS
      datb[a] = (datb[a] - cala[a]) / calb[a]
    end  end
  end

  if flag == 1
    for a = 1:dirn + 1  push!(data, datb[a])  end
  end
end
numb = div(length(data), dirn + 1)                                            # shape the 2D data array (6 rows, numb cols)
datc = reshape(data, (dirn + 1, numb))
println("found $numb values valid across analyses")

for a = 1:numb, b = 1:dirn                                                    # subtract the reference in situ from the rest
  datc[b,a] -= datc[dirn+1,a]                                                 # and report the mean and stdev of the diff
end
resa = mean(datc, 2)
resb =  std(datc, 2)
@printf("diff  mean %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %d %s\n", resa[1], resa[2], resa[3], resa[4], resa[5], resa[6], numb, ARGS[2])
@printf("diff stdev %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %d %s\n", resb[1], resb[2], resb[3], resb[4], resb[5], resb[6], numb, ARGS[2])

close(fpa)
for a = 1:dirn                                                                # then close this set of files
  close(fpn[a])
end

if CALIB == 0  tail = ".summ"  end
if CALIB == 1  tail = ".sumc"  end
fpb = My.ouvre(ARGS[1] * "." * ARGS[2] * tail, "w")
form = @sprintf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %d %s\n", resa[1], resa[2], resa[3], resa[4], resa[5], resa[6], numb, ARGS[2]) ; write(fpb, form)
form = @sprintf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %d %s\n", resb[1], resb[2], resb[3], resb[4], resb[5], resb[6], numb, ARGS[2]) ; write(fpb, form)
close(fpb)
exit(0)
