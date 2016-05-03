#=
 = Report on simple metrics of comparison (bias and RMSE) among eight analyses
 = (if available), computed with respect to a corresponding set of ICOADS in situ
 = (ship and buoy) estimates of the same variable.  If any one analysis is missing
 = then ICOADS is considered missing also.  Rescaling analyses prior to comparison
 = permits a simple evaluation of the rescaling - RD March 2016.
 =#

using My
const SHFX             = 1                              # indecies of all data variables
const LHFX             = 2
const WSPD             = 9
const AIRT             = 12
const SSTT             = 14
const SHUM             = 15

const VARN             = 22                             # number of output variables (besides date/lat/lon)
const MISS             = -9999.0                        # generic missing value
const CALIB            = 1                              # flag determining whether to calibrate (1 = yes)

if size(ARGS) != (2,)
  print("\nUsage: jjj $(basename(@__FILE__)) all.flux.daily_2.0_valid_remainder shfx/lhfx/wspd/airt/sstt/shum\n\n")
  exit(1)
end

vind = 0
if     ARGS[2] == "shfx"  vind = SHFX
  cala = [ 2.2,          4.1,     3.4,           3.4,      3.3,     1.9,      3.3,       5.6]
  calb = [0.96,         0.92,    0.82,          0.98,     0.91,    0.87,     0.95,      0.85]
elseif ARGS[2] == "lhfx"  vind = LHFX
  cala = [ 2.2,          4.1,     3.4,           3.4,      3.3,     1.9,      3.3,       5.6]
  calb = [0.96,         0.92,    0.82,          0.98,     0.91,    0.87,     0.95,      0.85]
elseif ARGS[2] == "wspd"  vind = WSPD
  cala = [ 2.2,          4.1,     3.4,           3.4,      3.3,     1.9,      3.3,       5.6]
  calb = [0.96,         0.92,    0.82,          0.98,     0.91,    0.87,     0.95,      0.85]
elseif ARGS[2] == "airt"  vind = AIRT
  cala = [ 2.2,          4.1,     3.4,           3.4,      3.3,     1.9,      3.3,       5.6]
  calb = [0.96,         0.92,    0.82,          0.98,     0.91,    0.87,     0.95,      0.85]
elseif ARGS[2] == "sstt"  vind = SSTT
  cala = [ 2.2,          4.1,     3.4,           3.4,      3.3,     1.9,      3.3,       5.6]
  calb = [0.96,         0.92,    0.82,          0.98,     0.91,    0.87,     0.95,      0.85]
elseif ARGS[2] == "shum"  vind = SHUM
  cala = [ 2.2,          4.1,     3.4,           3.4,      3.3,     1.9,      3.3,       5.6]
  calb = [0.96,         0.92,    0.82,          0.98,     0.91,    0.87,     0.95,      0.85]
end

dirs = ["cfsr", "erainterim", "hoaps", "ifremerflux", "jofuro", "merra", "oaflux", "seaflux"]
(dirn,) = size(dirs)

fpa = My.ouvre(ARGS[1], "r")                                                  # open the insitu and analysis files
fpn = Array(IOStream, 0)
for (a, dira) in enumerate(dirs)
  fpz = My.ouvre(ARGS[1] * "." * dira, "r") ; push!(fpn, fpz)
end

data = Array(Float64, 0)                                                      # read sets of valid data (allowing for
datb = Array(Float64, dirn + 1)                                               # missing cfsr LHFX and jofuro AIRT/SSTT)
for line in eachline(fpa)
  flag = 1
  datb[dirn+1] = float(split(line)[vind])
  if datb[dirn+1] < -333.0 || datb[dirn+1] > 3333.0  flag = 0  end
  for a = 1:dirn
    linz = readline(fpn[a])
    datb[a] = float(split(linz)[vind])
    if (datb[a] < -333.0 || datb[a] > 3333.0) && !(vind == LHFX && a == 1) &&
      !(vind == AIRT && a == 5) && !(vind == SSTT && a == 5)  flag = 0
    end
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
numb = div(length(data), dirn + 1)                                            # shape the 2D data array (9 rows, numb cols)
datc = reshape(data, (dirn + 1, numb))
println("found $numb values valid across analyses")

for a = 1:numb, b = 1:dirn                                                    # subtract the reference in situ from the rest
  datc[b,a] -= datc[dirn+1,a]                                                 # and report the mean and stdev of the diff
end
resa = mean(datc, 2)
resb =  std(datc, 2)
@printf("diff  mean %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %d %s\n", resa[1], resa[2], resa[3], resa[4], resa[5], resa[6], resa[7], resa[8], resa[9], numb, ARGS[2])
@printf("diff stdev %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %d %s\n", resb[1], resb[2], resb[3], resb[4], resb[5], resb[6], resb[7], resb[8], resb[9], numb, ARGS[2])

close(fpa)
for a = 1:dirn                                                                # then close this set of files
  close(fpn[a])
end

if CALIB == 0  tail = ".summ"  end
if CALIB == 1  tail = ".sumc"  end
fpb = My.ouvre(ARGS[1] * "." * ARGS[2] * tail, "w")
form = @sprintf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %d %s\n", resa[1], resa[2], resa[3], resa[4], resa[5], resa[6], resa[7], resa[8], resa[9], numb, ARGS[2]) ; write(fpb, form)
form = @sprintf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %d %s\n", resb[1], resb[2], resb[3], resb[4], resb[5], resb[6], resb[7], resb[8], resb[9], numb, ARGS[2]) ; write(fpb, form)
close(fpb)
exit(0)
