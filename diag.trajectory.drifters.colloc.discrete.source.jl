#=
 = Loop through all locations of interest and assemble valid collocations
 = from the various in situ and analysis subdirs - RD February, March 2016.
 =#

using My
const SHFX             = 1                              # indecies of all data variables
const LHFX             = 2
const WSPD             = 9
const AIRT             = 12
const SSTT             = 14
const SHUM             = 15
const MISS             = -9999.0                        # generic missing value

if size(ARGS) != (1,)
  print("\nUsage: jjj $(basename(@__FILE__)) all.flux.daily.locate_2.0_calib.airt.got2000_obs\n\n")
  exit(1)
end

vind = 0                                                                      # allow the number of output variables
contains(ARGS[1], "shfx") && (vind = SHFX)                                    # to be a function of the variable (i.e.,
contains(ARGS[1], "lhfx") && (vind = LHFX)                                    # omit cfsr LHFX and jofuro AIRT/SSTT)
contains(ARGS[1], "wspd") && (vind = WSPD)
contains(ARGS[1], "airt") && (vind = AIRT)
contains(ARGS[1], "sstt") && (vind = SSTT)
contains(ARGS[1], "shum") && (vind = SHUM)
if contains(ARGS[1], "lhfx")
  dirs = [        "erainterim", "hoaps", "ifremerflux", "jofuro", "merra", "oaflux", "seaflux"] ; fend = ".coml"
elseif contains(ARGS[1], "airt") || contains(ARGS[1], "sstt")
  dirs = ["cfsr", "erainterim", "hoaps", "ifremerflux",           "merra", "oaflux", "seaflux"] ; fend = ".comt"
else
  dirs = ["cfsr", "erainterim", "hoaps", "ifremerflux", "jofuro", "merra", "oaflux", "seaflux"] ; fend = ".comb"
end
(dirn,) = size(dirs)

fpa = My.ouvre(ARGS[1],        "r")
fpb = My.ouvre(ARGS[1] * fend, "w")

for line in eachline(fpa)                                                     # loop through the insitu locations
  vals = split(line)
  dat =       vals[ 4][1:8] ; datind = round(Int, My.datesous("19990930", dat, "dy"))
  lat = float(vals[ 5])
  lon = float(vals[ 6]) ; lon < -180 && (lon += 360) ; lon > 180 && (lon -= 360)
  flx = float(vals[vind])
  hum = float(vals[SHUM])
  spd = float(vals[WSPD])
  air = float(vals[AIRT])
  sst = float(vals[SSTT])
  out = @sprintf("%s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f", dat, lat, lon, flx, hum, spd, air, sst)
  tmp = @sprintf("%9.3f.%9.3f", lat, lon) ; tail = replace(tmp, " ", ".")

  bef = fill(MISS, dirn)                                                      # add analysis bef/aft to insitu data
  aft = fill(MISS, dirn)
  flag = true
  for (a, dira) in enumerate(dirs)
    fpc = My.ouvre("$dira/$dira.$tail.bef", "r", false) ; linb = readlines(fpc) ; close(fpc) ; fpc = 0
    fpd = My.ouvre("$dira/$dira.$tail.aft", "r", false) ; lina = readlines(fpd) ; close(fpd) ; fpd = 0

    tmp = split(linb[datind]) ; bef[a] = float(tmp[vind])
                                newdat = tmp[4][1:8] ; if dat != newdat  println("ERROR : $dat != newdat")  end
    tmp = split(lina[datind]) ; aft[a] = float(tmp[vind])
                                newdat = tmp[4][1:8] ; if dat != newdat  println("ERROR : $dat != newdat")  end
    if bef[a] < -333.0 || bef[a] > 333.0 || aft[a] < -333.0 || aft[a] > 333.0  flag = false  end
  end
  gc()

  if flag                                                                     # and store the line if all values exist
    for (a, dira) in enumerate(dirs)
      tmp = @sprintf(" %9.3f %9.3f", bef[a], aft[a]) ; out *= tmp
    end
    out *= "\n" ; write(fpb, out)
  end
end

close(fpa)
close(fpb)
exit(0)


#=
all/all.flux.daily.locate_2.0_calib.airt.got2000_obs
  -12.89    60.86      0000   199910010000 -30.500   15.250  1012.00  160.000   10.300   -3.523    9.679    17.00    13.50    16.50    9.559  1335.47    10.00    10.00
cfsr/cfsr...-30.500....15.250.aft
    3.62 -9999.00      0000     1999100112 -30.500   15.250 -9999.00 -9999.000    4.922 -9999.000 -9999.000    15.71 -9999.00    16.36    8.224 -9999.00 -9999.00 -9999.00
=#
