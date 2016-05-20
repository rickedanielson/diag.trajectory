#=
 = Loop through all locations of interest and assemble valid collocations
 = from the various in situ and analysis subdirs.  Ensure that both current
 = components are valid for inclusion - RD May 2016.
 =#

using My
const UCUR             = 10                             # indecies of all data variables
const VCUR             = 11
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs\n\n")
  exit(1)
end

vind = 0                                                                      # identify the output variable
contains(ARGS[1], "ucur") && (vind = UCUR ; vine = VCUR)
contains(ARGS[1], "vcur") && (vind = VCUR ; vine = UCUR)
dirs = ["v2.0_global_025_deg_ekman_15m", "v2.0_global_025_deg_ekman_hs", "v2.0_global_025_deg_geostrophic", "v2.0_global_025_deg_total_15m", "v2.0_global_025_deg_total_hs"]
dirn = length(dirs)

function read_nth_line(fn::AbstractString, ln::Int64)
  stream = open(fn, "r") 
  for i = 1:ln-1  readline(stream)  end
  line =          readline(stream)
  close(stream)
  return line
end

fpa = My.ouvre(ARGS[1],           "r")
fpb = My.ouvre(ARGS[1] * ".comb", "w")

for line in eachline(fpa)                                                     # loop through the insitu locations
  vals = split(line)
  dat =       vals[   1]  ; datind = round(Int, 4 * My.datesous("2012083118", dat, "dy"))
  lat = float(vals[   2])
  lon = float(vals[   3]) ; lon < -180 && (lon += 360) ; lon > 180 && (lon -= 360)
  cur = float(vals[vind])
  out = @sprintf("%s %9.3f %9.3f %9.3f", dat, lat, lon, cur)
  tmp = @sprintf("%9.3f.%9.3f", lat, lon) ; tail = replace(tmp, " ", ".")

  bef = fill(MISS, dirn)                                                      # add analysis bef/aft to insitu data
  aft = fill(MISS, dirn)
  flag = true
  for (a, dira) in enumerate(dirs)
    tmp = split(read_nth_line("$dira/$dira.$tail.bef", datind)) ; bef[a] = float(tmp[vind]) ; chkbef = float(tmp[vine])
    newdat = tmp[1] ; if dat != newdat  println("ERROR : $dat != $newdat") ; exit(-1)  end
    tmp = split(read_nth_line("$dira/$dira.$tail.aft", datind)) ; aft[a] = float(tmp[vind]) ; chkaft = float(tmp[vine])
    newdat = tmp[1] ; if dat != newdat  println("ERROR : $dat != $newdat") ; exit(-1)  end
    if bef[a] < -333.0 || bef[a] > 333.0 || aft[a] < -333.0 || aft[a] > 333.0 ||
       chkbef < -333.0 || chkbef > 333.0 || chkaft < -333.0 || chkaft > 333.0  flag = false  end
  end

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
all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs
2012090706   33.37500  161.87500      81961     52901  3 22895.25000   161.86800    33.30400    -0.50710     0.74060    -0.41013     0.50649     0.12458     0.03886     0.11398     0.02028    -0.26240    -0.12750    -0.09581     0.08826    -0.22333    -0.15568
v2.0_global_025_deg_ekman_15m/v2.0_global_025_deg_ekman_15m.....6.625..-108.125.bef
2013012600    6.62500 -108.12500   99999999     99999  9 -9999.00000 -9999.00000 -9999.00000    -0.04890     0.10378
=#
