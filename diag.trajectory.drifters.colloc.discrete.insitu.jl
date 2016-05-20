#=
 = Loop through a superset of insitu obs and, based on a subset of locations of
 = interest, write the corresponding subset of observations - RD March 2016.
 =#

using My

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000\n\n")
  exit(1)
end

locs = Set(Array(Tuple{Float64, Float64}, 0))                                 # initialize the set of
fpa = My.ouvre(ARGS[2], "r") ; lines = readlines(fpa) ; close(fpa)            # locations of interest
for line in lines  vals = split(line) ; push!(locs, (float(vals[1]), float(vals[2])))  end

fpa = My.ouvre(ARGS[1],          "r")
fpb = My.ouvre(ARGS[2] * "_obs", "w")
for line in eachline(fpa)
  vals = split(line)
  if in((float(vals[2]), float(vals[3])), locs)
    write(fpb, line)
  end
end
close(fpa)
close(fpb)
exit(0)
