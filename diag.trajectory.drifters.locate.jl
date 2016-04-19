#=
 = Identify and count the location of all observations in the input file,
 = where locations are defined at the resolution of a grid (for subsequent
 = collocation) - RD April 2016.
 =#

using My

if size(ARGS) != (1,)
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt\n")
  exit(1)
end

lats = collect( -90.0:0.25:89.75)                                             # then define the collocation grid
lons = collect(-180.0:0.25:179.75)                                            # and initialize the count
subs = Set(Array(Tuple{Int64, Int64}, 0))
numb = zeros(length(lons), length(lats))

fpa = My.ouvre(ARGS[1],             "r")
fpb = My.ouvre(ARGS[1] * ".locate", "w")

for line in readlines(fpa)                                                    # identify and count the collocations
  vals = split(line)
  indlat = findin(lats, float(vals[2]))[1]
  indlon = findin(lons, float(vals[3]))[1]
  push!(subs, (indlat, indlon))
  numb[indlon,indlat] += 1.0
end

for loc in subs
  (indlat, indlon) = loc
  line = @sprintf("%8.2f %8.2f %8.0f\n", lats[indlat], lons[indlon], numb[indlon,indlat])
  write(fpb, line)
end

close(fpa)
close(fpb)
exit(0)
