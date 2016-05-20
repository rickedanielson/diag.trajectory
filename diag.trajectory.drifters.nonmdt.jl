#=
 = Isolate the drifters not in the CNES-CLS-2013 MDT (so pass only data
 = from September 2012 onward), convert the CNES Julian day to a regular
 = date, and add the location of the nearest gridbox - RD April 2016.
 =#

using My

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2012_drogON.asc\n\n")
  exit(1)
end

lats = collect( -79.875:0.25:79.875)                                          # define the target grid
lons = collect(-179.875:0.25:179.875)

fpa = My.ouvre(ARGS[1],             "r")
fpb = My.ouvre(ARGS[1] * ".nonmdt", "w")

for line in eachline(fpa)
  vars = split(line)
  if float(vars[4]) >= 22889
    jday = My.dateadd("1950010100", float(vars[4]), "dy")
    lat = float(vars[6])
    lon = float(vars[5]) ; lon < -180 && (lon += 360) ; lon > 180 && (lon -= 360)
    dellat, indlat = findmin(abs(lats - lat))
    dellon, indlon = findmin(abs(lons - lon))
    repos = @sprintf("%10.5f %10.5f", lats[indlat], lons[indlon])
    write(fpb, jday * " " * repos * " " * line)
  end
end

close(fpa)
close(fpb)
exit(0)
