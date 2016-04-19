#=
 = Split many daily observations by location and store the resulting files
 = in an insitu dir (assuming that the observations have been sorted already).
 = First read the locations of interest, as given, for example, by a list of
 = calibration locations and only store daily obs groups at these locations
 = - RD March 2016.
 =#

using My
const LEN              = 100
const LOTS             = 1000
const TIMS             = 3745                           # number in timeseries
const START            = 2                              # make START-1 a valid array index

if size(ARGS) != (2,)
  print("\nUsage: jjj $(basename(@__FILE__)) all/all.flux.daily_2.0_locate.calib all.flux.daily.sort\n\n")
  exit(1)
end

locs = Set(Array(Tuple{Float64, Float64}, 0))
lins = Array(UTF8String, 1)

fpa = My.ouvre(ARGS[1], "r")                                                  # read the locations of interest
for line in eachline(fpa)
  vals = split(line)
  lat = float(vals[1])
  lon = float(vals[2])
  push!(locs, (lat, lon))
end
close(fpa)

n = 0 ; i = START                                                             # having initialized arrays with the first
fpb = My.ouvre(ARGS[2], "r")                                                  # entry undefined, starting with the second
for line in eachline(fpb)                                                     # entry, data for a new location is stored
  tmpa = line[1:30]
  tmpb = line[31:38]
  tmpc = line[43:end]
  line = tmpa * "  " * tmpb * "12" * tmpc
  push!(lins, line)
  if i != START && lins[i][44:59] != lins[i-1][44:59]                         # at end of one location and beginning of
    vals = split(lins[i-1])                                                   # next (where 44:59 cover lat and lon):
    lat = float(vals[5])
    lon = float(vals[6])
    if in((lat, lon), locs)                                                   # if the location is of interest then loop
      tmp = @sprintf("insitu/insitu.%9.3f.%9.3f", lat, lon)                   # through both TIMS and existing obs, and
      tmp = replace(tmp, " ", ".") ; fpc = My.ouvre(tmp, "w")                 # write either missing or valid data lines
      locind = START
      locdat = lins[locind][33:42]
      date = "1999100112"
      while parse(Int, date) < 2010010100
        if date == locdat
          formb = lins[locind]
          if locind < i - 1  locind += 1  end
          locdat = lins[locind][33:42]
        else
          formb = @sprintf("-9999.00 -9999.00      0000 %14s %7.3f %8.3f -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00\n",
            date, lat, lon)
        end
        write(fpc, formb)
        date = My.dateadd(date, 1, "dy")
      end
      close(fpc)
      n += 1
    end
    lins = Array(UTF8String, 1) ; push!(lins, line)                           # then reset arrays with the new starting line
    i = START
  end
  i += 1
end
close(fpb)

vals = split(lins[i-1])                                                       # write the last file, if also of interest
lat = float(vals[5])
lon = float(vals[6])
if in((lat, lon), locs)
  tmp = @sprintf("insitu/insitu.%9.3f.%9.3f", lat, lon)
  tmp = replace(tmp, " ", ".") ; fpc = My.ouvre(tmp, "w")
  locind = START
  locdat = lins[locind][33:42] 
  date = "1999100112"
  while parse(Int, date) < 2010010100
    if date == locdat
      formb = lins[locind]
      if locind < i - 1  locind += 1  end
      locdat = lins[locind][33:42] 
    else
      formb = @sprintf("-9999.00 -9999.00      0000 %14s %7.3f %8.3f -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00 -9999.00\n",
        date, lat, lon)
    end
    write(fpc, formb)
    date = My.dateadd(date, 1, "dy")
  end
  close(fpc)
  n += 1
end

print("wrote $n insitu files\n\n")
exit(0)


#=
for j = START:i-1  write(fpc, lins[j])  end
=#
