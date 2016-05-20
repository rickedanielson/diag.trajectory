#=
 = Split many drifter observations by location and store the resulting files
 = in an insitu dir (assuming that the observations have been sorted already).
 = First read the locations of interest, as given, for example, by a list of
 = calibration locations and only store drifter obs groups at these locations
 = - RD April 2016.
 =#

using My
const LEN              = 100
const LOTS             = 1000
const TIMS             = 3408                           # number in 6-h timeseries
const START            = 2                              # make START-1 a valid array index

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib buoydata_1993_2014_drogON.asc.nonmdt.sort\n\n")
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
  push!(lins, line)
  if i != START && lins[i][12:32] != lins[i-1][12:32]                         # at end of one location and beginning of
    vals = split(lins[i-1])                                                   # next (where 12:32 cover lat and lon):
    lat = float(vals[2])
    lon = float(vals[3])
    if in((lat, lon), locs)                                                   # if the location is of interest then loop
      tmp = @sprintf("insitu/insitu.%9.3f.%9.3f", lat, lon)                   # through both TIMS and existing obs, and
      tmp = replace(tmp, " ", ".") ; fpc = My.ouvre(tmp, "w")                 # write either missing or valid data lines
      locind = START
      locdat = lins[locind][1:10]
      date = "2012090100"
      while parse(Int, date) < 2015010100
        while locind < i - 1 && parse(Int, date) > parse(Int, locdat)         # note that date shouldn't ever pass locdat
          locind += 1 ; locdat = lins[locind][1:10]                           # unless there is more than one collocated
          println("WARNING : collocation detected")                           # drifter, in which case the first is used
        end
        if date == locdat
          formb = lins[locind]
          if locind < i - 1  locind += 1  end
          locdat = lins[locind][1:10]
        else
          formb = @sprintf("%10s %10.5f %10.5f   99999999     99999  9 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000\n", date, lat, lon)
        end
        write(fpc, formb)
        date = My.dateadd(date, 6, "hr")
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
lat = float(vals[2])
lon = float(vals[3])
if in((lat, lon), locs)
  tmp = @sprintf("insitu/insitu.%9.3f.%9.3f", lat, lon)
  tmp = replace(tmp, " ", ".") ; fpc = My.ouvre(tmp, "w")
  locind = START
  locdat = lins[locind][1:10] 
  date = "2012090100"
  while parse(Int, date) < 2015010100
    while locind < i - 1 && parse(Int, date) > parse(Int, locdat)             # note that date shouldn't ever pass locdat
      locind += 1 ; locdat = lins[locind][1:10]                               # unless there is more than one collocated
      println("WARNING : collocation detected")                               # drifter, in which case the first is used
    end
    if date == locdat
      formb = lins[locind]
      if locind < i - 1  locind += 1  end
      locdat = lins[locind][1:10] 
    else
      formb = @sprintf("%10s %10.5f %10.5f   99999999     99999  9 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000 -9999.00000\n", date, lat, lon)
    end
    write(fpc, formb)
    date = My.dateadd(date, 6, "hr")
  end
  close(fpc)
  n += 1
end

print("wrote $n insitu files\n\n")
exit(0)


#=
for j = START:i-1  write(fpc, lins[j])  end
2014100118    0.00000    2.25000     116130     31910  3 23649.75000     2.19200     0.11800     0.30970    -0.09700     0.20166     0.10230    -0.01952     0.01180     0.00068    -0.02360     0.01224     0.02165    -0.01290     0.06116  9999.00000  9999.00000
=#
