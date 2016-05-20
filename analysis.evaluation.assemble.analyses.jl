#=
 = Loop through all analyses and extract variables of interest at the position
 = of a set of ICOADS ship and buoy observations.  Include a conversion of units,
 = as required - RD March 2016.
 =#

using My, NetCDF
const MISS             = -9999.0                        # generic missing value
const START            = 2                              # make START-1 a valid array index

const SHFX             = 1                              # identify indecies of all data variables
const LHFX             = 2
const WSPD             = 3
const AIRT             = 4
const SSTT             = 5
const SHUM             = 6
const PARAMS           = 6

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) all.flux.daily_2.0_valid_remainder merra\n\n")
  exit(1)
end

tail = "120000-OHF-L4-global_daily_0.25x0.25-v0.7-f01.0.nc"
vars = ["surface_upward_sensible_heat_flux", "surface_upward_latent_heat_flux", "wind_speed",
        "air_temperature", "sea_surface_temperature", "air_surface_specific_humidity"]
varn = length(vars)

ARGS[2] ==        "cfsr" && (stem =        "cfsr-" ; vars[LHFX] = "" ; vars[AIRT] = "2t")
ARGS[2] ==  "erainterim" && (stem =  "erainterim-")
ARGS[2] ==       "hoaps" && (stem =       "hoaps-")
ARGS[2] == "ifremerflux" && (stem = "ifremerflux-")
ARGS[2] ==      "jofuro" && (stem =      "jofuro-" ; vars[AIRT] = "" ; vars[SSTT] = "")
ARGS[2] ==       "merra" && (stem =       "merra-")
ARGS[2] ==      "oaflux" && (stem =      "oaflux-")                           # define the available variables
ARGS[2] ==     "seaflux" && (stem =     "seaflux-")                           # and the filename stems/ending

lats =  -90.0:0.25:89.75  ; latn = length(lats)
lons = -180.0:0.25:179.75 ; lonn = length(lons)

fpa = My.ouvre(ARGS[1],                 "r")
fpb = My.ouvre(ARGS[1] * "." * ARGS[2], "w")

n = 0 ; i = START                                                             # having initialized arrays with the first
lins = Array(UTF8String, 1)                                                   # entry undefined, starting with the second
for line in eachline(fpa)                                                     # entry, accumulate data for the current date
  push!(lins, line)                                                           # and process these once a new date is read
  if i != START && lins[i][31:38] != lins[i-1][31:38]                         # (careful: try introduces new scope)
    stor = fill(MISS, lonn, latn, PARAMS)
    date = lins[i-1][31:38]
    file = ARGS[2] * "/" * stem * date * tail ; println(file)
    if isfile(file)
      for a = 1:PARAMS
        flag = dat3 = false
        if vars[a] != ""
          flag = nc = true
                   try    nc = NetCDF.open(file, mode=NC_NOWRITE, readdimvar=false)          catch;  flag = false  end
          if flag  try  dat3 = NetCDF.readvar(nc, vars[a], start=[1,1,1], count=[-1,-1,-1])  catch;  flag = false  end  end
                   try         NetCDF.close(nc)                                                                    end
        end
        if flag
          for b = 1:latn, c = 1:lonn
            stor[c,b,a] = dat3[c,b,1]
          end
        end
      end

      if ARGS[2] == "cfsr"
        for b = 1:latn, c = 1:lonn
          if 200.0 < stor[c,b,AIRT] < 350.0  stor[c,b,AIRT] -= 273.15  end
        end
      end
      for b = 1:latn, c = 1:lonn
        if 200.0 < stor[c,b,SSTT] < 350.0  stor[c,b,SSTT] -= 273.15  end
      end
    end

    for a = START:i - 1                                                       # echo the ICOADS lines, but with analysis
      vals = split(lins[a])                                                   # data replacing some variables
      lat = float(vals[5])
      lon = float(vals[6]) ; lon < -180 && (lon += 360) ; lon > 180 && (lon -= 360)
      dellat, indlat = findmin(abs(lats - lat))
      dellon, indlon = findmin(abs(lons - lon))
      form = @sprintf("%8.2f %8.2f %s %8.3f %s %8.2f %s %8.2f %8.3f %s\n",
        stor[indlon,indlat,SHFX], stor[indlon,indlat,LHFX], lins[a][ 19: 77], stor[indlon,indlat,WSPD], lins[a][ 88:104],
        stor[indlon,indlat,AIRT], lins[a][115:122], stor[indlon,indlat,SSTT], stor[indlon,indlat,SHUM], lins[a][142:167])
      write(fpb, form)
    end
    lins = Array(UTF8String, 1) ; push!(lins, line)                           # then reset arrays with the new starting line
    n += 1 ; i = START
  end
  i += 1
end

stor = fill(MISS, lonn, latn, PARAMS)                                         # write the last set of dates
date = lins[i-1][31:38]
file = ARGS[2] * "/" * stem * date * tail ; println(file)
if isfile(file)
  for a = 1:PARAMS
    flag = dat3 = false
    if vars[a] != ""
      flag = nc = true
               try    nc = NetCDF.open(file, mode=NC_NOWRITE, readdimvar=false)          catch;  flag = false  end
      if flag  try  dat3 = NetCDF.readvar(nc, vars[a], start=[1,1,1], count=[-1,-1,-1])  catch;  flag = false  end  end
               try         NetCDF.close(nc)                                                                    end
    end
    if flag
      for b = 1:latn, c = 1:lonn
        stor[c,b,a] = dat3[c,b,1]
      end
    end
  end

  if ARGS[2] == "cfsr"
    for b = 1:latn, c = 1:lonn
      if 200.0 < stor[c,b,AIRT] < 350.0  stor[c,b,AIRT] -= 273.15  end
    end
  end
  for b = 1:latn, c = 1:lonn
    if 200.0 < stor[c,b,SSTT] < 350.0  stor[c,b,SSTT] -= 273.15  end
  end
end

for a = START:i - 1                                                       # echo the ICOADS lines, but with analysis
  vals = split(lins[a])                                                   # data replacing some variables
  lat = float(vals[5])
  lon = float(vals[6]) ; lon < -180 && (lon += 360) ; lon > 180 && (lon -= 360)
  dellat, indlat = findmin(abs(lats - lat))
  dellon, indlon = findmin(abs(lons - lon))
  form = @sprintf("%8.2f %8.2f %s %8.3f %s %8.2f %s %8.2f %8.3f %s\n",
    stor[indlon,indlat,SHFX], stor[indlon,indlat,LHFX], lins[a][ 19: 77], stor[indlon,indlat,WSPD], lins[a][ 88:104],
    stor[indlon,indlat,AIRT], lins[a][115:122], stor[indlon,indlat,SSTT], stor[indlon,indlat,SHUM], lins[a][142:167])
  write(fpb, form)
end
n += 1
print("wrote $n dates\n\n")

close(fpa)
close(fpb)
exit(0)


#=
    2.50   107.13      0000   199910010000   1.000  -28.000  1014.00  180.000    9.000    0.000    9.000    26.50    23.80    27.30   18.282  1998.87    10.00    10.00
   17.73   151.63 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    7.700 xxxxxxxxxxxxxxxxx    24.00 xxxxxxxx    26.30   14.934 xxxxxxxxxxxxxxxxxxxxxxxxxx

form = @sprintf("%8.2f %8.2f %9s %14s %7.3f %8.3f %8.2f %8.3f %8.3f %8.3f %8.3f %8.2f %8.2f %8.2f %8.3f %8.2f %8.2f %8.2f\n",
        shf, lhf, iden, date, val[:LAT], val[:LON], val[:SLP], val[:D], val[:W], uwnd, vwnd, val[:AT], val[:DPT], val[:SST], shum, blay, hgtu, hgtt)
=#
