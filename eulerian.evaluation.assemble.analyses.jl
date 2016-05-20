#=
 = Loop through all current components and extract variables of interest
 = at the position of a set of drifter observations (accommodating CLS's
 = geostrophic grid difference) - RD April 2016.
 =#

using My, NetCDF
const MISS             = -9999.0                        # generic missing value
const START            = 2                              # make START-1 a valid array index

const UCUR             = 1                              # identify indecies of all data variables
const VCUR             = 2
const PARAMS           = 2

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs v2.0_global_025_deg_total_15m\n\n")
  exit(1)
end

vars = Array(UTF8String, 2)
ARGS[2] == "v2.0_global_025_deg_ekman_15m"   && (vars[UCUR] =       "eastward_ekman_current_velocity" ; vars[VCUR] =       "northward_ekman_current_velocity")
ARGS[2] == "v2.0_global_025_deg_ekman_hs"    && (vars[UCUR] =       "eastward_ekman_current_velocity" ; vars[VCUR] =       "northward_ekman_current_velocity")
ARGS[2] == "v2.0_global_025_deg_geostrophic" && (vars[UCUR] = "eastward_geostrophic_current_velocity" ; vars[VCUR] = "northward_geostrophic_current_velocity")
ARGS[2] == "v2.0_global_025_deg_total_15m"   && (vars[UCUR] =    "eastward_eulerian_current_velocity" ; vars[VCUR] =    "northward_eulerian_current_velocity")
ARGS[2] == "v2.0_global_025_deg_total_hs"    && (vars[UCUR] =    "eastward_eulerian_current_velocity" ; vars[VCUR] =    "northward_eulerian_current_velocity")

ARGS[2] == "v2.0_global_025_deg_ekman_15m"   && (tail = "0000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v02.0-fv01.0.nc")
ARGS[2] == "v2.0_global_025_deg_ekman_hs"    && (tail = "0000-GLOBCURRENT-L4-CURekm_hs-ERAWS_EEM-v02.0-fv01.0.nc")
ARGS[2] == "v2.0_global_025_deg_geostrophic" && (tail = "0000-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v02.0-fv01.0.nc")
ARGS[2] == "v2.0_global_025_deg_total_15m"   && (tail = "0000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v02.0-fv01.0.nc")
ARGS[2] == "v2.0_global_025_deg_total_hs"    && (tail = "0000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc")

lats =  -79.875:0.25:79.875  ; latn = length(lats)
lons = -179.875:0.25:179.875 ; lonn = length(lons)
if ARGS[2] == "v2.0_global_025_deg_geostrophic"
  lats = -89.875:0.25:89.875  ; latn = length(lats)                           # for CLS geostrophy
  lons =   0.125:0.25:359.875 ; lonn = length(lons)
end

fpa = My.ouvre(ARGS[1],                 "r")
fpb = My.ouvre(ARGS[1] * "." * ARGS[2], "w")

n = 0 ; i = START                                                             # having initialized arrays with the first
lins = Array(UTF8String, 1)                                                   # entry undefined, starting with the second
for line in eachline(fpa)                                                     # entry, accumulate data for the current date
  push!(lins, line)                                                           # and process these once a new date is read
  if i != START && lins[i][1:10] != lins[i-1][1:10]                           # (careful: try introduces new scope)
    stor = fill(MISS, lonn, latn, PARAMS)
    if ARGS[2] == "v2.0_global_025_deg_geostrophic"  date = lins[i-1][1:8] * "00"
    else                                             date = lins[i-1][1:10]  end
    file = ARGS[2] * "/" * date * tail ; println(file)
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
    end

    for a = START:i - 1                                                       # echo the drifter data lines, but with analysis
      vals = split(lins[a])                                                   # velocity replacing drifter velocity
      lat = float(vals[2])
      lon = float(vals[3])
      if ARGS[2] == "v2.0_global_025_deg_geostrophic"  lon <    0 && (lon += 360) ; lon > 360 && (lon -= 360)
      else                                             lon < -180 && (lon += 360) ; lon > 180 && (lon -= 360)  end
      dellat, indlat = findmin(abs(lats - lat))
      dellon, indlon = findmin(abs(lons - lon))
      form = @sprintf("%s %11.5f %11.5f %s", lins[a][1:92], stor[indlon,indlat,UCUR], stor[indlon,indlat,VCUR], lins[a][118:end])
      write(fpb, form)
    end
    lins = Array(UTF8String, 1) ; push!(lins, line)                           # then reset arrays with the new starting line
    n += 1 ; i = START
  end
  i += 1
end

stor = fill(MISS, lonn, latn, PARAMS)                                         # write the last set of dates
if ARGS[2] == "v2.0_global_025_deg_geostrophic"  date = lins[i-1][1:8] * "00"
else                                             date = lins[i-1][1:10]  end
file = ARGS[2] * "/" * date * tail ; println(file)
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
end

for a = START:i - 1                                                           # echo the drifter data lines, but with analysis
  vals = split(lins[a])                                                       # velocity replacing drifter velocity
  lat = float(vals[2])
  lon = float(vals[3])
  if ARGS[2] == "v2.0_global_025_deg_geostrophic"  lon <    0 && (lon += 360) ; lon > 360 && (lon -= 360)
  else                                             lon < -180 && (lon += 360) ; lon > 180 && (lon -= 360)  end
  dellat, indlat = findmin(abs(lats - lat))
  dellon, indlon = findmin(abs(lons - lon))
  form = @sprintf("%s %11.5f %11.5f %s", lins[a][1:92], stor[indlon,indlat,UCUR], stor[indlon,indlat,VCUR], lins[a][118:end])
  write(fpb, form)
end
n += 1
print("wrote $n dates\n\n")

close(fpa)
close(fpb)
exit(0)


#=
2012090100   31.75000  164.00000      81961     52901  3 22889.00000   163.93700    31.82100    -0.29050     0.31430    -0.15673     0.07769     0.05344    -0.00772     0.05938    -0.00146    -0.06511     0.00914    -0.00765     0.02985    -0.10288     0.03854

const EKLO             = 1                              # identify indecies of all data variables
const EKHI             = 2
const GEOS             = 3
const TTLO             = 4
const TTHI             = 5
=#
