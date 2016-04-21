#=
 = Loop through all analyses and extract variables of interest at
 = the position of a set of drifter observations (accommodating the
 = CLS geostrophic grid difference) - RD April 2016.
 =#

using My, NetCDF
const MISS             = -9999.0                        # generic missing value
const UCUR             = 1                              # identify indecies of all data variables
const VCUR             = 2
const PARAMS           = 2

if size(ARGS) != (2,)
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.sort v2.0_global_025_deg_ekman_15m\n\n")
  exit(1)
end

buoy = Set(Array(Tuple{Float64, Float64}, 0))                                 # initialize a grid mask that is
fpa = My.ouvre(ARGS[1], "r") ; lines = readlines(fpa) ; close(fpa)            # true wherever buoys are located
for line in lines  vals = split(line) ; push!(buoy, (float(vals[2]), float(vals[1])))  end

lats =  -79.875:0.25:79.875  ; latn = length(lats)
lons = -179.875:0.25:179.875 ; lonn = length(lons)
if ARGS[2] == "v2.0_global_025_deg_geostrophic"
  lats = -89.875:0.25:89.875  ; latn = length(lats)                           # (for CLS geostrophy)
  lons =   0.125:0.25:359.875 ; lonn = length(lons)                           # locs[   1,  1] = (-180.00,-90.00)
end                                                                           # locs[   1,720] = (-180.00, 89.75)
locs = [   (x, y)        for x=lons, y=lats]                                  # locs[1440,  1] = ( 179.75,-90.00)
mask = [in((x, y), buoy) for x=lons, y=lats]                                  # locs[1440,720] = ( 179.75, 89.75)
fend = locs[mask .== true] ; (locn,) = size(fend)                             # and define all filename endings

vars = Array(UTF8String, 2)
ARGS[2] == "v2.0_global_025_deg_ekman_15m"   && (vars[UCUR] =       "eastward_ekman_current_velocity" ; vars[VCUR] =       "northward_ekman_current_velocity")
ARGS[2] == "v2.0_global_025_deg_ekman_hs"    && (vars[UCUR] =       "eastward_ekman_current_velocity" ; vars[VCUR] =       "northward_ekman_current_velocity")
ARGS[2] == "v2.0_global_025_deg_geostrophic" && (vars[UCUR] = "eastward_geostrophic_current_velocity" ; vars[VCUR] = "northward_geostrophic_current_velocity")
ARGS[2] == "v2.0_global_025_deg_total_15m"   && (vars[UCUR] =    "eastward_eulerian_current_velocity" ; vars[VCUR] =    "northward_eulerian_current_velocity")
ARGS[2] == "v2.0_global_025_deg_total_hs"    && (vars[UCUR] =    "eastward_eulerian_current_velocity" ; vars[VCUR] =    "northward_eulerian_current_velocity")

ARGS[2] == "v2.0_global_025_deg_ekman_15m"   && (tail =   "0000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v02.0-fv01.0.nc")
ARGS[2] == "v2.0_global_025_deg_ekman_hs"    && (tail =   "0000-GLOBCURRENT-L4-CURekm_hs-ERAWS_EEM-v02.0-fv01.0.nc")
ARGS[2] == "v2.0_global_025_deg_geostrophic" && (tail = "000000-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v02.0-fv01.0.nc")
ARGS[2] == "v2.0_global_025_deg_total_15m"   && (tail =   "0000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v02.0-fv01.0.nc")
ARGS[2] == "v2.0_global_025_deg_total_hs"    && (tail =   "0000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc")

fpn = Array(IOStream, 0)                                                      # open a set of output files
for a = 1:locn                                                                # at all buoy locations
  (lon, lat) = fend[a]
  tmp = @sprintf("%9.3f.%9.3f", lat, lon)
  tmq = replace(tmp, " ", ".")
  tmr = "$(ARGS[2])/$(ARGS[2]).$tmq"
  fpa = My.ouvre(tmr, "w")
  push!(fpn, fpa)
end

date = "2012090100"                                                           # and fill the output files (all
while parse(Int, date) < 2015010100                                           # locs and vars) one date at a time
  stor = fill(MISS, locn, PARAMS)                                             # (careful: try introduces new scope)
  if ARGS[2] == "v2.0_global_025_deg_geostrophic"  file = ARGS[2] * "/" * date[1:8] * tail
  else                                             file = ARGS[2] * "/" * date      * tail  end
  println(file)

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
        dat2 = dat3[:,:,1]
        stor[:,a] = dat2[mask .== true]
      end
    end
  end

  for a = 1:locn
    (lon, lat) = fend[a]
    stor[a,:] = map(x -> begin
          if x > -MISS return MISS
      elseif x <  MISS return MISS
      else return x
      end
    end, stor[a,:])
    form = @sprintf("%10s %10.5f %10.5f   99999999     99999  9 -9999.00000 -9999.00000 -9999.00000 %11.5f %11.5f\n",
      date, lat, lon, stor[a,UCUR], stor[a,VCUR])
    write(fpn[a], form)
  end
  date = My.dateadd(date, 6, "hr")
end

print("closing...\n")
for a = 1:locn  close(fpn[a])  end
exit(0)
