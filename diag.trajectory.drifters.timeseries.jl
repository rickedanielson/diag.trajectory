#=
 = Loop through all analyses and extract variables of interest at the position
 = of a set of buoy observations.  Include a conversion of units, as required
 = - RD September 2015, February, March 2016.
 =#

using My, NetCDF
const MISS             = -9999.0                        # generic missing value
const SHFX             = 1                              # identify indecies of all data variables
const LHFX             = 2
const IDEN             = 3
const DATE             = 4
const LATT             = 5
const LONN             = 6
const SLPP             = 7
const WDIR             = 8
const WSPD             = 9
const UWND             = 10
const VWND             = 11
const AIRT             = 12
const DPTT             = 13
const SSTT             = 14
const SHUM             = 15
const BLAY             = 16
const HGTU             = 17
const HGTT             = 18
const PARAMS           = 18

if size(ARGS) != (2,)
  print("\nUsage: jjj $(basename(@__FILE__)) all.flux.daily_2.0_locate.calib.sortac merra\n\n")
  exit(1)
end

buoy = Set(Array(Tuple{Float64, Float64}, 0))                                 # initialize a grid mask that is
fpa = My.ouvre(ARGS[1], "r") ; lines = readlines(fpa) ; close(fpa)            # true wherever buoys are located
for line in lines  vals = split(line) ; push!(buoy, (float(vals[2]), float(vals[1])))  end

lats =  -90.0:0.25:89.75                                                      # locs[   1,  1] = (-180.00,-90.00)
lons = -180.0:0.25:179.75                                                     # locs[   1,720] = (-180.00, 89.75)
locs = [   (x, y)        for x=lons, y=lats]                                  # locs[1440,  1] = ( 179.75,-90.00)
mask = [in((x, y), buoy) for x=lons, y=lats]                                  # locs[1440,720] = ( 179.75, 89.75)
fend = locs[mask .== true] ; (locn,) = size(fend)                             # and define all filename endings

vars = ["surface_upward_sensible_heat_flux", "surface_upward_latent_heat_flux", "", "", "", "", "", "", "wind_speed",
        "", "", "air_temperature", "", "sea_surface_temperature", "air_surface_specific_humidity", "", "", ""]
tail = "-OHF-L4-global_daily_0.25x0.25-v0.7-f01.0.nc"

ARGS[2] ==        "cfsr" && (stem =        "cfsr-" ; vars[LHFX] = "" ; vars[AIRT] = "2t")
ARGS[2] ==  "erainterim" && (stem =  "erainterim-")
ARGS[2] ==       "hoaps" && (stem =       "hoaps-")
ARGS[2] == "ifremerflux" && (stem = "ifremerflux-")
ARGS[2] ==      "jofuro" && (stem =      "jofuro-" ; vars[AIRT] = "" ; vars[SSTT] = "")
ARGS[2] ==       "merra" && (stem =       "merra-")
ARGS[2] ==      "oaflux" && (stem =      "oaflux-")                           # define the available variables
ARGS[2] ==     "seaflux" && (stem =     "seaflux-")                           # and the filename stems/ending

fpn = Array(IOStream, 0)                                                      # then open a set of output files
for a = 1:locn                                                                # at all buoy locations
  (lon, lat) = fend[a]
  tmp = @sprintf("%9.3f.%9.3f", lat, lon)
  tmq = replace(tmp, " ", ".")
  tmr = "$(ARGS[2])/$(ARGS[2]).$tmq"
  fpa = My.ouvre(tmr, "w")
  push!(fpn, fpa)
end

date = "19991001120000"                                                       # and fill the output files (all
while parse(Int, date) < 20100101000000                                       # locs and vars) one date at a time
  stor = fill(MISS, locn, PARAMS)                                             # (careful: try introduces new scope)
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
    if ARGS[2] == "cfsr" && 200.0 < stor[a,AIRT] < 350.0  stor[a,AIRT] -= 273.15  end
    if                      200.0 < stor[a,SSTT] < 350.0  stor[a,SSTT] -= 273.15  end
    form = @sprintf("%8.2f %8.2f %9s %14s %7.3f %8.3f %8.2f %8.3f %8.3f %8.3f %8.3f %8.2f %8.2f %8.2f %8.3f %8.2f %8.2f %8.2f\n",
      stor[a,SHFX], stor[a,LHFX], "0000", date[1:10], lat, lon, stor[a,SLPP], stor[a,WDIR], stor[a,WSPD], stor[a,UWND], stor[a,VWND],
      stor[a,AIRT], stor[a,DPTT], stor[a,SSTT], stor[a,SHUM], stor[a,BLAY], stor[a,HGTU], stor[a,HGTT])
    write(fpn[a], form)
  end
  date = My.dateadd(date, 1, "dy")
end

print("closing...\n")
for a = 1:locn  close(fpn[a])  end
exit(0)


#=
        try    nc = NetCDF.open(file, mode=NC_NOWRITE, readdimvar=false)          end
        try  dat3 = NetCDF.readvar(nc, vars[a], start=[1,1,1], count=[-1,-1,-1])  end
        try         NetCDF.close(nc)                                              end

form = @sprintf("%8.2f %8.2f %9s %14s %7.3f %8.3f %8.2f %8.3f %8.3f %8.3f %8.3f %8.2f %8.2f %8.2f %8.3f %8.2f %8.2f %8.2f\n",
        shf, lhf, iden, date, val[:LAT], val[:LON], val[:SLP], val[:D], val[:W], uwnd, vwnd, val[:AT], val[:DPT], val[:SST], shum,
        blay, hgtu, hgtt

    if vars[3] == "dew_point_temperature"
    end
    tmp = @sprintf("%s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n", date[1:10],
      stor[a,1], stor[a,2], stor[a,3], stor[a,4], stor[a,5], stor[a,6])
    write(fpn[a], tmp)
=#
