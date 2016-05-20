#=
 = Loop through all analyses and visualize the availability of all
 = variables at the position of interest - RD April 2016.
 =#

using My, Winston
const UCUR             = 1                              # identify indecies of all data variables
const VCUR             = 2
const PARAMS           = 2

const UCUP             = 10                             # position of data variables on input data lines
const VCUP             = 11

const LOTS             = 1200                           # width of each analysis timeseries
const DAYS             = 3408                           # number of 6-h periods between 2012-09-01-00 and 2014-12-31-18
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) ....45.000...-45.500\n\n")
  exit(1)
end

dirs = ["v2.0_global_025_deg_ekman_15m", "v2.0_global_025_deg_ekman_hs", "v2.0_global_025_deg_geostrophic", "v2.0_global_025_deg_total_15m", "v2.0_global_025_deg_total_hs", "insitu"]
dirz = [                    "ekman_15m",                     "ekman_hs",                     "geostrophic",                     "total_15m",                     "total_hs", "insitu"]
data = Array(Float64, LOTS, DAYS, PARAMS)

for (a, dir) in enumerate(dirs)                                               # and read the current components
fila = "$dir/$dir$(ARGS[1])"
  fpa = My.ouvre(fila, "r")
  lines = readlines(fpa)
  close(fpa)
  for (b, linb) in enumerate(lines)
    vals = float(split(linb))
    data[a,b,UCUR] = -333 < vals[UCUP] < 333 ? 1.0 : 0.0
    data[a,b,VCUR] = -333 < vals[UCUP] < 333 ? 1.0 : 0.0
  end
end

for a = 1:PARAMS, b = 1:DAYS                                                  # map each dir to a thick line
  data[1001:1200,b,a] = data[6,b,a]
  data[ 801:1000,b,a] = data[5,b,a]
  data[ 601: 800,b,a] = data[4,b,a]
  data[ 401: 600,b,a] = data[3,b,a]
  data[ 201: 400,b,a] = data[2,b,a]
  data[   1: 200,b,a] = data[1,b,a]
end

xlab = Array(UTF8String, 0)                                                   # initialize the date label strings
xpos = Array(     Int64, 0)                                                   # (first date and then first day of
date = "2012-09-01-00"                                                        #  each subsequent year)
push!(xlab, "1 Oct\n$(date[1:4])") ; push!(xpos, 1)
for a = 2:DAYS
  date = dateadd(date, 6, "hr")
  if date[6:13] == "01-01-00" # && date[1:4] != "2000"
    push!(xlab, date[1:4]) ; push!(xpos, a)
  end
end

Colors.colormap("Blues", 3)
ppp  = Winston.Table(1,2) ; setattr(ppp, "cellpadding", -4.0)
for z = 1:PARAMS
  z == UCUR && (varname =      "a) Zonal Wind Component" ; tpos = (1,1))
  z == VCUR && (varname = "b) Meridional Wind Component" ; tpos = (1,2))

  tmp = Winston.imagesc(data[:,:,z])
  setattr(tmp,    "title",               varname) ; setattr(tmp,    "aspect_ratio",           0.35)
  setattr(tmp.x1, "ticks",                  xpos)#; setattr(tmp.x2, "draw_nothing",           true)
  setattr(tmp.x1, "tickdir",                   1)#; setattr(tmp.x2, "tickdir",                   1)
  setattr(tmp.x1, "ticklabels",             xlab) ; setattr(tmp.x2, "draw_ticks",            false)
  setattr(tmp.x1, "draw_subticks",         false) ; setattr(tmp.x2, "draw_subticks",         false)
  setattr(tmp.y1, "ticks", collect(100:200:1100))#; setattr(tmp.y2, "draw_nothing",           true)
  setattr(tmp.y1, "tickdir",                   1)#; setattr(tmp.y2, "tickdir",                   1)
  setattr(tmp.y1, "ticklabels",             dirz) ; setattr(tmp.y2, "draw_ticks",            false)
  setattr(tmp.y1, "draw_subticks",         false) ; setattr(tmp.y2, "draw_subticks",         false)
  ppp[tpos...] = Winston.add(tmp)

# @show getattr(tmp, :aspect_ratio)
  tpos[1] <= 2 && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:fontsize => 3, :color => "transparent"))
  tpos[2] == 1 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:fontsize => 3, :color => "black"))
  tpos[2] == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:fontsize => 3, :color => "transparent"))
end

xyzzy = "plot.avail$(ARGS[1]).png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
exit(0)
