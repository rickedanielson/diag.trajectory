#=
 = Loop through all analyses and plot the binned
 = sums of all available variables - RD May 2016.
 =#

using My, Winston
const UCUR             = 1                              # indecies of all data variables
const VCUR             = 2
const PARAMS           = 2

const TIMS             = 3408                           # number in timeseries
const MISS             = -9999.0                        # generic missing value

if size(ARGS) != (0,) 
  print("\nUsage: jjj $(basename(@__FILE__))\n\n")
  exit(1)
end

dirs = ["v2.0_global_025_deg_ekman_15m", "v2.0_global_025_deg_ekman_hs", "v2.0_global_025_deg_geostrophic", "v2.0_global_025_deg_total_15m", "v2.0_global_025_deg_total_hs"]
dirz = ["15-m Ekman", "Surface Ekman", "Geostrophic", "15-m Total", "Surface Total"]
dirn = length(dirs)

ucui = 0.002 ; ucus = collect(-5.0 : ucui : 5.0) ; ucun = zeros(length(ucus), length(dirs))
vcui = 0.002 ; vcus = collect(-5.0 : vcui : 5.0) ; vcun = zeros(length(vcus), length(dirs))

function restore(bound::Array{Float64,1}, grid::Array{Float64,2}, pname::UTF8String)
  fname = "histogr." * pname * ".dat"
  fpa = My.ouvre(fname, "r")
  for (a, vala) in enumerate(bound)
    line = readline(fpa)
    (grid[a,1], grid[a,2], grid[a,3], grid[a,4], grid[a,5]) = float(split(line))
  end
  close(fpa)
end

restore(ucus, ucun, utf8("ucur"))
restore(vcus, vcun, utf8("vcur"))

ppp = Winston.Table(1,2) ; setattr(ppp, "cellpadding", -0.5)                  # and then create the plots
for z = 1:PARAMS
  z == UCUR && (varname =      "a) Zonal Current (ms^{-1})" ; bound = ucus ; grid = ucun ; tpos = (1,1) ; delt = ucui)
  z == VCUR && (varname = "b) Meridional Current (ms^{-1})" ; bound = vcus ; grid = vcun ; tpos = (1,2) ; delt = vcui)

  bound += 0.5 * delt                                                         # make bound refer to grid midpoints
  z == UCUR && (xmin = -0.6 ; xmax = 0.6 ; ymin = 0 ; ymax = 300000)          # and locate the plot limits
  z == VCUR && (xmin = -0.6 ; xmax = 0.6 ; ymin = 0 ; ymax = 300000)

  ump = Array(Any, dirn)
  cols = [  "red",    "red", "green",  "blue",   "blue"]
  kynd = ["solid", "dashed", "solid", "solid", "dashed"]

  tmp = Winston.FramedPlot(title="$varname", xrange = (xmin,xmax), yrange = (ymin,ymax))
  ppp[tpos...] = Winston.add(tmp)

  for a = 1:dirn
    ump[a] = Winston.Curve(bound, grid[:,a], "color", parse(Winston.Colorant, cols[a]))
             style(ump[a], kind = kynd[a])
             setattr(ump[a], label = dirz[a])
             Winston.add(ppp[tpos...], ump[a])
  end
  if z == UCUR
    tmp = Winston.Legend(.07, .90, Any[ump[1], ump[2], ump[3], ump[4], ump[5]])
          Winston.add(ppp[tpos...], tmp)
  end
end

xyzzy = "histogr_z.list.png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
exit(0)


# if ( ARGS[1] ==  "erainterim" || ARGS[1] ==       "hoaps" || ARGS[1] == "ifremerflux" ||
#      ARGS[1] ==       "merra" || ARGS[1] ==      "oaflux" || ARGS[1] ==     "seaflux" ||
#     (ARGS[1] ==        "cfsr" &&              z != LHFX)   ||
#     (ARGS[1] ==      "jofuro" && z != SSTT && z != AIRT))
