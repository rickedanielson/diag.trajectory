#=
 = Loop through the timeseries and plot the corresponding forward and
 = backward extrapolated timeseries of all available variables, relative to the actual
 = (uninterpolated) values.  Note that BEF refers to an interpolation using analysis
 = data from before the extrapolation; AFT extrapolations use analysis data afterward.
 = Where one analysis is unavailable, all analyses are skipped - RD May 2016.
 =#

using My, Winston
const UCUR             = 1                              # indecies of all data variables
const VCUR             = 2
const PARAMS           = 2

const BEF              = 1                              # indecies of the source of extrapolations
const NOW              = 2
const AFT              = 3
const SRCS             = 3

const TIMS             = 3408                           # number in timeseries
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) v2.0_global_025_deg_geostrophic\n\n")
  exit(1)
end

dirs = ["v2.0_global_025_deg_ekman_15m", "v2.0_global_025_deg_ekman_hs", "v2.0_global_025_deg_geostrophic", "v2.0_global_025_deg_total_15m", "v2.0_global_025_deg_total_hs"]
dirn = length(dirs)

ucui = 0.002 ; ucus = collect(-5.0 : ucui : 5.0) ; ucun = zeros(length(ucus), length(ucus), length(dirs))
vcui = 0.002 ; vcus = collect(-5.0 : vcui : 5.0) ; vcun = zeros(length(vcus), length(vcus), length(dirs))

function restore(bound::Array{Float64,1}, grid::Array{Float64,3}, pname::UTF8String)
  fname = "extrapolated.histogr." * pname * ".dat"
  fpa = My.ouvre(fname, "r")
  for (a, vala) in enumerate(bound)
    for (b, valb) in enumerate(bound)
      line = readline(fpa)
      (grid[b,a,1], grid[b,a,2], grid[b,a,3], grid[b,a,4], grid[b,a,5]) = float(split(line))
    end
  end
  close(fpa)
end

restore(ucus, ucun, utf8("ucur"))
restore(vcus, vcun, utf8("vcur"))

function point(bound::Array{Float64,1}, grid::Array{Float64,3}, plotind::Int64)
  xpts = Array(Float64, 0)
  ypts = Array(Float64, 0)
  zpts = Array(Float64, 0)
  for (a, vala) in enumerate(bound)
    for (b, valb) in enumerate(bound)
      if grid[b,a,plotind] > 0
        push!(xpts, vala)
        push!(ypts, valb)
        push!(zpts, float(grid[b,a,plotind]))
      end
    end
  end
  return(xpts, ypts, zpts)
end

ARGS[1] == "v2.0_global_025_deg_ekman_15m"   && (plotind = 1 ; plotitle = "Ekman 15-m")
ARGS[1] == "v2.0_global_025_deg_ekman_hs"    && (plotind = 2 ; plotitle = "Ekman Surface")
ARGS[1] == "v2.0_global_025_deg_geostrophic" && (plotind = 3 ; plotitle = "Geostrophic")
ARGS[1] == "v2.0_global_025_deg_total_15m"   && (plotind = 4 ; plotitle = "Total 15-m")
ARGS[1] == "v2.0_global_025_deg_total_hs"    && (plotind = 5 ; plotitle = "Total Surface")

ppp = Winston.Table(1,2) ; setattr(ppp, "cellpadding", -0.5)                  # and then create the scatterplots
for z = 1:PARAMS
  z == UCUR && (varname =      "a) Zonal Current (ms^{-1})" ; (xpts, ypts, zpts) = point(ucus, ucun, plotind) ; tpos = (1,1) ; delt = ucui)
  z == VCUR && (varname = "b) Meridional Current (ms^{-1})" ; (xpts, ypts, zpts) = point(vcus, vcun, plotind) ; tpos = (1,2) ; delt = vcui)

  xpts += 0.5 * delt                                                          # make xpts and ypts refer to grid midpoints
  ypts += 0.5 * delt                                                          # and locate the plot limits
# xmin = minimum(xpts) - delt * 5
# xmax = maximum(xpts) + delt * 5
# ymin = minimum(ypts) - delt * 5
# ymax = maximum(ypts) + delt * 5
  z == UCUR && (xmin = -3.0 ; xmax = 3.0 ; ymin = -3.0 ; ymax = 3.0)
  z == VCUR && (xmin = -3.0 ; xmax = 3.0 ; ymin = -3.0 ; ymax = 3.0)

  cols = ["red",  "blue", "green", "orange", "black", "white"]
  lims = [    1,      10,     100,     1000,   10000,  100000]

  tmp = Winston.FramedPlot(title="$varname", xrange = (xmin,xmax), yrange = (ymin,ymax))
  ppp[tpos...] = Winston.add(tmp)

  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = 0.1)
          Winston.add(ppp[tpos...], tmp)
    if z == UCUR
      tmp = Winston.PlotLabel(0.05, 0.93, plotitle, "texthalign", "left", "size", 3.4)
            Winston.add(ppp[tpos...], tmp)
    end
    if z == VCUR
      tmp = Winston.PlotLabel(0.04, 1.05 - a * 0.05, "<span foreground=\"$(cols[length(cols) - a + 1])\">\\geq $(lims[length(cols) - a + 1])</span>", "texthalign", "left", "size", 2.0)
            Winston.add(ppp[tpos...], tmp)
    end
  end

  (intbef, slobef) = linreg(xpts, ypts, zpts)
  tmp = Winston.Slope(slobef, (0, intbef), kind = "solid", "linewidth", 5, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[tpos...], tmp)
  tmp = Winston.Slope(     1, (0,      0), kind = "solid")
        Winston.add(ppp[tpos...], tmp)
end

xyzzy = "extrapolated." * ARGS[1] * ".png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
exit(0)


#=
    ump = Array(Any, length(cols))
    ump[a] = Winston.Curve(specval[1:end,z], spectra[a,1:end,z], "color", parse(Winston.Colorant, cols[b]))
             style(ump[a], kind = kynd[b])
             setattr(ump[a], label = "($(specstr[a,z])) $(dirs[a])")
             Winston.add(ppp[tpos...], ump[a])
    b += 1
  end
  tmp = Winston.Legend(.23, .92, Any[ump[order[1]], ump[order[2]], ump[order[3]], ump[order[4]]])
        Winston.add(ppp[tpos...], tmp)
  tmp = Winston.Legend(.58, .92, Any[ump[order[5]], ump[order[6]], ump[order[7]], ump[order[8]]])
        Winston.add(ppp[tpos...], tmp)

          title="$varname Spectra (dB)", xlog = true,
          xlabel="Timescale (days)", xrange = (1/1000,1/2), yrange = (ymin,ymax))
          xlog = true, xrange = (1/1000,1/2), yrange = (ymin,ymax))
  setattr(tmp.x1, "ticks",          xposa) ; setattr(tmp.x2, "ticks",          xposb) ; setattr(tmp.y1, "ticks",          yposa)
  setattr(tmp.x1, "tickdir",            1) ; setattr(tmp.x2, "tickdir",           -1) ; setattr(tmp.y1, "tickdir",            1)
  setattr(tmp.x1, "ticklabels_offset",  0) ; setattr(tmp.x2, "ticklabels_offset", -7) ; setattr(tmp.y1, "ticklabels_offset",  0)
  setattr(tmp.x1, "ticklabels",     xlaba) ; setattr(tmp.x2, "ticklabels",     xlabb) ; setattr(tmp.y1, "ticklabels",     ylaba)
  setattr(tmp.x1, "draw_subticks",  false) ; setattr(tmp.x2, "draw_subticks",   true) ; setattr(tmp.y1, "draw_subticks",   true)
  tpos[1] <= 2 && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
  tpos[1] >= 2 && setattr(tmp.x2, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
  tpos[2] == 1 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "black"))
  tpos[2] == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
=#
