#=
 = Loop through the available collocations and plot the corresponding forward and
 = backward extrapolations relative to the actual (uninterpolated) values.  Note that
 = BEF refers to an extrapolation using analysis data from before the target value and
 = AFT refers to an extrapolation using data from afterward - RD May, November 2016.
 =#

using My, Winston
const ODAT             = 1                              # identify indecies of the input data:
const OLAT             = 2                              # date/lat/lon on the collocation grid
const OLON             = 3
const OCUR             = 4                              # then five buoy parameters

const MINSUM           = 10                             # minimum number of samples for an average
const LINWID           = 3                              # plotting line width
const DOTSIZE          = 0.2                            # size of scatterplot dots
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_extra.ucur.got2000_obs.comb v2.0_global_025_deg_total_15m\n\n")
  exit(1)
end
contains(ARGS[1], ".ucur") && (ARGS222 = replace(ARGS[1], "ucur", "vcur"))
contains(ARGS[1], ".vcur") && (ARGS222 = replace(ARGS[1], "vcur", "ucur"))

analysi = "GlobCurrent"
observr = "Drifter"
varname = ""
if contains(ARGS[1], "ucur")                                                  # define the plot dimensions and type
  varname = "U (ms^{-1})"
  xmin = -1.495 ; xmax = 1.495 ; ymin = -1.495 ; ymax = 1.495
elseif contains(ARGS[1], "vcur")
  varname = "V (ms^{-1})"
  xmin = -1.495 ; xmax = 1.495 ; ymin = -1.495 ; ymax = 1.495
elseif contains(ARGS[1], "wcur")
  varname = "Speed (ms^{-1})"
  xmin = -0.1 ; xmax = 2.0 ; ymin = -0.1 ; ymax = 2.0
end
step = 0.02 ; bound = collect(-4.0:step:4.0)
cols = ["red",  "blue", "green", "orange"]
lims = [    1,      10,     100,     1000]

gridbefaftnow = zeros(length(bound), length(bound), 2)                        # read the unadjusted zonal before and after grids
gridbefaftobs = zeros(length(bound), length(bound), 2)
gridbefaftone = zeros(length(bound), length(bound), 1)
gridnowobsone = zeros(length(bound), length(bound), 1)
fname = ARGS[1] * "." * ARGS[2] * ".extra.dat"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  for (b, valb) in enumerate(bound)
    line = readline(fpa)
    (gridbefaftnow[b,a,1], gridbefaftnow[b,a,2], gridbefaftobs[b,a,1], gridbefaftobs[b,a,2], gridbefaftone[b,a,1], gridnowobsone[b,a,1]) = float(split(line))
  end
end
close(fpa)

gricbefaftnow = zeros(length(bound), length(bound), 2)                        # and the unadjusted meridional before and after grids
gricbefaftobs = zeros(length(bound), length(bound), 2)
gricbefaftone = zeros(length(bound), length(bound), 1)
gricnowobsone = zeros(length(bound), length(bound), 1)
fname = ARGS222 * "." * ARGS[2] * ".extra.dat"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  for (b, valb) in enumerate(bound)
    line = readline(fpa)
    (gricbefaftnow[b,a,1], gricbefaftnow[b,a,2], gricbefaftobs[b,a,1], gricbefaftobs[b,a,2], gricbefaftone[b,a,1], gricnowobsone[b,a,1]) = float(split(line))
  end
end
close(fpa)

sumbefnow = zeros(length(bound)) ; conbefnow = zeros(length(bound))           # as well as the corresponding count and sum
sumaftnow = zeros(length(bound)) ; conaftnow = zeros(length(bound))           # of extrapolation values in each TOTN interval
sumbefobs = zeros(length(bound)) ; conbefobs = zeros(length(bound))
sumnowobs = zeros(length(bound)) ; connowobs = zeros(length(bound))
sumaftobs = zeros(length(bound)) ; conaftobs = zeros(length(bound))
sumobsbef = zeros(length(bound)) ; conobsbef = zeros(length(bound))
sumobsnow = zeros(length(bound)) ; conobsnow = zeros(length(bound))
sumobsaft = zeros(length(bound)) ; conobsaft = zeros(length(bound))
sumbefaft = zeros(length(bound)) ; conbefaft = zeros(length(bound))
sumaftbef = zeros(length(bound)) ; conaftbef = zeros(length(bound))
fname = ARGS[1] * "." * ARGS[2] * ".extra.sum"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  line = readline(fpa)
  (sumbefnow[a], conbefnow[a], sumaftnow[a], conaftnow[a],
   sumbefobs[a], conbefobs[a], sumaftobs[a], conaftobs[a],
   sumobsbef[a], conobsbef[a], sumobsaft[a], conobsaft[a],
   sumbefaft[a], conbefaft[a], sumaftbef[a], conaftbef[a],
   sumnowobs[a], connowobs[a], sumobsnow[a], conobsnow[a]) = float(split(line))
end
close(fpa)

sucbefnow = zeros(length(bound)) ; cocbefnow = zeros(length(bound))           # and the corresponding adjusted values
sucaftnow = zeros(length(bound)) ; cocaftnow = zeros(length(bound))
sucbefobs = zeros(length(bound)) ; cocbefobs = zeros(length(bound))
sucnowobs = zeros(length(bound)) ; cocnowobs = zeros(length(bound))
sucaftobs = zeros(length(bound)) ; cocaftobs = zeros(length(bound))
sucobsbef = zeros(length(bound)) ; cocobsbef = zeros(length(bound))
sucobsnow = zeros(length(bound)) ; cocobsnow = zeros(length(bound))
sucobsaft = zeros(length(bound)) ; cocobsaft = zeros(length(bound))
sucbefaft = zeros(length(bound)) ; cocbefaft = zeros(length(bound)) 
sucaftbef = zeros(length(bound)) ; cocaftbef = zeros(length(bound))
fname = ARGS222 * "." * ARGS[2] * ".extra.sum"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  line = readline(fpa)
  (sucbefnow[a], cocbefnow[a], sucaftnow[a], cocaftnow[a],
   sucbefobs[a], cocbefobs[a], sucaftobs[a], cocaftobs[a],
   sucobsbef[a], cocobsbef[a], sucobsaft[a], cocobsaft[a],
   sucbefaft[a], cocbefaft[a], sucaftbef[a], cocaftbef[a],
   sucnowobs[a], cocnowobs[a], sucobsnow[a], cocobsnow[a]) = float(split(line))
end
close(fpa)

function smooth(bnda::Array{Float64,1}, suma::Array{Float64,1}, avga::Array{Float64,1})
  suml = length(suma) ; sumb = zeros(suml) ; avgb = zeros(suml)
  for a = 1:suml
    if a == 1 || a == suml
      sumb[a] = -1
    else
      sumb[a] = suma[a-1] + suma[a] + suma[a+1]
      avgb[a] = avga[a-1] + avga[a] + avga[a+1]
    end
  end
  mask = sumb .>= MINSUM ; avgc = avgb[mask] ./ sumb[mask] ; bndc = bnda[mask]
  return(avgc, bndc)
end
(avgbefnow, bndbefnow) = smooth(bound, sumbefnow, conbefnow)                  # and allow a smoothing of these bin averages
(avgaftnow, bndaftnow) = smooth(bound, sumaftnow, conaftnow)                  # and the corresponding adjusted bin averages
(avgbefobs, bndbefobs) = smooth(bound, sumbefobs, conbefobs)
(avgnowobs, bndnowobs) = smooth(bound, sumnowobs, connowobs)
(avgaftobs, bndaftobs) = smooth(bound, sumaftobs, conaftobs)
(avgobsbef, bndobsbef) = smooth(bound, sumobsbef, conobsbef)
(avgobsnow, bndobsnow) = smooth(bound, sumobsnow, conobsnow)
(avgobsaft, bndobsaft) = smooth(bound, sumobsaft, conobsaft)
(avgbefaft, bndbefaft) = smooth(bound, sumbefaft, conbefaft)
(avgaftbef, bndaftbef) = smooth(bound, sumaftbef, conaftbef)

(avcbefnow, bncbefnow) = smooth(bound, sucbefnow, cocbefnow)
(avcaftnow, bncaftnow) = smooth(bound, sucaftnow, cocaftnow)
(avcbefobs, bncbefobs) = smooth(bound, sucbefobs, cocbefobs)
(avcnowobs, bncnowobs) = smooth(bound, sucnowobs, cocnowobs)
(avcaftobs, bncaftobs) = smooth(bound, sucaftobs, cocaftobs)
(avcobsbef, bncobsbef) = smooth(bound, sucobsbef, cocobsbef)
(avcobsnow, bncobsnow) = smooth(bound, sucobsnow, cocobsnow)
(avcobsaft, bncobsaft) = smooth(bound, sucobsaft, cocobsaft)
(avcbefaft, bncbefaft) = smooth(bound, sucbefaft, cocbefaft)
(avcaftbef, bncaftbef) = smooth(bound, sucaftbef, cocaftbef)

fname = ARGS[1] * "." * ARGS[2] * ".extra.reg"                                # finally read the zonal regression coefficient pairs
fpa = My.ouvre(fname, "r")
line = readline(fpa)
(intbefnow, slobefnow, intaftnow, sloaftnow,
 intbefobs, slobefobs, intaftobs, sloaftobs,
 intobsbef, sloobsbef, intobsaft, sloobsaft,
 intbefaft, slobefaft, intaftbef, sloaftbef,
 intnowobs, slonowobs, intobsnow, sloobsnow) = float(split(line))
close(fpa)

fname = ARGS222 * "." * ARGS[2] * ".extra.reg"                                # and the corresponding meridional pairs
fpa = My.ouvre(fname, "r")
line = readline(fpa)
(incbefnow, slcbefnow, incaftnow, slcaftnow,
 incbefobs, slcbefobs, incaftobs, slcaftobs,
 incobsbef, slcobsbef, incobsaft, slcobsaft,
 incbefaft, slcbefaft, incaftbef, slcaftbef,
 incnowobs, slcnowobs, incobsnow, slcobsnow) = float(split(line))
close(fpa)

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

ppp = Winston.Table(1,2) ; setattr(ppp, "cellpadding", -0.6)
for z = 1:2                                                                   # then plot NOW vs OBS scatterplots for
  z == 1 && ((xpts, ypts, zpts) = point(bound, gridnowobsone, 1))             # the zonal and meridional components
  z == 2 && ((xpts, ypts, zpts) = point(bound, gricnowobsone, 1))
  xpts += 0.5 * step
  ypts += 0.5 * step

  z == 1 && (plotnam = "Zonal GlobCurrent vs Drifter (ms^{-1})")
  z == 2 && (plotnam = "Meridional GlobCurrent vs Drifter (ms^{-1})")
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax), title = plotnam)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
# z == 1 && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  z == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[1,z] = Winston.add(tmp)

  z == 1 && (inter = intnowobs ; slope = slonowobs ; intes = intobsnow ; slopf = sloobsnow)
  z == 2 && (inter = incnowobs ; slope = slcnowobs ; intes = incobsnow ; slopf = slcobsnow)
  tmp = Winston.Slope(      slope, (0, inter), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[1,z], tmp)
  tmp = Winston.Slope(1.0 / slopf, (intes, 0), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[1,z], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[1,z], tmp)

  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = DOTSIZE)
          Winston.add(ppp[1,z], tmp)
    if z == 1
      tmp = Winston.PlotLabel(0.08, 1.00 - a * 0.07, "<span foreground=\"$(cols[length(cols) - a + 1])\">\\geq $(lims[length(cols) - a + 1])</span>", "texthalign", "left", "size", 3.0)
            Winston.add(ppp[1,z], tmp)
    end
  end
  if z == 1
    tmp = Winston.Curve(bndnowobs, avgnowobs, kind = "solid", "linewidth", LINWID) 
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.Curve(avgobsnow, bndobsnow, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.PlotLabel(0.80, 0.10, "a", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[1,z], tmp)
  else
    tmp = Winston.Curve(bncnowobs, avcnowobs, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.Curve(avcobsnow, bncobsnow, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.PlotLabel(0.80, 0.10, "b", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[1,z], tmp)
  end
end

xyzzy = ARGS[1] * "." * ARGS[2] * ".nowobs.png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 600, "height", 400)
exit(0)

#=
for z = 1:1                                                                   # add the unadjusted BEF vs AFT scatterplot
  (xpts, ypts, zpts) = point(bound, gridbefaftone, z)
   xpts += 0.5 * step
   ypts += 0.5 * step

  z == 1 && (plotnam = "Forecast vs Revcast")
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax), title = plotnam)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
# POSTCAL && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
             setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[1,z+4] = Winston.add(tmp)

  z == 1 && (inter = intbefaft ; slope = slobefaft ; intes = intaftbef ; slopf = sloaftbef)
  tmp = Winston.Slope(      slope, (0, inter), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[1,z+4], tmp)
  tmp = Winston.Slope(1.0 / slopf, (intes, 0), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[1,z+4], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[1,z+4], tmp)

  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = DOTSIZE)
          Winston.add(ppp[1,z+4], tmp)
  end
  if z == 1
    tmp = Winston.Curve(bndbefaft, avgbefaft, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z+4], tmp)
    tmp = Winston.Curve(avgaftbef, bndaftbef, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z+4], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "c", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[1,z+4], tmp)
  end
end

if POSTCAL
for z = 1:2                                                                   # add the adjusted BEF/AFT vs OBS scatterplots
  (xpts, ypts, zpts) = point(bound, gricbefaftobs, z)
   xpts += 0.5 * step
   ypts += 0.5 * step

  z == 1 && (plotnam = "Adj-forecast vs " * observr)# * "\n" * varname)
  z == 2 && (plotnam =  "Adj-revcast vs " * observr)# * "\n" * varname)
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax), title = plotnam)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
# z == 1 && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  z == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[2,z+2] = Winston.add(tmp)

  z == 1 && (inter = incbefobs ; slope = slcbefobs ; intes = incobsbef ; slopf = slcobsbef)
  z == 2 && (inter = incaftobs ; slope = slcaftobs ; intes = incobsaft ; slopf = slcobsaft)
  tmp = Winston.Slope(      slope, (0, inter), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[2,z+2], tmp)
  tmp = Winston.Slope(1.0 / slopf, (intes, 0), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[2,z+2], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[2,z+2], tmp)

  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = DOTSIZE)
          Winston.add(ppp[2,z+2], tmp)
  end
  if z == 1
    tmp = Winston.Curve(bncbefobs, avcbefobs, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[2,z+2], tmp)
    tmp = Winston.Curve(avcobsbef, bncobsbef, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[2,z+2], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "d", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[2,z+2], tmp)
  else
    tmp = Winston.Curve(bncaftobs, avcaftobs, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[2,z+2], tmp)
    tmp = Winston.Curve(avcobsaft, bncobsaft, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[2,z+2], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "e", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[2,z+2], tmp)
  end
end

for z = 1:1                                                                   # add the adjusted BEF vs AFT scatterplot
  (xpts, ypts, zpts) = point(bound, gricbefaftone, z)
   xpts += 0.5 * step
   ypts += 0.5 * step

  z == 1 && (plotnam = "Adj-forecast vs Adj-revcast")# vs " * observr)# * "\n" * varname)
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax), title = plotnam)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
# z == 1 && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
            setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[2,z+4] = Winston.add(tmp)

  z == 1 && (inter = incbefaft ; slope = slcbefaft ; intes = incaftbef ; slopf = slcaftbef)
  tmp = Winston.Slope(      slope, (0, inter), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[2,z+4], tmp)
  tmp = Winston.Slope(1.0 / slopf, (intes, 0), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[2,z+4], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[2,z+4], tmp)

  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = DOTSIZE)
          Winston.add(ppp[2,z+4], tmp)
  end
  if z == 1
    tmp = Winston.Curve(bncbefaft, avcbefaft, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[2,z+4], tmp)
    tmp = Winston.Curve(avcaftbef, bncaftbef, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[2,z+4], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "f", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[2,z+4], tmp)
  end
end
end

xyzzy = ARGS[1] * "." * ARGS[2] * ".nowobs.png"
print("writing $xyzzy\n")
POSTCAL && Winston.savefig(ppp, xyzzy, "width", 2100, "height", 1000)
POSTCAL || Winston.savefig(ppp, xyzzy, "width", 2100, "height",  500)
exit(0)
=#
