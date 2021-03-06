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
const POSTCAL          = true                           # include a bottom row of plots after adjustment
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_extra.ucur.got2000_obs.comb v2.0_global_025_deg_total_15m\n\n")
  exit(1)
end

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

gridbefaftnow = zeros(length(bound), length(bound), 2)                        # read the unadjusted before and after grids
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

gricbefaftnow = zeros(length(bound), length(bound), 2)                        # and the adjusted before and after grids
gricbefaftobs = zeros(length(bound), length(bound), 2)
gricbefaftone = zeros(length(bound), length(bound), 1)
gricnowobsone = zeros(length(bound), length(bound), 1)
fname = ARGS[1] * "." * ARGS[2] * ".extra.dau"
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
cutbefnow = zeros(length(bound)) ; cutaftnow = zeros(length(bound))
fname = ARGS[1] * "." * ARGS[2] * ".extra.sum"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  line = readline(fpa)
  (sumbefnow[a], conbefnow[a], sumaftnow[a], conaftnow[a],
   sumbefobs[a], conbefobs[a], sumaftobs[a], conaftobs[a],
   sumobsbef[a], conobsbef[a], sumobsaft[a], conobsaft[a],
   sumbefaft[a], conbefaft[a], sumaftbef[a], conaftbef[a],
   sumnowobs[a], connowobs[a], sumobsnow[a], conobsnow[a],
                 cutbefnow[a],               cutaftnow[a]) = float(split(line))
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
cuubefnow = zeros(length(bound)) ; cuuaftnow = zeros(length(bound))
fname = ARGS[1] * "." * ARGS[2] * ".extra.sun"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  line = readline(fpa)
  (sucbefnow[a], cocbefnow[a], sucaftnow[a], cocaftnow[a],
   sucbefobs[a], cocbefobs[a], sucaftobs[a], cocaftobs[a],
   sucobsbef[a], cocobsbef[a], sucobsaft[a], cocobsaft[a],
   sucbefaft[a], cocbefaft[a], sucaftbef[a], cocaftbef[a],
   sucnowobs[a], cocnowobs[a], sucobsnow[a], cocobsnow[a],
                 cuubefnow[a],               cuuaftnow[a]) = float(split(line))
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

fname = ARGS[1] * "." * ARGS[2] * ".extra.reg"                                # finally read the regression coefficient pairs
fpa = My.ouvre(fname, "r")
line = readline(fpa)
(intbefnow, slobefnow, intaftnow, sloaftnow,
 intbefobs, slobefobs, intaftobs, sloaftobs,
 intobsbef, sloobsbef, intobsaft, sloobsaft,
 intbefaft, slobefaft, intaftbef, sloaftbef,
 intnowobs, slonowobs, intobsnow, sloobsnow,
            covbefnow,            covaftnow,
            covbefobs,            covaftobs,
            covobsbef,            covobsaft,
            covbefaft,            covaftbef,
            covnowobs,            covobsnow) = float(split(line))
close(fpa)

fname = ARGS[1] * "." * ARGS[2] * ".extra.reh"                                # and the corresponding adjusted pairs
fpa = My.ouvre(fname, "r")
line = readline(fpa)
(incbefnow, slcbefnow, incaftnow, slcaftnow,
 incbefobs, slcbefobs, incaftobs, slcaftobs,
 incobsbef, slcobsbef, incobsaft, slcobsaft,
 incbefaft, slcbefaft, incaftbef, slcaftbef,
 incnowobs, slcnowobs, incobsnow, slcobsnow,
            cowbefnow,            cowaftnow,
            cowbefobs,            cowaftobs,
            cowobsbef,            cowobsaft,
            cowbefaft,            cowaftbef,
            cownowobs,            cowobsnow) = float(split(line))
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
#=
out = FramedPlot(xrange = (0,1), yrange = (0,1))
#setattr(out.x1, "draw_axis", false) ; setattr(out.x2, "draw_axis", false)
#setattr(out.y1, "draw_axis", false) ; setattr(out.y2, "draw_axis", false)
setattr(out.x1, "draw_ticks", false) ; setattr(out.x2, "draw_ticks", false)
setattr(out.y1, "draw_ticks", false) ; setattr(out.y2, "draw_ticks", false)
setattr(out.x1, "draw_ticklabels", false) ; setattr(out.x2, "draw_ticklabels", false)
setattr(out.y1, "draw_ticklabels", false) ; setattr(out.y2, "draw_ticklabels", false)
add(out, Points(0.5, 0.5, color="white"))
#x = linspace(0, 10)
##p1 = Points(0.5, 0.5, color="white")
#p2 = plot(x, sin(x), color="red")
#inset = PlotInset((0.1, 0.6), (0.7, 0.9), p2)
#add(out, inset)
#ppp = out
#if false
=#
POSTCAL && (ppp = Winston.Table(2,5) ; setattr(ppp, "cellpadding", -0.6))
POSTCAL || (ppp = Winston.Table(1,5)) #; setattr(ppp, "cellpadding", -0.5)    # and then create the unadjusted BEF/AFT
for z = 1:2                                                                   # vs NOW scatterplots, respectively (make
  (xpts, ypts, zpts) = point(bound, gridbefaftnow, z)                         # xpts and ypts refer to grid midpoints)
   xpts += 0.5 * step
   ypts += 0.5 * step

  z == 1 && (plotnam = "Forecast vs " * analysi)# * "\n" * varname)
  z == 2 && (plotnam =  "Revcast vs " * analysi)# * "\n" * varname)
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax), title = plotnam)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
# POSTCAL && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "lightgray"))
   z == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[1,z] = Winston.add(tmp)

  z == 1 && (inter = intbefnow ; slope = slobefnow)
  z == 2 && (inter = intaftnow ; slope = sloaftnow)
  tmp = Winston.Slope(slope, (0, inter), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
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
    tmp = Winston.Curve(bndbefnow, avgbefnow, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.Curve(bound, cutbefnow*100, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "a", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[1,z], tmp)
  else
    tmp = Winston.Curve(bndaftnow, avgaftnow, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.Curve(bound, cutaftnow*100, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "b", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[1,z], tmp)
    tmp = Winston.PlotLabel(0.20, 0.90, varname, "texthalign", "left", "size", 3.0)
          Winston.add(ppp[1,z], tmp)
  end
end

if POSTCAL
for z = 1:2                                                                   # add the adjusted BEF/AFT vs NOW scatterplots
  (xpts, ypts, zpts) = point(bound, gricbefaftnow, z)
   xpts += 0.5 * step
   ypts += 0.5 * step

  z == 1 && (plotnam = "Adj-forecast vs " * analysi)# * "\n" * varname)
  z == 2 && (plotnam =  "Adj-revcast vs " * analysi)# * "\n" * varname)
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax), title = plotnam)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
#           setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  z == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[2,z] = Winston.add(tmp)

  z == 1 && (inter = incbefnow ; slope = slcbefnow)
  z == 2 && (inter = incaftnow ; slope = slcaftnow)
  tmp = Winston.Slope(slope, (0, inter), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[2,z], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[2,z], tmp)

  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = DOTSIZE)
          Winston.add(ppp[2,z], tmp)
  end
  if z == 1
    tmp = Winston.Curve(bncbefnow, avcbefnow, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[2,z], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "c", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[2,z], tmp)
  else
    tmp = Winston.Curve(bncaftnow, avcaftnow, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[2,z], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "d", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[2,z], tmp)
  end
end
end

for z = 1:2                                                                   # add the unadjusted BEF/AFT vs OBS scatterplots
  (xpts, ypts, zpts) = point(bound, gridbefaftobs, z)
   xpts += 0.5 * step
   ypts += 0.5 * step

  z == 1 && (plotnam = "Forecast vs " * observr)# * "\n" * varname)
  z == 2 && (plotnam =  "Revcast vs " * observr)# * "\n" * varname)
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax), title = plotnam)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
# POSTCAL && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
   z == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[1,z+2] = Winston.add(tmp)

  z == 1 && (inter = intbefobs ; slope = slobefobs ; intes = intobsbef ; slopf = sloobsbef)
  z == 2 && (inter = intaftobs ; slope = sloaftobs ; intes = intobsaft ; slopf = sloobsaft)
  tmp = Winston.Slope(      slope, (0, inter), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[1,z+2], tmp)
  tmp = Winston.Slope(1.0 / slopf, (intes, 0), kind = "solid", "linewidth", LINWID, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[1,z+2], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[1,z+2], tmp)

  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = DOTSIZE)
          Winston.add(ppp[1,z+2], tmp)
    if z == 1
      tmp = Winston.PlotLabel(0.08, 1.00 - a * 0.07, "<span foreground=\"$(cols[length(cols) - a + 1])\">\\geq $(lims[length(cols) - a + 1])</span>", "texthalign", "left", "size", 3.0)
            Winston.add(ppp[1,z+2], tmp)
    end
  end
  if z == 1
    tmp = Winston.Curve(bndbefobs, avgbefobs, kind = "solid", "linewidth", LINWID) 
          Winston.add(ppp[1,z+2], tmp)
    tmp = Winston.Curve(avgobsbef, bndobsbef, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z+2], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "a", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[1,z+2], tmp)
  else
    tmp = Winston.Curve(bndaftobs, avgaftobs, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z+2], tmp)
    tmp = Winston.Curve(avgobsaft, bndobsaft, kind = "solid", "linewidth", LINWID)
          Winston.add(ppp[1,z+2], tmp)
    tmp = Winston.PlotLabel(0.85, 0.10, "b", "texthalign", "center", "size", 3.0)
          Winston.add(ppp[1,z+2], tmp)
    tmp = Winston.PlotLabel(0.20, 0.90, varname, "texthalign", "left", "size", 3.0)
          Winston.add(ppp[1,z+2], tmp)
  end
end

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

#end
#inset = PlotInset((0.1, 0.9), (0.1, 0.6), ppp)
#add(out, inset)

xyzzy = ARGS[1] * "." * ARGS[2] * ".extra.png"
print("writing $xyzzy\n")
POSTCAL && Winston.savefig(ppp, xyzzy, "width", 2100, "height", 1000)
POSTCAL || Winston.savefig(ppp, xyzzy, "width", 2100, "height",  500)
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

# setattr(tmp.x1, "tickdir",           -1) ; setattr(tmp.x2, "tickdir",           -1) ; setattr(tmp.y1, "tickdir",            1) ; setattr(tmp.y2, "tickdir",           -1)
# setattr(tmp.x1, "ticklabels_offset",  0) ; setattr(tmp.x2, "ticklabels_offset", -7) ; setattr(tmp.y1, "ticklabels_offset",  0) ; setattr(tmp.y2, "ticklabels_offset", -7)
# setattr(tmp.x1, "draw_subticks",  false) ; setattr(tmp.x2, "draw_subticks",   true) ; setattr(tmp.y1, "draw_subticks",  false) ; setattr(tmp.y1, "draw_subticks",   true)
#           setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "antiquewhite"))
# z == 1 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "black"))
# z == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "antiquewhite"))
#           setattr(tmp.x1, "draw_subticks", false)
# z == 1 && setattr(tmp.y1, "draw_subticks",  true)
# z == 2 && setattr(tmp.y1, "draw_subticks", false)
# if z == 2
# tmp = Winston.PlotLabel(0.50, 0.90, plotitle, "texthalign", "center", "size", 2.0)
#       Winston.add(ppp[1,z], tmp)
# end

  cols = ["red",  "blue", "green", "orange"]
  lims = [    1,      10,     100,     1000]
  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = 0.1)
          Winston.add(ppp[1,z], tmp)
    if z == 1
      tmp = Winston.Curve(bbndb, avgb, kind = "solid")
            Winston.add(ppp[1,z], tmp)
#     tmp = Winston.PlotLabel(0.08, 1.00 - a * 0.07, "<span foreground=\"$(cols[length(cols) - a + 1])\">\\geq $(lims[length(cols) - a + 1])</span>", "texthalign", "left", "size", 3.0)
#           Winston.add(ppp[1,z], tmp)
    end
    if z == 2
      tmp = Winston.Curve(bbnda, avga, kind = "solid")
            Winston.add(ppp[1,z], tmp)
    end
  end

  tmp = Winston.PlotLabel(0.50, 0.90, plotitle, "texthalign", "center", "size", 2.0)
        Winston.add(ppp[2,z], tmp)
=#
