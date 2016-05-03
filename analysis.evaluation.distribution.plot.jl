#=
 = Plot histograms of all parameters - RD April 2016.
 =#

using My, Winston
const SHFX             = 3                              # identify indecies of all data variables
const LHFX             = 2
const WSPD             = 6
const SHUM             = 4
const SSTT             = 5
const AIRT             = 1
const PARAMS           = 6

const TIMES            = 2556                           # number of days during 2001-2007
const MISS             = -9999.0                        # generic missing value

if size(ARGS) != (6,)
  print("\nUsage: jjj $(basename(@__FILE__)) all.flux.daily.locate_2.0_calib.????.got2000\n\n")
  exit(1)
end

pars = ["AIRT", "LHFX", "SHFX", "SHUM", "SSTT", "WSPD"]
dirs = ["cfsr", "erainterim", "hoaps", "ifremerflux", "jofuro", "merra", "oaflux", "seaflux"]
dirn = length(dirs)

ppp     = Winston.Table(3,2) ; setattr(ppp, "cellpadding", -0.5)
#specval = Array(Float64,          SPECS, PARAMS)
#spectra = Array(Float64,    dirn, SPECS, PARAMS)
#specvar = zeros(Float64,    dirn,        PARAMS)
#specstr = Array(UTF8String, dirn,        PARAMS)

for z = 1:PARAMS
  dataval = Array(Float64, 0, 0)
  fpa = My.ouvre(ARGS[z], "r") ; lines = readlines(fpa) ; close(fpa)
  for (a, line) in enumerate(lines)
    vals = float(split(line))
    specval[a,z] = vals[1]
    for b = 1:dirn
      if vals[b+1] == MISS
        spectra[b,a,z] = MISS
      else
        spectra[b,a,z] = 10.0 * log10(vals[b+1])
        if vals[1] > 1 / VARSUM
          specvar[b,z] += vals[b+1]
        end
      end
    end
  end
  order = sortperm(specvar[:,z], rev=true)

  for a in order
    contains(ARGS[z], "shfx") && (specstr[a,z] = @sprintf("%.0lf", specvar[a,z]))
    contains(ARGS[z], "lhfx") && (specstr[a,z] = @sprintf("%.0lf", specvar[a,z]))
    contains(ARGS[z], "wspd") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
    contains(ARGS[z], "airt") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
    contains(ARGS[z], "sstt") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
    contains(ARGS[z], "shum") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
    print("$(pars[z]) variance $(specstr[a,z]) in $(dirs[a])\n")
  end

  xlaba = ["1000", "365", "182", "100", "60", "30", "14", "7", "3.5", "2"]
  xlabb = ["0.002", "0.3"]
  xposa = float(xlaba).^-1
  xposb = float(xlabb)

  contains(ARGS[z], "shfx") && (varname = "a) Sensible Heat Flux" ; ymin = -15 ; ymax = 40 ; tpos = (1,1))
  contains(ARGS[z], "lhfx") && (varname = "b) Latent Heat Flux"   ; ymin = -15 ; ymax = 40 ; tpos = (1,2))
  contains(ARGS[z], "wspd") && (varname = "c) Wind Speed"         ; ymin = -30 ; ymax = 15 ; tpos = (2,1))
  contains(ARGS[z], "shum") && (varname = "d) Specific Humidity"  ; ymin = -30 ; ymax = 15 ; tpos = (2,2))
  contains(ARGS[z], "sstt") && (varname = "e) Sea Surface Temp"   ; ymin = -30 ; ymax = 15 ; tpos = (3,1))
  contains(ARGS[z], "airt") && (varname = "f) Air Temperature"    ; ymin = -30 ; ymax = 15 ; tpos = (3,2))
  tpos[1] == 1 && (ylaba = ["-10", "0", "10", "20", "30", "40"])
  tpos[1] != 1 && (ylaba = ["-30", "-20", "-10", "0", "10"])
  yposa = float(ylaba)

  tmp = Winston.FramedPlot(
#         title="$varname Spectra (dB)", xlog = true,
#         xlabel="Timescale (days)", xrange = (1/1000,1/2), yrange = (ymin,ymax))
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
  ppp[tpos...] = Winston.add(tmp)

  ump = Array(Any, 8)
  cols = [  "red",  "blue", "green", "orange",    "red",   "blue",  "green", "orange"]
  kynd = ["solid", "solid", "solid",  "solid", "dashed", "dashed", "dashed", "dashed"]

  b = 1
  for a in order
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
  tmp = Winston.DataLabel(0.0012, ymin + 0.12 * (ymax - ymin), varname, "texthalign", "left", "size", 1.4)
        Winston.add(ppp[tpos...], tmp)
end

xyzzy = ARGS[1]*".png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
exit(0)


#=
  setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:fontsize => 30))
  tpos[2] == 2 && setattr(tmp.y1, "draw_ticklabels", false)
  setattr(tmp, "aspect_ratio", 0.2)
  tpos[1] <= 2 && setattr(tmp.x1, "draw_ticklabels", false)
  @show getattr(tmp.x1, "label_style")
  setattr(tmp.x1, "label_style", Dict{Symbol,Any}(:fontsize => 8.0, :color => "white"))
=#
