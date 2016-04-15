#=
 = Plot residence time within an irregular domain by reading all
 = trajectory files in this dir and constructing an histogram - RD
 = September 2015
 =#

using My, LsqFit
const LEN              = 1000
const LOTS             = 10000

if size(ARGS) != (0,)
  print("\nUsage: jj $(basename(@__FILE__))\n\n")
  exit(1)
end

trajfiles = filter(x -> (Base.contains(x, ".traj.") &&                        # then identify all simulated
                        !Base.contains(x,    ".nc") &&                        # trajectory files in this dir
                         filesize(x) > 150), readdir("."))

presence = 0.0
transitions = 0.0
count = zeros(Int64, LOTS)                                                    # and for each, translate the column
for (a, fila) in enumerate(trajfiles)                                         # indicating presence in the domain to
  str = ""                                                                    # a string and then count word length
  fpa = My.ouvre(fila, "r", false)
  for (b, linb) in enumerate(eachline(fpa))
    tmp = split(linb)
    simind = int(float(tmp[7]))
    str *= simind > 1 ? "a" : " "
  end
  close(fpa)
  tmp = split(str)
  for word in tmp
    count[length(word)] += 1
  end

  strlen = length(str) < 700 ? length(str) : 700                              # also within the first six months or so
  for b = 2:strlen                                                            # count the presence in, and transition
    if                    str[b] == 'a'  presence    += 1.0  end              # (per day) into and out of, the domain
    if str[b-1] == ' ' && str[b] == 'a'  transitions += 1.0  end              # and then get an average residence time
    if str[b-1] == 'a' && str[b] == ' '  transitions += 1.0  end
  end
end
presence    /= 700.0
transitions /= 350.0
residence = presence / transitions
summary = @sprintf("Daily averages over\n the first six months\n\n%.2f residence time (days)\n%.2f drift presence\n   %.2f boundary flux (/day)\n",
  residence, presence, transitions)
print(summary)

xvals = zeros(Float64, 60)                                                    # then sum the 6-h residence times
yvals = zeros(Float64, 60)                                                    # over 2.5-day increments
a = 1
for b = 1:60 # LOTS / 10
  xvals[b] = b * 2.5
  yvals[b] = 0.0
  for c = 1:10
    yvals[b] += count[a]
    a += 1
  end
end

xyzzy = "lofoten_basin_residence_time.png"
print("writing $xyzzy\n")
using Winston
tmp = Winston.FramedPlot(
        title="Lofoten Basin simulated drift (2010-2012)",
        xlabel="Residence time (days)", xrange = (0,150),
        ylabel="Frequency",             yrange = (0,40))
ppp = Winston.add(tmp)
tmp = Winston.Curve(xvals, yvals, "color", parse(Winston.Colorant, "red"))
      Winston.add(ppp, tmp)
tmp = Winston.PlotLabel(0.3, 0.7, summary, "texthalign", "left", "size", 3.0)
      Winston.add(ppp, tmp)
      Winston.savefig(ppp, xyzzy, "width", 700, "height", 700)
exit(0)


#=
model(x, p) = p[1]./(x.*p[2]).^2.0
fit = curve_fit(model, xvals[1:60], yvals[1:60], [200.0, 5.0])
zvals = model(xvals, fit.param)
print("$(fit.param[1]) $(fit.param[2])\n")

model(x, p) = p[3] .* exp(-p[2] .* x) .* exp(-p[1] .* exp(-p[2] .* x))
model(x, p) = p[1]*exp(-x.*p[2])
avals = model(xvals, [1387.3688908760541 -1.0])

tmp = Winston.Curve(xvals, zvals, "color", parse(Winston.Colorant, "blue"))
      Winston.add(ppp, tmp)
tmp = Winston.Curve(xvals, avals, "color", parse(Winston.Colorant, "green"))
      Winston.add(ppp, tmp)
=#
