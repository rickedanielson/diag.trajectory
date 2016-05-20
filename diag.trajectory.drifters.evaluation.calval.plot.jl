#=
 = Loop through all analyses and plot the binned
 = sums of all available variables - RD April 2016.
 =#

using My, Winston
const TARG             = 3
const ALPH             = 4                              # error model x = ALPH + BETA * truth + error
const BETA             = 5                              # error model x = ALPH + BETA * truth + error
const SIGM             = 6                              # triple coll RMSE
const CORR             = 7                              # triple coll correlation coefficient
const PARAMS           = 4

const FNUM             = 4                              # number of files expected
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != FNUM
  print("\nUsage: jjj $(basename(@__FILE__)) *cali.locl\n\n")
  exit(1)
end

tars = collect(0.1:0.1:1.3)
tarn = zeros(length(tars), FNUM)
alpn = zeros(length(tars), FNUM) ; alpo = zeros(length(tars), FNUM)
betn = zeros(length(tars), FNUM) ; beto = zeros(length(tars), FNUM)
sign = zeros(length(tars), FNUM) ; sigo = zeros(length(tars), FNUM)
corn = zeros(length(tars), FNUM) ; coro = zeros(length(tars), FNUM)

for a = 1:FNUM
  fpa = My.ouvre(ARGS[a], "r")
  readline(fpa)
  for b = 1:length(tars)
    line = readline(fpa) ; vals = split(line)
    tarn[b,a] = float(vals[TARG])
    alpn[b,a] = float(vals[ALPH]) ; betn[b,a] = float(vals[BETA])
    sign[b,a] = float(vals[SIGM]) ; corn[b,a] = float(vals[CORR])
  end

# const VCLAV = [    -0.02881019      0.00653482      0.01137809]
  ln = readline(fpa) ; tmp = float(split(ln[16:end-2])) ; function fun(x)  tmp[1]x.^2 + tmp[2]x + tmp[3]  end ; alpo[:,a] = fun(tars)
  ln = readline(fpa) ; tmp = float(split(ln[16:end-2])) ; function fun(x)  tmp[1]x.^2 + tmp[2]x + tmp[3]  end ; beto[:,a] = fun(tars)
  ln = readline(fpa) ; tmp = float(split(ln[16:end-2])) ; function fun(x)  tmp[1]x.^2 + tmp[2]x + tmp[3]  end ; sigo[:,a] = fun(tars)
  ln = readline(fpa) ; tmp = float(split(ln[16:end-2])) ; function fun(x)  tmp[1]x.^2 + tmp[2]x + tmp[3]  end ; coro[:,a] = fun(tars)
  close(fpa)
end

ppp = Winston.Table(2,2) ; setattr(ppp, "cellpadding", -0.5)                  # and then create the plots
for z = ALPH:CORR
  z == ALPH && (varname = "a) Bias (ms^{-1})" ; bound = tarn ; grid = alpn ; tpos = (1,1) ; grie = alpo)
  z == BETA && (varname = "b) Slope"          ; bound = tarn ; grid = betn ; tpos = (1,2) ; grie = beto)
  z == SIGM && (varname = "c) RMSE (ms^{-1})" ; bound = tarn ; grid = sign ; tpos = (2,1) ; grie = sigo)
  z == CORR && (varname = "d) Correlation"    ; bound = tarn ; grid = corn ; tpos = (2,2) ; grie = coro)

  z == ALPH && (xmin = 0.05 ; xmax = 1.35 ; ymin = -0.1 ; ymax = 0.05)           # and locate the plot limits
  z == BETA && (xmin = 0.05 ; xmax = 1.35 ; ymin =  0.4 ; ymax = 3.5)
  z == SIGM && (xmin = 0.05 ; xmax = 1.35 ; ymin =  0.0 ; ymax = 0.2)
  z == CORR && (xmin = 0.05 ; xmax = 1.35 ; ymin =  0.8 ; ymax = 1.0)

  ump = Array(Any, 2*FNUM)
  cols = [    "red",    "blue",   "green",  "orange",     "red",    "blue",   "green",  "orange"]
  kynd = [  "solid",   "solid",   "solid",   "solid",  "dashed",  "dashed",  "dashed",  "dashed"]
  dirs = ["Grp-A U", "Grp-A V", "Grp-B U", "Grp-B V", "Est-A U", "Est-A V", "Est-B U", "Est-B V"]
# xmin = 0.0 ; xmax = 1.4 ; ymin = -0.5 ; ymax = 1.0

  tmp = Winston.FramedPlot(title="$varname", xrange = (xmin,xmax), yrange = (ymin,ymax))
  ppp[tpos...] = Winston.add(tmp)

  for a = 1:FNUM
    ump[a]      = Winston.Curve(bound[:,a], grid[:,a], "color", parse(Winston.Colorant, cols[a]))
                  style(ump[a], kind = kynd[a])
                  setattr(ump[a], label = dirs[a])
                  Winston.add(ppp[tpos...], ump[a])
    ump[a+FNUM] = Winston.Curve(      tars, grie[:,a], "color", parse(Winston.Colorant, cols[a+FNUM]))
                  style(ump[a+FNUM], kind = kynd[a+FNUM])
                  setattr(ump[a+FNUM], label = dirs[a+FNUM])
                  Winston.add(ppp[tpos...], ump[a+FNUM])
  end
  if z == BETA
    tmp = Winston.Legend(.45, .82, Any[ump[1], ump[2], ump[3], ump[4]])
          Winston.add(ppp[tpos...], tmp)
    tmp = Winston.Legend(.70, .82, Any[ump[5], ump[6], ump[7], ump[8]])
          Winston.add(ppp[tpos...], tmp)
  end
end

xyzzy = "calval_globcurrent.png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
exit(0)
