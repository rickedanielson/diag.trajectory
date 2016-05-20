#=
 = Loop through the six-hourly Ekman and geostrophic extrapolated
 = timeseries and add them where both are valid - RD May 2016.
 =#

using My
const TIMS             = 3408                           # number in timeseries
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 4
  print("\nUsage: jjj $(basename(@__FILE__)) v2.0_global_025_deg_geostrophic v2.0_global_025_deg_ekman_15m v2.0_global_025_deg_total_15m z.listaa\n\n")
  exit(1)
end

fpa = My.ouvre("$(ARGS[1])/$(ARGS[4])", "r")                                  # loop through the list of locations
files = readlines(fpa) ; close(fpa)                                           # and process each timeseries
for fila in files
  (z, tail) = split(strip(fila), ARGS[1])
  filb = ARGS[1] * "/" * ARGS[1] * tail * ".bef"
  filc = ARGS[2] * "/" * ARGS[2] * tail * ".bef"
  fild = ARGS[3] * "/" * ARGS[3] * tail * ".bef"
  fpb  = My.ouvre(filb, "r", false)
  fpc  = My.ouvre(filc, "r", false)
  fpd  = My.ouvre(fild, "w")
  for a = 1:TIMS
    line = readline(fpb) ; vals = split(line) ; ub = float(vals[10]) ; vb = float(vals[11])
    line = readline(fpc) ; vals = split(line) ; uc = float(vals[10]) ; vc = float(vals[11])
    ud = -333 < ub < 333 && -333 < uc < 333 ? ub + uc : MISS
    vd = -333 < vb < 333 && -333 < vc < 333 ? vb + vc : MISS
    form = @sprintf("%s %11.5f %11.5f\n", line[1:92], ud, vd)
    write(fpd, form)
  end
  close(fpb)
  close(fpc)
  close(fpd)

  filb = ARGS[1] * "/" * ARGS[1] * tail * ".aft"
  filc = ARGS[2] * "/" * ARGS[2] * tail * ".aft"
  fild = ARGS[3] * "/" * ARGS[3] * tail * ".aft"
  fpb  = My.ouvre(filb, "r", false)
  fpc  = My.ouvre(filc, "r", false)
  fpd  = My.ouvre(fild, "w")
  for a = 1:TIMS
    line = readline(fpb) ; vals = split(line) ; ub = float(vals[10]) ; vb = float(vals[11])
    line = readline(fpc) ; vals = split(line) ; uc = float(vals[10]) ; vc = float(vals[11])
    ud = -333 < ub < 333 && -333 < uc < 333 ? ub + uc : MISS
    vd = -333 < vb < 333 && -333 < vc < 333 ? vb + vc : MISS
    form = @sprintf("%s %11.5f %11.5f\n", line[1:92], ud, vd)
    write(fpd, form)
  end
  close(fpb)
  close(fpc)
  close(fpd)
end
exit(0)
