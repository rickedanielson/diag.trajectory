#=
 = Loop through the daily geostrophic timeseries and create the corresponding forward and
 = backward extrapolated timeseries for the two (u,v) current components.  First convert the
 = six-hourly data to daily and include a gap of a day or so between the end of the timeseries
 = and the extrapolated day of interest (at either end) - RD April 2016.
 =#

using My, Interpolations
const UCUR             = 1                              # indecies of all data variables
const VCUR             = 2
const PARS             = 2

const BEF              = 1                              # indecies of the source of extrapolations
const NOW              = 2
const AFT              = 3
const SRCS             = 3

const EXTRA            = 9                              # number of points used for extrapolation
const EXTRB            = 2                              # number of points in gap between EXTRA and target
const TIMS             = 3408                           # number in timeseries
const DAYS             = 852                            # days in timeseries (TIMS/4)
const MISS             = -9999.0                        # generic missing value

if size(ARGS) != (2,)
  print("\nUsage: jjj $(basename(@__FILE__)) v2.0_global_025_deg_geostrophic z.listah\n\n")
  exit(1)
end

inner = div(EXTRA - 1, 2)
outer = div(EXTRA + 1, 2)
dats = Array(UTF8String, TIMS)
datt = Array(Float64,    TIMS, SRCS, PARS)
datb = Array(Float64,    DAYS, SRCS, PARS)

fpa = My.ouvre("$(ARGS[1])/$(ARGS[2])", "r")                                  # loop through the list of locations
files = readlines(fpa) ; close(fpa)                                           # and process each timeseries
for fila in files
  fila = strip(fila)
  fpa = My.ouvre("$(ARGS[1])/$fila", "r", false)
  lines = readlines(fpa) ; close(fpa)
  for (a, line) in enumerate(lines)
    vals = split(line)
    datt[a,NOW,UCUR] = float(vals[10])
    datt[a,NOW,VCUR] = float(vals[11])
    dats[a]          =       vals[1]
  end

  for c = 1:PARS                                                              # retain the daily geostrophy timestep
    for b = 1:SRCS                                                            # (as in the source data)
      for a = 1:DAYS
        data[a,b,c] = datt[4*a,b,c]
      end
    end
  end

  for b = 1:PARS                                                              # set to missing the first few BEF
    for a = 1:EXTRA+1                                                         # extrapolations (not defined below)
      data[a,BEF,b] = MISS
    end
  end

#=for a = 1:PARS                                                              # simultaneously extrap from BEF and AFT
    for b = 1+outer:TIMS-outer
      tmp = vec(data[a,NOW,b-inner:b+inner])
      if all(-333 .< tmp .< 3333)
        itp = interpolate(tmp, BSpline(Quadratic(Line())), OnCell())
        data[a,BEF,b+outer] = itp[10]
        data[a,AFT,b-outer] = itp[0]
      else
        data[a,BEF,b+outer] = data[a,AFT,b-outer] = MISS
      end
    end
  end =#

  for a = 1:PARS                                                              # set to missing the last few AFT
    for b = 0:EXTRA                                                           # extrapolations (not defined above)
      data[a,AFT,TIMS-b] = MISS
    end
  end

  filb = "$fila.bef"                                                          # then save all extrapolations
  filc = "$fila.aft"
  fpb = My.ouvre("$(ARGS[1])/$filb", "w", false)
  fpc = My.ouvre("$(ARGS[1])/$filc", "w", false)
  (lll, lat, lon) = split(replace(fila, r"[\.]{2,}", " "))
  for a = 1:TIMS
    formb = @sprintf("%10s %10.5f %10.5f   99999999     99999  9 -9999.00000 -9999.00000 -9999.00000 %11.5f %11.5f\n",
      dats[a], float(lat), float(lon), datt[UCUR,BEF,a], datt[VCUR,BEF,a])
    formc = @sprintf("%10s %10.5f %10.5f   99999999     99999  9 -9999.00000 -9999.00000 -9999.00000 %11.5f %11.5f\n",
      dats[a], float(lat), float(lon), datt[UCUR,AFT,a], datt[VCUR,AFT,a])
    write(fpb, formb)
    write(fpc, formc)
  end
  close(fpb)
  close(fpc)
end
exit(0)


#=
2012090600   13.87500  141.37500   99999999     99999  9 -9999.00000 -9999.00000 -9999.00000    -0.11356    -0.02688
=#
