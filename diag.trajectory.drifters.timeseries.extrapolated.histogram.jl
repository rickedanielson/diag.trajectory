#=
 = Loop through the timeseries and grid the corresponding forward and
 = backward extrapolated timeseries of all available variables, relative to the actual
 = (uninterpolated) values.  Note that BEF refers to an interpolation using analysis
 = data from before the extrapolation; AFT extrapolations use analysis data afterward.
 = Where one analysis is unavailable, all analyses are skipped - RD May 2016.
 =#

using My
const UCUR             = 1                              # indecies of all data variables
const VCUR             = 2
const PARAMS           = 2

const BEF              = 1                              # indecies of the source of extrapolations
const NOW              = 2
const AFT              = 3
const SRCS             = 3

const TIMS             = 3408                           # number in timeseries
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1 && argc != 2
  print("\nUsage: jjj $(basename(@__FILE__)) z.list [30]\n\n")
  exit(1)
end
maxfiles = 9e9 ; argc == 2 && (maxfiles = parse(Int64, ARGS[2]))

dirs = ["v2.0_global_025_deg_ekman_15m", "v2.0_global_025_deg_ekman_hs", "v2.0_global_025_deg_geostrophic", "v2.0_global_025_deg_total_15m", "v2.0_global_025_deg_total_hs"]
dirn = length(dirs)

ucui = 0.002 ; ucus = collect(-5.0 : ucui : 5.0) ; ucun = zeros(length(ucus), length(ucus), length(dirs))
vcui = 0.002 ; vcus = collect(-5.0 : vcui : 5.0) ; vcun = zeros(length(vcus), length(vcus), length(dirs))

function count(bound::Array{Float64,1}, grid::Array{Float64,3}, bef::Array{Float64,1}, now::Array{Float64,1}, aft::Array{Float64,1})
  flag = false
  for a = 1:dirn
    (bef[a] < -333 || now[a] < -333 || aft[a] < -333 ||
     bef[a] >  333 || now[a] >  333 || aft[a] >  333) && (flag = true)
  end
  flag && return

  for a = 1:dirn
    delbef, indbef = findmin(abs(bound - bef[a])) ; bound[indbef] > bef[a] && indbef > 1 && (indbef -= 1)
    delnow, indnow = findmin(abs(bound - now[a])) ; bound[indnow] > now[a] && indnow > 1 && (indnow -= 1)
    delaft, indaft = findmin(abs(bound - aft[a])) ; bound[indaft] > aft[a] && indaft > 1 && (indaft -= 1)
    grid[indbef,indnow,a] += 1 ; grid[indaft,indnow,a] += 1
  end
end                                                                           # (grid boundaries refer to lower limits)

nfile = 0                                                                     # loop through the list of locations and
fpna = Array(IOStream, dirn)                                                  # grid the three timeseries, where valid
fpnb = Array(IOStream, dirn)
fpnc = Array(IOStream, dirn)
data = Array(Float64,  dirn, PARAMS, SRCS)
fpa = My.ouvre(dirs[1] * "/" * ARGS[1], "r")
files = readlines(fpa) ; close(fpa)
for (a, fila) in enumerate(files)
  (z, tail) = split(strip(fila), dirs[1])
  for b = 1:dirn
    fpna[b] = My.ouvre("$(dirs[b])/$(dirs[b])$tail.bef", "r", false)
    fpnb[b] = My.ouvre("$(dirs[b])/$(dirs[b])$tail",     "r")
    fpnc[b] = My.ouvre("$(dirs[b])/$(dirs[b])$tail.aft", "r", false)
  end

  for b = 1:TIMS
    for c = 1:dirn
      line = readline(fpna[c]) ; vala = split(line)
      line = readline(fpnb[c]) ; valb = split(line)
      line = readline(fpnc[c]) ; valc = split(line)
      data[c,UCUR,BEF] = float(vala[10]) ; data[c,UCUR,NOW] = float(valb[10]) ; data[c,UCUR,AFT] = float(valc[11])
      data[c,VCUR,BEF] = float(vala[11]) ; data[c,VCUR,NOW] = float(valb[11]) ; data[c,VCUR,AFT] = float(valc[11])
    end
    count(ucus, ucun, data[:,UCUR,BEF], data[:,UCUR,NOW], data[:,UCUR,AFT])
    count(vcus, vcun, data[:,VCUR,BEF], data[:,VCUR,NOW], data[:,VCUR,AFT])
  end

  for b = 1:dirn
    close(fpna[b])
    close(fpnb[b])
    close(fpnc[b])
  end
  nfile += 1 ; if nfile >= maxfiles  break  end
end
print("read $nfile files\n")

function store(bound::Array{Float64,1}, grid::Array{Float64,3}, pname::UTF8String)
  fname = "extrapolated.histogr." * pname * ".dat"
  fpa = My.ouvre(fname, "w")
  for (a, vala) in enumerate(bound)
    for (b, valb) in enumerate(bound)
      @printf(fpa, "%15.8f %15.8f %15.8f %15.8f %15.8f\n",
        grid[b,a,1], grid[b,a,2], grid[b,a,3], grid[b,a,4], grid[b,a,5])
    end
  end
  close(fpa)
end

store(ucus, ucun, utf8("ucur"))
store(vcus, vcun, utf8("vcur"))
exit(0)


#=
  CFSR = try  findin(dirs, [  "cfsr"])[1]  catch  0  end
  JOFU = try  findin(dirs, ["jofuro"])[1]  catch  0  end
    if CFSR > 0  data[CFSR,LHFX,:]                     = [3333 3333 3333] end  # set expected missing values to be
    if JOFU > 0  data[JOFU,AIRT,:] = data[JOFU,SSTT,:] = [3333 3333 3333] end  # outside the plotting range
=#
