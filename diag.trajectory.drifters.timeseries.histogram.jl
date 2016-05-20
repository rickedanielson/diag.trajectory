#=
 = Loop through the timeseries of all analyses and bin the values of all available
 = variables.  Where one analysis is unavailable, all analyses are skipped - RD May 2016.
 =#

using My
const UCUR             = 1                              # indecies of all data variables
const VCUR             = 2
const PARAMS           = 2

const TIMS             = 3408                           # number in timeseries
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1 && argc != 2
  print("\nUsage: jjj $(basename(@__FILE__)) z.list [30]\n\n")
  exit(1)
end
maxfiles = 9e9 ; argc == 2 && (maxfiles = parse(Int64, ARGS[2]))

dirs = ["v2.0_global_025_deg_ekman_15m", "v2.0_global_025_deg_ekman_hs", "v2.0_global_025_deg_geostrophic", "v2.0_global_025_deg_total_15m", "v2.0_global_025_deg_total_hs"]
dirn = length(dirs)

ucui = 0.002 ; ucus = collect(-5.0 : ucui : 5.0) ; ucun = zeros(length(ucus), length(dirs))
vcui = 0.002 ; vcus = collect(-5.0 : vcui : 5.0) ; vcun = zeros(length(vcus), length(dirs))

function count(bound::Array{Float64,1}, grid::Array{Float64,2}, now::Array{Float64,1})
  flag = false
  for a = 1:dirn
    (now[a] < -333 || now[a] > 333) && (flag = true)
  end
  flag && return

  for a = 1:dirn
    delnow, indnow = findmin(abs(bound - now[a])) ; bound[indnow] > now[a] && indnow > 1 && (indnow -= 1)
    grid[indnow,a] += 1
  end
end                                                                           # (grid boundaries refer to lower limits)

nfile = 0                                                                     # loop through the list of locations and
fpn  = Array(IOStream, dirn)                                                  # grid the three timeseries, where valid
data = Array(Float64,  dirn, PARAMS)
fpa = My.ouvre(dirs[1] * "/" * ARGS[1], "r")
files = readlines(fpa) ; close(fpa)
for (a, fila) in enumerate(files)
  (z, tail) = split(strip(fila), dirs[1])
  for b = 1:dirn
    fpn[b] = My.ouvre("$(dirs[b])/$(dirs[b])$tail", "r")
  end

  for b = 1:TIMS
    for c = 1:dirn
      line = readline(fpn[c]) ; vals = split(line)
      data[c,UCUR] = float(vals[10])
      data[c,VCUR] = float(vals[11])
    end
    count(ucus, ucun, data[:,UCUR])
    count(vcus, vcun, data[:,VCUR])
  end

  for b = 1:dirn
    close(fpn[b])
  end
  nfile += 1 ; if nfile >= maxfiles  break  end
end
print("read $nfile files\n")

function store(bound::Array{Float64,1}, grid::Array{Float64,2}, pname::UTF8String)
  fname = "histogr." * pname * ".dat"
  fpa = My.ouvre(fname, "w")
  for (a, vala) in enumerate(bound)
    @printf(fpa, "%15.8f %15.8f %15.8f %15.8f %15.8f\n",
      grid[a,1], grid[a,2], grid[a,3], grid[a,4], grid[a,5])
  end
  close(fpa)
end

store(ucus, ucun, utf8("ucur"))
store(vcus, vcun, utf8("vcur"))
exit(0)


#=
  CFSR = try  findin(dirs, [  "cfsr"])[1]  catch  0  end
  JOFU = try  findin(dirs, ["jofuro"])[1]  catch  0  end
    if CFSR > 0  data[CFSR,LHFX]                   = 333  end                 # set expected missing values to be
    if JOFU > 0  data[JOFU,AIRT] = data[JOFU,SSTT] = 333  end                 # outside the plotting range
    if float(valb[ 1]) < 50 && (float(vala[ 1]) > 300 || float(valc[ 1]) > 300)
      @show vala ; @show valb ; @show valc  end
    @printf(fpa, "%15.8f\n", grid[a,1])
=#
