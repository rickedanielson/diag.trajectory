#=
 = Identify all files like 2010-06-17-0000.00092861 in the current dir and for
 = each of the corresponding simulated trajectory files that exists in this dir
 = (and was launched from the position of the real drifter), identify the shapefile
 = domain that contains this initial position and verify that links to the trajectory
 = file (and perhaps links to all related files also) exist from a parallel dir that
 = is specific to that shapefile domain - RD September 2015
 =#

using My, LibGEOS
const LEN              = 101
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 0
  print("\nUsage: jj $(basename(@__FILE__))\n\n")
  exit(1)
end

polygons = Array(Polygon, 0)
##(polynames, polypoints) = My.obtenpoly("/home/ricani/data/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS")
##(polynames, polypoints) = My.obtenpoly("/home/ricani/data/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS_valid")
##(polynames, polypoints) = My.obtenpoly("/home/ricani/data/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS_valid_simple")
# (polynames, polypoints) = My.obtenpoly("/home/ricani/data/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS.gra/GlobCurrent_shapefiles_CNES-CLS13_RSMAS")
##(polynames, polypoints) = My.obtenpoly("/home/ricani/data/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS_valid.gra/GlobCurrent_shapefiles_CNES-CLS13_RSMAS_valid")
  (polynames, polypoints) = My.obtenpoly("/home/cercache/users/rdaniels/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS_valid")
for (a, line) in enumerate(polypoints)
# print("$a :$line: $a\n")
# if a == 13 || a == 16
#   tmp = LibGEOS.parseWKT("POLYGON((888 888,889 888,889 889,888 889,888 888))")
#   push!(polygons, tmp)
# else
    tmp = LibGEOS.parseWKT("POLYGON(($line))")
    push!(polygons, tmp)
# end
end

files = filter(x -> ismatch(r"^\d{4}-\d{2}-\d{2}-\d{4}\.\d{8}$", x), readdir("."))
(count,) = size(files) ; print("found $count files\n")
curdir = split(pwd(), "/")[end]

for (a, fila) in enumerate(files)
  if fila[1:3] == "201"
    fpa = My.ouvre(fila, "r", false)
    for (b, linb) in enumerate(eachline(fpa))
      tmp = split(linb)
      sim = fila*".traj."*tmp[1]
      if isfile(sim)
        lon = float(tmp[6]) ; lon > 180 && (lon -= 360)
        lat = float(tmp[7])
        loc = LibGEOS.parseWKT("POINT($lon $lat)")
        for (c, polyc) in enumerate(polygons)
          if LibGEOS.contains(polyc, loc)
            newdir = "../$curdir.$(polynames[c])"
            newlnk = "../$curdir.$(polynames[c])/$sim"
            tarlnk = "../$curdir/$sim"
             !isdir(newdir) &&           mkdir(newdir)
            !islink(newlnk) && symlink(tarlnk, newlnk)
            newlnk = "../$curdir.$(polynames[c])/$fila"
            tarlnk = "../$curdir/$fila"
            !islink(newlnk) && symlink(tarlnk, newlnk)
          end
        end
      end
    end
    close(fpa)
  end
end
