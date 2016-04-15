#=
 = Categorize positions along each of a set of simulated trajectories
 = according to whether positions fall inside or outside an irregular
 = domain, which is given as a shapefile.  A nearly identical set of
 = trajectory files is written in place - RD September 2015
 =#

using My, LibGEOS
const LEN              = 100

if size(ARGS) != (0,)
  print("\nUsage: jj $(basename(@__FILE__))\n\n")
  exit(1)
end

poly = ""                                                                     # read the domain of interest
fila = "/home/ricani/work/workt/norskhavet.lofoten"                           # (as a GEOS POLYGON)
fpa = My.ouvre(fila, "r") ; lines = readlines(fpa) ; close(fpa)
for (a, line) in enumerate(lines)
  (lon, lat) = split(line) ; poly *= "$lon $lat,"
end
lofoten = parseWKT("POLYGON(($(poly[1:end-1])))")

trajfiles = filter(x -> (Base.contains(x, ".traj.") &&                        # then identify all simulated
                        !Base.contains(x,    ".nc") &&                        # trajectory files in this dir
                         filesize(x) > 150), readdir("."))

for (a, fila) in enumerate(trajfiles)                                         # and update them with a nonzero
  filb = fila*".xyzzy"                                                        # collocation indicator (in place)
  fpa = My.ouvre(fila, "r")
  fpb = My.ouvre(filb, "w")
  for (b, linb) in enumerate(eachline(fpa))
    tmp = split(linb)
    simlat = float(tmp[2])
    simlon = float(tmp[3]) ; simlon > 180 && (simlon -= 360)
    locate = LibGEOS.parseWKT("POINT($simlon $simlat)")
    if LibGEOS.contains(lofoten, locate)
      outline = @sprintf("%s %13.5f %13.5f\n", linb[1:75], 2.0, 0.0)
    else
      outline = @sprintf("%s %13.5f %13.5f\n", linb[1:75], 1.0, 0.0)
    end
    write(fpb, outline)
  end
  close(fpb)
  close(fpa)

  print("mv $filb $fila\n")
  mv(filb, fila)
end


#    g1 = parseWKT("POLYGON((1 1,1 5,5 5,5 1,1 1))")
#    g2 = parseWKT("POINT(2 2)")
#    contains(g1, g2) => true
#    contains(g2, g1) => false
#geom1_ = LibGEOS.geomFromWKT("POLYGON((1 1,1 5,5 5,5 1,1 1))")
#geom2_ = LibGEOS.geomFromWKT("POINT(2 2)")
#@fact LibGEOS.within(geom1_, geom2_) => false
#@fact LibGEOS.within(geom2_, geom1_) => true
#LibGEOS.destroyGeom(geom1_)
#LibGEOS.destroyGeom(geom2_)
