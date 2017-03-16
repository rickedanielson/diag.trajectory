#=
 = Copy the landmask research.jisao.washington.edu/data_sets/elevation/fractional_land.0.25-deg.nc
 = (which is conveniently defined on the same grid as AVISO) and replace landmask with distance to
 = the coast - RD January 2017.
 =#

using My, NetCDF
const LATS             = 720                            # number of grid latitudes
const LONS             = 1440                           # number of grid longitudes

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) fractional_land.0.25-deg.nc\n\n")
  exit(1)
end

infile = ARGS[1]                                                              # create the output file
utfile = ARGS[1] * ".dist.to.coast.nc"
command = "cp $infile $utfile\n"
print(command) ; run(`$(split(command))`)

dist = Array(  Int16,     LONS, LATS)                                         # copy the rescaled landmask
land = Array(Float64, 3 * LONS, LATS)                                         # three times in longitude
lons = Array(Float64, 3 * LONS)
data = map(Float64, ncread(utfile, "data", start=[1,1,1], count=[LONS,LATS,1])[:,:,1])
scal =            ncgetatt(utfile, "data", "scale_factor")
lats = map(Float64, ncread(utfile,  "lat", start=[1],     count=[LATS]))
temp = map(Float64, ncread(utfile,  "lon", start=[1],     count=[LONS]))
for a = 1:LONS  lons[a] = lons[a+LONS] = lons[a+2*LONS] = temp[a]  end
for a = 1:LATS, b = 1:LONS
  land[b,a] = land[b+LONS,a] = land[b+2*LONS,a] = scal * data[b,a]
end

function mindis(lats::Array{Float64,1}, lons::Array{Float64,1}, land::Array{Float64,2}, dist::Array{Int16,2})
  for a = 1:LATS, b = 1:LONS
    if land[b,a] > 0.9
      dist[b,a] = 0
    else                                                                      # for locations over the ocean
      lata = lats[a] ; lona = lons[b] ; mindis = 32767000                     # get the minimum distance to land
      flag = true ; z = 1                                                     # (defined as landmask > 0.9) by a
      while flag                                                              # search over an increasingly large
        for c = -z:z-1, d = -z                                                # square perimeter (not ideal)
          e = a + c ; f = b + d + LONS                                        # starting on the west perimeter
          if 1 <= e <= LATS && land[f,e] > 0.9
            latb = lats[e] ; lonb = lons[f]
            tmp = My.latlondist(lata, lona, latb, lonb)
            if tmp < mindis  mindis = tmp ; flag = false  end
          end
        end
        for c = z, d = -z:z-1                                                 # south perimeter
          e = a + c ; f = b + d + LONS
          if 1 <= e <= LATS && land[f,e] > 0.9
            latb = lats[e] ; lonb = lons[f]
            tmp = My.latlondist(lata, lona, latb, lonb)
            if tmp < mindis  mindis = tmp ; flag = false  end
          end
        end
        for c = z:-1:1-z, d = z                                               # east perimeter
          e = a + c ; f = b + d + LONS
          if 1 <= e <= LATS && land[f,e] > 0.9
            latb = lats[e] ; lonb = lons[f]
            tmp = My.latlondist(lata, lona, latb, lonb)
            if tmp < mindis  mindis = tmp ; flag = false  end
          end
        end
        for c = -z, d = z:-1:1-z                                              # north perimeter
          e = a + c ; f = b + d + LONS
          if 1 <= e <= LATS && land[f,e] > 0.9
            latb = lats[e] ; lonb = lons[f]
            tmp = My.latlondist(lata, lona, latb, lonb)
            if tmp < mindis  mindis = tmp ; flag = false  end
          end
        end
        z += 1
      end
      dist[b,a] = round(Int16, mindis / 1000.0)
    end
  end
end

mindis(lats, lons, land, dist)
ncwrite(dist, utfile, "data", start=[1,1,1], count=[LONS,LATS,1])             # write the distance and close
ncclose(utfile)
exit(0)

#=
dist[b,a] = round(Int16, 32760)
@show lata, lona, dist[b,a]
=#
