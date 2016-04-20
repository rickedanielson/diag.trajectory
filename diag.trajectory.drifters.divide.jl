#=
 = Split the GlobCurrent drifter track file from M.-H. Rio
 = into a set of individual track files and convert the CNES
 = Julian day to a regular date - RD July 2015
 =#

using My
const LEN              = 100
const LOTS             = 1000
const START            = 2                              # make START-1 a valid array index

if size(ARGS) != (1,)
  print("\nUsage: julia $(basename(@__FILE__)) buoydata_1993_2012_drogON.asc\n\n")
  exit(1)
end

lins = Array(UTF8String, 1)                                                   # initialize some arrays with the first entry undefined
buoy = Array(UTF8String, 1)                                                   # and starting with the second entry, store new tracks
jday = Array(UTF8String, 1)

n = 0 ; i = START
fpa = My.ouvre(ARGS[1], "r")
for line in eachline(fpa)
  vars = split(line)
  push!(buoy, vars[1])                                                        # IDbuoy is the buoy identifier
# push!(drog, vars[3])                                                        # DROG is the buoy death code*
  push!(jday, My.dateadd("1950-01-01-0000", float(vars[4]), "dy"))            # Jcnes is a CNES Julian day (number of days since 1 Jan 1950)
# push!(lonn, vars[5])                                                        # LON is the drifter position longitude (0-360)
# push!(latt, vars[6])                                                        # LAT is the drifter position latitude (-90 90)
  push!(lins, line)

  if i != START && buoy[i] != buoy[i-1]                                       # at end of one track and beginning of next:
if (float(jday[START][1:4]) > 2009)
    fpb = My.ouvre(jday[START]*"."*("00000000"*buoy[i-1])[end-7:end], "w")    # write the previous track file (8-digit number)
    for j = START:i-1  write(fpb, "$(jday[j]) $(lins[j])")  end
    close(fpb)
end

    lins = Array(UTF8String, 1) ; push!(lins, line)                           # then reset arrays with the new starting track line
    buoy = Array(UTF8String, 1) ; push!(buoy, vars[1])
    jday = Array(UTF8String, 1) ; push!(jday, My.dateadd("1950-01-01-0000", float(vars[4]), "dy"))
    n += 1 ; i = START
  end
  i += 1
end
close(fpa)

fpb = My.ouvre(jday[START]*"."*("00000000"*buoy[i-1])[end-7:end], "w")        # write the last track file (8-digit number)
for j = START:i-1  write(fpb,"$(jday[j]) $(lins[j])")  end
close(fpb)

print("wrote $(n+1) track files\n\n")
exit(0)


#=       *DEATH CODES
0       = buoy still reporting as of 6-30-13
1       = buoy ran aground
2       = buoy was picked up
3       = buoy quit transmitting
4       = Unreliable transmissions at end of trajectory
5       = Bad battery voltage
6       = Place in inactive status while transmitting good position  =#
