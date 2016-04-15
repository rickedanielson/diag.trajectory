#=
 = Isolate the drifters not in the CNES-CLS-2013 MDT (so pass only data
 = from September 2012 onward) and convert the CNES Julian day to a regular
 = date - RD April 2016.
 =#

using My

if size(ARGS) != (1,)
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2012_drogON.asc\n\n")
  exit(1)
end

fpa = My.ouvre(ARGS[1],             "r")
fpb = My.ouvre(ARGS[1] * ".nonmdt", "w")

for line in eachline(fpa)
  vars = split(line)
  if float(vars[4]) >= 22889
    jday = My.dateadd("1950-01-01-0000", float(vars[4]), "dy")
    write(fpb, jday * " " * line)
  end
end

close(fpa)
close(fpb)
exit(0)
