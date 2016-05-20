#=
 = Reformat second-order polynomial coeff text files as grads
 = display commands and markdown tables.  It is assumed that
 = the files passed as arguments are the available group of
 = analyses for a given variable (e.g., wspd) - RD April 2016.
 =#

using My
const VARS             = 4                              # number of triple collocation metrics
const COEF             = 10                             # number of polynomial coefficients
const MISS             = -9999.0                        # generic missing value
const VARN             = ["Bias", "Slope", "RMSE", "Target Correlation"]

if (argc = length(ARGS)) == 0
  print("\nUsage: jjj $(basename(@__FILE__)) all.flux.daily.locate_2.0_calib.wspd.got2000_obs.comb.*\n\n")
  exit(1)
end

contains(ARGS[1], "shfx") && (varname = "Sensible Heat Flux")
contains(ARGS[1], "lhfx") && (varname = "Latent Heat Flux")
contains(ARGS[1], "wspd") && (varname = "Wind Speed")
contains(ARGS[1], "airt") && (varname = "Air Temperature")
contains(ARGS[1], "sstt") && (varname = "SST")
contains(ARGS[1], "shum") && (varname = "Specific Humidity")

data = Array(Float64, COEF, VARS, argc)

for a = 1:argc                                                                # read the four types of triple collocation
  fpa = My.ouvre(ARGS[a], "r") ; lines = readlines(fpa) ; close(fpa)          # polynomial coefficients for each analysis
  for b = 1:VARS
    vals = float(split(lines[b]))
    for c = 1:COEF
      data[c,b,a] = vals[c]
    end
  end
end

for a = 1:argc                                                                # format the coefficients for grads display
  tmp = replace(ARGS[a], "_obs.com", "_obs.gra")
  fpb = My.ouvre(tmp, "w")
  for b = 1:VARS
    form = @sprintf("\"d (%.8f)*air*air + (%.8f)*spd*spd + (%.8f)*(sst-273.15)*(sst-273.15) + (%.8f)*air*spd + (%.8f)*air*(sst-273.15) + (%.8f)*spd*(sst-273.15) + (%.8f)*air + (%.8f)*spd + (%.8f)*(sst-273.15) + (%.8f)\"\n", data[1,b,a], data[2,b,a], data[3,b,a], data[4,b,a], data[5,b,a], data[6,b,a], data[7,b,a], data[8,b,a], data[9,b,a], data[10,b,a])
    write(fpb, form)
  end
  close(fpb)
end

function go(val)                                                              # allow for tiny coefficients in the table
  if abs(val) < 0.00001  tmp = @sprintf("%.6f", val) ; return(tmp)  end
  if abs(val) < 0.0001   tmp = @sprintf("%.5f", val) ; return(tmp)  end
  if abs(val) < 0.001    tmp = @sprintf("%.4f", val) ; return(tmp)  end
  if abs(val) < 0.01     tmp = @sprintf("%.3f", val) ; return(tmp)  end
                         tmp = @sprintf("%.2f", val) ; return(tmp)
end

(tmp,) = split(replace(ARGS[1], "_obs.com", "_obs.trip.metricsEND"), "END")   # create tables for all analyses and each
fpc = My.ouvre(tmp, "w")                                                      # of the four triple collocation metrics
for a = 3:VARS
  form = "\n$varname $(VARN[a])\n\n"
  write(fpc, form)
  form = "\n| | air\\\*air | spd\\\*spd | sst\\\*sst | air\\\*spd | air\\\*sst | spd\\\*sst | air | spd | sst | C[10] |\n"
  write(fpc, form)
  form = "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |\n"
  write(fpc, form)
  for b = 1:argc
    contains(ARGS[b],        "cfsr") && (dirname =    "CFSR")
    contains(ARGS[b],  "erainterim") && (dirname =     "ERA")
    contains(ARGS[b],       "hoaps") && (dirname =   "HOAPS")
    contains(ARGS[b], "ifremerflux") && (dirname = "Ifremer")
    contains(ARGS[b],      "jofuro") && (dirname = "J-Ofuro")
    contains(ARGS[b],       "merra") && (dirname =   "Merra")
    contains(ARGS[b],      "oaflux") && (dirname =  "OAFlux")
    contains(ARGS[b],     "seaflux") && (dirname = "SeaFlux")
    form = @sprintf("| %7s | %8s |  %8s |  %8s |  %8s |  %8s |  %8s |  %8s |  %8s |  %8s |  %8s |\n",
      dirname, go(data[1,a,b]), go(data[2,a,b]), go(data[3,a,b]), go(data[4,a,b]), go(data[5,a,b]),
               go(data[6,a,b]), go(data[7,a,b]), go(data[8,a,b]), go(data[9,a,b]), go(data[10,a,b]))
    write(fpc, form)
  end
end
close(fpc)

exit(0)

#=
pars = ["AIRT", "LHFX", "SHFX", "SHUM", "SSTT", "WSPD"]
dirs = ["cfsr", "erainterim", "hoaps", "ifremerflux", "jofuro", "merra", "oaflux", "seaflux"]
dirn = length(dirs)

  contains(ARGS[a], "shfx") && (specstr[a,z] = @sprintf("%.0lf", specvar[a,z]))
  contains(ARGS[a], "lhfx") && (specstr[a,z] = @sprintf("%.0lf", specvar[a,z]))
  contains(ARGS[a], "wspd") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
  contains(ARGS[a], "airt") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
  contains(ARGS[a], "sstt") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
  contains(ARGS[a], "shum") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
=#
