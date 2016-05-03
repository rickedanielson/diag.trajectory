#=
 = Reformat triple collocation metric summaries as markdown
 = tables.  It is assumed that the files passed as arguments
 = are the available group of analyses for a given variable
 = (e.g., wspd) - RD April 2016.
 =#

using My
const VARS             = 4                              # number of triple collocation metrics
const COEF             = 10                             # number of polynomial coefficients
const MISS             = -9999.0                        # generic missing value

(argc,) = size(ARGS) ; if argc == 0
  print("\nUsage: jjj $(basename(@__FILE__)) all.flux.daily.locate_2.0_calib*\n\n")
  exit(1)
end

for a = 1:argc                                                                # convert each input file to markdown
  contains(ARGS[a], "shfx") && (varname = "Sensible Heat Flux")
  contains(ARGS[a], "lhfx") && (varname = "Latent Heat Flux")
  contains(ARGS[a], "wspd") && (varname = "Wind Speed")
  contains(ARGS[a], "airt") && (varname = "Air Temperature")
  contains(ARGS[a], "sstt") && (varname = "SST")
  contains(ARGS[a], "shum") && (varname = "Specific Humidity")

  fpa = My.ouvre(ARGS[a],         "r") ; lines = readlines(fpa) ; close(fpa)
  fpb = My.ouvre(ARGS[a] * ".md", "w")

  form = "\n$varname\n\n"
  write(fpb, form)
  form = "\n| | Bias | Slope | RMSE | Correlation |\n"
  write(fpb, form)
  form = "| --- | --- | --- | --- | --- |\n"
  write(fpb, form)

  (linc,) = size(lines)
  for b = 8:linc
    contains(lines[b],        "cfsr") && (dirname =    "CFSR")
    contains(lines[b],  "erainterim") && (dirname =     "ERA")
    contains(lines[b],       "hoaps") && (dirname =   "HOAPS")
    contains(lines[b], "ifremerflux") && (dirname = "Ifremer")
    contains(lines[b],      "jofuro") && (dirname = "J-Ofuro")
    contains(lines[b],       "merra") && (dirname =   "Merra")
    contains(lines[b],      "oaflux") && (dirname =  "OAFlux")
    contains(lines[b],     "seaflux") && (dirname = "SeaFlux")
    vals = split(lines[b])

    form = @sprintf("| %7s | %8s |  %8s |  %8s |  %8s |\n", dirname, vals[2], vals[3], vals[4], vals[5])
    write(fpb, form)
  end
  close(fpb)
end
exit(0)
