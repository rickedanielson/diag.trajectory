#=
 = Reformat triple collocation metric summaries as markdown
 = tables.  It is assumed that the files passed as arguments
 = are the available group of analyses for a given variable
 = (e.g., ucur) - RD June 2016.
 =#

using My
const VARS             = 4                              # number of triple collocation metrics
const COEF             = 10                             # number of polynomial coefficients
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) == 0
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comc.cali.pair.morph\n\n")
  exit(1)
end
dirname = "GlobCurrent"
contains(ARGS[1], "morph") && (dirname *= " with Gaussian anamorphosis")

fpa = My.ouvre(ARGS[1] * ".md", "w")                                          # use the first file to name the markdown file
form = "\n$dirname\n\n"
write(fpa, form)
form = "\n| | Bias | Slope | RMSE | Correlation |\n"
write(fpa, form)
form = "| --- | --- | --- | --- | --- |\n"
write(fpa, form)

for a = 1:argc                                                                # and convert each input file to markdown
  contains(ARGS[a], "ucur") && (varname = "Zonal Current")
  contains(ARGS[a], "vcur") && (varname = "Meridional Current")
  fpb = My.ouvre(ARGS[a], "r") ; lines = readlines(fpb) ; close(fpb)
  vals = split(lines[4])
  valz = split(lines[8])

  stra = abs(float(vals[2]) - 0) < abs(float(valz[2]) - 0) ? "**" * vals[2] * "/" * valz[2] * "**" : vals[2] * "/" * valz[2]
  strb = abs(float(vals[3]) - 1) < abs(float(valz[3]) - 1) ? "**" * vals[3] * "/" * valz[3] * "**" : vals[3] * "/" * valz[3]
  strc =     float(vals[4])      <     float(valz[4])      ? "**" * vals[4] * "/" * valz[4] * "**" : vals[4] * "/" * valz[4]
  strd =     float(vals[5])      >     float(valz[5])      ? "**" * vals[5] * "/" * valz[5] * "**" : vals[5] * "/" * valz[5]
  form = @sprintf("| %7s | %18s |  %18s |  %18s |  %18s |\n", varname, stra, strb, strc, strd)
  write(fpa, form)
end

close(fpa)
exit(0)
