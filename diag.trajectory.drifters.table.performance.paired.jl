#=
 = Reformat triple collocation metric summaries as markdown
 = tables.  It is assumed that the files passed as arguments
 = are the available group of analyses for a given variable
 = (e.g., ucur) - RD June, November 2016.
 =#

using My
const ALPH             = 2                              # error model x = ALPH + BETA * truth + error
const BETA             = 3                              # error model x = ALPH + BETA * truth + error
const SIGM             = 4                              # triple coll RMSE
const CORR             = 5                              # triple coll correlation coefficient
const VARS             = 4                              # number of triple collocation metrics
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) == 0
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comc.cali.pair.morph\n\n")
  exit(1)
end
dirname = "GlobCurrent"
contains(ARGS[1], "morph") && (dirname *= " with Gaussian anamorphosis")

fpa = My.ouvre(ARGS[1] * ".md", "w")                                          # use the first to name the markdown file
form = "\nBegin Table\n\n" ; write(fpa, form)
form = "|Group-1 Group-2 Recal-1 Recal-2|Obs RMSE|Obs Corr|Anlys Bias Addit|Anlys Bias Multi|Anlys RMSE|Anlys Corr|\n" ; write(fpa, form)
form = "|-|-|-|-|-|-|-|\n" ; write(fpa, form)
form  = @sprintf("|%s|||||||\n", dirname) ; write(fpa, form)

for a = 1:argc                                                                # and convert each input file to markdown
  contains(ARGS[a], "ucur") && (varname =      "Zonal Current")
  contains(ARGS[a], "vcur") && (varname = "Meridional Current")
  contains(ARGS[a], "wcur") && (varname =      "Speed Current")

  dirname = "GlobCurrent"
  contains(ARGS[a], "morph") && (dirname *= " with Gaussian anamorphosis")
  fpb = My.ouvre(ARGS[a], "r") ; lines = readlines(fpb) ; close(fpb)
  vala = split(lines[ 4]) ; valb = split(lines[ 5]) ; obsa = split(lines[ 6]) ; obsb = split(lines[ 7])
  valc = split(lines[14]) ; vald = split(lines[15]) ; obsc = split(lines[16]) ; obsd = split(lines[17])

  stra  =     float(obsa[SIGM])      <     float(obsc[SIGM])      ?  "**" * obsa[SIGM] * "**" :  "" * obsa[SIGM]
  strb  =     float(obsa[CORR])      >     float(obsc[CORR])      ?  "**" * obsa[CORR] * "**" :  "" * obsa[CORR]
  strc  = abs(float(vala[ALPH]) - 0) < abs(float(valc[ALPH]) - 0) ?  "**" * vala[ALPH] * "**" :  "" * vala[ALPH]
  strd  = abs(float(vala[BETA]) - 1) < abs(float(valc[BETA]) - 1) ?  "**" * vala[BETA] * "**" :  "" * vala[BETA]
  stre  =     float(vala[SIGM])      <     float(valc[SIGM])      ?  "**" * vala[SIGM] * "**" :  "" * vala[SIGM]
  strf  =     float(vala[CORR])      >     float(valc[CORR])      ?  "**" * vala[CORR] * "**" :  "" * vala[CORR]
  stra *=     float(obsb[SIGM])      <     float(obsd[SIGM])      ? " **" * obsb[SIGM] * "**" : " " * obsb[SIGM]
  strb *=     float(obsb[CORR])      >     float(obsd[CORR])      ? " **" * obsb[CORR] * "**" : " " * obsb[CORR]
  strc *= abs(float(valb[ALPH]) - 0) < abs(float(vald[ALPH]) - 0) ? " **" * valb[ALPH] * "**" : " " * valb[ALPH]
  strd *= abs(float(valb[BETA]) - 1) < abs(float(vald[BETA]) - 1) ? " **" * valb[BETA] * "**" : " " * valb[BETA]
  stre *=     float(valb[SIGM])      <     float(vald[SIGM])      ? " **" * valb[SIGM] * "**" : " " * valb[SIGM]
  strf *=     float(valb[CORR])      >     float(vald[CORR])      ? " **" * valb[CORR] * "**" : " " * valb[CORR]
  strc *= abs(float(vala[ALPH]) - 0) < abs(float(valc[ALPH]) - 0) ? " **" * valc[ALPH] * "**" : " " * valc[ALPH]
  strd *= abs(float(vala[BETA]) - 1) < abs(float(valc[BETA]) - 1) ? " **" * valc[BETA] * "**" : " " * valc[BETA]
  stre *=     float(vala[SIGM])      <     float(valc[SIGM])      ? " **" * valc[SIGM] * "**" : " " * valc[SIGM]
  strc *= abs(float(valb[ALPH]) - 0) < abs(float(vald[ALPH]) - 0) ? " **" * vald[ALPH] * "**" : " " * vald[ALPH]
  strd *= abs(float(valb[BETA]) - 1) < abs(float(vald[BETA]) - 1) ? " **" * vald[BETA] * "**" : " " * vald[BETA]
  stre *=     float(valb[SIGM])      <     float(vald[SIGM])      ? " **" * vald[SIGM] * "**" : " " * vald[SIGM]

  if float(obsa[SIGM]) != float(obsc[SIGM])
    stra *=   float(obsa[SIGM])      <     float(obsc[SIGM])      ? " **" * obsc[SIGM] * "**" : " " * obsc[SIGM]
  end
  if float(obsa[CORR]) != float(obsc[CORR])
    strb *=   float(obsa[CORR])      >     float(obsc[CORR])      ? " **" * obsc[CORR] * "**" : " " * obsc[CORR]
  end
  if float(vala[CORR]) != float(valc[CORR])
    strf *=   float(vala[CORR])      >     float(valc[CORR])      ? " **" * valc[CORR] * "**" : " " * valc[CORR]
  end
  if float(obsb[SIGM]) != float(obsd[SIGM])
    stra *=   float(obsb[SIGM])      <     float(obsd[SIGM])      ? " **" * obsd[SIGM] * "**" : " " * obsd[SIGM]
  end
  if float(obsb[CORR]) != float(obsd[CORR])
    strb *=   float(obsb[CORR])      >     float(obsd[CORR])      ? " **" * obsd[CORR] * "**" : " " * obsd[CORR]
  end
  if float(valb[CORR]) != float(vald[CORR])
    strf *=   float(valb[CORR])      >     float(vald[CORR])      ? " **" * vald[CORR] * "**" : " " * vald[CORR]
  end
  form  = @sprintf("|%s|%s|%s|%s|%s|%s|%s|\n", varname, stra, strb, strc, strd, stre, strf) ; write(fpa, form)
end

form = "\nEnd Table\n" ; write(fpa, form)
close(fpa)
exit(0)

#=
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
=#
