#=
 = Split the assembled triple collocation calibration and validation data by region,
 = where the regions are defined by the CLS MDT and the climatological ocean surface
 = current descriptions at http://oceancurrents.rsmas.miami.edu - RD February 2017.
 =#

using My, LibGEOS

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_extra.ucur.got2000_obs.comb\n\n")
  exit(1)
end

region = Array(UTF8String, 26)                                                # define regions that all subregions belong to
region[ 1] = "mediterranean"
region[ 2] = "hudson"
region[ 3] = "black"
region[ 4] = "arctic"
region[ 5] = "indiasubtropical"
region[ 6] = "kuroshio"
region[ 7] = "southpacific"
region[ 8] = "pacificequatorial"
region[ 9] = "northpacificsubtropical"
region[10] = "northpacific"
region[11] = "korea"
region[12] = "indochina"
region[13] = "arabiaindia"
region[14] = "eastafrica"
region[15] = "indiaequatorial"
region[16] = "australianewzealand"
region[17] = "agulhas"
region[18] = "brasilmalvinas"
region[19] = "antarctic"
region[20] = "southatlanticsubtropical"
region[21] = "amazoncaribbean"
region[22] = "atlanticequatorial"
region[23] = "northatlanticsubtropical"
region[24] = "gulfstream"
region[25] = "gulfofmexico"
region[26] = "northatlantic"

fpn = Array(IOStream, 0)                                                      # and open a set of output files for each region
for reg in region
  tmp = ARGS[1] * "." * reg
  fpa = My.ouvre(tmp, "w")
  push!(fpn, fpa)
end

polygons = Array(Polygon, 0)
(polynames, polypoints) = My.obtenpoly("/home/ricani/data/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS.gra/GlobCurrent_shapefiles_CNES-CLS13_RSMAS")
for (a, line) in enumerate(polypoints)
  tmp = LibGEOS.parseWKT("POLYGON(($line))")                                  # then get the subregion polygons
  push!(polygons, tmp)
end

fpi = zeros(Integer, length(polynames))                                       # and associate regions and subregions
for (a, poly) in enumerate(polynames)
  poly ==                     "gulfstream" && (fpi[a] = find(region .==               "gulfstream")[1])
  poly ==                       "kuroshio" && (fpi[a] = find(region .==                 "kuroshio")[1])
  poly ==               "agulhasextension" && (fpi[a] = find(region .==                  "agulhas")[1])
  poly ==                      "caribbean" && (fpi[a] = find(region .==          "amazoncaribbean")[1])
  poly ==                        "florida" && (fpi[a] = find(region .==               "gulfstream")[1])
  poly ==                           "loop" && (fpi[a] = find(region .==             "gulfofmexico")[1])
  poly ==                         "mexico" && (fpi[a] = find(region .==             "gulfofmexico")[1])
  poly ==                        "yucatan" && (fpi[a] = find(region .==             "gulfofmexico")[1])
  poly ==                         "canary" && (fpi[a] = find(region .== "northatlanticsubtropical")[1])
  poly ==                         "brasil" && (fpi[a] = find(region .==           "brasilmalvinas")[1])
  poly ==                         "azores" && (fpi[a] = find(region .== "northatlanticsubtropical")[1])
  poly ==                      "antillies" && (fpi[a] = find(region .== "northatlanticsubtropical")[1])
  poly ==                       "benguela" && (fpi[a] = find(region .== "southatlanticsubtropical")[1])
  poly ==                         "angola" && (fpi[a] = find(region .==       "atlanticequatorial")[1])
  poly == "atlanticnorthequatorialcounter" && (fpi[a] = find(region .==       "atlanticequatorial")[1])
  poly ==                  "northatlantic" && (fpi[a] = find(region .==               "gulfstream")[1])
  poly ==             "northatlanticdrift" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                       "labrador" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                         "guinea" && (fpi[a] = find(region .==       "atlanticequatorial")[1])
  poly ==                         "guiana" && (fpi[a] = find(region .==          "amazoncaribbean")[1])
  poly ==                    "northbrasil" && (fpi[a] = find(region .==          "amazoncaribbean")[1])
  poly ==       "northatlanticsubtropical" && (fpi[a] = find(region .== "northatlanticsubtropical")[1])
  poly ==                  "westgreenland" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==               "falklandmalvinas" && (fpi[a] = find(region .==           "brasilmalvinas")[1])
  poly ==                     "slopeshelf" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                       "portugal" && (fpi[a] = find(region .== "northatlanticsubtropical")[1])
  poly ==                         "biscay" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                         "celtic" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                          "irish" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                         "baltic" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                 "englishchannel" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                       "northsea" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                  "mediterranean" && (fpi[a] = find(region .==            "mediterranean")[1])
  poly ==                    "easticeland" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                      "norwegian" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                  "eastgreenland" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                 "agulhascoastal" && (fpi[a] = find(region .==                  "agulhas")[1])
  poly ==           "agulhasretroflection" && (fpi[a] = find(region .==                  "agulhas")[1])
  poly ==                         "hudson" && (fpi[a] = find(region .==                   "hudson")[1])
  poly ==                          "black" && (fpi[a] = find(region .==                    "black")[1])
  poly ==        "atlanticnorthequatorial" && (fpi[a] = find(region .== "northatlanticsubtropical")[1])
  poly ==        "atlanticsouthequatorial" && (fpi[a] = find(region .==       "atlanticequatorial")[1])
  poly ==             "atlanticequatorial" && (fpi[a] = find(region .==       "atlanticequatorial")[1])
  poly ==                         "arctic" && (fpi[a] = find(region .==                   "arctic")[1])
  poly ==       "southatlanticsubtropical" && (fpi[a] = find(region .== "southatlanticsubtropical")[1])
  poly ==                    "spitzbergen" && (fpi[a] = find(region .==            "northatlantic")[1])
  poly ==                 "eastaustralian" && (fpi[a] = find(region .==      "australianewzealand")[1])
  poly ==                        "leeuwin" && (fpi[a] = find(region .==      "australianewzealand")[1])
  poly ==                      "antarctic" && (fpi[a] = find(region .==                "antarctic")[1])
  poly ==                          "drake" && (fpi[a] = find(region .==                "antarctic")[1])
  poly ==                        "okhotsk" && (fpi[a] = find(region .==             "northpacific")[1])
  poly ==                        "oyashio" && (fpi[a] = find(region .==             "northpacific")[1])
  poly ==              "kuroshioextension" && (fpi[a] = find(region .==                 "kuroshio")[1])
  poly ==                  "eastkamchatka" && (fpi[a] = find(region .==             "northpacific")[1])
  poly ==               "eastindiacoastal" && (fpi[a] = find(region .==              "arabiaindia")[1])
  poly ==               "westindiacoastal" && (fpi[a] = find(region .==              "arabiaindia")[1])
  poly ==                           "oman" && (fpi[a] = find(region .==              "arabiaindia")[1])
  poly ==                        "persian" && (fpi[a] = find(region .==              "arabiaindia")[1])
  poly ==                         "redsea" && (fpi[a] = find(region .==               "eastafrica")[1])
  poly ==                         "somali" && (fpi[a] = find(region .==               "eastafrica")[1])
  poly ==                 "eastmadagascar" && (fpi[a] = find(region .==                  "agulhas")[1])
  poly ==                     "mozambique" && (fpi[a] = find(region .==                  "agulhas")[1])
  poly ==                    "eastafrican" && (fpi[a] = find(region .==               "eastafrica")[1])
  poly ==               "indiasubtropical" && (fpi[a] = find(region .==         "indiasubtropical")[1])
  poly ==                 "northaustralia" && (fpi[a] = find(region .==      "australianewzealand")[1])
  poly ==           "indiasouthequatorial" && (fpi[a] = find(region .==          "indiaequatorial")[1])
  poly ==         "indiaequatorialcounter" && (fpi[a] = find(region .==          "indiaequatorial")[1])
  poly ==                     "indonesian" && (fpi[a] = find(region .==                "indochina")[1])
  poly ==                     "southchina" && (fpi[a] = find(region .==                "indochina")[1])
  poly ==                      "yellowsea" && (fpi[a] = find(region .==                    "korea")[1])
  poly ==                     "seaofjapan" && (fpi[a] = find(region .==                    "korea")[1])
  poly ==                    "tehuantepec" && (fpi[a] = find(region .==        "pacificequatorial")[1])
  poly ==                     "california" && (fpi[a] = find(region .==  "northpacificsubtropical")[1])
  poly ==                         "bering" && (fpi[a] = find(region .==             "northpacific")[1])
  poly ==                         "alaska" && (fpi[a] = find(region .==             "northpacific")[1])
  poly ==                   "alaskastream" && (fpi[a] = find(region .==             "northpacific")[1])
  poly ==              "southpacificdrift" && (fpi[a] = find(region .==             "southpacific")[1])
  poly ==        "southpacificsubtropical" && (fpi[a] = find(region .==             "southpacific")[1])
  poly ==                         "tasman" && (fpi[a] = find(region .==      "australianewzealand")[1])
  poly ==                       "westland" && (fpi[a] = find(region .==      "australianewzealand")[1])
  poly ==                      "southland" && (fpi[a] = find(region .==      "australianewzealand")[1])
  poly ==                       "eastcape" && (fpi[a] = find(region .==      "australianewzealand")[1])
  poly ==                   "eastauckland" && (fpi[a] = find(region .==      "australianewzealand")[1])
  poly ==              "northpacificdrift" && (fpi[a] = find(region .==             "northpacific")[1])
  poly ==        "northpacificsubtropical" && (fpi[a] = find(region .==  "northpacificsubtropical")[1])
  poly ==                           "peru" && (fpi[a] = find(region .==             "southpacific")[1])
  poly ==         "pacificsouthequatorial" && (fpi[a] = find(region .==        "pacificequatorial")[1])
  poly ==              "pacificequatorial" && (fpi[a] = find(region .==        "pacificequatorial")[1])
  poly ==  "pacificnorthequatorialcounter" && (fpi[a] = find(region .==        "pacificequatorial")[1])
  poly ==         "pacificnorthequatorial" && (fpi[a] = find(region .==        "pacificequatorial")[1])
  poly ==                         "bengal" && (fpi[a] = find(region .==              "arabiaindia")[1])
end

fpa = My.ouvre(ARGS[1], "r")                                                  # read all triple collocations in one gulp
lines = readlines(fpa) ; close(fpa)
for (a, line) in enumerate(lines)
  tmp = split(line)
  lat = float(tmp[2])
  lon = float(tmp[3]) ; lon > 180 && (lon -= 360)
  loc = LibGEOS.parseWKT("POINT($lon $lat)")
  ind = 0
  for (c, polyc) in enumerate(polygons)                                       # then find the bounding polygon for each line
    if LibGEOS.contains(polyc, loc)
      ind = c
      break
    end
  end
  if ind == 0                                                                 # (allowing for a 180-deg shift in longitude)
    lon += 360
    loc = LibGEOS.parseWKT("POINT($lon $lat)")
    for (c, polyc) in enumerate(polygons)
      if LibGEOS.contains(polyc, loc)
        ind = c
        break
      end
    end
  end
  write(fpn[fpi[ind]], line)                                                  # and write each line to the regional files
end
exit(0)
