```shell
`
# get a copy
git clone git@github.com:rickedanielson/diag.trajectory.git

# requirements on ubuntu 14.04 (local) and at Ifremer (12.04)
julia (http://julialang.org/)
GNU parallel (http://www.gnu.org/software/parallel/)
alias wrkt 'cd ~/work/workt; ls'

# copy data then isolate and plot the location of drifters not in the CNES-CLS-2013 MDT
wrkt ; cd all ; cp /home/cercache/project/globcurrent/data/third-party/insitu/drifters-rio/* .
       mkdir plot.available plot.histogr plot.locate plot.scatter
       jjj diag.trajectory.drifters.nonmdt.jl buoydata_1993_2014_drogON.asc
       jjj diag.trajectory.drifters.locate.jl buoydata_1993_2014_drogON.asc.nonmdt
       xvfb-run -a grads -blc "diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate"
       mv plot.trajectory.drifters.dots*locate*png plot.locate

# split the drifter data by location into calibration and validation groups
wrkt ; cd all
       jjj diag.trajectory.drifters.split.jl buoydata_1993_2014_drogON.asc.nonmdt.locate
       jjj diag.trajectory.drifters.split.jl buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid
       mv buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid           buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib_remainder
       mv buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_2.0_calib buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid
       mv buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_2.0_valid buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder
       xvfb-run -a grads -blc "diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib"
       xvfb-run -a grads -blc "diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid"
       xvfb-run -a grads -blc "diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib_remainder"
       xvfb-run -a grads -blc "diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder"
       mv plot.trajectory.drifters.dots*locate*png plot.locate

# further split the in situ cal/val observations by location and store files in an insitu dir
wrkt ; mkdir insitu
       sort -k2,2 -k3,3 -k1,1 all/buoydata_1993_2014_drogON.asc.nonmdt > buoydata_1993_2014_drogON.asc.nonmdt.sort
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.split.location.jl ::: all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_?ali? ::: buoydata_1993_2014_drogON.asc.nonmdt.sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       cd insitu ; ls -1 insitu..????.???..????.??? > z.list ; cd .. ; wc insitu/z.list
       rm commands buoydata_1993_2014_drogON.asc.nonmdt.sort

# create local links to all analysis data files and example ncdumps too
wrkt ; mkdir v2.0_global_025_deg_ekman_15m v2.0_global_025_deg_ekman_hs v2.0_global_025_deg_geostrophic v2.0_global_025_deg_total_15m v2.0_global_025_deg_total_hs
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_ekman_15m   ; jjj diag.trajectory.drifters.links.jl /home/cercache/project/globcurrent/data/globcurrent/v2.0/global_025_deg/ekman_15m
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_ekman_hs    ; jjj diag.trajectory.drifters.links.jl /home/cercache/project/globcurrent/data/globcurrent/v2.0/global_025_deg/ekman_hs
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_geostrophic ; jjj diag.trajectory.drifters.links.jl /home/cercache/project/globcurrent/data/globcurrent/v2.0/global_025_deg/geostrophic
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_total_15m   ; jjj diag.trajectory.drifters.links.jl /home/cercache/project/globcurrent/data/globcurrent/v2.0/global_025_deg/total_15m
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_total_hs    ; jjj diag.trajectory.drifters.links.jl /home/cercache/project/globcurrent/data/globcurrent/v2.0/global_025_deg/total_hs
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_ekman_15m   ; mkdir old ; mv 200* old ; mv 201[01]* old ; mv 20120[1-8]* old
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_ekman_hs    ; mkdir old ; mv 200* old ; mv 201[01]* old ; mv 20120[1-8]* old
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_geostrophic ; mkdir old ; mv 200* old ; mv 201[01]* old ; mv 20120[1-8]* old
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_total_15m   ; mkdir old ; mv 200* old ; mv 201[01]* old ; mv 20120[1-8]* old
       cd /home/cercache/users/rdaniels/work/workt/v2.0_global_025_deg_total_hs    ; mkdir old ; mv 200* old ; mv 201[01]* old ; mv 20120[1-8]* old
wrkt ; mkdir ncdump
       ncdump   v2.0_global_025_deg_ekman_15m/20120901000000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v02.0-fv01.0.nc > ncdump/v2.0_global_025_deg_ekman_15m
       ncdump    v2.0_global_025_deg_ekman_hs/20120901000000-GLOBCURRENT-L4-CURekm_hs-ERAWS_EEM-v02.0-fv01.0.nc  > ncdump/v2.0_global_025_deg_ekman_hs
       ncdump v2.0_global_025_deg_geostrophic/20120901000000-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v02.0-fv01.0.nc     > ncdump/v2.0_global_025_deg_geostrophic
       ncdump   v2.0_global_025_deg_total_15m/20120901000000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v02.0-fv01.0.nc   > ncdump/v2.0_global_025_deg_total_15m
       ncdump    v2.0_global_025_deg_total_hs/20120901000000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc    > ncdump/v2.0_global_025_deg_total_hs

# assemble a large third dataset for analysis evaulation (with respect to the 2.0_valid_remainder in situ obs)
wrkt ; cd all ; jjj eulerian.evaluation.assemble.insitu.jl buoydata_1993_2014_drogON.asc.nonmdt  buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder
       cd .. ; sort all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs    > buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs.sort
       split -l 171000  buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs.sort buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/eulerian.evaluation.assemble.analyses.jl ::: buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs?? ::: v2.0_global_025_deg_ekman_15m v2.0_global_025_deg_ekman_hs v2.0_global_025_deg_geostrophic v2.0_global_025_deg_total_15m v2.0_global_025_deg_total_hs | grep buoy | sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs??.v2.0_global_025_deg_ekman_15m   | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs.v2.0_global_025_deg_ekman_15m
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs??.v2.0_global_025_deg_ekman_hs    | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs.v2.0_global_025_deg_ekman_hs
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs??.v2.0_global_025_deg_geostrophic | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs.v2.0_global_025_deg_geostrophic
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs??.v2.0_global_025_deg_total_15m   | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs.v2.0_global_025_deg_total_15m
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs??.v2.0_global_025_deg_total_hs    | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs.v2.0_global_025_deg_total_hs
       wc all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs.*
       rm commands buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs*

# perform an initial eight-analysis evaulation, but without calibration (versus the 2.0_valid_remainder obs)
wrkt ; parallel --dry-run /home1/homedir1/perso/rdaniels/bin/eulerian.evaluation.versus.insitu.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid_remainder_obs ::: ucur vcur | grep buoy | sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       rm commands ; cd all ; cat *ucur.summ *vcur.summ

# return to the cal/val locations to create analysis timeseries (some locations missing too much of this timeseries are later ignored)
wrkt ; sort      all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib    > buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.sort
       split -l 1000 buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.sort buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.sort
       sort      all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid    > buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.sort
       split -l 1000 buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.sort buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.sort
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.jl ::: buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_?ali?.sort?? ::: v2.0_global_025_deg_ekman_15m v2.0_global_025_deg_ekman_hs v2.0_global_025_deg_geostrophic v2.0_global_025_deg_total_15m v2.0_global_025_deg_total_hs | grep buoy | sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       rm commands buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_?ali?.sor*

# verify that each subdir contains the expected number of files (e.g., 4339 + 3619 = 7958 files with 3408 dates)
wrkt ; cd v2.0_global_025_deg_ekman_15m   ; ls -1   v2.0_global_025_deg_ekman_15m..????.???..????.??? > z.list ; split -l 1000 z.list z.list ; cd ..
       cd v2.0_global_025_deg_ekman_hs    ; ls -1    v2.0_global_025_deg_ekman_hs..????.???..????.??? > z.list ; split -l 1000 z.list z.list ; cd ..
       cd v2.0_global_025_deg_geostrophic ; ls -1 v2.0_global_025_deg_geostrophic..????.???..????.??? > z.list ; split -l 1000 z.list z.list ; cd ..
       cd v2.0_global_025_deg_total_15m   ; ls -1   v2.0_global_025_deg_total_15m..????.???..????.??? > z.list ; split -l 1000 z.list z.list ; cd ..
       cd v2.0_global_025_deg_total_hs    ; ls -1    v2.0_global_025_deg_total_hs..????.???..????.??? > z.list ; split -l 1000 z.list z.list ; cd ..
       wc *[a-z]/z.list

# plot examples of temporal coverage by all analyses (include in situ) at a few locations
wrkt ; xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.available.jl ....44.875...-45.125
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.available.jl ....55.625...-10.875
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.available.jl ....47.375....-9.375
       mv plot.avail.* all/plot.available

# plot histograms
wrkt ; echo /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.histogram.jl z.list > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia
       rm commands
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.histoplot.jl
       gzip histogr*dat ; mv histogr* all/plot.histogr

# create the forward and backward extrapolated timeseries, but just for the geostrophic component
wrkt ; cd v2.0_global_025_deg_geostrophic ; ls z.list?? ; cd ..
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.jl v2.0_global_025_deg_geostrophic ::: z.listaa z.listab z.listac z.listad z.listae z.listaf z.listag z.listah | grep drift | sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       rm commands

# plot extrapolation histograms (forward and backward versus the actual values for assessment of bias in the extrapolation method)
wrkt ; nohup julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histogram.jl z.list > xcom &
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_ekman_15m
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_ekman_hs
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_geostrophic
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_total_15m
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_total_hs
       gzip extrapolated.histogr*dat ; mv extrapolated.histogr* all/plot.histogr



# store drifters (lasting more than a day) as individual files
wrkt
mkdir buoydata_1993_2012_drogON buoydata_1993_2012_drogOFF
jj diag.trajectory.drifters.divide.jl buoydata_1993_2012_drogON.asc
jj diag.trajectory.drifters.check.jl
mv [12]* buoydata_1993_2012_drogON
jj diag.trajectory.drifters.divide.jl buoydata_1993_2012_drogOFF.asc
jj diag.trajectory.drifters.check.jl
mv [12]* buoydata_1993_2012_drogOFF

# compute simulated trajectories (and check for any missing)
cd buoydata_1993_2012_drogON
cd buoydata_1993_2012_drogOFF
ls -1 201?-??-??-????.???????? | /home5/begmeil/tools/gogolist/bin/gogolist.py -e diag.trajectory.drifters --mem=2000mb
jj diag.trajectory.drifters.check.jl

# animate the drifters
cd buoydata_1993_2012_drogON
cd buoydata_1993_2012_drogOFF
ls -1 ????-??-??-????.???????? | /home5/begmeil/tools/gogolist/bin/gogolist.py -e "julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.plot.jl" --mem=2000mb

# and divide by region
jjo diag.trajectory.drifters.divide.again.jl buoydata_1993_2012_drogON  ~/data/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS_valid.shp
jjo diag.trajectory.drifters.divide.again.jl buoydata_1993_2012_drogOFF ~/data/mdt/GlobCurrent_shapefiles_CNES-CLS13_RSMAS_valid.shp

# then compile validation stats
jj 





# create the subdomain files (at Ifremer)
wrkt ; julia
STEMA = "/home1/homedir1/perso/rdaniels"
try ; global STEMA ; STEMA = ENV["STEMA"] ; end
push!(LOAD_PATH, STEMA * "/prog/code.julia.modules/")
gcfiles = filter(x -> (contains(x, "ALT_OI-v01.0-fv01.0.nc")), readdir("linkt"))
gcfiles = filter(x -> (contains(x, "ERAWS_EEM-v01.0-fv01.0.nc")), readdir("linkv"))
fpa = open("exec", "w")
for file in gcfiles
  write(fpa, "ncks -d lat,65.,75. -d lon,-10.,20. /home/cercache/users/rdaniels/work/workt/linkv/$file -o norskhavet_$file\n")
end
close(fpa)

# create template for Lofoten basin GlobCurrent timeseries
date.extend 2010-01-01-00 24 1095 > norskhavet_2010-2012.geostrophic
nc.template.norskhavet.surf.current norskhavet_2010-2012.geostrophic
nc.transfer.norskhavet.surf.current norskhavet_2010-2012.geostrophic norskhavet_20121221000000-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v01.0-fv01.0.nc norskhavet_2010-2012.geostrophic.surf.current.nc
cp norskhavet_2010-2012.geostrophic.surf.current.nc norskhavet_2010-2012.geostrophic.surf.current.low.030.nc
cp norskhavet_2010-2012.geostrophic.surf.current.nc norskhavet_2010-2012.geostrophic.surf.current.hig.030.nc
cp norskhavet_2010-2012.geostrophic.surf.current.nc norskhavet_2010-2012.geostrophic.surf.current.low.060.nc
cp norskhavet_2010-2012.geostrophic.surf.current.nc norskhavet_2010-2012.geostrophic.surf.current.hig.060.nc
cp norskhavet_2010-2012.geostrophic.surf.current.nc norskhavet_2010-2012.geostrophic.surf.current.low.090.nc
cp norskhavet_2010-2012.geostrophic.surf.current.nc norskhavet_2010-2012.geostrophic.surf.current.hig.090.nc
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012.geostrophic .surf.current 030 u
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012.geostrophic .surf.current 030 v
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012.geostrophic .surf.current 060 u
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012.geostrophic .surf.current 060 v
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012.geostrophic .surf.current 090 u
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012.geostrophic .surf.current 090 v

date.extend 2010-01-01-03  3 8766 > norskhavet_2010-2012
nc.template.norskhavet.surf.current norskhavet_2010-2012
nc.transfer.norskhavet.surf.current norskhavet_2010-2012 norskhavet_20100101030000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v01.0-fv01.0.nc norskhavet_2010-2012.surf.current.nc
cp norskhavet_2010-2012.surf.current.nc norskhavet_2010-2012.surf.current.low.030.nc
cp norskhavet_2010-2012.surf.current.nc norskhavet_2010-2012.surf.current.hig.030.nc
cp norskhavet_2010-2012.surf.current.nc norskhavet_2010-2012.surf.current.low.060.nc
cp norskhavet_2010-2012.surf.current.nc norskhavet_2010-2012.surf.current.hig.060.nc
cp norskhavet_2010-2012.surf.current.nc norskhavet_2010-2012.surf.current.low.090.nc
cp norskhavet_2010-2012.surf.current.nc norskhavet_2010-2012.surf.current.hig.090.nc
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012 .surf.current 030 ucomb
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012 .surf.current 030 vcomb
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012 .surf.current 060 ucomb
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012 .surf.current 060 vcomb
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012 .surf.current 090 ucomb
diag.trajectory.drifters.norskhavet.fft norskhavet_2010-2012 .surf.current 090 vcomb
cp norskhavet_2010-2012.surf*.nc /net/sverdrup-2/vol/Service1/FTPRoot/pub/rick

date.extend 2010-01-01-03  3 8766 > norskhavet_2010-2012.ekman15
nc.template.norskhavet.surf.current norskhavet_2010-2012.ekman15
nc.transfer.norskhavet.surf.current norskhavet_2010-2012.ekman15 norskhavet_20121214150000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v01.0-fv01.0.nc norskhavet_2010-2012.ekman15.surf.current.nc
date.extend 2010-01-01-03  3 8766 > norskhavet_2010-2012.ekman00
nc.template.norskhavet.surf.current norskhavet_2010-2012.ekman00
nc.transfer.norskhavet.surf.current norskhavet_2010-2012.ekman00 norskhavet_20110701060000-GLOBCURRENT-L4-CURekm_hs-ERAWS_EEM-v01.0-fv01.0.nc  norskhavet_2010-2012.ekman00.surf.current.nc

# restore Lofoten drifters (lasting more than a day) as individual files
wrkt ; mkdir lofoten ; cd lofoten
mkdir buoydata_1993_2012_drogON        buoydata_1993_2012_drogOFF
cp ../buoydata_1993_2012_drogON.asc ../buoydata_1993_2012_drogOFF.asc .
jj diag.trajectory.drifters.divide.jl buoydata_1993_2012_drogON.asc
jj diag.trajectory.drifters.check.jl
mv [12]* buoydata_1993_2012_drogON
jj diag.trajectory.drifters.divide.jl buoydata_1993_2012_drogOFF.asc
jj diag.trajectory.drifters.check.jl
mv [12]* buoydata_1993_2012_drogOFF

# divide the drifter files by region
wrkt ; cd lofoten/buoydata_1993_2012_drogON
jjo diag.trajectory.drifters.divide.lofoten.jl
wrkt ; cd lofoten/buoydata_1993_2012_drogOFF
jjo diag.trajectory.drifters.divide.lofoten.jl

# compute simulated trajectories (and check for any missing)
cd buoydata_1993_2012_drogON
cd buoydata_1993_2012_drogOFF
ls -1 201?-??-??-????.???????? | /home5/begmeil/tools/gogolist/bin/gogolist.py -e diag.trajectory.drifters --mem=2000mb
jj diag.trajectory.drifters.check.jl






convert -delay 0 -loop 0 plot.globcurrent.agulhas.* GlobCurrent_2012-09_Agulhas.gif
convert -delay 0 -loop 0 plot.globcurrent.northatl.* GlobCurrent_2012-09_NorthAtl.gif

date.extend 2012-05-05-00 3 3111 > datelist ; tail datelist


finf "*20120920*" > zxc ; vi zxc
:1,33s/\.\//cp \.\//
:1,33s/\.nc/\.nc lum.test/


wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_010_deg/ekman_15m/'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_010_deg/ekman_hs/'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_010_deg/geostrophic/'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_010_deg/stokes_drift/'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_010_deg/tidal/'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_010_deg/total_15m/'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_010_deg/total_hs/'



wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_010_deg/stokes_drift/2012/245/*'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_025_deg/v0.1/ekman_15m/2012/245/*'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_025_deg/v0.1/ekman_hs/2012/245/*'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_025_deg/v0.1/geostrophic/2012/245/*'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_025_deg/v0.1/total_15m/2012/245/*'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_025_deg/v0.1/total_15m/2012/245/*'
wget 'ftp://gg1f3e8:xG3jZhT9@eftp.ifremer.fr:21/data/globcurrent/global_025_deg/v0.1/total_hs/2012/245/*'


jj diag.trajectory.drifters.divide.jl buoydata_1993_2012_drogON.asc
mv [12]* buoydata_1993_2012_drogON

jj diag.trajectory.drifters.divide.jl buoydata_1993_2012_drogOFF.asc
mv [12]* buoydata_1993_2012_drogOFF.asc

cd buoydata_1993_2012_drogON
doit2 2012-04-12-0000.00109460.traj*      (nc.template.track.drifters)
jj diag.trajectory.drifters.plot.jl 2012-04-12-0000.00109460


~/work/workt/buoydata_1993_2012_drogON > jj diag.trajectory.drifters.Agulhas.jl | grep -v reading
jj diag.trajectory.drifters.plot.jl 1994-09-21-0000.07722594
jj diag.trajectory.drifters.plot.jl 1994-09-26-0600.07722587
jj diag.trajectory.drifters.plot.jl 1994-10-01-0600.09422042
jj diag.trajectory.drifters.plot.jl 1994-10-01-1800.09422041
jj diag.trajectory.drifters.plot.jl 1994-10-01-1800.09422098
jj diag.trajectory.drifters.plot.jl 1994-10-02-0000.07722589
jj diag.trajectory.drifters.plot.jl 1994-10-03-0000.07722583
jj diag.trajectory.drifters.plot.jl 1994-10-03-0000.09422095
jj diag.trajectory.drifters.plot.jl 1994-10-03-1200.07722593
jj diag.trajectory.drifters.plot.jl 1994-11-09-0000.09422091
jj diag.trajectory.drifters.plot.jl 1994-11-20-0600.09422100
jj diag.trajectory.drifters.plot.jl 1995-03-02-1200.09421934
jj diag.trajectory.drifters.plot.jl 1995-03-04-1800.09421910
jj diag.trajectory.drifters.plot.jl 1995-03-23-1800.09421943
jj diag.trajectory.drifters.plot.jl 1995-05-12-0600.09422069
jj diag.trajectory.drifters.plot.jl 1995-05-13-1200.09422070
jj diag.trajectory.drifters.plot.jl 1995-06-06-1800.09422076
jj diag.trajectory.drifters.plot.jl 1995-07-02-1800.09421897
jj diag.trajectory.drifters.plot.jl 1995-11-13-0000.09422084
jj diag.trajectory.drifters.plot.jl 1995-11-14-0600.09422074
jj diag.trajectory.drifters.plot.jl 1995-11-15-1200.09422083
jj diag.trajectory.drifters.plot.jl 1995-12-01-0600.09525851
jj diag.trajectory.drifters.plot.jl 1995-12-01-0600.09525852
jj diag.trajectory.drifters.plot.jl 1995-12-21-1800.09525862
jj diag.trajectory.drifters.plot.jl 1995-12-27-0600.09525861
jj diag.trajectory.drifters.plot.jl 1996-03-01-1800.09525867
jj diag.trajectory.drifters.plot.jl 1996-04-01-0600.09422068
jj diag.trajectory.drifters.plot.jl 1996-04-02-0000.09422067
jj diag.trajectory.drifters.plot.jl 1996-04-07-1800.09421924
jj diag.trajectory.drifters.plot.jl 1996-04-21-1200.09421931
jj diag.trajectory.drifters.plot.jl 1996-04-26-0000.09525826
jj diag.trajectory.drifters.plot.jl 1996-04-27-0600.09525833
jj diag.trajectory.drifters.plot.jl 1996-04-27-1800.09525831
jj diag.trajectory.drifters.plot.jl 1996-08-09-0000.09526393
jj diag.trajectory.drifters.plot.jl 1996-11-01-0600.09616587
jj diag.trajectory.drifters.plot.jl 1996-12-08-0000.09524493
jj diag.trajectory.drifters.plot.jl 1996-12-08-0000.09524501
jj diag.trajectory.drifters.plot.jl 1996-12-10-1200.09524502
jj diag.trajectory.drifters.plot.jl 1996-12-11-0600.09616583
jj diag.trajectory.drifters.plot.jl 1996-12-12-0600.09525837
jj diag.trajectory.drifters.plot.jl 1996-12-14-1800.09616596
jj diag.trajectory.drifters.plot.jl 1997-01-10-1800.09616676
jj diag.trajectory.drifters.plot.jl 1997-01-11-0000.09616672
jj diag.trajectory.drifters.plot.jl 1997-01-11-1800.09524459
jj diag.trajectory.drifters.plot.jl 1997-01-11-1800.09616677
jj diag.trajectory.drifters.plot.jl 1997-01-12-0600.09616678
jj diag.trajectory.drifters.plot.jl 1997-01-12-1800.09524469
jj diag.trajectory.drifters.plot.jl 1997-02-01-1800.09619431
jj diag.trajectory.drifters.plot.jl 1997-02-06-1800.09619412
jj diag.trajectory.drifters.plot.jl 1997-02-07-0000.09619429
jj diag.trajectory.drifters.plot.jl 1997-02-08-0600.09619426
jj diag.trajectory.drifters.plot.jl 1997-02-19-0600.09619407
jj diag.trajectory.drifters.plot.jl 1997-04-27-0000.09619423
jj diag.trajectory.drifters.plot.jl 1997-04-29-0600.09619428
jj diag.trajectory.drifters.plot.jl 1997-05-26-1800.09619415
jj diag.trajectory.drifters.plot.jl 1997-05-27-1800.09619410
jj diag.trajectory.drifters.plot.jl 1997-05-30-0600.09619409
jj diag.trajectory.drifters.plot.jl 1997-10-30-0000.09525824
jj diag.trajectory.drifters.plot.jl 1997-12-06-0600.09619816
jj diag.trajectory.drifters.plot.jl 1997-12-10-1800.09729821
jj diag.trajectory.drifters.plot.jl 1997-12-11-1800.09729825
jj diag.trajectory.drifters.plot.jl 1998-01-07-0000.09619818
jj diag.trajectory.drifters.plot.jl 1998-01-23-0600.09721481
jj diag.trajectory.drifters.plot.jl 1998-11-10-1800.09729785
jj diag.trajectory.drifters.plot.jl 1998-12-07-1800.09729717
jj diag.trajectory.drifters.plot.jl 1999-01-27-0000.09628902
jj diag.trajectory.drifters.plot.jl 1999-01-30-0000.09628911
jj diag.trajectory.drifters.plot.jl 1999-02-09-1800.09619815
jj diag.trajectory.drifters.plot.jl 1999-02-22-1200.09730553
jj diag.trajectory.drifters.plot.jl 1999-03-17-0000.09825463
jj diag.trajectory.drifters.plot.jl 1999-03-17-1200.09826019
jj diag.trajectory.drifters.plot.jl 1999-04-19-1200.09915131
jj diag.trajectory.drifters.plot.jl 1999-05-05-1800.09810599
jj diag.trajectory.drifters.plot.jl 1999-07-01-0600.09912352
jj diag.trajectory.drifters.plot.jl 1999-07-01-0600.09912367
jj diag.trajectory.drifters.plot.jl 1999-08-01-0600.00018323
jj diag.trajectory.drifters.plot.jl 1999-09-21-1200.09714829
jj diag.trajectory.drifters.plot.jl 1999-09-23-0000.09714825
jj diag.trajectory.drifters.plot.jl 1999-11-28-1800.00008523
jj diag.trajectory.drifters.plot.jl 2000-01-18-1800.00026399
jj diag.trajectory.drifters.plot.jl 2000-02-03-0000.00009379
jj diag.trajectory.drifters.plot.jl 2000-03-09-0600.09820361
jj diag.trajectory.drifters.plot.jl 2000-03-10-1800.00020364
jj diag.trajectory.drifters.plot.jl 2000-04-12-1200.00018679
jj diag.trajectory.drifters.plot.jl 2000-04-29-1800.00018677
jj diag.trajectory.drifters.plot.jl 2000-07-17-1200.00025018
jj diag.trajectory.drifters.plot.jl 2000-07-18-0600.00025020
jj diag.trajectory.drifters.plot.jl 2000-07-26-0600.00016983
jj diag.trajectory.drifters.plot.jl 2000-07-27-0600.00012221
jj diag.trajectory.drifters.plot.jl 2000-07-27-0600.00016866
jj diag.trajectory.drifters.plot.jl 2000-10-06-0000.00007402
jj diag.trajectory.drifters.plot.jl 2000-10-06-0600.00007574
jj diag.trajectory.drifters.plot.jl 2000-12-07-1800.00009394
jj diag.trajectory.drifters.plot.jl 2000-12-10-1800.00018742
jj diag.trajectory.drifters.plot.jl 2000-12-19-1800.00018755
jj diag.trajectory.drifters.plot.jl 2001-01-10-0000.00030717
jj diag.trajectory.drifters.plot.jl 2001-04-07-1800.00018691
jj diag.trajectory.drifters.plot.jl 2001-04-08-0600.00018692
jj diag.trajectory.drifters.plot.jl 2001-04-23-0600.00030052
jj diag.trajectory.drifters.plot.jl 2001-08-02-0600.00021606
jj diag.trajectory.drifters.plot.jl 2001-08-02-1800.00024057
jj diag.trajectory.drifters.plot.jl 2001-09-18-1800.00025567
jj diag.trajectory.drifters.plot.jl 2001-09-20-0600.00025722
jj diag.trajectory.drifters.plot.jl 2001-09-21-0600.00024828
jj diag.trajectory.drifters.plot.jl 2001-11-08-1200.00029495
jj diag.trajectory.drifters.plot.jl 2001-12-09-1200.00030724
jj diag.trajectory.drifters.plot.jl 2001-12-10-1800.00030723
jj diag.trajectory.drifters.plot.jl 2002-02-19-1800.00012371
jj diag.trajectory.drifters.plot.jl 2002-03-10-1800.00030439
jj diag.trajectory.drifters.plot.jl 2002-04-01-0600.00034160
jj diag.trajectory.drifters.plot.jl 2002-04-05-0600.00034161
jj diag.trajectory.drifters.plot.jl 2002-04-29-0000.00034354
jj diag.trajectory.drifters.plot.jl 2002-09-20-0000.02234182
jj diag.trajectory.drifters.plot.jl 2002-09-20-1800.00034183
jj diag.trajectory.drifters.plot.jl 2002-11-06-0600.00013572
jj diag.trajectory.drifters.plot.jl 2002-12-06-0600.00036972
jj diag.trajectory.drifters.plot.jl 2002-12-08-0000.00036973
jj diag.trajectory.drifters.plot.jl 2003-01-05-0600.02236978
jj diag.trajectory.drifters.plot.jl 2003-02-10-1200.00025516
jj diag.trajectory.drifters.plot.jl 2003-03-27-1200.00039662
jj diag.trajectory.drifters.plot.jl 2003-04-28-0600.00039663
jj diag.trajectory.drifters.plot.jl 2003-08-15-1200.00039671
jj diag.trajectory.drifters.plot.jl 2003-08-17-0000.02339668
jj diag.trajectory.drifters.plot.jl 2003-09-15-0600.00039670
jj diag.trajectory.drifters.plot.jl 2003-09-17-0000.00039669
jj diag.trajectory.drifters.plot.jl 2003-12-11-1200.00044297
jj diag.trajectory.drifters.plot.jl 2004-01-13-1800.02339285
jj diag.trajectory.drifters.plot.jl 2004-02-04-1800.02339294
jj diag.trajectory.drifters.plot.jl 2004-02-05-0600.00039293
jj diag.trajectory.drifters.plot.jl 2004-02-15-0600.02340007
jj diag.trajectory.drifters.plot.jl 2004-02-27-1200.02339255
jj diag.trajectory.drifters.plot.jl 2004-04-04-1800.02339263
jj diag.trajectory.drifters.plot.jl 2004-04-06-1200.02239182
jj diag.trajectory.drifters.plot.jl 2004-05-23-1200.00036929
jj diag.trajectory.drifters.plot.jl 2004-08-24-0600.02339258
jj diag.trajectory.drifters.plot.jl 2004-09-04-0000.00045973
jj diag.trajectory.drifters.plot.jl 2004-11-11-0000.00045969
jj diag.trajectory.drifters.plot.jl 2004-11-13-1800.00045962
jj diag.trajectory.drifters.plot.jl 2004-11-16-0000.00045964
jj diag.trajectory.drifters.plot.jl 2004-11-18-0000.00045966
jj diag.trajectory.drifters.plot.jl 2004-12-18-0000.00053390
jj diag.trajectory.drifters.plot.jl 2005-04-12-0000.00053341
jj diag.trajectory.drifters.plot.jl 2005-05-15-1800.00045959
jj diag.trajectory.drifters.plot.jl 2005-05-16-0000.00045958
jj diag.trajectory.drifters.plot.jl 2005-05-24-0600.00052828
jj diag.trajectory.drifters.plot.jl 2005-11-04-0000.00060298
jj diag.trajectory.drifters.plot.jl 2005-12-05-1800.00060309
jj diag.trajectory.drifters.plot.jl 2006-02-19-0000.00060287
jj diag.trajectory.drifters.plot.jl 2006-03-05-0600.00036165
jj diag.trajectory.drifters.plot.jl 2006-12-11-0600.00063875
jj diag.trajectory.drifters.plot.jl 2007-03-12-1200.00071094
jj diag.trajectory.drifters.plot.jl 2007-03-12-1800.00071114
jj diag.trajectory.drifters.plot.jl 2007-03-13-1800.00071124
jj diag.trajectory.drifters.plot.jl 2007-03-17-0000.00062077
jj diag.trajectory.drifters.plot.jl 2007-09-18-1800.00072196
jj diag.trajectory.drifters.plot.jl 2007-09-20-0600.00072199
jj diag.trajectory.drifters.plot.jl 2007-11-01-0000.00060328
jj diag.trajectory.drifters.plot.jl 2007-12-08-0000.00072251
jj diag.trajectory.drifters.plot.jl 2007-12-08-1800.00072320
jj diag.trajectory.drifters.plot.jl 2007-12-10-1800.00072247
jj diag.trajectory.drifters.plot.jl 2008-01-02-0600.00072259
jj diag.trajectory.drifters.plot.jl 2008-01-20-0600.00072325
jj diag.trajectory.drifters.plot.jl 2008-01-20-1800.00072217
jj diag.trajectory.drifters.plot.jl 2008-02-14-1800.00072332
jj diag.trajectory.drifters.plot.jl 2008-02-15-0600.00072246
jj diag.trajectory.drifters.plot.jl 2008-03-05-1800.00044647
jj diag.trajectory.drifters.plot.jl 2008-03-13-0000.00072238
jj diag.trajectory.drifters.plot.jl 2008-04-10-0600.00072338
jj diag.trajectory.drifters.plot.jl 2008-04-28-1800.00041534
jj diag.trajectory.drifters.plot.jl 2008-05-28-1200.00071201
jj diag.trajectory.drifters.plot.jl 2008-11-06-0000.00083344
jj diag.trajectory.drifters.plot.jl 2008-11-07-0000.00083341
jj diag.trajectory.drifters.plot.jl 2008-12-15-1200.00083358
jj diag.trajectory.drifters.plot.jl 2008-12-18-1200.00083356
jj diag.trajectory.drifters.plot.jl 2008-12-27-1800.00072182
jj diag.trajectory.drifters.plot.jl 2009-04-11-0000.00060291
jj diag.trajectory.drifters.plot.jl 2009-05-12-1800.00041577
jj diag.trajectory.drifters.plot.jl 2009-06-26-1800.00083506
jj diag.trajectory.drifters.plot.jl 2009-09-21-0000.00060362
jj diag.trajectory.drifters.plot.jl 2009-09-21-1800.00060348
jj diag.trajectory.drifters.plot.jl 2009-09-23-1800.00060344
jj diag.trajectory.drifters.plot.jl 2009-11-01-0600.00072213
jj diag.trajectory.drifters.plot.jl 2009-12-13-1200.00041536
jj diag.trajectory.drifters.plot.jl 2010-01-23-0600.00041547
jj diag.trajectory.drifters.plot.jl 2010-01-26-1200.00075285
jj diag.trajectory.drifters.plot.jl 2010-01-27-1200.00090504
jj diag.trajectory.drifters.plot.jl 2010-02-14-0600.00090490
jj diag.trajectory.drifters.plot.jl 2010-02-15-0600.00090501
jj diag.trajectory.drifters.plot.jl 2010-02-19-0000.00090492
jj diag.trajectory.drifters.plot.jl 2010-02-21-0600.00075297
jj diag.trajectory.drifters.plot.jl 2010-03-08-0600.00075281
jj diag.trajectory.drifters.plot.jl 2010-04-10-0000.00079093
jj diag.trajectory.drifters.plot.jl 2010-04-12-1200.00075288
jj diag.trajectory.drifters.plot.jl 2010-05-16-1800.00079089
jj diag.trajectory.drifters.plot.jl 2010-08-10-0000.00079094
jj diag.trajectory.drifters.plot.jl 2010-09-17-0000.00071500
jj diag.trajectory.drifters.plot.jl 2010-09-17-1200.00071499
jj diag.trajectory.drifters.plot.jl 2010-09-21-0600.00071486
jj diag.trajectory.drifters.plot.jl 2010-10-04-1800.00070972
jj diag.trajectory.drifters.plot.jl 2010-10-10-1800.00070969
jj diag.trajectory.drifters.plot.jl 2010-11-11-0000.00071029
jj diag.trajectory.drifters.plot.jl 2010-11-11-0600.00071032
jj diag.trajectory.drifters.plot.jl 2010-11-11-1800.00071030
jj diag.trajectory.drifters.plot.jl 2010-11-12-0000.00071031
jj diag.trajectory.drifters.plot.jl 2010-11-13-0000.00075300
jj diag.trajectory.drifters.plot.jl 2010-11-13-0600.00075303
jj diag.trajectory.drifters.plot.jl 2010-11-13-1800.00071033
jj diag.trajectory.drifters.plot.jl 2010-11-17-1200.00075293
jj diag.trajectory.drifters.plot.jl 2010-11-23-1200.00079090
jj diag.trajectory.drifters.plot.jl 2010-11-28-1800.00075302
jj diag.trajectory.drifters.plot.jl 2011-01-04-0600.00039273
jj diag.trajectory.drifters.plot.jl 2011-01-25-0000.00092967
jj diag.trajectory.drifters.plot.jl 2011-02-21-0600.00072249
jj diag.trajectory.drifters.plot.jl 2011-02-21-1200.00075279
jj diag.trajectory.drifters.plot.jl 2011-03-10-1200.00092989
jj diag.trajectory.drifters.plot.jl 2011-04-20-0000.00100962
jj diag.trajectory.drifters.plot.jl 2011-04-20-0000.00100966
jj diag.trajectory.drifters.plot.jl 2011-05-09-0000.00090163
jj diag.trajectory.drifters.plot.jl 2011-09-25-1200.00060346
jj diag.trajectory.drifters.plot.jl 2011-09-26-0000.00075294
jj diag.trajectory.drifters.plot.jl 2011-09-26-0600.00075290
jj diag.trajectory.drifters.plot.jl 2011-09-26-1200.00075292
jj diag.trajectory.drifters.plot.jl 2011-09-26-1800.00075287
jj diag.trajectory.drifters.plot.jl 2011-12-11-0000.00090156
jj diag.trajectory.drifters.plot.jl 2012-04-15-0000.00039801
