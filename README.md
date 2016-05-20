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
       mkdir plot.available plot.histogr plot.locate
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

# perform an initial analysis evaulation, but without calibration (and including all components, versus the 2.0_valid_remainder obs)
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

# verify that each subdir contains the expected number of files (e.g., 4357 + 3640 = 7997 files with 3408 dates)
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

# create the forward and backward extrapolated timeseries (separately for the geostrophic and Ekman components and combined for the total)
wrkt ; cd v2.0_global_025_deg_geostrophic ; ls z.list?? ; cd ..
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.ekman.jl ::: v2.0_global_025_deg_ekman_15m v2.0_global_025_deg_ekman_hs ::: z.listaa z.listab z.listac z.listad z.listae z.listaf z.listag z.listah | grep drift | sort > commands
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.geostrophic.jl v2.0_global_025_deg_geostrophic ::: z.listaa z.listab z.listac z.listad z.listae z.listaf z.listag z.listah | grep drift | sort >> commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       rm commands
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.total.jl v2.0_global_025_deg_geostrophic v2.0_global_025_deg_ekman_15m v2.0_global_025_deg_total_15m ::: z.listaa z.listab z.listac z.listad z.listae z.listaf z.listag z.listah | grep drift | sort > commands
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.total.jl v2.0_global_025_deg_geostrophic v2.0_global_025_deg_ekman_hs v2.0_global_025_deg_total_hs ::: z.listaa z.listab z.listac z.listad z.listae z.listaf z.listag z.listah | grep drift | sort >> commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       rm commands

# plot extrapolation histograms (forward and backward versus the actual values for assessment of bias in the extrapolation method)
wrkt ; nohup julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histogram.jl z.list > xcom &
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_ekman_15m
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_ekman_hs
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_geostrophic
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_total_15m
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.extrapolated.histoplot.jl v2.0_global_025_deg_total_hs
       gzip extrapolated.histogr*dat ; mv extrapolated.* all/plot.histogr

# identify the subset of drifter cal/val locations for which analyses are also available for much of 2001-2007 (call these the collocations)
wrkt ; mkdir fft
       split -l 400 all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib
       split -l 400 all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.nfft.jl ::: buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_?ali??? ::: ucur vcur | grep drift | sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib??.ucur.got2000 | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib??.ucur.not2000 | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.not2000
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid??.ucur.got2000 | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.got2000
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid??.ucur.not2000 | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.not2000
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib??.vcur.got2000 | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.vcur.got2000
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib??.vcur.not2000 | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.vcur.not2000
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid??.vcur.got2000 | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.vcur.got2000
       cat buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid??.vcur.not2000 | sort > all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.vcur.not2000
#      sort all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib > aa
#      cat  all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000 all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.not2000 | sort > bb ; diff aa bb ; rm aa bb
#      sort all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid > aa
#      cat  all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.got2000 all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.not2000 | sort > bb ; diff aa bb ; rm aa bb
       mv commands buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_?ali?* all/limbo

# plot the location and distribution of the collocations
wrkt ; cd all
       echo julia /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.evaluation.distribution.plot.jl buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.????.got2000 > commands
       echo grads -blc \"diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000\" >> commands
       echo grads -blc \"diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.vcur.got2000\" >> commands
       echo grads -blc \"diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.got2000\" >> commands
       echo grads -blc \"diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.vcur.got2000\" >> commands
       echo grads -blc \"diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.not2000\" >> commands
       echo grads -blc \"diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.vcur.not2000\" >> commands
       echo grads -blc \"diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.not2000\" >> commands
       echo grads -blc \"diag.trajectory.drifters.locate buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.vcur.not2000\" >> commands ; echo "
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e xvfb-run --mem=2000mb
       rm commands ; mv plot.trajectory.drifters.dots*locate*png plot.locate

# plot the all-collocation-averaged (weighted by daily obs number) one-sided NFFT spectra for each variable and identify submonthly variance
wrkt ; parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.heat.flux.timeseries.nfft.avg.jl  ::: all/all.flux.daily.locate_2.0_?ali?.????.got2000      | grep flux | sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e              julia  --mem=2000mb
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.heat.flux.timeseries.nfft.plot.jl all/all.flux.daily.locate_2.0_calib.????.got2000.spec
       xvfb-run -a julia /home1/homedir1/perso/rdaniels/bin/diag.heat.flux.timeseries.nfft.plot.jl all/all.flux.daily.locate_2.0_valid.????.got2000.spec
       cd all ; cat *calib.shfx.got2000.spec.stat *calib.lhfx.got2000.spec.stat *calib.wspd.got2000.spec.stat *calib.airt.got2000.spec.stat *calib.sstt.got2000.spec.stat *calib.shum.got2000.spec.stat
                cat *valid.shfx.got2000.spec.stat *valid.lhfx.got2000.spec.stat *valid.wspd.got2000.spec.stat *valid.airt.got2000.spec.stat *valid.sstt.got2000.spec.stat *valid.shum.got2000.spec.stat
       cd .. ; rm commands ; vi coads.gts.ncepnrt.heat.flux.colloc.discrete.triple.jl

# assemble the insitu and analysis data for a triple collocation cal/val
wrkt ; parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.colloc.discrete.insitu.jl all/buoydata_1993_2014_drogON.asc.nonmdt ::: all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_?ali?.????.got2000     | grep drift | sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia --mem=2000mb
       parallel --dry-run /home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.colloc.discrete.source.jl                                          ::: all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_?ali?.????.got2000_obs | grep drift | sort > commands
       cat commands | /home5/begmeil/tools/gogolist/bin/gogolist.py -e julia
       rm commands

# perform a global triple collocation cal/val (GLOBAL = true; RECALIB = false then true) and compare calibrated and uncalibrated RMSE
wrkt ; mkdir all/zali.recalib.false all/zali.recalib.true
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.vcur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.vcur.got2000_obs.comb
       mv all/*cali.glob all/zali.recalib.false
       cat               all/zali.recalib.false/*cali.glob | grep const
       vi  diag.trajectory.drifters.colloc.discrete.triple.jl                 (and inject the various alpha and beta)
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.vcur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.vcur.got2000_obs.comb
       mv all/*cali.glob all/zali.recalib.true
       cat               all/zali.recalib.false/*cali.glob | grep total
       cat               all/zali.recalib.true/*cali.glob  | grep total

# perform a local triple collocation cal/val (GLOBAL = false; RECALIB = false then true) and compare calibrated and uncalibrated RMSE
wrkt ; mkdir all/zali.recalib.local
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.vcur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.vcur.got2000_obs.comb
       mv all/*cali.locl all/zali.recalib.false
       cat               all/zali.recalib.false/*cali.locl | grep const
       vi  diag.trajectory.drifters.colloc.discrete.triple.jl                 (and inject the various alpha and beta)
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.ucur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib.vcur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.ucur.got2000_obs.comb
       jjj diag.trajectory.drifters.colloc.discrete.triple.jl all/buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_valid.vcur.got2000_obs.comb
       mv all/*cali.glob all/zali.recalib.local
       cat               all/zali.recalib.false/*cali.glob | grep total
       cat               all/zali.recalib.true/*cali.glob  | grep total
       cat               all/zali.recalib.local/*cali.glob | grep total

# plot local calibration and performance as a function of the variable of interest
wrkt ; cd all
       jjj diag.trajectory.drifters.evaluation.calval.plot.jl all/zali.recalib.false/*cali.locl
       mv 
