#!/usr/bin/perl -w

# Convert 6-hourly netcdf files to monthly netcdf files and
# perform a regional subsetting - RD February 2012, September 2015
# followed by tar cvfz norskhavet_2011.taz norskhavet_2011*nc

$outste = "norskhavet_";
$srcdir = "/home/cercache/users/rdaniels/work/workt/links";
$latrng =  "65.,75.";
$lonrng = "-10.,20.";

printf("reading the list of files in %s\n",$srcdir);                          # read the list of all source files
opendir(WRK, $srcdir) or die("ERROR reading $srcdir : $!\n");                 # (e.g., 20110703090000-GLOBCURRENT-L4
@srcfile = grep /^201.*fv01.0.nc$/, readdir(WRK) ; closedir(WRK);             #  -CUReul_15m-ALT_SUM-v01.0-fv01.0.nc)
@srcfile = sort(@srcfile);

foreach $tmpfila (@srcfile) {                                                 # for each global file of interest
  $tmpfilb = substr($tmpfila,0,6);                                            # perform the regional subsetting and
  if ($tmpfilb > 201107 && $tmpfilb < 201213) {                               # add the result to the monthly file
    $srcfil = "$srcdir/$tmpfila";
    $outfil = $outste.$tmpfila;
    $command = "ncks -d lat,$latrng -d lon,$lonrng $srcfil -o $outfil\n";
    print $command ; system($command);
#   $command = "ncrcat --rec_apn tmp$outfil $outfil\n";
#   print $command ; system($command);
#   $command = "rm tmp$outfil\n";
#   print $command ; system($command);
  }
}
