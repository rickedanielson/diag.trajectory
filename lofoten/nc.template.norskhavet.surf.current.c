/*
 * This program is designed to create a netcdf template (i.e. define the
 * fundamental dimensions, variables, and attributes, including data for
 * the coordinate variables) for surface-level spherical coordinate data.
 * The input dates are assumed to be in increasing order - RD November 1999.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "netcdf.h"
#include "udunits.h"

#define LEN        100
#define LOTS       10000
#define PARAMS     5                         /* number of parameters */
#define LEVELS     1                         /* number of pressure levels */
#define LATS       100                       /* number of domain latitudes */
#define LONS       300                       /* number of domain longitudes */
#define FIRSTLAT   65.05                     /* first latitude (deg), last is FIRST - DEL * (LATS-1) */
#define DELLAT     0.1                       /* difference in adjacent latitudes (deg) */
#define FIRSTLON  -9.95                      /* first longitude (deg), last is FIRST + DEL * (LONS-1) */
#define DELLON     0.1                       /* difference in adjacent longitudes (deg) */

main(int argc, char *argv[])
{
    time_t time_val;                                                          /* time declarations (first!) */
    struct tm tma;
    char *month[12] = {"January","February","March","April","May","June","July",
                       "August","September","October","November","December"};

    FILE *fpa, *fopen();                                                      /* basic declarations */
    char infila[LEN], outfila[LEN];
    char date[LEN], tag[LEN], line[5*LOTS];
    int a, b, c, d, count;

    int status, ncid;                                                         /* netCDF declarations */
    int time_dim, level_dim, lat_dim, lon_dim;
    int time_id,  level_id,  lat_id,  lon_id,  par_id[PARAMS], dims[4];
    short short_val;
    float float_val, float_range[2];
    double double_range[2];
    char attrib[LEN], lina[LEN], linb[LEN];
    static size_t begin[3] = {0,0,0}, level_end = LEVELS, lat_end = LATS, lon_end = LONS, len_end = LEN;
    size_t time_end, indic_end[3];

    char udunita[LEN], udunitb[LEN];                                          /* udunits declarations */
    utUnit utunita, utunitb;
    double slope, intercept, times[LOTS];
    char predate[] = "hours since ";
    char postdate[] = ":00:0.0";

/*  char *long_name[PARAMS] = {"eastward_ekman_current_velocity",
                               "northward_ekman_current_velocity",  */
    char *long_name[PARAMS] = {"eastward_eulerian_current_velocity",
                               "northward_eulerian_current_velocity",
                               "eulerian_current_velocity_magnitude",
                               "eulerian_current_velocity_uncertainty",
                               "eulerian_current_velocity_flag"};
    char *parname[PARAMS] = {          "u",     "v",     "w", "error", "flag"};
    char *units[PARAMS] = {          "m/s",   "m/s",   "m/s",   "m/s",    "-"};
    float valid_range_low[PARAMS] = { -32.,    -32.,    -32.,    -32.,   -32.};
    float valid_range_high[PARAMS] = { 32.,     32.,     32.,     32.,    32.};
    float add_offset[PARAMS] = {        0.,      0.,      0.,      0.,     0.};
    float scale_factor[PARAMS] = {   0.001,   0.001,   0.001,   0.001,  0.001};
    short missing_value[PARAMS] = {  32766,   32766,   32766,   32766,  32766};
    short precision[PARAMS] = {          3,       3,       3,       3,      3};
    int numdims[PARAMS] = {              4,       4,       4,       4,      4};

    static float level[LEVELS] = {15};                                        /* define dimension values */
    float lat[LOTS], lon[LOTS];

    for (a = 0; a < LATS; a++)
      lat[a] = FIRSTLAT + a * DELLAT;
    for (a = 0; a < LONS; a++)
      lon[a] = FIRSTLON + a * DELLON;

    if (argc != 2) {
      printf("Usage: %s datelist\n",argv[0]);
      exit(1);
    }

    strcpy(infila,argv[1]);                                                   /* open datelist */
    if ((fpa = fopen(infila,"r")) == NULL) {
      printf("couldn't open %s\n",argv[1]);
      exit(1);
    }
    printf("using %s",infila);

    if (utInit("") != 0) {
      printf("couldn't initialize Unidata units library\n");
      exit(2);
    }

    a = 0;
    while (fgets(line,LOTS,fpa) != NULL) {                                    /* get dates */
      sscanf(line, "%s %s", tag, date);                                       /* form the desired date */
      strcpy(udunita,predate);
      strcat(udunita,date);
      strcat(udunita,postdate);
      udunita[22] = ' ';
      strcpy(udunitb,"hours since 1-1-1 00:00:0.0");                          /* and the reference date */
      if (utScan(udunita,&utunita) != 0 || utScan(udunitb,&utunitb) != 0)     /* to obtain the number of hours */
        printf("either %s or %s is incomprehensible\n",udunita,udunitb);      /* from the reference date */
      if (utConvert(&utunita,&utunitb,&slope,&intercept) != 0)                /* to the desired date */
        printf("couldn't convert %s to %s\n",udunita,udunitb);
      times[a] = intercept;                                                   /* store and count all referenced dates */
      a++;
    }
    fclose(fpa);
    utTerm();

    if (a == 0) {
      printf(" but no tracks were found...\n");
      exit(1);
    }
    count = a;

    strcpy(outfila,argv[1]);                                                  /* form the template file name */
    strcat(outfila,".surf.current.nc");
    status = nc_create(outfila, NC_CLOBBER, &ncid);                           /* create file and enter define mode */
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    printf(" to create %s\n",outfila);

    status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);               /* define dimensions */
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_def_dim(ncid, "level", LEVELS, &level_dim);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_def_dim(ncid, "lat", LATS, &lat_dim);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_def_dim(ncid, "lon", LONS, &lon_dim);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    printf("with dimensions of time %d, level %d, lat %d, and lon %d\n",count,LEVELS,LATS,LONS);

    dims[0] = time_dim;                                                       /* define dimension variables */
    status = nc_def_var(ncid, "time", NC_DOUBLE, 1, dims, &time_id);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    dims[0] = level_dim;
    status = nc_def_var(ncid, "level", NC_FLOAT, 1, dims, &level_id);
    dims[0] = lat_dim;
    status = nc_def_var(ncid, "lat", NC_FLOAT, 1, dims, &lat_id);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    dims[0] = lon_dim;
    status = nc_def_var(ncid, "lon", NC_FLOAT, 1, dims, &lon_id);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    printf("including variables");
    dims[0] = time_dim;                                                       /* define variables */
    for (a = 0; a < PARAMS; a++) {
      if (numdims[a] == 3) {
        dims[1] = lat_dim;
        dims[2] = lon_dim;
      }
      else {
        dims[1] = level_dim;
        dims[2] = lat_dim;
        dims[3] = lon_dim;
      }
      status = nc_def_var(ncid, parname[a], NC_SHORT, numdims[a], dims, &par_id[a]);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      printf(" %s",parname[a]);
    }
    printf("\n");

    strcpy(attrib,"long_name");                                               /* assign time attributes */
    strcpy(line,"Time");
    status = nc_put_att_text(ncid, time_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"units");
    strcpy(line,"hours since 1-1-1 00:00:0.0");
    status = nc_put_att_text(ncid, time_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"actual_range");
    double_range[0] = times[0];
    double_range[1] = times[count-1];
    status = nc_put_att_double(ncid, time_id, attrib, NC_DOUBLE, 2, double_range);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"delta_t");
    strcpy(line,"0000-00-00 06:00:00");
    status = nc_put_att_text(ncid, time_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    strcpy(attrib,"long_name");                                               /* assign level attributes */
    strcpy(line,"Depth below ocean surface");
    status = nc_put_att_text(ncid, level_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"units");
    strcpy(line,"m");
    status = nc_put_att_text(ncid, level_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"actual_range");
    float_range[0] = level[LEVELS-1];
    float_range[1] = level[0];
    status = nc_put_att_float(ncid, level_id, attrib, NC_FLOAT, 2, float_range);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"positive");
    strcpy(line,"down");
    status = nc_put_att_text(ncid, level_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"GRIB_id");
    short_val = 100;
    status = nc_put_att_short(ncid, level_id, attrib, NC_SHORT, 1, &short_val);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"GRIB_name");
    strcpy(line,"hPa");
    status = nc_put_att_text(ncid, level_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    strcpy(attrib,"long_name");                                               /* assign lat attributes */
    strcpy(line,"Latitude");
    status = nc_put_att_text(ncid, lat_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"units");
    strcpy(line,"degrees_north");
    status = nc_put_att_text(ncid, lat_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"actual_range");
    float_range[0] = lat[0];
    float_range[1] = lat[LATS-1];
    status = nc_put_att_float(ncid, lat_id, attrib, NC_FLOAT, 2, float_range);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    strcpy(attrib,"long_name");                                               /* assign lon attributes */
    strcpy(line,"Longitude");
    status = nc_put_att_text(ncid, lon_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"units");
    strcpy(line,"degrees_east");
    status = nc_put_att_text(ncid, lon_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"actual_range");
    float_range[0] = lon[0];
    float_range[1] = lon[LONS-1];
    status = nc_put_att_float(ncid, lon_id, attrib, NC_FLOAT, 2, float_range);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    for (a = 0; a < PARAMS; a++) {                                            /* assign parameter attributes */
      strcpy(attrib,"long_name");
      strcpy(line,long_name[a]);
      status = nc_put_att_text(ncid, par_id[a], attrib, strlen(line)+1, line);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      strcpy(attrib,"units");
      strcpy(line,units[a]);
      status = nc_put_att_text(ncid, par_id[a], attrib, strlen(line)+1, line);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      strcpy(attrib,"valid_range");
      float_range[0] = valid_range_low[a];
      float_range[1] = valid_range_high[a];
      status = nc_put_att_float(ncid, par_id[a], attrib, NC_FLOAT, 2, float_range);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      strcpy(attrib,"add_offset");
      float_val = add_offset[a];
      status = nc_put_att_float(ncid, par_id[a], attrib, NC_FLOAT, 1, &float_val);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      strcpy(attrib,"scale_factor");
      float_val = scale_factor[a];
      status = nc_put_att_float(ncid, par_id[a], attrib, NC_FLOAT, 1, &float_val);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      strcpy(attrib,"missing_value");
      short_val = missing_value[a];
      status = nc_put_att_short(ncid, par_id[a], attrib, NC_SHORT, 1, &short_val);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      strcpy(attrib,"precision");
      short_val = precision[a];
      status = nc_put_att_short(ncid, par_id[a], attrib, NC_SHORT, 1, &short_val);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    }

/*
strcpy(attrib,"title") ; strcpy(line,"Stokes component of near surface current velocity");  ** assign global attributes **
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"abstract") ; strcpy(line,"This component is updated every 3-6 hours but new information is\navailable only locally, where wind and predicted or observed\nwave data have been acquired during this interval.");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"topiccategory") ; strcpy(line,"Oceans ClimatologyMeteorologyAtmosphere");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
** strcpy(attrib,"keywords") ; strcpy(line,""EARTH SCIENCE","OCEANS","OCEAN CIRCULATION","OCEAN CURRENTS","",""");
"CONSORTIA/INSTITUTIONS","","","","NERSC","Nansen Environmental and Remote Sensing Centre","http://www.nersc.no/main/index2.php"
"GOVERNMENT AGENCIES-NON-US","FRANCE","","","FR/IFREMER/CERSAT","Centre ERS d'Archivage et de Traitement, French Research Institute for Exploitation of the Sea, France","http://www.ifremer.fr/cersat/index.html"
"G - I","GHRSST-PP","GODAE High Resolution Sea Surface Temperature Pilot Project"
"Solar/Space Observing Instruments","","","","",""
strcpy(attrib,"gcmd_keywords") ; strcpy(line,"");  **
strcpy(attrib,"activity_type") ; strcpy(line,"L3 analysis product");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"Conventions") ; strcpy(line,"CF-1.0");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"product_name") ; strcpy(line,"GlobCurrent Version 2");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"history") ; strcpy(line,"2007-05-12 creation\n2007-05-14 update");
    time_val = time(NULL);
    tma = *localtime(&time_val);
    strcpy(line,"2007-05-12 creation\n");
    sprintf(lina,"%d-%d-%d %s %d.",1900+tma.tm_year,tma.tm_mon,tma.tm_mday,month[tma.tm_mon],1900+tma.tm_year);
    strcat(line,lina);
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"area") ; strcpy(line,"Global");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"southernmost_latitude") ; strcpy(line,"90S");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"northernmost_latitude") ; strcpy(line,"90N");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"westernmost_longitude") ; strcpy(line,"0E");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"easternmost_longitude") ; strcpy(line,"359.75E");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"start_date") ; strcpy(line,"2012-01-01 00:00:00 UTC");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"stop_date") ; strcpy(line,"2012-12-31 18:00:00 UTC");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"institution") ; strcpy(line,"Nansen Environmental and Remote Sensing Center");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"PI_name") ; strcpy(line,"Johnny Johannessen");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"contact") ; strcpy(line,"Johnny.Johannessen@nersc.no");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"distribution_statement") ; strcpy(line,"Free");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
strcpy(attrib,"project_name") ; strcpy(line,"GlobCurrent");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}  */

    strcpy(attrib,"title");
    sprintf(line,"GlobCurrent data for the %s datelist",argv[1]);
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"history");
    time_val = time(NULL);
    tma = *localtime(&time_val);
    strcpy(line,"This surface ocean current data was cropped from an\n");
    sprintf(lina,"archive at Ifremer - RD %s %d.",month[tma.tm_mon],1900+tma.tm_year);
    strcat(line,lina);
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"description");
    strcpy(line,"Data is either a combination such as geostrophic\n");
    strcat(line,"and Ekman current components, estimated separately\n");
    strcat(line,"and combined linearly, or is just one component.");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    status = nc_enddef(ncid);                                                 /* leave define mode */
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    time_end = (long)count;                                                   /* store dimension variables */
    status = nc_put_vara_double(ncid, time_id, begin, &time_end, times);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, level_id, begin, &level_end, level);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, lat_id, begin, &lat_end, lat);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, lon_id, begin, &lon_end, lon);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    status = nc_close(ncid);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
}

/*
topiccategory	A blank separated list of topic keywords describing the dataset. See below for applicable keywords.
keywords	A blank separated list of keywords describing the dataset. See below for list of applicable keywords.
gcmd_keywords	Newline separated list of GCMD scientific keywords describing the various variables. This will be used to categorize the datasets according to the "Topics and variables" menu selection in the metadata search facility. If proper standard names have been used for the variables, data will be mapped for search under Topics and variables even without the gcmd_keywords attribute. For more information see below.
activity_type	Comma separated list of activity types. See list below for applicable descriptions.
Conventions	The metadata convention used, should be "CF-1.0"
product_name	A product name of the dataset.
history	Modification history of the dataset. Should be of the form:
2007-05-12 creation
2007-06-10 revision and separated by newlines.
area	Area name describing the geographical area being studied. If several area names are used, separate them using comma. See below
southernmost_latitude	Elements to describe a geographical bounding box for the data. Should be a floating point value (decimal degrees).
northernmost_latitude	Elements to describe a geographical bounding box for the data. Should be a floating point value (decimal degrees).
westernmost_longitude	Elements to describe a geographical bounding box for the data. Should be a floating point value (decimal degrees).
easternmost_longitude	Elements to describe a geographical bounding box for the data. Should be a floating point value (decimal degrees).
start_date	Start date and time of the dataset in the form "2007-06-12 12:30:00 UTC"
stop_date	Stop date and time of the dataset in the form "2007-06-12 12:30:00 UTC"
institution	Name of the institution responsible for the dataset. Please use one of the standardised names below (short or long name).
PI_name	Name of the person responsible for the data set.
contact	email address to responsible user support or principal investigator. If the email address of the principal investigator is used, the variable "PI_name" should be set accordingly.
distribution_statement	A distribution statement, see below for a applicable list.
project_name 	Name of the project within which the data were collected
*/
