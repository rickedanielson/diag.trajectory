/*
 * This program is designed to create a netcdf template for (i.e. define
 * the fundamental dimensions, variables, and attributes, including data
 * for the coordinate variables) and transfer data on the track of some
 * event (e.g. a cyclone).  The input dates are assumed to be in increasing
 * order.  A time index is included to facilitate areal averaging over a
 * non-stationary domain, in GrADS - RD November 1999, January 2000.
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <stdlib.h>
#include "netcdf.h"
#include "udunits.h"

#define LEN        100
#define LOTS       4500
#define PARAMS     12                        /* number of parameters */
#define LATS       1                         /* number of latitudes used to define a track position */
#define LONS       1                         /* number of longitudes used to define a track position */

main(int argc, char *argv[])
{
    time_t time_val;                                                          /* time declarations (first!) */
    struct tm tma;
    char *month[] = {"January","February","March","April","May","June","July",
                     "August","September","October","November","December"};

    FILE *fpa, *fopen();                                                      /* basic declarations */
    int a, b, c, d, count;
    char infila[LEN], outfila[LEN];
    char date[LOTS][LEN], tag[LOTS][LEN], notes[LOTS*LEN], missnotes[LOTS*LEN], line[5*LOTS];
    float lats[LOTS], lons[LOTS], latd[LOTS], lond[LOTS], hour[LOTS], dist[LOTS], bear[LOTS],
          miss[LOTS], tindexvals[LOTS];

    int status, ncid;                                                         /* netCDF declarations */
    int time_dim, lat_dim, lon_dim, len_dim;
    int time_id,  lat_id,  lon_id,  len_id,  par_id[PARAMS];
    int dims[4], int_range[2];
    short short_val;
    float float_val, float_range[2];
    double double_range[2];
    char attrib[LEN], lina[LEN], linb[LEN];
    static size_t begin[] = {0, 0, 0, 0}, end[] = {0, 0, 0, 0};

    char udunita[LEN], udunitb[LEN];                                          /* udunits declarations */
    utUnit utunita, utunitb;
    double slope, intercept, times[LOTS];
    char predate[] = "hours since ";
    char postdate[] = ":00:0.0";

    char *long_name[PARAMS] = {"Latitude of simulated drifter",               /* define parameters */
                               "Longitude of simulated drifter",
                               "Time index of simulated drifter",
                               "Time index for GrADs",
                               "Latitude of actual drifter",
                               "Longitude of actual drifter",
                               "Distance from actual to simulated drifter",
                               "Bearing from actual to simulated drifter",
                               "Collocated variable a (e.g., bathymetric depth)",
                               "Collocated variable b (e.g., distance to coast)",
                               "Collocated variable c",
                               "Collocated variable d"};
    char *parname[PARAMS] = {     "simlat",        "simlon",     "simtim",  "tindex",
                                  "drilat",        "drilon",       "dist",    "bear",
                                   "colla",         "collb",      "collc",   "colld"};
    char *units[PARAMS] = {"degrees_north",  "degrees_east",          "h",        "",
                           "degrees_north",  "degrees_east",          "m", "degrees",
                                        "",              "",           "",        ""};
    float missing_value[PARAMS] = { -9999.,          -9999.,       -9999.,    -9999.,
                                    -9999.,          -9999.,       -9999.,    -9999.,
                                    -9999.,          -9999.,       -9999.,    -9999.};
    nc_type partype[PARAMS] = {   NC_FLOAT,        NC_FLOAT,     NC_FLOAT,  NC_FLOAT,
                                  NC_FLOAT,        NC_FLOAT,     NC_FLOAT,  NC_FLOAT,
                                  NC_FLOAT,        NC_FLOAT,     NC_FLOAT,  NC_FLOAT};
    int pardims[PARAMS] = {              3,               3,            3,         3,
                                         3,               3,            3,         3,
                                         3,               3,            3,         3};

    float lat = 45.0, lon = 180.0;                                            /* define (bogus) dimension values */
    int lens[LOTS];
    for (a = 0; a < LEN; a++)
      lens[a] = a;

    if (argc != 2) {
      printf("Usage: %s 2012-10-20-0000.00000000.traj.2012-05-22-0000\n",argv[0]);
      exit(1);
    }

    if ((fpa = fopen(argv[1],"r")) == NULL) {
      fprintf(stderr, "ERROR : couldn't open %s\n", argv[1]);
      exit(1);
    }
    printf("reading %s\n",argv[1]);

    if (utInit("") != 0) {
      printf("couldn't initialize Unidata units library\n");
      exit(2);
    }

    a = 0;
    while (fgets(line,LOTS,fpa) != NULL) {                                    /* get lines from datelist until none left */
      sscanf(line, "%s %f %f %f %f %f %f %f",
        date[a], &latd[a], &lond[a], &hour[a], &lats[a], &lons[a], &dist[a], &bear[a]);
      if (lons[a] > 180) lons[a] -= 360.0;
      if (lond[a] > 180) lond[a] -= 360.0;
      date[a][13] = '\0';
      miss[a] = -9999.;
      strcpy(udunita,predate);
      strcat(udunita,date[a]);
      strcat(udunita,postdate);
      udunita[22] = ' ';
      strcpy(udunitb,"hours since 1-1-1 00:00:0.0");                          /* form the udunits date from the original */
      if (utScan(udunita,&utunita) != 0 || utScan(udunitb,&utunitb) != 0)     /* obtain the number of hours */
        printf("either %s or %s is incomprehensible\n",udunita,udunitb);      /* from the reference date */
      if (utConvert(&utunita,&utunitb,&slope,&intercept) != 0)                /* to the desired date */
        printf("couldn't convert %s to %s\n",udunita,udunitb);
      times[a] = intercept;                                                   /* store and count all referenced dates */
      a++;
    }
    fclose(fpa);
    utTerm();

    if (a == 0) {
      fprintf(stderr, "ERROR : no tracks were found...\n");
      exit(1);
    }
    count = a;

    strcpy(outfila,argv[1]);                                                  /* form the template file name */
    strcat(outfila,".nc");
    status = nc_create(outfila, NC_CLOBBER, &ncid);                           /* create file and enter define mode */
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    printf("writing %s\n",outfila);

    status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);               /* define dimensions */
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_def_dim(ncid, "lat", LATS, &lat_dim);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_def_dim(ncid, "lon", LONS, &lon_dim);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
/*  status = nc_def_dim(ncid, "len", LEN, &len_dim);
    printf("with dimensions of time %d, lat %d, lon %d, and len %d\n",count,LATS,LONS,LEN);  */
    printf("with dimensions of time %d, lat %d, and lon %d\n",count,LATS,LONS);

    dims[0] = time_dim;                                                       /* define dimension variables */
    status = nc_def_var(ncid, "time", NC_DOUBLE, 1, dims, &time_id);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    dims[0] = lat_dim;
    status = nc_def_var(ncid, "lat", NC_FLOAT, 1, dims, &lat_id);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    dims[0] = lon_dim;
    status = nc_def_var(ncid, "lon", NC_FLOAT, 1, dims, &lon_id);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
/*  dims[0] = len_dim;
    status = nc_def_var(ncid, "len", NC_INT, 1, dims, &len_id);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}  */

    printf("including variables");
    dims[0] = time_dim;                                                       /* define variables */
    dims[1] = lat_dim;
    dims[2] = lon_dim;
    dims[3] = len_dim;
    for (a = 0; a < PARAMS; a++) {
      status = nc_def_var(ncid, parname[a], partype[a], pardims[a], dims, &par_id[a]);
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

    strcpy(attrib,"long_name");                                               /* assign lat attributes */
    strcpy(line,"Bogus Latitude");
    status = nc_put_att_text(ncid, lat_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"units");
    strcpy(line,"degrees_north");
    status = nc_put_att_text(ncid, lat_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    strcpy(attrib,"long_name");                                               /* assign lon attributes */
    strcpy(line,"Bogus Longitude");
    status = nc_put_att_text(ncid, lon_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"units");
    strcpy(line,"degrees_east");
    status = nc_put_att_text(ncid, lon_id, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    for (a = 0; a < PARAMS; a++) {                                            /* assign parameter attributes */
      strcpy(attrib,"long_name");
      strcpy(line,long_name[a]);
      status = nc_put_att_text(ncid, par_id[a], attrib, strlen(line)+1, line);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      if (partype[a] == NC_FLOAT) {
        strcpy(attrib,"units");
        strcpy(line,units[a]);
        status = nc_put_att_text(ncid, par_id[a], attrib, strlen(line)+1, line);
        if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
        strcpy(attrib,"missing_value");
        float_val = missing_value[a];
        status = nc_put_att_float(ncid, par_id[a], attrib, NC_FLOAT, 1, &float_val);
        if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      }
    }

    strcpy(attrib,"Conventions");                                             /* assign global attributes */
    strcpy(line,"COARDS");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"title");
    sprintf(line,"Trajectory data for the %s datelist",argv[1]);
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"history");
    time_val = time(NULL);
    tma = *localtime(&time_val);
    strcpy(line,"Simulated trajectories have been launched from the position\n");
    strcat(line,"of each six-hourly drifter location and advected forward until either the\n");
    sprintf(lina,"24th day or the end of the actual drifter is reached - RD %s %d.",month[tma.tm_mon],1900+tma.tm_year);
    strcat(line,lina);
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    strcpy(attrib,"description");
    strcpy(line,"Interpolation is based on the three-hourly GlobCurrent\n");
    strcat(line,"analyses of combined current.  Spatially the nearest neighbour is employed\n");
    strcat(line,"while temporal interpolation is linear.  The position update is taken from\n");
    strcat(line,"www.movable-type.co.uk/scripts/latlong.html (destination given distance and\n");
    strcat(line,"bearing from start).  It is expected that trajectories also depend slightly\n");
    strcat(line,"on wind and short waves.  Depth and distance to the coast may be employed\n");
    strcat(line,"as proxies for simulation quality");
    status = nc_put_att_text(ncid, NC_GLOBAL, attrib, strlen(line)+1, line);

    status = nc_enddef(ncid);                                                 /* leave define mode */
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    end[0] = (long)count;                                                     /* store dimension variables */
    status = nc_put_vara_double(ncid, time_id, begin, end, times);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    end[0] = LATS;
    status = nc_put_vara_float(ncid, lat_id, begin, end, &lat);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    end[0] = LONS;
    status = nc_put_vara_float(ncid, lon_id, begin, end, &lon);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
/*  end[0] = LEN;
    status = nc_put_vara_int(ncid, len_id, begin, end, lens);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}  */

    end[0] = (long)count;                                                     /* store variables */
    end[1] = LATS;
    end[2] = LONS;
    end[3] = LEN;
    status = nc_put_vara_float(ncid, par_id[0], begin, end, lats);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[1], begin, end, lons);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[2], begin, end, hour);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    for (a = 1; a <= LOTS; a++)
      tindexvals[a-1] = a;
    status = nc_put_vara_float(ncid, par_id[3], begin, end, tindexvals);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    status = nc_put_vara_float(ncid, par_id[4], begin, end, latd);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[5], begin, end, lond);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[6], begin, end, dist);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[7], begin, end, bear);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[8], begin, end, miss);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[9], begin, end, miss);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[10], begin, end, miss);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_put_vara_float(ncid, par_id[11], begin, end, miss);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    status = nc_close(ncid);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
}
