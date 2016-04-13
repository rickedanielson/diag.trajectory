/*
 * This program is designed to transfer TIME-LEVEL-LAT-LON time-mean
 * grids from source netcdf files (e.g. from the CDC) to TIME-LEVEL-LAT-LON
 * grids of a destination netcdf file using a list of dates as input.
 * We assume that space has already been allocated for the new parameter
 * in the destination file.  A sub-domain can be defined - RD October 1999
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "netcdf.h"
#include "/home/ricani/prog/include.netcdf/include/cdcnetcdf.h"

#define LEN        100
#define LOTS       3000
#define PARAMS     2
#define LEVSA      1                         /* number of input domain pressure levels */
#define LATSA      100                       /* number of input domain latitudes */
#define LONSA      300                       /* number of input domain longitudes */
#define LEVSB      1                         /* number of output domain pressure levels */
#define LATSB      100                       /* number of output domain latitudes */
#define LONSB      300                       /* number of output domain longitudes */
#define MISSING    -8e9    

void transfer(float in[LEVSA*LATSA*LONSA], float out[LEVSB*LATSB*LONSB]);
void intread(char file[], char param[], char date[], float *array, int arraylen);

float in[LEVSA*LATSA*LONSA], outa[LEVSB*LATSB*LONSB], outb[LEVSB*LATSB*LONSB];
float indata[LEVSA][LATSA][LONSA], outdata[LEVSB][LATSB][LONSB];

main(int argc, char *argv[])
{
    FILE *fpa, *fpb, *fopen();
    int a, b, c, d, param;
    char infila[LEN], outfila[LEN], date[LEN], tag[LEN], year[LEN], line[LOTS];

    if (argc != 4) {
      printf("Usage: %s norskhavet_2010-2012 norskhavet_20100101030000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v01.0-fv01.0.nc norskhavet_2010-2012.surf.current.nc\n",argv[0]);
      exit(1);
    }
    strcpy(infila,argv[2]);
    strcpy(outfila,argv[3]);

    if ((fpa = fopen(argv[1], "r")) == NULL) {                      /* try opening datelist in the current dir */
      fprintf(stderr, "ERROR : couldn't open %s\n", argv[1]);
      exit(1);
    }
    printf("using datelist %s\n",argv[1]);

    while (fgets(line,LOTS,fpa) != NULL) {                          /* for each date in datelist */
      sscanf(line, "%*s %s", date);
      infila[11] = date[ 0] ; infila[12] = date[ 1] ; infila[13] = date[ 2];
      infila[14] = date[ 3] ; infila[15] = date[ 5] ; infila[16] = date[ 6];
      infila[17] = date[ 8] ; infila[18] = date[ 9] ; infila[19] = date[11];
      infila[20] = date[12];
/* 0123 56 89 12    0123456789012345678901234
   2010-01-01-03 -> norskhavet_20121229120000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v01.0-fv01.0.nc  */
/*    intread(  infila,  "eastward_eulerian_current_velocity", date, in, LEVSA*LATSA*LONSA);
      intread(  infila,  "eastward_geostrophic_current_velocity", date, in, LEVSA*LATSA*LONSA);  */
      intread(  infila,  "eastward_ekman_current_velocity",       date, in, LEVSA*LATSA*LONSA);
      cdcwrite(outfila,                                      "u", date, in, LEVSB*LATSB*LONSB);
/*    intread(  infila, "northward_geostrophic_current_velocity", date, in, LEVSA*LATSA*LONSA);  */
      intread(  infila, "northward_ekman_current_velocity",       date, in, LEVSA*LATSA*LONSA);
      cdcwrite(outfila,                                      "v", date, in, LEVSB*LATSB*LONSB);
/*    for (a = 0; a < LATSB * LONSB; a++)
        outa[a] = pow(outa[a] * outa[a] + outb[a] * outb[a], 0.5);
      cdcwrite(outfila, "w_Stokes", date, outa, LEVSB*LATSB*LONSB);  */
    }
}

void transfer(float in[LEVSA*LATSA*LONSA], float out[LEVSB*LATSB*LONSB])
{
    int a, b, c, d;

    d = 0;
    for (a = 0; a < LEVSA; a++)
      for (b = 0; b < LATSA; b++)
        for (c = 0; c < LONSA; c++)
          indata[a][b][c] = in[d++];

    for (a = 0; a < LEVSB; a++)
      for (b = 0; b < LATSB; b++) {
        for (c = 0; c < LONSA / 2; c++)
          outdata[a][b][c] = indata[a][b][c+LONSA/2];
        for (c = LONSA / 2; c < LONSB; c++)
          outdata[a][b][c] = indata[a][b][c-LONSA/2];
      }

    d = 0;
    for (a = 0; a < LEVSB; a++)
      for (b = 0; b < LATSB; b++)
        for (c = 0; c < LONSB; c++)
          out[d++] = outdata[a][b][c];
}

void intread(char file[], char param[], char date[], float *array, int arraylen)
{
    int match, loopa, loopb, datalen, mesgflag;
    double sarray[arraylen];
    double packedval;
    int status, ncid, dim_id, par_id, par_ndims, par_dims[NC_MAX_VAR_DIMS], par_natts;
    char dim_name[NC_MAX_NAME], dim_specs[NC_MAX_VAR_DIMS*NC_MAX_NAME];
    nc_type par_type;
    float offset, scale;
    size_t dim_len, origin[NC_MAX_VAR_DIMS], interval[NC_MAX_VAR_DIMS];

    if (strcmp(param,"") == 0) {                                              /* just return zeros if the */
      for (loopa = 0; loopa < arraylen; loopa++)                              /* parameter name is null */
        array[loopa] = 0.0;
      return;
    }

    status = nc_open(file, NC_NOWRITE, &ncid);                                /* open the data file */
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    status = nc_inq_varid(ncid, param, &par_id);                              /* get parameter information */
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_inq_var(ncid, par_id, 0, &par_type, &par_ndims, par_dims, &par_natts);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

    datalen = 1;                                                              /* get dimension information */
    strcpy(dim_specs,"");
    for (loopa = 0; loopa < par_ndims; loopa++) {
      status = nc_inq_dim(ncid, par_dims[loopa], dim_name, &dim_len);
      if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
      strcat(dim_specs,dim_name);
      strcat(dim_specs,"=");
      if (strcmp(dim_name,"time") == 0) {
        timematch(ncid, date, &match);
        origin[loopa] = match;
        interval[loopa] = 1;
        sprintf(dim_name,"%d",match+1);
        strcat(dim_specs,dim_name);
      }
      else {
        origin[loopa] = 0;
        interval[loopa] = dim_len;
        datalen *= (int)dim_len;
        sprintf(dim_name,"1-%d",(int)dim_len);
        strcat(dim_specs,dim_name);
      }
      if (loopa != par_ndims-1)
        strcat(dim_specs," ");
    }

    if (datalen > arraylen) {                                                 /* exit if the number of available */
      printf("ERROR: ");                                                      /* data values exceeds expectations */
      printf("%d %s values requested from %s ",datalen,param,file);
      printf("but array can only accommodate %d values\n",arraylen);
      exit(-1);
    }
    else if (datalen < arraylen) {                                            /* warn if the number of available */
      printf("WARNING: ");                                                    /* data values is less than expected */
      printf("%d %s values requested from %s ",datalen,param,file);
      printf("but array is expected to accommodate %d values\n",arraylen);
    }

    printf("reading %s %s %-12s %s\n",file,date,param,dim_specs);             /* read the packed grid */
    status = nc_get_vara_double(ncid, par_id, origin, interval, sarray);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}

/*  status = nc_get_att_float(ncid, par_id, "add_offset", &offset);           /* read the offset and scale factor *
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
    status = nc_get_att_float(ncid, par_id, "scale_factor", &scale);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}  */

/*  for (loopa = 0; loopa < datalen; loopa++)                                 /* and unpack the grid *
      array[loopa] = sarray[loopa] * scale + offset;  */

    mesgflag = 0;
    for (loopa = 0; loopa < datalen; loopa++) {                               /* and unpack the grid */
      if (sarray[loopa] < -32766) {                                           /* but if any packed data values */
        array[loopa] = MISSING;                                               /* exceeded the allowable range in */
        mesgflag++;                                                           /* values then set them to missing */
      }
      else
        array[loopa] = sarray[loopa];
    }

    if (mesgflag == 1) {                                                      /* and notify the user */
      printf("WARNING: ");
      printf("one %s value exceeded the valid packing range ",param);
      printf("so the corresponding unpacked value has been reset to MISSING (%f)\n",MISSING);
    }
    else if (mesgflag > 1) {
      printf("WARNING: ");
      printf("%d %s values exceeded the valid packing range ",mesgflag,param);
      printf("so all corresponding unpacked values have been reset to MISSING (%f)\n",MISSING);
    }

    status = nc_close(ncid);
    if (status != NC_NOERR) {printf("%s\n",nc_strerror(status));exit(-1);}
}
