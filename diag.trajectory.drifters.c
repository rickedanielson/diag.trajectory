/*
 * This program is designed to compute trajectories using a spherical coordinate,
 * near-surface timeseries.  Trajectories are launched every six hours beginning at
 * the location of an actual drifter, are limited to about 25 days in length (or to
 * the end of the drifter track, whichever is first), and are stored in individual
 * ascii files.  Position update from www.movable-type.co.uk/scripts/latlong.html
 * (destination given distance and bearing from start) is employed and following
 * Liu et al. [2014] (http://ocg6.marine.usf.edu/~liu/skill_score_drifter.html),
 * a skill score is included - RD November 2004, July, September 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "prog/include.netcdf/include/cdcnetcdf.h"

#define LEN        100
#define LOTS       15000
#define D2R        (3.141592654/180.0)       /* degrees to radians conversion */
#define EARTH      6.371e6                   /* mean earth radius (m) */
#define ONEDAY     24                        /* length of one day (hours) for calculating dates */
#define H2S        3600.0                    /* hours to seconds conversion */
#define MISS      -9999.0                    /* generic missing value */
#define INTVALS    6                         /* number of lats and lons within which to interpolate */

#define LATS       1600                      /* number of latitudes  in the calculation domain */
#define LONS       3600                      /* number of longitudes in the calculation domain */
#define FIRSTLAT  -79.95                     /* first latitude       in the calculation domain (degrees) */
#define FIRSTLON  -179.95                    /* first longitude      in the calculation domain (degrees) */
#define DELLAT     0.1                       /* difference in adjacent latitudes (degrees) */
#define DELLON     0.1                       /* difference in adjacent longitudes (degrees) */

#define UVEL       0                         /* zonal velocity */
#define VVEL       1                         /* meridional velocity */
#define PARAMS     2                         /* number of parameters */

#define INI        0                         /* initial (0-hr) time */
#define MID        1                         /* middle  (3-hr) time */
#define FIN        2                         /* final   (6-hr) time */
#define TIMES      3                         /* number of analyses used to interpolate in time */
#define DELTIM     3                         /* difference in adjacent analysis times (hours) */
#define TIMSTEPS   11                        /* number of times to interpolate between analyses */

#define DELDRIFT   6                         /* difference in adjacent drifter times (hours) */
#define TIMDRIFT   101                       /* maximum number of drifter times in one trajectory */

void traj(int a, int b);
void interp(double array[LATS][LONS], double lat, double lon, double *val);
void spline(double x[], double y[], int n, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void dateval2str(int  yr, int  mo, int  dy, int  hr, int  mi, char date[]);
void datestr2val(int *yr, int *mo, int *dy, int *hr, int *mi, char date[]);

int count;
double data[PARAMS][TIMES][LATS][LONS], array[LATS*LONS];                     /* declare globally to use all RAM */
double terp[PARAMS][LATS][LONS], lat[LOTS][TIMDRIFT], lon[LOTS][TIMDRIFT];
double obsdis[LOTS][TIMDRIFT], simdis[LOTS][TIMDRIFT], deldis[LOTS][TIMDRIFT];
double obssum[LOTS],           simsum[LOTS],           delsum[LOTS];

main(int argc, char *argv[])
{
    FILE *fpa, *fpb, *fopen();
    int a, b, c, d, e, f, outside, yr, mo, dy, hr, mi;
    long int datenow, datebef[LOTS], datemid[LOTS], dateaft[LOTS];
    char date[LOTS][LEN], tmpdate[LEN], line[LOTS], infila[LEN], infilb[LEN], infilc[LEN], outfile[LEN];
    char *parname[PARAMS] = {"eastward_eulerian_current_velocity", "northward_eulerian_current_velocity"};
    double deldate, dist, bear, skill;

    if (argc != 2) {
      printf("Usage: %s track\n", argv[0]);
      exit(1);
    }

    if ((fpa = fopen(argv[1], "r")) == NULL) {
      fprintf(stderr, "ERROR : couldn't open %s\n", argv[1]);
      exit(-1);
    }
    printf("reading %s\n", argv[1]);

    for (a = 0; a < LOTS; a++)                                                /* initialize all positions as valid and */
      for (b = 0; b < TIMDRIFT; b++)                                          /* define the start of each trajectory and */
        lat[a][b] = lon[a][b] = -361.0;                                       /* their (inclusive) maximum date bracket */

    count = 0;
    while (fgets(line, LOTS, fpa) != NULL) {
      sscanf(line, "%s %*s %*s %*s %*s %lf %lf", date[count], &lon[count][0], &lat[count][0]);
      strcpy(tmpdate, date[count]);
      datestr2val(&yr, &mo, &dy, &hr, &mi, tmpdate);
      datebef[count] = yr * 100000000L + mo * 1000000L + dy * 10000L + hr * 100L + mi;
      strcpy(line, tmpdate+13);
      tmpdate[13] = '\0' ; dateshift(tmpdate,                             3) ; strcat(tmpdate, line);
      datestr2val(&yr, &mo, &dy, &hr, &mi, tmpdate);
      datemid[count] = yr * 100000000L + mo * 1000000L + dy * 10000L + hr * 100L + mi;
      tmpdate[13] = '\0' ; dateshift(tmpdate, DELDRIFT * (TIMDRIFT - 1) - 3) ; strcat(tmpdate, line);
      datestr2val(&yr, &mo, &dy, &hr, &mi, tmpdate);
      dateaft[count] = yr * 100000000L + mo * 1000000L + dy * 10000L + hr * 100L + mi;
      printf("%s %ld %ld %ld\n", date[count], datebef[count], datemid[count], dateaft[count]);
      count++;
    }
    fclose(fpa);
    printf("%d track entries were found\n",count);
    if (count > LOTS) {
      fprintf(stderr, "ERROR : %d track entries exceeds the maximum expected (%d)\n", count, LOTS);
      exit(-1);
    }

    for (a = 1; a < count; a++) {                                             /* first check that the drifter */
      strcpy(infila, date[a-1]) ; infila[13] = '\0';                          /* dates are strictly six-hourly */
      strcpy(infilb, date[a  ]) ; infilb[13] = '\0';
      getdeldate(infila, infilb, &deldate);
      if (deldate != DELDRIFT * 60) {
        fprintf(stderr, "ERROR : %s and %s are not %d hours apart for drifter %s\n\n", date[a-1], date[a], DELDRIFT, argv[1]);
        exit(-1);
      }
    }

    datenow = datebef[count-1];
    for (a = 0; a < count; a++) {                                             /* then constrain the end bracket date */
      for (b = count; b < TIMDRIFT + a; b++)                                  /* (set any positions beyond this date */
        lat[a][b-a] = lon[a][b-a] = MISS;                                     /*  to missing) */
      if (count - a < TIMDRIFT) dateaft[a] = datenow;
      printf("%s %ld %ld %ld\n", date[a], datebef[a], datemid[a], dateaft[a]);
    }

    sprintf(infila, "../links/%ld", datebef[0]);                              /* read the first INI grids into FIN */
    strcat(infila, "00-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v01.0-fv01.0.nc");
    strcpy(tmpdate, date[0]) ; tmpdate[13] = '\0';
    for (b = 0; b < PARAMS; b++) {
      cgcread(infila, parname[b], tmpdate, array, LATS*LONS);
      e = 0 ; for (c = 0; c < LATS; c++) for (d = 0; d < LONS; d++) data[b][FIN][c][d] = array[e++];
    }

    for (a = 0; a < count - 1; a++) {                                         /* then for each date interval shift FIN */
      sprintf(infilb, "../links/%ld", datemid[a]);                            /* to INI and read the new MID/FIN data */
      sprintf(infilc, "../links/%ld", datebef[a+1]);                          /* and advance all valid trajectories */
      strcat(infilb, "00-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v01.0-fv01.0.nc");
      strcat(infilc, "00-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v01.0-fv01.0.nc");
      for (b = 0; b < PARAMS; b++) {
        e = 0 ; for (c = 0; c < LATS; c++) for (d = 0; d < LONS; d++) data[b][INI][c][d] = data[b][FIN][c][d];
        strcpy(tmpdate, date[a+1]) ; tmpdate[13] = '\0' ; dateshift(tmpdate, -3);
        cgcread(infilb, parname[b], tmpdate, array, LATS*LONS);
        e = 0 ; for (c = 0; c < LATS; c++) for (d = 0; d < LONS; d++) data[b][MID][c][d] = array[e++];
        strcpy(tmpdate, date[a+1]) ; tmpdate[13] = '\0';
        cgcread(infilc, parname[b], tmpdate, array, LATS*LONS);
        e = 0 ; for (c = 0; c < LATS; c++) for (d = 0; d < LONS; d++) data[b][FIN][c][d] = array[e++];
      }

      for (b = 0; b < a + 1; b++)                                             /* by looping through each trajectory that */
        if (a + 1 - b < TIMDRIFT) {                                           /* remains valid and advancing one drifter */
          if (lat[b][a+1-b-1] > -362 && lat[b][a+1-b] > -362)                 /* timestep (any position that is invalid */
            traj(b, a+1-b);                                                   /* initially remains so; only positions */
          else                                                                /* equal to -361 are updated) */
            lat[b][a+1-b] = lon[b][a+1-b] = MISS;
        }
    }

    for (a = 0; a < count; a++) {                                             /* calculate incremental distance measures */
      obsdis[a][0] = 0.0;
      simdis[a][0] = 0.0;
      deldis[a][0] = 0.0;
      for (b = a + 1; b < count && b < TIMDRIFT + a; b++) {
        bear = cos(lat[b-1][    0]*D2R) * cos(lat[b][  0]*D2R) * cos(lon[b-1][    0]*D2R - lon[b][  0]*D2R) + sin(lat[b-1][    0]*D2R) * sin(lat[b][  0]*D2R) ; if (bear > 1) bear = 1.0 ; obsdis[a][b-a] = EARTH * acos(bear);
        bear = cos(lat[a  ][b-a-1]*D2R) * cos(lat[a][b-a]*D2R) * cos(lon[a  ][b-a-1]*D2R - lon[a][b-a]*D2R) + sin(lat[a  ][b-a-1]*D2R) * sin(lat[a][b-a]*D2R) ; if (bear > 1) bear = 1.0 ; simdis[a][b-a] = EARTH * acos(bear);
        bear = cos(lat[a  ][b-a  ]*D2R) * cos(lat[b][  0]*D2R) * cos(lon[a  ][b-a  ]*D2R - lon[b][  0]*D2R) + sin(lat[a  ][b-a  ]*D2R) * sin(lat[b][  0]*D2R) ; if (bear > 1) bear = 1.0 ; deldis[a][b-a] = EARTH * acos(bear);
      }
    }

    for (a = 0; a < count; a++) {                                             /* include distance and bearing in output */
      strcpy(outfile, argv[1]);                                               /* (i.e., from the drifter at index [b][0] */
      strcat(outfile, ".traj.");                                              /*  to the calculated trajectory position */
      strcat(outfile, date[a]);                                               /*  at [a][b-a]) as well as a skill score */
      if ((fpb = fopen(outfile, "w")) == NULL) {
        fprintf(stderr, "ERROR : couldn't open %s\n", outfile);
        exit(-1);
      }
      printf("writing %s\n", outfile);
/*    obssum[0] = 0.0;
      simsum[0] = 0.0;
      delsum[0] = 0.0;
      for (b = a + 1; b < count && b < TIMDRIFT + a; b++) {
        obssum[b-a] = obssum[b-a-1] + obsdis[a][b-a];
        simsum[b-a] = simsum[b-a-1] + simdis[a][b-a];
        delsum[b-a] = delsum[b-a-1] + deldis[a][b-a];  */
      for (b = a; b < count - 11 && b < TIMDRIFT + a - 11; b++) {
/*      obssum[b-a] = 3.0 * (obsdis[a][b-a  ] + obsdis[a][b-a+1] + obsdis[a][b-a+ 2] + obsdis[a][b-a+ 3]) +
                      2.0 * (obsdis[a][b-a+4] + obsdis[a][b-a+5] + obsdis[a][b-a+ 6] + obsdis[a][b-a+ 7]) +
                             obsdis[a][b-a+8] + obsdis[a][b-a+9] + obsdis[a][b-a+10] + obsdis[a][b-a+11];
        simsum[b-a] = 3.0 * (simdis[a][b-a  ] + simdis[a][b-a+1] + simdis[a][b-a+ 2] + simdis[a][b-a+ 3]) +
                      2.0 * (simdis[a][b-a+4] + simdis[a][b-a+5] + simdis[a][b-a+ 6] + simdis[a][b-a+ 7]) +
                             simdis[a][b-a+8] + simdis[a][b-a+9] + simdis[a][b-a+10] + simdis[a][b-a+11];
        delsum[b-a] =        deldis[a][b-a  ] + deldis[a][b-a+1] + deldis[a][b-a+ 2] + deldis[a][b-a+ 3] +
                             deldis[a][b-a+4] + deldis[a][b-a+5] + deldis[a][b-a+ 6] + deldis[a][b-a+ 7] +
                             deldis[a][b-a+8] + deldis[a][b-a+9] + deldis[a][b-a+10] + deldis[a][b-a+11];  */
        obssum[b-a] = 12.0 * obsdis[a][b-a  ] + 11.0 * obsdis[a][b-a+1] + 10.0 * obsdis[a][b-a+ 2] + 9.0 * obsdis[a][b-a+ 3] +
                       8.0 * obsdis[a][b-a+4] +  7.0 * obsdis[a][b-a+5] +  6.0 * obsdis[a][b-a+ 6] + 5.0 * obsdis[a][b-a+ 7] +
                       4.0 * obsdis[a][b-a+8] +  3.0 * obsdis[a][b-a+9] +  2.0 * obsdis[a][b-a+10] +       obsdis[a][b-a+11];
        simsum[b-a] = 12.0 * simdis[a][b-a  ] + 11.0 * simdis[a][b-a+1] + 10.0 * simdis[a][b-a+ 2] + 9.0 * simdis[a][b-a+ 3] +
                       8.0 * simdis[a][b-a+4] +  7.0 * simdis[a][b-a+5] +  6.0 * simdis[a][b-a+ 6] + 5.0 * simdis[a][b-a+ 7] +
                       4.0 * simdis[a][b-a+8] +  3.0 * simdis[a][b-a+9] +  2.0 * simdis[a][b-a+10] +       simdis[a][b-a+11];
        delsum[b-a] =        deldis[a][b-a  ] +        deldis[a][b-a+1] +        deldis[a][b-a+ 2] +       deldis[a][b-a+ 3] +
                             deldis[a][b-a+4] +        deldis[a][b-a+5] +        deldis[a][b-a+ 6] +       deldis[a][b-a+ 7] +
                             deldis[a][b-a+8] +        deldis[a][b-a+9] +        deldis[a][b-a+10] +       deldis[a][b-a+11];
      }
      for (     ; b < count      && b < TIMDRIFT + a     ; b++)
        obssum[b-a] = simsum[b-a] = delsum[b-a] = 0.0;
/*    for (b = a + 1; b < count && b < TIMDRIFT + a; b++) {
        obssum[b-a] += obssum[b-a-1];
        simsum[b-a] += simsum[b-a-1];
      }
      obssum[0] = 1e-9;
      simsum[0] = 1e-9;  */

      for (b = a; b < count && b < TIMDRIFT + a; b++) {
        bear = cos(lat[a][b-a]*D2R) * cos(lat[b][0]*D2R) * cos(lon[a][b-a]*D2R - lon[b][0]*D2R) + sin(lat[a][b-a]*D2R) * sin(lat[b][0]*D2R);
        if (bear > 1) bear = 1.0;
        dist = EARTH * acos(bear);
/*      dist = EARTH * acos(cos(lat[a][b-a]*D2R) * cos(lat[b][0]*D2R) *
                            cos(lon[a][b-a]*D2R  -     lon[b][0]*D2R) +
                            sin(lat[a][b-a]*D2R) * sin(lat[b][0]*D2R));  */
        bear = atan2(sin(lon[a][b-a]*D2R - lon[b][0]*D2R) * cos(lat[a][b-a]*D2R), cos(lat[b][0]*D2R) * sin(lat[a][b-a]*D2R) -
                     cos(lon[a][b-a]*D2R - lon[b][0]*D2R)                       * sin(lat[b][0]*D2R) * cos(lat[a][b-a]*D2R)) / D2R;
        skill = 1.0 - 2.0 * delsum[b-a] / (obssum[b-a] + simsum[b-a]);
        if ((obssum[b-a] + simsum[b-a]) < 1e-9 || skill < 0) skill = 0.0;
/*      skill = 1.0 - delsum[b-a] / obssum[b-a] ; if (skill < 0) skill = 0.0;  */
        fprintf(fpb, "%s %13.5f %13.5f %5d %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f\n",
          date[b], lat[b][0], lon[b][0], (b - a) * DELDRIFT, lat[a][b-a], lon[a][b-a], dist, bear, skill, obssum[b-a], simsum[b-a], delsum[b-a]);
      }
      fclose(fpb);
    }
    exit(0);
}

void traj(int aa, int bb)
{
    int a, b, c, d, e, f, latind, lonind;
    double latt, lonn, weight, duration, vala, valb, dist, bear;

    latt = lat[aa][bb-1];                                                     /* start with initial position and */
    lonn = lon[aa][bb-1] ; if (lonn > 180) lonn -= 360.0;                     /* its nearest gridpoint position */
    latind = lonind = -1;
    for (a = 0; a < LATS; a++)
      if (FIRSTLAT + (a - 0.5) * DELLAT < latt && FIRSTLAT + (a + 0.5) * DELLAT >= latt) {
        latind = a ; a = LATS;
      }
    for (a = 0; a < LONS; a++)
      if (FIRSTLON + (a - 0.5) * DELLON < lonn && FIRSTLON + (a + 0.5) * DELLON >= lonn) {
        lonind = a ; a = LONS;
      }
    if (latind == -1 || lonind == -1) {
      fprintf(stderr, "ERROR: lat %f or lon %f is outside the interpolation range\n", latt, lonn);
      exit(-2);
    }

    for (e = 0; e < 2; e++)                                                   /* loop through two 3-h periods; */
      for (f = 0; f < TIMSTEPS + 2; f++) {                                    /* for each time step interpolate */
        weight = (double)f / (double)(TIMSTEPS + 1);                          /* linearly in time (on the first */
        duration = (double)DELTIM * H2S / (double)(TIMSTEPS + 1);             /* and last steps divide duration */
        if (f == 0 || f == TIMSTEPS + 1)                                      /* in half) but in space simply use */
          duration *= 0.5;                                                    /* the nearest gridpoint position */
/*      printf("advancing by %13.5lf hours ", duration / H2S);  */

        if (e == 0) {
          vala = (1.0 - weight) * data[UVEL][INI][latind][lonind] + weight * data[UVEL][MID][latind][lonind];
          valb = (1.0 - weight) * data[VVEL][INI][latind][lonind] + weight * data[VVEL][MID][latind][lonind];
        }
        else {
          vala = (1.0 - weight) * data[UVEL][MID][latind][lonind] + weight * data[UVEL][FIN][latind][lonind];
          valb = (1.0 - weight) * data[VVEL][MID][latind][lonind] + weight * data[VVEL][FIN][latind][lonind];
        }
        dist = pow(vala * vala + valb * valb, 0.5) * duration;
/*      printf("and %13.5lf meters... %lf %lf %lf\n", dist, vala, valb, weight);  */

        latt *= D2R;
        dist /= EARTH;
        bear  = atan2(vala, valb);
        vala  =  asin(cos(bear) * sin(dist) * cos(latt) + cos(dist) * sin(latt));
        valb  = atan2(sin(bear) * sin(dist) * cos(latt),  cos(dist) - sin(latt) * sin(vala));
        latt  = vala / D2R;
        lonn += valb / D2R;

        latind = lonind = -1;                                                 /* update nearest gridpoint position */
        for (a = 0; a < LATS; a++)
          if (FIRSTLAT + (a - 0.5) * DELLAT < latt && FIRSTLAT + (a + 0.5) * DELLAT >= latt) {
            latind = a ; a = LATS;
          }
        for (a = 0; a < LONS; a++)
          if (FIRSTLON + (a - 0.5) * DELLON < lonn && FIRSTLON + (a + 0.5) * DELLON >= lonn) {
            lonind = a ; a = LONS;
          }

        if (latind == -1 || lonind == -1) {                                   /* but invalidate the rest of the */
          latt = lonn = MISS ; e = f = LOTS;                                  /* trajectory where it leads out */
        }                                                                     /* of the ocean */
        else if (data[UVEL][INI][latind][lonind] < MISS || data[VVEL][INI][latind][lonind] < MISS ||
                 data[UVEL][MID][latind][lonind] < MISS || data[VVEL][MID][latind][lonind] < MISS ||
                 data[UVEL][FIN][latind][lonind] < MISS || data[VVEL][FIN][latind][lonind] < MISS) {
          latt = lonn = MISS ; e = f = LOTS;
        }
      }

                                                  lat[aa][bb] = latt;         /* and update final position */
    if (lonn != MISS && lonn < 0) lonn += 360.0 ; lon[aa][bb] = lonn;
}

void interp(double array[LATS][LONS], double lat, double lon, double *val)
{
    int a, b, latind, lonind;
    double latcen, loncen, lonints[INTVALS];
    double x[LOTS], y[LOTS], y2[LOTS], misslat, misslon, missval;

    latind = lonind = 0;                                                      /* locate the position in the grid */
    for (a = INTVALS/2; a < LATS-INTVALS/2; a++)
      if (FIRSTLAT - (a-1) * DELLAT > lat && FIRSTLAT - a * DELLAT <= lat) {
        latcen = FIRSTLAT - a * DELLAT;
        latind = a;
      }
    for (a = INTVALS/2; a < LONS-INTVALS/2; a++)
      if (FIRSTLON + (a-1) * DELLON < lon && FIRSTLON + a * DELLON >= lon) {
        loncen = FIRSTLON + a * DELLON;
        lonind = a;
      }

    if (latind == 0 || lonind == 0) {                                         /* quit if interpolation isn't needed */
      printf("ERROR: lat %f or lon %f is outside the interpolation range\n",lat,lon);
      exit(-2);
    }
    if (lat == FIRSTLAT - latind * DELLAT && lon == FIRSTLON + lonind * DELLON) {
      *val = array[latind][lonind];
      return;
    }

    misslon = lon;                                                            /* start interpolating longitudinally */
    for (a = 0; a < INTVALS; a++)                                             /* set the local grid longitudes */
      x[a] = FIRSTLON + (lonind+a-INTVALS/2) * DELLON;

    for (a = 0; a < INTVALS; a++) {                                           /* and loop through the latitudes */
      for (b = 0; b < INTVALS; b++)
        y[b] = array[latind+a-INTVALS/2][lonind+b-INTVALS/2];
      spline(x,y,INTVALS,y2);
      splint(x,y,y2,INTVALS,misslon,&missval);
      lonints[a] = missval;
    }

    misslat = lat;                                                            /* and finish with the latitudinal */
    for (a = 0; a < INTVALS; a++) {                                           /* interpolation */
      x[INTVALS-1-a] = FIRSTLAT - (latind+a-INTVALS/2) * DELLAT;
      y[INTVALS-1-a] = lonints[a];
    }
    spline(x,y,INTVALS,y2);
    splint(x,y,y2,INTVALS,misslat,&missval);
    *val = missval;
}

void spline(double x[LOTS], double y[LOTS], int n, double y2[LOTS])
{
    int i, k;
    double p, sig, u[LOTS];

    y2[0] = y2[n-1] = 0.0;
    u[0] = 0.0;

    for (i = 1; i < n-1; i++) {
      sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
      p = sig * y2[i-1] + 2.0;
      y2[i] = (sig - 1.0) / p;
      u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
      u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
    }

    for (k = n-2; k >=0; k--)
      y2[k] = y2[k] * y2[k+1] + u[k];
}

void splint(double xa[LOTS], double ya[LOTS], double y2a[LOTS], int n, double x, double *y)
{
    int k, klo, khi;
    double a, b, h;

    klo = 0;
    khi = n-1;

    while (khi - klo > 1) {
      k = (khi + klo) >> 1;
      if (xa[k] > x) khi = k;
      else klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0) {
      printf("ERROR: bad splint...\n");
      exit(1);
    }

    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    *y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6.0;
}

void dateval2str(int  yr, int  mo, int  dy, int  hr, int  mi, char date[])
{
    char cyr[5], cmo[3], cdy[3], chr[3], cmi[3];

                sprintf(cyr,"%d",yr);
    if (mo > 9) sprintf(cmo,"%d",mo) ; else sprintf(cmo,"0%d",mo);
    if (dy > 9) sprintf(cdy,"%d",dy) ; else sprintf(cdy,"0%d",dy);
    if (hr > 9) sprintf(chr,"%d",hr) ; else sprintf(chr,"0%d",hr);
    if (mi > 9) sprintf(cmi,"%d",mi) ; else sprintf(cmi,"0%d",mi);
    sprintf(date,"%s-%s-%s-%s%s",cyr,cmo,cdy,chr,cmi);
}

void datestr2val(int *yr, int *mo, int *dy, int *hr, int *mi, char date[])
{
    char tmp[LEN];
    int pos, chs;

    *yr = *mo = *dy = *hr = *mi = 0;
    pos =  0 ; chs = 4 ; strncpy(tmp,date+pos,chs) ; tmp[chs] = '\0' ; *yr = atoi(tmp);
    pos =  5 ; chs = 2 ; strncpy(tmp,date+pos,chs) ; tmp[chs] = '\0' ; *mo = atoi(tmp);
    pos =  8 ; chs = 2 ; strncpy(tmp,date+pos,chs) ; tmp[chs] = '\0' ; *dy = atoi(tmp);
    pos = 11 ; chs = 2 ; strncpy(tmp,date+pos,chs) ; tmp[chs] = '\0' ; *hr = atoi(tmp);
    pos = 13 ; chs = 2 ; strncpy(tmp,date+pos,chs) ; tmp[chs] = '\0' ; *mi = atoi(tmp);
}
