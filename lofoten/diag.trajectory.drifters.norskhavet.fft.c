/*
 * This program is designed to compute a Fourier decomposition and
 * filter it to obtain a lowpass time series and highpass anomalies
 * - RD October, November 2002, October 2007, January 2008, October
 * 2015
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "prog/include.netcdf/include/genalloc.h"
#include "prog/include.netcdf/include/cdcnetcdf.h"
#include "prog/include.netcdf/include/dateshift.h"
#include "netcdf.h"

#define LEN        100
#define LOTS       10000
#define PI         3.141592654
#define D2R        (3.141592654/180.0)       /* degrees to radians conversion */
#define SHORTLIM   32766                     /* packed data value limit */

#define LATS       100                       /* number of domain latitudes */
#define LONS       300                       /* number of domain longitudes */
#define FIRSTLAT   65.05                     /* first latitude  in the domain (degrees) */
#define FIRSTLON  -9.95                      /* first longitude in the domain (degrees) */
#define DELLAT     0.1                       /* difference in adjacent latitudes (degrees) */
#define DELLON     0.1                       /* difference in adjacent longitudes (degrees) */
#define MISS       -9999.0                   /* generic missing data value */
#define DELMIS     -8888.0                   /* generic missing data value comparison */
#define BIGMIS     -9e9                      /* generic missing data value for output */
#define BUTORD     5.0                       /* Butterworth filter order */

#define RAW        0                         /* timeseries of unfiltered data */
#define LOW        1                         /* timeseries and spectral coefficients of lowpass filtered data */
#define HIG        2                         /* timeseries and spectral coefficients of higpass filtered data */
#define PARAMS     3                         /* number of parameters */

void fft(int freqs, int isign, int var);

float ****grid;
double **data;

main(int argc, char *argv[])
{
    FILE *fpa, *fpb, *fopen();
    int a, b, c, d, e, count, latind, lonind, times, freqs, var;
    char date[LOTS][LEN], infile[LEN], outfila[LEN], outfilb[LEN], parname[LEN], line[LOTS], command[LOTS];
    double timean, cutlow, cuthig, filtlow, filthig;
    float array[LATS*LONS], lata, lona;

    int mesgflag, time_dim;  short *sarray;  double packedval, badval;
    int status, ncid, ncida, ncidb, dim_id, par_id, par_ndims, par_dims[3], par_natts;
    char dim_name[NC_MAX_NAME];  nc_type par_type;  float offset, scale;
    size_t dim_len, origin[3], interval[3];

/*  double x[2*CONS], y[2*CONS], y2[2*CONS], ans, rrr[LEN];
    fftw_plan f10, f15, f20, f24, f30, f36, f40, f45;
    fftw_plan b10, b15, b20, b24, b30, b36, b40, b45;
    fftw_complex ccc[LEN];  */

    if (argc != 5) {
      fprintf(stderr,"Usage: %s datelist data_end cutoff_low-hig parname\n",argv[0]);
      fprintf(stderr," e.g.: %s norskhavet_2010-2012 .surf.current 030 ucomb\n",argv[0]);
      fprintf(stderr,"       and cutoff timescales 090 010 are for 90-d and 10-d cutoffs\n");
      exit(1);
    }
    sscanf(argv[3],"%lf",&cutlow) ; /* cutlow *= 8.0; */
    strcpy(parname, argv[4]);

    sprintf(infile,"%s%s.nc",argv[1],argv[2]);                                /* create the output files */
    sprintf(outfila,"%s%s.low.%s.nc",argv[1],argv[2],argv[3]);
    sprintf(outfilb,"%s%s.hig.%s.nc",argv[1],argv[2],argv[3]);
/*  sprintf(command,"cp %s %s\n",infile,outfila) ; printf("%s",command) ; system(command);
    sprintf(command,"cp %s %s\n",infile,outfilb) ; printf("%s",command) ; system(command);  */

    if ((fpa = fopen(argv[1],"r")) == NULL) {                                 /* read the datelist */
      fprintf(stderr, "ERROR: couldn't open %s\n",argv[1]);
      exit(-1);
    }
    count = 0;
    while (fgets(line,LOTS,fpa) != NULL) {
      sscanf(line, "%*s %s", date[count]);
      count++;
    }
    fclose(fpa);
    printf("%d dates were found\n",count);
    times = count;

    freqs = pow(2.0, ceil(log((double)times) / log(2.0)));                    /* compute the next higher power of */
    printf("times = %d and freqs = %d\n", times, freqs);                      /* two as the number of frequencies */
    get_mem4Dfloat(&grid, PARAMS, times, LATS, LONS);
    get_mem2Ddouble(&data, PARAMS, 2 * freqs);

    d = 0;
    for (a = 0; a < times; a++) {
      cdcread(infile, parname, date[a], array, LATS*LONS);
      d = 0 ; for (b = 0; b < LATS; b++)  for (c = 0; c < LONS; c++)  grid[RAW][a][b][c] = array[d++];
    }

/*  printf("filtering at");  */
    for (a = 0; a < LATS; a++)                                                /* filter the entire series */
      for (b = 0; b < LONS; b++) {                                            /* at each valid location */
        lata = FIRSTLAT + a * DELLAT;
        lona = FIRSTLON + b * DELLON;
        timean = 0;
        mesgflag = 0;
        for (c = 0; c < times; c++) {
          if (grid[RAW][c][a][b] < DELMIS) {
            mesgflag = 1;
            c = times;
          }
          else {
            data[RAW][2*c] = grid[RAW][c][a][b];
            data[RAW][2*c+1] = 0.0;
            timean += data[RAW][2*c];
          }
        }

        if (mesgflag == 0) {
/*        printf(" (%8.2f, %8.2f)",lata,lona);  */                            /* if the timeseries is */
          timean /= (double)times;                                            /* valid, then remove its */
          for (c = 0; c < times; c++) {                                       /* mean, multiply it by */
            data[LOW][2*c] = data[RAW][2*c] - timean;                         /* a Bartlett window, and */
            data[LOW][2*c+1] = 0.0;                                           /* pad to freqs with zeros */
          }
          for (c = 0; c < times; c++)
            data[LOW][2*c] *= 1.0 - fabs(2.0 / ((double)times + 1.0) * ((double)c + 1.0) - 1.0);
          for (c = times; c < freqs; c++)
            data[LOW][2*c] = data[LOW][2*c+1] = 0.0;

          fft(freqs, 1, LOW);                                                 /* compute the spectrum */

          data[LOW][0] *= 1.0;                                                /* filter high frequencies */
          data[LOW][1] *= 1.0;
          for (c = 1; c < freqs/2; c++) {
            filtlow = pow(1.0 + pow((double)c * cutlow / (double)freqs, 2.0 * BUTORD), -0.5);
            data[LOW][2*c]           *= filtlow;
            data[LOW][2*c+1]         *= filtlow;
            data[LOW][2*freqs-2*c]   *= filtlow;
            data[LOW][2*freqs-2*c+1] *= filtlow;
          }
          c = freqs/2;
          filtlow = pow(1.0 + pow((double)c * cutlow / (double)freqs, 2.0 * BUTORD), -0.5);
          data[LOW][freqs]   *= filtlow;
          data[LOW][freqs+1] *= filtlow;

          fft(freqs, -1, LOW);                                                /* get the filtered timeseries */

          for (c = 0; c < times; c++) {                                       /* un-window, and add mean */
            data[LOW][2*c] /= 1.0 - fabs(2.0 / ((double)times + 1.0) * ((double)c + 1.0) - 1.0);
            data[LOW][2*c] += timean;
          }

          for (c = 0; c < times; c++)                                         /* compute highpass as anomaly */
            data[HIG][2*c] = data[RAW][2*c] - data[LOW][2*c];

          for (c = 0; c < times; c++) {
            grid[LOW][c][a][b] = data[LOW][2*c];
            grid[HIG][c][a][b] = data[HIG][2*c];
          }
        }
        else
          for (c = 0; c < times; c++)
            grid[LOW][c][a][b] = grid[HIG][c][a][b] = BIGMIS;
      }
    printf("\n");

    for (a = 0; a < times; a++) {
      d = 0 ; for (b = 0; b < LATS; b++)  for (c = 0; c < LONS; c++)  array[d++] = grid[LOW][a][b][c];
      cdcwrite(outfila, parname, date[a], array, LATS*LONS);
      d = 0 ; for (b = 0; b < LATS; b++)  for (c = 0; c < LONS; c++)  array[d++] = grid[HIG][a][b][c];
      cdcwrite(outfilb, parname, date[a], array, LATS*LONS);
    }

    free_mem4Dfloat(grid,PARAMS,times);
    free_mem2Ddouble(data);
    exit(0);
}

void fft(int freqs, int isign, int var)
{
    unsigned long nn, n, mmax, m, j, istep, i;
    double vals[2*freqs+1], wr, wi, wpr, wpi, wtemp, theta, tempr, tempi;

    nn = freqs;
    for (i = 0; i < 2*freqs; i++)
      vals[i+1] = data[var][i];

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {                                              /* perform bit reversal */
      if (j > i) {                                                            /* and exchange complex */
        tempr = vals[j];    vals[j]   = vals[i];    vals[i]   = tempr;        /* numbers */
        tempr = vals[j+1];  vals[j+1] = vals[i+1];  vals[i+1] = tempr;
      }
      m = nn;
      while (m >= 2 && j > m) {
        j -= m;
        m >>= 1;
      }
      j += m;
    }

    mmax = 2;
    while (n > mmax) {                                                        /* do outer loop log2(nn) times */
      istep = mmax << 1;                                                      /* initialize the trigonometric */
      theta = isign * (6.28318530717959 / mmax);                              /* recurrence */
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (m = 1; m < mmax; m += 2) {
        for (i = m; i <= n; i += istep) {
          j = i + mmax;
          tempr = wr * vals[j]   - wi * vals[j+1];
          tempi = wr * vals[j+1] + wi * vals[j];
          vals[j]    = vals[i]   - tempr;
          vals[j+1]  = vals[i+1] - tempi;
          vals[i]   += tempr;
          vals[i+1] += tempi;
        }
        wr = (wtemp=wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
      }
      mmax = istep;
    }

    if (isign == 1)  for (i = 0; i < 2*freqs; i++)  data[var][i] = vals[i+1];
    else             for (i = 0; i < 2*freqs; i++)  data[var][i] = vals[i+1] / (double)freqs;
}

/*
 * This subroutine is from Numerical Recipes and either replaces
 * vals[1..2*nn] by its discrete Fourier transform (if isign is 1)
 * or replaces vals[1. .2*nn] by nn times its inverse discrete Fourier
 * transform (if isign is -1).  Data is a complex array of length nn,
 * or equivalently, a double precision array of length 2*nn, where nn
 * MUST be an integer power of two (and this is not checked for).
 */


/*  e = 40;  printf("%f %f %f %f\n", grid[RAW][a][e][e], grid[LOW][a][e][e] + grid[HIG][a][e][e], grid[LOW][a][e][e], grid[HIG][a][e][e]);  */

/*  strcpy(outfila, "norskhavet_2010-2012.surf.current.low");
    if ((fpa = fopen(outfila, "w")) == NULL) {
      fprintf(stderr, "ERROR: couldn't open %s\n", outfila);
      exit(-1);
    }
    printf("writing %s\n", outfila);
    for (a = 0; a < times; a++)
      for (b = 0; b < LATS; b++) {
        for (c = 0; c < LONS; c++)
          printf("%f ", grid[RAW][a][b][c]);
        printf("\n");
      }
    fclose(fpa);  */
