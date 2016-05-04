/*
 * Create NFFT power spectral coefficients (one-sided) at locations where
 * analyses provide good coverage of the 2001-2007 period - RD October 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include "nfft3util.h"
#include "nfft3.h"

#define LEN        100
#define LOTS       500
#define DAYS       3408                      /* number of days between 1999-10-01 and 2009-12-31 */
#define TIMES      2920                      /* number of days between 2001-01-01 and 2007-12-31 (and spectral coefficients) */
#define MISS      -9999.0                    /* generic missing value */

main (int argc, char *argv[])
{
    FILE *fpa, *fpb, *fpc, *fpd, *fopen();
    int a, b, c, d, count, lines, obsnum;
    char infile[LEN], outfile[LEN], line[LOTS];
    float tmp, tor[DAYS], flx[DAYS], suma, sumb;
    nfft_plan p;

    if (argc != 2) {
      printf("Usage: %s seaflux....65.500....-5.000.fft\n",argv[0]);
      exit(1);
    }
    strcpy( infile, argv[1]);                                                 /* define the in/out files */
    strcpy(outfile, argv[1]);
    strcat(outfile, "est");

    suma = sumb = 0.0;
/*  printf("reading %s\n", infile);  */                                       /* then read the flux timeseries */
    if ((fpa = fopen(infile,"r")) == NULL) {
      fprintf(stderr, "ERROR : couldn't open %s\n", infile);
      exit(1);
    }
    a = 0;
    while (fgets(line,LOTS,fpa) != NULL) {
      sscanf(line, "%f %f", &tor[a], &flx[a]);
      suma += powf(flx[a], 2.0);
      a++;
    }
    fclose(fpa);
    obsnum = a;
/*  printf("%d %f\n", obsnum, suma / (float)obsnum); */

    nfft_init_1d(&p, TIMES, obsnum);                                          /* initialize a one-dimensional plan */
    for (a = 0; a < obsnum; a++)                                              /* and pass the [-0.5,0.5) torus points */
      p.x[a] = tor[a];
    if (p.nfft_flags & PRE_ONE_PSI)                                           /* precompute psi, the entries of the matrix B */
      nfft_precompute_one_psi(&p);
    for (a = 0; a < obsnum; a++) {                                            /* and pass the [-0.5,0.5) torus points */
      p.f[a][0] = flx[a];
      p.f[a][1] = 0.0;
    }
    nfft_adjoint(&p);

/*  for (a = 0; a < obsnum; a++)
      printf("%f %f\n", p.f_hat[a][0], p.f_hat[a][1]);  */

/*  printf("writing %s\n", outfile);  */                                      /* and store the power spectrum (in dB) */
    if ((fpa = fopen(outfile,"w")) == NULL) {                                 /* (note that FFT requires normalization by */
      fprintf(stderr, "ERROR : couldn't open %s\n", outfile);                 /*  TIMES^2 but NFFT employs TIMES*obsnum */
      exit(1);                                                                /*  to satisfy Parseval's equation) */
    }
    for (a = 0; a <= TIMES/2; a++) {
      if (a == 0)
        tmp = powf(p.f_hat[TIMES/2  ][0], 2.0) + powf(p.f_hat[TIMES/2  ][1], 2.0);
      else if (a == TIMES/2)
        tmp = powf(p.f_hat[0        ][0], 2.0) + powf(p.f_hat[0        ][1], 2.0);
      else
        tmp = powf(p.f_hat[TIMES/2+a][0], 2.0) + powf(p.f_hat[TIMES/2+a][1], 2.0) +
              powf(p.f_hat[TIMES/2-a][0], 2.0) + powf(p.f_hat[TIMES/2-a][1], 2.0);
/*    tmp /= powf((float)TIMES, 2.0);  */
      tmp /= (float)TIMES * (float)obsnum;
      fprintf(fpa, "%15.8f %15.8f\n", (float)a / (float)TIMES, tmp);
/*    sumb += tmp;
      fprintf(fpa, "%15.8f %15.8f\n", (float)a / (float)TIMES, 10.0 * log10(tmp));  */
    }
    fclose(fpa);
/*  printf("%d %f\n", TIMES / 2 + 1, sumb); */
    exit(0);
}


/*  nfft_vpr_complex(p.f_hat, p.N_total, "adjoint nfft, vector f_hat");

     ** init pseudo random Fourier coefficients and show them *
    nfft_vrand_unit_complex(p.f_hat, p.N_total);
    nfft_vpr_complex(p.f_hat, p.N_total, "given Fourier coefficients, vector f_hat");

     ** direct trafo and show the result *
    nfft_trafo_direct(&p);
    nfft_vpr_complex(p.f, p.M_total, "ndft, vector f");

     ** approx. trafo and show the result *
    nfft_trafo(&p);
    nfft_vpr_complex(p.f, p.M_total, "nfft, vector f");

     ** approx. adjoint and show the result *
    nfft_adjoint_direct(&p);
    nfft_vpr_complex(p.f_hat, p.N_total, "adjoint ndft, vector f_hat");

     ** approx. adjoint and show the result *
    nfft_adjoint(&p);
    nfft_vpr_complex(p.f_hat, p.N_total, "adjoint nfft, vector f_hat");

     ** finalise the one dimensional plan *
    nfft_finalize(&p); */
