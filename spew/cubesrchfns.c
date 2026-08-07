#include <stdio.h>
#include <stdlib.h>	
#include <spewhead.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <gsl/gsl_sort_float.h>
#include <gsl/gsl_statistics_float.h>

/* ----------------------------------------------------------------------------------------------------

    Functions to seacrh outliers data cube

                                                            Last updated: AB (6 August 2026)

-------------------------------------------------------------------------------------------------------*/


int cubecln (float *datac, float *datap, int datadim, int *dimlens, float *noise, int thrds, float sigthresh) {

//  Function to clean cubes in place ******** under development **************

//      datac   = datacube
//      datadim = number of dimensions
//      dimlens = Array of dimension lengths (Time first)
//      noise   = Noise map  

    int     t, i, j, i0, j0, l, m, n, p, x, y, kmax, qmax;
    float   maxval;

    printf("\nCleaning a %d dimensional array (",datadim);
    for(i = 0; i < datadim; i++)
        printf(" %d ",dimlens[i]);
    printf(")\n\n");

    //#pragma omp parallel for num_threads(thrds) private(j,p)
    for(t = 0; t < dimlens[0]; t++) {

        p       = dimlens[2] * dimlens[1] * t ;
        maxval  = gsl_stats_float_max(datac + p, 1, dimlens[2] * dimlens[1]);
        kmax    = gsl_stats_float_max_index(datac + p, 1, dimlens[2] * dimlens[1]);
        i0      = kmax / dimlens[2] ;
        j0      = kmax % dimlens[2] ;

        if (maxval > noise[kmax]*sigthresh) {
            printf("Cleaning at %d  %f  %d (%d %d) \n",t,(maxval/noise[kmax]),kmax,i0,j0);
            
            qmax= gsl_stats_float_max_index(datap + p, 1, dimlens[2] * dimlens[1]);
            l   = qmax / dimlens[2] ;
            m   = qmax % dimlens[2] ;

            printf("PSF max at %d %d \n",l,m);

            for(i = MAX(0, (i0-l)); i < MIN((i0-l+dimlens[1]), dimlens[1]); i++) {
                for(j = MAX(0, (j0-m)); j < MIN((j0-m+dimlens[2]), dimlens[2]); j++) {

                    x   = MAX(0, (l-i0)) ;
                    y   = MAX(0, (m-j0)) ;
                    
                    datac[i*dimlens[2] + j] = datac[i*dimlens[2] + j] - datap[x*dimlens[2] + y];
                }
            }      
        }    
    }
    
    printf("\n   Returning noise map \n");

    return 0;
} 
//  -------------------------------------------------------------------------------------------------------









