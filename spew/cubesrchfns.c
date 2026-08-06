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

    int     t, i, j, p, kmax;
    float   maxval;

    printf("\nCleaning a %d dimensional array (",datadim);
    for(i = 0; i < datadim; i++)
        printf(" %d ",dimlens[i]);
    printf(")\n");

    //#pragma omp parallel for num_threads(thrds) private(j,p)
    for(t = 0; t < dimlens[0]; t++) {

        p   = dimlens[2] * dimlens[1] * t ;
        maxval  = gsl_stats_float_max(datac + p, 1, dimlens[2] * dimlens[1]);
        kmax    = gsl_stats_float_max_index(datac + p, 1, dimlens[2] * dimlens[1]);

        if (maxval > noise[kmax]*sigthresh) {
            printf("Cleaning at %d  %f  %d \n",t,(maxval/noise[kmax]),kmax);
            !!!!!!!!!!!!!!!! Now SUBTRACT !!!!!!!!!!!!!!!!!!
            **********  Or DO IN PYTHON *******************
        }    
    }
    
    printf("\n   Returning noise map \n");

    return 0;
} 
//  -------------------------------------------------------------------------------------------------------









