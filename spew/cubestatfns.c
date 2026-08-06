#include <stdio.h>
#include <stdlib.h>	
#include <spewhead.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <gsl/gsl_sort_float.h>
#include <gsl/gsl_statistics_float.h>

/* ----------------------------------------------------------------------------------------------------

    Functions to calculate staistics of a (multi dimensional) data cube

                                                            Last updated: AB (6 August 2026)

-------------------------------------------------------------------------------------------------------*/


int calcubenoise (float *datac, int datadim, int *dimlens, float *noise, int thrds) {

//  Function to calcluate noise in a multi dimensional data cube

//      datac   = datacube
//      datadim = number of dimensions
//      dimlens = Array of dimension lengths
//                  The first two are spatial dimentions
//      noise   = Array to store the computed noise   

    int i,j,p;

    printf("\nComputing noise in a %d dimensional array (",datadim);
    for(i = 0; i < datadim; i++)
        printf(" %d ",dimlens[i]);
    printf(")\n");

    #pragma omp parallel for num_threads(thrds) private(j,p)
    for(i = 0; i < dimlens[0]; i++) {
        for(j = 0; j < dimlens[1]; j++) {
            
            double *work = (double *) malloc( dimlens[2]*sizeof(double));
            p   = dimlens[2] * ( i*dimlens[1] + j ) ;

            noise[ (i*dimlens[1]) + j ] = (float) gsl_stats_float_mad(datac + p, 1, dimlens[2], work) ;

            free(work);
        }
    }
    
    printf("\n   Returning noise map \n");

    return 0;
} 
//  -------------------------------------------------------------------------------------------------------









