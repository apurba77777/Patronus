#include <stdio.h>
#include <stdlib.h>	
#include <spewhead.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <gsl/gsl_sort_float.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>

/* ----------------------------------------------------------------------------------------------------

    Functions to seacrh outliers data cube

                                                            Last updated: AB (6 August 2026)

-------------------------------------------------------------------------------------------------------*/


int cubecln (float *datac, float *datap, int datadim, int *dimlens, float *noise, int thrds, float sigthresh, float restbeam, int spikemax) {

//  Function to clean cubes in place 

//      datac       = datacube
//      datap       = PSF cube
//      datadim     = number of dimensions
//      dimlens     = Array of dimension lengths (Time first)
//      noise       = Noise map
//      thrds       = Number of threads
//      sigthresh   = Threshold in unit of noise RMS 
//      restbeam    = FWHM of restoring beam in pixels
//      spikemax    = Maximum spikes to clean per time   

    gsl_set_error_handler_off();

    int     t, i, j, i0, j0, l, m, n, p, x, y, spi, kmax, qmax;
    float   maxval,psfmax,drr;

    printf("\nCleaning a %d dimensional array (",datadim);
    for(i = 0; i < datadim; i++)
        printf(" %d ",dimlens[i]);
    printf(")\n\n");

    #pragma omp parallel for num_threads(thrds) private(i,j,i0,j0,l,m,n,p,x,y,spi,kmax,qmax,maxval,psfmax,drr)
    for(t = 0; t < dimlens[0]; t++) {

        int     spikex[spikemax], spikey[spikemax] ;
        float   sflux[spikemax] ;

        spi     = 0;

        p       = dimlens[2] * dimlens[1] * t ;
        maxval  = gsl_stats_float_max(datac + p, 1, dimlens[2] * dimlens[1]);
        kmax    = gsl_stats_float_max_index(datac + p, 1, dimlens[2] * dimlens[1]);
        i0      = kmax / dimlens[2] ;
        j0      = kmax % dimlens[2] ;

        /*.......................... Clean spikes above threshold .......................*/

        while (maxval > noise[kmax]*sigthresh*0.8) {
            printf("Cleaning at %d  %f  %d (%d %d) \n",t,(maxval/noise[kmax]),kmax,i0,j0);
            
            psfmax  = gsl_stats_float_max(datap + p, 1, dimlens[2] * dimlens[1]);
            qmax    = gsl_stats_float_max_index(datap + p, 1, dimlens[2] * dimlens[1]);
            l       = qmax / dimlens[2] ;
            m       = qmax % dimlens[2] ;
            
            spikex[spi] = i0;
            spikey[spi] = j0;
            sflux[spi]  = maxval;
            spi++;

            if ( spi >= spikemax ) {
                printf("Triggers exceed limit (%d) at t = %d \n",spikemax,t);
                break;
            } 

            for(i = MAX(0, (i0-l)); i < MIN((i0-l+dimlens[1]), dimlens[1]); i++) {
                for(j = MAX(0, (j0-m)); j < MIN((j0-m+dimlens[2]), dimlens[2]); j++) {

                    x   = MAX(0, (l-i0)) + (i - MAX(0, (i0-l))) ; 
                    y   = MAX(0, (m-j0)) + (j - MAX(0, (j0-m))) ;
                    
                    datac[p + i*dimlens[2] + j] = datac[p + i*dimlens[2] + j] 
                                    - (maxval/psfmax) * datap[p + x*dimlens[2] + y];
                }
            }   
            
            maxval  = gsl_stats_float_max(datac + p, 1, dimlens[2] * dimlens[1]);
            kmax    = gsl_stats_float_max_index(datac + p, 1, dimlens[2] * dimlens[1]);
            i0      = kmax / dimlens[2] ;
            j0      = kmax % dimlens[2] ;
            //printf("Now max at %d  %f  %d (%d %d) \n",t,(maxval/noise[kmax]),kmax,i0,j0);
        }   

        /*------------------- Restore spikes as circular Gaussians ----------------------*/        
        
        for ( n = 0; n < spi; n++) {

            i0      = spikex[n] ;
            j0      = spikey[n] ;
            maxval  = sflux[n] ;

            for(i = MAX(0, (i0-l)); i < MIN((i0-l+dimlens[1]), dimlens[1]); i++) {
                for(j = MAX(0, (j0-m)); j < MIN((j0-m+dimlens[2]), dimlens[2]); j++) {

                    x   = MAX(0, (l-i0)) + (i - MAX(0, (i0-l))) ; 
                    y   = MAX(0, (m-j0)) + (j - MAX(0, (j0-m))) ;

                    drr = ((float) (i - i0))*((float) (i - i0)) +
                            ((float) (j - j0))*((float) (j - j0)) ; 
                    
                    datac[p + i*dimlens[2] + j] = datac[p + i*dimlens[2] + j] 
                                    + (float) maxval*gsl_sf_exp( - 4 * M_LN2 * drr / (restbeam*restbeam));
                }
            } 
        }
    }
    
    printf("\n   Cleaning done!!! \n");

    return 0;
} 
//  -------------------------------------------------------------------------------------------------------



int subgaussian (float *imarr, float maxval, float restbeam, int i0, int j0, int dlen1, int dlen2) {

    //  Funtction to subtract a Gaussian from an image 

    int     i, j;
    float   drr;

    for(i = 0; i < dlen1; i++) {
        for(j = 0; j < dlen2; j++) {

            drr = ((float) (i - i0))*((float) (i - i0)) + ((float) (j - j0))*((float) (j - j0)) ; 
        
            imarr[i*dlen2 + j] = imarr[i*dlen2 + j] - (float) maxval*gsl_sf_exp( - 4 * M_LN2 * drr / (restbeam*restbeam));
        }
    }

    return(0);
}
//  -------------------------------------------------------------------------------------------------------



float calclocalnoise (float *imarr, int locnoise, int i0, int j0, int dlen1, int dlen2) {

    //  Calculate local noise in a plane

    int     i, l, m, iloc, rsize;
    float   rmsloc;

    float *locimg = (float *) malloc( locnoise*locnoise*sizeof(float));            

    l   = locnoise / 2 ;
    m   = locnoise / 2 ;
    iloc= 0 ;

    rsize   = MIN((j0-m+locnoise), dlen2) -  MAX(0, (j0-m)) ;

    for (i = MAX(0, (i0-l)); i < MIN((i0-l+locnoise), dlen1); i++) {
            
        memcpy(locimg + iloc*rsize, imarr + i*dlen2 + MAX(0, (j0-m)), rsize*sizeof(float)) ;
        iloc++;
    }

    double *work = (double *) malloc( iloc*rsize*sizeof(double));
    rmsloc       = (float) gsl_stats_float_mad(locimg, 1, iloc*rsize, work) ;
    
    free(locimg) ;
    free(work) ;

    return(rmsloc);
}
//  -------------------------------------------------------------------------------------------------------



int srchspike (float *datac, float *spikes, int datadim, int *dimlens, float *noise, int thrds, float sigthresh, float imgthresh, float restbeam, int spikemax, int locnoise) {

    //  Function to search for spikes in a cube

    //  datac       = datacube
    //  spikes      = Array of spikes
    //  datadim     = number of dimensions
    //  dimlens     = Array of dimension lengths (Time first)
    //  noise       = Noise map
    //  thrds       = Number of threads
    //  sigthresh   = Threshold in unit of noise RMS 
    //  restbeam    = FWHM of restoring beam in pixels
    //  spikemax    = Maximum spikes to clean per time   

    int     t, i, j, i0, j0, l, m, p, spi, kmax, qmax;
    float   maxval,drr;

    printf("\nSearching in a %d dimensional array (",datadim);
    for(i = 0; i < datadim; i++)
        printf(" %d ",dimlens[i]);
    printf(")\n\n");

    #pragma omp parallel for num_threads(thrds) private(i,j,i0,j0,l,m,p,spi,kmax,qmax,maxval,drr)
    for(t = 0; t < dimlens[0]; t++) {

        p       = dimlens[2] * dimlens[1] * t ;
        float *plnim = (float *) malloc(dimlens[1]*dimlens[2]*sizeof(float));

        spi     = 0 ;
        memcpy( plnim, datac + p, dimlens[2]*dimlens[1]*sizeof(float));

        maxval  = gsl_stats_float_max( (float *) plnim, 1, dimlens[2]*dimlens[1]);
        kmax    = gsl_stats_float_max_index( (float *) plnim, 1, dimlens[2]*dimlens[1]);
        i0      = kmax / dimlens[2] ;
        j0      = kmax % dimlens[2] ;

        /*.......................... Search for spikes above threshold .......................*/

        while ( maxval > noise[kmax]*sigthresh ) { 
            float   rmsloc;
            rmsloc  = calclocalnoise (plnim, locnoise, i0, j0, dimlens[1], dimlens[2]);

            if ( maxval > rmsloc*imgthresh ) {
                printf("Spike at %d  %f/%f  %d (%d %d) \n",t,(maxval/noise[kmax]),(maxval/rmsloc),kmax,i0,j0);

                spikes[t*spikemax*6 + spi*6]    = (float) t ;
                spikes[t*spikemax*6 + spi*6 + 1]= (float) i0 ;
                spikes[t*spikemax*6 + spi*6 + 2]= (float) j0 ;
                spikes[t*spikemax*6 + spi*6 + 3]= maxval ;
                spikes[t*spikemax*6 + spi*6 + 4]= maxval / noise[kmax] ;
                spikes[t*spikemax*6 + spi*6 + 5]= maxval / rmsloc ;

                spi++;

                if ( spi >= spikemax ) {
                    printf("Triggers exceed limit (%d) at t = %d \n",spikemax,t);
                    break;
                }
            }

            subgaussian (plnim, maxval, restbeam, i0, j0, dimlens[1], dimlens[2]) ;

            maxval  = gsl_stats_float_max( (float *) plnim, 1, dimlens[2]*dimlens[1]);
            kmax    = gsl_stats_float_max_index( (float *) plnim, 1, dimlens[2]*dimlens[1]);
            i0      = kmax / dimlens[2] ;
            j0      = kmax % dimlens[2] ;
            //printf("Now max at %d  %f  %d (%d %d) \n",t,(maxval/noise[kmax]),kmax,i0,j0);        
        }

        free(plnim);
    }

    return 0;
} 
//  -------------------------------------------------------------------------------------------------------









