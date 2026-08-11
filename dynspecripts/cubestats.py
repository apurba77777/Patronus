import os,sys
import numpy as np
import ctypes as ctp
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

#   --------------------------------------------------------------------------------------------------------------------------
#   Calculate statistics of the data cube
#
#                                              AB  [last updated: 5 August 2026] 
#
#   --------------------------------------------------------------------------------------------------------------------------


def noisemap (cubedata, spewdir, pars=None):

    #   Generate noise map for the cube

    spew        = ctp.CDLL( os.path.abspath(f"{spewdir}/cubestatfns.so") ) 
    #int calcubenoise (double *datac, int datadim, int *dimlens, double *noise)

    tdata       = np.ascontiguousarray(cubedata, dtype='float32')    
    spatnoise   = np.ascontiguousarray(np.zeros((tdata.shape[0], tdata.shape[1]), dtype='float32'))
    datadims    = np.ascontiguousarray(tdata.shape, dtype='intc')

    print(f"\nCube dimensions {tdata.shape} type {tdata.dtype}")
    
    spew.calcubenoise.argtypes = [
        ctp.POINTER(ctp.c_float),
        ctp.c_int,
        ctp.POINTER(ctp.c_int),
        ctp.POINTER(ctp.c_float),
        ctp.c_int
    ]

    spew.calcubenoise.restype = ctp.c_int

    tdataptr    = tdata.ctypes.data_as(ctp.POINTER(ctp.c_float))
    noiseptr    = spatnoise.ctypes.data_as(ctp.POINTER(ctp.c_float))
    dimptr      = datadims.ctypes.data_as(ctp.POINTER(ctp.c_int))

    retval      = spew.calcubenoise(tdataptr, np.intc(tdata.ndim), dimptr, noiseptr, np.intc(pars['Threads']))

    spatnoise   = np.reshape(spatnoise, (tdata.shape[0], tdata.shape[1]))
    
    fig     = plt.figure(figsize=(5,4))     
    ax5     = fig.add_subplot(111)

    plt.imshow(spatnoise, origin='lower', interpolation='none', aspect='auto', cmap='plasma')
    plt.colorbar()
    plt.tight_layout()
    plt.show()   

    return(spatnoise)
#   -----------------------------------------------------------------------------------------------------



def timeavg (cubedata, psfdata, mjdsecs, dtsec, tavgfac=1, tshift=0, pars=None):

    #   Average datacubes in time

    if (tavgfac <= 1):
        print("No averaging required...")
        return(cubedata, mjdsecs)

    print(f"Original time resolution = {dtsec: .2f} seconds ")

    avgmjdsecs  = []
    bindices    = []
    tstart      = tshift

    while (tstart < len(mjdsecs)):
        mjd0        = mjdsecs[tstart] - 0.5*dtsec
        tindices    = [tstart]
        ti          = 1

        while ( (tstart + ti < len(mjdsecs)) and ((mjdsecs[tstart + ti] - mjd0) <= dtsec*tavgfac)):
            tindices.append(tstart + ti)
            ti      = ti + 1

        if (len(tindices) == tavgfac):
            bindices.append(tindices)
            avgmjdsecs.append(np.nanmedian(mjdsecs[tindices]))
        tstart      = tstart + ti

    avgmjdsecs  = np.array(avgmjdsecs)
    print(f"After averaging length = {len(avgmjdsecs)} \
                \n  resolution = {np.nanmedian(avgmjdsecs[1:] - avgmjdsecs[:-1]): .2f} sec")   

    avgcube     = np.zeros((len(avgmjdsecs), cubedata.shape[1], cubedata.shape[2]), dtype='float32')
    avgpsf      = np.zeros((len(avgmjdsecs), cubedata.shape[1], cubedata.shape[2]), dtype='float32')

    for i in range(0, len(avgmjdsecs)):
        avgcube[i]  = np.nanmean(cubedata[bindices[i]], axis=0)
        avgpsf[i]   = np.nanmean(psfdata[bindices[i]], axis=0)        

    return(avgcube, avgpsf, avgmjdsecs)
#   -----------------------------------------------------------------------------------------------------

