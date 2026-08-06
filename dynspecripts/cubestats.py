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

def noisemap (fitsfile, spewdir, pars=None):

    #   Generate noise map for the cube

    spew        = ctp.CDLL( os.path.abspath(f"{spewdir}/cubestatfns.so") ) 
    #int calcubenoise (double *datac, int datadim, int *dimlens, double *noise)
    
    hdulist     = fits.open(fitsfile+".fits")
    tfdata      = hdulist[0].data
    hdulist.close()

    cubedata    = np.transpose(np.nanmean(tfdata, axis=1), axes=(2,1,0))

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

    print(retval)

    spatnoise   = np.reshape(spatnoise, (tdata.shape[0], tdata.shape[1]))
    
    fig     = plt.figure(figsize=(4,4))     
    ax5     = fig.add_subplot(111)

    plt.imshow(spatnoise, origin='lower', interpolation='none', aspect='auto', cmap='plasma')
    plt.tight_layout()
    plt.show()   

    return(0)
#   -----------------------------------------------------------------------------------------------------




