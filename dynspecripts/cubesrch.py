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


def cleancube (cubedata, cubepsf, spewdir, nsmap, pars=None):

    #   Clean the data cube

    spew        = ctp.CDLL( os.path.abspath(f"{spewdir}/cubesrchfns.so") ) 
    #int cubecln (float *datac, float *datap, int datadim, int *dimlens, float *noise, 
    #                                    int thrds, float sigthresh, float restbeam, int spikemax)

    tdata       = np.ascontiguousarray(cubedata, dtype='float32') 
    pdata       = np.ascontiguousarray(cubepsf, dtype='float32')    
    spatnoise   = np.ascontiguousarray(nsmap, dtype='float32')
    datadims    = np.ascontiguousarray(tdata.shape, dtype='intc')

    print(f"\nCube dimensions {tdata.shape} type {tdata.dtype}")

    fig     = plt.figure(figsize=(5,4))     
    ax5     = fig.add_subplot(111)

    #plt.imshow(spatnoise, origin='lower', interpolation='none', aspect='auto', cmap='plasma')
    plt.imshow(tdata[201], origin='lower', interpolation='none', aspect='auto', cmap='plasma')
    plt.colorbar()
    plt.tight_layout()
    plt.show()  
    
    spew.cubecln.argtypes = [
        ctp.POINTER(ctp.c_float),
        ctp.POINTER(ctp.c_float),
        ctp.c_int,
        ctp.POINTER(ctp.c_int),
        ctp.POINTER(ctp.c_float),
        ctp.c_int,
        ctp.c_float,
        ctp.c_float,
        ctp.c_int
    ]

    spew.cubecln.restype = ctp.c_int

    tdataptr    = tdata.ctypes.data_as(ctp.POINTER(ctp.c_float))
    pdataptr    = pdata.ctypes.data_as(ctp.POINTER(ctp.c_float))
    noiseptr    = spatnoise.ctypes.data_as(ctp.POINTER(ctp.c_float))
    dimptr      = datadims.ctypes.data_as(ctp.POINTER(ctp.c_int))

    retval      = spew.cubecln(tdataptr, pdataptr, np.intc(tdata.ndim), dimptr, noiseptr, np.intc(pars['Threads']), \
                            np.single(pars['SigThresh']), np.single(pars['RestBeam']), np.intc(pars['MaxSrc']))
    
    fig     = plt.figure(figsize=(5,4))     
    ax5     = fig.add_subplot(111)

    #plt.imshow(spatnoise, origin='lower', interpolation='none', aspect='auto', cmap='plasma')
    plt.imshow(tdata[202], origin='lower', interpolation='none', aspect='auto', cmap='plasma')
    plt.colorbar()
    plt.tight_layout()
    plt.show()   

    return(0)
#   -----------------------------------------------------------------------------------------------------


def searchcube (cubedata, spewdir, nsmap, pars=None):

    #   Search for spikes in the data cube

    spew        = ctp.CDLL( os.path.abspath(f"{spewdir}/cubesrchfns.so") ) 
    #int srchspike (float *datac, float *spikes, int datadim, int *dimlens, float *noise, 
    #                                    int thrds, float sigthresh, float restbeam, int spikemax)

    tdata       = np.ascontiguousarray(cubedata, dtype='float32') 
    spikes      = np.ascontiguousarray(np.zeros( pars['MaxSrc']*cubedata.shape[0]*4, dtype='float32'))    
    spatnoise   = np.ascontiguousarray(nsmap, dtype='float32')
    datadims    = np.ascontiguousarray(tdata.shape, dtype='intc')

    print(f"\nCube dimensions {tdata.shape} type {tdata.dtype}")
    
    spew.srchspike.argtypes = [
        ctp.POINTER(ctp.c_float),
        ctp.POINTER(ctp.c_float),
        ctp.c_int,
        ctp.POINTER(ctp.c_int),
        ctp.POINTER(ctp.c_float),
        ctp.c_int,
        ctp.c_float,
        ctp.c_float,
        ctp.c_int
    ]

    spew.srchspike.restype = ctp.c_int

    tdataptr    = tdata.ctypes.data_as(ctp.POINTER(ctp.c_float))
    spikeptr    = spikes.ctypes.data_as(ctp.POINTER(ctp.c_float))
    noiseptr    = spatnoise.ctypes.data_as(ctp.POINTER(ctp.c_float))
    dimptr      = datadims.ctypes.data_as(ctp.POINTER(ctp.c_int))

    retval      = spew.srchspike(tdataptr, spikeptr, np.intc(tdata.ndim), dimptr, noiseptr, np.intc(pars['Threads']), \
                            np.single(pars['SigThresh']), np.single(pars['RestBeam']), np.intc(pars['MaxClnSrc']))
    
    # fig     = plt.figure(figsize=(5,4))     
    # ax5     = fig.add_subplot(111)

    # #plt.imshow(spatnoise, origin='lower', interpolation='none', aspect='auto', cmap='plasma')
    # plt.imshow(tdata[202], origin='lower', interpolation='none', aspect='auto', cmap='plasma')
    # plt.colorbar()
    # plt.tight_layout()
    # plt.show()   

    return(0)
#   -----------------------------------------------------------------------------------------------------




