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
    
    return(0)
#   -----------------------------------------------------------------------------------------------------


def searchcube (cubedata, spewdir, nsmap, pars=None):

    #   Search for spikes in the data cube

    spew        = ctp.CDLL( os.path.abspath(f"{spewdir}/cubesrchfns.so") ) 
    #int srchspike (float *datac, float *spikes, int datadim, int *dimlens, float *noise, 
    #                int thrds, float sigthresh, float imgthresh, float restbeam, int spikemax, int locnoise)

    tdata       = np.ascontiguousarray(cubedata, dtype='float32') 
    spikes      = np.ascontiguousarray(np.zeros( pars['MaxSrc']*cubedata.shape[0]*6, dtype='float32'))    
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
        ctp.c_float,
        ctp.c_int,
        ctp.c_int
    ]

    spew.srchspike.restype = ctp.c_int

    tdataptr    = tdata.ctypes.data_as(ctp.POINTER(ctp.c_float))
    spikeptr    = spikes.ctypes.data_as(ctp.POINTER(ctp.c_float))
    noiseptr    = spatnoise.ctypes.data_as(ctp.POINTER(ctp.c_float))
    dimptr      = datadims.ctypes.data_as(ctp.POINTER(ctp.c_int))

    retval      = spew.srchspike(tdataptr, spikeptr, np.intc(tdata.ndim), dimptr, noiseptr, np.intc(pars['Threads']), \
                                np.single(pars['SigThresh']), np.single(pars['ImgThresh']), np.single(pars['RestBeam']), \
                                    np.intc(pars['MaxSrc']), np.intc(pars['LocNoise']))
    

    spikes      = np.reshape(spikes,(cubedata.shape[0]*pars['MaxSrc'], 6))
    #   t   px  py  maxv    tsnr    psnr
    #   0   1   2   3       4       5

    spikes      = spikes[spikes[:,4] > pars['SigThresh']]
    print(f"\n  Found {len(spikes)} objects")

    for i in range(0, len(spikes)):
        fig     = plt.figure(figsize=(5,4))     
        ax5     = fig.add_subplot(111)
        
        plt.imshow(tdata[int(spikes[i,0])], origin='lower', interpolation='none', aspect='auto', cmap='plasma')
        plt.colorbar()
        plt.tight_layout()
        plt.show()   

    return(0)
#   -----------------------------------------------------------------------------------------------------




