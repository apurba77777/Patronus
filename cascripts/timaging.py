import os,sys
import numpy as np
import casatasks as ct
import casatools
from casaplotms import plotms
from astropy.io import fits



def maketime (visfile, pars=None):

    #   Make a list of time stamps in the given visibility dataset

    mjdfile = visfile+"_mjds.txt"
    
    msmd    = casatools.msmetadata()
    msmd.open(visfile+".ms")
    times   = msmd.timesforfield(0)
    msmd.done()    

    mjds    = times/86400

    mjds    = mjds[: : pars['TimeBins']]
    times   = times[: : pars['TimeBins']]

    print(f"\n Binning {pars['TimeBins']} timestamps...\n")

    print(f"Writing to {mjdfile}\n")    
    np.savetxt(mjdfile, np.array([times,mjds]).T, fmt="%.2f   %.6f")

    return(0)
#   -----------------------------------------------------------------------------------------------------



def makehdr (fitshdr, imdhr, tmjds):

    #   Build a FITS header from image header


    dtsec   = np.nanmedian(tmjds[1:,0] - tmjds[:-1,0])
    print(f"Time resolution = {dtsec} seconds")


    fitshdr['TIMESYS']  = imdhr['TIMESYS']
    fitshdr['RADESYS']  = imdhr['RADESYS']

    fitshdr['CTYPE1']   = imdhr['CTYPE1']
    fitshdr['CTYPE2']   = imdhr['CTYPE2']
    fitshdr['CTYPE3']   = imdhr['CTYPE3']
    fitshdr['CTYPE4']   = ('TIME', 'The time values are in an extension table')

    fitshdr['CUNIT1']   = imdhr['CUNIT1']
    fitshdr['CUNIT2']   = imdhr['CUNIT2']
    fitshdr['CUNIT3']   = imdhr['CUNIT3']
    fitshdr['CUNIT4']   = 'MJDSEC'

    fitshdr['CRPIX1']   = imdhr['CRPIX1']
    fitshdr['CRPIX2']   = imdhr['CRPIX2']
    fitshdr['CRPIX3']   = imdhr['CRPIX3']
    fitshdr['CRPIX4']   = 1.0

    fitshdr['CRVAL1']   = imdhr['CRVAL1']
    fitshdr['CRVAL2']   = imdhr['CRVAL2']
    fitshdr['CRVAL3']   = imdhr['CRVAL3']
    fitshdr['CRVAL4']   = 0.0

    fitshdr['CDELT1']   = imdhr['CDELT1']
    fitshdr['CDELT2']   = imdhr['CDELT2']
    fitshdr['CDELT3']   = imdhr['CDELT3']
    fitshdr['CDELT4']   = dtsec

    return(0)
#   -----------------------------------------------------------------------------------------------------



def timager (visfile, tmjds, pars=None, ntime=-1):

    #   Make a time series of image cubes 

    casq    = casatools.quanta()

    dtsec   = np.nanmedian(tmjds[1:,0] - tmjds[:-1,0])
    print(f"\n Time resolution = {dtsec} seconds\n")

    mjdlims = np.array([ (tmjds[:,0] - (dtsec/2.0)) / 86400 , (tmjds[:,0] + (dtsec/2.0)) / 86400 ]).T

    if (ntime >= 3):
        mjdlims = mjdlims[:ntime]

    print(f"\n Imaging {mjdlims.shape[0]} time intervals...\n")

    for ki in range(0,len(mjdlims)):

        os.system("rm -rf tcube_"+str(ki)+"*")

        mjdlim  = mjdlims[ki]
        
        startm  = casq.time(casq.quantity(mjdlim[0], 'd'),form=["ymd"]) 
        stoptm  = casq.time(casq.quantity(mjdlim[1], 'd'),form=["ymd"]) 

        tstring = startm[0]+"~"+stoptm[0]
        print(f"    Imaging time range {ki} / {len(mjdlims)}: {tstring}")

        ct.tclean(
            vis=visfile+".ms", \
            imagename="tcube_"+str(ki), \
            datacolumn="corrected", \
            imsize=pars['TimgSize'], \
            cell=pars['TcellSize'], \
            phasecenter=pars['TphaseCen'], \
            restoration=True, \
            selectdata=True, \
            field=pars['TargetName'],\
            uvrange=pars['FinUvLim'], \
            timerange=tstring, \
            specmode='mfs', \
            gridder='widefield', \
            wprojplanes=pars['TwprojPln'], \
            pblimit=0.2, \
            deconvolver='hogbom', \
            weighting='briggs',\
            robust=pars['TwtRobust'], \
            nsigma=pars['TinSigma'], \
            niter=pars['TimNiter']
        )
        

    casq.done()

    return(0)
#   -----------------------------------------------------------------------------------------------------



def tfimager (visfile, tmjds, pars=None, ntime=-1):

    #   Make a time series of image cubes 

    casq    = casatools.quanta()

    dtsec   = np.nanmedian(tmjds[1:,0] - tmjds[:-1,0])
    print(f"\n Time resolution = {dtsec} seconds\n")

    mjdlims = np.array([ (tmjds[:,0] - (dtsec/2.0)) / 86400 , (tmjds[:,0] + (dtsec/2.0)) / 86400 ]).T

    if (ntime >= 3):
        mjdlims = mjdlims[:ntime]

    print(f"\n Imaging {mjdlims.shape[0]} time intervals...\n")

    for ki in range(0,len(mjdlims)):

        os.system("rm -rf tfcube_"+str(ki)+"*")

        mjdlim  = mjdlims[ki]
        
        startm  = casq.time(casq.quantity(mjdlim[0], 'd'),form=["ymd"]) 
        stoptm  = casq.time(casq.quantity(mjdlim[1], 'd'),form=["ymd"]) 

        tstring = startm[0]+"~"+stoptm[0]
        print(f"     Imaging time range {ki} / {len(mjdlims)}: {tstring}")

        ct.tclean(
            vis=visfile+".ms", \
            imagename="tfcube_"+str(ki), \
            datacolumn="corrected", \
            imsize=pars['TimgSize'], \
            cell=pars['TcellSize'], \
            phasecenter=pars['TphaseCen'], \
            restoration=True, \
            selectdata=True, \
            field=pars['TargetName'],\
            uvrange=pars['FinUvLim'], \
            timerange=tstring, \
            specmode='cubedata', \
            width=pars['TavgChan'], \
            gridder='widefield', \
            wprojplanes=pars['TwprojPln'], \
            pblimit=0.2, \
            deconvolver='hogbom', \
            weighting='briggs',\
            robust=pars['TwtRobust'], \
            nsigma=pars['TinSigma'], \
            niter=pars['TimNiter']
        )
        

    casq.done()

    return(0)
#   -----------------------------------------------------------------------------------------------------



def makefits (tmjds, cubename, ntime=-1, nchan=-1):

    #   Combine image time series into a single FITS



    iman    = casatools.image()

    if (ntime >= 3):
        mjdarr  = tmjds[:ntime]
    else:
        mjdarr  = tmjds

    if (nchan < 1):
        ct.exportfits(imagename="tcube_0.image", fitsimage="tcube_0.fits", dropstokes=True, overwrite=True)
        refcube = fits.open("tcube_0.fits")
    else:
        ct.exportfits(imagename="tfcube_0.image", fitsimage="tfcube_0.fits", dropstokes=True, overwrite=True)
        refcube = fits.open("tfcube_0.fits")
        
    imhdr   = refcube[0].header
    refcube.close() 

    tfcube  = np.zeros((len(mjdarr), imhdr['NAXIS3'], imhdr['NAXIS2'], imhdr['NAXIS1']), dtype='float32')
    pfcube  = np.zeros((len(mjdarr), imhdr['NAXIS3'], imhdr['NAXIS2'], imhdr['NAXIS1']), dtype='float32')

    print(f"Reading {len(mjdarr)} images into cube {cubename}.fits...\n")

    for ki in range(0,len(mjdarr)):
        
        if (nchan < 1):
            iman.open("tcube_"+str(ki)+".image")
            timecube    = np.array([iman.getchunk(dropdeg=True)])
            #print(timecube.shape)
            tfcube[ki] = np.transpose(timecube, (0,2,1))
            iman.done()

            iman.open("tcube_"+str(ki)+".psf")
            psfcube    = np.array([iman.getchunk(dropdeg=True)])
            #print(psfcube.shape)
            pfcube[ki] = np.transpose(psfcube, (0,2,1))
            iman.done()

        else:
            iman.open("tfcube_"+str(ki)+".image")
            timecube    = iman.getchunk(dropdeg=True)
            #print(timecube.shape)
            tfcube[ki] = np.transpose(timecube, (2,1,0))
            iman.done()     

            iman.open("tfcube_"+str(ki)+".psf")
            psfcube    = iman.getchunk(dropdeg=True)
            #print(psfcube.shape)
            pfcube[ki] = np.transpose(psfcube, (2,1,0))
            iman.done()   

    iman.close()
    
    print(f"\n  Cube dimensions = {tfcube.shape} \n")

    #   Image 4d cube
    tfhdu       = fits.PrimaryHDU()
    tfhdu.data  = tfcube    
    makehdr(tfhdu.header, imhdr, tmjds)

    hdulist     = fits.HDUList([tfhdu])

    #   Time table
    tcol    = fits.Column(name='MJDSEC', format='D', unit='MJDSEC', array=mjdarr[:,0])
    coldefs = fits.ColDefs([tcol])
    tabhdu  = fits.BinTableHDU.from_columns(coldefs)
    tabhdu.header['EXTNAME']   = "TIME"
    tabhdu.header['TIMESYS']  = imhdr['TIMESYS']

    hdulist.append(tabhdu)    
    hdulist.writeto(cubename+".fits", overwrite=True)

    #   PSF 4d cube
    pfhdu       = fits.PrimaryHDU()
    pfhdu.data  = pfcube    
    makehdr(pfhdu.header, imhdr, tmjds)
    phdulist    = fits.HDUList([pfhdu])
    phdulist.append(tabhdu)    
    phdulist.writeto(cubename+"_psf.fits", overwrite=True)

    print("\n    Clearing images...\n")
    for ki in range(0,len(mjdarr)):
        if (nchan < 1):
            os.system("rm -rf tcube_"+str(ki)+"*")
        else:
            os.system("rm -rf tfcube_"+str(ki)+"*")


    return(0)
#   -----------------------------------------------------------------------------------------------------
