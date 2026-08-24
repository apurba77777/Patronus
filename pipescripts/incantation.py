import os,sys
import argparse as ap
import yaml as ym
from cascripts.utils import *
from dynspecripts.transearch import *
from dynspecripts.cubestats import *
from dynspecripts.cubesrch import *

#	---------------------------------------------------------------------------------------------------------
#
#	Incatations for spells and charms  
#                                                   AB  [last updated: 1 July 2026]
#
#   This programme can be used to identify transient events in a time-resolved image cube
#
#   To run this programme, use the following command with a python executor
#
#   incantation.py      --[spell(s)]                //  processing step(s) -- see description below
#                       --infile [param_YAML]       //  YAML file containing input parameters
#                       --pipedir [pipe_direcory]   //  Path to the pipeline itself
#
#   Muggle-friendly spells              
#                       getdspec    //  Generate dynamic spectrum at a specific sky position
#                       mapnoise    //  Generate spatial map of noise
#                       cleensweep  //  Clean and search for transients at the original time resolution
#                       acleensweep //  Clean and search for transients after time averaging
#
#   Simple & convenient charms          
#                       obliviate   //  Clear existing files 
#                       lumos       //  List Usable Modes On Screen
#                       revelio     //  Reveal configuration parameters
#
#   Advanced spells and charms (Should NOT be attempted before passing O. W. L.s)
#                        
#                       accio       //  Accumulate Continuum Components in Image Output
#                       scourgify   //  Scrutinize Calibration Outputs and Ultimate Robustness of Gains with Image Files Yielded
#                       incendio    //  Image Normal Continuum Emission using Nice Data from Interferometric Observations
#  
#   Dangerous spells and curses ( Extreme caution recommended !! Should NOT be attempted before passing N.E.W.T.s)
#
#                       crucio      //  Calibrate Response for an Uncorrupted Channel Isolated from Observation 
#
#	--------------------------------------------------------------------------------------------------------


#   Get command line arguments
argus   = incan_args()

print("\n-------------------------------------------------------\n")
print("        Concentrate And Search Transients (CAST) ")
print("\n-------------------------------------------------------\n")

#   Read and complain about input parameters
if (argus.infile == None):
    print(" Missing input YAML file! Please provide one...")
    sys.exit()
else:
    with open(argus.infile+'.yml', 'r') as infl:
        pars = ym.load(infl, Loader=ym.SafeLoader)
        if (argus.revelio):
            print(" Inputs provided -- \n")
            print(ym.dump(pars, sort_keys=False))

#   Make data directories
niffler(pars)

#   List supported modes
if (argus.lumos):  
    incan_spells()

#   --------------------------- Spells   ----------------------------------------------

fitslist   = [ pars['OutDir']+pars['CubeDir']+fname for fname in pars['FitsNames'] ]

#   Generate a spatial map of noise
if (argus.mapnoise):      
    for fitsname in fitslist:
        # tfdata      = fits.getdata(fitsname+".fits", ext=0)
        # cubedata    = np.transpose(np.nanmean(tfdata, axis=1), axes=(2,1,0))
        # nsmap       = noisemap (cubedata, argus.pipedir+"/spew/", pars=pars)
        # np.save(fitsname+'_noisemap.npy', nsmap)
        with fits.open(fitsname+".fits") as tfdhu:
            tfdata      = tfdhu[0].data
            spatpix     = (tfdata.shape[2], tfdata.shape[3])
            spatnoise   = np.zeros((tfdata.shape[3], tfdata.shape[2]), dtype='float32')
            subpix      = 1 + (np.array(spatpix) // pars['SubIms'])
            print(f"Spatial size {spatpix}")
            print(f"Dividing in {pars['SubIms']} sub-images of {subpix}")

            for sx in range (0, pars['SubIms']):
                xstart = sx * subpix[0]
                xstop  = min(xstart + subpix[0], spatpix[0])
                for sy in range (0, pars['SubIms']):
                    ystart = sy * subpix[1]
                    ystop  = min(ystart + subpix[1], spatpix[1])
                    print(xstart,xstop,ystart,ystop)
                    tfcrop      = tfdata[:, :, xstart:xstop, ystart:ystop]
                    cubedata    = np.transpose(np.nanmean(tfcrop, axis=1), axes=(2,1,0))
                    nsmap       = noisemap (cubedata, argus.pipedir+"/spew/", pars=pars)
                    spatnoise[ystart:ystop, xstart:xstop]   = nsmap
                    del nsmap
                    del cubedata

            np.save(fitsname+'_noisemap.npy', spatnoise)

            fig     = plt.figure(figsize=(5,4))     
            ax5     = fig.add_subplot(111)
        
            plt.imshow(spatnoise, origin='lower', interpolation='none', aspect='auto', cmap='plasma')
            plt.colorbar()
            plt.tight_layout()
            plt.show()



#   Clean and search for transients at the original time resolution
if (argus.cleansweep):      
    for fitsname in fitslist:    
        tfdata      = fits.getdata(fitsname+".fits", ext=0)
        psfdata     = fits.getdata(fitsname+"_psf.fits", ext=0)
        if (tfdata.shape[1] > 1):
            cubedata    = np.nanmean(tfdata, axis=1)
            cubepsf     = np.nanmean(psfdata, axis=1)
        else:
            cubedata    = tfdata[:,0]
            cubepsf     = psfdata[:,0]
        nsmap       = np.load(fitsname+'_noisemap.npy')
        print("Attempting to clean cube...")
        cleancube (cubedata, cubepsf, argus.pipedir+"/spew/", nsmap, pars=pars)
        searchcube (cubedata, argus.pipedir+"/spew/", nsmap, pars=pars)


#   Clean and search for transients after time averaging
if (argus.acleansweep):      
    for fitsname in fitslist:    
        hdulist     = fits.open(fitsname+".fits")
        tfdata      = hdulist[0].data
        mjdtime     = hdulist[1].data['MJDSEC']
        dtsec       = hdulist[0].header['CDELT4']
        hdulisp     = fits.open(fitsname+"_psf.fits")
        psfdata     = hdulisp[0].data
        if (tfdata.shape[1] > 1):
            cubedata    = np.nanmean(tfdata, axis=1)
            cubepsf     = np.nanmean(psfdata, axis=1)
        else:
            cubedata    = tfdata[:,0]
            cubepsf     = psfdata[:,0]
        nsmap       = np.load(fitsname+'_noisemap.npy')
        print("Attempting to clean cube...")
        hdulist.close()
        hdulisp.close()

        for tavg in pars['TavgFacs']:
            print(f"Time averaging by a factor of {tavg}")
            for tmov in [0, int(tavg/2)]:
                nsmap       = np.load(fitsname+'_noisemap.npy') / np.sqrt(float(tavg))
                avgcube, avgpsf, avgmjdsec  = timeavg (cubedata, psfdata, mjdtime, dtsec, tavgfac=tavg, tshift=tmov, pars=None)
                cleancube (avgcube, avgpsf, argus.pipedir+"/spew/", nsmap, pars=pars)
                searchcube (avgcube, argus.pipedir+"/spew/", nsmap, pars=pars)
        
    
#   Search at a specific sky position
if (argus.getdspec):      
    fitslist   = [ pars['OutDir']+pars['CubeDir']+fname for fname in pars['FitsNames'] ]
    for fitsname in fitslist:
        getdynspec (fitsname, pars=pars)


#   Show patronus
if (argus.expecto_patronum or argus.lumos or argus.revelio):
    patronus_charm (argus)


print("\n----------------------------------------------------------------------------")
print("----------------------------------------------------------------------------")

