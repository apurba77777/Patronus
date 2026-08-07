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

#   --------------------------- Tasks   ------------

#   Generate a spatial map of noise
if (argus.mapnoise):      
    fitslist   = [ pars['OutDir']+pars['CubeDir']+fname for fname in pars['FitsNames'] ]
    for fitsname in fitslist:
        tfdata      = fits.getdata(fitsname+".fits", ext=0)
        psfdata     = fits.getdata(fitsname+"_psf.fits", ext=0)
        
        #cubedata    = np.transpose(np.nanmean(tfdata, axis=1), axes=(2,1,0))
        #nsmap       = noisemap (cubedata, argus.pipedir+"/spew/", pars=pars)
        #np.save(fitsname+'_noisemap.npy', nsmap)
        cubedata    = np.nanmean(tfdata, axis=1)
        cubepsf     = np.nanmean(psfdata, axis=1)
        nsmap       = np.load(fitsname+'_noisemap.npy')
        cleancube (cubedata, cubepsf, argus.pipedir+"/spew/", nsmap, pars=pars)
        




        


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

