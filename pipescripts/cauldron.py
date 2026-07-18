import os,sys
import argparse as ap
import yaml as ym
from specscripts.ingredients import *
from specscripts.auxfns import *
from specscripts.compilecats import *
from specscripts.processubcubes import *
from specscripts.constructcats import *
from specscripts.spectralstack import *

#	---------------------------------------------------------------------------------------------------------
#
#	Cauldron to bre potions for spectral line analysis 
#                                                   AB  [last updated: 10 July 2026]
#
#   This programme can be used for radio spectroscopy 
#
#   To run this programme, use the following command with a python executor
#
#   cauldron.py         --[potion(s)]               //  processing step(s) -- see description below
#                       --infile [param_YAML]       //  YAML file containing input parameters
#                       --pipedir [pipe_direcory]   //  Path to the pipeline itself
#
#   Muggle-friendly potions (can be run without any prior knowledge of potion making)  
#             
#                           excat           //  Extract catalogue for spectroscopic analysis
#                           exubcubes       //  Extract subcubes
#                           smoothcubes     //  Spatially smooth subcubes
#                           intercubes      //  Spectrally interpolate subcubes
#                           regcubes        //  Regrid subcubes to spatial pixels
#                           redcubes        //  Reduce subcubes
#                           getmastercat    //  Construct master catalogue
#                           getfinecat      //  Construct catalogue with well-behaved spectra
#                           stackcats       //  Construct sample for stacking
#                           stackcubes      //  Stack spectral cubes
#                          
#
#   Simple & convenient charms          
#                           obliviate       //  Clear existing files 
#                           lumos           //  List Usable Modes On Screen
#                           revelio         //  Reveal configuration parameters
#
#   Advanced spells and charms (Should NOT be attempted before passing O. W. L.s)
#                        
#                           accio       //  Accumulate Continuum Components in Image Output
#                           scourgify   //  Scrutinize Calibration Outputs and Ultimate Robustness of Gains with Image Files Yielded
#                           incendio    //  Image Normal Continuum Emission using Nice Data from Interferometric Observations
#  
#   Dangerous spells and curses (Extreme caution recommended!! Should NOT be attempted before passing N.E.W.T.s)
#
#                           crucio      //  Calibrate Response for an Uncorrupted Channel Isolated from Observation 
#
#	--------------------------------------------------------------------------------------------------------


#   Get command line arguments
argus   = potion_args()

print("\n-------------------------------------------------------\n")
print("               Brewing potions                    ")
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
phials(pars)

#   List supported modes
if (argus.lumos):  
    potion_ingredients()

#   --------------------------- Potion brewing tasks   ----------------------------------

#   Extract COSMOS catalogue
if (argus.excat):  
    
    srcs    = exkhastocat(pars)
    if (pars['SrcType']=="galaxy"):
        srcs    = srcs[srcs[:,14]==0]
        print("Selecting galaxies only")
    elif (pars['SrcType']=="agn"):
        srcs    = srcs[srcs[:,14]==2]
        print("Selecting AGNs only")
    print('Extracted sources		= %d'%len(srcs))
    np.savetxt(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['CatFile'], srcs, \
               fmt = '%d  %d  %f  %f  %f  %d  %f  %d  %d  %d  %f  %f  %f  %f  %d')

#	i	id	ra	dec	z	zq  d   px	py	ochan   SFR	lsm	NUV R   Type
#	0	1	2	3	4	5	6	7	8	9	    10	11	12  13  14

#   Extract subcubes
if (argus.exubcubes):  
    dets    = np.loadtxt(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['CatFile'])
    (detid, detz, posra, posdec, obschan) = (dets[:,1], dets[:,4], dets[:,2], dets[:,3], dets[:,9]) 
    detid   = detid.astype(int)
    exsrcs  = exubcubes(detid, detz, posra, posdec, obschan, istart=0, pars=pars, ovrt=argus.obliviate, nobj=-1)
    dets    = dets[exsrcs]
    print('Extracted sources		= %d'%len(dets))
    np.savetxt(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['CubeCatFile'], dets, \
                fmt = '%d  %d  %f  %f  %f  %d  %f  %d  %d  %d  %f  %f  %f  %f  %d')


#   Spatially smooth subcubes
if (argus.smoothcubes):  
    if (pars['SmKpc'] > 0):
        print("\n Smoothing has not yet been implemented...\n")
    else:
        print("\n No smoothing required...\n")


#   Spectrally interpolate subcubes
if (argus.intercubes):  
    dets    = np.loadtxt(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['CubeCatFile'])
    (detid, detz, posra, posdec, obschan) = (dets[:,1], dets[:,4], dets[:,2], dets[:,3], dets[:,9]) 
    detid   = detid.astype(int)
    intercubes(detid, detz, posra, posdec, obschan, istart=0, pars=pars, ovrt=argus.obliviate, nobj=-1)


#   Regrid subcubes to spatial pixels
if (argus.regcubes):  
    if (pars['ReGrid']):        
        print("\n Regridding has not yet been implemented...\n")
    else:        
        print("\n No regridding required...\n")


#   Reduce subcubes and extract spectra
if (argus.redcubes):  
    dets    = np.loadtxt(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['CubeCatFile'])
    (detid, detz, posra, posdec, obschan) = (dets[:,1], dets[:,4], dets[:,2], dets[:,3], dets[:,9]) 
    detid   = detid.astype(int)
    if (pars['ReGrid']):
        print("\n Regridding has not yet been implemented...\n")
        sys.exit()
    else:
        getcubedir  = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['IscubeDir']
        getnoisedir = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['InoiseDir']
        redscubedir = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedCubes']
        outspecdir  = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedSpecs']
        coskip      = int(pars['BeamSec']/pars['BeamSec'])

    reducecubes(getcubedir, getnoisedir, redscubedir, outspecdir, detid, istart=0, pars=pars, cskip=coskip, nobj=-1)
    

#   Construct master catalogue
if (argus.getmastercat):  
    dets    = np.loadtxt(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['CubeCatFile'])
    gdets   = conmastercat(dets, pars=pars)
    print('\n     Sources	= %d'%len(gdets))
    np.savetxt(pars['WorkDir']+'/'+pars['StacatDir']+'/mastercat_'+pars['ColType']+'.cat', gdets, \
                fmt = '%d  %d  %f  %f  %f  %d  %f  %d  %d  %d  %f  %f  %f  %f  %d')
    

#   Construct catalogue with well-behaved spectra
if (argus.getfinecat):  
    dets    = np.loadtxt(pars['WorkDir']+'/'+pars['StacatDir']+'/mastercat_'+pars['ColType']+'.cat')
    if (pars['ReGrid']):
        print("\n Regridding has not yet been implemented...\n")
        sys.exit()
    else:
        redscubedir = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedCubes']
        outspecdir  = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedSpecs']

    gdets   = confinecat(redscubedir, outspecdir, dets, pars=pars)
    print('\n     Sources	= %d'%len(gdets))
    np.savetxt(pars['WorkDir']+'/'+pars['StacatDir']+'/finecat_'+pars['ColType']+'.cat', gdets, \
                fmt = '%d  %d  %f  %f  %f  %d  %f  %d  %d  %d  %f  %f  %f  %f  %d')


#   Construct sample for stacking
if (argus.stackcats):  
    refcat  = np.loadtxt(pars['WorkDir']+'/'+pars['StacatDir']+'/mastercat_'+pars['NbrType']+'.cat')
    detcat  = np.loadtxt(pars['WorkDir']+'/'+pars['StacatDir']+'/finecat_'+pars['ColType']+'.cat')

    constackcat(detcat, refcat, pars=pars)


#   Stack spectral cubes
if (argus.stackcubes): 
    if (pars['ReGrid']):
        print("\n Regridding has not yet been implemented...\n")
        sys.exit()
    else:
        redscubedir = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedCubes']
        outspecdir  = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedSpecs']
    
    if ((pars['BinEdges']==None) or (pars['BinParam']==None) or (pars['BinParam']=="")):
        sampfile    = open(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+".pkl", 'rb')	
        gsamp		= pkl.load(sampfile)
        sampfile.close()
        stackubes(gsamp, outspecdir, redscubedir, pars=pars)
    else:
        for i in range (0, len(pars['BinEdges'])-1):
            sampfile    = open(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars["BinParam"]+"bin_"+str(i)+".pkl", 'rb')	
            gsamp		= pkl.load(sampfile)
            sampfile.close()
            stackubes(gsamp, outspecdir, redscubedir, pars=pars)
    
    
    
    









    