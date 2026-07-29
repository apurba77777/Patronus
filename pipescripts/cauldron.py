import os,sys
import argparse as ap
import yaml as ym
import pandas as pd
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
#                           samplot         //  Plot sample statistics
#                           stackcubes      //  Stack spectral cubes
#                           stackerrs       //  Calculate errors
#                           stackmass       //  Calculate average mass
#
#   Simple & convenient charms          
#                           obliviate       //  Clear existing files 
#                           lumos           //  List Usable Modes On Screen
#                           veritaserum     //  Reveal configuration parameters
#
#   Advanced spells and charms (Should NOT be attempted before passing O. W. L.s)
#                        
#                           !accio          //  Accumulate Continuum Components in Image Output
#                           !scourgify      //  Scrutinize Calibration Outputs and Ultimate Robustness of Gains with Image Files Yielded
#                           !incendio       //  Image Normal Continuum Emission using Nice Data from Interferometric Observations
#  
#   Dangerous spells and curses (Extreme caution recommended!! Should NOT be attempted before passing N.E.W.T.s)
#
#                           polyjuice       //
#                           amortentia      //  
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
    np.savetxt(pars['CatDir']+'/'+pars['CatFile'], srcs, \
               fmt = '%d  %d  %f  %f  %f  %d  %f  %d  %d  %d  %f  %f  %f  %f  %d  %f  %f  %f  %f  %f')

#	i	id	ra	dec	z	zq  d   px	py	ochan   SFR	lsm	NUV-r r-J  Type  U-V  V-J   SFR10  SFR100  lgm
#	0	1	2	3	4	5	6	7	8	9	    10	11	12    13   14    15   16    17     18      19

#   Extract subcubes
if (argus.exubcubes):  
    dets    = np.loadtxt(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['CatFile'])
    (detid, detz, posra, posdec, obschan) = (dets[:,1], dets[:,4], dets[:,2], dets[:,3], dets[:,9]) 
    detid   = detid.astype(int)
    exsrcs  = exubcubes(detid, detz, posra, posdec, obschan, istart=0, pars=pars, ovrt=argus.obliviate, nobj=-1)
    dets    = dets[exsrcs]
    print('Extracted sources		= %d'%len(dets))
    np.savetxt(pars['CatDir']+'/'+pars['CubeCatFile'], dets, \
                fmt = '%d  %d  %f  %f  %f  %d  %f  %d  %d  %d  %f  %f  %f  %f  %d  %f  %f  %f  %f  %f')


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
    intercubes(detid, detz, posra, posdec, obschan, istart=argus.istart, pars=pars, ovrt=argus.obliviate, nobj=-1)


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

    reducecubes(getcubedir, getnoisedir, redscubedir, outspecdir, detid, istart=argus.istart, pars=pars, cskip=coskip, nobj=-1)
    

#   Construct master catalogue
if (argus.getmastercat):  
    dets    = np.loadtxt(f"{pars['CatDir']}/{pars['CubeCatFile']}")
    gdets   = conmastercat(dets, pars=pars)
    print('\n     Sources	= %d'%len(gdets))
    np.savetxt(f"{pars['WorkDir']}/{pars['StacatDir']}/mastercat_{pars['StackName']}_{pars['ColType']}.cat", gdets, \
                fmt = '%d  %d  %f  %f  %f  %d  %f  %d  %d  %d  %f  %f  %f  %f  %d  %f  %f  %f  %f  %f')
    

#   Construct catalogue with well-behaved spectra
if (argus.getfinecat):  
    dets    = np.loadtxt(f"{pars['WorkDir']}/{pars['StacatDir']}/mastercat_{pars['StackName']}_{pars['ColType']}.cat")
    if (pars['ReGrid']):
        print("\n Regridding has not yet been implemented...\n")
        sys.exit()
    else:
        redscubedir = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedCubes']
        outspecdir  = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedSpecs']

    gdets   = confinecat(redscubedir, outspecdir, dets, pars=pars)
    print('\n     Sources	= %d'%len(gdets))
    np.savetxt(f"{pars['WorkDir']}/{pars['StacatDir']}/finecat_{pars['StackName']}_{pars['ColType']}.cat", gdets, \
                fmt = '%d  %d  %f  %f  %f  %d  %f  %d  %d  %d  %f  %f  %f  %f  %d  %f  %f  %f  %f  %f')


#   Construct sample for stacking
if (argus.stackcats or argus.amortentia):  
    refcat  = np.loadtxt(f"{pars['WorkDir']}/{pars['StacatDir']}/mastercat_{pars['StackName']}_{pars['NbrType']}.cat")
    detcat  = np.loadtxt(f"{pars['WorkDir']}/{pars['StacatDir']}/finecat_{pars['StackName']}_{pars['ColType']}.cat")
    constackcat(detcat, refcat, pars=pars)


#   Plot sample statistica
if (argus.samplot):
    if ((pars['BinEdges']==None) or (pars['BinParam']==None) or (pars['BinParam']=="")):
        samfilename = pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+".pkl"
        with open(samfilename, 'rb') as sampfile:	
            gsamp	= pkl.load(sampfile)
        zlsmplt(gsamp.sampcat, pars=pars)
        mainseqplt(gsamp.sampcat, pars=pars)
        zdist([gsamp.sampcat], f'{pars['StackName']}_{pars['ColType']}', pars=pars)       
        
    else:
        detcats = []
        for i in range (0, len(pars['BinEdges'])-1):
            samfilename = pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+"_"+pars["BinParam"]+"bin_"+str(i)+".pkl"
            with open(samfilename, 'rb') as sampfile:	
                gsamp	= pkl.load(sampfile)
            detcats.append(gsamp.sampcat)
        zdist(detcats, f"{pars['StackName']}_{pars['ColType']}_{pars['BinParam']}", pars=pars)             


#   Stack spectral cubes and/or calculate errors
if (argus.stackcubes or argus.stackerrs): 
    if (pars['ReGrid']):
        print("\n Regridding has not yet been implemented...\n")
        sys.exit()
    else:
        redscubedir = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedCubes']
        outspecdir  = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedSpecs']
    
    if ((pars['BinEdges']==None) or (pars['BinParam']==None) or (pars['BinParam']=="")):
        samfilename = pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+".pkl"
        with open(samfilename, 'rb') as sampfile:	
            gsamp	= pkl.load(sampfile)
        
        if (argus.stackcubes):
            mgsamp  = stackubes(gsamp, outspecdir, redscubedir, pars=pars)
            with open(samfilename, 'wb') as sampfile:	
                pkl.dump(mgsamp, sampfile)
            
        if (argus.stackerrs):
            mgsamp  = stackerrors(gsamp, outspecdir, pars=pars)
            with open(samfilename, 'wb') as sampfile:	
                pkl.dump(mgsamp, sampfile)
        
    else:
        for i in range (0, len(pars['BinEdges'])-1):
            samfilename = pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+"_"+pars["BinParam"]+"bin_"+str(i)+".pkl"
            with open(samfilename, 'rb') as sampfile:	
                gsamp	= pkl.load(sampfile)
            
            if (argus.stackcubes):
                mgsamp  = stackubes(gsamp, outspecdir, redscubedir, pars=pars)
                with open(samfilename, 'wb') as sampfile:
                    pkl.dump(mgsamp, sampfile)
                
            if (argus.stackerrs):
                mgsamp  = stackerrors(gsamp, outspecdir, pars=pars)
                with open(samfilename, 'wb') as sampfile:
                    pkl.dump(mgsamp, sampfile)
            
    
#   Calculate average mass
if (argus.stackmass):   
    if ((pars['BinEdges']==None) or (pars['BinParam']==None) or (pars['BinParam']=="")):
        samfilename = pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+".pkl"
        with open(samfilename, 'rb') as sampfile:	
            gsamp	= pkl.load(sampfile)
        
        mgsamp      = stackmasses(gsamp, pars=pars)
        results     = resultable(mgsamp, pars=pars)
        with open(samfilename, 'wb') as sampfile:
            pkl.dump(mgsamp, sampfile)
        
        results.to_csv(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+".txt")
    else:
        reslist = []
        for i in range (0, len(pars['BinEdges'])-1):
            samfilename = pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+"_"+pars["BinParam"]+"bin_"+str(i)+".pkl"
            with open(samfilename, 'rb') as sampfile:	
                gsamp	= pkl.load(sampfile)

            mgsamp      = stackmasses(gsamp, pars=pars)
            reslist.append(resultable(mgsamp, pars=pars))
            with open(samfilename, 'wb') as sampfile:	
                pkl.dump(mgsamp, sampfile)
            
        results = pd.concat(reslist)
        results.to_csv(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+"_"+pars["BinParam"]+".txt")
    
    
#   Stack spectral cubes, calculate errors and average mass
if (argus.amortentia): 
    if (pars['ReGrid']):
        print("\n Regridding has not yet been implemented...\n")
        sys.exit()
    else:
        redscubedir = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedCubes']
        outspecdir  = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['RedSpecs']
    
    if ((pars['BinEdges']==None) or (pars['BinParam']==None) or (pars['BinParam']=="")):
        samfilename = pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+".pkl"
        with open(samfilename, 'rb') as sampfile:	
            gsamp	= pkl.load(sampfile)
        
        mgsamp  = stackubes(gsamp, outspecdir, redscubedir, pars=pars)
        mgsamp  = stackerrors(mgsamp, outspecdir, pars=pars)
        mgsamp  = stackmasses(mgsamp, pars=pars)
        results = resultable(mgsamp, pars=pars)
        with open(samfilename, 'wb') as sampfile:	
            pkl.dump(mgsamp, sampfile)

        results.to_csv(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+".txt")
        
    else:
        reslist = []
        for i in range (0, len(pars['BinEdges'])-1):
            samfilename = pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+"_"+pars["BinParam"]+"bin_"+str(i)+".pkl"
            with open(samfilename, 'rb') as sampfile:	
                gsamp	= pkl.load(sampfile)
                        
            mgsamp  = stackubes(gsamp, outspecdir, redscubedir, pars=pars)       
            mgsamp  = stackerrors(mgsamp, outspecdir, pars=pars)
            mgsamp  = stackmasses(mgsamp, pars=pars)
            reslist.append(resultable(mgsamp, pars=pars))
            with open(samfilename, 'wb') as sampfile:
                pkl.dump(mgsamp, sampfile)  
        results = pd.concat(reslist)
        results.to_csv(pars['WorkDir']+'/'+pars['StacatDir']+'/'+pars['StackName']+"_"+pars['ColType']+"_"+pars["BinParam"]+".txt")  









    