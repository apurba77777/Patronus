import os,sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import pickle as pkl
from collections import namedtuple
import matplotlib.pyplot as plt
from specscripts.auxfns import *
from specscripts.stackfns import *

#   ---------------------------------------------------------------------------------------------------
#
#       Process and analyze stacked cubes 
#
#   ---------------------------------------------------------------------------------------------------

def stackubes(gsamp, getspecdir, getcubedir, pars=None):

    #   Stack spectral cubes

    domedian	= (pars['StackStat']=="median")
    if (domedian):
        print("Stacking for median cube ...")
    else:
        print("Stacking for mean cube ...")

    gdets   = gsamp.sampcat
    print('Total cubes to be stacked	=	%d'%(len(gdets)))		

    allnoiseful	= []
    allcubefull	= []
    allspecfull	= []
    ngal		= 0

    for i in range (0, len(gdets)):
        detid	= int(gdets[i,1])
        detz	= gdets[i,4]
        detpb	= 1.0
        dlgpc	= d_lum_gpc(detz)

        ngal	= ngal+1
        ispecful= np.loadtxt(getspecdir+f'spec_{detid}_{pars['SmKpc']}.txt')	
        ivelarr	= ispecful[:,0]
        ivelreg	= np.where(np.abs(ivelarr) < pars['HalfLen'])[0]
        ivelarr	= ivelarr[ivelreg]
        inoise	= ispecful[ivelreg,2]
        inoise	= inoise/detpb		        #	Flux density noise correcting for primary beam

        allnoiseful.append(inoise)					
	
        icube	= fits.open(getcubedir+f'xubcube_{detid}_{pars['SmKpc']}.fits')	
        idata	= icube[0].data	
        idata	= idata[ivelreg]	
        icarr	= idata[:, int(idata.shape[1]/2) - pars['BSize'] : int(idata.shape[1]/2) + pars['BSize'] , \
                      int(idata.shape[2]/2) - pars['BSize'] : int(idata.shape[2]/2) + pars['BSize'] ]

        icarr	= icarr/detpb									    #	Flux density correcting for primary beam
        icarr	= icarr*(4*np.pi*dlgpc*dlgpc)/(1.0+detz)		    #	Luminosity density
        icarr0	= icarr[:, pars['BSize'], pars['BSize']]

        allcubefull.append(icarr)
        allspecfull.append(icarr0)
        icube.close()
    
    allnoiseful	= np.array(allnoiseful)
    allcubefull	= np.array(allcubefull)
    allspecfull	= np.array(allspecfull)    
        
    #print(allnoiseful.shape, allcubefull.shape)	
    print(f"Cubes stacked	= {ngal}")	

    velarr,stackarr,velavg,stackavg	=	stackcube_lumz_fluxwt(allnoiseful, pars['VelRes'], allcubefull, \
                                                        pars['RmsPow'], domedian)	

    stackarr	= subase(stackarr, pars['FitOrder'], pars['ExclChans'], pars['VelRes'])
    stackavg	= subase(stackavg, pars['FitOrder'], pars['ExAvgChans'], 2*pars['VelRes'])
    rmsarr		= planerms(stackarr, pars['ExRad'], pars['ExEdge'])
    rmsavg		= planerms(stackavg, pars['ExRad'], pars['ExEdge'])

    randrms,randrmsavg  = specrms_lum_fluxwt_random_fullspec(allspecfull, allnoiseful, pars['Realn'], pars['RmsPow'])
    
    velarr0, stackarr0, rmsnoisenormal, velarravg0, stackarravg0, avgrmsnoisenormal	= \
        stackspecs_int_lumz_fluxwt_singpix(allnoiseful, pars['VelRes'], allspecfull, pars['FitOrder'], \
                                           pars['ExclChans'], pars['RmsPow'], domedian)









    if (domedian):
        print("Saving median cube...")
        gsamp	= gsamp._replace(medstkcube = stackarr, medstkcubeavg = stackavg, medplnrmsarr = rmsarr, \
                               medplnrmsavg = rmsavg, medrandarr = randrms, medrandarravg = randrmsavg)
    else:
        print("Saving mean cube...")
        gsamp	= gsamp._replace(meanstkcube = stackarr, meanstkcubeavg = stackavg,	meanplnrmsarr = rmsarr, \
                               meanplnrmsavg = rmsavg, meanrandarr = randrms, meanrandarravg = randrmsavg)

    return (gsamp)
#   ---------------------------------------------------------------------------------------------------


