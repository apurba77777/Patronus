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
                      int(idata.shape[2]/2) + pars['BSize'] : int(idata.shape[2]/2) + pars['BSize'] ]

        icarr	= icarr/detpb									#	Flux density correcting for primary beam
        icarr	= icarr*(4*np.pi*dlgpc*dlgpc)/(1.0+detz)		#	Luminosity density
        
        allcubefull.append(icarr)
        icube.close()
    
    allnoiseful	=	np.array(allnoiseful)
    allcubefull	=	np.array(allcubefull)
        
    print(allnoiseful.shape, allcubefull.shape)	
    print(f"Cubes stacked	= {ngal}")	

    velarr,stackarr,velavg,stackavg	=	#stackspecs_int_lumz_fluxwt_cube (allnoiseful, velres, allcubefull, rmspow, domedian)	

    #stackarr	=	subase(stackarr,fitorder,exclchans,velres)
    #stackavg	=	subase(stackavg,fitorder,np.array(exclchans)/2,2*velres)

    #rmsarr		=	planerms(stackarr, exrad, exedge)
    #rmsavg		=	planerms(stackavg, exrad, exedge)

    return
#   ---------------------------------------------------------------------------------------------------