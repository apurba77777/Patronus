import os,sys
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from scipy.stats import kstest
from scipy.stats import anderson
import matplotlib.pyplot as plt
from specscripts.auxfns import *
from specscripts.samplotfns import *

#   ---------------------------------------------------------------------------------------------------
#
#       Functions to construct final catalogues 
#
#   ---------------------------------------------------------------------------------------------------

def confinecat(redscubedir, outspecdir, gdets, pars=None):
    
    #   catalogue with well-behaved spectra	

    #	i	id	ra	dec	z	zq  d   px	py	ochan   SFR	lsm	NUV R   Type
    #	0	1	2	3	4	5	6	7	8	9	    10	11	12  13  14

    print('Sample size full		=	%d'%len(gdets))

    gdets00		=	np.copy(gdets)
    kindices	=	[]
    gindices	=	[]
    dindices	=	[]
    oindices	=	[]

    for i in range (0,len(gdets)):
		
        ispecraw	=	np.loadtxt(outspecdir+'/spec_'+str(int(gdets[i,1]))+'_'+str(pars['SmKpc'])+'.txt')	
        #	relvel	fluxrms_interpolated(Jy)	fluxrms_regridded(Jy)	flux(Jy)		

        ispecraw	=	ispecraw[np.abs(ispecraw[:,0]) < pars['HalfLen']]
        ispecfull	=	ispecraw[np.isfinite(ispecraw[:,2])]
        inoise		=	ispecfull[:,2]
        ispec		=	ispecfull[:,3]
        fullnoise	=	np.nanstd(ispec)

        basespecfull=	ispecfull[np.abs(ispecfull[:,0]) > pars['HalfWidth']]
        baseinoise	=	basespecfull[:,2]
        basespec	=	basespecfull[:,3]

        #	Maximum allowed flagging criterion
        if(len(ispec) >= pars['MincFrac']*len(ispecraw)):								
            kindices.append(i)
            ispecsnr	=	ispec/inoise

            #	Maximum allowed S/N criterion
            if((np.nanmax(np.abs(ispecsnr)) < pars['MaxOut']) and (np.nanmax(np.abs(ispec/fullnoise)) < pars['MaxOut'])):									
                oindices.append(i)
                basesnr		=	basespec/baseinoise
                snrnoise	=	np.nanstd(basesnr)
                ksd,ksp		=	kstest(basesnr,'norm',args=(0.0,snrnoise))				
                if(ksp > pars['KspLim']):					#	KS test criterion
                    gindices.append(i)					
                    ads,adcv,adsl	=	anderson(basesnr,dist='norm')	
                    
                    if(ads < pars['AdsLim']):				#	AD test criterion
                        dindices.append(i)

    print('Excluding >min flagged		= %d'%len(gdets00[kindices]))	
    print('Excluding n-sigma fetures	= %d'%len(gdets00[oindices]))			
    print('Excluding KS failed		= %d'%len(gdets00[gindices]))	
    print('Excluding AD failed		= %d'%len(gdets00[dindices]))	

    gdets00		= gdets00[dindices]		

    #	How does the noise change with velocity averaging ? 

    cenchan		=	int(2*pars['HalfLen'] / pars['VelRes'])
    nichan		=	pars['CavgFac'] * ( cenchan // pars['CavgFac'])
    nujy		=	[]
    gdetsout	=	[]
    gdets		=	gdets00

    for i in range (0,len(gdets)):
        #print("src	%d (%d / %d)"%(int(gdets[i,0]), i, len(gdets)))	
        detid	=	int(gdets[i,1])		
        icube	=	fits.open(redscubedir+'/xubcube_'+str(detid)+'_'+str(pars['SmKpc'])+'.fits')
        idata	=	icube[0].data
        avgcube	=	np.nanmean(np.reshape(idata[ np.arange( (idata.shape[0]-nichan)//2 , \
                            (idata.shape[0]+nichan)//2, 1, dtype=int)], (-1,pars['CavgFac'],idata.shape[1],idata.shape[2])), axis=1)

        irms	=	np.nanstd(idata) * 1.0e6
        avgrms	=	np.nanstd(avgcube) * 1.0e6
        nujy.append([irms, avgrms])
        icube.close()

        if (((avgrms/irms) - pars['RatLims'][0])*((avgrms/irms) - pars['RatLims'][1]) <= 0.0):
            gdetsout.append(gdets[i])

    nujy	= np.array(nujy)
    setrats	= nujy[:,1]/nujy[:,0]
    gdetsout= np.array(gdetsout)

    noiseratplt(setrats, pars=pars)

    return (gdetsout) 
#	--------------------------------------------------------------------------------------------------














