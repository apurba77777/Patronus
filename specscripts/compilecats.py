import os,sys
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from specscripts.auxfns import *

#   ---------------------------------------------------------------------------------------------------
#
#       Functions to compile catalogues for analysis from parent catalogues
#
#   ---------------------------------------------------------------------------------------------------

	
def exkhastocat(pars):
    
    #   Function to extract sources from COSMOS compilation catalogue by Khastobhan et al.
    
    reffile =   fits.open(pars['CubeDir']+'/'+pars['ImgCube']+".fits")    
    testhd	= 	reffile[0].header							#	image header
    w		=	WCS(testhd,naxis=2)

    if (pars['CenRa']==None):
        ra0	=	testhd['CRVAL1']							#	Central RA
    else:
        ra0	=	pars['CenRa']
    
    if (pars['CenDec']==None):
        dec0=	testhd['CRVAL2']							# 	Central DEC
    else:
        dec0=	pars['CenDec']

    refreq	=	testhd['CRVAL3']/1.0e6						#	Reference frequency in MHz
    refloc	=	testhd['CRPIX3'] - pars['PiXys']			#	Reference location in frequency axis
    dfreq	=	testhd['CDELT3']/1.0e6						#	Frequency increment in MHz	
    freqlen	=	testhd['NAXIS3']							#	Length of frequency axis
    freq	=	refreq + dfreq*((freqlen/2)-refloc)			#	Central frequency in MHz
    print('Central frequency		= %f MHz'%freq)
    print('Channel width			= %f MHz'%dfreq)
    reffile.close()	
    
    srcs	=	[]
    spzfile	=	fits.open(pars['CatDir']+'/'+pars['SpzCat'])
    sedfile	=	fits.open(pars['CatDir']+'/'+pars['SedCat'])	
    spdata	=	spzfile[1].data
    sedata	=	sedfile[1].data
    spzid	=	spdata['Id_COS20_Classic']
    sedid	=	sedata['Id']
    srcra	=	spdata['ra_corrected']							#	source RAs
    srcdec	=	spdata['dec_corrected']							#	source DECs
    srcz	=	spdata['specz']								    #	source redshifts
    srczq	=	spdata['flag']								    #	source redshift quality, Good = 3,4,13,14
    sedz	=	sedata['zs']								    #	source redshifts
    srctype	=	spdata['photoz_type']						    #	source type, 0 = galaxy, 1 = star, 2 = AGN

    srcnuv	=	sedata['MNUV']							        #	Abs NUV mag
    srcr	=	sedata['MR']								    #	Abs R mag
    srcsfr  =   sedata['SFR_best']                              #   SFR
    srclsm  =   sedata['mass_best']                             #   log(stellar mass)

    srcn	=	len(srcra)									    #	number of sources 

    if (len(srcz) != len(sedz)):
        print("Catalogue mismatch !!!")
        sys.exit()

    for i in range (0,srcn):
        if (np.isin(srczq[i], np.array([3,4,13,14]))):
            if (np.abs(srcz[i] - sedz[i]) > pars['ZTol']):
                print("Redshift mismatch !!!", i, srcz[i], sedz[i], srczq[i])
                sys.exit()
            
            if (spzid[i] != sedid[i]):
                print("ID mismatch !!!", i, spzid[i], sedid[i])
                sys.exit()
            
            if ((srcz[i] - pars['ZLim'][0])*(srcz[i] - pars['ZLim'][1]) <= 0.0):	
                if (srclsm[i] > 0.0):
                    dtheta	=	distang(srcra[i],ra0,srcdec[i],dec0)/3600.0
                    obsfreq	=	pars['RestFreq']/(1.0+srcz[i])
                    obschan	=	int(round(((obsfreq-refreq)/dfreq) - refloc))
                    
                    if (dtheta <= pars['MaxRadeg']):
                        if ((obschan > 0) and (obschan < (freqlen-1))):
                            tmp	=	np.rint(w.wcs_world2pix(srcra[i],srcdec[i], pars['PiXys']))
                            px	=	int(tmp[0])
                            py	=	int(tmp[1])				
                                            
                            srcs.append([i, spzid[i], srcra[i], srcdec[i], srcz[i], srczq[i], dtheta, px, py, obschan, \
                                         srcsfr[i], srclsm[i], srcnuv[i], srcr[i], srctype[i]])
                            #	i	id	ra	dec	z	zq  d   px	py	ochan   SFR	lsm	NUV R   Type
                            #	0	1	2	3	4	5	6	7	8	9	    10	11	12  13  14          
    srcs	=	np.array(srcs)
    srcs	=	srcs[np.argsort(srcs[:,4])]	
    spzfile.close()	
    sedfile.close()		
    
    return (srcs)			
#	---------------------------------------------------------------------------------------------------


def conmastercat(gdets, pars=None):
    
    #   Function to construct the master sample		

    #	i	id	ra	dec	z	zq  d   px	py	ochan   SFR	lsm	NUV R   Type
    #	0	1	2	3	4	5	6	7	8	9	    10	11	12  13  14

    if (pars['ColType']=="blue"):
        stackid	= 'mastercat_blue'
        nuvrlim	= [-10000.0,4.0]
    elif (pars['ColType']=="red"):
        stackid	= 'mastercat_red'
        nuvrlim	= [4.0,100000.0]
    else:
        stackid	= 'mastercat_all'
        nuvrlim	= [-100000.0,100000.0]

    print('Sample size full		=	%d'%len(gdets))

    gdets	= gdets[(gdets[:,12] - gdets[:,13]) > nuvrlim[0]]
    gdets	= gdets[(gdets[:,12] - gdets[:,13]) < nuvrlim[1]]
    print('Within NUV-r limit		=	%d'%len(gdets))

    gdets	= gdets[gdets[:,4] >= pars['ZLim'][0]]
    gdets	= gdets[gdets[:,4] <= pars['ZLim'][1]]
    print('Within redshift range		=	%d'%len(gdets))

    gdets	= gdets[gdets[:,11] >  pars['LsmLim'][0]]
    gdets	= gdets[gdets[:,11] <= pars['LsmLim'][1]]
    print('Within lSM limit		=	%d'%len(gdets))

    return (gdets, stackid) 
#	--------------------------------------------------------------------------------------------------













