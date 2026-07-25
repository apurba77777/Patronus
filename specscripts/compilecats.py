import os,sys
from astropy.io import fits
from astropy.table import Table
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
    
    reffile = fits.open(pars['CubeDir']+'/'+pars['ImgCube']+".fits")    
    testhd	= reffile[0].header							#	image header
    w		= WCS(testhd,naxis=2)

    if (pars['CenRa']==None):
        ra0	= testhd['CRVAL1']							#	Central RA
    else:
        ra0	= pars['CenRa']
    
    if (pars['CenDec']==None):
        dec0= testhd['CRVAL2']							# 	Central DEC
    else:
        dec0= pars['CenDec']

    refreq	= testhd['CRVAL3']/1.0e6					#	Reference frequency in MHz
    refloc	= testhd['CRPIX3'] - pars['PiXys']			#	Reference location in frequency axis
    dfreq	= testhd['CDELT3']/1.0e6					#	Frequency increment in MHz	
    freqlen	= testhd['NAXIS3']							#	Length of frequency axis
    freq	= refreq + dfreq*((freqlen/2)-refloc)		#	Central frequency in MHz
    print('Central frequency	= %f MHz'%freq)
    print('Channel width		= %f MHz'%dfreq)
    reffile.close()	
    
    srcs	= []
    spdata	= Table.read(pars['CatDir']+'/'+pars['SpzCat'], format='fits')
    sedata	= Table.read(pars['CatDir']+'/'+pars['SedCat'], format='fits')	

    matching= np.isin(spdata['Id_specz'], sedata['Id_specz'])
    spdata  = spdata[matching]

    spdata.sort(['Id_specz'])
    sedata.sort(['Id_specz'])

    spzid	= spdata['Id_COS20_Classic']
    srcra	= spdata['ra_corrected']					#	source RAs
    srcdec	= spdata['dec_corrected']					#	source DECs
    srcz	= spdata['specz']							#	source redshifts
    srczq	= spdata['flag']							#	source redshift quality, Good = 3,4,13,14
    srctype	= spdata['photoz_type']						#	source type, 0 = galaxy, 1 = star, 2 = AGN
    srcsur  = spdata['survey']                          #   Survey IDs
    
    sedid	= sedata['Id_specz']    
    sedz	= sedata['specz']						    #	source redshifts
    srcnuvr	= sedata['best.param.restframe_galex.NUV-subaru.suprime.r']		    #	NUV - r
    srcrj   = sedata['best.param.restframe_subaru.suprime.r-paranal.vircam.J']  #   r - j
    srcuv   = sedata['best.param.restframe_cfht.megacam.u-subaru.suprime.V']    #   U - V
    srcvj   = sedata['best.param.restframe_subaru.suprime.V-paranal.vircam.J']  #   V - J

    srcsfr      = sedata['best.sfh.sfr']                    #   Instanteneous SFR
    srcsfr10    = sedata['best.sfh.sfr10Myrs']              #   Instanteneous SFR
    srcsfr100   = sedata['best.sfh.sfr100Myrs']             #   Instanteneous SFR
    srclsm      = np.log10(sedata['best.stellar.m_star'])   #   log(stellar mass)
    srclgm      = np.log10(sedata['best.stellar.m_gas'])    #   log(gas mass)

    srcn	= len(srcra)	
    print(f"Starting with {srcn} objects")							    #	number of sources 

    if (len(srcz) != len(sedz)):
        print("Catalogue mismatch !!!")
        print(f"{len(srcz)} vs {len(sedz)}")
        sys.exit()

    for i in range (0,srcn):
        if (np.isin(srczq[i], np.array(pars['ZQuality'])) and np.isin(srcsur[i], np.array(pars['SurveyIds']))):
            if (np.abs(srcz[i] - sedz[i]) > pars['ZTol']):
                print("Redshift mismatch !!!", i, srcz[i], sedz[i], srczq[i])
                sys.exit()
            
            if (spdata['Id_specz'][i] != sedid[i]):
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
                                            
                            srcs.append([i, spzid[i], srcra[i], srcdec[i], srcz[i], srczq[i], dtheta, px, py, \
                                         obschan, srcsfr[i], srclsm[i], srcnuvr[i], srcrj[i], srctype[i], \
                                            srcuv[i], srcvj[i], srcsfr10[i], srcsfr100[i], srclgm[i]])

        #	i	id	ra	dec	z	zq  d   px	py	ochan   SFR	lsm	NUV-r r-J  Type  U-V  V-J   SFR10  SFR100  lgm
        #	0	1	2	3	4	5	6	7	8	9	    10	11	12    13   14    15   16    17     18      19
                                 
    srcs	=	np.array(srcs)
    srcs	=	srcs[np.argsort(srcs[:,4])]		
    
    return (srcs)			
#	---------------------------------------------------------------------------------------------------


def conmastercat(gdets, pars=None):
    
    #   Function to construct the master sample		

    #	i	id	ra	dec	z	zq  d   px	py	ochan   SFR	lsm	NUV-r r-J  Type  U-V  V-J   SFR10  SFR100  lgm
    #	0	1	2	3	4	5	6	7	8	9	    10	11	12    13   14    15   16    17     18      19

    if (pars['ColType']=="blue"):
        nuvrlim	= [-10000.0,4.0]
    elif (pars['ColType']=="red"):
        nuvrlim	= [4.0,100000.0]
    else:
        nuvrlim	= [-100000.0,100000.0]

    print('Sample size full		=	%d'%len(gdets))

    gdets	= gdets[gdets[:,4] >= pars['ZLim'][0]]
    gdets	= gdets[gdets[:,4] <= pars['ZLim'][1]]
    print('Within redshift range		=	%d'%len(gdets))

    gdets	= gdets[gdets[:,11] >  pars['LsmLim'][0]]
    gdets	= gdets[gdets[:,11] <= pars['LsmLim'][1]]
    print('Within lSM limit		=	%d'%len(gdets))

    return (gdets) 
#	--------------------------------------------------------------------------------------------------













