import os,sys
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import casatasks as ct
from specscripts.auxfns import *

#   ---------------------------------------------------------------------------------------------------
#
#       Functions to process subcubes 
#
#   ---------------------------------------------------------------------------------------------------


def exubcubes(detid, detz, posra, posdec, obschan, istart=0, pars=None, ovrt=True, nobj=-1):
    
    #   Extract subcubes for a given catalogue
    
    detkept	= []	
    ntotal  = len(detid)
    if (nobj > 0):
        ntotal  = nobj
        
    for i in range (istart, ntotal):

        bchan	= int(obschan[i] - pars['FSize'])
        echan	= int(obschan[i] + pars['FSize'])
        ra      = posra[i]
        dec     = posdec[i]

        if ((bchan >= 0) and (echan < pars['FLen'])):
            detkept.append(i)
        
            print('Extracting source =	%d / %d	(%ld, RA = %f, dec = %f, z = %f)'%(i, ntotal, detid[i], ra, dec, detz[i]))

            ct.imsubimage(
                imagename = pars['CubeDir']+'/'+pars['ImgCube']+'.image', \
                outfile = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['ScubeDir']+'/subcube_'+str(detid[i])+'_0.im', \
                region = 'centerbox[['+str(ra)+'deg ,'+str(dec)+'deg],['+str(pars['BxSize'])+'pix,'+str(pars['BxSize'])+'pix]]', \
                chans = str(bchan)+'~'+str(echan), \
                overwrite = ovrt
            ) 

            print("Exporting FITS image")
            ct.exportfits(
                imagename = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['ScubeDir']+'/subcube_'+str(detid[i])+'_0.im', \
                fitsimage = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['ScubeDir']+'/subcube_'+str(detid[i])+'_0.fits', \
                overwrite = ovrt
            )

            q1	= (1+detz[i])*(pars['N0Kpc']/2)*180*3600/(d_com_gpc(detz[i])*1.0e6*np.pi)
            q2	= (1+detz[i])*((pars['N0Kpc']/2) + pars['NKpc'])*180*3600 / (d_com_gpc(detz[i])*1.0e6*np.pi)

            xstat = ct.imstat(
                        imagename = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['ScubeDir']+'/subcube_'+str(detid[i])+'_0.im', \
                        axes = [0,1], \
                        region = 'annulus[ [ '+str(ra)+'deg , '+str(dec)+'deg], ['+str(q1)+'arcsec , '+str(q2)+'arcsec] ]'
                    )

            noisespec	= xstat['sigma']	
            print ("Noise length	= %d"%len(noisespec))
            np.savetxt(pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['NoiseDir']+'/spec_'+str(detid[i])+'_noise_0.txt', noisespec)

            ct.specflux(
                imagename = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['ScubeDir']+'/subcube_'+str(detid[i])+'_0.im', \
			    region = 'circle[ [ '+str(ra)+'deg , '+str(dec)+'deg] ,'+str(0.5)+'pix]', \
			    unit = 'MHz', \
                function = 'mean', \
                logfile = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['SpecDir']+'/spec_'+str(detid[i])+'_0.txt', \
			    major = "20.0arcsec", \
                minor = "20.0arcsec", \
                overwrite = ovrt
            )

        else:
            print(f"Rejected source	{i} ({bchan} / {echan})")

    detkept	= np.array(detkept)

    return (detkept) 
#	----------------------------------------------------------------------------------------------------