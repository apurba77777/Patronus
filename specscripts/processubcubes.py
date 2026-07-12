import os,sys
import subprocess
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import casatasks as ct
import casatools as ctools
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
			    region = 'circle[ [ '+str(ra)+'deg , '+str(dec)+'deg] ,'+str(1.0)+'pix]', \
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


def intercubes(detid, detz, posra, posdec, obschan, istart=0, pars=None, ovrt=True, nobj=-1):
    
    #   Spectrally interpolate subcubes to a common velocity grid

    chanstat    = np.loadtxt(pars['CubeDir']+"/"+pars['ChanFile'])
    ntotal  = len(detid)
    if (nobj > 0):
        ntotal  = nobj

    iman    = ctools.image()  

    for i in range (istart, ntotal):
        ra      = posra[i]
        dec     = posdec[i]

        if (pars['SmKpc'] > 0):
            print("Spatial smoothing has not yet been implemented...\n")
        else:
            pmax	= np.ones(2*pars['FSize']+1, dtype=float)
        
        print('Interpolating source =	%d / %d	(%ld, RA = %f, dec = %f, z = %f)'%(i, ntotal, detid[i], ra, dec, detz[i]))
        

        subprocess.run(["cp", "-r", pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['ScubeDir']+"/subcube_"+str(detid[i])+"_"+str(pars['SmKpc'])+".im",
                        pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['IscubeDir']+"/subcube_"+str(detid[i])+"_"+str(pars['SmKpc'])+".im"])


        specfile	= np.loadtxt(pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['SpecDir']+'/spec_'+str(detid[i])+'_'+str(pars['SmKpc'])+'.txt')

        rfreqarr	=	specfile[:,2]
        rvelarr		=	relvel(pars['RestFreq'], (1.0 + detz[i])*rfreqarr)
        ivelarr0	=	np.arange( -pars['FSize']*pars['VelRes'], (pars['FSize']+1)*pars['VelRes'], pars['VelRes'], dtype=float)
        ivelarr		=	ivelarr0[::-1]

        iman.open(pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['IscubeDir']+"/subcube_"+str(detid[i])+"_"+str(pars['SmKpc'])+".im")
        imarr		=	iman.getchunk()		
        carr		=	np.copy(imarr)
        ccarr		=	np.zeros(carr.shape, dtype='float32')

        gchans		=	np.ones(ivelarr.shape, dtype=int)

        for k in range (0, len(gchans)):		
	
            bc		= False			
            dvelarr	= rvelarr - ivelarr[k]
            nearest	= np.argmin(np.abs(dvelarr))
			
            if (dvelarr[nearest]<0):
                lowind	=	nearest-1
            else:
                lowind	=	nearest
			
            if (lowind < 0):
                bc		=	True
            elif (lowind >= 2*pars['FSize']):
                bc		=	True
            elif (np.isnan(pmax[k])):
                bc		=	True
			
            if (bc):
                gchans[k]	=	0				
            else:				
                if (chanstat[int(obschan[i])-pars['FSize']+lowind]==0 or chanstat[int(obschan[i])-pars['FSize']+lowind+1]==0):
                    gchans[k]	=	0	
			
                else:
                    wl				=	(rvelarr[lowind]-ivelarr[k])/(rvelarr[lowind]-rvelarr[lowind+1])	
                    wr				=	(ivelarr[k]-rvelarr[lowind+1])/(rvelarr[lowind]-rvelarr[lowind+1])
                    ccarr[:,:,k]	=	(wl * carr[:,:,lowind]/pmax[lowind]) + (wr * carr[:,:,lowind+1]/pmax[lowind+1])				
		
        iman.putchunk(ccarr)
        iman.done()

        #	Extract spectra		.................................	
            
        ct.specflux(
            imagename = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['IscubeDir']+"/subcube_"+str(detid[i])+"_"+str(pars['SmKpc'])+".im", \
            region = 'circle[ [ '+str(ra)+'deg , '+str(dec)+'deg] ,'+str(1.0)+'pix]', \
            unit = 'MHz', \
            function = 'mean', \
            logfile = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['IntSpecs']+'/spec_'+str(detid[i])+'_'+str(pars['SmKpc'])+'.txt', \
            major = "20.0arcsec", \
            minor = "20.0arcsec", \
            overwrite = ovrt
        )

        #	Extract noise spectra	..........................	
					
        q1	= (1+detz[i])*(pars['N0Kpc']/2)*180*3600/(d_com_gpc(detz[i])*1.0e6*np.pi)
        q2	= (1+detz[i])*((pars['N0Kpc']/2) + pars['NKpc'])*180*3600 / (d_com_gpc(detz[i])*1.0e6*np.pi)

        xstat	=   ct.imstat(
                        imagename = pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['IscubeDir']+"/subcube_"+str(detid[i])+"_"+str(pars['SmKpc'])+".im", \
                        axes = [0,1], \
				        region = 'annulus[ [ '+str(ra)+'deg , '+str(dec)+'deg], ['+str(q1)+'arcsec , '+str(q2)+'arcsec] ]'
                    )

        noisespec	=	xstat['sigma']	
        noisespec[gchans==0] = 0.0
        np.savetxt(pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['InoiseDir']+'/spec_'+str(detid[i])+'_'+str(pars['SmKpc'])+'_noise.txt',noisespec)

        ct.exportfits(
            imagename=pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['IscubeDir']+"/subcube_"+str(detid[i])+"_"+str(pars['SmKpc'])+".im", \
            fitsimage=pars['WorkDir']+'/'+pars['SubDir']+'/'+pars['IscubeDir']+"/xubcube_"+str(detid[i])+"_"+str(pars['SmKpc'])+".fits", \
            overwrite = ovrt
        )

    iman.close()

    return 
#	----------------------------------------------------------------------------------------------------