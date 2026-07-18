import os,sys
import numpy as np
import pickle as pkl
from collections import namedtuple
import copy
from specscripts.auxfns import *

#   --------------------------------------------------------------------------------------------------------------
#
#       Functions to stack spectral cubes 
#
#   --------------------------------------------------------------------------------------------------------------


def stackcube_lumz_fluxwt(allnoiseful, velres, allcubefull, rmspow, medstack):	
	
    #   Stack spectral luminosity cubes weighting by flux density errors

	fsize	= int(allnoiseful.shape[1]/2)
	
	if(allnoiseful.shape[1]%2 == 0):
		print("fsize trouble !!!")
		return 0
	
	bsize	= allcubefull.shape[2]	
	carr0	= np.zeros((2*fsize+1, bsize, bsize),dtype=float)	
	velarr	= np.linspace(fsize*velres, -fsize*velres, 2*fsize+1)
	wtarr	= np.zeros(2*fsize+1,dtype=float)
		
	for i in range (0,2*fsize+1):	
		stackfluxrms	= allnoiseful[:,i] 
		ginds			= np.isfinite(stackfluxrms)		
		fluxrms			= stackfluxrms[ginds]
		islice			= allcubefull[:,i]
				
		if (len(fluxrms) > 0):			
			fluxwt		= 1.0/(fluxrms**rmspow)
			totalwt		= np.sum(fluxwt)
			wtarr[i]	= totalwt
			
			for mm in range (0, bsize):
				for nn in range (0, bsize):
					fluxarr	= islice[:,mm,nn]
					fluxarr	= fluxarr[ginds]
					if (medstack):
						carr0[i,mm,nn]	= np.nanmedian(fluxarr)
					else:
						carr0[i,mm,nn]	= np.nansum(fluxarr*fluxwt)/totalwt			
		else:
			carr0[i]	= np.nan
		
	#	Hanning smoothing and resampling
		
	f0			=	0	
	avgspeclen	=	int((2*fsize+1)/2)		
	if (avgspeclen%2 == 0):
		avgspeclen 	= avgspeclen - 1	
		f0			= 1
		
	halflen		= int(avgspeclen/2)
	avgcarr0	= np.zeros((avgspeclen,bsize,bsize),dtype=float)
	velarravg	= np.zeros(avgspeclen,dtype=float)
	
	for j in range (0,avgspeclen):
		velarravg[j]	= velarr[f0+2*j+1]					
		for mm in range (0,bsize):
			for nn in range (0,bsize):
				avgcarr0[j, mm,nn]	= 0.25*carr0[f0+2*j,mm,nn] + 0.5*carr0[f0+2*j+1,mm,nn] + \
					                    0.25*carr0[f0+2*j+2,mm,nn]		
							
	return (velarr, carr0, velarravg, avgcarr0)
#   --------------------------------------------------------------------------------------------------------------


def subase (carr0, fitpoly, linechans, velres):	
	
    #   Subtract spectral baseline

	fsize	= int(carr0.shape[0]/2)
	
	if(carr0.shape[0]%2 == 0):
		print("fsize trouble !!!")
		return 0
	
	bsize	= carr0.shape[1]
	velarr	= np.linspace(fsize*velres, -fsize*velres, 2*fsize+1)
	wtarr	= np.ones(2*fsize+1,dtype=float)
		
	wtarr[fsize+int(round(linechans[0])):fsize+int(round(linechans[1]))+1]	=	0.0																			
	for mm in range (0,bsize):
		for nn in range (0,bsize):
			basefitline		= np.poly1d(np.polyfit(velarr,carr0[:,mm,nn],deg=fitpoly, w=wtarr))	
			carr0[:,mm,nn]	= carr0[:,mm,nn] - basefitline(velarr)
								
	return (carr0)
#   --------------------------------------------------------------------------------------------------------------


def planerms(carr0, exrad, exedge):	
	
    #   Calculate RMS noise on the frequency planes of the stacked cube
	
	rmsarr	=	np.zeros(carr0.shape[0],dtype=float)
	
	for c in range (0,carr0.shape[0]):
		planearr	=	[]
		for mm in range (exedge,carr0.shape[1]-exedge):
			for nn in range (exedge,carr0.shape[2]-exedge):
				if(np.sqrt((mm-carr0.shape[1]/2)**2 + (nn-carr0.shape[2]/2)**2) > float(exrad)):
					planearr.append(carr0[c,mm,nn])		
									
		rmsarr[c]	=	np.nanstd(np.array(planearr))									

	return (rmsarr)
#   --------------------------------------------------------------------------------------------------------------


def singleplanerms(carr0, exrad, exedge):	
	
    #   Calculate RMS noise on a single frequency plane

	rmsval		=	0.0	
	planearr	=	[]
	for mm in range (exedge,carr0.shape[0]-exedge):
		for nn in range (exedge,carr0.shape[1]-exedge):
			if(np.sqrt((mm-carr0.shape[0]/2)**2 + (nn-carr0.shape[1]/2)**2) > float(exrad)):
				planearr.append(carr0[mm,nn])	
									
	rmsval		=	1.48*np.nanmedian(np.abs(np.array(planearr) - np.nanmedian(np.array(planearr))))									

	return(rmsval)
#   --------------------------------------------------------------------------------------------------------------


def cuberms_lum_fluxwt_random(allcubefull, allnoiseful, realn, rmspow):
	
    #   Calculate RMS noise in the stacked cube by random stacking

	nspec	=	allnoiseful.shape[0]
	speclen	=	allnoiseful.shape[1]
	bsize	=	allcubefull.shape[2]
	indarr	=	np.random.randint(0,speclen,size=(realn,nspec))	
	rmsarr	=	np.zeros((realn,bsize,bsize),dtype=float)
		
	for i in range (0,realn):		
		temp	=	np.zeros(np.shape(allcubefull[0,0]), dtype=float)
		wtotal	=	0.0
		
		for k in range (0,nspec):			
			wt		=	0.0
			indx	=	indarr[i,k]
			
			if(np.isfinite(allnoiseful[k,indx]) and (allnoiseful[k,indx] > 0.0)):
				wt	=	1.0/(allnoiseful[k,indx]**rmspow)			
				temp	+= 	wt*allcubefull[k,indx]
				wtotal	+=	wt
		
		if(wtotal>0.0):
			rmsarr[i]	=	temp/wtotal

	ensrms		=	np.std(rmsarr)

	return(ensrms)
#   --------------------------------------------------------------------------------------------------------------


def	find_h1(lumspec,lumrms,thresh):
	
    #   Find the HI 21 cm line
	
	lumsnr	=	lumspec/lumrms
	k0		=	int(len(lumspec)/2)
	klow	=	k0 - 1
	
	while (klow >=0 and lumsnr[klow]>thresh):
		klow	=	klow - 1
	
	klow	=	klow + 1
	
	if (klow==k0):
		while (klow<len(lumspec) and lumsnr[klow]<thresh):
			klow	=	klow + 1
	
	khi		=	klow + 1
	while (khi<len(lumspec) and lumsnr[khi]>thresh):
		khi 	=	khi + 1
	
	w0		=	khi - klow
	if (khi>=klow):
		lumax	=	max(lumspec[klow-1:khi+2])
	else:
		lumax	=	0.0
		
	low50	=	klow
	while (low50<len(lumspec) and lumspec[low50]<lumax*0.5):
		low50	=	low50 + 1
		
	hi50	=	low50 + 1
	while (hi50<len(lumspec) and lumspec[hi50]>0.5*lumax):
		hi50	=	hi50 + 1

	w50		=	hi50 - low50
	
	lumint	=	np.sum(lumspec[klow:khi])	
	rmsint	=	np.sqrt(np.sum(lumrms[klow:khi]**2))

	return (lumint, rmsint, w50, w0, low50, hi50, klow, khi)
#   --------------------------------------------------------------------------------------------------------------


def stackspecs_int_lumz_fluxwt_singpix (allnoiseful, velres, allspecfull, fitpoly, linechans, rmspow, medstack):	
	
    #   Stack spectra from a single spatial pixel of the cubes

	fsize	=	int(allnoiseful.shape[1]/2)
	if (allnoiseful.shape[1]%2 == 0):
		print("fsize trouble !!!")
		return 0
			
	carr0		=	np.zeros(2*fsize+1,dtype=float)	
	velarr		=	np.linspace(fsize*velres, -fsize*velres, 2*fsize+1)
	wtarr		=	np.zeros(2*fsize+1,dtype=float)
	rmsnoisearr	=	np.zeros(2*fsize+1,dtype=float)	
			
	for i in range (0,2*fsize+1):	
		stackfluxrms	=	allnoiseful[:,i] 
		ginds			=	np.isfinite(stackfluxrms)		
		fluxrms			=	stackfluxrms[ginds]
		islice			=	allspecfull[:,i]
				
		if (len(fluxrms)):			
			fluxwt			=	1.0/(fluxrms**rmspow)
			totalwt			=	np.sum(fluxwt)
			wtarr[i]		=	totalwt		
			fluxarr			=	islice[ginds]
			if (medstack):
				carr0[i]		=	np.nanmedian(fluxarr)
			else:
				carr0[i]		=	np.nansum(fluxarr*fluxwt)/totalwt	
			rmsnoisearr[i]	=	np.sqrt(np.sum((fluxrms**2)*fluxwt))/totalwt		
		else:
			carr0[i]		=	np.nan
			rmsnoisearr[i]	=	np.nan	
	
	wtarr			=	np.ones(2*fsize+1,dtype=float)		
	wtarr[fsize+linechans[0]:fsize+linechans[1]+1]	=	0.0																			
	basefitline		=	np.poly1d(np.polyfit(velarr,carr0,deg=fitpoly, w=wtarr))	
	carr0			=	carr0 - basefitline(velarr)
	
	#	Hanning smoothing and resampling
	f0			=	0	
	avgspeclen	=	int((2*fsize+1)/2)		
	if (avgspeclen%2 == 0):
		avgspeclen 	= 	avgspeclen - 1	
		f0			=	1
	halflen		=	int(avgspeclen/2)
		
	avgcarr0		=	np.zeros(avgspeclen,dtype=float)
	velarravg		=	np.zeros(avgspeclen,dtype=float)
	avgrmsnoisearr	=	np.zeros(avgspeclen,dtype=float)
	
	for j in range (0,avgspeclen):
		velarravg[j]		=	velarr[f0+2*j+1]					
		avgcarr0[j]			=	0.25*carr0[f0+2*j] + 0.5*carr0[f0+2*j+1] + 0.25*carr0[f0+2*j+2]
		avgrmsnoisearr[j]	=	np.sqrt(0.25*(rmsnoisearr[f0+2*j]**2) + \
							        0.5*(rmsnoisearr[f0+2*j+1]**2) + 0.25*(rmsnoisearr[f0+2*j+2]**2))

	return(velarr,carr0,rmsnoisearr,velarravg,avgcarr0,avgrmsnoisearr)
#   --------------------------------------------------------------------------------------------------------------


def specrms_lum_fluxwt_random_fullspec(allspecfull, allnoiseful, realn, rmspow):
	
    #   Calculate RMS noise on the stacked spectrum by random stacking
	
	nspec	=	allnoiseful.shape[0]
	speclen	=	allnoiseful.shape[1]
	fsize	=	int(speclen/2)
	indarr	=	np.random.randint(0,speclen,size=(realn,nspec))	
	rmsarr	=	np.zeros((realn,speclen),dtype=float)
		
	for i in range (0,realn):		
		tempspecarr		=	np.zeros((nspec,speclen),dtype=float)
		wtotalarr		=	np.zeros(speclen,dtype=float)
		wtarr			=	np.zeros((nspec,speclen),dtype=float)
		
		for k in range (0,nspec):			
			indx0	=	indarr[i,k]
			
			rolnoise	=	np.roll(allnoiseful[k],indx0)
			rolspec		=	np.roll(allspecfull[k],indx0)
			
			tempspecarr[k]	=	rolspec
			wtarr[k]		=	1.0/(rolnoise**rmspow)
			
		wtotalarr	=	np.nansum(wtarr,axis=0)
		tempspecarr	=	tempspecarr*wtarr
		rmsarr[i]	=	np.nansum(tempspecarr,axis=0)/wtotalarr
		
	ensrmsarr		=	np.nanstd(rmsarr, axis=0)
	
	#	Hanning smoothing and resampling
	f0			=	0	
	avgspeclen	=	int((2*fsize+1)/2)		
	if (avgspeclen%2 == 0):
		avgspeclen 	= 	avgspeclen - 1	
		f0			=	1
	halflen		=	int(avgspeclen/2)	
	avgrmsarr	=	np.zeros(avgspeclen,dtype=float)
	
	for j in range (0,avgspeclen):				
		avgrmsarr[j]		=	np.sqrt((0.25*ensrmsarr[f0+2*j])**2 + \
						            (0.5*ensrmsarr[f0+2*j+1])**2 + (0.25*ensrmsarr[f0+2*j+2])**2)
				
	return(ensrmsarr,avgrmsarr)
#   --------------------------------------------------------------------------------------------------------------

