import os,sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import pickle as pkl
from collections import namedtuple
import matplotlib.pyplot as plt
from specscripts.auxfns import *
from specscripts.stackfns import *
from specscripts.stackplotfns import *

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
                      int(idata.shape[2]/2) - pars['BSize'] : int(idata.shape[2]/2) + pars['BSize'] ]

        icarr	= icarr/detpb									    #	Flux density correcting for primary beam
        icarr	= icarr*(4*np.pi*dlgpc*dlgpc)/(1.0+detz)		    #	Luminosity density

        allcubefull.append(icarr)
        icube.close()
    
    allnoiseful	= np.array(allnoiseful)
    allcubefull	= np.array(allcubefull)
        
    #print(allnoiseful.shape, allcubefull.shape)	
    print(f"Cubes stacked	= {ngal}")	

    velarr,stackarr,velavg,stackavg	=	stackcube_lumz_fluxwt(allnoiseful, pars['VelRes'], allcubefull, \
                                                        pars['RmsPow'], domedian)	

    stackarr	= subase(stackarr, pars['FitOrder'], pars['ExclChans'], pars['VelRes'])
    stackavg	= subase(stackavg, pars['FitOrder'], pars['ExAvgChans'], 2*pars['VelRes'])
    rmsarr		= planerms(stackarr, pars['ExRad'], pars['ExEdge'])
    rmsavg		= planerms(stackavg, pars['ExRad'], pars['ExEdge'])

    if (domedian):
        print("Saving median cube...")
        gsamp	= gsamp._replace(medstkcube = stackarr, medstkcubeavg = stackavg, medplnrmsarr = rmsarr, \
                               medplnrmsavg = rmsavg)
    else:
        print("Saving mean cube...")
        gsamp	= gsamp._replace(meanstkcube = stackarr, meanstkcubeavg = stackavg,	meanplnrmsarr = rmsarr, \
                               meanplnrmsavg = rmsavg)

    return (gsamp)
#   ---------------------------------------------------------------------------------------------------


def stackerrors(gsamp, getspecdir, pars=None):

    #   Calculate errors

    domedian	= (pars['StackStat']=="median")
    if (domedian):
        print("Stacking for median cube ...")
    else:
        print("Stacking for mean cube ...")

    gdets   = gsamp.sampcat
    print('Total cubes to be stacked	=	%d'%(len(gdets)))		

    allnoiseful	= []
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
	
        icarr0  = ispecful[ivelreg,3]
        icarr0	= icarr0/detpb									    #	Flux density correcting for primary beam
        icarr0	= icarr0*(4*np.pi*dlgpc*dlgpc)/(1.0+detz)		    #	Luminosity density

        allspecfull.append(icarr0)
    
    allnoiseful	= np.array(allnoiseful)
    allspecfull	= np.array(allspecfull)    
        
    #print(allnoiseful.shape, allcubefull.shape)	
    print(f"Spectra stacked	= {ngal}")	

    if (domedian):
        print("Using median cube...")
        planermsarr		=	gsamp.medplnrmsarr
        planermsavg		=	gsamp.medplnrmsavg 
    else:
        print("Using mean cube...")
        planermsarr		=	gsamp.meanplnrmsarr
        planermsavg		=	gsamp.meanplnrmsavg 

    velarr0, stackarr0, rmsnoisenormal, velarravg0, stackarravg0, avgrmsnoisenormal	= \
        stackspecs_int_lumz_fluxwt_singpix(allnoiseful, pars['VelRes'], allspecfull, pars['FitOrder'], \
                                           pars['ExclChans'], pars['RmsPow'], domedian)

    lind	= int(len(velarr0)/2)
    rind	= int(len(velarr0)/2)
    psnrarr	= stackarr0/planermsarr

    while(psnrarr[lind-1] >= pars['Thresh']):
        lind	= lind-1
    while(psnrarr[rind+1] >= pars['Thresh']):
        rind	= rind+1	
    velmin		= velarr0[lind] - pars['VelRes']/2
    velmax		= velarr0[rind] + pars['VelRes']/2
    print("Integrating from %.1f to %.1f km/s"%(velmin,velmax))	
    lumint0		= np.sum(stackarr0[lind:rind+1])
    masslumv0	= masslumvel_noz(lumint0*1.0e6/(4*np.pi), pars['VelRes'])

    randmislum	= intrms_lum_fluxwt_random(allspecfull,allnoiseful,pars['Realn'],pars['RmsPow'],[lind,rind],domedian)

    avglind		= int(len(velarravg0)/2)
    avgrind		= int(len(velarravg0)/2)
    psnravg		= stackarravg0/planermsavg

    while(psnravg[avglind-1] >= pars['Thresh']):
        avglind	= avglind-1
    while(psnravg[avgrind+1] >= pars['Thresh']):
        avgrind	= avgrind+1
    velminavg	= velarravg0[avglind] - pars['VelRes']
    velmaxavg	= velarravg0[avgrind] + pars['VelRes']
    print("Integrating from %.1f to %.1f km/s"%(velminavg,velmaxavg))		
    lumintavg0	= np.sum(stackarravg0[avglind:avgrind+1])
    masslumvavg0= masslumvel_noz(lumintavg0*1.0e6/(4*np.pi), 2*pars['VelRes'])

    chanw		= rind+1-lind
    randmass	= masslumvel_noz(randmislum*1.0e6/(4*np.pi), pars['VelRes'])	
    print("Average HI mass =  %.2e	Random err = %.2e"%(masslumv0,randmass))	
    print("Average HI mass =  %.2e	Random err = %.2e"%(masslumvavg0,randmass))

    # Jackknife resampling

    templumarr		= []
    templumarravg	= []

    for rn in range(0, ngal):
        specfulljack= np.delete(allspecfull, rn, axis=0)
        noisefuljack= np.delete(allnoiseful, rn, axis=0)		
        velarr, stackarr, rmsnoisenormal, velarravg, stackarravg, avgrmsnoisenormal	= \
            stackspecs_int_lumz_fluxwt_singpix(noisefuljack, pars['VelRes'], specfulljack, pars['FitOrder'], \
                                               pars['ExclChans'], pars['RmsPow'], domedian)
        lumint	    = np.sum(stackarr[lind:rind+1])		
        templumarr.append(lumint)	
        lumintavg	= np.sum(stackarravg[avglind:avgrind+1])		
        templumarravg.append(lumintavg)	

    templumarr		= np.array(templumarr)
    templumarravg	= np.array(templumarravg)	
    jkavglum	    = np.mean(templumarr)	
    jkstderr	    = np.std(templumarr)*np.sqrt(float(ngal-1))
    jkavglumavg	    = np.mean(templumarravg)	
    jkstderravg	    = np.std(templumarravg)*np.sqrt(float(ngal-1))
    jkmassavg	    = masslumvel_noz(jkavglum*1.0e6/(4*np.pi),pars['VelRes'])
    jkmasserr	    = masslumvel_noz(jkstderr*1.0e6/(4*np.pi),pars['VelRes'])
    jkmassavgavg	= masslumvel_noz(jkavglumavg*1.0e6/(4*np.pi),2*pars['VelRes'])
    jkmasserravg	= masslumvel_noz(jkstderravg*1.0e6/(4*np.pi),2*pars['VelRes'])
    logmassarr		= np.log10(masslumvel_noz(templumarr*1.0e6/(4*np.pi),pars['VelRes']))
    logmassarravg	= np.log10(masslumvel_noz(templumarravg*1.0e6/(4*np.pi),2*pars['VelRes']))
    jkmassavglog	= np.mean(logmassarr)
    jkmasserrlog	= np.std(logmassarr)*np.sqrt(float(ngal-1))
    jkmassavglogavg	= np.mean(logmassarravg)
    jkmasserrlogavg	= np.std(logmassarravg)*np.sqrt(float(ngal-1))

    print("Jacknife avg =  %.2e	Err = %.2e"%(jkmassavg, jkmasserr))
    print("Jacknife log =  %.2f	Err = %.2f"%(jkmassavglog, jkmasserrlog))
    print("Jacknife avg =  %.2e	Err = %.2e"%(jkmassavgavg, jkmasserravg))
    print("Jacknife log =  %.2f	Err = %.2f"%(jkmassavglogavg, jkmasserrlogavg))

    if (domedian):
        print("Saving median values...")
        gsamp	= gsamp._replace(emedmhrand = randmass, emedmhjk = jkmasserr, emedmhjkavg = jkmasserravg)
    else:
        print("Saving mean values...")
        gsamp	= gsamp._replace(emeanmhrand = randmass, emeanmhjk = jkmasserr, emeanmhjkavg = jkmasserravg)

    return (gsamp)
#   ---------------------------------------------------------------------------------------------------


def stackmasses(gsamp, pars=None):

    #   Calculate stacked masses and plot

    resplotdir  = f"{pars["WorkDir"]}/{pars["ResplotDir"]}/"
    domedian	= (pars['StackStat']=="median")
    if (domedian):
        print("Stacking for median cube ...")
    else:
        print("Stacking for mean cube ...")

    gdets		= gsamp.sampcat
    print('Total cubes stacked	= %d'%(len(gdets)))

    meanz   = np.median(gdets[:,4])
    meancol	= 0.0
    meannuvr= 0.0
    meanlb	= 0.0
    meanmb	= 0.0
    meansm	= np.median(10.0**gdets[:,11])
    medsfr	= np.median(gdets[:,10])

    print("z = %.3f	Col = %.3f	NUV-r = %.3f"%(meanz, meancol, meannuvr))
    print("MB = %.3f	SM = %.2e   SFR = %.2e"%(meanmb, meansm, medsfr))    

    if (domedian):
        print("Fetching median cube ...")
        stackarr	= gsamp.medstkcube	
        stackavg	= gsamp.medstkcubeavg
        planermsarr	= gsamp.medplnrmsarr
        planermsavg	= gsamp.medplnrmsavg
        jkmasserr	= gsamp.emedmhjk
        jkmasserravg= gsamp.emedmhjkavg
        randmasserr	= gsamp.emedmhrand
    else:
        print("Fetching mean cube ...")
        stackarr	= gsamp.meanstkcube	
        stackavg	= gsamp.meanstkcubeavg
        planermsarr	= gsamp.meanplnrmsarr
        planermsavg	= gsamp.meanplnrmsavg
        jkmasserr	= gsamp.emeanmhjk
        jkmasserravg= gsamp.emeanmhjkavg
        randmasserr	= gsamp.emeanmhrand

    stackarr	= 1.0e6*subase(stackarr, pars["FitOrder"], pars["ExclChans"], pars["VelRes"])
    stackavg	= 1.0e6*subase(stackavg, pars["FitOrder"], np.array(pars["ExclChans"])/2, 2*pars["VelRes"])
    planermsarr	= planermsarr*1.0e6
    planermsavg	= planermsavg*1.0e6
    lumspec		= stackarr[:, pars["BSize"], pars["BSize"]]
    lumspecavg	= stackavg[:, pars["BSize"], pars["BSize"]]

    fsize		= int(stackarr.shape[0]/2)	
    velarr		= np.linspace(fsize*pars["VelRes"], -fsize*pars["VelRes"], 2*fsize+1)
    fsizeavg	= int(stackavg.shape[0]/2)	
    velavg		= np.linspace(fsizeavg*2*pars["VelRes"], -fsizeavg*2*pars["VelRes"], 2*fsizeavg+1)	

    lind		= int(len(velarr)/2)
    rind		= int(len(velarr)/2)	
    psnrarr		= lumspec/planermsarr

    while(psnrarr[lind-1] >= pars["Thresh"]):
        lind	= lind-1
    while(psnrarr[rind+1] >= pars["Thresh"]):
        rind	= rind+1	
        
    avglind		= int(len(velavg)/2)
    avgrind		= int(len(velavg)/2)	
    psnravg		= lumspecavg/planermsavg

    while(psnravg[avglind-1] >= pars["Thresh"]):
        avglind	= avglind-1
    while(psnravg[avgrind+1] >= pars["Thresh"]):
        avgrind	= avgrind+1	

    hmap		= np.nansum(stackarr[lind : rind+1], axis=0)
    hmapavg		= np.nansum(stackavg[avglind : avgrind+1], axis=0)

    velmin		= velarr[lind] + pars["VelRes"]/2
    velmax		= velarr[rind] - pars["VelRes"]/2
    print("\n%.1f to %.1f km/s"%(velmin,velmax))		
    lumint		= np.sum(lumspec[lind:rind+1])
    rmsplane	= np.sqrt(np.sum(planermsarr[lind:rind+1]**2))															
    masslumv	= masslumvel_noz(lumint/(4*np.pi), pars["VelRes"])
    planerr		= masslumvel_noz(rmsplane/(4*np.pi), pars["VelRes"])
    print("\nAverage mass 	=  %.2f	+/- %.2f (random) +/- %.2f (plane) +/- %.2f (jk)"\
                                %(masslumv/1.0e9,randmasserr/1.0e9,planerr/1.0e9,jkmasserr/1.0e9))
    print("Plane noise = ",rmsplane/(1+rind-lind))

    velminavg	= velavg[avglind] + pars["VelRes"]
    velmaxavg	= velavg[avgrind] - pars["VelRes"]
    print("\n%.1f to %.1f km/s"%(velminavg,velmaxavg))		
    lumintavg	= np.sum(lumspecavg[avglind:avgrind+1])
    rmsplaneavg	= np.sqrt(np.sum(planermsavg[avglind:avgrind+1]**2))	
    masslumvavg	= masslumvel_noz(lumintavg/(4*np.pi),2*pars["VelRes"])
    planerravg	= masslumvel_noz(rmsplaneavg/(4*np.pi),2*pars["VelRes"])
    print("\nAverage mass 	=  %.2f	+/- %.2f (random) +/- %.2f (plane) +/- %.2f (jk)"\
                                %(masslumvavg/1.0e9,randmasserr/1.0e9,planerravg/1.0e9,jkmasserravg/1.0e9))
    print("Averange plane noise = ",rmsplaneavg/(1+avgrind-avglind))

    if (domedian):
        print("Saving median values...")
        gsamp	= gsamp._replace(medlumint = lumint/(4*np.pi), medlumintavg = lumintavg/(4*np.pi), 
                               medmh9 = masslumv, medmh9avg = masslumvavg, emedpln = planerr, \
                                emedplnavg = planerravg, medstkspec = np.array([velarr,lumspec,planermsarr]).T, \
                                    medstkspecavg = np.array([velavg,lumspecavg,planermsavg]).T)
        
        stackspecplt(gsamp.medstkspec, pars["VelRes"], pars=pars)
        stackspecplt(gsamp.medstkspecavg, 2 * pars["VelRes"], pars=pars)
    else:
        print("Saving mean values...")
        gsamp	= gsamp._replace(meanlumint = lumint/(4*np.pi), meanlumintavg = lumintavg/(4*np.pi), \
                               meanmh9 = masslumv, meanmh9avg = masslumvavg, emeanpln = planerr, \
                                emeanplnavg = planerravg, meanstkspec = np.array([velarr,lumspec,planermsarr]).T, \
                                    meanstkspecavg = np.array([velavg,lumspecavg,planermsavg]).T)

        stackspecplt(gsamp.meanstkspec, pars["VelRes"], pars=pars)
        stackspecplt(gsamp.meanstkspecavg, 2 * pars["VelRes"], pars=pars)

    stackmaplt(hmap, pars["VelRes"], pars=pars)
    stackmaplt(hmapavg, 2*pars["VelRes"], pars=pars)

    return 
#   ---------------------------------------------------------------------------------------------------