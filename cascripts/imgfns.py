import os,sys
import numpy as np
from astropy.io import fits
import casatasks as ct
import casatools
from casaplotms import plotms
import bdsf as sf
import matplotlib.pyplot as plt



def avgtarget (targetvis, pars=None):
    
    #   Channel average calibrated target visibilities


    avgvis  = targetvis+"_avg"

    try:
        os.system("rm -rf "+pars['WorkDir']+pars['ImgUvDir']+'/'+avgvis+".m*")
    except:
        print("Nothing to delete...")


    print("Channel averaging...\n")

    tavgdo = False 
    if (pars['CalTimeAvg'] != ""):
        print(f"Time-averaging to {pars['CalTimeAvg']}")
        tavgdo = True

    ct.mstransform(
        vis=pars['WorkDir']+pars['UvMsDir']+'/'+targetvis+".ms", \
        outputvis=pars['WorkDir']+pars['ImgUvDir']+'/'+avgvis+".ms", \
        datacolumn="data", \
        keepflags=False, \
        hanning=True, \
        chanaverage=True, \
        chanbin=pars['TarChanAvg'], \
        timeaverage=tavgdo, \
        timebin=pars['CalTimeAvg']
    )


    print("Generating diagnostic plots...")

    plotms(vis=pars['WorkDir']+pars['ImgUvDir']+'/'+avgvis+".ms", xaxis="frequency", yaxis="amp", gridrows=2, \
           height=1000, width=1000, gridcols=1, showgui=False)

    plotms(vis=pars['WorkDir']+pars['ImgUvDir']+'/'+avgvis+".ms", xaxis="frequency", yaxis="phase", gridrows=2, gridcols=1, rowindex=1, \
            plotindex=1, clearplots=False, plotfile=pars['WorkDir']+pars['LogDir']+"/target_"+avgvis+"_freq.png", \
                height=1000, width=1000, overwrite=True, showgui=False)

    plotms(vis=pars['WorkDir']+pars['ImgUvDir']+'/'+avgvis+".ms", xaxis="uvwave", yaxis="amp", gridrows=2, \
           height=1000, width=1000, gridcols=1, showgui=False)

    plotms(vis=pars['WorkDir']+pars['ImgUvDir']+'/'+avgvis+".ms", xaxis="uvwave", yaxis="phase", gridrows=2, gridcols=1, rowindex=1, \
            plotindex=1, clearplots=False, plotfile=pars['WorkDir']+pars['LogDir']+"/target_"+avgvis+"_uv.png", \
                height=1000, width=1000, overwrite=True, showgui=False)
    
    print("\n Done!\n")

    return (0)
#   -----------------------------------------------------------------------------------------------------



def imgtarget (targetvislist, imgname, dosavemodel=True, dointeractive=False, pars=None, clnmask=None):
    
    #   Image calibrated target  
    
    imgpre  = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+'_'+imgname

    print("\nVisibilities -- ",targetvislist)
    print("Imagename -- ",imgpre)    

    print("Clearing existing image components...\n")
    os.system("rm -rf "+imgpre+".*")

    savemod = 'none'
    if (dosavemodel):
        savemod = 'modelcolumn'

    hrad    = float(pars['ImgSize'][0])/2
    rgnmask = 'Circle[['+str(hrad)+'pix, '+str(hrad)+'pix],'+str(hrad)+'pix]'

    if (clnmask != None):
        print(f"Clean mask = {clnmask}")
        rgnmask = clnmask
    else:
        print(f"No cleanmask provided...")

    print("\nMaking image...\n")
    ct.tclean(
        vis=targetvislist, \
        imagename=imgpre, \
        datacolumn="corrected", \
        imsize=pars['ImgSize'], \
        cell=pars['CellSize'], \
        selectdata=True, \
        field=pars['TargetName'],\
        uvrange=pars['ImgUvLim'], \
        startmodel=pars['PreModel'], \
        specmode='mfs', \
        gridder='widefield', \
        wprojplanes=pars['WprojPln'], \
        pblimit=-0.1, \
        deconvolver='mtmfs', \
        scales=pars['DeconScls'], \
        smallscalebias=pars['SclBias'], \
        weighting='briggs',\
        robust=pars['WtRobust'], \
        uvtaper=pars['ImUvTaper'], \
        niter=pars['ImNiter'], \
        nsigma=pars['ClnSigma'], \
        interactive=dointeractive, \
        usemask='user', \
        mask=rgnmask, \
        pbmask=0.2, \
        savemodel=savemod        
    )
    
    print(" Done!\n")

    return (0)
#   -----------------------------------------------------------------------------------------------------



def selfcal (targetvis, calfile, scalint="100s", gcmode=None, pars = None):
    
    #   Self-calibrate target visibilities

    if (os.path.exists(calfile)):
        os.system("rm -rf "+calfile)
    
    scmode = "ap"
    apmode = ""
    if (gcmode!=None and gcmode=="p"):
        scmode="p"
        apmode='calonly'

    print("Self-calibrating in "+scmode+" mode...\n")
    
    ct.gaincal(
        vis=targetvis+".ms", \
        caltable=calfile, \
        uvrange=pars['CalUvLim'], \
        solint=scalint, \
        refant=pars['RefAntenna'], \
        minsnr=pars['MinSNR'], \
        calmode=scmode, \
        solmode='R', \
        rmsthresh=pars['OutThresh']
    )

    print("\nApplying calibration...\n")
    ct.applycal(vis=targetvis+".ms", gaintable=[calfile], applymode=apmode)


    print("\nPlotting solutions...\n")

    plotms(
        vis=calfile, \
        xaxis="time", \
        yaxis="amp", \
        gridrows=2, \
        height=1000, \
        width=1000, \
        gridcols=1, \
        coloraxis='antenna1', \
        showgui=False
    )

    plotms(
        vis=calfile, \
        xaxis="time", \
        yaxis="phase", \
        gridrows=2, \
        gridcols=1, \
        rowindex=1, \
        plotindex=1, \
        clearplots=False, \
        plotfile=pars['WorkDir']+pars['LogDir']+"/selfcal.png", \
        height=1000, \
        width=1000, \
        coloraxis='antenna1', \
        overwrite=True, \
        showgui=False
    )

    print("\n Done!\n")

    return (0)
#   -----------------------------------------------------------------------------------------------------



def flagcaltarget (tarfile, pars=None, pyx="python", ankdir=None, ankin=None, ovrt=False):
    
    #   Flag calibrated target data 
    
    
    if (os.path.exists(tarfile+".fits") and ovrt):
        os.system("rm -rf "+tarfile+".fits")
    

    #   ----------------------------------
    #   Flag target with aNKflag
    #   ----------------------------------
    
    if ( (ankin != None) and (os.path.exists(ankin+'.yml')) ):

        print("Exporting target FITS...\n")
        ct.exportuvfits(vis=tarfile+".ms", fitsfile=tarfile+".fits", datacolumn='corrected', overwrite=ovrt)

        
        print("Flagging the target aNKflag...\n")        
    
        if (not os.path.exists("glogout.dat")):
            os.system("cp "+ankdir+"/glogout.dat .")

        tarcmd = pyx +" " + ankdir + "/runank.py --ankdir " + ankdir + " --scratchdir ankscratch/ " + \
                    " --parfile " + ankin + " --infilename " + tarfile + " --outfilename " + tarfile+"_f" + \
                    " --logfile " +pars['WorkDir']+pars['LogDir']+"/sc_target_"+pars['TargetName'] + \
                    " --flagmode uvbin --targetype=normal --clearscratch --nthreads " + str(pars['FlgThreads'])

        print("Running \n" + tarcmd)
        os.system(tarcmd)

        if (os.path.exists(tarfile+"_f.ms")):
            os.system("rm -rf "+tarfile+"_f.ms")

        print("Converting FITS back to MS...\n")
        ct.importuvfits(vis=tarfile+"_f.ms", fitsfile=tarfile+"_f.fits")
        
        print("Generating diagnostic plots...")

        plotms(vis=tarfile+"_f.ms", xaxis="frequency", yaxis="amp", gridrows=2, \
               height=1000, width=1000, gridcols=1, showgui=False)

        plotms(vis=tarfile+"_f.ms", xaxis="frequency", yaxis="phase", gridrows=2, gridcols=1, rowindex=1, \
                plotindex=1, clearplots=False, plotfile=pars['WorkDir']+pars['LogDir']+"/sc_target_"+pars['TargetName']+"_freq.png", \
                    height=1000, width=1000, overwrite=True, showgui=False)

        plotms(vis=tarfile+"_f.ms", xaxis="row", yaxis="amp", gridrows=2, \
               height=1000, width=1000, gridcols=1, showgui=False)

        plotms(vis=tarfile+"_f.ms", xaxis="row", yaxis="phase", gridrows=2, gridcols=1, rowindex=1, \
                plotindex=1, clearplots=False, plotfile=pars['WorkDir']+pars['LogDir']+"/sc_target_"+pars['TargetName']+"_row.png", \
                    height=1000, width=1000, overwrite=True, showgui=False)
        


        #   Open and edit the flags in the MS files
        print("Copying flags to the original file...\n")
        wms     = casatools.ms()

        wms.open(tarfile+"_f.ms")
        wms.selectinit(datadescid=0)
        tarcalflgs  = wms.getdata(["flag"])["flag"]
        #print(bpcalflgs.shape)
        wms.close()        

        wms.open(tarfile+".ms", nomodify=False)
        wms.selectinit(datadescid=0)
        fldflg  = wms.getdata(["flag"])
        fldflg["flag"] = tarcalflgs
        #print(fldflg["flag"].shape)
        wms.putdata(fldflg)
        wms.close()

    
    print(" Done!\n")


    return (0)
#   -----------------------------------------------------------------------------------------------------



def finalimg (targetvislist, dosavemodel=True, pars=None):
    
    #   Attempt to produce the *final* continuum image 
    
    imgpre  = pars['WorkDir']+pars['ImgDir']+'/'+pars['FinImage']

    print("\nVisibilities -- ",targetvislist)
    print("Imagename -- ",imgpre)    

    print("Clearing existing image components...\n")
    os.system("rm -rf "+imgpre+".*")

    savemod = 'none'
    if (dosavemodel):
        savemod = 'modelcolumn'

    hrad    = float(pars['ImgSize'][0])/2
    finmask = 'Circle[['+str(hrad)+'pix, '+str(hrad)+'pix],'+str(hrad)+'pix]'

    print("\nMaking final continuum image...\n")
    ct.tclean(
        vis=targetvislist, \
        imagename=imgpre, \
        datacolumn="corrected", \
        imsize=pars['ImgSize'], \
        cell=pars['CellSize'], \
        selectdata=True, \
        field=pars['TargetName'],\
        uvrange=pars['FinUvLim'], \
        specmode='mfs', \
        gridder='widefield', \
        wprojplanes=pars['WprojPln'], \
        pblimit=0.1, \
        deconvolver='mtmfs', \
        scales=pars['DeconScls'], \
        smallscalebias=pars['SclBias'], \
        weighting='briggs',\
        robust=pars['ImRobust'], \
        uvtaper=pars['ImUvTaper'], \
        niter=pars['FiNiter'], \
        nsigma=pars['FinSigma'], \
        interactive=False, \
        usemask='user', \
        mask=finmask, \
        pbmask=0.2, \
        savemodel=savemod        
    )

    print("Exporting final image to Output Directory...")

    ct.exportfits(
        imagename=imgpre+".image.tt0", \
        fitsimage=pars['OutDir']+'/'+pars['FinImage']+"_continuum.fits", \
        overwrite=True
    )

    ct.exportfits(
        imagename=imgpre+".alpha", \
        fitsimage=pars['OutDir']+'/'+pars['FinImage']+"_alpha.fits", \
        overwrite=True
    )

    ct.exportfits(
        imagename=imgpre+".alpha.error", \
        fitsimage=pars['OutDir']+'/'+pars['FinImage']+"_alpha_error.fits", \
        overwrite=True
    )

    sfimg   = sf.process_image(
                pars['OutDir']+'/'+pars['FinImage']+"_continuum.fits", \
                adaptive_rms_box = True, \
                advanced_opts = True, \
                group_by_isl = False, \
                interactive = False, \
                thresh_isl = 5.0, \
                thresh_pix = 3.0
            )

    sfimg.write_catalog(
        outfile=pars['OutDir']+'/'+pars['FinImage']+"_src_catalogue.fits", \
        catalog_type='srl', \
        format='fits', \
        clobber=True
    )
    
    print(" Done!\n")

    return (0)
#   -----------------------------------------------------------------------------------------------------



def getuvsub (ivis, calfile, pars=None): 

    #   Prepare calibrated and model subtracted visibilities 

    finmodel0   = pars['WorkDir']+pars['ImgDir']+'/'+pars['FinImage']+'.model.tt0'
    finmodel1   = pars['WorkDir']+pars['ImgDir']+'/'+pars['FinImage']+'.model.tt1'

    targetvis   = pars['WorkDir']+pars['UvMsDir']+'/'+ivis+".ms"
    outfits     = pars['WorkDir']+pars['ImgUvDir']+'/'+ivis+"_uvsub.fits"

    print("\nApplying calibration...\n")
    ct.applycal(vis=targetvis, gaintable=[pars['WorkDir']+pars['ImgUvDir']+calfile], applymode='')

    imgpre  = "junk"

    print("\nVisibilities -- ",targetvis)
    print("Imagename -- ",imgpre)    

    os.system("rm -rf "+imgpre+".*")

    hrad    = float(pars['ImgSize'][0])/2
    finmask = 'Circle[['+str(hrad)+'pix, '+str(hrad)+'pix],'+str(hrad)+'pix]'

    print("\nMaking the junk image...\n")
    ct.tclean(
        vis=targetvis, \
        imagename=imgpre, \
        startmodel=[finmodel0, finmodel1], \
        datacolumn="corrected", \
        imsize=pars['ImgSize'], \
        cell=pars['CellSize'], \
        selectdata=True, \
        field=pars['TargetName'],\
        specmode='mfs', \
        gridder='widefield', \
        wprojplanes=pars['WprojPln'], \
        pblimit=0.1, \
        deconvolver='mtmfs', \
        scales=pars['DeconScls'], \
        smallscalebias=pars['SclBias'], \
        weighting='briggs',\
        robust=pars['ImRobust'], \
        uvtaper=pars['ImUvTaper'], \
        niter=0, \
        nsigma=10.0, \
        interactive=False, \
        usemask='user', \
        mask=finmask, \
        pbmask=0.2, \
        savemodel='modelcolumn'        
    )

    os.system("rm -rf "+imgpre+".*")

    print("\n Subtracting model visibilities...\n")
    ct.uvsub(vis=targetvis)

    print(" Exporting target FITS...\n")
    ct.exportuvfits(vis=targetvis, fitsfile=outfits, datacolumn='corrected', overwrite=True)
    
    print(" Done!\n")

    return (0)
#   -----------------------------------------------------------------------------------------------------



def flagavguvsub (tarfile, pars=None, pyx="python", ankdir=None, ankin=None, ovrt=False):
    
    #   Flag calibrated continuum subtracted visibilities and average in channel
    
    
    if (os.path.exists(tarfile+"_f.fits") and ovrt):
        os.system("rm -rf "+tarfile+"_f.fits")
    
    filetoavg   = tarfile

    #   ----------------------------------
    #   Flag target with aNKflag
    #   ----------------------------------
    
    if ( (ankin != None) and (os.path.exists(ankin+'.yml')) ):
        
        print("Flagging the target aNKflag...\n")        
    
        if (not os.path.exists("glogout.dat")):
            os.system("cp "+ankdir+"/glogout.dat .")

        tarcmd = pyx +" " + ankdir + "/runank.py --ankdir " + ankdir + " --scratchdir ankscratch/ " + \
                    " --parfile " + ankin + " --infilename " + tarfile + " --outfilename " + tarfile+"_f" + \
                    " --logfile " +pars['WorkDir']+pars['LogDir']+"/uvsub_"+pars['TargetName'] + \
                    " --flagmode uvbin --targetype=uvsub --clearscratch --nthreads " + str(pars['FlgThreads'])

        print("Running \n" + tarcmd)
        os.system(tarcmd)

        filetoavg   = tarfile+"_f"
    

    if (os.path.exists(filetoavg+".ms")):
        os.system("rm -rf "+filetoavg+".ms")

    print("Converting flagged FITS back to MS...\n")
    ct.importuvfits(vis=filetoavg+".ms", fitsfile=filetoavg+".fits")


    if (os.path.exists(filetoavg+"_avg.ms")):
        os.system("rm -rf "+filetoavg+"_avg.ms")


    print("Channel averaging...\n")
    ct.mstransform(
        vis=filetoavg+".ms", \
        outputvis=filetoavg+"_avg.ms", \
        datacolumn="data", \
        keepflags=False, \
        hanning=True, \
        chanaverage=True, \
        chanbin=pars['FinChanAvg']
    )



    print("Generating diagnostic plots...")

    plotms(vis=filetoavg+"_avg.ms", xaxis="frequency", yaxis="amp", gridrows=2, \
           height=1000, width=1000, gridcols=1, showgui=False)

    plotms(vis=filetoavg+"_avg.ms", xaxis="frequency", yaxis="phase", gridrows=2, gridcols=1, rowindex=1, \
            plotindex=1, clearplots=False, plotfile=pars['WorkDir']+pars['LogDir']+"/uvsub_"+pars['TargetName']+"_freq.png", \
                height=1000, width=1000, overwrite=True, showgui=False)

    plotms(vis=filetoavg+"_avg.ms", xaxis="row", yaxis="amp", gridrows=2, \
           height=1000, width=1000, gridcols=1, showgui=False)

    plotms(vis=filetoavg+"_avg.ms", xaxis="row", yaxis="phase", gridrows=2, gridcols=1, rowindex=1, \
            plotindex=1, clearplots=False, plotfile=pars['WorkDir']+pars['LogDir']+"/uvsub_"+pars['TargetName']+"_row.png", \
                height=1000, width=1000, overwrite=True, showgui=False)


    
    print(" Done!\n")


    return (0)
#   -----------------------------------------------------------------------------------------------------



def findsrcs (imgname, pars=None):
    
    #   Find sources in an image and make a catalogue

    ct.exportfits(
        imagename=imgname+".image.tt0", \
        fitsimage=imgname+".fits", \
        overwrite=True
    )

    sfimg   = sf.process_image(
                imgname+".fits", \
                rms_box=(int(min(pars['ImgSize'])/10), int(min(pars['ImgSize'])/40)), \
                rms_map=True, \
                thresh='fdr', \
                group_by_isl = False, \
                interactive = False, \
                thresh_isl = pars['IslThresh'], \
                thresh_pix = pars['PeakThresh']
            )
    
    sfimg.export_image(
        outfile = imgname+"_src_mask.fits", \
        clobber = True, \
        img_format = 'fits', \
        img_type = 'island_mask', \
        pad_image = True
    )

    refits  = fits.open(imgname+".fits")
    refhd   = refits[0].header

    tarfits = fits.open(imgname+"_src_mask.fits", mode="update")
    tarfits[0].header['BMAJ'] = refhd['BMAJ']
    tarfits[0].header['BMIN'] = refhd['BMIN']
    tarfits[0].header['BPA'] = refhd['BPA']
    tarfits[0].header['SPECSYS'] = refhd['SPECSYS']
    tarfits[0].header['VELREF'] = refhd['VELREF']
    tarfits[0].header['DATE-OBS'] = refhd['DATE-OBS']
    tarfits[0].header['TIMESYS'] = refhd['TIMESYS']
    tarfits.flush()

    refits.close()
    tarfits.close()

    sfimg.write_catalog(
        outfile=imgname+"_srcat.fits", \
        catalog_type='gaul', \
        format='fits', \
        clobber=True
    )

    ct.importfits(
        fitsimage=imgname+"_src_mask.fits", \
        imagename=imgname+"_src_mask.mask", \
        overwrite=True
    )

    return (0)
#   -----------------------------------------------------------------------------------------------------



def readsfcat (imgname, pars=None):
    
    #   Read a source catalogue

    fitscat = fits.open(imgname+"_srcat.fits")

    sfcat   = fitscat[1].data    
    sfcat   = sfcat[sfcat['S_Code']=='S']
    #print(sfcat['Peak_flux'])
    fitscat.close()

    return (sfcat)
#   -----------------------------------------------------------------------------------------------------



def checkselfcal (imgold, imgnew, pars=None):
    
    #   Check if self calibration has converged

    print(f"\n  Checking for convergence... \n")

    hasconverged    = True

    iman    = casatools.image()

    print(f"Reading image {imgold}...")
    iman.open(imgold+".image.tt0")
    imgoldata    = iman.getchunk(dropdeg=True)
    iman.done()

    print(f"Reading image {imgnew}...")
    iman.open(imgnew+".image.tt0")
    imgnewdata    = iman.getchunk(dropdeg=True)
    iman.done()

    iman.close()  

    boxize  = min(imgoldata.shape) * pars['CenFrac'] // 2  
    
    cropold = imgoldata[int(imgoldata.shape[0]//2 - boxize): int(imgoldata.shape[0]//2 + boxize), \
                      int(imgoldata.shape[1]//2 - boxize): int(imgoldata.shape[1]//2 + boxize)]
    
    cropnew = imgnewdata[int(imgnewdata.shape[0]//2 - boxize): int(imgnewdata.shape[0]//2 + boxize), \
                      int(imgnewdata.shape[1]//2 - boxize): int(imgnewdata.shape[1]//2 + boxize)]
    
    cropdiff= cropnew - cropold

    diffrms = 1.48 * np.nanmedian(np.abs(cropdiff))
    oldrms  = 1.48 * np.nanmedian(np.abs(cropold))
    newrms  = 1.48 * np.nanmedian(np.abs(cropnew))
    print(f"\n  Central noise: {oldrms:.1e} (old), {newrms: .1e} (new), {diffrms: .1e} (difference)")

    #   Convergence based on difference image

    if (diffrms / newrms > pars['TolRms']):
        hasconverged = False
        print(f"  Residual / new = {diffrms / newrms} \n")
        return (hasconverged)
    
    if (np.abs(newrms-oldrms)/oldrms > pars['TolFrac']):
        hasconverged = False
        print(f"  Fractional change = {np.abs(newrms-oldrms)/oldrms} \n")
        return (hasconverged)
    
    #plt.imshow(cropdiff.T / newrms, origin='lower', vmin=-3, vmax=5)
    #plt.show()

    oldcat  = readsfcat (imgold, pars)
    newcat  = readsfcat (imgnew, pars)

    oldcat  = oldcat[np.argsort(oldcat['Peak_flux'])]
    newcat  = newcat[np.argsort(newcat['Peak_flux'])]

    nsrcmin  = min(len(oldcat), len(newcat))

    if (nsrcmin > pars['MaxSrcs']):   
        oldcat  = oldcat[:pars['MaxSrcs']]
        newcat  = newcat[:pars['MaxSrcs']]
    else:
        oldcat  = oldcat[:nsrcmin]
        newcat  = newcat[:nsrcmin]

    #print(oldcat['Peak_flux'])
    #print(oldcat['Xposn'])
    #print(newcat['Peak_flux'])
    #print(newcat['Xposn'])

    xshift  = np.abs(oldcat['Xposn'] - newcat['Xposn']) 
    yshift  = np.abs(oldcat['Yposn'] - newcat['Yposn'])

    #   Convergence in source position

    if (max(max(xshift), max(yshift)) > pars['TolPos'] ):
        hasconverged = False
        print(f"  Maximum position offset = {max(xshift)}, {max(yshift)} \n")
        return (hasconverged)
    
    dpeak       = np.abs(oldcat['Peak_flux'] - newcat['Peak_flux']) 
    dpeakrel    = dpeak / newcat['E_Peak_flux']
    dpeakfrac   = dpeak / newcat['Peak_flux']
    
    if (max(dpeakfrac) > pars['TolPeak'] and max(dpeakrel) > 2.0 ):
        hasconverged = False
        print(f"  Maximum difference in peak = {max(dpeakrel)} (relative), {max(dpeakfrac)} (fractional) \n")
        
    if (hasconverged):
        print(f"\n -- Calibration process converged -- \n")

    return (hasconverged)
#   -----------------------------------------------------------------------------------------------------



def subandimg (targetvislist, pars=None):
    
    #   Make sub-band images 
    
    if ((pars['SubRanges']==None) or (len(pars['SubRanges']) < 3)):
        print("Not making sub-band imaging. It's a waste of time...")
        return (0)

    nbands  = len(pars['SubRanges']) - 1
    print(f"\n  Making images for {nbands} sub-bands...")

    print("******************************************************")
    print(" Assuming same frequency settings for all visibilities")
    print(" If not, good luck...")
    print("******************************************************")

    wmsmd   = casatools.msmetadata()
    wmsmd.open(targetvislist[0])
    chan_freqs  = wmsmd.chanfreqs(0)/1.0e6
    wmsmd.done()

    hrad    = float(pars['ImgSize'][0])/2
    finmask = 'Circle[['+str(hrad)+'pix, '+str(hrad)+'pix],'+str(hrad)+'pix]'

    print("\nVisibilities -- ",targetvislist)

    for i in range (0, nbands):
        imgpre  = pars['WorkDir']+pars['ImgDir']+'/'+pars['FinImage']+"_sub_"+str(pars['SubRanges'][i])+"_"+str(pars['SubRanges'][i+1])        
        print("Imagename -- ",imgpre)    
    
        print("Clearing existing image components...\n")
        os.system("rm -rf "+imgpre+".*")

        cl      = max(np.argmin(np.abs(chan_freqs - pars['SubRanges'][i])), 0) 
        cr      = min(np.argmin(np.abs(chan_freqs - pars['SubRanges'][i+1])), len(chan_freqs)) 
        chanstr = "0:"+str(min(cl,cr))+"~"+str(max(cl,cr) )

        print(f"  Imaging frequence range [{pars['SubRanges'][i]} , {pars['SubRanges'][i+1]}], channels {chanstr}")

        ct.tclean(
            vis=targetvislist, \
            imagename=imgpre, \
            datacolumn="corrected", \
            spw=chanstr, \
            imsize=pars['ImgSize'], \
            cell=pars['CellSize'], \
            selectdata=True, \
            field=pars['TargetName'],\
            uvrange=pars['FinUvLim'], \
            specmode='mfs', \
            gridder='widefield', \
            wprojplanes=pars['WprojPln'], \
            pblimit=0.1, \
            deconvolver='mtmfs', \
            scales=pars['DeconScls'], \
            smallscalebias=pars['SclBias'], \
            weighting='briggs',\
            robust=pars['ImRobust'], \
            uvtaper=pars['ImUvTaper'], \
            niter=pars['FiNiter'], \
            nsigma=pars['FinSigma'], \
            interactive=False, \
            usemask='user', \
            mask=finmask, \
            pbmask=0.2       
        )

        print("Exporting final image to Output Directory...")

        ct.exportfits(
            imagename=imgpre+".image.tt0", \
            fitsimage=pars['OutDir']+'/'+pars['FinImage']+"_sub_"+str(pars['SubRanges'][i])+"_"+str(pars['SubRanges'][i+1])+".fits", \
            overwrite=True
        )
    
    print(" Done!\n")

    return (0)
#   -----------------------------------------------------------------------------------------------------



def splitsubands (visfile, pars):
    
    #   Split individual sub-bands 

    if ((pars['ScalFreq']==None) or (len(pars['ScalFreq']) < 2)):
        print("Simply copying the visibilities...")
        os.system(f"rm -rf {visfile}_b0.ms")
        os.system(f"cp -r {visfile}.ms {visfile}_b0.ms")
        return (0)

    wmsmd   = casatools.msmetadata()
    wmsmd.open(visfile+".ms")
    chan_freqs  = wmsmd.chanfreqs(0)/1.0e6
    wmsmd.done()

    os.system(f"rm -rf {visfile}_b{0}.ms")    

    cl      = max(np.argmin(np.abs(chan_freqs - pars['ScalFreq'][0])), 0) 
    cr      = min(np.argmin(np.abs(chan_freqs - pars['ScalFreq'][1])), len(chan_freqs)) 
    chanstr = "0:"+str(min(cl,cr))+"~"+str(max(cl,cr) )

    print(f"  Self-cal frequency range {pars['ScalFreq']}, channels {chanstr}")

    ct.mstransform(
        vis=visfile+".ms", \
        outputvis=visfile+"_b0.ms", \
        datacolumn="DATA", \
        keepflags=False, \
        spw=chanstr, \
        correlation=pars['CorrProds']
    )        
    
    print(" Done!\n")

    return (0)
#   -----------------------------------------------------------------------------------------------------
