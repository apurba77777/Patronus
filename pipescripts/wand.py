import os,sys
import argparse as ap
import yaml as ym
from cascripts.calfns import *
from cascripts.imgfns import *
from cascripts.timaging import *
from cascripts.utils import *

#	---------------------------------------------------------------------------------------------------------
#
#	The wand used to cast spells 
#                                                   AB  [last updated: 1 July 2026]
#
#   This programme can be used to calibrate and image (GMRT) visibility data
#
#   To run this programme, use the following command with a python executor
#
#   wand.py             --[spell(s)]                //  processing step(s) -- see description below
#                       --infile [param_YAML]       //  YAML file containing input parameters
#                       --flgin [ankflag_YAML]      //  YAML file containing aNKflag parameters
#                       --rfifile [RFI_text]        //  TEXT file with RFI frequency ranges
#                       --pipedir [pipe_direcory]   //  Path to the pipeline itself
#                       --imgname [image_name]      //  Name of the image (only for imaging)
#                       --oldimg [old_image]        //  Old image (for checking selfcal)
#                       --savemodel                 //  Save model column?
#                       --intmask                   //  Interactive masking?
#                       --calmode [p/ap]            //  Self-calibration mode (p/ap)
#
#
#   Muggle-friendly spells              fitstoms    //  Convert FITS to MS
#                                       initrawms   //  Initialize raw MS
#                                       makech0     //  Create single channel file
#                                       fluxch0     //  Set flux density of single channel file
#                                       calch0      //  Calibrate single channel file
#                                       flagch0     //  Flag single channel file
#                                       exbpcal     //  Extract bandpass calibrator file
#                                       calbpcal    //  Calibrate bandpass
#                                       flagbpcal   //  Flag bandpass calibrator file
#                                       extarget    //  Extract calibrated target file
#                                       flagtarget  //  Flag calibrated target file
#
#                                       avgtarget   //  Channel average target visibilities
#                                       imgtarget   //  Image the calibrated target
#                                       selfcal     //  Self calibrate
#                                       flagselfcal //  Flag calibrated visibilities
#                                       getuvsub    //  Subtract the final continuum model
#                                       flaguvsub   //  Flag continuum subtracted visibilities
#                                       
#                                       metronome   //  Make list of timestamps in MJD
#                                       snapshot    //  Make snapshot images/cubes
#
#   Simple & convenient charms          obliviate   //  Clear existing files 
#                                       lumos       //  List Usable Modes On Screen
#                                       revelio     //  Reveal configuration parameters
#
#   Advanced spells and charms (Should NOT be attempted before passing O. W. L.s)
#                        
#                                       accio       //  Accumulate Continuum Components in Image Output
#                                       scourgify   //  Scrutinize Calibration Outputs and Ultimate Robustness of Gains with Image Files Yielded
#                                       incendio    //  Image Normal Continuum Emission using Nice Data from Interferometric Observations
#  
#   Dangerous spells and curses ( Extreme caution recommended !! Should NOT be attempted before passing N.E.W.T.s)
#
#                                       crucio      //  Calibrate Response for an Uncorrupted Channel Isolated from Observation 
#                                       defodio     //  Determine Effects of Frequency Ousting Detected Interference in Observation
#                                       confringo   //  Calibrate Observation for Normal and Frequency Response of the Instrument with Natural Good Objects
#                                       imperio     //  Iterative Mapping of Persistent Emission in Radio using Interferometric Observations
#                                       reducto     //  Reduce entire dataset to usable calibrated target outputs
#                                       rictusempra //  Remove Image Components Through Uv Subtraction and Endeavour Mitigation of Persistent Radio Aberrations 
#                                       petrificus  //  Produce and Encapsulate Time Resolved Images into a Fits Image by Combining Unique Snapshots
#
#	--------------------------------------------------------------------------------------------------------


#   Get command line arguments
argus   = get_args()

print("\n---------------------------------------\n")
print("            Wand at ready !             ")
print("\n---------------------------------------\n")

#   Read and complain about input parameters
if (argus.infile == None):
    print(" Missing input YAML file! Please provide one...")
    sys.exit()
else:
    with open(argus.infile+'.yml', 'r') as infl:
        pars = ym.load(infl, Loader=ym.SafeLoader)
        if (argus.revelio):
            print(" Inputs provided -- \n")
            print(ym.dump(pars, sort_keys=False))


if (argus.flgin == None):
    print(" Missing YAML file for flagging! aNKflag won't run...\n")


#   Make data directories
conjure_boxes(pars)

#   List supported modes
if (argus.lumos):  
    print_spells()


#   --------------------------- Execute spells, charms and curses(!)   --------------------------------------------------

#   Convert fits to MS
if (argus.fitstoms):  
    importrawuvfile(pars['RawUvFile'], pars['RawFlagFiles'], pars, ovrt=argus.obliviate)


#   Initialize raw MS file
if (argus.initrawms):   
    initrawuvfile(pars['WorkDir']+pars['UvMsDir']+pars['ReducedName'], pars, rfifreq=argus.rfifile, ovrt=argus.obliviate)


#   Create single channel file
if (argus.makech0):  
    makesinglechan(pars['WorkDir']+pars['UvMsDir']+pars['ReducedName'], pars, ovrt=argus.obliviate)


#   Set flux scale of single channel file
if (argus.fluxch0):  
    setfluxsinglechan(pars)


#   Calibrate single channel file
if (argus.calch0):  
    calsinglechan(pars)


#   Flag single channel file
if (argus.flagch0):  
    flagsinglechan(pars, ankdir=argus.pipedir+"ankflag_3/", ankin=argus.flgin, ovrt=argus.obliviate)


#   Extract bpcal file
if (argus.exbpcal):  
    exbpcal(pars['WorkDir']+pars['UvMsDir']+pars['ReducedName'], pars['FluxCal'], pars)


#   Calibrate bandpass
if (argus.calbpcal):  
    calbpcal( pars['FluxCal'], pars)


#   Flag bandpass calibrator file
if (argus.flagbpcal):  
    flagbpcal(pars['FluxCal'], pars, ankdir=argus.pipedir+"ankflag_3/", ankin=argus.flgin, ovrt=argus.obliviate)


#   Difficult magic ********* Calibrate time dependent gain in single channel 
#   Calibrtation - flagging loop
if (argus.crucio or argus.confringo or argus.reducto):

    makesinglechan(pars['WorkDir']+pars['UvMsDir']+pars['ReducedName'], pars, ovrt=argus.obliviate)
    setfluxsinglechan(pars)
    
    for i in range(0, pars['CalIter']):
        print(f'\n Gain calibration iteration {i}...\n')
        flagsinglechan(pars, ankdir=argus.pipedir+"ankflag_3/", ankin=argus.flgin, ovrt=True)
        calsinglechan(pars)


#   Difficult magic ********* Calibrate time dependent gain in single channel
#   Calibrtation - flagging loop
if (argus.defodio or argus.confringo or argus.reducto):

    exbpcal(pars['WorkDir']+pars['UvMsDir']+pars['ReducedName'], pars['FluxCal'], pars)

    for i in range(0, pars['BpIter']):
        print(f'\n Bandpass calibration iteration {i}...\n')
        flagbpcal(pars['FluxCal'], pars, ankdir=argus.pipedir+"ankflag_3/", ankin=argus.flgin, ovrt=True)
        calbpcal( pars['FluxCal'], pars)


#   Extract calibrated target file
if (argus.extarget or argus.defodio or argus.confringo or argus.reducto):  
    extarget(pars['WorkDir']+pars['UvMsDir']+pars['ReducedName'], pars['FluxCal'], pars)


#   Flag calibrated target file
if (argus.flagtarget or argus.defodio or argus.confringo or argus.reducto): 
    flagtarget(pars['TargetName'], pars, ankdir=argus.pipedir+"ankflag_3/", ankin=argus.flgin, ovrt=argus.obliviate)


#   ------------------------    Imaging and selfcal tasks ---------------------------------------------------------

if (pars['VisList'] == None):   
    #print("\n ---- No vislist provided ----")
    vislist   = [ pars['ReducedName']+"_"+pars['TargetName'] ]
else:
    #print("\n Found list of visibilities")
    vislist   = pars['VisList']


#   Channel average target visibilities
if (argus.avgtarget or argus.imperio or argus.reducto):   
    print(f"Working with {vislist}\n")
    for ivis in vislist:
        avgtarget (ivis, pars)


#   Image the target field
if (argus.imgtarget):      
    listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_avg.ms" for vis in vislist ]
    imgtarget(listofvis, argus.imgname, argus.savemodel, argus.intmask, pars, clnmask=pars['MaskFile'])


#   Self-calibrate
if (argus.selfcal):      
    listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_avg" for vis in vislist ]
    for ivis in listofvis:
        selfcal (ivis, ivis+".scal", argus.calmode, pars)
    

#   Find sources and make a catalogue
if (argus.accio):  
    imgfile  = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+'_'+argus.imgname   
    findsrcs (imgfile, pars)


#   Check if self calibration converged
if (argus.scourgify):  
    img1  = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+'_'+argus.imgname
    img2  = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+'_'+argus.oldimg    
    cstatus  = checkselfcal (img1, img2, pars)
    if(cstatus):
        print("  Selfcal has converged !! \n")
    else:
        print("  Nope! Need to continue... \n")


#   Flag calibrated visibilities
if (argus.flagselfcal):      

    listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_avg" for vis in vislist ]
    for ivis in listofvis:
        flagcaltarget (ivis, pars, ankdir=argus.pipedir+"ankflag_3/", ankin=argus.flgin, ovrt=argus.obliviate)


#   Difficult magic ********* Self-calibrate until it converges
if (argus.imperio or argus.reducto):
    
    #   Make the zeroth image
    listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_avg.ms" for vis in vislist ]
    imgtarget(listofvis, "fscal_0", dosavemodel=False, dointeractive=argus.intmask, pars=pars, clnmask=None)

    #   Find sources and make a clean mask
    imgfile  = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+'_fscal_0'  
    findsrcs (imgfile, pars)

    atmpt   = 1
    cstatus = False
    scalmode= "p"
    if (pars['scapp']=="ap"):
        scalmode    = "ap"

    while ((atmpt <= pars['MaxIter']) and (not cstatus)):
        print(f"\n\n Self calibration: iteration {atmpt} Mode {scalmode}")

        #   Make image for this iteration
        listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_avg.ms" for vis in vislist ]
        oldmask     = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+"_fscal_"+str(atmpt-1)+"_src_mask.mask"
        
        if ( (pars['MaskFile'] != None) and (pars['MaskFile'] != "") ):
            if ( os.path.exists(pars['MaskFile'])):
                oldmask     = pars['MaskFile']
                print("Using given Mask file for cleaning...")

        imgtarget(listofvis, "fscal_"+str(atmpt), dosavemodel=True, dointeractive=argus.intmask, pars=pars, clnmask=oldmask)

        #   Find sources and make a clean mask
        imgfile  = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+"_fscal_"+str(atmpt)
        findsrcs (imgfile, pars)

        #   Check if selfcal has converged
        if (atmpt > 1):
            oldim   = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+"_fscal_"+str(atmpt-1)
            newim   = pars['WorkDir']+pars['ImgDir']+'/'+pars['TargetName']+"_fscal_"+str(atmpt)   
            cstatus = checkselfcal (oldim, newim, pars)
            if (cstatus):
                print(f"  Selfcal has converged in {scalmode} mode \n")
                if (scalmode=="p" and pars['scapp']!="p"):
                    cstatus     = False
                    scalmode    = "ap"
            else:
                print("  Nope! Need to continue... \n")

        #   Find gain solutions from self-calibration
        listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_avg" for vis in vislist ]
        for ivis in listofvis:
            selfcal (ivis, ivis+".scal", scalmode, pars)

        #   Flag calibrated visibilities
        for ivis in listofvis:
            flagcaltarget (ivis, pars, ankdir=argus.pipedir+"ankflag_3/", ankin=argus.flgin, ovrt=True)

        atmpt += 1


#   Attempt to produce the *final* continuum image
if (argus.incendio or argus.reducto):      
    listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_avg.ms" for vis in vislist ]
    finalimg(listofvis, dosavemodel=True, pars=pars)


if (pars['VisUvSub'] == None):   
    #print("\n ---- No vislist provided ----")
    visuvsublist   = [ pars['ReducedName']+"_"+pars['TargetName'] ]
else:
    #print("\n Found list of visibilities")
    visuvsublist   = pars['VisUvSub']


#   Channel average target visibilities
if (argus.getuvsub or argus.rictusempra):      
    for ivis in visuvsublist:
        getuvsub (ivis, ivis+"_avg.scal", pars)


#   Flag continuum subtracted visibilities
if (argus.flaguvsub or argus.rictusempra):      

    listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_uvsub" for vis in visuvsublist ]

    for ivis in listofvis:
        flagavguvsub (ivis, pars, ankdir=argus.pipedir+"ankflag_3/", ankin=argus.flgin, ovrt=argus.obliviate)


#   Make list of timestamps in MJD
if (argus.metronome or argus.petrificus):      
    listofvis   = [ pars['WorkDir']+pars['ImgUvDir']+vis+"_uvsub_f_avg" for vis in visuvsublist ]

    for ivis in listofvis:
        maketime (ivis)


#   Make snapshot images
if (argus.snapshot or argus.petrificus):

    listofvis   = [ vis+"_uvsub_f_avg" for vis in visuvsublist ]

    for ivis in listofvis:
        times = np.loadtxt(pars['WorkDir']+pars['ImgUvDir']+ivis+"_mjds.txt")
        
        if (pars['TavgChan'] < 1):
            print('\n Averaging over all frequencies... \n')
            timager (pars['WorkDir']+pars['ImgUvDir']+ivis, times, pars, ntime=pars['TNImg'])
            cubename    = pars['OutDir']+pars['CubeDir']+"/"+ivis+"_tcube"
        else:
            tfimager (pars['WorkDir']+pars['ImgUvDir']+ivis, times, pars, ntime=pars['TNImg'])
            cubename    = pars['OutDir']+pars['CubeDir']+"/"+ivis+"_tfcube"
        
        makefits (times, cubename, ntime=pars['TNImg'], nchan=pars['TavgChan'])
    

#   Show patronus
if (argus.expecto_patronum or argus.lumos or argus.revelio or argus.accio or argus.scourgify or argus.crucio or argus.defodio or \
    argus.confringo or argus.reducto or argus.imperio or argus.incendio or argus.rictusempra or argus.petrificus):

    patronus_charm (argus)


print("\n----------------------------------------------------------------------------")
print("----------------------------------------------------------------------------")
