import os,sys
import socket
import platform
import glob
import numpy as np
import argparse as ap
import yaml as ym


#   --------------------------------------------------------------------------------------------------------------------------
#   Utility functions for wand.py
#
#                                              AB  [last updated: 24 June 2026] 
#
#   --------------------------------------------------------------------------------------------------------------------------

def get_args():
    
    #   Read command-line arguments

    parser = ap.ArgumentParser(
        description = "Wand for casting spells and charms"
    )

    #   Input files, identifiers & parameters
    parser.add_argument("--infile", help = "YAML file with input params", type = str, default = None)
    parser.add_argument("--flgin", help = "YAML file with flagging params", type = str, default = None)
    parser.add_argument("--rfifile", help = "File with list of RFI frequencies", type = str, default = None)
    parser.add_argument("--pipedir", help = "Directory to the pipeline", type = str, default = None)
    parser.add_argument("--pypath", help = "Path to python executable", type = str, default = "python")
    parser.add_argument("--imgname", help = "Name of the image (only for imaging)", type = str, default = "random")
    parser.add_argument("--oldimg", help = "Old image (for checking selfcal)", type = str, default = "random")
    parser.add_argument("--savemodel", help = "Save model column?", action='store_true')
    parser.add_argument("--intmask", help = "Interactive masking?", action='store_true')
    parser.add_argument("--calmode", help = "Self-calibration mode (p/ap)", type = str, default = "p")

    #   Calibration tasks
    parser.add_argument("--fitstoms", help = "Convert FITS to MS", action='store_true')
    parser.add_argument("--initrawms", help = "Initialize raw MS", action='store_true')
    parser.add_argument("--makech0", help = "Create single channel file", action='store_true')
    parser.add_argument("--fluxch0", help = "Set flux density of single channel file", action='store_true')
    parser.add_argument("--calch0", help = "Calibrate single channel file", action='store_true')
    parser.add_argument("--flagch0", help = "Flag single channel file", action='store_true')
    parser.add_argument("--exbpcal", help = "Extract bandpass calibrator file", action='store_true')
    parser.add_argument("--calbpcal", help = "Calibrate bandpass", action='store_true')
    parser.add_argument("--flagbpcal", help = "Flag bandpass calibrator file", action='store_true')
    parser.add_argument("--phcalbpcal", help = "Calibrate bandpass with phase cal", action='store_true')
    parser.add_argument("--phflagbpcal", help = "Flag phase cal bandpass calibrator file", action='store_true')
    parser.add_argument("--modpcalspec", help = "Model the phase cal spectrum", action='store_true')
    parser.add_argument("--extarget", help = "Extract calibrated target file", action='store_true')
    parser.add_argument("--flagtarget", help = "Flag calibrated target file", action='store_true')
    #   Imaging tasks
    parser.add_argument("--avgtarget", help = "Channel average target visibilities", action='store_true')
    parser.add_argument("--imgtarget", help = "Image the calibrated target", action='store_true')
    parser.add_argument("--selfcal", help = "Self calibrate", action='store_true')    
    parser.add_argument("--flagselfcal", help = "Flag self-calibrated visibilities", action='store_true')
    parser.add_argument("--getuvsub", help = "Subtract the final continuum model", action='store_true')
    parser.add_argument("--flaguvsub", help = "Flag continuum subtracted visibilities", action='store_true')
    #   Snapshot imaging tasks
    parser.add_argument("--metronome", help = "Make list of timestamps in MJD", action='store_true')
    parser.add_argument("--snapshot", help = "Make snapshot images/cubes", action='store_true')    

    #   Options/utilities
    parser.add_argument("--obliviate", help = "Clear existing files ?", action='store_true')
    parser.add_argument("--lumos", help = "List usable modes on screen", action='store_true')
    parser.add_argument("--revelio", help = "Reveal configuration parameters", action='store_true')
    
    parser.add_argument("--accio", help = "Accumulate Continuum Components in Image Output", action='store_true')
    parser.add_argument("--scourgify", help = "Scrutinize Calibration Outputs and Ultimate Robustness of Gains with Image Files Yielded", action='store_true')
    parser.add_argument("--incendio", help = "Image Normal Continuum Emission using Nice Data from Interferometric Observations", action='store_true')

    #   Combined tasks
    parser.add_argument("--crucio", help = "Calibrate Response for an Uncorrupted Channel Isolated from Observation", action='store_true')
    parser.add_argument("--defodio", help = "Determine Effects of Frequency Ousting Detected Interference in Observation", action='store_true')
    parser.add_argument("--confringo", help = "Calibrate Observation for Normal and Frequency Response of the Instrument with Natural Good Objects", action='store_true')
    parser.add_argument("--imperio", help = "Iterative Mapping of Persistent Emission in Radio using Interferometric Observations", action='store_true')
    parser.add_argument("--reducto", help = "Reduce Entire Dataset to Useful Calibrated Target Outputs", action='store_true')
    parser.add_argument("--rictusempra", help = "Remove Image Components Through Uv Subtraction and Endeavour Mitigation of Persistent Radio Altercations", action='store_true')
    parser.add_argument("--petrificus", help = "Produce and Encapsulate Time Resolved Images into a Fits Image by Combining Unique Snapshots", action='store_true')

    parser.add_argument("--expecto_patronum", help = "What it says...", action='store_true')

    args = parser.parse_args()

    return args
#   --------------------------------------------------------------------------------------------------------------------------


def print_spells():

    #   Print charms and spells

    print("\n*** List of spells ***\n")

    print(" Muggle-friendly spells \n")

    print("      --fitstoms    - Convert FITS to MS")
    print("      --initrawms   - Initialize raw MS")
    print("      --makech0     - Create single channel file")
    print("      --flagch0     - Flag single channel file")
    print("      --fluxch0     - Set flux density of single channel file")
    print("      --calch0      - Calibrate single channel file")
    print("      --exbpcal     - Extract bandpass calibrator file")
    print("      --flagbpcal   - Flag bandpass calibrator file")
    print("      --calbpcal    - Calibrate bandpass")
    print("      --extarget    - Extract calibrated target file")
    print("      --flagtarget  - Flag calibrated target file\n")
    
    print("      --avgtarget   - Channel average target visibilities")
    print("      --imgtarget   - Image the calibrated target")
    print("      --selfcal     - Self calibrate")
    print("      --flagselfcal - Flag self-calibrated visibilities")
    print("      --getuvsub    - Subtract the final continuum model")
    print("      --flaguvsub   - Flag continuum subtracted visibilities\n")    

    print("      --metronome   - Make list of timestamps in MJD")
    print("      --snapshot    - Make snapshot images/cubes")

    print("\n Simple & convenient charms \n")

    print("      --obliviate   - Clear existing files")  
    print("      --lumos       - List Usable Modes On Screen")
    print("      --revelio     - Reveal configuration parameters\n")

    print("\n Advanced spells and charms (Do not attempt before passing O W L) \n")

    print("      --accio       - Accumulate Continuum Components in Image Output")
    print("      --scourgify   - Scrutinize Calibration Outputs and Ultimate Robustness of Gains with Image Files Yielded")    
    print("      --incendio    - Image Normal Continuum Emission using Nice Data from Interferometric Observations\n")  

    print("\n Dangeroous spells and curses ( Extreme caution recommended !!) \n")  

    print("      --crucio      - Calibrate Response for an Uncorrupted Channel Isolated from Observation")
    print("      --defodio     - Determine Effects of Frequency Ousting Detected Interference in Observation")
    print("      --confringo   - Calibrate Observation for Normal and Frequency Response of the Instrument with Natural Good Objects")  
    print("      --imperio     - Iterative Mapping of Persistent Emission in Radio using Interferometric Observations")
    print("      --reducto     - Reduce Entire Dataset to Useful Calibrated Target Outputs")
    print("      --rictusempra - Remove Image Components Through Uv Subtraction and Endeavour Mitigation of Persistent Radio Altercations")
    print("      --petrificus  - Produce and Encapsulate Time Resolved Images into a Fits Image by Combining Unique Snapshots\n")
    
    return (0)
#   --------------------------------------------------------------------------------------------------------------------------


def conjure_boxes(pars):

    #   Conjure boxes to keep intermediate data products

    if (os.path.exists(pars['WorkDir'])):
        print("Found ",pars['WorkDir'])
    else:
        print("Creating ",pars['WorkDir'])
        os.system("mkdir "+pars['WorkDir'])

    if (os.path.exists(pars['WorkDir']+pars['UvMsDir'])):
        print("Found ",pars['WorkDir']+pars['UvMsDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['UvMsDir'])
        os.system("mkdir "+pars['WorkDir']+pars['UvMsDir'])

    if (os.path.exists(pars['WorkDir']+pars['LogDir'])):
        print("Found ",pars['WorkDir']+pars['LogDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['LogDir'])
        os.system("mkdir "+pars['WorkDir']+pars['LogDir'])

    if (os.path.exists(pars['WorkDir']+pars['ImgUvDir'])):
        print("Found ",pars['WorkDir']+pars['ImgUvDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['ImgUvDir'])
        os.system("mkdir "+pars['WorkDir']+pars['ImgUvDir'])

    if (os.path.exists(pars['WorkDir']+pars['ImgDir'])):
        print("Found ",pars['WorkDir']+pars['ImgDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['ImgDir'])
        os.system("mkdir "+pars['WorkDir']+pars['ImgDir'])
    
    if (os.path.exists(pars['OutDir'])):
        print("Found ",pars['OutDir'])
    else:
        print("Creating ",pars['OutDir'])
        os.system("mkdir "+pars['OutDir'])
    
    if (os.path.exists(pars['OutDir']+pars['CubeDir'])):
        print("Found ",pars['OutDir']+pars['CubeDir'])
    else:
        print("Creating ",pars['OutDir']+pars['CubeDir'])
        os.system("mkdir "+pars['OutDir']+pars['CubeDir'])

    return (0)
#   --------------------------------------------------------------------------------------------------------------------------


def patronus_charm (args):

    print("\n ***************************************************** ")
    print(" **                                                 ** ")
    print(" **             Conjuring PATRONUS                  ** ")
    print(" **                                                 ** ")
    print(" ***************************************************** ")

    patid   = "random"
    caster  = os.getlogin().lower()
    wands   = [os.uname().sysname.lower(), platform.uname().system.lower(), platform.system().lower()] 
    spls    = [os.uname().nodename.lower(), socket.gethostname().lower(), platform.uname().node.lower()]

    ptype   = ( (sum(ord(c) for c in caster) + sum(ord(c) for c in spls[0])) % 9 ) + 1

    if ( ("darwin" in ".".join(wands)) or ("mac" in ".".join(wands)) or ("mac" in ".".join(spls)) ):
        patid   = f"{args.pipedir}/patlib/ipat.txt"

    elif ( ("jit" in caster) or ("sur" in caster) or ("dal" in caster) or ("jit" in ".".join(spls)) or ("sur" in ".".join(spls)) ):
        patid   = f"{args.pipedir}/patlib/c_1_0.txt"
    
    elif ( ("rna" in caster) or ("das" in caster) or ("budo" in caster) or ("rna" in ".".join(spls)) or ("das" in ".".join(spls)) ):
        patlist = glob.glob(f"{args.pipedir}/patlib/f_1_*.txt")
        patid   = patlist[np.random.randint(0,high=len(patlist))]

    elif ("charizard" in ".".join(spls)):
        patlist = glob.glob(f"{args.pipedir}/patlib/p_1_*.txt")
        patid   = patlist[np.random.randint(0,high=len(patlist))]

    elif ("garfield" in ".".join(spls)):
            patlist = glob.glob(f"{args.pipedir}/patlib/p_7_*.txt")
            patid   = patlist[np.random.randint(0,high=len(patlist))]
    
    elif ("paradox" in ".".join(spls)):
            patlist = glob.glob(f"{args.pipedir}/patlib/p_1_*.txt")
            patid   = patlist[np.random.randint(0,high=len(patlist))]
    
    else:
        patlist = glob.glob(f"{args.pipedir}/patlib/p_{ptype}_*.txt")
        patid   = patlist[np.random.randint(0,high=len(patlist))]

    if (os.path.exists(patid)):
        os.system(f"cat {patid}")
    else:
        print("\n\n                Wand broken. Couldn't conjure PATRONUS.\n\n")

    return (0)
#   --------------------------------------------------------------------------------------------------------------------------


def incan_args():
    
    #   Read command-line arguments for incantations

    parser = ap.ArgumentParser(
        description = "Incantations for casting spells and charms"
    )

    #   Input files, identifiers & parameters
    parser.add_argument("--infile", help = "YAML file with input params", type = str, default = None)
    parser.add_argument("--pipedir", help = "Directory to the pipeline", type = str, default = None)

    #   Options/utilities
    parser.add_argument("--obliviate", help = "Clear existing files ?", action='store_true')
    parser.add_argument("--lumos", help = "List usable modes on screen", action='store_true')
    parser.add_argument("--revelio", help = "Reveal configuration parameters", action='store_true')

    #   Search at a specific sky position
    parser.add_argument("--getdspec", help = "Generate dynamic spectrum at a specific sky position", action='store_true')
    parser.add_argument("--mapnoise", help = "Generate spatial map of noise", action='store_true')
    parser.add_argument("--cleansweep", help = "Clean and search for transients at the original time resolution", action='store_true')
    parser.add_argument("--acleansweep", help = "Clean and search for transients after time averaging", action='store_true')


    parser.add_argument("--expecto_patronum", help = "What it says...", action='store_true')
    
    args = parser.parse_args()

    return args
#   --------------------------------------------------------------------------------------------------------------------------


def incan_spells():

    #   Print charms and spells for incantations

    print("\n*** List of spells ***\n")

    print(" Muggle-friendly spells \n")

    print("      --getdspec    - Generate dynamic spectrum at a specific sky position")

    print("\n Simple & convenient charms \n")

    print("      --obliviate   - Clear existing files")  
    print("      --lumos       - List Usable Modes On Screen")
    print("      --revelio     - Reveal configuration parameters\n")

    print("\n Advanced spells and charms (Do not attempt before passing O W L) \n")

    print("      --accio       - Accumulate Continuum Components in Image Output")
    
    return (0)
#   --------------------------------------------------------------------------------------------------------------------------


def niffler (pars):

    #   Create directories for Nifflers to store their treasures

    if (os.path.exists(pars['OutDir'])):
        print("Found ",pars['OutDir'])
    else:
        print("Creating ",pars['OutDir'])
        os.system("mkdir "+pars['OutDir'])

    if (os.path.exists(pars['OutDir']+pars['CubeDir'])):
        print("Found ",pars['OutDir']+pars['CubeDir'])
    else:
        print("Creating ",pars['OutDir']+pars['CubeDir'])
        os.system("mkdir "+pars['OutDir']+pars['CubeDir'])

    if (os.path.exists(pars['OutDir']+pars['PlotDir'])):
        print("Found ",pars['OutDir']+pars['PlotDir'])
    else:
        print("Creating ",pars['OutDir']+pars['PlotDir'])
        os.system("mkdir "+pars['OutDir']+pars['PlotDir'])

    if (os.path.exists(pars['OutDir']+pars['CanDir'])):
        print("Found ",pars['OutDir']+pars['CanDir'])
    else:
        print("Creating ",pars['OutDir']+pars['CanDir'])
        os.system("mkdir "+pars['OutDir']+pars['CanDir'])    

    return (0)
#   --------------------------------------------------------------------------------------------------------------------------

