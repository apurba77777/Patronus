import os,sys
import argparse as ap
import yaml as ym
from specscripts.ingredients import *
from specscripts.auxfns import *

#	---------------------------------------------------------------------------------------------------------
#
#	Cauldron to bre potions for spectral line analysis 
#                                                   AB  [last updated: 10 July 2026]
#
#   This programme can be used for radio spectroscopy 
#
#   To run this programme, use the following command with a python executor
#
#   cauldron.py         --[potion(s)]               //  processing step(s) -- see description below
#                       --infile [param_YAML]       //  YAML file containing input parameters
#                       --pipedir [pipe_direcory]   //  Path to the pipeline itself
#
#   Muggle-friendly potions (can be run without any prior knowledge of potion making)              
#                           getdspec    //  Generate dynamic spectrum at a specific sky position
#                           initrawms   //  Initialize raw MS
#
#   Simple & convenient charms          
#                           obliviate   //  Clear existing files 
#                           lumos       //  List Usable Modes On Screen
#                           revelio     //  Reveal configuration parameters
#
#   Advanced spells and charms (Should NOT be attempted before passing O. W. L.s)
#                        
#                           accio       //  Accumulate Continuum Components in Image Output
#                           scourgify   //  Scrutinize Calibration Outputs and Ultimate Robustness of Gains with Image Files Yielded
#                           incendio    //  Image Normal Continuum Emission using Nice Data from Interferometric Observations
#  
#   Dangerous spells and curses (Extreme caution recommended!! Should NOT be attempted before passing N.E.W.T.s)
#
#                           crucio      //  Calibrate Response for an Uncorrupted Channel Isolated from Observation 
#
#	--------------------------------------------------------------------------------------------------------


#   Get command line arguments
argus   = potion_args()

print("\n-------------------------------------------------------\n")
print("               Brewing potions                    ")
print("\n-------------------------------------------------------\n")

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

#   Make data directories
phials(pars)

#   List supported modes
if (argus.lumos):  
    potion_ingredients()

#   --------------------------- Tasks   ------------


#   Search at a specific sky position

# if (argus.getdspec):  
    
#     fitslist   = [ pars['OutDir']+pars['CubeDir']+fname for fname in pars['FitsNames'] ]

#     for fitsname in fitslist:
#         getdynspec (fitsname, pars)



