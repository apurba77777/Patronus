import os,sys
import numpy as np
import argparse as ap
import yaml as ym


#   --------------------------------------------------------------------------------------------------------------------------
#   Ingredients for potion making
#
#                                              AB  [last updated: 10 July 2026] 
#
#   --------------------------------------------------------------------------------------------------------------------------


def potion_args():
    
    #   Read command-line arguments for potion making

    parser = ap.ArgumentParser(
        description = "Ingredients for potion making"
    )

    #   Input files, identifiers & parameters
    parser.add_argument("--infile", help = "YAML file with input params", type = str, default = None)
    parser.add_argument("--pipedir", help = "Directory to the pipeline", type = str, default = None)

    #   Options/utilities
    parser.add_argument("--obliviate", help = "Clear existing files ?", action='store_true')
    parser.add_argument("--lumos", help = "List usable modes on screen", action='store_true')
    parser.add_argument("--revelio", help = "Reveal configuration parameters", action='store_true')

    #   Search at a specific sky position
    parser.add_argument("--excat", help = "Extract catalogue for spectroscopic analysis", action='store_true')
    parser.add_argument("--exubcubes", help = "Extract subcubes", action='store_true')
    parser.add_argument("--smoothcubes", help = "Spatially smooth subcubes", action='store_true')
        
    args = parser.parse_args()

    return args
#   --------------------------------------------------------------------------------------------------------------------------


def potion_ingredients():

    #   Print ingredients for potion making

    print("\n*** List of potions ***\n")

    print(" Muggle-friendly potions \n")

    print("      --excat       - Extract catalogue for spectroscopic analysis")

    print("\n Simple & convenient charms \n")

    print("      --obliviate   - Clear existing files")  
    print("      --lumos       - List Usable Modes On Screen")
    print("      --revelio     - Reveal configuration parameters\n")

    print("\n Advanced potions (Do not attempt before passing O W L) \n")

    print("      --accio       - Accumulate Continuum Components in Image Output")
    
    return (0)
#   ----------------------------------------------------------------------------------------------------------------


def phials (pars):

    #   Create directories for Nifflers to store their treasures

    if (os.path.exists(pars['WorkDir'])):
        print("Found ",pars['WorkDir'])
    else:
        print("Creating ",pars['WorkDir'])
        os.mkdir(pars['WorkDir'])

    if (os.path.exists(pars['WorkDir']+pars['StacatDir'])):
        print("Found ",pars['WorkDir']+pars['StacatDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['StacatDir'])
        os.mkdir(pars['WorkDir']+pars['StacatDir'])

    if (os.path.exists(pars['WorkDir']+pars['SubDir'])):
        print("Found ",pars['WorkDir']+pars['SubDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir'])
        os.mkdir(pars['WorkDir']+pars['SubDir'])

    if (os.path.exists(pars['WorkDir']+pars['ResplotDir'])):
        print("Found ",pars['WorkDir']+pars['ResplotDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['ResplotDir'])
        os.mkdir(pars['WorkDir']+pars['ResplotDir'])    

    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['ScubeDir'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['ScubeDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['ScubeDir'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['ScubeDir'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['IscubeDir'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['IscubeDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['IscubeDir'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['IscubeDir'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['OscubeDir'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['OscubeDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['OscubeDir'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['OscubeDir'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['RedregCubes'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['RedregCubes'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['RedregCubes'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['RedregCubes'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['RedintCubes'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['RedintCubes'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['RedintCubes'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['RedintCubes'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['SpecDir'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['SpecDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['SpecDir'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['SpecDir'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['IntSpecs'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['IntSpecs'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['IntSpecs'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['IntSpecs'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['RedregSpecs'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['RedregSpecs'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['RedregSpecs'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['RedregSpecs'])

    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['RedintSpecs'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['RedintSpecs'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['RedintSpecs'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['RedintSpecs'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['NoiseDir'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['NoiseDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['NoiseDir'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['NoiseDir'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['InoiseDir'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['InoiseDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['InoiseDir'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['InoiseDir'])
    
    if (os.path.exists(pars['WorkDir']+pars['SubDir']+pars['RenoiseDir'])):
        print("Found ",pars['WorkDir']+pars['SubDir']+pars['RenoiseDir'])
    else:
        print("Creating ",pars['WorkDir']+pars['SubDir']+pars['RenoiseDir'])
        os.mkdir(pars['WorkDir']+pars['SubDir']+pars['RenoiseDir'])

    return (0)
#   ----------------------------------------------------------------------------------------------------------------

