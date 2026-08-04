#!/bin/bash

#   ***********************************************************
#   
#   Conjure spells to calibrate and image visibility data
#
#   For basic instrctions, run with 
#                               --help 
#
#                                           AB 
#                                           Last updated 23/6/26
#   
#   ***********************************************************

#   --------------------------------------------------------
#   Default paths and input files
confile=""
pyex="python3"
Pipe_Dir="../../patronus/"
parfile="../inputfiles/calibration_b4_exam"
flgpar="../inputfiles/ankpar_b4"
rfifile="../inputfiles/rfifreqb4.txt"
searchpars="../inputfiles/search_b4"
logfile="patronus.log"

imgext=""
imgext0=""
savemodel=""
intmask=""
scalmode=""
ovr=""
sln=""
spl=""

#   --------------------------------------------------------                                                    
howtoconjure()
{
   # Print Instructions on terminal 
   echo
   echo "  Instructions to conjure "
   echo
   echo "options:"
   echo "--help             Print these options"
   echo "-c / --conf [C]    Configuration file (without .con)"
   echo "-s / --spell [S]   Spell"
   echo "-i / --imnew [I]   Name extension of the (new) image"
   echo "-j / --imold [I]   Name extension of the old image"
   echo "-p / --pcal        Phase-only (selfcal)"
   echo "--save             Save model column (imaging)"
   echo "--interactive      Interactive masking (imaging)"
   echo "--obliviate        Clear existing files if required"
   echo "--flagrate         Write terminal messages to logfile"
   echo 
}

#   ---------------------------------------------------------

# Process the input options       
for arg in "$@"; do
    case $arg in
        --help)
            howtoconjure
            exit;;
        -c|--conf)
            confile=$2
            source ${confile}.con
            echo "Using configuration from ${confile}.con"
            shift 2;;
        -s|--spell)
            spl="--$2"
            echo "Casting spell ${spl}"
            shift 2;;
        -i|--imnew)
            imgext="--imgname $2"
            echo "New image ${imgext}"
            shift 2;;
        -j|--imold)
            imgext0="--oldimg $2"
            echo "Old image ${imgext0}"
            shift 2;;
        -p|--pcal)
            echo "Phase-only self-cal chosen"
            scalmode="--calmode p"
            shift;;
        --save)
            echo "Saving model column"
            savemodel="--savemodel"
            shift;;
        --interactive)
            echo "Interactive masking chosen"
            intmask="--intmask"
            shift;;
        --obliviate)
            echo "Existing knowledge will be forgotten..."
            ovr="--obliviate"
            shift;;
        --flagrate)
            echo "Messages will be written to ${logfile}"
            sln="dumb"
            shift;;
    esac
done

# echo "${pyex}"
# echo "${Pipe_Dir}"
# echo "${parfile}"
# echo "${flgpar}"
# echo "${rfifile}"
# echo "${searchpars}"

#   --------------------------------------------------------

#   Construct commands

ankimgcmd="${pyex} -u ${Pipe_Dir}/pipescripts/wand.py \
    --pipedir ${Pipe_Dir} \
    --infile ${parfile} \
    --flgin ${flgpar} \
    --rfifile ${rfifile} \
    --pypath ${pyex}
    ${imgext} \
    ${imgext0} \
    ${savemodel} ${intmask} ${scalmode} \
    ${spl} ${ovr}"


export PYTHONPATH=${PYTHONPATH}:${Pipe_Dir} 

export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

echo
echo "Running command"    
echo ${ankimgcmd} 
if [ "$sln" != "" ]
then
    ${ankimgcmd} |& tee -a ${logfile}
else
    ${ankimgcmd}
fi

echo




