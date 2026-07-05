#!/bin/bash

#   ***********************************************************
#   
#   Cast spells to find transients
#
#   For basic instrctions, run with 
#                               --help 
#
#                                           AB 
#                                           Last updated 1/7/26
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

incmd="${pyex} -u ${Pipe_Dir}/pipescripts/incantation.py \
    --pipedir ${Pipe_Dir} \
    --infile ${searchpars} \
    ${spl} ${ovr}"


export PYTHONPATH=${PYTHONPATH}:${Pipe_Dir} 

export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

echo
echo "Running command"    
echo ${incmd} 
if [ "$sln" != "" ]
then
    ${incmd} |& tee -a ${logfile}
else
    ${incmd}
fi

echo




