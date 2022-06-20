#!/bin/bash

## INSTRUCTIONS: 
## 1. COPY THIS FILE TO YOUR LAB DIRECTORY:

# cp  lab/labcore/start_lab.sh  lab/start.sh

## 2. CHANGE USER VARIABLES BELOW TO POINT TO YOUR OWN DIRECTORIES


####  USER VARIABLES  ####


_LAB_MAIN_DIR=/net/user/mfbaker/testlab/branches/multi_time

_LAB_LIB_DIR="/net/user/mfbaker/testlab/branches/multi_time/lib/${OS_ARCH}"
# _LAB_LIB_DIR="/net/user/mfbaker/proto_multi/lib/${OS_ARCH}"
# To build the code from the same source directory on multiple machines
# (e.g. 32- and 64-bit cobalt machines), run the script on a different 
# machine with LIB pointing to a new directory

#ICETRAY_ENV="/net/user/mfbaker/std-processing/releases/${stdproc_version}/${OS_ARCH}/env-shell.sh"
ICETRAY_ENV="/net/user/mfbaker/offline-software/releases/V11-01-00/${OS_ARCH}/env-shell.sh"
#ICETRAY_ENV="/net/user/mfbaker/std-processing/releases/10-01-00/RHEL_4.0_amd64/env-shell.sh"
# Or, if no icetray environment is needed:
# ICETRAY_ENV=$SHELL

_LAB_DATA_DIR=/net/user/mfbaker/testlab/branches/multi_time/data
_LAB_RESULTS_DIR=/net/user/mfbaker/testlab/branches/multi_time/results
# These data and result Dirs are not essential to source code, 
# but may help to improve portability of scripts


####  MAIN SCRIPT HERE  ####


_LAB_CORE_DIR=$_LAB_MAIN_DIR/labcore
_LAB_PRELOAD=$_LAB_CORE_DIR/preload.C
# these enviroment variables usually do not need to be further adjusted


if [[ -z "$LAB_MAIN_DIR" ]]
    then
    
    printf "\nWelcome to the lab.\n"
    printf "   LAB_MAIN_DIR = %s\n" $_LAB_MAIN_DIR
    printf "   LAB_LIB_DIR  = %s\n" $_LAB_LIB_DIR
    printf "   LAB_CORE_DIR = %s\n" $_LAB_CORE_DIR
    printf "   LAB_PRELOAD  = %s\n" $_LAB_PRELOAD
    printf "   LAB_DATA_DIR = %s\n" $_LAB_DATA_DIR
    printf "   ICETRAY_ENV  = %s\n" $ICETRAY_ENV
    printf "\n"

    LAB_MAIN_DIR=$_LAB_MAIN_DIR \
	LAB_LIB_DIR=$_LAB_LIB_DIR \
	LAB_CORE_DIR=$_LAB_CORE_DIR \
	LAB_PRELOAD=$_LAB_PRELOAD \
        LAB_DATA_DIR=$_LAB_DATA_DIR \
	$ICETRAY_ENV $@ 
          # $@ allows passing arguments to the new shell

    printf "You are exiting the lab.\n"
    exit
else
    printf "You are alread in the lab.  You must exit before reloading.\n"
fi
