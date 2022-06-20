#!/bin/bash

## INSTRUCTIONS: 
## 1. COPY THIS FILE TO YOUR LAB DIRECTORY:

# cp  lab/labcore/start_lab.sh  lab/start.sh

## 2. CHANGE USER VARIABLES BELOW TO POINT TO YOUR OWN DIRECTORIES

## 3. CHOOSE ADDITIONAL PROJECTS (E.G. ASTRO, IF INSTALLED)

####  USER VARIABLES  ####


_LAB_MAIN_DIR=/home/cfinley/lab

_LAB_LIB_DIR=$_LAB_MAIN_DIR/lib64
# To build the code from the same source directory on multiple machines
# (e.g. 32- and 64-bit cobalt machines), run the script on a different 
# machine with LIB pointing to a new directory

ICETRAY_ENV=/home/cfinley/offline-software/releases/V10-03-00/build/env-shell.sh
# Or, if no icetray environment is needed:
# ICETRAY_ENV=$SHELL

# If you switch coord. projects, best to remove your LAB_LIB_DIR and rebuild!
_LAB_COORDINATE_PROJECT=coord_interface
#_LAB_COORDINATE_PROJECT=astro_interface

_LAB_DATA_LOCATION=Madison
# 'location' directory in root macro areas should contain a file with this name,
# identifying data directories on local system.
# With this, the primary scripts do not have to be modified for different users

_LAB_DATA_DIR=/home/cfinley/data
_LAB_RESULTS_DIR=/home/cfinley/results
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
    printf "   LAB_DATA_LOCATION = %s\n" $_LAB_DATA_LOCATION
    printf "   LAB_DATA_DIR = %s\n" $_LAB_DATA_DIR
    printf "   ICETRAY_ENV  = %s\n" $ICETRAY_ENV
    printf "   LAB_COORDINATE_PROJECT = %s\n" $_LAB_COORDINATE_PROJECT
    printf "\n"

    LAB_MAIN_DIR=$_LAB_MAIN_DIR \
	LAB_LIB_DIR=$_LAB_LIB_DIR \
	LAB_CORE_DIR=$_LAB_CORE_DIR \
	LAB_PRELOAD=$_LAB_PRELOAD \
	LAB_DATA_LOCATION=$_LAB_DATA_LOCATION \
        LAB_DATA_DIR=$_LAB_DATA_DIR \
	LAB_COORDINATE_PROJECT=$_LAB_COORDINATE_PROJECT \
	$ICETRAY_ENV $@ 
          # $@ allows passing arguments to the new shell

    printf "You are exiting the lab.\n"
    exit
else
    printf "You are alread in the lab.  You must exit before reloading.\n"
fi
