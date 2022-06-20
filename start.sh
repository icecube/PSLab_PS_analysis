#!/bin/bash

## INSTRUCTIONS: CHANGE USER VARIABLES BELOW TO POINT TO YOUR OWN PSLAB DIRECTORY

####  USER VARIABLES  ####
_LAB_MAIN_DIR=<full_path_to_your_psLab_folder>
if [ ! -d "$_LAB_MAIN_DIR" ]
then
    printf "ERROR: the path to psLab does not exist. Please, correct the path and run the script again\n"
    exit
fi

####  MAIN SCRIPT HERE  ####
_LAB_LIB_DIR=$_LAB_MAIN_DIR/lib64
ENV=$SHELL
_LAB_CORE_DIR=$_LAB_MAIN_DIR/labcore
_LAB_PRELOAD=$_LAB_CORE_DIR/preload.C

if [[ -z "$LAB_MAIN_DIR" ]]
    then
    
    printf "\nWelcome to the lab.\n"
    printf "   LAB_MAIN_DIR = %s\n" $_LAB_MAIN_DIR
    printf "   LAB_LIB_DIR  = %s\n" $_LAB_LIB_DIR
    printf "   LAB_CORE_DIR = %s\n" $_LAB_CORE_DIR
    printf "   LAB_PRELOAD  = %s\n" $_LAB_PRELOAD
    printf "   ENV  = %s\n" $ENV
    printf "\n" 

    LAB_MAIN_DIR=$_LAB_MAIN_DIR \
    LAB_LIB_DIR=$_LAB_LIB_DIR \
    LAB_CORE_DIR=$_LAB_CORE_DIR \
    LAB_PRELOAD=$_LAB_PRELOAD \
    $ENV $@

    err=$?
    if [ $err -ne 0 ]
    then
        echo "PYTHON EXECUTABLE FAILED"
        printf "You are exiting the lab.\n"
        exit $err
    fi

    printf "You are exiting the lab.\n"
    exit
else
    printf "You are alread in the lab.  You must exit before reloading.\n"
fi

