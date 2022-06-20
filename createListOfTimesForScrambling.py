#!/usr/bin/env python

"""
This script creates the list of times used for time scrambling in time-dep analyses and gives basic information about the files, such as start time, stop time, livetime. It takes as input good run lists in ROOT format.
The files have to contain the following columns:

start    - float64 - start MJD in days
stop     - float64 - stop MJD in days
livetime - float64 - livetime of run in days

"""

import numpy as np
import pandas as pd
import os
import random
import argparse
import ROOT as root

def tree2dict(tree):
    tstart, tstop, tlive = [], [], []
    for entry in tree:
        tstart.append(entry.start)
        tstop.append(entry.stop)
        tlive.append(entry.livetime)
    dictionary = {'start': tstart, 'stop': tstop, 'livetime': tlive}
    return dictionary

def getInfosGRL(GRL, mask=None):
    
    if mask!=None:
        start = GRL['start'][mask]
        stop = GRL['stop'][mask]
        livetime = GRL['livetime'][mask]
    
    else:        
        start = GRL['start']
        stop = GRL['stop']
        livetime = GRL['livetime']
    
    print('tmin [MJD]:', min(start), ', tmax [MJD]:', max(stop), ', total livetime [days]:', np.sum(livetime))
    
def createListOfTimes(listGRL, scaleNumberTimesPerRun, random_seed, outputFile):
    
    random.seed(random_seed)
   
    frame = pd.DataFrame()
    list_ = []
    
    #create a single dataframe with all the good run lists we want to use
    for file_ in list_GRL:
        df = pd.DataFrame(data=file_)
        list_.append(df)
        
    frame = pd.concat(list_)
    frame = frame.reset_index(drop=True) #reset indices of df so each entry has a single index
    
    start = frame['start']
    stop = frame['stop']
    livetime = frame['livetime']

    times_array = np.array([])

    for i in range (len(start)):
        #sampling with random uniform distribution, scaling number of times according to the livetime of each run
        #i.e. if run is longer, we will draw more times from this run than from a run which would be shorter
        times = np.random.uniform(start[i], stop[i],int(scaleNumberTimesPerRun*livetime[i]/np.amax(livetime)))
        times_array = np.append(times_array, times)

    np.set_printoptions(suppress=True)
    np.savetxt(outputFile, times_array, fmt='%.18f')  
    print('Saved list of', len(times_array), 'times in file', outputFile)

def generateInputFileName(season):
    if season=='IC40' or season=='IC59' or season=='IC79': return [ season+'_exp.root' ]
    elif season=='IC86_I': return [ 'IC86_2011_exp.root' ]
    elif season=='IC86_II': return [ 'IC86_2012_exp.root' ]
    elif season=='IC86_III': return [ 'IC86_2013_exp.root' ]
    elif season=='IC86_IV': return [ 'IC86_2014_exp.root' ]
    elif season=='IC86_V': return [ 'IC86_2015_exp.root' ]
    elif season=='IC86_VI': return [ 'IC86_2016_exp.root' ]
    elif season=='IC86_VII': return [ 'IC86_2017_exp.root' ]
    elif season=='IC86_II_VII': return ['IC86_2012_exp.root', 'IC86_2013_exp.root', 'IC86_2014_exp.root', 'IC86_2015_exp.root', 'IC86_2016_exp.root', 'IC86_2017_exp.root']

def generateOutputFileName(season):
    if season=='IC40' or season=='IC59' or season=='IC79':
        return 'HugeListOfTimes_'+season+'_GRL_version-003-p02.txt'
    elif season=='IC86_I': return 'HugeListOfTimes_IC86-I_GRL_version-003-p02.txt'
    elif season=='IC86_II': return 'HugeListOfTimes_IC86-II_GRL_version-003-p02.txt'
    elif season=='IC86_III': return 'HugeListOfTimes_IC86-III_GRL_version-003-p02.txt'
    elif season=='IC86_IV': return 'HugeListOfTimes_IC86-IV_GRL_version-003-p02.txt'
    elif season=='IC86_V': return 'HugeListOfTimes_IC86-V_GRL_version-003-p02.txt'
    elif season=='IC86_VI': return 'HugeListOfTimes_IC86-VI_GRL_version-003-p02.txt'
    elif season=='IC86_VII': return 'HugeListOfTimes_IC86-VII_GRL_version-003-p02.txt'
    elif season=='IC86_II_VII': return 'HugeListOfTimes_IC86-II-III-IV-V-VI-VII_GRL_version-003-p02.txt'

if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description='Gives infos and create a list of times taken from the good run list of one or many data taking periods.')
    parser.add_argument('-seed', required=True,  type=int, help='Seed for random uniform algorithm.')
    parser.add_argument('-season', required=False, type=str, help='The name of the IC season.')
    parser.add_argument('-scale', required=False, default=1000,  type=int, help='The scaling factor for extracting more or less times from each good run. IMPORTANT: you need to get a lot more times than entries in analysis data.')
    args = parser.parse_args()
    
    scaleNumberTimesPerRun = args.__dict__['scale']
    random_seed = args.__dict__['seed']

    if args.__dict__['season']: seasons = [args.__dict__['season']]
    else: seasons = [ 'IC40', 'IC59', 'IC79', 'IC86_I', 'IC86_II', 'IC86_III', 'IC86_IV', 'IC86_V', 'IC86_VI', 'IC86_VII', 'IC86_II_VII' ] 

    for season in seasons:
        listOfFiles = generateInputFileName(season)
        outputFile  = generateOutputFileName(season)
 
        list_GRL = []
        inputDir = 'data/GRL/'
        outputDir = 'ListOfTimes/'
        if not os.path.exists(outputDir):
            os.makedirs(outputDir) 

        for file in listOfFiles:
            f = root.TFile(inputDir+file)
            tree = f.Get("tree")
            GRL_file = tree2dict(tree)

            list_GRL.append(GRL_file)

            print('')
            print('Infos for file', file, ':')
            getInfosGRL(GRL_file)
            print('')
        
        print('Creating list of times with good run lists...')
        createListOfTimes(list_GRL, scaleNumberTimesPerRun, random_seed, outputDir+outputFile)
        print('')
