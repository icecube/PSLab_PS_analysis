import numpy as np
import pandas as pd
import ROOT as root
from array import array
import argparse
import os

def season2treename(year):
    treename = ''
    if   year=='IC40': treename = 'IC40_exp.root'
    elif year =='IC59': treename = 'IC59_exp.root'
    elif year =='IC79': treename = 'IC79_exp.root'
    elif year =='IC86_I': treename = 'IC86_2011_exp.root'
    elif year =='IC86_II': treename = 'IC86_2012_exp.root'
    elif year =='IC86_III': treename = 'IC86_2013_exp.root'
    elif year =='IC86_IV': treename = 'IC86_2014_exp.root'
    elif year =='IC86_V': treename = 'IC86_2015_exp.root'
    elif year =='IC86_VI': treename = 'IC86_2016_exp.root'
    elif year =='IC86_VII': treename = 'IC86_2017_exp.root'
    return treename

def create_roottree_datafile(data, season):

    tree = root.TTree("tree", "tree")

    time   = array('d',[0])
    ra     = array('d',[0])
    dec    = array('d',[0])
    sinDec = array('d',[0])
    azi    = array('d',[0])
    zen    = array('d',[0])
    angErr = array('d',[0])
    logE   = array('d',[0])
    run    = array('i',[0])
    event  = array('i',[0])

    tree.Branch("time",  time,  'time/D')
    tree.Branch("ra",  ra,  'ra/D')
    tree.Branch("dec",  dec,  'dec/D')
    tree.Branch("sinDec",  sinDec,  'sinDec/D')
    tree.Branch("azi",  azi,  'azi/D')
    tree.Branch("zen",  zen,  'zen/D')
    tree.Branch("angErr",  angErr,  'angErr/D')
    tree.Branch("logE",  logE,  'logE/D')
    tree.Branch("run", run, 'run/I')
    tree.Branch("event", event, 'event/I')

    for i in range(len(data)):
        time[0]   = data['MJD[days]'][i]
        ra[0]     = np.radians(data['RA[deg]'][i])
        dec[0]    = np.radians(data['Dec[deg]'][i])
        sinDec[0] = np.sin(np.radians(data['Dec[deg]'][i]))
        azi[0]    = np.radians(data['Azimuth[deg]'][i])
        zen[0]    = np.radians(data['Zenith[deg]'][i])
        angErr[0] = np.radians(data['AngErr[deg]'][i])
        logE[0]   = data['log10(E/GeV)'][i]
        run[0]    = 0
        event[0]  = i
        tree.Fill() 

    return tree 

def create_roottree_GRLfile(data, season):

    tree = root.TTree("tree", "tree")

    start     = array('d',[0])
    stop      = array('d',[0])
    livetime  = array('d',[0])

    tree.Branch("start",  start,  'start/D')
    tree.Branch("stop",  stop,  'stop/D')
    tree.Branch("livetime",  livetime,  'livetime/D')

    for i in range(len(data)):
        start[0]    = data['MJD_start[days]'][i]
        stop[0]     = data['MJD_stop[days]'][i]
        livetime[0] = stop[0] - start[0]
        tree.Fill() 

    return tree 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run psLab',prefix_chars='@')
    parser.add_argument('@inDir', required=True, type=str, help='directory containing data release')
    parser.add_argument('@season', required=False, type=str, help='IceCube season: IC40, IC59, IC79, IC86_I, IC86_II, IC86_III, IC86_IV, IC86_V, IC86_VI, IC86_VII')
    args = parser.parse_args()

    inputDir = args.__dict__["inDir"]
    if args.__dict__["season"]: seasons = [ args.__dict__["season"] ]
    else: seasons = ["IC40", "IC59", "IC79", "IC86_I", "IC86_II", "IC86_III", "IC86_IV", "IC86_V", "IC86_VI", "IC86_VII"]

    for season in seasons:
        #creating root tree with data events
        data = pd.read_csv(inputDir+'/events/'+season+'_exp.csv', delimiter=r"\s+")
        tree = create_roottree_datafile(data, season)
      
        directory = 'data/'
        if not os.path.exists(directory): os.makedirs(directory)
        outfile = root.TFile(directory+season2treename(season), 'recreate')
        tree.Write()
        outfile.Close()

        print('Saving data season '+season+'...')

        #creating root tree with Good Run List
        data    = pd.read_csv(inputDir+'/uptime/'+season+'_exp.csv', delimiter=r"\s+")
        treeGRL = create_roottree_GRLfile(data, season)

        directory = 'data/GRL/'
        if not os.path.exists(directory): os.makedirs(directory)
        outfile = root.TFile(directory+season2treename(season), 'recreate')
        treeGRL.Write()
        outfile.Close()

        print('... and GRL tree')
