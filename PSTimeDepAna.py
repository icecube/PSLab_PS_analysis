#!/usr/bin/env python
import argparse,os,sys
from math import sqrt,log10,radians
from ctypes import c_double
import copy
import numpy as np
import math
import time
from detector_response import EffectiveArea, ResponseDistribution

import ROOT,sys
import os
ROOT.PyConfig.IgnoreCommandLineOptions = True
gROOT = ROOT.gROOT
TCanvas = ROOT.TCanvas
TH2D = ROOT.TH2D
TH2F = ROOT.TH2F
TH1D = ROOT.TH1D
TGraph =ROOT.TGraph
TFile = ROOT.TFile

def set_likelihood(llh, ark, normFromGRL):
    llh.SetUseEnergy(True)
    llh.SetOptimizeTolerance(0.01)
    llh.SetMonitorLevel(0)
    llh.SetEMaxRatioWarnOnlyOnce(1)
    llh.close_ = 10.
    llh.JimsTerm_ = False;
    llh.SpectralPenalty_ = False
    llh.ndof_ = 3
    llh.SetLivetime(ark.livetime/86400.) #to convert in seconds
    #llh.useLCBkgProb_ = False
    llh.SetLocalCoordBkgProb(ark.lcBkgProb)
    llh.SetNormFromGRL(normFromGRL)
    llh.SetMissingRuns(ark.startMissRuns, ark.stopMissRuns)
    #llh.SetOptStoreRatios(True) 
    return llh

def load_detector_response(ark, season):
    if season=="IC86_II" or season=="IC86_III" or season=="IC86_IV" or season=="IC86_V" or season=="IC86_VI" or season=="IC86_VII" or season=="IC86_II_VII":
        season_data="IC86_II"
    else: season_data=season

    #read effective area
    datapath = "data_release/irfs/"+season_data+"_effectiveArea.csv"
    Aeff = EffectiveArea(datapath, season_data)

    #read detector reasponse distribution
    datapath = "data_release/irfs/"+season_data+"_smearing.csv"
    angErr = ResponseDistribution("AngErr", datapath, season_data)
    logEproxy = ResponseDistribution("logE", datapath, season_data)
    angSeparation = ResponseDistribution("PSF", datapath, season_data)

    #link detector response distribution to psark
    ark.SetRecoEnergyDistribution(logEproxy.get_histo())
    ark.SetPSFDistribution(angSeparation.get_histo())
    ark.SetRecoAngErrDistribution(angErr.get_histo())
    ark.SetEffectiveAreaDistribution(Aeff.get_histo())

    return ark

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run psLab',prefix_chars='@')
    parser.add_argument('@srcra',required=True,  type=float, help='RA where to inj/test [deg]')
    parser.add_argument('@srcdec',required=True,  type=float, help='DEC where to inj/test [deg]')
    parser.add_argument('@ns',required=False, default=0, type=int, help='number of neutrinos to inject: 0->background, >0->signal')
    parser.add_argument('@gamma',required=False, default=2, type=float, help='gamma, default 2')
    parser.add_argument('@sigmaT',required=False, default=10, type=float, help='flare duration (in days)')
    parser.add_argument('@upLimSigma',required=False,  default=200.0, type=float, help='upper limit for sigma_time (in days)')
    parser.add_argument('@lowLimSigma',required=False, type=float, help='low limit for sigma_time (in days)')
    parser.add_argument('@nInj',required=False, default=10, type=int, help='number of signal injections to perform in this job, def 10')
    parser.add_argument('@tInj',required=False, default=[57000], nargs="*", type=float, help='Time of injection (central value of Gauss time pdf)')
    parser.add_argument('@seed',required=True, type=int, help='seed')
    parser.add_argument('@nSigTrunc', required=False, default=4, type=float, help='number of sigma for truncated Gauss PDF')
    parser.add_argument('@TSthr',required=False, default=2, type=float, help='TS threshold for multi-flare selection')
    parser.add_argument('@unbl', required=False, type=int, default=0, help='run on real data (1) or generate pseudo-experiment (0)')
    args = parser.parse_args()
 
    ##############
    # parameters #
    ##############
    ##### fixed parameters #####
    ROOT.OPT_USEREALDATA=False
    ROOT.OPT_USENEWDATASTRUCTURE=True
    ROOT.OPT_USEANGFLOOR=True
    ROOT.OPT_NORMONLYOVERACTIVERUNS=True
    
    ##### input parameters #####
    ns=0
    sigmaT      = args.__dict__["sigmaT"]
    TSthr       = args.__dict__["TSthr"]
    seed        = args.__dict__["seed"]
    ns          = args.__dict__["ns"]
    nInjections = args.__dict__["nInj"]
    setUpLimSigma=False
    if args.__dict__["upLimSigma"] is not None:
        upLimSigma=args.__dict__["upLimSigma"]
        setUpLimSigma=True
    setLowLimSigma=False
    if args.__dict__["lowLimSigma"] is not None:
        lowLimSigma=args.__dict__["lowLimSigma"]
        setLowLimSigma=True
    nSigmaTrunc = args.__dict__["nSigTrunc"]
    gamma  = args.__dict__["gamma"]
    srcRA  = args.__dict__["srcra"]
    srcDEC = args.__dict__["srcdec"]
    if args.__dict__["unbl"] == 1:
        ROOT.OPT_USEREALDATA=True
        nInjections = 1
    print("central source coordinates: ra = %.2f, dec=%.2f" %(srcRA,srcDEC))

    ##################
    # load libraries #
    ################## 
    gROOT.Macro("$LAB_MAIN_DIR/llhTimeDep/loadlibs.C")
    ROOT.initialize_ran1(-seed)
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/macro_llh/ArkTimeExt.C+") 

    ######################
    # periods to analyse #
    ######################
    periods=["IC40", "IC59", "IC79", "IC86_I", "IC86_II_VII"] #IC40, IC59, IC79, IC86_I, IC86_II_VII
    ROOT.sample=""
    if len(periods) == 5:
      print("Analysing 10 years (IC40-IC59-IC79-IC86_I-IC86_II_VII)")
      ROOT.sample="_IC10years"
    else:
      for per in periods:
        print("Analysing ", per)
        ROOT.sample+="_"
        ROOT.sample+=per

    #create arks
    psarks=[]
    gROOT.ProcessLine(".L tree_loader.C")
    gROOT.ProcessLine(".L load_ark.C")
    for iark, season in enumerate(periods):
      if season == "IC40":
        gROOT.ProcessLine("I3Ark psark_40;")
        ROOT.psark_40 = load_detector_response(ROOT.psark_40, season)
        gROOT.ProcessLine("load_ark(psark_40, OPT_USEREALDATA, \"IC40\",\"SplineMPE\",OPT_USEANGFLOOR,OPT_USENEWDATASTRUCTURE)")
        psarks.append( ROOT.psark_40 )
      elif season == "IC59":
        gROOT.ProcessLine("I3Ark psark_59;")
        ROOT.psark_59 = load_detector_response(ROOT.psark_59, season)
        gROOT.ProcessLine("load_ark(psark_59, OPT_USEREALDATA, \"IC59\",\"SplineMPE\",OPT_USEANGFLOOR,OPT_USENEWDATASTRUCTURE)")
        psarks.append( ROOT.psark_59 )
      elif season == "IC79":
        gROOT.ProcessLine("I3Ark psark_79;")
        ROOT.psark_79 = load_detector_response(ROOT.psark_79, season)
        gROOT.ProcessLine("load_ark(psark_79, OPT_USEREALDATA, \"IC79\",\"SplineMPE\",OPT_USEANGFLOOR,OPT_USENEWDATASTRUCTURE)")
        psarks.append( ROOT.psark_79 )
      elif season == "IC86_I":
        gROOT.ProcessLine("I3Ark psark_86_I;")
        ROOT.psark_86_I = load_detector_response(ROOT.psark_86_I, season)
        gROOT.ProcessLine("load_ark(psark_86_I, OPT_USEREALDATA,\"IC86_I\",\"SplineMPE\",OPT_USEANGFLOOR,OPT_USENEWDATASTRUCTURE)")
        psarks.append(ROOT.psark_86_I)
      elif season == "IC86_II_VII":
        gROOT.ProcessLine("I3Ark psark_86_II_VII;")
        ROOT.psark_86_II_VII = load_detector_response(ROOT.psark_86_II_VII, season)
        gROOT.ProcessLine("load_ark(psark_86_II_VII, OPT_USEREALDATA,\"IC86_II_VII\",\"SplineMPE\",OPT_USEANGFLOOR,OPT_USENEWDATASTRUCTURE)")
        psarks.append(ROOT.psark_86_II_VII)

    ###########################################################
    # Create lists/vectors of detector response distributions #
    ###########################################################
    gROOT.ProcessLine("vector<TH2D*> effectiveAreas;")
    recoEnergyDistributions, angErrDistributions, PSFDistributions = [], [], []
    for iark in range(len(psarks)):
        recoEnergyDistributions.append( psarks[iark].GetRecoEnergyDistribution() )
        angErrDistributions.append( psarks[iark].GetRecoAngErrDistribution() )
        PSFDistributions.append( psarks[iark].GetPSFDistribution() )
        ROOT.effectiveAreas.push_back( psarks[iark].GetEffectiveAreaDistribution() )
    
    #creating likelihood single periods
    likelihoods=[]
    for per in periods:
      if per == "IC40":
        llh_IC40=ROOT.NewLlhGausTime()
        llh_IC40=set_likelihood(llh_IC40, ROOT.psark_40, ROOT.OPT_NORMONLYOVERACTIVERUNS)
        likelihoods.append(llh_IC40)
      elif per == "IC59":
        llh_IC59=ROOT.NewLlhGausTime()
        llh_IC59=set_likelihood(llh_IC59, ROOT.psark_59, ROOT.OPT_NORMONLYOVERACTIVERUNS)
        likelihoods.append(llh_IC59)
      elif per == "IC79":
        llh_IC79=ROOT.NewLlhGausTime()
        llh_IC79=set_likelihood(llh_IC79, ROOT.psark_79, ROOT.OPT_NORMONLYOVERACTIVERUNS)
        likelihoods.append(llh_IC79)
      elif per == "IC86_I": 
        llh_IC86_I=ROOT.NewLlhGausTime()
        llh_IC86_I=set_likelihood(llh_IC86_I, ROOT.psark_86_I, ROOT.OPT_NORMONLYOVERACTIVERUNS)
        likelihoods.append(llh_IC86_I)
      elif per == "IC86_II_VII":
        llh_IC86_II_VII=ROOT.NewLlhGausTime()
        llh_IC86_II_VII=set_likelihood(llh_IC86_II_VII, ROOT.psark_86_II_VII, ROOT.OPT_NORMONLYOVERACTIVERUNS)
        likelihoods.append(llh_IC86_II_VII)
   
    ########################
    # Creating PS injector #
    ########################
    injSrcCoord=ROOT.EquatorialDeg(srcRA, srcDEC)
    searchCoord=ROOT.EquatorialDeg(srcRA, srcDEC)

    #setting time PDF
    tmin=psarks[0].tmin
    tmax=psarks[-1].tmax
    print("overall time: tmin=%.1f, tmax=%.1f" %(tmin, tmax))
    timeSuff=""
    '''
    timePdfVector = ROOT.std.vector("GaussianTimePdf*")
    tPdf = timePdfVector()
    for i in range(0, len(args.__dict__["tInj"])):
      tmean=args.__dict__["tInj"][i]
      tPdf.push_back(ROOT.GaussianTimePdf(tmin, tmax, tmean, sigmaT, 1.))
      timeSuff+= ("_"+str(tmean))
    '''
    tmean=args.__dict__["tInj"][0]
    gROOT.ProcessLine("GaussianTimePdf* tPdf = new GaussianTimePdf("+str(tmin)+", "+str(tmax)+", "+str(tmean)+", "+str(sigmaT)+", 1.);")
 
    #setting power law
    gROOT.ProcessLine("PowerLawFlux *flux = new PowerLawFlux(1,-"+str(gamma)+");")

    ######## creating multiArk ########
    gROOT.ProcessLine("MultiArk mark;")
    for ark in psarks:
      ROOT.mark.AddArk(ark) 
    ROOT.mark.SetPointSource(injSrcCoord, ROOT.flux, ROOT.tPdf, ROOT.OPT_NORMONLYOVERACTIVERUNS)  

    ######## Linking the ark and point source to the PDF  ########
    for per in range(0, len(periods)):
      likelihoods[per].SetAnalysis(psarks[per].psData, searchCoord)                                                                                                                   
    maf=ROOT.MultiGaussAnalysisFn()
    for llh in likelihoods:
      maf.AddAnalysisFn(llh)

    maf.SetTimeBounds(tmin,tmax)
    maf.SetSearchCoord(searchCoord)
    maf.JimsTerm_ = True
    maf.SetNSigmaTrunc(nSigmaTrunc)
    maf.SetTSthr(TSthr)
    maf.SetOptStoreRatios(True)

    sigmaLim=""
    if setUpLimSigma == True:
        maf.SetUpLimitSigma(upLimSigma)
        print("setting upper limit sigma_time = ", upLimSigma)
        sigmaLim+="_upLimSigma_"+str(upLimSigma)
    if setLowLimSigma == True:
        maf.SetLowLimitSigma(lowLimSigma)
        print("setting low limit sigma_time = ", lowLimSigma)
        sigmaLim+="_lowLimSigma_"+str(lowLimSigma)

    ######## Setting parTranslator
    pt=ROOT.NewLlhGausTime_ParTranslator()
    pt.SetRange(1,4,31) #gamma_min, gamma_max, nBins
    gROOT.ProcessLine("MultiAnalysisSet* mas = dynamic_cast<MultiAnalysisSet*>(mark.psData);")
    pt.SetTranslator(ROOT.mas)
    maf.SetParTranslator(pt)

    gROOT.ProcessLine(".L SimpleAnalysisPS.C+")
    gROOT.ProcessLine("SimpleAnalysis_multiSet mps;") 
    ROOT.mps.SetRndOnlyTimes(False)
    print("initializing histos..")
    ROOT.mps.Initialize(ns, srcRA, srcDEC, 5)
    if ROOT.OPT_USEREALDATA == True:
        ROOT.mps.SetDoUnblind()

    fit_time = []
    print("processing sigma", sigmaT)
    for i in range(nInjections):
        fitStart_time = time.time()
        print("------ TRIAL ", i+1)
        for ark in psarks:
          nAzResetPS=ark.lcBkgProb.nbAz
          if sigmaT < 8e-2:
            nAzResetPS=int(360./sqrt(sigmaT))
          ark.psData.GetSource().SetTimeAzBins(nAzResetPS)
        print("Executing SimpleAnalysis_multiSet..")
        ROOT.mps.Execute(ROOT.mark, maf)
        fitEnd_time = time.time()
        print('fit time = ', fitEnd_time - fitStart_time)
        fit_time.append(fitEnd_time - fitStart_time)

    print('average fit time = ', sum(fit_time)/len(fit_time))

   
    directory="output/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory+="dec"+str(srcDEC)+"_ra"+str(srcRA)
    if not os.path.exists(directory):
        os.makedirs(directory)

    ROOT.mps.Write(""+directory+"/TS"+str(ROOT.sample)+"_RA_"+str(srcRA)+"_DEC_"+str(srcDEC)+"_tsigma_"+str(sigmaT)+"_time"+timeSuff+"_ns_"+str(ns)+sigmaLim+"_gamma_"+str(gamma)+"_nSigmaTrunc_"+str(nSigmaTrunc)+"_TSthr_"+str(TSthr)+"_seed_"+str(seed)+".root", "RECREATE")
   
    
    
    
