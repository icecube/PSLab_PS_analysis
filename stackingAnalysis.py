#!/usr/bin/env python
import argparse,os,sys
import copy
from math import sqrt,log10,radians
import copy
import numpy as np
import math
import time
from detector_response import EffectiveArea, ResponseDistribution

import ROOT
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
    llh.ndof_ = 3.
    llh.SetLivetime(ark.livetime/86400.) #to convert in seconds
    llh.SetLocalCoordBkgProb(ark.lcBkgProb)
    return llh

def calculate_fluencescale(fluxscale, tminvect, tmaxvect, weights):
    nSrc = len(weights)
    WT = 0
    for isrc in range(nSrc): WT += ( weights[isrc] * (tmaxvect[isrc]-tminvect[isrc]) )
    return WT*fluxscale

def load_detector_response(ark, season_data): 

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
    parser.add_argument('@ns',required=True,  type=int, help='number of signal-like neutrinos to inject')
    parser.add_argument('@gamma',required=False, default=2, type=float, help='gamma, default 2')
    parser.add_argument('@nInj',required=False,  type=int, default=10, help='number of signal injections to perform in this job, def 10')
    parser.add_argument('@s',required=True,  type=int, help='seed')
    parser.add_argument('@hemi', required=True,  type=str, help='north or south')
    parser.add_argument('@unbl', required=False, type=int, default=0, help='run on real data (1) or scramble (0)')
    args = parser.parse_args()

    ROOT.OPT_USEREALDATA=False
    ##############
    # parameters #
    ##############
    ##### input parameters #####
    seed        = args.__dict__["s"]
    nInjections = args.__dict__["nInj"]
    gamma       = args.__dict__["gamma"]
    ns          = args.__dict__["ns"]
    hemi        = args.__dict__["hemi"]
    if args.__dict__["unbl"] == 1: ROOT.OPT_USEREALDATA=True
    elif args.__dict__["unbl"] == 0: ROOT.OPT_USEREALDATA=False
    else: print("ERROR: UNBL VALUE NOT CORRECT. SCRAMBLING DATA")

    ##### fixed parameters #####
    ROOT.OPT_USENEWDATASTRUCTURE=True
    ROOT.OPT_USEANGFLOOR=True
    ROOT.OPT_NORMONLYOVERACTIVERUNS=False
    gROOT.SetBatch()

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

    #########################
    # Set detector response #
    ######################### 
    ##### create arks #####
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
        ROOT.psark_86_II_VII = load_detector_response(ROOT.psark_86_II_VII, "IC86_II")
        gROOT.ProcessLine("load_ark(psark_86_II_VII, OPT_USEREALDATA,\"IC86_II_VII\",\"SplineMPE\",OPT_USEANGFLOOR,OPT_USENEWDATASTRUCTURE)")
        psarks.append(ROOT.psark_86_II_VII)

    gROOT.ProcessLine("vector<double> perTminVect;")
    gROOT.ProcessLine("vector<double> perTmaxVect;")
    for ark in psarks:
        ROOT.perTminVect.push_back(ark.tmin)
        ROOT.perTmaxVect.push_back(ark.tmax)  
    
    ###############
    # Source info #
    ###############
    filename='catalogue.txt'
    data = np.loadtxt(filename, unpack=True,
                     dtype={'names': ('name', 'ra', 'dec', 'tmin', 'tmax', 'angErr', 'fluence'), 'formats': ('|S15', float, float, float, float, float, float)})
    src_name = np.array(data[0])
    ra       = np.array(data[1])
    dec      = np.array(data[2])
    tmin     = np.array(data[3])
    tmax     = np.array(data[4])
    sigma    = np.array(data[5])
    fluence  = np.array(data[6])

    mask = sigma<0.001
    sigma[mask] = 0.001

    if hemi=='north':   mask_hemi = (dec>=-5)
    elif hemi=='south': mask_hemi = (dec<-5)
    else: sys.exit('hemisphere ({}) not spelled correctly'.format(hemi))
    src_name = src_name[mask_hemi]
    ra       = ra[mask_hemi]
    dec      = dec[mask_hemi]
    tmin     = tmin[mask_hemi]
    tmax     = tmax[mask_hemi]
    sigma    = sigma[mask_hemi]
    fluence  = fluence[mask_hemi]
    nSrc = len(dec)
    print('sources selected in hemisphere {}: {}'.format(hemi, nSrc)) 

    gROOT.ProcessLine("vector<EquatorialDeg> srcLocations;")
    gROOT.ProcessLine("vector<double> srcSigmas;")
    gROOT.ProcessLine("vector<double> srcTminVect;")
    gROOT.ProcessLine("vector<double> srcTmaxVect;")
    gROOT.ProcessLine("vector<double> enhancementFactor;")
    WeightTot=0
    for i in range(nSrc):
        ROOT.srcLocations.push_back( ROOT.EquatorialDeg(ra[i], dec[i]) )
        ROOT.srcSigmas.push_back( sigma[i] )
        ROOT.srcTminVect.push_back(tmin[i])
        ROOT.srcTmaxVect.push_back(tmax[i])
        weight = fluence[i]
        WeightTot += weight
        ROOT.enhancementFactor.push_back(weight)

    for i in range(len(ROOT.enhancementFactor)):
        ROOT.enhancementFactor[i] /= WeightTot
 
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

    #########################
    # Creating PS injectors #
    #########################
    fluxScale = 0
    gROOT.ProcessLine("vector<I3Event> sourceEvents;")
    gROOT.ProcessLine("vector<I3MultiSignalGenerator*> msg;")
    gROOT.ProcessLine("I3MultiSignalGenerator* msg_tmp;")
    gROOT.ProcessLine("TH2D* effectiveArea;")
    gROOT.ProcessLine("BoxTimePdf *tPdf;")
    gROOT.ProcessLine("PowerLawFlux* pflux;")
    gROOT.ProcessLine("I3PointGenerator* i3point;")
    gROOT.ProcessLine("double nsrcnev;")
    iflux=0
    for iark, ark in enumerate(psarks):
        gROOT.ProcessLine("msg_tmp = new I3MultiSignalGenerator();")
        gROOT.ProcessLine("effectiveArea = new TH2D();")
        for isrc in range(nSrc):
            gROOT.ProcessLine("tPdf = new BoxTimePdf("+str(ark.tmin)+", "+str(ark.tmax)+", "+str(ROOT.srcTminVect[isrc])+", "+str(ROOT.srcTmaxVect[isrc])+", 1. , 0.);")
            if ROOT.tPdf.GetNorm()==0: continue
            gROOT.ProcessLine("pflux = new PowerLawFlux(1,-"+str(gamma)+");")
            gROOT.ProcessLine("i3point = new I3PointGenerator(*pflux, srcLocations["+str(isrc)+"], "+str(ark.livetime)+", tPdf, effectiveAreas["+str(iark)+"]);")
            #ROOT.i3point.SetEffectiveAreaDistribution(ROOT.effectiveAreas[iark])
            ROOT.i3point.SetSourceSigmas(ROOT.srcSigmas[isrc])
            gROOT.ProcessLine("nsrcnev = i3point->GetMeanSrcNev();")
            gROOT.ProcessLine("msg_tmp->AddSignal(*i3point, "+str(ROOT.enhancementFactor[isrc])+");")
            fluxScale += ROOT.nsrcnev * ROOT.tPdf.GetNorm() * ROOT.enhancementFactor[isrc]
            iflux+=1
        ROOT.msg_tmp.SetRecoLogEproxyDistribution( recoEnergyDistributions[iark] )
        ROOT.msg_tmp.SetRecoAngErrDistribution( angErrDistributions[iark] )
        ROOT.msg_tmp.SetPSFDistribution( PSFDistributions[iark] )
        gROOT.ProcessLine("msg.push_back(msg_tmp);")
    if fluxScale!=0: fluxScale = 1./fluxScale
    
    #######################
    # Adding sourceModule #
    #######################
    gROOT.ProcessLine("vector<SourceModule*> srcMods;")
    for iark in range(len(psarks)): gROOT.ProcessLine("srcMods.push_back(msg["+str(iark)+"]);")
    
    #############################
    # Stacking weights VS gamma #
    #############################
    gammaMin, gammaMax, nGammaBins = -4, -1, 30
    srcWeightsArray = []
    gROOT.ProcessLine("vector<vector<double>> srcWeightsArray_tmp;")
    gROOT.ProcessLine("srcWeightsArray_tmp.resize( "+str(len(ROOT.srcLocations))+");")
    gROOT.ProcessLine("FluxBase *fluxPtr;")
    for iark, ark in enumerate(psarks):
        gROOT.ProcessLine("effectiveArea = new TH2D();")
        ROOT.effectiveArea = psarks[iark].GetEffectiveAreaDistribution()
        for isrc in range(len(ROOT.srcLocations)): ROOT.srcWeightsArray_tmp[isrc].clear()
        testGamma = 0
        for isrc in range(len(ROOT.srcLocations)):
            gROOT.ProcessLine("tPdf = new BoxTimePdf("+str(ark.tmin)+", "+str(ark.tmax)+", "+str(ROOT.srcTminVect[isrc])+", "+str(ROOT.srcTmaxVect[isrc])+", 1. , 0.);")
            for gidx in range(nGammaBins):
                if ROOT.tPdf.GetNorm()==0: srcWeight = 0
                else:
                    testGamma = (gidx+0.5)*(gammaMax-gammaMin)/nGammaBins + gammaMin
                    gROOT.ProcessLine("fluxPtr = new PowerLawFlux(1, "+str(testGamma)+");")
                    gROOT.ProcessLine("i3point = new I3PointGenerator(*fluxPtr, srcLocations["+str(isrc)+"], "+str(ark.livetime)+", tPdf, effectiveAreas["+str(iark)+"]);") 
                    ROOT.i3point.SetEffectiveAreaDistribution(ROOT.effectiveArea)
                    gROOT.ProcessLine("nsrcnev = i3point->GetMeanSrcNev();")
                    srcWeight = ROOT.nsrcnev * ROOT.tPdf.GetNorm() * ROOT.enhancementFactor[isrc]
                gROOT.ProcessLine("srcWeightsArray_tmp["+str(isrc)+"].push_back("+str(srcWeight)+");")

        ##### Normalize the weights #####
        for gidx in range(nGammaBins):
            sumRow = 0
            for isrc in range(len(ROOT.srcLocations)): sumRow += ROOT.srcWeightsArray_tmp[isrc][gidx]
            for isrc in range(len(ROOT.srcLocations)):
                ROOT.srcWeightsArray_tmp[isrc][gidx] /= sumRow
        srcWeightsArray.append( copy.deepcopy(ROOT.srcWeightsArray_tmp) ) 

    ##################################
    # Initializeing stacking classes #
    ##################################
    maf = ROOT.MultiBoxAnalysisFnStack()
    llh_stack_tmp = []
    for iark, ark in enumerate(psarks):
        psarks[iark].SetSource(ROOT.srcMods[iark])
        llh_stack_tmp.append( ROOT.NewLlhBoxStack() )
        llh_stack_tmp[iark].SetUseEnergy(True)
        llh_stack_tmp[iark].SetOptimizeTolerance(0.00)
        llh_stack_tmp[iark].SetMonitorLevel(0)
        llh_stack_tmp[iark].SetAnalysisSet(psarks[iark].psData)
        llh_stack_tmp[iark].SetSourceCoords(ROOT.srcLocations)
        llh_stack_tmp[iark].SetStackedWeightTable(srcWeightsArray[iark])
        llh_stack_tmp[iark].SetSourceSigmas(ROOT.srcSigmas)
        llh_stack_tmp[iark].SetSourceTminTmax(ROOT.srcTminVect, ROOT.srcTmaxVect)
        llh_stack_tmp[iark].SetPeriodTminTmax(ark.tmin, ark.tmax)
        maf.AddAnalysisFn(llh_stack_tmp[iark])

    ##### create multiark #####
    gROOT.ProcessLine("MultiArk mark;")
    for ark in psarks: ROOT.mark.AddArk(ark)

    #########################
    # Setting parTranslator #
    #########################
    pt=ROOT.NewLlhBoxTime_ParTranslator()
    pt.SetRange(1,4,31) #gamma_min, gamma_max, nBins
    pt.SetPeriodsTimeBounds(ROOT.perTminVect, ROOT.perTmaxVect)
    pt.SetSourceTimeBounds(ROOT.srcTminVect, ROOT.srcTmaxVect)
    pt.SetEnhanceFactor(ROOT.enhancementFactor)
    gROOT.ProcessLine("MultiAnalysisSet* mas = dynamic_cast<MultiAnalysisSet*>(mark.psData);")
    pt.SetTranslatorStacking(ROOT.mas)
    maf.SetParTranslator(pt)

    ##############################################
    # Calculating flux and fluence scale factors #
    ##############################################
    #NOTE: the fluenceScale/fluxScale must be MULTIPLIED by ns to obtain the corresponding fluence/flux
    fluenceScale = calculate_fluencescale(fluxScale, ROOT.srcTminVect, ROOT.srcTmaxVect, ROOT.enhancementFactor) # GeV cm^-2, normalization at 1 GeV
    fluenceScale = fluenceScale/1000. #TeV cm^-2, normalization at 1 GeV

    ###########################
    # Creating SimpleAnalysis #
    ###########################
    gROOT.ProcessLine(".L SimpleAnalysisStack.C+")
    gROOT.ProcessLine("SimpleAnalysis_multiSet mps;")
    ROOT.mps.Initialize(ns, len(periods), fluenceScale, fluxScale)
    if ROOT.OPT_USEREALDATA == True:
        ROOT.mps.SetDoUnblind()

    fit_time = []
    for i in range(nInjections):
        fitStart_time = time.time()
        print("------ TRIAL ", i+1)
        print("Executing SimpleAnalysis_multiSet..")
        print("fluxscale = {:.2e}, fluencescale = {:.2e}".format(fluxScale, fluenceScale))
        maf.SetGammaFixed(gamma)
        ROOT.mps.Execute(ROOT.mark, maf)
        fitEnd_time = time.time()
        print('fit time = ', fitEnd_time - fitStart_time)
        fit_time.append(fitEnd_time - fitStart_time)

    print('average fit time = ', sum(fit_time)/len(fit_time))

    directory="output/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory+="stacking/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory+=(hemi+"/")
    if not os.path.exists(directory):
        os.makedirs(directory)

    if ROOT.OPT_USEREALDATA == False: ROOT.mps.Write(""+directory+"TS_stacking_"+hemi+"_ns_"+str(ns)+"_gamma_"+str(gamma)+"_seed_"+str(seed)+".root", "RECREATE")
    else: ROOT.mps.Write(""+directory+"TS_stacking_"+hemi+"_gamma_"+str(gamma)+".root", "RECREATE")
        
 
