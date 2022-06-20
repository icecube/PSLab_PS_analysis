import pandas as pd
import numpy as np
import ROOT as root

TH1DVector = root.std.vector("TH1D")
TH1DVectorOfVector = root.std.vector( root.std.vector("TH1D") )

class DetectorResponse(object):
    def __init__(self, datapath, season):
        self.season     = season
        self.data       = pd.read_csv(datapath, delimiter=r"\s+")
        mask = (self.data['logEmin_nu']<9) #Since IC40 has effective area values up to 1e10
        self.logEbins   = self.generate_bins(self.data['logEmin_nu'][mask], self.data['logEmax_nu'][mask])
        self.decbins    = self.generate_bins(self.data['decmin'], self.data['decmax'])
        self.sindecbins = np.sin(np.deg2rad(self.decbins))

    def generate_bins(self, binmin, binmax):
        binmin_unique  = np.unique(binmin)
        binmax_unique  = np.unique(binmax)
        binedge_unique = np.append(binmin_unique, binmax_unique[-1])
        return binedge_unique

    def get_histo(self):
        return self.histo

class EffectiveArea(DetectorResponse):
    def __init__(self, datapath, season):
        super(EffectiveArea, self).__init__(datapath, season)
        self.histo = root.TH2D('Aeff'+season, 'Aeff'+season, len(self.logEbins)-1, self.logEbins, len(self.sindecbins)-1,  self.sindecbins)
        self.fill_histo()

    def set_bin_content(self, logE, dec):
        mask = ( (self.data['logEmin_nu']==logE) & (self.data['decmin']==dec) )
        Aeff = self.data['Aeff'][mask]

        binx = self.histo.GetXaxis().FindBin(logE)
        biny = self.histo.GetYaxis().FindBin(np.sin(np.deg2rad(dec)))
        self.histo.SetBinContent(binx, biny, Aeff)

    def fill_histo(self):  
        for logE in self.logEbins[:-1]:
            if logE>9: break
            for dec in self.decbins[:-1]:
                self.set_bin_content(logE, dec)

    def save_histo(self, outfilepath):
        canvas = root.TCanvas("c_Aeff_"+self.season, "c_Aeff_"+self.season, 800, 600)
        self.histo.SetMinimum(0.999e1)
        self.histo.SetMaximum(1e7)
        self.histo.GetXaxis().SetRange(1, 35)

        self.histo.GetXaxis().SetTitle("log_{10}(E/GeV)")
        self.histo.GetYaxis().SetTitle("\sin(\delta)")
        self.histo.GetZaxis().SetTitle("Effective area [cm^{2}]")

        self.histo.SetContour(100)
        self.histo.Draw('colz0')

        canvas.SetLogz()
        root.gStyle.SetOptStat(0)
        canvas.Draw()

        fileOutput =  root.TFile(outfilepath, 'recreate')
        canvas.Write()
        fileOutput.Close()

class ResponseDistribution(DetectorResponse):
    def __init__(self, parameter, datapath, season):
        super(ResponseDistribution, self).__init__(datapath, season) 
        self.histo = root.std.vector( root.std.vector("TH1D*") )() #TH1DVectorOfVector() #root.std.vector( root.std.vector("TH1D*") )() 
        self.par   = parameter
        self.fill_response_distribution()

    def get_frac_counts(self, parmin, counts):
        frac_counts = []
        for par in np.unique(parmin):
            mask = (parmin==par)
            frac_counts.append( np.sum(counts[mask]) )
        return frac_counts 
 
    def make_histo(self, name, bins, counts):
        hist = root.TH1D(name, name, len(bins)-1, bins)
        for i in range(len(bins[:-1])):
            binx = hist.GetXaxis().FindBin(bins[i])
            hist.SetBinContent(binx, counts[i])
        return hist

    def fill_response_distribution(self):
        for dec in self.decbins[:-1]:
            histo_decBand = root.std.vector("TH1D*")() #TH1DVector()
            for logE in self.logEbins[:-1]:
                mask    = ((self.data['logEmin_nu']==logE) & (self.data['decmin']==dec))
                parmin      = self.data[self.par+'min'][mask]
                parmax      = self.data[self.par+'max'][mask]
                parbins     = self.generate_bins(parmin, parmax)
                frac_counts = self.get_frac_counts(parmin, self.data['Fractional_Counts'][mask])
              
                histo_name = self.par+'_'+self.season+'_logE_'+str(logE)+'_sindec_'+str(dec) 
                histo_decBand.push_back( self.make_histo(histo_name, parbins, frac_counts) )

            self.histo.push_back(histo_decBand)
        #self.histo = root.std.vector(root.std.vector("TH1D*"))(self.histo)

    def save_histo(self, outfilepath, hemi):
        fileOutput =  root.TFile(outfilepath, 'recreate')

        for i in range(len(self.histo[hemi])): self.histo[hemi][i].Write()

        fileOutput.Close()
