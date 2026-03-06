import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities,PlotTools
from config.SignalDef import SIGNAL_DEFINATION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS


mnvplotter = PlotUtils.MnvPlotter()
#Subtracts background from difefrent sources from signal data and/or mc prediction of data 
def BackgroundSubtraction(data_hists, mc_hists):
    data_hists.POTScale(False)   #pot scaled mc, not bin width normalized
    mc_hists.POTScale(False)     
    out_data = data_hists.GetHist().Clone()
    out_data.AddMissingErrorBandsAndFillWithCV(mc_hists.GetHist())
    first_subtraction = True
    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINATION:  #everything other than NCCohPi0
            out_data = SubtractPoissonHistograms(out_data,mc_hists.hists[group])   #mc_hists.hists[group] should be weighted! #nimmy
    for i in range(out_data.GetSize()): 
        if out_data.GetBinContent(i) < 0: 
            out_data.SetBinContent(i,0.01)
    print("out_data=",out_data.Print())
    return out_data

#Function that do the subtraction do the error propagation
def SubtractPoissonHistograms(signal,background_band,first_subtraction = False): # Sets error of signal where signal - background_band is being calculat
    errors = []
    for i in range(signal.GetSize()): # Adds errors in quadrature
        errors.append(math.sqrt(signal.GetBinError(i)**2 + background_band.GetBinError(i)**2))
        #if first_subtraction: errors.append(math.sqrt(pot_scale*signal.GetBinError(i)**2 + pot_scale*background_band.GetBinError(i)**2)) # Scales the mc total error by the POT scale just once
        #else: errors.append(math.sqrt(signal.GetBinError(i)**2 + pot_scale*background_band.GetBinError(i)**2))
    signal.Add(background_band,-1) # Subtracts the band from signal
    for i in range(signal.GetSize()): # Sets the error of signal to be what was calculated before
        signal.SetBinError(i,errors[i])
    return signal

def GetBins(hist,print_results = True):
    bins = []
    for i in range(hist.GetSize()):
        bins.append(hist.GetBinContent(i))
    if print_results: print(hist,bins)
    return bins

def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None

def GetPOTScale(data_path,mc_path):
    pots = [None,None]
    for i,t in enumerate(["data","mc"]):
        path = [data_path,mc_path][i]
        try:
            pots[i]= getPOTFromFile(path) or getPOT(playlist,t,ntuple_tag)
        except KeyError:
            pots[i]=None
    pot_scale = pots[0]/pots[1] if pots.count(None) == 0 else 1.0
    print(pot_scale)
    return pot_scale


playlist = AnalysisConfig.playlist
print("playlist=" ,playlist)
#pot_mc_path = "/exp/minerva/app/users/nalex/nccoh/AL9-NCCohPi0/tools/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
#pot_mc_file = ROOT.TFile.Open(pot_mc_path)

#mc_path = "/exp/minerva/app/users/nalex/nccoh/AL9-NCCohPi0/tools/kin_dist_"+str(playlist)+"_weightedmc.root"
mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(playlist)+"_final_weightedmc.root"
mc_file = ROOT.TFile.Open(mc_path)

data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)

pot_scale = GetPOTScale(data_path,mc_path)

mcEel_hists = HistHolder("Lepton Energy",mc_file,"Signal",is_mc = True,pot = pot_scale) # Gets the two electron energy histograms
dataEel_hist = HistHolder("Lepton Energy",data_file,"Signal",is_mc = False) # I'm using HistHolders so I can access the backgrounds


Eel_hist_data = dataEel_hist.GetHist().Clone()

new_dataEel_hist = BackgroundSubtraction(dataEel_hist,mcEel_hists)

GetBins(new_dataEel_hist,True)

mc_file.Close()
data_file.Close()

output_datafile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_subtracted_collab1_fspline.root","RECREATE")
new_dataEel_hist.Write()
output_datafile.Close()



new_dataEel_hist.Scale(1,"width")  #now these are potscaled and bin width noirmalized, only for plotting

#Plotting BG subtracted data(new_dataEel_hist) and data (Eel_hist) and make bin width normalized plots
new_dataEel_hist.GetXaxis().SetTitle("Electron Cone Energy (GeV)")
new_dataEel_hist.GetXaxis().SetTitleFont(62)
new_dataEel_hist.GetXaxis().SetTitleSize(0.045)
new_dataEel_hist.GetYaxis().SetTitle("Background Subtracted Bin Counts")
new_dataEel_hist.GetYaxis().SetTitleFont(62)
new_dataEel_hist.GetYaxis().SetTitleSize(0.045)
new_dataEel_hist.SetTitle("Eel Data after background subtraction")
new_dataEel_hist.SetTitleFont(62)
new_dataEel_hist.SetLabelFont(62,"xyz")
new_dataEel_hist.SetMarkerStyle(21)
new_dataEel_hist.SetLineColor(4)
new_dataEel_hist.SetMarkerColor(4)
new_dataEel_hist.SetLineWidth(4)
new_dataEel_hist.SetStats(0)

Eel_hist_data.Scale(1,"width")

Eel_hist_data.SetLineWidth(4)
Eel_hist_data.SetLineColor(2)
Eel_hist_data.SetMarkerColor(2)
Eel_hist_data.SetMarkerStyle(21)
Eel_hist_data.SetTitle("Electron Cone Energy Data")
Eel_canvas = ROOT.TCanvas("c1","c1",1600,1200)
Eel_hist_data.Draw()

new_dataEel_hist.Draw("SAME")
legend = ROOT.gPad.BuildLegend(0.5,0.5,0.7,0.7)
legend.SetTextSize(0.02) 
legend.SetBorderSize(0)
legend.Draw()
#Eel_canvas.SetLogy()
Eel_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/BWN_Eel_subtracted_"+str(playlist)+".png","png")

Eel_canvas.SetLogy()
Eel_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/BWN_Eel_subtracted_"+str(playlist)+"_log.png","png")



Eel_error_data_canvas = ROOT.TCanvas("c1","c1",1600,1200)
PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
PlotTools.MakeErrPlot(new_dataEel_hist)
PlotTools.MNVPLOTTER.axis_maximum = 0.2
#mnvplotter.DrawErrorSummary(new_dataEel_hist)
Eel_error_data_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/Eel_Error_data_subtracted_"+str(playlist)+".png","png")



