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
from config.DrawingConfig import Categories


mnvplotter = PlotUtils.MnvPlotter()
MNVPLOTTER = PlotUtils.MnvPlotter()
#CANVAS = ROOT.TCanvas("c2","c2",1600,1000)
SELECTED_SIDEBANDS = AnalysisConfig.sidebands
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
            pots[i]= getPOTFromFile(path) #or getPOT(playlist,t,ntuple_tag)
        except KeyError:
            pots[i]=None
    pot_scale = pots[0]/pots[1] if pots.count(None) == 0 else 1.0
    print(pot_scale)
    return pot_scale

CANVAS = ROOT.TCanvas("c2","c2",1600,1000)


def MakePlot(data_hists, mc_hists):
    hists = data_hists#.GetHist()
    mc_list,color,title = mc_hists.GetCateList(Categories,)
    #for hist in mc_list:
    #    hist.PopVertErrorBand("LowQ2Pi_None")
    #    hist.PopVertErrorBand("MK_model")
    #    hist.PopVertErrorBand("SuSA_Valencia_Weight")

    print(title)
    mc_hist_signal=[]
    color_signal=[]
    title_signal=[]

    mc_hist_signal.append(mc_list[10])
    color_signal.append(color[10])
    title_signal.append(title[10])
    Eel_error_mc_canvas = ROOT.TCanvas("c1","c1",1600,1200)
    #PlotTools.MakeDataMCStackedPlot(hists, mc_list,color,title,"TR")
    PlotTools.MakeDataMCStackedPlot(hists, mc_hist_signal,color_signal,title_signal,"TR")
    legend = ROOT.gPad.GetPrimitive("TPave")  # Get the existing legend
    #legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    legend.SetTextSize(0.02)
    if legend:
        legend.AddEntry(hists,"Data after background subtraction", "p")  # Add new entry
        legend.Draw()  # Redraw the legend
    Eel_error_mc_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/BWN_Bkg_subtracted_Data_mcstack"+str(playlist)+".png","png")
    print(f"Plot {mc_hists.plot_name} made.")
    return True

playlist = AnalysisConfig.playlist
print("playlist=" ,playlist)
pot_mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
pot_mc_file = ROOT.TFile.Open(pot_mc_path)
mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(playlist)+"_final_weightedmc.root"
mc_file = ROOT.TFile.Open(mc_path)
data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)

subtracted_data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_subtracted_collab1_fspline.root"
subtracted_data_file = ROOT.TFile.Open(subtracted_data_path)
subtracted_data_Eel = subtracted_data_file.Eel.Clone()
pot_scale = GetPOTScale(data_path,pot_mc_path)
print("pot_scale=", pot_scale)
mc_pot = getPOTFromFile(pot_mc_path)
print("mc pot = " ,mc_pot)
data_pot = getPOTFromFile(data_path)
print("data pot = " , data_pot)
standPOT = data_pot if data_pot is not None else mc_pot
print("stand_pot = ", standPOT)

mcEel_hists = HistHolder("Lepton Energy",mc_file,"Signal",True,pot_scale)#mc_pot,standPOT)
dataEel_hist = HistHolder("Lepton Energy",data_file,"Signal",False)#,data_pot,standPOT)  #, data_pot,standPOT)
mcEel_hists.POTScale(True)
dataEel_hist.POTScale(True)
subtracted_data_Eel.Scale(1,"width")
#MakePlot(dataEel_hist,mcEel_hists)
MakePlot(subtracted_data_Eel,mcEel_hists)
mc_file.Close()
data_file.Close()
pot_mc_file.Close()

