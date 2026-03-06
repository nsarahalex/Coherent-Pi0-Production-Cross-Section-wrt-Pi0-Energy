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
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,GENIE_ERROR_GROUPS #,DETAILED_ERROR_GROUPS
from config.DrawingConfig import Categories

ROOT.gStyle.SetPalette(ROOT.kBird)
mnvplotter = PlotUtils.MnvPlotter()
MNVPLOTTER = PlotUtils.MnvPlotter()
#print(mnvplotter.error_summary_group_color_map)
#print(mnvplotter.error_band_color_map)
#CANVAS = ROOT.TCanvas("c2","c2",1600,1000)
SELECTED_SIDEBANDS = AnalysisConfig.sidebands
def GetBins(hist,print_results = True):
    bins = []
    for i in range(hist.GetSize()):
        bins.append(hist.GetBinContent(i))
    if print_results: print(hist,bins)
    return bins


CANVAS = ROOT.TCanvas("c2","c2",1600,1000)



playlist = AnalysisConfig.playlist
print("playlist=" ,playlist)
#mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
#mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(playlist)+"_eavail_weightedmc.root"
mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(playlist)+"_final_weightedmc.root"
mc_file = ROOT.TFile.Open(mc_path)

data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)

pot_scale = Utilities.GetPOTScale(data_path,mc_path)
Sidebands = ["Signal","Diffractive", "Electron_Neutrino", "Pi0","Pi0_muonexits"]
for sides in Sidebands:
    mcEel_hists = HistHolder("Lepton Energy",mc_file,sides,True,pot_scale) #mc_pot,standPOT)
    dataEel_hist = HistHolder("Lepton Energy",data_file,sides,False)  #, data_pot,standPOT)
    
    #mcEel_hists.ResumTotal()
    mcEel_hists.POTScale(True)
    dataEel_hist.POTScale(True)

    mc_list,color,title = mcEel_hists.GetCateList(Categories)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    #PlotTools.MNVPLOTTER.axis_maximum = 1200
    PlotTools.MakeDataMCStackedPlot(dataEel_hist.GetHist(),mc_list,color,title,"TR")
    c1.Update()
    ROOT.gStyle.SetTextSize(0.003)
    #c1.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/Weighed_BWN_Eel_0_"+str(sides)+"_"+str(playlist)+".png","png")
    #c1.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/Weighed_BWN_Eel_1_"+str(sides)+"_"+str(playlist)+".png","png")
    c1.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/Weighed_BWN_Eel_2_"+str(sides)+"_"+str(playlist)+".png","png")
    #PlotTools.MNVPLOTTER.axis_maximum = 0.25
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
    #PlotTools.updatePlotterErrorGroup(GENIE_ERROR_GROUPS)
    PlotTools.MakeErrPlot(mcEel_hists.GetHist())
    #c1.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/Weighed_Eel_error_0_"+str(sides)+"_"+str(playlist)+".png","png")
    #c1.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/Weighed_Eel_error_1_"+str(sides)+"_"+str(playlist)+".png","png")
    c1.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/Weighed_Eel_error_2_"+str(sides)+"_"+str(playlist)+".png","png")
    c1.Close()

