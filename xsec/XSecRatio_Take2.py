import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities, PlotTools
from tools.PlotLibrary import PLOT_SETTINGS
from config.SystematicsConfig import USE_NUE_CONSTRAINT,CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS,AnaNuPDG
from config import DrawingConfig
from config.CutConfig import NEUTRINO_ENERGY_RANGE,FIDUCIAL_Z_RANGE,FIDUCIAL_APOTHEM
from functools import partial
import numpy as np
mnvplotter = PlotUtils.MnvPlotter()

ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch(True)

playlist = AnalysisConfig.playlist
#FHCFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e/NCCoh_04_16_FHC_xsec.root","READ")
#RHCFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e/NCCoh_04_10_RHC_xsec.root","READ")
RHCFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/Coh_RHC_01_26_xsec.root","READ")
FHCFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/Coh_FHC_01_26_xsec.root","READ")

PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)

RHC_data_hist = RHCFile.pionE_recovered.Clone()
RHC_sim_hist = RHCFile.PiZeroE_NCCohPi0.Clone()#NCCohPi0.Clone()
RHC_sim_hist.AddMissingErrorBandsAndFillWithCV(RHC_data_hist)

FHC_data_hist = FHCFile.pionE_recovered.Clone()
FHC_sim_hist = FHCFile.PiZeroE_NCCohPi0.Clone()
FHC_data_hist.AddMissingErrorBandsAndFillWithCV(RHC_data_hist)
FHC_sim_hist.AddMissingErrorBandsAndFillWithCV(RHC_data_hist)

RHC_FHC_ratio_hist_data = RHCFile.PiZeroE_NCCohPi0.Clone()#_NCCohPi0.Clone()
RHC_FHC_ratio_hist_sim = RHCFile.PiZeroE_NCCohPi0.Clone()#_NCCohPi0.Clone()
RHC_ratio = RHCFile.pionE_recovered.Clone()
FHC_ratio = FHCFile.pionE_recovered.Clone()

print("RHC_data_hist:", list(RHC_data_hist.GetVertErrorBandNames()))
print("RHC_sim_hist:", list(RHC_sim_hist.GetVertErrorBandNames()))
print("FHC_data_hist:", list(FHC_data_hist.GetVertErrorBandNames()))
print("FHC_sim_hist:", list(FHC_sim_hist.GetVertErrorBandNames()))


RHC_canvas = ROOT.TCanvas("c2","c2",1600,1200)
mnvplotter.DrawDataMCRatio(RHC_data_hist, RHC_sim_hist, 1.0 ,True,True,0,5)
RHC_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/RHC_ratio_MNV.png","png")

FHC_canvas = ROOT.TCanvas("c2","c2",1600,1200)
mnvplotter.DrawDataMCRatio(FHC_data_hist, FHC_sim_hist, 1.0 ,True,True,0,5)
FHC_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/FHC_ratio_MNV.png","png")


RHCFile.Close()
FHCFile.Close()

print(RHC_data_hist.GetVertErrorBandNames())

RHC_ratio.Divide(RHC_data_hist,RHC_sim_hist)
FHC_ratio.Divide(FHC_data_hist,FHC_sim_hist)

RHC_FHC_ratio_hist_sim.Divide(RHC_sim_hist,FHC_sim_hist)

RHC_FHC_ratio_hist_data.Divide(RHC_data_hist,FHC_data_hist)

outPath = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/Nimmy_XSec_ratio.root"
outFile = ROOT.TFile.Open(outPath,"RECREATE")

RHC_canvas = ROOT.TCanvas("c2","c2",1600,1200)

RHCratio = RHC_ratio.GetCVHistoWithError()
RHCratio.Print("All")
RHCratio.SetMarkerStyle(21)
RHCratio.SetMarkerSize(1.5)
RHCratio.SetMarkerColor(4)
RHCratio.SetLineColor(9)
RHCratio.SetMaximum(max([RHCratio.GetMaximum(),RHCratio.GetMaximum()])*1.25)
RHCratio.GetXaxis().SetTitle("Pion Energy(GeV)")
RHCratio.GetXaxis().SetTitleFont(62)
RHCratio.GetXaxis().SetRangeUser(0,0.3)
RHCratio.GetYaxis().SetRangeUser(0,6)
RHCratio.GetXaxis().SetTitleSize(0.045)
RHCratio.GetYaxis().SetTitle("Data/Simulation, XSecRatio")
RHCratio.GetYaxis().SetTitleFont(62)
RHCratio.GetYaxis().SetTitleSize(0.045)
RHCratio.SetStats(0)
RHCratio.SetTitle("Data")
RHCratio.SetLabelFont(62,"xyz")
RHCratio.SetLineWidth(4)

RHCratio.Draw("MINO")
#pion_legend = ROOT.gPad.BuildLegend(0.6,0.7,0.85,0.85)
#pion_legend.SetTextFont(62)
#pion_legend.SetTextSize(0.035)
RHC_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/RHC_ratio.png","png")
PlotTools.MakeErrPlot(RHC_data_hist)
RHC_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/RHC_ratio_error.png","png")



FHC_canvas = ROOT.TCanvas("c2","c2",1600,1200)

FHCratio = FHC_ratio.GetCVHistoWithError()
FHCratio.Print("All")
FHCratio.SetMarkerStyle(21)
FHCratio.SetMarkerSize(1.5)
FHCratio.SetMarkerColor(4)
FHCratio.SetLineColor(9)
FHCratio.SetMaximum(max([FHCratio.GetMaximum(),FHCratio.GetMaximum()])*1.25)
FHCratio.GetXaxis().SetTitle("Pion Energy(GeV)")
FHCratio.GetXaxis().SetTitleFont(62)
FHCratio.GetXaxis().SetRangeUser(0,0.3)
FHCratio.GetYaxis().SetRangeUser(0,6)
FHCratio.GetXaxis().SetTitleSize(0.045)
FHCratio.GetYaxis().SetTitle("Data/Simulation, XSecRatio")
FHCratio.GetYaxis().SetTitleFont(62)
FHCratio.GetYaxis().SetTitleSize(0.045)
FHCratio.SetStats(0)
#FHCratio.SetTitle("Data")
FHCratio.SetLabelFont(62,"xyz")
FHCratio.SetLineWidth(4)

#FHCratio.Scale(1,"width")
#pion_fitted_err.Draw("MIN0 E1")
FHCratio.Draw("MINO")
#pion_legend = ROOT.gPad.BuildLegend(0.6,0.7,0.85,0.85)
#pion_legend.SetTextFont(62)
#pion_legend.SetTextSize(0.035)
FHC_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/FHC_ratio.png","png")
PlotTools.MakeErrPlot(FHC_data_hist)
FHC_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/FHC_ratio_error.png","png")


RHC_FHC_ratio_data_canvas = ROOT.TCanvas("c2","c2",1600,1200)

dataRHCFHCratio = RHC_FHC_ratio_hist_data#.GetCVHistoWithError()
dataRHCFHCratio.Print("All")
dataRHCFHCratio.SetMarkerStyle(21)
dataRHCFHCratio.SetMarkerSize(1.5)
dataRHCFHCratio.SetMarkerColor(4)
dataRHCFHCratio.SetLineColor(9)
dataRHCFHCratio.SetMaximum(max([dataRHCFHCratio.GetMaximum(),dataRHCFHCratio.GetMaximum()])*1.25)
dataRHCFHCratio.GetXaxis().SetTitle("Pion Energy(GeV)")
dataRHCFHCratio.GetXaxis().SetTitleFont(62)
#dataRHCFHCratio.GetXaxis().SetRangeUser(0,0.3)
dataRHCFHCratio.GetYaxis().SetRangeUser(0,2)
dataRHCFHCratio.GetXaxis().SetTitleSize(0.045)
dataRHCFHCratio.GetYaxis().SetTitle("RHC/FHC DATA, XSecRatio")
dataRHCFHCratio.GetYaxis().SetTitleFont(62)
dataRHCFHCratio.GetYaxis().SetTitleSize(0.045)
dataRHCFHCratio.SetStats(0)
#FHCratio.SetTitle("Data")
dataRHCFHCratio.SetLabelFont(62,"xyz")
dataRHCFHCratio.SetLineWidth(4)

#dataRHCFHCratio.Scale(1,"width")
#pion_fitted_err.Draw("MIN0 E1")
dataRHCFHCratio.Draw("MINO")
#pion_legend = ROOT.gPad.BuildLegend(0.6,0.7,0.85,0.85)
#pion_legend.SetTextFont(62)
#pion_legend.SetTextSize(0.035)
RHC_FHC_ratio_data_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/RHC_FHC_ratio_data.png","png")

PlotTools.MakeErrPlot(RHC_FHC_ratio_hist_data)
RHC_FHC_ratio_data_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/RHC_FHC_ratio_data_error.png","png")


RHC_FHC_ratio_sim_canvas = ROOT.TCanvas("c2","c2",1600,1200)

simRHCFHCratio = RHC_FHC_ratio_hist_sim#.GetCVHistoWithError()
simRHCFHCratio.Print("All")
simRHCFHCratio.SetMarkerStyle(21)
simRHCFHCratio.SetMarkerSize(1.5)
simRHCFHCratio.SetMarkerColor(4)
simRHCFHCratio.SetLineColor(9)
simRHCFHCratio.SetMaximum(max([simRHCFHCratio.GetMaximum(),simRHCFHCratio.GetMaximum()])*1.25)
simRHCFHCratio.GetXaxis().SetTitle("Pion Energy(GeV)")
simRHCFHCratio.GetXaxis().SetTitleFont(62)
#simRHCFHCratio.GetXaxis().SetRangeUser(0,0.3)
simRHCFHCratio.GetYaxis().SetRangeUser(0,2)
simRHCFHCratio.GetXaxis().SetTitleSize(0.045)
simRHCFHCratio.GetYaxis().SetTitle("RHC/FHC SIM, XSecRatio")
simRHCFHCratio.GetYaxis().SetTitleFont(62)
simRHCFHCratio.GetYaxis().SetTitleSize(0.045)
simRHCFHCratio.SetStats(0)
#FHCratio.SetTitle("Data")
simRHCFHCratio.SetLabelFont(62,"xyz")
simRHCFHCratio.SetLineWidth(4)

#simRHCFHCratio.Scale(1,"width")
#pion_fitted_err.Draw("MIN0 E1")
simRHCFHCratio.Draw("MINO")
#pion_legend = ROOT.gPad.BuildLegend(0.6,0.7,0.85,0.85)
#pion_legend.SetTextFont(62)
#pion_legend.SetTextSize(0.035)
RHC_FHC_ratio_sim_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/RHC_FHC_ratio_sim.png","png")
PlotTools.MakeErrPlot(RHC_FHC_ratio_hist_sim)
RHC_FHC_ratio_sim_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/RHC_FHC_ratio_sim_error.png","png")

RHC_FHC_ratio_hist_sim.Write()
RHC_FHC_ratio_hist_data.Write()
RHC_data_hist.Write()
FHC_data_hist.Write()
outFile.Close()
