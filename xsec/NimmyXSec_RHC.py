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

ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch(True)

USE_BIGNUE = True
threshold = 100 if USE_BIGNUE else 1
TARGET_UTILS = PlotUtils.TargetUtils.Get()
warping_errorband = ["fsi_weight","SuSA_Valencia_Weight","MK_model","LowQ2Pi_Joint","LowQ2Pi_NUPI0","LowQ2Pi_None"]
#warping_errorband = ["SuSA_Valencia_Weight"]
FLUX5a="minervame5a"
FLUX6a="minervame6a"

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

def GetXSectionHistogram(unfolded,efficiency,is_mc):
    #divide by efficiency
    efficiency.AddMissingErrorBandsAndFillWithCV(unfolded)
    unfolded.Divide(unfolded,efficiency)
    #divide by flux
    DivideFlux(unfolded,is_mc)
    #divide by N Hydrogen
    nplanes = 2*(80-27+1) # Fiducial Volume -> Modules 27-80
    mass = TARGET_UTILS.GetTrackerMass(nplanes,is_mc,FIDUCIAL_APOTHEM)
    massFrac = TARGET_UTILS.GetTrackerElementMassFraction(6,is_mc)
    targets = mass*massFrac*TARGET_UTILS.GetTrackerElementAtomsPerGram(6)#*6 # That last *1(*6) is because there's 1(12) proton/atom for H(C)
    h_nucleon = GetNnucleonError(unfolded,targets)
    unfolded.Divide(unfolded,h_nucleon)
    unfolded = BinNormalize(unfolded)
    return unfolded

def GetNnucleonError(hist,ntargets):
    hist_target = hist.Clone("number_of_targets")
    hist_target.ClearAllErrorBands()
    hist_target.Reset()
    errband_name = "Target_Mass_CH"
    band_err = 0.014
    hist_target.AddVertErrorBand(errband_name,2)

    for i in range(hist_target.GetSize()):
        hist_target.SetBinContent(i,ntargets)
        hist_target.SetBinError(i,0)
        hist_target.GetVertErrorBand(errband_name).SetBinContent(i,ntargets)
        hist_target.GetVertErrorBand(errband_name).GetHist(0).SetBinContent(i,ntargets*(1-band_err))
        hist_target.GetVertErrorBand(errband_name).GetHist(1).SetBinContent(i,ntargets*(1+band_err))
    hist_target.AddMissingErrorBandsAndFillWithCV(hist)
    print("Target Normalization: {:.4E},{:.4E}".format(ntargets,ntargets*0.014))
    return hist_target

def DivideFlux(unfolded,is_mc):
    nuPDGs = [-14,14,-12,12]
    nuPDGs.remove(AnaNuPDG)
    print(nuPDGs)
    print("I'm here")
    frw5a= PlotUtils.flux_reweighter(FLUX5a,AnaNuPDG,USE_NUE_CONSTRAINT) #playlist is dummy for now
    frw6a= PlotUtils.flux_reweighter(FLUX6a,AnaNuPDG,USE_NUE_CONSTRAINT)
    flux5a = frw5a.GetIntegratedFluxReweighted(AnaNuPDG,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
    flux6a = frw6a.GetIntegratedFluxReweighted(AnaNuPDG,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
    flux5a.PopVertErrorBand("Flux_BeamFocus")
    flux6a.PopVertErrorBand("Flux_BeamFocus")
    flux5a.PopVertErrorBand("ppfx1_Total")
    flux6a.PopVertErrorBand("ppfx1_Total")
    flux5a.Scale(frac5A)    
    flux5a.Scale(1e-4 *(mc_pot if is_mc else data_pot)) #change unit to nu/cm^2
    flux6a.Scale(frac6A)    
    flux6a.Scale(1e-4 * (mc_pot if is_mc else data_pot)) #change unit to nu/cm^2    
    
    for pdg in nuPDGs:
        frw5a= PlotUtils.flux_reweighter(FLUX5a,pdg,USE_NUE_CONSTRAINT) #playlist is dummy for now
        addedflux5a = frw5a.GetIntegratedFluxReweighted(pdg,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
        addedflux5a.Scale(frac5A)
        addedflux5a.Scale(1e-4 *(mc_pot if is_mc else data_pot)) #change unit to nu/cm^2
        flux5a.Add(addedflux5a)
        frw6a= PlotUtils.flux_reweighter(FLUX6a,pdg,USE_NUE_CONSTRAINT) #playlist is dummy for now
        addedflux6a = frw6a.GetIntegratedFluxReweighted(pdg,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
        addedflux6a.Scale(frac6A)
        addedflux6a.Scale(1e-4 *(mc_pot if is_mc else data_pot)) #change unit to nu/cm^2
        flux6a.Add(addedflux6a)
    flux5a.Add(flux6a)
        
    print("flux5a=",flux5a.Print("All"))
    print ("Flux Normalization: {:.4E},{:.4E}".format(flux5a.GetBinContent(1,1),flux5a.GetTotalError(False).GetBinContent(1,1)))
    unfolded.Divide(unfolded,flux5a)
def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None

def BinNormalize(hist):
    bin_sizes = hist.Clone()
    bin_sizes.ClearAllErrorBands()
    bin_sizes.Reset()
    
    xaxis = hist.GetXaxis()
    for i in range(hist.GetNbinsX()):
        bin_sizes.SetBinContent(i+1,xaxis.GetBinLowEdge(i+2)-xaxis.GetBinLowEdge(i+1))
    bin_sizes.AddMissingErrorBandsAndFillWithCV(hist)
    #bin_sizes.Print("all")
    hist.Divide(hist,bin_sizes)
    return hist

def Integral(hist):
    integral = 0
    errors = np.zeros(hist.GetNbinsX())
    for i in range(hist.GetNbinsX()):
        integral += hist.GetBinContent(i+1)*(hist.GetXaxis().GetBinLowEdge(i+2)-hist.GetXaxis().GetBinLowEdge(i+1))
        errors[i] = hist.GetBinError(i+1)*(hist.GetXaxis().GetBinLowEdge(i+2)-hist.GetXaxis().GetBinLowEdge(i+1))
    error = np.sqrt(np.sum(np.square(errors)))
    return integral,error

def TotalXSecPrinter(stat_hist,sys_hist,name):
    integral,stat_error = Integral(stat_hist)
    sys_integral,total_error = Integral(sys_hist)
    print("The integral of the",name,"cross section is",integral,"the stat error is",stat_error,"the sys error is",np.sqrt(total_error**2 - stat_error**2))


playlist = AnalysisConfig.playlist
unfoldedFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/"+str(playlist)+"_unfolded.root","READ")
pion_matrix_hist = unfoldedFile.pionE_recovered.Clone()
mcTruthFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(playlist)+"_collab1_fspline.root","READ")
pion_sim_hist = mcTruthFile.PiZeroE_NCCohPi0
mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
mc_file = ROOT.TFile.Open(mc_path)
data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)

pot_scale = GetPOTScale(data_path,mc_path)

pion_sim_hist.Scale(pot_scale)


efficiencyFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/"+str(playlist)+"_calculatedEfficiency.root","READ")
pion_efficiency = efficiencyFile.pion_efficiency.Clone()
unfoldedFile.Close()
efficiencyFile.Close()

data_pot = getPOTFromFile("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_collab1_fspline.root")
mc_pot = getPOTFromFile("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(playlist)+"_collab1_fspline.root")

pot_5A = getPOTFromFile("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me5A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00123000_Playlist.root")
print("mc_pot=",mc_pot)
print("data_pot=",data_pot)
print("5Apot=",pot_5A)
frac5A = pot_5A/mc_pot #if is_mc else data_pot)
print("frac5A=",frac5A)
frac6A = 1- frac5A
print("frac6A=",frac6A)


pion_matrix_xsec = GetXSectionHistogram(pion_matrix_hist,pion_efficiency,False)
pion_sim_xsec = GetXSectionHistogram(pion_sim_hist,pion_efficiency,False)


outPath = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/"+str(playlist)+"_xsec.root"
outFile = ROOT.TFile.Open(outPath,"RECREATE")
pion_matrix_xsec.Write()
pion_sim_xsec.Write()
outFile.Close()

pion_canvas = ROOT.TCanvas("c2","c2",1600,1200)

pion_recovered_err = pion_matrix_xsec.GetCVHistoWithError()
TotalXSecPrinter(pion_matrix_xsec,pion_recovered_err,"Pion Matrix")

pion_recovered_err.SetMarkerStyle(21)
pion_recovered_err.SetMarkerSize(1.5)
pion_recovered_err.SetMarkerColor(4)
pion_recovered_err.SetLineColor(4)
#pion_recovered_err.SetMaximum(max([pion_matrix_xsec.GetMaximum(),pion_sim_xsec.GetMaximum()])*1.25)
#pion_recovered_err.SetMinimum(min([pion_matrix_xsec.GetMaximum(),pion_sim_xsec.GetMinimum()])/3)
#pion_recovered_err.Scale(1,"width")
pion_recovered_err.GetXaxis().SetTitle("Pion Energy (GeV)")
pion_recovered_err.GetXaxis().SetTitleFont(62)
pion_recovered_err.GetXaxis().SetTitleSize(0.045)
pion_recovered_err.GetYaxis().SetTitle("d#sigma/dE_{#pi^{0}} (cm^{2}/GeV/Carbon)")
pion_recovered_err.GetYaxis().SetTitleFont(62)
pion_recovered_err.GetYaxis().SetTitleSize(0.045)
pion_recovered_err.SetTitle("Unfolded data")
pion_recovered_err.SetLabelFont(62,"xyz")
pion_recovered_err.SetLineWidth(4)

pion_sim_err = pion_sim_xsec.GetCVHistoWithError()
pion_sim_err.SetMarkerStyle(20)
pion_sim_err.SetMarkerSize(2.5)
pion_sim_err.SetMarkerColor(8)
pion_sim_err.SetLineColor(8)
#pion_sim_err.Scale(1,"width")
pion_sim_err.SetTitle("Simulation")
pion_sim_err.SetLineWidth(4)
pion_recovered_err.SetMaximum(1.3e-39)
pion_recovered_err.Draw("MIN0 E0")
pion_sim_err.Draw("SAME E0")
pion_legend = ROOT.gPad.BuildLegend(0.6,0.7,0.8,0.85)
pion_legend.SetTextFont(62)
pion_legend.SetTextSize(0.035)
pion_legend.SetBorderSize(0)
pion_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_xsec"+str(playlist)+".png","png")

pion_canvas.SetLogy(1)
pion_recovered_err.GetYaxis().SetRangeUser(1e-43,3e-39)
#pion_recovered_err.GetYaxis().SetRangeUser(pion_sim_err.GetMinimum()/3,pion_recovered_err.GetMaximum()*3)
pion_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_xsec_log"+str(playlist)+".png","png")



X_Sec_error_canvas = ROOT.TCanvas("c1","c1",1600,1200)
PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
PlotTools.MakeErrPlot(pion_matrix_xsec)
PlotTools.MNVPLOTTER.axis_maximum = 0.2
X_Sec_error_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/X_Sec_error"+str(playlist)+".png","png")


