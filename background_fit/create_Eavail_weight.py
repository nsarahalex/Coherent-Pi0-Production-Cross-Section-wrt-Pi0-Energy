import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
import array


weighting_playlist = AnalysisConfig.playlist



#Kevin suggested using FHC value or avg value to scale. Add new systematics corresponding to stat uncertainity on weight

def MakeWeightHist(hist):
    #Make a new MnvH1D similar to hist
    nbins = hist.GetNbinsX()
    edges = [hist.GetBinLowEdge(i+1) for i in range(nbins)]
    edges.append(hist.GetBinLowEdge(nbins+1))
    edges_array = array.array('d',edges)
    weight_hist = ROOT.PlotUtils.MnvH1D("Eavail_weight","Eavail_weight",nbins,edges_array)
    val = 1.172  #1.261
    #Fill with val, set stat uncertainity to 0
    for bin_idx in range(1,nbins +1):
        weight_hist.SetBinContent(bin_idx,val)
        weight_hist.SetBinError(bin_idx, 0.0)
    #Create a new Sys Universe to include stat fluction in Eavail weight
    n_universes = 2
    weight_hist.AddVertErrorBand("Eavail_Sys",n_universes)
    new_band = weight_hist.GetVertErrorBand("Eavail_Sys")
    for bin_idx in range(1,nbins +1):
        new_band.GetHist(0).SetBinContent(bin_idx, val+0.081)
        new_band.GetHist(0).SetBinError(bin_idx, 0.0)
        new_band.GetHist(1).SetBinContent(bin_idx, val-0.081)   #0.143
        new_band.GetHist(1).SetBinError(bin_idx, 0.0)
    weight_hist.AddMissingErrorBandsAndFillWithCV(hist)
    return weight_hist

mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(weighting_playlist)+"_collab1_fspline.root"

mc_file = ROOT.TFile.Open(mc_path)


Eavail_weight_file = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(weighting_playlist)+"_Eavail_weight.root","RECREATE")

Lepton_Energy_hist = mc_file.Eel.Clone()
Eavail_weight_hist = MakeWeightHist(Lepton_Energy_hist)
Eavail_weight_file.cd()
Eavail_weight_hist.Write()
print("Done")
#Eavail_weight_file.Close()

mc_file.Close()

