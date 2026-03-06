



import ROOT
import PlotUtils
from config.SystematicsConfig import USE_NUE_CONSTRAINT
from config.AnalysisConfig import AnalysisConfig
from config.CutConfig import NEUTRINO_ENERGY_RANGE


playlist = AnalysisConfig.playlist
energy_binning = [1.5, 2, 2.5, 3.25, 4, 5, 6, 7, 10, 14, 24]
def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None



data_pot = getPOTFromFile("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_collab1_fspline.root")
mc_pot = getPOTFromFile("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(playlist)+"_collab1_fspline.root")

pot_5A = getPOTFromFile("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me5A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00123000_Playlist.root")
print("5Apot=",pot_5A)
frac5A = pot_5A/mc_pot #if is_mc else data_pot)
print("frac5A=",frac5A)
frac6A = 1- frac5A
print("frac6A=",frac6A)



#FLUX = [("minervame1d1m1nweightedave",1)]

FLUX = [("minervame5a",frac5A),("minervame6a",frac6A)]
def flux_val(FLUX,pdg,USE_NUE_CONSTRAINT):
    frw = PlotUtils.flux_reweighter(FLUX, pdg, USE_NUE_CONSTRAINT)
    return frw

energy_binning = [1.5, 2, 2.5, 3.25, 4, 5, 6, 7, 10, 14, 24]

flux_by_pdg = {}


PDG =  [-14,14,-12,12]


unfoldedFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/"+str(playlist)+"_unfolded.root","READ")
unfolded = unfoldedFile.pionE_recovered.Clone()

for pdg in PDG:
    total_flux_hist = None
    for (fluxs,scalefactors) in FLUX:
        frw = flux_val(fluxs,pdg,USE_NUE_CONSTRAINT)
        hist = frw.GetIntegratedFluxReweighted(pdg, unfolded, NEUTRINO_ENERGY_RANGE[0], NEUTRINO_ENERGY_RANGE[1], False)
        hist.Scale(scalefactors)
    
        if total_flux_hist is None:
            total_flux_hist = hist.Clone(f"total_flux_{pdg}")
        else:
            total_flux_hist.Add(hist)
    total_flux = total_flux_hist.Integral()
    flux_by_pdg[pdg] = []
    
    for i in range(len(energy_binning) - 1):
        low, high = energy_binning[i], energy_binning[i+1]
        
        bin_low =  total_flux_hist.FindBin(low)
        bin_high =  total_flux_hist.FindBin(high)
        
        bin_val =  total_flux_hist.Integral(bin_low, bin_high)
        frac = bin_val/total_flux
        flux_by_pdg[pdg].append(frac)
        
    print(f"Flux contribution for PDG {pdg}: {flux_by_pdg[pdg]}")


