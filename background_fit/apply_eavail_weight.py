import ROOT
import re
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
import cppyy
from tools import Utilities
import math
import array
from tools.PlotLibrary import HistHolder


weighting_playlist = AnalysisConfig.playlist

mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(weighting_playlist)+"_collab1_fspline.root"
mc_file = ROOT.TFile.Open(mc_path)
data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(weighting_playlist)+"_collab1_fspline.root"

pot_scale = Utilities.GetPOTScale(data_path,mc_path)
print("pot_scale=",pot_scale)



Eavail_weight_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(weighting_playlist)+"_Eavail_weight.root"
Eavail_weight_file = ROOT.TFile.Open(Eavail_weight_path)
Eavail_weight_hist = Eavail_weight_file.Eavail_weight.Clone()

#mc_file.Eel_Diffractive_CCPi0.GetVertErrorBand("Flux").GetHist(1).Print("All")


output_file = ROOT.TFile("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(weighting_playlist)+"_eavail_weightedmc.root", "RECREATE")





all_mc_keys = [key.GetName() for key in mc_file.GetListOfKeys()]
mc_keys = [key for key in all_mc_keys if re.match(r"^Eel", key)]

def CopyCleanMetaTreeFromFile(source_file, output_file):
    if not source_file.Get("Meta"):
        print("No Meta tree found!")
        return
    
    old_tree = source_file.Get("Meta")
    new_tree = ROOT.TTree("Meta", "Cleaned Meta tree")

    pot_val = array.array('d', [0.0])  # POT_Used is a double
    new_tree.Branch("POT_Used", pot_val, "POT_Used/D")

    n_entries = old_tree.GetEntries()
    old_tree.SetBranchStatus("*", 0)
    old_tree.SetBranchStatus("POT_Used", 1)
    old_tree.SetBranchAddress("POT_Used", pot_val)

    for i in range(n_entries):
        if old_tree.GetEntry(i) > 0:
            new_tree.Fill()
        else:
            print(f"Skipping corrupted entry {i}")

    output_file.cd()
    new_tree.Write()


def should_weight_pi0(hist_name):
    valid_prefixes = ["Eel_", "Eel_Electron_Neutrino_", "Eel_Diffractive_"]
    valid_suffixes = ["NCPi0", "CCPi0"]
    
    for prefix in valid_prefixes:
        if hist_name.startswith(prefix):
            for suffix in valid_suffixes:
                if hist_name == f"{prefix}{suffix}":
                    return True
    return False

#Step 1: Copy Everything from mc_file to output file:

for key in mc_file.GetListOfKeys():
    name = key.GetName()
    if name == "Meta":
        continue
    obj = key.ReadObj()
    output_file.cd()
    obj.Write()

#Step2: Weight the selected Eel ones with Eavail weight

processed_hists = set()
for mc_hist_name in mc_keys:
    mc_hist = mc_file.Get(mc_hist_name).Clone()
    mc_hist.AddMissingErrorBandsAndFillWithCV(Eavail_weight_hist)
    if isinstance(mc_hist, cppyy.gbl.PlotUtils.MnvH1D):
        mc_hist = ROOT.PlotUtils.MnvH1D(mc_hist)
        if should_weight_pi0(mc_hist_name):
            mc_hist.Multiply(mc_hist,Eavail_weight_hist)
        processed_hists.add(mc_hist_name) 
        output_file.cd()
        mc_hist.Write(mc_hist_name, ROOT.TObject.kOverwrite)

sides = ["Signal","Diffractive","Electron_Neutrino","Pi0","Pi0_muonexits"]
for side in sides:
    Cone_hist = HistHolder("Lepton Energy",output_file,side,is_mc = True,pot = "pot_scale")
    Cone_hist.ResumTotal()
    output_file.cd()
    Cone_hist.GetHist().Write(Cone_hist.GetHist().GetName(),ROOT.TObject.kOverwrite)

CopyCleanMetaTreeFromFile(mc_file, output_file)

print("Done")
#output_file.Close()
#weight_file.Close()
