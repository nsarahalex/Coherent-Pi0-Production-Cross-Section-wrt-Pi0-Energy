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

mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(weighting_playlist)+"_eavail_weightedmc.root"
mc_file = ROOT.TFile.Open(mc_path)
data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(weighting_playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)

pot_scale = Utilities.GetPOTScale(data_path,mc_path)
print("pot_scale=",pot_scale)
weight_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(weighting_playlist)+"_weights_with_staterror.root"
weight_file = ROOT.TFile.Open(weight_path)

#mc_file.Eel_Diffractive_CCPi0.GetVertErrorBand("Flux").GetHist(1).Print("All")


output_file = ROOT.TFile("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(weighting_playlist)+"_final_weightedmc.root", "RECREATE")

all_mc_keys = [key.GetName() for key in mc_file.GetListOfKeys()]
mc_keys = [key for key in all_mc_keys if re.match(r"^Eel", key)]
weights_keys = [key.GetName() for key in weight_file.GetListOfKeys()]
match_rules = {
    "coherent":r"NCCohPi0|lowNCPi0",
    "ncdiff": r"NCDiff",
    "CCpizero": r"CCPi0",
    "NCpizero": r"(?<!low)NCPi0",
    "electron_neutrino":r"notCCQElike|CCQElike",
}


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



#Step 1: Copy Everything from mc_file to output file:

for key in mc_file.GetListOfKeys():
    name = key.GetName()
    if name == "Meta":
        continue
    obj = key.ReadObj()
    output_file.cd()
    obj.Write()

#Step2: Weigh the Eel
hist_pairs = []
for weight_key, mc_string in match_rules.items():
    for mc_key in mc_keys:
        if re.search(mc_string,mc_key):
            result_name = f"{mc_key}"
            hist_pairs.append((mc_key,weight_key,result_name))

for mc_hist_name, weight_hist_name, result in hist_pairs:
     weight_hist = weight_file.Get(weight_hist_name)
     result_hist = mc_file.Get(mc_hist_name).Clone()
     weight_hist.AddMissingErrorBandsAndFillWithCV(result_hist)
     if isinstance(result_hist, cppyy.gbl.PlotUtils.MnvH1D):
         weight_hist = ROOT.PlotUtils.MnvH1D(weight_hist)
         result_hist = ROOT.PlotUtils.MnvH1D(result_hist)
         result_hist.Multiply(result_hist,weight_hist)
         if weight_hist_name == "ncdiff":
            result_hist.Scale(6.0)
         output_file.cd()
         result_hist.Write(mc_hist_name, ROOT.TObject.kOverwrite)

#Step3: Resum Total
sides = ["Signal","Diffractive","Electron_Neutrino","Pi0","Pi0_muonexits"]
for side in sides:
    Cone_hist = HistHolder("Lepton Energy",output_file,side,is_mc = True,pot = "pot_scale")
    Cone_hist.ResumTotal()
    output_file.cd()
    Cone_hist.GetHist().Write(Cone_hist.GetHist().GetName(),ROOT.TObject.kOverwrite)

#Step4: Copy Meta Tree
CopyCleanMetaTreeFromFile(mc_file, output_file)
print("Done")
#output_file.Close()
#weight_file.Close()
