import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
from config.SystematicsConfig import GENIE_ERROR_GROUPS,CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS
from tools import Utilities,PlotTools

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.15)



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

def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None
def format_hist(hist, title=None):
    """Apply consistent styling to histograms."""
    hist.SetTitle(title or hist.GetTitle())
    hist.SetTitleFont(62)
    hist.SetLabelFont(62, "xyz")
    hist.SetMarkerStyle(21)
    hist.SetMarkerSize(2)
    hist.SetLineWidth(4)
    hist.SetStats(0)

    for axis in (hist.GetXaxis(), hist.GetYaxis()):
        axis.SetTitleFont(62)
        axis.SetTitleSize(0.045)

    return hist

def save_canvas(canvas, path, logy=False):
    if logy:
        canvas.SetLogy()
    canvas.Print(str(path), "png")


playlist = AnalysisConfig.playlist
#mc_file = "/exp/minerva/app/users/nalex/nccoh/AL9-NCCohPi0/ROOTFILES/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
mc_file = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(playlist)+"_final_weightedmc.root"
rootFile = ROOT.TFile.Open(mc_file,"READ")

data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_collab1_fspline.root"
pot_scale = GetPOTScale(data_path,mc_file)

print("pot_scale=",pot_scale)

Eel = rootFile.Eel.Clone()
pion_numerator = rootFile.PiZeroE_NCCohPi0.Clone()   # PiZeroE_NCCohPi0
pion_denominator = rootFile.true_pizeroE.Clone()


pion_numerator.AddMissingErrorBandsAndFillWithCV(Eel)
pion_denominator.AddMissingErrorBandsAndFillWithCV(Eel)

pion_numerator.Scale(pot_scale)
pion_denominator.Scale(pot_scale)


pion_efficiency = pion_numerator.Clone()

pion_efficiency.Divide(pion_numerator,pion_denominator,1.0,1.0)

pion_efficiency.GetYaxis().SetTitle("Efficiency")
pion_efficiency.SetName("pion_efficiency")
pion_efficiency.GetXaxis().SetTitle("Pion Energy (GeV)")
format_hist(pion_efficiency, "pion_efficiency")



efficiencyFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/"+str(playlist)+"_calculatedEfficiency.root","RECREATE")
pion_efficiency.Write()
efficiencyFile.Close()


pion_canvas = ROOT.TCanvas("c2","c2",1600,1200)
pion_efficiency.Draw("MINO")  #MIN0
save_canvas(pion_canvas,"/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_efficiency_"+str(playlist)+".png")



Efficiency_numerator_denominator_canvas = ROOT.TCanvas("c1","c1",1600,1200)
numerator = pion_numerator.GetCVHistoWithError()
denominator = pion_denominator.GetCVHistoWithError()
format_hist(numerator, "Numerator")
format_hist(denominator, "Denominator")
numerator.Scale(1, "width")
denominator.Scale(1, "width")
denominator.SetLineColor(ROOT.kRed)

denominator.Draw("E1")
numerator.Draw("E1 SAME")
ROOT.gPad.BuildLegend(0.5,0.5,0.8,0.8)
save_canvas(Efficiency_numerator_denominator_canvas,"/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/BWN_Efficiency_numerator_denominator_"+str(playlist)+"_.png")

Efficiency_numerator_denominator_canvas.SetLogy()
save_canvas(Efficiency_numerator_denominator_canvas,"/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/BWN_Efficiency_numerator_denominator_"+str(playlist)+"_log_.png")

error_canvas = ROOT.TCanvas("c1","c1",1600,1200)
PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)

def make_error_plot(hist, output):
    canvas = ROOT.TCanvas("", "", 1600, 1200)
    PlotTools.MakeErrPlot(hist)
    save_canvas(canvas, output)

make_error_plot(pion_numerator, f"/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/numerator_error_{playlist}.png")
make_error_plot(pion_denominator, f"/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/denominator_error_{playlist}.png")
make_error_plot(pion_efficiency, f"/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/Efficiency_Error_{playlist}.png")


efficiencyFile.Close()
rootFile.Close()
