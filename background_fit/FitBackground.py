import ROOT
import PlotUtils
import math
from array import array
import numpy as np
import os

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from config.PlotConfig import ELECTRON_ENERGY_BINNING
from tools.Utilities import getPOT
from tools import Utilities,PlotTools
from config.SignalDef import SIGNAL_DEFINATION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS
MNVPLOTTER = PlotUtils.MnvPlotter()
ROOT.TH1.AddDirectory(True) # If you make this false, the error messages go away but the fitter just stops working because it depends on weird ROOT behavior with the POT scaling applying in ways I don't understand but it works as is
# If you really want to make this false, you gotta go back and rethink how the fitter grabs the root files and turn on the POT scaling. If you don't fix how it grabs ROOT files after making it false, this code will take up all the RAM on your system. 

def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None

def GetPOTScale(data_path,mc_path): # This is to correct for the fact that the mc and data have different simulated numbers of target protons in the neutrino maker
    pots = [None,None]
    for i,t in enumerate(["data","mc"]):
        path = [data_path,mc_path][i]
        try:
            pots[i]= getPOTFromFile(path) #or getPOT(playlist,t,ntuple_tag)
            print("I'm here")
            print(pots)
        except KeyError:
            pots[i]=None
    pot_scale = pots[0]/pots[1] if pots.count(None) == 0 else 1.0
    print("POT Scale:",pot_scale, "Data POT:",pots[0],"Sim POT:",pots[1])
    return pot_scale

def AddTwoHistograms(band1,band2,pot_scale): # Adds band2 to band 1 and adds the pot scaled errors in quadrature
    errors = []
    for i in range(band1.GetSize()): # Adds errors in quadrature
        #POTScaling scales the bin errors also, so there is no need to accoun for the pot_scale again. Just sum in quadrature
        errors.append(math.sqrt(band1.GetBinError(i)**2 + band2.GetBinError(i)**2))
    band1.Add(band2,1) # Adds band2 to band1
    for i in range(band1.GetSize()): # Sets the error of signal to be what was calculated before
        band1.SetBinError(i,errors[i])
    return band1


def SubtractPoissonHistograms(signal,background_band,pot_scale,first_subtraction = False): # Sets error of signal where signal - background_band is being calculat
    errors = []
    for i in range(signal.GetSize()): 
        #POTScaling scales the bin errors also, so there is no need to accoun for the pot_scale again. Just sum in quadrature
        errors.append(math.sqrt(signal.GetBinError(i)**2 + background_band.GetBinError(i)**2))
    signal.Add(background_band,-1) # Subtracts the band from signal
    for i in range(signal.GetSize()): # Sets the error of signal to be what was calculated before
        signal.SetBinError(i,errors[i])
    return signal


def PrimeExtractHists(data_hists, mc_hists,sideband):
    pot_scale = mc_hists.pot_scale
    data = data_hists.GetHist().Clone()
    bands = ["NCCohPi0","NCDiff","CCQElike","notCCQElike","NCPi0","CCPi0","lowNCPi0","NCDiff_NPS"]
    #coherent = mc_hists.hists["NCCohPi0"].Clone()
    coherent = AddTwoHistograms(mc_hists.hists["NCCohPi0"].Clone(),mc_hists.hists["lowNCPi0"].Clone(),pot_scale)
    #ncdiff = mc_hists.hists["NCDiff"].Clone()
    ncdiff = AddTwoHistograms(mc_hists.hists["NCDiff"].Clone(),mc_hists.hists["NCDiff_NPS"].Clone(),pot_scale)
    #ncdiff.Scale(6.0)
    nu_e = AddTwoHistograms(mc_hists.hists["CCQElike"].Clone(),mc_hists.hists["notCCQElike"].Clone(),pot_scale)
    CCpiZero = mc_hists.hists["CCPi0"].Clone()
    NCpiZero = mc_hists.hists["NCPi0"].Clone()
    total_mc = mc_hists.hists["Total"].Clone()
    
    for i,band in enumerate(bands): # Here, each of the subbands that are being calculated are removed from the total
        total_mc = SubtractPoissonHistograms(total_mc,mc_hists.hists[band],pot_scale,not bool(i)) 
    augment = SubtractPoissonHistograms(data,total_mc,pot_scale) # This subtracts the LHS const from the data to give the RHS of the equation
    #print("augment=",augment)
    return [coherent,nu_e,NCpiZero,CCpiZero,ncdiff,augment]

def ErrExtractHists(data_hists, mc_hists, sideband, err_type, universe_num):
    pot_scale = mc_hists.pot_scale
    #print("err_pot=",pot_scale)
    data = data_hists.GetHist().Clone()
    coherent = AddTwoHistograms(mc_hists.hists["NCCohPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["lowNCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),pot_scale)
    #ncdiff = mc_hists.hists["NCDiff"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    ncdiff= AddTwoHistograms(mc_hists.hists["NCDiff"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["NCDiff_NPS"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),pot_scale)
    #ncdiff.Scale(6)
    nu_e = AddTwoHistograms(mc_hists.hists["CCQElike"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["notCCQElike"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),pot_scale)
    CCpiZero= mc_hists.hists["CCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone() #
    NCpiZero= mc_hists.hists["NCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone() #
    total_mc = mc_hists.hists["Total"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    bands = ["NCCohPi0","NCDiff","CCQElike","notCCQElike","NCPi0","CCPi0","lowNCPi0","NCDiff_NPS"]
    for i,band in enumerate(bands): # Here, each of the subbands that are being calculated are removed from the total
        total_mc = SubtractPoissonHistograms(total_mc,mc_hists.hists[band].GetVertErrorBand(err_type).GetHist(universe_num),pot_scale,not bool(i))
    other = PlotUtils.MnvH1D(total_mc) # I have no idea why the add command can't do this itself but whatever
    augment = SubtractPoissonHistograms(data,other,pot_scale) # This subtracts the LHS const from the data to give the RHS of the equation
    return [coherent,nu_e,NCpiZero,CCpiZero,ncdiff,augment]

def HistListtoArray(histlist):
#    """Convert list of ROOT histograms to a numpy array"""
    histarray = np.array([[hist.GetBinContent(i+1) for i in range(hist.GetSize()-1)]
                          for hist in histlist])
    if len(histlist) == 1:
        histarray = histarray.flatten()
    return histarray


def SetupMatrices(dataarray):
    dataarray = np.array(dataarray)
    num_matrices = dataarray.shape[2]

    matrices = dataarray[:, 0:5, :].transpose(2, 0, 1)   # (num_matrices, 5, 5)
    vectors = dataarray[:, 5, :].T                        # (num_matrices, 5)
    vectors -= np.sum(matrices, axis=2)

    return matrices, vectors


def DoLinearAlgebra(inversesArray, vectorArray):
    solutions = np.matmul(inversesArray, vectorArray[..., None]).squeeze(-1)
    return 1.0 + solutions


def NormalInverse(matrices):
    matrices = np.asarray(matrices)
    return np.linalg.inv(matrices)

def MPInverse(matrices):
    matrices = np.asarray(matrices)
    N = matrices.shape[0]
    U, S, Vt = np.linalg.svd(matrices)
    
        #RHC
    discarded = np.array([[False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,True],
                         [False,False,False,False,True],
                         [False,False,False, True,True],
                         [False,False,False,True,True]],dtype = bool) #Picks singular values to be discarded
    """
    #FHC
    discarded = np.array([[False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,False],
                         [False,False,False,False,True],
                         [False,False,False,True,True],
                         [False,False,False,False,True],
                         [False,False,False,False,True]],dtype = bool) #Picks singular values to be discarded
    """
    if discarded.shape != S.shape:
        raise ValueError(f"Discard mask shape {discarded.shape} "
                         f"does not match singular values shape {S.shape}")
    S_inv = np.zeros_like(S)
    keep = ~discarded
    S_inv[keep] = 1.0 / S[keep]
    
    Sigma_plus = np.zeros_like(matrices)
    idx = np.arange(S.shape[1])
    Sigma_plus[:, idx, idx] = S_inv
     
    mpinverse = Vt.transpose(0,2,1) @ Sigma_plus @ U.transpose(0,2,1)
    return mpinverse 

    
def GetBins(hist,print_results = False): # If you want to print out the binds in a hist, use this
    bins = []
    for i in range(hist.GetSize()):
        bins.append(hist.GetBinContent(i))
    if print_results: print(hist,bins)
    return bins

def GetPreweightedROOTFiles(weighting_playlist):
    mc_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(weighting_playlist)+"_eavail_weightedmc.root"
    mc_file = ROOT.TFile.Open(mc_path)
    data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(weighting_playlist)+"_collab1_fspline.root"
    data_file = ROOT.TFile.Open(data_path)
    pot_scale = GetPOTScale(data_path,mc_path)
    return mc_file,data_file,pot_scale

def GetCVWeights(mc_file,data_file,pot_scale):
    sidebands = ["Signal","Electron_Neutrino","Pi0","Pi0_muonexits","Diffractive"]
    # First going through the weights for the prime universe because I want to have this be its own thing
    sideband_list = []
    for sideband in sidebands:
        mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale) # Gets the two electron energy histograms
        dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False) # I'm using HistHolders to access the backgrounds
        mcEel_hists.POTScale(False)
        dataEel_hist.POTScale(False)
        listofhists = PrimeExtractHists(dataEel_hist, mcEel_hists,sideband)
        arrayofhists = HistListtoArray(listofhists)
        sideband_list.append(arrayofhists)
        del mcEel_hists, dataEel_hist
    bigArray = np.array(sideband_list)
    matrices,vectors = SetupMatrices(bigArray)
    Chi2test(matrices,vectors)
    invmatrices = MPInverse(matrices)
    primeweights = DoLinearAlgebra(invmatrices,vectors)
    WasIRight(matrices,primeweights,vectors)
    
    cov_matrix_data = [np.diag(np.abs(vector)) for vector in vectors]
    cov_matrix_data = np.array(cov_matrix_data)
    weight_covariances = [inv @ cov @ inv.T for inv, cov in zip(invmatrices, cov_matrix_data)]
    reco_weight_errors = [np.sqrt(np.diag(cov)) for cov in weight_covariances]
    return primeweights,reco_weight_errors


def Chi2test(matrices,vectors):
    for i in range(len(matrices)):
        D  = vectors[i]
        err = np.sqrt(D)
        U, s, Vt = np.linalg.svd(matrices[i], full_matrices=False)

        N = len(D)
        print("eigenvalues=",s)
        for k in range(1, len(s)+1):

            s_inv = np.zeros_like(s)
            s_inv[:k] = 1.0 / s[:k]
            print("s_inv=",s_inv)
            MC_pinv_k = Vt.T @ np.diag(s_inv) @ U.T
            w_k = MC_pinv_k @ D
            Dpred = matrices[i] @ w_k

            chi2 = np.sum(((D - Dpred) / err)**2)
            ndf = N - k

            print("i=",i,"k=",k, "chi2=",chi2, "chi2/ndf=",chi2 / ndf if ndf > 0 else None)

def WasIRight(matrices,weights,vectors):
    for i in range(len(matrices)):
        print("MC: ",np.dot(matrices[i],weights[i]-np.ones(5)))
        print("Data: ",vectors[i])


def RunAlexWeighter():
    #print("Starting Alex's Reweighting Weighter")
    mc_file,data_file,pot_scale = GetPreweightedROOTFiles(None)
    cvweights = GetCVWeights(mc_file,data_file,pot_scale)
    universe_weights_dict = GetSystematicUniverseWeights(mc_file,data_file,pot_scale)
    #print("Weighter Ran Successfully")
    return cvweights,universe_weights_dict




def GetSystematicUniverseWeights(mc_file,data_file,pot_scale):
    # Now for the Systematic Universes
    sidebands = ["Signal","Electron_Neutrino","Pi0","Pi0_muonexits","Diffractive"]
    err_weights_dict = {}
    err_weights_stat_err_dict = {}
    universe_list = [("Flux",200),("GENIE_AGKYxF1pi",2),("GENIE_AhtBY",2),("GENIE_BhtBY",2),("GENIE_CCQEPauliSupViaKF",2),("GENIE_CV1uBY",2),("GENIE_CV2uBY",2),("GENIE_EtaNCEL",2),("GENIE_FrAbs_N",2),("GENIE_FrAbs_pi",2),("GENIE_FrCEx_N",2),("GENIE_FrCEx_pi",2),("GENIE_FrElas_N",2),("GENIE_FrElas_pi",2),("GENIE_FrInel_N",2),("GENIE_FrPiProd_N",2),("GENIE_FrPiProd_pi",2),("GENIE_MFP_N",2),("GENIE_MFP_pi",2),("GENIE_MaNCEL",2),("GENIE_MaRES",2),("GENIE_MvRES",2),("GENIE_NormCCRES",2),("GENIE_NormDISCC",2),("GENIE_NormNCRES",2),("GENIE_RDecBR1gamma",2),("GENIE_Rvn2pi",2),("GENIE_Rvp2pi",2),("GENIE_Theta_Delta2Npi",2),("GENIE_VecFFCCQEshape",2),("GENIE_Rvn1pi",2),("GENIE_Rvp1pi",2),("Low_Recoil_2p2h_Tune",3),("RPA_HighQ2",2),("RPA_LowQ2",2),("LowQ2Pi",2),("MK_model",1),("GENIE_MaCCQE",2),("fsi_weight",3),("SuSA_Valencia_Weight",1),("GEANT_Proton",2),("GEANT_Neutron",2),("GEANT_Pion",2),("Target_Mass_CH",2),("Muon_Energy_MINERvA",2),("Muon_Energy_MINOS",2),("beam_angle",4),("response_p",2),("response_meson",2),("response_em",2),("response_other",2),("response_xtalk",2),("Leakage_Uncertainty",2),("Eavail_Sys",2)]
    for error_type,error_universes in universe_list:
        err_type_weights = []
        err_type_weight_errors = []
        for universe_num in range(error_universes):
            sideband_list = []
            for sideband in sidebands:
                mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale)
                dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False)
                mcEel_hists.POTScale(False)   #If i turn these on my vertical uncertainity goes very very very high??
                dataEel_hist.POTScale(False)
                listofhists = ErrExtractHists(dataEel_hist, mcEel_hists, sideband, error_type, universe_num)
                arrayofhists = HistListtoArray(listofhists)
                sideband_list.append(arrayofhists)
                del mcEel_hists, dataEel_hist
            bigArray = np.array(sideband_list)
            matrices,vectors = SetupMatrices(bigArray)
            invmatrices = MPInverse(matrices)
            #invmatrices = NormalInverse(matrices)
            weights = DoLinearAlgebra(invmatrices, vectors)
            err_type_weights.append(weights)
            
            #Calculate stat uncertainity
            cov_matrix_data = [np.diag(np.abs(vector)) for vector in vectors]
            weight_covariances = [inv @ cov @ inv.T for inv, cov in zip(invmatrices, cov_matrix_data)]
            reco_y_error = [np.sqrt(np.diag(cov)) for cov in weight_covariances]
            err_type_weight_errors.append(reco_y_error)

        err_weights_dict[error_type] = np.array(err_type_weights)
        err_weights_stat_err_dict[error_type] = np.array(err_type_weight_errors)
    return err_weights_dict,err_weights_stat_err_dict


def plot(hist):
    c1_canvas = ROOT.TCanvas(f"{hist.GetName()}_cv","",1600,1200)
    h_err = hist.GetCVHistoWithError()
    h_err.SetMarkerStyle(21)
    h_err.SetMarkerColor(4)
    h_err.SetLineColor(9)
    h_err.SetMarkerSize(2)
    h_err.SetLineWidth(3)
    h_err.GetXaxis().SetTitle("Lepton Energy (GeV)")
    h_err.SetStats(0)
    h_err.Draw("MIN0 E1")
    c1_canvas.Print(f"/exp/minerva/data/users/nalex/nu_e2/results/NCDiff/{hist.GetName()}_err_2_{preweighted_playlist}.png","png")

    c1_error_canvas = ROOT.TCanvas("c1","c1",1600,1200)
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
    PlotTools.MakeErrPlot(hist)
    c1_error_canvas.Print(f"/exp/minerva/data/users/nalex/nu_e2/results/NCDiff/{hist.GetName()}_wt_error_2_{preweighted_playlist}.png","png")


if __name__ == "__main__":
    preweighted_playlist = AnalysisConfig.playlist
    mc_file,data_file,pot_scale = GetPreweightedROOTFiles(preweighted_playlist)
    cvweights,reco_y_error = GetCVWeights(mc_file,data_file,pot_scale)
    print("cvweights=",cvweights)
    universe_weights_dict,sys_stat_error = GetSystematicUniverseWeights(mc_file,data_file,pot_scale)
    np.set_printoptions(precision=6,suppress = True,threshold = np.inf)
    np.savetxt("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/weights.txt",cvweights,delimiter=",")
    dict_str = str(universe_weights_dict)
    with open("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/systematic_universes_weights.txt", "w") as file:
        file.write(dict_str)
    Eel_hist = mc_file.Eel.Clone()
    

    mc_file.Close()
    data_file.Close()

    weightFile = ROOT.TFile.Open("/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_"+str(preweighted_playlist)+"_weights_with_staterror.root","RECREATE")

    hist_defs = {
    "coherent": "Weights on Coherent Events",
    "electron_neutrino": "Weights on Electron Neutrino Events",
    "NCpizero": "Weights on NCPiZero Events",
    "CCpizero": "Weights on CCPiZero Events",
    "ncdiff": "Weights on NCDiff Events"
    }
    hists = {}
    for name, title in hist_defs.items():
        hists[name] = PlotUtils.MnvH1D(name, title, len(ELECTRON_ENERGY_BINNING)-1, array("d", ELECTRON_ENERGY_BINNING))

    # Fill CV and bin errors
    for i in range(len(ELECTRON_ENERGY_BINNING)-1):
        for j, name in enumerate(hist_defs.keys()):
            hists[name].SetBinContent(i+1, cvweights[i, j])
            hists[name].SetBinError(i+1, reco_y_error[i][j])

    # Add missing error bands using Eel_hist
    for hist in hists.values():
        hist.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    # Fill universe weights and errors
    for hist in hists.values():
        for error_type in hist.GetErrorBandNames():
            error_type = str(error_type)
            n_hists = hist.GetVertErrorBand(error_type).GetNHists()
            for k in range(n_hists):
                vert_hist = hist.GetVertErrorBand(error_type).GetHist(k)
                for i in range(len(ELECTRON_ENERGY_BINNING)-1):
                    vert_hist.SetBinContent(i+1, universe_weights_dict[error_type][k][i, list(hist_defs.keys()).index(hist.GetName())])
                    vert_hist.SetBinError(i+1, sys_stat_error[error_type][k][i, list(hist_defs.keys()).index(hist.GetName())])
    for hist in hists.values():
        hist.GetXaxis().SetTitle("Lepton Energy")
        hist.SetMarkerStyle(21)
        hist.Write()
    plot(hists["coherent"])
    plot(hists["ncdiff"])
    plot(hists["NCpizero"])
    plot(hists["CCpizero"])
    plot(hists["electron_neutrino"])
    #del weightFile, hists, Eel_hist
    print("Done")

