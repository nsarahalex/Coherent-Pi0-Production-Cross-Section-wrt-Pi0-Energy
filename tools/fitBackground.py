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
            #print("I'm here")
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
        #errors.append(math.sqrt(math.sqrt(pot_scale)*band1.GetBinError(i)**2 + math.sqrt(pot_scale)*band2.GetBinError(i)**2))
    band1.Add(band2,1) # Adds band2 to band1
    for i in range(band1.GetSize()): # Sets the error of signal to be what was calculated before
        band1.SetBinError(i,errors[i])
    return band1

def AddThreeHistograms(band1,band2,band3,pot_scale): # Adds band2 to band 1 and adds the pot scaled errors in quadrature
    errors = []
    for i in range(band1.GetSize()): # Adds errors in quadrature
        #POTScaling scales the bin errors also, so there is no need to accoun for the pot_scale again. Just sum in quadrature
        errors.append(math.sqrt(band1.GetBinError(i)**2 + band2.GetBinError(i)**2 + band3.GetBinError(i)**2))
        #errors.append(math.sqrt(pot_scale*band1.GetBinError(i)**2 + pot_scale*band2.GetBinError(i)**2+pot_scale*band3.GetBinError(i)**2))
    #GetBins(band1,True)

    band1.Add(band2,1) # Adds band2 to band1
    band1.Add(band3,1)
    #GetBins(band1,True)
    for i in range(band1.GetSize()): # Sets the error of signal to be what was calculated before
        band1.SetBinError(i,errors[i])
    return band1


def SubtractPoissonHistograms(signal,background_band,pot_scale,first_subtraction = False): # Sets error of signal where signal - background_band is being calculat
    errors = []
    for i in range(signal.GetSize()): 
        #POTScaling scales the bin errors also, so there is no need to accoun for the pot_scale again. Just sum in quadrature
        errors.append(math.sqrt(signal.GetBinError(i)**2 + background_band.GetBinError(i)**2))
        #if first_subtraction: errors.append(math.sqrt(pot_scale*signal.GetBinError(i)**2 + pot_scale*background_band.GetBinError(i)**2)) # Scales the mc total error by the POT scale just once
        #else: errors.append(math.sqrt(signal.GetBinError(i)**2 + pot_scale*background_band.GetBinError(i)**2))
    signal.Add(background_band,-1) # Subtracts the band from signal
    for i in range(signal.GetSize()): # Sets the error of signal to be what was calculated before
        signal.SetBinError(i,errors[i])
    return signal

def Scalepizero(hist):
    ####################
    ##RHC
#    val = 1.248  # Scale factor 
#    val_err = 0.260  # Uncertainty on scale factor
    ####################
    #FHC
    val = 1.274  # Scale factor
    val_err = 0.121  # Uncertainty on scale factor
    #####################
    nbins = hist.GetNbinsX()
    for i in range(1, nbins + 1):  # Bins go from 1 to nbins
        content = hist.GetBinContent(i)
        error = hist.GetBinError(i)
        new_content = content * val
        new_error = math.sqrt((content * val_err) ** 2 + (val * error) ** 2)
        hist.SetBinContent(i, new_content)
        hist.SetBinError(i, new_error)
    return hist

def PrimeExtractHists(data_hists, mc_hists,sideband):
    pot_scale = mc_hists.pot_scale
    data = data_hists.GetHist().Clone()
    bands = ["NCCohPi0","NCDiff","CCQElike","notCCQElike","NCPi0","CCPi0","lowNCPi0"]
    #coherent = mc_hists.hists["NCCohPi0"].Clone()
    coherent = AddTwoHistograms(mc_hists.hists["NCCohPi0"].Clone(),mc_hists.hists["lowNCPi0"].Clone(),pot_scale)
    ncdiff = mc_hists.hists["NCDiff"].Clone()
    nu_e = AddTwoHistograms(mc_hists.hists["CCQElike"].Clone(),mc_hists.hists["notCCQElike"].Clone(),pot_scale)
    piZero = AddTwoHistograms(mc_hists.hists["NCPi0"].Clone(),mc_hists.hists["CCPi0"].Clone(),pot_scale)
    if sideband !="Pi0":
        piZero = Scalepizero(piZero)
    total_mc = mc_hists.hists["Total"].Clone()
    for i,band in enumerate(bands): # Here, each of the subbands that are being calculated are removed from the total
        total_mc = SubtractPoissonHistograms(total_mc,mc_hists.hists[band],pot_scale,not bool(i)) 
    augment = SubtractPoissonHistograms(data,total_mc,pot_scale) # This subtracts the LHS const from the data to give the RHS of the equation
    return [coherent,nu_e,piZero,ncdiff,augment]

def ErrExtractHists(data_hists, mc_hists, sideband, err_type, universe_num):
    pot_scale = mc_hists.pot_scale
    data = data_hists.GetHist().Clone()
    coherent = AddTwoHistograms(mc_hists.hists["NCCohPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["lowNCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),pot_scale)
    ncdiff = mc_hists.hists["NCDiff"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    nu_e = AddTwoHistograms(mc_hists.hists["CCQElike"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["notCCQElike"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),pot_scale)
    piZero= AddTwoHistograms(mc_hists.hists["NCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["CCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),pot_scale)
    if sideband !="Pi0":
        piZero = Scalepizero(piZero)
    total_mc = mc_hists.hists["Total"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    bands = ["NCCohPi0","NCDiff","CCQElike","notCCQElike","NCPi0","CCPi0","lowNCPi0"]
    for i,band in enumerate(bands): # Here, each of the subbands that are being calculated are removed from the total
        total_mc = SubtractPoissonHistograms(total_mc,mc_hists.hists[band].GetVertErrorBand(err_type).GetHist(universe_num),pot_scale,not bool(i))
    other = PlotUtils.MnvH1D(total_mc) # I have no idea why the add command can't do this itself but whatever
    augment = SubtractPoissonHistograms(data,other,pot_scale) # This subtracts the LHS const from the data to give the RHS of the equation
    return [coherent,nu_e,piZero,ncdiff,augment]


def HistListtoArray(histlist): # This converts the root histograms being passed to it into an array
    histarray = np.zeros((len(histlist),histlist[0].GetSize()-1))
    for i,hist in enumerate(histlist):
        for j in range(hist.GetSize()-1):
            histarray[i,j] = hist.GetBinContent(j+1)
    if len(histlist) == 1:
        histarray = histarray.flatten()
    return histarray

def SetupMatrices(dataarray): # This converts the weird data array that has a weird data structure into matrices that get the weights for each bin
    numofmatrices = len(dataarray[0,0])
    matrices = np.zeros((numofmatrices,4,4))
    vectors = np.zeros((numofmatrices,4))
    for i in range(numofmatrices):
        for j in range(4):
            matrices[i,j] = dataarray[j,0:4,i]
            vectors[i,j] = dataarray[j,4,i]
    for i in range(numofmatrices):
        for j in range(len(matrices[0])):
            vectors[i,j] -= np.sum(matrices[i,j,0:4]) # This makes it so that the matrices are finding (weight-1) 
    return matrices,vectors

def Weightuncertainity(dataarray): # This converts the weird data array that has a weird data structure into matrices that get the weights for each bin
    numofmatrices = len(dataarray[0,0])
    matrices = np.zeros((numofmatrices,4,4))
    vectors = np.zeros((numofmatrices,4))
    for i in range(numofmatrices):
        for j in range(4):
            vectors[i,j] = dataarray[j,4,i]
    return vectors


def DoLinearAlgebra(inversesArray,vectorArray):
    solutions = np.zeros(vectorArray.shape)
    for i in range(len(inversesArray)):
        #invMatrix = np.linalg.inv(matrixArray[i])
        solutions[i] = np.dot(inversesArray[i],vectorArray[i])
    final_weights = np.add(np.ones(solutions.shape),solutions) # Corrects for the fact that the matrices are solving for weight-1
    return final_weights

def NormalInverse(matrices): #This shouldn't be called anywhere it's useful for testing how much the SVD made things weird
    inverses = []
    for matrix in matrices:
        inverse = np.linalg.inv(matrix)
        inverses.append(inverse)
    return inverses

def MPInverse(matrices):
    svds = []
    for matrix in matrices: # Calculates all of the SVDs and puts them in a list
        decomposition = np.linalg.svd(matrix)
        svds.append(decomposition)
        #print("decomposition[1]=",decomposition[1])

    discarded = np.array([[False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False, False],
                          [False, False, False, False],
                          [False,False, False,False]],dtype = bool) #Picks singular values to be discarded
    # These may look different for you and you may need to rethink these
    for i,svd in enumerate(svds):
        for j in range(len(svd[1])):
            if discarded[i,j]:
                svd[1][j] = 0
    inverses = []
    for svd in svds: # Constructs the MP Pseudoinverse using the SVD
        Ustar = np.transpose(svd[0])
        V = np.transpose(svd[2])
        sigma = svd[1]
        sigmaPlus = np.zeros((len(sigma),len(sigma)))
        for i in range(len(sigma)):
            if sigma[i] == 0:
                sigmaPlus[i,i] = 0
            else: sigmaPlus[i,i] = 1/sigma[i]
        mpinverse = np.matmul(np.matmul(V,sigmaPlus),Ustar)
        inverses.append(mpinverse)
    return inverses

def GetBins(hist,print_results = False): # If you want to print out the binds in a hist, use this
    bins = []
    for i in range(hist.GetSize()):
        bins.append(hist.GetBinContent(i))
    if print_results: print(hist,bins)
    return bins

def GetPreweightedROOTFiles(weighting_playlist):
    mc_path = "/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/tools/kin_dist_mc"+str(weighting_playlist)+"_collab1_fspline.root"
    mc_file = ROOT.TFile.Open(mc_path)
    data_path = "/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/tools/kin_dist_data"+str(weighting_playlist)+"_collab1_fspline.root"
    data_file = ROOT.TFile.Open(data_path)
    pot_scale = GetPOTScale(data_path,mc_path)
    return mc_file,data_file,pot_scale

def GetCVWeights(mc_file,data_file,pot_scale):
    sidebands = ["Signal","Diffractive","Electron_Neutrino","Pi0"]
    # First going through the weights for the prime universe because I want to have this be its own thing
    sideband_list = []
    for sideband in sidebands:
        mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale) # Gets the two electron energy histograms
        dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False) # I'm using HistHolders to access the backgrounds
        #mcEel_hists.hists["NCCohPi0"].Print("All")
        mcEel_hists.POTScale(False)
        dataEel_hist.POTScale(False)
        #mcEel_hists.hists["NCCohPi0"].Print("All" )
        listofhists = PrimeExtractHists(dataEel_hist, mcEel_hists,sideband)
        arrayofhists = HistListtoArray(listofhists)
        sideband_list.append(arrayofhists)
    bigArray = np.array(sideband_list)
    matrices,vectors = SetupMatrices(bigArray)
    invmatrices = MPInverse(matrices)
    primeweights = DoLinearAlgebra(invmatrices,vectors)
    WasIRight(matrices,primeweights,vectors)
    
    weight_uncertainity = Weightuncertainity(bigArray)
    cov_matrix_data = [np.diag(np.sqrt(np.abs(weight_uncertainity[i]))) for i in range(len(vectors))]
    cov_matrix_data = np.array(cov_matrix_data)
    weight_covariances = [invmatrices[i] @ cov_matrix_data[i] @ invmatrices[i].T for i in range(len(vectors))]
    reco_weight_errors = [np.sqrt(np.diag(cov)) for cov in weight_covariances]
    return primeweights,reco_weight_errors

def WasIRight(matrices,weights,vectors):
    for i in range(len(matrices)):
        print("MC: ",np.dot(matrices[i],weights[i]-np.ones(4)))
        print("Data: ",vectors[i])


def RunAlexWeighter():
    #print("Starting Alex's Reweighting Weighter")
    mc_file,data_file,pot_scale = GetPreweightedROOTFiles(None)
    cvweights = GetCVWeights(mc_file,data_file,pot_scale)
    universe_weights_dict = GetSystematicUniverseWeights(mc_file,data_file,pot_scale)
    #print("Weighter Ran Successfully")
    return cvweights,universe_weights_dict

#def GetSystematicUniverseWeights(mc_file,data_file,pot_scale):
#    # Now for the Systematic Universes
#    sidebands = ["Signal","Diffractive","Electron_Neutrino","Pi0"]
#    #sidebands = ["Diffractive","Electron_Neutrino","Pi0","Signal"]
#    err_weights_dict = {}
#    universe_list = [["Flux",200],["GENIE_AGKYxF1pi",2],["GENIE_AhtBY",2],["GENIE_BhtBY",2],["GENIE_CCQEPauliSupViaKF",2],["GENIE_CV1uBY",2],["GENIE_CV2uBY",2],["GENIE_EtaNCEL",2],["GENIE_FrAbs_N",2],["GENIE_FrAbs_pi",2],["GENIE_FrCEx_N",2],["GENIE_FrCEx_pi",2],["GENIE_FrElas_N",2],["GENIE_FrElas_pi",2],["GENIE_FrInel_N",2],["GENIE_FrPiProd_N",2],["GENIE_FrPiProd_pi",2],["GENIE_MFP_N",2],["GENIE_MFP_pi",2],["GENIE_MaNCEL",2],["GENIE_MaRES",2],["GENIE_MvRES",2],["GENIE_NormCCRES",2],["GENIE_NormDISCC",2],["GENIE_NormNCRES",2],["GENIE_RDecBR1gamma",2],["GENIE_Rvn2pi",2],["GENIE_Rvp2pi",2],["GENIE_Theta_Delta2Npi",2],["GENIE_VecFFCCQEshape",2],["GENIE_Rvn1pi",2],["GENIE_Rvp1pi",2],["GENIE_MaCCQE",2],["Low_Recoil_2p2h_Tune",3],["RPA_HighQ2",2],["RPA_LowQ2",2],["LowQ2Pi",2],["LowQ2Pi_None",1],["MK_model",1],["fsi_weight",3],["SuSA_Valencia_Weight",1],["GEANT_Proton",2],["GEANT_Neutron",2],["GEANT_Pion",2],["Target_Mass_CH",2],["Muon_Energy_MINERvA",2],["Muon_Energy_MINOS",2],["beam_angle",4],["response_p",2],["response_meson",2],["response_em",2],["response_other",2],["response_xtalk",2],["Leakage_Uncertainty",2]]
#
#    for error_type,error_universes in universe_list:
#        err_type_weights = []
#        for universe_num in range(error_universes):
#            sideband_list = []
#            for sideband in sidebands:
#                mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale)
#                dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False)
#                #mcEel_hists.POTScale(False)   If i turn these on my vertical uncertainity goes very very very high??
#                #dataEel_hist.POTScale(False)
#                listofhists = ErrExtractHists(dataEel_hist, mcEel_hists, error_type, universe_num)
#                arrayofhists = HistListtoArray(listofhists)
#                sideband_list.append(arrayofhists)
#            bigArray = np.array(sideband_list)
#            matrices,vectors = SetupMatrices(bigArray)
#            inverse_matrices = np.array(MPInverse(matrices))
#            err_type_weights.append(DoLinearAlgebra(MPInverse(matrices),vectors))
#            del mcEel_hists, dataEel_hist
#            #*****************************
#        #********************************
#        err_weights_dict[error_type] = np.array(err_type_weights)
#    return err_weights_dict



def GetSystematicUniverseWeights(mc_file,data_file,pot_scale):
    # Now for the Systematic Universes
    sidebands = ["Signal","Diffractive","Electron_Neutrino","Pi0"]
    err_weights_dict = {}
    err_weights_stat_err_dict = {}
    universe_list = [["Flux",200],["GENIE_AGKYxF1pi",2],["GENIE_AhtBY",2],["GENIE_BhtBY",2],["GENIE_CCQEPauliSupViaKF",2],["GENIE_CV1uBY",2],["GENIE_CV2uBY",2],["GENIE_EtaNCEL",2],["GENIE_FrAbs_N",2],["GENIE_FrAbs_pi",2],["GENIE_FrCEx_N",2],["GENIE_FrCEx_pi",2],["GENIE_FrElas_N",2],["GENIE_FrElas_pi",2],["GENIE_FrInel_N",2],["GENIE_FrPiProd_N",2],["GENIE_FrPiProd_pi",2],["GENIE_MFP_N",2],["GENIE_MFP_pi",2],["GENIE_MaNCEL",2],["GENIE_MaRES",2],["GENIE_MvRES",2],["GENIE_NormCCRES",2],["GENIE_NormDISCC",2],["GENIE_NormNCRES",2],["GENIE_RDecBR1gamma",2],["GENIE_Rvn2pi",2],["GENIE_Rvp2pi",2],["GENIE_Theta_Delta2Npi",2],["GENIE_VecFFCCQEshape",2],["GENIE_Rvn1pi",2],["GENIE_Rvp1pi",2],["GENIE_MaCCQE",2],["Low_Recoil_2p2h_Tune",3],["RPA_HighQ2",2],["RPA_LowQ2",2],["LowQ2Pi",2],["LowQ2Pi_None",1],["MK_model",1],["fsi_weight",3],["SuSA_Valencia_Weight",1],["GEANT_Proton",2],["GEANT_Neutron",2],["GEANT_Pion",2],["Target_Mass_CH",2],["Muon_Energy_MINERvA",2],["Muon_Energy_MINOS",2],["beam_angle",4],["response_p",2],["response_meson",2],["response_em",2],["response_other",2],["response_xtalk",2],["Leakage_Uncertainty",2]]

    for error_type,error_universes in universe_list:
        err_type_weights = []
        err_type_weight_errors = []
        for universe_num in range(error_universes):
            sideband_list = []
            for sideband in sidebands:
                mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale)
                dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False)
#                #mcEel_hists.POTScale(False)   If i turn these on my vertical uncertainity goes very very very high??
#                #dataEel_hist.POTScale(False)
                listofhists = ErrExtractHists(dataEel_hist, mcEel_hists, sideband, error_type, universe_num)
                arrayofhists = HistListtoArray(listofhists)
                sideband_list.append(arrayofhists)
            bigArray = np.array(sideband_list)
            matrices,vectors = SetupMatrices(bigArray)
            invmatrices = MPInverse(matrices)
            weights = DoLinearAlgebra(invmatrices, vectors)
            err_type_weights.append(weights)
            
            #Calculate stat uncertainity
            weight_uncertainty = Weightuncertainity(bigArray)   #N_data
            cov_matrix_data = [np.diag(np.sqrt(np.abs(weight_uncertainty[i]))) for i in range(len(vectors))]
            weight_covariance = [invmatrices[i] @ cov_matrix_data[i] @ invmatrices[i].T for i in range(len(vectors))]
            reco_y_error = [np.sqrt(np.diag(cov)) for cov in weight_covariance]
            err_type_weight_errors.append(reco_y_error)

            del mcEel_hists, dataEel_hist
        err_weights_dict[error_type] = np.array(err_type_weights)
        err_weights_stat_err_dict[error_type] = np.array(err_type_weight_errors)
    return err_weights_dict,err_weights_stat_err_dict



if __name__ == "__main__":
    preweighted_playlist = AnalysisConfig.playlist
    mc_file,data_file,pot_scale = GetPreweightedROOTFiles(preweighted_playlist)
    cvweights,reco_y_error = GetCVWeights(mc_file,data_file,pot_scale)
    print("cvweights=",cvweights)
    universe_weights_dict,sys_stat_error = GetSystematicUniverseWeights(mc_file,data_file,pot_scale)
    np.set_printoptions(precision=6,suppress = True,threshold = np.inf)
    np.savetxt("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/tools/weights.txt",cvweights,delimiter=",")

    dict_str = str(universe_weights_dict)
    with open("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/tools/systematic_universes_weights.txt", "w") as file:
        file.write(dict_str)
    Eel_hist = mc_file.Eel.Clone()
    weightFile = ROOT.TFile.Open("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/tools/kin_dist_"+str(preweighted_playlist)+"_weights_with_staterror.root","RECREATE")
    coherent = PlotUtils.MnvH1D("coherent","Weights on Coherent Events",len(ELECTRON_ENERGY_BINNING)-1,array("d",ELECTRON_ENERGY_BINNING))
    electron_neutrino = PlotUtils.MnvH1D("electron_neutrino","Weights on Electron Neutrino Events",len(ELECTRON_ENERGY_BINNING)-1,array("d",ELECTRON_ENERGY_BINNING))
    pizero = PlotUtils.MnvH1D("pizero","Weights on PiZero Events",len(ELECTRON_ENERGY_BINNING)-1,array("d",ELECTRON_ENERGY_BINNING))
    ncdiff = PlotUtils.MnvH1D("ncdiff","Weights on NCDiff Events",len(ELECTRON_ENERGY_BINNING)-1,array("d",ELECTRON_ENERGY_BINNING))
    for i in range(len(ELECTRON_ENERGY_BINNING)-1):
        coherent.SetBinContent(i+1,cvweights[i,0])
        electron_neutrino.SetBinContent(i+1,cvweights[i,1])
        pizero.SetBinContent(i+1,cvweights[i,2])
        ncdiff.SetBinContent(i+1,cvweights[i,3])
        coherent.SetBinError(i+1,reco_y_error[i][0])
        electron_neutrino.SetBinError(i+1,reco_y_error[i][1])
        pizero.SetBinError(i+1,reco_y_error[i][2])
        ncdiff.SetBinError(i+1,reco_y_error[i][3])
    coherent.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    electron_neutrino.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    pizero.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    ncdiff.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    for error_type in coherent.GetErrorBandNames():
        error_type = str(error_type) # I really hate that this line has to be here
        for k in range(coherent.GetVertErrorBand(error_type).GetNHists()):
            for i in range(len(ELECTRON_ENERGY_BINNING)-1):
                coherent.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,universe_weights_dict[error_type][k][i,0])
                coherent.GetVertErrorBand(error_type).GetHist(k).SetBinError(i+1,sys_stat_error[error_type][k][i,0])
                electron_neutrino.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,universe_weights_dict[error_type][k][i,1])
                electron_neutrino.GetVertErrorBand(error_type).GetHist(k).SetBinError(i+1,sys_stat_error[error_type][k][i,1])
                pizero.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,universe_weights_dict[error_type][k][i,2])
                pizero.GetVertErrorBand(error_type).GetHist(k).SetBinError(i+1,sys_stat_error[error_type][k][i,2])
                ncdiff.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,universe_weights_dict[error_type][k][i,3])
                ncdiff.GetVertErrorBand(error_type).GetHist(k).SetBinError(i+1,sys_stat_error[error_type][k][i,3])
    coherent.GetXaxis().SetTitle("Lepton Energy")
    coherent.SetMarkerStyle(21)
    electron_neutrino.GetXaxis().SetTitle("Lepton Energy")
    electron_neutrino.SetMarkerStyle(21)
    pizero.GetXaxis().SetTitle("Lepton Energy")
    pizero.SetMarkerStyle(21)
    ncdiff.GetXaxis().SetTitle("Lepton Energy")
    ncdiff.SetMarkerStyle(21)
    coherent.Write()
    pizero.Write()
    ncdiff.Write()
    electron_neutrino.Write()
    #weightFile.Close()


    coherent_canvas = ROOT.TCanvas("c2","c2",1600,1200)
    coherent_err = coherent.GetCVHistoWithError()
    coherent_err.SetMarkerStyle(21)
    coherent_err.SetMarkerColor(4)
    coherent_err.SetLineColor(9)
    coherent_err.SetMarkerSize(2)
    coherent_err.SetLineWidth(3)
    coherent_err.GetXaxis().SetTitle("Lepton Energy (GeV)")
    coherent_err.SetStats(0)
    coherent_err.Draw("MIN0 E1")
    #coherent_err.Print("All")
    #ROOT.gPad.BuildLegend(0.85,0.8,0.85,0.6)
    coherent_canvas.Print("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/results_storage/coherent_err_2"+str(preweighted_playlist)+".png","png")
    
    Coh_wt_error_canvas = ROOT.TCanvas("c1","c1",1600,1200)
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
    PlotTools.MakeErrPlot(coherent)
    #PlotTools.MNVPLOTTER.axis_maximum = 0.2
    #mnvplotter.DrawErrorSummary(new_dataEel_hist)
    Coh_wt_error_canvas.Print("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/results_storage/Coh_wt_error_2"+str(preweighted_playlist)+".png","png")


    pizero_canvas = ROOT.TCanvas("c2","c2",1600,1200)
    pizero_err = pizero.GetCVHistoWithError()
    pizero_err.SetMarkerStyle(21)
    pizero_err.SetMarkerColor(4)
    pizero_err.SetLineColor(9)
    pizero_err.SetMarkerSize(2)
    pizero_err.SetLineWidth(3)
    pizero_err.GetXaxis().SetTitle("Lepton Energy (GeV)")
    pizero_err.SetStats(0)
    pizero_err.Draw("MIN0 E1")
    #pizero_err.Print("All")
    #ROOT.gPad.BuildLegend(0.85,0.8,0.85,0.6)
    pizero_canvas.Print("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/results_storage/pizero_err_2"+str(preweighted_playlist)+".png","png")
    
    pizero_wt_error_canvas = ROOT.TCanvas("c1","c1",1600,1200)
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
    PlotTools.MakeErrPlot(pizero)
    #PlotTools.MNVPLOTTER.axis_maximum = 0.2
    #mnvplotter.DrawErrorSummary(new_dataEel_hist)
    pizero_wt_error_canvas.Print("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/results_storage/PiZero_wt_error_2"+str(preweighted_playlist)+".png","png")



    ncdiff_canvas = ROOT.TCanvas("c2","c2",1600,1200)
    ncdiff_err = ncdiff.GetCVHistoWithError()
    ncdiff_err.SetMarkerStyle(21)
    ncdiff_err.SetMarkerColor(4)
    ncdiff_err.SetLineColor(9)
    ncdiff_err.SetMarkerSize(2)
    ncdiff_err.SetLineWidth(3)
    ncdiff_err.GetXaxis().SetTitle("Lepton Energy (GeV)")
    ncdiff_err.SetStats(0)
    ncdiff_err.Draw("MIN0 E1")
    #ncdiff_err.Print("All")
    #ROOT.gPad.BuildLegend(0.85,0.8,0.85,0.6)
    ncdiff_canvas.Print("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/results_storage/ncdiff_err_2"+str(preweighted_playlist)+".png","png")

    ncdiff_wt_error_canvas = ROOT.TCanvas("c1","c1",1600,1200)
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
    PlotTools.MakeErrPlot(ncdiff)
    #PlotTools.MNVPLOTTER.axis_maximum = 0.2
    #mnvplotter.DrawErrorSummary(new_dataEel_hist)
    ncdiff_wt_error_canvas.Print("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/results_storage/Ncdiff_wt_error_2"+str(preweighted_playlist)+".png","png")


    electron_neutrino_canvas = ROOT.TCanvas("c2","c2",1600,1200)
    electron_neutrino_err = electron_neutrino.GetCVHistoWithError()
    electron_neutrino_err.SetMarkerStyle(21)
    electron_neutrino_err.SetMarkerColor(4)
    electron_neutrino_err.SetLineColor(9)
    electron_neutrino_err.SetMarkerSize(2)
    electron_neutrino_err.SetLineWidth(3)
    electron_neutrino_err.GetXaxis().SetTitle("Lepton Energy (GeV)")
    electron_neutrino_err.SetStats(0)
    electron_neutrino_err.Draw("MIN0 E1")
    #electron_neutrino_err.Print("All")
    #ROOT.gPad.BuildLegend(0.85,0.8,0.85,0.6)
    electron_neutrino_canvas.Print("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/results_storage/electron_neutrino_err_2"+str(preweighted_playlist)+".png","png")
    
    electron_neutrino_wt_error_canvas = ROOT.TCanvas("c1","c1",1600,1200)
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
    PlotTools.MakeErrPlot(electron_neutrino)
    #PlotTools.MNVPLOTTER.axis_maximum = 0.2
    #mnvplotter.DrawErrorSummary(new_dataEel_hist)
    electron_neutrino_wt_error_canvas.Print("/exp/minerva/app/users/nalex/nccoh/NCCoh-X-Sec/results_storage/electron_neutrino_wt_error_2"+str(preweighted_playlist)+".png","png")


#weightFile.Close()
#mc_file.Close()
#data_file.Close()
