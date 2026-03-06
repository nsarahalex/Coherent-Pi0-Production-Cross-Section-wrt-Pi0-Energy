# Nimmy Sar\h
# Matrix Unfolding using SVD for Coherent Pion Production

import ROOT
import numpy as np
import math
import PlotUtils
from collections import OrderedDict
from array import array

from config.AnalysisConfig import AnalysisConfig
from tools import Utilities,PlotTools
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS

ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# Classes to convert ROOT histograms into numpy arrays with values, error, and binning
class Variable1D:
    def __init__(self,name,filePath,skipFirst = False,isEfficiency = False,univ_type = None,univ_num = None):
        self.name = name
        self.isCV = not bool(univ_type)
        self.isEfficiency = isEfficiency
        if self.isCV:
            self.values,self.error,self.bins = self.GetCVValsErrorandBins(name,filePath,skipFirst)
            self.universes = self.UniverseListAssembler(name,filePath,skipFirst)
        else:
            self.univ_type = univ_type
            self.univ_num = univ_num
            self.values,self.error,self.bins = self.GetUnivValsErrorandBins(name,filePath,skipFirst,univ_type,univ_num)
        #if self.isEfficiency: self.values = np.ones(len(self.values)) # For testing with no efficiency

    def GetCVValsErrorandBins(self,name,filePath,skipFirst):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        self.hist = dataHist.Clone()
        dataArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
        errorArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
        dataBins = np.zeros(len(dataArray)+1)
        for i in range(len(dataArray)):
            dataArray[i] = dataHist.GetBinContent(i+1+skipFirst)
            errorArray[i] = dataHist.GetBinError(i+1+skipFirst)
            dataBins[i] = dataHist.GetBinLowEdge(i+1+skipFirst)
        dataBins[i+1] = dataHist.GetBinLowEdge(i+2+skipFirst) # Gets upper edge of last bin
        return dataArray,errorArray,dataBins

    def UniverseListAssembler(self,name,filePath,skip1st):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        universedict = {}
        for error_type in dataHist.GetErrorBandNames():
            error_type_univ_list = []
            for k in range(dataHist.GetVertErrorBand(error_type).GetNHists()):
                error_type_univ_list.append(Variable1D(name,filePath,skipFirst = skip1st, isEfficiency = self.isEfficiency,univ_type = error_type,univ_num = k))
            universedict[str(error_type)] = error_type_univ_list
        return universedict

    def GetUnivValsErrorandBins(self,name,filePath,skipFirst,error_type,k):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        self.hist = dataHist.Clone()
        dataArray = np.zeros(dataHist.GetVertErrorBand(error_type).GetHist(k).GetNbinsX()-skipFirst)
        errorArray = np.zeros(dataHist.GetVertErrorBand(error_type).GetHist(k).GetNbinsX()-skipFirst)
        dataBins = np.zeros(len(dataArray)+1)
        for i in range(len(dataArray)):
            dataArray[i] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinContent(i+1+skipFirst)
            errorArray[i] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinError(i+1+skipFirst)
            dataBins[i] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinLowEdge(i+1+skipFirst)
        dataBins[i+1] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinLowEdge(i+2+skipFirst) # Gets upper edge of last bin
        return dataArray,errorArray,dataBins



class MigrationMatrix:
    def __init__(self,name,filePath,skipFirstX = False, univ_type = None,univ_num = None):
        self.name = name
        self.isCV = not bool(univ_type)
        if self.isCV:
            self.values,self.error,self.xbins,self.ybins = self.GetCVValsErrorandBins(name,filePath,skipFirstX)
            self.universes = self.UniverseListAssembler(name,filePath,skipFirstX)
        else:
            self.univ_type = univ_type
            self.univ_num = univ_num
            self.values,self.error,self.xbins,self.ybins = self.GetUnivValsErrorandBins(name,filePath,skipFirstX,univ_type,univ_num)

    def GetCVValsErrorandBins(self,name,filePath,skipFirstX):
        dataFile = ROOT.TFile.Open(filePath,"READ")

        dataHist = eval("dataFile."+str(name)+"_NCCohPi0.Clone()")
        self.hist = dataHist.Clone()
        dataArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        errorArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        for i in range(len(dataArray)):
            for j in range(len(dataArray[0])):
                dataArray[i,j] = dataHist.GetBinContent(j+1+skipFirstX,i+1)
                errorArray[i,j] = dataHist.GetBinError(j+1+skipFirstX,i+1)
        for i in range(len(dataArray)):
            row_sum = np.sum(dataArray[i])
            dataArray[i] /= row_sum
            errorArray[i] /=row_sum
            #errorArray[i] = np.sqrt(dataArray[i] * (1 - dataArray[i]) / row_sum)
            #errorArray[i] /= np.sum(dataArray[i])  #nimmy changed this from /=np.sum(dataArray[i])
            #dataArray[i] /= np.sum(dataArray[i])
        dataArray = np.transpose(dataArray)
        errorArray = np.transpose(errorArray)
        ybins = np.zeros(len(dataArray)+1)
        xbins = np.zeros(len(dataArray[0])+1)
        yaxis = dataHist.GetYaxis()
        xaxis = dataHist.GetXaxis()
        for i in range(len(ybins)):
            ybins[i] = yaxis.GetBinLowEdge(i+1)
        for j in range(len(xbins)):
            xbins[j] = xaxis.GetBinLowEdge(j+1+skipFirstX)
        #for i in range(len(xbins)-1): # this chunk is for xbin normalized migration
        #    dataArray[:,i] /= (xbins[i+1]-xbins[i])
        #for i in range(len(dataArray)): dataArray[i] /= np.sum(dataArray[i])
        return dataArray,errorArray,xbins,ybins

    def UniverseListAssembler(self,name,filePath,skip1stX):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+"_NCCohPi0.Clone()")
        dataHist.AddMissingErrorBandsAndFillWithCV(Eel.hist)
        universedict = {}
        for error_type in dataHist.GetErrorBandNames():
            error_type_univ_list = []
            for k in range(dataHist.GetVertErrorBand(error_type).GetNHists()):
                error_type_univ_list.append(MigrationMatrix(name,filePath,skipFirstX = skip1stX,univ_type = error_type,univ_num = k))
            universedict[str(error_type)] = error_type_univ_list
        return universedict

    def GetUnivValsErrorandBins(self,name,filePath,skipFirstX,error_type,k):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+"_NCCohPi0.Clone()")
        dataHist.AddMissingErrorBandsAndFillWithCV(Eel.hist)
        dataArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        errorArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        for i in range(len(dataArray)):
            for j in range(len(dataArray[0])):
                dataArray[i,j] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinContent(j+1+skipFirstX,i+1)
                errorArray[i,j] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinError(j+1+skipFirstX,i+1)
        for i in range(len(dataArray)):
            row_sum = np.sum(dataArray[i])
            dataArray[i] /= row_sum
            errorArray[i] /=row_sum
            #errorArray[i] = np.sqrt(dataArray[i] * (1 - dataArray[i]) / row_sum)

             #dataArray[i] /= np.sum(dataArray[i])
        dataArray = np.transpose(dataArray)
        errorArray = np.transpose(errorArray)
        ybins = np.zeros(len(dataArray)+1)
        xbins = np.zeros(len(dataArray[0])+1)
        yaxis = dataHist.GetYaxis()
        xaxis = dataHist.GetXaxis()
        for i in range(len(ybins)):
            ybins[i] = yaxis.GetBinLowEdge(i+1)
        for j in range(len(xbins)):
            xbins[j] = xaxis.GetBinLowEdge(j+1+skipFirstX)
        #for i in range(len(xbins)-1): # this chunk is for xbin normalized migration
        #    dataArray[:,i] /= (xbins[i+1]-xbins[i])
        #for i in range(len(dataArray)): dataArray[i] /= np.sum(dataArray[i])

        return dataArray,errorArray,xbins,ybins





def CovtoCorrMatrix(cov_matrix,yname,printMatrix = True,isCV=True):
    if not isCV: return None
    diag = np.sqrt(np.diag(cov_matrix))
    denom = np.outer(diag, diag)
    corr_matrix = cov_matrix / denom
    corr_matrix[denom == 0] = 0
    return corr_matrix

def MPInverse(matrix):
    threshold = 0.4
    U,sigma,Vt = np.linalg.svd(matrix,full_matrices=False)
    print(sigma) 
    discarded = np.array([False,False,False,False,False,False,False,False,True,True],dtype = bool)
    #discarded = np.array([False,False,False,False,False,False,False,True,True,True],dtype = bool)
    
    if len(discarded) != len(sigma):
        raise ValueError("Discarded mask length must match number of singular values")
    

    sigma_inv = np.zeros_like(sigma)
    sigma_inv[~discarded] = 1.0 / sigma[~discarded]
   
    mpinverse = (Vt.T * sigma_inv) @ U.T
    return mpinverse,threshold
def get_cov_matrix(hist):
    """Return a clean 2D NumPy array of the covariance matrix from a ROOT MnvH1D/MnvH2D."""
    cov_root = hist.GetStatErrorMatrix()
    nrows, ncols = cov_root.GetNrows(), cov_root.GetNcols()
    cov_matrix = np.zeros((nrows, ncols), dtype=float)
    for i in range(nrows):
        for j in range(ncols):
            cov_matrix[i, j] = cov_root[i, j]
    # Optionally trim underflow/overflow bins
    if np.sum(cov_matrix[0]) == 0 and nrows > 2:
        cov_matrix = cov_matrix[1:-1, 1:-1]
    return cov_matrix



def InverseMigration(data, migration, efficiency):
    invmatrix = MPInverse(migration.values)[0]
    cov_x = get_cov_matrix(data.hist) #np.diag(data.values)
    # 3. Propagate covariance
    #print("invmatrix.shape=",invmatrix.shape)
    #print("cov_x.shape=", cov_x.shape)
    cov_y = invmatrix @ cov_x @ invmatrix.T
    # 4. True bin reconstruction
    reco_y_dist = invmatrix @ data.values
    reco_y_error = np.sqrt(np.diag(cov_y))
    return reco_y_dist, reco_y_error, cov_y

def MatrixTruthRecovery(data,migration,efficiency,xname,yname,PiZeroE_NCCOhPi0_dist): # These have to be the CV classes that include the universes
    mmatrix = migration.values 
    invmatrix,threshold = MPInverse(mmatrix)
    reco_y_dist = invmatrix @ data.values #truth
    x_predicted_dist = mmatrix @ reco_y_dist #refolding truth to get reco
    nbins = PiZeroE_NCCOhPi0_dist.GetNbinsX() 
    refold_pizeroE_dist = mmatrix @ np.array([PiZeroE_NCCOhPi0_dist.GetBinContent(bin_idx)for bin_idx in range(1, nbins + 1)])
    
    cov_root_x = data.hist.GetStatErrorMatrix()
    cov_matrix_x = np.zeros((cov_root_x.GetNrows(),cov_root_x.GetNcols()))
    for i in range(cov_root_x.GetNrows()):
        for j in range(cov_root_x.GetNcols()):
            cov_matrix_x[i,j] = cov_root_x[i,j]

    if np.all(cov_matrix_x[:, 1] == 0):
        cov_matrix_x = cov_matrix_x[2:-1, 2:-1]
    else:
        cov_matrix_x = cov_matrix_x[1:-1, 1:-1]
    cov_matrix_x = np.array(cov_matrix_x, dtype=float)
    cov_matrix_y = invmatrix @ cov_matrix_x @ invmatrix.T
    #cov_matrix_y = invmatrix @ cov_matrix_x @ np.transpose(invmatrix) #cov_matrix_y = variance_in_true_coh
    Corr_matrix = CovtoCorrMatrix(cov_matrix_y,yname)
    new_cov_matrix_x = mmatrix @ cov_matrix_y @ mmatrix.T #correlation matrix
    reco_y_error = np.sqrt(np.diag(cov_matrix_y)) #reco_y_error = cov0
    x_predicted_error = np.sqrt(np.diag(new_cov_matrix_x))
    reco_y_error_hist =  PlotUtils.MnvH1D(str(yname)+"_reco_y_error_hist","True "+str(yname)+" From Inverse Migration Matrix",len(mmatrix),array("d",migration.ybins))
    correlation_hist = PlotUtils.MnvH2D(str(yname)+"_correlation_hist","True "+str(yname)+" From Inverse Migration Matrix",len(mmatrix),array("d",migration.xbins),len(mmatrix),array("d",migration.ybins)) 
    y_recovered_hist = PlotUtils.MnvH1D(str(yname)+"_recovered","True "+str(yname)+" From Inverse Migration Matrix",len(mmatrix),array("d",migration.ybins))
    x_predicted_hist = PlotUtils.MnvH1D(str(xname)+"_matrixpredicted","Predicted "+str(xname)+" From MPInverse Matrix",len(mmatrix[0]),array("d",migration.xbins))
    refold_pizeroE_hist = PlotUtils.MnvH1D(str(xname)+"_piZeroE_NCCoherent_refolded_Eel","piZeroE_refolded_Eel "+str(xname)+" From MPInverse Matrix",len(mmatrix[0]),array("d",migration.xbins))
    for i in range(len(reco_y_dist)):
        y_recovered_hist.SetBinContent(i+1,reco_y_dist[i])
        y_recovered_hist.SetBinError(i+1,reco_y_error[i])
        reco_y_error_hist.SetBinContent(i+1,reco_y_error[i])
    for j in range(len(x_predicted_dist)):
        x_predicted_hist.SetBinContent(j+1,x_predicted_dist[j])
        x_predicted_hist.SetBinError(j+1,x_predicted_error[j])
        refold_pizeroE_hist.SetBinContent(j+1,refold_pizeroE_dist[j])
    x_predicted_hist.AddMissingErrorBandsAndFillWithCV(data.hist)
    refold_pizeroE_hist.AddMissingErrorBandsAndFillWithCV(data.hist)
    y_recovered_hist.AddMissingErrorBandsAndFillWithCV(data.hist)
    chi2 = CalChi2(x_predicted_dist,data.error,data.values)
    print("done untl here")
    for ix in range(np.shape(Corr_matrix)[0]):
        for iy in range(np.shape(Corr_matrix)[1]):
            correlation_hist.SetBinContent(ix+1,iy+1,Corr_matrix[ix][iy])
    for error_type in data.hist.GetErrorBandNames():
        #error_type = str(error_type) # I really hate that this line has to be here
        for k in range(data.hist.GetVertErrorBand(error_type).GetNHists()):
            data_univ = data.universes[error_type][k]
            mig_univ = migration.universes[error_type][k]
            eff_univ = efficiency.universes[error_type][k]
            y_recovered_dist, y_recovered_error, cov_y_univ = InverseMigration(data_univ, mig_univ, eff_univ)
            x_predicted_dist = mig_univ.values @ y_recovered_dist
            x_predicted_error = np.sqrt(np.diag(mig_univ.values @ cov_y_univ @ mig_univ.values.T))
            for i in range(len(y_recovered_dist)): 
                y_recovered_hist.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,y_recovered_dist[i])
                y_recovered_hist.GetVertErrorBand(error_type).GetHist(k).SetBinError(i+1,y_recovered_error[i])
            for j in range(len(x_predicted_dist)):
                x_predicted_hist.GetVertErrorBand(error_type).GetHist(k).SetBinContent(j+1,x_predicted_dist[j])
                x_predicted_hist.GetVertErrorBand(error_type).GetHist(k).SetBinError(j+1,x_predicted_error[j])
                refold_pizeroE_hist.GetVertErrorBand(error_type).GetHist(k).SetBinContent(j+1,refold_pizeroE_dist[j])
    print("The Chi2 for the prediction found from the matrix in",xname,"is",str(chi2))
    return y_recovered_hist,x_predicted_hist,reco_y_error_hist,correlation_hist,chi2,threshold,refold_pizeroE_hist

def CalChi2(prediction,error,data): # These should all be numpy arrays
    prediction = np.asarray(prediction)
    data = np.asarray(data)
    error = np.asarray(error)
    chi2 = np.sum((prediction - data)**2 / error**2)
    return chi2

def HistCalChi2(prediction,data):
    data = data.GetCVHistoWithError()
    chi2 = 0
    for i in range(data.GetNbinsX()):
        chi2 += (prediction.GetBinContent(i+1)-data.GetBinContent(i+1))**2/(data.GetBinError(i+1)**2)
    return chi2


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
        except KeyError:
            pots[i]=None
    pot_scale = pots[0]/pots[1] if pots.count(None) == 0 else 1.0
    print("POT Scale:",pot_scale, "Data POT:",pots[0],"Sim POT:",pots[1])
    return pot_scale


playlist = AnalysisConfig.playlist

if __name__ == "__main__":
    subtractedDataPath = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_subtracted_collab1_fspline.root"
    migrationPath = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
    efficiencyPath = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/"+str(playlist)+"_calculatedEfficiency.root"
    pion_efficiency = Variable1D("pion_efficiency",efficiencyPath,isEfficiency = True)
    Eel = Variable1D("Eel",subtractedDataPath)
    pion_migration = MigrationMatrix("ee_vs_pionE",migrationPath)
    #pion_migration.hist.AddMissingErrorBandsAndFillWithCV(Eel.hist)
    migrationFile = ROOT.TFile.Open(migrationPath,"READ")
    PiZeroE_NCCOhPi0_dist = migrationFile.PiZeroE_NCCohPi0.Clone()
    PiZeroE_NCCOhPi0_dist.AddMissingErrorBandsAndFillWithCV(Eel.hist)
    true_pionE = migrationFile.true_pizeroE.Clone()
    true_pionE.AddMissingErrorBandsAndFillWithCV(Eel.hist)
    data_path = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/kin_dist_data"+str(playlist)+"_collab1_fspline.root"

    pot_scale = GetPOTScale(data_path,migrationPath)
    PiZeroE_NCCOhPi0_dist.Scale(pot_scale)
    true_pionE.Scale(pot_scale)
    pion_recovered_hist, Eel_matrixpredicted_hist,sqrt_diag_cov,correlation_hist,chi2,threshold,refold_pizeroE_hist = MatrixTruthRecovery(Eel,pion_migration,pion_efficiency,"Eel","pionE",PiZeroE_NCCOhPi0_dist)


    '''Eel_chi2 = HistCalChi2(pion_fitted_hist,pion_recovered_hist)
    print("Agreement Chi2s:"Eel:",str(Eel_chi2))'''

    #outPath = "/minerva/data/users/nalex/nu_e/Nimmy_04_18_2024_nx_unfolded.root"
    outPath = "/exp/minerva/data/users/nalex/nu_e2/ROOTFILES/Coh/"+str(playlist)+"_unfolded.root"
    outFile = ROOT.TFile.Open(outPath,"RECREATE")
    pion_recovered_hist.Write()
    Eel_matrixpredicted_hist.Write()
    outFile.Close()


    pion_canvas = ROOT.TCanvas("c2","c2",1600,1200)
    pion_recovered_err = pion_recovered_hist.GetCVHistoWithError()
    pion_recovered_err.Scale(1,"width")
    pion_recovered_err.SetMarkerStyle(21)
    pion_recovered_err.SetMarkerColor(4)
    pion_recovered_err.SetMarkerSize(2)
    pion_recovered_err.SetLineWidth(2)
    pion_recovered_err.SetLineColor(4)
    pion_recovered_err.SetTitle("Unfolded data")
    pion_recovered_err.GetYaxis().SetTitle("dNEvents/dE_{#pi^{0}}")
    pion_recovered_err.GetXaxis().SetTitle("Recovered Pion Energy (GeV)")

    truePiZeroE = true_pionE.GetCVHistoWithError()
    truePiZeroE.Scale(1,"width")
    truePiZeroE.SetMarkerStyle(21)
    truePiZeroE.SetMarkerColor(20)
    truePiZeroE.SetMarkerSize(2)
    truePiZeroE.SetLineWidth(2)
    truePiZeroE.SetLineColor(20)
    truePiZeroE.GetXaxis().SetTitle("True Pion Energy (GeV)")
    truePiZeroE.SetTitle("Simulation")


    PiZeroE_NCCohPi0_err = PiZeroE_NCCOhPi0_dist.GetCVHistoWithError()
    PiZeroE_NCCohPi0_err.Scale(1,"width")
    PiZeroE_NCCohPi0_err.SetMarkerStyle(21)
    PiZeroE_NCCohPi0_err.SetMarkerColor(8)
    PiZeroE_NCCohPi0_err.SetMarkerSize(2)
    PiZeroE_NCCohPi0_err.SetLineWidth(2)
    PiZeroE_NCCohPi0_err.SetLineColor(8)
    PiZeroE_NCCohPi0_err.SetTitle("dNEvents/dE_{#pi^{0}}")
    PiZeroE_NCCohPi0_err.GetXaxis().SetTitle("Pion Energy (GeV)")
    PiZeroE_NCCohPi0_err.SetTitle("Simulation")

    pion_recovered_err.Draw("MIN0 E1")
    PiZeroE_NCCohPi0_err.Draw("MIN0 E1 SAME")
    #truePiZeroE.Draw("MIN0 E1 SAME")
    legend = ROOT.gPad.BuildLegend(0.35,0.7,0.55,0.9)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    pion_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_recovered_"+str(playlist)+".png","png")

    pion_canvas.SetLogy()
    pion_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_recovered_"+str(playlist)+"_log.png","png")


    pion_fit_accuracy_canvas = ROOT.TCanvas("c4","c4",1600,1200)
    Eel_matrixpredicted_err = Eel_matrixpredicted_hist.GetCVHistoWithError()
    Eel_matrixpredicted_err.Scale(1,"width")
    Eel_matrixpredicted_err.SetMarkerStyle(21)
    Eel_matrixpredicted_err.SetMarkerColor(2)
    Eel_matrixpredicted_err.SetLineColor(2)
    Eel_matrixpredicted_err.SetMarkerSize(2)
    Eel_matrixpredicted_err.SetLineWidth(2)
    Eel_matrixpredicted_err.SetTitle("Refolded obtained pionE (GeV)")

    Eel_data_err = Eel.hist.GetCVHistoWithError()
    Eel_data_err.Scale(1,"width")
    Eel_data_err.SetMarkerStyle(21)
    Eel_data_err.SetMarkerColor(8)
    Eel_data_err.SetLineColor(8)
    Eel_data_err.SetMarkerSize(2)
    Eel_data_err.SetLineWidth(2)
    Eel_data_err.GetXaxis().SetTitleSize(0.04)
    Eel_data_err.GetYaxis().SetTitle("NEvents/GeV")
    Eel_data_err.SetTitle("Subtracted Electron Cone Energy (GeV)")
    Eel_data_err.GetXaxis().SetTitle("GeV")

    refold_pizeroE_err = refold_pizeroE_hist.GetCVHistoWithError()
    refold_pizeroE_err.Scale(1,"width")
    refold_pizeroE_err.SetMarkerStyle(21)
    refold_pizeroE_err.SetMarkerColor(4)
    refold_pizeroE_err.SetLineColor(4)
    refold_pizeroE_err.SetMarkerSize(2)
    refold_pizeroE_err.SetLineWidth(2)
    refold_pizeroE_err.GetXaxis().SetTitleSize(0.04)
    refold_pizeroE_err.GetXaxis().SetTitle("GeV")
    refold_pizeroE_err.GetYaxis().SetTitle("NEvents/GeV")
    refold_pizeroE_err.SetTitle("Refolded PiZeroE_NCCohPi0 from mc(GeV)")


    Eel_data_err.Draw("MIN0 E1")
    Eel_matrixpredicted_err.Draw("SAME MIN0 E1")
    refold_pizeroE_err.Draw("SAME MIN0 E1")
    legend = ROOT.gPad.BuildLegend(0.78,0.7,0.9,0.8)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    pion_fit_accuracy_canvas.SetTitle("Accuracy")
    pion_fit_accuracy_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_accuracy_"+str(playlist)+".png","png")


    pion_fit_accuracy_canvas.SetLogy()
    pion_fit_accuracy_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_accuracy_"+str(playlist)+"_log.png","png")


    pion_Eel_comp_canvas = ROOT.TCanvas("c4","c4",1600,1200)
    pion_recovered_err.Draw("MIN0 E1")
    Eel_matrixpredicted_err.Draw("SAME MIN0 E1 ")
    legend =ROOT.gPad.BuildLegend(0.35,0.7,0.55,0.8)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    pion_Eel_comp_canvas.SetTitle("Comparison")
    pion_Eel_comp_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_Eel_comparison_"+str(playlist)+".png","png")

    pion_Eel_comp_canvas.SetLogy()
    pion_Eel_comp_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/pion_Eel_comparison_"+str(playlist)+"_log.png","png")


    sqrt_diag_cov_canvas = ROOT.TCanvas("c4","c4",1600,1200)
    sqrt_diag_cov = sqrt_diag_cov#.GetCVHistoWithError()
    #sqrt_diag_cov.Scale(1,"width")
    sqrt_diag_cov.SetTitle(f"Threshold {threshold}, chi2 {chi2}")
    sqrt_diag_cov.SetMarkerStyle(21)
    sqrt_diag_cov.SetMarkerColor(9)
    sqrt_diag_cov.SetLineColor(9)
    sqrt_diag_cov.SetLineWidth(6)
    sqrt_diag_cov.Draw("HIST")
    title = ROOT.TLatex()
    title.SetNDC()  # Use normalized device coordinates
    title.SetTextSize(0.035)
    title.DrawLatex(0.33, 0.94, f"Threshold {threshold}, chi2 {chi2:.3g}")
    sqrt_diag_cov_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/sqrt_diag_cov_threshold_"+str(threshold)+"_"+str(playlist)+".png","png")


    correlation_hist_canvas = ROOT.TCanvas("c4","c4",1600,1200)
    correlation_hist = correlation_hist#.GetCVHistoWithError()
    correlation_hist.SetTitle(f"Threshold {threshold}, chi2 {chi2:.3g}")
    correlation_hist.SetMarkerStyle(21)
    correlation_hist.SetMarkerColor(9)
    correlation_hist.SetMinimum(-1)
    correlation_hist.SetMaximum(1)
    correlation_hist.SetLineColor(9)
    #ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetPalette(ROOT.kCubehelix)
    correlation_hist.Draw("COLZ")
    #ROOT.gPad.BuildLegend(0.45,0.4,0.85,0.6)
    title = ROOT.TLatex()
    title.SetNDC()  # Use normalized device coordinates
    title.SetTextSize(0.035)
    title.DrawLatex(0.33, 0.94, f"Threshold {threshold}, chi2 {chi2:.3g}") 

    correlation_hist_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/correlation_hist_threshold_"+str(threshold)+"_"+str(playlist)+".png","png")


    Unfolded_error_canvas = ROOT.TCanvas("c1","c1",1600,1200)
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
    PlotTools.MakeErrPlot(pion_recovered_hist)
    PlotTools.MNVPLOTTER.axis_maximum = 0.2
    Unfolded_error_canvas.Print("/exp/minerva/data/users/nalex/nu_e2/results/NCCoh/UnfoldedpionE_error_"+str(playlist)+".png","png")


