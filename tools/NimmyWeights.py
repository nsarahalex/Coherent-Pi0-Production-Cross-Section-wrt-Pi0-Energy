import os
import sys
import ROOT
from ROOT import *
import PlotUtils
from functools import partial
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.PlotLibrary import HistHolder
from PlotUtils.HistWrapper import HistWrapper
import tools.TruthTools as TruthTools
from config.PlotConfig import ELECTRON_ENERGY_BINNING as bins
import numpy as np
from tools import fitBackground

IsCC = lambda event: event.mc_current == 1
IsNC = lambda event: event.mc_current == 2
IsCoherent = lambda event: event.mc_intType == 4
IsNuE = lambda event: abs(event.mc_incoming) == 12
IsPi0InFinalState = lambda event: 111 in event.mc_FSPartPDG

IsHeavyBaryon = lambda event: 3112 in event.mc_FSPartPDG or 3122 in event.mc_FSPartPDG or 3212 in event.mc_FSPartPDG or 3222 in event.mc_FSPartPDG or 4112 in event.mc_FSPartPDG or 4122 in event.mc_FSPartPDG or 4212 in event.mc_FSPartPDG or 4222 in event.mc_FSPartPDG
IsMeson = lambda event: 211 in event.mc_FSPartPDG or -211 in event.mc_FSPartPDG or 321 in event.mc_FSPartPDG or -321 in event.mc_FSPartPDG or 323 in event.mc_FSPartPDG or -323 in event.mc_FSPartPDG  or 111 in event.mc_FSPartPDG or 130 in event.mc_FSPartPDG or 310 in event.mc_FSPartPDG or 311 in event.mc_FSPartPDG
IsPhoton = lambda event: 22 in event.mc_FSPartPDG and event.mc_FSPartPDG[0] != 22

IsCoherentEvent = lambda event: IsCoherent(event) and IsNC(event) and IsPi0InFinalState(event)
IsnonCCQEnu_e = lambda event: IsNuE(event) and (IsPi0InFinalState(event) or IsMeson(event) or IsPhoton(event) or IsHeavyBaryon(event))
IsCCQEnu_e = lambda event: IsCC(event) and IsNuE(event) and not IsPi0InFinalState(event) and not IsMeson(event) and not IsHeavyBaryon(event) and not IsPhoton(event)
IsNCPi0 = lambda event: IsNC(event) and IsPi0InFinalState(event)
IsCCPi0 = lambda event: IsCC(event) and IsPi0InFinalState(event)
IsNCDiff = lambda event: event.mc_intType == 10

if not AnalysisConfig.unweighted:
    mainweights,univweights = fitBackground.RunAlexWeighter()

def GetAlexWeights(event):
    doWeights = not AnalysisConfig.unweighted
    weight = 1
    Elep = event.kin_cal.reco_E_lep
    pionE = TruthTools.PiZeroE(event)
    if (IsCCQEnu_e(event) or IsnonCCQEnu_e(event)) and not IsNCDiff(event):
        eventType = 1
    elif IsCoherentEvent(event): # Coherent event
        eventType = 3  
    elif IsNCDiff(event): # Diffractive pion production
        eventType = 0
    elif IsNCPi0(event) or IsCCPi0(event): # Pi0 events
        eventType = 2
    else: return weight # Returns the default of weight is 1 if it doesn't match any of the event types

    if not doWeights: return weight # This needs to be here to escape if weights aren't needed
    if event.ShortName() == "cv":
        weights = mainweights
    else: weights = univweights[event.ShortName()][event.ithInWrapper]
    if eventType == 3: # NCCoh 
        pionE = TruthTools.PiZeroE(event)
        for i in range(len(bins)-1):
            if bins[i] <= pionE < bins[i+1]:
                weight *= weights[i,eventType]
                break
    else: #Non-signal uses Elep
        Elep = event.kin_cal.reco_E_lep
        for i in range(len(bins)-1):
            if Elep is not None:
               if bins[i] <= Elep < bins[i+1]:
                  weight *= weights[i,eventType]
                  break
    return weight

