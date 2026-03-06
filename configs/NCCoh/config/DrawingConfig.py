import ROOT
from collections import OrderedDict
from config import PlotConfig
from tools.PlotLibrary import HistHolder
from tools import PlotTools
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS,GENIE_ERROR_GROUPS



Default_Plot_Type="stacked"
Default_Scale = lambda histHolder:histHolder.POTScale(True)
DefaultSlicer = PlotTools.PrepareSlicer


# the order and color of each band in stacked plots defined here.
COLORS=ROOT.MnvColors.GetColors()
Categories = OrderedDict()
Categories["NCDiff"] = {
    "title":"NC Diffractive #pi^{0}",
    "color" : COLORS[2]
}
Categories["NCCohPi0"]= {
    "title": "NC Coherent #pi^{0}",
    "color" : COLORS[4]
}
Categories["lowNCPi0"]= {
    #"title": "NC Coh #pi^{0} , E_{#pi^{0}}<1.5GeV, inside cone",
    "title": "NC Coh #pi^{0} , E_{#pi^{0}}<1.5GeV",
    "color" : ROOT.kRed 
}
Categories["NCPi0"]= {
    "title": "NC Non-Coherent#pi^{0}",
    "color" : COLORS[5]
}
Categories["NCDiff_NPS"]= {
    "title": "NCDiff_NPS",
    "color" : ROOT.kRed+3
}


Categories["CCPi0"]= {
    "title": "Not #nu_{e} CC #pi^{0}",
    "color" : COLORS[3]
}
##Categories["CCPi0mid"]= {    
##    "title": "Not #nu_{e} CC #pi^{0} Pi0",
##    "color" : COLORS[6]
##}
##Categories["CCPi0high"]= {  #added for Eavail
##    "title": "Not #nu_{e} CC #pi^{0} Excluded",
##    "color" : COLORS[2]
##}
Categories["CCQElike"]={
    "title" : "CCQElike #nu_{e}",
    "color": COLORS[0]
}
Categories["notCCQElike"]={
    "title" : "notCCQElike #nu_{e}",
    "color": COLORS[1]
}
Categories["NuEElastic"] = {
    "title":"#nu+e elastic",
    "color" : COLORS[6]
}
#Categories["Other"] = {
#    "title":"Others",
#    #"cate" : {"NuEElastic"},
#    
#    ##    "cate" : {"Other" ,"NonPhaseSpace","NonFiducial","NuEElastic"},
#    "color" : COLORS[7]
#}
Categories["NonPhaseSpace"] = {
    "title":"NC Coh #pi^{0},E_{#pi^{0}}>1.5GeV,#theta_{#pi^{0}}>30^{o}",
    "color" : ROOT.kMagenta
}
Categories["lowNonPhaseSpace"] = {
    "title":"NC Coh #pi^{0},E_{#pi^{0}}<1.5GeV,#theta_{#pi^{0}}>30^{o}",
    "color" : ROOT.kTeal
}
Categories["NonFiducial"] = {
    "title":"NonFiducial",
    "color" : ROOT.kGray
}
Categories["Other"] = {
    "title":"Others",
    "cate" : {"Other"},# ,"NonPhaseSpace","lowNonPhaseSpace"},
    "color" : COLORS[7]
}

SignalDecomposition = {
    "CCNuEQE" :
    {
        "title" : "CC #nu_{e}-QE",
        "color": COLORS[0]
    },
    "CCNuEDelta" : {
        "title" : "CC #nu_{e}-Delta",
        "color": COLORS[1]
    },
    "CCNuEDIS" : {
        "title" : "CC #nu_{e}-DIS",
        "color": COLORS[2]
    },
    "CCNuE2p2h" : {
        "title" : "CC #nu_{e}-2p2h",
        "color": COLORS[5]
    },
    "CCNuE": {
        "title" : "CC #nu_{e}-Other",
        "color": COLORS[3]
    },
    "Background" : {
        "title" : "Backgrounds",
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCCOH","NCRES","NCDIS","NCOther","NuEElastic","Other","CCNuEAntiNu","NonPhaseSpace"},
        "color": COLORS[4]
    }
    # "CCNuEAntiNu" : {
    #     "title" : "CC #bar{#nu_{e}}",
    #     "color": COLORS[4]
    # }
}

SignalBackground = {
    "Signal" :
    {
        "title" : "CC #nu_{e}",
        "cate": {
            "CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"
        },
        "color": COLORS[0]
    },
    "Background" : {
        "title" : "Backgrounds",
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCCOH","NCRES","NCDIS","NCOther","NuEElastic","Other","CCNuEAntiNu","NonPhaseSpace"},
        "color": COLORS[4]
    }
}

NuEElasticCategory = {
    "Nu+e" :{
        "title": "nu+e",
        "cate" :{"NuEElastic"},
        "color":COLORS[1]
    },
    "CCNuE":
    {
        "title" : "CC #nu_{e}",
        "cate": {
            "CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"
        },
        "color": COLORS[0]
    },
    "NCCOH":{
        "title" : "NCCoh",
        "color": COLORS[3]
    },
    "ExcessModel": {
        "title": "NC Diffractive",
        "color":COLORS[2]
    },
    "Other":{
        "title" : "Others",
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCRES","NCDIS","NCOther","Other","CCNuEAntiNu","NonPhaseSpace"},
        "color": COLORS[4]
    }
}

MODELS = {
    "MnvTune v1": {
        "errorband":(None,None),
        "color":COLORS[0]
    },
    "2p2h Tune (QE)": {
        "errorband":("Low_Recoil_2p2h_Tune",2),
        "color":COLORS[1]
    },
    "SuSA 2p2h" : {
        "errorband":("SuSA_Valencia_Weight",0),
        "color":COLORS[2]
    },
    "MK model": {
        "errorband": ("MK_model",0),
        "color":COLORS[7]
    },
    "Low Q2 Pion Joint": {
        "errorband" : ("LowQ2Pi_Joint",0),
        "color": COLORS[5]
    },
    "Low Q2 Pion NuPi0": {
        "errorband" : ("LowQ2Pi_NUPI0",0),
        "color": COLORS[6]
    }
}

DefaultPlotters={
    "comp":{"func": PlotTools.PrepareComp},
    "ratio":{"func": PlotTools.PrepareRatio},
    "err": {"func": PlotTools.PrepareErr,
            #"args":(True,False,GENIE_ERROR_GROUPS)},    
            "args":(True,False,CONSOLIDATED_ERROR_GROUPS)},
    "stacked":{"func":PlotTools.PrepareStack,
               "args": (Categories,)},
    "diff":{"func":PlotTools.PrepareDiff},
    "migration":{"func":PlotTools.PrepareMigration},
    "sigdep":{"func":PlotTools.PrepareSignalDecompose,
              "args": (SignalDecomposition,True,False)},
    "errband":{"func":PlotTools.PrepareErrorBand},
    "model":{},
    "model_ratio" :{},
    "sigdep":{},

    "mcstack":{"func":PlotTools.PrepareMCStackforNimmy,
               "args": (Categories,)},
    
    "mconlystack":{"func":PlotTools.PrepareMCStackforNimmy_v2,
               "args": (Categories,)},
    }

PLOTS_TO_MAKE = [
#     {"name":"Lackner Migration",
#      "plot_type": "migration"},
#     {"name":"IUE vs Lackner Reco",
#      "plot_type":"migration"},
#     {"name":"t vs Lackner Truth",
#      "plot_type":"migration"},
#     {"name":"Inline Upstream Energy vs Inline Upstream Energy Weighted Position",
#      "slicer" :lambda hist: PlotTools.Make2DSlice(hist,False,bin_start = 2,bin_end = 0,interval = 1),
#      "plot_type":"stacked",
      #"canvasconfig" : PlotTools.Logx,
      #"canvasconfig" : PlotTools.Logy,
#      },
     {"name":"Inline Upstream Energy Weighted Position",
     "plot_type":"stacked",
     "canvasconfig" : PlotTools.Logx,
      },
     {"name":"Electron Energy vs Pion Energy",
     "plot_type":"migration"},
#    {"name":"Electron Energy vs Pion Energy",
#     "plot_type":"stacked",
##     "slicer" :lambda hist:PlotTools.Make2DSlice(hist,True,bin_start = 0,interval = 1),
##     "tag": "stack"},
##    {"name":"Neutrino Energy",
##      "plot_type":"stacked",
##      "scale" : lambda histHolder:histHolder.POTScale(True)},
##    {"name":"Lepton Energy",
##      "plot_type":"ratio"},
    {"name":"Lepton Energy",
      "plot_type":"err"},
##    {"name":"Lepton Energy",
##      "plot_type":"comp"},
##     {"name":"PiZero E NCCoherent",
##      "plot_type":"err"},
##     {"name":"PiZero E NCCoherent",
##      "plot_type":"stacked",
##      #"canvasconfig" : PlotTools.Logx},
##      },
##    {"name":"PiZeroE",
##     "plot_type":"stacked"},
##    {"name":"Lepton Pt",
##      "plot_type": "stacked",
##      "scale" : lambda histHolder:histHolder.POTScale(True)},
##    {"name":"Lepton Energy",
##      "plot_type":"diff"}, 
##    {"name":"Lepton Energy",
##      "plot_type":"sigdep"},
##     {"name":"Lepton Energy",
##      "plot_type":"stacked",
##      "scale" : lambda histHolder:histHolder.POTScale(False),
##      "tag" : "Not BWN"},
     {"name":"Lepton Energy",
      "plot_type":"stacked",
      "scale" : lambda histHolder:histHolder.POTScale(True)},
##
##    # {"variables":["Lepton Visible Energy","Visible Energy"]},

    {"name":"Epi(1-cos(pi))",
     "plot_type":"stacked",
     "scale" : lambda histHolder:histHolder.POTScale(False)},
##    {"name":"Weighted_Epi(1-cos(pi))",
#     "plot_type":"stacked",
#     "scale" : lambda histHolder:histHolder.POTScale(True),
     #"slicer" :lambda hist:PlotTools.Make2DSlice(hist,True,bin_start = 0,interval = 1),
#     },
     #     "scale" : lambda histHolder:histHolder.POTScale(True)},
#    {"name":"PiZero E NC_Coherent_diff_bin",
#     "plot_type":"stacked",
#     "scale" : lambda histHolder:histHolder.POTScale(False),
#     "tag" : "not_BWN"},
#    {"name":"PiZero E NC_Coherent_diff_bin",
#     "plot_type":"stacked",
#     "scale" : lambda histHolder:histHolder.POTScale(True),
#     "tag" : "BWN"},
#    {"name":"Lepton Energy vs pionE",
#     "plot_type": "migration",
#     #"slicer":lambda hist:PlotTools.Make2DSlice(hist,True,bin_start = -1,interval = 1),
#     #"scale" : lambda histHolder:histHolder.POTScale(False)},
#     },
#     {"name": "Electron Energy vs Pion Energy New",
#     "plot_type": "migration",
#      "tag" : "updated"},
#     {"name": "Electron Energy vs Pion Energy New",
#     "plot_type": "stacked"},
#     {"name": "Lepton Energy vs true_Elep",
#     "plot_type": "migration"},
#    {"name": "True Signal Pi Zero Energy",
#     "plot_type": "stacked"}, 
    {"name": "True Signal Pi Zero Energy",
     "plot_type": "err"}, 
#    {"name": "True Signal Pi Zero Energy_diff_bin",
#     "plot_type": "stacked"},
#    {"name": "Reco vs True PionE",
#     "plot_type": "migration"},
#    {"name": "Reco vs True Pion Theta",
#     "plot_type": "migration"},
#    {"name": "True Pi0",
#     "plot_type": "stacked"},
#    {"name": "Reco Lepton Energy",
#     "plot_type": "stacked"},
#    {"name": "Pi0 True theta",
#     "plot_type": "stacked"},
#    {"name": "Reco Lep theta",
#     "plot_type": "stacked"},
    ]




