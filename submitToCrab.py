import os,sys
import string, re
from time import gmtime, localtime, strftime


channels  = ["Summer09-7TeV-ZeeJets-Pt_0_15",
            "Summer09-7TeV-ZeeJets-Pt_15_20",
            "Summer09-7TeV-ZeeJets-Pt_20_30",
            "Summer09-7TeV-ZeeJets-Pt_30_50",
            "Summer09-7TeV-ZeeJets-Pt_50_80",
            "Summer09-7TeV-ZeeJets-Pt_80_120",
            "Summer09-7TeV-ZeeJets-Pt_120_170",
            "Summer09-7TeV-ZeeJets-Pt_170_230",
            "Summer09-7TeV-ZeeJets-Pt_230_300",
            "Summer09-7TeV-ZeeJets-Pt_300_Inf"]


## channels  = ["Summer09-10TeV-ZeeJets-Pt_0_15",
##              "Summer09-10TeV-ZeeJets-Pt_15_20",
##              "Summer09-10TeV-ZeeJets-Pt_20_30",
##              "Summer09-10TeV-ZeeJets-Pt_30_50",
##              "Summer09-10TeV-ZeeJets-Pt_50_80",
##              "Summer09-10TeV-ZeeJets-Pt_80_120",
##              "Summer09-10TeV-ZeeJets-Pt_120_170",
##              "Summer09-10TeV-ZeeJets-Pt_170_230",
##              "Summer09-10TeV-ZeeJets-Pt_230_300",
##              "Summer09-10TeV-ZeeJets-Pt_300_Inf"]





dataset  = ["/ZeeJet_Pt0to15/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt15to20/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt20to30/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt30to50/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt50to80/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt80to120/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt120to170/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt170to230/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt230to300/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO",
           "/ZeeJet_Pt300toInf/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO"]


## dataset  = ["/ZeeJet_Pt0to15/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt15to20/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt20to30/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt30to50/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt50to80/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt80to120/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt120to170/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt170to230/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt230to300/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
##             "/ZeeJet_Pt300toInf/Summer09-MC_31X_V3-v1/GEN-SIM-RECO"]





condor  = [1,1,1,1,1,1,1,1,1,1]

MyResilientArea = "/dmehmet/ZeeJet_Summer09"


def changeMainConfigFile(outfile):
    fin  = open("Electron_EDM_Ntuple_cfg.py")
    pset_cfg      = "py_" + outfile + ".py"
    outfile_root  = "Summer09_7TeV_" + outfile + ".root"
    fout = open(pset_cfg,"w")
    for line in fin.readlines():
        if  line.find("demo.root")!=-1:
            line=line.replace("demo.root",outfile_root)
        fout.write(line)
    print pset_cfg + " has been written.\n"


def changeCrabTemplateFile(outfile, index):
    fin  = open("crabTemplate.cfg")
    pset_cfg      = "py_" + outfile + ".py"
    pset_crab     = "crabjob_" + outfile + ".cfg"
    outfile_root  = "Summer09_7TeV_" + outfile + ".root"
    fout = open(pset_crab,"w")
    for line in fin.readlines():
        if  line.find("mydataset")!=-1:
            line=line.replace("mydataset",dataset[index])
        if line.find("myanalysis")!=-1:
            line=line.replace("myanalysis",pset_cfg)    
        if  line.find("myrootfile")!=-1:
            line=line.replace("myrootfile",outfile_root)
        if  line.find("myresilient")!=-1:
            line=line.replace("myresilient",MyResilientArea)    
        if line.find("glite")!=-1 and condor[index]!=0:
            line=line.replace("glite", "condor")        
        fout.write(line)        
    if condor[index]!=0:
        fout.write("ce_white_list = cmssrm.fnal.gov")
      
    print pset_crab + " has been written.\n"


    
###################
for i in range(len(channels)):
    changeMainConfigFile(channels[i])
    changeCrabTemplateFile(channels[i],i)

for i in range(len(channels)):
    #if i<9: continue
    submitcommand = "crab -create -cfg " + "crabjob_" + channels[i] + ".cfg"
    child   = os.system(submitcommand)
    child2   = os.system("crab -submit")
