import os,sys
import string, re
from time import gmtime, localtime, strftime


channels  = ["0_15","15_20","20_30","30_50",
             "50_80","80_120","120_170",
             "170_230","230_300","300_Inf"]



for i in range(len(channels)):
    submitcommand3 = "hadd Spring10_7TeV_Spring10-7TeV-ZeeJets-Pt_" + channels[i] + ".root" + " " + "Spring10_7TeV_Spring10-7TeV-ZeeJets-Pt_" + channels[i] + "_*.root"
    child   = os.system(submitcommand3)
    submitcommand4 = "rm Spring10_7TeV_Spring10-7TeV-ZeeJets-Pt_" + channels[i] + "_*.root"
    child2   = os.system(submitcommand4)
