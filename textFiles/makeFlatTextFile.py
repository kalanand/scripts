import os,sys
import string, re
from time import gmtime, localtime, strftime
from ROOT import gROOT

dataFile2010 = "Data_WenuJets_2010_36invpb.root"
outputTextFile2010 = "Data_WenuJets_2010"

dataFile2011 = "Data_WenuJets_2011_22invpb.root"
outputTextFile2011 = "Data_WenuJets_2011"

mcFileWJets = "MC_WJets_enujj.root"
outputTextFileWJets = "MC_WJets_WenuJets"


mcFileWW = "MC_WW_enujj.root"
outputTextFileWW = "MC_WW_WenuJets"

mcFileWZ = "MC_WZ_enujj.root"
outputTextFileWZ = "MC_WZ_WenuJets"

mcFileTT = "MC_TT_enujj.root"
outputTextFileTT = "MC_TT_WenuJets"


varsToSave = "event_met_pfmet:event_met_pfmetPhi:event_met_pfmetsignificance:event_met_pfsumet:event_nPV"
cuts = "W_electron_et>25. && event_met_pfmet>25. && W_mt>50. && JetPFCor_Pt[0]>20. && JetPFCor_Pt[1]>20."


####################################################
gROOT.Reset()
gROOT.ProcessLine("TFile::Open(\"" + dataFile2010 +"\")")
gROOT.ProcessLine("WJet->SetScanField(-1)")
gROOT.ProcessLine("WJet->Scan(\"" + varsToSave + "\",\"" + cuts + "\"); >" + outputTextFile2010)
gROOT.ProcessLine(".q")
os.system("./cleanscan.sh " + outputTextFile2010)
os.system("rm " + outputTextFile2010)
####################################################
gROOT.Reset()
gROOT.ProcessLine("TFile::Open(\"" + dataFile2011 +"\")")
gROOT.ProcessLine("WJet->SetScanField(-1)")
gROOT.ProcessLine("WJet->Scan(\"" + varsToSave + "\",\"" + cuts + "\"); >" + outputTextFile2011)
gROOT.ProcessLine(".q")
os.system("./cleanscan.sh " + outputTextFile2011)
os.system("rm " + outputTextFile2011)
####################################################
gROOT.Reset()
gROOT.ProcessLine("TFile::Open(\"" + mcFileWJets +"\")")
gROOT.ProcessLine("WJet->SetScanField(-1)")
gROOT.ProcessLine("WJet->Scan(\"" + varsToSave + "\",\"" + cuts + "\"); >" + outputTextFileWJets)
gROOT.ProcessLine(".q")
os.system("./cleanscan.sh " + outputTextFileWJets)
os.system("rm " + outputTextFileWJets)
####################################################
gROOT.Reset()
gROOT.ProcessLine("TFile::Open(\"" + mcFileWW +"\")")
gROOT.ProcessLine("WJet->SetScanField(-1)")
gROOT.ProcessLine("WJet->Scan(\"" + varsToSave + "\",\"" + cuts + "\"); >" + outputTextFileWW)
gROOT.ProcessLine(".q")
os.system("./cleanscan.sh " + outputTextFileWW)
os.system("rm " + outputTextFileWW)
####################################################
gROOT.Reset()
gROOT.ProcessLine("TFile::Open(\"" + mcFileWZ +"\")")
gROOT.ProcessLine("WJet->SetScanField(-1)")
gROOT.ProcessLine("WJet->Scan(\"" + varsToSave + "\",\"" + cuts + "\"); >" + outputTextFileWZ)
gROOT.ProcessLine(".q")
os.system("./cleanscan.sh " + outputTextFileWZ)
os.system("rm " + outputTextFileWZ)
####################################################
gROOT.Reset()
gROOT.ProcessLine("TFile::Open(\"" + mcFileTT +"\")")
gROOT.ProcessLine("WJet->SetScanField(-1)")
gROOT.ProcessLine("WJet->Scan(\"" + varsToSave + "\",\"" + cuts + "\"); >" + outputTextFileTT)
gROOT.ProcessLine(".q")
os.system("./cleanscan.sh " + outputTextFileTT)
os.system("rm " + outputTextFileTT)
####################################################
