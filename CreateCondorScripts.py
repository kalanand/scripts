#! /usr/bin/env python

#****************************************************************
#               Dr. Konstantinos Kousouris                      *
#                 Fermilab, Batavia, IL                         *
#****************************************************************
import os
#******************   definitions  **********************************
NALG = 6
SHORT_JET_ALG_LIST = ["IC5","AKT5","KT4","KT6","SC5","SC7"]
JET_ALG_LIST =       ["iterativeCone5","antikt5","kt4","kt6","sisCone5","sisCone7"]
JET_TYPE =           ["Calo","PF"] 

PREFIX =             "MCTruthTree"
SAMPLE_LIST=["QCDDiJet_Pt0to15",
             "QCDDiJet_Pt15to20",
             "QCDDiJet_Pt20to30",
	     "QCDDiJet_Pt30to50",
	     "QCDDiJet_Pt50to80",
	     "QCDDiJet_Pt80to120",
             "QCDDiJet_Pt120to170",
             "QCDDiJet_Pt170to230",
	     "QCDDiJet_Pt230to300",
	     "QCDDiJet_Pt300to380",
	     "QCDDiJet_Pt380to470",
             "QCDDiJet_Pt470to600",
             "QCDDiJet_Pt600to800",
	     "QCDDiJet_Pt800to1000",
	     "QCDDiJet_Pt1000to1400",
	     "QCDDiJet_Pt1400to1800",
             "QCDDiJet_Pt1800to2200",
             "QCDDiJet_Pt2200to2600",
	     "QCDDiJet_Pt2600to3000",
	     "QCDDiJet_Pt3000to3500",
	     "QCDDiJet_Pt3500toInf"]
	     
DATA_LIST = ["/QCDDiJet_Pt0to15/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
             "/QCDDiJet_Pt15to20/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt20to30/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt30to50/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt50to80/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt80to120/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt120to170/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt170to230/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt230to300/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt300to380/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt380to470/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt470to600/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt600to800/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt800to1000/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt1000to1400/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt1400to1800/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt1800to2200/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt2200to2600/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt2600to3000/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt3000to3500/Summer09-MC_31X_V3-v1/GEN-SIM-RECO",
	     "/QCDDiJet_Pt3500toInf/Summer09-MC_31X_V3-v1/GEN-SIM-RECO"]	     
#*********************************************************************
nfiles = 21
sample = 0
while sample < nfiles: #Starting the loop for different samples
    #************* Configuration file *********************
    cfg_body = "MCTruthTree_"+SAMPLE_LIST[sample]
    cfg_name = cfg_body+"_cfg.py" 
    file = open(cfg_name,'w')
    file.write("process = cms.Process(\"Ana\") \n")
    file.write("process.load(\"FWCore.MessageService.MessageLogger_cfi\") \n")
    file.write("process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1)) \n")
    file.write("process.load(\"JetMETCorrections.Configuration.JetPlusTrackCorrections_cff\") \n")
    file.write("process.source = cms.Source(\"PoolSource\",fileNames = cms.untracked.vstring()) \n")
    file.write("\n")
    i=0
    path = "process.ZSPJetCorrections*process.JetPlusTrackCorrections*"
    #################  JPT #####################################
    module = SHORT_JET_ALG_LIST[i]+"JPTMctruthTree"
    path+="process."+module+"*"
    phrase = "process."+module+" = cms.EDAnalyzer(\"CaloMCTruthTreeProducer\", \n"
    file.write(phrase)
    phrase = "     jets               = cms.string(\'JetPlusTrackZSPCorJetIcone5\'), \n"
    file.write(phrase)
    phrase = "     genjets            = cms.string(\'iterativeCone5GenJets\'), \n"
    file.write(phrase)
    histoname = module+"_"+SAMPLE_LIST[sample]+".root"
    phrase = "     histogramFile      = cms.string(\'"+module+"_"+SAMPLE_LIST[sample]+".root\') \n"
    file.write(phrase)  
    file.write(")\n")
    file.write("\n")
    while i<NALG:
       for jet in JET_TYPE:
          module = SHORT_JET_ALG_LIST[i]+jet+"MctruthTree"
          path+="process."+module+"*"
          phrase = "process."+module+" = cms.EDAnalyzer(\""+jet+"MCTruthTreeProducer\", \n"
          file.write(phrase)
          phrase = "     jets               = cms.string(\'"+JET_ALG_LIST[i]+jet+"Jets\'), \n"
          file.write(phrase)
          phrase = "     genjets            = cms.string(\'"+JET_ALG_LIST[i]+"GenJets\'), \n"
          file.write(phrase)
	  histoname+ = module+"_"+SAMPLE_LIST[sample]+".root"
          phrase = "     histogramFile      = cms.string(\'"+module+"_"+SAMPLE_LIST[sample]+".root\') \n"
          file.write(phrase)  
          file.write(")\n")
          file.write("\n") 
       i+=1
    phrase = "process.p = cms.Path("+path[:-1]+") \n"
    file.write(phrase)
    file.write("process.MessageLogger.cerr.FwkReport.reportEvery = 1000 \n")
    file.close()
    #####################  CRAB FILE ###################################### 
    crab_body = "crab_"+SAMPLE_LIST[sample]
    crab_name = crab_body+".cfg"
    print crab_name
    file = open(crab_name,'w')
    file.write("[CRAB]\n")
    file.write("jobtype                 = cmssw \n")
    file.write("scheduler               = glite \n")
    file.write("[CMSSW] \n")
    phrase = "datasetpath             = "+DATA_LIST[sample]+"\n"
    file.write(phrase)
    phrase = "pset                    = "+cfg_name+"\n" 
    file.write(phrase)
    file.write("total_number_of_events  = -1\n")
    file.write("events_per_job          = 5000\n")
    file.write("output_file             = ")
    histo_name = "DBTree"+"_"+FILE_LIST[sample]+".root \n"
    file.write(histo_name)
    file.write("[USER]\n")
    file.write("eMail                   = Konstantinos.Kousouris@cern.ch \n")
    file.write("return_data             = 0\n")
    file.write("copy_data               = 1\n")
    file.write("storage_element         = cmssrm.fnal.gov\n")
    file.write("storage_path            = /srm/managerv2?SFN=/resilient/kkousour/Summer09/\n")
    phrase = "user_remote_dir         = "+FILE_LIST[sample]+"\n"
    file.write(phrase)
    file.write("srm_version             = 2\n")
    file.write("use_central_bossDB      = 0\n")
    file.write("use_boss_rt             = 0\n")
    file.write("[EDG]\n")
    file.write("rb                      = CERN\n") 
    file.write("proxy_server            = myproxy.cern.ch\n") 
    file.write("virtual_organization    = cms\n")
    file.write("retry_count             = 0\n")
    file.write("lcg_catalog_type        = lfc\n")
    file.write("lfc_host                = lfc-cms-test.cern.ch\n")
    file.write("lfc_home                = /grid/cms\n")
    file.close()
    sample +=1
print "Done........" #End of sample loop 
