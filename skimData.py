import FWCore.ParameterSet.Config as cms
process = cms.Process("SKIM")

#------------------------------------------
# Load standard sequences.
#------------------------------------------
process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR10_P_V4::All' 


## process.load('Configuration.StandardSequences.MagneticField_38T_cff')



## Technical Trigger bits, beam coincidence/ BPTX, etc ...
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.L1T1coll=process.hltLevel1GTSeed.clone()
process.L1T1coll.L1TechTriggerSeeding = cms.bool(True)
process.L1T1coll.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')

process.l1tcollpath = cms.Path(process.L1T1coll)

process.primaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
    filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.noscraping = cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.25)
)



# Physics declared 
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'


###############################################################################
############## We should require all the above for every path ################
process.goodCollSequence=cms.sequence( process.L1T1coll +
                                       process.primaryVertexFilter +
                                       process.noscraping +
                                       process.hltPhysicsDeclared )
###############################################################################



##-------- Muon events of interest --------
process.HLTMu =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring("HLT_L1MuOpen", "HLT_L1Mu", "HLT_Mu3", "HLT_Mu5",
                "HLT_Mu9", "HLT_Mu11", "HLT_Mu13", "HLT_Mu15", "HLT_Mu15_L1Mu7",
                "HLT_L2Mu9", "HLT_IsoMu11", "HLT_IsoMu13", "HLT_IsoMu15"),
                # provide list of HLT paths (or patterns) you want
     eventSetupPathsKey = cms.string(''),
      # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
     andOr = cms.bool(True),
     # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
     throw = cms.bool(False)    # throw exception on unknown path names
 )
process.pathHLTMu = cms.Path(process.goodCollSequence+process.HLTMu)




##-------- Electron events of interest --------
process.HLTEle =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring("HLT_L1SingleEG5","HLT_L1SingleEG8",
                            "HLT_L1DoubleEG5","HLT_Photon10_L1R",
                            "HLT_Photon15_L1R","HLT_Photon15_LooseEcalIso_L1R",
                            "HLT_Photon20_L1R"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
 )
process.pathHLTEle = cms.Path(process.goodCollSequence+process.noscraping+process.HLTEle)





##-------- Jet events of interest --------
process.HLTJets = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths           = cms.vstring('HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30','HLT_Jet50','HLT_Jet50'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(True)  # throw exception on unknown path names
)

process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_ReReco332_cff")
process.antikt5CaloJetsCor  = cms.EDFilter("CaloJetRefSelector",  
    src = cms.InputTag("L2L3CorJetAK5Calo"),
    cut = cms.string('pt > 20.0 & abs( eta ) < 3.0')
)
process.pathHLJet = cms.Path(process.HLTJets*process.L2L3CorJetAK5Calo*process.antikt5CaloJetsCor)





#############  output module if just want to skim by HLT path ##############
process.out = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('process.pathHLTMu','process.pathHLTEle','process.pathHLTJet' )
       ), 
    fileName = cms.untracked.string('mySkim.root')
)
process.p = cms.EndPath(process.out)


#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000







#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     'dcache:/pnfs/cms/WAX/resilient/kalanand/data/BeamCommissioning/data_2010_run_132601_1_1.root')
)










