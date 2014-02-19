import FWCore.ParameterSet.Config as cms
process = cms.Process("SKIM")
process.load('FWCore.MessageService.MessageLogger_cfi')


#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/Run2010A/JetMET/RECO/v4/000/144/114/F4E53E34-31B4-DF11-9A2C-0030487CD6D8.root'
    )
)

#############   Path       ###########################
process.ak5CaloJetsSel = cms.EDFilter("CaloJetSelector",  
    src = cms.InputTag("ak5CaloJets"),
    cut = cms.string('pt > 100.0 && eta<2.6 && eta>-2.6')
)

process.dijetMass = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("ak5CaloJetsSel ak5CaloJetsSel"), 
    checkCharge = cms.bool(False),                           
    cut   = cms.string("360 < mass < 7000"),
    filter = cms.bool(True),
    minN    = cms.int32(1)      
)


##-------- Jets Triggers --------------
process.HLTJets = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths           = cms.vstring('HLT_Jet100U'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(True)
)
process.skimPath = cms.Path(process.ak5CaloJetsSel*process.HLTJets*process.dijetMass)

#############   output module ########################
process.out = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('skimPath')), 
    fileName = cms.untracked.string('JetSkim.root')
)

process.p = cms.EndPath(process.out)


#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

