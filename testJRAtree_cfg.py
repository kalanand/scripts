import FWCore.ParameterSet.Config as cms
import JetMETAnalysis.JetAnalyzers.JRA_Defaults_cff as Defaults;
import JetMETAnalysis.JetAnalyzers.JRA_TreeDefaults_cff as Tree;

Defaults.JetPtEta = cms.PSet(
    etaMin = cms.double(-5.0),
    etaMax = cms.double(5.0),
    ptMin  = cms.double(1.0)
)
Defaults.RefPtEta = cms.PSet(
    etaMin = cms.double(-5.0),
    etaMax = cms.double(5.0),
    ptMin = cms.double(1.0)
)
Defaults.JetResponseParameters = Tree.JetResponseParameters


process = cms.Process("JRA")
process.load("JetMETAnalysis.JetAnalyzers.JRA_TreeDefaults_cff")
process.load("JetMETAnalysis.JetAnalyzers.JRA_Paths_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("JetMETCorrections.Configuration.JetPlusTrackCorrections_cff")
process.load("JetMETCorrections.Configuration.ZSPJetCorrections219_cff")
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Summer09/QCD_Pt300/GEN-SIM-RECO/MC_31X_V3-v1/0026/CC1E3FD7-EE89-DE11-A34E-001EC9AA9E78.root'
)
)
process.TFileService = cms.Service("TFileService",
    fileName      = cms.string('JRAt.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.recoJets = cms.Path(process.ZSPJetCorrections+process.JetPlusTrackCorrections)

process.schedule = cms.Schedule(
    process.recoJets,
    # uncorrected jets
    process.ic5jptJRA, 
    process.kt4caloJRA,
    process.kt6caloJRA,
    process.sc5caloJRA,
    process.sc7caloJRA,
    process.ic5caloJRA,
    process.ak5caloJRA,
    process.kt4pfJRA,
    process.kt6pfJRA,
    process.sc5pfJRA,
    process.sc7pfJRA,
    process.ic5pfJRA,
    process.ak5pfJRA
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10
