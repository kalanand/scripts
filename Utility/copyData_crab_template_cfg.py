import FWCore.ParameterSet.Config as cms
process = cms.Process("Transfer")
process.load('FWCore.MessageService.MessageLogger_cfi')


#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)


#############   output module ########################
process.output = cms.OutputModule("PoolOutputModule", 
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('exotica-EXODiPhoSkimOct09_7TeV.root')
)


process.p = cms.EndPath(process.output)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

