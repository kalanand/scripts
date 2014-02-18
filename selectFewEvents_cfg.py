import FWCore.ParameterSet.Config as cms
process = cms.Process("SKIM")
process.load('FWCore.MessageService.MessageLogger_cfi')


nEventsToKeep = 2000

inputFile = 'file:/uscms_data/d2/kalanand/ttbsm_53x_data_1_1_8sY.root'
outputFile = 'ttbsm_53x_data_1_1_8sY.root'


## inputFile = 'file:/uscms_data/d2/kalanand/ttbsm_52x_mc_128_1_lmR.root'
## outputFile = 'ttbsm_52x_mc_128_1_lmR.root'



#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nEventsToKeep)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFile)
)
#############   output module ########################
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(outputFile)
)
process.p = cms.EndPath(process.out)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

