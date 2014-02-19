import FWCore.ParameterSet.Config as cms
process = cms.Process("COPY")

#------------------------------------------
# Load standard sequences.
#------------------------------------------
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000



#############  output module ##############
process.out = cms.OutputModule("PoolOutputModule",
   outputCommands = cms.untracked.vstring('keep *'),                           
   fileName = cms.untracked.string('data_copied.root')
)
process.p = cms.EndPath(process.out)




#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20000)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)










