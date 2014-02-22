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
process.L1T1coll.L1SeedsLogicalExpression = cms.string('0 AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')

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





###############################################################################
############## We should require all the above for every path ################
process.goodCollSequence=cms.sequence( process.L1T1coll +
                                       process.primaryVertexFilter +
                                       process.noscraping )
###############################################################################




#  SuperClusters  ################

process.superClusters = cms.EDFilter("SuperClusterMerger",
   src = cms.VInputTag(cms.InputTag("hybridSuperClusters","", "RECO"),
                       cms.InputTag("multi5x5SuperClustersWithPreshower","", "RECO"))  
)



process.mySuperClusters = cms.EDFilter("SuperClusterSelector",
    src = cms.InputTag("superClusters"),
    cut = cms.string('energy*sin(position.theta)>15.0')
)



process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
   src = cms.InputTag("mySuperClusters"),
   particleType = cms.int32(11),
)




process.ZToEE = cms.EDFilter("NamedCandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 120'),
    name = cms.string('ZToEE'),
    roles = cms.vstring('electron1','electron2'),
    decay = cms.string('superClusterCands superClusterCands'),
    checkCharge = cms.bool(False)                  
)


process.ZeeSequence = cms.Sequence( process.superClusters * process.mySuperClusters *
                                    process.superClusterCands * process.ZToEE
                                    )




process.mypath = cms.Path(process.goodCollSequence*process.noscraping*
                              process.ZeeSequence)




#############  output module if just want to skim by HLT path ##############
process.out = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('process.')
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










