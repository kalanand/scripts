# Auto generated configuration file
# using: 
# Revision: 1.123 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/PYTHIA6_QCDplus3TeVcontact_pt_80_120_7TeV_cff.py -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,HLT --conditions FrontierConditions_GlobalTag,STARTUP31X_V2::All --fileout GenHLT_8E29.root --number 100 --mc --no_exec --datatier GEN-SIM-RAW --eventcontent RAWSIM --processName HLT8E29
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT8E29')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('HLTrigger/Configuration/HLT_8E29_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('QCD+3TeVcontact-pt-80-120 at 7TeV'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/GenProduction/python/PYTHIA6_QCDplus3TeVcontact_pt_80_120_7TeV_cff.py,v $')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("EmptySource")

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('GenHLT_8E29.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'STARTUP31X_V2::All'
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    crossSection = cms.untracked.double(783700),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(51)=10042     ! CTEQ6L1 structure function chosen', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model', 
            'MSTU(21)=1     ! Check on possible errors during program execution', 
            'PARP(82)=1.8387   ! pt cutoff for multiparton interactions', 
            'PARP(89)=1960. ! sqrts for which PARP82 is set', 
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter', 
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter', 
            'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
            'PARP(67)=2.5    ! amount of initial-state radiation', 
            'PARP(85)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(86)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(62)=1.25   ! ', 
            'PARP(64)=0.2    ! ', 
            'MSTP(91)=1     !', 
            'PARP(91)=2.1   ! kt distribution', 
            'PARP(93)=15.0  ! '),
        processParameters = cms.vstring('MSEL=0    !(D=1) to select between full user control (0, then use MSUB) and some preprogrammed alternative', 
            'ITCM(5)=2      ! Switch on contact inteactions for all quarks', 
            'RTCM(41)=5000  ! Set Contact Scale Lambda to 3 TeV', 
            'RTCM(42)=1     ! Sign of contact interaction is +', 
            'MSUB(381)=1    ! qi qj -> qi qj via QCD plus a contact interaction', 
            'MSUB(382)=1    ! qi qiBar -> qk qkBar via QCD plus a contact interaction', 
            'MSUB(13)=1     ! qi qiBar -> g g via normal QCD', 
            'MSUB(28)=1     ! qi g -> qi g  via normal QCD', 
            'MSUB(53)=1     ! g g -> qk qkbar via normal QCD', 
            'MSUB(68)=1     ! g g -> g g via normal QCD', 
            'CKIN(3)=230  ! minimum pt hat for hard interactions', 
            'CKIN(4)=300  ! maximum pt hat for hard interactions'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)
process.ProductionFilterSequence = cms.Sequence(process.generator)

process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.JetProducers.ak5GenJets_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')


# Path and EndPath definitions
#process.generation_step = cms.Path(process.pgen)
process.generation_step = cms.Path(process.genParticles *process.genParticlesForJets *process.ak5GenJets *process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step)
#process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.out_step])
# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq


### Next chained config file ###

# Auto generated configuration file
# using: 
# Revision: 1.123 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step2 --step=DIGI,L1,DIGI2RAW,HLT --conditions=FrontierConditions_GlobalTag,MC_31X_V3::All --number=100 --mc --datatier GEN-SIM-RAW --eventcontent=RAWSIM --processName=HLT --no_exec --filein file:GenHLT_8E29.root --fileout GenHLT_8E29_1E31.root
#import FWCore.ParameterSet.Config as cms
#
#process = cms.Process('HLT')
#
## import of standard configurations
#process.load('Configuration/StandardSequences/Services_cff')
#process.load('FWCore/MessageService/MessageLogger_cfi')
#process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
#process.load('Configuration/StandardSequences/GeometryIdeal_cff')
#process.load('Configuration/StandardSequences/MagneticField_38T_cff')
#process.load('Configuration/StandardSequences/Digi_cff')
#process.load('Configuration/StandardSequences/SimL1Emulator_cff')
#process.load('Configuration/StandardSequences/DigiToRaw_cff')
#process.load('HLTrigger/Configuration/HLT_1E31_cff')
#process.load('Configuration/StandardSequences/EndOfProcess_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.load('Configuration/EventContent/EventContent_cff')
#
#process.configurationMetadata = cms.untracked.PSet(
#    version = cms.untracked.string('$Revision: 1.123 $'),
#    annotation = cms.untracked.string('step2 nevts:100'),
#    name = cms.untracked.string('PyReleaseValidation')
#)
#process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(10)
#)
#process.options = cms.untracked.PSet(
#    Rethrow = cms.untracked.vstring('ProductNotFound')
#)
## Input source
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:GenHLT_8E29.root')
#)
#
## Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    outputCommands = process.RAWSIMEventContent.outputCommands,
#    fileName = cms.untracked.string('GenHLT_8E29_1E31.root'),
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('GEN-SIM-RAW'),
#        filterName = cms.untracked.string('')
#    )
#)
#
## Additional output definition
#
## Other statements
#process.GlobalTag.globaltag = 'MC_31X_V3::All'
#
## Path and EndPath definitions
#process.digitisation_step = cms.Path(process.pdigi)
#process.L1simulation_step = cms.Path(process.SimL1Emulator)
#process.digi2raw_step = cms.Path(process.DigiToRaw)
#process.endjob_step = cms.Path(process.endOfProcess)
#process.out_step = cms.EndPath(process.output)
#
## Schedule definition
##process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
##process.schedule.extend(process.HLTSchedule)
##process.schedule.extend([process.endjob_step,process.out_step])

