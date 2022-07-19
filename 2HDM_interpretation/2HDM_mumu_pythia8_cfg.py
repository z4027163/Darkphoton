###########################################################################################################################################################################
### cmsDriver command
### cmsDriver.py Configuration/GenProduction/python/HIG-RunIISummer20UL18wmLHEGEN-00789-fragment.py --python_filename HIG-RunIISummer20UL18wmLHEGEN-00789_1_cfg.py --eventcontent RAWSIM,LHE --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN,LHE --fileout file:HIG-RunIISummer20UL18wmLHEGEN-00789.root --conditions 106X_upgrade2018_realistic_v4 --beamspot Realistic25ns13TeVEarly2018Collision --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(200)" --step LHE,GEN --geometry DB:Extended --era Run2_2018 --no_exec --mc -n $EVENTS || exit $? ;
###########################################################################################################################################################################

mass = 1.6 #sets dark photon mass

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process('SIM', eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2018Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(10000)
)

from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')
params.register('seed',1234567,VarParsing.multiplicity.singleton,VarParsing.varType.int,'the random seed')
params.register('suffix','theOutName',VarParsing.multiplicity.singleton,VarParsing.varType.string, 'the file qill contain this substring')
params.parseArguments()

process.RandomNumberGeneratorService.generator.initialSeed = params.seed

# Input source
process.source = cms.Source("EmptySource")
process.options = cms.untracked.PSet()

# Output definition
process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    
    fileName = cms.untracked.string('file:lightPseudoscalarH_mumu_filter_10k.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v4', '')

# Pythia settings
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
                                 pythiaPylistVerbosity = cms.untracked.int32(0),
                                 filterEfficiency = cms.untracked.double(1.0),
                                 pythiaHepMCVerbosity = cms.untracked.bool(False),
                                 comEnergy = cms.double(13000.0),
                                 crossSection = cms.untracked.double(0),
                                 maxEventsToPrint = cms.untracked.int32(0),
                                 PythiaParameters = cms.PSet(
                                     pythia8CommonSettings = cms.vstring(#'Tune:preferLHAPDF = 2',
                                                                         'Main:timesAllowErrors = 10000', 
                                                                         'Check:epTolErr = 0.01', 
                                                                         'Beams:setProductionScalesFromLHEF = off', 
                                                                         'SLHA:keepSM = on', 
                                                                         'ParticleDecays:limitTau0 = off', 
                                                                         'ParticleDecays:allowPhotonRadiation = on',
                                                                         'ParticleDecays:limitTau0 = off',
                                                                         #'PhaseSpace:mHatMin = 0.05',
                                                                         #'PhaseSpace:mHatMax = 20',
                                                                     ),
                                     pythia8CP5Settings = cms.vstring(#'Tune:pp 14', 
                                                                      #'Tune:ee 7', 
                                                                  ),
                                     processParameters = cms.vstring('Higgs:useBSM = on',
                                                                     #'HiggsBSM:gg2A3 = on', 
                                                                     'HiggsBSM:qg2A3q = on',
                                                                     'HiggsBSM:gg2A3bbbar = on',
                                                                     'HiggsBSM:qqbar2A3bbbar = on',
                                                                     'HiggsBSM:gg2A3g(l:t) = on',
                                                                     'HiggsBSM:qg2A3q(l:t) = on',
                                                                     'HiggsBSM:qqbar2A3g(l:t) = on',
                                                                     '36:onMode = off',
                                                                     '36:onIfAny = 13',
                                                                     '36:mMin = '+str(mass-0.5),
                                                                     '36:mMax = '+str(mass+0.5),
                                                                     'Higgs:clipWings = off',
                                                                     '36:m0 = '+str(mass),
                                                                     #'36:mWidth = 0.01',
                                                                 ),
                                     parameterSets = cms.vstring(
                                         'pythia8CommonSettings', 
                                         'pythia8CP5Settings', 
                                         'processParameters')
                                 )
                             )

#process.ProductionFilterSequence = cms.Sequence(process.generator)

process.mumugenfilter = cms.EDFilter("MCParticlePairFilter",
                                     Status = cms.untracked.vint32(1, 1),
                                     MinPt = cms.untracked.vdouble(4, 4),
                                     MaxEta = cms.untracked.vdouble(1.9, 1.9),
                                     MinEta = cms.untracked.vdouble(-1.9, -1.9),
                                     #MinInvMass = cms.untracked.double(2.0),
                                     #MaxInvMass = cms.untracked.double(4.0),
                                     ParticleCharge = cms.untracked.int32(-1),
                                     ParticleID1 = cms.untracked.vint32(13),
                                     ParticleID2 = cms.untracked.vint32(13)
                                 )

process.ProductionFilterSequence = cms.Sequence(process.generator * process.mumugenfilter)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

# End of customisation functions
