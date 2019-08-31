# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --filein dbs:/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM --fileout file:mini.root --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions 94X_mc2017_realistic_v14 --step PAT --nThreads 4 --scenario pp --era Run2_2017,run2_miniAOD_94XFall17 --python_filename mini_1_cfg.py --no_exec -n 100
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2017_cff import Run2_2017
from Configuration.Eras.Modifier_run2_miniAOD_94XFall17_cff import run2_miniAOD_94XFall17

process = cms.Process('PAT',Run2_2017,run2_miniAOD_94XFall17)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/40000/02CF836E-AB60-E811-A5A7-1866DA85DC7F.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/40000/4E7AD41E-6761-E811-A91B-008CFA1CB470.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/40000/6A70ACD7-8162-E811-B933-008CFAF292AC.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/60000/805EACDA-A064-E811-AECD-44A842B2D631.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/60000/B4153EF4-A064-E811-A2E9-0025901DF5EC.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/6AE398D5-CC6E-E811-A235-3417EBE465FE.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/608599C5-CB6E-E811-B254-10983627C3CE.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/F684C323-F46E-E811-A37A-7CD30AC03006.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/18A216BB-1F6F-E811-9D59-7CD30AC0311A.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/C0FC2675-0C6F-E811-9C9B-7CD30ABD2EE8.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/C0B58583-186F-E811-A53E-3417EBE465FE.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/F0EB2853-3F6F-E811-9B6C-7CD30ABD2EE8.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/B08AEFA6-4D6F-E811-B996-3417EBE465FE.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/D634ADD9-6E6F-E811-8D7F-000E1EB004E0.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/B8E50AC7-6D6F-E811-B21F-F04DA27541B7.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/E6616CF2-F06E-E811-95C2-008CFAF554D0.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/0A93D330-196F-E811-861E-3417EBE5264B.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/603ACD1A-426F-E811-9082-3417EBE5264B.root', 
        '/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/70000/F81E7AA5-476F-E811-A38C-A0369FC51EE4.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:mini.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14', '')

# Path and EndPath definitions
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,process.MINIAODSIMoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(1)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

# cusomization for muon cleaned taus
import PhysicsTools.PatAlgos.tools.helpers as configtools
from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

# lower the pt threshold
def lowerTauPt(process,postfix='',tauPt=8, jetPt=5):
    from FWCore.ParameterSet.MassReplace import massSearchReplaceParam
    massSearchReplaceParam(getattr(process,'PATTauSequence'+postfix),'minJetPt',14,jetPt)
    getattr(process,'selectedPatTaus'+postfix).cut = cms.string("pt > {} && tauID(\'decayModeFindingNewDMs\')> 0.5".format(tauPt))

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
patAlgosToolsTask = configtools.getPatAlgosToolsTask(process)

process.PATTauSequence = cms.Sequence(process.PFTau+process.makePatTaus+process.selectedPatTaus)

#########################
### Muon Cleaned Taus ###
#########################
process.recoMuonsForJetCleaning = cms.EDFilter('MuonRefSelector',
    src = cms.InputTag('muons'),
    cut = cms.string('pt > 3.0 && isPFMuon && (isGlobalMuon || isTrackerMuon)'),
)

process.ak4PFJetsMuonCleaned = cms.EDProducer(
    'MuonCleanedJetProducer',
    jetSrc = cms.InputTag("ak4PFJets"),
    muonSrc = cms.InputTag("recoMuonsForJetCleaning"),
    pfCandSrc = cms.InputTag("particleFlow"),
)

process.muonCleanedHPSPFTausTask = cms.Task(
    process.recoMuonsForJetCleaning,
    process.ak4PFJetsMuonCleaned
)
patAlgosToolsTask.add(process.muonCleanedHPSPFTausTask)

jetSrc = cms.InputTag('ak4PFJetsMuonCleaned')
#pfCandSrc = cms.InputTag('ak4PFJetsMuonCleaned','particleFlowMuonCleaned')
#pfCandSrc = cms.InputTag('particleFlow')
pfAssocMap = cms.InputTag('ak4PFJetsMuonCleaned','pfCandAssocMapForIsolation')

process.PATTauSequenceMuonCleaned = cloneProcessingSnippet(process,process.PATTauSequence, 'MuonCleaned', addToTask = True)
massSearchReplaceAnyInputTag(process.PATTauSequenceMuonCleaned,cms.InputTag('ak4PFJets'),jetSrc)
#massSearchReplaceAnyInputTag(process.PATTauSequenceMuonCleaned,cms.InputTag('particleFlow'),pfCandSrc)
#process.combinatoricRecoTausMuonCleaned.builders[0].pfCandSrc = pfCandSrc
# alternative approach
process.recoTauAK4PFJets08RegionMuonCleaned.pfCandAssocMapSrc = pfAssocMap

process.slimmedTausMuonCleaned = process.slimmedTausNoDeepIDs.clone(
    src = cms.InputTag('selectedPatTausMuonCleaned'),
)
patAlgosToolsTask.add(process.slimmedTausMuonCleaned)

lowerTauPt(process,'MuonCleaned')

# MC customization
process.tauGenJetsMuonCleaned.GenParticles = "prunedGenParticles"
process.patTausMuonCleaned.embedGenMatch = False

# additional changes to standard MiniAOD content
process.MINIAODSIMoutput.outputCommands += [
    'keep *_slimmedTausMuonCleaned_*_*',
]
