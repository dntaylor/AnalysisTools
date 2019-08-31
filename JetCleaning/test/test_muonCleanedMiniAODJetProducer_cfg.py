import FWCore.ParameterSet.Config as cms

maxEvents = 10000
doLower = True
doMuTau = True
doETau = True
doTauTau = False
nCore = 8

from Configuration.StandardSequences.Eras import eras
era = eras.Run2_2018
process = cms.Process("BoostedTauReco", era)
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource", 
    fileNames=readFiles, 
    secondaryFileNames=secFiles,
)

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(maxEvents)
)

readFiles.extend([
    'file:mini.root',
    #'/store/mc/RunIIFall17MiniAODv2/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/40000/8683D81E-BC62-E811-8C66-D4AE52901DF0.root',
])

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')


import Configuration.EventContent.EventContent_cff as evtContent
process.output = cms.OutputModule(
    'PoolOutputModule',
    fileName=cms.untracked.string('miniAOD_BoostedDiTauReco.root'),
    fastCloning=cms.untracked.bool(False),
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string('MINIAODSIM'),
        filterName=cms.untracked.string('')
    ),
    outputCommands = evtContent.MINIAODSIMEventContent.outputCommands,
    SelectEvents=cms.untracked.PSet(
        SelectEvents=cms.vstring('*',)
        )
    )

process.out = cms.EndPath(process.output)

######################
### rerun tau reco ###
######################
#import RecoTauTag.Configuration.tools.adaptToRunAtMiniAOD as tauAtMiniTools
import AnalysisTools.JetCleaning.adaptToRunAtMiniAOD as tauAtMiniTools
tauAtMiniTools.addTauReReco(process)

# lower the pt threshold
def lowerTauPt(process,postfix='',tauPt=8, jetPt=5):
    from FWCore.ParameterSet.MassReplace import massSearchReplaceParam
    massSearchReplaceParam(getattr(process,'miniAODTausTask'+postfix),'minJetPt',14,jetPt)
    getattr(process,'selectedPatTaus'+postfix).cut = cms.string("pt > {} && tauID(\'decayModeFindingNewDMs\')> 0.5".format(tauPt))

#########################
### muon cleaned taus ###
#########################


# alternative approach, full recluster jets without muons
jetSrc = 'patJetsMuonCleaned'
recoJetSrc = 'ak4PFJetsMuonCleaned'
pfCandSrc = "muonCleanedPackedPFCandidates"

process.recoMuonsForJetCleaning = cms.EDFilter('PATMuonRefSelector',
    src = cms.InputTag('slimmedMuons'),
    cut = cms.string('pt > 3.0 && isPFMuon && (isGlobalMuon || isTrackerMuon)'),
)

process.muonCleanedPackedPFCandidates = cms.EDProducer("MuonCleanedPackedCandidateProducer",
    src = cms.InputTag("recoMuonsForJetCleaning"),
    pfCandSrc = cms.InputTag("packedPFCandidates"),
)

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4PFJetsMuonCleaned = ak4PFJets.clone(
    src=cms.InputTag('muonCleanedPackedPFCandidates')
)

from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
patJetsNew = _patJets.clone(
    jetSource            = cms.InputTag(recoJetSrc),
    addJetCorrFactors    = cms.bool(False),
    jetCorrFactorsSource = cms.VInputTag(),
    addBTagInfo          = cms.bool(False),
    addDiscriminators    = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    addAssociatedTracks  = cms.bool(False),
    addJetCharge         = cms.bool(False),
    addGenPartonMatch    = cms.bool(False),
    embedGenPartonMatch  = cms.bool(False),
    addGenJetMatch       = cms.bool(False),
    getJetMCFlavour      = cms.bool(False),
    addJetFlavourInfo    = cms.bool(False),
)
setattr(process,jetSrc,patJetsNew)


process.muonCleanedHPSPFTausTask = cms.Task(
    process.recoMuonsForJetCleaning,
    process.muonCleanedPackedPFCandidates,
    process.ak4PFJetsMuonCleaned,
    process.patJetsMuonCleaned,
)

tauAtMiniTools.adaptTauToMiniAODReReco(process, jetSrc=jetSrc, pfCandSrc=pfCandSrc, postfix='MuonCleaned')
process.miniAODTausTaskMuonCleaned.add(process.muonCleanedHPSPFTausTask)
lowerTauPt(process,postfix='MuonCleaned')
process.output.outputCommands.append('keep *_selectedPatTausMuonCleaned_*_*')



###############
### options ###
###############
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
)
process.options.numberOfThreads = cms.untracked.uint32(nCore)
process.options.numberOfStreams = cms.untracked.uint32(0)

process.options = cms.untracked.PSet(
    process.options,
    wantSummary=cms.untracked.bool(False)
)
