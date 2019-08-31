import FWCore.ParameterSet.Config as cms

######
# Tools to adapt Tau sequences to run tau ReReco+PAT at MiniAOD samples
# M. Bluj, NCBJ Warsaw
# based on work of J. Steggemann, CERN
# Created: 9 Nov. 2017
######

from PhysicsTools.PatAlgos.tools.helpers import *

#####
def addTauReReco(process):
    #PAT
    process.load('PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff')
    process.load('PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi')
    process.selectedPatTaus.cut="pt > 18. && tauID(\'decayModeFindingNewDMs\')> 0.5" #Cut as in MiniAOD
    #Tau RECO
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

#####
def adaptTauToMiniAODReReco(process,
                            reclusterJets=False,
                            jetSrc='slimmedJets',
                            pfCandSrc='packedPFCandidates',
                            vertexSrc='offlineSlimmedPrimaryVertices',
                            postfix=''):
    # TRYING TO MAKE THINGS MINIAOD COMPATIBLE, FROM THE START, TO THE END, 1 BY 1
    #print '[adaptTauToMiniAODReReco]: Start'

    if not isinstance(pfCandSrc,cms.InputTag): pfCandSrc = cms.InputTag(pfCandSrc)

    # Add new jet collections if reclustering is demanded
    if reclusterJets:
        jetPostfix = 'PAT'
        jetSrc = 'patJets'+jetPostfix
        recoJetSrc = 'ak4PFJets'+jetPostfix
        from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
        ak4PFJetsNew = ak4PFJets.clone(
            src=pfCandSrc
        )
        setattr(process,recoJetSrc,ak4PFJetsNew)
        # trivial PATJets
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
 
    # so this adds all tracks to jet in some deltaR region. we don't have tracks so don't need it :D
    # process.ak4PFJetTracksAssociatorAtVertex.jets = cms.InputTag(jetSrc)
    
    # Remove ak4PFJetTracksAssociatorAtVertex from recoTauCommonSequence
    # Remove pfRecoTauTagInfoProducer from recoTauCommonSequence since it uses the jet-track association
    # HOWEVER, may use https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Isolated_Tracks
    # probably needs recovery of the two modules above

    process.recoTauAK4Jets08Region = cms.EDProducer("RecoTauPatJetRegionProducer",
        deltaR = process.recoTauAK4PFJets08Region.deltaR,
        maxJetAbsEta = process.recoTauAK4PFJets08Region.maxJetAbsEta,
        minJetPt = process.recoTauAK4PFJets08Region.minJetPt,
        pfCandAssocMapSrc = cms.InputTag(""),
        pfCandSrc = cms.InputTag('packedPFCandidates'),
        src = cms.InputTag(jetSrc)
    )

    process.recoTauPileUpVertices.src = cms.InputTag("offlineSlimmedPrimaryVertices")
    # Redefine recoTauCommonTask 
    # with redefined region and PU vertices, and w/o track-to-vertex associator and tauTagInfo (the two latter are probably obsolete and not needed at all)
    process.recoTauCommonTask = cms.Task(
        process.recoTauAK4Jets08Region,
        process.recoTauPileUpVertices
    )

    # Adapt TauPiZeros producer
    process.ak4PFJetsLegacyHPSPiZeros.builders[0].qualityCuts.primaryVertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.ak4PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag(jetSrc)

    # Adapt TauChargedHadrons producer
    for builder in process.ak4PFJetsRecoTauChargedHadrons.builders:
        builder.qualityCuts.primaryVertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices")
        if builder.name.value() == 'tracks': #replace plugin based on generalTracks by one based on lostTracks
            builder.name = 'lostTracks'
            builder.plugin = 'PFRecoTauChargedHadronFromLostTrackPlugin'
            builder.srcTracks = cms.InputTag("lostTracks")
    process.ak4PFJetsRecoTauChargedHadrons.jetSrc = cms.InputTag(jetSrc)

    # Adapt combinatoricRecoTau producer
    # this was not in the official pull request, not needed
    #process.combinatoricRecoTaus._TypedParameterizable__type = 'RecoTauPatProducer'
    process.combinatoricRecoTaus.jetRegionSrc = 'recoTauAK4Jets08Region'
    process.combinatoricRecoTaus.jetSrc = jetSrc
    # Adapt builders
    for builder in process.combinatoricRecoTaus.builders:
        for name,value in builder.parameters_().iteritems():
            if name == 'qualityCuts':
                builder.qualityCuts.primaryVertexSrc = 'offlineSlimmedPrimaryVertices'
            elif name == 'pfCandSrc':
                builder.pfCandSrc = pfCandSrc
    # Adapt supported modifiers and remove unsupported ones 
    modifiersToRemove_ = cms.VPSet()
    for mod in process.combinatoricRecoTaus.modifiers:
        if mod.name.value() == 'elec_rej':
            modifiersToRemove_.append(mod)
            continue
        elif mod.name.value() == 'TTIworkaround':
            modifiersToRemove_.append(mod)
            continue
        for name,value in mod.parameters_().iteritems():
            if name == 'qualityCuts':
                mod.qualityCuts.primaryVertexSrc = 'offlineSlimmedPrimaryVertices'
    for mod in modifiersToRemove_:
        process.combinatoricRecoTaus.modifiers.remove(mod)
        #print "\t\t Removing '%s' modifier from 'combinatoricRecoTaus'" %mod.name.value()

    # Redefine tau PV producer
    process.hpsPFTauPrimaryVertexProducer.__dict__['_TypedParameterizable__type'] = 'PFTauMiniAODPrimaryVertexProducer'
    process.hpsPFTauPrimaryVertexProducer.PVTag = 'offlineSlimmedPrimaryVertices'
    process.hpsPFTauPrimaryVertexProducer.packedCandidatesTag = pfCandSrc
    process.hpsPFTauPrimaryVertexProducer.lostCandidatesTag = cms.InputTag("lostTracks")

    # Redefine tau SV producer
    process.hpsPFTauSecondaryVertexProducer = cms.EDProducer("PFTauSecondaryVertexProducer",
        PFTauTag = cms.InputTag("hpsPFTauProducer")
    )
    
    # Instead add against-mu discriminants which are MiniAOD compatible
    from RecoTauTag.RecoTau.hpsPFTauDiscriminationByAMuonRejectionSimple_cff import hpsPFTauDiscriminationByLooseMuonRejectionSimple, hpsPFTauDiscriminationByTightMuonRejectionSimple
    
    process.hpsPFTauDiscriminationByLooseMuonRejectionSimple = hpsPFTauDiscriminationByLooseMuonRejectionSimple
    process.hpsPFTauDiscriminationByTightMuonRejectionSimple = hpsPFTauDiscriminationByTightMuonRejectionSimple

    #####
    # PAT part in the following

    process.tauGenJets.GenParticles = cms.InputTag("prunedGenParticles")
    process.tauMatch.matched = cms.InputTag("prunedGenParticles")

    # Remove unsupported tauIDs
    for name, src in process.patTaus.tauIDSources.parameters_().iteritems():
        if name.find('againstElectron') > -1 or name.find('againstMuon') > -1:
            delattr(process.patTaus.tauIDSources,name)
    # Add MiniAOD specific ones
    setattr(process.patTaus.tauIDSources,'againstMuonLooseSimple',cms.InputTag('hpsPFTauDiscriminationByLooseMuonRejectionSimple'))
    setattr(process.patTaus.tauIDSources,'againstMuonTightSimple',cms.InputTag('hpsPFTauDiscriminationByTightMuonRejectionSimple'))

    #Task/Sequence for tau rereco
    process.miniAODTausTask = cms.Task()
    if reclusterJets:
            process.miniAODTausTask.add(getattr(process,recoJetSrc))
            process.miniAODTausTask.add(getattr(process,jetSrc))
    process.miniAODTausTask.add(process.PFTauTask)
    process.miniAODTausTask.add(process.hpsPFTauDiscriminationByLooseMuonRejectionSimple)
    process.miniAODTausTask.add(process.hpsPFTauDiscriminationByTightMuonRejectionSimple)
    process.miniAODTausTask.add(process.makePatTausTask)
    process.miniAODTausTask.add(process.selectedPatTaus)

    # Remove RecoTau producers which are not supported (yet?), i.e. against-e/mu discriminats
    for moduleName in process.miniAODTausTask.moduleNames(): 
        if 'ElectronRejection' in moduleName or ('MuonRejection' in moduleName and 'Simple' not in moduleName):
            process.miniAODTausTask.remove(getattr(process, moduleName))

    massSearchReplaceAnyInputTag(process.miniAODTausTask,cms.InputTag('ak4PFJets'),cms.InputTag(jetSrc))
    massSearchReplaceAnyInputTag(process.miniAODTausTask,cms.InputTag('particleFlow'),pfCandSrc)
    massSearchReplaceAnyInputTag(process.miniAODTausTask,cms.InputTag('packedPFCandidates'),pfCandSrc)
    massSearchReplaceAnyInputTag(process.miniAODTausTask,cms.InputTag('offlinePrimaryVertices'),cms.InputTag(vertexSrc))
    # currently, the tau reco requires same pf coll for both the jets and the following:
    #process.combinatoricRecoTaus.builders[0].pfCandSrc = cms.InputTag('packedPFCandidates')

    # manually change 1 at a time
    # this one breaks things... but its the one we definitely need!
    #process.recoTauAK4Jets08Region.pfCandSrc = cms.InputTag('packedPFCandidates')



    # now clone the task and add the postfix
    if postfix:
        setattr(process,'miniAODTausTask'+postfix,cloneTask(process,process.miniAODTausTask,postfix))
        setattr(process,'miniAODTausSequence'+postfix,cms.Sequence(getattr(process,'miniAODTausTask'+postfix)))
        #Path with tau rereco (Needed?)
        setattr(process,'TauReco'+postfix,cms.Path(getattr(process,'miniAODTausSequence'+postfix)))
    else:
        process.miniAODTausSequence = cms.Sequence(process.miniAODTausTask)
        #Path with tau rereco (Needed?)
        process.TauReco = cms.Path(process.miniAODTausSequence)
    


    #print '[adaptTauToMiniAODReReco]: Done!'

#####
def setOutputModule(mode=0,postfix=''):
    #mode = 0: store original MiniAOD and new selectedPatTaus 
    #mode = 1: store original MiniAOD, new selectedPatTaus, and all PFTau products as in AOD (except of unsuported ones), plus a few additional collections (charged hadrons, pi zeros, combinatoric reco taus)

    import Configuration.EventContent.EventContent_cff as evtContent
    output = cms.OutputModule(
        'PoolOutputModule',
        fileName=cms.untracked.string('miniAOD_TauReco.root'),
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
    output.outputCommands.append('keep *_selectedPatTaus{}_*_*'.format(postfix))
    if mode==1:
        for prod in evtContent.RecoTauTagAOD.outputCommands:
            if prod.find('ElectronRejection') > -1:
                continue
            if prod.find('MuonRejection') > -1:
                continue
            output.outputCommands.append(prod)
        jetPostfix = 'PAT'+postfix
        output.outputCommands.append('keep *_hpsPFTauDiscriminationByLooseMuonRejectionSimple{}_*_*'.format(postfix))
        output.outputCommands.append('keep *_hpsPFTauDiscriminationByTightMuonRejectionSimple{}_*_*'.format(postfix))
        output.outputCommands.append('keep *_combinatoricReco*_*_*')
        output.outputCommands.append('keep *_ak4PFJetsRecoTauChargedHadrons{}_*_*'.format(postfix))
        output.outputCommands.append('keep *_ak4PFJetsLegacyHPSPiZeros{}_*_*'.format(postfix))
        output.outputCommands.append('keep *_patJets{}_*_*'.format(jetPostfix))

    return output

#####

def cloneTask(process, task, postfix, removePostfix="", noClones = [], addToTask = False):
    """
    ------------------------------------------------------------------
    copy a task plus the modules and tasks therein
    both are renamed by getting a postfix
    input tags are automatically adjusted
    ------------------------------------------------------------------
    """
    result = task
    if not postfix == "":
        visitor = CloneTaskVisitor(process, task.label(), postfix, removePostfix, noClones, addToTask)
        task.visit(visitor)
        result = visitor.clonedTask()
    return result

class CloneTaskVisitor(object):
    """Visitor that travels within a cms.Task, and returns a cloned version of the Task.
    All modules and tasks are cloned and a postfix is added"""
    def __init__(self, process, label, postfix, removePostfix="", noClones = [], addToTask = False):
        self._process = process
        self._postfix = postfix
        self._removePostfix = removePostfix
        self._noClones = noClones
        self._addToTask = addToTask
        self._moduleLabels = []
        self._clonedTask = cms.Task()
        setattr(process, self._newLabel(label), self._clonedTask)
        if addToTask:
            self._patAlgosToolsTask = getPatAlgosToolsTask(process)

    def enter(self, visitee):
        if isinstance(visitee, cms._Module):
            label = visitee.label()
            newModule = None
            if label in self._noClones: #keep unchanged
                newModule = getattr(self._process, label)
            elif label in self._moduleLabels: # has the module already been cloned ?
                newModule = getattr(self._process, self._newLabel(label))
            else:
                self._moduleLabels.append(label)
                newModule = visitee.clone()
                setattr(self._process, self._newLabel(label), newModule)
                if self._addToTask:
                    self._patAlgosToolsTask.add(getattr(self._process, self._newLabel(label)))
            self.__addToTopTask(newModule)

    def leave(self, visitee):
        pass

    def clonedTask(self):
        for label in self._moduleLabels:
            massSearchReplaceAnyInputTag(self._clonedTask, label, self._newLabel(label), moduleLabelOnly=True, verbose=False)
        self._moduleLabels = [] # prevent the InputTag replacement next time the 'clonedTask' function is called.
        return self._clonedTask

    def _newLabel(self, label):
        if self._removePostfix != "":
            if label[-len(self._removePostfix):] == self._removePostfix:
                label = label[0:-len(self._removePostfix)]
            else:
                raise Exception("Tried to remove postfix %s from label %s, but it wasn't there" % (self._removePostfix, label))
        return label + self._postfix

    def __addToTopTask(self, visitee):
        self._clonedTask.add(visitee)
