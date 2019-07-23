// -*- C++ -*-
//
// Package:    MuonCleanedJetProducer
// Class:      MuonCleanedJetProducer
// 
/**\class MuonCleanedJetProducer MuonCleanedJetProducer.cc

 Description: Removes PF muons from PFJet candidates and reconstructs the jets
	          Associates those muons to the jets from which they were removed

 Implementation:
     [Notes on implementation]
*/
//
//          Author:  Devin Taylor
//
//


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "TLorentzVector.h"
#include "TMath.h"

//
// class declaration
//

class MuonCleanedJetProducer : public edm::stream::EDProducer<>
{
   public:
      explicit MuonCleanedJetProducer(const edm::ParameterSet&);
      ~MuonCleanedJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      typedef edm::AssociationMap<edm::OneToMany<std::vector<reco::PFJet>, std::vector<reco::PFCandidate>, unsigned int> >
            JetToPFCandidateAssociation;
      
      // ----------member data ---------------------------

      // source of the jets to be cleaned of muons
      edm::EDGetTokenT<reco::PFJetCollection> jetSrc_;

      // source of muons that, if found within jet, should be removed
      edm::EDGetTokenT<reco::MuonRefVector> muonSrc_;

      // source of PF candidates
      edm::EDGetTokenT<reco::PFCandidateCollection> pfCandSrc_;

      edm::ParameterSet* cfg_;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MuonCleanedJetProducer::MuonCleanedJetProducer(const edm::ParameterSet& iConfig):
  jetSrc_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  muonSrc_(consumes<reco::MuonRefVector>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  pfCandSrc_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSrc")))
{
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);

  //register your products
  produces<reco::PFJetCollection>();
  produces<edm::RefVector<reco::PFCandidateCollection> >("particleFlowMuonCleaned");
  produces<JetToPFCandidateAssociation>("pfCandAssocMapForIsolation");
}


MuonCleanedJetProducer::~MuonCleanedJetProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonCleanedJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByToken(jetSrc_, pfJets);
  std::unique_ptr<reco::PFJetCollection> cleanedJets = std::make_unique<reco::PFJetCollection>();
  edm::RefProd<reco::PFJetCollection> selectedJetRefProd = iEvent.getRefBeforePut<reco::PFJetCollection>();

  auto selectedJetPFCandidateAssociationForIsolation =
        std::make_unique<JetToPFCandidateAssociation>(&iEvent.productGetter());

  edm::Handle<reco::MuonRefVector> muons;
  iEvent.getByToken(muonSrc_, muons);

  edm::Handle<reco::PFCandidateCollection> pfCands;
  iEvent.getByToken(pfCandSrc_, pfCands);
  std::unique_ptr<edm::RefVector<reco::PFCandidateCollection> > pfCandsExcludingMuons(new edm::RefVector<reco::PFCandidateCollection>);

  //fill an STL container with muon ref keys
  std::vector<unsigned int> muonRefKeys;
  if (muons.isValid()) 
  {
    for (reco::MuonRefVector::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon)
    {
      muonRefKeys.push_back(iMuon->key());
    }
  }

  // Do cleaning
  for (reco::PFJetCollection::const_iterator iJet = pfJets->begin(); iJet != pfJets->end(); ++iJet)
  {
    std::vector<reco::PFCandidatePtr> jetPFCands = iJet->getPFConstituents();
    reco::PFJet::Specific specs = iJet->getSpecific();
    math::XYZTLorentzVector pfmomentum;
    std::vector<edm::Ptr<reco::Candidate> > jetConstituents;
    jetConstituents.clear();

    for (std::vector<edm::Ptr<reco::PFCandidate> >::iterator i = jetPFCands.begin(); i != jetPFCands.end(); ++i)
    {
      reco::PFCandidate pfCand = *i;
      
      // Is the PF Candidate a muon?
      if (pfCand.particleId() == 3) //Reference: https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_1_17/doc/html/d8/d17/PFCandidate_8h_source.html
      {
        reco::MuonRef theRecoMuon = pfCand.muonRef();

        //does this muon pass the desired muon ID?
        std::vector<unsigned int>::const_iterator iMuon = std::find(muonRefKeys.begin(), muonRefKeys.end(), theRecoMuon.key());
   
        if (iMuon != muonRefKeys.end()) 
   	    {
          specs.mMuonEnergy -= pfCand.p4().e();
          specs.mMuonMultiplicity -= 1;
          specs.mChargedMuEnergy -= pfCand.p4().e();
          specs.mChargedMultiplicity -= 1;
        }
        else
        {
          pfmomentum += pfCand.p4(); // total p4()
          jetConstituents.push_back((*i));
        }
      }
      else // if it's not a muon
      {
        pfmomentum += pfCand.p4(); // total p4()
        jetConstituents.push_back((*i));
      }
    } // loop over PF candidates

    // Build a new jet without the muon
    reco::PFJet muonfreePFJet(pfmomentum, specs, jetConstituents);
    cleanedJets->push_back( muonfreePFJet );

    // get ref to output jet for isolation
    edm::Ref<reco::PFJetCollection> jetRef(selectedJetRefProd, cleanedJets->size() - 1);

    // get muons not in jet and add to isolation association
    for (size_t i = 0; i < pfCands->size(); ++i) {
      if ((*pfCands)[i].particleId() == 3) {
        reco::MuonRef theRecoMuon = (*pfCands)[i].muonRef();
        std::vector<unsigned int>::const_iterator iMuon = std::find(muonRefKeys.begin(), muonRefKeys.end(), theRecoMuon.key());
        if (iMuon != muonRefKeys.end()) {
          reco::PFCandidateRef pfCandRef(pfCands,i);
          selectedJetPFCandidateAssociationForIsolation->insert(jetRef, pfCandRef);
        }
      }
    }
  } // loop over jets
  
  // build a collection of PF candidates excluding muons
  // we will still tag the jet as signal-like by the presence of a muon IN the jet, but this 
  // ensures that such jets also cannot have the muon enter the isolation candidate collection
  for (unsigned int i = 0; i < pfCands->size(); ++i)
  {
    reco::MuonRef removedMuRef = (*pfCands)[i].muonRef();
    if ((removedMuRef.isNonnull() && (std::find(muonRefKeys.begin(), muonRefKeys.end(), removedMuRef.key()) == muonRefKeys.end())) || removedMuRef.isNull()) 
    {
      pfCandsExcludingMuons->push_back(edm::RefVector<reco::PFCandidateCollection>::value_type(pfCands, i));
    }
  }

  iEvent.put(std::move(cleanedJets));

  //put the soft-muon-free PF cands into the event
  iEvent.put(std::move(pfCandsExcludingMuons), "particleFlowMuonCleaned");
  iEvent.put(std::move(selectedJetPFCandidateAssociationForIsolation), "pfCandAssocMapForIsolation");

}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonCleanedJetProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonCleanedJetProducer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonCleanedJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonCleanedJetProducer);
