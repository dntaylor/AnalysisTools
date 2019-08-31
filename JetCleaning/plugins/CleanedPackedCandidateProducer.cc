// -*- C++ -*-
//
// Package:    CleanedPackedCandidateProducer
// Class:      CleanedPackedCandidateProducer
// 
/**\class CleanedPackedCandidateProducer CleanedPackedCandidateProducer.cc

 Description: Removes PF cands from PackedCandidates

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Devin Taylor
//         Created:  26-03-19
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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"
#include "TMath.h"

//
// class declaration
//
template<class PatType>
class CleanedPackedCandidateProducer : public edm::stream::EDProducer<>
{
   public:
      explicit CleanedPackedCandidateProducer(const edm::ParameterSet&);
      ~CleanedPackedCandidateProducer() override {};

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      
      // source of objs to remove
      edm::EDGetTokenT<edm::RefVector<std::vector<PatType> > > src_;

      // source of PF candidates
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandSrc_;

      edm::ParameterSet* cfg_;

};

//
// constructors and destructor
//
template<class PatType>
CleanedPackedCandidateProducer<PatType>::CleanedPackedCandidateProducer(const edm::ParameterSet& iConfig):
  src_(consumes<edm::RefVector<std::vector<PatType> > >(iConfig.getParameter<edm::InputTag>("src"))),
  pfCandSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSrc")))
{
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);

  //register your products
  // Produce refs, doesnt work
  //produces<edm::RefVector<pat::PackedCandidateCollection> >();
  // produces objs, only works during runtime
  produces<pat::PackedCandidateCollection>();

}


//
// member functions
//
template<class PatType>
void
CleanedPackedCandidateProducer<PatType>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::RefVector<std::vector<PatType> > > objs;
  iEvent.getByToken(src_, objs);

  edm::Handle<pat::PackedCandidateCollection> pfCands;
  iEvent.getByToken(pfCandSrc_, pfCands);
  //std::unique_ptr<edm::RefVector<pat::PackedCandidateCollection> > pfCandsExcluding(new edm::RefVector<pat::PackedCandidateCollection>);
  std::unique_ptr<pat::PackedCandidateCollection> pfCandsExcluding(new pat::PackedCandidateCollection);

  std::vector<reco::CandidatePtr> PFs;
  if (objs.isValid()) 
  {
    for (typename edm::RefVector<std::vector<PatType> >::const_iterator iObj = objs->begin(); iObj != objs->end(); ++iObj)
    {
      for (unsigned int j = 0; j<(*iObj)->numberOfSourceCandidatePtrs(); j++) {
        PFs.push_back((*iObj)->sourceCandidatePtr(j));
      }
    }
  }

  // fill the vector of particleFlow without objs
  for (unsigned int i = 0; i != pfCands->size(); i++) {
    if (std::find(PFs.begin(), PFs.end(), reco::CandidatePtr(pfCands,i)) == PFs.end()) {
      //pfCandsExcluding->push_back(edm::RefVector<pat::PackedCandidateCollection>::value_type(pfCands, i));
      const pat::PackedCandidate &pf = (*pfCands)[i];
      pfCandsExcluding->push_back(pf);
    }
  }

  //put the cleaned PF cands into the event
  iEvent.put(std::move(pfCandsExcluding));

}

//define this as a plug-in
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
typedef CleanedPackedCandidateProducer<pat::Muon> MuonCleanedPackedCandidateProducer;
typedef CleanedPackedCandidateProducer<pat::Electron> ElectronCleanedPackedCandidateProducer;
typedef CleanedPackedCandidateProducer<pat::Photon> PhotonCleanedPackedCandidateProducer;

DEFINE_FWK_MODULE(MuonCleanedPackedCandidateProducer);
DEFINE_FWK_MODULE(ElectronCleanedPackedCandidateProducer);
DEFINE_FWK_MODULE(PhotonCleanedPackedCandidateProducer);
