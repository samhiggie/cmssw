// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/RecoMuonFilter
// Class:      RecoMuonFilter
// 
/**\class RecoMuonFilter RecoMuonFilter.cc PhaseTwoAnalysis/RecoMuonFilter/plugins/RecoMuonFilter.cc

Description: adds a vector of reco muons

Implementation:
- muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/UPGTrackerTDRStudies#Muon_identification
*/
//
// Original Author:  Elvire Bouvier
//         Created:  Sun, 02 Jul 2017 20:54:15 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include "Math/GenVector/VectorUtil.h"

//
// class declaration
//

class RecoMuonFilter : public edm::stream::EDProducer<> {
  public:
    explicit RecoMuonFilter(const edm::ParameterSet&);
    ~RecoMuonFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    bool isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
    edm::EDGetTokenT<std::vector<reco::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandsNoLepToken_;
    
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
RecoMuonFilter::RecoMuonFilter(const edm::ParameterSet& iConfig):
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  pfCandsNoLepToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandsNoLep")))    
{
  produces<std::vector<reco::Muon>>("LooseMuons");
  produces<std::vector<double>>("LooseMuonRelIso");
  produces<std::vector<reco::Muon>>("MediumMuons");
  produces<std::vector<double>>("MediumMuonRelIso");
  produces<std::vector<reco::Muon>>("TightMuons");
  produces<std::vector<double>>("TightMuonRelIso");

}


RecoMuonFilter::~RecoMuonFilter()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
  void
RecoMuonFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);
  // Vertices
  int prVtx = -1;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
  }

  Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  Handle<std::vector<reco::PFCandidate>> pfCandsNoLep;
  iEvent.getByToken(pfCandsNoLepToken_, pfCandsNoLep);  
  std::unique_ptr<std::vector<reco::Muon>> filteredLooseMuons;
  std::unique_ptr<std::vector<double>> filteredLooseMuonRelIso;
  std::unique_ptr<std::vector<reco::Muon>> filteredMediumMuons;
  std::unique_ptr<std::vector<double>> filteredMediumMuonRelIso;
  std::unique_ptr<std::vector<reco::Muon>> filteredTightMuons;
  std::unique_ptr<std::vector<double>> filteredTightMuonRelIso;
  std::vector<reco::Muon> looseVec;
  std::vector<double> looseIsoVec;
  std::vector<reco::Muon> mediumVec;
  std::vector<double> mediumIsoVec;
  std::vector<reco::Muon> tightVec;
  std::vector<double> tightIsoVec;
  for (size_t i = 0; i < muons->size(); i++) {
    if (muons->at(i).pt() < 10.) continue;
    if (fabs(muons->at(i).eta()) > 3.) continue;

    bool isLoose = (fabs(muons->at(i).eta()) < 2.4 && muon::isLooseMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.5));
    bool isMedium = (fabs(muons->at(i).eta()) < 2.4 && muon::isMediumMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.3));
    bool isTight = (fabs(muons->at(i).eta()) < 2.4 && prVtx > -0.5 && muon::isTightMuon(muons->at(i),vertices->at(prVtx))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.1));
    double relIso = 0.;
    for (size_t k = 0; k < pfCandsNoLep->size(); k++) {
      if (ROOT::Math::VectorUtil::DeltaR(muons->at(i).p4(),pfCandsNoLep->at(k).p4()) > 0.3) continue;
      relIso += pfCandsNoLep->at(k).pt();
    }
    if (muons->at(i).pt() > 0.) relIso = relIso / muons->at(i).pt();
    else relIso = -1.;

    if (!isLoose) continue;
    looseVec.push_back(muons->at(i));
    looseIsoVec.push_back(relIso);

    if (!isMedium) continue;
    mediumVec.push_back(muons->at(i));
    mediumIsoVec.push_back(relIso);

    if (!isTight) continue;
    tightVec.push_back(muons->at(i));
    tightIsoVec.push_back(relIso);

  }

  filteredLooseMuons.reset(new std::vector<reco::Muon>(looseVec));
  filteredLooseMuonRelIso.reset(new std::vector<double>(looseIsoVec));
  filteredMediumMuons.reset(new std::vector<reco::Muon>(mediumVec));
  filteredMediumMuonRelIso.reset(new std::vector<double>(mediumIsoVec));
  filteredTightMuons.reset(new std::vector<reco::Muon>(tightVec));
  filteredTightMuonRelIso.reset(new std::vector<double>(tightIsoVec));

  iEvent.put(std::move(filteredLooseMuons), "LooseMuons");
  iEvent.put(std::move(filteredLooseMuonRelIso), "LooseMuonRelIso");
  iEvent.put(std::move(filteredMediumMuons), "MediumMuons");
  iEvent.put(std::move(filteredMediumMuonRelIso), "MediumMuonRelIso");
  iEvent.put(std::move(filteredTightMuons), "TightMuons");
  iEvent.put(std::move(filteredTightMuonRelIso), "TightMuonRelIso");

  return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
  void
RecoMuonFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
RecoMuonFilter::endStream() {
}

  bool 
RecoMuonFilter::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaX = 999;
    double deltaY = 999;
    double pullX = 999;
    double pullY = 999;
    double deltaPhi = 999;

    bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for(std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber){

      for (std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment){

        if (chamber->detector() == 5){

          deltaX   = fabs(chamber->x - segment->x);
          deltaY   = fabs(chamber->y - segment->y);
          pullX    = fabs(chamber->x - segment->x) / std::sqrt(chamber->xErr + segment->xErr);
          pullY    = fabs(chamber->y - segment->y) / std::sqrt(chamber->yErr + segment->yErr);
          deltaPhi = fabs(atan(chamber->dXdZ) - atan(segment->dXdZ));

        }
      }
    }

    if ((pullX < pullXCut) || (deltaX < dXCut)) X_MatchFound = true;
    if ((pullY < pullYCut) || (deltaY < dYCut)) Y_MatchFound = true;
    if (deltaPhi < dPhi) Dir_MatchFound = true;

    result = X_MatchFound && Y_MatchFound && Dir_MatchFound;

  }

  return result;

}

// ------------ method called when starting to processes a run  ------------
/*
   void
   RecoMuonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   RecoMuonFilter::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   RecoMuonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   RecoMuonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecoMuonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoMuonFilter);
