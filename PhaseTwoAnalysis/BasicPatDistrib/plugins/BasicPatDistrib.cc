// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/BasicPatDistrib
// Class:      BasicPatDistrib
// 
/**\class BasicPatDistrib BasicPatDistrib.cc PhaseTwoAnalysis/BasicPatDistrib/plugins/BasicPatDistrib.cc

Description: produces histograms of basic quantities from PAT collections

Implementation:
   - lepton isolation might need to be refined
   - muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/UPGTrackerTDRStudies#Muon_identification
   - electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
      /!\ no ID is implemented for forward electrons as:
      - PFClusterProducer does not run on miniAOD
      - jurassic isolation needs tracks
   - PF jet ID comes from Run-2 https://github.com/cms-sw/cmssw/blob/CMSSW_9_1_1_patch1/PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
   - no JEC applied
   - b-tagging WP come from Run-2 https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco#Supported_Algorithms_and_Operati 
      - for pfCombinedInclusiveSecondaryVertexV2BJetTags: L = 0.5426, M = 0.8484, T = 0.9535)
      - for deepCSV: L = 0.2219, M = 0.6324, T = 0.8958
*/
/*
    Modifier: Sam Higginbotham 
    Date: July 6th,2017. 
    - Photon ID:Selection based on https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details 
    - Starting template taken DIRECTLY from RecoEgamma/Examples/plugins/PatPhotonSimpleAnalyzer.cc  (Author:  M.B. Anderson based on simple photon analyzer by:  J. Stilley, A. Askew)


*/
//
// Original Author:  Elvire Bouvier
//         Created:  Wed, 14 Jun 2017 14:16:22 GMT
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"//
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"

//for photons
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "TMath.h"
#include "TTree.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class BasicPatDistrib : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit BasicPatDistrib(const edm::ParameterSet&);
    ~BasicPatDistrib();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0,
      TRUE_PROMPT_ELECTRON,
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON};  


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    bool isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isME0MuonSel(reco::Muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);

    // ----------member data ---------------------------
    edm::Service<TFileService> fs_;

    bool useDeepCSV_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
    edm::EDGetTokenT<std::vector<pat::Electron>> elecsToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
    PFJetIDSelectionFunctor jetIDLoose_;
    PFJetIDSelectionFunctor jetIDTight_;
    edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
    edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    edm::EDGetTokenT<std::vector<pat::Photon>> photonsToken_;

    // MC truth in fiducial phase space
    TH1D* h_genMuons_n_;
    TH1D* h_genMuons_pt_;
    TH1D* h_genMuons_phi_;
    TH1D* h_genMuons_eta_;
    TH1D* h_genMuons_iso_;  
    TH1D* h_genElecs_n_;
    TH1D* h_genElecs_pt_;
    TH1D* h_genElecs_phi_;
    TH1D* h_genElecs_eta_;
    TH1D* h_genElecs_iso_; 
    TH1D* h_genJets_n_;
    TH1D* h_genJets_pt_;
    TH1D* h_genJets_phi_;
    TH1D* h_genJets_eta_;

    // Vertices
    TH1D* h_allVertices_n_;
    // ... that pass ID
    TH1D* h_goodVertices_n_;

    // Muons
    TH1D* h_allMuons_n_;
    TH1D* h_allMuons_pt_;
    TH1D* h_allMuons_phi_;
    TH1D* h_allMuons_eta_;
    TH1D* h_allMuons_iso_;
    TH1D* h_allMuons_id_;
    // ... that pass kin cuts, tight ID, and are isolated
    TH1D* h_goodMuons_n_;
    TH1D* h_goodMuons_pt_;
    TH1D* h_goodMuons_phi_;
    TH1D* h_goodMuons_eta_;
    TH1D* h_goodMuons_iso_;

    // Elecs
    TH1D* h_allElecs_n_;
    TH1D* h_allElecs_pt_;
    TH1D* h_allElecs_phi_;
    TH1D* h_allElecs_eta_;
    TH1D* h_allElecs_iso_;
    TH1D* h_allElecs_id_;
    // ... that pass kin cuts, tight ID, and are isolated
    TH1D* h_goodElecs_n_;
    TH1D* h_goodElecs_pt_;
    TH1D* h_goodElecs_phi_;
    TH1D* h_goodElecs_eta_;
    TH1D* h_goodElecs_iso_;

    // Jets
    TH1D* h_allJets_n_;
    TH1D* h_allJets_pt_;
    TH1D* h_allJets_phi_;
    TH1D* h_allJets_eta_;
    TH1D* h_allJets_csv_;
    TH1D* h_allJets_id_;
    // ... that pass kin cuts, loose ID
    TH1D* h_goodJets_n_;
    TH1D* h_goodJets_nb_;
    TH1D* h_goodJets_pt_;
    TH1D* h_goodJets_phi_;
    TH1D* h_goodJets_eta_;
    TH1D* h_goodJets_csv_;
    TH1D* h_goodLJets_n_;
    TH1D* h_goodLJets_nb_;
    TH1D* h_goodLJets_pt_;
    TH1D* h_goodLJets_phi_;
    TH1D* h_goodLJets_eta_;
    TH1D* h_goodLJets_csv_;
    TH1D* h_goodBJets_n_;
    TH1D* h_goodBJets_nb_;
    TH1D* h_goodBJets_pt_;
    TH1D* h_goodBJets_phi_;
    TH1D* h_goodBJets_eta_;
    TH1D* h_goodBJets_csv_;

    TH1D* h_goodMET_pt_;
    TH1D* h_goodMET_phi_;      

    //Photons! Incorporated from PatPhotonSimpleAnalyzer.h (see .cc file in the intro)    
    //std::string outputFile_;   // output file
    double minPhotonEt_;       // minimum photon Et
    double minPhotonAbsEta_;   // min and
    double maxPhotonAbsEta_;   // max abs(eta)
    double minPhotonR9_;       // minimum R9 = E(3x3)/E(SuperCluster)
    double maxPhotonHoverE_;   // maximum HCAL / ECAL
    double minIetaIeta_;       // min sigma Ieta Ieta... ECAL crystal location
    bool   createPhotonTTree_; // Create a TTree of photon variables

    // Will be used for creating TTree of photons.
    // These names did not have to match those from a phtn->...
    // but do match for clarity.
    struct struct_recPhoton {
        float isolationEcalRecHit;
        float isolationHcalRecHit;
        float isolationSolidTrkCone;
        float isolationHollowTrkCone;
        float nTrkSolidCone;
        float nTrkHollowCone;
        float isEBGap;
        float isEEGap;
        float isEBEEGap;
        float r9;
        float et;
        float eta;
        float phi;
        float hadronicOverEm;
        float ecalIso;
        float hcalIso;
        float trackIso;
    } ;
    struct_recPhoton recPhoton;

    // root file to store histograms
    //TFile*  rootFile_;

    // data members for histograms to be filled

    // PhotonID Histograms
    TH1F* h_isoEcalRecHit_;
    TH1F* h_isoHcalRecHit_;
    TH1F* h_trk_pt_solid_;
    TH1F* h_trk_pt_hollow_;
    TH1F* h_ntrk_solid_;
    TH1F* h_ntrk_hollow_;
    TH1F* h_ebgap_;
    TH1F* h_eeGap_;
    TH1F* h_ebeeGap_;
    TH1F* h_r9_;

    // Photon Histograms
    TH1F* h_photonEt_;
    TH1F* h_photonEta_;
    TH1F* h_photonPhi_;
    TH1F* h_hadoverem_;
    TH1F* h_photonIetaIeta_;

    // Photon's SuperCluster Histograms
    TH1F* h_photonScEt_;
    TH1F* h_photonScEta_;
    TH1F* h_photonScPhi_;
    TH1F* h_photonScEtaWidth_;

    // Composite or Other Histograms
    TH1F* h_photonInAnyGap_;
    TH1F* h_nPassingPho_;
    TH1F* h_nPho_;

    // TTree
    TTree* tree_PhotonAll_;

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
BasicPatDistrib::BasicPatDistrib(const edm::ParameterSet& iConfig):
  useDeepCSV_(iConfig.getParameter<bool>("useDeepCSV")),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),  
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetIDLoose_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE), 
  jetIDTight_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT), 
  metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
  genPartsToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
  photonsToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
  minPhotonEt_(iConfig.getParameter<double>("minPhotonEt")),
  minPhotonAbsEta_(iConfig.getParameter<double>("minPhotonAbsEta")),
  maxPhotonAbsEta_(iConfig.getParameter<double>("maxPhotonAbsEta")),
  minPhotonR9_(iConfig.getParameter<double>("minPhotonR9")),
  maxPhotonHoverE_(iConfig.getParameter<double>("maxPhotonHoverE")),
  minIetaIeta_(iConfig.getParameter<double>("minIetaIeta"))
  //outputFile_(iConfig.getParameter<std::string>("outputFile"))
{
  //Photon PATS 
  // Read Parameters from configuration file

  // output filename
  // Read variables that must be passed to allow a
  //  supercluster to be placed in histograms as a photon.
  
  // Read variable to that decidedes whether
  // a TTree of photons is created or not
  createPhotonTTree_ = iConfig.getParameter<bool>("createPhotonTTree");
  std::cout<<"Hello User I am running!"<<std::endl;
  // open output file to store histograms
  //rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");



  //now do what ever initialization is needed
  usesResource("TFileService");

  // MC truth in fiducial phase space
  h_genMuons_n_ = fs_->make<TH1D>("GenMuonsN",";Muon multiplicity;Events / 1", 4, 0., 4.);
  h_genMuons_pt_ = fs_->make<TH1D>("GenMuonsPt",";p_{T}(#mu) (GeV);Events / (5 GeV)", 26, 20., 150.);
  h_genMuons_phi_ = fs_->make<TH1D>("GenMuonsPhi",";#phi(#mu);Events / 0.2", 30, -3., 3.);
  h_genMuons_eta_ = fs_->make<TH1D>("GenMuonsEta",";#eta(#mu);Events / 0.2", 30, -3., 3.);
  h_genMuons_iso_ = fs_->make<TH1D>("GenMuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.01", 20, 0., 0.2); 
  h_genElecs_n_ = fs_->make<TH1D>("GenElecsN",";Electron multiplicity;Events / 1", 4, 0., 4.);
  h_genElecs_pt_ = fs_->make<TH1D>("GenElecsPt",";p_{T}(e) (GeV);Events / (5 GeV)", 26, 20., 150.);
  h_genElecs_phi_ = fs_->make<TH1D>("GenElecsPhi",";#phi(e);Events / 0.2", 30, -3., 3.);
  h_genElecs_eta_ = fs_->make<TH1D>("GenElecsEta",";#eta(e);Events / 0.2", 30, -3., 3.);
  h_genElecs_iso_ = fs_->make<TH1D>("GenElecsIso",";I_{rel}^{PUPPI}(e);Events / 0.01", 20, 0., 0.2); 
  h_genJets_n_ = fs_->make<TH1D>("GenJetsN",";Jet multiplicity;Events / 1", 14, 0., 14.);
  h_genJets_pt_ = fs_->make<TH1D>("GenJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_genJets_phi_ = fs_->make<TH1D>("GenJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_genJets_eta_ = fs_->make<TH1D>("GenJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);

  // Vertices
  h_allVertices_n_ = fs_->make<TH1D>("AllVertices",";Vertex multiplicity;Events / 1", 7, 0., 7.);
  // ... that pass ID
  h_goodVertices_n_ = fs_->make<TH1D>("GoodVertices",";Vertex multiplicity;Events / 1", 7, 0., 7.);

  // Muons
  h_allMuons_n_ = fs_->make<TH1D>("AllMuonsN",";Muon multiplicity;Events / 1", 6, 0., 6.);
  h_allMuons_pt_ = fs_->make<TH1D>("AllMuonsPt",";p_{T}(#mu) (GeV);Events / (2 GeV)", 75, 0., 150.);
  h_allMuons_phi_ = fs_->make<TH1D>("AllMuonsPhi",";#phi(#mu);Events / 0.1", 60, -3., 3.);
  h_allMuons_eta_ = fs_->make<TH1D>("AllMuonsEta",";#eta(#mu);Events / 0.1", 60, -3., 3.);
  h_allMuons_iso_ = fs_->make<TH1D>("AllMuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.01", 40, 0., 0.4);
  h_allMuons_id_ = fs_->make<TH1D>("AllMuonsID",";;Muons / 1", 4, 0., 4.);
  h_allMuons_id_->SetOption("bar");
  h_allMuons_id_->SetBarWidth(0.75);
  h_allMuons_id_->SetBarOffset(0.125);
  h_allMuons_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allMuons_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allMuons_id_->GetXaxis()->SetBinLabel(3,"Medium");
  h_allMuons_id_->GetXaxis()->SetBinLabel(4,"Tight");
  // ... that pass kin cuts, tight ID, and are isolated
  h_goodMuons_n_ = fs_->make<TH1D>("GoodMuonsN",";Muon multiplicity;Events / 1", 4, 0., 4.);
  h_goodMuons_pt_ = fs_->make<TH1D>("GoodMuonsPt",";p_{T}(#mu) (GeV);Events / (5 GeV)", 26, 20., 150.);
  h_goodMuons_phi_ = fs_->make<TH1D>("GoodMuonsPhi",";#phi(#mu);Events / 0.2", 30, -3., 3.);
  h_goodMuons_eta_ = fs_->make<TH1D>("GoodMuonsEta",";#eta(#mu);Events / 0.2", 30, -3., 3.);
  h_goodMuons_iso_ = fs_->make<TH1D>("GoodMuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.01", 20, 0., 0.2);

  // Elecs
  h_allElecs_n_ = fs_->make<TH1D>("AllElecsN",";Electron multiplicity;Events / 1", 6, 0., 6.);
  h_allElecs_pt_ = fs_->make<TH1D>("AllElecsPt",";p_{T}(e) (GeV);Events / (2 GeV)", 75, 0., 150.);
  h_allElecs_phi_ = fs_->make<TH1D>("AllElecsPhi",";#phi(e);Events / 0.1", 60, -3., 3.);
  h_allElecs_eta_ = fs_->make<TH1D>("AllElecsEta",";#eta(e);Events / 0.1", 60, -3., 3.);
  h_allElecs_iso_ = fs_->make<TH1D>("AllElecsIso",";I_{rel}^{PUPPI}(e);Events / 0.01", 40, 0., 0.4);
  h_allElecs_id_ = fs_->make<TH1D>("AllElecsID",";;Electrons / 1", 4, 0., 4.);
  h_allElecs_id_->SetOption("bar");
  h_allElecs_id_->SetBarWidth(0.75);
  h_allElecs_id_->SetBarOffset(0.125);
  h_allElecs_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allElecs_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allElecs_id_->GetXaxis()->SetBinLabel(3,"Medium");
  h_allElecs_id_->GetXaxis()->SetBinLabel(4,"Tight");
  // ... that pass kin cuts, tight ID, and are isolated
  h_goodElecs_n_ = fs_->make<TH1D>("GoodElecsN",";Electron multiplicity;Events / 1", 4, 0., 4.);
  h_goodElecs_pt_ = fs_->make<TH1D>("GoodElecsPt",";p_{T}(e) (GeV);Events / (5 GeV)", 26, 20., 150.);
  h_goodElecs_phi_ = fs_->make<TH1D>("GoodElecsPhi",";#phi(e);Events / 0.2", 30, -3., 3.);
  h_goodElecs_eta_ = fs_->make<TH1D>("GoodElecsEta",";#eta(e);Events / 0.2", 30, -3., 3.);
  h_goodElecs_iso_ = fs_->make<TH1D>("GoodElecsIso",";I_{rel}^{PUPPI}(e);Events / 0.01", 20, 0., 0.2);

  // Jets
  h_allJets_n_ = fs_->make<TH1D>("AllJetsN",";Jet multiplicity;Events / 1", 15, 0., 15.);
  h_allJets_pt_ = fs_->make<TH1D>("AllJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 100, 0., 200.);
  h_allJets_phi_ = fs_->make<TH1D>("AllJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_allJets_eta_ = fs_->make<TH1D>("AllJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_allJets_csv_ = fs_->make<TH1D>("AllJetsCSV",";CSV discriminant;Events / 0.02", 50, 0., 1.);
  h_allJets_id_ = fs_->make<TH1D>("AllJetsID",";;Jets / 1", 3, 0., 3.);
  h_allJets_id_->SetOption("bar");
  h_allJets_id_->SetBarWidth(0.75);
  h_allJets_id_->SetBarOffset(0.125);
  h_allJets_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allJets_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allJets_id_->GetXaxis()->SetBinLabel(3,"Tight");
  // ... that pass kin cuts, loose ID
  h_goodJets_n_ = fs_->make<TH1D>("GoodJetsN",";Jet multiplicity;Events / 1", 14, 0., 14.);
  h_goodJets_nb_ = fs_->make<TH1D>("GoodJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodJets_pt_ = fs_->make<TH1D>("GoodJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_goodJets_phi_ = fs_->make<TH1D>("GoodJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_goodJets_eta_ = fs_->make<TH1D>("GoodJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_goodJets_csv_ = fs_->make<TH1D>("GoodJetsCSV",";CSV discriminant;Events / 0.02", 50, 0., 1.);
  h_goodLJets_n_ = fs_->make<TH1D>("GoodLightJetsN",";Jet multiplicity;Events / 1", 12, 0., 12.);
  h_goodLJets_nb_ = fs_->make<TH1D>("GoodLightJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodLJets_pt_ = fs_->make<TH1D>("GoodLightJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_goodLJets_phi_ = fs_->make<TH1D>("GoodLightJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_goodLJets_eta_ = fs_->make<TH1D>("GoodLightJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_goodLJets_csv_ = fs_->make<TH1D>("GoodLightJetsCSV",";CSV discriminant;Events / 0.02", 50, 0., 1.);
  h_goodBJets_n_ = fs_->make<TH1D>("GoodBtaggedJetsN",";Jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodBJets_nb_ = fs_->make<TH1D>("GoodBtaggedJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodBJets_pt_ = fs_->make<TH1D>("GoodBtaggedJetsPt",";p_{T}(jet) (GeV);Events / (5 GeV)", 36, 20., 200.);
  h_goodBJets_phi_ = fs_->make<TH1D>("GoodBtaggedJetsPhi",";#phi(jet);Events / 0.2", 30, -3., 3.);
  h_goodBJets_eta_ = fs_->make<TH1D>("GoodBtaggedJetsEta",";#eta(jet);Events / 0.2", 50, -5., 5.);
  h_goodBJets_csv_ = fs_->make<TH1D>("GoodBtaggedJetsCSV",";CSV discriminant;Events / 0.01", 20, 0.8, 1.);

  // MET
  h_goodMET_pt_ = fs_->make<TH1D>("GoodMETPt",";p_{T}(MET) (GeV);Events / (5 GeV)", 60, 0., 300.);
  h_goodMET_phi_ = fs_->make<TH1D>("GoodMETPhi",";#phi(MET);Events / 0.2", 30, -3., 3.);

  // PhotonID Histograms
  h_isoEcalRecHit_ = fs_->make<TH1F>("photonEcalIso",          "Ecal Rec Hit Isolation", 100, 0, 100);
  h_isoHcalRecHit_ = fs_->make<TH1F>("photonHcalIso",          "Hcal Rec Hit Isolation", 100, 0, 100);
  h_trk_pt_solid_  = fs_->make<TH1F>("photonTrackSolidIso",    "Sum of track pT in a cone of #DeltaR" , 100, 0, 100);
  h_trk_pt_hollow_ = fs_->make<TH1F>("photonTrackHollowIso",   "Sum of track pT in a hollow cone" ,     100, 0, 100);
  h_ntrk_solid_    = fs_->make<TH1F>("photonTrackCountSolid",  "Number of tracks in a cone of #DeltaR", 100, 0, 100);
  h_ntrk_hollow_   = fs_->make<TH1F>("photonTrackCountHollow", "Number of tracks in a hollow cone",     100, 0, 100);
  h_ebgap_         = fs_->make<TH1F>("photonInEBgap",          "Ecal Barrel gap flag",  2, -0.5, 1.5);
  h_eeGap_         = fs_->make<TH1F>("photonInEEgap",          "Ecal Endcap gap flag",  2, -0.5, 1.5);
  h_ebeeGap_       = fs_->make<TH1F>("photonInEEgap",          "Ecal Barrel/Endcap gap flag",  2, -0.5, 1.5);
  h_r9_            = fs_->make<TH1F>("photonR9",               "R9 = E(3x3) / E(SuperCluster)", 300, 0, 3);

  // Photon Histograms
  h_photonEt_      = fs_->make<TH1F>("photonEt",     "Photon E_{T}",  200,  0, 200);
  h_photonEta_     = fs_->make<TH1F>("photonEta",    "Photon #eta",   200, -4,   4);
  h_photonPhi_     = fs_->make<TH1F>("photonPhi",    "Photon #phi",   200, -1.*TMath::Pi(), TMath::Pi());
  h_hadoverem_     = fs_->make<TH1F>("photonHoverE", "Hadronic over EM", 200, 0, 1);
  h_photonIetaIeta_ = fs_->make<TH1F>("photonSigmaIetaIeta","Photon #sigma_{i#etai#eta}",1500,0.0,0.3),

  // Photon's SuperCluster Histograms
  h_photonScEt_       = fs_->make<TH1F>("photonScEt",  "Photon SuperCluster E_{T}", 200,  0, 200);
  h_photonScEta_      = fs_->make<TH1F>("photonScEta", "Photon #eta",               200, -4,   4);
  h_photonScPhi_      = fs_->make<TH1F>("photonScPhi", "Photon #phi", 200, -1.*TMath::Pi(), TMath::Pi());
  h_photonScEtaWidth_ = fs_->make<TH1F>("photonScEtaWidth","#eta-width",            100,  0,  .1);

  // Composite or Other Histograms
  h_photonInAnyGap_   = fs_->make<TH1F>("photonInAnyGap",     "Photon in any gap flag",  2, -0.5, 1.5);
  h_nPassingPho_      = fs_->make<TH1F>("photonPassingCount", "Total number photons (0=NotPassing, 1=Passing)", 2, -0.5, 1.5);
  h_nPho_             = fs_->make<TH1F>("photonCount",        "Number of photons passing cuts in event",  10,  0,  10);

  // Create a TTree of photons if set to 'True' in config file
  if ( createPhotonTTree_ ) {
    tree_PhotonAll_     = fs_->make<TTree>("TreePhotonAll", "Reconstructed Photon");
    tree_PhotonAll_->Branch("recPhoton", &recPhoton.isolationEcalRecHit, "isolationEcalRecHit/F:isolationHcalRecHit:isolationSolidTrkCone:isolationHollowTrkCone:nTrkSolidCone:nTrkHollowCone:isEBGap:isEEGap:isEBEEGap:r9:et:eta:phi:hadronicOverEm:ecalIso:hcalIso:trackIso");
  }
}


BasicPatDistrib::~BasicPatDistrib()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //delete rootFile_;

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
BasicPatDistrib::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //scope   
  using namespace edm;
  using namespace std;
  //Gathering Physics Objects 
  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);

  Handle<std::vector<pat::Electron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(convToken_, conversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();  

  Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);

  Handle<std::vector<pat::MET>> mets;
  iEvent.getByToken(metsToken_, mets);

  Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  Handle<std::vector<pat::PackedGenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);

  Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genJetsToken_, genJets);

  //Handle< View<pat::Photon> >  photonHandle;
  //iEvent.getByLabel("slimmedPhotons", photonHandle);
  //View<pat::Photon> photons = *photonHandle; 
  
  Handle<std::vector<pat::Photon>>  photonHandle;
  iEvent.getByToken(photonsToken_, photonHandle);
  std::vector<pat::Photon> photons = *photonHandle; 

  // Vertices
  int prVtx = -1;
  size_t nVtx = 0;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
    ++ nVtx;
  }
  if (prVtx < 0) return;
  h_goodVertices_n_->Fill(nVtx);
  h_allVertices_n_->Fill(vertices->size());
   
  // MC truth in fiducial phase space
  std::vector<size_t> jGenJets;
  size_t nGenJets = 0;
  for (size_t i = 0; i < genJets->size(); i++) {
    bool overlaps = false;
    for (size_t j = 0; j < genParts->size(); j++) {
      if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
      if (fabs(genJets->at(i).pt()-genParts->at(j).pt()) < 0.01*genParts->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    jGenJets.push_back(i);

    if (genJets->at(i).pt() < 30.) continue;
    if (fabs(genJets->at(i).eta()) > 4.7) continue;
    h_genJets_pt_->Fill(genJets->at(i).pt());
    h_genJets_phi_->Fill(genJets->at(i).phi());
    h_genJets_eta_->Fill(genJets->at(i).eta());
    ++nGenJets;
  }
  h_genJets_n_->Fill(nGenJets);

  size_t nGenMuons = 0;
  size_t nGenElecs = 0;
  for (size_t i = 0; i < genParts->size(); i++) {
    if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
    if (fabs(genParts->at(i).eta()) > 2.8) continue;
    double genIso = 0.;
    for (size_t j = 0; j < jGenJets.size(); j++) {
      if (ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),genJets->at(jGenJets[j]).p4()) > 0.7) continue; 
      std::vector<const reco::Candidate *> jconst = genJets->at(jGenJets[j]).getJetConstituentsQuick();
      for (size_t k = 0; k < jconst.size(); k++) {
        double deltaR = ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),jconst[k]->p4());
        if (deltaR < 0.01) continue;
        if (abs(genParts->at(i).pdgId()) == 13 && deltaR > 0.4) continue;
        if (abs(genParts->at(i).pdgId()) == 11 && deltaR > 0.3) continue;
        genIso = genIso + jconst[k]->pt();
      }
    }
    genIso = genIso / genParts->at(i).pt();
    if (genIso > 0.15) continue;
    if (abs(genParts->at(i).pdgId()) == 13) {
      if (genParts->at(i).pt() < 26.) continue;
      h_genMuons_pt_->Fill(genParts->at(i).pt());
      h_genMuons_phi_->Fill(genParts->at(i).phi());
      h_genMuons_eta_->Fill(genParts->at(i).eta());
      h_genMuons_iso_->Fill(genIso); 
      ++nGenMuons;
    }
    if (abs(genParts->at(i).pdgId()) == 11) {
      if (genParts->at(i).pt() < 30.) continue;
      h_genElecs_pt_->Fill(genParts->at(i).pt());
      h_genElecs_phi_->Fill(genParts->at(i).phi());
      h_genElecs_eta_->Fill(genParts->at(i).eta());
      h_genElecs_iso_->Fill(genIso); 
      ++nGenElecs;
    }
  }
  h_genMuons_n_->Fill(nGenMuons);
  h_genElecs_n_->Fill(nGenElecs);

  // Muons
  size_t nGoodMuons = 0;
  for (size_t i = 0; i < muons->size(); i++) {
    h_allMuons_pt_->Fill(muons->at(i).pt());
    h_allMuons_phi_->Fill(muons->at(i).phi());
    h_allMuons_eta_->Fill(muons->at(i).eta());
    h_allMuons_iso_->Fill((muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt());
    h_allMuons_id_->Fill(0.);
    if ((fabs(muons->at(i).eta()) < 2.4 && muon::isLooseMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.5))) h_allMuons_id_->Fill(1.);
    if ((fabs(muons->at(i).eta()) < 2.4 && muon::isMediumMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.3))) h_allMuons_id_->Fill(2.);
    if ((fabs(muons->at(i).eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(muons->at(i),vertices->at(prVtx))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.1))) h_allMuons_id_->Fill(3.);

    if (muons->at(i).pt() < 26.) continue;
    if (fabs(muons->at(i).eta()) > 2.8) continue;
    if ((muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt() > 0.15) continue;
    if ((fabs(muons->at(i).eta()) < 2.4 && (vertices->size() <= 0 || !muon::isTightMuon(muons->at(i),vertices->at(prVtx)))) || (fabs(muons->at(i).eta()) > 2.4 && !isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.1))) continue;
    h_goodMuons_pt_->Fill(muons->at(i).pt());
    h_goodMuons_phi_->Fill(muons->at(i).phi());
    h_goodMuons_eta_->Fill(muons->at(i).eta());
    h_goodMuons_iso_->Fill((muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt());
    ++nGoodMuons;
  }
  h_goodMuons_n_->Fill(nGoodMuons);
  h_allMuons_n_->Fill(muons->size());
 
  // Electrons
  size_t nGoodElecs = 0;
  for (size_t i = 0; i < elecs->size(); i++) {
    h_allElecs_pt_->Fill(elecs->at(i).pt());
    h_allElecs_phi_->Fill(elecs->at(i).phi());
    h_allElecs_eta_->Fill(elecs->at(i).eta());
    h_allElecs_iso_->Fill((elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt());
    h_allElecs_id_->Fill(0.);
    if (isLooseElec(elecs->at(i),conversions,beamspot)) h_allElecs_id_->Fill(1.);    
    if (isMediumElec(elecs->at(i),conversions,beamspot)) h_allElecs_id_->Fill(2.);    
    if (isTightElec(elecs->at(i),conversions,beamspot)) h_allElecs_id_->Fill(3.);    

    if (elecs->at(i).pt() < 30.) continue;
    if (fabs(elecs->at(i).eta()) > 2.8) continue;
    if ((elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt() > 0.15) continue;
    if (!isTightElec(elecs->at(i),conversions,beamspot)) continue;    
    h_goodElecs_pt_->Fill(elecs->at(i).pt());
    h_goodElecs_phi_->Fill(elecs->at(i).phi());
    h_goodElecs_eta_->Fill(elecs->at(i).eta());
    h_goodElecs_iso_->Fill((elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt());
    ++nGoodElecs;
  }
  h_goodElecs_n_->Fill(nGoodElecs);
  h_allElecs_n_->Fill(elecs->size());
  
  // Jets
  size_t nGoodJets = 0;
  size_t nbGoodJets = 0;
  size_t nGoodLightJets = 0;
  size_t nbGoodLightJets = 0;
  size_t nGoodBtaggedJets = 0;
  size_t nbGoodBtaggedJets = 0;
  for (size_t i =0; i < jets->size(); i++) {
    bool overlaps = false;
    for (size_t j = 0; j < elecs->size(); j++) {
      if (fabs(jets->at(i).pt()-elecs->at(j).pt()) < 0.01*elecs->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(elecs->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    for (size_t j = 0; j < muons->size(); j++) {
      if (fabs(jets->at(i).pt()-muons->at(j).pt()) < 0.01*muons->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(muons->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;

    double btagDisc = -1.;
    if (useDeepCSV_)
        btagDisc = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
            jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
    else
        btagDisc = jets->at(i).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    
    h_allJets_pt_->Fill(jets->at(i).pt());
    h_allJets_phi_->Fill(jets->at(i).phi());
    h_allJets_eta_->Fill(jets->at(i).eta());
    h_allJets_csv_->Fill(btagDisc); 
    h_allJets_id_->Fill(0.);
    pat::strbitset retLoose = jetIDLoose_.getBitTemplate();
    retLoose.set(false);
    if (jetIDLoose_(jets->at(i), retLoose)) h_allJets_id_->Fill(1.);
    pat::strbitset retTight = jetIDTight_.getBitTemplate();
    retTight.set(false);
    if (jetIDTight_(jets->at(i), retTight)) h_allJets_id_->Fill(2.);

    if (jets->at(i).pt() < 30.) continue;
    if (fabs(jets->at(i).eta()) > 4.7) continue;
    if (!jetIDLoose_(jets->at(i), retLoose)) continue;
    h_goodJets_pt_->Fill(jets->at(i).pt());
    h_goodJets_phi_->Fill(jets->at(i).phi());
    h_goodJets_eta_->Fill(jets->at(i).eta());
    h_goodJets_csv_->Fill(btagDisc); 
    ++nGoodJets;
    if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodJets;
    if ((useDeepCSV_ && btagDisc > 0.6324)
            || (!useDeepCSV_ && btagDisc > 0.8484)){  
      h_goodBJets_pt_->Fill(jets->at(i).pt());
      h_goodBJets_phi_->Fill(jets->at(i).phi());
      h_goodBJets_eta_->Fill(jets->at(i).eta());
      h_goodBJets_csv_->Fill(btagDisc);
      ++nGoodBtaggedJets;
      if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodBtaggedJets;
    } else {
      h_goodLJets_pt_->Fill(jets->at(i).pt());
      h_goodLJets_phi_->Fill(jets->at(i).phi());
      h_goodLJets_eta_->Fill(jets->at(i).eta());
      h_goodLJets_csv_->Fill(btagDisc); 
      ++nGoodLightJets;
      if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodLightJets;
    }
  }
  h_goodLJets_n_->Fill(nGoodLightJets);
  h_goodLJets_nb_->Fill(nbGoodLightJets);
  h_goodBJets_n_->Fill(nGoodBtaggedJets);
  h_goodBJets_nb_->Fill(nbGoodBtaggedJets);
  h_goodJets_n_->Fill(nGoodJets);
  h_goodJets_nb_->Fill(nbGoodJets);
  h_allJets_n_->Fill(jets->size());
  
  // MET
  if (mets->size() > 0) {
    h_goodMET_pt_->Fill(mets->at(0).pt());
    h_goodMET_phi_->Fill(mets->at(0).phi());
  }

  // Photons
  int photonCounter = 0;


  for (int i=0; i<int(photons.size()); i++)
  {

    pat::Photon currentPhoton = photons.at(i);

    float photonEt       = currentPhoton.et();
    float superClusterEt = (currentPhoton.superCluster()->energy())/(cosh(currentPhoton.superCluster()->position().eta()));

    // Only store photon candidates (SuperClusters) that pass some simple cuts
    bool passCuts = (              photonEt < (minPhotonEt_ + 0.0047*currentPhoton.pt()    )) &&
       //             (      fabs(currentPhoton.eta()) > minPhotonAbsEta_ ) &&
         //           (      fabs(currentPhoton.eta()) < maxPhotonAbsEta_ ) &&
           //         (          currentPhoton.r9() > minPhotonR9_     ) &&
                    ( currentPhoton.hadronicOverEm() < maxPhotonHoverE_ ) &&
                    ( currentPhoton.sigmaIetaIeta() > minIetaIeta_ ) ;

    if ( passCuts )
    {
      ///////////////////////////////////////////////////////
      //                fill histograms                    //
      ///////////////////////////////////////////////////////
      // PhotonID Variables
      h_isoEcalRecHit_->Fill(currentPhoton.ecalRecHitSumEtConeDR04());
      h_isoHcalRecHit_->Fill(currentPhoton.hcalTowerSumEtConeDR04());
      h_trk_pt_solid_ ->Fill(currentPhoton.trkSumPtSolidConeDR04());
      h_trk_pt_hollow_->Fill(currentPhoton.trkSumPtHollowConeDR04());
      h_ntrk_solid_->   Fill(currentPhoton.nTrkSolidConeDR04());
      h_ntrk_hollow_->  Fill(currentPhoton.nTrkHollowConeDR04());
      h_ebgap_->        Fill(currentPhoton.isEBGap());
      h_eeGap_->        Fill(currentPhoton.isEEGap());
      h_ebeeGap_->      Fill(currentPhoton.isEBEEGap());
      h_r9_->           Fill(currentPhoton.r9());

      // Photon Variables
      h_photonEt_->  Fill(photonEt);
      h_photonEta_-> Fill(currentPhoton.eta());
      h_photonPhi_-> Fill(currentPhoton.phi());
      h_hadoverem_-> Fill(currentPhoton.hadronicOverEm());
      h_photonIetaIeta_->Fill(currentPhoton.sigmaIetaIeta());

      // Photon's SuperCluster Variables
      // eta is with respect to detector (not physics) vertex,
      // thus Et and eta are different from photon.
      h_photonScEt_->      Fill(superClusterEt);
      h_photonScEta_->     Fill(currentPhoton.superCluster()->position().eta());
      h_photonScPhi_->     Fill(currentPhoton.superCluster()->position().phi());
      h_photonScEtaWidth_->Fill(currentPhoton.superCluster()->etaWidth());

      // It passed photon cuts, mark it
      h_nPassingPho_->Fill(1.0);

      ///////////////////////////////////////////////////////
      //                fill TTree (optional)              //
      ///////////////////////////////////////////////////////
      if ( createPhotonTTree_ ) {
        recPhoton.isolationEcalRecHit    = currentPhoton.ecalRecHitSumEtConeDR04();
        recPhoton.isolationHcalRecHit    = currentPhoton.hcalTowerSumEtConeDR04();
        recPhoton.isolationSolidTrkCone  = currentPhoton.trkSumPtSolidConeDR04();
        recPhoton.isolationHollowTrkCone = currentPhoton.trkSumPtHollowConeDR04();
        recPhoton.nTrkSolidCone          = currentPhoton.nTrkSolidConeDR04();
        recPhoton.nTrkHollowCone         = currentPhoton.nTrkHollowConeDR04();
        recPhoton.isEBGap                = currentPhoton.isEBGap();
        recPhoton.isEEGap                = currentPhoton.isEEGap();
        recPhoton.isEBEEGap              = currentPhoton.isEBEEGap();
        recPhoton.r9                     = currentPhoton.r9();
        recPhoton.et                     = currentPhoton.et();
        recPhoton.eta                    = currentPhoton.eta();
        recPhoton.phi                    = currentPhoton.phi();
        recPhoton.hadronicOverEm         = currentPhoton.hadronicOverEm();
        recPhoton.ecalIso                = currentPhoton.ecalIso();
        recPhoton.hcalIso                = currentPhoton.hcalIso();
        recPhoton.trackIso               = currentPhoton.trackIso();

        // Fill the tree (this records all the recPhoton.* since
        // tree_PhotonAll_ was set to point at that.
        tree_PhotonAll_->Fill();
      }

      // Record whether it was near any module gap.
      // Very convoluted at the moment.
      bool inAnyGap = currentPhoton.isEBEEGap() || (currentPhoton.isEB()&&currentPhoton.isEBGap()) || (currentPhoton.isEE()&&currentPhoton.isEEGap());
      if (inAnyGap) {
        h_photonInAnyGap_->Fill(1.0);
      } else {
        h_photonInAnyGap_->Fill(0.0);
      }

      photonCounter++;
    }
    else
    {
      // This didn't pass photon cuts, mark it
      h_nPassingPho_->Fill(0.0);
    }

  } // End Loop over photons
  h_nPho_->Fill(photonCounter);




}


// ------------ method check that an e passes loose ID ----------------------------------
  bool
BasicPatDistrib::isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.02992) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.004119) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.05176) return false;
  if (patEl.hcalOverEcal() > 6.741) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 2.5) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 73.76) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes medium ID ----------------------------------
  bool
BasicPatDistrib::isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01609) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001766) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.03130) return false;
  if (patEl.hcalOverEcal() > 7.371) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.325) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 22.6) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes tight ID ----------------------------------
  bool
BasicPatDistrib::isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01614) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001322) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.06129) return false;
  if (patEl.hcalOverEcal() > 4.492) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.255) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 18.26) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method to improve ME0 muon ID ----------------
  bool 
BasicPatDistrib::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
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

// ------------ method called once each job just before starting event loop  ------------
  void 
BasicPatDistrib::beginJob()
{
  // go to *OUR* rootfile
  //rootFile_->cd();

  // Book Histograms
  

}

// ------------ method called once each job just after ending the event loop  ------------
  void 
BasicPatDistrib::endJob() 
{

  // go to *OUR* root file and store histograms
  ////rootFile_->cd();

  //// PhotonID Histograms
  //h_isoEcalRecHit_->Write();
  //h_isoHcalRecHit_->Write();
  //h_trk_pt_solid_-> Write();
  //h_trk_pt_hollow_->Write();
  //h_ntrk_solid_->   Write();
  //h_ntrk_hollow_->  Write();
  //h_ebgap_->     Write();
  //h_eeGap_->     Write();
  //h_ebeeGap_->   Write();
  //h_r9_->        Write();

  //// Photon Histograms
  //h_photonEt_->  Write();
  //h_photonEta_-> Write();
  //h_photonPhi_-> Write();
  //h_hadoverem_-> Write();

  //// Photon's SuperCluster Histograms
  //h_photonScEt_->      Write();
  //h_photonScEta_->     Write();
  //h_photonScPhi_->     Write();
  //h_photonScEtaWidth_->Write();

  //// Composite or Other Histograms
  //h_photonInAnyGap_->Write();
  //h_nPassingPho_->   Write();
  //h_nPho_->          Write();

  //// Write the root file (really writes the TTree)
  //rootFile_->Write();
  //rootFile_->Close();


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BasicPatDistrib::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicPatDistrib);
