#ifndef PCCAlCaRecoProducers_AlcaPCCProducer_h
#define PCCAlCaRecoProducers_AlcaPCCProducer_h

/**_________________________________________________________________
   class:   AlcaPCCProducer.h
   package: Calibration/PCCAlCaRecoProducers
   


 author: Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch)

________________________________________________________________**/


// C++ standard
#include <string>
#include <map>
// CMS
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "DataFormats/PCC/interface/PCC.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TObject.h"
#include "TH1F.h"


using namespace reco;

//Trees for files manipulation and storage
class TTree;
class DetId;


class AlcaPCCProducer : public edm::one::EDProducer<edm::EndLuminosityBlockProducer,edm::one::WatchLuminosityBlocks> {
     public:
      explicit AlcaPCCProducer(const edm::ParameterSet&);
      ~AlcaPCCProducer();

     private:
      virtual void beginLuminosityBlock     (edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup) override final;
      virtual void endLuminosityBlock       (edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup) override final;
      virtual void endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup) override final;
      virtual void produce                  (edm::Event& iEvent, const edm::EventSetup& iSetup) override final;
  
      int ftotalevents;
      int fitNLumi_;
      int resetFitNLumi_;
      int countEvt_;       //counter
      int countLumi_;      //counter
      int ftmprun0, ftmprun;
      int beginLumiOfBSFit_;
      int endLumiOfBSFit_;
      std::time_t refBStime[2];

      //bool write2DB_;

      reco::PCC * thePCCob;
      //std::unique_ptr<reco::PCC> thePCCob;
      
      edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> >  pixelToken;
      edm::EDGetTokenT<reco::VertexCollection> recoVtxToken;
      edm::EDGetTokenT<std::vector< PileupSummaryInfo> > pileUpToken;
      edm::EDGetTokenT<reco::CaloJetCollection>  hltjetsToken_;
      float *jhcalpt, *jhcalphi, *jhcaleta, *jhcale, *jhcalemf, *jhcaln90, *jhcaln90hits;
      int nhjetcal;

      edm::InputTag   fPrimaryVertexCollectionLabel;
      edm::InputTag   fPixelClusterLabel;
      edm::InputTag   fPileUpInfoLabel;

      static const int MAX_VERTICES=200;

      // saving events per LS, LN or event
      std::string saveType = "LumiSect"; // LumiSect or LumiNib or Event
      std::string sampleType="MC"; // MC or DATA
      bool saveAndReset;
      bool sameEvent;
      bool sameLumiNib;
      bool sameLumiSect;
      bool firstEvent;

          
       // Lumi stuff
      TTree * tree;
      int run;
      int LS=-99;    // set to indicate first pass of analyze method
      int LN=-99;    // set to indicate first pass of analyze method
      int event=-99; // set to indicate first pass of analyze method
      int bunchCrossing=-99;    // local variable only
      int orbit=-99;
      
      std::pair<int,int> bxModKey;    // local variable only
     
      int eventCounter=0;
      int totalEvents;
      
      bool includeVertexInformation;
      bool includePixels;
      bool includeJets;
      bool splitByBX;

      std::map<int,int> nGoodVtx;
      std::map<int,int> nValidVtx;
      std::map<std::pair<int,int>,int> nPixelClusters;
      std::map<std::pair<int,int>,int> nClusters;
      std::map<int,int> layers;

      std::map<std::pair<int,int>,float> meanPixelClusters;
      std::map<std::pair<int,int>,float> meanPixelClustersError;
      

     

      UInt_t timeStamp_begin;
      UInt_t timeStamp_local;
      UInt_t timeStamp_end;
      std::map<int,int> BXNo;


};

#endif
