/**_________________________________________________________________
   class:   AlcaPCCProducer.cc
   


   authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 

   ________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "DataFormats/PCC/interface/PCC.h"
#include "Calibration/PCCAlCaRecoProducers/interface/AlcaPCCProducer.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
//
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//ROOT
#include "TMath.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
//--------------------------------------------------------------------------------------------------
AlcaPCCProducer::AlcaPCCProducer(const edm::ParameterSet& iConfig){
  // get parameter
  //write2DB_        = iConfig.getParameter<edm::ParameterSet>("AlcaPCCProducerParameters").getParameter<bool>("WriteToDB");
  //Adding the branches to tree for the number of pixel clusters.

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree","Pixel Cluster Counters");
  tree->Branch("run",&run,"run/I");
  tree->Branch("LS",&LS,"LS/I");
  tree->Branch("LN",&LN,"LN/I");
  tree->Branch("timeStamp_begin",&timeStamp_begin,"timeStamp_begin/i");
  tree->Branch("timeStamp_end",&timeStamp_end,"timeStamp_end/i");
  tree->Branch("eventCounter",&eventCounter,"eventCounter/I");
  tree->Branch("BXNo","map<int,int>",&BXNo);
  tree->Branch("nPixelClusters","map<std::pair<int,int>,int>",&nPixelClusters);
  tree->Branch("nClusters",     "map<std::pair<int,int>,int>",&nClusters);
  //tree->Branch("nPixelClusters","map<int,int>",&nPixelClusters);
  //tree->Branch("nClusters","map<int,int>",&nClusters);
  tree->Branch("layers","map<int,int>",&layers);
  //    
  tree->Branch("event",&event,"event/i");
  tree->Branch("orbit",&orbit,"orbit/I");
  tree->Branch("bunchCrossing",&bunchCrossing,"bunchCrossing/I");


  thePCCob = new reco::PCC; 
 
  ftotalevents = 0;
  ftmprun0 = ftmprun = -1;
  countLumi_ = 0;
  beginLumiOfBSFit_ = endLumiOfBSFit_ = -1;
  
  produces<reco::PCC, edm::InLumi>("alcaPCC");
}

//--------------------------------------------------------------------------------------------------
AlcaPCCProducer::~AlcaPCCProducer(){
  delete thePCCob;
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  //Scope
  using namespace edm;
  using reco::VertexCollection;



  ftotalevents++;
  // 
  unsigned int bx=iEvent.bunchCrossing();
  std::cout<<"The Bunch Crossing"<<bx<<std::endl;
  std::cout<<"The Bunch Crossing"<<std::endl;
  thePCCob->eventCounter(bx);

  // Get the Run, Lumi Section, and Event numbers, etc.
  run   = iEvent.id().run();
  LS    = iEvent.getLuminosityBlock().luminosityBlock();
  //LN    = -99; // FIXME need the luminibble
  event = iEvent.id().event();
  bunchCrossing   = iEvent.bunchCrossing();

  //Looping over the clusters and adding the counts up  
  edm::Handle< edmNew::DetSetVector<SiPixelCluster> > hClusterColl;
  iEvent.getByToken(pixelToken,hClusterColl);
        
  const edmNew::DetSetVector<SiPixelCluster>& clustColl = *(hClusterColl.product()); 
        // ----------------------------------------------------------------------
        // -- Clusters without tracks

for (TrackerGeometry::DetContainer::const_iterator it = TG->dets().begin(); it != TG->dets().end(); it++){
    //if (dynamic_cast<PixelGeomDetUnit*>((*it)) != 0){ 
        DetId detId = (*it)->geographicalId();

        bxModKey.second=detId();

        // -- clusters on this det
        edmNew::DetSetVector<SiPixelCluster>::const_iterator isearch = clustColl.find(detId);
        if (isearch != clustColl.end()) {  // Not an empty iterator
            edmNew::DetSet<SiPixelCluster>::const_iterator  di;
            for (di = isearch->begin(); di != isearch->end(); ++di) {
                if(nPixelClusters.count(bxModKey)==0){
                    nPixelClusters[bxModKey]=0;
                }
                nPixelClusters[bxModKey] = nPixelClusters[bxModKey]+1;
            }

            int nCluster = isearch->size();
            if(nClusters.count(bxModKey)==0){
                nClusters[bxModKey]=0;
            }
            nClusters[bxModKey] += nCluster;

            if (detId.subdetId() == PixelSubdetector::PixelBarrel) {
                PixelBarrelName detName = PixelBarrelName(detId);
                int layer = detName.layerName();
                if(layers.count(detId())==0){
                    layers[detId()]=layer;
                }
            } else {
                assert(detId.subdetId() == PixelSubdetector::PixelEndcap);
                PixelEndcapName detName = PixelEndcapName(detId);
                int disk = detName.diskName();
                if(layers.count(detId())==0){
                    layers[detId()]=disk+3;
                }
            }
        }
    //}
}

  
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
  const edm::TimeValue_t fbegintimestamp = lumiSeg.beginTime().value();
  const std::time_t ftmptime = fbegintimestamp >> 32;

  if ( countLumi_ == 0 || (resetFitNLumi_ > 0 && countLumi_%resetFitNLumi_ == 0) ) {
    ftmprun0 = lumiSeg.run();
    ftmprun = ftmprun0;
    beginLumiOfBSFit_ = lumiSeg.luminosityBlock();
    refBStime[0] = ftmptime;
  }
    
  countLumi_++;
  
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup){
  const edm::TimeValue_t fendtimestamp = lumiSeg.endTime().value();
  const std::time_t fendtime = fendtimestamp >> 32;
  refBStime[1] = fendtime;
    
  endLumiOfBSFit_ = lumiSeg.luminosityBlock();
    
  if ( fitNLumi_ == -1 && resetFitNLumi_ == -1 ) return;
	
  if (fitNLumi_ > 0 && countLumi_%fitNLumi_!=0) return;

    
//Instantiation of Pixel Object 
  auto result = std::make_unique<reco::PCC>();
  *result = *thePCCob;
  lumiSeg.put(std::move(result), std::string("alcaPCC")); 
  delete thePCCob;
  
  thePCCob = new reco::PCC; 
}

DEFINE_FWK_MODULE(AlcaPCCProducer);
