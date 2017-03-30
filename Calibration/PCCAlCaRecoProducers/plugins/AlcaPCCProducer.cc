/**_________________________________________________________________
   class:   AlcaPCCProducer.cc
   


   authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 

   ________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "DataFormats/Luminosity/interface/PCC.h"
#include "Calibration/PCCAlCaRecoProducers/interface/AlcaPCCProducer.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "TMath.h"

//--------------------------------------------------------------------------------------------------
AlcaPCCProducer::AlcaPCCProducer(const edm::ParameterSet& iConfig)
{
    fPixelClusterLabel = iConfig.getParameter<edm::ParameterSet>("AlcaPCCProducerParameters").getParameter<edm::InputTag>("pixelClusterLabel");
    //thePCCob = new reco::PCC; 

    //std::cout<<"A Print Statement"<<std::endl;
    ftotalevents = 0;
    ftmprun0 = ftmprun = -1;
    countLumi_ = 0;
    beginLumiOfBSFit_ = endLumiOfBSFit_ = -1;
    
    produces<reco::PCC, edm::InLumi>("alcaPCC");
    pixelToken=consumes<edmNew::DetSetVector<SiPixelCluster> >(fPixelClusterLabel);
}

//--------------------------------------------------------------------------------------------------
AlcaPCCProducer::~AlcaPCCProducer(){
  //delete thePCCob;
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  ftotalevents++;
   
  unsigned int bx=iEvent.bunchCrossing();
  //std::cout<<"The Bunch Crossing"<<bx<<std::endl;
  thePCCob->eventCounter(bx);

  //Looping over the clusters and adding the counts up  
  edm::Handle< edmNew::DetSetVector<SiPixelCluster> > hClusterColl;
  iEvent.getByToken(pixelToken,hClusterColl);
        
  const edmNew::DetSetVector<SiPixelCluster>& clustColl = *(hClusterColl.product()); 
        // ----------------------------------------------------------------------
        // -- Clusters without tracks
  for (edmNew::DetSetVector<SiPixelCluster>::const_iterator isearch = clustColl.begin();  isearch != clustColl.end(); ++isearch){
        edmNew::DetSet<SiPixelCluster>  mod = *isearch;
        if(mod.empty()) { continue; }
        DetId detId = mod.id();
        
        // -- clusters on this det
        edmNew::DetSet<SiPixelCluster>::const_iterator  di;
        int nClusterCount=0;
        for (di = isearch->begin(); di != isearch->end(); ++di) {
            nClusterCount++;
        }

        int nCluster = isearch->size();
        if(nCluster!=nClusterCount) {
            std::cout<<"counting yields "<<nClusterCount<<" but the size is "<<nCluster<<"; they should match."<<std::endl;
        }
        thePCCob->Increment(detId(), bx, nCluster);
    }
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
  const edm::TimeValue_t fbegintimestamp = lumiSeg.beginTime().value();
  const std::time_t ftmptime = fbegintimestamp >> 32;
  //std::cout<<"Print begin Lumi Block"<<std::endl;
  //New PCC object at the beginning of each lumi section
  thePCCob = std::make_unique<reco::PCC>();

  if ( countLumi_ == 0 || (resetFitNLumi_ > 0 && countLumi_%resetFitNLumi_ == 0) ) {
    ftmprun0 = lumiSeg.run();
    ftmprun = ftmprun0;
    beginLumiOfBSFit_ = lumiSeg.luminosityBlock();
    refBStime[0] = ftmptime;
  }
    
  countLumi_++;
  //std::cout<<"The Count Lumi "<<countLumi_<<std::endl;
  //std::cout<<"The ftotal "<<ftotalevents<<std::endl;
  
  
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup){
  const edm::TimeValue_t fendtimestamp = lumiSeg.endTime().value();
  const std::time_t fendtime = fendtimestamp >> 32;
  refBStime[1] = fendtime;
  //std::cout<<"Print end Lumi Block"<<std::endl;
    
  endLumiOfBSFit_ = lumiSeg.luminosityBlock();
    
  if ( fitNLumi_ == -1 && resetFitNLumi_ == -1 ) return;
	
  if (fitNLumi_ > 0 && countLumi_%fitNLumi_!=0) return;

    
//Instantiation of Pixel Object 
  thePCCob->printVector();
  lumiSeg.put(std::move(thePCCob), std::string("alcaPCC")); 
  
}

DEFINE_FWK_MODULE(AlcaPCCProducer);
