/**_________________________________________________________________
   class:   AlcaPCCProducer.cc
   


   authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 

   ________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "DataFormats/PCC/interface/PCC.h"
#include "Calibration/PCCAlCaRecoProducers/interface/AlcaPCCProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TMath.h"

//--------------------------------------------------------------------------------------------------
AlcaPCCProducer::AlcaPCCProducer(const edm::ParameterSet& iConfig){
  // get parameter
  //write2DB_        = iConfig.getParameter<edm::ParameterSet>("AlcaPCCProducerParameters").getParameter<bool>("WriteToDB");
	
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
  ftotalevents++;
   
  unsigned int bx=iEvent.bunchCrossing();
  std::cout<<"The Bunch Crossing"<<bx<<std::endl;
  //thePCCob->Increment(iEvent);
  thePCCob->eventCounter(bx);
  //ftmprun = iEvent.id().run();
  //alcaPCC.Increment(iEvent)
  
  //The addition of incrementing the counts and collecting vertices
  
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
