/**_________________________________________________________________
   class:   AlcaPCCProducer.cc
   package: RecoVertex/PCCProducer
   


   author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
   Geng-Yuan Jeng, UC Riverside (Geng-Yuan.Jeng@cern.ch)


   ________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "DataFormats/PCC/interface/PCC.h"
#include "Calibration/TkAlCaRecoProducers/interface/AlcaPCCProducer.h"
#include "RecoVertex/PCCProducer/interface/BSFitter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "CondFormats/DataRecord/interface/PCCObjectsRcd.h"
#include "CondFormats/PCCObjects/interface/PCCObjects.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TMath.h"

//--------------------------------------------------------------------------------------------------
AlcaPCCProducer::AlcaPCCProducer(const edm::ParameterSet& iConfig){
  // get parameter
  write2DB_        = iConfig.getParameter<edm::ParameterSet>("AlcaPCCProducerParameters").getParameter<bool>("WriteToDB");

  ftotalevents = 0;
  ftmprun0 = ftmprun = -1;
  countLumi_ = 0;
  beginLumiOfBSFit_ = endLumiOfBSFit_ = -1;
  
  produces<reco::PCC, edm::InLumi>("alcaPCC");
}

//--------------------------------------------------------------------------------------------------
AlcaPCCProducer::~AlcaPCCProducer(){
  delete theBeamFitter;
}

//--------------------------------------------------------------------------------------------------
void AlcaPCCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  ftotalevents++;
  //theBeamFitter->readEvent(iEvent);
  //ftmprun = iEvent.id().run();
  alcaPCC.Increment(iEvent)
  
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

  theBeamFitter->setFitLSRange(beginLumiOfBSFit_,endLumiOfBSFit_);
  theBeamFitter->setRefTime(refBStime[0],refBStime[1]);
  theBeamFitter->setRun(ftmprun0);
    
  std::pair<int,int> LSRange = theBeamFitter->getFitLSRange();

  reco::PCC bs;
  if (theBeamFitter->runPVandTrkFitter()){
    bs = theBeamFitter->getPCC();
    edm::LogInfo("AlcaPCCProducer")
        << "\n RESULTS OF DEFAULT FIT " << std::endl
        << " for runs: " << ftmprun0 << " - " << ftmprun << std::endl
        << " for lumi blocks : " << LSRange.first << " - " << LSRange.second << std::endl
        << " lumi counter # " << countLumi_ << std::endl
        << bs << std::endl
        << "fit done. \n" << std::endl;	
  }
  else { // Fill in empty beam spot if beamfit fails
    bs.setType(reco::PCC::Fake);
    edm::LogInfo("AlcaPCCProducer")
        << "\n Empty Beam spot fit" << std::endl
        << " for runs: " << ftmprun0 << " - " << ftmprun << std::endl
        << " for lumi blocks : " << LSRange.first << " - " << LSRange.second << std::endl
        << " lumi counter # " << countLumi_ << std::endl
        << bs << std::endl
        << "fit failed \n" << std::endl;
  }

  auto result = std::make_unique<reco::PCC>();
  *result = bs;
  lumiSeg.put(std::move(result), std::string("alcaPCC"));
	
  if (resetFitNLumi_ > 0 && countLumi_%resetFitNLumi_ == 0) {
    std::vector<BSTrkParameters> theBSvector = theBeamFitter->getBSvector();
    edm::LogInfo("AlcaPCCProducer")
        << "Total number of tracks accumulated = " << theBSvector.size() << std::endl
        << "Reset track collection for beam fit" <<std::endl;
    theBeamFitter->resetTrkVector();
    theBeamFitter->resetLSRange();
    theBeamFitter->resetCutFlow();
    theBeamFitter->resetRefTime();
    theBeamFitter->resetPVFitter();
    countLumi_=0;
  }
}

DEFINE_FWK_MODULE(AlcaPCCProducer);
