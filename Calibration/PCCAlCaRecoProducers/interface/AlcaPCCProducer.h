#ifndef PCCAlCaRecoProducer_AlcaPCCProducer_h
#define PCCAlCaRecoProducer_AlcaPCCProducer_h

/**_________________________________________________________________
   class:   AlcaPCCProducer.h
   package: Calibration/TkAlCaRecoProducers
   


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)


________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PCC/interface/PCC.h"

    class AlcaPCCProducer : public edm::one::EDProducer<edm::EndLuminosityBlockProducer,
                                                             edm::one::WatchLuminosityBlocks> {
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

  PCC * thePCCob;
};

#endif
