#ifndef PCCAlCaRecoProducers_AlcaPCCProducer_h
#define PCCAlCaRecoProducers_AlcaPCCProducer_h

/**_________________________________________________________________
   class:   AlcaPCCProducer.h
   package: Calibration/PCCAlCaRecoProducers
   

 authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 

________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Luminosity/interface/PCC.h"

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
 
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> >  pixelToken;
    edm::InputTag   fPixelClusterLabel;

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

    //reco::PCC * thePCCob;
    std::unique_ptr<reco::PCC> thePCCob;

};

#endif
