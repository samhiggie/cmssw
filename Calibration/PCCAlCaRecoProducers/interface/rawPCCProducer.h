#ifndef PCCAlCaRecoProducers_rawPCCProducer_h
#define PCCAlCaRecoProducers_rawPCCProducer_h

/**_________________________________________________________________
   class:   rawPCCProducer.h
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
#include "DataFormats/Luminosity/interface/LumiInfo.h"

class rawPCCProducer : public edm::one::EDProducer<edm::EndLuminosityBlockProducer,
                                                             edm::one::WatchLuminosityBlocks> {
  public:
    explicit rawPCCProducer(const edm::ParameterSet&);
    ~rawPCCProducer();

  private:
    virtual void beginLuminosityBlock     (edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void endLuminosityBlock       (edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void produce                  (edm::Event& iEvent, const edm::EventSetup& iSetup) override final;
 
    edm::EDGetTokenT<reco::PixelClusterCounts>  PCCToken;
    //edm::EDGetTokenT< std::unique_ptr<reco::PixelClusterCounts> >  PCCToken;
    edm::InputTag   PCCsrc_;
    edm::InputTag   ProdInst_;

    std::string trigstring_; //specifies the trigger Rand or ZeroBias 
    std::vector<int> clusters_;//Will fill this with content from PCC
    int ftotalevents;
    int resetNLumi_;
    int countEvt_;       //counter
    int countLumi_;      //counter

    //New output object
    std::unique_ptr<LumiInfo> theLumiOb;

    //Input Information 
    //typedef std::vector<std::vector<int> > clusters_;
    std::unique_ptr<reco::PixelClusterCounts> thePCCob;

};

#endif
