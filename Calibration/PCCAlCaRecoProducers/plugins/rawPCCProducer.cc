/**_________________________________________________________________
class:   rawPCCProducer.cc

description: Takes reco 


authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 

________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "DataFormats/Luminosity/interface/PCC.h"
#include "DataFormats/Luminosity/interface/LumiInfo.h"
#include "Calibration/PCCAlCaRecoProducers/interface/rawPCCProducer.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "TMath.h"

//--------------------------------------------------------------------------------------------------
rawPCCProducer::rawPCCProducer(const edm::ParameterSet& iConfig)
{
    //Config Parameters from the config file
    PCCsrc_ = iConfig.getParameter<edm::ParameterSet>("rawPCCProducerParameters").getParameter<std::string>("PCCobLabel");
    ProdInst_ = iConfig.getParameter<edm::ParameterSet>("rawPCCProducerParameters").getParameter<std::string>("ProdInst");
    resetNLumi_ = iConfig.getParameter<edm::ParameterSet>("rawPCCProducerParameters").getUntrackedParameter<int>("resetEveryNLumi",-1);
    trigstring_ = iConfig.getParameter<edm::ParameterSet>("rawPCCProducerParameters").getUntrackedParameter<std::string>("trigstring","alcaLumi");
    //Initialization of Params in rawPCC
    ftotalevents = 0;
    countLumi_ = 0;

    //edm::InputTag PCCInputTag_(PCCsrc_, ProdInst_);
    edm::InputTag PCCInputTag_("alcaPCCProducer", "");
    //recoPixelClusterCounts_ alcaPCCProducer __ALCARECO
    
    //produces<reco::PixelClusterCounts, edm::InLumi>("alcaLumi");
    PCCToken=consumes<reco::PixelClusterCounts, edm::InLumi>(PCCInputTag_);
    //PCCToken=consumes<std::unique_ptr<reco::PixelClusterCounts> >(PCCsrc_);
    produces<LumiInfo, edm::InLumi>(trigstring_);
}

//--------------------------------------------------------------------------------------------------
rawPCCProducer::~rawPCCProducer(){
    //delete thePCCob;
}

//--------------------------------------------------------------------------------------------------
void rawPCCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
    //namespaces

    std::cout<<"A Print Statement"<<std::endl;

      
}

//--------------------------------------------------------------------------------------------------
void rawPCCProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    std::cout<<"Begin Lumi-Block"<<std::endl;

    //Grabbing the info from the PCC ob
    edm::Handle<reco::PixelClusterCounts> PCCHandle; 
    //edm::Handle< std::unique_ptr<reco::PixelClusterCounts> > PCCHandle; 
    lumiSeg.getByToken(PCCToken,PCCHandle);
    
    const reco::PixelClusterCounts& PCCob = *(PCCHandle.product()); 
    //const std::unique_ptr<reco::PixelClusterCounts>& PCCob = *(PCCHandle.product()); 

    //Making the vectors to loop over our end products   
    //auto_ptr<LumiInfo> alcaLumi( new LumiInfoi ); 

    clusters_ = PCCob.readCounts();
    //clusters_ = PCCob->readCounts();
    
    //std::cout<<"Print begin Lumi Block"<<std::endl;
    //New PCC object at the beginning of each lumi section
    //thePCCob = std::make_unique<reco::PixelClusterCounts>();
    countLumi_++;
    //std::cout<<"The Count Lumi "<<countLumi_<<std::endl;
    //std::cout<<"The ftotal "<<ftotalevents<<std::endl;
    std::cout<<"End of begin Lumi-Block"<<std::endl;


}

//--------------------------------------------------------------------------------------------------
void rawPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
}

//--------------------------------------------------------------------------------------------------
void rawPCCProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup){
    //std::cout<<"Print end Lumi Block"<<std::endl;

    //if ( fitNLumi_ == -1 && resetFitNLumi_ == -1 ) return;

    //if (fitNLumi_ > 0 && countLumi_%fitNLumi_!=0) return;


    //Instantiation of Pixel Object 
    lumiSeg.put(std::move(theLumiOb), std::string("alcaLumi")); 

}

DEFINE_FWK_MODULE(rawPCCProducer);
