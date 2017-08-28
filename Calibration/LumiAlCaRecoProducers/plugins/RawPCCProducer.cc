/**_________________________________________________________________
class:   RawPCCProducer.cc

description: Creates a LumiInfo object that will contain the luminosity per bunch crossing
             Along with the total lumi and the statistical error... (standard sqrt(N) stats) 


authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 

________________________________________________________________**/


// C++ standard
#include <string>
#include <vector>
// CMS
#include "DataFormats/Luminosity/interface/PixelClusterCounts.h"
#include "DataFormats/Luminosity/interface/LumiInfo.h"
#include "DataFormats/Luminosity/interface/LumiConstants.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
// CMS
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "TMath.h"

class RawPCCProducer : public edm::one::EDProducer<edm::EndLuminosityBlockProducer,
                                                             edm::one::WatchLuminosityBlocks> {
  public:
    explicit RawPCCProducer(const edm::ParameterSet&);
    ~RawPCCProducer();

  private:
    virtual void beginLuminosityBlock     (edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void endLuminosityBlock       (edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void produce                  (edm::Event& iEvent, const edm::EventSetup& iSetup) override final;
 
    edm::EDGetTokenT<reco::PixelClusterCounts>  PCCToken;
    std::string   PCCsrc_;//input file EDproducer module label 
    std::string   ProdInst_;//input file product instance 

    std::vector<int>   modVeto_;//The list of modules to skip in the lumi calc. 
    std::string trigstring_; //specifies the trigger Rand or ZeroBias 
    std::vector<int> clusters_;//Will fill this with content from PCC
    std::vector<int> modID_;//vector with Module IDs 1-1 map to bunch x-ing in clusers_
    std::vector<int> events_;//vector with total events at each bxid.
    std::vector<int> clusterPerBX_;//new vector containing clusters per bxid 
    std::vector<float> rawlumiBX_;//new vector containing clusters per bxid 
    std::vector<float> errOnLumiByBX_;//standard error per bx
    std::vector<int> goodMods_;//The indicies of all the good modules - not vetoed
    float totalLumi_;//The total raw luminosity from the pixel clusters - not scaled
    float statErrOnLumi_;//the statistical error on the lumi - large num ie sqrt(N)

    
    //New output object
    std::unique_ptr<LumiInfo> theLumiOb;


};

//--------------------------------------------------------------------------------------------------
RawPCCProducer::RawPCCProducer(const edm::ParameterSet& iConfig)
{
    //Config Parameters from the config file
    PCCsrc_ = iConfig.getParameter<edm::ParameterSet>("RawPCCProducerParameters").getParameter<std::string>("PCCobLabel");
    ProdInst_ = iConfig.getParameter<edm::ParameterSet>("RawPCCProducerParameters").getParameter<std::string>("ProdInst");
    trigstring_ = iConfig.getParameter<edm::ParameterSet>("RawPCCProducerParameters").getUntrackedParameter<std::string>("trigstring","alcaLumi");
    modVeto_ = iConfig.getParameter<edm::ParameterSet>("RawPCCProducerParameters").getParameter<std::vector<int>>("modVeto");
    //Initialization of Params in rawPCC

    edm::InputTag PCCInputTag_(PCCsrc_, ProdInst_);
    clusterPerBX_.resize(LumiConstants::numBX,0);//new vector containing clusters per bxid 
    rawlumiBX_.resize(LumiConstants::numBX,0);//new vector containing clusters per bxid 
    errOnLumiByBX_.resize(LumiConstants::numBX,0);
    PCCToken=consumes<reco::PixelClusterCounts, edm::InLumi>(PCCInputTag_);
    produces<LumiInfo, edm::InLumi>(trigstring_);
}

//--------------------------------------------------------------------------------------------------
RawPCCProducer::~RawPCCProducer(){
}

//--------------------------------------------------------------------------------------------------
void RawPCCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

      
}

//--------------------------------------------------------------------------------------------------
void RawPCCProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    theLumiOb = std::make_unique<LumiInfo>(); 

}

//--------------------------------------------------------------------------------------------------
void RawPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    totalLumi_=0.0;
    statErrOnLumi_=0.0;
    clusterPerBX_.resize(LumiConstants::numBX,0);//new vector containing clusters per bxid 
    rawlumiBX_.resize(LumiConstants::numBX,0);//new vector containing clusters per bxid 
    errOnLumiByBX_.resize(LumiConstants::numBX,0);
    goodMods_.clear();

    edm::Handle<reco::PixelClusterCounts> PCCHandle; 
    lumiSeg.getByToken(PCCToken,PCCHandle);
    
    const reco::PixelClusterCounts& PCCob = *(PCCHandle.product()); 

    //Making the vectors to loop over our end products   

    modID_ = PCCob.readModID();
    events_= PCCob.readEvents();
    clusters_ = PCCob.readCounts();
    

    for (unsigned int i=0;i<modID_.size();i++){
        if (std::find(modVeto_.begin(),modVeto_.end(), modID_.at(i)) == modVeto_.end()){
            goodMods_.push_back(i);
        }
    }

    //summing over modules
    for (int bx=0;bx<int(LumiConstants::numBX);bx++){
        for (unsigned int i=0;i<goodMods_.size();i++){
            clusterPerBX_.at(bx)+=clusters_.at(goodMods_.at(i)*int(LumiConstants::numBX)+bx);
           
        }
        if (clusterPerBX_.at(bx)!=0){
            errOnLumiByBX_.at(bx)=1/TMath::Sqrt(clusterPerBX_.at(bx));
        }
        else{
            errOnLumiByBX_.at(bx)=0.0;
        }
    }
    std::cout<<"Print end Lumi Block"<<std::endl;
    for (unsigned int i=0;i<clusterPerBX_.size();i++){
        if (events_.at(i)!=0){
            //rawlumiBX_.at(i)=clusterPerBX_.at(i)/float(events_.at(i));
            rawlumiBX_.at(i)=clusterPerBX_.at(i);
        }
        else{
            rawlumiBX_.at(i)=0.0;//Set the lumi to 0 if there is no events...
        }
        totalLumi_+=rawlumiBX_.at(i);        
        statErrOnLumi_+=float(events_.at(i));       
      
    }
    if (statErrOnLumi_!=0){statErrOnLumi_=1/TMath::Sqrt(statErrOnLumi_);}


    theLumiOb->setTotalLumi(totalLumi_);
    theLumiOb->setStatErrorOnLumi(statErrOnLumi_);

    //Note! Setting the error on the lumi per bx as the number of events... we can compute the true luminosity and stats later... using as place holders for time being. 
    std::vector<float> events(events_.begin(), events_.end());
    theLumiOb->setErrLumiBX(events);
    theLumiOb->setInstLumi(rawlumiBX_);
 
    
 

}

//--------------------------------------------------------------------------------------------------
void RawPCCProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup){
   lumiSeg.put(std::move(theLumiOb), std::string(trigstring_)); 

}

DEFINE_FWK_MODULE(RawPCCProducer);
