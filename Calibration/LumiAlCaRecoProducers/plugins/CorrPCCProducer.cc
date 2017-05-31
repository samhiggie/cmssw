/**_________________________________________________________________
class:   CorrPCCProducer.cc

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
#include "FWCore/Framework/interface/Run.h"
#include "TMath.h"

class CorrPCCProducer : public edm::one::EDProducer<edm::EndRunProducer,edm::one::WatchRuns,edm::EndLuminosityBlockProducer,edm::one::WatchLuminosityBlocks> {
  public:
    explicit CorrPCCProducer(const edm::ParameterSet&);
    ~CorrPCCProducer();

  private:
    virtual void MakeCorrections ();
    virtual void beginRun(edm::Run const& runSeg, const edm::EventSetup& iSetup) override final;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup);
    virtual void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup);
    virtual void endRun(edm::Run const& runSeg, const edm::EventSetup& iSetup)override final;
    virtual void endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void endRunProduce(edm::Run& runSeg, const edm::EventSetup& iSetup);
    virtual void produce                  (edm::Event& iEvent, const edm::EventSetup& iSetup) override final;
 
    //Old input object
    edm::EDGetTokenT<LumiInfo>  LumiToken;
    std::string   PCCsrc_;//input file EDproducer module label 
    std::string   ProdInst_;//input file product instance 

    std::string trigstring_; //specifies the trigger Rand or ZeroBias 
    std::vector<float> rawlumiBX_;//new vector containing clusters per bxid 
    std::vector<float> errOnLumiByBX_;//standard error per bx
    std::vector<float> correctionList_;//list of scale factors to apply.
    float totalLumi_;//The total raw luminosity from the pixel clusters - not scaled
    float statErrOnLumi_;//the statistical error on the lumi - large num ie sqrt(N)
    int countLumi_;//The lumisection count... the size of the lumiblock
    int resetNLumi_;//The number of lumisections per block. 

    
    //New output object
    std::unique_ptr<LumiInfo> outLumiOb;


};

//--------------------------------------------------------------------------------------------------
CorrPCCProducer::CorrPCCProducer(const edm::ParameterSet& iConfig)
{
    //Config Parameters from the config file
    PCCsrc_ = iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<std::string>("inLumiObLabel");
    ProdInst_ = iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<std::string>("ProdInst");
    trigstring_ = iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getUntrackedParameter<std::string>("trigstring","alcaLumi");
    resetNLumi_=iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<int>("resetEveryNLumi");
    //Initialization of Params
    countLumi_=0;

    //Input tag for raw lumi
    edm::InputTag PCCInputTag_(PCCsrc_, ProdInst_);

    LumiToken=consumes<LumiInfo, edm::InLumi>(PCCInputTag_);
    produces<LumiInfo, edm::InLumi>(trigstring_);
}

//--------------------------------------------------------------------------------------------------
CorrPCCProducer::~CorrPCCProducer(){
}
//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::MakeCorrections(){
    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        correctionList_.push_back(1.0);
    }
}


//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
    //namespaces

    //std::cout<<"A Print Statement"<<std::endl;

      
}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    std::cout<<"Begin Lumi-Block"<<std::endl;
    outLumiOb = std::make_unique<LumiInfo>(); 
    //LumiInfo outLumiOb; 
    countLumi_++;

}
//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::beginRun(edm::Run const& runSeg, const edm::EventSetup& iSetup){
    std::cout<<"Begin Run"<<std::endl;
    outLumiOb = std::make_unique<LumiInfo>(); 
    //LumiInfo outLumiOb; 

}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    //check to see if end of lumisection

    if (resetNLumi_ > 0 && countLumi_%resetNLumi_!=0) return;



    edm::Handle<LumiInfo> PCCHandle; 
    lumiSeg.getByToken(LumiToken,PCCHandle);
    
    const LumiInfo& inLumiOb = *(PCCHandle.product()); 

    //Making the vectors to loop over the lumisections for each run    
    rawlumiBX_= inLumiOb.getInstLumiAllBX();
    std::cout<<"The total Luminosity "<<inLumiOb.totalrawLuminosity()<<std::endl;
    //Example of forloop
    //for (unsigned int i=0;i<modID_.size();i++){
        
    //summing over modules
    
    std::cout<<"Print end Lumi Block"<<std::endl;


}
//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endRun(edm::Run const& runSeg, const edm::EventSetup& iSetup){
    //std::cout<<"End Run"<<std::endl;
    //outLumiOb = std::make_unique<LumiInfo>(); 
    //LumiInfo outLumiOb; 

}
//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup){

    lumiSeg.put(std::move(outLumiOb), std::string(trigstring_)); 

}
//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endRunProduce(edm::Run& runSeg, const edm::EventSetup& iSetup){
    std::cout<<"End Run"<<std::endl;
    //outLumiOb = std::make_unique<LumiInfo>(); 
    //LumiInfo outLumiOb; 
    //place a save here.
    runSeg.put(std::move(outLumiOb), std::string(trigstring_)); 

}


DEFINE_FWK_MODULE(CorrPCCProducer);
