/**_________________________________________________________________
class:   LumiPCCProducer.cc

description: Creates a LumiInfo object that will contain the luminosity per bunch crossing
             Along with the total lumi and the statistical error... (standard sqrt(N) stats) 


authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 

________________________________________________________________**/


// C++ standard
#include <string>
#include <iostream> 
#include <fstream> 
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
#include "FWCore/Framework/interface/IOVSyncValue.h"
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
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/LumiCorrections/interface/LumiCorrections.h"
#include "CondFormats/DataRecord/interface/LumiCorrectionsRcd.h"

class LumiPCCProducer : public edm::one::EDProducer<edm::EndRunProducer,edm::one::WatchRuns,edm::EndLuminosityBlockProducer,
                                                             edm::one::WatchLuminosityBlocks> {
  public:
    explicit LumiPCCProducer(const edm::ParameterSet&);
    ~LumiPCCProducer();

  private:
    virtual void beginRun(edm::Run const& runSeg, const edm::EventSetup& iSetup) override final;
    virtual void beginLuminosityBlock     (edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void endLuminosityBlock       (edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup) override final;
    virtual void endRun(edm::Run const& runSeg, const edm::EventSetup& iSetup)override final;
    virtual void endRunProduce(edm::Run& runSeg, const edm::EventSetup& iSetup);
    virtual void produce                  (edm::Event& iEvent, const edm::EventSetup& iSetup) override final;
 
    edm::EDGetTokenT<LumiInfo>  LumiToken;
    std::string   PCCsrc_;//input file EDproducer module label 
    std::string   ProdInst_;//input file product instance 

    std::string trigstring_; //specifies the trigger Rand or ZeroBias 
    std::string label_;
    std::vector<float> rawlumiBX_;//new vector containing clusters per bxid 
    std::vector<float> correctedLumiBX_;
    std::vector<float> corrlist_;

    float totalLumi;
     
    //New output object
    std::unique_ptr<LumiInfo> theLumiOb;

    std::ofstream csvfile;

};

//--------------------------------------------------------------------------------------------------
LumiPCCProducer::LumiPCCProducer(const edm::ParameterSet& iConfig)
{
    //Config Parameters from the config file
    PCCsrc_ = iConfig.getParameter<edm::ParameterSet>("LumiPCCProducerParameters").getParameter<std::string>("PCCobLabel");
    ProdInst_ = iConfig.getParameter<edm::ParameterSet>("LumiPCCProducerParameters").getParameter<std::string>("ProdInst");
    trigstring_ = iConfig.getParameter<edm::ParameterSet>("LumiPCCProducerParameters").getUntrackedParameter<std::string>("trigstring","alcaLumi");
    label_ = iConfig.getParameter<edm::ParameterSet>("LumiPCCProducerParameters").getParameter<std::string>("label");
    //Initialization of Params in rawPCC

    edm::InputTag LumiInputTag_(PCCsrc_, ProdInst_);
    produces<LumiInfo, edm::InLumi>(trigstring_);
    LumiToken=consumes<LumiInfo, edm::InLumi>(LumiInputTag_);

    totalLumi = 0.0;

    std::ofstream csvfile;//
}

//--------------------------------------------------------------------------------------------------
LumiPCCProducer::~LumiPCCProducer(){
}
//--------------------------------------------------------------------------------------------------
void LumiPCCProducer::beginRun(edm::Run const& runSeg, const edm::EventSetup& iSetup){
    std::cout<<"Begin Run"<<std::endl;
    csvfile<<std::to_string(runSeg.run())<<",";
    //LumiInfo outLumiOb; 
}
//--------------------------------------------------------------------------------------------------
void LumiPCCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
    //namespaces

    //std::cout<<"A Print Statement"<<std::endl;

      
}

//--------------------------------------------------------------------------------------------------
void LumiPCCProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    std::cout<<"Begin Lumi-Block"<<std::endl;
    theLumiOb = std::make_unique<LumiInfo>(); 
    //LumiInfo theLumiOb; 

    csvfile.open(label_, std::ios_base::app);
    if (csvfile.is_open()) { std::cout<<"My File is open!!"<<std::endl;}
    csvfile<<std::to_string(lumiSeg.luminosityBlock())<<",";
}

//--------------------------------------------------------------------------------------------------
void LumiPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
}

//--------------------------------------------------------------------------------------------------
void LumiPCCProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup){
    //reset parameters 
    rawlumiBX_.resize(LumiConstants::numBX,0);//new vector containing clusters per bxid 
    correctedLumiBX_.resize(LumiConstants::numBX,0);

    //Reading the info from the Lumi Info Object  
    edm::Handle<LumiInfo> PCCHandle; 
    lumiSeg.getByToken(LumiToken,PCCHandle);
    
    const LumiInfo& inLumiOb = *(PCCHandle.product()); 

    rawlumiBX_= inLumiOb.getInstLumiAllBX();
    //Reading the info from the data base
    edm::ESHandle< LumiCorrections > corrHandle;//The corrections stored per bunch crossing 
    iSetup.get<LumiCorrectionsRcd>().get(corrHandle);
    const LumiCorrections *mycorrections = corrHandle.product();
    corrlist_ = mycorrections->GetCorrectionsBX();

    std::cout<<"End Luminosity Block"<<std::endl;

    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        correctedLumiBX_[bx] = rawlumiBX_[bx]*corrlist_[bx];
        totalLumi+= correctedLumiBX_[bx];
        //output to csv file:
    }

    
    //theLumiOb->setTotalLumi(correctedLumiBX_);
    theLumiOb->setInstLumi(correctedLumiBX_);
 
    csvfile<<std::to_string(totalLumi);

    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
       csvfile<<","<<std::to_string(correctedLumiBX_[bx]);
    }
    csvfile<<std::endl;   
    lumiSeg.put(std::move(theLumiOb), std::string(trigstring_)); 
    totalLumi = 0.0;
    
    std::cout<<"writing the file"<<std::endl;
    csvfile.close();
}
//--------------------------------------------------------------------------------------------------
void LumiPCCProducer::endRun(edm::Run const& runSeg, const edm::EventSetup& iSetup){
    //std::cout<<"End Run"<<std::endl;
    //std::cout<<"The end time of the run: "<<std::to_string(runSeg.endTime())<<std::endl;
    //csvfile<<std::endl;   
 
}

//--------------------------------------------------------------------------------------------------
void LumiPCCProducer::endRunProduce(edm::Run& runSeg, const edm::EventSetup& iSetup){
    

}
DEFINE_FWK_MODULE(LumiPCCProducer);
