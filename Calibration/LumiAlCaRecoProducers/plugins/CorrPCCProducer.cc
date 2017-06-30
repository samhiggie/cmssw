/**_________________________________________________________________
class:   CorrPCCProducer.cc
description: Computes the type1 corrections to the luminosity
authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 
________________________________________________________________**/

// C++ standard
#include <string>
#include <vector>
// CMS
#include "DataFormats/Luminosity/interface/PixelClusterCounts.h"
#include "DataFormats/Luminosity/interface/LumiInfo.h"
#include "DataFormats/Luminosity/interface/LumiConstants.h"
//Test for intellegent lumiblocks
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
    void MakeCorrectionTemplate ();
    float GetMaximum(std::vector<float>);
    void EstimateType1Frac(std::vector<float>, float& );
    void CalculateCorrections (std::vector<float>, std::vector<float>&, float&);
    std::vector<float>& MakeCorrections (std::vector<float>&);
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

    std::string trigstring_; //specifices the iov LSs for the object that is saved.  
    std::vector<float> rawlumiBX_;//new vector containing clusters per bxid 
    std::vector<float> errOnLumiByBX_;//standard error per bx
    std::vector<float> totalLumiByBX_;//summed lumi
    std::vector<float> correctionTemplate_;
    std::vector<float> corr_list_;//list of scale factors to apply.
    float Overall_corr;//The Overall correction to the integrated luminosity
    float totalLumi_;//The total raw luminosity from the pixel clusters - not scaled
    float statErrOnLumi_;//the statistical error on the lumi - large num ie sqrt(N)
    float type1Frac;
    float mean_type1;
    float mean_type2; 
    int countLumi_;//The lumisection count... the size of the lumiblock
    int resetNLumi_;//The number of lumisections per block.
    int LSrun_;//The parameters that will save beginning and end LS for saving corr to runs
    int startLS;//Starting lumisection for the iov that we save with the lumiInfo object
    int endLS;//Ending lumisection for the iov that we save with the lumiInfo object.

    //double type1_; //Initial type 1 correction factor
    double type2_a_;//amplitude for the type 2 correction 
    double type2_b_;//decay width for the type 2 correction

    
    //output 
    std::unique_ptr<LumiInfo> outLumiOb;//lumi object with corrections per BX
    std::unique_ptr<float> Type1frac;
    std::unique_ptr<float> Type1res;
    std::unique_ptr<float> Type2res;
    std::unique_ptr<float> T1fUnc;
    std::unique_ptr<float> T1rUnc;
    std::unique_ptr<float> T2rUnc;
};

//--------------------------------------------------------------------------------------------------
CorrPCCProducer::CorrPCCProducer(const edm::ParameterSet& iConfig)
{
    //Config Parameters from the config file
    PCCsrc_ = iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<std::string>("inLumiObLabel");
    ProdInst_ = iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<std::string>("ProdInst");
    trigstring_ = iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getUntrackedParameter<std::string>("trigstring","alcaLumi");
    resetNLumi_=iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<int>("resetEveryNLumi");
    //type1_= iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<double>("type1");
    type2_a_ = iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<double>("type2_a");
    type2_b_ = iConfig.getParameter<edm::ParameterSet>("CorrPCCProducerParameters").getParameter<double>("type2_b");
    //Initialization of Params
    countLumi_=0;
    //resetNLumi_=10;
    //the LS range for the product instance

    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        totalLumiByBX_.push_back(0);
    }
  
    //Initialization of Temparory Corrected PCC
     
    //for (size_t bx=0; bx<LumiConstants::numBX; bx++){
    //    corrected_tmp_.push_back(0);    
   
    //} 
    //Generate a pseudo correction list:Add function to created list here?
    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        corr_list_.push_back(1.0);
    }
   
    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        correctionTemplate_.push_back(1.0);
    
    }

    //Input tag for raw lumi
    edm::InputTag PCCInputTag_(PCCsrc_, ProdInst_);

    LumiToken=consumes<LumiInfo, edm::InLumi>(PCCInputTag_);

    //producers 
    produces<LumiInfo, edm::InRun>(trigstring_);   
    produces<float, edm::InRun>( "Type1frac" );
    produces<float, edm::InRun>( "Type1res" );
    produces<float, edm::InRun>( "Type2res" );
    produces<float, edm::InRun>( "T1fUnc" );
    produces<float, edm::InRun>( "T1rUnc" );
    produces<float, edm::InRun>( "T2rUnc" );
}

//--------------------------------------------------------------------------------------------------
CorrPCCProducer::~CorrPCCProducer(){
}
//--------------------------------------------------------------------------------------------------
std::vector<float>& CorrPCCProducer::MakeCorrections(std::vector<float>& corrected_){
        
    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        corrected_.at(bx)=corrected_.at(bx)*corr_list_.at(bx);//Applying the corrections
    } 
    return corrected_;
}

//--------------------------------------------------------------------------------------------------
void  CorrPCCProducer::MakeCorrectionTemplate (){
    for(unsigned int bx=1;bx<LumiConstants::numBX;bx++){
       correctionTemplate_.at(bx)=type2_a_*exp(-(bx-1)*type2_b_);
    }

}

//--------------------------------------------------------------------------------------------------
float CorrPCCProducer::GetMaximum(std::vector<float> lumi_vector){
    float max_lumi=0;
    for(size_t i=0; i<lumi_vector.size(); i++){
        if(lumi_vector.at(i)>max_lumi)  max_lumi = lumi_vector.at(i);
    }   
    return max_lumi;
}



//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::EstimateType1Frac(std::vector<float> uncorrPCCPerBX, float &type1Frac){
    
    std::vector<float> corrected_tmp_; 
    //for(size_t i=0; i<corrected_tmp_.size(); i++){
    for(size_t i=0; i<uncorrPCCPerBX.size(); i++){
        corrected_tmp_.push_back(uncorrPCCPerBX.at(i));
    }
    bool gap=false;
    int idl=0;
    int num_cut=20;
    float noise=0;
//Finds the gap and the noise 
    for(int l=0; l<500; l++){
        if (corrected_tmp_.at(l)==0 && corrected_tmp_.at(l+1)==0 && corrected_tmp_.at(l+2)==0){
            gap=true;
        }
        if(gap && corrected_tmp_.at(l)!=0 && idl<num_cut){
//is this full flat noise subtraction?
            noise+=corrected_tmp_.at(l);
            idl+=1;
        }

    }

    if(idl!=0){
        noise=noise/idl;
    }

    else{
        noise=0;

    }
   

    //Apply initial type 1 correction
   for(size_t k=0;k<LumiConstants::numBX-1; k++){ 
       float bin_k = corrected_tmp_.at(k);
       corrected_tmp_.at(k+1)=corrected_tmp_.at(k+1)-type1Frac*bin_k;
  
   }
 
   mean_type2 = 0;  //Calculate the mean value of the type 2 residual 

   //Apply type 2 correction
   for(size_t i=0; i<LumiConstants::numBX-1; i++){
       for(size_t j=i+1; j<i+LumiConstants::numBX-1; j++){
           float bin_i = corrected_tmp_.at(i);
           if (j<LumiConstants::numBX){
               corrected_tmp_.at(j)=corrected_tmp_.at(j)-bin_i*correctionTemplate_.at(j-i);
           } 
           
           else{
               corrected_tmp_.at(j-LumiConstants::numBX) = corrected_tmp_.at(j-LumiConstants::numBX)-bin_i*correctionTemplate_.at(j-i);
           }
       }

    }
    //mean_type2=mean_type2/(LumiConstants::numBX);
  
    //Apply additional iteration for type 1 correction
    float lumiMax = GetMaximum(corrected_tmp_);
    float threshold = lumiMax*0.2;  //need to be changed to GetMaximum()*0.2
   
    mean_type1 = 0;  //Calculate the mean value of the type 1 residual 
    int nTrain = 0;
    for(size_t ibx=2; ibx<LumiConstants::numBX-5; ibx++){
        //float lumiM1 = corrected_tmp_.at(ibx-1);
        float lumi   = corrected_tmp_.at(ibx);
        float lumiP1 = corrected_tmp_.at(ibx+1);
        float lumiP2 = corrected_tmp_.at(ibx+2);
        float lumiP3 = corrected_tmp_.at(ibx+3);
        float lumiP4 = corrected_tmp_.at(ibx+4);
        float lumiP5 = corrected_tmp_.at(ibx+5);
   
        if(lumi>threshold && lumiP1<threshold && lumiP2<threshold){
            mean_type1+=(lumiP1-(lumiP3+lumiP4+lumiP5)/3)/lumi;
            mean_type2+=(lumiP3+lumiP4+lumiP5)/(3*lumi);
            nTrain+=1;
        }
    } 
    //Added condition in case nTrain is 0 and we get the nan(s)  
    if (nTrain!=0){
    mean_type1 = mean_type1/nTrain;
    mean_type2 = mean_type2/nTrain;
    }
    else{mean_type1=0;}

    type1Frac+=mean_type1;

}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::CalculateCorrections (std::vector<float> uncorrected, std::vector<float>& corr_list_, float& Overall_corr){
    std::cout<<"Making Corrections"<<std::endl;
    float type1frac = 0;
    EstimateType1Frac(uncorrected, type1frac);
    std::cout<<type1frac<<std::endl;
    EstimateType1Frac(uncorrected, type1frac);
    std::cout<<type1frac<<std::endl;

    //Find the abort gap and calculate the noise
    //std::vector<float> corrected_tmp_;
    std::vector<float> corrected_tmp_;
    for(size_t i=0; i<corrected_tmp_.size(); i++){
        corrected_tmp_.push_back(uncorrected.at(i));
    }

    //Apply initial type 1 correction
    for(size_t k=0;k<LumiConstants::numBX-1; k++){ 
       float bin_k = corrected_tmp_.at(k);
       corrected_tmp_.at(k+1)=corrected_tmp_.at(k+1)-type1frac*bin_k;
  
    }
 

   //Apply type 2 correction
   for(size_t i=0; i<LumiConstants::numBX-1; i++){
       for(size_t j=i+1; j<i+LumiConstants::numBX-1; j++){
           float bin_i = corrected_tmp_.at(i);
           if (j<LumiConstants::numBX){
               corrected_tmp_.at(j)=corrected_tmp_.at(j)-bin_i*correctionTemplate_.at(j-i);
           } 
           
           else{
               corrected_tmp_.at(j-LumiConstants::numBX) = corrected_tmp_.at(j-LumiConstants::numBX)-bin_i*correctionTemplate_.at(j-i);
           }
       }

    }
  
    //Apply additional iteration for type 1 correction
    float lumiMax = GetMaximum(corrected_tmp_);
    float threshold = lumiMax*0.2;  //need to be changed to GetMaximum()*0.2
   
    //float mean_type1 = 0;  //Calculate the mean value of the type 1 residual 
    //int nTrain = 0;
    //for(size_t ibx=2; ibx<LumiConstants::numBX-5; ibx++){
    //    //float lumiM1 = corrected_tmp_.at(ibx-1);
    //    float lumi   = corrected_tmp_.at(ibx);
    //    float lumiP1 = corrected_tmp_.at(ibx+1);
    //    float lumiP2 = corrected_tmp_.at(ibx+2);
    //    float lumiP3 = corrected_tmp_.at(ibx+3);
    //    float lumiP4 = corrected_tmp_.at(ibx+4);
    //    float lumiP5 = corrected_tmp_.at(ibx+5);
   
    //    if(lumi>threshold && lumiP1<threshold && lumiP2<threshold){
    //        mean_type1+=(lumiP1-(lumiP3+lumiP4+lumiP5)/3)/lumi;
    //        nTrain+=1;
    //    }
    //} 
  
    //mean_type1 = mean_type1/nTrain;
  
    //for (size_t ibx=0; ibx<LumiConstants::numBX-1; ibx++){
    //    float bin_i = corrected_tmp_.at(ibx);
    //    corrected_tmp_.at(ibx+1) = corrected_tmp_.at(ibx+1)-bin_i*mean_type1;

    //}
   
   float Integral_Uncorr=0;
   float Integral_Corr = 0;
   //Calculate Per-BX correction factor and overall correction factor
   for (size_t ibx=0; ibx<corrected_tmp_.size(); ibx++){
       if(corrected_tmp_.at(ibx)>threshold){
           Integral_Uncorr+=uncorrected.at(ibx);
           Integral_Corr+=corrected_tmp_.at(ibx);
       }
       
       corr_list_.at(ibx) = corrected_tmp_.at(ibx)/uncorrected.at(ibx);

   }
 
   Overall_corr = Integral_Corr/Integral_Uncorr;  
   //std::pair <float, std::vector<float>> corr_pair_ = std::make_pair(Integral_Corr/Integral_Uncorr, corrected_tmp_); 
   //return corr_pair_;
}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
    //std::cout<<"A Print Statement"<<std::endl;

      
}
//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::beginRun(edm::Run const& runSeg, const edm::EventSetup& iSetup){
    std::cout<<"Begin Run"<<std::endl;
    //LumiInfo outLumiOb; 
    LSrun_=1;
}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    std::cout<<"Begin Lumi-Block"<<std::endl;
    //outLumiOb = std::make_unique<LumiInfo>(); 
    //LumiInfo outLumiOb; 
    countLumi_++;
    
    //Lumi-Range setters 
    if (LSrun_==1){
    startLS=lumiSeg.luminosityBlock();
    LSrun_=0;
    }
   

}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    endLS = lumiSeg.luminosityBlock();
    //Try to sum over the lumisections and map to pair std::vector 
    edm::Handle<LumiInfo> PCCHandle; 
    lumiSeg.getByToken(LumiToken,PCCHandle);
    
    const LumiInfo& inLumiOb = *(PCCHandle.product()); 

    //Making the vectors to loop over the lumisections for each run    
    rawlumiBX_= inLumiOb.getInstLumiAllBX();

            
    //summing over lumisections
    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        totalLumiByBX_[bx]+=rawlumiBX_[bx];
    }
    
    //std::cout<<"Print end Lumi Block"<<std::endl;


}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endRun(edm::Run const& runSeg, const edm::EventSetup& iSetup){
    //std::cout<<"End Run"<<std::endl;
    //std::cout<<"The end time of the run: "<<std::to_string(runSeg.endTime())<<std::endl;
}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup){
    //Save if # of lumisections are reached
    //std::cout<<"The end time of the block: "<<std::to_string(lumiSeg.endTime())<<std::endl;
    //We want an object that is LS iov and the object that is saved. Somehting like
    //std::map<std::pair<int,int>, LumiInfo*>
    
    //LuminosityBlockRange lbr;
    //std::cout<<"lbr :"<<lbr.endLumi()<<std::endl;

    //Derive correction every resetNLumi 
    //if (resetNLumi_ > 0 && countLumi_%resetNLumi_!=0) return;
   
   


}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endRunProduce(edm::Run& runSeg, const edm::EventSetup& iSetup){
    std::cout<<"End Run"<<std::endl;
    //Setting the corrections
    //loop over all the pointers, identifying the iov and object. Combine the corrections?
    outLumiOb = std::make_unique<LumiInfo>(); 
    Type1frac = std::make_unique<float>(); 
    Type1res = std::make_unique<float>(); 
    Type2res = std::make_unique<float>(); 
    T1fUnc = std::make_unique<float>(); 
    T1rUnc = std::make_unique<float>(); 
    T2rUnc = std::make_unique<float>(); 
    //std::auto_ptr<float> Type1frac ( new float );
       
    //filling the map to identify the correction interval to the lumiObject 
    //std::pair<int,int> iov;//LS interval of validity 
    //iov = std::make_pair(iov1,iov2);
    
    //std::map<std::pair<int,int>, LumiInfo*> 
    //setting the map of the iov to the lumiInfo object.
    //lsmap_[iov]=&inLumiOb; 
    

    //Deriving the PCC corrections
    //CalculateCorrections(totalLumiByBX_,corr_list_, Overall_corr);

    //Setting the values for output
    outLumiOb->setInstLumi(corr_list_);   
    outLumiOb->setTotalLumi(Overall_corr);   
    outLumiOb->setStatErrorOnLumi(float(startLS));
    outLumiOb->setDeadFraction(float(endLS));
    *Type1frac=type1Frac;
    *Type1res=mean_type1;
    *Type2res=mean_type2;
    *T1fUnc=3.14159;
    *T1rUnc=3.14159;
    *T2rUnc=3.14159;

    //reset! 
    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        totalLumiByBX_[bx]=0;//reset the total after the save
    }

    runSeg.put(std::move(outLumiOb),std::string(trigstring_)); 
    runSeg.put(std::move(Type1frac),"Type1frac"); 
    runSeg.put(std::move(Type1res),"Type1res"); 
    runSeg.put(std::move(Type2res),"Type2res"); 
    runSeg.put(std::move(T1fUnc),"T1fUnc"); 
    runSeg.put(std::move(T1rUnc),"T1rUnc"); 
    runSeg.put(std::move(T2rUnc),"T2rUnc"); 
    
 


}

DEFINE_FWK_MODULE(CorrPCCProducer);
