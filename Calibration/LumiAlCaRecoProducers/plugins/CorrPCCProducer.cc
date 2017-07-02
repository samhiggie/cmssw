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

    std::string lumistring = "LumiInfo";   
    std::string type1facstring = "Type1or2FraOrRes";   

    std::string trigstring_; //specifices the iov LSs for the object that is saved.  
    std::vector<float> rawlumiBX_;//new vector containing clusters per bxid 
    std::vector<float> errOnLumiByBX_;//standard error per bx
    std::vector<float> totalLumiByBX_;//summed lumi
    std::vector<float> correctionTemplate_;
    std::vector<float> corr_list_;//list of scale factors to apply.
    float Overall_corr;//The Overall correction to the integrated luminosity

    std::map<std::pair<int,int>, std::unique_ptr<LumiInfo>>::iterator it; 
    std::map<std::pair<int,int>, std::unique_ptr<LumiInfo>> myInfoPointers;//map to obtain iov for lumiOb corrections to the luminosity. 

    float type1frac;
    std::vector<float> t1fUncVect;
    float mean_type1;//Type 1 residual 
    float mean_type2;//Type 2 residual 
    float t1fUnc;//Type 1 fraction uncertainty rms
    float t1rUnc;//Type 1 residual uncertainty rms 
    float t2rUnc;//Type 2 residual uncertainty rms 
    int nTrain;//Number of bunch trains used in calc type 1 and 2 res, frac.
    int countLumi_;//The lumisection count... the size of the lumiblock
    int resetNLumi_;//The number of lumisections per block.
    int LSrun_;//The parameters that will save beginning and end LS for saving corr to runs
    int iov1;//beginning lumisection for iov
    int iov2;//end lumisection for iov
    int startLS;//Starting lumisection for the iov that we save with the lumiInfo object
    int thisLS;//Ending lumisection for the iov that we save with the lumiInfo object.

    //double type1_; //Initial type 1 correction factor
    double type2_a_;//amplitude for the type 2 correction 
    double type2_b_;//decay width for the type 2 correction

    
    //output 
    std::unique_ptr<LumiInfo> outLumiOb;//lumi object with corrections per BX
    std::unique_ptr<LumiInfo> thisLSLumiInfo;//lumi object with corrections per BX
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
    for ( int iLum=0; iLum<5; iLum++){
        produces<LumiInfo, edm::InRun>(lumistring+std::to_string(iLum));
        produces<float, edm::InRun>(type1facstring+std::to_string(iLum));
    }
        produces<LumiInfo, edm::InRun>(trigstring_);   
    //produces<LumiInfo, edm::InRun>(trigstring_);   
    //produces<float, edm::InRun>( "Type1frac" );
    //produces<float, edm::InRun>( "Type1res" );
    //produces<float, edm::InRun>( "Type2res" );
    //produces<float, edm::InRun>( "T1fUnc" );
    //produces<float, edm::InRun>( "T1rUnc" );
    //produces<float, edm::InRun>( "T2rUnc" );
    
    
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
void CorrPCCProducer::EstimateType1Frac(std::vector<float> uncorrPCCPerBX, float &type1frac){
    
    std::vector<float> corrected_tmp_; 
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
    //mean_type2=mean_type2/(LumiConstants::numBX);
  
    //Apply additional iteration for type 1 correction
    float lumiMax = GetMaximum(corrected_tmp_);
    float threshold = lumiMax*0.2;  //need to be changed to GetMaximum()*0.2
   
    mean_type1 = 0;  //Calculate the mean value of the type 1 residual 
    //mean_type2 = 0;  //Calculate the mean value of the type 2 residual 
    nTrain = 0;
    for(size_t ibx=2; ibx<LumiConstants::numBX-5; ibx++){
        //float lumiM1 = corrected_tmp_.at(ibx-1);
        float lumi   = corrected_tmp_.at(ibx);
        float lumiP1 = corrected_tmp_.at(ibx+1);
        float lumiP2 = corrected_tmp_.at(ibx+2);
        float lumiP3 = corrected_tmp_.at(ibx+3);
        float lumiP4 = corrected_tmp_.at(ibx+4);
        float lumiP5 = corrected_tmp_.at(ibx+5);
   
        //Where type 1 and type 2 residuals are computed
        if(lumi>threshold && lumiP1<threshold && lumiP2<threshold){
            mean_type1+=(lumiP1-(lumiP3+lumiP4+lumiP5)/3)/lumi;
            //mean_type2+=(lumiP3+lumiP4+lumiP5)/(3*lumi);
            nTrain+=1;
        }
    } 
    //Added condition in case nTrain is 0 and we get the nan(s)  
    if (nTrain!=0){
    mean_type1 = mean_type1/nTrain;
    }
    else{mean_type1=0;}

    type1frac+=mean_type1;
    //t1fUncVect.push_back(mean_type1);
}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::CalculateCorrections (std::vector<float> uncorrected, std::vector<float>& corr_list_, float& Overall_corr){
    std::cout<<"Making Corrections"<<std::endl;
    type1frac = 0;
    EstimateType1Frac(uncorrected, type1frac);
    std::cout<<type1frac<<std::endl;
    EstimateType1Frac(uncorrected, type1frac);
    std::cout<<type1frac<<std::endl;

    //Find the abort gap and calculate the noise
    std::vector<float> corrected_tmp_;
    for(size_t i=0; i<uncorrected.size(); i++){
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

   float lumiMax = GetMaximum(corrected_tmp_);
   float threshold = lumiMax*0.2;  //need to be changed to GetMaximum()*0.2

    //change naming convention here.... call them residuals
    mean_type1 = 0;  //Calculate the mean value of the type 1 residual 
    mean_type2 = 0;  //Calculate the mean value of the type 2 residual
    float mean_sq_type1 = 0; 
    nTrain = 0;
    for(size_t ibx=2; ibx<LumiConstants::numBX-5; ibx++){
            //float lumiM1 = corrected_tmp_.at(ibx-1);
            float lumi   = corrected_tmp_.at(ibx);
            float lumiP1 = corrected_tmp_.at(ibx+1);
            float lumiP2 = corrected_tmp_.at(ibx+2);
            float lumiP3 = corrected_tmp_.at(ibx+3);
            float lumiP4 = corrected_tmp_.at(ibx+4);
            float lumiP5 = corrected_tmp_.at(ibx+5);
       
            //Where type 1 and type 2 residuals are computed
            if(lumi>threshold && lumiP1<threshold && lumiP2<threshold){
                mean_type1+=(lumiP1-(lumiP3+lumiP4+lumiP5)/3)/lumi;
                mean_sq_type1+=TMath::Power((lumiP1-(lumiP3+lumiP4+lumiP5)/3)/lumi, 2);
                mean_type2+=(lumiP3+lumiP4+lumiP5)/(3*lumi);
                //mean_type2+=TMath::Power((lumiP3+lumiP4+lumiP5)/(3*lumi)i,2);
                nTrain+=1;
            }
        } 
    //Added condition in case nTrain is 0 and we get the nan(s)  
    if (nTrain!=0){
        mean_type1 = mean_type1/nTrain;
        mean_type2 = mean_type2/nTrain;
        mean_sq_type1 = mean_sq_type1/nTrain;
        t1rUnc = sqrt(mean_sq_type1-mean_type1*mean_type1)/sqrt(nTrain);
        t2rUnc =TMath::Sqrt(t2rUnc/nTrain);
        t2rUnc = 1/sqrt(3*nTrain); 
    }
    else{mean_type1=0;}

    //type1frac+=mean_type1;
    t1fUncVect.push_back(mean_type1);

           
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

    //beginning of lumiblock
    iov2 = lumiSeg.luminosityBlock();
        
    //Lumi-Range setters 
    if (LSrun_==1){
    startLS=lumiSeg.luminosityBlock();
    LSrun_=0;
    }
   

}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    thisLS = lumiSeg.luminosityBlock();
    iov1 = lumiSeg.luminosityBlock();
    //Try to sum over the lumisections and map to pair std::vector 
    edm::Handle<LumiInfo> PCCHandle; 
    lumiSeg.getByToken(LumiToken,PCCHandle);
    
    const LumiInfo& inLumiOb = *(PCCHandle.product()); 

    bool found=false;//Going to match iovs to lumiInfo objects


    //Making the vectors to loop over the lumisections for each run    
    rawlumiBX_= inLumiOb.getInstLumiAllBX();

        
    //filling the map to identify the correction interval to the lumiObject 
    //std::pair<int,int> iov;//LS interval of validity 
    //iov = std::make_pair(iov1,iov2);
    
    //setting the map of the iov to the lumiInfo object.
    //myInfoPointers[iov] = &inLumiOb; 
    
    //const LumiInfo* thisLSLumiInfo; //lumioutput object but will be used with map and saved to Run
    //thisLSLumiInfo = std::make_unique<LumiInfo>(); 
    //thisLSLumiInfo = std::unique_ptr<LumiInfo>(); 


    for(it=myInfoPointers.begin(); (it!=myInfoPointers.end()); ++it) {
        if( (thisLS >= it->first.first) && (thisLS<= it->first.second) ){
          thisLSLumiInfo = it->second;
          found=true;
          break;
        }
    }

    if(!found) {
        //insert new pair
        std::pair<int,int> lsKey;
        //find appropriate window
        // this is a function of the nLS block length and thisLS
        // e.g. if thisLS=56 and nLS=50, then lsKey=(51,100)
        lsKey=std::make_pair(thisLS - thisLS%resetNLumi_+1,thisLS - thisLS%resetNLumi_+resetNLumi_);
        //myInfoPointers[lsKey]= new LumiInfo();
        myInfoPointers[lsKey]= std::unique_ptr<LumiInfo>();
        thisLSLumiInfo = myInfoPointers[lsKey];
    }


    //Deriving the PCC corrections every N lumi
    if(resetNLumi_!=0 && countLumi_%resetNLumi_==0){
        CalculateCorrections(totalLumiByBX_,corr_list_, Overall_corr);
    }
           
    //summing over lumisections
    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        //thisLSLumiInfo->instantLumiByBX_.at(bx)+=rawlumiBX_[bx];
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
    //Setting the values for output
    outLumiOb->setInstLumi(corr_list_);   
    outLumiOb->setTotalLumi(Overall_corr);   
    outLumiOb->setStatErrorOnLumi(float(startLS));
    outLumiOb->setDeadFraction(float(thisLS));
    *Type1frac=type1frac;
    *Type1res=mean_type1;
    *Type2res=mean_type2;
    t1fUnc=(TMath::RMS(t1fUncVect.begin(),t1fUncVect.end()))/TMath::Sqrt(nTrain);
    *T1fUnc=t1fUnc;
    *T1rUnc=t1rUnc;
    *T2rUnc=t2rUnc;

    //reset! 
    for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
        totalLumiByBX_[bx]=0;//reset the total after the save
    }
    
    runSeg.put(std::move(outLumiOb),std::string(lumistring+std::to_string(1))); 
    runSeg.put(std::move(Type1frac),type1facstring+std::to_string(0)); 
    runSeg.put(std::move(Type1res),type1facstring+std::to_string(1)); 
    runSeg.put(std::move(Type2res),type1facstring+std::to_string(2)); 
    runSeg.put(std::move(T1fUnc),type1facstring+std::to_string(3)); 
    runSeg.put(std::move(T1rUnc),type1facstring+std::to_string(4)); 
    

    for(it=myInfoPointers.begin(); (it!=myInfoPointers.end()); ++it) {
        runSeg.put(std::move(it->second),std::string(trigstring_)); 

    }


    //runSeg.put(std::move(outLumiOb),std::string(trigstring_)); 
    //runSeg.put(std::move(Type1frac),"Type1frac"); 
    //runSeg.put(std::move(Type1res),"Type1res"); 
    //runSeg.put(std::move(Type2res),"Type2res"); 
    //runSeg.put(std::move(T1fUnc),"T1fUnc"); 
    //runSeg.put(std::move(T1rUnc),"T1rUnc"); 
    //runSeg.put(std::move(T2rUnc),"T2rUnc"); 
    
 

}

DEFINE_FWK_MODULE(CorrPCCProducer);
