/**_________________________________________________________________
class:   CorrPCCProducer.cc
description: Computes the type1 corrections to the luminosity
authors:Sam Higginbotham (shigginb@cern.ch) and Chris Palmer (capalmer@cern.ch) 
Notes:
1.The producers must be initiated in the Default constructor
2.Must use a std::unique_ptr to save objects to an EDM Producer submethod 
3.using .put(std::move(object),std::to_string(prodInst)):
    a.Cannot save 2 objects to the same product instance
    b.Cannot save same object to different product instance 
Solutions:
1.Save the map<std::pair<int>,*LumiInfo> to Run along with the object. Then in next producer then match the LS range to the object in the runs, gathering the objects at the lumiseciton level. ->Downside could be needing a master normtag-like file with all valid LS.
2.Save at the LS 
________________________________________________________________**/

// C++ standard
#include <memory>
#include <string>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <iostream> 
#include <map>
#include <utility>
// CMS
#include "DataFormats/Luminosity/interface/PixelClusterCounts.h"
#include "DataFormats/Luminosity/interface/LumiInfo.h"
#include "DataFormats/Luminosity/interface/LumiConstants.h"
//Test for intellegent lumiblocks
#include "DataFormats/Luminosity/interface/LumiConstants.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/FileBlock.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/IOVSyncValue.h"
// CMS
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "TMath.h"
#include "TH1.h"
#include "TFile.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/MyLumiCorrections/interface/MyLumiCorrections.h"
#include "CondFormats/DataRecord/interface/MyLumiCorrectionsRcd.h"
#include "CondFormats/Serialization/interface/Serializable.h"

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
    std::string type1facString = "Type1fraction";   
    std::string type1resString = "Type1residual";   
    std::string type2resString = "Type2residual"; 

    std::string trigstring_; //specifices the iov LSs for the object that is saved.  
    std::vector<float> rawlumiBX_;//new vector containing clusters per bxid 
    std::vector<float> errOnLumiByBX_;//standard error per bx
    std::vector<float> totalLumiByBX_;//summed lumi
    std::vector<float> correctionTemplate_;
    std::vector<float> corr_list_;//list of scale factors to apply.
    float Overall_corr;//The Overall correction to the integrated luminosity

    int totalLS=3000;//Max number of Lumisections in a run! Change this later to something more robust
    //mergedLS=100;
    float nBlocks=0;
    //std::vector<std::pair<int,int>> iovs;//move this to global later 
    std::map<std::pair<int,int>, LumiInfo*>::iterator it; 
    std::map<std::pair<int,int>, LumiInfo*> myInfoPointers;//map to obtain iov for lumiOb corrections to the luminosity. 
    //The Penultimate solution....
    //std::vector<std::pair<std::pair<int,int>,LumiInfo*>> myInfoVector;//Will save this object to the run

    //Why not some histos at this point?!
    TH1D  *clustersBxHist;
    TH1D  *corrLumiHist;
    TH1F *myhist;
    TH1D *type1fracHist;
    TH1D *type1resHist;
    TH1D *type2resHist;
    TList *hlist;//list for the clusters and corrections 

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
    int iov1;//beginning lumisection for iov
    int iov2;//end lumisection for iov
    int thisLS;//Ending lumisection for the iov that we save with the lumiInfo object.

    //double type1_; //Initial type 1 correction factor
    double type2_a_;//amplitude for the type 2 correction 
    double type2_b_;//decay width for the type 2 correction

    
    //output 
    std::unique_ptr<LumiInfo> outLumiOb;//lumi object with corrections per BX
    //std::unique_ptr<LumiInfo> thisLSLumiInfo;//lumi object with corrections per BX
    std::unique_ptr<float> Type1frac;
    std::unique_ptr<float> Type1res;
    std::unique_ptr<float> Type2res;
    std::unique_ptr<float> T1fUnc;
    std::unique_ptr<float> T1rUnc;
    std::unique_ptr<float> T2rUnc;
    
    MyLumiCorrections* pMyLumiCorrections;
    

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

    //Making the template for the corrections 
    MakeCorrectionTemplate(); 

    //Input tag for raw lumi
    edm::InputTag PCCInputTag_(PCCsrc_, ProdInst_);

    LumiToken=consumes<LumiInfo, edm::InLumi>(PCCInputTag_);
    

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
void  CorrPCCProducer::MakeCorrectionTemplate(){
    for(unsigned int bx=1;bx<LumiConstants::numBX;bx++){
       correctionTemplate_.at(bx)=type2_a_*exp(-(float(bx)-1)*type2_b_);
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


    //Make Histos for uncorrected and uncrorrected*corr_list_


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

    t1fUncVect.push_back(mean_type1);

           
   float Integral_Uncorr=0;
   float Integral_Corr = 0;
   //Calculate Per-BX correction factor and overall correction factor
   for (size_t ibx=0; ibx<corrected_tmp_.size(); ibx++){
       if(corrected_tmp_.at(ibx)>threshold){
           Integral_Uncorr+=uncorrected.at(ibx);
           Integral_Corr+=corrected_tmp_.at(ibx);
       }
       if(corrected_tmp_.at(ibx)!=0.0&&uncorrected.at(ibx)!=0.0){ 
       corr_list_.at(ibx) = corrected_tmp_.at(ibx)/uncorrected.at(ibx);
       }
       else{corr_list_.at(ibx) = 0.0;}
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
}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){
    std::cout<<"Begin Lumi-Block"<<std::endl;
    countLumi_++;

    //beginning of lumiblock
    iov1 = lumiSeg.luminosityBlock();
        
    //Lumi-Range setters 
    if(resetNLumi_!=0&&countLumi_%resetNLumi_!=0){return;}
    nBlocks++;

}


//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, const edm::EventSetup& iSetup){

    thisLS = lumiSeg.luminosityBlock();

    //Try to sum over the lumisections and map to pair std::vector 
    edm::Handle<LumiInfo> PCCHandle; 
    lumiSeg.getByToken(LumiToken,PCCHandle);
    
    const LumiInfo& inLumiOb = *(PCCHandle.product()); 

    bool found=false;//Going to match iovs to lumiInfo objects


    //Making the vectors to loop over the lumisections for each run    
    rawlumiBX_= inLumiOb.getInstLumiAllBX();

        
    LumiInfo* thisLSLumiInfo; //lumioutput object but will be used with map and saved to Run
    std::memcpy(&thisLSLumiInfo,&inLumiOb,sizeof(thisLSLumiInfo));

    for(it=myInfoPointers.begin(); (it!=myInfoPointers.end()); ++it) {
        if( (thisLS >= it->first.first) && (thisLS<= it->first.second) ){
            thisLSLumiInfo = it->second;
            found=true;
            //summing over lumisections
            for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
                totalLumiByBX_[bx]+=rawlumiBX_[bx];
            }
            thisLSLumiInfo->setInstLumi(totalLumiByBX_);
            //std::cout<<"Total Lumi By BX "<<totalLumiByBX_.at(200)<<std::endl;
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
        std::cout<<"Found a new range and object "<<lsKey.first<<":"<<lsKey.second<<std::endl;
        myInfoPointers[lsKey]= new LumiInfo();
        thisLSLumiInfo = myInfoPointers[lsKey];
        for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
                totalLumiByBX_[bx]=0.0;
        }
    }
}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endRun(edm::Run const& runSeg, const edm::EventSetup& iSetup){
    //std::cout<<"End Run"<<std::endl;
    //std::cout<<"The end time of the run: "<<std::to_string(runSeg.endTime())<<std::endl;
 
}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumiSeg, const edm::EventSetup& iSetup){


}

//--------------------------------------------------------------------------------------------------
void CorrPCCProducer::endRunProduce(edm::Run& runSeg, const edm::EventSetup& iSetup){
    std::cout<<"End Run"<<std::endl;
    //Setting the corrections
    outLumiOb = std::make_unique<LumiInfo>(); 

    //Setting up database parameters
    edm::Service<cond::service::PoolDBOutputService> poolDbService;
    cond::Time_t thisIOV = 1;


    char *filename = new char[50];
    sprintf(filename,"Hist_Run_%d.root",runSeg.run());
    TFile * outf = new TFile(filename, "RECREATE");
    outf->cd();
    TH1F *myhist[int(nBlocks)+1];
    char *histname = new char[50];

    type1fracHist = new TH1D("Type 1 Fraction","Type 1 Fraction",1000,-0.5,0.5);
    type1resHist = new TH1D("Type 1 Residual","Type 1 Residual",4000,-0.2,0.2);
    type2resHist = new TH1D("Type 2 Residual","Type 2 Residual",4000,-0.2,0.2);

    //Setting the data in the run but in a lumisection range
    for(it=myInfoPointers.begin(); (it!=myInfoPointers.end()); ++it) {

        totalLumiByBX_=it->second->getInstLumiAllBX();

        CalculateCorrections(totalLumiByBX_,corr_list_, Overall_corr); 

        //std::vector<float> *pointer;
        //float *tp1;
        //float *tp2;
        //float *tp3;
        //float *tp4;

        //pointer = &corr_list_;
        //tp1 = &type1frac;
        //tp2 = &mean_type1;
        //tp3 = &mean_type2;
        //tp4 = &Overall_corr;

        edm::LuminosityBlockID lu(runSeg.id().run(),edm::LuminosityBlockNumber_t (it->first.first));
        thisIOV = (cond::Time_t)(lu.value()); 
        std::cout<<"This IOV "<<thisIOV<<std::endl;
        //thisIOV = (cond::Time_t)(edm::LuminosityBlockNumber_t (it->first.first));

        //Writing the corrections to SQL lite file for db. 
        MyLumiCorrections* pMyLumiCorrections = new MyLumiCorrections();
        pMyLumiCorrections->SetOverallCorrection(Overall_corr);
        pMyLumiCorrections->SetType1Fraction(type1frac);
        pMyLumiCorrections->SetType1Residual(mean_type1);
        pMyLumiCorrections->SetType2Residual(mean_type2);
        



        poolDbService->writeOne<MyLumiCorrections>(pMyLumiCorrections,thisIOV,"MyLumiCorrectionsRcd"); 
        //poolDbService->writeOne<std::vector<float>>(pointer,thisIOV,"Corrections"); 
        //poolDbService->writeOne<float>(tp1,thisIOV,"Type1Frac"); 
        //poolDbService->writeOne<float>(tp2,thisIOV,"Type1Residual"); 
        //poolDbService->writeOne<float>(tp3,thisIOV,"Type2Residual"); 
        //poolDbService->writeOne<float>(tp4,thisIOV,"OverallCorr"); 

         
        //histos
        int block = (it->first.first-1)*nBlocks/totalLS;  //iBlock*totalLS/nBlocks+1
        sprintf(histname, "CorrectedLumi_%d_%d_%d",block,it->first.first,it->first.second);
        std::cout<<"Histogram Name "<<histname<<std::endl;
        myhist[block]=new TH1F(histname,"",LumiConstants::numBX,1,LumiConstants::numBX);
        for(unsigned int bx=0;bx<LumiConstants::numBX;bx++){
            Double_t value = totalLumiByBX_[bx]*corr_list_[bx];
            myhist[block]->SetBinContent(bx,value);
            //std::cout<<"Bx Number"<<bx<<"Value "<<value<<std::endl;
        }
        myhist[block]->Write(); 
        type1fracHist->Fill(type1frac);
        type1resHist->Fill(mean_type1);
        type2resHist->Fill(mean_type2);

        //reset just in case
        type1frac=0.0;
        mean_type1=0.0;
        mean_type2=0.0;
        
    }
     
     type1fracHist->Write();
     type1resHist->Write();
     type2resHist->Write();
     outf->Write();
    //reset! 
    myInfoPointers.clear();


 

}

DEFINE_FWK_MODULE(CorrPCCProducer);
