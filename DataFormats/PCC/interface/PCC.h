#ifndef PCC_h
#define PCC_h
/** \class reco::PCC
 *  
 * Reconstructed PCC object that will contain the moduleID, BX.
 *
 * \authors: Chris Palmer capalmer@cern.ch
 * \Sam Higginbotham shiggib@cern.ch
 *
 *
 */

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

using namespace reco;
class PCC {
    public:
        PCC() : events(3564){}

        void Increment(int mD,int BXid,int count){
            std::vector<int>::iterator it;
            it = std::find(ModID.begin(), ModID.end(), mD);
            size_t modIndex = it - ModID.begin();

        if(it == ModID.end()){
                std::vector<int> empBX(3564);//make 0s
                ModID.push_back(mD);
                //events.push_back(BXid);
                counts.push_back(empBX); 
                //std::cout<<"Mod ID entry "<<mD<<" BXid "<<BXid<<" count "<<count<<std::endl;
                }
            //cout<<"Iterator :"<<*it<<endl;   
            //counts[modIndex].assign(BXid,BXid,count);
            //counts[modIndex].push_back(count);
            counts[modIndex][BXid]=counts[modIndex][BXid]+count;
            //std::replace (counts[modIndex].at(BXid),counts[modIndex].at(BXid),0,count);
            }
        ////////////////////////////////////////////
        void eventCounter(int BXid){
            events[BXid]++;//make sure this is at the right index!
}
        std::vector<std::vector<int> > read_counts(){
                return(counts);
            }
        ////////////////////////////////////////////
        void printVector()
        { 
            int irow = 0;
            std::vector< std::vector<int> >::const_iterator row; 
            std::vector<int>::const_iterator col; 

            for (row = counts.begin(); row != counts.end(); ++row)
            { 
            
               std::cout << ModID[irow] << " :";
               for (col = row->begin(); col != row->end(); ++col)
               { 
                  std::cout << *col << " "; 
               } 
               irow++;
                std::cout<<std::endl;
            } 

        }

      private:
        std::vector<std::vector<int> > counts;
        std::vector<int> events;
        std::vector<int> ModID;
        std::vector<int> empBX;

};

  ///
  //std::ostream& operator<< ( std::ostream&, PCC beam );

#endif
