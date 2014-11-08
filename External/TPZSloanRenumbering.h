//---------------------------------------------------------------------------

#ifndef TPZSloanRenumberingH
#define TPZSloanRenumberingH

#include "pzrenumbering.h"
#include <map>
#include "pzmanvector.h"
#include <list>
#include "pzstack.h"

class TPZSloanRenumbering : public TPZRenumbering {

  private:

    //Auxiliar structure: a list without duplicates and with control of history
    struct SList{
      //it stores the elements of the list
      TPZVec<long> fList;

      //actual size of the list
      long fSize;

      //history: stores every element that has been in the list ever
      std::set<long> fHistory;

      SList(long nnodes){
        fList.Resize(nnodes);
        //fHistory.Resize(nnodes);
        //fHistory.Fill(0);
        fSize = 0;
      }
      ~SList(){
        //nothing here
      }


      void Reset(){
        fHistory.clear();
        fSize = 0;
      }

      //returns true if object is inserted or false otherwise (i.e. already existed)
      bool push_back(long val){
        if( (fHistory.insert(val)).second == true){
          fList[fSize] = val;
          fSize++;
          return true;
        }
        else return false;
      }

      void remove(long index){
        fList[index] = fList[fSize-1];
        fSize--;
      }

      //Find the element and remove it. Returns true if element val was found or false otherwise
      bool FindAndRemove(long val){
        for(long i = 0; i < fSize; i++){
          if(fList[i] == val){
            this->remove(i);
            return true;
          }
        }
        return false;
      }//method
    };


    long W1() const{ return 1; }
    long W2() const{ return 2; }

    enum ENodeStatus
    {
      EInactive = 0,
      EPreActive = 1,
      EActive = 2,
      EPostActive = 3
    };

    long FindHighestPriority(const SList &Q,const TPZVec<long> &priority, long &Qindex) const;



  public:


    virtual void Resequence(TPZVec<long> &permGather, TPZVec<long> &permScatter);

    TPZSloanRenumbering(long NElements, long NNodes);

    virtual ~TPZSloanRenumbering();

};
#endif
