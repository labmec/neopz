//---------------------------------------------------------------------------

#ifndef TPZSloanRenumberingH
#define TPZSloanRenumberingH

#include "pzrenumbering.h"
#include <map>
#include "pzmanvector.h"
#include <list>
#include "pzstack.h"

//#define _OPTSET4History_
class TPZSloanRenumbering : public TPZRenumbering

{
    
private:
    
    //Auxiliar structure: a list without duplicates and with control of history
    struct SList{
        //it stores the elements of the list
        TPZVec<long> fList;
        
        //actual size of the list
        long fSize;
        
        //history: stores every element that has been in the list ever
#ifdef _OPTSET4History_
        std::set<long> fHistory;
#else
        TPZVec<long> fHistory;
#endif
        
        SList(long nnodes){
            fList.Resize(nnodes);
#ifndef _OPTSET4History_
            fHistory.Resize(nnodes);
            fHistory.Fill(0);
#endif
            fSize = 0;
        }
        ~SList(){
            //nothing here
        }
        
        
        void Reset(){
#ifdef _OPTSET4History_
            fHistory.clear();
#else
            fHistory.Fill(0);
#endif
            fSize = 0;
        }
        
        //returns true if object is inserted or false otherwise (i.e. already existed)
        bool push_back(long val){
#ifdef _OPTSET4History_
            if( (fHistory.insert(val)).second == true){
#else
                if( fHistory[val] == 0 ){
                    fHistory[val] = 1;
#endif
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
        
        struct TNo
        {
            TNo *fprev;
            TNo *fnext;
            
            long fIndex;
            long fSequenceNumber;
            
            int fPriority;
            int fStatus;
            
            TNo() : fprev(0),fnext(0),fIndex(-1),fSequenceNumber(-1),fPriority(-1),fStatus(EInactive)
            {
                
            }
        };
        
        std::map<int, TNo *> fActive;
        
        TPZVec<TNo> fAllNodes;
        
        void InsertNode(TNo *node)
        {
            if( fActive.find(node->fPriority) == fActive.end())
            {
                node->fprev = 0;
                node->fnext = 0;
                fActive[node->fPriority] = node;
            }
            else
            {
                TNo *first = fActive[node->fPriority];
                TNo *second = first->fnext;
                first->fnext = node;
                node->fprev = first;
                node->fnext = second;
                if(second) second->fprev = node;
            }
        }
        
        TNo *PopHighestPriorityNode()
        {
#ifdef DEBUG
            if(fActive.size() == 0) DebugStop();
#endif
            TNo* result = fActive.rbegin()->second;
            int priority = fActive.rbegin()->first;
            TNo* next = result->fnext;
            if(! next)
            {
                fActive.erase(priority);
            } else{
                next->fprev = result->fprev;
                if(next->fprev) next->fprev->fnext = next;
                fActive[priority] = next;
            }
            result->fprev = 0;
            result->fnext = 0;
            return result;
        }
        
        void TransferPriority(TNo *no, int newpriority)
        {
            int priority = no->fPriority;
#ifdef DEBUG
            if(fActive.find(priority) == fActive.end()) DebugStop();
            if(no->fPriority == newpriority) DebugStop();
#endif
            TNo *first = fActive[priority];
            if(no == first)
            {
                first = no->fnext;
                if(!first )
                {
                    fActive.erase(priority);
                }
                else
                {
                    first->fprev = no->fprev;
                    if(no->fprev) no->fprev->fnext = first;
                    fActive[priority] = first;
                }
            }
            else
            {
                TNo *prev = no->fprev;
                TNo *next = no->fnext;
                if(prev) prev->fnext = next;
                if(next) next->fprev = prev;
            }
            no->fprev = 0;
            no->fnext = 0;
            no->fPriority = newpriority;
            InsertNode(no);
        }
        
    public:
        
        
        virtual void Resequence(TPZVec<long> &permGather, TPZVec<long> &permScatter);

        virtual void Resequence2(TPZVec<long> &permGather, TPZVec<long> &permScatter);

        
        TPZSloanRenumbering(long NElements, long NNodes);
        
        virtual ~TPZSloanRenumbering();
        
    };
#endif
