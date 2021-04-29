//---------------------------------------------------------------------------

#ifndef TPZSloanRenumberingH
#define TPZSloanRenumberingH

#include "TPZRenumbering.h"
#include <map>
#include "pzmanvector.h"
#include <list>
#include "pzstack.h"

//#define _OPTSET4History_
class TPZSloanRenumbering : public TPZRenumbering

{
    
private:
    
    //Auxiliary structure: a list without duplicates and with control of history
    struct SList{
        //it stores the elements of the list
        TPZVec<int64_t> fList;
        
        //actual size of the list
        int64_t fSize;
        
        //history: stores every element that has been in the list ever
#ifdef _OPTSET4History_
        std::set<int64_t> fHistory;
#else
        TPZVec<int64_t> fHistory;
#endif
        
        SList(int64_t nnodes){
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
        bool push_back(int64_t val){
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
            
            void remove(int64_t index){
                fList[index] = fList[fSize-1];
                fSize--;
            }
            
            //Find the element and remove it. Returns true if element val was found or false otherwise
            bool FindAndRemove(int64_t val){
                for(int64_t i = 0; i < fSize; i++){
                    if(fList[i] == val){
                        this->remove(i);
                        return true;
                    }
                }
                return false;
            }//method
        };
        
        
        int64_t W1() const{ return 1; }
        int64_t W2() const{ return 2; }
        
        enum ENodeStatus
        {
            EInactive = 0,
            EPreActive = 1,
            EActive = 2,
            EPostActive = 3
        };
        
        int64_t FindHighestPriority(const SList &Q,const TPZVec<int64_t> &priority, int64_t &Qindex) const;
        
        struct TNo
        {
            TNo *fprev;
            TNo *fnext;
            
            int64_t fIndex;
            int64_t fSequenceNumber;
            
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
#ifdef PZDEBUG
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
#ifdef PZDEBUG
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
        
        
        virtual void Resequence(TPZVec<int64_t> &permGather, TPZVec<int64_t> &permScatter) override;

        virtual void Resequence2(TPZVec<int64_t> &permGather, TPZVec<int64_t> &permScatter);

        
        TPZSloanRenumbering(int64_t NElements, int64_t NNodes);
        
        TPZSloanRenumbering();
        
        virtual ~TPZSloanRenumbering();
        
    };
#endif
