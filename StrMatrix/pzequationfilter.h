#ifndef PZEQUATIONFILTERHPP
#define PZEQUATIONFILTERHPP

#include "pzvec.h"
#include <set>

class TPZEquationFilter
{
public:
    TPZEquationFilter(int numeq) : fNumEq(numeq), fMinEq(0), fMaxEq(numeq), fActiveEqs(), fDestIndices()
    {
        
    }
    
    void SetMinMaxEq(int mineq, int maxeq)
    {
        if (mineq <0 || mineq>=fNumEq || maxeq < mineq || maxeq > fNumEq) {
            // unhandled case - inconsistent data
            DebugStop();
        }
        
        fMinEq = mineq;
        fMaxEq = maxeq;
        
        if (fNumEq && fDestIndices.size() == fNumEq) {
            // please implement me!
            DebugStop();
        }
    }
    
    void SetActiveEquations(TPZVec<int> &active)
    {
        if (fMinEq != 0 || fMaxEq != fNumEq) {
            // Please implement me
            DebugStop();
        }
        std::set<int> activeset;
        int neq = active.size();
        if (neq) {
            activeset.insert(&active[0], &active[neq-1]+1);
        }

        fDestIndices.Resize(fNumEq);
        fDestIndices.Fill(-1);
        int count = 0;
        fActiveEqs.Resize(activeset.size());
        for (std::set<int>::iterator it=activeset.begin(); it != activeset.end(); it++) {
            fActiveEqs[count] = *it;
            fDestIndices[*it] = count++;
        }
    }
    
    void Reset()
    {
        fMinEq = 0;
        fMaxEq = fNumEq;
        fActiveEqs.Resize(0);
        fDestIndices.Resize(0);
    }
    
    void Filter(TPZVec<long> &orig, TPZVec<long> &dest) const
    {
        if (fMinEq == 0 && fMaxEq == fNumEq && fDestIndices.size() == 0) {
            return;
        }
        if (fMinEq != 0 || fMaxEq != fNumEq) 
        {
            int count = 0;
            int numeq = dest.size();
            for (int i=0; i<numeq; i++) {
                if (dest[i] >= fMinEq && dest[i] < fMaxEq) {
                    orig[count] = orig[count];
                    dest[count] = dest[i]-fMinEq;
                    count++;
                }
            }
            orig.Resize(count);
            dest.Resize(count);
            if (fNumEq && fDestIndices.size()) {
                DebugStop();
            }
        }
        if (fDestIndices.size()) {
            if (fMinEq != 0 || fMaxEq != fNumEq) {
                DebugStop();
            }
            int count = 0;
            int numeq = dest.size();
            for (int i=0; i<numeq; i++) {
                if (fDestIndices[dest[i]] != -1) {
                    orig[count] = orig[i];
                    dest[count] = fDestIndices[dest[i]];
                    count++;
                }
            }
            orig.Resize(count);
            dest.Resize(count);
        }
    }
    
    void Filter(TPZVec<long> &dest) const
    {
        if (fMinEq == 0 && fMaxEq == fNumEq && fDestIndices.size() == 0) {
            return;
        }
        if (fMinEq != 0 || fMaxEq != fNumEq) 
        {
            int count = 0;
            int numeq = dest.size();
            for (int i=0; i<numeq; i++) {
                if (dest[i] >= fMinEq && dest[i] < fMaxEq) {
                    dest[count] = dest[i]-fMinEq;
                    count++;
                }
            }
            dest.Resize(count);
            if (fNumEq && fDestIndices.size()) {
                DebugStop();
            }
        }
        if (fDestIndices.size()) {
            if (fMinEq != 0 || fMaxEq != fNumEq) {
                DebugStop();
            }
            int count = 0;
            int numeq = dest.size();
            for (int i=0; i<numeq; i++) {
                if (fDestIndices[dest[i]] != -1) {
                    dest[count] = fDestIndices[dest[i]];
                    count++;
                }
            }
            dest.Resize(count);
        }
    }
    
    int NEq() const
    {
        if (!Consistent()) {
            DebugStop();
        }
        if (fMinEq != 0 || fMaxEq != fNumEq) {
            return fMaxEq-fMinEq;
        }
        if (fActiveEqs.size()) {
            return fActiveEqs.size();
        }
        return fNumEq;
    }
    
    bool IsActive() const
    {
        if (fMinEq == 0 && fMaxEq == fNumEq && fDestIndices.size() == 0) {
            return false;
        }
        else {
            return true;
        }

    }
    
    int NumActive(int minindex, int maxindex)
    {
        if (maxindex < minindex || !Consistent()) {
            DebugStop();
        }
        if (!IsActive()) {
            return maxindex-minindex;
        }
        if (fMinEq != 0 || fMaxEq != fNumEq) {
            if (maxindex < fMinEq) {
                maxindex = fMinEq;
            }
            if (minindex < fMinEq) {
                minindex = fMinEq;
            }
            if (minindex > fMaxEq) {
                minindex = fMaxEq;
            }
            if (maxindex > fMaxEq) {
                maxindex = fMaxEq;
            }
            return maxindex-minindex;
        }
        int numactive = 0;
        for (int i=minindex; i<maxindex; i++) {
            if (fDestIndices[i] != -1) {
                numactive++;
            }
        }
        return numactive;
    }
    
    /*!
     \fn FilterSkyline()
     */
    void FilterSkyline(TPZVec<int> &skyline)
    {
        if (!IsActive()) {
            return;
        }
        if (!Consistent()) {
            DebugStop();
        }
        if (!fActiveEqs.size()) 
        {
            //  int neq = skyline.NElements();
            int ieq;
            for(ieq = fMinEq; ieq < fMaxEq; ieq++)
            {
                skyline[ieq-fMinEq] = skyline[ieq]-fMinEq;
            }
            skyline.Resize(fMaxEq-fMinEq);
        }
        else {
            for (int ieq = 0; ieq<fActiveEqs.size(); ieq++) 
            {
                int skyl = skyline[fActiveEqs[ieq]];
                while (fDestIndices[skyl] == -1 && skyl < fNumEq) {
                    skyl++;
                }
#ifdef DEBUG
                // all active equations should have a destination
                if (skyl > fActiveEqs[ieq] || fDestIndices[skyl] < 0) {
                    DebugStop();
                }
#endif
                skyline[ieq] = fDestIndices[skyl];
            }
            skyline.resize(fActiveEqs.size());
        }
    }
    

private:
    
    int fNumEq;
    int fMinEq;
    int fMaxEq;
    
    TPZVec<int> fActiveEqs;
    TPZVec<int> fDestIndices;
    
    bool Consistent() const
    {
        bool filterone = false;
        if (fMinEq != 0 || fMaxEq != fNumEq) {
            filterone = true;
        }
        bool filtertwo = false;
        if (fNumEq && fDestIndices.size() == fNumEq) {
            filtertwo = true;
        }
        if (filterone && filtertwo) {
            return false;
        }
        return true;
    }
    
};

#endif