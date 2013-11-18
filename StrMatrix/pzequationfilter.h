#ifndef PZEQUATIONFILTERHPP
#define PZEQUATIONFILTERHPP

#include "pzvec.h"
#include "pzfmatrix.h"

#include <set>

class TPZEquationFilter
{
public:

    TPZEquationFilter(long numeq) : fNumEq(numeq), fActiveEqs(), fDestIndices()
    {

    }

    TPZEquationFilter(const TPZEquationFilter&cp):fNumEq(cp.fNumEq),
           fActiveEqs(cp.fActiveEqs),fDestIndices(cp.fDestIndices)
    {
      ///nothing here
    }

    ~TPZEquationFilter()
    {
      ///nothing here
    }

    TPZEquationFilter & operator=(const TPZEquationFilter&cp)
    {
        this->fNumEq = cp.fNumEq;
        this->fActiveEqs = cp.fActiveEqs;
        this->fDestIndices = cp.fDestIndices;
        return *this;
    }

    ///Define as equacoes ativas de [mineq, maxeq)
    void SetMinMaxEq(long mineq, long maxeq)
    {
      if (mineq < 0 || mineq > fNumEq ||
          maxeq < 0 || maxeq > fNumEq ||
          mineq > maxeq){
          DebugStop();
      }

      const long n = maxeq-mineq;
      TPZVec<long> activeEquations(n);
      for(long i = 0; i < n; i++){
        activeEquations[i] = i + mineq;
      }
      this->SetActiveEquations( activeEquations );
    }

    ///Define as equacoes ativas
    void SetActiveEquations(TPZVec<long> &active)
    {
        if(fActiveEqs.NElements()) DebugStop();///oops, call reset first

        ///removendo duplicados e reordenando
        std::set<long> activeset;
        long neq = active.size();
        if (neq) {
            activeset.insert(&active[0], &active[neq-1]+1);
        }

        fDestIndices.Resize(fNumEq);
        fDestIndices.Fill(-1);
        long count = 0;
        fActiveEqs.Resize(activeset.size());
        for (std::set<long>::iterator it=activeset.begin(); it != activeset.end(); it++) {
            fActiveEqs[count] = *it;
            fDestIndices[*it] = count++;
        }
    }

    /// Reset method
    void Reset()
    {
        fActiveEqs.Resize(0);
        fDestIndices.Resize(0);
    }

    /** Filtra as equações:
     * @param orig [in][out] - remove de orig equacoes nao ativas
     * @param dest [in][out] - remove de dest as equcoes nao ativas
     */
    void Filter(TPZVec<long> &orig, TPZVec<long> &dest) const
    {
        if (fDestIndices.size() == 0) {
            return;
        }
        else {
            long count = 0;
            long numeq = dest.size();
            for (long i=0; i<numeq; i++) {
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

    /** Filtra as equações:
      * @param dest [in][out] - remove de dest as equacoes nao ativas
      */
    void Filter(TPZVec<long> &dest) const
    {
        if (fDestIndices.size() == 0) {
            return;
        }
        else{
            long count = 0;
            long numeq = dest.size();
            for (long i=0; i<numeq; i++) {
                if (fDestIndices[dest[i]] != -1) {
                    dest[count] = fDestIndices[dest[i]];
                    count++;
                }
            }
            dest.Resize(count);
        }
    }

    ///Retorna o numero de equacoes ativas do sistema
    long NActiveEquations() const
    {
        if (IsActive()) {
            return fActiveEqs.size();
        }
        else
        {
            return fNumEq;
        }
    }

    ///Retorna o numero de equacoes do sistema original
    long NEqExpand() const
    {
        return fNumEq;
    }

	/**
	 * Returns true if the filter is active
	 */
    bool IsActive() const
    {
        if (fActiveEqs.size() == fNumEq || fActiveEqs.size() == 0) {
            return false;
        }
        else {
            return true;
        }

    }

	/**
	 * Expands the vector small to a original system, fill zeros into the no active equations.
	 */
	template<class TVar>
    void Scatter(const TPZFMatrix<TVar> &vsmall, TPZFMatrix<TVar> &vexpand) const
    {
        long neqcondense = this->NActiveEquations();
        if(vsmall.Rows() != neqcondense || vexpand.Rows() != fNumEq)
        {
            DebugStop();
        }
        if(! IsActive())
        {
            vexpand = vsmall;
            return;
        }
        vexpand.Zero();

#ifdef DEBUG
        {
            for(long i=0; i<neqcondense; i++)
            {
                if(fActiveEqs[i] >= fNumEq)
                {
                    DebugStop();
                }
            }
        }
#endif
        for(long i=0; i<neqcondense; i++) vexpand(fActiveEqs[i],0) = vsmall.GetVal(i,0);
    }

    /**
     * @brief Reduce the vector to the number of active equations.
     */
    template<class T>
    void Gather(const TPZFMatrix<T> &large, TPZFMatrix<T> &gathered) const
    {
        long neqcondense = this->NActiveEquations();
        if(gathered.Rows() != neqcondense || large.Rows() != fNumEq)
        {
            DebugStop();
        }
        if(! IsActive())
        {
            gathered = large;
            return;
        }
        gathered.Zero();
        for(long i=0; i<neqcondense; i++) gathered(i,0) = large.GetVal(fActiveEqs[i],0);
    }

    /**
     * @brief Returns the number of active equations between [minindex,maxindex]
     */
    long NumActive(long minindex, long maxindex) const
    {
        if (minindex < 0 || maxindex < 0 || minindex > fNumEq || maxindex > fNumEq ||
            maxindex < minindex) {
            DebugStop();
        }
        if (!IsActive()) {
            return maxindex-minindex;
        }
        int numactive = 0;
        for (long i=minindex; i<maxindex; i++) {
            if (fDestIndices[i] != -1) {
                numactive++;
            }
        }
        return numactive;
    }
    
    void FilterSkyline(TPZVec<long> &skyline) const
    {
        if (!IsActive()) {
            return;
        }

        for (long ieq = 0; ieq<fActiveEqs.size(); ieq++)
        {
            long skyl = skyline[fActiveEqs[ieq]];
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
        skyline.Resize(fActiveEqs.size());

    }


private:

    /// Numero de equacoes do sistema original
    long fNumEq;
    
    /// Equacoes ativas
    TPZVec<long> fActiveEqs;

    /// Posicao das equacoes originais no sistema reduzido
    TPZVec<long> fDestIndices;
    
};

#endif