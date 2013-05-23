#ifndef PZEQUATIONFILTERHPP
#define PZEQUATIONFILTERHPP

#include "pzvec.h"
#include <set>

class TPZEquationFilter
{
public:

    TPZEquationFilter(int numeq) : fNumEq(numeq), fNumActive(numeq), fActiveEqs(), fDestIndices()
    {

    }

    TPZEquationFilter(const TPZEquationFilter&cp):fNumEq(cp.fNumEq),fNumActive(cp.fNumActive),
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
        this->fNumActive = cp.fNumActive;
        this->fActiveEqs = cp.fActiveEqs;
        this->fDestIndices = cp.fDestIndices;
        return *this;
    }

    ///Define as equacoes ativas de [mineq, maxeq)
    void SetMinMaxEq(int mineq, int maxeq)
    {
      if (mineq < 0 || mineq > fNumEq ||
          maxeq < 0 || maxeq > fNumEq ||
          mineq > maxeq){
          DebugStop();
      }

      const int n = maxeq-mineq;
      TPZVec<int> activeEquations(n);
      for(int i = 0; i < n; i++){
        activeEquations[i] = i + mineq;
      }
      this->SetActiveEquations( activeEquations );
        fNumActive = activeEquations.size();
    }

    ///Define as equacoes ativas
    void SetActiveEquations(TPZVec<int> &active)
    {
        if(fActiveEqs.NElements()) DebugStop();///oops, call reset first

        ///removendo duplicados e reordenando
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

    /// Reset method
    void Reset()
    {
        fNumActive = fNumEq;
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

    /** Filtra as equações:
      * @param dest [in][out] - remove de dest as equacoes nao ativas
      */
    void Filter(TPZVec<long> &dest) const
    {
        if (fDestIndices.size() == 0) {
            return;
        }
        else{
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

    ///Retorna o numero de equacoes ativas do sistema
    int NActiveEquations() const
    {
        return fNumActive;
    }

    ///Retorna o numero de equacoes do sistema original
    int NEqExpand() const
    {
        return fNumEq;
    }

    ///Retorna se o filtro esta ativo
    bool IsActive() const
    {
        if (fNumActive == fNumEq) {
            return false;
        }
        else {
            return true;
        }

    }

    ///Expande o vetor para o sistema original, zerando posicao de equacoes nao ativas
    template<class T>
    void Scatter(const TPZFMatrix<T> &small, TPZFMatrix<T> &expand) const
    {
        int neqcondense = this->NActiveEquations();
        if(small.Rows() != neqcondense || expand.Rows() != fNumEq)
        {
            DebugStop();
        }
        if(! IsActive())
        {
            expand = small;
            return;
        }
        expand.Zero();


#ifdef DEBUG
        {
            for(int i=0; i<neqcondense; i++)
            {
                if(fActiveEqs[i] >= fNumEq)
                {
                    DebugStop();
                }
            }
        }
#endif
        for(int i=0; i<neqcondense; i++) expand(fActiveEqs[i],0) = small.GetVal(i,0);

    }

    /**
     * @brief Reduce the vector to the number of active equations.
     */
    template<class T>
    void Gather(const TPZFMatrix<T> &large, TPZFMatrix<T> &gathered) const
    {
        int neqcondense = this->NActiveEquations();
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
        for(int i=0; i<neqcondense; i++) gathered(i,0) = large.GetVal(fActiveEqs[i],0);
    }

    /**
     * @brief Returns the number of active equations between [minindex,maxindex]
     */
    int NumActive(int minindex, int maxindex) const
    {
        if (minindex < 0 || maxindex < 0 || minindex >= fNumEq || maxindex >= fNumEq ||
            maxindex < minindex) {
            DebugStop();
        }
        if (!IsActive()) {
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
    void FilterSkyline(TPZVec<int> &skyline) const
    {
        if (!IsActive()) {
            return;
        }

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
        skyline.Resize(fActiveEqs.size());

    }


private:

    /// Numero de equacoes do sistema original
    int fNumEq;
    
    /// Numero de equacoes ativas
    int fNumActive;

    /// Equacoes ativas
    TPZVec<int> fActiveEqs;

    /// Posicao das equacoes originais no sistema reduzido
    TPZVec<int> fDestIndices;
    
};

#endif