#ifndef PZEQUATIONFILTERHPP
#define PZEQUATIONFILTERHPP

#include "pzvec.h"
#include "pzfmatrix.h"

#include <set>

class TPZEquationFilter : public TPZSavable {
public:

    TPZEquationFilter() {} 
    
    TPZEquationFilter(int64_t numeq) : fNumEq(numeq), fIsActive(false), fActiveEqs(), fDestIndices()
    {

    }

    TPZEquationFilter(const TPZEquationFilter&cp):fNumEq(cp.fNumEq),fIsActive(cp.fIsActive),
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
        this->fIsActive = cp.fIsActive;
        this->fActiveEqs = cp.fActiveEqs;
        this->fDestIndices = cp.fDestIndices;
        return *this;
    }

    int ClassId() const override{
        return Hash("TPZEquationFilter");
    }
    
    void Read(TPZStream& buf, void* context) override {
        buf.Read(&fNumEq);
        buf.Read(fIsActive);
        buf.Read(fActiveEqs);
        buf.Read(fDestIndices);
    }
    
    void Write(TPZStream &buf, int withclassid) const override{
        buf.Write(&fNumEq);
        buf.Write(fIsActive);
        buf.Write(fActiveEqs);
        buf.Write(fDestIndices);
    }

    ///Define as equacoes ativas de [mineq, maxeq)
    void SetMinMaxEq(int64_t mineq, int64_t maxeq)
    {
      if (mineq < 0 || mineq > fNumEq ||
          maxeq < 0 || maxeq > fNumEq ||
          mineq > maxeq){
          DebugStop();
      }

      const int64_t n = maxeq-mineq;
      TPZVec<int64_t> activeEquations(n);
      for(int64_t i = 0; i < n; i++){
        activeEquations[i] = i + mineq;
      }
      this->SetActiveEquations( activeEquations );
    }

    ///Define as equacoes ativas
    void SetActiveEquations(TPZVec<int64_t> &active)
    {
        if(fActiveEqs.NElements()) DebugStop();///oops, call reset first

        fIsActive = true;
        ///removendo duplicados e reordenando
        std::set<int64_t> activeset;
        int64_t neq = active.size();
        if (neq) {
            activeset.insert(&active[0], &active[neq-1]+1);
        }
#ifdef PZDEBUG       
        if (activeset.size() > fNumEq){
            std::cout << "active set with size = " << activeset.size() << 
                         " should be greater than the number of equations " << fNumEq << std::endl;
            DebugStop();
        }
        if (*activeset.rbegin() > fNumEq){
            std::cout << "active set rbegin = " << *activeset.rbegin() << 
                         " cannot be greater than the number of equations " << fNumEq << std::endl;
            DebugStop();
        }
#endif
        fDestIndices.Resize(fNumEq);
        fDestIndices.Fill(-1);
        int64_t count = 0;
        fActiveEqs.Resize(activeset.size());
        for (std::set<int64_t>::iterator it=activeset.begin(); it != activeset.end(); it++) {
            fActiveEqs[count] = *it;
            fDestIndices[*it] = count++;
        }
    }

    /// Reset method
    void Reset()
    {
        fIsActive = false;
        fActiveEqs.Resize(0);
        fDestIndices.Resize(0);
    }

    /** Filtra as equações:
     * @param orig [in][out] - remove de orig equacoes nao ativas
     * @param dest [in][out] - remove de dest as equcoes nao ativas
     */
    void Filter(TPZVec<int64_t> &orig, TPZVec<int64_t> &dest) const
    {
        if (fDestIndices.size() == 0) {
            return;
        }
        else {
            int64_t count = 0;
            int64_t numeq = dest.size();
            for (int64_t i=0; i<numeq; i++) {
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
    void Filter(TPZVec<int64_t> &dest) const
    {
        if (fDestIndices.size() == 0) {
            return;
        }
        else{
            int64_t count = 0;
            int64_t numeq = dest.size();
            for (int64_t i=0; i<numeq; i++) {
                if (fDestIndices[dest[i]] != -1) {
                    dest[count] = fDestIndices[dest[i]];
                    count++;
                }
            }
            dest.Resize(count);
        }
    }

    ///Retorna o numero de equacoes ativas do sistema
    int64_t NActiveEquations() const
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
    int64_t NEqExpand() const
    {
        return fNumEq;
    }

	/**
	 * Returns true if the filter is active
	 */
    bool IsActive() const
    {
        return fIsActive;

    }

	/**
	 * Expands the vector small to a original system, fill zeros into the no active equations.
	 */
    void Scatter(const TPZBaseMatrix &vsmall, TPZBaseMatrix &vexpand) const;

    /**
     * @brief Reduce the vector to the number of active equations.
     */
    void Gather(const TPZBaseMatrix &large, TPZBaseMatrix &gathered) const;
    /**
     * @brief Returns the number of active equations between [minindex,maxindex]
     */
    int64_t NumActive(int64_t minindex, int64_t maxindex) const
    {
        if (minindex < 0 || maxindex < 0 || minindex > fNumEq || maxindex > fNumEq ||
            maxindex < minindex) {
            DebugStop();
        }
        if (!IsActive()) {
            return maxindex-minindex;
        }
        int numactive = 0;
        for (int64_t i=minindex; i<maxindex; i++) {
            if (fDestIndices[i] != -1) {
                numactive++;
            }
        }
        return numactive;
    }
    
    void FilterSkyline(TPZVec<int64_t> &skyline) const
    {
        if (!IsActive()) {
            return;
        }

        for (int64_t ieq = 0; ieq<fActiveEqs.size(); ieq++)
        {
            int64_t skyl = skyline[fActiveEqs[ieq]];
            while (fDestIndices[skyl] == -1 && skyl < fNumEq) {
                skyl++;
            }
#ifdef PZDEBUG
            // all active equations should have a destination
            if (skyl > fActiveEqs[ieq] || fDestIndices[skyl] < 0) {
                DebugStop();
            }
#endif
            skyline[ieq] = fDestIndices[skyl];
        }
        skyline.Resize(fActiveEqs.size());

    }

    void SetNumEq(const int64_t numEq){
        fNumEq = numEq;
    }

private:

    /// Numero de equacoes do sistema original
    int64_t fNumEq{0};
    
    /// Flag indicating whether the filter is active
    bool fIsActive{false};
    
    /// Equacoes ativas
    TPZVec<int64_t> fActiveEqs;

    /// Posicao das equacoes originais no sistema reduzido
    TPZVec<int64_t> fDestIndices;

    template <class TVar>
    void ScatterInternal(const TPZFMatrix<TVar> &vsmall,
                         TPZFMatrix<TVar> &vexpand) const;

    template <class TVar>
    void GatherInternal(const TPZFMatrix<TVar> &large,
                        TPZFMatrix<TVar> &gathered) const;
};

#endif
