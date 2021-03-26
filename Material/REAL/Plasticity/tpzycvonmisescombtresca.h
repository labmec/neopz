/**
 * @file
 */

#ifndef TPZYCVONMISESCOMBTRESCA_H
#define TPZYCVONMISESCOMBTRESCA_H

#include "pzlog.h"

#include "TPZYCVonMises.h"
#include "tpzyctrescaregularized.h"

/**
 * @author LabMeC
 */
class TPZYCVonMisesCombTresca : public TPZPlasticCriterion {
	
public:
		
    public:
virtual int ClassId() const override;

    
	const char * Name() const
    {
	   return "TPZYCVonMisesCombTresca";	
    }
	
    void Print(std::ostream & out) const override
    {
		out << "\n" << this->Name();
		out << "\n fVonMises :\n";
		fVonMises.Print(out);
		out << "\n fTresca   :\n";
		fTresca.Print(out);
    }
	
	int GetForceYield()
	{
		return 0; // nothing to be done in this yield criterium
	}
	
	void SetForceYield(const int forceYield)
	{
		// nothing to be done in this yield criterium
	}
	
	/**
	 * Checks if the proposed yield state leads to post-peak material behaviour. If so, the material
	 * is forced to behave in post-peak in order to avoid equation switching during Newton's method
	 * in the PlasticLoop routines.
	 * @param[in] sigma stress state
	 * @param[in] A Thermo Force
	 */
	void SetYieldStatusMode(const TPZTensor<REAL> & sigma, const REAL & A)
	{
		// nothing to be done in this yield criterium
	}
    
    void Write(TPZStream& buf, int withclassid) const override{
        fVonMises.Write(buf, withclassid);
        fTresca.Write(buf, withclassid);
    }

    void Read(TPZStream& buf, void* context) override{
        fVonMises.Read(buf, context);
        fTresca.Read(buf, context);
    }
    
protected:
  /** @brief Pointer to Von Mises's yield criteria object */
  TPZYCVonMises fVonMises;

  /** @brief Pointer to Tresca's yield criteria object */
  TPZYCTrescaRegularized   fTresca;

public:

  enum { NYield = TPZYCVonMises::NYield + TPZYCTrescaRegularized::NYield };

  /** @brief Constructor */
  TPZYCVonMisesCombTresca();

  /**
   *Calculo do criterio de plastificacao
   * @param[in] sigma tensao atual
   * @param[in] A forca thermodinamica atual
   * @param[out] res Derivative
   * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
   */  
  template < class T>
  void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield = 0) const;

  /**
   *Derivada da funcao de plastificacao
   * @param[in] sigma tensao atual
   * @param[in] A forca termodinamica atual
   * @param[out] NDir Derivada com respeito a tensao
   * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
   */
  template <class T>
  void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & NDir, int checkForcedYield = 0) const;

  /**
   *Derivada da funcao de plastificacao com respeito a forca termodinamica
   *@param[in] sigma tensao atual
   *@param[in] A forca termodinamica atual
   *@param[out] h Derivada com respeito a forca termodinamica
   *@param[in] checkForcedYield indicates wether to force post-peak failure behavior
   */
  template <class T>
  void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield = 0) const;

    /**
     * Multiplicador para o caso onde utilizamos uma variavel de dano modificada
     */
    template <class T>
    void AlphaMultiplier(const T &A, T &multiplier) const
    {
        multiplier = T(1.);
    }
    
    void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override{
        TPZTensor<STATE> sigmaTensor;
        sigmaTensor.XX() = sigma[0];
        sigmaTensor.YY() = sigma[1];
        sigmaTensor.ZZ() = sigma[2];
        Compute(sigmaTensor, kprev, yield, 0);
    }
    
    virtual int GetNYield() const override{
        return as_integer(NYield);
    }
    
};



inline TPZYCVonMisesCombTresca::TPZYCVonMisesCombTresca() : fVonMises(), fTresca()
{
}


template < class T>
void TPZYCVonMisesCombTresca::Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield) const
{
#ifdef PZ_LOG
  TPZLogger logger("plasticity.ycvonmisescombtresca");
#endif
  if (res.NElements()!=2)
  {
#ifdef PZ_LOG
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << " wrong vector size: res has dimension = " << res.NElements()
          << " waiting for size = 2";
    LOGPZ_ERROR(logger,sout.str().c_str());
#endif
    return;
  }
  TPZVec<T> aux (1);
  const T aux1 = A*T(cos(3.141592/6.)/cos(3.141592/12.));
//  const T aux1 = A;
  fVonMises.Compute(sigma,aux1,aux,checkForcedYield);
  res[0] = aux[0];
  fTresca.Compute(sigma,A,aux,checkForcedYield);
  res[1] = aux[0];
#ifdef PZ_LOG
  {
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << " valores dos phi " << res;
    sout << "\nvalor aux1 " << aux1;
    sout << "\nvalor A " << A;
    LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif
}


template <class T>
void TPZYCVonMisesCombTresca::N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & NDir, int checkForcedYield) const
{
#ifdef PZ_LOG
  TPZLogger logger("plasticity.ycvonmisescombtresca");
#endif
  if (NDir.NElements()!=2)
  {
#ifdef PZ_LOG
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << " wrong vector size: Ndir has dimension = " << NDir.NElements()
          << " waiting for size = 2";
    LOGPZ_ERROR(logger,sout.str().c_str());
#endif
    return;
  }
  TPZVec< TPZTensor<T> > aux (1);
  fVonMises.N(sigma,A,aux,checkForcedYield);
  NDir[0] = aux[0];
  fTresca.N(sigma,A,aux,checkForcedYield);
  NDir[1] = aux[0];
}


template <class T>
void TPZYCVonMisesCombTresca::H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield) const
{
#ifdef PZ_LOG
  TPZLogger logger("plasticity.ycvonmisescombtresca");
#endif
  if (h.NElements()!=2)
  {
#ifdef PZ_LOG
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << " wrong vector size: h has dimension = " << h.NElements()
          << " waiting for size = 2";
    LOGPZ_ERROR(logger,sout.str().c_str());
#endif
    return;
  }
  TPZVec<T> aux (1);
  fVonMises.H(sigma,A,aux,checkForcedYield);
  h[0] = aux[0];
  fTresca.H(sigma,A,aux,checkForcedYield);
  h[1] = aux [0];
}

#endif
