/**
 * @file
 */

#ifndef TPZYCTRESCA_H
#define TPZYCTRESCA_H

#include "TPZTensor.h"
#include "pzlog.h"
#include "TPZSavable.h"
#include "TPZPlasticCriterion.h"


class TPZYCTresca : public TPZPlasticCriterion {

public:
  
  enum {NYield = 1};
    
  public:
virtual int ClassId() const override;

  
    const char * Name() const
    {
	   return "TPZYCTresca";	
    }
	
    void Print(std::ostream & out) const override
    {
       out << Name();
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
	 * @param [in] sigma stress state
	 * @param [in] A Thermo Force
	 */
	void SetYieldStatusMode(const TPZTensor<REAL> & sigma, const REAL & A)
	{
		// nothing to be done in this yield criterium
	}
	
  /**
   * Evaluate the yield criteria
   * @param [in] sigma current stress tensor
   * @param [in] A current thermodynamical force
   * @param [out] result
   * @param [in] checkForcedYield indicates wether to force post-peak failure behavior
   */
  template < class T>
  void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &result, int checkForcedYield = 0) const;

  /**
   * Derivative of the yield function
   * @param [in] sigma current stress tensor
   * @param [in] A current thermodynamical force
   * @param [out] Ndir Stress derivative
   * @param [in] checkForcedYield indicates wether to force post-peak failure behavior
   */
  template <class T>
  void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield = 0) const;

  /**
   * Derivative of the yield function with respect to the thermodynamical force
   * @param [in] sigma current stress tensor
   * @param [in] A current thermodynamical force
   * @param [out] h Derivative with respect to thermodynamical force
   * @param [in] checkForcedYield indicates wether to force post-peak failure behavior
   */
  template <class T>
  void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield = 0) const
  {
    h[0] = 1.;
  }

    /**
     * Multiplicador para o caso onde utilizamos uma variavel de dano modificada
     */
    template <class T>
    void AlphaMultiplier(const T &A, T &multiplier) const
    {
        multiplier = T(1.);
    }
    
    void Read(TPZStream& buf, void* context) override {
        
    }
    
    void Write(TPZStream& buf, int withclassid) const override {
        
    }
    void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override{
        TPZTensor<STATE> sigmaTensor;
        sigmaTensor.XX() = sigma[0];
        sigmaTensor.YY() = sigma[1];
        sigmaTensor.ZZ() = sigma[2];
        Compute(sigmaTensor, kprev, yield, 0);
    }
    
    virtual int GetNYield() const override {
        return as_integer(NYield);
    }



protected:
  /**
   * @brief Compute the inverse angle of the tresca yield criterium formula and the related data
   * @param sigma [in] stress tensor
   * @param theta [out] one third of the asin of the inverse angle
   * @param gradtheta [out] gradient of the inverse angle
   */
  template <class T>
  void GradTheta(const TPZTensor<T> & sigma,T & theta,TPZTensor<T> & gradtheta) const;

  /** @brief Computes the inverse angle of the tresca yield criterium formula */
  template <class T>
  T InverseAngle(const TPZTensor<T> &deviatoric) const;

  /** @brief Computes the derivative of the inverse angle of the tresca yield criterium formula */
  template <class T>
  void GradInverseAngle(const TPZTensor<T> &sigma, TPZTensor<T> &grad) const;
 
  public:

//////////////////CheckConv related methods/////////////////////

    /** @brief number of types of residuals */
    int NumCases() 
    {
      return 5;
    }
    static TPZTensor<REAL> gRefTension;

    /** @brief LoadState will keep a given state as static variable of the class */
    void LoadState(TPZFMatrix<REAL> &state)
    {
      int i;
      for(i=0; i<6; i++) gRefTension.fData[i] = state(i,0);
    }
    void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase)
    {
      switch(icase)
      {
        case 0:
        {
          TPZTensor<REAL> grad;
          this->gRefTension.dJ2(grad);
          tangent.Redim(1,6);
          int i;
          for(i=0; i<6; i++)
          {
	    tangent(0,i) = grad.fData[i];
          }
        }
        break;
        case 1:
        {
          TPZTensor<REAL> grad;
          this->gRefTension.dJ3(grad);
          tangent.Redim(1,6);
          int i;
          for(i=0; i<6; i++)
          {
    	    tangent(0,i) = grad.fData[i];
          }
        }
        break;
        case 2:
        {
          TPZTensor<REAL> grad;
          this->GradInverseAngle(gRefTension,grad);
          tangent.Redim(1,6);
          int i;
          for(i=0; i<6; i++)
          {
    	    tangent(0,i) = grad.fData[i];
          }
        }
        break;
        case 3:
        {
          REAL theta;
          TPZTensor<REAL> grad;
          this->GradTheta<REAL>(gRefTension,theta,grad);
          tangent.Redim(1,6);
          int i;
          for(i=0; i<6; i++)
          {
            tangent(0,i) = grad.fData[i];
          }
        }
        break;
        case 4:
        {
          REAL yield = 1.e6;
          TPZVec<TPZTensor<REAL> > grad(1);
          this->N<REAL>(gRefTension,yield,grad,0);
          tangent.Redim(1,6);
          int i;
          for(i=0; i<6; i++)
          {
            tangent(0,i) = grad[0].fData[i];
          }
        }
      }
    }

    void Residual(TPZFMatrix<REAL> &res,int icase)
    {
      REAL inv = this->InverseAngle(gRefTension);
      res.Redim(1,1);
      TPZTensor<REAL> grad;
      switch(icase)
      {
        case 0:
          res(0,0) = gRefTension.J2();
        break;
        case 1:
          res(0,0) = gRefTension.J3();
        break;
        case 2:
          res(0,0) = inv;
        break;
        case 3:
          REAL theta;
          GradTheta<REAL>(gRefTension,theta,grad);
          res(0,0) = theta;
        break;
        case 4:
          TPZVec<REAL> phi(1);
          REAL yield = 1.e6;
          this->Compute(gRefTension,yield,phi,0);
          res(0,0) = phi[0];
        break;
      }
    }

//////////////////CheckConv related methods/////////////////////
public:


};


template < class T>
void TPZYCTresca::Compute(const TPZTensor<T> & sigma, const T & A,TPZVec<T> &result, int checkForcedYield) const
{
  T inverseangle = InverseAngle(sigma);
  T theta = asin(inverseangle)*(1./3.);
  result[0] = sqrt(sigma.J2())*cos(theta)*2.-A;
}

template <class T> 
void TPZYCTresca::N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{
#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.yctresca"));
#endif

/*  T invangle = InverseAngle(sigma);
  T theta = (1./3.)*asin(invangle);
  TPZTensor<T> gradinvangle;
  GradInverseAngle(sigma,gradinvangle);
  T derivasin = 1./sqrt(1.-invangle*invangle);*/
//  T invangle;
  T theta;
  TPZTensor<T> gradtheta;
  //T derivasin;
  //GradASinInvAngle (sigma,invangle,theta,gradinvangle,derivasin);
#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout;
    sout << "A  " << A << "    N: sigma " << sigma;
    LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif
  GradTheta<T> (sigma,theta,gradtheta);
#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout;
    sout << "N: theta " << theta;
    sout << "gradtheta " << gradtheta;
    LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif

  T sqj2 = sqrt(sigma.J2());
  TPZTensor<T> dj2;
  sigma.dJ2(dj2);
  
  TPZTensor<T> resultdtheta(gradtheta);
  resultdtheta.Multiply(sqj2*sin(theta),T(-2.));
  TPZTensor<T> resultdj2(dj2);
  resultdj2.Multiply(cos(theta)/sqj2,T(1.));
  TPZTensor<T> result(resultdtheta);
  result.Add(resultdj2,1.);
#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout;
    sout << "sqj2 " << sqj2;
    sout << "   theta" << theta;
    T coco = sin (theta);
    sout << "  sintheta " << coco;
    sout << "dj2" << dj2;
    sout << "N: resultdtheta " << theta;
    LOGPZ_DEBUG(logger,sout.str().c_str());
  }
  { 
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << std::endl;
    sout << "gradtheta " << gradtheta << std::endl;
    sout << "Result " << result;
    LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif
  Ndir[0] = result;
}

template <class T>
T TPZYCTresca::InverseAngle(const TPZTensor<T> &deviatoric) const
{
  T j2 = deviatoric.J2();
  T j3 = deviatoric.J3();
  if(TPZExtractVal::val(j2) < 1.e-6) return T(0.);
  return j3/j2/sqrt(j2)*(-3.*sqrt(3.)/2.); 
}

template <class T>
void TPZYCTresca::GradInverseAngle(const TPZTensor<T> &sigma, TPZTensor<T> &grad) const
{
  T j2 = sigma.J2();
  if(TPZExtractVal::val(j2) < 1.e-6) return;
  T j3 = sigma.J3();
  TPZTensor<T> dj2,dj3;
  sigma.dJ2(dj2);
  sigma.dJ3(dj3);
  //if(j2 < 1.e-6) return;
  dj3.Multiply(1./(j2*sqrt(j2)),-3.*sqrt(3.)/2.);
  dj2.Multiply(j3/(j2*j2*sqrt(j2)),(9./4.)*sqrt(3.));
  dj2.Add(dj3,1.);
  grad = dj2;
  return;
}

template <class T>
void TPZYCTresca::GradTheta(const TPZTensor<T> & sigma,T & theta, TPZTensor<T> & gradtheta) const
{

#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.yctresca"));

  {
    LOGPZ_DEBUG(logger,__PRETTY_FUNCTION__);
  }
#endif
  T invangle = InverseAngle(sigma);
  theta = (1./3.)*asin(invangle);
  GradInverseAngle(sigma,gradtheta);
  T derivasin = 1./sqrt(1.-invangle*invangle);
  gradtheta.Multiply(derivasin,1./3.);
}

#endif //TPZYCTRESCA_H
