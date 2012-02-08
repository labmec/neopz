// $Id: tpzyctrescaregularized.h,v 1.8 2010-06-11 22:12:14 diogo Exp $
/***************************************************************************
 *   Copyright (C) 2004 by N� os �dios                                   *
 *   labmec@labmec.fec.unicamp.br                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TPZYCTRESCAREGULARIZED_H
#define TPZYCTRESCAREGULARIZED_H

#include "TPZYCTresca.h"
#include "pzlog.h"
#include <math.h>
#include <iostream> 

/**
This class implements a tresca yield criteria whrere the gradient of the inverse angle is regularized.
@author LabMeC
*/
class TPZYCTrescaRegularized : public TPZYCTresca
{
public:
		
    const char * Name() const
    {
	   return "TPZYCTrescaRegularized";	
    }
	
    void Print(std::ostream & out) const
    {
       out << Name();
    }
		
  /**
   * Calculo do criterio de plastificacao
   * @param [in] sigma tensao atual
   * @param [in] A forca thermodinamica atual
   * @param [in] checkForcedYield indicates wether to force post-peak failure behavior
   */	
  template < class T>
  void Compute(const TPZTensor<T> & sigma, const T & A,TPZVec<T> &result, int checkForcedYield = 0) const;

protected:
  /**
   * Compute the inverse angle of the tresca yield criterium formula and
   * the related data
   * @param sigma [in] stress tensor
   * @param invangle [out] inverse angle
   * @param theta [out] one third of the asin of the inverse angle
   * @param grainvangle [out] gradient of the inverse angle
   * @param derivasin [out] derivative of the inverse angle
   */
  template <class T>
  void GradTheta(const TPZTensor<T> & sigma,T & theta, TPZTensor<T> & gradasin) const;

public:
/**
 * Derivada da funcao de plastificacao
 * @param [in] sigma tensao atual
 * @param [in] A forca termodinamica atual
 * @param [out] Derivida com respeito a tensao
 * @param [in] checkForcedYield indicates wether to force post-peak failure behavior
 */
template <class T> 
void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield = 0) const;

//////////////////CheckConv related methods/////////////////////
public:

/**
number of types of residuals
*/
int NumCases()
{
  return 2;
}

static TPZTensor<REAL> gRefTension;
/**
LoadState will keep a given state as static variable of the class
*/
void LoadState(TPZFMatrix &state)
{
  int i;
  for(i=0; i<6; i++) gRefTension.fData[i] = state(i,0);
}

void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &, int icase)
{
  TPZTensor<double> gradtheta;
  double theta,A(1.e7);
  TPZVec<double> phivec(1,0.);
  TPZVec<TPZTensor<double> > ndir(1);
  int i;
  switch(icase)
  {
    case 0:
      GradTheta(gRefTension,theta,gradtheta);
      tangent.Redim(1,6);
      for(i=0; i<6; i++) tangent(0,i) = gradtheta.fData[i];
      break;
    case 1:
      N(gRefTension,A,ndir,0);
      tangent.Redim(1,6);
      for(i=0; i<6; i++) tangent(0,i) = ndir[0].fData[i];
  }
}

void Residual(TPZFMatrix &res,int icase)
{
  TPZTensor<double> gradtheta;
  double theta,A(1.e7);
  TPZVec<double> phivec(1,0.);
  switch(icase)
  {
    case 0:
      GradTheta(gRefTension,theta,gradtheta);
      res.Redim(1,1);
      res(0,0) = theta;
      break;
    case 1:
      Compute(gRefTension,A,phivec,0);
      res.Redim(1,1);
      res(0,0) = phivec[0];
      break;
  }

}
//////////////////CheckConv related methods/////////////////////

};


/**
 * Calculo do criterio de plastificacao
 * @param [in] sigma tensao atual
 * @param [in] A forca thermodinamica atual
 */

template < class T>
void TPZYCTrescaRegularized::Compute(const TPZTensor<T> & sigma, const T & A,TPZVec<T> &result, int checkForcedYield) const
{
#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.yctresca"));
#endif
	
//  result[0] = sqrt(sigma.J2()) - A;

//  return;
  const double tol = 1.e-4;
  T invangle = InverseAngle(sigma);
  TPZTensor <T> s;
  sigma.S(s);
  double aux = (1. - tol);
  if (fabs(invangle) < aux)
  {
#ifdef LOG4CXX_PLASTICITY
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << " O angulo nao corresponde " << invangle;
    LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
    TPZYCTresca::Compute(sigma,A,result,checkForcedYield);
  } else
  {
    if (fabs(shapeFAD::val(invangle)) > 1.) invangle = 1.;
    REAL asineps = asin(1. - tol);
    T tayext = (fabs(invangle) - T(1.-tol))*T(1./(sqrt(1.- (1. - tol)*(1. - tol))));
    invangle = (tayext + T(asineps)) / 3.;
    //invangle = alpha * (1./3.) * (asin(1. - 1.e-6) + ( fabs(alpha)-(1.-1.e-6))*(1./(sqrt(1.-(1.-1.e-6)*(1.-1.e-6)))));
#ifdef LOG4CXX_PLASTICITY
    std::stringstream sout;
    T small = cos(invangle)*2.-T(sqrt(3.));
    
    sout << __PRETTY_FUNCTION__ << "2.*cos(invangle)-sqrt(3.) " << small;
    sout << " invangle " << invangle;
    LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
    result[0] = sqrt(sigma.J2())*cos(invangle)*2. - A;
    //result[0] = sqrt(3.)* sqrt(sigma.J2()) - A;
  }
}


/**
 * Computes the inverse angle of the tresca yield criterium formula and
 * the related data
 * @param sigma [in] stress tensor
 * @param invangle [out] inverse angle
 * @param theta [out] one third of the asin of the inverse angle
 * @param grainvangle [out] gradient of the inverse angle
 * @param derivasin [out] derivative of the inverse angle
 */
template <class T>
void TPZYCTrescaRegularized::GradTheta(const TPZTensor<T> & sigma,T & theta, TPZTensor<T> & gradtheta) const
{
#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.yctresca"));
#endif

  const double tol = 1.e-4;

  T invangle = InverseAngle(sigma);

  {
/*    std::stringstream sout;
    sout << "GradTheta:Invangle " << invangle;
    LOGPZ_DEBUG(logger,sout.str().c_str());*/
  }

  if (fabs(invangle) < (1. - 1.e-6))
  {
	#ifdef LOG4CXX_PLASTICITY
    {
      std::stringstream sout;
      sout << "calling GradTheta from father... invangle " << invangle;
      LOGPZ_DEBUG(logger,sout.str().c_str());
    }
    #endif

    TPZYCTresca::GradTheta(sigma,theta,gradtheta);
    #ifdef LOG4CXX_PLASTICITY
    {
      std::stringstream sout;
      sout << "sigma "  << sigma << std::endl;
      sout << "theta " << theta << std::endl;
      sout << "gratheta " << gradtheta << std::endl;
      LOGPZ_DEBUG(logger,sout.str().c_str());
    }
    #endif


  } else
  {
    //double alpha;
    //if (fabs(invangle) > 1.)
//     alpha = (invangle > 0.) ? 1. : -1.;
//     theta = (alpha/3.)* (asin(1. - 1.e-6) + (fabs(alpha) - (1.-1.e-6))*(1./(sqrt(1.-(1.-1.e-6)*(1.-1.e-6)))));
//     GradInverseAngle(sigma,gradtheta);
//     T derivasin = 1./sqrt(1.-(1. - 1.e-6)*(1. - 1.e-6));
//     gradtheta.Multiply(derivasin,1./3.);
    REAL alpha = (shapeFAD::val(invangle) > 0.) ? 1. : -1.;
    if (fabs(shapeFAD::val(invangle) ) > 1.) invangle = 1.;
    REAL asineps = asin(1. - tol);
    T tayext = (fabs(invangle) - T(1.- tol))*T(1./(sqrt(1.- (1. - tol)*(1. - tol))));
    theta = (tayext + T(asineps)) * T(alpha / 3.);
    GradInverseAngle(sigma,gradtheta);
    REAL derivasin = 1./sqrt(1.-(1. - tol)*(1. - tol));
    gradtheta.Multiply(T(derivasin),T(1./3.));
    {
/*      std::stringstream sout;
      sout << "Invangle " << invangle << std::endl;
      sout << "alpha    " << alpha << std::endl;
      sout << "asineps  " << asineps << std::endl;
      sout << "tayext   " << tayext << std::endl;
      sout << "theta    " << theta << std::endl;
      sout << "gradtheta" << gradtheta << std::endl;
      sout << "derivsin " << derivasin << std::endl;
      sout << "gradtheta" << gradtheta << std::endl;
      LOGPZ_DEBUG(logger,sout.str().c_str());*/
    }
  }
}

/**
 * Derivada da funcao de plastificacao
 * @param [in] sigma tensao atual
 * @param [in] A forca termodinamica atual
 * @param [out] Derivida com respeito a tensao
 * @param [in] forceYield indicates wether to force post-peak failure behavior
 */
template <class T> 
void TPZYCTrescaRegularized::N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{
#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.yctresca"));
#endif

//   sigma.dJ2(Ndir[0]);
//    Ndir[0].Multiply(sqrt(3.)*0.5*(1./sqrt(sigma.J2())),1.);
//  return;



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
  {
/*    std::stringstream sout;
    sout << "A  " << A << "    N: sigma " << sigma;
    LOGPZ_DEBUG(logger,sout.str().c_str());*/
  }

  GradTheta<T> (sigma,theta,gradtheta);

  {
/*    std::stringstream sout;
    sout << "N: theta " << theta;
    sout << "gradtheta " << gradtheta;
    LOGPZ_DEBUG(logger,sout.str().c_str());*/
  }

  T sqj2 = sqrt(sigma.J2());
  TPZTensor<T> dj2;
  sigma.dJ2(dj2);
  
  TPZTensor<T> resultdtheta(gradtheta);
  resultdtheta.Multiply(sqj2*sin(theta),-2.);
  TPZTensor<T> resultdj2(dj2);
   resultdj2.Multiply(cos(theta)/sqj2,1.);
   TPZTensor<T> result(resultdj2);
//   result.Add(resultdtheta,1.);

//  resultdj2.Multiply(1./sqj2,sqrt(3.)/2.);


  {
/*    std::stringstream sout;
    sout << "sqj2 " << sqj2;
    sout << "   theta" << theta;
    T coco = sin (theta);
    sout << "  sintheta " << coco;
    sout << "dj2" << dj2;
    sout << "N: resultdtheta " << theta;
    LOGPZ_DEBUG(logger,sout.str().c_str());*/
  }
  { 
/*    std::stringstream sout;
    sout << "gradtheta " << gradtheta;
    sout << "Result N = " << result;
    LOGPZ_DEBUG(logger,sout.str().c_str());*/
  }

  Ndir[0] = result;
//  Ndir[0] = resultdj2;
}




#endif
