/*
 *  TPZYCWillamWarnke.h
 *  ElastoPlasticModels
 *
 *  Created by Diogo Cecilio on 12/13/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

//$Id: TPZYCDruckerPrager.h,v 1.5 2009-12-14 21:59:52 erick Exp $
#ifndef TPZYCWILLAMWARNKE_H
#define TPZYCWILLAMWARNKE_H

#include "TPZTensor.h"
#include "pzfmatrix.h"

#include "pzsave.h"





/**
 Implementa  a plastificacao do criterio de Von Mises
 */
class TPZYCWillamWarnke: public TPZSaveable {
    
	
public:
	
	
    TPZYCWillamWarnke():fa(1.),fb(1.),fcConcrete(20.),fcoesion(9.2376),fa0(0.1025),fa1(-0.8403),fa2(-0.0910),fb0(0.1025),fb1(-0.4507),fb2(-0.1018){}
	
    TPZYCWillamWarnke(const TPZYCWillamWarnke & source)
    {
		fa  = source.fa;
   		fb  = source.fb;
		fcConcrete = source.fcConcrete;
		fa0 = source.fa0;
		fa1 = source.fa1;
		fa2 = source.fa2;
	    fb0 = source.fb0;
		fb1 = source.fb1;
		fb2 = source.fb2;
    }
	
	enum {NYield = 1};
	
    const char * Name() const
    {
		return "TPZYCDruckerPrager";	
    }
	
    void Print(std::ostream & out) const
    {
		out << Name();
    }
	
	/**
	 * Setup of material parameters
	 * @param [in] phi Mohr Coulomb's internal friction angle
	 * @param [in] innerMCFit If one, Drucker Prager model is inscribed in a referred Mohr Coulomb envelope. If zero, circumscribed.
	 * VERY IMPORTANT!! The ThermoForceA parameters should be set as:
	 *      fk: Herdening slope for the cohesion
	 *      fYield0 : equivalent Mohr Coulomb cohesion C
	 */
	void SetUp(const REAL a, const REAL b,const REAL fConcrete) 
	{
		fa = a;
		fb = b;
		fcConcrete = fConcrete;
				
	}
	
	int GetForceYield()
	{
		return 0; //nothing to be done in this yield criterium
	}
	
	void SetForceYield(const int forceYield)
	{
		//nothing to be done in this yield criterium
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
	 Calculo do criterio de plastificacao 
	 @param [in] sigma tensao atual
	 @param [in] A forca thermodinamica atual
	 @param [in] checkForcedYield indicates wether to force post-peak failure behavior
	 */  
    template < class T>
    void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield = 0) const;
    
    /**
	 Derivada da funcao de plastificacao
	 @param [in] sigma tensao atual
	 @param [in] A forca termodinamica atual
	 @param [out] Derivida com respeito a tensao
	 @param [in] checkForcedYield indicates wether to force post-peak failure behavior
	 */
    template <class T> 
    void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield = 0) const;
	
    /**
	 Derivada da funcao de plastificacao com respeito a forca termodinamica
	 @param [in] sigma tensao atual
	 @param [in] A forca termodinamica atual
	 @param [out] Derivida com respeito a forca termodinamica
	 @param [in] checkForcedYield indicates wether to force post-peak failure behavior
	 */
    template <class T> 
    void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield = 0) const;
	
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
	virtual int ClassId() const;
	
	TPZTensor<REAL> gRefTension;
	
	int NumCases() 
    {
		return 1;
    }
    /**
	 LoadState will keep a given state as static variable of the class
	 */
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
				TPZVec<TPZTensor<REAL> > Ndir(1);
				REAL yield = 1.e6;
				this->N<REAL>(gRefTension,yield, Ndir, 0); 
				tangent.Redim(1,6);
				for(int i=0; i<6; i++)
				{
					tangent(0,i) = Ndir[0].fData[i];
				}
				break;
			}

		}
    }
	
    void Residual(TPZFMatrix<REAL> &res,int icase)
    {
		
		res.Redim(1,1);
		TPZTensor<REAL> grad;
		
		switch(icase)
		{
			case 0:
			{
				TPZVec<REAL> phi(1);
				REAL yield = 1.e6;
				this->Compute(gRefTension,yield,phi,0);
				res(0,0) = phi[0];
				break;
			}
		}
		
    }
	
public:
	REAL fPhi;
	REAL fa, fb, fcConcrete,fcoesion,fa0,fa1,fa2,fb0,fb1,fb2;
protected:
	
	
	//////////////////CheckConv related methods/////////////////////
};

template <class T>
void TPZYCWillamWarnke::Compute(const TPZTensor<T> & sigma,const T & A, TPZVec<T> &res, int checkForcedYield) const
{

	TPZTensor<T> gradlode,dJ2,dJ3;
	T lode,I1,J2,J3;
	TPZTensor<T> eigenval,dSigma1,dSigma2,dSigma3;
	sigma.Eigenvalue(eigenval, dSigma1, dSigma2, dSigma3);
	
	I1 = sigma.I1();
	J2 = sigma.J2();
	J3 = sigma.J3();
	sigma.dJ2(dJ2);
	sigma.dJ3(dJ3);
/*	
	T I13=I1/T(3.);
	T rhoc = -(1./(2.*fb2))*(fb1+ sqrt((fb1*fb1)-T(4.)*T(fb2)*(T(fb0)-I13)));
	T rhot = -(1./(2.*fa2))*(fa1+ sqrt((fa1*fa1)-T(4.)*T(fa2)*(T(fa0)-I13)));
	
//	std::cout << " rhoc " << std::endl;
//	std::cout << rhoc<<std::endl;
//	
//	std::cout << " rhot " << std::endl;
//	std::cout << rhot <<std::endl;
	
	sigma.Lodeangle(gradlode,lode);
	
	T s = (T(2.)*rhoc)*((rhoc*rhoc)-(rhot*rhot));
	T u = T(4.)*(rhoc*rhoc-rhot*rhot)*cos(lode)*cos(lode)+T(5.)*rhot*rhot-T(4.)*rhot*rhoc;
	T t = rhoc*(T(2.)*rhot-rhoc)*sqrt(u);
	T v = T(4.)*(rhoc*rhoc)*cos(lode)*cos(lode)+(rhoc-T(2.)*rhot)*(rhoc-T(2.)*rhot);
//	std::cout << " s " << std::endl;
//	std::cout << s <<std::endl;
//	std::cout << " u " << std::endl;
//	std::cout << u <<std::endl;
//	std::cout << " v " << std::endl;
//	std::cout << v <<std::endl;
//	std::cout << " t " << std::endl;
//	std::cout << t <<std::endl;
	res[0] = (s+t)/v;
	*/
//	std::cout << " PHI " << std::endl;
//	std::cout << res[0] <<std::endl;
	
	//-1 + (a*I1(x))/fc + ((b - c*Cos(3*theta(x)))*Sqrt(J2(x)))/fc


	
	
	sigma.Lodeangle(gradlode,lode);
//	if(fabs( shapeFAD::val(J2) ) < 1.e-6)J2 += 1.e-6;
//	T sqrtJ2 = sqrt(J2);
//	T lode1 =-asin( ( T( 3.) * sqrt( T( 3.) ) * J3 ) /( T( 2.) *  sqrt(J2*J2*J2) ) )/T( 3.);
//	cout << " sigma " << endl;
//	cout << sigma << endl;
//	cout << " gradlode " << endl;
//	cout << gradlode << endl;
//	cout << " lode " << endl;
//	cout << lode << endl;
	T temp1;
	temp1 = ((I1*T(fa))/(T(fcConcrete)));
//	cout << "temp1" <<endl;
//	cout << temp1<<endl;
	T cos3lode = cos(T(3.)*lode);
	T temp4 = ( T(fb) - ( A * cos3lode ) );
//	cout << "fb" <<endl;
//	cout << fb<<endl;
	
//	cout << "A" <<endl;
//	cout << A <<endl;
	
//	cout << "cos3lode" <<endl;
//	cout << cos3lode <<endl;
	
//	cout << "temp4" <<endl;
//	cout << temp4<<endl;
	T temp2 = ( ( temp4 ) * sqrt(J2) )/(T(fcConcrete));
	T resp = (temp1 + temp2) - 1.;
//	std::cout << " PHI " << std::endl;
//	std::cout << resp <<std::endl;

	res[0] = resp;


}

template <class T> 
void TPZYCWillamWarnke::N(const TPZTensor<T> & sigma,const T & A, TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{
	TPZTensor<T> gradlode,dJ2,dJ3,dI1,S;
	T lode,I1,J2,J3;
	
	I1 = sigma.I1();
	J2 = sigma.J2();
	J3 = sigma.J3();
	sigma.dJ2(dJ2);
	sigma.dJ3(dJ3);
	sigma.S(S);
	dI1.Identity();
	sigma.Lodeangle(gradlode,lode);
	
/*	
	T I13=I1/T(3.);

	T rhoc = -(1./(2.*fb2))*(fb1+ sqrt((fb1*fb1)-T(4.)*T(fb2)*(T(fb0)-I13)));
	T rhot = -(1./(2.*fa2))*(fa1+sqrt((fa1*fa1)-T(4.)*T(fa2)*(T(fa0)-I13)));
	
	
	T rhotL = T(1.)/(fa1+T(2.)*T(fa2)*rhot);
	T rhocL = T(1.)/(fb1+T(2.)*T(fb2)*rhoc);
	
	T s = T(2.)*rhoc*(rhoc*rhoc-rhot*rhot);
	T u = T(4.)*(rhoc*rhoc-rhot*rhot)*cos(lode)*cos(lode)+T(5.)*rhot*rhot-T(4.)*rhot*rhoc;
	T t = rhoc*(T(2.)*rhot-rhoc)*sqrt(u);
	T v = T(4.)*(rhoc*rhoc)*cos(lode)*cos(lode)+(rhoc-T(2.)*rhot)*(rhoc-T(2.)*rhot);
	
	T rho = sqrt(T(2.)*J2);
	T temp1 = -(rhotL/(T(3.)*v))*(T(4.)*rho*(rhoc-T(2.)*rhot)+ T(8.)*  rho * rhot * cos(lode) * cos(lode) - ( T(4.)*rhoc*rhot *cos(lode) - T(2.)*rhoc*sqrt(u)+rhoc*(T(2.)*rhot-rhoc)*(T(1.)/sqrt(u))*(T(4.)*rhot*cos(lode)*cos(lode)-T(5.)*rhot+T(2.)*rhoc)));
	T temp21 = T(2.)*rho*(rhoc-T(2.)*rhot)+ T(8.) * rho * rhoc * cos(lode) * cos(lode);//rhocL/(T(3.)*v))
	T temp22 = -( T(2.)*(T(3.)*rhoc*rhoc-rhot*rhot)*cos(lode)+T(2.)*(rhot-rhoc)*sqrt(u)+T(2.)*rhoc*(T(2.)*rhot-rhoc)*(T(1.)/sqrt(u))*((T(2.)*rhoc*cos(lode)*cos(lode))-rhot));
	T temp2 = (rhocL/(T(3.)*v))*(temp21 + temp22);
//	T temp2 =  (rhocL/(T(3.)*v))*(T(2.)*rho*(rhoc-T(2.)*rhot)+ T(8.) * rho * rhoc * cos(lode) * cos(lode) - ( T(2.)*(T(3.)*rhoc*rhoc-rhot*rhot)*cos(lode)+T(2.)*(rhot-rhoc)*sqrt(u)+T(2.)*rhoc*(T(2.)*rhot-rhoc)*(T(1.)/sqrt(u))*(T(2.)*rhoc*cos(lode)*cos(lode)-rhot));
	//A0 = GRAD f/I1
	T A0 = temp1+temp2;
	T temp3 = (sqrt(T(3.))*(rhoc*rhoc-rhot*rhot))/(v*(T(4.)*cos(lode)*cos(lode)-1)*sqrt(J2*J2*J2));
	T tt = ((T(2.)*(T(1.)/sqrt(u))))* cos(lode);
	T ttt = T(4.)*rho*cos(lode);
	T temp4 = ( ttt - rhoc * ( T(1.) + tt * ( T(2.) * rhot - rhoc ) ) );
	//A2 = GRAD f/J3
	T A2 = temp3*temp4;
	//A1 = GRAD f/J2
	T A1 = -((T(3.)*J3)/(T(2.)*J2))*A2;
	
	
	dI1.Multiply(A0, T(1.));
	S.Multiply(A1, T(1.));
	dJ3.Multiply(A2, T(1.));
//	
//	cout << "dI1*A0" <<endl;
//	cout << dI1 << endl;
//	
//	cout << "dJ2*A1" <<endl;
//	cout << dJ2 << endl;
//	
//	cout << "dJ3*A2" <<endl;
//	cout << dJ3 << endl;
	
	dI1.Add(dJ2, T(1));
	dI1.Add(dJ3, T(1));
	
	Ndir[0] = dI1;
*/
	
/*	
	T sqrtJ2 = sqrt(J2);
	T lode1 =-asin( ( T( 3.) * sqrt( T( 3.) ) * J3 ) /( T( 2.) *  sqrt(J2*J2*J2) ) )/T( 3.);
	T temp11= -( T(9.*sqrt(3.)) * J3) / ( T(4.)*sqrt(J2*J2*J2*J2*J2) );
	TPZTensor<T> resultdj2(dJ2);
	resultdj2.Multiply(temp11,T(1.));
	T temp2 = ( T(3.)*sqrt( T( 3. ) ) )/ ( T( 2. ) * sqrt( J2 * J2 * J2 ) );
	TPZTensor<T> resultdj3(dJ3);
	resultdj3.Multiply(temp2,T(1.));
	T temp33 = ( T( 27. ) * ( J3 * J3 ) )/( T( 4. ) *  J2 * J2 * J2);
	T temp3 = ( T( 3. )  *  sqrt( T( 1. )- temp33 ));
	TPZTensor<T>  RESULT(resultdj2);
	RESULT.Add(resultdj3, T(1.));
	RESULT.Multiply( T(1.) / temp3 ,T(-1.));
	
*/	

//	(a*DI1)/fc + (DJ2*(b - c*Cos(3*theta(x) ) ) )/(2.*fc*Sqrt(J2(x))) +(3*c*Dtheta*Sqrt(J2(x))*Sin(3*theta(x)))/fc
	T temp111 = T(fa)/T(fcConcrete);
	dI1.Identity();
	dI1.Multiply(temp111,T(1.));
	T temp222 = ( T(fb) - (A * cos( T(3.) * lode )) )/(T(2.) * T(fcConcrete) * sqrt(J2));
	dJ2.Multiply(temp222, T(1.));
	
	T temp333 = ( T(3.) * A *sqrt(J2) * sin(T(3.) * lode) )/T(fcConcrete);
	gradlode.Multiply(temp333, T(1.));
	dI1.Add(dJ2,T(1));
	dI1.Add(gradlode,T(1));
	Ndir[0] = dI1;

	
}

template <class T> 
void TPZYCWillamWarnke::H(const TPZTensor<T> & sigma,const T & A, TPZVec<T> & h, int checkForcedYield) const
{

    h[0] = 1.;
}


inline void TPZYCWillamWarnke::Write(TPZStream &buf, int withclassid)
{
}
inline void TPZYCWillamWarnke::Read(TPZStream &buf, void *context)
{
}
inline int TPZYCWillamWarnke::ClassId() const
{
	return 888888;
}
inline template class TPZRestoreClass<TPZYCWillamWarnke, 888888> ;

#endif//TPZYCWillamWarnke