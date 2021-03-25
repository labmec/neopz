// $Id: TPZYCSandlerDimaggio.h,v 1.11 2009-06-29 22:54:01 erick Exp $

#ifndef TPZYCSANDLERDIMAGGIOL_H
#define TPZYCSANDLERDIMAGGIOL_H

#include "TPZYCSandlerDimaggio.h"
#include "pzlog.h"

#ifndef CHECKCONV
#define CHECKCONV
#include "checkconv.h"
#endif

#include "fadType.h"

#ifdef PZ_LOG
#include "pzlog.h"

static TPZLogger loggerSML("material.plasticity.SML");

#endif

/**
Implementa as funções de potencial plástico e yield criterium do 
modelo constitutivo associativo de Sandler e Dimaggio (1971), desenvolvido
inicialmente para arenitos (Ranch McCormic Sand)
*/
class TPZYCSandlerDimaggioL : public TPZYCSandlerDimaggio {
    

public:

  enum {NYield = 2};
  
int ClassId() const override;

	
    TPZYCSandlerDimaggioL():TPZYCSandlerDimaggio() { }
	
    TPZYCSandlerDimaggioL(const TPZYCSandlerDimaggioL & source) : TPZYCSandlerDimaggio(source)
    {
    }

    TPZYCSandlerDimaggioL & operator=(const TPZYCSandlerDimaggioL & source)
    {
        TPZYCSandlerDimaggio::operator=(source);
		return *this;
    }
    
    virtual ~TPZYCSandlerDimaggioL()
    {
        
    }


	const char * Name() const
    {
	   return "TPZYCSandlerDimaggioL";	
    }
	
    void Print(std::ostream & out) const override
    {
		out << "\n" << this->Name();
        TPZYCSandlerDimaggio::Print(out);
    }
	
	
    /**
	 * @brief Calculo do criterio de plastificacao 
	 * @param[in] sigma tensao atual
	 * @param[in] A forca thermodinamica atual
	 * @param[out] res Result
	 * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */  
    template < class T>
    void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield) const;

    /**
    Derivada da derivada da funcao de potencial plastico (direção de plastificação)
    @param[in] sigma tensao atual
    @param[in] A forca termodinamica atual
    @param[out] Ndir Derivada com respeito a tensao
	@param[in] checkForcedYield indicates wether to force post-peak failure behavior
    */
    template <class T> 
    void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const;

    /**
    Derivada da funcao de plastificacao com respeito a forca termodinamica
    @param[in] sigma tensao atual
    @param[in] A forca termodinamica atual
    @param[out] h Derivada com respeito a forca termodinamica
	@param[in] checkForcedYield indicates wether to force post-peak failure behavior
    */
    template <class T> 
    void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield) const;

    /**
     * Multiplicador para o caso onde utilizamos uma variavel de dano modificada
     */
    template <class T>
    void AlphaMultiplier(const T &A, T &multiplier) const;
    
    /**
     Projeto o ponto sigtrial sobre a superficie de plastificacao (se precisar) e atualiza a variavel de dano
     Este metodo utiliza backtracking
     @param[in] ER resposta elastica
     @param[in] variavel de dano atual
     @param[in] tensao da resposta elastica
     @param[out] epspproj dano apos a projecao
     @param[out] delgamma fatores multiplicadores para projetar o ponto sobre a superficie
     @param[sigproj] tensor de tensao projetado
     */
    void InitialGuess(const TPZElasticResponse &ER, REAL L, TPZTensor<REAL> &sigtrial, REAL &epspproj, 
                      TPZVec<REAL> &delgamma, TPZTensor<REAL> &sigproj);

    inline REAL InitialDamage()
    {
        // for this particular case!!! toto
        REAL X,L,alpha(0.);
        ComputeX(alpha, X);
        SolveL(X, L);
//        L -= 0.01;
        return L;

    }
    // 3 K Depsp/DL (L-Lextern) = (I1(sigtrial)-L)
    void ComputeLatIntersectionLeft(const TPZElasticResponse &ER, REAL &L, TPZVec<REAL> &sigtrialIJ) const;
    
    // 3 K Depsp/DL (LIntersect-Lextern) = (I1(sigtrial)-I1(Lintersect) )
    void ComputeLatIntersectionRight(const TPZElasticResponse &ER, REAL &L, TPZVec<REAL> &sigtrialIJ) const;
private:
    
	void NewtonF2(const TPZElasticResponse &ER, REAL &L, TPZVec<REAL> &sigtrialIJ)
    {
        DebugStop();
    }
public:

    virtual int GetNYield() const  override {
        return as_integer(NYield);
    }
		
//////////////////CheckConv related methods/////////////////////

    /**
    number of types of residuals
    */
    inline int NumCases();

    static TPZTensor<REAL> gRefTension;
    /**
    LoadState will keep a given state as static variable of the class
    */
    inline void LoadState(TPZFMatrix<REAL> &state);

    inline void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &, int icase);

    inline void Residual(TPZFMatrix<REAL> &res,int icase);

    static void CheckConv();
	
//////////////////CheckConv related methods/////////////////////

//////////////////Internal routines verification/////////////////

	static void TestSolveL();
//////////////////Internal routines verification/////////////////

	static void McCormicRanchSand(TPZYCSandlerDimaggio & material);
public:
  
};



template <class T>
inline void TPZYCSandlerDimaggioL::Compute(const TPZTensor<T> & sigma,const T & A, TPZVec<T> &res, int checkForcedYield) const
{
	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	// The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.
#ifdef PZ_LOG
    {
        TPZLogger logger("plasticity.SandlerDimaggioL");
        std::stringstream sout;
        sout << ">>> TPZYCSandlerDimaggio::Compute *** - Plastic Potential / Yield - associative";
        LOGPZ_INFO(logger,sout.str().c_str());
    }
#endif
    
//	REAL L_REAL = TPZExtractVal::val(A);
//    REAL X_REAL;
//	ComputeX((REAL)TPZExtractVal::val(A), X_REAL);
//	LInitialGuess(X_REAL, L_REAL);
//	SolveL(X_REAL, L_REAL);
    
    T I1 = sigma.I1();
    T J2 = sigma.J2();
	
    if (fIsonCap == false) 
    {
        // f1 - Modified Drucker-Prager as shear Yield Criterium
        T FI1;
        ComputeF(I1, FI1);
        if(fabs((REAL)TPZExtractVal::val(J2)) < 1.e-6)
        {
            res[0] = - FI1;	// avoiding nan derivatives
        }else{
            res[0] = sqrt(J2) - FI1;
        }
        
        // f2 - ellipsoidal hardening/softening cap
        T L = A;
        T FL;
    //    T L(L_REAL), X;
    //   	ComputeX(A, X);
    //	SolveL(X, L); // evaluating the derivatives of L
        ComputeF(L, FL);
        
        if(fabs( (REAL)TPZExtractVal::val(FL) ) < 0.00001)
        {
    #ifdef PZ_LOG
            {
                TPZLogger logger("plasticity.SandlerDimaggio");
                std::stringstream sout;
                sout << "*** TPZYCSandlerDimaggio::ComputePlasticPotential ***";
                sout << "\nDivision by F=" << TPZExtractVal::val(L) << " at f2 - ellipsoidal hardening/softening cap";
                LOGPZ_WARN(logger,sout.str().c_str());
            }
    #endif
        }
        
        T Temp1( (L - I1)/(FL * T(fR) ) );
        Temp1 *= Temp1;
        T Temp2 = J2 / FL / FL;
        
        res[1] = Temp1 + Temp2 - T(1.);
    }
	else {
        REAL lmax = LMax();
        REAL FI1;
        ComputeF(lmax , FI1);
        if(fabs((REAL)TPZExtractVal::val(J2)) < 1.e-6)
        {
            res[1] = - FI1;	// avoiding nan derivatives
        }else{
            res[1] = sqrt(J2) - FI1;
        }
        res[0] = I1-T(lmax);
#ifdef PZDEBUG
        {
            std::stringstream sout;
            T sqj2 = J2;
            if(fabs(TPZExtractVal::val(J2)) > 1.e-6)
            {
                sqj2 = sqrt(J2);
            }
                
            sout << "Computing the distance from cap entry I1 " << I1 << " sqJ2 " << sqj2 << " lmax " << lmax << " F(lmax) " << FI1;
            LOGPZ_DEBUG(loggerSML, sout.str())
        }
#endif
    }
	return;
}

template <class T> 
inline void TPZYCSandlerDimaggioL::N(const TPZTensor<T> & sigma, const T & A, TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{

	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	//The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.

	
    #ifdef PZ_LOG
        {
          TPZLogger logger("plasticity.SandlerDimaggioL");
          std::stringstream sout;
          sout << ">>> TPZYCSandlerDimaggio::N *** - Plastification direction - associative";
          LOGPZ_INFO(logger,sout.str().c_str());
        }
    #endif
	
	REAL ResTol = 1.e-6;
	
	REAL L_REAL;
//  REAL X_REAL;
//	ComputeX((REAL)TPZExtractVal::val(A), X_REAL);
//	LInitialGuess(X_REAL, L_REAL);
//	SolveL(X_REAL, L_REAL, ResTol);

    L_REAL = TPZExtractVal::val(A);
    T I1 = sigma.I1();
    T J2 = sigma.J2();
    T SQRTJ2;
    if (TPZExtractVal::val(J2) > 1.e-12) {
        SQRTJ2 = sqrt(J2);
    }
    else {
        SQRTJ2 = sqrt(T(1.e-6)+J2);
    }
	
    if(fIsonCap == false)
	{
		//f1 - Modified Drucker-Prager as shear Yield Criterium / Plastic Potential
        REAL fz = FZero();
        T Temp1(0.);
        if (TPZExtractVal::val(I1) < fz ) 
        {
            Temp1 = I1 * T(fB);
            Temp1 = exp( Temp1 ) * T (fB * fC);
		
            if((REAL)TPZExtractVal::val(SQRTJ2) < 1.e-6) // just for robustness. f1 shouldn't be reached when J2 = 0.
            {
                #ifdef PZ_LOG
                {
                   TPZLogger logger("plasticity.SandlerDimaggio");
                   std::stringstream sout;
                   sout << "*** TPZYCSandlerDimaggio::N *** - SQRT(J2) = " << TPZExtractVal::val(SQRTJ2) <<  " < 1.e-6 causes error in 0-th yield function. Imposing J2 = 1.e-6 instead";
                   LOGPZ_WARN(logger,sout.str().c_str());
                }
                #endif
                SQRTJ2 = T(1.e-6)+SQRTJ2;
            }
        }		
		Temp1 = Temp1 - I1 / SQRTJ2 / T(6.);
		T Temp2 = T(1.) / SQRTJ2;
		T Temp3 = Temp2 / T(2.);
		
		Ndir[0].XX() = Temp1 + sigma.XX() * Temp3;
		Ndir[0].YY() = Temp1 + sigma.YY() * Temp3;
		Ndir[0].ZZ() = Temp1 + sigma.ZZ() * Temp3;
//	    Ndir[0].YZ() = sigma.YZ() * Temp2;
//	    Ndir[0].XZ() = sigma.XZ() * Temp2;
//	    Ndir[0].XY() = sigma.XY() * Temp2;
	    Ndir[0].YZ() = sigma.YZ() * Temp3;
	    Ndir[0].XZ() = sigma.XZ() * Temp3;
	    Ndir[0].XY() = sigma.XY() * Temp3;
#ifdef PZ_LOG
        {
            std::stringstream sout;
            Ndir[0].Print(sout);
            LOGPZ_DEBUG(loggerSML, sout.str())
        }
#endif
    }
    else {
		Ndir[0].XX() = 1.;
		Ndir[0].YY() = 1.;
		Ndir[0].ZZ() = 1.;
        //	    Ndir[0].YZ() = sigma.YZ() * Temp2;
        //	    Ndir[0].XZ() = sigma.XZ() * Temp2;
        //	    Ndir[0].XY() = sigma.XY() * Temp2;
	    Ndir[0].YZ() = 0.;
	    Ndir[0].XZ() = 0.;
	    Ndir[0].XY() = 0.;
        
    }
	
    if (fIsonCap == false)
    {//f2 - ellipsoidal hardening/softening cap

    T FL;
    T L = A;
    T X;
//        T L(L_REAL * 1.- ResTol); // guaranteeing that the function will be evaluated
    ComputeX(A, X);
//		SolveL(X, L, ResTol); // evaluating the derivatives of L
		
		ComputeF(L, FL);
        // the radius of the ellips needs to be positive
        // this should be taken care of by the computation of L which is limited by LMax() 
        if (TPZExtractVal::val(FL) <= 0.) {
            DebugStop();
        }
		T FL2 = FL * FL;
		T FL3 = FL2;// / T(2.);
	
		T Temp = (I1-L)/ T(fR * fR) - I1 / T(6.);
		Temp = Temp / FL2 * T(2.);

#ifdef PZ_LOG
            {
               TPZLogger logger("plasticity.SandlerDimaggio");
               std::stringstream sout;
               sout << "*** TPZYCSandlerDimaggio::N *** X = " << X
					<< "\n L = " << L << " L_REAL = " << L_REAL
					<< "\n FL = " << FL
					<< "\n Temp = " << Temp
                << "\n 3 Temp + I1/FL2 " << T(3.)*Temp+I1/FL2;
                
               LOGPZ_DEBUG(logger,sout.str().c_str());
            }
#endif
		
		Ndir[1].XX() = Temp + sigma.XX() / FL2;
		Ndir[1].YY() = Temp + sigma.YY() / FL2;
		Ndir[1].ZZ() = Temp + sigma.ZZ() / FL2;
        Ndir[1].YZ() = sigma.YZ() / FL3;
        Ndir[1].XZ() = sigma.XZ() / FL3;
        Ndir[1].XY() = sigma.XY() / FL3;
	}
    else {
        TPZTensor<T> Deviatoric;
        sigma.S(Deviatoric);
        Ndir[1] = Deviatoric;
    }
	
    #ifdef PZ_LOG
    TPZLogger logger("pz.plasticity.SandlerDimaggio.main");
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "<< TPZYCSandlerDimaggio::N *** \n sigma = \n" << sigma
			 << "\nI1 = " << I1 
			 << "\nJ2 = " << J2
			 << "\nSQRTJ2 = " << SQRTJ2
			 << "\nNdir = \n" << Ndir;
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
    #endif
	
	return;
}

/**
 * Multiplicador para o caso onde utilizamos uma variavel de dano modificada
 */
template <class T>
inline void TPZYCSandlerDimaggioL::AlphaMultiplier(const T &A, T &multiplier) const
{
    this->DEpspDL(A, multiplier);
}

template <class T> 
inline void TPZYCSandlerDimaggioL::H(const TPZTensor<T> & sigma,const T & A, TPZVec<T> & h, int checkForcedYield) const
{

 	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	//The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.
	
    #ifdef PZ_LOG
        {
          TPZLogger logger("plasticity.SandlerDimaggio");
          std::stringstream sout;
          sout << ">>> TPZYCSandlerDimaggio::H *** - Hardening modulus";
          LOGPZ_INFO(logger,sout.str().c_str());
        }
    #endif
	
//	REAL L_REAL = TPZExtractVal::val(A);
//    REAL X_REAL;
//	ComputeX((REAL)TPZExtractVal::val(A), X_REAL);
//	LInitialGuess(X_REAL, L_REAL);
//	SolveL(X_REAL, L_REAL);
    
    if (fIsonCap == false) 
    {
        T I1 = sigma.I1();
        
        {//f1 - Modified Drucker-Prager as shear Yield Criterium / Plastic Potential

            h[0] =  exp (I1 * T(fB) ) * T(3. * fB * fC)  ;
            // h > 0 because plastic deformation is dilatant
        }	

        {//f2 - ellipsoidal hardening/softening cap
            T FL;
            T L = A;
    //        T X, L(L_REAL);
    //    	ComputeX(A, X);
    //		SolveL(X, L); // evaluating the derivatives of L
            ComputeF(L, FL);
            T FL2 = FL * FL;
    //		DebugStop();
            h[1] =  (I1 - L) / FL2 * T ( 6. / fR / fR);
            // h <=0 because plastic deformation is compactant
        }
    }
	else
    {
        h[0] = 0.;
        h[1] = 0.;
    }
	return;
	
}



//////////////////CheckConv related methods/////////////////////

//////////////////Internal routines verification/////////////////

/**
 Projeto o ponto sigtrial sobre a superficie de plastificacao (se precisar) e atualiza a variavel de dano
 Este metodo utiliza backtracking
 @param[in] ER resposta elastica
 @param[in] variavel de dano atual
 @param[in] tensao da resposta elastica
 @param[out] epspproj dano apos a projecao
 @param[out] delgamma fatores multiplicadores para projetar o ponto sobre a superficie
 @param[sigproj] tensor de tensao projetado
 */
inline void TPZYCSandlerDimaggioL::InitialGuess(const TPZElasticResponse &ER, REAL Lextern, TPZTensor<REAL> &sigtrial, REAL &Lproj, 
                  TPZVec<REAL> &delgamma, TPZTensor<REAL> &sigproj)
{
    TPZManVector<REAL,2> yield(2,0.);
    REAL L = Lextern;
    TPZTensor<REAL> S;
    sigtrial.S(S);
    REAL J2 = S.J2();
    REAL sqJ2 = sqrt(fabs(J2));
    REAL I1 = sigtrial.I1();
    TPZManVector<REAL,2> sigtrialIJ(2),sigtrialIJF1(2),sigtrialIJF2(2),sigtrialIJkeep;
    sigtrialIJ[0] = I1;
    sigtrialIJ[1] = sqJ2;
    sigtrialIJkeep = sigtrialIJ;
    Compute(sigtrial, L, yield, 0);
    int surfaceprojected = -1;
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "Value of fIsonCap " << fIsonCap;
        if (fIsonCap == true) {
            sout << "\nI should stop";
        }
        LOGPZ_DEBUG(loggerSML, sout.str())
    }
#endif
    fIsonCap = false;

    if (yield[1] <= 0.) {
        sigtrialIJF1 = sigtrialIJ;
        if (yield[0] > 0.) 
        {
            NewtonF1(ER, L, sigtrialIJF1);
            sigtrial.Adjust(sigtrialIJF1, sigproj);
            surfaceprojected = 0;

        }
        else {
            // fully elastic
            Lproj = Lextern;
            sigproj = sigtrial;
            delgamma.Fill(0.);
#ifdef PZ_LOG
            if(loggerSML.isDebugEnabled())
            {
                std::stringstream sout;
                if(yield[0] <= 0.)
                {
                    sout << "ELASTIC YIELD\n";
                }
                sout << "Deformation elastic yield = " << yield;
                sout << "delgamma condition " << (delgamma[0] > 0.) << " " << (delgamma[1] > 0.) << std::endl;
                LOGPZ_DEBUG(loggerSML, sout.str())
            }
#endif
            return;

        }
        if (sigtrialIJF1[0] <= LMax()) 
        {
            // after the projection the point is left from the cap surface
            sigtrialIJ = sigtrialIJF1;
            Lproj = L;
        }
        else {
            fIsonCap = true;
            sigtrialIJ[0] = LMax();
            REAL F;
            ComputeF(sigtrialIJ[0],F);
            if (sqJ2 > F) {
                sigtrialIJ[1] = F;
            }
            else
            {
                sigtrialIJ[1] = sqJ2;
            }
            sigtrial.Adjust(sigtrialIJ, sigproj);
            L = LMax();
            Lproj = Lextern;
#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "Projecting on the cap after projection on F1 : sigtrialIJ " << sigtrialIJ << " L " << L << " Lproj " << Lproj;
                LOGPZ_DEBUG(loggerSML, sout.str())
            }
#endif
        }
    }
    else {
        surfaceprojected = 1;
        NewtonF2L(ER, L, sigtrialIJ);
        
        sigtrial.Adjust(sigtrialIJ, sigproj);
        if(sigtrialIJ[0] > LMax())
        {
            fIsonCap = true;
            sigtrialIJ[0] = LMax();
            REAL F;
            ComputeF(sigtrialIJ[0],F);
            if (sqJ2 > F) {
                sigtrialIJ[1] = F;
            }
            else
            {
                sigtrialIJ[1] = sqJ2;
            }
            sigtrial.Adjust(sigtrialIJ, sigproj);
            L = LMax();
            Lproj = Lextern;
#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "Projecting on the cap after projection on F2 : sigtrialIJ " << sigtrialIJ << " L " << L << " Lproj " << Lproj;
                LOGPZ_DEBUG(loggerSML, sout.str())
            }
#endif
        }
        else
        {
            Lproj = L;
        }
    }
//    REAL LF1 = L;
//    REAL LF2 = L;
//    if(I1 < Lextern)
//    {
//        fIsonCap = false;
//        // we can only project on the surface F2
//        if (yield[1] > 0.) {
//            // project the point on surface F2
//            sigtrialIJF2 = sigtrialIJ;
//            TPZTensor<REAL> sigtrialF2;
//            NewtonF2L(ER, LF2, sigtrialIJF2);
//        }
//        sigtrial.Adjust(sigtrialIJF2, sigproj);
////        Lproj = LF2;
//        sigtrialIJ = sigtrialIJF2;
//        L = LF2;
//    }
//    else {
//        // the point can project on F1 or F2 or can be projected on the cap
//        if (yield[0] > 0.) {
//            sigtrialIJF1 = sigtrialIJ;
//            NewtonF1(ER, LF1, sigtrialIJF1);
//
//            // verify if the projected point plastifies F2
//            if (sigtrialIJF1[0] < LF1) {
//                // the projected point is left of L
//                // this means F2 must be active
//                // the increase in L must be at the intersection of both surfaces this intersection point is L itself!
//                // 3 K Depsp/DL (L-Lextern) = (I1(sigtrial)-L)
//                TPZManVector<REAL> sigtrialIntersectIJ(sigtrialIJ);
//                REAL LFIntersect = Lextern;
//                ComputeLatIntersectionLeft(ER, LFIntersect, sigtrialIJ);
//                
//            }
//            else {
//                // verify if the projected point is within F2 
//                TPZManVector<REAL,2> yieldProj(2);
//                TPZTensor<REAL> sigtrialF1;
//                sigtrial.Adjust(sigtrialIJF1,sigtrialF1);
//                Compute(sigtrialF1,LF1,yieldProj,0);
//                if (yieldProj[1] > 0.) {
//                    // the projection on F1 is within F2
//                    // the increase in L is such that the projection point is at the intersection of F1 and F2 to the right
//                    // unlikely case!
//                    // 3 K Depsp/DL (L-Lextern) = (I1(sigtrial) - I1intersect(L) )
//                    TPZManVector<REAL> sigtrialIntersectIJ(sigtrialIJ);
//                    REAL LFIntersect = Lextern;
//                    ComputeLatIntersectionRight(ER, LFIntersect, sigtrialIJ);  
//                    sigtrial.Adjust(sigtrialIJ, sigproj);
//                    L = LFIntersect;
//                }
//                else {
//                    L = LF1;
//                    sigtrialIJ = sigtrialIJF1;
//                    sigtrial.Adjust(sigtrialIJF1, sigproj);
//                }
//            }
//        }
//        else if(yield[1] > 0.)
//        {
//            sigtrialIJF2 = sigtrialIJ;
//            LF2 = Lextern;
//            NewtonF2(ER, LF2, sigtrialIJF2);
//            TPZTensor<REAL> sigtrialF2;
//            sigtrial.Adjust(sigtrialIJF2, sigtrialF2);
//            TPZManVector<REAL,2> yieldproj(2,0.);
//            Compute(sigtrialF2,LF2,yieldproj,0);
//            if (yieldproj[0] > 0.) {
//                // Don't know what to do??
//                DebugStop();
//            }
//            else {
//                L = LF2;
//                sigtrialIJ = sigtrialIJF2;
//            }
//        }
//        if (LF1 > LMax()) {
//            fIsonCap = true;
//        }
//        
//    }
    
    TPZManVector<TPZTensor<REAL>,2> Ndir(2);
    this->N(sigproj, Lproj, Ndir, 1);
    TPZTensor<REAL> sigPlast(sigtrial),epsPlast;
    sigPlast.Add(sigproj, -1.);
    TPZManVector<REAL,2> sigPlastIJ(2);
    sigPlastIJ[0] = sigPlast.I1();
    sigPlastIJ[1] = sqrt(sigPlast.J2());
    ER.ComputeStrain(sigPlast, epsPlast);
    TPZManVector<REAL,2> epsplastIJ(2);
    epsplastIJ[0] = epsPlast.I1();
    epsplastIJ[1] = sqrt(epsPlast.J2());
    REAL I1err = sigPlastIJ[0]-epsplastIJ[0]*3.*ER.K();
    REAL J2err = sigPlastIJ[1]-epsplastIJ[1]*2.*ER.G();
    if (fabs(I1err) > 1.e-10 || fabs(J2err) > 1.e-10) {
        DebugStop();
    }
    if (fIsonCap == false) 
    {
        if (surfaceprojected != 0 && surfaceprojected != 1) {
            DebugStop();
        }
        REAL i1Ndir = Ndir[surfaceprojected].I1();
        REAL sqj2Ndir = sqrt(Ndir[surfaceprojected].J2());
        REAL theta1 = atan2(sqj2Ndir,i1Ndir);
        REAL theta2 = atan2(epsplastIJ[1],epsplastIJ[0]);
        REAL difftheta = theta1-theta2;
        REAL sigplnorm = sigPlast.Norm();
        if (fabs(difftheta) > 1.e-4 && sigplnorm > 1.e-9) {
            std::cout << "difftheta " << difftheta << " theta1 " << theta1 << " theta2 " << theta2 << std::endl;
        }
        REAL theta = 0.;
        if (surfaceprojected == 1) {
            REAL FL;
            ComputeF(Lproj, FL);
            REAL cst = i1Ndir*(FL*fR)/6.;
            REAL sst = sqj2Ndir*FL;
            REAL check = 1.-cst*cst-sst*sst;
            if (fabs(check) > 1.e-6) {
                std::cout << "check = " << check << " cst " << cst << " sst " << sst << std::endl;
            }
            theta = atan2(sst,cst);
            
            REAL xdist = sigtrialIJ[0]-Lproj;
            REAL ydist = sigtrialIJ[1];
            REAL theta3 = atan2(ydist, xdist/fR);
            
            REAL thetadiff = theta-theta3;
            if (fabs(thetadiff) > 1.e-6) {
                std::cout  << " thetadiff " << thetadiff << " theta " << theta << " theta3 " << theta3 << std::endl;
            }
                // verificar pelo angulo
            REAL verify = FuncTheta2L(ER, theta, Lproj, sigtrialIJkeep);
            if (fabs(verify) > 1.e-9) {
                std::cout << "Validity of functheta " << verify << std::endl;
            }
        }
        Lproj = L;
        REAL scale = epsPlast.Norm()/Ndir[surfaceprojected].Norm();
        for (int i=0; i<6; i++) {
            REAL diff = fabs(scale*Ndir[surfaceprojected][i]-epsPlast[i]);
            if (diff > 1.e-6) {
                DebugStop();
            }
        }
        delgamma[0] = 0.;
        delgamma[1] = 0.;
        delgamma[surfaceprojected] = scale;
    }
    else
    {
        REAL lmax = LMax();
        REAL F;
        ComputeF(lmax , F);
        REAL i1Ndir = Ndir[0].I1();
        REAL sqj2Ndir = sqrt(Ndir[1].J2());
        REAL i1epsplast = epsPlast.I1();
        REAL sqj2epsplast = sqrt(epsPlast.J2());
        if (I1 > lmax ) {
            delgamma[0] = i1epsplast/i1Ndir;
        }
        else
        {
            // there would be no reason to fall on the cap
            DebugStop();
        }
        if (sqJ2 > F) {
            delgamma[1] = sqj2epsplast/sqj2Ndir;
        }
        else
        {
            delgamma[1] = 0.;
        }
    }
    Compute(sigproj, Lproj, yield, 0);
#ifdef PZ_LOG
    if(loggerSM.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "After projecting the point yield = " << yield;
        sout << "\ndelgamma = " << delgamma;
        LOGPZ_DEBUG(loggerSM, sout.str())
    }
#endif
    if (yield[0] > 1.e-8 && fIsonCap == true) {
        DebugStop();
    }
    if (yield[1] > 1.e-6) {
        DebugStop();
    }
}

// 3 K Depsp/DL (L-Lextern) = (I1(sigtrial)-L)
inline void TPZYCSandlerDimaggioL::ComputeLatIntersectionLeft(const TPZElasticResponse &ER, REAL &Lextern, TPZVec<REAL> &sigtrialIJ) const
{
    REAL depspdl;
    REAL L = Lextern;
    DEpspDL(L, depspdl);
    REAL I1 = sigtrialIJ[0];
    REAL res = 3.*ER.K()*depspdl*(L-Lextern)-(I1-L);
    while (fabs(res) > 1.e-10) {
        REAL dres;
        REAL d2epspdl2;
        D2EpspDL2(L, d2epspdl2 );
        dres = 3.*ER.K()*depspdl + 1 + 3.*ER.K()*(L-Lextern)*d2epspdl2;
        L -= res/dres;
        DEpspDL(L, depspdl);
        res = 3.*ER.K()*depspdl*(L-Lextern)-(I1-L);
    }
}

// 3 K Depsp/DL (LIntersect-Lextern) = (I1(sigtrial)-I1(Lintersect) )
inline void TPZYCSandlerDimaggioL::ComputeLatIntersectionRight(const TPZElasticResponse &ER, REAL &L, TPZVec<REAL> &sigtrialIJ) const
{
    DebugStop();
}


#endif //TPZYCSANDLERDIMAGGIO_H
