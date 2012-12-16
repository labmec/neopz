// $Id: TPZYCSandlerDimaggio.h,v 1.11 2009-06-29 22:54:01 erick Exp $

#ifndef TPZYCSANDLERDIMAGGIO_H
#define TPZYCSANDLERDIMAGGIO_H

#include "TPZTensor.h"
#include "pzfmatrix.h"
#include "TPZElasticResponse.h"
#include "pzlog.h"

#ifndef CHECKCONV
#define CHECKCONV
#include "checkconv.h"
#endif

#include "fadType.h"

#ifdef LOG4CXX
#include "pzlog.h"

static LoggerPtr loggerSM(Logger::getLogger("plasticity.SM"));

#endif

/**
Implementa as funções de potencial plástico e yield criterium do 
modelo constitutivo associativo de Sandler e Dimaggio (1971), desenvolvido
inicialmente para arenitos (Ranch McCormic Sand)
*/
class TPZYCSandlerDimaggio {
    

public:

  enum {NYield = 2};
	
    TPZYCSandlerDimaggio():fA(0.),fB(0.),fC(0.),fD(0.),fW(0.),fR(0.){ }
	
    TPZYCSandlerDimaggio(const TPZYCSandlerDimaggio & source)
    {
		fA = source.fA;
		fB = source.fB;
		fC = source.fC;
		fD = source.fD;
		fW = source.fW;
		fR = source.fR;
    }

    TPZYCSandlerDimaggio & operator=(const TPZYCSandlerDimaggio & source)
    {
		fA = source.fA;
		fB = source.fB;
		fC = source.fC;
		fD = source.fD;
		fW = source.fW;
		fR = source.fR;
		return *this;
    }

	const char * Name() const
    {
	   return "TPZYCSandlerDimaggio";	
    }
	
    void Print(std::ostream & out) const
    {
		out << "\n" << this->Name();
		out << "\n fA = " << fA;
		out << "\n fB = " << fB;
		out << "\n fC = " << fC;
		out << "\n fD = " << fD;
		out << "\n fR = " << fR;
		out << "\n fW = " << fW;
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

    inline void SetUp(REAL A, REAL B, REAL C, REAL D, REAL R, REAL W)
    {
	   fA = A;
	   fB = B;
	   fC = C;
	   fD = D;
	   fR = R;
	   fW = W;
    }
    
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
    void InitialGuess(const TPZElasticResponse &ER, REAL epsp, TPZTensor<REAL> &sigtrial, REAL &epspproj, 
                      TPZVec<REAL> &delgamma, TPZTensor<REAL> &sigproj);

    /**
     * value of x for which F(x)=0
     */
    REAL FZero() const
    {
        return log(fA/fC)/fB;
    }
    
    /**
     * maximum value of L allowed
     */
    REAL LMax() const
    {
        return FZero()-0.001;
    }
private:
   /**
    Solves for the invariant I1 value at the intersection
	of shear and hardening cap yield criteria. 
	The current value of L when entering the funcition is assumed
	to be the initial guess. 
	In this implementation L should be negative in compression.
   */
   template <class T>
   void SolveL(const T & X, T & L, REAL relTol = 1.e-6) const;

    /** compute the derivative of the L function as a function of epsp (A)
     * @param[in] L value of the L function corresponding to A
     * @param[in] value of the volumetric plastic strain
     * @param[out] derivative of L with respect to A
     */
    template <class T>
    void ComputeDL(const T &L, const T &A, T &DL) const;
   /**
    Might be a reasonable initial guess for L when no better data is available. 
	In this implementation L should be negative in compression.
   */
   template <class T>
   void LInitialGuess(const T & X, T & L) const;
	
   /**
    Evaluates the F(L) grouping
   */
   template <class T>
   void ComputeF(const T & L, T & F) const;
	
   /**
    Evaluates the F(L) grouping total derivative
   */
   template <class T>
   void ComputedF(const T & L, T & dF) const;

   /**
    Evaluates X(EpsilonPvol), the value of the first invariant
	of an hidrostatic stress tensor at the cap yield surface.
	In this implementation X should negative in compression.
   */
   template <class T>
   void ComputeX(const T & A, T & X) const;
    
    /**
     * Computes the value of volumetric plastic strain as a function of L
     */
    REAL ComputeEpsp(const REAL L) const
    {
        REAL FL;
        ComputeF(L, FL);
        REAL X = L-fR*FL;
        REAL eps = fW*(exp(fD*X)-1);
        return eps;
    }
    /**
     Compute the value of the equation which equates the evolution of the plastic deformation
     */
    REAL FuncEpsp(const TPZElasticResponse &ER, REAL theta, REAL epsp, REAL delepsp, TPZVec<REAL> &sigtrialIJ)
    {
        REAL X,L;
        ComputeX(epsp+delepsp, X);
        LInitialGuess(X, L);
        SolveL(X, L);
        REAL F;
        ComputeF(L, F);
        const REAL K = ER.fLambda+2.*ER.fMu/3.;
        REAL funcepsp = 3.*K*delepsp-(sigtrialIJ[0]-(L+F*fR*cos(theta)));
        return funcepsp;
    }
    
    /**
     * compute the value of the equation which determines the orthogonality of the projection
     */
    REAL FuncTheta(const TPZElasticResponse &ER, REAL theta, REAL epsp, REAL delepsp, TPZVec<REAL> &sigtrialIJ) const
    {
        REAL I1 = sigtrialIJ[0];
        REAL sqJ2 = sigtrialIJ[1];
        REAL X,L;
        ComputeX(epsp+delepsp, X);
        SolveL(X, L);
        REAL F;
        ComputeF(L, F);
        const REAL K = ER.fLambda+2.*ER.fMu/3.;
        const REAL G = ER.fMu;
        REAL y = 9.*K*(sqJ2-F*sin(theta));
        REAL x = G*fR*(I1-(L+F*fR*cos(theta)));
        REAL res = theta-atan2(y, x);
        return res;
    }
    /**
     * compute the value of the equation which determines the orthogonality of the projection
     */
    REAL FuncTheta2(const TPZElasticResponse &ER, REAL theta, REAL epsp, REAL delepsp, TPZVec<REAL> &sigtrialIJ) const
    {
        REAL I1 = sigtrialIJ[0];
        REAL sqJ2 = sigtrialIJ[1];
        REAL X,L;
        ComputeX(epsp+delepsp, X);
        SolveL(X, L);
        REAL F;
        ComputeF(L, F);
        const REAL K = ER.fLambda+2.*ER.fMu/3.;
        const REAL G = ER.fMu;
        REAL y = (sqJ2-F*sin(theta));
        REAL x = G*fR/(9.*K)*(I1-(L+F*fR*cos(theta)));
        REAL res = x*sin(theta)-y*cos(theta);
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "x = " << x << " y = " << y << " theta = " << theta << " res = " << res;
            LOGPZ_DEBUG(loggerSM, sout.str())
        }
#endif
        return res;
    }
    /**
     * compute the value of the distance which determines the orthogonality of the projection
     */
    REAL DistTheta(const TPZElasticResponse &ER, REAL theta, REAL epsp, REAL delepsp, TPZVec<REAL> &sigtrialIJ) const
    {
        REAL I1 = sigtrialIJ[0];
        REAL sqJ2 = sigtrialIJ[1];
        REAL X,L;
        ComputeX(epsp+delepsp, X);
        SolveL(X, L);
        REAL F;
        ComputeF(L, F);
        const REAL K = ER.fLambda+2.*ER.fMu/3.;
        const REAL G = ER.fMu;
        REAL y = (sqJ2-F*sin(theta));
        REAL x = (I1-(L+F*fR*cos(theta)));
        REAL dist = G*x*x/2.+9.*K*y*y/2.;
        return dist;
    }
    
    /**
     * compute the derivative of the function FuncEpsp with respect to theta and delepsp
     */
    void DFuncEpsp(const TPZElasticResponse &ER, REAL theta, REAL epsp, REAL delepsp, TPZVec<REAL> &result) const
    {
        REAL X,L,DL,F,DF,Dtheta,Depsp;
        ComputeX(epsp+delepsp, X);
        LInitialGuess(X, L);
        SolveL(X, L);
        ComputeDL(L, epsp+delepsp, DL);
        ComputeF(L, F);
        ComputedF(L, DF);
        DF *= DL;
        Dtheta = -F*fR*sin(theta);
        const REAL K = ER.fLambda+2.*ER.fMu/3.;
        Depsp = 3.*K+(DL+DF*fR*cos(theta));
        result[0] = Dtheta;
        result[1] = Depsp;        
    }
    
    /**
     * compute the derivative of the function FuncTheta with respect to theta and delepsp
     */
    void DFuncTheta(const TPZElasticResponse &ER, REAL theta, REAL epsp, REAL delepsp, TPZVec<REAL> &sigtrialIJ, TPZVec<REAL> &result) const
    {
        REAL I1 = sigtrialIJ[0];
        REAL sqJ2 = sigtrialIJ[1];
        REAL X,L,DL,F,DF,Dtheta,Depsp;
        ComputeX(epsp+delepsp, X);
        LInitialGuess(X, L);
        SolveL(X, L);
        ComputeDL(L, epsp+delepsp, DL);
        ComputeF(L, F);
        ComputedF(L, DF);
        DF *= DL;
        const REAL K = ER.fLambda+2.*ER.fMu/3.;
        const REAL G = ER.fMu;
        REAL x = G*fR*(I1-(L+F*fR*cos(theta)));
        REAL y = 9.*K*(sqJ2-F*sin(theta));
        REAL dxepsp = -G*fR*(DL+DF*fR*cos(theta));
        REAL dyepsp = -9.*K*DF*sin(theta);
        REAL dxtheta = G*fR*F*fR*sin(theta);
        REAL dytheta = -9.*K*F*cos(theta);
        REAL denom = x*x+y*y;
        Dtheta = 1-(-y*dxtheta+x*dytheta)/denom;
        Depsp = -(-y*dxepsp+x*dyepsp)/denom;
        result[0] = Dtheta;
        result[1] = Depsp;
    }
    
    /**
     * compute the derivative of the function FuncTheta with respect to theta and delepsp
     */
    void DFuncTheta2(const TPZElasticResponse &ER, REAL theta, REAL epsp, REAL delepsp, TPZVec<REAL> &sigtrialIJ, TPZVec<REAL> &result) const
    {
        REAL I1 = sigtrialIJ[0];
        REAL sqJ2 = sigtrialIJ[1];
        REAL X,L,DL,F,DF,Dtheta,Depsp;
        ComputeX(epsp+delepsp, X);
        LInitialGuess(X, L);
        SolveL(X, L);
        ComputeDL(L, epsp+delepsp, DL);
        ComputeF(L, F);
        ComputedF(L, DF);
        DF *= DL;
        const REAL K = ER.fLambda+2.*ER.fMu/3.;
        const REAL G = ER.fMu;
        REAL x = G*fR/(9.*K)*(I1-(L+F*fR*cos(theta)));
        REAL y = (sqJ2-F*sin(theta));
        REAL dxepsp = -G*fR/(9.*K)*(DL+DF*fR*cos(theta));
        REAL dyepsp = -DF*sin(theta);
        REAL dxtheta = G*fR/(9.*K)*F*fR*sin(theta);
        REAL dytheta = -F*cos(theta);
        Dtheta = dxtheta*sin(theta)-dytheta*cos(theta)+x*cos(theta)+y*sin(theta);
        Depsp = dxepsp*sin(theta)-dyepsp*cos(theta);
        result[0] = Dtheta;
        result[1] = Depsp;
    }
    
    void UpdateSigtrialIJ(const TPZElasticResponse &ER, REAL epsp, REAL theta, TPZVec<REAL> &sigtrialIJ)
    {
        REAL X,L,F;
        ComputeX(epsp, X);
        LInitialGuess(X, L);
        SolveL(X, L);
        ComputeF(L,F);
        sigtrialIJ[0] = L+F*fR*cos(theta);
        sigtrialIJ[1] = F*sin(theta);
    }
    
	void NewtonF2(const TPZElasticResponse &ER, REAL &epsp, TPZVec<REAL> &sigtrialIJ)
    {
        REAL restheta, resdelepsp,disttheta;
        REAL theta = 0.;
        REAL delepsp = 0.;
        disttheta = DistTheta(ER, theta, epsp, delepsp, sigtrialIJ);
        // Look for a best guess for theta
        for (REAL thetaguess=0.; thetaguess <= M_PI; thetaguess += M_PI/20.) {
            REAL distnew = DistTheta(ER, thetaguess, epsp, delepsp, sigtrialIJ);
            if (fabs(distnew) < fabs(disttheta)) {
                theta = thetaguess;
                disttheta = distnew;
            }
        }
        REAL Xini,Lini,Fini;
        ComputeX(epsp, Xini);
        LInitialGuess(Xini, Lini);
        SolveL(Xini, Lini);
        ComputeF(Lini, Fini);
        bool secondquadrant = false;
        if (sigtrialIJ[0] < Lini) {
            secondquadrant = true;
        }
        if (theta < M_PI/2 && secondquadrant) {
            theta = M_PI/2+M_PI/20.;
        }
        if (Lini < -5.) {
            REAL arcs = sigtrialIJ[1]/Fini;
            if (arcs < 1.) {
                theta = asin(arcs);
                if (secondquadrant) {
                    theta = M_PI-theta;
                }
            }
            REAL Lnew = sigtrialIJ[0]-Fini*fR*cos(theta);
            REAL epsnew = ComputeEpsp(Lnew);
            delepsp = epsnew-epsp;
        }
        
        // perform Newton iterations
        restheta = FuncTheta2(ER, theta, epsp, delepsp, sigtrialIJ);
        resdelepsp = FuncEpsp(ER, theta, epsp, delepsp, sigtrialIJ);
        REAL error = sqrt(restheta*restheta+resdelepsp*resdelepsp);
        int count = 0;
        while((fabs(restheta) > 1.e-10 || fabs(resdelepsp) > 1.e-10) && count < 100)
        {
            REAL errprev = error;
            TPZFNMatrix<4,REAL> tangent(2,2);
            TPZFNMatrix<2,REAL> resmat(2,1);
            TPZManVector<REAL,2> tantheta(2,0.), tandelepsp(2,0.);
            DFuncEpsp(ER, theta, epsp, delepsp, tandelepsp);
            DFuncTheta2(ER, theta, epsp, delepsp, sigtrialIJ, tantheta);
            for (int i=0; i<2; i++) {
                tangent(1,i) = tandelepsp[i];
                tangent(0,i) = tantheta[i];
            }
            resmat(0,0) = restheta;
            resmat(1,0) = resdelepsp;
#ifdef LOG4CXX
            {
                std::stringstream sout;
                tangent.Print("tangent matrix",sout);
                resmat.Print("residual",sout);
                LOGPZ_DEBUG(loggerSM, sout.str())
            }
#endif
            std::list<int> singular;
            tangent.SolveDirect(resmat, ELU, singular);
            REAL scale = 1.;
            if (epsp+delepsp-resmat(1,0) < -fW) {
                scale = 0.9999*(fW+epsp+delepsp)/resmat(1,0);
                resmat *= scale;
            }
            REAL thetaprev = theta;
            REAL delepspprev = delepsp;
            theta = thetaprev-resmat(0,0);
            delepsp = delepspprev-resmat(1,0);
            restheta = FuncTheta2(ER, theta, epsp, delepsp, sigtrialIJ);
            resdelepsp = FuncEpsp(ER, theta, epsp, delepsp, sigtrialIJ);
            error = sqrt(restheta*restheta+resdelepsp*resdelepsp);
            int iline = 0;
            while(error > errprev && iline < 5)
            {
                resmat *= 0.5;
                theta = thetaprev-resmat(0,0);
                delepsp = delepspprev-resmat(1,0);
                restheta = FuncTheta2(ER, theta, epsp, delepsp, sigtrialIJ);
                resdelepsp = FuncEpsp(ER, theta, epsp, delepsp, sigtrialIJ);
                error = sqrt(restheta*restheta+resdelepsp*resdelepsp);
                iline++;
            }
            count++;
        }
#ifdef LOG4CXX
        if(count > 10)
        {
            std::stringstream sout;
            sout << "interation count " << count;
            LOGPZ_DEBUG(loggerSM, sout.str())
        }
#endif
        
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "sigtrialIJ input " << sigtrialIJ << " epsp input " << epsp << "\ndelepsp = " << delepsp
            << " theta = " << theta;
            LOGPZ_DEBUG(loggerSM, sout.str())
        }
#endif
        epsp += delepsp;
        UpdateSigtrialIJ(ER, epsp, theta, sigtrialIJ);
    }
public:

		
   /**
     Parameter related to the YC and Plastic Potential
     A, B and C are constants in the modified Drucker-Prager
	 shear yield criteria.
   */
   REAL fA, fB, fC;

   /**
     Parameters related to the YC and Plastic Potential
	 The D and W parameters correlates the total plastic strain to the
	 hydrostatic loading level. It is thus related to the cap
	 hardening/softening.
   */
   REAL fD, fW;

   /**
     Parameter related to the YC and Plastic Potential
     half-axis ratio for the ellipsoidal hardening/softening
	 cap.
   */
   REAL fR;


public:
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
inline void TPZYCSandlerDimaggio::Compute(const TPZTensor<T> & sigma,const T & A, TPZVec<T> &res, int checkForcedYield) const
{
	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	// The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.
	
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << ">>> TPZYCSandlerDimaggio::Compute *** - Plastic Potential / Yield - associative";
          LOGPZ_INFO(logger,sout.str().c_str());
        }
    #endif

	REAL L_REAL, X_REAL;
	ComputeX((REAL)shapeFAD::val(A), X_REAL);
	LInitialGuess(X_REAL, L_REAL);
	SolveL(X_REAL, L_REAL);
    
    T I1 = sigma.I1();
    T J2 = sigma.J2();
	
	// f1 - Modified Drucker-Prager as shear Yield Criterium
	T FI1;
	ComputeF(I1, FI1);
	if(fabs((REAL)shapeFAD::val(J2)) < 1.e-6)
	{
		res[0] = - FI1;	// avoiding nan derivatives
	}else{
		res[0] = sqrt(J2) - FI1;
	}
	
	// f2 - ellipsoidal hardening/softening cap
	T FL, L(L_REAL), X;
   	ComputeX(A, X);
	SolveL(X, L); // evaluating the derivatives of L
	ComputeF(L, FL);
		
	if(fabs( (REAL)shapeFAD::val(FL) ) < 0.00001)
	{
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << "*** TPZYCSandlerDimaggio::ComputePlasticPotential ***";
          sout << "\nDivision by F=" << shapeFAD::val(L) << " at f2 - ellipsoidal hardening/softening cap";
          LOGPZ_WARN(logger,sout.str().c_str());
        }
    #endif
	}
	
	T Temp1( (L - I1)/(FL * T(fR) ) );
	Temp1 *= Temp1;
	T Temp2 = J2 / FL / FL;
		
	res[1] = Temp1 + Temp2 - T(1.);
	
	return;
}

template <class T> 
inline void TPZYCSandlerDimaggio::N(const TPZTensor<T> & sigma, const T & A, TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{

	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	//The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.

	
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << ">>> TPZYCSandlerDimaggio::N *** - Plastification direction - associative";
          LOGPZ_INFO(logger,sout.str().c_str());
        }
    #endif
	
	REAL ResTol = 1.e-6;
	
	REAL L_REAL, X_REAL;
	ComputeX((REAL)shapeFAD::val(A), X_REAL);
	LInitialGuess(X_REAL, L_REAL);
	SolveL(X_REAL, L_REAL, ResTol);
    
    T I1 = sigma.I1();
    T J2 = sigma.J2();
    T SQRTJ2 = sqrt(J2);
	
	{
		//f1 - Modified Drucker-Prager as shear Yield Criterium / Plastic Potential
        REAL fz = FZero();
        T Temp1(0.);
        if (shapeFAD::val(I1) < fz ) 
        {
            Temp1 = I1 * T(fB);
            Temp1 = exp( Temp1 ) * T (fB * fC);
		
            if((REAL)shapeFAD::val(SQRTJ2) < 1.e-6) // just for robustness. f1 shouldn't be reached when J2 = 0.
            {
                #ifdef LOG4CXX_PLASTICITY
                {
                   LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
                   std::stringstream sout;
                   sout << "*** TPZYCSandlerDimaggio::N *** - SQRT(J2) = " << shapeFAD::val(SQRTJ2) <<  " < 1.e-6 causes error in 0-th yield function. Imposing J2 = 1.e-6 instead";
                   LOGPZ_WARN(logger,sout.str().c_str());
                }
                #endif
                SQRTJ2 = T(1.e-6);
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
	}
	
	{//f2 - ellipsoidal hardening/softening cap

		T FL, X, L(L_REAL * 1.- ResTol); // guaranteeing that the function will be evaluated
	   	ComputeX(A, X);
		SolveL(X, L, ResTol); // evaluating the derivatives of L
		
		ComputeF(L, FL);
        // the radius of the ellips needs to be positive
        // this should be taken care of by the computation of L which is limited by LMax() 
        if (shapeFAD::val(FL) <= 0.) {
            DebugStop();
        }
		T FL2 = FL * FL;
		T FL3 = FL2;// / T(2.);
	
		T Temp = (I1-L)/ T(fR * fR) - I1 / T(6.);
		Temp = Temp / FL2 * T(2.);

			#ifdef LOG4CXX_PLASTICITY
            {
               LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
               std::stringstream sout;
               sout << "*** TPZYCSandlerDimaggio::N *** X = " << X
					<< "\n L = " << L << " L_REAL = " << L_REAL
					<< "\n FL = " << FL
					<< "\n Temp = " << Temp;
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
	
    #ifdef LOG4CXX
    {
        LoggerPtr logger(Logger::getLogger("pz.plasticity.SandlerDimaggio.main"));
        std::stringstream sout;
        sout << "<< TPZYCSandlerDimaggio::N *** \n sigma = \n" << sigma
			 << "\nI1 = " << I1 
			 << "\nJ2 = " << J2
			 << "\nSQRTJ2 = " << SQRTJ2
			 << "\nNdir = \n" << Ndir;
        //LOGPZ_DEBUG(logger,sout.str().c_str());
    }
    #endif
	
	return;
}

template <class T> 
inline void TPZYCSandlerDimaggio::H(const TPZTensor<T> & sigma,const T & A, TPZVec<T> & h, int checkForcedYield) const
{

 	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	//The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.
	
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << ">>> TPZYCSandlerDimaggio::H *** - Hardening modulus";
          LOGPZ_INFO(logger,sout.str().c_str());
        }
    #endif
	
	REAL L_REAL, X_REAL;
	ComputeX((REAL)shapeFAD::val(A), X_REAL);
	LInitialGuess(X_REAL, L_REAL);
	SolveL(X_REAL, L_REAL);
    
    T I1 = sigma.I1();
	
	{//f1 - Modified Drucker-Prager as shear Yield Criterium / Plastic Potential

		h[0] = exp (I1 * T(fB) ) * T(3. * fB * fC);
		// h > 0 because plastic deformation is dilatant
	}	

	{//f2 - ellipsoidal hardening/softening cap
		T FL, X, L(L_REAL);
    	ComputeX(A, X);
		SolveL(X, L); // evaluating the derivatives of L
		ComputeF(L, FL);
		T FL2 = FL * FL;
		
		h[1] = (I1 - L) / FL2 * T ( 6. / fR / fR);
		// h <=0 because plastic deformation is compactant
	}
	
	return;
	
}

template <class T>
inline void TPZYCSandlerDimaggio::SolveL(const T & X, T & L, REAL relTol) const
{
	T F, dF, res, dRes;
	
    ComputeF(L, F);
	res = F * T(fR) + X - L;
	
	int i = 0; // ensuring evaluating at least once the Newton method
	// so that the evaluation with FAD Type and converged L evaluates the
	// derivatives
	while(
		  fabs( (REAL)shapeFAD::val(res) ) /
		  max( (REAL)fabs(shapeFAD::val(X) ), (REAL).000001 ) // avoiding division by zero
		  > relTol ||
		  i<1)
	{
        i++;
		
        ComputedF(L,dF);
        dRes = dF * T(fR) - T(1.);
		
        L -= res/dRes;
																	  
	    ComputeF(L, F);
	    res = F * T(fR) + X - L;
	}
    if (shapeFAD::val(L) > LMax()) {
        L = T(LMax());
    }
	
}

template <class T>
inline void TPZYCSandlerDimaggio::LInitialGuess(const T & X, T & L) const
{
	T FAprox;
	ComputeF(X/T(2.), FAprox);	
	L = X + FAprox * T(fR);
}

/** compute the derivative of the L function as a function of epsp (A)
 * @param[in] L value of the L function corresponding to A
 * @param[in] value of the volumetric plastic strain
 * @param[out] derivative of L with respect to A
 */
template <class T>
inline void TPZYCSandlerDimaggio::ComputeDL(const T &L, const T &A, T &DL) const
{
    REAL LMx = LMax();
    if (shapeFAD::val(L) >= LMx) {
        DL = T(0.);
    }
    else 
    {
        DL = T(1.)/(T(fD)*(A+T(fW))*(T(1.)+T(fR*fB*fC)*exp(T(fB)*L)));
    }
}

template <class T>
inline void TPZYCSandlerDimaggio::ComputeF(const T & L, T & F) const
{
	F = L * T(fB);
	F = exp(F) * T(fC);
	F = T(fA) - F;	
}

template <class T>
inline void TPZYCSandlerDimaggio::ComputedF(const T & L, T & dF) const
{
	dF = L * T(fB);
	dF = exp(dF) * T(-fC*fB);	
}

template <class T>
inline void TPZYCSandlerDimaggio::ComputeX(const T & A, T & X) const
{
	REAL ep = - 0.999999 * fW;
	
	if(shapeFAD::val(A) < ep)
	{
		T dXdep = T(fW / (ep + fW) / fD);
		X = T( log( ep/fW + 1. ) / fD );
		X = X + dXdep * (A - T(ep));
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << "*** TPZYCSandlerDimaggio::ComputeX *** ##### Changing hardening equation to adjusted linear - Excessive volumetric plastic strain - Please check aftwerwards if alpha = epsP.I1() to verify if results are consistent! ######";
          LOGPZ_WARN(logger,sout.str().c_str());
        }
    #endif
	}else{
		X = A / T(fW) + T(1.);
        if (shapeFAD::val(X) <= 0.) {
            DebugStop();
        }
		X = log( X ) / T(fD);
	}
	
}

//////////////////CheckConv related methods/////////////////////

inline int TPZYCSandlerDimaggio::NumCases()
{
    return 4;
}

inline void TPZYCSandlerDimaggio::LoadState(TPZFMatrix<REAL> &state)
{

  int i;
  for(i=0; i<6; i++) gRefTension.fData[i] = state(i,0);

#ifdef LOG4CXX_PLASTICITY
  {
  LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));

  std::stringstream sout;
  sout << ">>> LoadState *** Tension " << state;
  LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif
}

inline void TPZYCSandlerDimaggio::ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &, int icase)
{

  const int nVars = 6;
  typedef TFad<nVars,REAL> TFAD;

  int i, j;
  TPZVec< TPZTensor < REAL > > N_Dir(2);
  REAL A = -0.05; // example epsVP value to be used with checkconv
  TPZTensor < TFAD > Sigma_FAD;
  TFAD A_FAD(A);
  TPZVec< TPZTensor < TFAD > > N_Dir_FAD(2);

  switch(icase)
  {
    case 0:
	case 1:
      //Compute N
      tangent.Redim(1,nVars);
      N(gRefTension, A, N_Dir, 0);
      for(i=0; i<nVars; i++)
          tangent(0,i) = N_Dir[icase].fData[i];
      break;
    case 2:
	case 3:
      //Compute derivatives of N
      tangent.Redim(nVars,nVars);
      gRefTension.CopyTo(Sigma_FAD);
      for(i = 0;i < nVars; i++)
         Sigma_FAD.fData[i].diff(i,nVars);
      N(Sigma_FAD, A_FAD, N_Dir_FAD, 0);
      for(i = 0; i < nVars; i++)
        for(j = 0; j < nVars; j++)
          tangent(i,j) = N_Dir_FAD[icase-2].fData[i].dx(j);
    break;

  }
#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
  std::stringstream sout;
  sout << ">>> ComputeTangent *** " << tangent;
  LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
}

inline void TPZYCSandlerDimaggio::Residual(TPZFMatrix<REAL> &res,int icase)
{

  int i;
  const int nVars = 6;
  typedef TFad<nVars,REAL> TFAD;

  TPZVec< REAL > PlasticPot(2);
  REAL A = -0.05; // example epsVP value to be used with checkconv
  TPZVec< TPZTensor<REAL> > N_Dir(2);
  switch(icase)
  {
    case 0:
	case 1:
      //Compute PlasticPotential
      res.Redim(1,1);
      Compute(gRefTension, A, PlasticPot, 0);
      res(0,0) = PlasticPot[icase];
    break;
    case 2:
	case 3:
      //Compute Ndir
      res.Redim(nVars,1);
      N(gRefTension, A, N_Dir,0);
      for(i = 0; i < nVars; i++)
         res(i,0) = N_Dir[icase-2].fData[i];
    break;
  }
	
#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
  std::stringstream sout;
  sout << ">>> Residual *** " << res;
  LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
  
}

inline void TPZYCSandlerDimaggio::CheckConv()
{

   const int nVars = 6;

   // Creating the Sandler Dimaggio obejct
   TPZYCSandlerDimaggio YCSandlerDimaggio;

   TPZYCSandlerDimaggio::McCormicRanchSand(YCSandlerDimaggio);

   REAL multipl = 1.;
	
   TPZFNMatrix<nVars> sigma(nVars,1), Range(nVars,1);
   sigma(_XX_,0) = -0.17*multipl;
   sigma(_YY_,0) = -0.13*multipl;
   sigma(_ZZ_,0) = -0.11*multipl;
   sigma(_XY_,0) = -0.7 *multipl;
   sigma(_XZ_,0) = -0.5 *multipl;
   sigma(_ZZ_,0) = -0.3 *multipl;

   Range = sigma * (1./19.);
   TPZVec< REAL > Coefs(1,1.);
   CheckConvergence(YCSandlerDimaggio, sigma, Range, Coefs);
   
}

//////////////////CheckConv related methods/////////////////////

//////////////////Internal routines verification/////////////////

inline void TPZYCSandlerDimaggio::TestSolveL()
{
   #ifdef LOG4CXX_PLASTICITY
   LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
   {
     std::stringstream sout;
     sout << "<<< TPZYCSandlerDimaggio::TestSolveL ***";
     LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
   }
   #endif
	
   const int oneVar = 1;
   typedef TFad<oneVar, REAL> TFAD_ONE;
	
   // Creating the Sandler Dimaggio obejct
   TPZYCSandlerDimaggio YCSandlerDimaggio;

   TPZYCSandlerDimaggio::McCormicRanchSand(YCSandlerDimaggio);

   // Verifying if the SolveL routines are working
   // It needs a correct evaluation of both F and dF
	
   //Supposing a plastic volumetric strain
   REAL X, L, dF, epsVP = -.05;
   TFAD_ONE L_FAD, F_FAD;
   YCSandlerDimaggio.ComputeX(epsVP, X);

   YCSandlerDimaggio.LInitialGuess(X, L);
   L_FAD = L;
   L_FAD.diff(0, 0);
	
   YCSandlerDimaggio.ComputeF(L_FAD, F_FAD);
   YCSandlerDimaggio.ComputedF((REAL)shapeFAD::val(L_FAD), dF);
	
   #ifdef LOG4CXX_PLASTICITY
   {
     std::stringstream sout;
     sout << "*** TPZYCSandlerDimaggio::TestSolveL ***";
     sout << "\nDerivative evaluated with FAD (dF/dL)FAD = " << F_FAD.dx(0);
     sout << "\nDerivative evaluated explicitly (dF/dL) =  " << dF;
     sout << "\nProposed L = " << L;
     LOGPZ_DEBUG(loggerSandlerDimaggio,sout.str().c_str());
   }
   #endif
  
   YCSandlerDimaggio.SolveL(X, L);
   #ifdef LOG4CXX_PLASTICITY
   {
     std::stringstream sout;
     sout << "*** TPZYCSandlerDimaggio::TestSolveL ***";
     sout << "\nSolved   L = " << L;
     LOGPZ_DEBUG(loggerSandlerDimaggio,sout.str().c_str());
   }
   #endif
		
   //REAL multipl = 1.; // testing the shear yield criterium (f1)
   REAL multipl = 10.; // testing the elipsoidal compressive cap (f2)
	
   // Checking if NDir:(1,1,1,0,0,0) equals H
   // verifying if the TFAD derivatives of N equal the Ndir vector.
   const int sixVars = 6;
   typedef TFad<sixVars, REAL> TFAD_SIX;

   TPZTensor<REAL> sigma;
   TPZVec<TPZTensor<REAL> > Ndir(2);
   TPZVec<REAL> h(2);
   TPZTensor<TFAD_SIX> sigma_FAD;
   TFAD_SIX epsVP_FAD(epsVP);
   TPZVec< TFAD_SIX > PlasticPot(2);
	
   sigma.XX() = -0.17*multipl;
   sigma.YY() = -0.13*multipl;
   sigma.ZZ() = -0.11*multipl;
   sigma.XY() = -0.7 *multipl;
   sigma.YZ() = -0.5 *multipl;
   sigma.XZ() = -0.3 *multipl;
	
   sigma_FAD.XX() = sigma.XX();
   sigma_FAD.XX().diff(0,sixVars);
	
   sigma_FAD.YY() = sigma.YY();
   sigma_FAD.YY().diff(3,sixVars);
	
   sigma_FAD.ZZ() = sigma.ZZ();
   sigma_FAD.ZZ().diff(5,sixVars);
	
   sigma_FAD.XY() = sigma.XY();
   sigma_FAD.XY().diff(1,sixVars);
	
   sigma_FAD.YZ() = sigma.YZ();
   sigma_FAD.YZ().diff(4,sixVars);
	
   sigma_FAD.XZ() = sigma.XZ();
   sigma_FAD.XZ().diff(2,sixVars);
	
   YCSandlerDimaggio.N(sigma, epsVP, Ndir, 0);
   YCSandlerDimaggio.H(sigma, epsVP, h, 0);
   YCSandlerDimaggio.Compute(sigma_FAD, epsVP_FAD, PlasticPot, 0);

   #ifdef LOG4CXX_PLASTICITY
   {
	 std::stringstream sout;
     sout << "*** TPZYCSandlerDimaggio::TestSolveL ***";
     sout << "\nVerifying if H equals depsVP/dSigma = Ndir:(1 1 1 0 0 0)" ;
	
     for(int i = 0; i < NYield; i++)
     {

	    sout << "\n" << i << "-th Yield Criterium";
        sout << "\nepsVP calculated from Ndir = " << Ndir[i].XX()+Ndir[i].YY()+Ndir[i].ZZ();
        sout << "\nepsVP calculated from H =    " << h[i];
        sout << "\nVerifying if dPlasticPot/dSigma equals Ndir" ;
	    sout << "\nNdir evaluated explicitly with function N():" << Ndir[i];
        sout << "\nN evaluated through FAD evaluations: dPlasticPot/dsigma:" << PlasticPot[i];

	 }

	 LOGPZ_DEBUG(loggerSandlerDimaggio,sout.str().c_str()); 
   }
   #endif
	
}
//////////////////Internal routines verification/////////////////


inline void TPZYCSandlerDimaggio::McCormicRanchSand(TPZYCSandlerDimaggio & material)
{
   #ifdef LOG4CXX_PLASTICITY
   LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
   {
      std::stringstream sout;
      sout << ">>> TPZYCSandlerDimaggio::McCormicRanchSand ***";
      LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
   }
   #endif
   // setup with data from McCormic Ranch Sand
   // Dimaggio, Frank L. Sandler, Ivan S. Material model for granular soils
   // J. of the Eng. Mech. Div. vol. 97 n0 EM3 
   // pp 935-949,1971
   // OBS: stresses in ksi
   REAL A = 0.25,
        B = 0.67,
        C = 0.18,
        D = 0.67,
        R = 2.5,
        W = 0.066;
	
   material.SetUp(A, B, C, D, R, W);
}

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
inline void TPZYCSandlerDimaggio::InitialGuess(const TPZElasticResponse &ER, REAL epsp, TPZTensor<REAL> &sigtrial, REAL &epspproj, 
                  TPZVec<REAL> &delgamma, TPZTensor<REAL> &sigproj)
{
    TPZManVector<REAL,2> yield(2,0.);
    Compute(sigtrial, epsp, yield, 0);
    if (yield[0] <= 0. && yield[1] <= 0.) {
        epspproj = epsp;
        sigproj = sigtrial;
        delgamma.Fill(0.);
#ifdef LOG4CXX
        if(loggerSM->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Deformation elastic yield = " << yield;
            sout << "delgamma condition " << (delgamma[0] > 0.) << " " << (delgamma[1] > 0.) << std::endl;
            LOGPZ_DEBUG(loggerSM, sout.str())
        }
#endif
        return;
    }
    TPZTensor<REAL> S;
    sigtrial.S(S);
    REAL J2 = S.J2();
    REAL sqJ2 = sqrt(fabs(J2));
    REAL I1 = sigtrial.I1();
    if (yield[1] > 0.) {
        // project the point on surface F2
        TPZManVector<REAL,2> sigtrialIJ(2);
        sigtrialIJ[0] = I1;
        sigtrialIJ[1] = sqJ2;
        NewtonF2(ER, epsp, sigtrialIJ);
        REAL j2Scale = 1;
        if (sqJ2 > 1.e-6) {
            j2Scale = sigtrialIJ[1]/sqJ2;
        }
        sigproj.Identity();
        sigproj.Multiply(1./3., sigtrialIJ[0]);
        sigproj.Add(S, j2Scale);
        epspproj = epsp;
        TPZManVector<TPZTensor<REAL>,2> Ndir(2);
        this->N(sigproj, epspproj, Ndir, 1);
        TPZTensor<REAL> sigPlast(sigtrial),epsPlast;
        sigPlast.Add(sigproj, -1.);
        ER.ComputeDeformation(sigPlast, epsPlast);
        REAL scale = epsPlast.Norm()/Ndir[1].Norm();
        for (int i=0; i<6; i++) {
            REAL diff = fabs(scale*Ndir[1][i]-epsPlast[i]);
            if (diff > 1.e-6) {
                DebugStop();
            }
        }
        delgamma[0] = 0.;
        delgamma[1] = scale;
    }
    Compute(sigproj, epspproj, yield, 0);
#ifdef LOG4CXX
    if(loggerSM->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "After projecting the point yield = " << yield;
        sout << "\ndelgamma = " << delgamma;
        LOGPZ_DEBUG(loggerSM, sout.str())
    }
#endif
    if (yield[0] > 0.) {
        ;
//        DebugStop();
    }
    if (yield[1] > 1.e-6) {
        DebugStop();
    }
    // falta calcular delgamma
}


#endif //TPZYCSANDLERDIMAGGIO_H
