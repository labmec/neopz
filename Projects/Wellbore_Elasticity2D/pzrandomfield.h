#ifndef PZRANDOMFIELD_H
#define PZRANDOMFIELD_H

#include "pzfunction.h"
#include <iostream>
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzerror.h"
#include "tpzverysparsematrix.h"
#include "pzsfulmat.h"

#include <math.h>
#include <complex>
#include <string>

template<class TVar>
class TPZRandomField : public TPZFunction<TVar>
{
	
    void (*fFunc)(const TPZVec<REAL> &x, TPZVec<TVar> &f);
    void (*fFunc2)(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    void (*fFunc3)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    
    int fPorder;
    int nSquareElements;
    TPZGeoMesh* gmesh;
    
public:
	
	/** @brief Class constructor */
	TPZRandomField(TPZGeoMesh* geometricMesh, int numSquareElems) : TPZFunction<TVar>(), fPorder(-1)
    {
        fFunc = 0;
		fFunc2 = 0;
		fFunc3 = 0;
        gmesh = geometricMesh;
        nSquareElements = numSquareElems;
    }
	
	/** @brief Class destructor */
	virtual ~TPZRandomField()
    {
        
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val))
    {
        fFunc = FuncPtr;
		fFunc2 = 0;
		fFunc3 = 0;
		fPorder = -1;
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf))
    {
		fFunc = 0;
        fFunc2 = FuncPtr;		
		fFunc3 = 0;
		fPorder = -1;
    }
	
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf))
    {
		fFunc = 0;
        fFunc2 = 0;		
		fFunc3 = FuncPtr;
		fPorder = -1;
    }	
    
    TPZRandomField(const TPZRandomField &cp) : fFunc(cp.fFunc), fFunc2(cp.fFunc2), fFunc3(cp.fFunc3), fPorder(cp.fPorder)
    {
        
    }
    

    TPZRandomField &operator=(const TPZRandomField &cp)
    {
        fFunc = cp.fFunc;
		fFunc2 = cp.fFunc2;
		fFunc3 = cp.fFunc3;
        fPorder = cp.fPorder;
        return *this;
    }
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
    {
        if (!fFunc2) {
			DebugStop();
		}
        fFunc2(x, f, df);
    }
	
	/**
	 * @brief Performs time dependent function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param ftime  time to evaluate	 
	 * @param f function values
	 * @param gradf function derivatives
	 */	
	virtual void Execute(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf)
    {
        if (!fFunc3) {
			DebugStop();
		}
        fFunc3(x, ftime, f, gradf);
    }	
    
	/**
	 * @brief Execute method receiving axes. It is used in shape functions
	 * @note NOT IMPLEMENTED
	 */
	virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<TVar> &f, TPZFMatrix<TVar> &df){
        DebugStop();
    }
    
    /**
	 * @brief Simpler version of Execute method which does not compute function derivatives 
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f) {
        /* The arc4random() function uses the key stream generator employed by the arc4 cipher,
         * which uses 8*8 8 bit S-Boxes.  The S-Boxes can be in about (2**1700) states.  The
         * arc4random() function returns pseudo-randomnumbers in the range of 0 to (2**32)-1, and
         * therefore has twice the range of rand(3) and random(3).
         * arc4random_uniform() will return a uniformly distributed random number less than
         * upper_bound. arc4random_uniform() is recommended over constructions like
         * ``arc4random() % upper_bound'' as it avoids "modulo bias" when the upper bound is not a
         * power of two.
         */
		
        //REAL E = rand() % 3000 + 15300 + 1; // not uniform
        REAL E = arc4random_uniform(3001) + 15300; //uniform
		f[0] = E;
        
        // TODO - chamar calcStochasticField para obter vetor E[] correlacionado
        
    }
    
    TPZVec<TVar> calcStochasticField(){
        
        // Stochastic Field
        TPZFMatrix<REAL> K = calcCorrelationMatrix();
    
        // TODO - Implementar SVD
        // Correlacionar vetor randomico E[]    
        return NULL;
    }
    
    TPZFMatrix<REAL> calcCorrelationMatrix() {
        std::cout << "\nCria matriz da norma entre os centroides (para a matriz de correlacao)" << std::endl;
        
        // Refinamento de elementos selecionados
        REAL e = M_E; // Numero de Euler
        REAL scale = 1.; // Valor de alpha, escala normalizada
        
        TPZFMatrix<REAL> CenterNorm(nSquareElements, nSquareElements, 0.0);
        
        TPZManVector<REAL, 3> CenterPoint1, CenterPoint2;
        
        // Matriz da distancia entre os centroides
        TPZGeoEl *gel;
        TPZVec<TPZGeoEl *> sub;
        TPZManVector<REAL> centerpsi(3), center(3);
        
        // Matriz de correlacao
        TPZFMatrix<REAL> KCorr(nSquareElements, nSquareElements, 0.0);
        
        std::cout << "\nMatriz de correlacao" << std::endl;
        for (int i = 0; i < nSquareElements; i++) {
            for (int j = 0; j < nSquareElements; j++) {
                gel = gmesh->ElementVec()[i];
                gel->CenterPoint(8, centerpsi);
                gel->X(centerpsi, center);
                
                CenterPoint1 = center;
                
                gel = gmesh->ElementVec()[j];
                gel->CenterPoint(8, centerpsi);
                gel->X(centerpsi, center);
                
                CenterPoint2 = center;
                
                REAL dx = pow((CenterPoint2[0]-CenterPoint1[0]), 2);
                REAL dy = pow((CenterPoint2[1]-CenterPoint1[1]), 2);
                REAL dz = pow((CenterPoint2[2]-CenterPoint1[2]), 2);
                
                CenterNorm(i,j) = sqrt(dx + dy + dz);
                
                REAL r = CenterNorm(i,j);
                REAL r2 = pow(r, 2);
                KCorr(i,j) = pow(e, (-scale * r2));
            }
        }
        
        return KCorr;
    }
    
    /** @brief Returns number of functions. */
	virtual int NFunctions()
    {
        return 1;
    }
    
    void SetPolynomialOrder(int porder)
    {
        fPorder = porder;
    }
	
	/** @brief Polynomial order of this function. */
	/** In case of non-polynomial function it can be a reasonable approximation order. */
	virtual int PolynomialOrder() 
    {
#ifdef DEBUG
        if (fPorder == -1) {
            DebugStop();
        }
#endif
        return fPorder;
    }
	
	/** @brief Unique identifier for serialization purposes */
	virtual int ClassId() const
    {
        return -1;
    }
	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid)
    {
//        DebugStop();
        TPZSaveable::Write(buf,withclassid);
    }
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context)
    {
        DebugStop();
    }
};

#endif
