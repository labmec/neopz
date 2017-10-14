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
    void (*fFunc4)(const TPZVec<REAL> &f, int id);
    
    int fPorder;
    int fnSquareElements;
    TPZGeoMesh* fgmesh;
    TPZFMatrix<TVar> fK;
    TPZFMatrix<TVar> fRand_E;
    TPZFMatrix<TVar> fE;
    
public:
	
	/** @brief Class constructor */
	TPZRandomField(TPZGeoMesh* geometricMesh, int numSquareElems) : TPZFunction<TVar>(), fPorder(-1)
    {
        fFunc = 0;
		fFunc2 = 0;
		fFunc3 = 0;
        fFunc4 = 0;
        
        fgmesh = geometricMesh;
        fnSquareElements = numSquareElems; // number of Square Elements
        fK = calcCorrelationMatrix();      // Correlation matrix K
        PrintCorrelation(); // Chama para que possa imprimir automaticamente
        
        // Random Vector
        TPZFMatrix<TVar> Rand_E (1, fnSquareElements, 0.);
        for (int i = 0; i < fnSquareElements; i++) {
            Rand_E(0,i) = arc4random_uniform(3001) + 15300; //uniform
        }
        // std::cout << Rand_E << std::endl; //teste
        fRand_E = Rand_E;
        
        //std::ofstream out_VecE("/Users/batalha/Desktop/fE.txt");
        //fK.Print("E = ", out_VecE, EMathematicaInput);
        
        // Multiplying decomposed Matrix M (U*Sqrt(S)) and random vector fRand_E
        TPZFMatrix<TVar> M = readDecomposedMatrixFromFile();
        fE = fRand_E * M; // Obtem valores de E correlacionados
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
        fFunc4 = 0;
		fPorder = -1;
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf))
    {
		fFunc = 0;
        fFunc2 = FuncPtr;		
        fFunc3 = 0;
        fFunc4 = 0;
		fPorder = -1;
    }
	
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf))
    {
		fFunc = 0;
        fFunc2 = 0;		
        fFunc3 = FuncPtr;
        fFunc4 = 0;
		fPorder = -1;
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &f, int id))
    {
        fFunc = 0;
        fFunc2 = 0;
        fFunc3 = 0;
        fFunc4 = FuncPtr;
        fPorder = -1;
    }
    
    TPZRandomField(const TPZRandomField &cp) : fFunc(cp.fFunc), fFunc2(cp.fFunc2), fFunc3(cp.fFunc3), fFunc4(cp.fFunc4), fPorder(cp.fPorder)
    {
        
    }
    

    TPZRandomField &operator=(const TPZRandomField &cp)
    {
        fFunc = cp.fFunc;
		fFunc2 = cp.fFunc2;
		fFunc3 = cp.fFunc3;
        fFunc4 = cp.fFunc4;
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
    
    // Call this method in the PZMatElasticity 2D (Gaussian Field)
    virtual void Execute(const TPZVec<TVar> &f, int id) {
        f[0] = fE(0, id);
        //std::cout << f[0];
    }
    
    TPZVec<TVar> calcStochasticField(){
        
        // Stochastic Field
        TPZFMatrix<REAL> K = calcCorrelationMatrix();
    
        // TODO - Implementar SVD
        // Correlacionar vetor randomico E[]
        
        return NULL;
    }
    
    
    // Print Correlation Matrix to be decomposed at Mathematica   NANANANANAN pedreiro
    virtual void PrintCorrelation() {
    
        std::ofstream out_kmatrix("KCorr.txt");
        fK.Print("KCorr = ",out_kmatrix,EMathematicaInput);
    }
    
    TPZFMatrix<TVar> readDecomposedMatrixFromFile() {
        TPZFMatrix<TVar> M (fnSquareElements, fnSquareElements, 0.);
        
        // Setar valores de M obtidos do Mathematica (Decomposed Matrix)
        std::ifstream DecMatFile("/Users/batalha/Desktop/decomposed_matrix.tbl");
        std::string line;
        
        int i = 0;
        
        while (std::getline(DecMatFile, line)) {
            REAL value;
            int j = 0;
            std::stringstream ss(line);
            
            while (ss >> value) {
                M(i, j) = value;
                j++;
            }
            i++;
        }
        
        return M;
    }
    
    TPZFMatrix<REAL> calcCorrelationMatrix() {
        
        std::cout << "\nCria matriz da norma entre os centroides (para a matriz de correlacao)" << std::endl;
        
        // Refinamento de elementos selecionados
        REAL e = M_E; // Numero de Euler
        REAL scale = 1.; // Valor de alpha, escala normalizada // variar: 1/4; 1.0; 4.0 // NANANANA pedreiro
        
        TPZFMatrix<REAL> CenterNorm(fnSquareElements, fnSquareElements, 0.0);
        
        TPZManVector<REAL, 3> CenterPoint1, CenterPoint2;
        
        // Elemento analizado
        TPZGeoEl *gel1;
        TPZVec<TPZGeoEl *> sub1;
        TPZManVector<REAL> centerpsi1(3), center1(3);
        
        // Outros elementos
        TPZGeoEl *gel2;
        TPZVec<TPZGeoEl *> sub2;
        TPZManVector<REAL> centerpsi2(3), center2(3);
        
        // Matriz de correlacao
        TPZFMatrix<REAL> KCorr(fnSquareElements, fnSquareElements, 0.0);
        
        // Matriz da distancia entre os centroides
        for (int i = 0; i < fnSquareElements; i++) {
            for (int j = 0; j < fnSquareElements; j++) {
                gel1 = fgmesh->ElementVec()[i];
                gel1->CenterPoint(8, centerpsi1);
                gel1->X(centerpsi1, center1);
                
                CenterPoint1 = center1;
                
                gel2 = fgmesh->ElementVec()[j];
                gel2->CenterPoint(8, centerpsi2);
                gel2->X(centerpsi2, center2);
                
                CenterPoint2 = center2;
                
                //	/*3*/	EQuadrilateral
                if (gel1->Type() == 3 && gel2->Type() == 3) {
                    
                    REAL dx = pow((CenterPoint2[0]-CenterPoint1[0]), 2);
                    REAL dy = pow((CenterPoint2[1]-CenterPoint1[1]), 2);
                    REAL dz = pow((CenterPoint2[2]-CenterPoint1[2]), 2);
                    
                    CenterNorm(i,j) = sqrt(dx + dy + dz);
                    
                    REAL r = CenterNorm(i,j);
                    REAL r2 = pow(r, 2);
                    KCorr(i,j) = pow(e, (-scale * r2));
                }
                
                else {
                    
                    // Verifica se el atual eh quadrilatero
                    std::cout<< "Element Type Error" << std::endl;
    
                }
            }
        }
//        std::cout<< KCorr << std::endl; // teste verifica tam de K
        
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
