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
#include <random>

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
    TPZFMatrix<TVar> fRand_U;
    TPZFMatrix<TVar> fU;
    
public:
	
	/** @brief Class constructor */
	TPZRandomField(TPZGeoMesh* geometricMesh, int numSquareElems) : TPZFunction<TVar>(), fPorder(-1)
    {
        fFunc  = 0;
		fFunc2 = 0;
		fFunc3 = 0;
        fFunc4 = 0;
        
        fgmesh = geometricMesh;
        fnSquareElements = numSquareElems;  // number of Square Elements
        fK = calcCorrelationMatrix();       // Correlation matrix K
        PrintCorrelation();                 // Exporta KCoor .txt
        
        // Random Vector U - Normal Distribution
        TPZFMatrix<TVar> Rand_U (fnSquareElements, 1, 0.);
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.,1.0);
        for (int i = 0; i < fnSquareElements; i++) {
            Rand_U(i,0) = distribution(generator); //uniform
        }
        fRand_U = Rand_U;
        
        // Multiplying decomposed Matrix M (U*Sqrt(S)) and random normal vector fRand_U
        TPZFMatrix<TVar> M = readDecomposedMatrixFromFile();
        fU = M * fRand_U; // Obtem valores correlacionados
        
        
//        // Exporta vetor randomico para validar no mathematica
//        std::ofstream out_VecRand_U("/Users/batalha/Desktop/Rand_U.txt");
//        Rand_U.Print("ERand = ", out_VecRand_U, EMathematicaInput);
//        
//        // Exporta vetor correlacionado para validar no mathematica
//        std::ofstream out_M("/Users/batalha/Desktop/M.txt");
//        M.Print("M = ", out_M, EMathematicaInput);
//        
//        // Exporta vetor correlacionado para validar no mathematica
//        std::ofstream out_VecUCorr("/Users/batalha/Desktop/fUCorr.txt");
//        fU.Print("ECorr = ", out_VecUCorr, EMathematicaInput);
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
    // This method gives random values of Young Modulus from a specific range
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f) {
        //REAL E = rand() % 3000 + 15300 + 1; // not uniform
        REAL E = arc4random_uniform(3001) + 15300; //uniform distribution
		f[0]   = E;
        
    }
    
    // Call this method in the PZMatElasticity 2D (Gaussian Field)
    virtual void Execute(const TPZVec<TVar> &f, int id) {
        f[0] = fU(id, 0); // gives to the TPZMaterial the correlated variable of the element
    }
    
    // Calc Correlation Matrix
    TPZVec<TVar> calcStochasticField(){
        
        // Stochastic Field
        TPZFMatrix<REAL> K = calcCorrelationMatrix();
        
        return NULL;
    }
    
    
    // Print Correlation Matrix to be decomposed at Mathematica   
    virtual void PrintCorrelation() {
    
        std::ofstream out_kmatrix("KCorr.txt");
        fK.Print("KCorr = ",out_kmatrix,EMathematicaInput);
    }
    
    // Read Decomposed Matrix from File
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
    
    // Calcula Correlation Matrix
    TPZFMatrix<REAL> calcCorrelationMatrix() {
        
        std::cout << "\nCria matriz da norma entre os centroides (para a matriz de correlacao)" << std::endl;
        
        // Refinamento de elementos selecionados
        REAL e = M_E; // Numero de Euler
        REAL scale = 0.0; // Valor de alpha, escala normalizada // variar: 1/4; 1.0; 4.0
        
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
        
        /* Checking correlation Matrix at Matlab */
        // Matriz de coordenadas
//        TPZFMatrix<REAL> Coordinates(fnSquareElements, 4, 0.0);
//        TPZGeoEl *gel;
//        TPZManVector<REAL> centerpsi(3), center(3);
//        TPZManVector<REAL, 3> CenterPoint;
//        
//        for (int i = 0; i < fnSquareElements; i++) {
//                gel = fgmesh->ElementVec()[i];
//                gel->CenterPoint(8, centerpsi);
//                gel->X(centerpsi, center);
//                
//                CenterPoint = center;
//            
//            //Coordinates
//            REAL xx = CenterPoint[0];
//            REAL yy = CenterPoint[1];
//            REAL zz = CenterPoint[2];
//            
//            Coordinates(i, 0) = i+1;
//            Coordinates(i, 1) = xx;
//            Coordinates(i, 2) = yy;
//            Coordinates(i, 3) = zz;
//
//        }
//        //std::cout << Coordinates << std::endl;
//        std::ofstream out_Coordinates("Coordinates.txt");
//        Coordinates.Print("XYZ = ",out_Coordinates,EMathematicaInput);
        
        
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
