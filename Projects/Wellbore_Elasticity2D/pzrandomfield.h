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
#include <chrono>


template<class TVar>
class TPZRandomField : public TPZFunction<TVar>
{
	
    void (*fFunc)(const TPZVec<REAL> &x, TPZVec<TVar> &f);
    void (*fFunc2)(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    void (*fFunc3)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    void (*fFunc4)(const TPZVec<REAL> &f, int id);
    
    int fPorder;
    int fnSquareElements;
    int fmatsize;
    int fstochasticInclined;
    REAL fdirection;
    REAL finclination;
    REAL frw; // wellbore radius
    REAL frext; // external radius
    REAL fH; // cylinder total height
    REAL fh; // elements height (square elements for now)
    TPZGeoMesh* fgmesh; //geometric mesh
    TPZFMatrix<TVar> fK; // correlation matrix
    TPZFMatrix<TVar> fRand_U; //random distribution
    TPZFMatrix<TVar> fU; // random correlated* distribution
    TPZFMatrix<TVar> fM; //Decomposed matrix from Mathematica
    bool fexponential; //exponential function
    bool fspheric; //spheric function
    bool fnormDistribution; // normal distribution
    bool flognormDistribution; // lognormal distribution
    
public:
	
    /** @brief Class constructor */
    TPZRandomField(TPZGeoMesh* geometricMesh, int numSquareElems, int stochasticInclined, REAL direction,
    REAL inclination, REAL rw, REAL rext, const TPZFMatrix<TVar> &M) : TPZFunction<TVar>(), fPorder(-1)
    {
        fFunc  = 0;
        fFunc2 = 0;
        fFunc3 = 0;
        fFunc4 = 0;
        
        fgmesh = geometricMesh;
        fnSquareElements = numSquareElems;  // number of Square Elements
        fstochasticInclined = stochasticInclined;
        fdirection = direction;
        finclination = inclination;
        fnormDistribution = true;
        flognormDistribution = false;
        
        // this should not be here, try to get from fgmesh
        frw = rw;
        frext = rext;
        fM = M;
        
        CalculateStochasticField(); //This should be called by the object

    }
    
    virtual void CalculateStochasticField(){
        
        SetFieldDistribution(fM, fnormDistribution, flognormDistribution); //This should be called by the object
        
        SetFieldGeometry(fgmesh,fnSquareElements,frw, frext); //This should be called by the object
        
        if (fstochasticInclined == 1) {
            
            fK = calcCorrelationMatrixInclined();  // Correlation matrix K
            
            PrintCorrelation();                    // Exporta KCoor .txt
            
            // Create function to decompose fK using SVD decomposition
            
            GetDistribution(fmatsize);             // Get random distribution
            
            // Multiplying decomposed Matrix M (U*Sqrt(S)) and random normal vector fRand_U
            fU = fM * fRand_U; // Get correlated random distribution
            // In this function fM should be replaced by the left singular vetor U and the square root of the diagonal matrix S, then multiply by fRand_U
        }
        else{
            
            fK = calcCorrelationMatrix();       // Correlation matrix K
            
            PrintCorrelation();                 // Exporta KCoor .txt
            
            // Create function to decompose fK using SVD decomposition
            
            GetDistribution(fnSquareElements);  // Get random distribution
            
            // Multiplying decomposed Matrix M (U*Sqrt(S)) and random normal vector fRand_U
            fU = fM * fRand_U; // Get correlated random distribution
            // In this function fM should be replaced by the left singular vetor U and the square root of the diagonal matrix S, then multiply by fRand_U
        }
        
    }
    
    virtual void SetInclinedField(int stochasticInclined,REAL direction, REAL inclination)
    {
        fstochasticInclined = stochasticInclined;
        fdirection = direction;
        finclination = inclination;
    }
    
    virtual void SetFieldGeometry(TPZGeoMesh* geometricMesh,int numSquareElems,REAL rw, REAL rext)
    {
        
        fgmesh = geometricMesh;
        fnSquareElements = numSquareElems;  // number of Square Elements
        frw = rw;
        frext = rext;

        if (fstochasticInclined==true){
            InclinedFieldGeometry();
        }
    }
    
    virtual void InclinedFieldGeometry(){
        int nLayers = 8;
        fH = 2 * frext; // altura total do cilindro em metros
        fh = fH / nLayers; // altura de cada cubo (elemento) em metros
        fmatsize = fnSquareElements * (fH/fh) + fnSquareElements;
        
    }
    
    virtual void SetFieldDistribution(const TPZFMatrix<TVar> &M, bool normDistrib, bool lognormDistrib){
        fnormDistribution = normDistrib;
        flognormDistribution = lognormDistrib;
        fM = M; //This should be removed later
    }
    
    
    virtual void EvaluateCorrelation(int stochasticInclined)
    {
        if (stochasticInclined == 1) {
            
            int nLayers = 8;
            fH = 2 * frext; // altura total do cilindro em metros
            fh = fH / nLayers; // altura de cada cubo (elemento) em metros
            fmatsize = fnSquareElements * (fH/fh) + fnSquareElements;
            
            fK = calcCorrelationMatrixInclined();  // Correlation matrix K
            
            PrintCorrelation();                    // Exporta KCoor .txt
            
            // Create function to decompose fK using SVD decomposition
            
            GetDistribution(fmatsize);             // Get random distribution
            
            // Multiplying decomposed Matrix M (U*Sqrt(S)) and random normal vector fRand_U
            fU = fM * fRand_U; // Get correlated random distribution
            // In this function fM should be replaced by the left singular vetor U and the square root of the diagonal matrix S, then multiply by fRand_U
        }
        else{
            
            fK = calcCorrelationMatrix();       // Correlation matrix K
            
            PrintCorrelation();                 // Exporta KCoor .txt
            
            // Create function to decompose fK using SVD decomposition
            
            GetDistribution(fnSquareElements);  // Get random distribution
            
            // Multiplying decomposed Matrix M (U*Sqrt(S)) and random normal vector fRand_U
            fU = fM * fRand_U; // Get correlated random distribution
            // In this function fM should be replaced by the left singular vetor U and the square root of the diagonal matrix S, then multiply by fRand_U
        }
    }
    
    TPZFMatrix<TVar> GetDistribution(int matrixSize)
    {
        //if(fnormDistribution==true) {
        
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution(0.,1.0);
        
        // Random Vector U
        TPZFMatrix<TVar> Rand_U (matrixSize, 1, 0.);
        
        for (int i = 0; i < matrixSize; i++) {
            Rand_U(i,0) = distribution(generator);
            distribution.reset();
            fRand_U = Rand_U;
        }
        //}
        
//        else if(flognormDistribution==true){
//
//            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//            std::default_random_engine generator (seed);
//            std::lognormal_distribution<double> distribution(0.,1.0);
//
//            // Random Vector U
//            TPZFMatrix<TVar> Rand_U (matrixSize, 1, 0.);
//
//            for (int i = 0; i < matrixSize; i++) {
//                Rand_U(i,0) = distribution(generator);
//                distribution.reset();
//                fRand_U = Rand_U;
//            }
//        }
        return fRand_U;
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
        //REAL E = arc4random_uniform(3001) + 15300; //uniform distribution
		//f[0]   = E;
        
    }
    
    // Call this method in the PZMatElasticity 2D (Gaussian Field)
    virtual void Execute(const TPZVec<TVar> &f, int id) {
        f[0] = fU(id, 0); // gives to the TPZMaterial the correlated variable of the element
    }
    
    // Call this method in the PZMatElasticity for Element Failure Area
    virtual void ExecuteArea(const TPZVec<TVar> &f, int id) {
        TPZGeoEl *gelFail;
        gelFail = fgmesh->ElementVec()[id];
        if (gelFail->Type() == 3) {
            f[0] = gelFail->SideArea(8);
        }
        else {
        f[0] = 0.0;
        }
        
         //
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
    
    
    // Calcula Correlation Matrix
    TPZFMatrix<REAL> calcCorrelationMatrix() {
        
        std::cout << "\nCria matriz da norma entre os centroides (para a matriz de correlacao)" << std::endl;
        
        // Refinamento de elementos selecionados
        REAL e = M_E; // Numero de Euler
        REAL scale = 0.5; //40.0 * frw; // Valor de alpha, escala normalizada // variar: 1/4; 1.0; 4.0
        
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
                    
                    // if (fexponetial==true){
                    KCorr(i,j) = pow(e, -((r2*r2)/(scale*scale)));
                    //}
                    
                    //else if (fspherical==true){
                    //insert function
                    //}
                }
                
                else {
                    // Verifica se el atual eh quadrilatero
                    std::cout<< "Element Type Error" << std::endl;
                }
            }
        }
        return KCorr;
    }
    
    
    // Calcula Correlation Matrix para Poços Inclinados
    TPZFMatrix<REAL> calcCorrelationMatrixInclined() {
        
        std::cout << "\nCria matriz dos centroides dos elementos " << std::endl;
        
        // Refinamento de elementos selecionados
        REAL e = M_E; // Numero de Euler
        REAL scale = 0.5; //40.0 * frw; // Valor de alpha, escala normalizada // variar: 1/4; 1.0; 4.0
        // scale simulado: 40* frext, errado corrigir
        
        TPZFMatrix<REAL> CenterNorm(fmatsize, fmatsize, 0.0);
        
        // Matriz de correlacao
        TPZFMatrix<REAL> KCorr(fmatsize, fmatsize, 0.0);
        
        REAL Pi = M_PI;
        
        //******* angulos COLOCADOS A MAO para fazer teste *********
        REAL alpha = 0.; // azimuth
        REAL beta = 0.; // inclination
        alpha = (fdirection*(Pi/180)); // azimuth
        beta = (finclination*(Pi/180)); // inclination
        
        //  Geeting all coordinates
        TPZGeoEl *gel;
        TPZFMatrix<REAL> Coordinates(fmatsize, 4, 0.0); //nanana
        TPZFMatrix<REAL> rotCoordinates(fnSquareElements, 4, 0.0);
        TPZManVector<REAL> centerpsi(3), center(3);
        TPZManVector<REAL, 3> CenterPoint;
        
        for (int i = 0; i < fnSquareElements; i++) {
            gel = fgmesh->ElementVec()[i];
            gel->CenterPoint(8, centerpsi);
            gel->X(centerpsi, center);
            
            CenterPoint = center;
            
            //	/*3*/	EQuadrilateral
            if (gel->Type() == 3) {
                //Coordinates
                REAL xx = CenterPoint[0];
                REAL yy = CenterPoint[1];
                REAL zz = CenterPoint[2];
                
                Coordinates(i, 0) = i;
                Coordinates(i, 1) = xx;
                Coordinates(i, 2) = yy;
                Coordinates(i, 3) = zz;
            }
        }
        
        int z = 0; // z <= (rext/h);
        int signal = 1;
        REAL altura = (z+(z-1))*(fh/2);
        for (int k = fnSquareElements; k < fmatsize; k += fnSquareElements) {
            if (k >= fmatsize/2 && signal > 0) {
                z = 0;
                signal = -1;
            }
            
            if (k % fnSquareElements == 0) {
                z++;
                altura = signal * (z + (z-1)) * (fh/2);
            }
            
            //std::cout << k << std::endl;
            for (int j = 0; j < fnSquareElements; j++) {
                Coordinates(k+j, 0) = k+j;
                Coordinates(k+j, 1) = Coordinates(j, 1);
                Coordinates(k+j, 2) = Coordinates(j, 2);
                Coordinates(k+j, 3) = altura;
            }
        }
        
        
        // Rotate fnSquareElements Coordinates and alocate in rotCoordinates
        for (int i = 0; i < fnSquareElements; i++) {
            rotCoordinates(i, 0) = i;
            rotCoordinates(i, 1) = Coordinates(i,1)*cos(alpha)*cos(beta) + Coordinates(i,2)*
            cos(beta)*sin(alpha) - Coordinates(i,3)*sin(beta);
            rotCoordinates(i, 2) = Coordinates(i,2)*cos(alpha) - Coordinates(i,1)*sin(alpha);
            rotCoordinates(i, 3) = Coordinates(i,3)*cos(beta) + Coordinates(i,1)*cos(alpha)*
            sin(beta) + Coordinates(i,2)*sin(alpha)*sin(beta);
        }
        
        // Getting rotCoordinates in Coordinates
        for (int i = 0; i < fnSquareElements; i++) {
            Coordinates(i, 0) = i;
            Coordinates(i, 1) = rotCoordinates(i, 1);
            Coordinates(i, 2) = rotCoordinates(i, 2);
            Coordinates(i, 3) = rotCoordinates(i, 3);
        }
        
        
        //std::cout << Coordinates << std::endl;
        std::ofstream out_Coordinates("Coordinates.txt");
        Coordinates.Print("XYZ = ",out_Coordinates,EMathematicaInput);
        
        //std::cout << rotCoordinates << std::endl;
        std::ofstream out_rotCoordinates("rotCoordinates.txt");
        rotCoordinates.Print("XYZ = ",out_rotCoordinates,EMathematicaInput);
        
        std::cout << "\nCria matriz da norma entre os centroides e Matriz de Correlacao" << std::endl;
        
        // Matriz da distancia entre os centroides
        for (int i = 0; i < fmatsize; i++) {
            for (int j = 0; j < fmatsize; j++) {
                
                REAL dx = pow((Coordinates(i,1)-Coordinates(j,1)), 2);
                REAL dy = pow((Coordinates(i,2)-Coordinates(j,2)), 2);
                REAL dz = pow((Coordinates(i,3)-Coordinates(j,3)), 2);
                
                CenterNorm(i,j) = sqrt(dx + dy + dz);
                
                REAL r = CenterNorm(i,j);
                REAL r2 = pow(r, 2);
                
                // if (fexponetial==true){
                KCorr(i,j) = pow(e, -((r2*r2)/(scale*scale)));
                //}
                
                //else if (fspherical==true){
                //insert function
                //}
            }
        }
        
        //        std::cout << "Numero colunas: " << KCorr.Cols() << std::endl;
        //        std::cout << "Numero linhas: " << KCorr.Rows() << std::endl;
        //        std::cout << "Penultimo valor " << KCorr(fmatsize-1,fmatsize-2) << std::endl;
        //        std::cout << "Ultimo valor " << KCorr(fmatsize-1,fmatsize-1) << std::endl;
        
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
