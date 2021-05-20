/**
 * \file
 * @brief Contains the TPZMulticamadaOrthotropic class.
 */

#ifndef MULTICAMADAORTH
#define MULTICAMADAORTH

#include <iostream>
#include "pzvec.h"
#include "pzstack.h"
#include "TPZPlacaOrthotropic.h"
#include "tpzautopointer.h"

class TPZGeoMesh;
class TPZCompMesh;
class TPZPlacaOrthotropic;
class TPZMatOrthotropic;
class TPZLinearAnalysis;
class TPZMaterial;

/**
 * @ingroup material
 * @brief Gerencia um conjunto de placas
 * dispostas em forma multicamada
 */

class TPZMulticamadaOrthotropic {
	
	/** @brief Geometric mesh with shells */
	TPZGeoMesh             *fGeoMesh;
	/** @brief Computational mesh to calculations */
	TPZCompMesh            *fCompMesh;
	/** @brief Shells vector */
	TPZStack<TPZPlacaOrthotropic> fPlacaOrth;

	/** @brief Dimension of the shells (must to be constant for all shells) */
	REAL fDx,fDy;
	/** @brief Number of elementos at x and y axes : fNelx, fNely */
	int64_t fNelx, fNely;
	REAL fZMin, fZMax;

	REAL fMX[3],fMY[3],fMXY[3],fQX[3],fQY[3],fNX[3],fNY[3],fNXY[3];
	REAL fdMXdX[3],fdMXdY[3],fdMYdX[3],fdMYdY[3],fdMXYdX[3],fdMXYdY[3],
	fdQXdX[3],fdQXdY[3],fdQYdX[3],fdQYdY[3],fdNXdX[3],fdNXdY[3],fdNYdX[3],fdNYdY[3],fdNXYdX[3],fdNXYdY[3];
	int fLinearX,fLinearY;
	TPZManVector<REAL,3> fDirx, fDiry;
	
	/**
	 * @brief Relaxation factor to correct resulting forces.
	 * @since Feb 10, 2004
	 */
	REAL fCorrect;
	
public:
	/** @brief Construtor */
	TPZMulticamadaOrthotropic(REAL z,REAL dx,REAL dy, int64_t nelx, int64_t nely, REAL Correct = 1.0);
	/** @brief Destrutor */
	~TPZMulticamadaOrthotropic(){}
	
	/** @brief Adds shells */
	void AddPlacaOrtho(TPZMaterial * material, REAL height);
	/** @brief Creates a computational mesh to all the shells */
	void GenerateMesh();

	void Print(std::ostream &out = std::cout);
	/*criando método para retornar a altura da multicamada*/
	REAL Height();
	/*criando método para retornar fPlacaOrth*/
	TPZVec<TPZPlacaOrthotropic> &RPlacaOrtho(){return fPlacaOrth;}
	/*criando método para contar quant de placas*/
	int NPlacas();
	/**
	 * @brief Compute a tension state corresponding to the difference between the target state \n
	 * and tension state loaded in the solution
	 */
	void AnalyticTensor(TPZVec<REAL> &co, TPZFMatrix<REAL> &tensor);
	
	/** @brief Tensor which needs to be applied at the given coordinate */
	void Tensor(TPZVec<REAL> &x, int placa, TPZFMatrix<REAL> &tensor);
	/** @brief Computes the global efforts of the finite element solution */
	void ComputeCenterForces();
	
	void ComputeSolution(std::ostream &out = std::cout,int print = 0);
	
	void ComputeSolution(TPZMaterial *mat,std::ofstream &out,int64_t numiter);

	/**
	 * @name Set data methods
	 * @{
	 */
	 
	void SetMX(REAL MX) { 
		fMX[0] = MX;
		fMX[2] = MX;
	}
	
	void SetNX(REAL NX) {
		fNX[0] = NX;
		fNX[2] = NX;
	}
	
	void SetNY(REAL NY) {
		fNY[0] = NY;
		fNY[2] = NY;
	}
	
	void SetMY(REAL MX) { 
		fMY[0] = MX;
		fMY[2] = MX;
	}
	
	void SetNXY(REAL NX) {
		fNXY[0] = NX;
		fNXY[2] = NX;
	}
	
	void SetMXY(REAL MXY) {
		fMXY[0] = MXY;
		fMXY[2] = MXY;
	}
	
	void SetQX(REAL QX) {
		fQX[0] = QX;
		fQX[2] = QX;
		fdMXdX[0] = QX;
		fdMXdX[2] = QX;
		if(QX != 0.) fLinearX = 1;
		else fLinearX = 0;
	}
	
	void SetQY(REAL QY) {
		fQY[0] = QY;
		fQY[2] = QY;
		fdMYdY[0] = QY;
		fdMYdY[2] = QY;
		if(QY != 0.) fLinearY = 1;
		else fLinearY = 0;
		
	}

	/** @} */
	
	TPZGeoMesh *GeoMesh(){return fGeoMesh;}
	
	TPZCompMesh *CompMesh(){return fCompMesh;}
	
	void SetCorrect(REAL Correct){ fCorrect = Correct;}
	
	REAL CorrectFactor(){return fCorrect;}
	
	void PrintTensors(std::ostream &out);
	
	void PrintTensors(std::ostream &out,TPZFMatrix<REAL> &tensorin,TPZFMatrix<REAL> &tensorout);
	
	void PrintCenterForces(std::ostream &out);
};

#endif
