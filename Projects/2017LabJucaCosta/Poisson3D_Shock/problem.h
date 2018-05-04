#ifndef MYPROBLEMCLASS
#define MYPROBLEMCLASS

#include <fstream>
#include <cmath>

#include "pzfmatrix.h"

/**
 * Class ProblemClass implements a special problem based on Poisson equation, given the rigth function and the exact solution if it exists.
 * @ingroup tutorial
 */
class ProblemClass {
    /**
     * External log file
     */
    static std::ofstream fLogFile;
    
    
public:
    
    /**
     * Constructor
     * Creates an object from a number of state variables
     */
    ProblemClass(int nstate) {
	}
    
    /**
     * Print oject data
     */
    void Print(char *msg = 0, std::ostream &out = std::cout) {
	}
        
    /**
     * Validation routines
     */
    static int main() {
		return 1;
	}
    
};

/** PROBLEM WITH HIGH GRADIENT ON CIRCUNFERENCE  ---  DATA */
extern STATE ValueK;
extern REAL GlobalMaxError;
extern STATE F;
// Circunference with high gradient - data
extern TPZManVector<REAL,3> CCircle;
extern REAL RCircle;
extern int ModelDimension;

REAL VarTimesOneMinusVar(int var,int dim,const TPZVec<REAL> &x);
REAL VarTimesVarMinusOne(int var,int dim,const TPZVec<REAL> &x);

// Computes ||x-x0||^2
REAL SquareNorm(int dim,const TPZManVector<REAL> &x,TPZVec<REAL> &x0);

REAL ArgumentArc(int dim,REAL Radius,const TPZVec<REAL> &x,TPZVec<REAL> &x0);

REAL PolinomicValue(int var,int dim,const TPZVec<REAL> &x,TPZVec<REAL> &x0);

STATE GradientNorm(TPZInterpolatedElement *el);
STATE Laplacian(TPZInterpolatedElement *el);
bool GradientAndLaplacian(TPZInterpolatedElement *el,REAL &Grad,REAL &Laplacian);
void ComputingMaxGradientAndLaplacian(TPZCompMesh *cmesh,REAL &MaxGrad,REAL &MaxLaplacian);


bool LaplacianValue(TPZCompEl *el,REAL &Laplacian);
void ComputingMaxLaplacian(TPZCompMesh *cmesh,REAL &MaxLaplacian,REAL &Min);


/* To high gradient generate by arc tangent term - and considering this term joined with Plus Pi/2 */
void ExactSolutionArcTangent2(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol, TPZVec<STATE> &ddsol);

void ExactSolutionArcTangent(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
bool GradientAndLaplacianOnCorners(TPZInterpolatedElement *el,REAL &Grad,REAL &Laplacian);
/** We are considering - f, because is as TPZMatPoisson3d was implemented in Contribute method */
void RightTermArcTangentBad(const TPZVec<REAL> &x, TPZVec<STATE> &force, TPZFMatrix<STATE> &dforce);
void RightTermArcTangent(const TPZVec<REAL> &x, TPZVec<STATE> &force);

/* To high gradient generate by arc tangent term - and a term Pi was factored as coefficient of the funcion */
bool GradientAndLaplacianToSphereOnCorners(TPZInterpolatedElement *el,REAL &Grad,REAL &Laplacian);
STATE LaplacianToSphere(TPZInterpolatedElement *el);
STATE GradientNormToSphere(TPZInterpolatedElement *el);
void ExactSolutionSphere(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol, TPZVec<STATE> &ddsol);
void ExactSolutionSphere(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
// We are considering - f, because is as TPZMatPoisson3d was implemented in Contribute method *
void RightTermSphere(const TPZVec<REAL> &x, TPZVec<STATE> &force, TPZFMatrix<STATE> &dforce);


/** Laplace problem on L-Shape domain */
void ExactSolLaplace(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
void ExactSolLaplaceBC(const TPZVec<REAL> &x, TPZVec<STATE> &sol);


/** UTILITIES */
bool CheckFunctionAndDerivatives(std::ifstream &inpoint,std::ifstream &math);

#endif
