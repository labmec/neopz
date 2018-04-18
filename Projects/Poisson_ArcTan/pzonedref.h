#ifndef TPZONEDREF
#define TPZONEDREF

#include <fstream>
#include <cmath>

#include "pzfmatrix.h"

/**
 * Class TPZOneDRef implements one dimensional
 * refinement pattern
 * @ingroup analysis
 */
class TPZOneDRef {
    /**
     * External log file
     */
    static std::ofstream fLogFile;
    
    /**
     * \f[fMS_{1}S_{1}=\int _{\widetilde{\Omega _{1}}}\psi _{s_{i}.}.\psi _{s_{j}.}.d\widetilde{\Omega }_{1}\f]
     */
    TPZFMatrix<REAL>	fMS1S1;
    
    /**
     * \f[fMS_{1}B=\int _{\widetilde{\Omega _{1}}}\psi _{s_{i}.}.\psi _{b_{j}.}.d\widetilde{\Omega }_{1}\f]
     */
    TPZFMatrix<REAL>	fMS1B;
    
    /**
     * \f[fMS_{2}B=\int _{\widetilde{\Omega _{1}}}\psi _{s_{i}.}.\psi _{b_{j}.}.d\widetilde{\Omega }_{2}\f]
     */
    TPZFMatrix<REAL>	fMS2B;
    
    /**
     * \f[fMBB=\int _{\widetilde{\Omega }_{1}+\widetilde{\Omega }_{2}}\psi _{s_{i}.}.\psi _{s_{j}.}.d\widetilde{\Omega }\f]
     */
    TPZFMatrix<REAL>	fMBB;
    
    /**
     * \f[fKS_{1}S_{1}=\int _{\widetilde{\Omega _{1}}}\psi '_{s_{i}.}.\psi '_{s_{j}.}.d\widetilde{\Omega }_{1}\f]
     */
    
    TPZFMatrix<REAL>	fKS1S1;
    
    /**
     * \f[fKS_{1}B=\int _{\widetilde{\Omega _{1}}}\psi '_{s_{i}.}.\psi '_{b_{j}.}.d\widetilde{\Omega }_{1}\f]
     */
    TPZFMatrix<REAL>	fKS1B;
    
    /**
     *\f[fKS_{2}B=\int _{\widetilde{\Omega _{2}}}\psi '_{s_{i}.}.\psi '_{b_{j}.}.d\widetilde{\Omega }_{2}\f]
     */
    TPZFMatrix<REAL>	fKS2B;
    
    /**
     * \f[fKBB=\int _{\widetilde{\Omega }_{1}+\widetilde{\Omega }_{2}}\psi '_{b_{i}.}.\psi '_{b_{j}.}.d\widetilde{\Omega }\f]
     */
    TPZFMatrix<REAL>	fKBB;
    
    /**
     * fU : externally provided solution, adapted for the local id sequence and with homogeneous boundary values
     * fRhs* : inner product of the externally provided solution with the shape functions
     */
    TPZFMatrix<REAL>fU,fRhs1,fRhs2,fRhsb,fStiffU;
    
    /**
     * fM : contains the scaled sum of fMS1S1 and fKS1S1
     */
    TPZFMatrix<REAL>fM;
    /**
     * Size of the small elements
     */
    REAL fDelx;
    /**
     * polynomial order of the externally provided solution
     */
    int fp1,fp2,fpb;
    
    /**
     * fNState : number of state variables
     */
    int fNState;
    
public:
    
    /**
     * Constructor
     * Creates an object from a number of state variables
     */
    TPZOneDRef(int nstate);
    
    /**
     * Computes the refinement pattern which best will approximate the given solution
     * @param U computed one dimensional solution
     * @param id ids which determine the orientation of the shape functions
     * @param p1 [in] interpolation order of the first element, on output, best approximation to U
     * @param p2 [in] interpolation order of the second element, on output, best approximation to U
     * @param hp1
     * @param hp2
     * @param hperror
     * @param delx
     */
    REAL BestPattern(TPZFMatrix<REAL> &U, TPZVec<int64_t> &id, int &p1, int &p2, int &hp1, int &hp2, REAL &hperror, REAL delx);
    
    /**
     * Print oject data
     */
    void Print(char *msg = 0, std::ostream &out = std::cout);
    
    static int gMaxP;
    
    /**
     * Validation routines
     */
    static int main();
private:
    
    /**
     * Verify the compatibility odd functions
     * the function must start at smaller id ant end at greater id
     */
    void TransformU(TPZFMatrix<REAL>&U, TPZVec<int64_t> &id, int p1, int p2);
    
    /**
     * Integrate Matrices
     */
    void IntegrateMatrices();
    
    /**
     * Load and transform small elements solution
     * to coarse element interpolation order
     * @param U - coarse solution matrix
     * @param p1 - small element one refinement order
     * @param p2 - small element two refinement order
     * @param delx
     */
    void LoadU(TPZFMatrix<REAL>&U, int p1, int p2, REAL delx);
    
    /**
     * Store element vectors compabilization
     */
    void BuildIndexVectors(int p1, int p2, int pb);
    
    /**
     * Returns the minimum obtained error using
     * p1 and p2 orders
     */
    REAL Error(int p1, int p2);
    
    /**
     * Return the minimum error obtained
     * using non ?h? refinement
     */
    REAL Error(int pb);
    
    /**
     * Generates de stiffness matrix using p1 and p2 order
     */
    void BuildStiffness(int p1, int p2, TPZFMatrix<REAL>&stiff);
    
    /**
     * Flag to activate the p refinement. By default is set to true
     * If the maximum p order is verified then this flag is set to false 
     */
    int fTryP;
    
};


#endif
