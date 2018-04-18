#ifndef CEDRICTEST_HH
#define CEDRICTEST_HH

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"

class TPZGeoMesh;

/** To numerical test from Cedric+Renato publication */
class TCedricTest {

public:
    static TPZManVector<REAL,3> fX0, fEps;
    
    TPZGeoMesh fDeformed;
public:
    TCedricTest();

    void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem);

    /** Constructing geometrical mesh depends on type of element wished. */
    TPZGeoMesh *HexahedralMesh(int64_t nelem,int MaterialId);
    TPZGeoMesh *PyramidalAndTetrahedralMesh(int64_t nelem,int MaterialId);
    TPZGeoMesh *TetrahedralMesh(int64_t nelem,int MaterialId);
    TPZGeoMesh *TetrahedralMeshUsingRefinement(int64_t nelem,int MaterialId);

    int AddBoundaryElements(TPZGeoMesh *gmesh);

    TPZCompMesh *GenerateCompMesh(TPZGeoMesh *gmesh);
    
    void CreateCondensedElements(TPZCompMesh *cmesh);
    
    void UnwrapElements(TPZCompMesh *cmesh);

    static REAL fx(REAL x, REAL x0, REAL eps) {
        REAL a = x*x*(1-x)*(1-x)*exp(-(x-x0)*(x-x0)/eps);
//        REAL a = 4.*x*(1-x);
        return a;
    }

    static REAL dfx(REAL x, REAL x0, REAL eps) {
        REAL a = 2*x*(1-x)*(1-x)*exp(-(x-x0)*(x-x0)/eps);
        REAL b = -2*(1-x)*x*x*exp(-(x-x0)*(x-x0)/eps);
        REAL c = -2.*(x-x0)*x*x*(1-x)*(1-x)*exp(-(x-x0)*(x-x0)/eps)/eps;
        REAL result = a+b+c; 
//        REAL result = 4.*(1-2.*x);
        return result;
    }

    static REAL d2fx(REAL x, REAL x0, REAL eps) {

        REAL result = 2*pow(M_E,((x - x0)*(-x + x0))/eps)*(1 - x)*(1-x) -
        8*pow(M_E,((x - x0)*(-x + x0))/eps)*(1 - x)*x +
        2*pow(M_E,((x - x0)*(-x + x0))/eps)*(x*x) -
        (2*pow(M_E,((x - x0)*(-x + x0))/eps)*(1 - x)*(1-x)*x*x)/eps +
        4*pow(M_E,((x - x0)*(-x + x0))/eps)*(1 - x)*(1-x)*x*
        (-((x - x0)/eps) + (-x + x0)/eps) -
        4*pow(M_E,((x - x0)*(-x + x0))/eps)*(1 - x)*(x*x)*
        (-((x - x0)/eps) + (-x + x0)/eps) +
        pow(M_E,((x - x0)*(-x + x0))/eps)*(1 - x)*(1-x)*(x*x)*
        (-((x - x0)/eps) + (-x + x0)/eps)*(-((x - x0)/eps) + (-x + x0)/eps);
//        REAL result = -8.;
        return result;
    }
    
    static void Exact(const TPZVec<REAL> &x, TPZVec<STATE> &func, TPZFMatrix<STATE> &deriv) {
        REAL v[3] = {fx(x[0],fX0[0],fEps[0]), fx(x[1],fX0[1],fEps[1]), fx(x[2],fX0[2],fEps[2])};
        func[0] = v[0]*v[1]*v[2];
        for(int i=0;i<3;i++) {
            REAL dvz = dfx(x[i],fX0[i],fEps[i]);
            deriv(i,0) = dvz*v[(i+1)%3]*v[(i+2)%3];
        }
    }
    
    /// Deform the geometric mesh according to the coordinates of fDeformed
    void DeformGMesh(TPZGeoMesh &gmesh);
    
    /// verify if the faces without neighbour should be orthogonal to the main planes
    void CheckConsistency(TPZGeoMesh *mesh);
    
    void Run(int nsubdivisions,int geocase,int POrder,int MaterialId,std::ostream &out=std::cout);
    
    void LoadInterpolation(TPZCompMesh *cmesh);
    
    void InterpolationError(int nsubdivisions,int geocase, int MaterialId,std::ostream &out);
};

#endif
