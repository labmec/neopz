#ifndef CEDRICTEST_HH
#define CEDRICTEST_HH

#include "pzmanvector.h"

class TPZGeoMesh;

/** To numerical test from Cedric+Renato publication */
class TCedricTest {

public:
    static TPZManVector<REAL,3> fX0, fEps;
    
    TPZGeoMesh fDeformed;
public:
    TCedricTest();

    void GenerateNodes(TPZGeoMesh *gmesh, long nelem);

    /** Constructing geometrical mesh depends on type of element wished. */
    TPZGeoMesh *HexahedralMesh(long nelem,int MaterialId);
    TPZGeoMesh *PyramidalAndTetrahedralMesh(long nelem,int MaterialId);
    TPZGeoMesh *TetrahedralMesh(long nelem,int MaterialId);
    TPZGeoMesh *TetrahedralMeshUsingRefinement(long nelem,int MaterialId);

    int AddBoundaryElements(TPZGeoMesh *gmesh);

    TPZCompMesh *GenerateCompMesh(TPZGeoMesh *gmesh);

    static REAL fx(REAL x, REAL x0, REAL eps) {
        REAL result = 0.;
        REAL a = exp(-(x-x0)*(x-x0)/eps);
        REAL b = exp(-x0*x0/eps)*(1.-x);
        REAL c = exp(-(1.-x0)*(1.-x0)/eps)*x;
        result = a-b-c;
        return result;
    }

    static REAL dfx(REAL x, REAL x0, REAL eps) {
        REAL a = -exp(-(1-x0)*(1.-x0)/eps);
        REAL b = exp(-x0*x0/eps);
        REAL c = -2.*(x-x0)*exp(-(x-x0)*(x-x0)/eps)/eps;
        REAL result = a+b+c; 
        return result;
    }

    static REAL d2fx(REAL x, REAL x0, REAL eps) {
        REAL a = 2.*exp(-(x-x0)*(x-x0)/eps)/eps;
        REAL b = -4.*(x-x0)*(x-x0)*exp(-(x-x0)*(x-x0)/eps)/eps/eps;
        REAL result = a+b;
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
};

#endif
