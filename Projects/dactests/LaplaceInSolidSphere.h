//
//  LaplaceInSolidSphere.hpp
//  PZ
//
//  Created by Omar on 5/9/16.
//
//

#ifndef LaplaceInSolidSphere_hpp
#define LaplaceInSolidSphere_hpp

#include <stdio.h>



#include <iostream>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"

#include "tpzchangeel.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPBSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "TPZReadGIDGrid.h"
#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

//#include "pzhdivfull.h"
#include "pzelchdiv.h"

#include "pzgeopyramid.h"

#include "pznumeric.h"

#include "TPZExtendGridDimension.h"
#include "pzelchdivbound2.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"

#include "TPZLagrangeMultiplier.h"
#include "pzmatmixedpoisson3d.h"
#include "PZMatPoissonD3.h"

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

#include "pzcheckgeom.h"


#include "pyramidalmesh.h"

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;
using namespace pzshape;

//static int matId = 1;
//
//static int dirichlet = 0;
//static int neumann = 1;
//
//static int bc0 = -1;
//static int bc1 = -2;
//static int bc2 = -3;
//static int bc3 = -4;
//static int bc4 = -5;
//static int bc5 = -6;
//static int matskeleton = -7;

class LaplaceInSolidSphere {
private:
    int fDim;
    
    int fmatId;
    
    int fdirichlet;
    int fneumann;
    
    int fbc0;
    int fbc1;
    int fbc2;
    int fbc3;
    int fbc4;
    int fbc5;
    int fmatskeleton;
    
    bool fisH1;
    
    bool fIsNonLinearMeshQ;
    
    
public:
    
    LaplaceInSolidSphere( );
    
    ~LaplaceInSolidSphere();
    
    void Run(int ordemP, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv, std::ofstream &saidaErro, bool HdivMaisMais);
    
    int getDimension() const {return fDim;};
    
    TPZGeoMesh *GMeshSphericalShell( int ndiv);
    
    //---------------
    
    TPZGeoMesh *GMeshSphericalRingQuarter(int dimensao, bool triang, int ndiv);
    
    TPZGeoMesh *GMeshTropicodeCancer(int ndiv , TPZVec<bool>  &CurvesSides, bool isPlane, int plane);
    
    TPZGeoMesh *GMeshCirculoPolarArtico(int ndiv , TPZVec<bool>  &CurvesSides, bool isPlane, int plane);
    
    // theta (0,pi) angulo que se inicia no polo norte. phi (0,2pi) o angulo no plano xy
    TPZVec<REAL> SphereToKartesian(REAL r, REAL theta, REAL phi);
    TPZVec<REAL> SphereToKartesian(TPZManVector<REAL> xc, REAL r, REAL theta, REAL phi);
    TPZManVector<REAL,3> ParametricSphere(REAL radius, REAL phi,REAL theta);
    
    TPZGeoMesh *MakeSphereFromQuadrilateralFaces(int ndiv);
    
    TPZGeoMesh *MakeSphereFromLinearQuadrilateralFaces(int ndiv);
    
    //--------------------
    
    /* Malhas computacionais */
    TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
    
    void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);
    void SetupDisconnectedHdivboud(const int left,const int rigth, TPZCompMesh * cmesh);
    
    //solucao exata
    static void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
    
    //lado direito da equacao
    static void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
    static void ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux);
    
    //Para condicao de contorno de Dirichlet
    static void ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    
    //Para condicao de contorno de Neumann
    static void ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    
    static void ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  REAL &error_primal , REAL & error_dual);
    
    void ChangeExternalOrderConnects(TPZCompMesh *mesh);
    
    void SetNonLinearMesh(bool nonlinear){
     
        fIsNonLinearMeshQ = nonlinear;
        
    }
    
};


#endif /* LaplaceInSolidSphere_hpp */
