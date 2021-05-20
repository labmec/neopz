//
//  hdivCurvedJCompAppMath.cpp
//  PZ
//
//  Created by Douglas Castro on 18/06/15.
//
//

#ifndef __PZ__hdivCurvedJCompAppMath__
#define __PZ__hdivCurvedJCompAppMath__

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
#include "tpzquadratictrig.h"

#include "tpzchangeel.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZBSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParFrontStructMatrix.h"

#include "TPZLinearAnalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "TPZReadGIDGrid.h"
#include "TPZLinearAnalysis.h"

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

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;
using namespace pzshape;

class hdivCurvedJCompAppMath{

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
    
    static bool probAtCircle;
    static bool probAtCylinder;
    static bool probAtSphere;

    int fmatskeleton;
    
    bool ftriang;
    
    bool isgeoblend;

public:
    
    enum ApproximationSpace { EH1, EHDiv, EHDivStar, EHDivStarStar };
    
    enum Eltype { EQuad, ETriang };
    
    enum geomDomain {ECircle, ECylinder, ESphere};
    
    hdivCurvedJCompAppMath();
    hdivCurvedJCompAppMath(geomDomain geodomain);
    
    ~hdivCurvedJCompAppMath();
    
    void Run(geomDomain geodomain, ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZFMatrix< REAL > &errors);

    void PrintErrors(geomDomain geodomain, ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZVec<REAL> &errors, std::ostream &output);

private:
    
    /*  Geometrical meshes */
    TPZGeoMesh *MakeCircle( int ndiv);
    TPZManVector<REAL,3> ParametricCircle(REAL radius,REAL theta);
    TPZManVector<REAL,3> ParametricSphere(REAL radius,REAL phi,REAL theta);
    TPZGeoMesh *MakeSphereFromQuadrilateral(int dimensao, bool triang, int ndiv);
    void RotateNode(TPZVec<REAL> &iCoords, REAL CounterClockwiseAngle, int &Axis);
    TPZGeoMesh *GMeshCilindricalMesh( int ndiv);
    void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);
    
    /* Computational meshes */
    TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
    void SetupDisconnectedHdivboud(const int left,const int rigth, TPZCompMesh * cmesh);
    
    //Exact Solution
    static void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
    static void SolExataH1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
    
    //Force function
    static void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
    static void ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux);
    
    //Dirichlet B. Conditions
    static void ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    static void ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp);
    
    
    void ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, int pos, TPZFMatrix< REAL > &errors );
    void ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, std::ostream &out, int DoFT, int DofCond);
    
    void ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, int pos, TPZFMatrix< REAL > &errors);
    
    void ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, std::ostream &out, int DoFT, int DofCond);
    
    void ChangeExternalOrderConnects(TPZCompMesh *mesh);
    
    void SolveSyst(TPZLinearAnalysis &an, TPZCompMesh *fCmesh);
    
    void setTriangTrue()
    {
        ftriang = true;
    }

    void setDimension(geomDomain geodomain)
    {
      switch(geodomain){
	case ECircle:
	{
	  fDim = 2;
	}
	  break;
	case ECylinder:
	{
	  fDim = 2;
	}
	break;
	case ESphere:
	{
	  fDim = 3;
	}
	  break;
	default:
            break;
      }
	

    }
    

};


#endif /* defined(__PZ__hdivCurvedJCompAppMath__) */
