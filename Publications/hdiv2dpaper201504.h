//
//  Hdiv2dPaper201504.cpp
//  PZ
//
//  Created by Douglas Castro on 18/03/15.
//
//

#ifndef __PZ__hdiv2dpaper201504__
#define __PZ__hdiv2dpaper201504__

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

class Hdiv2dPaper201504{

private:
    int fDim;
    
    int fmatId;
    
    int fdirichlet;
    int fneumann;
    
    int fbc1;
    int fbc2;
    int fbc3;
    int fbc4;

    int fmatskeleton;
    
    bool ftriang;

public:
    
    enum ApproximationSpace { EH1, EHDiv, EHDivStar, EHDivStarStar };
    
    enum Eltype { EQuad, ETriangle };
    
    Hdiv2dPaper201504();
    
    ~Hdiv2dPaper201504();
    
    void Run(ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZFMatrix< REAL > &errors);

    void PrintErrors(ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZVec<REAL> &errors, std::ostream &output);

private:
    
    /*  Geometrical meshes */
    TPZGeoMesh *GMesh(int dim, bool ftriang, int ndiv);
    
    /* Computational meshes */
    TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
    
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

    
};


#endif /* defined(__PZ__LaplaceInQuadrilateral__) */
