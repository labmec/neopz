//
//  hdiv3dpaper201504.h
//  PZ
//
//  Created by Douglas Castro on 22/04/15.
//
//

#ifndef __PZ__hdiv3dpaper201504__
#define __PZ__hdiv3dpaper201504__

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



/**
 * Classe que define o problema do laplaceano no cubo
 */

class Hdiv3dPaper201504{
private:
    int fDim;
    
    int fmatId;
    
    int fdirichlet;
    int fneumann ;
    
    int fbc0;
    int fbc1;
    int fbc2;
    int fbc3;
    int fbc4;
    int fbc5 ;
    int fmatskeleton;
    
    bool fisH1;
    
    bool ftetra;
    
    bool fprisma;
    
    bool isgeoblend;
    
    TPZFMatrix< int > tetraedra_2;
    //int tetraedra_2[6][4];
    
public:
    
    enum ApproximationSpace { EH1, EHDiv, EHDivStar, EHDivStarStar };
    
    enum Eltype { ECub, EPrism, ETetra };
    
    Hdiv3dPaper201504();
    
    ~Hdiv3dPaper201504();
    
    void Run(ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZFMatrix< REAL > &errors);
    
    void PrintErrors(ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZVec<REAL> &errors, std::ostream &output);
    
private:
    
    TPZGeoMesh *GMeshWithPrism( int ndiv);
    
    TPZGeoMesh *CreateOneCuboWithTetraedrons(int64_t nelem);
    
    void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem);
    
    TPZGeoMesh *CreateOneCubo(int nref);
    TPZGeoMesh *CreateOneQuadraticCube(int nref);
    
    /* Malhas computacionais */
    TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
    
    //solucao exata
    static void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
    static void SolExataH1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
    
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
    
    void ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, int pos, TPZFMatrix< REAL > &errors );
    void ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, std::ostream &out, int DoFT, int DofCond);
    
    void ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, int pos, TPZFMatrix< REAL > &errors);
    
    void ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, std::ostream &out, int DoFT, int DofCond);
    
    void ChangeExternalOrderConnects(TPZCompMesh *mesh);
    
    void SolveSyst(TPZLinearAnalysis &an, TPZCompMesh *fCmesh);
    
    bool MyDoubleComparer(REAL a, REAL b)
    {
	if (IsZero(a-b)){
	    return true;
	}
	else{
	    return false;
	}
    }
    
    void setTetraTrue()
    {
        ftetra = true;
    }
    void setPrismaTrue()
    {
        fprisma = true;
    }
    void setH1True()
    {
        fisH1 = true;
    }
    bool getIsH1(bool &EH1){
        EH1 = fisH1;
        return EH1;
    }
    
    
};

#endif /* defined(__PZ__hdiv3dpaper201504__) */
