//
//  LaplaceInCircle.h
//  PZ
//
//  Created by Douglas Castro on 1/29/15.
//
//

#ifndef __PZ__LaplaceInCircle__
#define __PZ__LaplaceInCircle__

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
#include "tpzquadratictrig.h"
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
#include "pzmatmixedpoisson3d.h"
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



#include "pyramidalmesh.h"

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;
using namespace pzshape;



/** 
 * Classe que define o problema do laplaceano no circulo
 */

class LaplaceInCircle{
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
    
    bool ftriang;
    
    bool isgeoblend;
    
    
public:
    
    
    LaplaceInCircle();
    
    ~LaplaceInCircle();
    
    void Run(int ordemP, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv, std::ofstream &saidaErro, bool HdivMaisMais);
    
    /**
     * Cria uma malha geometrica formada por elementos do tipo Geoblend 
     */
    TPZGeoMesh *GMeshCirculoGeob( int ndiv);
    TPZGeoMesh *GMeshCirculoGeobQuart( int ndiv);
    TPZGeoMesh *GMeshCirculoGeobQuartT( int ndiv);
    TPZGeoMesh *GMeshCirculoTriangGeob( int ndiv);
    TPZGeoMesh *GMeshQuartoCirculoGeobAT( int ndiv);
    TPZGeoMesh *GMeshCirculoQuadraticQuartT( int ndiv);
    TPZGeoMesh *MakeQuarterOfCircle( int ndiv);
    TPZGeoMesh *MakeCircle( int ndiv);    
    /**
     * Cria uma malha geometrica formada por elementos do tipo Geoblend
     */
    TPZGeoMesh *GMeshCirculoQuad( int ndiv);
    
    TPZManVector<REAL,3> ParametricCircle(REAL radius,REAL theta);
    
    
    TPZGeoMesh *GmeshCirculoPorElementosRetos( int ndiv);
    TPZVec<REAL> PolarToKartesian(REAL r, REAL theta, TPZManVector<REAL> xc);
    
    
    
    /* Malhas computacionais */
    TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim);
    TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
    TPZCompMesh *CMeshMixedWrap(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
    
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

    static void ErrorL2(TPZCompMesh *l2mesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv);
    
    static void ErrorHDiv(TPZCompMesh *hdivmesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv);
    
    static void ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, std::ostream &out, int DoFT, int DofCond);
    
    void ChangeExternalOrderConnects(TPZCompMesh *mesh);
    
    void setTriangTrue()
    {
        ftriang = true;
    }
    void setH1True()
    {
        fisH1 = true;
    }
    bool getIsH1(bool &EH1){
        EH1 = fisH1;
		return fisH1;
    }
    
};

#endif /* defined(__PZ__LaplaceInCircle__) */
