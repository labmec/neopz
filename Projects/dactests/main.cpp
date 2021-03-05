
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

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

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


#include "pyramidalmesh.h"

#include "tools.h"
#include "LaplaceInCylinder.h"
#include "LaplaceInCircle.h"
#include "LaplaceInSphere.h"
#include "LaplaceInQuadrilateral.h"
#include "LaplaceInCube.h"
#include "LaplaceInSolidSphere.h"

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;
using namespace pzshape;

//int matId = 1;
//
//int dirichlet = 0;
//int neumann = 1;
//
//int bc0 = -1;
//int bc1 = -2;
//int bc2 = -3;
//int bc3 = -4;
//int bc4 = -5;
//int bc5 = -6;
//int matskeleton = -7;

//// just for print data
///** @brief Map used norms */
std::map<REAL,REAL> fDebugMapL2, fDebugMapHdiv;
///** @brief Map used Degrees of Freedom */
std::map<int,int> fDebugDoF;

//int tetraedra_2[6][4]=
//{
//    {1,2,5,4},
//    {4,7,3,2},
//    {0,1,2,4},
//    {0,2,3,4},
//    {4,5,6,2},
//    {4,6,7,2}
//};
//
//int piramide_2[6][5]=
//{
//    {0,1,2,3,8},
//    {0,1,5,4,8},
//    {1,2,6,5,8},
//    {3,2,6,7,8},
//    {0,3,7,4,8},
//    {4,5,6,7,8}
//};


//bool MyDoubleComparer(REAL a, REAL b);
//
//void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem);
//void GenerateNodesPyramid(TPZGeoMesh *gmesh, int64_t nelem);
//
//TPZGeoMesh *GMesh(int dimensao, bool ftriang, int ndiv);
//
//TPZGeoMesh *GMeshCirculoGeob(int dimensao, int ndiv);
//TPZGeoMesh *GMeshCirculoQuad(int dimensao, int ndiv);
//
//TPZGeoMesh *CreateOneCubo(int nref=0);
//TPZGeoMesh *CreateOneCuboWithTetraedrons(int64_t nelem=1, int MaterialId=1);
//TPZGeoMesh *GMeshCubeWithPyramids(int64_t nelem=1, int MaterialId=1);
//
//
//TPZGeoMesh * BasicForm(int n, REAL t, REAL dt);
//void Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X);
//void Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X);
//void Parametricfunction3(const TPZVec<STATE> &par, TPZVec<STATE> &X);

int dim = 2;
REAL aa = 0.0;
REAL bb = 0.0;
REAL cc = 0.0;
REAL Epsilon = 0.4;


//REAL const Pi = M_PI;//4.*atan(1.);



//bool ftriang = false;//true;//
bool IsCubedomain = false;
bool IsPrism = false;
bool IsTetra = false;
bool IsPiram = false;
bool IsSphere = true;

//bool isspheredomain = true, iscircledomain = false, iscylinderdomain = false, isquaddomain = false;
//bool iscircledomain = true, isspheredomain = false, iscylinderdomain = false, isquaddomain = false;
// bool iscylinderdomain = true, iscircledomain = false, isspheredomain = false, isquaddomain = false;
bool isquaddomain = true, iscircledomain = false, isspheredomain = false, iscylinderdomain = false;


#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

#include "pztransfer.h"

int main(int argc, char *argv[])
{
//#ifdef LOG4CXX
//    InitializePZLOG();
//#endif
    
//  gRefDBase.InitializeAllUniformRefPatterns();
//	gRefDBase.InitializeRefPatterns();

    int p = 1;
    int ndiv = 0;
    ofstream saidaerros("ErroNormas.txt",ios::app);
    
    for(p=1;p<2;p++)
    {
        saidaerros << "\nPARA p = " << p << endl;
//        saidaerros << "ndiv " << setw(6) << "DoFT" << setw(20) << "DofCond" << setw(20) << "ErroL2Primal" << setw(20) << "ErroL2Dual" << setw(20) << "ErroL2Div" << setw(20) << "ErroHDivDual"  << endl;
        saidaerros << "ndiv" << setw(10) <<"NDoF"<< setw(20)<<"NDoFCond"<< setw(20)<< "Assemble"<< setw(20) << "Solve" << setw(20) <<"Ttotal" << setw(20) << "Error primal" << setw(20) <<"Error dual \n";
        
        for (ndiv=0; ndiv<1; ndiv++)
        {
            
            if (dim==2)
            {
                //TPZGeoMesh *gmesh2d = GMesh(2, ftriang, ndiv);
                if (iscircledomain) {
                    LaplaceInCircle  * circ = new LaplaceInCircle();
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    circ->Run( k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else if (iscylinderdomain)
                {
                    LaplaceInCylinder * cilind = new LaplaceInCylinder( );
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    cilind->Run( k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else  if(isspheredomain)
                {
                    LaplaceInSphere * sphere = new LaplaceInSphere( );
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    sphere->Run( k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else if(isquaddomain)
                {
                    LaplaceInQuadrilateral * quad = new LaplaceInQuadrilateral();
                    //quad->setTriangTrue();
                    //quad->setH1True();
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    quad->Run(k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else
                {
                    //Ops! Nao esta fazendo nada.
                    DebugStop();
                }
            }
            else // dim == 3
            {
                if (IsCubedomain)
                {
                    LaplaceInCube * cubo = new LaplaceInCube();
                    //cubo->setTetraTrue();
                    //cubo->setPrismaTrue();
//                    cubo->setH1True();
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    cubo->Run(k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);

                }
                else if(IsTetra)
                {
                    DebugStop();
//                    REAL dndiv = ndiv;
//                    int nref = (int) pow(2., dndiv);
//                    
//                    gmesh = CreateOneCuboWithTetraedrons(nref, matId);
                }
                else if (IsPrism)
                {
                    // Tem que escolher a malha certa.
                    DebugStop();
                }
                else if(IsPiram)
                {
                    DebugStop();
//                    REAL dndiv = ndiv;
//                    int nref = (int) pow(2., dndiv);
//                    gmesh = GMeshCubeWithPyramids( nref);
//                    PyramidalMesh(gmesh,ndiv);
                }
                else if(IsSphere)
                {
                    LaplaceInSolidSphere * sphere = new LaplaceInSolidSphere();
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    sphere->SetNonLinearMesh(false);
                    sphere->Run( k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else
                {
                    // Nenhuma malha foi escolhida
                    DebugStop();
                }
            }
            

//            std::cout<< " grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
//            std::string filename("InputDataMeta");
//            std::string L2("L2.txt");
//            std::string Hdiv("Hdiv.txt");
//            std::string HdivData,L2Data;
//            HdivData = filename+str+Hdiv;
//            L2Data = filename+str+L2;
//            
//            PrintDebugMapForMathematica(HdivData,L2Data);
            
            std::cout<< " FIM - grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
            
        }
        
        saidaerros << "\n ----------------------------------------------------------------------------- " << endl;
//        fDebugMapHdiv.clear();
//        fDebugMapL2.clear();
    }
    
    std::cout<< " fim " << std::endl;
    
	return EXIT_SUCCESS;
}



