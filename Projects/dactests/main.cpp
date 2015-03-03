
#ifdef HAVE_CONFIG_H
#include <config.h>
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
#include "tpzquadratictrig.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"

#include "tpzchangeel.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
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

#include "pzhdivfull.h"
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


//#include "pyramidalmesh.h"

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;
using namespace pzshape;

int matId = 1;

int dirichlet = 0;
int neumann = 1;

int bc0 = -1;
int bc1 = -2;
int bc2 = -3;
int bc3 = -4;
int bc4 = -5;
int bc5 = -6;
int matskeleton = -7;

// just for print data
/** @brief Map used norms */
std::map<REAL,REAL> fDebugMapL2, fDebugMapHdiv;

int tetraedra_2[6][4]=
{
    {1,2,5,4},
    {4,7,3,2},
    {0,1,2,4},
    {0,2,3,4},
    {4,5,6,2},
    {4,6,7,2}
};

int piramide_2[6][5]=
{
    {0,1,2,3,8},
    {0,1,5,4,8},
    {1,2,6,5,8},
    {3,2,6,7,8},
    {0,3,7,4,8},
    {4,5,6,7,8}
};


bool MyDoubleComparer(REAL a, REAL b);

void GenerateNodes(TPZGeoMesh *gmesh, long nelem);
void GenerateNodesPyramid(TPZGeoMesh *gmesh, long nelem);

void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);
void RotateNode(TPZVec<STATE> &iCoords, REAL CounterClockwiseAngle, int &Axis);

TPZGeoMesh *GMesh(int dimensao, bool ftriang, int ndiv);

TPZGeoMesh *GMeshCirculoGeob(int dimensao, int ndiv);
TPZGeoMesh *GMeshCirculoGeobTriang(int dimensao, int ndiv);
TPZGeoMesh *GMeshCirculoQuad(int dimensao, int ndiv);

TPZGeoMesh *GMeshSphericalShell(int dimensao, bool triang, int ndiv);
TPZGeoMesh *GMeshSphericalShellGeob(int dimensao, int ndiv);
TPZGeoMesh *GMeshCilindricalMesh( int ndiv);
TPZGeoMesh *GMeshCilindricalMeshF( int ndiv);
TPZGeoMesh *GMeshCilindricalMeshR( int ndiv);

TPZGeoMesh *CreateOneCubo(int nref=0);
TPZGeoMesh *CreateOneCuboWithTetraedrons(long nelem=1, int MaterialId=1);
TPZGeoMesh *GMeshCubeWithPyramids(long nelem=1, int MaterialId=1);

TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);

void UniformRefinement(TPZGeoMesh* gmesh, int nDiv);
void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcess(TPZAnalysis &an, std::string plotfile);

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv);
void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);
/** @brief Prints debug map in Mathematica style */
void PrintDebugMapForMathematica(std::string filenameHdiv, std::string filenameL2);

//solucao exata
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);

//lado direito da equacao
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
void ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux);

//Para condicao de contorno de Dirichlet
void ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno de Neumann
void ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//criar elementos esqueleto
void AddWrap(TPZMultiphysicsElement *mfcel, int matskeleton, TPZStack< TPZStack< TPZMultiphysicsElement *, 7> > &ListGroupEl);

TPZGeoMesh * BasicForm(int n, REAL t, REAL dt);
void Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X);
void Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X);
void Parametricfunction3(const TPZVec<STATE> &par, TPZVec<STATE> &X);

int dim = 2;
REAL aa = 0.0;
REAL bb = 0.0;
REAL cc = 0.0;
REAL Epsilon = 0.4;
// tensor de permutacao
TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);

REAL const Pi = M_PI;//4.*atan(1.);



bool ftriang = false;//true;//
bool IsCube = false;
bool IsPrism = false;
bool IsTetra = false;
bool IsPiram = true;

//bool isspherical = true, isgeoblend = false, iscilindro = false;
bool isgeoblend = true, isspherical = false, iscilindro = false;
//bool iscilindro = true, isgeoblend = false, isspherical = false;

bool isH1 = false;


#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

#include "pztransfer.h"

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
//  gRefDBase.InitializeAllUniformRefPatterns();
//	gRefDBase.InitializeRefPatterns();

    int p = 1;
    int ndiv = 0;

    // Hard coded
    int sizeTP = dim;//isspherical ? 3 : dim;
    if (isspherical) {
        TP.Resize(sizeTP, sizeTP);
        InvTP.Resize(sizeTP, sizeTP);
    }
    for (int id = 0; id < sizeTP; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    ofstream saidaerro("../ErroPoissonHdivMalhaQuad.txt",ios::app);
    //int tipo = 1;
    //ofstream saidaerro("../ErroPoissonHdivMalhaTriang.txt",ios::app);
    
    for(p=1;p<5;p++)
    {
        int pq = p;
        int pp = p;
        
        for (ndiv=0; ndiv<6; ndiv++)
        {
            std::cout<< " INICIO - grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
            std::cout<< " Dimensao == " << dim << std::endl;

            TPZGeoMesh *gmesh;
            
            if (dim==2)
            {
                //TPZGeoMesh *gmesh2d = GMesh(2, ftriang, ndiv);
                if (isgeoblend) {

                    gmesh = GMeshCirculoGeob(2, ndiv);
                    //gmesh = GMeshCirculoGeobTriang(2, ndiv);
//                    {
//                        //  Print Geometrical Base Mesh
//                        std::ofstream Dummyfile("GeometricMesh2DH0.vtk");
//                        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
//                    }
//                    gmesh = GMeshCirculoGeob(2, 3);
//                    {
//                        //  Print Geometrical Base Mesh
//                        std::ofstream Dummyfile("GeometricMesh2DH3.vtk");
//                        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
//                    }
                    //return EXIT_SUCCESS;
                }
                else if(isspherical)
                {
                    std::string GridFileName;
                    GridFileName += "esferacoarse.dump";
                    
                    TPZReadGIDGrid GeometryInfo;
                    GeometryInfo.SetfDimensionlessL(1.0);
                    //gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
                    gmesh = GMeshSphericalShell(2, ftriang, ndiv);
                }
                else if (iscilindro)
                {
                    //gmesh = GMeshCilindricalMesh( ndiv);
                    //gmesh = GMeshCilindricalMeshF( ndiv);
                    gmesh = GMeshCilindricalMeshR( ndiv);
                }
                else
                {
                    gmesh = GMesh(2, ftriang, ndiv);
                }
                
                gmesh->SetDimension(dim);
//                ofstream arg1("gmesh2d.txt");
//                gmesh->Print(arg1);
                
                {
                    //  Print Geometrical Base Mesh
                    
                    // rotacao da malha geometrica
                    int Axis;
                    REAL theta, dump = 0.0;
                    
                    theta = 90.0;
                    Axis = 2;
                    RotateGeomesh(gmesh, theta*dump, Axis);
                    
                    std::ofstream Dummyfile("GeometrysphereGID.vtk");
                    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
                }
 
            }
            else // dim == 3
            {
                if (IsCube)
                {
                    // Tem que escolher a malha certa.
                    DebugStop();
                    //TPZGeoMesh *gmesh2d = GMesh(2, ftriang, ndiv);
                    
                    //REAL layerthickness = 1.;
                    //TPZExtendGridDimension extend(gmesh2d,layerthickness);
                    //TPZGeoMesh *gmesh3d = extend.ExtendedMesh(1,bc0,bc5);
                    //gmesh = gmesh3d;

                }
                else if(IsTetra)
                {
                    REAL dndiv = ndiv;
                    int nref = (int) pow(2., dndiv);
                    
                    gmesh = CreateOneCuboWithTetraedrons(nref, matId);
                }
                else if (IsPrism)
                {
                    // Tem que escolher a malha certa.
                    DebugStop();
                }
                else if(IsPiram)
                {
                    REAL dndiv = ndiv;
                    int nref = (int) pow(2., dndiv);
                    gmesh = GMeshCubeWithPyramids( nref);
//                    PyramidalMesh(gmesh,ndiv);
                }
                else
                {
                    // Nenhuma malha foi escolhida
                    DebugStop();
                }
                gmesh->SetDimension(dim);
                ofstream arg1("gmesh3d.txt");
                gmesh->Print(arg1);
                {
                    //  Print Geometrical Base Mesh
                    std::ofstream Dummyfile("GeometricMesh3D.vtk");
                    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
                }

            }
            
            //            int n = ndiv;
            //            REAL t= 0.0;
            //            REAL dt = 1./(pow(2.0, 1.0*n));
            //            int nt= int(1/dt);
            //            TPZGeoMesh *gmesh =  BasicForm(nt,t,dt);
            
            //            ofstream arg1("gmesh.txt");
            //            gmesh->Print(arg1);
            
            
            
            // rotacao da malha geometrica
            int Axis;
            REAL theta, dump = 0.0;
            
            theta = 45.0;
            Axis = 1;
            RotateGeomesh(gmesh, theta*dump, Axis);
            
            theta = -45.0;
            Axis = 2;
            RotateGeomesh(gmesh, theta*dump, Axis);
            
            theta = 45.0;
            Axis = 3;
            RotateGeomesh(gmesh, theta*dump, Axis);
            
            TPZCompMesh *cmesh2 = CMeshPressure(gmesh, pp, dim);
            TPZCompMesh *cmesh1 = CMeshFlux(gmesh, pq, dim);

// Um teste para a solucao via H1, sem hdiv
            if (isH1) {
                TPZCompMesh *cmeshH1 = CMeshH1(gmesh, pq, dim);
                TPZAnalysis anh1(cmeshH1, true);
                
//                TPZFMatrix<REAL> IntegraloverS;
//                
//                anh1.Assemble();
//                IntegraloverS=anh1.Rhs();
//                std::stringstream fileout;
//                IntegraloverS.Print("IntegraloverS = ",std::cout,EMathematicaInput);
                
                
                
                
                SolveSyst(anh1, cmeshH1);
                
                stringstream refh1,grauh1;
                grauh1 << p;
                refh1 << ndiv;
                string strgh1 = grauh1.str();
                string strrh1 = refh1.str();
                std::string plotnameh1("OurSolutionH1");
                std::string Grauh1("P");
                std::string Refh1("H");
                std::string VTKh1(".vtk");
                std::string plotDatah1;
                plotDatah1 = plotnameh1+Grauh1+strgh1+Refh1+strrh1+VTKh1;
                std::string plotfileh1(plotDatah1);
                
                PosProcess(anh1, plotfileh1);
                
                std::cout<< " FIM - H1 - grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
                
                return EXIT_SUCCESS;
            }
            
// exit
            
            
//            ofstream arg1("cmeshflux.txt");
//            cmesh1->Print(arg1);
//
//            ofstream arg2("cmeshpressure.txt");
//            cmesh2->Print(arg2);
//            
//            ofstream arg4("gmesh2.txt");
//            gmesh->Print(arg4);
//            
            
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            
            
//            int Prows =cmesh2->Solution().Rows();
//            TPZFMatrix<STATE> alphaP(Prows,1,0.0);
//            alphaP(26,0)=1.0;
//            cmesh2->LoadSolution(alphaP);
//            TPZAnalysis anP(cmesh2,false);
//            std::string plotfile2("AlphasPressure.vtk");
//            PosProcess(anP, plotfile2);
            
            
            
            TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
            
            //TestMesh(mphysics);
//            ofstream arg5("cmeshmultiphysics.txt");
//            mphysics->Print(arg5);
            
            
//            TPZCompEl *cel = mphysics->Element(0);
//            TPZElementMatrix ek,ef;
//            cel->CalcStiff(ek, ef);
            
            TPZAnalysis an(mphysics, true);

            SolveSyst(an, mphysics);
//            std::string plotfile("OurSolution.vtk");
//            PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
            
//            TPZAutoPointer< TPZMatrix<REAL> > matK;
//            TPZFMatrix<STATE> fvec;
//            matK=an.Solver().Matrix();
//            fvec = an.Rhs();
//            
//            std::stringstream sout;
//            matK->Print("matK = ", std::cout,EMathematicaInput);
//            fvec.Print("fvec = ", std::cout,EMathematicaInput);
            
            stringstream ref,grau;
            grau << p;
            ref << ndiv;
            string strg = grau.str();
            string strr = ref.str();
            std::string plotname("OurSolutionMeta");
            std::string Grau("P");
            std::string Ref("H");
            std::string VTK(".vtk");
            std::string plotData;
            plotData = plotname+Grau+strg+Ref+strr+VTK;
            std::string plotfile(plotData);
            
            PosProcessMultphysics(meshvec,  mphysics, an, plotfile);



            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            TPZVec<REAL> erros;
            
            //            TPZManVector<REAL,3> myerrors(3,0.);
            //            an.SetExact(SolExata);
            //            an.PostProcessError(myerrors);
            
            //            saidaerro << "Valor de epsilone " << EPSILON << std::endl;
            //            saidaerro << "Numero de threads " << numthreads << std::endl;
            saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
            ErrorHDiv(cmesh1, saidaerro, p, ndiv);
            
            saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
            ErrorL2(cmesh2, saidaerro, p, ndiv);

            std::cout << "Postprocessed\n";
            
            //Plot da solucao aproximada
            // para estudo do ero nao precisa gerar dados de saida
            //PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
            
            
            stringstream ss;
            ss << p;
            string str = ss.str();
            
            std::cout<< " grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
            std::string filename("InputDataMeta");
            std::string L2("L2.txt");
            std::string Hdiv("Hdiv.txt");
            std::string HdivData,L2Data;
            HdivData = filename+str+Hdiv;
            L2Data = filename+str+L2;
            
            PrintDebugMapForMathematica(HdivData,L2Data);
            
            std::cout<< " FIM - grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
            
        }
        fDebugMapHdiv.clear();
        fDebugMapL2.clear();
    }
    
    std::cout<< " fim " << std::endl;
    
	return EXIT_SUCCESS;
}

TPZGeoMesh * BasicForm(int n, REAL t, REAL dt)
{
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
    GeoMesh1->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh1->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
    GeoMesh1->BuildConnectivity();
    GeoMesh1->SetDimension(0);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew1.txt");
        GeoMesh1->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh1,Dummyfile, true);
    }
    
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh1);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc = new TPZDummyFunction<STATE>(Parametricfunction);
    CreateGridFrom.SetParametricFunction(ParFunc);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh2 = CreateGridFrom.ComputeExtrusion(t, dt, n);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew2.txt");
        GeoMesh2->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew2.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh2,Dummyfile, true);
    }
    
    
    
    TPZHierarquicalGrid CreateGridFrom2(GeoMesh2);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc2 = new TPZDummyFunction<STATE>(Parametricfunction2);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dt, n);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew3.txt");
        GeoMesh3->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew3.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh3,Dummyfile, true);
    }
    
    
    
    TPZHierarquicalGrid CreateGridFrom3(GeoMesh3);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc3 = new TPZDummyFunction<STATE>(Parametricfunction3);
    CreateGridFrom3.SetParametricFunction(ParFunc3);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh4 = CreateGridFrom3.ComputeExtrusion(t, dt, n);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew4.txt");
        GeoMesh4->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew4.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh4,Dummyfile, true);
    }
    return GeoMesh4;
}

TPZGeoMesh *GMesh(int d, bool ftriang, int ndiv)
{
    
    if(dim!=2)
    {
        std::cout << "dimensao errada" << std::endl;
        dim = 2;
        DebugStop();
    }
    
    int Qnodes =  4;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    gmesh->SetDimension(dim);
    
    TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolTriang(3);
    TPZVec <long> TopolLine(2);
    TPZVec <long> TopolPoint(1);
    
    //indice dos nos
    long id = 0;
    //    REAL valx;
    //    for(int xi = 0; xi < Qnodes/2; xi++)
    //    {
    //        valx = xi*Lx;
    //        Node[id].SetNodeId(id);
    //        Node[id].SetCoord(0 ,valx );//coord X
    //        Node[id].SetCoord(1 ,0. );//coord Y
    //        gmesh->NodeVec()[id] = Node[id];
    //        id++;
    //    }
    //
    //    for(int xi = 0; xi < Qnodes/2; xi++)
    //    {
    //        valx = Lx - xi*Lx;
    //        Node[id].SetNodeId(id);
    //        Node[id].SetCoord(0 ,valx );//coord X
    //        Node[id].SetCoord(1 ,Ly);//coord Y
    //        gmesh->NodeVec()[id] = Node[id];
    //        id++;
    //    }
    //
    TPZManVector<REAL,3> coord(2,0.);
    int in = 0;
    //c0
    coord[0] = -1.0;//0.0; //
    coord[1] = -1.0;//0.0; //
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c1
    coord[0] =   1.0;//1.0;//
    coord[1] =  -1.0;//0.0;//
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c2
    coord[0] =   1.0;//0.0;//
    coord[1] =   1.0;//1.0;//
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c3
    coord[0] =  -1.0;//1.0;//
    coord[1] =   1.0;//1.0;//
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    //indice dos elementos
    id = 0;
    
    if(ftriang==true) // triangulo
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 0;
        TopolTriang[1] = 2;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
        id++;
    }
    else{
        
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
    }
    
    gmesh->BuildConnectivity();
    
    /// gmesh para aqui
    
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < ndiv; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Geometrica Inicial\n ";
    //        gmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    return gmesh;
}

TPZGeoMesh *CreateOneCuboWithTetraedrons(long nelem, int MaterialId)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh,nelem);
    
    for (long i=0; i<nelem; i++) {
        for (long j=0; j<nelem; j++) {
            for (long k=0; k<nelem; k++) {
                TPZManVector<long,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef LOG4CXX
                if(logdata->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logdata, sout.str())
                }
#endif
                for (int el=0; el<6; el++)
                {
                    TPZManVector<long,4> elnodes(4);
                    long index;
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra_2[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    
    // Boundary Conditions
    const int numelements = gmesh->NElements();
    //    const int bczMinus = -3, bczplus = -2, bcids = -1;
    //    const int bczMinus = -1, bczplus = -1, bcids = -1;
    
    for(int el=0; el<numelements; el++)
    {
        TPZManVector <TPZGeoNode,4> Nodefinder(4);
        TPZManVector <REAL,3> nodecoord(3);
        TPZGeoEl *tetra = gmesh->ElementVec()[el];
        TPZVec<long> ncoordVec(0); long sizeOfVec = 0;
        
        // na face z = 0
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bc0);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 0
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bc1);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 1
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bc2);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 1
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bc3);
        }
        
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 0
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bc4);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face z = 1
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bc5);
        }
        
        
        
    }
    
    return gmesh;
}

bool MyDoubleComparer(REAL a, REAL b)
{
    if (IsZero(a-b)){
        return true;
    }
    else{
        return false;
    }
}

void GenerateNodes(TPZGeoMesh *gmesh, long nelem)
{
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (long i=0; i<=nelem; i++) {
        for (long j=0; j<=nelem; j++) {
            for (long k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
}

void GenerateNodesPyramid(TPZGeoMesh *gmesh, long nelem)
{
    long sizenodevec = (nelem+1)*(nelem+1)*(nelem+1)+(nelem*nelem*nelem);
    gmesh->NodeVec().Resize(sizenodevec);
    int posicao = 0;
    for (long i=0; i<=nelem; i++) {
        for (long j=0; j<=nelem; j++) {
            for (long k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                posicao = i*(nelem+1)*(nelem+1)+j*(nelem+1)+k;
                gmesh->NodeVec()[posicao].Initialize(x, *gmesh);
                 
            }
        }
    }
    for (long i=0; i<nelem; i++) {
        for (long j=0; j<nelem; j++) {
            for (long k=0; k<nelem; k++) {
                TPZManVector<REAL,3> xc(3);
                xc[0] = k*1./nelem + (0.5)*1./nelem;
                xc[1] = j*1./nelem + (0.5)*1./nelem;
                xc[2] = i*1./nelem + (0.5)*1./nelem;
                posicao++;
                gmesh->NodeVec()[posicao].Initialize(xc, *gmesh);
            }
        }
    }
}


TPZGeoMesh *CreateOneCubo(int nref)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int nnodes = 8;
    int idf0=-1;
    int idf1=-2;
    int idf2=-3;
    int idf3=-4;
    int idf4=-5;
    int idf5=-6;
    
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    
    //cubo [0,1]ˆ3
    //    //c0
    //    coord[0] = 0.0;
    //    coord[1] = 0.0;
    //    coord[2] = 0.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in++);
    //    //c1
    //    coord[0] =  1.0;
    //    coord[1] = 0.0;
    //    coord[2] = 0.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in++);
    //    //c2
    //    coord[0] =  1.0;
    //    coord[1] =  1.0;
    //    coord[2] = 0.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in++);
    //    //c3
    //    coord[0] = 0.0;
    //    coord[1] =  1.0;
    //    coord[2] = 0.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in++);
    //
    //    //c4
    //    coord[0] = 0.0;
    //    coord[1] = 0.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in++);
    //    //c5
    //    coord[0] =  1.0;
    //    coord[1] = 0.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in++);
    //    //c6
    //    coord[0] =  1.0;
    //    coord[1] =  1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in++);
    //    //c7
    //    coord[0] = 0.0;
    //    coord[1] =  1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in++);
    
    // cubo [-1,1]^3
    //c0
    coord[0] = -1.0;
    coord[1] = -1.0;
    coord[2] = -1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c1
    coord[0] =  1.0;
    coord[1] = -1.0;
    coord[2] = -1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c2
    coord[0] =  1.0;
    coord[1] =  1.0;
    coord[2] = -1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c3
    coord[0] = -1.0;
    coord[1] =  1.0;
    coord[2] = -1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    
    //c4
    coord[0] = -1.0;
    coord[1] = -1.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c5
    coord[0] =  1.0;
    coord[1] = -1.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c6
    coord[0] =  1.0;
    coord[1] =  1.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c7
    coord[0] = -1.0;
    coord[1] =  1.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    
    
    int index = 0;
    
    TPZVec<long> TopologyQuad(4);
    
    // bottom
    TopologyQuad[0]=0;
    TopologyQuad[1]=1;
    TopologyQuad[2]=2;
    TopologyQuad[3]=3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf0,*gmesh);
    index++;
    
    // Front
    TopologyQuad[0]=0;
    TopologyQuad[1]=1;
    TopologyQuad[2]=5;
    TopologyQuad[3]=4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf1,*gmesh);
    index++;
    
    // Rigth
    TopologyQuad[0]=1;
    TopologyQuad[1]=2;
    TopologyQuad[2]=6;
    TopologyQuad[3]=5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf2,*gmesh);
    index++;
    // Back
    TopologyQuad[0]=3;
    TopologyQuad[1]=2;
    TopologyQuad[2]=6;
    TopologyQuad[3]=7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf3,*gmesh);
    index++;
    
    // Left
    TopologyQuad[0]=0;
    TopologyQuad[1]=3;
    TopologyQuad[2]=7;
    TopologyQuad[3]=4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf4,*gmesh);
    index++;
    
    // Top
    TopologyQuad[0]=4;
    TopologyQuad[1]=5;
    TopologyQuad[2]=6;
    TopologyQuad[3]=7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf5,*gmesh);
    index++;
    
    TPZManVector<long,8> TopolCubo(8,0);
    TopolCubo[0] = 0;
    TopolCubo[1] = 1;
    TopolCubo[2] = 2;
    TopolCubo[3] = 3;
    TopolCubo[4] = 4;
    TopolCubo[5] = 5;
    TopolCubo[6] = 6;
    TopolCubo[7] = 7;
    
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index, TopolCubo, matId, *gmesh);
    index++;
    
    gmesh->BuildConnectivity();
    
    /// gmesh para aqui
    
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("SingleCubeWithBcs.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZGeoMesh *GMeshCubeWithPyramids(long nelem, int MaterialId)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodesPyramid(gmesh,nelem);
    
    int posicao = 0;
    for (long i=0; i<=nelem; i++) {
        for (long j=0; j<=nelem; j++) {
            for (long k=0; k<=nelem; k++) {
                posicao = i*(nelem+1)*(nelem+1)+j*(nelem+1)+k;
            }
        }
    }

    
    for (long i=0; i<nelem; i++) {
        for (long j=0; j<nelem; j++) {
            for (long k=0; k<nelem; k++) {
                TPZManVector<long,9> nodes(9,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[8] = posicao + (k)*(nelem)*(nelem)+(j)*(nelem)+i+1;
#ifdef LOG4CXX
                if(logdata->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Pyramidal nodes " << nodes;
                    LOGPZ_DEBUG(logdata, sout.str())
                }
#endif
                for (int el=0; el<6; el++)
                {
                    TPZManVector<long,5> elnodes(5);
                    long index;
                    for (int il=0; il<5; il++) {
                        elnodes[il] = nodes[piramide_2[el][il]];
                    }
                    gmesh->CreateGeoElement(EPiramide, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    
    // Boundary Conditions
    const int numelements = gmesh->NElements();
    //    const int bczMinus = -3, bczplus = -2, bcids = -1;
    //    const int bczMinus = -1, bczplus = -1, bcids = -1;
    
    for(int el=0; el<numelements; el++)
    {
        TPZManVector <TPZGeoNode,5> Nodefinder(5);
        TPZManVector <REAL,3> nodecoord(3);
        TPZGeoEl *piramide = gmesh->ElementVec()[el];
        TPZVec<long> ncoordVec(0); long sizeOfVec = 0;
        
        // na face z = 0
        for (int i = 0; i < 5; i++)
        {
            long pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc0);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 0
        for (int i = 0; i < 5; i++)
        {
            long pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc1);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 1
        for (int i = 0; i < 5; i++)
        {
            long pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc2);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 1
        for (int i = 0; i < 5; i++)
        {
            long pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc3);
        }
        
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 0
        for (int i = 0; i < 5; i++)
        {
            long pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc4);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face z = 1
        for (int i = 0; i < 5; i++)
        {
            long pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc5);
        }
        
        
        
    }
    
    
    std::ofstream out("CubeWithBcs.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

void UniformRefinement(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    //	gmesh->BuildConnectivity();
}

TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
//    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
//    TPZAutoPointer<TPZFunction<STATE> > forcef;
//    forcef = new TPZDummyFunction<STATE>(Forcing);
//    material->SetForcingFunction(forcef);
    
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingH1);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    if( dim == 3 ) { BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);}
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    if( dim == 3 ) { BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);}
    
    cmesh->InsertMaterialObject(mat);
 
    cmesh->SetAllCreateFunctionsContinuous(); 
    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond0); }
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond5); }
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}


TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    //  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    if( dim == 3 ) { BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);}
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    if( dim == 3 ) { BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);}
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond0); }
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond5); }
    
    cmesh->SetDefaultOrder(pOrder);
    
    
    if (!isgeoblend) {
        TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(matskeleton, dim-1, 1);
        TPZMaterial * mat2(matskelet);
        cmesh->InsertMaterialObject(mat2);
    }
    
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    return cmesh;

}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    material->NStateVariables();
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    //    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    //    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    //    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    //    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    //    TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    //    TPZMaterial * BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
    //
    //    cmesh->InsertMaterialObject(BCond0);
    //    cmesh->InsertMaterialObject(BCond1);
    //    cmesh->InsertMaterialObject(BCond2);
    //    cmesh->InsertMaterialObject(BCond3);
    //    cmesh->InsertMaterialObject(BCond4);
    //    cmesh->InsertMaterialObject(BCond5);
    
    cmesh->SetDefaultOrder(pOrder);
    //cmesh->SetDimModel(dim);
    
    bool h1function = true;//false em triangulo
    if(h1function){
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    else{
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
        newnod.SetLagrangeMultiplier(1);
    }
    
    if(!h1function)
    {
        
        int nel = cmesh->NElements();
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            celdisc->SetConstC(1.);
            celdisc->SetCenterPoint(0, 0.);
            celdisc->SetCenterPoint(1, 0.);
            celdisc->SetCenterPoint(2, 0.);
            celdisc->SetTrueUseQsiEta();
            //celdisc->SetFalseUseQsiEta();
            
            //            TPZVec<REAL> qsi(3,0.);
            //            qsi[0] = 0.5;
            //            qsi[1] = 0.5;
            //            TPZFMatrix<REAL> phi;
            //            TPZFMatrix<REAL> dphi;
            //            celdisc->Shape(qsi, phi,dphi);
            //            phi.Print("phi = ");
            
            
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                if(ftriang==true) celdisc->SetTotalOrderShape();
                else celdisc->SetTensorialShape();
            }
            
        }
    }
    
#ifdef DEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_2 pressure\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
    return cmesh;

}


#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
{
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim = gmesh->Dimension();
    
    bool hdivantigo;
    
    //TPZMatPoissonD3 *material = new TPZMatPoissonD3(matId,dim); hdivantigo = true; // nesse material tem que ser true
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(matId,dim); hdivantigo = false; // nesse material tem que ser false
    //TPZMixedPoisson *material = new TPZMixedPoisson(matId,dim); hdivantigo = true;; // nesse material tem que ser true
    
    //incluindo os dados do problema
    if (hdivantigo) {
        TPZFNMatrix<2,REAL> PermTensor(dim,dim,0.);
        TPZFNMatrix<2,REAL> InvPermTensor(dim,dim,0.);
        
        for (int i=0; i<dim; i++)
        {
            PermTensor(i,i) = 1.0;
        }
        InvPermTensor=PermTensor;
        material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    }
    
    //incluindo os dados do problema
    TPZFNMatrix<2,REAL> PermTensor(dim,dim,0.);
    TPZFNMatrix<2,REAL> InvPermTensor(dim,dim,0.);
    
    PermTensor = TP;
    InvPermTensor = InvTP;
    
    
    material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    
    solexata = new TPZDummyFunction<STATE>(SolExata);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    if (dim==3)
    {
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZAutoPointer<TPZFunction<STATE> > FBCond0 = new TPZDummyFunction<STATE>(ForcingBC0D);
        BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
        BCond0->SetForcingFunction(FBCond0);
    }

    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1D);
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
//    BCond1 = material->CreateBC(mat, bc1,neumann, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2D);
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    BCond2->SetForcingFunction(FBCond2);
//    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2N);
//    BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
//    BCond2->SetForcingFunction(FBCond2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3D);
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    BCond3->SetForcingFunction(FBCond3);
//    BCond3 = material->CreateBC(mat, bc3,neumann, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4D);
    BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    BCond4->SetForcingFunction(FBCond4);
//    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4N);
//    BCond4 = material->CreateBC(mat, bc4,neumann, val1, val2);
//    BCond4->SetForcingFunction(FBCond4);
    
    if (dim==3)
    {
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZAutoPointer<TPZFunction<STATE> > FBCond5 = new TPZDummyFunction<STATE>(ForcingBC5D);
        BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
        BCond5->SetForcingFunction(FBCond5);
    }
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    if( dim == 3 ) { mphysics->InsertMaterialObject(BCond0); }
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    if( dim == 3 ) { mphysics->InsertMaterialObject(BCond5); }
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    if(hdivantigo){ // ! colocada para testar material que fiz, o antigo
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    }
    //Creating multiphysic elements containing skeletal elements.
    else
    {
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        
        //        TPZMaterial * skeletonEl = material->CreateBC(mat, matskeleton, 3, val1, val2);
        //        mphysics->InsertMaterialObject(skeletonEl);
        
        //        TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(matskeleton, dim-1, 1);
        //        TPZMaterial * mat2(matskelet);
        //        mphysics->InsertMaterialObject(mat2);
        
        int nel = mphysics->ElementVec().NElements();
        TPZStack< TPZStack< TPZMultiphysicsElement *,7> > wrapEl;
        for(int el = 0; el < nel; el++)
        {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(mphysics->Element(el));
            if(mfcel->Dimension()==dim) TPZBuildMultiphysicsMesh::AddWrap(mfcel, matId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
        }
        meshvec[0]->CleanUpUnconnectedNodes();
        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
        //------- Create and add group elements -------
        long index, nenvel;
        nenvel = wrapEl.NElements();
        for(int ienv=0; ienv<nenvel; ienv++){
            TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
            nel = wrapEl[ienv].NElements();
            for(int jel=0; jel<nel; jel++){
                elgr->AddElement(wrapEl[ienv][jel]);
            }
        }
    }
    
    return mphysics;

}

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    cout <<"Numero de equacoes "<< fCmesh->NEquations()<< endl;
    
	bool isdirect = true;
    bool simetrico = true;
    if (isdirect)
    {
        if (simetrico)
        {
            //TPZBandStructMatrix full(fCmesh);
            TPZSkylineStructMatrix skylstr(fCmesh); //caso simetrico
            //    TPZSkylineNSymStructMatrix full(fCmesh);
            an.SetStructuralMatrix(skylstr);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt); //caso simetrico
            an.SetSolver(step);
            an.Run();
        }
        else
        {
            TPZBandStructMatrix full(fCmesh);
            an.SetStructuralMatrix(full);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELU);
            an.SetSolver(step);
            an.Run();
        }
                
    }
    else
    {
        TPZSkylineStructMatrix skylstr(fCmesh); //caso simetrico
        skylstr.SetNumThreads(10);
        an.SetStructuralMatrix(skylstr);
        
        TPZAutoPointer<TPZMatrix<STATE> > matbeingcopied = skylstr.Create();
        TPZAutoPointer<TPZMatrix<STATE> > matClone = matbeingcopied->Clone();
        
        TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>(matClone);
        TPZStepSolver<STATE> *Solver = new TPZStepSolver<STATE>(matbeingcopied);
        precond->SetReferenceMatrix(matbeingcopied);
        precond->SetDirect(ELDLt);
        Solver->SetGMRES(20, 20, *precond, 1.e-18, 0);
        //        Solver->SetCG(10, *precond, 1.0e-10, 0);
        an.SetSolver(*Solver);
        an.Run();
    }
    
    

}

void PosProcess(TPZAnalysis &an, std::string plotfile){
    TPZManVector<std::string,10> scalnames(1), vecnames(0);
    scalnames[0] = "Solution";
    //const int dim = 3;
    int div = 2;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
//    std::ofstream out("malhaNormal.txt");
//    an.Print("nothing",out);
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZManVector<std::string,10> scalnames(2), vecnames(2);
    vecnames[0]  = "Flux";
    vecnames[1]  = "ExactFlux";
    //    vecnames[1]  = "GradFluxX";
    //    vecnames[2]  = "GradFluxY";
    //    vecnames[2]  = "GradFluxz";
    scalnames[0] = "Pressure";
    scalnames[1] = "ExactPressure";
    //scalnames[2] = "Rhs";
    
    //const int dim = 3;
    int div = 2;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    
//    std::ofstream out("malha.txt");
//    an.Print("nothing",out);

    
    //    mphysics->Solution().Print("Solucao");
    
}

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    
    
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica
    
    
    REAL x,y,z;
    // Defining data
    REAL a=M_PI/2.0;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    // Hard coded
    if (isspherical) {
        flux.Resize(3, 1);
    }
    
    
    REAL r = sqrt(x*x+y*y+z*z);
//    //REAL theta = (atan2(sqrt(x*x+y*y),z));
//    REAL theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
//    REAL phi = atan2(y,x);
//    
//    solp[0] = (a-theta)*sin(theta)*sin(theta);
//    flux(0,0)= -(cos(phi)*cos(theta)*(2.0*(a - theta)*cos(theta) - sin(theta))*sin(theta))/r;
//    flux(1,0)= -(cos(theta)*sin(phi)*(2.0*(a - theta)*cos(theta) - sin(theta))*sin(theta))/r;
//    flux(2,0)= -(sin(theta)*sin(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta)))/r;
    
    /* Problem with phi */
//    solp[0] = (a-theta)*phi*phi*sin(theta)*sin(theta);
//    flux(0,0)= (phi*sin(theta)*(2.0*(-a + theta)*sin(phi)-phi*cos(phi)*cos(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta))))/r;
//    flux(1,0)= (phi*sin(theta)*(2.0*(a - theta)*(cos(phi) + phi*cos(theta)*cos(theta)*sin(phi)) -phi*cos(theta)*sin(phi)*sin(theta)))/r;
//    flux(2,0)= (phi*phi*sin(theta)*sin(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta)))/r;
    
    // Problema na circunferencia
    r = sqrt(x*x+y*y);
    if (r<1.0) {
        solp[0] = 3.0*exp(1.0/(x*x+y*y-1.0));
        flux(0,0) = 6.0*x*exp(1.0/(x*x+y*y-1.0))/((x*x+y*y-1.0)*(x*x+y*y-1.0));
        flux(1,0) = 6.0*y*exp(1.0/(x*x+y*y-1.0))/((x*x+y*y-1.0)*(x*x+y*y-1.0));
        //flux(2,0) = 0.0;
    }
    else
    {
        solp[0] = 0.0;
        flux(0,0) = 0.0;
        flux(1,0) = 0.0;
      //  flux(2,0) = 0.0;
    }

    
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){

    
    REAL x,y,z;
    // Defining data
    REAL a=M_PI/2.0;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
   
    REAL r = sqrt(x*x+y*y+z*z);
    //REAL theta = (atan2(sqrt(x*x+y*y),z));
    REAL theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    REAL phi = atan2(y,x);
    
    ff[0] = (-1.0/(2.0*r*r))*(2.0*(a-theta)*(1.0+3.0*cos(2.0*theta))- 5.0*sin(2.0*theta));
//
//    std::cout << " x =" << x << "y =" << y << "z =" << z << std::endl;
//    std::cout << " Theta =" << theta << "phi =" << phi << "r =" << r << std::endl;
//    std::cout << " ff[0] =" << ff[0] << "a =" << a << std::endl;
    
//    /* Problem with phi */
//      ff[0] = (-1.0/(r*r))*((a-theta)*(2.0-phi*phi+3.0*phi*phi*cos(2.0*theta))- 5.0*phi*phi*cos(theta)*sin(theta));
    
    // Problema na circunferencia
    r = sqrt(x*x+y*y);
    if (r<1.0) {
        ff[0] = -12.0*(-1.0 + x*x*x*x + y*y + y*y*y*y + x*x*(1.0+2.0*y*y))*exp(1.0/(x*x+y*y-1.0))/((x*x+y*y-1.0)*(x*x+y*y-1.0)*(x*x+y*y-1.0)*(x*x+y*y-1.0));
    }
    else
    {
        ff[0] = 0.0;
    }


}

void ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
{
    REAL x,y,z;
    // Defining data
    REAL a=M_PI/2.0;
    flux.Resize(3, 1);
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = (atan2(sqrt(x*x+y*y),z));
    REAL phi = atan2(y,x);
    
//    ff[0] = (1.0/(2.0*r*r))*(2.0*(a-theta)*(1.0+3.0*cos(2.0*theta))- 5.0*sin(2.0*theta));
//    flux(0,0)= -(cos(phi)*cos(theta)*(2.0*(a - theta)*cos(theta) - sin(theta))*sin(theta))/r;
//    flux(1,0)= -(cos(theta)*sin(phi)*(2.0*(a - theta)*cos(theta) - sin(theta))*sin(theta))/r;
//    flux(2,0)= -(sin(theta)*sin(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta)))/r;
    
    /* Problem with phi */
    ff[0] = (1.0/(r*r))*((a-theta)*(2.0-phi*phi+3.0*phi*phi*cos(2.0*theta)- 5.0*phi*phi*cos(theta)*sin(theta)));
    flux(0,0)= -(phi*sin(theta)*(2.0*(-a + theta)*sin(phi)-phi*cos(phi)*cos(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta))))/r;
    flux(1,0)= -(phi*sin(theta)*(2.0*(a - theta)*(cos(phi) + phi*cos(theta)*cos(theta)*sin(phi)) -phi*cos(theta)*sin(phi)*sin(theta)))/r;
    flux(2,0)= -(phi*phi*sin(theta)*sin(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta)))/r;
    
}

void ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    solp[0]=0.0; // When theta-> pi/2;
}

void ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0]=0.0;//0.5; // When theta-> pi/2;
}


void ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0]=0.0;//0.5; // When theta-> pi/2;
}


void ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0]=0.0; // When theta-> pi/2;
}


void ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0]=0.0;//0.5; // When theta-> pi/2;
}


void ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0]=0.0; // When theta-> pi/2;
}



void ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){

    DebugStop();
}

void ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
    
}

void ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
        DebugStop();
}

void ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    DebugStop();
}

void ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
}

void ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
}

void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv)
{
    long nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<STATE,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    out << "L2 Norm for flux - "<< endl; //L2 Norm for divergence - Hdiv Norm for flux " << endl;
    out <<  setw(16) << sqrt(globalerrors[1]) <<endl;// setw(25)  << sqrt(globalerrors[2]) << setw(21)  << sqrt(globalerrors[3]) << endl;
    //
    //    out << "L2 Norm for flux = "    << sqrt(globalerrors[1]) << endl;
    //    out << "L2 Norm for divergence = "    << sqrt(globalerrors[2])  <<endl;
    //    out << "Hdiv Norm for flux = "    << sqrt(globalerrors[3])  <<endl;
    //
    fDebugMapHdiv.insert(std::pair<REAL, REAL> (ndiv,sqrt(globalerrors[1])));
}

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv)
{
    long nel = l2mesh->NElements();
    //int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    globalerrors.Fill(0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<STATE,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, NULL);
        int nerr = elerror.size();
        // globalerrors.resize(nerr);
#ifdef LOG4CXX
        if (logdata->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logdata, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with L2 space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    out << "L2 Norm = "    << sqrt(globalerrors[1]) << endl;
    
    fDebugMapL2.insert(std::pair<REAL, REAL> (ndiv,sqrt(globalerrors[1])));

}

void PrintDebugMapForMathematica(std::string filenameHdiv, std::string filenameL2)
{
    std::ofstream outHdiv(filenameHdiv.c_str());
    std::ofstream outL2(filenameL2.c_str());
    
    for (std::map<REAL,REAL>::const_iterator it = fDebugMapL2.begin(); it != fDebugMapL2.end(); it++) {
        outL2 << it->first << "   " << it->second << std::endl;
    }
    outL2.close();
    
    
    for (std::map<REAL,REAL>::const_iterator it = fDebugMapHdiv.begin(); it != fDebugMapHdiv.end(); it++) {
        outHdiv <<  it->first << "   " << it->second << std::endl;
    }
    outHdiv.close();
}

// malha que representa um circulo de raio R formado por 4 quadrilateros com geoblends
TPZGeoMesh *GMeshCirculoGeob(int dimensao, int ndiv)
{
    if(dim != 2)
    {
        DebugStop();
    }
    
    TPZGeoMesh * gmesh;// = new TPZGeoMesh;
    
    // description of Geometry and application
    // 2D Cylindrical Domain boundaries
    //int matId = 1;
    int arc1 = bc1; // -1;
    int arc2 = bc2; // -2;
    int arc3 = bc3; // -3;
    int arc4 = bc4; // -4;
    
    int nodenumber = 13;
    REAL ModelRadius = 1.5;
    REAL ModelRadiusInt2 = 0.5;//ModelRadius/2.;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordantes for Arc3D 1
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //6
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //7
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //8
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //9
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadiusInt2);//coord Y
    id++;
    //10
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //11
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadiusInt2);//coord Y
    id++;
    //12
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    
    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(4);
    
    // Create Geometrical Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 9;
    nodeindex[3] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 10;
    nodeindex[3] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 11;
    nodeindex[3] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 8;
    nodeindex[3] = 11;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #5
    nodeindex[0] = 8;
    nodeindex[1] = 9;
    nodeindex[2] = 10;
    nodeindex[3] = 11;
    new TPZGeoElRefPattern < pzgeom::TPZGeoQuad >  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
//    // Create Geometrical Quad #5
//    nodeindex.resize(3);
//    nodeindex[0] = 8;
//    nodeindex[1] = 9;
//    nodeindex[2] = 12;
//    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, matId,*gmesh);
//    elementid++;
//    // Create Geometrical Quad #5
//    nodeindex[0] = 9;
//    nodeindex[1] = 10;
//    nodeindex[2] = 12;
//    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, matId,*gmesh);
//    elementid++;
//    // Create Geometrical Quad #5
//    nodeindex[0] = 10;
//    nodeindex[1] = 11;
//    nodeindex[2] = 12;
//    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, matId,*gmesh);
//    elementid++;
//    // Create Geometrical Quad #5
//    nodeindex[0] = 11;
//    nodeindex[1] = 8;
//    nodeindex[2] = 12;
//    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, matId,*gmesh);
//    elementid++;
    
    // Definition of Arc coordenates
    nodeindex.resize(3);
    // Create Geometrical Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2, *gmesh);
    
    
    
    gmesh->BuildConnectivity();
    
    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("CurvoCirculo.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZGeoMesh *GMeshCirculoGeobTriang(int dimensao, int ndiv)
{
    if(dim != 2)
    {
        DebugStop();
    }
    
    TPZGeoMesh * gmesh;// = new TPZGeoMesh;
    
    // description of Geometry and application
    // 2D Cylindrical Domain boundaries
    //int matId = 1;
    int arc1 = bc1; // -1;
    int arc2 = bc2; // -2;
    int arc3 = bc3; // -3;
    int arc4 = bc4; // -4;
    
    int nodenumber = 13;
    REAL ModelRadius = 1.5;
    REAL ModelRadiusInt2 = 0.5;//ModelRadius/2.;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordantes for Arc3D 1
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //6
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //7
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //8
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //9
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadiusInt2);//coord Y
    id++;
    //10
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //11
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadiusInt2);//coord Y
    id++;
    //12
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    
    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(3);
    
    // Create Geometrical Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #2
    nodeindex[0] = 0;
    nodeindex[1] = 9;
    nodeindex[2] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #3
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #4
    nodeindex[0] = 1;
    nodeindex[1] = 10;
    nodeindex[2] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #4
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 11;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    // Create Geometrical Quad #4
    nodeindex[0] = 2;
    nodeindex[1] = 11;
    nodeindex[2] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    // Create Geometrical Quad #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    // Create Geometrical Quad #4
    nodeindex[0] = 3;
    nodeindex[1] = 8;
    nodeindex[2] = 11;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #5
    nodeindex.resize(3);
    nodeindex[0] = 8;
    nodeindex[1] = 9;
    nodeindex[2] = 12;
    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    // Create Geometrical Quad #5
    nodeindex[0] = 9;
    nodeindex[1] = 10;
    nodeindex[2] = 12;
    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    // Create Geometrical Quad #5
    nodeindex[0] = 10;
    nodeindex[1] = 11;
    nodeindex[2] = 12;
    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    // Create Geometrical Quad #5
    nodeindex[0] = 11;
    nodeindex[1] = 8;
    nodeindex[2] = 12;
    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Definition of Arc coordenates
    nodeindex.resize(3);
    // Create Geometrical Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2, *gmesh);
    
    
    
    gmesh->BuildConnectivity();
    
    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("CurvoCirculoCTriang.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

// malha que representa um circulo de raio R formado por 5 quadrilateros com quadraticos
TPZGeoMesh *GMeshCirculoQuad(int dimensao, int ndiv)
{
    if(dim != 2)
    {
        DebugStop();
    }

    TPZGeoMesh * gmesh;// = new TPZGeoMesh;

    // description of Geometry and application
    // 2D Cylindrical Domain boundaries
    //int matId = 1;
    int arc1 = bc1; // -1;
    int arc2 = bc2; // -2;
    int arc3 = bc3; // -3;
    int arc4 = bc4; // -4;

    int nodenumber = 20;
    REAL ModelRadius = 1.5;
    REAL ModelRadiusInt2 = 0.5;//ModelRadius/2.;

    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);

    // Setting node coordantes for Arc3D 1
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //6
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //7
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //8
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //9
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadiusInt2);//coord Y
    id++;
    //10
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //11
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadiusInt2);//coord Y
    id++;
    //12
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadiusInt2/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadiusInt2/2.0);//coord Y
    id++;
    //13
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadiusInt2/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadiusInt2/2.0);//coord Y
    id++;
    //14
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadiusInt2/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadiusInt2/2.0);//coord Y
    id++;
    //15
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadiusInt2/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadiusInt2/2.0);//coord Y
    id++;
    //16
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,(ModelRadiusInt2+ModelRadius)/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //17
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,(ModelRadiusInt2+ModelRadius)/2.0);//coord Y
    id++;
    //18
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-(ModelRadiusInt2+ModelRadius)/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //19
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-(ModelRadiusInt2+ModelRadius)/2.0);//coord Y



    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(8);

    // Create Quadratic Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 9;
    nodeindex[3] = 8;
    nodeindex[4] = 4;
    nodeindex[5] = 17;
    nodeindex[6] = 12;
    nodeindex[7] = 16;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, matId,*gmesh);
    //new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;

    // Create Quadratic Quad #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 10;
    nodeindex[3] = 9;
    nodeindex[4] = 5;
    nodeindex[5] = 18;
    nodeindex[6] = 13;
    nodeindex[7] = 17;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, matId,*gmesh);
    //new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;

    // Create Quadratic Quad #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 11;
    nodeindex[3] = 10;
    nodeindex[4] = 6;
    nodeindex[5] = 19;
    nodeindex[6] = 14;
    nodeindex[7] = 18;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, matId,*gmesh);
    //new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;

    // Create Quadratic Quad #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 8;
    nodeindex[3] = 11;
    nodeindex[4] = 7;
    nodeindex[5] = 16;
    nodeindex[6] = 15;
    nodeindex[7] = 19;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, matId,*gmesh);
    //new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;

    // Create Quadratic Quad #5
    nodeindex[0] = 8;
    nodeindex[1] = 9;
    nodeindex[2] = 10;
    nodeindex[3] = 11;
    nodeindex[4] = 12;
    nodeindex[5] = 13;
    nodeindex[6] = 14;
    nodeindex[7] = 15;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, matId,*gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZGeoQuad >  (elementid,nodeindex, matId,*gmesh);
    elementid++;

    // Definition of Arc coordenates
    nodeindex.resize(3);
    
    // Create Quadratic Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc3, *gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1, *gmesh);
    elementid++;

    // Create Quadratic Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc4, *gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2, *gmesh);
    elementid++;

    // Create Quadratic Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc1, *gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3, *gmesh);
    elementid++;

    // Create Quadratic Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc2, *gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4, *gmesh);



    gmesh->BuildConnectivity();

    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }


    std::ofstream out("CinrculoQuadratico.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);

    return gmesh;
}




TPZGeoMesh *GMeshSphericalShell(int dimensao, bool triang, int ndiv)
{
    
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    long materialId = matId;
    long arc1 = bc1; // -1;
    long arc2 = bc2; // -2;
    long arc3 = bc3; // -3;
    long arc4 = bc4; // -4;
    
    int nnodes = 9;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(dim);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 2.;
    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    int Axis = 3;
    REAL theta = 0.0;
    
    int id = 0;
    //no 0
    coord[0] = xc[0] + r;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + 0.;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 1
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + r;
    coord[2] = xc[2] + 0.;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 2
    coord[0] = xc[0] + -r;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + 0.;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 3
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + -r;
    coord[2] = xc[2] + 0.;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 4
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + r;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 5
    coord[0] = xc[0] + r*sqrt(2.)/2.;
    coord[1] = xc[1] + r*sqrt(2.)/2.;
    coord[2] = xc[2] + 0.;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 6
    coord[0] = xc[0] + -r*sqrt(2.)/2.;
    coord[1] = xc[1] + r*sqrt(2.)/2.;
    coord[2] = xc[2] + 0.;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 7
    coord[0] = xc[0] + -r*sqrt(2.)/2.;
    coord[1] = xc[1] + -r*sqrt(2.)/2.;
    coord[2] = xc[2] + 0.;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 8
    coord[0] = xc[0] + r*sqrt(2.)/2.;
    coord[1] = xc[1] + -r*sqrt(2.)/2.;
    coord[2] = xc[2] + 0.;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;

    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<long> topology(3);
    
    // El 0
    topology[0] = 0;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 4;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereEighth1 =
    new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
    SphereEighth1->Geom().SetData(r,xc);
    elementid++;

    // El 1
    topology[0] = 1;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 2;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 4;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereEighth2 =
    new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
    SphereEighth2->Geom().SetData(r,xc);
    elementid++;//
    
    // El 2
    topology[0] = 2;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 3;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 4;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereEighth3 =
    new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
    SphereEighth3->Geom().SetData(r,xc);
    elementid++;//
    
    // El 3
    topology[0] = 3;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 0;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 4;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereEighth4 =
    new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
    SphereEighth4->Geom().SetData(r,xc);
    elementid++;//
    
    // El linha
    // Definition of Arc coordenates
    topology.resize(3);
    // Create Geometrical Arc #1
    topology[0] = 0;
    topology[1] = 1;
    topology[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc3, *geomesh);
    elementid++;

    // Create Geometrical Arc #2
    topology[0] = 1;
    topology[1] = 2;
    topology[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc4, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #3
    topology[0] = 2;
    topology[1] = 3;
    topology[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc1, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #4
    topology[0] = 3;
    topology[1] = 0;
    topology[2] = 8;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc2, *geomesh);
    
    
//    // El linha
//    // Definition of Arc coordenates
//    topology.resize(2);
//    // Create Geometrical Arc #1
//    topology[0] = 0;
//    topology[1] = 1;
//    
//    new TPZGeoElRefPattern < pzgeom::TPZGeoBlend< pzgeom::TPZGeoLinear > > (elementid,topology, arc3, *geomesh);
//    elementid++;
//    
//    // Create Geometrical Arc #2
//    topology[0] = 1;
//    topology[1] = 2;
//    new TPZGeoElRefPattern < pzgeom::TPZGeoBlend< pzgeom::TPZGeoLinear > > (elementid,topology, arc4, *geomesh);
//    elementid++;
//    
//    // Create Geometrical Arc #3
//    topology[0] = 2;
//    topology[1] = 3;
//    new TPZGeoElRefPattern < pzgeom::TPZGeoBlend< pzgeom::TPZGeoLinear > > (elementid,topology, arc1, *geomesh);
//    elementid++;
//    
//    // Create Geometrical Arc #4
//    topology[0] = 3;
//    topology[1] = 0;
//    new TPZGeoElRefPattern < pzgeom::TPZGeoBlend< pzgeom::TPZGeoLinear > > (elementid,topology, arc2, *geomesh);
//    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = 2;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement())
            {
                gel->Divide(sons);
            }
        }
    }
    
    std::ofstream outfile("malhaCascaesfera.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    

}


TPZGeoMesh *GMeshSphericalShellGeob(int dimensao, int ndiv)
{
    
    if(dim != 2)
    {
        DebugStop();
    }
    
    TPZGeoMesh * gmesh;// = new TPZGeoMesh;
    
    // description of Geometry and application
    // 2D Cylindrical Domain boundaries
    //int matId = 1;
    int arc1 = bc1; // -1;
    int arc2 = bc2; // -2;
    int arc3 = bc3; // -3;
    int arc4 = bc4; // -4;
    
    int nodenumber = 13;
    REAL ModelRadius = 2;
    //REAL ModelRadiusInt2 = ModelRadius/2.;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordantes for Arc3D 1
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //6
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //7
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //8
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,sqrt(2.)*ModelRadius/2.);//coord Z
    id++;
    //9
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,sqrt(2.)*ModelRadius/2.);//coord Z
    id++;
    //10
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,sqrt(2.)*ModelRadius/2.);//coord Z
    id++;
    //11
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,sqrt(2.)*ModelRadius/2.);//coord Z
    id++;
    //12
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,ModelRadius);//coord Z
    //id++;
    
    
    
    int elementid = 0;
    TPZVec < long > nodeindex(6,0.0);
    nodeindex.resize(6);
    
    // Create Quadratic Triang #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 12;
    nodeindex[3] = 4;
    nodeindex[4] = 9;
    nodeindex[5] = 8;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticTrig>  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Quadratic Triang #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 12;
    nodeindex[3] = 5;
    nodeindex[4] = 10;
    nodeindex[5] = 9;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticTrig>  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Quadratic Triang #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 12;
    nodeindex[3] = 6;
    nodeindex[4] = 11;
    nodeindex[5] = 10;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticTrig>  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Quadratic Triang #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 12;
    nodeindex[3] = 7;
    nodeindex[4] = 8;
    nodeindex[5] = 11;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticTrig>  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    
    // Definition of Arc coordenates
    nodeindex.resize(3);
    
    // Create Quadratic Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc3, *gmesh);
    elementid++;
    
    // Create Quadratic Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc4, *gmesh);
    elementid++;
    
    // Create Quadratic Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    // Create Quadratic Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc2, *gmesh);
    
    
    
    gmesh->BuildConnectivity();
    
    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("EsferaQuadratica.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
    
}


TPZGeoMesh *GMeshCilindricalMesh( int ndiv)
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    /// Materiais
    long materialId = matId;
    long arc1 = bc1; // -1;
    long arc2 = bc2; // -2;
    long arc3 = bc3; // -3;
    long arc4 = bc4; // -4;
    
    int nodenumber = 6;
    REAL ModelRadius = 1.5;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordinates
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius/2.);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,ModelRadius/2.);//coord Z
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius/2.);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,ModelRadius/2.);//coord Z
    
    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(4);
    
    // Create Geometrical Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 2;
    nodeindex[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    
    // Definition of Arc coordenates
    nodeindex.resize(3);
    // Create Geometrical Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3, *gmesh);
    elementid++;
    
    nodeindex.resize(2);
    
    // Create Geometrical Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc4, *gmesh);
    elementid++;
    
    nodeindex.resize(3);
    
    // Create Geometrical Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    nodeindex.resize(2);
    
    // Create Geometrical Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc2, *gmesh);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    gmesh->BuildConnectivity();
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = 1;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malhaCilindrica.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outfile, true);
    
    return gmesh;

}

TPZGeoMesh *GMeshCilindricalMeshR( int ndiv)
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    /// Materiais
    long materialId = matId;
    long arc1 = bc1; // -1;
    long arc2 = bc2; // -2;
    long arc3 = bc3; // -3;
    long arc4 = bc4; // -4;
    
    int nodenumber = 6;
    REAL ModelRadius = 1.5;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordinates
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,ModelRadius/2.);//coord Z
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,ModelRadius/2.);//coord Z
    
    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(4);
    
    // Create Geometrical Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 2;
    nodeindex[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    
    // Definition of Arc coordenates
    
    
    nodeindex.resize(2);
    
    // Create Geometrical Arc #2
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    nodeindex.resize(3);
    // Create Geometrical Arc #1
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2, *gmesh);
    elementid++;
    
    nodeindex.resize(2);
    // Create Geometrical Arc #4
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc3, *gmesh);

    
    nodeindex.resize(3);
    
    // Create Geometrical Arc #3
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4, *gmesh);
    elementid++;
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    gmesh->BuildConnectivity();
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = 1;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malhaCilindrica.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outfile, true);
    
    return gmesh;
    
}


TPZGeoMesh *GMeshCilindricalMeshF( int ndiv)
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    /// Materiais
    long materialId = matId;
    long arc1 = bc1; // -1;
    long arc2 = bc2; // -2;
    long arc3 = bc3; // -3;
    long arc4 = bc4; // -4;
    
    int nodenumber = 4;
    REAL ModelRadius = 1.5;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordinates
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,0.0);//coord Z

    
    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(4);
    
    // Create Geometrical Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 2;
    nodeindex[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    
    // Definition of Arc coordenates
    nodeindex.resize(2);
    // Create Geometrical Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc3, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    nodeindex.resize(2);
    
    // Create Geometrical Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc4, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc2, *gmesh);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    gmesh->BuildConnectivity();
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = 1;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malhaCilindricaF.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outfile, true);
    
    return gmesh;
    
}


#include "pzstack.h"
void AddWrap(TPZMultiphysicsElement *mfcel, int matskeleton, TPZStack< TPZStack<TPZMultiphysicsElement *,7> > &ListGroupEl)
{
    TPZCompMesh *multiMesh = mfcel->Mesh();
    TPZInterpolationSpace *hdivel = dynamic_cast<TPZInterpolationSpace *> (mfcel->Element(0));
    TPZCompElDisc *discel = dynamic_cast<TPZCompElDisc *>(mfcel->Element(1));
    TPZGeoEl *gel = mfcel->Reference();
    
    int dimMesh = mfcel->Mesh()->Dimension();
    if (!hdivel || !discel || gel->Dimension() != dimMesh) {
        DebugStop();
    }
    
    //wrap element
    TPZStack<TPZMultiphysicsElement *, 7> wrapEl;
    wrapEl.push_back(mfcel);
    
    for (int side = 0; side < gel->NSides(); side++)
    {
        if (gel->SideDimension(side) != gel->Dimension()-1) {
            continue;
        }
        TPZGeoEl *gelbound = gel->CreateBCGeoEl(side, matskeleton);
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(hdivel);
        int loccon = intel->SideConnectLocId(0,side);
        long index;
        
        TPZInterpolationSpace *bound;
        MElementType elType = gel->Type(side);
        switch(elType)
        {
            case(EOned)://line
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(bound);
                hdivbound->SetSideOrient(sideorient);
                break;
            }
            case(ETriangle)://triangle
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeTriang>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeTriang> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeTriang> *>(bound);
                hdivbound->SetSideOrient(sideorient);
                break;
            }
            case(EQuadrilateral)://quadrilateral
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeQuad>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(bound);
                hdivbound->SetSideOrient(sideorient);
                break;
            }
                
            default:
            {
                bound=0;
                std::cout << "ElementType not found!";
                DebugStop();
                break;
            }
        }
        
        long sideconnectindex = intel->ConnectIndex(loccon);
        bound->SetConnectIndex(0, sideconnectindex);
        //bound->Print(std::cout);
        
        TPZCompEl *newMFBound = multiMesh->CreateCompEl(gelbound, index);
        TPZMultiphysicsElement *locMF = dynamic_cast<TPZMultiphysicsElement *>(newMFBound);
        
        locMF->AddElement(bound, 0);
        locMF->AddElement(TPZCompElSide(discel,side), 1);
        
        wrapEl.push_back(locMF);
    }
    
    ListGroupEl.push_back(wrapEl);
}


void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    
    switch (Axis) {
        case 1:
        {
            RotationMatrix(0,0) = 1.0;
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(1,2) =   -sin(theta);
            RotationMatrix(2,1) =   +sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 2:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,2) =   +sin(theta);
            RotationMatrix(1,1) = 1.0;
            RotationMatrix(2,0) =   -sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 3:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
        default:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
    }
    
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    //RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

void RotateNode(TPZVec<STATE> &iCoords, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    
    switch (Axis) {
        case 1:
        {
            RotationMatrix(0,0) = 1.0;
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(1,2) =   -sin(theta);
            RotationMatrix(2,1) =   +sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 2:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,2) =   +sin(theta);
            RotationMatrix(1,1) = 1.0;
            RotationMatrix(2,0) =   -sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 3:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
        default:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
    }
    
    TPZVec<STATE> iCoordsRotated(3,0.0);
    // Apply rotation
    iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
    iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
    iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
    iCoords = iCoordsRotated;
}

void Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0];
    X[1] = 0.0;
    X[2] = 0.0;
}

void Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void Parametricfunction3(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}


