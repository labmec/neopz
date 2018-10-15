
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgmesh.h"
#include "pzcheckgeom.h"
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
#include "pzvec_extras.h"
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

void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);
void RotateNode(TPZVec<REAL> &iCoords, REAL CounterClockwiseAngle, int &Axis);


TPZGeoMesh *GMeshSphericalShell(int dimensao, bool triang, int ndiv);
TPZGeoMesh *GMeshSphericalShell2(int dimensao, bool triang, int ndiv);
TPZGeoMesh *GMeshSliceSphericalShell(int dimensao, bool triang, int ndiv);
TPZGeoMesh *GMeshSphericalRing(int dimensao, bool triang,int ndiv);
TPZGeoMesh *GMeshSphericalRingQuarter(int dimensao, bool triang,int ndiv);
TPZVec<REAL> SphereToKartesian(REAL r, REAL theta, REAL phi);
TPZGeoMesh *GMeshTropicodeCancer(int ndiv , TPZVec<bool>  &CurvesSides, bool isPlane, int plane);
void PrintLS(TPZAnalysis *an);


TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);

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

int dim = 2;
REAL aa = 0.0;
REAL bb = 0.0;
REAL cc = 0.0;
REAL Epsilon = 0.4;
// tensor de permutacao
TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);

REAL const Pi = M_PI;//4.*atan(1.);

bool isH1 = false;

bool isgeoblend = false;

bool ftriang = false;

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.douglas"));
#endif

#define RING

#include "pztransfer.h"

int main(int argc, char *argv[])
{
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    std::string FileName = dirname;
    FileName = dirname + "/Projects/dactests/";
    FileName += "Douglas.cfg";
    //InitializePZLOG(FileName);
    InitializePZLOG();
#endif


    int p = 1;
    int ndiv = 0;

    // Hard coded
    int sizeTP = 3;
    
    TP.Resize(sizeTP, sizeTP);
    InvTP.Resize(sizeTP, sizeTP);
    
    for (int id = 0; id < sizeTP; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    
    for(p=4;p<5;p++)
    {
        int pq = p;
        int pp = p;
        
        for (ndiv=2; ndiv<3; ndiv++)
        {
            std::cout<< " INICIO - grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
            std::cout<< " Dimensao == " << dim << std::endl;
#ifdef RING
            //TPZGeoMesh *gmesh = GMeshSphericalRing(2, ftriang, ndiv);
            //TPZGeoMesh *gmesh = GMeshSphericalRingQuarter(2, ftriang, ndiv);
            //TPZGeoMesh *gmesh = GMeshSliceSphericalShell(2, ftriang, ndiv);
            
            // caso horizontal arestas retas
            TPZVec<bool> CurvesSides(4,true);
//            CurvesSides[1] = false;
//            CurvesSides[3] = false;
            TPZGeoMesh *gmesh = GMeshTropicodeCancer(ndiv,CurvesSides,false,3);
            
            // caso inclinado arestas retas
//            TPZVec<bool> CurvesSides(4,false);
//            CurvesSides[1] = false;
//            CurvesSides[3] = false;
//            TPZGeoMesh *gmesh = GMeshTropicodeCancer(ndiv,CurvesSides,false);
            
            // arestas tortas
//            TPZVec<bool> CurvesSides(4,true);
//            TPZGeoMesh *gmesh = GMeshTropicodeCancer(ndiv,CurvesSides,false);
            
            
            //TPZGeoMesh *gmesh = GMeshTropicodeCancer(ndiv,CurvesSides,false);
#else
            
            //TPZGeoMesh *gmesh = GMeshSphericalShell(2, ftriang, ndiv);
            //TPZGeoMesh *gmesh = GMeshSphericalShell2(2, ftriang, ndiv);
            TPZGeoMesh *gmesh = GMeshSliceSphericalShell(2, ftriang, ndiv);
#endif
           
            {
                gmesh->SetDimension(dim);
                ofstream arg1("gmesh2d.txt");
                gmesh->Print(arg1);
            }

            
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
            
            
//            ofstream argf("cmeshflux.txt");
//            cmesh1->Print(argf);
//
//            ofstream argp("cmeshpressure.txt");
//            cmesh2->Print(argp);
            
            
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
    
            
            TPZAnalysis an(mphysics, false);

            SolveSyst(an, mphysics);
            
//            PrintLS(&an);
            
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

            std::cout << "Postprocessed\n";
            {
                ofstream argp("cmeshpressure.txt");
                cmesh2->Print(argp);
                
                ofstream arg5("cmeshmultiphysics.txt");
                mphysics->Print(arg5);
            }
            
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


void PrintLS(TPZAnalysis *an)
{
    an->Assemble();
    TPZAutoPointer< TPZMatrix<STATE> > Kglob;
    TPZFMatrix<STATE> Fglob;
    Kglob=an->Solver().Matrix();
    Fglob = an->Rhs();
    
#ifdef LOG4CXX
    if(logdata->isDebugEnabled())
    {
        std::stringstream sout;
        Kglob->Print("matK = ", sout,EMathematicaInput);
        Fglob.Print("fvec = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logdata,sout.str())
    }
#endif
    
}

TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
//    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
//    TPZAutoPointer<TPZFunction<STATE> > forcef;
//    forcef = new TPZDummyFunction<STATE>(Forcing, 5);
//    material->SetForcingFunction(forcef);
    
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingH1, 5);
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
    
    
    if( dim == 3 ) {
        cmesh->InsertMaterialObject(BCond0);
    }
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    if(dim == 3 )
    {
        cmesh->InsertMaterialObject(BCond5);
    }
    
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
    
#ifdef PZDEBUG
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
    TPZFNMatrix<9,REAL> PermTensor(3,3,0.);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3,0.);
    
    PermTensor = TP;
    InvPermTensor = InvTP;
    
    
    material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    
    solexata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing, 5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(0);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    if (!hdivantigo)
    {
        material->UseSecondIntegrationByParts();
    }
    
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
        TPZAutoPointer<TPZFunction<STATE> > FBCond0 = new TPZDummyFunction<STATE>(ForcingBC0D, 5);
        BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
        BCond0->SetForcingFunction(FBCond0);
    }

    
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1D, 5);
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
//    BCond1 = material->CreateBC(mat, bc1,neumann, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2N, 5);
    BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    BCond2->SetForcingFunction(FBCond2);
//    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2N, 5);
//    BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
//    BCond2->SetForcingFunction(FBCond2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3D, 5);
//    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    BCond3->SetForcingFunction(FBCond3);
//
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4N, 5);
    BCond4 = material->CreateBC(mat, bc4,neumann, val1, val2);
    BCond4->SetForcingFunction(FBCond4);
//    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4N, 5);
//    BCond4 = material->CreateBC(mat, bc4,neumann, val1, val2);
//    BCond4->SetForcingFunction(FBCond4);
    
    
    
    if (dim==3)
    {
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZAutoPointer<TPZFunction<STATE> > FBCond5 = new TPZDummyFunction<STATE>(ForcingBC5D, 5);
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
            if(mfcel->Dimension()==dim) {
                TPZBuildMultiphysicsMesh::AddWrap(mfcel, matId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
            }
        }
        meshvec[0]->CleanUpUnconnectedNodes();
        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
        //------- Create and add group elements -------
        int64_t index, nenvel;
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
    TPZManVector<std::string,10> scalnames(3), vecnames(3);
    vecnames[0]  = "Flux";
    vecnames[1]  = "ExactFlux";
    vecnames[2]  = "GradP";
    //    vecnames[1]  = "GradFluxX";
    //    vecnames[2]  = "GradFluxY";
    //    vecnames[2]  = "GradFluxz";
    scalnames[0] = "Pressure";
    scalnames[1] = "ExactPressure";
    scalnames[2] = "Rhs";
    
    //const int dim = 3;
    int div = 4;
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
    REAL a=5.0*M_PI/16.0;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    // Hard coded

    flux.Resize(3, 1);

    
    
    REAL r = sqrt(x*x+y*y+z*z);
    //REAL theta = (atan2(sqrt(x*x+y*y),z));
    REAL theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    REAL phi = atan2(y,x);
    
    
    

#ifdef RING
    // Anel
    solp[0] = (a-theta);
    flux(0,0)= (cos(phi)*cos(theta))/r;
    flux(1,0)= (cos(theta)*sin(phi))/r;
    flux(2,0)= -(sin(theta))/r;
    
//    // teste 1
//    solp[0] = 0.5*(x*x);
//    flux(0,0) = -1.0;
//    flux(1,0) = -1.0;
//    flux(2,0) = 0.0;
    
#else
    
    // problema na casca
    
//    solp[0] = (a-theta)*sin(theta)*sin(theta);
//    flux(0,0)= (cos(phi)*cos(theta)*(2.0*(a - theta)*cos(theta) - sin(theta))*sin(theta))/r;
//    flux(1,0)= (cos(theta)*sin(phi)*(2.0*(a - theta)*cos(theta) - sin(theta))*sin(theta))/r;
//    flux(2,0)= (sin(theta)*sin(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta)))/r;
    
    /* Problem with phi */
//    solp[0] = (a-theta)*phi*phi*sin(theta)*sin(theta);
//    flux(0,0)= (phi*sin(theta)*(2.0*(-a + theta)*sin(phi)-phi*cos(phi)*cos(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta))))/r;
//    flux(1,0)= (phi*sin(theta)*(2.0*(a - theta)*(cos(phi) + phi*cos(theta)*cos(theta)*sin(phi)) -phi*cos(theta)*sin(phi)*sin(theta)))/r;
//    flux(2,0)= (phi*phi*sin(theta)*sin(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta)))/r;
    
    
    
#endif
    
    // Problema na circunferencia
//    r = sqrt(x*x+y*y);
//    if (r<1) {
//        solp[0] = 3.0*exp(1.0/(x*x+y*y-1.0));
//        flux(0,0) = 6.0*x*exp(1.0/(x*x+y*y-1.0))/((x*x+y*y-1.0)*(x*x+y*y-1.0));
//        flux(1,0) = 6.0*y*exp(1.0/(x*x+y*y-1.0))/((x*x+y*y-1.0)*(x*x+y*y-1.0));
//        //flux(2,0) = 0.0;
//    }
//    else
//    {
//        solp[0] = 0.0;
//        flux(0,0) = 0.0;
//        flux(1,0) = 0.0;
//      //  flux(2,0) = 0.0;
//    }
//
//    
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){

    
    REAL x,y,z;
    // Defining data
//    REAL a=M_PI/2.0;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
   
    REAL r = sqrt(x*x+y*y+z*z);
    //REAL theta = (atan2(sqrt(x*x+y*y),z));
    REAL theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
//    REAL phi = atan2(y,x);
    
#ifdef RING
    // anel
    ff[0] = (1.0/(r*r))*(1.0/tan(theta));
    ff[0] = 0.0;
    
#else
    // casca
    ff[0] = (-1.0/(2.0*r*r))*(2.0*(a-theta)*(1.0+3.0*cos(2.0*theta))- 5.0*sin(2.0*theta));
//
//    std::cout << " x =" << x << "y =" << y << "z =" << z << std::endl;
//    std::cout << " Theta =" << theta << "phi =" << phi << "r =" << r << std::endl;
//    std::cout << " ff[0] =" << ff[0] << "a =" << a << std::endl;
    
//    /* Problem with phi */
//      ff[0] = (-1.0/(r*r))*((a-theta)*(2.0-phi*phi+3.0*phi*phi*cos(2.0*theta))- 5.0*phi*phi*cos(theta)*sin(theta));
    
#endif
    
    // Problema na circunferencia
//    r = sqrt(x*x+y*y);
//    if (r<1) {
//        ff[0] = -12.0*(-1.0+x*x*x*x+y*y+y*y*y*y+x*x*(1.0+2.0*y*y))*exp(1.0/(x*x+y*y-1.0))/((x*x+y*y-1.0)*(x*x+y*y-1.0)*(x*x+y*y-1.0)*(x*x+y*y-1.0));
//    }
//    else
//    {
//        ff[0] = 0.0;
//    }


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
    
    REAL x,y,z;
    // Defining data
    REAL a=M_PI/2.0;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    
    
#ifdef RING
    // Anel
    solp[0] = (a-theta);
    solp[0] = -4.0;
    
#else
    
    // problema na casca
    REAL t = (a-theta)*sin(theta)*sin(theta);
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
    
    /* Problem with phi */
    //    solp[0] = (a-theta)*phi*phi*sin(theta)*sin(theta);
    
#endif
}


void ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x,y,z;
    // Defining data
    REAL a=M_PI/2.0;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    
#ifdef RING
    // Anel
    solp[0] = (a-theta);
    
#else
    
    // problema na casca
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
    
    /* Problem with phi */
    //    solp[0] = (a-theta)*phi*phi*sin(theta)*sin(theta);
    
#endif
    
}


void ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x,y,z;
    // Defining data
    REAL a=M_PI/2.0;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    
    
#ifdef RING
    // Anel
    solp[0] = (a-theta);
    solp[0] = 4.0;
    
#else
    
    // problema na casca
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
    
    /* Problem with phi */
    //    solp[0] = (a-theta)*phi*phi*sin(theta)*sin(theta);
    
#endif
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
    
    disp[0] = 0.0;
}

void ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    disp[0] = -1.0;
}

void ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    disp[0] = 0.0;
}

void ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
}

void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, false);
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
}

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv)
{
    int64_t nel = l2mesh->NElements();
    //int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, false);
        int nerr = elerror.size();
        globalerrors.resize(nerr);
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


TPZGeoMesh *GMeshSphericalRingQuarter(int dimensao, bool triang, int ndiv)
{
    
    bool ftriangulo = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int materialId = matId;
    int64_t arc1 = bc1; // -1;
    int64_t arc2 = bc2; // -2;
    int64_t arc3 = bc3; // -3;
    int64_t arc4 = bc4; // -4;
    
    int nnodes = 4;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(dim);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 1.;
    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    REAL theta = Pi/2.0 ;
    REAL phi = Pi/6.0;
    
    int Axis = 3;
    REAL angulo = 0.0;
    
    int id = 0;
    //no 0
    coord[0] = xc[0] + r*sin(theta + Pi/6.0)*cos(-phi);
    coord[1] = xc[1] + r*sin(theta + Pi/6.0)*sin(-phi);
    coord[2] = xc[2] + r*cos(theta + Pi/6.0);
    
    RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 1
    coord[0] = xc[0] + r*sin(theta + Pi/6.0)*cos(phi);
    coord[1] = xc[1] + r*sin(theta + Pi/6.0)*sin(phi);
    coord[2] = xc[2] + r*cos(theta + Pi/6.0);
    RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 2
    coord[0] = xc[0] + r*sin(theta - Pi/6.0)*cos(phi);
    coord[1] = xc[1] + r*sin(theta - Pi/6.0)*sin(phi);
    coord[2] = xc[2] + r*cos(theta - Pi/6.0);
    RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    
    //no 3
    coord[0] = xc[0] + r*sin(theta - Pi/6.0)*cos(-phi);
    coord[1] = xc[1] + r*sin(theta - Pi/6.0)*sin(-phi);
    coord[2] = xc[2] + r*cos(theta - Pi/6.0);
    RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    //id++;
    
    int64_t elementid = 0;

    // Using triangle to sphere special map
    TPZVec<int64_t> topology;
    
    if (ftriangulo)
    {
        
        
        topology.Resize(3);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 2;
        
        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereEighth1 =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
        SphereEighth1->Geom().SetData(r,xc);
        elementid++;
        
        // El 0
        topology[0] = 0;
        topology[1] = 2;
        topology[2] = 3;
        

        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereEighth2 =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
        SphereEighth2->Geom().SetData(r,xc);
        //elementid++;

        
    }
    else
    {
        
        topology.Resize(4);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 2;
        topology[3] = 3;
        
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereEighth1 =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
        SphereEighth1->Geom().SetData(r,xc);
        elementid++;
        
    }
    
    // El linha
    // Definition of Arc coordenates
    topology.resize(2);
    // Create Geometrical Arc #1
    topology[0] = 0;
    topology[1] = 1;
    {
        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc3, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #2
    topology[0] = 1;
    topology[1] = 2;
    {
        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc2, *geomesh);
        // arc->Geom().Print(std::cout);
    }
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc4, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #3
    topology[0] = 2;
    topology[1] = 3;
    {
        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc3, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc1, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #4
    topology[0] = 3;
    topology[1] = 0;
    {
        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc4, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc2, *geomesh);
    //    elementid++;
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    

    const unsigned int nel = 1;//geomesh->NElements();
    for (int i = 0 ; i < nel ; i++){
        TPZGeoEl *gel = geomesh->Element(i);
        if (!gel) {
            continue;
        }
        const int npt = 5;
        TPZManVector<REAL,3> qsi(3,0.), x(3,0.);
        for (int iqsi = 0; iqsi < npt; iqsi++) {
            for (int ieta = 0; ieta < npt; ieta++) {
                qsi[0] = -1 + iqsi * 2./(npt-1);
                qsi[1] = -1 + ieta * 2./(npt-1);
                gel->X(qsi, x);
                //calcula a funcao e ve se zero
                cout << " qsi = " << qsi << " x = " << x << endl;
                cout << " r = " << sqrt(x[0]*x[0]+x[1]*x[1]) << endl;
                cout << " R = " << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << endl;
                //cout << " phi
            }
        }
    }
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
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
    
    std::ofstream outfile("malhaCascaesferaQuarto.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
    
}


TPZGeoMesh *GMeshTropicodeCancer(int ndiv , TPZVec<bool>  &CurvesSides, bool isPlane, int plane)
{
    
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int materialId = matId;
    int64_t arc1 = bc1; // -1;
    int64_t arc2 = bc2; // -2;
    int64_t arc3 = bc3; // -3;
    int64_t arc4 = bc4; // -4;
    
    int nnodes = 8;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(dim);
    
    if (CurvesSides.size() != 4) {
        DebugStop();
    }
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
//    REAL beta = M_PI/16.0;
//    REAL alpha = M_PI/4.0;
//    REAL h = M_PI/4.0;
//    REAL k = M_PI/4.0;
    
    REAL beta = M_PI/3.0;
    REAL alpha = M_PI/2.0;
    REAL h = M_PI/2.0;
    REAL k = 0.0;
    
    int id=0;
    TPZVec<REAL> coord(3,0.);

    if (isPlane) {
        
        switch (plane) {
            case 1:
            {
                //no 0
                coord = SphereToKartesian(1.0, h-beta, k-alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 1
                coord = SphereToKartesian(1.0, h-beta, k+alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 2
                coord = SphereToKartesian(1.0, h-beta, k+alpha + 2.0*(h-alpha));
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 3
                coord = SphereToKartesian(1.0, h-beta, k-alpha - 2.0*(h-alpha));
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
            }
                break;
                
                case 2:
            {
                //no 0
                coord = SphereToKartesian(1.0, h-beta, k-alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 1
                coord = SphereToKartesian(1.0, h-beta, k+alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 2
                coord = SphereToKartesian(1.0, h+beta, k+alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 3
                coord = SphereToKartesian(1.0, h+beta, k-alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
            }
                break;
                
            case 3:
            {
                //no 0
                coord = SphereToKartesian(1.0, h-beta, k+alpha + 2.0*(h-alpha));
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 1
                coord = SphereToKartesian(1.0, h-beta, k + alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 2
                coord = SphereToKartesian(1.0, h+beta, k + alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 3
                coord = SphereToKartesian(1.0, h+beta, k + alpha + 2.0*(h-alpha));
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
            }
                break;
                
            case 4:
            {
                // domain with equal lentgh
                //no 0
                coord[0] = -M_PI/8.0;
                coord[1] = -sqrt(3.0)/2.0;
                coord[2] = 0.0;
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;

                //no 1
                coord[0] = M_PI/8.0;
                coord[1] = -sqrt(3.0)/2.0;
                coord[2] = 0.0;
                
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                
                //no 2
                coord[0] = M_PI/8.0;
                coord[1] = sqrt(3.0)/2.0;
                coord[2] = 0.0;
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;

                //no 3
                coord[0] = -M_PI/8.0;
                coord[1] = sqrt(3.0)/2.0;
                coord[2] = 0.0;
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;

            }
                break;
                
            default:
            {
                DebugStop();
            }
                break;
        }
        
    }
    else{
        //no 0
        coord = SphereToKartesian(1.0, h + beta , k  ); //k + M_PI/22.0
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 1
        coord = SphereToKartesian(1.0, h + beta , k + alpha  );//k + alpha - M_PI/22.0
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 2
        coord = SphereToKartesian(1.0, h - beta , k + alpha );
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 3
        coord = SphereToKartesian(1.0, h - beta , k);
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 4
        coord = SphereToKartesian(1.0, h + beta , k + alpha/2.0);
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 5
        coord = SphereToKartesian(1.0, h, k + alpha);
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 6
        coord = SphereToKartesian(1.0, h - beta , k + alpha/2.0);
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 7
        coord = SphereToKartesian(1.0, h, k);
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
    }
    
    geomesh->SetMaxNodeId(id);
    int64_t elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(3);
    TPZVec<int64_t> topologyLine(2);
    

    // Side 0
    if (CurvesSides[0]) {
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 4;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc1, *geomesh);
        elementid++;
    }
    else{
        topologyLine[0] = 0;
        topologyLine[1] = 1;
        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,topologyLine, arc1, *geomesh);
        elementid++;
    }
    
    // Side 1
    if (CurvesSides[1]) {
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 5;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc2, *geomesh);
        elementid++;
    }
    else{
        topologyLine[0] = 1;
        topologyLine[1] = 2;
        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,topologyLine, arc2, *geomesh);
        elementid++;
    }
    
    // Side 2
    if (CurvesSides[2]) {
        topology[0] = 2;
        topology[1] = 3;
        topology[2] = 6;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc3, *geomesh);
        elementid++;
    }
    else{
        topologyLine[0] = 2;
        topologyLine[1] = 3;
        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,topologyLine, arc3, *geomesh);
        elementid++;
    }
    
    // Side 3
    if (CurvesSides[3]) {
        topology[0] = 3;
        topology[1] = 0;
        topology[2] = 7;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc4, *geomesh);
        elementid++;
    }
    else{
        topologyLine[0] = 3;
        topologyLine[1] = 0;
        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,topologyLine, arc4, *geomesh);
        elementid++;
    }
    
    REAL r = 1.0;
    TPZManVector<REAL,3> xc(3,0.0);
    
    // Create Geometrical Quad #1
    topology.Resize(4);    
    topology[0] = 0;
    topology[1] = 1;
    topology[2] = 2;
    topology[3] = 3;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere<pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > > * SphereEighth1 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > > (elementid, topology,materialId,*geomesh);
    SphereEighth1->Geom().SetData(r,xc);
    elementid++;
    geomesh->SetMaxElementId(elementid);
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    

    
    
//    const unsigned int nel = geomesh->NElements();
//    for (int i = 0 ; i < nel ; i++){
//        TPZGeoEl *gel = geomesh->Element(i);
//        if (!gel) {
//            continue;
//        }
//        const int npt = 3;
//        TPZManVector<REAL,3> qsi(3,0.), x(3,0.);
//        for (int iqsi = 0; iqsi < npt; iqsi++) {
//            for (int ieta = 0; ieta < npt; ieta++) {
//                qsi[0] = -1 + iqsi * 2./(npt-1);
//                qsi[1] = -1 + ieta * 2./(npt-1);
//                gel->X(qsi, x);
//                //calcula a funcao e ve se zero
//                cout << " qsi = " << qsi << " x = " << x << endl;
//                cout << " r = " << sqrt(x[0]*x[0]+x[1]*x[1]) << endl;
//                cout << " R = " << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << endl;
//                
//            }
//        }
//    }
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
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
    
    std::ofstream outfileAnel("TropicodeCancer.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfileAnel, true);
    
#ifdef PZDEBUG
    TPZCheckGeom check(geomesh);
    check.CheckUniqueId();
#endif
    return geomesh;
    

}

TPZVec<REAL> SphereToKartesian(REAL r, REAL theta, REAL phi)
{
    TPZVec<REAL> xyz(3,0.0);
    xyz[0] = r*cos(phi)*sin(theta);
    xyz[1] = r*sin(theta)*sin(phi);
    xyz[2] = r*cos(theta);
    return xyz;
}


TPZGeoMesh *GMeshSphericalShell(int dimensao, bool triang, int ndiv)
{
    
    
    bool ftriangulo = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = matId;
    int64_t arc1 = bc1; // -1;
    
    int nnodes = 9;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(dim);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 1.;
    REAL z = r/2.;
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
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + z;
    coord[0] = xc[0] + sqrt(r*r - coord[2]*coord[2] - coord[1]*coord[1]);
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 5
    coord[0] = xc[0] + 0.;
    coord[2] = xc[2] + z;
    coord[1] = xc[1] + sqrt(r*r - coord[0]*coord[0] - coord[2]*coord[2]);
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 6
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + z;
    coord[0] = xc[0] + -sqrt(r*r - coord[2]*coord[2] - coord[1]*coord[1]);
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 7
    coord[0] = xc[0] + 0.;
    coord[2] = xc[2] + z;
    coord[1] = xc[1] + -sqrt(r*r - coord[0]*coord[0] - coord[2]*coord[2]);
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 8
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + r;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    //id++;
    
    int64_t elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(4);
    
    
    if (ftriangulo)
    {
        
        topology.resize(3);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 0;
        topology[1] = 5;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <>  > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 1;
        topology[1] = 6;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 4
        topology[0] = 2;
        topology[1] = 3;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 5
        topology[0] = 2;
        topology[1] = 7;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 6
        topology[0] = 3;
        topology[1] = 0;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 7
        topology[0] = 3;
        topology[1] = 4;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
    }
    else
    {
        
        topology.resize(4);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 5;
        topology[3] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingQ->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 6;
        topology[3] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingQ->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 2;
        topology[1] = 3;
        topology[2] = 7;
        topology[3] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingQ->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 3;
        topology[1] = 0;
        topology[2] = 4;
        topology[3] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingQ->Geom().SetData(r,xc);
            elementid++;
        }
        
    }
    
    topology.resize(3);
    // El 4
    topology[0] = 4;
    topology[1] = 5;
    topology[2] = 8;
    {
        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereTriangQ =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
        SphereTriangQ->Geom().SetData(r,xc);
        elementid++;
    }
    // El 5
    topology[0] = 5;
    topology[1] = 6;
    topology[2] = 8;
    {
        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereTriangQ =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
        SphereTriangQ->Geom().SetData(r,xc);
        elementid++;
    }
    // El 6
    topology[0] = 6;
    topology[1] = 7;
    topology[2] = 8;
    {
        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereTriangQ =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
        SphereTriangQ->Geom().SetData(r,xc);
        elementid++;
    }
    // El 4
    topology[0] = 7;
    topology[1] = 4;
    topology[2] = 8;
    {
        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereTriangQ =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
        SphereTriangQ->Geom().SetData(r,xc);
        elementid++;
    }
    
    // El linha
    // Definition of Arc coordenates
    topology.resize(2);
    // Create Geometrical Arc #0
    topology[0] = 0;
    topology[1] = 1;
    {
        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc3 */, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    elementid++;
    
    // Create Geometrical Arc #1
    topology[0] = 1;
    topology[1] = 2;
    {
        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc4 */, *geomesh);
        // arc->Geom().Print(std::cout);
    }
    elementid++;
    
    // Create Geometrical Arc #2
    topology[0] = 2;
    topology[1] = 3;
    {
        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc1 */, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    elementid++;
    
    // Create Geometrical Arc #3
    topology[0] = 3;
    topology[1] = 0;
    {
        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc2 */, *geomesh);
        //arc->Geom().Print(std::cout);
    }

    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
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



TPZGeoMesh *GMeshSphericalShell2(int dimensao, bool triang, int ndiv)
{
    
    
    bool ftriangulo = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int materialId = matId;
    int64_t arc1 = bc1; // -1;
    
    int nnodes = 37;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(dim);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 1.;

    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    int Axis = 3;
    REAL theta = 0.0;
    
    int sind = (nnodes+1)/3;
    
    TPZManVector<int> indinf(sind);
    TPZManVector<int> indmid(sind);
    TPZManVector<int> indsup(sind);
    int cont = 0;
    
    int id = 0;
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 12; i++) {
            //no id
            coord[0] = xc[0] + r*cos(i*Pi/6.0)*sin(Pi/2.0 - j*Pi/6.0);
            coord[1] = xc[1] + r*sin(i*Pi/6.0)*sin(Pi/2.0 - j*Pi/6.0);
            coord[2] = xc[2] + r*cos(Pi/2.0 - j*Pi/6.0);
            RotateNode(coord, theta, Axis);
            node.SetNodeId(id);
            node.SetCoord(coord);
            geomesh->NodeVec()[id] = node;
            if (j==0) {
                indinf[cont] = id;
                cont++;
            }
            else if (j==1)
            {
                indmid[cont] = id;
                cont++;
            }
            else
            {
                indsup[cont] = id;
                cont++;
            }
            id++;
            
        }
        cont = 0;
    }
    
    //no polo norte
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + r;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    //id++;
    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(4);
    
    
    if (ftriangulo)
    {
        
        topology.resize(3);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 0;
        topology[1] = 5;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 1;
        topology[1] = 6;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 4
        topology[0] = 2;
        topology[1] = 3;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 5
        topology[0] = 2;
        topology[1] = 7;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 6
        topology[0] = 3;
        topology[1] = 0;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 7
        topology[0] = 3;
        topology[1] = 4;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
    }
    else
    {
        
        topology.resize(4);
        
        for (int nelinf = 0; nelinf < 12; nelinf++) {
            // El nel
            
            topology[0] = indinf[nelinf];
            topology[1] = indinf[(nelinf+1)%12];
            topology[2] = indmid[(nelinf+1)%12];
            topology[3] = indmid[nelinf];
            {
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
                new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
                SphereRingQ->Geom().SetData(r,xc);
                elementid++;
            }
        }
        
        for (int nelmid = 0; nelmid < 12; nelmid++) {
            // El nel
            topology[0] = indmid[nelmid];
            topology[1] = indmid[(nelmid+1)%12];
            topology[2] = indsup[(nelmid+1)%12];
            topology[3] = indsup[nelmid];
            {
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
                new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
                SphereRingQ->Geom().SetData(r,xc);
                elementid++;
            }
        }
        
        topology.resize(3);
        
//        for (int nelsup = 0; nelsup < 12; nelsup++) {
//            // El nel
////            int a = indsup[nelsup];
////            int b = indsup[(nelsup+1)%12];
////            int c = polo;
//            topology[0] = indsup[nelsup];
//            topology[1] = indsup[(nelsup+1)%12];
//            topology[2] = polo;
//            {
//                TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereRingQ =
//                new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
//                SphereRingQ->Geom().SetData(r,xc);
//                elementid++;
//            }
//        }
        
        
    }
    
    // El linha
    // Definition of Arc
    topology.resize(2);
    
    for (int nelinf = 0; nelinf < 12; nelinf++) {
        // El nel
        topology[0] = indinf[nelinf];
        topology[1] = indinf[(nelinf+1)%12];
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    // arcos do topo
    for (int nelsup = 0; nelsup < 12; nelsup++) {
        // El nel
        topology[0] = indsup[nelsup];
        topology[1] = indsup[(nelsup+1)%12];
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    
    
//    // Create Geometrical Arc #0
//    topology[0] = 0;
//    topology[1] = 1;
//    {
//        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc3 */, *geomesh);
//        //arc->Geom().Print(std::cout);
//    }
//    elementid++;
//    
//    // Create Geometrical Arc #1
//    topology[0] = 1;
//    topology[1] = 2;
//    {
//        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc4 */, *geomesh);
//        // arc->Geom().Print(std::cout);
//    }
//    elementid++;
//    
//    // Create Geometrical Arc #2
//    topology[0] = 2;
//    topology[1] = 3;
//    {
//        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc1 */, *geomesh);
//        //arc->Geom().Print(std::cout);
//    }
//    elementid++;
//    
//    // Create Geometrical Arc #3
//    topology[0] = 3;
//    topology[1] = 0;
//    {
//        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc2 */, *geomesh);
//        //arc->Geom().Print(std::cout);
//    }
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
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

TPZGeoMesh *GMeshSliceSphericalShell(int dimensao, bool triang, int ndiv)
{
    
    int nfatias = 1;
    const REAL r = 1.;
    // centro da esfera
    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    // para rotacao da malha, se quiser
    int Axis = 3;
    REAL theta = 0.0;
    
    
    bool ftriangulo = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int materialId = matId;
    int64_t arc1 = bc1; // -1;
    int64_t arc3 = bc3; // -3;
    
    int nnodes = nfatias == 12 ? 37 :(nfatias+1)*3+1;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(dim);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    
    // quantidade de indices
    int sind = (nnodes+1)/3;
    
    TPZManVector<int> indinf(sind);
    TPZManVector<int> indmid(sind);
    TPZManVector<int> indsup(sind);
    int cont = 0;
    int npts = nfatias == 12 ? 12 : nfatias + 1;
    
    int id = 0;
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < npts; i++) {
            //no id
            coord[0] = xc[0] + r*cos(i*Pi/6.0)*sin(Pi/2.0 - j*Pi/6.0);
            coord[1] = xc[1] + r*sin(i*Pi/6.0)*sin(Pi/2.0 - j*Pi/6.0);
            coord[2] = xc[2] + r*cos(Pi/2.0 - j*Pi/6.0);
            RotateNode(coord, theta, Axis);
            node.SetNodeId(id);
            node.SetCoord(coord);
            geomesh->NodeVec()[id] = node;
            if (j==0) {
                indinf[cont] = id;
                cont++;
            }
            else if (j==1)
            {
                indmid[cont] = id;
                cont++;
            }
            else
            {
                indsup[cont] = id;
                cont++;
            }
            id++;
            
        }
        cont = 0;
    }
    
    //no polo norte
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + r;
    RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(4);
    
    
    if (ftriangulo)
    {
        DebugStop();
        topology.resize(3);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 0;
        topology[1] = 5;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 1;
        topology[1] = 6;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 4
        topology[0] = 2;
        topology[1] = 3;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 5
        topology[0] = 2;
        topology[1] = 7;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 6
        topology[0] = 3;
        topology[1] = 0;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 7
        topology[0] = 3;
        topology[1] = 4;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere <> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
    }
    else
    {
        
        topology.resize(4);
        
        for (int nelinf = 0; nelinf < nfatias; nelinf++) {
            // El nel
            
            topology[0] = indinf[nelinf];
            topology[1] = indinf[(nelinf+1)%npts];
            topology[2] = indmid[(nelinf+1)%npts];
            topology[3] = indmid[nelinf];
            {
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
                new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
                SphereRingQ->Geom().SetData(r,xc);
                elementid++;
            }
        }
        
        for (int nelmid = 0; nelmid < nfatias; nelmid++) {
            // El nel
            topology[0] = indmid[nelmid];
            topology[1] = indmid[(nelmid+1)%npts];
            topology[2] = indsup[(nelmid+1)%npts];
            topology[3] = indsup[nelmid];
            {
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
                new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
                SphereRingQ->Geom().SetData(r,xc);
                elementid++;
            }
        }
        
        topology.resize(3);
        
        //        for (int nelsup = 0; nelsup < nfatias; nelsup++) {
        //            // El nel
        ////            int a = indsup[nelsup];
        ////            int b = indsup[(nelsup+1)%npts];
        ////            int c = polo;
        //            topology[0] = indsup[nelsup];
        //            topology[1] = indsup[(nelsup+1)%npts];
        //            topology[2] = polo;
        //            {
        //                TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereRingQ =
        //                new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
        //                SphereRingQ->Geom().SetData(r,xc);
        //                elementid++;
        //            }
        //        }
        
        
    }
    
    // El linha
    // Definition of Arc
    topology.resize(2);
    
    for (int nelinf = 0; nelinf < nfatias; nelinf++) {
        // El nel
        topology[0] = indinf[nelinf];
        topology[1] = indinf[(nelinf+1)%npts];
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    // arcos do topo
    for (int nelsup = 0; nelsup < nfatias; nelsup++) {
        // El nel
        topology[0] = indsup[nelsup];
        topology[1] = indsup[(nelsup+1)%npts];
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    
    // condicoes lado
    if (nfatias<12) {
        topology[0] = 1;
        topology[1] = 3;
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc3 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
        topology[0] = 3;
        topology[1] = 5;
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc3 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
        topology[0] = 0;
        topology[1] = 2;
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc3 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
        topology[0] = 2 ;
        topology[1] = 4;
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc3 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    
    
    //    // Create Geometrical Arc #0
    //    topology[0] = 0;
    //    topology[1] = 1;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc3 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #1
    //    topology[0] = 1;
    //    topology[1] = 2;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc4 */, *geomesh);
    //        // arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #2
    //    topology[0] = 2;
    //    topology[1] = 3;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc1 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #3
    //    topology[0] = 3;
    //    topology[1] = 0;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, arc1 /** arc2 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
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
    
    std::ofstream outfile("malhaFatiaCascaEsfera.vtk");
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
    TPZVec < int64_t > nodeindex(6,0.0);
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
    int64_t arc1 = bc1; // -1;
    int64_t arc2 = bc2; // -2;
    int64_t arc3 = bc3; // -3;
    int64_t arc4 = bc4; // -4;
    
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
    TPZVec < int64_t > nodeindex(3,0.0);
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
    int64_t arc1 = bc1; // -1;
    int64_t arc2 = bc2; // -2;
    int64_t arc3 = bc3; // -3;
    int64_t arc4 = bc4; // -4;
    
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
    TPZVec < int64_t > nodeindex(3,0.0);
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
    int64_t arc1 = bc1; // -1;
    int64_t arc2 = bc2; // -2;
    int64_t arc3 = bc3; // -3;
    int64_t arc4 = bc4; // -4;
    
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
    TPZVec < int64_t > nodeindex(3,0.0);
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
        int64_t index;
        
        TPZInterpolationSpace *bound;
        MElementType elType = gel->Type(side);
        switch(elType)
        {
            case(EOned)://line
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(bound);
                hdivbound->SetSideOrient(3,sideorient);
                break;
            }
            case(ETriangle)://triangle
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeTriang>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeTriang> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeTriang> *>(bound);
                hdivbound->SetSideOrient(6,sideorient);
                break;
            }
            case(EQuadrilateral)://quadrilateral
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeQuad>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(bound);
                hdivbound->SetSideOrient(8,sideorient);
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
        
        int64_t sideconnectindex = intel->ConnectIndex(loccon);
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
    
    TPZVec<REAL> iCoords(3,0.0);
    TPZVec<REAL> iCoordsRotated(3,0.0);
    
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

void RotateNode(TPZVec<REAL> &iCoords, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<REAL> RotationMatrix(3,3,0.0);
    
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
    
    TPZVec<REAL> iCoordsRotated(3,0.0);
    // Apply rotation
    iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
    iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
    iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
    iCoords = iCoordsRotated;
}




