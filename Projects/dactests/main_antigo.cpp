
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

#include "PZMatPoissonD3.h"

#include "pznumeric.h"

#include "TPZExtendGridDimension.h"

#include <iostream>
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


TPZGeoMesh *GMesh(int dimensao, int tipo, REAL Lx, REAL Ly, REAL Lz);
TPZGeoMesh *GMeshDeformed();

TPZGeoMesh *CreateOneCubo();
TPZGeoMesh *CreateGMeshCubo();
TPZGeoMesh *CreateGMeshTetra();
TPZGeoMesh *CreateGMeshPiramide();
TPZGeoMesh *CreateGMeshPrisma();

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, int dim);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcess(TPZAnalysis &an, std::string plotfile);

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv);
void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);

void TestVec(TPZGeoEl *gel, ostream &out);
TPZGeoMesh *MalhaCubo();
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);

//solucao exata
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);

//lado direito da equacao
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno de Neumann
void ForcingBC0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

REAL const Pi = 4.*atan(1.);

bool ftriang = false;//true
bool iscontinuous = true;
int Maps();

/*
 bool ftriang = true;//seta polinomios completos ou totais
 bool isStab = false;//ativa ou nao estabilizacao
 bool iscontinuou = false;//ativa h1 ou l2 para pressao
 bool useh2 = false;//ativa o termo h2 para penalizacao
 REAL delta1 = 0.5;//seta os valores de delta 1 e 2 para a penalizacao
 REAL delta2 = 0.5;
 bool isFullHdiv=false;//seta espaco completo ou nao para o fluxo
 */

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif



//#define MAIN3D
#define MAIN2D

#include "pztransfer.h"

#ifdef MAIN3D
int main3d(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    int dim = 3;
    //TPZGeoMesh *gmesh = GMesh(dim, 2, 1, 1, 1);
    ofstream saidavets("../vetoresTeB.txt",ios::app);
    
    
    //TPZGeoMesh *gmesh = CreateOneCubo();
    //TPZGeoMesh *gmesh = CreateGMeshCubo();
    
    TPZGeoMesh *gmesh = MalhaCubo(); //cubo em tetraedros
    //TPZGeoMesh *gmesh = CreateGMeshPiramide();
    //TPZGeoMesh *gmesh = CreateGMeshPrisma();
    //TPZGeoMesh *gmesh = GMeshDeformed();
    
    ofstream arg200("gmesh1.txt");
    gmesh->Print(arg200);

    //TPZGeoEl *gel;
    //gel =  gmesh->ElementVec()[0];
//    int nele = gmesh->NElements();
//    
//    for (int iel = 0 ; iel < nele; iel++)
//    {
//        TPZGeoEl *gel;
//        gel =  gmesh->ElementVec()[iel]; // vejo a dimensao do elemento (lado) aqui para saber se esta devolvendo os vetores corretamente.
//        int nsides = gel->NSides();
//        for (int is = 0 ; is < nsides; is++)
//        {
//            TPZGeoElSide igeoside(gel,is);
//            
//            if(igeoside.Dimension() < 2)
//            {
//                continue;
//            }
//            int iside = igeoside.Side();
//            TPZAutoPointer<TPZIntPoints> IntegrationRule = gel->CreateSideIntegrationRule(iside, 2);
//            TPZManVector<int,3> order(gel->Dimension(),2);
//            IntegrationRule->SetOrder(order);
//            TPZFMatrix<REAL> directions;
//            TPZVec<int> vectorsides;
//            TPZManVector<REAL,3> qsi(gel->Dimension());
//            int npoints = IntegrationRule->NPoints();
//            for (int ipoint =0; ipoint < npoints; ipoint++) {
//                REAL weitgh;
//                IntegrationRule->Point(ipoint, qsi, weitgh);
//                gel->Directions(iside, qsi, directions, vectorsides);
//            }
//        }
//        //TestVec(gel, saidavets);
//    }
    int nele = gmesh->NElements();
    
    
    TPZFMatrix<REAL> normais;
    TPZVec<int> ladosdasnormais;
    

    for (int iel = 0 ; iel < nele; iel++)
    {
    
        TPZGeoEl *gel;
        gel =  gmesh->ElementVec()[iel]; 
        int nsides = gel->NSides();
        
        normais.Resize(3, nsides*dim);
        ladosdasnormais.Resize(nsides*dim);
        int contlado = 0;
        
        int geldim = gel->Dimension();
        
        for (int is = 0 ; is < nsides; is++)
        {
            TPZGeoElSide igeoside(gel,is);  
            int isdim = igeoside.Dimension();
            if (isdim<2)
            {
                continue;
            }
            TPZVec<int> permutegather;
            if (isdim==2/*spacedim-1*/ && geldim==3/*spacedim*/ )
            {
                gel->HDivPermutation(is, permutegather);
//                cout << "Elemento " << iel << " Lado " << is << endl;
//                cout << permutegather << endl;
            }
            

            
            TPZManVector<REAL,3> center(3); // por no ponto de integracao
            gel->CenterPoint(is, center);
//            igeoside.CenterPoint(center);
            TPZVec<REAL> ptx = center ;
//            cout << "ponto ptx " << ptx << endl;
            TPZFMatrix<REAL> directions;
            TPZVec<int> sidenormals;
            gel->Directions(is, ptx, directions,sidenormals);
            
            TPZFMatrix<REAL> vetores;
            int nvec = directions.Cols();
            vetores.Resize(directions.Rows(), nvec);
            TPZVec<int> lados(nvec,0);
            
            
            if (is<(nsides-1))
            {
                for (int ip=0; ip < nvec; ip++)
                {
                    for(int id = 0; id<3; id++)
                    {
                        vetores(id,ip) = directions(id,permutegather[ip]);
                    }
                    lados[ip] = sidenormals[permutegather[ip]];
                }
                for (int ip = 0 ; ip < nvec; ip++)
                {
                    for(int id = 0; id<3; id++)
                    {
                        normais(id,contlado) = vetores(id,ip);
                    }
                    ladosdasnormais[contlado] = lados[ip];
                    contlado++;
                }
                
            }
            else
            {
                for (int ip = 0 ; ip < nvec; ip++)
                {
                    for(int id = 0; id<3; id++)
                    {
                        normais(id,contlado) = directions(id,ip);
                    }
                    ladosdasnormais[contlado] = sidenormals[ip];
                    contlado++;
                }
            }
            
            
            //determina o indice dos vizinhos de menor dimensao.
//            TPZStack<int> lowdim;
//            gel->LowerDimensionSides(is,lowdim);
            //anexa o indice do lado side
//            lowdim.Push(is);
            
//            for (int nn = 0; nn<4; nn++)
//            {
//                cout << "no " << gel->SideNodeLocIndex(is, nn) << endl;
//            }
            //retornou os mesmos ids dos nos, dados por LowerDimensionSides
            
//            cout << "indices lado " << is << endl;
//            cout << lowdim << endl;
            
            
            if ((is==22 && iel==0)||(is==24 && iel==1)) {
#ifdef LOG4CXX
                std::stringstream sout;
                sout << std::endl;
                sout << "LADO " << is << endl;
                //            directions.Print("direcoes orig ", sout);
                sout << std::endl;
                sout << "sidenormals  orig " << sidenormals << std::endl;
                sout << std::endl;
                sout << "permutacao " << permutegather << std::endl;
                sout << "PERMUTADO LADO " << is << endl;
                //            vetores.Print("direcoes perm ", sout);
                sout << std::endl;
                sout << "sidenormals perm " << lados << std::endl;
                
                LOGPZ_DEBUG(logdata, sout.str());
#endif

            }
            
            
        }
#ifdef LOG4CXX
        std::stringstream sout;
        
        sout << "LadosN " << ladosdasnormais << std::endl;
        sout << "vetores " << normais << std::endl;
        
        LOGPZ_DEBUG(logdata, sout.str());
#endif
        
        
    }
    
    
    // Inicio da relacao entre shape e vetor
    TPZGeoEl *gel;
    gel =  gmesh->ElementVec()[0];
    //MElementType elt = gel::Type(26);
    int64_t index;
    
    MElementType jssj =  gel->Type();
    
    //TPZCompElHDiv< TPZShapeCube > elemhdiv(cmesh1,gel,index);
    
    //TPZVec<std::pair<int,int64_t> >  ShapeAndVec;
    //elemhdiv.IndexShapeToVec(ladosdasnormais,ShapeAndVec,1);

    
    return 0;
    
    
    int p = 1;
    
    int pq = p;
    int pp;
    if(ftriang==true){
        pp = pq-1;
    }else{
        pp = pq;
    }
    
    ofstream arg("gmesh1.txt");
    gmesh->Print(arg);
    
    TPZCompMesh *cmesh1 = CMeshFlux(gmesh, pq, dim);
    TPZCompMesh *cmesh2 = CMeshPressure(gmesh, pp, dim);
    
    
    ofstream arg1("cmeshflux.txt");
    cmesh1->Print(arg1);
    
    ofstream arg2("cmeshpressure.txt");
    cmesh2->Print(arg2);
    
    ofstream arg4("gmesh2.txt");
    gmesh->Print(arg4);
    
    
    //malha multifisica
    TPZVec<TPZCompMesh *> meshvec(2);
    meshvec[0] = cmesh1;
    meshvec[1] = cmesh2;
    
    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec, dim);
    ofstream arg5("cmeshmultiphysics.txt");
    mphysics->Print(arg5);

    return 0;
    
    
    
    
    TPZAnalysis an(mphysics);
    string plotfile("Solution_mphysics.vtk");
    PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
    
    SolveSyst(an, mphysics);
    
    return 0;
}
#endif



#ifdef MAIN2D
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    //int outcome = Maps();
    
    REAL Lx = 1.; // limite inferior e superior do intervalo em x
    REAL Ly = 1.; // limite inferior e superior do intervalo em y
    REAL Lz = 1.; // limite inferior e superior do intervalo em z
    int p = 1;
    int ndiv = 0;
    
    
    // Para dimensao 2
    // tipo 1 triangulo
    // tipo 2 quadrilatero
    
    int dim = 2;
    int tipo = 2;
    ofstream saidaerro("../ErroPoissonHdivMalhaQuad.txt",ios::app);
    //int tipo = 1;
    //ofstream saidaerro("../ErroPoissonHdivMalhaTriang.txt",ios::app);
    ofstream saidavets("../vetoresTeB.txt",ios::app);
    
    for(p=1;p<2;p++)
    {
        int pq = p;
        int pp = p;
        if(ftriang==true){
            pp = pq-1;
        }else{
            pp = pq;
        }
        
        for (ndiv=0; ndiv<1; ndiv++)
        {
//            //  Reading mesh
//            std::string GridFileName;
//            GridFileName += "Crazy.dump";
//            
//            TPZReadGIDGrid GeometryInfo;
//            GeometryInfo.SetfDimensionlessL(1.0);
//            TPZGeoMesh * gmesh2dC = GeometryInfo.GeometricGIDMesh(GridFileName);
//            {
//                //  Print Geometrical Base Mesh
//                std::ofstream argument("GeometicMeshC.txt");
//                gmesh2dC->Print(argument);
//                std::ofstream Dummyfile("GeometricMeshC.vtk");
//                TPZVTKGeoMesh::PrintGMeshVTK(gmesh2dC,Dummyfile, true);
//            }
//            
//            TPZExtendGridDimension extend(gmesh2dC,layerthickness);
//            TPZGeoMesh *gmesh3d = extend.ExtendedMesh(1,-7,-8);
//            
//            ofstream arq3("gmesh3d.txt");
//            gmesh3d->Print(arq3);
//            TPZFMatrix<STATE> Tr(3,3,0.0);
//            //            Tr(0,0)=1.0;
//            //            Tr(1,1)=1.0;
//            //            Tr(2,2)=1.0;
//            
//            Tr(0,0)=0.5;
//            Tr(0,1)=2.1;
//            Tr(0,2)=0.0;
//            
//            Tr(1,0)=1.0;
//            Tr(1,1)=1.8;
//            Tr(1,2)=0.0;
//            
//            Tr(2,0)=1.0;
//            Tr(2,1)=2.0;
//            Tr(2,2)=1.0;
//            TPZExtendGridDimension::DeformMesh(Tr, gmesh3d);
//            TPZGeoMesh *gmesh = gmesh3d;
//            dim = 3;
//            {
//                //  Print Geometrical Base Mesh
//                std::ofstream Dummyfile("GeometricMesh.vtk");
//                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
//            }


            
             
            //TPZGeoMesh *gmesh = GMeshDeformed();
            TPZGeoMesh *gmesh2d = GMesh(dim, tipo, Lx, Ly, Lz);
            REAL layerthickness = 1.;
            UniformRefine(gmesh2d, ndiv);
            
            {
                //  Print Geometrical Base Mesh
                std::ofstream Dummyfile("GeometricMesh2D.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh2d,Dummyfile, true);
            }
            
            TPZExtendGridDimension extend(gmesh2d,layerthickness);
            TPZGeoMesh *gmesh3d = extend.ExtendedMesh(1,-5,-6);
            ofstream arg("gmesh1.txt");
            gmesh2d->Print(arg);
            
            ofstream arq3("gmesh3d.txt");
            gmesh3d->Print(arq3);
            TPZFMatrix<STATE> Tr(3,3,0.0);

            TPZGeoMesh *gmesh = gmesh3d;
            dim = 3;
            
            
//            TPZGeoMesh *gmesh = gmesh2d;
//            dim = 2;
            {
                //  Print Geometrical Base Mesh
                std::ofstream Dummyfile("GeometricMesh.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
            }
            
//            
//            {
//                //	Print Geometrical Base Mesh
//                std::ofstream Dummyfile("GeometricMesh.vtk");
//                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
//            }
//            
            TPZCompMesh *cmesh2 = CMeshPressure(gmesh, pp, dim);
            TPZCompMesh *cmesh1 = CMeshFlux(gmesh, pq, dim);


            
            ofstream arg1("cmeshflux.txt");
            cmesh1->Print(arg1);
//
//            ofstream arg2("cmeshpressure.txt");
//            cmesh2->Print(arg2);
//            
//            ofstream arg4("gmesh2.txt");
//            gmesh->Print(arg4);
            
            
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            
            
//            int Prows =cmesh2->Solution().Rows();
//            TPZFMatrix<STATE> alphaP(Prows,1,0.0);
//            alphaP(26,0)=1.0;
//            cmesh2->LoadSolution(alphaP);
//            TPZAnalysis anP(cmesh2);
//            std::string plotfile2("AlphasPressure.vtk");
//            PosProcess(anP, plotfile2);
            
            
            
            TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec, dim);
            ofstream arg5("cmeshmultiphysics.txt");
            mphysics->Print(arg5);
            TPZCompEl *cel = mphysics->Element(0);
            TPZElementMatrix ek,ef;
            cel->CalcStiff(ek, ef);
            
            TPZAnalysis an(mphysics);

            SolveSyst(an, mphysics);
            std::string plotfile("GSaida.vtk");
            PosProcessMultphysics(meshvec,  mphysics, an, plotfile);

            
            
          //
            //    TPZAutoPointer< TPZMatrix<REAL> > matKdoAna;
            //    matKdoAna = an.Solver().Matrix();
            //
            //#ifdef LOG4CXX
            //    if(logdata->isDebugEnabled())
            //    {
            //        std::stringstream sout;
            //        matKdoAna->Print("Eglobal = ", sout,EMathematicaInput);
            //        an.Rhs().Print("Fglobal = ", sout,EMathematicaInput);
            //        an.Solution().Print("Sglobal = ", sout,EMathematicaInput);
            //        LOGPZ_DEBUG(logdata,sout.str())
            //    }
            //#endif
            
            
            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            TPZVec<REAL> erros;
            
            //            TPZManVector<REAL,3> myerrors(3,0.);
            //            an.SetExact(SolExata);
            //            an.PostProcessError(myerrors);
            
            //            saidaerro << "Valor de epsilone " << EPSILON << std::endl;
            //            saidaerro << "Numero de threads " << numthreads << std::endl;
//            saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
//            ErrorHDiv(cmesh1, saidaerro, p, ndiv);
//            
//            saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
//            ErrorL2(cmesh2, saidaerro, p, ndiv);
//            
            std::cout << "Postprocessed\n";
            
            //Plot da solucao aproximada
            // para estudo do ero nao precisa gerar dados de saida
            //PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
            
            std::cout<< " grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
            
            
        }
        
    }
    
    
    
    
    std::cout<< " fim " << std::endl;
    
	return EXIT_SUCCESS;
}
#endif

TPZGeoMesh *GMesh(int d, int tipo, REAL Lx, REAL Ly, REAL Lz){
    
    int dim = d;
    if(d<2 || d>3)
    {
        std::cout << "dimensao errada" << std::endl;
        dim = 2;
        DebugStop();
    }
    
    int Qnodes = dim == 2 ? 4 : 8;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
    TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
    
    //indice dos nos
    int64_t id = 0;
    REAL valx;
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = xi*Lx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,0. );//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = Lx - xi*Lx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,Ly);//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    //indice dos elementos
    id = 0;
    
    if(tipo==1) // triangulo
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
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
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    
    gmesh->BuildConnectivity();
    
    
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

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
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


TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
	// int dim = 2;
	TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
	TPZMaterial * mat(material);
	material->NStateVariables();
    
    //  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, -5,neumann , val1, val2);
    TPZMaterial * BCond5 = material->CreateBC(mat, -6,neumann, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
	
    cmesh->SetDimModel(dim);
    
    //cmesh->SetAllCreateFunctionsHDivFull();
    cmesh->SetAllCreateFunctionsHDiv();
    
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    
	cmesh->SetDefaultOrder(pOrder);
    
	
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
	//int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
	material->NStateVariables();
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, -5,neumann, val1, val2);
    TPZMaterial * BCond5 = material->CreateBC(mat, -6,neumann, val1, val2);
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    if(iscontinuous == true){
        cmesh->SetAllCreateFunctionsContinuous();
    }
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    if (iscontinuous==false) {
        int ncon = cmesh->NConnects();
        for(int i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            //newnod.SetPressure(true);
            newnod.SetLagrangeMultiplier(1);
        }
        
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

TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, int dim)
{
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    //int dim =2;
    bool intface;
    TPZMatPoissonD3 *material = new TPZMatPoissonD3(matId,dim); intface = true; // nesse material tem que ser true
    //TPZMixedPoisson *material = new TPZMixedPoisson(matId,dim); intface = false; // nesse material tem que ser false
    
    //incluindo os dados do problema
    //incluindo os dados do problema
    if (!intface) {
        TPZFNMatrix<2,REAL> PermTensor(dim,dim,0.);
        TPZFNMatrix<2,REAL> InvPermTensor(dim,dim,0.);
        for (int i=0; i<dim; i++)
        {
            PermTensor(i,i) = 1.0;
        }
        InvPermTensor=PermTensor;
        material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    }
    
    
    //    REAL coefk = 1.;
    //    REAL coefvisc = 1.;
    //    material->SetPermeability(coefk);
    //    material->SetViscosity(coefvisc);
    //
    //    if(isStab==true){
    //        material->SetStabilizedMethod();
    //        material->SetStabilizationCoeficients(delta1,delta2);
    //	}
    //    if(isStab==true && useh2==true) material->SetHdois();
    //
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing, 5);
    dum->SetPolynomialOrder(4);
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
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond4 = material->CreateBC(mat, -5,neumann, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond5 = material->CreateBC(mat, -6,neumann, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);
    
    //Fazendo auto build
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    
	// Creation of interface elements
    if (intface)
    {
        int nel = mphysics->ElementVec().NElements();
        for(int el = 0; el < nel; el++)
        {
            TPZCompEl * compEl = mphysics->ElementVec()[el];
            if(!compEl) continue;
            int index = compEl ->Index();
            if(compEl->Dimension() == mphysics->Dimension())
            {
                TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(mphysics->ElementVec()[index]);
                if(!InterpEl) continue;
                InterpEl->CreateInterfaces();
            }
        }

    }
	   
    return mphysics;
}

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
    //    TPZSkylineNSymStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt); //caso simetrico
    //	step.SetDirect(ELU);
	an.SetSolver(step);
    //    an.Assemble();
	an.Run();
    
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void PosProcess(TPZAnalysis &an, std::string plotfile){
    TPZManVector<std::string,10> scalnames(1), vecnames(0);
    scalnames[0]= "Solution";
    
    const int dim = 3;
    int div = 2;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    std::ofstream out("malhaNormal.txt");
    an.Print("nothing",out);
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(3), vecnames(3);
	vecnames[0]  = "Flux";
    vecnames[1]  = "GradFluxX";
    vecnames[2]  = "GradFluxY";
    scalnames[0] = "Pressure";
    scalnames[1] = "H1ErrorPerArea";
    scalnames[2] = "ExactPressure";
    
	const int dim = mphysics->Dimension();
	int div =2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
    
    //    mphysics->Solution().Print("Solucao");
    
}

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    // (dimensao espacial == 2)
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    solp[0] = sin(Pi*x)*sin(Pi*y);
    flux(0,0)=-Pi*cos(Pi*x)*sin(Pi*y);
    flux(1,0)=-Pi*cos(Pi*y)*sin(Pi*x);
    flux(2,0)=2*Pi*Pi*sin(Pi*y)*sin(Pi*x);
   
    // (dimensao espacial == 3)
//    solp.Resize(1, 0.);
//    flux.Resize(4, 1.);
//    flux(0,0)=flux(1,0)=flux(2,0)=flux(3,0)=0.;
//    double x = pt[0];
//    double y = pt[1];
//    double z = pt[2];
//    solp[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
//    flux(0,0)=-Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
//    flux(1,0)=-Pi*cos(Pi*y)*sin(Pi*x)*sin(Pi*z);
//    flux(2,0)=-Pi*cos(Pi*z)*sin(Pi*x)*sin(Pi*y);
//    flux(3,0)=3.*Pi*Pi*sin(Pi*y)*sin(Pi*x)*sin(Pi*z);
    
    
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    // (dimensao espacial == 2)
    double x = pt[0];
    double y = pt[1];
    disp[0] = 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
    
    // (dimensao espacial == 3)
//    double x = pt[0];
//    double y = pt[1];
//    double z = pt[2];
//    disp[0] = 3.*Pi*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
    
	
}

void ForcingBC0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    disp[0]= Pi*cos(Pi*y)*sin(Pi*x);
}

void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    disp[0]= -Pi*cos(Pi*y)*sin(Pi*x);
}

int Maps()
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(2);
    
    // Create nodes
    int MatID = 1;
    int nodeId=0;
    
    TPZVec<REAL> coord(3,0.);
    TPZVec<int64_t> TopologyPoint(1);
    TPZVec<int64_t> TopologyLinear(2);
    
    coord[0]=0.0; // xcoord
    coord[1]=0.0; // ycoord
    coord[2]=0.0; // zcoord
    gmesh->NodeVec()[nodeId].SetCoord(coord);
    gmesh->NodeVec()[nodeId].SetNodeId(nodeId);
    TopologyPoint[0] = nodeId;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > *Mypoint1=new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (nodeId++, TopologyPoint, MatID,*gmesh);
    
    coord[0]=1.0; // xcoord
    coord[1]=1.0; // ycoord
    coord[2]=1.0; // zcoord
    gmesh->NodeVec()[nodeId].SetCoord(coord);
    gmesh->NodeVec()[nodeId].SetNodeId(nodeId);
    TopologyPoint[0] = nodeId;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > *Mypoint2=new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (nodeId++, TopologyPoint, MatID,*gmesh);
    
    // Creating line
    TopologyLinear[0] = 0;
    TopologyLinear[1] = 1;
    TPZGeoElRefPattern < pzgeom::TPZGeoLinear > *Myline=new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (nodeId++, TopologyLinear, 2,*gmesh);
    
    
    
    
    TPZVec<REAL> qsi(1,1.0);
    TPZVec<REAL> coordX(3,1.0);
    
    Myline->X(qsi, coordX);
    std::cout << "Xcoord = " << coordX[0] << std::endl;
    std::cout << "Ycoord = " << coordX[1] << std::endl;
    std::cout << "Zcoord = " << coordX[2] << std::endl;
    TPZFMatrix<REAL> jac,axes,jacinverse;
    REAL detjac;
    
    Myline->Jacobian(qsi, jac, axes, detjac, jacinverse);
    
    std::cout << "jacx = " << jac(0,0) << std::endl;
    std::cout << "jacinversex = " << jacinverse(0,0) << std::endl;
    
    //	std::string GeoGridFile;
    //    //	GeoGridFile = "SQDomain.dump";
    //	GeoGridFile = "MeshTest.dump";
    //    TPZReadGIDGrid GIDMesh;
    //    TPZGeoMesh *gmesh = GIDMesh.GeometricGIDMesh(GeoGridFile);
    //
    //	{
    //		//	Print Geometrical Base Mesh
    //		std::ofstream Dummyfile("GeometricGID.vtk");
    //		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    //	}
    
    gmesh->BuildConnectivity();
    gmesh->Print();
    
    return 0;
    
}


void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv)
{
    int64_t nel = hdivmesh->NElements();
    //int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    out << "L2 Norm for flux - L2 Norm for divergence - Hdiv Norm for flux " << endl;
    out <<  setw(16) << sqrt(globalerrors[1]) << setw(25)  << sqrt(globalerrors[2]) << setw(21)  << sqrt(globalerrors[3]) << endl;
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
        cel->EvaluateError(SolExata, elerror, 0);
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
//REAL determinante(TPZVec<REAL> &v1,TPZVec<REAL> &v2, TPZVec<REAL> &v3 );
//REAL determinante(TPZVec<REAL> &v1,TPZVec<REAL> &v2, TPZVec<REAL> &v3 )
//{
//    return (v1[0]*v2[1]*v3[2]+v1[1]*v2[2]*v3[0]+v1[2]*v2[0]*v3[1]-v1[2]*v2[1]*v3[0]-v1[0]*v2[2]*v3[1]-v1[1]*v2[0]*v3[2]);
//}

void TestVec(TPZGeoEl *gel, ostream &out)
{
    
    
    
    int nsides = gel->NSides();
    for(int side=0; side<nsides; side++)
    {
        int numbernormals = 0;
        
        int dimension = gel->Dimension();
        TPZManVector<REAL,3> center(dimension);
        
        int sidedimension = gel->SideDimension(side);
        
        //determina o indice dos vizinhos de menor dimensao.
        TPZStack<int> lowdim;
        gel->LowerDimensionSides(side,lowdim);
        //anexa o indice do lado side
        lowdim.Push(side);
        
        // the normals corresponding to the internal shape functions
        // Compute the number of normals we need to compute
        if(sidedimension == dimension-1)
        {
            numbernormals = lowdim.NElements();
        }
        else if(sidedimension == dimension)
        {
            numbernormals = nsides*dimension;
        }
        else
        {
            numbernormals = 0;
        }
        
        if(!numbernormals)
        {
            continue;
        }


        gel->CenterPoint(nsides-1, center);
        TPZVec<REAL> ptx = center ;
        TPZFMatrix<REAL> directions;
        TPZVec<int> sidenormals;
        gel->Directions(side, ptx, directions,sidenormals);
        
#ifdef LOG4CXX
        std::stringstream sout;
        sout << std::endl;
        directions.Print("direcoes ", sout);
        sout << std::endl;
        sout << "sidenormals " << sidenormals << std::endl;
        sout << std::endl;
        sout << "lowdim " << lowdim << std::endl;
        LOGPZ_DEBUG(logdata, sout.str());
#endif
        int ind = TPZGeoElSide(gel,side).Id();
        std::cout << " indice do elemento " << gel->Dimension() << "D  eh " << ind << std::endl;
        
    }
    
    

    
    
    

    
}

TPZGeoMesh *MalhaCubo()
{
	int64_t numnodes=-1;
	int64_t numelements=-1;
	
	string FileName, dirname = PZSOURCEDIR;
	FileName = dirname + "/Projects/dactests/";
	FileName += "cube1.txt";
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		ifstream read (FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements") countelements = false;
			if(countelements) numelements++;
		}
	}
	
	TPZGeoMesh * gMesh = new TPZGeoMesh;
	
	gMesh -> NodeVec().Resize(numnodes);
	
	TPZManVector <int64_t> TopolTetra(4);
	
	const int64_t Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int64_t nodeId = 0, elementId = 0, matElId = 1;
	
	ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int64_t in;
	for(in=0; in<numnodes; in++)
	{
		read >> nodeId;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		Node[nodeId-1].SetNodeId(nodeId);
		Node[nodeId-1].SetCoord(0,nodecoordX);
		Node[nodeId-1].SetCoord(1,nodecoordY);
		Node[nodeId-1].SetCoord(2,nodecoordZ);
		gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
	}
	
	{
		read.close();
		read.open(FileName.c_str());
		
		int64_t l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		
		int64_t el;
		int neumann1 = -4, neumann2 = -5;
		//std::set<int> ncoordz; //jeitoCaju
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> TopolTetra[0]; //node 1
			read >> TopolTetra[1]; //node 2
			read >> TopolTetra[2]; //node 3
			read >> TopolTetra[3]; //node 4
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
			TopolTetra[0]--;
			TopolTetra[1]--;
			TopolTetra[2]--;
			TopolTetra[3]--;
			
			int64_t index = el;
			
			new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		}
		
		gMesh->BuildConnectivity();
		
		// Colocando as condicoes de contorno
		for(el=0; el<numelements; el++)
		{
			TPZManVector <TPZGeoNode,4> Nodefinder(4);
			TPZManVector <REAL,3> nodecoord(3);
			TPZGeoEl *tetra = gMesh->ElementVec()[el];
			
			// na face x = 1
			TPZVec<int64_t> ncoordzVec(0); int64_t sizeOfVec = 0;
			for (int i = 0; i < 4; i++)
			{
				int64_t pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[0] == 1.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = pos;
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,neumann1);
			}
			
			// Na face x = -1
			ncoordzVec.Resize(0);
			sizeOfVec = 0;
			for (int i = 0; i < 4; i++)
			{
				int64_t pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];
				
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[0] == -1.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = pos;
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,neumann2);
			}
			
		}
		
		TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
		yz[0] = 1.;
		z[2] = -1;
		int bcidxyz = -1, bcidyz = -2, bcidz = -3;
		SetPointBC(gMesh, xyz, bcidxyz);
		SetPointBC(gMesh, yz, bcidyz);
		SetPointBC(gMesh, z, bcidz);
		
	}
	
	ofstream arg("malhaPZ1BC.txt");
	gMesh->Print(arg);
	
	std::ofstream out("Cube.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gMesh, out, true);
	
	return gMesh;
}

/// Generate a boundary geometric element at the indicated node
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
{
	// look for an element/corner node whose distance is close to start
	TPZGeoNode *gn1 = gr->FindNode(x);
	int64_t iel;
	int64_t nelem = gr->ElementVec().NElements();
	TPZGeoEl *gel;
	for (iel = 0; iel<nelem; iel++) {
		gel = gr->ElementVec()[iel];
		if(!gel) continue;
		int nc = gel->NCornerNodes();
		int c;
		for (c=0; c<nc; c++) {
			TPZGeoNode *gn = gel->NodePtr(c);
			if (gn == gn1) {
				break;
			}
		}
		if (c<nc) {
			TPZGeoElBC(gel, c, bc);
			return;
		}
	}
}

TPZGeoMesh *CreateOneCubo()
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int nnodes = 8;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    
    //c0
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c1
    coord[0] = 0.5;
    coord[1] = 0;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c2
    coord[0] = 0.5;
    coord[1] = 1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c3
    coord[0] = 0;
    coord[1] = 1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c4
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c5
    coord[0] = 0.5;
    coord[1] = 0;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c6
    coord[0] = 0.5;
    coord[1] = 1;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c7
    coord[0] = 0;
    coord[1] = 1;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    
     
    
    int index = 0;
    TPZManVector<int64_t,8> TopolCubo(8,0);
    
    
    TopolCubo[0] = 0;
    TopolCubo[1] = 1;
    TopolCubo[2] = 2;
    TopolCubo[3] = 3;
    TopolCubo[4] = 4;
    TopolCubo[5] = 5;
    TopolCubo[6] = 6;
    TopolCubo[7] = 7;
    
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index, TopolCubo, matId, *gmesh);
    
    
    gmesh->BuildConnectivity();
    
    /// gmesh para aqui
    
    int nref = 0;
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
    
    
    std::ofstream out("SingleCubeDAC.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZGeoMesh *CreateGMeshCubo()
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int nnodes = 16;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
   
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    
    //c0
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c1
    coord[0] = 0.5;
    coord[1] = 0;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c2
    coord[0] = 0.5;
    coord[1] = 1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c3
    coord[0] = 0;
    coord[1] = 1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c4
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c5
    coord[0] = 0.5;
    coord[1] = 0;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c6
    coord[0] = 0.5;
    coord[1] = 1;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c7
    coord[0] = 0;
    coord[1] = 1;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    
    //c1
    coord[0] = 1;
    coord[1] = 0;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c2
    coord[0] = 1;
    coord[1] = 1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c5
    coord[0] = 1;
    coord[1] = 0;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c6
    coord[0] = 1;
    coord[1] = 1;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;

    
    
    int index = 0;
    TPZManVector<int64_t,8> TopolCubo(8,0);

    
    TopolCubo[0] = 0;
    TopolCubo[1] = 1;
    TopolCubo[2] = 2;
    TopolCubo[3] = 3;
    TopolCubo[4] = 4;
    TopolCubo[5] = 5;
    TopolCubo[6] = 6;
    TopolCubo[7] = 7;
    
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index, TopolCubo, matId, *gmesh);
    
    
    TopolCubo[0] = 1;
    TopolCubo[1] = 8;
    TopolCubo[2] = 9;
    TopolCubo[3] = 2;
    TopolCubo[4] = 5;
    TopolCubo[5] = 10;
    TopolCubo[6] = 11;
    TopolCubo[7] = 6;
    
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index, TopolCubo, matId, *gmesh);
    
    
    gmesh->BuildConnectivity();
    
    /// gmesh para aqui
    
    int nref = 0;
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
    
    
    std::ofstream out("CubeDAC.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZGeoMesh *CreateGMeshPiramide()
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int matid = 1;
    int nnodes = 5;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    //REAL l = 1;
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    //c0
    coord[0] = -1;
    coord[1] = -1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c1
    coord[0] = 1;
    coord[1] = -1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c2
    coord[0] = 1;
    coord[1] = 1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c3
    coord[0] = -1;
    coord[1] = 1;
    coord[2] = 0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c4
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    

    int index = 0;
    TPZManVector<int64_t,5> TopolPiramide(5,0);
    
    TopolPiramide[0] = 0;
    TopolPiramide[1] = 1;
    TopolPiramide[2] = 2;
    TopolPiramide[3] = 3;
    TopolPiramide[4] = 4;

    
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid> (index, TopolPiramide, matid, *gmesh);
    
    gmesh->BuildConnectivity();
    
    int nref = 0;
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
    
    
    std::ofstream out("PiramideDAC.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZGeoMesh *CreateGMeshPrisma()
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int matid = 1;
    int nnodes = 6;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    //REAL l = 1;
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    //c0
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;//-1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c1
    coord[0] = 1;
    coord[1] = 0;
    coord[2] = 0;//-1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c2
    coord[0] = 0;
    coord[1] = 1;
    coord[2] = 0;//-1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c3
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c4
    coord[0] = 1;
    coord[1] = 0;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c5
    coord[0] = 0;
    coord[1] = 1;
    coord[2] = 1;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    
    int index = 0;
    TPZManVector<int64_t,6> TopolPrisma(6,0);
    
    TopolPrisma[0] = 0;
    TopolPrisma[1] = 1;
    TopolPrisma[2] = 2;
    TopolPrisma[3] = 3;
    TopolPrisma[4] = 4;
    TopolPrisma[5] = 5;
    
    
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (index, TopolPrisma, matid, *gmesh);
    
    gmesh->BuildConnectivity();
    
    int nref = 0;
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
    
    
    std::ofstream out("PiramideDAC.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

/*TPZGeoMesh *GMeshDeformed(int d, int tipo, REAL Lx, REAL Ly, REAL Lz){
    
    int dim = d;
    if(d<2 || d>3)
    {
        std::cout << "dimensao errada" << std::endl;
        dim = 2;
        DebugStop();
    }
    
    int Qnodes = dim == 2 ? 4 : 8;
    REAL alpha = 0.125;
    REAL beta = 0.5;
    int totalnodes = Qnodes + 4;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(totalnodes);
    TPZVec<TPZGeoNode> Node(totalnodes);
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
    TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolArc(3);
    TPZVec <int64_t> TopolPoint(1);
    
    //indice dos nos
    int64_t id = 0;
    REAL valx;
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = xi*Lx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,0. );//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = Lx - xi*Lx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,Ly);//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,-alpha*Lx );//coord X
    Node[id].SetCoord(1 ,beta*Ly);//coord Y
    gmesh->NodeVec()[id] = Node[id];
    id++; //id 4
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 , beta*Lx );//coord X
    Node[id].SetCoord(1 ,-alpha*Ly);//coord Y
    gmesh->NodeVec()[id] = Node[id];
    id++; //id 5
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 , Lx - alpha*Lx );//coord X
    Node[id].SetCoord(1 ,beta*Ly);//coord Y
    gmesh->NodeVec()[id] = Node[id];
    id++; //id 6
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 , beta*Lx );//coord X
    Node[id].SetCoord(1 ,Ly - alpha*Ly);//coord Y
    gmesh->NodeVec()[id] = Node[id];
    id++; //id 7
    
    //indice dos elementos
    id = 0;
    
    if(tipo==1) // triangulo
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolArc[0] = 0;
        TopolArc[1] = 5;
        TopolArc[2] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,bc0,*gmesh);
        id++;
        
        TopolArc[0] = 2;
        TopolArc[1] = 6;
        TopolArc[2] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,bc1,*gmesh);
        id++;
        
        TopolArc[0] = 3;
        TopolArc[1] = 7;
        TopolArc[2] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,bc2,*gmesh);
        id++;
        
        TopolArc[0] = 3;
        TopolArc[1] = 4;
        TopolArc[2] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,bc3,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (id,TopolQuad,matId,*gmesh);
        id++;
        
        TopolArc[0] = 5;
        TopolArc[1] = 0;
        TopolArc[2] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,bc0,*gmesh);
        id++;
        
        TopolArc[0] = 2;
        TopolArc[1] = 6;
        TopolArc[2] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,bc1,*gmesh);
        id++;
        
        TopolArc[0] = 3;
        TopolArc[1] = 7;
        TopolArc[2] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,bc2,*gmesh);
        id++;
        
        TopolArc[0] = 3;
        TopolArc[1] = 4;
        TopolArc[2] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,bc3,*gmesh);
    }
    
    gmesh->BuildConnectivity();
    
    
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
}*/

TPZGeoMesh *GMeshDeformed(){
    
    // description of Geometry and application
    // 2D Cylindrical Domain boundaries
    int matId = 1;
    int arc1 = -1;
    int arc2 = -2;
    int arc3 = -3;
    int arc4 = -4;
    int Point1 = -11;
    int Point2 = -22;
    int Point3 = -33;
    int Point4 = -44;
    int WellPoint = -55;
    
    int nodenumber = 9;
    REAL ModelRadius = 1.0;
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordantes for Arc3D 1
    int id = 0;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    id++;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    id++;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2)*ModelRadius/2.);//coord Y
    id++;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2)*ModelRadius/2.);//coord Y
    id++;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y
    id++;
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y
    id++;
    
    int elementid = 0;
    // Create Geometrical Arc #1
    // Definition of Arc coordenates
    TPZVec < int64_t > nodeindex(3,0.0);
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 5;
    TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc1 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1,*gmesh);
    elementid++;
    
    // Create Geometrical Arc #2
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 6;
    TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc2 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2,*gmesh);
    elementid++;
    
    // Create Geometrical Arc #3
    nodeindex[0] = 3;
    nodeindex[1] = 4;
    nodeindex[2] = 7;
    TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc3 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3,*gmesh);
    elementid++;
    
    // Create Geometrical Arc #4
    nodeindex[0] = 4;
    nodeindex[1] = 1;
    nodeindex[2] = 8;
    TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc4 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4,*gmesh);
    elementid++;
    
    // Create Geometrical Point #1
    nodeindex.resize(1);
    nodeindex[0] = 1;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint1 = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, Point1,*gmesh);
    elementid++;
    
    // Create Geometrical Point #2
    nodeindex.resize(1);
    nodeindex[0] = 3;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint2 = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, Point2,*gmesh);
    elementid++;
    
    // Create Geometrical Point #3
    nodeindex.resize(1);
    nodeindex[0] = 2;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint3 = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, Point3,*gmesh);
    elementid++;
    
    // Create Geometrical Point #4
    nodeindex.resize(1);
    nodeindex[0] = 4;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint4 = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, Point4,*gmesh);
    elementid++;
    
    // Create Geometrical Point for fluid injection or Production #1
    nodeindex.resize(1);
    nodeindex[0] = 0;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, WellPoint,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #1
    nodeindex.resize(4);
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 3;
    nodeindex[3] = 4;
    TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > * Quad1 = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
//    // Create Geometrical triangle #1
//    nodeindex.resize(3);
//    nodeindex[0] = 0;
//    nodeindex[1] = 1;
//    nodeindex[2] = 2;
//    TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > *triangle1 = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
//    elementid++;
//    
//    // Create Geometrical triangle #2
//    nodeindex[0] = 0;
//    nodeindex[1] = 2;
//    nodeindex[2] = 3;		
//    TPZGeoElRefPattern<pzgeom::TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle2 = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend< pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
//    elementid++;
//    
//    // Create Geometrical triangle #3		
//    nodeindex[0] = 0;
//    nodeindex[1] = 3;
//    nodeindex[2] = 4;		
//    TPZGeoElRefPattern<pzgeom::TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle3 = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
//    elementid++;
//    
//    // Create Geometrical triangle #4		
//    nodeindex[0] = 0;
//    nodeindex[1] = 4;
//    nodeindex[2] = 1;		
//    TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle4 = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend< pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
//    elementid++;
    

        

    gmesh->BuildConnectivity();
    
    int nref = 0;
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
    
    
    std::ofstream out("CurvoDAC.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}
