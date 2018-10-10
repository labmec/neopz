
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

// just for print data
/** @brief Map used norms */
std::map<REAL,REAL> fDebugMapL2, fDebugMapHdiv;


TPZGeoMesh *GMesh2dpol(int dimensao, int tipo, int ndiv);
TPZGeoMesh *GMeshDeformed2dpol();
TPZGeoMesh *CreateOneCubo2dpol(int nref=0);

TPZCompMesh *CMeshFlux2dpol(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshPressure2dpol(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshMixed2dpol(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);

void UniformRefine2dpol(TPZGeoMesh* gmesh, int nDiv);
void SolveSyst2dpol(TPZAnalysis &an, TPZCompMesh *fCmesh);
void PosProcessMultphysics2dpol(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcess2dpol(TPZAnalysis &an, std::string plotfile);

void ErrorL22dpol(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv);
void ErrorHDiv2dpol(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);

/** @brief Prints debug map in Mathematica style */
void PrintDebugMapForMathematica2dpol(std::string filenameHdiv, std::string filenameL2);


//solucao exata
void SolExata2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);

//lado direito da equacao
void Forcing2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno de Dirichlet
void ForcingBC0D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC1D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC2D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC3D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC4D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC5D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno de Neumann
void ForcingBC0N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC1N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC2N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC3N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC4N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC5N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

int dim = 2;
REAL aa = 0.0;
REAL bb = 0.0;
REAL cc = 0.0;
REAL Epsilon = 0.4;
// tensor de permutacao
TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);

// Para dimensao 2
// tipo 1 triangulo
// tipo 2 quadrilatero
int tipo = 2;
bool ftriang = false;//true

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

#include "pztransfer.h"

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    int p = 1;
    int ndiv = 0;
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    
    
    
    
    ofstream saidaerro("../ErroPoissonHdivMalhaQuad.txt",ios::app);
    //int tipo = 1;
    //ofstream saidaerro("../ErroPoissonHdivMalhaTriang.txt",ios::app);
    
    for(p=1;p<3;p++)
    {
        int pq = p;
        int pp = p;
        if(ftriang==true){
            pp = pq-1;
        }else{
            pp = pq;
        }
        
        for (ndiv=0; ndiv<3; ndiv++)
        {
            
            //            TPZGeoMesh *gmesh = CreateOneCubo(ndiv);
            
            TPZGeoMesh *gmesh = GMesh2dpol( dim, tipo, ndiv);
//            ofstream arg1("gmesh.txt");
//            gmesh->Print(arg1);
            {
                //  Print Geometrical Base Mesh
                std::ofstream Dummyfile("GeometricMesh.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
            }
            
            TPZCompMesh *cmesh2 = CMeshPressure2dpol(gmesh, pp, dim);
            TPZCompMesh *cmesh1 = CMeshFlux2dpol(gmesh, pq, dim);
            
            
            
            //            ofstream arg1("cmeshflux.txt");
            //            cmesh1->Print(arg1);
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
            
            
            
            TPZCompMesh * mphysics = CMeshMixed2dpol(gmesh,meshvec);
            //            ofstream arg5("cmeshmultiphysics.txt");
            //            mphysics->Print(arg5);
            //            TPZCompEl *cel = mphysics->Element(0);
            //            TPZElementMatrix ek,ef;
            //            cel->CalcStiff(ek, ef);
            
            TPZAnalysis an(mphysics);
            
            SolveSyst2dpol(an, mphysics);
            std::string plotfile("OurSolution.vtk");
            PosProcessMultphysics2dpol(meshvec,  mphysics, an, plotfile);
            
            
            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            TPZVec<REAL> erros;
            
            //            TPZManVector<REAL,3> myerrors(3,0.);
            //            an.SetExact(SolExata);
            //            an.PostProcessError(myerrors);
            
            //            saidaerro << "Valor de epsilone " << EPSILON << std::endl;
            //            saidaerro << "Numero de threads " << numthreads << std::endl;
            saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
            ErrorHDiv2dpol(cmesh1, saidaerro, p, ndiv);
            
            saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
            ErrorL22dpol(cmesh2, saidaerro, p, ndiv);
            //
            std::cout << "Postprocessed\n";
            
            //Plot da solucao aproximada
            // para estudo do ero nao precisa gerar dados de saida
            //PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
            
            stringstream ss;
            ss << p;
            string str = ss.str();
            
            std::cout<< " grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
            std::string filename("InputData");
            std::string L2("L2.txt");
            std::string Hdiv("Hdiv.txt");
            std::string HdivData,L2Data;
            HdivData = filename+str+Hdiv;
            L2Data = filename+str+L2;
            
            PrintDebugMapForMathematica2dpol(HdivData,L2Data);
            
        }
        
    }
    
    
    std::cout<< " fim " << std::endl;
    
    return EXIT_SUCCESS;
}

TPZGeoMesh *GMesh2dpol(int d, int tipo, int ndiv)
{
    
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
    
    gmesh->SetDimension(dim);
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
    TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
    
    //indice dos nos
    int64_t id = 0;
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
    coord[0] = -1.0;
    coord[1] = -1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c1
    coord[0] =  1.0;
    coord[1] = -1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c2
    coord[0] =  1.0;
    coord[1] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c3
    coord[0] = -1.0;
    coord[1] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
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

TPZGeoMesh *CreateOneCubo2dpol(int nref)
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
    //c0
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c1
    coord[0] =  1.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c2
    coord[0] =  1.0;
    coord[1] =  1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c3
    coord[0] = 0.0;
    coord[1] =  1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;

    //c4
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c5
    coord[0] =  1.0;
    coord[1] = 0.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c6
    coord[0] =  1.0;
    coord[1] =  1.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    //c7
    coord[0] = 0.0;
    coord[1] =  1.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in); in++;
    
    // cubo [-1,1]^3
//    //c0
//    coord[0] = -1.0;
//    coord[1] = -1.0;
//    coord[2] = -1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in++);
//    //c1
//    coord[0] =  1.0;
//    coord[1] = -1.0;
//    coord[2] = -1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in++);
//    //c2
//    coord[0] =  1.0;
//    coord[1] =  1.0;
//    coord[2] = -1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in++);
//    //c3
//    coord[0] = -1.0;
//    coord[1] =  1.0;
//    coord[2] = -1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in++);
//    
//    //c4
//    coord[0] = -1.0;
//    coord[1] = -1.0;
//    coord[2] =  1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in++);
//    //c5
//    coord[0] =  1.0;
//    coord[1] = -1.0;
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
//    coord[0] = -1.0;
//    coord[1] =  1.0;
//    coord[2] =  1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in++);
    
    
    
    int index = 0;
    
    TPZVec<int64_t> TopologyQuad(4);
    
    // bottom
    TopologyQuad[0]=0;
    TopologyQuad[1]=1;
    TopologyQuad[2]=2;
    TopologyQuad[3]=3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index++,TopologyQuad,idf0,*gmesh);
    
    // Front
    TopologyQuad[0]=0;
    TopologyQuad[1]=1;
    TopologyQuad[2]=5;
    TopologyQuad[3]=4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index++,TopologyQuad,idf1,*gmesh);
    
    // Rigth
    TopologyQuad[0]=1;
    TopologyQuad[1]=2;
    TopologyQuad[2]=6;
    TopologyQuad[3]=5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index++,TopologyQuad,idf2,*gmesh);
    
    // Back
    TopologyQuad[0]=3;
    TopologyQuad[1]=2;
    TopologyQuad[2]=6;
    TopologyQuad[3]=7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index++,TopologyQuad,idf3,*gmesh);
    
    // Left
    TopologyQuad[0]=0;
    TopologyQuad[1]=3;
    TopologyQuad[2]=7;
    TopologyQuad[3]=4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index++,TopologyQuad,idf4,*gmesh);
    
    // Top
    TopologyQuad[0]=4;
    TopologyQuad[1]=5;
    TopologyQuad[2]=6;
    TopologyQuad[3]=7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index++,TopologyQuad,idf5,*gmesh);
    
    TPZManVector<int64_t,8> TopolCubo(8,0);
    TopolCubo[0] = 0;
    TopolCubo[1] = 1;
    TopolCubo[2] = 2;
    TopolCubo[3] = 3;
    TopolCubo[4] = 4;
    TopolCubo[5] = 5;
    TopolCubo[6] = 6;
    TopolCubo[7] = 7;
    
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index++, TopolCubo, matId, *gmesh);
    
    
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

void UniformRefine2dpol(TPZGeoMesh* gmesh, int nDiv)
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


TPZCompMesh *CMeshFlux2dpol(TPZGeoMesh *gmesh, int pOrder, int dim)
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
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    //    TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    //    TPZMaterial * BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    //    cmesh->InsertMaterialObject(BCond4);
    //    cmesh->InsertMaterialObject(BCond5);
    
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

TPZCompMesh *CMeshPressure2dpol(TPZGeoMesh *gmesh, int pOrder, int dim)
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
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    //    TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    //    TPZMaterial * BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    //    cmesh->InsertMaterialObject(BCond4);
    //    cmesh->InsertMaterialObject(BCond5);
    
    cmesh->SetDefaultOrder(pOrder);
    //cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
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

TPZCompMesh *CMeshMixed2dpol(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
{
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim = gmesh->Dimension();
    bool intface;
    TPZMatPoissonD3 *material = new TPZMatPoissonD3(matId,dim); intface = true; // nesse material tem que ser true
    //TPZMixedPoisson *material = new TPZMixedPoisson(matId,dim); intface = false; // nesse material tem que ser false
    
    //incluindo os dados do problema
    //    if (!intface) {
    //        TPZFNMatrix<2,REAL> PermTensor(dim,dim,0.);
    //        TPZFNMatrix<2,REAL> InvPermTensor(dim,dim,0.);
    //
    //        for (int i=0; i<dim; i++)
    //        {
    //            PermTensor(i,i) = 1.0;
    //        }
    //        InvPermTensor=PermTensor;
    //        material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    //    }
    
    //incluindo os dados do problema
    TPZFNMatrix<2,REAL> PermTensor(dim,dim,0.);
    TPZFNMatrix<2,REAL> InvPermTensor(dim,dim,0.);
    
    PermTensor = TP;
    InvPermTensor = InvTP;
    
    
    material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    
    //solucao exata
    
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    
    solexata = new TPZDummyFunction<STATE>(SolExata2dpol, 5);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing2dpol, 5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(10);
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
//    TPZMaterial * BCond4;
//    TPZMaterial * BCond5;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    //#ifdef PROBARCTAN
//    TPZAutoPointer<TPZFunction<STATE> > FBCond0 = new TPZDummyFunction<STATE>(ForcingBC0D2dpol, 5);
//    BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
//    BCond0->SetForcingFunction(FBCond0);
//    //#else
    TPZAutoPointer<TPZFunction<STATE> > FBCond0 = new TPZDummyFunction<STATE>(ForcingBC0N2dpol, 5);
    BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    BCond0->SetForcingFunction(FBCond0);
    //#endif
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1D2dpol, 5);
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2D2dpol, 5);
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    BCond2->SetForcingFunction(FBCond2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3D2dpol, 5);
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    BCond3->SetForcingFunction(FBCond3);
    
    //    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4D2dpol, 5);
    //    BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    //    BCond4->SetForcingFunction(FBCond4);
    //
    //    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    //#ifdef PROBARCTAN
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond5 = new TPZDummyFunction<STATE>(ForcingBC5D2dpol, 5); //diriclet
    //    BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2); //dirichlet
    //    BCond5->SetForcingFunction(FBCond5);
    //#else
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond5 = new TPZDummyFunction<STATE>(ForcingBC5N2dpol, 5);
    //    BCond5 = material->CreateBC(mat, bc5,neumann, val1, val2);
    //    BCond5->SetForcingFunction(FBCond5);
    //#endif
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    //    mphysics->InsertMaterialObject(BCond4);
    //    mphysics->InsertMaterialObject(BCond5);
    
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

void SolveSyst2dpol(TPZAnalysis &an, TPZCompMesh *fCmesh)
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
    cout <<"Numero de equacoes "<< fCmesh->NEquations()<< endl;
    
    //Saida de Dados: solucao e  grafico no VT
    //	ofstream file("Solutout");
    //	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void PosProcess2dpol(TPZAnalysis &an, std::string plotfile){
    TPZManVector<std::string,10> scalnames(1), vecnames(0);
    scalnames[0]= "Solution";
    //const int dim = 3;
    int div = 2;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    //    std::ofstream out("malhaNormal.txt");
    //    an.Print("nothing",out);
}

void PosProcessMultphysics2dpol(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZManVector<std::string,10> scalnames(2), vecnames(2);
    vecnames[0]  = "Flux";
    vecnames[1]  = "ExactFlux";
    //    vecnames[1]  = "GradFluxX";
    //    vecnames[2]  = "GradFluxY";
    //    vecnames[2]  = "GradFluxz";
    scalnames[0] = "Pressure";
    scalnames[1] = "ExactPressure";
    
    //const int dim = 3;
    int div = 2;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    std::ofstream out("malha.txt");
    an.Print("nothing",out);
    
    //    mphysics->Solution().Print("Solucao");
    
}

void SolExata2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    solp[0]=0.;
    flux.Resize(dim, 1);
    

    solp.Resize(1, 0.);
    flux.Resize(dim, 1.);
    double x = pt[0];
    double y = pt[1];
    for(int d=0; d<dim;d++) flux(d,0)=0.;
    solp[0] = (-1.0+x)*(1.0+x)*(-1.0+y)*(1.0+y);
    flux(0,0)= -2.0*x*(-1.0+y*y);
    flux(1,0)= -2.0*y*(-1.0+x*x);
    
}

void Forcing2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    // TP e o tensor de pemearbilidade.
    double x = pt[0];
    double y = pt[1];
    
    disp[0] = -2.0*(-2.0+x*x+y*y);
    
}

void ForcingBC0D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    double y = pt[1];
    disp[0] = (-1.0+x)*(1.0+x)*(-1.0+y)*(1.0+y);

}

void ForcingBC1D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];

    disp[0] = (-1.0+x)*(1.0+x)*(-1.0+y)*(1.0+y);;

}

void ForcingBC2D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];

    
    disp[0] = (-1.0+x)*(1.0+x)*(-1.0+y)*(1.0+y);

}

void ForcingBC3D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    
    disp[0] = (-1.0+x)*(1.0+x)*(-1.0+y)*(1.0+y);
    
}

void ForcingBC4D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];

    disp[0] = (-1.0+x)*(1.0+x)*(-1.0+y)*(1.0+y);
    
}

void ForcingBC5D2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];

    disp[0] = (-1.0+x)*(1.0+x)*(-1.0+y)*(1.0+y);

}


void ForcingBC0N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    
    disp[0] = 2.0 - 2.0*x*x;  
    
}

void ForcingBC1N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double y = pt[2];

    disp[0] = 2.0-2.0*y*y;
    
}

void ForcingBC2N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[1];
    double y = pt[2];
    disp[0] =  (  exp((1.+aa+bb+cc+x+y)/Epsilon)*(TP(1,0)+TP(1,1)+TP(1,2))  )/( (exp((2.*(aa+bb+cc))/Epsilon)+exp((2.*(1.+x+y))/Epsilon))*Epsilon );
}

void ForcingBC3N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[2];
    disp[0] =  (  exp((1.+aa+bb+cc+x+y)/Epsilon)*(TP(0,0)+TP(0,1)+TP(0,2))  )/( (exp((2.*(aa+bb+cc))/Epsilon)+exp((2.*(1.+x+y))/Epsilon))*Epsilon );
}

void ForcingBC4N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    disp[0] = 0.;
}

void ForcingBC5N2dpol(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    disp[0] = 0.;
}

void ErrorHDiv2dpol(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue; // Filtering lower dimension elements
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata2dpol, elerror, 0);
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

void ErrorL22dpol(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv)
{
    int64_t nel = l2mesh->NElements();
    //int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata2dpol, elerror, 0);
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
    
    fDebugMapL2.insert(std::pair<REAL, REAL> (ndiv,sqrt(globalerrors[1])));
    
    
}


void PrintDebugMapForMathematica2dpol(std::string filenameHdiv, std::string filenameL2)
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

TPZGeoMesh *GMeshDeformed2dpol(){
    
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
