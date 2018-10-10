//
//  main-ZigZagMesh.cpp
//  PZ
//
//  Created by Agnaldo Farias on 24/01/17.
//
//

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "TPZGeoCube.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"
#include "TPZVTKGeoMesh.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "TPZLagrangeMultiplier.h"
#include "pzmatmixedpoisson3d.h"

#include "mixedpoisson.h"

#include "pzgengrid.h"

#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <cmath>
#include <set>

#include "TPZParFrontStructMatrix.h"

// To compute jump of pressure over all faces - Jorge 2018
#include "ErrorOnFaces.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
//#endif

using namespace std;

const int matId = 1;

const int bcdirichlet = 0;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;
int const bc4 = -5;
int const bc5 = -6;


bool fTriang = true;

TPZGeoMesh *MalhaGeom(int nelx, int nely, REAL Lx, REAL Ly,bool ftriang, bool zigzag);
TPZGeoMesh *CreateOneCubo(int ndiv);
TPZGeoMesh * CreateOneCuboWithTetraedrons(int ndiv);
void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem);
bool MyDoubleComparer(REAL a, REAL b);


TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh, bool secondIntegration);

TPZCompMesh *MalhaCompH1(TPZGeoMesh * gmesh,int ordem);


void UniformRefine2(TPZGeoMesh* gmesh, int nDiv);

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh,int numthreads);
void SaidaSolucao(TPZAnalysis &an, std::string plotfile);

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv);
void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out, STATE &errorL2, STATE &errordu);

void SolProblema(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &flux);
void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &res);
void ForcingTang3(const TPZVec<REAL> &pt, TPZVec<STATE> &res,TPZFMatrix<STATE> &disp);
void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);


void ChangeSideConnectOrderConnects(TPZCompMesh *mesh, int sideorder);
void ChangeInternalOrderH1(TPZCompMesh *mesh, int neworder);
void NEquationsCondensed(TPZCompMesh *cmesh, int64_t &neqglob,int64_t &neqcond, bool ismisto);
void ChangeInternalConnectOrder(TPZCompMesh *mesh, int addToOrder);

void PermeabilityTensor(const TPZVec<REAL> &pt,TPZVec<STATE> &kabs,TPZFMatrix<STATE> &tensorK);
void ReactionTerm(const TPZVec<REAL> &pt,TPZVec<STATE> &alpha,TPZFMatrix<STATE> &disp);

bool metodomisto;
bool solArcTan;
bool problema3D;

int main()
{
    HDivPiola = 0;
    bool SecondIntegration = false;
    
    //InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    TPZVec<REAL> erros;
    
    metodomisto = true; //false --> formulacao H1
    solArcTan = true; //false --> ploblema 2D tese Margui
    problema3D = false;
    bool zigzag = false;
    bool HDivMaisMais = false;
    int nmais = 0;

    
    ofstream saidaerrosHdiv("../Erro-Misto.txt");
    ofstream saidaerrosH1("../Erro-H1.txt");
    
    int maxp = 11;
    int maxhref = 7;
    TPZFMatrix<STATE> L2ErrorPrimal(maxhref,maxp-1,0.);
    TPZFMatrix<STATE> L2ErrorDual(maxhref,maxp-1,0.);
    TPZFMatrix<STATE> L2ErrorDiv(maxhref,maxp-1,0.);
    TPZFMatrix<STATE> HDivErrorDual(maxhref,maxp-1,0.);
    TPZFMatrix<STATE> JumpPressure(maxhref,maxp-1,0.);
    
    TPZFMatrix<STATE> L2ConvergPrimal(maxhref-1,maxp-1,0.);
    TPZFMatrix<STATE> L2ConvergDual(maxhref-1,maxp-1,0.);
    TPZFMatrix<STATE> L2ConvergDiv(maxhref-1,maxp-1,0.);
    TPZFMatrix<STATE> HDivConverg(maxhref-1,maxp-1,0.);
    
    TPZFMatrix<int> porders(maxhref,maxp-1,0.);
    TPZFMatrix<int> numhref(maxhref,maxp-1,0.);
    TPZFMatrix<int> DofTotal(maxhref,maxp-1,0.);
    TPZFMatrix<int> DofCond(maxhref,maxp-1,0.);
    
    int nelx = 1;
    int nely = 1;
    
    TPZGeoMesh * gmesh;
    for(int p = 1; p<maxp; p++)
    {
        int pq = p;
        int pp = p;
        if(HDivMaisMais){
            pp = p + nmais;//Aqui = comeca com 1
        }

        //refinamentos h adptativos
        for(int nref = 0; nref<maxhref; nref++){
            
            if(metodomisto==true){
                saidaerrosHdiv<< "Saida do erro para formulacão hdiv, com ordem p = " << p << "\n";
                if (zigzag) {
                    saidaerrosHdiv << "Malha irregular - Trapezoidal\n";
                }
            }
            else{
                saidaerrosH1<< "Saida do erro para formulacão H1, com ordem p  = " << p << "\n";
                if (zigzag) {
                    saidaerrosH1 << "Malha irregular - Trapezoidal\n";
                }
            }
            
            nelx = 1 << nref;
            nely = 1 << nref;
            if(problema3D)
            {
                //gmesh = CreateOneCubo(nref);
                gmesh = CreateOneCuboWithTetraedrons(nref);
                //UniformRefine2(gmesh, nref);
            }
            else
            {
                gmesh = MalhaGeom(nelx, nely, 1,1,fTriang,zigzag);
            }
           
            {
//                ofstream filemesh1("malhageometrica.txt");
//                gmesh->Print(filemesh1);
//                std::ofstream filemesh("MalhaGeometricaInicial.vtk");
//                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
            }
            
            //rodar formulacao mista
            if(metodomisto==true){
                
                // Criando a primeira malha computacional
                TPZCompMesh * cmesh1= CMeshFlux(pq, gmesh);
//                {
//                    ofstream filemesh2("MalhaFlux.txt");
//                    filemesh2<<"\nDOF HDiv: "<< cmesh1->NEquations()<<std::endl;
//                    cmesh1->Print(filemesh2);
//                }
                
                if(HDivMaisMais)
                {
                    ChangeInternalConnectOrder(cmesh1, nmais);
//                    ofstream filemesh2("MalhaFlux.txt");
//                    filemesh2<<"\nDOF HDiv++: "<< cmesh1->NEquations()<<std::endl;
//                    cmesh1->Print(filemesh2);
                }
                
                // Criando a segunda malha computacional
                TPZCompMesh * cmesh2 = CMeshPressure(pp, gmesh);
                gmesh->ResetReference();
                cmesh2->LoadReferences();

//                {
//                    ofstream filemesh3("MalhaPressao.txt");
//                    filemesh3<<"\nDOF Pressao: "<< cmesh2->NEquations()<<std::endl;
//                    cmesh2->Print(filemesh3);
//                }
                
                // Identifying all the boundary and inner faces on triangulation   - Jorge 2018
     //           TPZStack<TPZCompElSide> faces;
       //         TPZStack<TPZCompElSide> facefromneigh;
         //       IdentifyingFaces(cmesh2,faces,facefromneigh);

                
                // Criando a malha computacional multifísica
                //malha multifisica
                TPZManVector<TPZCompMesh *,2> meshvec(2);
                meshvec[0] = cmesh1;
                meshvec[1] = cmesh2;
                TPZCompMesh * mphysics = MalhaCompMultifisica(meshvec, gmesh, SecondIntegration);
                
                int nDofTotal;
                nDofTotal = meshvec[0]->NEquations() + meshvec[1]->NEquations();
                
                int nDofCondensed;
                nDofCondensed = mphysics->NEquations();
                
//                ofstream filemesh4("MalhaMultifisica.txt");
//                filemesh4<<"\nDOF Total Multifisica: "<< nDofTotal<<std::endl;
//                filemesh4<<"DOF Condensados Multifisica: "<< nDofCondensed <<std::endl;
//                mphysics->Print(filemesh4);
                
                saidaerrosHdiv<< "NRefinamento h  = "<< nref <<std::endl;
                saidaerrosHdiv<< "Grau de Liberdade Total = " << nDofTotal<<std::endl;
                saidaerrosHdiv<< "Grau de Liberdade Condensado = " << nDofCondensed<<std::endl;
                
                // Resolvendo o sistema linear
                TPZAnalysis an(mphysics);
                ResolverSistema(an, mphysics,6);
                
                //            ofstream filemesh4("MalhaMultifisica.txt");
                //            mphysics->Print(filemesh4);
                
                TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                
                //Arquivo de saida para plotar a solução
                if(p==1)
                {
                    string plotfile("../Solution_Hdiv.vtk");
                    SaidaSolucao(an, plotfile);
                }
                
                STATE errorPrimalL2;
                STATE errorDuL2;
                TPZVec<STATE> errorsHDiv;
                STATE JumpAsError = 0.;
                
                saidaerrosHdiv<<"Erro da simulacao multifisica  para o Fluxo\n";
                ErrorHDiv2(cmesh1,saidaerrosHdiv,errorsHDiv);
                
                saidaerrosHdiv<<"Erro da simulacao multifisica  para a Pressao\n";
                ErrorH1(cmesh2, saidaerrosHdiv,errorPrimalL2,errorDuL2);
                saidaerrosHdiv<<"Salto de pressao como Erro\n";
                if(ComputePressureJumpOnFaces(cmesh2,matId, JumpAsError)) {
                    saidaerrosHdiv << "Jump of pressure = "    << JumpAsError << std::endl;
                    saidaerrosHdiv << std::endl;
                }
                else
                    saidaerrosHdiv << "Jump of pressure couldn't to be computed."    << std::endl << std::endl;

                
                L2ErrorPrimal(nref,p-1) = errorPrimalL2;
                L2ErrorDual(nref,p-1) = errorsHDiv[0];
                L2ErrorDiv(nref,p-1) = errorsHDiv[1];
                HDivErrorDual(nref,p-1) = errorsHDiv[2];
                
                JumpPressure(nref,p-1) = JumpAsError;   /// Jorge
                
                porders(nref,p-1) = p;
                numhref(nref,p-1) = nref;
                DofTotal(nref,p-1) = nDofTotal;
                DofCond(nref,p-1) = nDofCondensed;
            }
            else{
                
                TPZCompMesh * cmesh= MalhaCompH1(gmesh,p);
                
                ofstream arg("../MalhaCompH1.txt");
                cmesh->Print(arg);
                
                //                ChangeInternalOrderH1(cmesh, p+4);
                //            ofstream arg2("MalhaCompH1-2.txt");
                //            cmesh->Print(arg2);
                
                saidaerrosH1<< "\nNRefinamento h  = "<< nref <<std::endl;
                
                int64_t neq, neqcond;
                NEquationsCondensed(cmesh, neq, neqcond, metodomisto);
                
                saidaerrosH1 << "Numero total de equacoes: " <<neq<< "\n";
                saidaerrosH1 << "Numero de equacoes condensaveis: " <<neqcond<< "\n";
                saidaerrosH1 << "Numero de equacoes final: " << neq - neqcond<< "\n\n";
                
                
                // Resolvendo o sistema linear
                TPZAnalysis an(cmesh);
                ResolverSistema(an, cmesh,8);
                
                //Arquivo de saida para plotar a solução
                //                string plotfile("Solution_H1.vtk");
                //                SaidaSolucao(an, plotfile);
                
                //Pos-processamento calculo do erro
                an.SetExact(*SolProblema);
                bool store_errors = false;
                an.PostProcessError(erros, store_errors, saidaerrosH1);
                saidaerrosH1<<"==========================\n\n";
                
            }
        }
    }
    
    if(metodomisto){
        for(int j=0; j<maxp-1; j++){
            for(int i=0; i<maxhref-1; i++){
                L2ConvergPrimal(i,j) = log(L2ErrorPrimal(i+1,j)/L2ErrorPrimal(i,j))/log(1./2.);
                L2ConvergDual(i,j) = log(L2ErrorDual(i+1,j)/L2ErrorDual(i,j))/log(1./2.);
                L2ConvergDiv(i,j) = log(L2ErrorDiv(i+1,j)/L2ErrorDiv(i,j))/log(1./2.);
                HDivConverg(i,j) = log(HDivErrorDual(i+1,j)/HDivErrorDual(i,j))/log(1./2.);
            }
        }
        
        std::ofstream errtable("../ErrosConvergencia.txt");
        if(HDivPiola==0){
            errtable << "\nSolucão HDiv com transformacao e normais Normalizadas";
        }else if(HDivPiola==1){
            errtable << "\nSolucão HDiv com transformada de Piola";
        }else{
             errtable << "\nSolucão HDiv com vetores constantes";
        }
        if(zigzag){
            errtable << ": MALHA TRAPEZOIDAL";
        }else{
            errtable << ": MALHA REGULAR";
        }
        errtable <<"\n\n";
        
        L2ErrorPrimal.Print("Erro na variavel Primal: norma L2 = ",errtable);
        L2ConvergPrimal.Print("Convergencia Primal: norma L2 = ",errtable);
        
        L2ErrorDual.Print("Erro na variavel Dual: norma L2 = ",errtable);
        L2ConvergDual.Print("Convergencia Dual: norma L2 = ",errtable);
        
        L2ErrorDiv.Print("Erro no Divergente: norma L2 = ",errtable);
        L2ConvergDiv.Print("Convergencia Divergente: norma L2 = ",errtable);
        
        HDivErrorDual.Print("Erro na variavel Dual: norma HDiv = ",errtable);
        HDivConverg.Print("Convergencia Dual: norma HDiv  = ",errtable);
        
        JumpPressure.Print("Jump of Pressure: = ",errtable);   /// Jorge
        
        porders.Print("porder = ",errtable);
        numhref.Print("numhref = ",errtable);
        DofTotal.Print("DOF Total = ",errtable);
        DofCond.Print("DOF Condesed = ",errtable);
    }
    std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}



TPZGeoMesh *MalhaGeom(int nelx, int nely, REAL Lx, REAL Ly,bool ftriang, bool zigzag){
    
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    TPZManVector<int,3> nx(2);
    nx[0] = nelx;
    nx[1] = nely;
    x1[0] = Lx;
    x1[1] = Ly;
    TPZGenGrid gengrid(nx,x0,x1);
    
    if (ftriang && zigzag) {
        std::cout << "Zigzag meshes cannot be created with triangular meshes\n";
        DebugStop();
    }
    
    if (ftriang) {
        gengrid.SetElementType(ETriangle);
    }
    if (zigzag) {
        gengrid.SetDistortion(0.25);
        //        gengrid.SetZigZagPattern();
    }
//    gengrid.SetDistortion(0.75);
    gengrid.Read(gmesh);
    x1[0] = Lx;
    x1[1] = 0.;
    gengrid.SetBC(gmesh, x0, x1, bc1);
    x0 = x1;
    x1[0] = Lx;
    x1[1] = Ly;
    gengrid.SetBC(gmesh, x0, x1, bc2);
    x0 = x1;
    x1[0] = 0.;
    x1[1] = Ly;
    gengrid.SetBC(gmesh, x0, x1, bc3);
    x0 = x1;
    x1[0] = 0.;
    x1[1] = 0.;
    gengrid.SetBC(gmesh, x0, x1, bc4);
    
    
//    REAL theta = 0.0;
//    // It represents a 3D rotation around the z axis.
//    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
//    RotationMatrix(0,0) =   +cos(theta);
//    RotationMatrix(0,1) =   -sin(theta);
//    RotationMatrix(1,0) =   +sin(theta);
//    RotationMatrix(1,1) =   +cos(theta);
//    RotationMatrix(2,2) = 1.0;
//    TPZVec<STATE> iCoords(3,0.0);
//    TPZVec<STATE> iCoordsRotated(3,0.0);
//    
//    int NumberofGeoNodes = gmesh->NNodes();
//    for (int inode = 0; inode < NumberofGeoNodes; inode++)
//    {
//        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
//        GeoNode.GetCoordinates(iCoords);
//        // Apply rotation
//        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
//        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
//        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
//        GeoNode.SetCoord(iCoordsRotated);
//        gmesh->NodeVec()[inode] = GeoNode;
//    }
    
/*
    int Qnodes = 4;
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
    
    
    if(ftriang==true)
    {
        gmesh->SetMaxNodeId(Qnodes);
        gmesh->NodeVec().Resize(Qnodes+1);
        Node.Resize(Qnodes+1);
        
        //No central
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0, Lx/2.);//coord X
        Node[id].SetCoord(1, Ly/2.);//coord Y
        gmesh->NodeVec()[id] = Node[id];
        
        //        id = 0;
        //        TopolTriang[0] = 0;
        //        TopolTriang[1] = 1;
        //        TopolTriang[2] = 3;
        //        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
        //        id++;
        //
        //        TopolTriang[0] = 2;
        //        TopolTriang[1] = 1;
        //        TopolTriang[2] = 3;
        //        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        //        id++;
        
        id = 0;
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 1;
        TopolTriang[1] = 2;
        TopolTriang[2] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 3;
        TopolTriang[2] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        
        TopolTriang[0] = 3;
        TopolTriang[1] = 0;
        TopolTriang[2] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        //        TopolLine[0] = 2;
        //        TopolLine[1] = 1;
        //        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        //        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        //        TopolLine[0] = 3;
        //        TopolLine[1] = 2;
        //        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        //        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
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
 
*/
    
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout<<"Malha Geometrica Inicial\n";
//        gmesh->Print(sout);
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
    
    return gmesh;
}


TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
    int dim =gmesh->Dimension();
    
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(matId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    //cmesh->SetAllCreateFunctionsHDivFull();
    
    cmesh->InsertMaterialObject(mat);
    
    ///Criar condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1, bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2, bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3, bcdirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bc4, bcdirichlet, val1, val2);
    
    ///Inserir condicoes de contorno
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    if(problema3D)
    {
        //bc0
        TPZMaterial * BCond4 = material->CreateBC(mat, bc0, bcdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond4);
        
        //bc5
        TPZMaterial * BCond5 = material->CreateBC(mat, bc5, bcdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond5);
    }

    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    //    TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(matskeleton, dim-1, 1);
    //    TPZMaterial * mat2(matskelet);
    //    cmesh->InsertMaterialObject(mat2);
    
    cmesh->AutoBuild();//Ajuste da estrutura de dados computacional
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
    int dim =gmesh->Dimension();
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(matId,dim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    
    //        TPZMaterial * BCond0 = material->CreateBC(mat, bc0,bcdirichlet, val1, val2);
    //        TPZMaterial * BCond1 = material->CreateBC(mat, bc1,bcdirichlet, val1, val2);
    //        TPZMaterial * BCond2 = material->CreateBC(mat, bc2,bcdirichlet, val1, val2);
    //        TPZMaterial * BCond3 = material->CreateBC(mat, bc3,bcdirichlet, val1, val2);
    //        //para malha com 4 elementos usar mais estes contornos
    //        //    TPZMaterial *BCond4 = mat->CreateBC (mat,-5,0,val1,val2);
    //        //    TPZMaterial *BCond5 = mat->CreateBC (mat,-6,0,val1,val2);
    //        //    TPZMaterial *BCond6 = mat->CreateBC (mat,-7,0,val1,val2);
    //        //    TPZMaterial *BCond7 = mat->CreateBC (mat,-8,0,val1,val2);
    //
    //
    //
    //
    //        ///Inserir condicoes de contorno
    //        cmesh->InsertMaterialObject(BCond0);
    //        cmesh->InsertMaterialObject(BCond1);
    //        cmesh->InsertMaterialObject(BCond2);
    //        cmesh->InsertMaterialObject(BCond3);
    //        //    cmesh->InsertMaterialObject(BCond4);
    //        //    cmesh->InsertMaterialObject(BCond5);
    //        //    cmesh->InsertMaterialObject(BCond6);
    //        //    cmesh->InsertMaterialObject(BCond7);
    
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    //cmesh->SetAllCreateFunctionsDiscontinuous();
    
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
        if(!celdisc) continue;
        celdisc->SetConstC(1.);
        celdisc->SetCenterPoint(0, 0.);
        celdisc->SetCenterPoint(1, 0.);
        celdisc->SetCenterPoint(2, 0.);
        celdisc->SetTrueUseQsiEta();
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(fTriang==true) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

TPZCompMesh *MalhaCompH1(TPZGeoMesh * gmesh,int ordem){
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =2;
    cmesh->SetDefaultOrder(ordem);
    cmesh->SetDimModel(dim);
    
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForcingTang3, 5);
    dum->SetPolynomialOrder(20);
    force = dum;
    material->SetForcingFunction(force);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolProblema, 5);
    material->SetForcingFunctionExact(solexata);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    cmesh->InsertMaterialObject(mat);
    
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,bcdirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,bcdirichlet, val1, val2);
    
    
    ///Inserir condicoes de contorno
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Fazendo auto build
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    return cmesh;
}

#include "pzintel.h"
TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh, bool secondIntegration){
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =gmesh->Dimension();
    
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(matId,dim);
    
    if(secondIntegration==true){
        material->UseSecondIntegrationByParts();
    }
    
    //incluindo os dados do problema
    //    REAL coefk = 1.;
    //    material->SetPermeability(coefk);
    REAL coefvisc = 1.;
    material->SetViscosity(coefvisc);
    
    //permeabilidade
    TPZFMatrix<REAL> Ktensor(dim,dim,0.);
    TPZFMatrix<REAL> InvK(dim,dim,0.);
    
    for(int i=0; i<dim; i++){
        Ktensor(i,i)=1.;
    }
    
    InvK=Ktensor;
    material->SetPermeabilityTensor(Ktensor,InvK);
    
    if(problema3D==false && solArcTan==false)
    {
        //permeabilidade
        TPZAutoPointer<TPZFunction<STATE> > tensorK;
        tensorK = new TPZDummyFunction<STATE>(PermeabilityTensor, 5);
        material->SetPermeabilityFunction(tensorK);
        
        //Termo de reacao: alpha*p
        material->SetReactionTerm(1.);
        TPZAutoPointer<TPZFunction<STATE> > reaction;
        reaction = new TPZDummyFunction<STATE>(ReactionTerm, 5);
        material->SetfReactionTermFunction(reaction);
    }
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolProblema, 5);
    material->SetForcingFunctionExact(solexata);
    
    int int_order = 10;
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForcingF, 5);
    dum->SetPolynomialOrder(int_order);
    force = dum;
    material->SetForcingFunction(force);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    mphysics->SetDimModel(dim);
    
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);

    if(problema3D)
    {
        //bc1
        TPZMaterial *bnd = mat->CreateBC(mat, bc1, 0, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc1exata;
        TPZDummyFunction<STATE> *dum1;
        dum1 = new TPZDummyFunction<STATE>(Dirichlet, 5);//NeumannBC1
        dum1->SetPolynomialOrder(int_order);
        bc1exata = dum1;
        bnd->SetForcingFunction(bc1exata);
        mphysics->InsertMaterialObject(bnd);
        
        //bc2
        bnd = mat->CreateBC (mat, bc2, 1, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc2exata;
        TPZDummyFunction<STATE> *dum2;
        dum2 = new TPZDummyFunction<STATE>(NeumannBC2, 5);//NeumannBC2
        dum2->SetPolynomialOrder(int_order);
        bc2exata = dum2;
        bnd->SetForcingFunction(bc2exata);
        mphysics->InsertMaterialObject(bnd);
        
        //bc3
        bnd = mat->CreateBC (mat, bc3, 0, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc3exata;
        TPZDummyFunction<STATE> *dum3;
        dum3 = new TPZDummyFunction<STATE>(Dirichlet, 5);
        dum3->SetPolynomialOrder(int_order);
        bc3exata = dum3;
        bnd->SetForcingFunction(bc3exata);
        mphysics->InsertMaterialObject(bnd);
        
        //bc4
        bnd = mat->CreateBC (mat, bc4, 0, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc4exata;
        TPZDummyFunction<STATE> *dum4;
        dum4 = new TPZDummyFunction<STATE>(Dirichlet, 5);
        dum4->SetPolynomialOrder(int_order);
        bc4exata = dum4;
        bnd->SetForcingFunction(bc4exata);
        mphysics->InsertMaterialObject(bnd);
        
        //bc0
        bnd = mat->CreateBC (mat, bc0, 0, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc0exata;
        TPZDummyFunction<STATE> *dum0;
        dum0 = new TPZDummyFunction<STATE>(Dirichlet, 5);
        dum0->SetPolynomialOrder(int_order);
        bc0exata = dum0;
        bnd->SetForcingFunction(bc0exata);
        mphysics->InsertMaterialObject(bnd);
        
        //bc5
        bnd = mat->CreateBC (mat, bc5, 1, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc5exata;
        TPZDummyFunction<STATE> *dum5;
        dum5 = new TPZDummyFunction<STATE>(NeumannAcima, 5);//NeumannAcima
        dum5->SetPolynomialOrder(int_order);
        bc5exata = dum5;
        bnd->SetForcingFunction(bc5exata);
        mphysics->InsertMaterialObject(bnd);
        
    }
    else if (solArcTan)
    {
        //bc1
        TPZMaterial *bnd = mat->CreateBC(mat, bc1, 0, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc1exata;
        TPZDummyFunction<STATE> *dum1;
        dum1 = new TPZDummyFunction<STATE>(Dirichlet, 5);//NeumannBC1
        dum1->SetPolynomialOrder(int_order);
        bc1exata = dum1;
        bnd->SetForcingFunction(bc1exata);
        mphysics->InsertMaterialObject(bnd);
        
        //bc2
        bnd = mat->CreateBC (mat, bc2, 1, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc2exata;
        TPZDummyFunction<STATE> *dum2;
        dum2 = new TPZDummyFunction<STATE>(NeumannBC2, 5);//NeumannBC2
        dum2->SetPolynomialOrder(int_order);
        bc2exata = dum2;
        bnd->SetForcingFunction(bc2exata);
        mphysics->InsertMaterialObject(bnd);
        
        //bc3
        bnd = mat->CreateBC (mat, bc3, 0, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc3exata;
        TPZDummyFunction<STATE> *dum3;
        dum3 = new TPZDummyFunction<STATE>(Dirichlet, 5);
        dum3->SetPolynomialOrder(int_order);
        bc3exata = dum3;
        bnd->SetForcingFunction(bc3exata);
        mphysics->InsertMaterialObject(bnd);
        
        //bc4
        bnd = mat->CreateBC (mat, bc4, 0, val1, val2);
        TPZAutoPointer<TPZFunction<STATE> > bc4exata;
        TPZDummyFunction<STATE> *dum4;
        dum4 = new TPZDummyFunction<STATE>(Dirichlet, 5);
        dum4->SetPolynomialOrder(int_order);
        bc4exata = dum4;
        bnd->SetForcingFunction(bc4exata);
        mphysics->InsertMaterialObject(bnd);
    }
    else
    {
        TPZMaterial * BCond0 = material->CreateBC(mat, bc1, bcdirichlet, val1, val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, bc2, bcdirichlet, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, bc3, bcdirichlet, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, bc4, bcdirichlet, val1, val2);
    
        ///Inserir condicoes de contorno
        mphysics->InsertMaterialObject(BCond0);
        mphysics->InsertMaterialObject(BCond1);
        mphysics->InsertMaterialObject(BCond2);
        mphysics->InsertMaterialObject(BCond3);
    }
          
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    
    if(material->IsUsedSecondIntegration()==false){
        // Creating multiphysic elements into mphysics computational mesh
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
        mphysics->Reference()->ResetReference();
        mphysics->LoadReferences();
        
        // create condensed elements
        // increase the NumElConnected of one pressure connects in order to prevent condensation
        mphysics->ComputeNodElCon();
        for (int64_t icel=0; icel < mphysics->NElements(); icel++) {
            TPZCompEl  * cel = mphysics->Element(icel);
            if(!cel) continue;
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    break;
                }
            }
            new TPZCondensedCompEl(cel);
        }
        
        
        mphysics->CleanUpUnconnectedNodes();
        mphysics->ExpandSolution();
    }
    
    //Creating multiphysic elements containing skeletal elements.
    else
    {
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        mphysics->Reference()->ResetReference();
        mphysics->LoadReferences();
        
        int64_t nel = mphysics->ElementVec().NElements();
        
        std::map<int64_t, int64_t> bctoel, eltowrap;
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = mphysics->Element(el);
            TPZGeoEl *gel = cel->Reference();
            int matid = gel->MaterialId();
            if (matid < 0) {
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->Dimension() == dim && neighbour.Element()->Reference()) {
                        // got you!!
                        bctoel[el] = neighbour.Element()->Reference()->Index();
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if (neighbour == gelside) {
                    DebugStop();
                }
            }
        }
        
        TPZStack< TPZStack< TPZMultiphysicsElement *,7> > wrapEl;
        for(int64_t el = 0; el < nel; el++)
        {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(mphysics->Element(el));
            if(mfcel->Dimension()==dim) TPZBuildMultiphysicsMesh::AddWrap(mfcel, matId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
        }
        
        for (int64_t el =0; el < wrapEl.size(); el++) {
            TPZCompEl *cel = wrapEl[el][0];
            int64_t index = cel->Index();
            eltowrap[index] = el;
        }
        
        meshvec[0]->CleanUpUnconnectedNodes();
        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
        std::map<int64_t, int64_t>::iterator it;
        for (it = bctoel.begin(); it != bctoel.end(); it++) {
            int64_t bcindex = it->first;
            int64_t elindex = it->second;
            if (eltowrap.find(elindex) == eltowrap.end()) {
                DebugStop();
            }
            int64_t wrapindex = eltowrap[elindex];
            TPZCompEl *bcel = mphysics->Element(bcindex);
            TPZMultiphysicsElement *bcmf = dynamic_cast<TPZMultiphysicsElement *>(bcel);
            if (!bcmf) {
                DebugStop();
            }
            wrapEl[wrapindex].Push(bcmf);
            
        }
        
        //------- Create and add group elements -------
        int64_t index, nenvel;
        nenvel = wrapEl.NElements();
        TPZStack<TPZElementGroup *> elgroups;
        for(int64_t ienv=0; ienv<nenvel; ienv++){
            TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
            elgroups.Push(elgr);
            nel = wrapEl[ienv].NElements();
            for(int jel=0; jel<nel; jel++){
                elgr->AddElement(wrapEl[ienv][jel]);
            }
        }
        
        mphysics->ComputeNodElCon();
        // create condensed elements
        // increase the NumElConnected of one pressure connects in order to prevent condensation
        for (int64_t ienv=0; ienv<nenvel; ienv++) {
            TPZElementGroup *elgr = elgroups[ienv];
            int nc = elgr->NConnects();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = elgr->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    break;
                }
            }
            TPZCondensedCompEl *condense = new TPZCondensedCompEl(elgr);
        }
    }
    
    mphysics->CleanUpUnconnectedNodes();
    mphysics->ExpandSolution();
    
    return mphysics;
}


#define VTK
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads)
{
   //TPZSkylineStructMatrix strmat(fCmesh);
    
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
    strmat.SetDecomposeType(ELDLt);
    
    
    if(numthreads>0){
        strmat.SetNumThreads(numthreads);
    }
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); //caso simetrico
    //	step.SetDirect(ELU);
    an.SetSolver(step);
    //    an.Assemble();
    an.Run();
    
    //Saida de Dados: solucao e  grafico no VTK
    //ofstream file("Solution.out");
    //an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}


void SaidaSolucao(TPZAnalysis &an, std::string plotfile){
    
    if(metodomisto==true)
    {
        TPZManVector<std::string,10> scalnames(5), vecnames(2);
        scalnames[0] = "Pressure";
        scalnames[1] = "ExactPressure";
        scalnames[2]="POrder";
        vecnames[0]= "Flux";
        vecnames[1]= "ExactFlux";
        scalnames[3] = "Divergence";
        scalnames[4] = "ExactDiv";
      
        const int dim = an.Mesh()->Dimension();

        int div = 0;
        an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
        an.PostProcess(div,dim);
    }
    
    //Rodar H1
    else{
        TPZManVector<std::string,10> scalnames(3), vecnames(1);
        scalnames[0] = "Solution";
        scalnames[1]="POrder";
        scalnames[2] = "ExactSolution";
        vecnames[0]= "Derivative";
        
        const int dim = an.Mesh()->Dimension();
        int div = 0;
        an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
        an.PostProcess(div,dim);
    }
}


void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    TPZStack<REAL> vech;
    
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        //--metodo para pegar tamanho da malha
        //        TPZMaterialData data;
        //        data.fNeedsHSize=true;
        //        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace * >(cel);
        //        TPZVec<REAL> qsi(2,0.);
        //        sp->ComputeRequiredData(data, qsi);
        //        REAL &hsize = data.HSize;
        //        vech.push_back(hsize);
        
        //----
        
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolProblema, elerror, 0);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    
    
    //    int nh = vech.size();
    //    REAL hmax=0;
    //    for(int i=0; i<nh; i++){
    //        if(vech[i]>hmax){
    //            hmax=vech[i];
    //        }
    //    }
    // out << "Errors associated with HDiv space\n";
    //out << " Hmax = "    << hmax << std::endl;
    
    
    //out << "Errors associated with HDiv space\n";
    out << "L2 Error Norm for flux = "    << sqrt(globerrors[1]) << endl;
    errorHDiv.Resize(3,0.);
    errorHDiv[0] = sqrt(globerrors[1]);
    errorHDiv[1] = sqrt(globerrors[2]);
    errorHDiv[2] = sqrt(globerrors[3]);
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<endl;
    
}

void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out, STATE &errorL2, STATE &errorDu)
{
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolProblema, elerror, 0);
        
        int nerr = elerror.size();
        globerrors.resize(nerr);
//#ifdef LOG4CXX
//        if (logger->isDebugEnabled()) {
//            std::stringstream sout;
//            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
//            LOGPZ_DEBUG(logger, sout.str())
//        }
//#endif
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
    }
    //out << "Errors associated with L2 or H1 space\n";
    out << "H1 Error Norm = "    << sqrt(globerrors[0]) << endl;
    out << "L2 Error Norm = "    << sqrt(globerrors[1]) << endl;
    out << "Semi H1 Norm = "    << sqrt(globerrors[2]) << endl;
    out << "=============================\n"<<endl;
    errorL2 = sqrt(globerrors[1]);
}


void UniformRefine2(TPZGeoMesh* gmesh, int nDiv)
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
    	gmesh->BuildConnectivity();
}

#define Power pow
#define ArcTan atan
#define Sqrt sqrt

void SolProblema(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &flux){
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    u.Resize(1, 0.);
    flux.Resize(3, 1);
    flux(0,0)=0., flux(1,0)=0., flux(2,0)=0.;
    
    if(problema3D)
    {
        flux.Resize(4, 1);
        flux(3,0)=0.;
        
        REAL x0 = 1.25, y0 = -0.25, z0 = -0.25;
        REAL r0 = M_PI/3.;
        REAL alpha = 5.;
        
        REAL temp1;
        REAL temp2;
        REAL r;
        REAL grad;
        temp1 = 0.,temp2=0., r=0., grad=0.;
        
        //solucao u
        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
        r = sqrt(temp1);
        temp2 = (r - r0)*alpha;
        u[0] = M_PI/2. - atan(temp2);
        
        //fluxo em x
        temp2 =  r*(1./alpha + (r - r0)*(r - r0)*alpha);
        temp1 = (x0-x);
        grad=temp1/temp2;
        if(metodomisto) grad *=-1.;
        flux(0,0) = grad;
        
        //fluxo em y
        temp1 = (y0-y);
        grad=temp1/temp2;
        if(metodomisto) grad *=-1.;
        flux(1,0) = grad;
        
        
        //fluxo em z
        temp1 = (z0-z);
        grad=temp1/temp2;
        if(metodomisto) grad *=-1.;
        flux(2,0) = grad;
        
        //Solucao do divergente de u
        REAL temp3=0., temp4 = 0., div=0.;
        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
        r = sqrt(temp1);
        temp2 = (2./alpha) + 2.*r0*(r0-r)*alpha;
        temp3 = r*((r - r0)*(r - r0) + 1./(alpha*alpha));
        temp4 = 1. + (r - r0)*(r - r0)*alpha*alpha;
        div = temp2/(temp3*temp4);
        flux(3,0) = div; //valor do divergente
    }
    else if(solArcTan)
    {
        REAL x0 = 1.25, y0 = -0.25;
        REAL r0 = M_PI/3.;
        REAL alpha = 5.;
        
        REAL temp1;
        REAL temp2;
        REAL r;
        REAL grad;
        temp1 = 0.,temp2=0., r=0., grad=0.;
        
        //solucao u
        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0);
        r = sqrt(temp1);
        temp2 = (r - r0)*alpha;
        u[0] = M_PI/2. - atan(temp2);
        
        //fluxo em x
        temp2 =  r*(1./alpha + (r - r0)*(r - r0)*alpha);
        temp1 = (x0-x);
        grad=temp1/temp2;
        if(metodomisto) grad *= -1.;
        flux(0,0) = grad;
        
        //fluxo em y
        temp1 = (y0-y);
        grad=temp1/temp2;
        if(metodomisto) grad *= -1.;
        flux(1,0) = grad;
        
        //Solucao do divergente de u
        REAL temp3=0., temp4 = 0., div=0.;
        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0);
        r = sqrt(temp1);
        temp2 = (1./alpha) + (r0*r0 - r*r)*alpha;
        temp3 = r*((r - r0)*(r - r0) + 1./(alpha*alpha));
        temp4 = 1. + (r - r0)*(r - r0)*alpha*alpha;
        div = temp2/(temp3*temp4);
        flux(2,0) = div; //valor do divergente

//------- Solucao sin(pix)sin(piy) -----------
//        u[0] = sin(M_PI*x)*sin(M_PI*y);
//        flux(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y);
//        flux(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x);
//        flux(2,0) = 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y); //valor do divergente
    }
    else
    {
        //Solucao u
        REAL solp = sin(M_PI*x)*sin(M_PI*y);
        u[0] = solp;
        
        REAL temp = 1. + 10.*x;
        
        //fluxo em x
        flux(0,0) = M_PI*temp*cos(M_PI*x)*sin(M_PI*y);
        if(metodomisto) flux(0,0) *=-1.;
        
        //fluxo em y
        flux(1,0) =  M_PI*temp*cos(M_PI*y)*sin(M_PI*x);
        if(metodomisto) flux(1,0) *=-1.;
        
        //divergente: -(dux/dx + duy/dy)
        flux(2,0) = -10.*M_PI*cos(M_PI*x)*sin(M_PI*y) + 2.*M_PI*M_PI*temp*sin(M_PI*x)*sin(M_PI*y);
    }
}


void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &res){
    
    
    double x = pt[0];
    double y = pt[1];
    res[0] = 0.;
    
    if(problema3D)
    {
        REAL z = pt[2];
        
        REAL x0 = 1.25, y0 = -0.25, z0 = -0.25;
        REAL r0 = M_PI/3.;
        REAL alpha = 5.;
        
        REAL temp1=0., temp2=0.,temp3=0., temp4=0., r=0., sol=0.;
        
        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
        r = sqrt(temp1);
        
        temp2 = (2./alpha) + 2.*r0*(r0-r)*alpha;
        temp3 = r*((r - r0)*(r - r0) + 1./(alpha*alpha));
        temp4 = 1. + (r - r0)*(r - r0)*alpha*alpha;
        
        sol = temp2/(temp3*temp4);
        
        res[0] = sol;
    }
    else if(solArcTan)
    {
        REAL x0 = 1.25, y0 = -0.25;
        REAL r0 = M_PI/3.;
        REAL alpha = 5.;
        
        REAL temp1=0., temp2=0.,temp3=0., temp4=0., r=0., sol=0.;
        
        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0);
        r = sqrt(temp1);
        
        temp2 = (1./alpha) + (r0*r0 - r*r)*alpha;
        temp3 = r*((r - r0)*(r - r0) + 1./(alpha*alpha));
        temp4 = 1. + (r - r0)*(r - r0)*alpha*alpha;
        
        sol = temp2/(temp3*temp4);
        res[0] = sol;

//------- Solucao sin(pix)sin(piy) -----------
//        res[0] = 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    }
    else
    {
        REAL solp =  sin(M_PI*x)*sin(M_PI*y);
        
        REAL temp1 = 1. + 10.*x;
        REAL temp2 = 1. - x*x - y*y;
        REAL temp3 = exp(temp2);

        res[0] = -10.*M_PI*cos(M_PI*x)*sin(M_PI*y) + (temp3 + 2.*M_PI*M_PI*temp1)*solp;
    }
}

void ForcingTang3(const TPZVec<REAL> &pt, TPZVec<STATE> &res,TPZFMatrix<STATE> &disp)
{
//    double x = pt[0];
//    double y = pt[1];
//    res[0] = 0.;
//    
//    if(solSinSin)
//    {
//        res[0] = 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
//    }
//    else
//    {
//        REAL solp =  sin(M_PI*x)*sin(M_PI*y);
//        
//        REAL temp1 = 1. + 10.*x;
//        REAL temp2 = 1. - x*x - y*y;
//        REAL temp3 = exp(temp2);
//        
//        res[0] = -10.*M_PI*cos(M_PI*x)*sin(M_PI*y) + (temp3 + 2.*M_PI*M_PI*temp1)*solp;
//    }
    
    ForcingF(pt, res);

    res[0] *=-1.;
}

void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result)
{
    TPZFMatrix<STATE> du(3,1);
    SolProblema(loc,result,du);
}

void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result)
{
    REAL normal[3] = {1.,0.,0.};
    TPZFNMatrix<5,STATE> du(3,1);
    
    SolProblema(loc,result,du);
    
    result.Resize(1);
    if(problema3D)
    {
        result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
    }
    else
    {
        result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
    }
}

void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result)
{
    REAL normal[3] = {0.,0.,1.};
//    TPZManVector<REAL> u(1);
    TPZFNMatrix<5,STATE> du(3,1);
    
    SolProblema(loc,result,du);
    
    result.Resize(1);
    if(problema3D)
    {
        result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
    }
    else
    {
        result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
    }
}


void PermeabilityTensor(const TPZVec<REAL> &pt,TPZVec<STATE> &kabs,TPZFMatrix<STATE> &tensorK)
{
    
    REAL x = pt[0];
    //REAL y = pt[1];
    tensorK.Resize(4,2);
    kabs.Resize(1, 0.);
    
    //K
    REAL temp = 1.+10.*x;
    tensorK(0,0) = temp;     tensorK(0,1) = 0.0;
    tensorK(1,0) = 0.0;      tensorK(1,1) = temp;
    
    //Kinv
    tensorK(2,0) = 1.0/temp;     tensorK(2,1) = 0.0;
    tensorK(3,0) = 0.0;          tensorK(3,1) = 1.0/temp;
}

void ReactionTerm(const TPZVec<REAL> &pt,TPZVec<STATE> &alpha,TPZFMatrix<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    disp.Resize(2,2);
    alpha.Resize(1, 0.);
    
    //Termo de reacao
    REAL temp = 1. - x*x - y*y;
    alpha = exp(temp);
}


//void ChangeExternalOrderConnects(TPZCompMesh *mesh){
//    
//    int nEl= mesh-> NElements();
//    int dim = mesh->Dimension();
//    
//    int cordermin = -1;
//    for (int iel=0; iel<nEl; iel++) {
//        TPZCompEl *cel = mesh->ElementVec()[iel];
//        if (!cel) continue;
//        int corder = 0;
//        int nshape = 0;
//        
//        TPZGeoEl *gel = cel->Reference();
//        int nsides = gel->NSides();
//        int nnodes = gel->NCornerNodes();
//        
//        if(cel->Dimension()== dim)
//        {
//            for(int is = nnodes; is<nsides-1; is++)
//            {
//                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
//                TPZConnect &conside = intel->MidSideConnect(is);
//                
//                corder = conside.Order();
//                if(corder!=cordermin)
//                {
//                    cordermin = corder-1;
//                    conside.SetOrder(cordermin);
//                    
//                    TPZGeoEl *gelface = gel->CreateBCGeoEl(is, -100);
//                    
//                    if(gelface->Type()==EQuadrilateral)
//                    {
//                        nshape = (cordermin+1)*(cordermin+1);
//                    }
//                    else if (gelface->Type()==ETriangle)
//                    {
//                        nshape = (cordermin+1)*(cordermin+2)/2;
//                    }
//                    else
//                    {
//                        //caso linera
//                        if(gelface->Type()!= EOned) DebugStop();
//                        nshape = cordermin+1;
//                    }
//                    delete gelface;
//                    
//                    conside.SetNShape(nshape);
//                    mesh->Block().Set(conside.SequenceNumber(), nshape);
//                }
//
//                
//                
////                TPZConnect &conel  = cel->Connect(icon);
////                corder = conel.Order();
////                nshape = conel.NShape();
////                if(corder!=cordermin)
////                {
////                    cordermin = corder-1;
////                    conel.SetOrder(cordermin);
////                    
////                                        conel.SetNShape(nshape-1);
////                    mesh->Block().Set(conel.SequenceNumber(),nshape-1);
////                }
//            }
//        }
//    }
//    mesh->ExpandSolution();
//    mesh->CleanUpUnconnectedNodes();
//}


void ChangeSideConnectOrderConnects(TPZCompMesh *mesh, int order){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        
        if(cel->Dimension()== dim)
        {
            for (int icon=0; icon<ncon-1; icon++)
            {
                TPZConnect &conel  = cel->Connect(icon);
               
                corder = conel.Order();
                nshape = conel.NShape();
                int64_t cindex = cel->ConnectIndex(icon);
                if(corder!=order)
                {
                    conel.SetOrder(order,cindex);
                    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                    nshape = intel->NConnectShapeF(icon,order);
                    conel.SetNShape(nshape);
                    mesh->Block().Set(conel.SequenceNumber(),nshape);
                }
            }
        }
    }
    mesh->ExpandSolution();
    mesh->CleanUpUnconnectedNodes();
}


void ChangeInternalOrderH1(TPZCompMesh *mesh, int neworder){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        
        int ninternalshape = 0;
        if(cel->Dimension()== dim)
        {
            TPZConnect &co  = cel->Connect(ncon-1);
            int64_t cindex = cel->ConnectIndex(ncon-1);
            ninternalshape = (neworder-1)*(neworder-1);
            co.SetOrder(neworder,cindex);
            co.SetNShape(ninternalshape);
            mesh->Block().Set(co.SequenceNumber(),ninternalshape);
        }
    }
    mesh->ExpandSolution();
    mesh->CleanUpUnconnectedNodes();
}



void NEquationsCondensed(TPZCompMesh *cmesh, int64_t &neqglob,int64_t &neqcond, bool ismisto){
    
    int64_t ncon = cmesh->NConnects();
    neqglob = 0;
    neqcond = 0;
    for(int i = 0; i< ncon; i++){
        TPZConnect &co  = cmesh->ConnectVec()[i];
        //if(co.HasDependency()) continue;
        if(co.HasDependency()  || co.IsCondensed() || !co.NElConnected() || co.SequenceNumber() == -1) continue;
        
        int nelc = co.NElConnected();
        if (nelc==0) DebugStop();
        
        int dofsize = co.NShape()*co.NState();
        neqglob += dofsize;
        
        //equacoes condensaveis
        if (nelc == 1){
            neqcond += dofsize;
        }
    }
    
    if(ismisto){
        int nel2D =0;
        for(int i=0; i<cmesh->NElements(); i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            if(!cel) continue;
            if(cel->Reference()->Dimension() == cmesh->Dimension()){
                nel2D++;
            }
        }
        neqcond = neqcond - nel2D;
    }
}

void ChangeInternalConnectOrder(TPZCompMesh *mesh, int addToOrder){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        int nshape2 = 0;
        
        if(cel->Dimension()== dim)
        {
            TPZConnect &conel = cel->Connect(ncon-1);
            corder = conel.Order();
            nshape = conel.NShape();
            
            int neworder = corder + addToOrder;//Aqui = +1
            int64_t cindex = cel->ConnectIndex(ncon-1);
            conel.SetOrder(neworder,cindex);
            
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            intel->SetPreferredOrder(neworder);
            nshape = intel->NConnectShapeF(ncon-1,neworder);
            
            if(dim==2 && addToOrder==1)
            {
                if(fTriang){
                    nshape2 = (corder + 2)*(corder + 2)-1;
                }else{//Quadrilateral
                    nshape2 = 2*(corder + 1)*(corder + 2);
                }
                if(nshape2!=nshape)
                {
                    DebugStop();
                }
            }
            
            conel.SetNShape(nshape);
            mesh->Block().Set(conel.SequenceNumber(),nshape);
        }
    }
    mesh->CleanUpUnconnectedNodes();
    mesh->ExpandSolution();
}


TPZGeoMesh *CreateOneCubo(int ndiv)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int nnodes = 8;
    
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    gmesh->SetMaxNodeId(8);
    
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    
    //cubo [0,1]ˆ3
    //c0
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c1
    coord[0] = 1.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c2
    coord[0] = 1.0;
    coord[1] = 1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c3
    coord[0] = 0.0;
    coord[1] = 1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    
    //c4
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c5
    coord[0] = 1.0;
    coord[1] = 0.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c6
    coord[0] = 1.0;
    coord[1] = 1.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c7
    coord[0] = 0.0;
    coord[1] = 1.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    gmesh->SetMaxElementId(in);
    
    int index = 0;
    
    TPZVec<int64_t> TopologyQuad(4);
    
    // bottom
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 1;
    TopologyQuad[2] = 2;
    TopologyQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc0,*gmesh);
    index++;
    
    // Front
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 1;
    TopologyQuad[2] = 5;
    TopologyQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc1,*gmesh);
    index++;
    
    // Rigth
    TopologyQuad[0] = 1;
    TopologyQuad[1] = 2;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc2,*gmesh);
    index++;
    
    // Back
    TopologyQuad[0] = 3;
    TopologyQuad[1] = 2;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc3,*gmesh);
    index++;
    
    // Left
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 3;
    TopologyQuad[2] = 7;
    TopologyQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc4,*gmesh);
    index++;
    
    // Top
    TopologyQuad[0] = 4;
    TopologyQuad[1] = 5;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc5,*gmesh);
    index++;
    
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
    index++;
    
    
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
    
    
    //    std::ofstream out("SingleCubeWithBcs.vtk");
    //    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZGeoMesh * CreateOneCuboWithTetraedrons(int ndiv)
{
    double dndiv = ndiv;
    int nelem = (int) pow(2., dndiv);
    
    int tetraedra_2[6][4] = { {1,2,5,4}, {4,7,3,2}, {0,1,2,4}, {0,2,3,4}, {4,5,6,2},{4,6,7,2} };
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh,nelem);
    gmesh->SetDimension(3);
    
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                
                for (int el=0; el<6; el++)
                {
                    TPZManVector<int64_t,4> elnodes(4);
                    int64_t index;
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra_2[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, matId, index);
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
        TPZVec<int64_t> ncoordVec(0);
        int64_t sizeOfVec = 0;
        
        // na face z = 0
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
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
            int64_t pos = tetra->NodeIndex(i);
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
            int64_t pos = tetra->NodeIndex(i);
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
            int64_t pos = tetra->NodeIndex(i);
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
            int64_t pos = tetra->NodeIndex(i);
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
            int64_t pos = tetra->NodeIndex(i);
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

void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem)
{
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (int64_t i=0; i<=nelem; i++) {
        for (int64_t j=0; j<=nelem; j++) {
            for (int64_t k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
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

