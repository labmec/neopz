//
//  main_mixed.cpp
//  PZ
//
//  Created by Agnaldo Farias on 3/05/15.
//  Copyright (c) 2015 LabMec-Unicamp. All rights reserved.
//

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzinterpolationspace.h"
#include "TPZCompElDisc.h"
#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzgengrid.h"

#include "pzlog.h"

#include "pzl2projection.h"
#include "pzmultiphysicselement.h"
#include "pzintel.h"
#include "TPZVTKGeoMesh.h"

#include "TPZMatDualHybridPoisson.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"
#include "pzintel.h"
#include "pzmatmixedpoisson3d.h"

#include "pzbstrmatrix.h"
#include "pzvisualmatrix.h"
#include "pzskylmat.h"

#include <iostream>
#include <math.h>

#include "TPZSSpStructMatrix.h"
#include "TPZCompMeshTools.h"
#include "TPZMatLaplacian.h"
#include "TPZFrontNonSym.h"
#include "TPZSkylineNSymStructMatrix.h"

using namespace std;
#include <algorithm>
#include <iterator>

const int matId = 1;

const int bcdirichlet = 0;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;
int const bc4 = -5;
int const bc5 = -6;

bool fTriang = false;

TPZGeoMesh *CreateOneCubo(int ndiv);
TPZGeoMesh * CreateOneCuboWithTetraedrons(int ndiv);
void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem);
bool MyDoubleComparer(REAL a, REAL b);


TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh, bool secondIntegration);
TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder, bool ismultiplierH1, REAL penaltyConst);
TPZCompMesh *MalhaCompH1(TPZGeoMesh * gmesh,int ordem, REAL penaltyConst);

void GroupElements(TPZCompMesh *cmesh, int dimproblema);


void UniformRefine2(TPZGeoMesh* gmesh, int nDiv);

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh,int numthreads);
void SaidaSolucao(TPZAnalysis &an, std::string plotfile);

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv);
void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out, TPZVec<STATE> &errorL2);

void SolProblema(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &flux);
void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &res);
void ForcingTang(const TPZVec<REAL> &pt, TPZVec<STATE> &res,TPZFMatrix<STATE> &disp);
void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannBC1(const TPZVec<REAL> &loc, TPZVec<STATE> &result);

void ChangeSideConnectOrderConnects(TPZCompMesh *mesh, int sideorder);
void ChangeInternalOrderH1(TPZCompMesh *mesh, int neworder);
void NEquationsCondensed(TPZCompMesh *cmesh, int64_t &neqglob,int64_t &neqcond, bool ismisto);
void ChangeInternalConnectOrder(TPZCompMesh *mesh, int addToOrder);


#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

int mainHybrid();
int mainH1();
int mainMixed();

bool hybrid = true;
bool metodomisto = false;//Rodar Metodo Misto
bool rodarSIPGD = false;//Rodar Metodo DG
bool rodarH1 = false;//Rodar Metodo H1

bool fTetra = false;//false: malha com Hexaedros;

int main()
{
    if(hybrid)
    {
        return mainHybrid();
    }
    if(metodomisto)
    {
        return mainMixed();
    }
    if(rodarSIPGD || rodarH1)
    {
        return mainH1();
    }
    return 0;
}

//Rodar Metodo Hibrido
int mainHybrid()
{
    //InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    ofstream saidaerrosHibrido("../Erro-Hibrido.txt");
    
    bool multiplicadorH1 = false;//true: Hibrido Continuo
    
    //bool MudarOrdemPdoMultiplicador = false;
    
    int maxp = 5;
    int ndiv = 0;

    std::ofstream myerrorfile("Resultado-Hibrido.txt");
    
    TPZGeoMesh * gmesh;
    TPZCompMesh * cmesh;
    for(int p = 4; p<maxp; p++)
    {
        myerrorfile<<"\nORDEM p = "<<p <<"\n\n";
        myerrorfile << "ndiv" << setw(10) <<"NDoF"<< setw(12)<<"NDoFCond" << "     Entradas" <<"       NumZeros" <<
        "       Razao" <<setw(19)<< "Assemble"<< setw(20)<<"Solve" << setw(20) <<"Ttotal" <<setw(12) <<"Error u" << setw(16)<<"Error gradU\n"<<std::endl;
        
        for(ndiv = 0; ndiv<5; ndiv++)
        {
            int64_t NDoF=0, NDoFCond=0;
            int64_t nNzeros=0;
            
            //Geo Mesh
            if(!fTetra){
                gmesh = CreateOneCubo(ndiv);
            }
            else{
                gmesh = CreateOneCuboWithTetraedrons(ndiv);
            }
            
            REAL Beta0 = 6.;
            cmesh = CreateHybridCompMesh(*gmesh, p, multiplicadorH1, Beta0);
            
//            ofstream out1("gmesh1.txt");
//            gmesh->Print(out1);
//            
//            ofstream out11("cmesh1.txt");
//            cmesh->Print(out11);
            
            //------- Criar elementos de Lagrange (Ribs)--------
            //materiais do problema
            std::set<int>matids;
            matids.insert(matId);
            int dim = gmesh->Dimension();
            int nbc = dim*2;
            for(int i=1; i<=nbc; i++)
            {
                if(dim==2){
                    matids.insert(-(i+1));
                }else{
                    matids.insert(-i);
                }
            }
            
            cmesh->ApproxSpace().Hybridize(*cmesh, matids, multiplicadorH1);
            
//            ofstream out2("gmesh2.txt");
//            gmesh->Print(out2);
//            
//            ofstream out22("cmesh2.txt");
//            cmesh->Print(out22);

            
            NDoF = cmesh->NEquations();
            
            //condesacao estatica
            GroupElements(cmesh, dim);
            cmesh->LoadReferences();//mapeia para a malha geometrica lo
//            ofstream out3("gmesh3.txt");
//            gmesh->Print(out3);
//            
//            ofstream out33("cmesh3.txt");
//            cmesh->Print(out33);
            
            NDoFCond = cmesh->NEquations();

            
            // Resolvendo o sistema linear
            TPZAnalysis analysis(cmesh);
            int64_t neq = NDoFCond;
            TPZVec<int64_t> skyline;
            cmesh->Skyline(skyline);
            TPZSkylMatrix<STATE> matsky(neq,skyline);
            nNzeros = matsky.GetNelemts();
            
            //MKL Pardiso
            TPZSymetricSpStructMatrix strmat(cmesh);
            int nthreads = 12;
            strmat.SetNumThreads(nthreads);
            analysis.SetStructuralMatrix(strmat);
            
//            TPZSkylineStructMatrix  strmat(cmesh);
//            strmat.SetNumThreads(6);
//            analysis.SetStructuralMatrix(strmat);
            
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt); //caso simetrico
            analysis.SetSolver(step);
            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            analysis.Assemble();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
//            analysis.Solve();
//            
#ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
#endif
            
//            ofstream out4("gmesh4.txt");
//            gmesh->Print(out4);
//            
//            ofstream out44("cmesh4.txt");
//            cmesh->Print(out44);
//
            
            //Arquivo de saida para plotar a solução
            if(p==100)
            {
                string plotfile("../Solution_Hibrido.vtk");
                SaidaSolucao(analysis, plotfile);
            }
            
            //visualizar matriz no vtk
            if(ndiv==100){
                TPZFMatrix<REAL> vismat(100,100);
                cmesh->ComputeFillIn(100,vismat);
                VisualMatrixVTK(vismat,"matrixstruct.vtk");
            }
            
            
            TPZVec<STATE> errorPrimalL2;
            errorPrimalL2.Resize(3, 0.);
            
//            saidaerrosHibrido<<"Erro da simulacao multifisica  para a Pressao\n";
//            TPZCompMeshTools::UnCondensedElements(cmesh);
//            TPZCompMeshTools::UnGroupElements(cmesh);
//            ErrorH1(cmesh, saidaerrosHibrido,errorPrimalL2);
            
//            analysis.SetExact(SolProblema);
//            TPZVec<REAL> erros(3);
//            analysis.PostProcessError(erros);
            
            REAL totalbanda = NDoFCond*NDoFCond;
            REAL NumZeros = totalbanda - nNzeros;
            REAL razao = NumZeros/totalbanda;
            
            myerrorfile << ndiv <<  setw(13) << NDoF << setw(12) << NDoFCond << setw(13)<< NDoFCond*NDoFCond
            << setw(15) << NumZeros << setw(12) << razao;
#ifdef USING_BOOST
            myerrorfile << "    " << (t2-t1) << "     " << (t3-t2) << "     " << (t2-t1)+(t3-t2);
#endif
            myerrorfile << setw(12) << errorPrimalL2[1] << setw(15) << errorPrimalL2[2] <<std::endl;
        }
        
        cmesh->CleanUp();
        delete cmesh;
        delete gmesh;
    }
    
    return EXIT_SUCCESS;
}

//Rodar H1 ou DG
int mainH1()
{
    //InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    ofstream saidaerrosH1("../Erro-H1.txt");
    
    int maxp = 5;
    int ndiv =0;
    
    std::ofstream myerrorfile("Resultados-H1-ou-DG.txt");
    
    TPZGeoMesh * gmesh;
    for(int p = 4; p<maxp; p++)
    {
        myerrorfile<<"\nORDEM p = "<<p <<"\n\n";
        myerrorfile << "ndiv" << setw(10) <<"NDoF"<< setw(12)<<"NDoFCond" << "     Entradas" <<"       NumZeros" <<
        "       Razao" <<setw(19)<< "Assemble"<< setw(20)<<"Solve" << setw(20) <<"Ttotal" <<setw(12) <<"Error u" << setw(16)<<"Error gradU\n"<<std::endl;
        
        //refinamentos h adptativos
        for(ndiv = 0; ndiv<5; ndiv++)
        {
            int64_t NDoF=0, NDoFCond=0;
            int64_t nNzeros=0;
            
            //Geo Mesh
            if(!fTetra){
                gmesh = CreateOneCubo(ndiv);
            }
            else{
                gmesh = CreateOneCuboWithTetraedrons(ndiv);
            }
            
            REAL Beta0 = 6.;
            TPZCompMesh * cmesh= MalhaCompH1(gmesh,p,Beta0);//Malha Computacional H1 ou DG
            
            if(rodarSIPGD){
                cmesh->LoadReferences();
                TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
                cmesh->ExpandSolution();
                cmesh->CleanUpUnconnectedNodes();
            }

            NDoF = cmesh->NEquations();
            
            if(!rodarSIPGD)
            {
                cmesh->Reference()->ResetReference();
                cmesh->LoadReferences();
                if(rodarSIPGD) DebugStop();
                for (int64_t iel=0; iel<cmesh->NElements(); iel++) {
                    TPZCompEl *cel = cmesh->Element(iel);
                    if(!cel) continue;
                    TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel);
                }
                cmesh->ExpandSolution();
                cmesh->CleanUpUnconnectedNodes();
                
            }
            NDoFCond = cmesh->NEquations();
            
                
            // Resolvendo o sistema linear
            TPZAnalysis analysis(cmesh);
            int64_t neq = NDoFCond;
            TPZVec<int64_t> skyline;
            cmesh->Skyline(skyline);
            TPZSkylMatrix<STATE> matsky(neq,skyline);
            nNzeros = matsky.GetNelemts();
            
            //MKL Partdiso
            TPZSymetricSpStructMatrix strmat(cmesh);
            int nthreads = 12;
            strmat.SetNumThreads(nthreads);
            analysis.SetStructuralMatrix(strmat);
            
            //Par Frontal
//            TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh);
//            strmat.SetDecomposeType(ELDLt);
//            strmat.SetNumThreads(6);
//            analysis.SetStructuralMatrix(strmat);
            
            
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt); //caso simetrico
            analysis.SetSolver(step);
            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            analysis.Assemble();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
//            analysis.Solve();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
#endif
        
            //Arquivo de saida para plotar a solução
            if(p==100)
            {
                string plotfile("../Solution_H1.vtk");
                SaidaSolucao(analysis, plotfile);
            }
            
            if(ndiv==100)
            {
                //visualizar matriz no vtk
                TPZFMatrix<REAL> vismat(100,100);
                cmesh->ComputeFillIn(100,vismat);
                VisualMatrixVTK(vismat,"matrixstruct.vtk");
            }
            
            
            TPZVec<STATE> errorPrimalL2;
            errorPrimalL2.Resize(3, 0.);
            
//            saidaerrosH1<<"Erro da simulacao multifisica  para a Pressao\n";
//            ErrorH1(cmesh, saidaerrosH1,errorPrimalL2);
            
            REAL totalbanda = NDoFCond*NDoFCond;
            REAL NumZeros = totalbanda - nNzeros;
            REAL razao = NumZeros/totalbanda;
            
            myerrorfile << ndiv <<  setw(13) << NDoF << setw(12) << NDoFCond << setw(13)<< NDoFCond*NDoFCond
            << setw(15) << NumZeros << setw(12) << razao;
#ifdef USING_BOOST
            myerrorfile << "    " << (t2-t1) << "     " << (t3-t2) << "     " << (t2-t1)+(t3-t2);
#endif
            myerrorfile << setw(12) << errorPrimalL2[1] << setw(15) << errorPrimalL2[2] <<std::endl;
        }
        
    }
    return EXIT_SUCCESS;
}

//Rodar Metodo Misto
int mainMixed()
{
    HDivPiola = 1;
    bool SecondIntegration = false;//Rodar false
    
    //InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    TPZVec<REAL> erros;
    
    bool HDivMaisMais = false;//true: formulacao Pk**
    int nmais = 1;

    ofstream saidaerrosHdiv("../Erro-Misto.txt");
    ofstream saidaerrosH1("../Erro-H1.txt");
    
    int maxp = 5;
    int ndiv =0;
    
    
    std::ofstream myerrorfile("Resultados-HDiv.txt");
    
    TPZGeoMesh * gmesh;
    for(int p = 4; p<maxp; p++)
    {
        int pq = p;
        int pp = p;
        if(HDivMaisMais){
            pp = p + nmais;//Aqui = comeca com 1
        }
        
        myerrorfile<<"\nORDEM p = "<<p <<"\n\n";
        myerrorfile << "ndiv" << setw(10) <<"NDoF"<< setw(12)<<"NDoFCond" << "     Entradas" <<"       NumZeros" <<
        "       Razao" <<setw(19)<< "Assemble"<< setw(20)<<"Solve" << setw(20) <<"Ttotal" <<setw(12) <<"Error u" << setw(16)<<"Error gradU\n"<<std::endl;
        
        //refinamentos h adptativos
        for(ndiv = 0; ndiv<5; ndiv++)
        {
            int64_t NDoF=0, NDoFCond=0;
            int64_t nNzeros=0;
            
            //Geo Mesh
            if(!fTetra){
                gmesh = CreateOneCubo(ndiv);
            }
            else{
                gmesh = CreateOneCuboWithTetraedrons(ndiv);
            }
            
                
            // Criando a primeira malha computacional
            TPZCompMesh * cmesh1= CMeshFlux(pq, gmesh);
            if(HDivMaisMais)
            {
                ChangeInternalConnectOrder(cmesh1, nmais);
            }
            
            // Criando a segunda malha computacional
            TPZCompMesh * cmesh2 = CMeshPressure(pp, gmesh);
            
            // Criando a malha computacional multifísica
            //malha multifisica
            TPZManVector<TPZCompMesh *,2> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            TPZCompMesh * mphysics = MalhaCompMultifisica(meshvec, gmesh, SecondIntegration);
            
            NDoF = meshvec[0]->NEquations() + meshvec[1]->NEquations();
            
            NDoFCond = mphysics->NEquations();
            
            
            // Resolvendo o sistema linear
            TPZAnalysis analysis(mphysics);
            int64_t neq = NDoFCond;
            TPZVec<int64_t> skyline;
            mphysics->Skyline(skyline);
            TPZSkylMatrix<STATE> matsky(neq,skyline);
            nNzeros = matsky.GetNelemts();
            
            //MKL Partdiso
            TPZSymetricSpStructMatrix strmat(mphysics);
            int nthreads = 12;
            strmat.SetNumThreads(nthreads);
            analysis.SetStructuralMatrix(strmat);
            
            //Par Frontal
            // TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(mphysics);
            // strmat.SetDecomposeType(ELDLt);
            // strmat.SetNumThreads(6);
            // analysis.SetStructuralMatrix(strmat);
            
            
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt); //caso simetrico
            analysis.SetSolver(step);
            
    #ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
    #endif
            analysis.Assemble();
            
    #ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
    #endif
//            analysis.Solve();
            
    #ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    #endif
            
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            
            //Arquivo de saida para plotar a solução
            if(p==100)
            {
                string plotfile("../Solution_Hdiv.vtk");
                SaidaSolucao(analysis, plotfile);
            }
            
            if(ndiv==100)
            {
              //visualizar matriz no vtk
              TPZFMatrix<REAL> vismat(100,100);
              mphysics->ComputeFillIn(100,vismat);
              VisualMatrixVTK(vismat,"matrixstruct.vtk");
            }
            
            TPZVec<STATE> errorPrimalL2;
            TPZVec<STATE> errorsHDiv;
            errorPrimalL2.Resize(3, 0.);
            errorsHDiv.Resize(4, 0.);

            
//            saidaerrosHdiv<<"Erro da simulacao multifisica  para o Fluxo\n";
//            ErrorHDiv2(cmesh1,saidaerrosHdiv,errorsHDiv);
//            
//            saidaerrosHdiv<<"Erro da simulacao multifisica  para a Pressao\n";
//            ErrorH1(cmesh2, saidaerrosHdiv,errorPrimalL2);
            
            REAL totalbanda = NDoFCond*NDoFCond;
            REAL NumZeros = totalbanda - nNzeros;
            REAL razao = NumZeros/totalbanda;
            
            myerrorfile << ndiv <<  setw(13) << NDoF << setw(12) << NDoFCond << setw(13)<< NDoFCond*NDoFCond
            << setw(15) << NumZeros << setw(12) << razao;
#ifdef USING_BOOST
            myerrorfile << "    " << (t2-t1) << "     " << (t3-t2) << "     " << (t2-t1)+(t3-t2);
#endif
            myerrorfile << setw(12) << errorPrimalL2[1] << setw(15) << errorsHDiv[1] <<std::endl;
        }
    }
    return EXIT_SUCCESS;
}


void GroupElements(TPZCompMesh *cmesh, int dimproblema)
{
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();

    int nel = cmesh->NElements();
    std::set<TPZCompEl *> celset;
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();
        if (dim ==dimproblema) {
            celset.insert(cel);
        }
    }
    std::set<int> elgroupindices;
    
    for (std::set<TPZCompEl *>::iterator it = celset.begin(); it != celset.end(); it++) {
        
        std::list<TPZCompEl *> group;
        group.push_back(*it);
        TPZCompEl *cel = *it;
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != dimproblema-1) {
                continue;
            }
            TPZStack<TPZCompElSide> connected;
            TPZCompElSide celside(cel,is);
            celside.EqualLevelElementList(connected, false, false);
            int neq = connected.NElements();
            for (int eq=0; eq<neq; eq++) {
                TPZCompElSide eqside = connected[eq];
                TPZCompEl *celeq = eqside.Element();
                TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(celeq);
                if (!intface) {
                    continue;
                }
                TPZCompEl *left = intface->LeftElement();
                if (left == cel) {
                    //put in the group
                    group.push_back(intface);
                }
            }
        }
        
        //#ifdef LOG4CXX
        //        if (logger->isDebugEnabled())
        //        {
        //            std::stringstream sout;
        //            for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
        //                (*it)->Print(sout);
        //            }
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //#endif
        
        int64_t index;
        TPZElementGroup *celgroup = new TPZElementGroup(*cmesh,index);
        elgroupindices.insert(index);
        for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
            celgroup->AddElement(*it);
        }
    }
    cmesh->ComputeNodElCon();
    
    for (std::set<int>::iterator it = elgroupindices.begin(); it!=elgroupindices.end(); it++) {
        TPZCompEl *cel = cmesh->ElementVec()[*it];
        TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel);
    }
    
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        cmesh->Print(sout);
//        LOGPZ_DEBUG(logger, sout.str())
//    }
//#endif
    
}


TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder, bool ismultiplierH1, REAL penaltyConst){
    
    //TPZCompEl::SetgOrder(porder);
    TPZCompMesh *cmesh = new TPZCompMesh(&gmesh);
    int dim = gmesh.Dimension();
    
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    
    // Criar e inserir os materiais na malha
    REAL beta = penaltyConst;
    TPZMatDualHybridPoisson *material = new TPZMatDualHybridPoisson(matId,0.,beta);
    TPZMaterial * automat(material);
    cmesh->InsertMaterialObject(automat);
    material->SetDimension(dim);
    
    //Condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.),val2(1,1,0.);
    
    int int_order = 4;
   
    //sol exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolProblema, 5);
    material->SetForcingFunctionExact(solexata);
    
    //vetor de carga: lada direita da equacao
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForcingTang, 5);
    dum->SetPolynomialOrder(int_order);
    force = dum;
    material->SetForcingFunction(force);

    //bc1
    TPZMaterial *bnd = automat->CreateBC(automat, bc1, 0, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc1exata;
    TPZDummyFunction<STATE> *dum1;
    dum1 = new TPZDummyFunction<STATE>(Dirichlet, 5);//NeumannBC1
    dum1->SetPolynomialOrder(int_order);
    bc1exata = dum1;
    bnd->SetForcingFunction(bc1exata);
    cmesh->InsertMaterialObject(bnd);

    
    //bc2
    bnd = automat->CreateBC (automat, bc2, 1, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc2exata;
    TPZDummyFunction<STATE> *dum2;
    dum2 = new TPZDummyFunction<STATE>(NeumannBC2, 5);//NeumannBC2
    dum2->SetPolynomialOrder(int_order);
    bc2exata = dum2;
    bnd->SetForcingFunction(bc2exata);
    cmesh->InsertMaterialObject(bnd);
    
    //bc3
    bnd = automat->CreateBC (automat, bc3, 0, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc3exata;
    TPZDummyFunction<STATE> *dum3;
    dum3 = new TPZDummyFunction<STATE>(Dirichlet, 5);
    dum3->SetPolynomialOrder(int_order);
    bc3exata = dum3;
    bnd->SetForcingFunction(bc3exata);
    cmesh->InsertMaterialObject(bnd);
    
    //bc4
    bnd = automat->CreateBC (automat, bc4, 0, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc4exata;
    TPZDummyFunction<STATE> *dum4;
    dum4 = new TPZDummyFunction<STATE>(Dirichlet, 5);
    dum4->SetPolynomialOrder(int_order);
    bc4exata = dum4;
    bnd->SetForcingFunction(bc4exata);
    cmesh->InsertMaterialObject(bnd);

    
    //bc0
    bnd = automat->CreateBC (automat, bc0, 0, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc0exata;
    TPZDummyFunction<STATE> *dum0;
    dum0 = new TPZDummyFunction<STATE>(Dirichlet, 5);
    dum0->SetPolynomialOrder(int_order);
    bc0exata = dum0;
    bnd->SetForcingFunction(bc0exata);
    cmesh->InsertMaterialObject(bnd);
    
    //bc5
    bnd = automat->CreateBC (automat, bc5, 1, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc5exata;
    TPZDummyFunction<STATE> *dum5;
    dum5 = new TPZDummyFunction<STATE>(NeumannAcima, 5);//NeumannAcima
    dum5->SetPolynomialOrder(int_order);
    bc5exata = dum5;
    bnd->SetForcingFunction(bc5exata);
    cmesh->InsertMaterialObject(bnd);

    
    // Ajuste da estrutura de dados computacional
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    std::set<int> matids;
    matids.insert(matId);
    cmesh->AutoBuild(matids);
    
    cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    matids.clear();
    int nbc = dim*2;
    for (int i=1; i<=nbc; i++)
    {
        if(dim==2){
            matids.insert(-(i+1));
        }else{
            matids.insert(-i);
        }
    }
    
    cmesh->SetDefaultOrder(porder);
    cmesh->AutoBuild(matids);
    cmesh->SetDimModel(dim);
    
    cmesh->AdjustBoundaryElements();//ajusta as condicoes de contorno
    cmesh->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateInterfaceElements(cmesh,true);
    
    //
    //    matids.insert(1);
    //    if(!fsolsuave){
    //        for (int i=1; i<=nbc; i++) {
    //            matids.insert(-i);
    //        }
    //    }else{
    //        matids.insert(-1);
    //    }
    //
    //    comp->ApproxSpace().Hybridize(*comp, matids, ismultiplierH1);
    cmesh->AdjustBoundaryElements();
    return cmesh;
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
    TPZMaterial * BCond5 = material->CreateBC(mat, bc0, bcdirichlet, val1, val2);
    TPZMaterial * BCond6 = material->CreateBC(mat, bc5, bcdirichlet, val1, val2);
    
    ///Inserir condicoes de contorno
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    cmesh->InsertMaterialObject(BCond6);
    
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
            if(fTriang==true || fTetra==true) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

TPZCompMesh *MalhaCompH1(TPZGeoMesh * gmesh,int ordem, REAL penaltyConst){
    
    int dim =gmesh->Dimension();
    TPZMatLaplacian *material = new TPZMatLaplacian(matId,dim);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    if(rodarSIPGD)
    {
        material->SetSymmetric();
        material->SetSolutionPenalty();
        material->SetValPenaltyConstant(penaltyConst);
    }
    
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolProblema, 5);
    material->SetForcingFunctionExact(solexata);
    
    int int_order = 4;
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForcingF, 5);
    dum->SetPolynomialOrder(int_order);
    force = dum;
    material->SetForcingFunction(force);
    
    //inserindo o material na malha computacional
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(ordem);
    
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(1,1,0.);
    
    //bc1
    TPZMaterial *bnd = mat->CreateBC(mat, bc1, 0, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc1exata;
    TPZDummyFunction<STATE> *dum1;
//    dum1 = new TPZDummyFunction<STATE>(NeumannBC1, 5);//NeumannBC1
    dum1 = new TPZDummyFunction<STATE>(Dirichlet, 5);
    dum1->SetPolynomialOrder(int_order);
    bc1exata = dum1;
    bnd->SetForcingFunction(bc1exata);
    cmesh->InsertMaterialObject(bnd);
    
    //bc2
    bnd = mat->CreateBC (mat, bc2, 1, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc2exata;
    TPZDummyFunction<STATE> *dum2;
    dum2 = new TPZDummyFunction<STATE>(NeumannBC2, 5);//NeumannBC2
    dum2->SetPolynomialOrder(int_order);
    bc2exata = dum2;
    bnd->SetForcingFunction(bc2exata);
    cmesh->InsertMaterialObject(bnd);
    
    //bc3
    bnd = mat->CreateBC (mat, bc3, 0, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc3exata;
    TPZDummyFunction<STATE> *dum3;
    dum3 = new TPZDummyFunction<STATE>(Dirichlet, 5);
    dum3->SetPolynomialOrder(int_order);
    bc3exata = dum3;
    bnd->SetForcingFunction(bc3exata);
    cmesh->InsertMaterialObject(bnd);
    
    //bc4
    bnd = mat->CreateBC (mat, bc4, 0, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc4exata;
    TPZDummyFunction<STATE> *dum4;
    dum4 = new TPZDummyFunction<STATE>(Dirichlet, 5);
    dum4->SetPolynomialOrder(int_order);
    bc4exata = dum4;
    bnd->SetForcingFunction(bc4exata);
    cmesh->InsertMaterialObject(bnd);
    
    //bc0
    bnd = mat->CreateBC (mat, bc0, 0, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc0exata;
    TPZDummyFunction<STATE> *dum0;
    dum0 = new TPZDummyFunction<STATE>(Dirichlet, 5);
    dum0->SetPolynomialOrder(int_order);
    bc0exata = dum0;
    bnd->SetForcingFunction(bc0exata);
    cmesh->InsertMaterialObject(bnd);
    
    //bc5
    bnd = mat->CreateBC (mat, bc5, 1, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bc5exata;
    TPZDummyFunction<STATE> *dum5;
    dum5 = new TPZDummyFunction<STATE>(NeumannAcima, 5);//NeumannAcima
    dum5->SetPolynomialOrder(int_order);
    bc5exata = dum5;
    bnd->SetForcingFunction(bc5exata);
    cmesh->InsertMaterialObject(bnd);
    
    if(rodarSIPGD){
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }else{
        cmesh->SetAllCreateFunctionsContinuous();
    }

    
    //Fazendo auto build
    cmesh->AutoBuild();
    
    
    if(rodarSIPGD)
    {
        int64_t nel = cmesh->NElements();
        for(int64_t i=0; i<nel; i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            if(!celdisc) continue;
            //            celdisc->SetConstC(1.);
            //            celdisc->SetCenterPoint(0, 0.);
            //            celdisc->SetCenterPoint(1, 0.);
            //            celdisc->SetCenterPoint(2, 0.);
            celdisc->SetFalseUseQsiEta();
            //if (dim_el != dim) continue;
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                if(fTetra==true || fTriang==true) celdisc->SetTotalOrderShape();
                else celdisc->SetTensorialShape();
            }
        }
    }

    
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
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolProblema, 5);
    material->SetForcingFunctionExact(solexata);
    
    int int_order = 4;
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
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
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
    else if (rodarH1 || rodarSIPGD)
    {
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
    else
    {
        TPZManVector<std::string,10> scalnames(3), vecnames(2);
        scalnames[0] = "Solution";
        scalnames[1]="POrder";
        scalnames[2] = "ExactSolution";
        vecnames[0]= "Grad";
        vecnames[1] = "ExactGrad";
        
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
    TPZManVector<REAL,10> globerrors(10,0.);
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

void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out, TPZVec<STATE> &errorL2)
{
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globerrors(10,0.);
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
    
    errorL2.Resize(3,0.);
    errorL2[0] = sqrt(globerrors[0]);//H1
    errorL2[1] = sqrt(globerrors[1]);//L2
    errorL2[2] = sqrt(globerrors[2]);//SemiH1
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


void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &res){
    
    
    double x = pt[0];
    double y = pt[1];
    res[0] = 0.;
    
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

void ForcingTang(const TPZVec<REAL> &pt, TPZVec<STATE> &res,TPZFMatrix<STATE> &disp)
{
    ForcingF(pt, res);
    //res[0] *=-1.;
}

void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result)
{
    TPZFMatrix<STATE> du(3,1);
    SolProblema(loc,result,du);
}

void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result)
{
    REAL normal[3] = {1.,0.,0.};
    TPZManVector<REAL> u(1);
    TPZFNMatrix<5,STATE> du(3,1);
    
    SolProblema(loc,result,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
}

void NeumannBC1(const TPZVec<REAL> &loc, TPZVec<STATE> &result){

    REAL normal[3] = {0.,-1.,0.};
    TPZManVector<STATE> u(1);
    TPZFNMatrix<5,STATE> du(3,1);

    SolProblema(loc,result,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
}

void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result)
{
    REAL normal[3] = {0.,0.,1.};
    TPZManVector<REAL> u(1);
    TPZFNMatrix<5,STATE> du(3,1);
    
    SolProblema(loc,result,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
    
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



////#define SolutionPoly
//#define SolutionShock
//
////with hybrid method
//
////nref: numero de refinamento
//TPZGeoMesh *CreateOneCubo(int ndiv);
//TPZGeoMesh * CreateOneCuboWithTetraedrons(int ndiv);
//TPZGeoMesh *GMesh2D(bool ftriang);
//TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder,bool ismultiplierH1, int ndiv);
//
//TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim, bool rodarSIPGD, int ndiv);
//
//void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem);
//bool MyDoubleComparer(REAL a, REAL b);
//
//void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);
//void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL,3> > &Points);
//void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined);
//void RegularizeMesh(TPZGeoMesh *gmesh);
//
//
//void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide);
//void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
//void Prefinamento(TPZCompMesh * cmesh, int ndiv,int porder);
//void PrefinamentoRibsHybridMesh(TPZCompMesh *cmesh);
//void SetPOrderRibsHybridMesh(TPZCompMesh *cmesh, int porder);
//
//void GroupElements(TPZCompMesh *cmesh, int dimproblema);
//
//// Exact functions
//void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannBC1(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannAbaixo(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void ForcingShockProblem2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
//
//void ExactSolution(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
//void ForcingFunction(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &df);
//void ForcingFunctionII(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
//
////shock problem
//void SolShockProblem(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
//void ForcingShockProblem(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &df);
//
//
//// Polynomial
//void PolyProblem(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
//void ForcingPolyProblem(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &df);
//
////problema Suave
//void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &df);
//void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
//
////malha Hdiv
//TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh);
//TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh);
//TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh, bool hdivskeleton, int ndiv);
//void ErrorH1(TPZCompMesh *l2mesh, TPZVec<STATE> &Error /*,std::ostream &out*/);
//void ErrorHDiv(TPZCompMesh *hdivmesh, TPZVec<STATE> &Error /*,std::ostream &out*/);
//
//void AjustarContorno(TPZGeoMesh *gmesh);
//
//#ifdef USING_BOOST
//#include "boost/date_time/posix_time/posix_time.hpp"
//#endif
//
//int dim_problema = 3;
//int nbc = dim_problema*2;
//bool fTriang = false;
//int flevel=2;
//
//bool rodarH1 = true;
//bool rodarSIPGD = false;
//bool rodarHdiv = false;
//
//REAL alpha_param = 200.;
//
//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.hibridoprimal"));
//#endif
//void ChangeInternalConnectOrder(TPZCompMesh *mesh);
//
//int main(int argc, char *argv[])
//{
//#ifdef LOG4CXX
//    InitializePZLOG();
//#endif
//
//    ///Refinamento
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//
//
//    bool multiplicadorH1 = true;
//    bool MudarOrdemPdoMultiplicador = false;
//    std::ofstream myerrorfile("Simulacao-Primal.txt");
//
//    if(rodarH1){
//        myerrorfile<<"\nDADOS PARA O REFINAMENTO hp E FORMULAÇÃO H1"<<std::endl;
//        myerrorfile<<std::endl;
//
//    }else if(rodarSIPGD){
//        myerrorfile<<"\nDADOS PARA O REFINAMENTO hp E FORMULAÇÃO DG"<<std::endl;
//        myerrorfile<<std::endl;
//
//    }
//    else if(multiplicadorH1){
//        myerrorfile<<"\nDADOS PARA O REFINAMENTO hp E MULTIPLICADOR CONTINUO"<<std::endl;
//        myerrorfile<<std::endl;
//
//    }else {
//        myerrorfile<<"\nDADOS PARA O REFINAMENTO hp E MULTIPLICADOR DESCONTINUO"<<std::endl;
//        myerrorfile<<std::endl;
//    }
//
//    TPZCompMesh *cmesh;
//    TPZGeoMesh *gmesh;
//
//    int pini =1;
//    for(int p = pini; p<2; p++)
//    {
//
//        myerrorfile<<"\nORDEM p = "<<p <<"\n\n";
//        if(dim_problema==0){
//            myerrorfile << "ndiv" << setw(10) <<"NDoF"<< setw(15)<<"NDoFCond" << setw(19)<< "Assemble"<< setw(20)<<
//            "Solve" << setw(20) <<"Ttotal" << setw(18) <<"Error u" << setw(20)<<"Error gradU\n";
//        }else{
//            myerrorfile << "ndiv" << setw(10) <<"NDoF"<< setw(12)<<"NDoFCond" << "     Entradas" <<"       NumZeros" <<
//            "       Razao" <<setw(19)<< "Assemble"<< setw(20)<<"Solve" << setw(20) <<"Ttotal" << setw(12) <<"Error u" << setw(16)<<"Error gradU\n";
//        }
//        for(int ndiv=1; ndiv<2; ndiv++){
//
//
//            if(dim_problema==2){
//                gmesh = GMesh2D(fTriang);//malha geometrica
//                UniformRefine(gmesh, flevel+ndiv);
////                RefiningNearCircunference(dim_problema, gmesh,ndiv,1);
////                //DirectionalRef(gmesh, 1, ndiv);
////                AjustarContorno(gmesh);
//
//
//            }else{
//                gmesh = CreateOneCubo(ndiv);
//                //gmesh = CreateOneCuboWithTetraedrons(ndiv);
//            }
//
////            {
////                std::ofstream out("gmesh.txt");
////                gmesh->Print(out);
////                std::ofstream filemesh("MalhaGeometricaInicial.vtk");
////                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
////            }
//
//
//            int64_t NDoF=0, NDoFCond=0;
//            int64_t nNzeros=0;
//
//            if(!rodarH1 && !rodarSIPGD){
//                cmesh= CreateHybridCompMesh(*gmesh, p, multiplicadorH1, ndiv);//malha computacional
//
//                //------- Criar elementos de Lagrange (Ribs)--------
//                //materiais do problema
//                std::set<int>matids;
//                matids.insert(matId);
//                for(int i=1; i<=nbc; i++)
//                {
//                    if(dim_problema==2){
//                        matids.insert(-(i+1));
//                    }else{
//                        matids.insert(-i);
//                    }
//                }
//                //if(dim_problema==2) Prefinamento(cmesh, ndiv, p);
//                cmesh->ApproxSpace().Hybridize(*cmesh, matids, multiplicadorH1);
//
//
//                if(MudarOrdemPdoMultiplicador){
//                    SetPOrderRibsHybridMesh(cmesh, 2/*p-(pini-1)*/);
//                    cmesh->CleanUpUnconnectedNodes();
//                    cmesh->ExpandSolution();
//                }
//
////                {
////                    std::ofstream out("cmeshHib2.txt");
////                    cmesh->Print(out);
////                }
//
//
//                //myerrorfile << "\nRefinamento h = "<< ndiv <<"\n";
//                NDoF = cmesh->NEquations();
//                //myerrorfile << "\nDOF Total = "<< cmesh->NEquations() << "\n";
//
//                //condesacao estatica
//                GroupElements(cmesh, dim_problema);
//                cmesh->LoadReferences();//mapeia para a malha geometrica lo
//                NDoFCond = cmesh->NEquations();
//
////                {
////                    std::ofstream out("cmeshHib2-AfterGroup.txt");
////                    cmesh->Print(out);
////                }
//            }
//            else {//Malha H1
//                {
//                    std::ofstream gout("gmeshH1-1.txt");
//                    gmesh->Print(gout);
//                }
//                cmesh = CMeshH1(gmesh, p, dim_problema, rodarSIPGD, ndiv);
//                cmesh->ExpandSolution();
//                cmesh->CleanUpUnconnectedNodes();
//                
//                {
//                    std::ofstream out("cmeshH1-1.txt");
//                    cmesh->Print(out);
//                    std::ofstream gout("gmeshH1-1.txt");
//                    cmesh->Reference()->Print(gout);
//                }
//
//                if(dim_problema==2) Prefinamento(cmesh, ndiv, p);
//
//                if(rodarSIPGD){
//                    cmesh->LoadReferences();
//                    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
//                    cmesh->ExpandSolution();
//                    cmesh->CleanUpUnconnectedNodes();
//                }
////
////                {
////                    std::ofstream out("cmeshH1-2.txt");
////                    cmesh->Print(out);
////
////                    std::ofstream out2("gmesh-2.txt");
////                    gmesh->Print(out2);
////                }
//
//                NDoF = cmesh->NEquations();
//                //condensar
//                if(rodarH1){
//                    cmesh->Reference()->ResetReference();
//                    cmesh->LoadReferences();
//                    if(rodarSIPGD) DebugStop();
//                    for (int64_t iel=0; iel<cmesh->NElements(); iel++) {
//                        TPZCompEl *cel = cmesh->Element(iel);
//                        if(!cel) continue;
//                        TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel);
//                }
//                    cmesh->ExpandSolution();
//                    cmesh->CleanUpUnconnectedNodes();
//                    
//                }
//                
//
//                NDoFCond = cmesh->NEquations();
//            }
//
//
//            if(1)
//            {
//                std::ofstream gout("../gmesh.txt");
//                cmesh->Reference()->Print(gout);
//                std::ofstream out("../cmesh.txt");
//                cmesh->Print(out);
//            }
//            //Resolver problema
//            TPZAnalysis analysis(cmesh);
//            if(dim_problema==2){
//
////                TPZSkylineStructMatrix skylstr(cmesh); //caso simetrico
////                //TPZSkylineNSymStructMatrix skylstr(cmesh); //caso nao simetrico
////                //skylstr.SetNumThreads(8);
////                analysis.SetStructuralMatrix(skylstr);
//
//                int64_t neq = NDoFCond;
//                TPZVec<int64_t> skyline;
//                cmesh->Skyline(skyline);
//                TPZSkylMatrix<STATE> matsky(neq,skyline);
//                nNzeros = matsky.GetNelemts();
//                
//
//                TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh);
//                strmat.SetDecomposeType(ELDLt);
//                strmat.SetNumThreads(6);
//                analysis.SetStructuralMatrix(strmat);
//
//
////                TPZBandStructMatrix bdmat(cmesh);
////                //bdmat.SetNumThreads(8);
////                analysis.SetStructuralMatrix(bdmat);
//
//
////                TPZParFrontStructMatrix<TPZFrontNonSym<STATE> > strmat(cmesh);
////                strmat.SetDecomposeType(ELU);
////                strmat.SetNumThreads(8);
////                analysis.SetStructuralMatrix(strmat);
//
//            }else{
//
//                int64_t neq = NDoFCond;
//                TPZVec<int64_t> skyline;
//                cmesh->Skyline(skyline);
//                TPZSkylMatrix<STATE> matsky(neq,skyline);
//                nNzeros = matsky.GetNelemts();
//
//                TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh);
//                strmat.SetDecomposeType(ELDLt);
//                strmat.SetNumThreads(8);
//                analysis.SetStructuralMatrix(strmat);
//            }
//
//            TPZStepSolver<STATE> step;
//            step.SetDirect(ELDLt); //caso simetrico
//            //step.SetDirect(ELU);
//            analysis.SetSolver(step);
//
//#ifdef USING_BOOST
//            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
//#else
//            REAL t1=0.;
//#endif
//            analysis.Assemble();
//
//#ifdef USING_BOOST
//            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
//#else
//            REAL t2 = 0.;
//#endif
//            analysis.Solve();
//
//#ifdef USING_BOOST
//            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
//#else
//            REAL t3 = 0.;
//#endif
//
//            //REAL t1=0., t2=0., t3=0.;
//
//            //            std::ofstream out("cmeshHib22.txt");
//            //            cmesh->Print(out);
//
//            if(p!=4){
//                TPZVec<std::string> scalnames(3), vecnames;
//                scalnames[0] = "Solution";
//                scalnames[1] = "POrder";
//                scalnames[2] = "ExactSolution";
//                if(rodarH1 || rodarSIPGD){
//                    vecnames.Resize(1);
//                    vecnames[0] = "Derivative";
//                }else{
//                    vecnames.Resize(2);
//                    vecnames[0] = "Grad";
//                    vecnames[1] = "ExactGrad";
//                }
//
//                std::stringstream name;
//                name << "Solution_bima" <<ndiv<< ".vtk";
//                std::string paraviewfile(name.str());
//                analysis.DefineGraphMesh(dim_problema,scalnames,vecnames,paraviewfile);
//                analysis.PostProcess(0);
//
//                //visualizar matriz no vtk
////                TPZFMatrix<REAL> vismat(100,100);
////                cmesh->ComputeFillIn(100,vismat);
////                VisualMatrixVTK(vismat,"matrixstruct.vtk");
//            }
//
//            analysis.SetExact(SolShockProblem);
//            TPZVec<REAL> erros(3);
//            analysis.PostProcessError(erros);
//
//            if(dim_problema==0){
//                //            myerrorfile << ndiv <<  setw(13) << NDoF << setw(15)<< NDoFCond <<"     "<< (t2-t1) << "     "<< (t3-t2) <<"     "<<(t2-t1)+(t3-t2) << setw(18) << erros[1]<< setw(19)<< erros[2]<<std::endl;
//            }else{
//
//                REAL totalbanda = NDoFCond*NDoFCond;
//                REAL NZeros = totalbanda - nNzeros;
//                REAL razao = NZeros/totalbanda;
//                myerrorfile << ndiv <<  setw(13) << NDoF << setw(12)<< NDoFCond <<setw(13)<< NDoFCond*NDoFCond <<setw(15)<<NZeros <<setw(12)<< razao <<"    "<< (t2-t1) << "     "<< (t3-t2) <<"     "<<(t2-t1)+(t3-t2) <<setw(12) << erros[1]<< setw(15)<< erros[2]<<std::endl;
//            }
//            delete cmesh;
//            delete gmesh;
//        }
//        myerrorfile<<"\n-------------------------------------------------------------------------"<<std::endl;
//    }
//
//    return 0;
//}
//
////Malha Hdiv
//bool HDivMaisMais = false;
//bool hp_method = false;
//
//int main2(int argc, char *argv[])
//{
//    //#ifdef LOG4CXX
//    //    InitializePZLOG();
//    //#endif
//    
//    ///Refinamento
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    
//    std::ofstream myerrorfile("Simulacao-Hdiv.txt");
//    myerrorfile<<"\nDADOS PARA O REFINAMENTO hp: Simulacao Hdiv"<<std::endl;
//    
//    
//    HDivPiola = 1;
//    bool hdivskeleton = false;
//    if(HDivPiola != 0)
//    {
//        hdivskeleton = false;
//    }
//    
//    
//    TPZGeoMesh *gmesh;
//    TPZCompMesh * cmesh1;
//    TPZCompMesh * cmesh2;
//    TPZCompMesh * mphysics;
//    int pini = 1;
//    for(int p = pini; p<2; p++)
//    {
//        int pp = p;
//        if(HDivMaisMais){
//            pp = p+1;
//        }
//        
//        myerrorfile<<"\nORDEM p = "<<p <<"\n\n";
//        
//        myerrorfile << "ndiv" << setw(10) <<"NDoF"<< setw(12)<<"NDoFCond" << "     Entradas" <<"       NumZeros" <<
//        "       Razao" <<setw(19)<< "Assemble"<< setw(20)<<"Solve" << setw(20) <<"Ttotal" <<setw(12) <<"Error u" << setw(16)<<"Error gradU\n"<<std::endl;
//        
//        for(int ndiv=2; ndiv<3; ndiv++){
//            
//            if(hp_method){
//                if(dim_problema==2){
//                    gmesh = GMesh2D(fTriang);//malha geometrica
//                    UniformRefine(gmesh, flevel);
//                    RefiningNearCircunference(dim_problema, gmesh,ndiv,1);
//                    AjustarContorno(gmesh);
//                }
//                else{
//                    flevel = 3;
//                    gmesh = CreateOneCubo(flevel);
//                    //gmesh = CreateOneCuboWithTetraedrons(ndiv);
//                    RefiningNearCircunference(dim_problema, gmesh,ndiv,1);
//                    AjustarContorno(gmesh);
//                }
//            }else{
//                gmesh = CreateOneCubo(ndiv);
//            }
//            
//            
//            //            {
//            //                std::ofstream out2("gmesh.txt");
//            //                gmesh->Print(out2);
//            //                std::ofstream filemesh("MalhaGeometricaInicial.vtk");
//            //                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
//            //            }
//            
//            
//            int64_t NDoF=0, NDoFCond=0;
//            int64_t nNzeros=0;
//            
//            cmesh1= CMeshFlux(p, gmesh);
//            if(hp_method) {
//                Prefinamento(cmesh1, ndiv, p);
//            }
//            if(HDivMaisMais){
//                ChangeInternalConnectOrder(cmesh1);
//            }
//            
//            cmesh2 = CMeshPressure(pp, gmesh);
//            if(hp_method) {
//                Prefinamento(cmesh2, ndiv, pp);
//            }
//            
//            NDoF = cmesh1->NEquations() + cmesh2->NEquations();
//            
//            //malha multifisica
//            TPZManVector<TPZCompMesh *,2> meshvec(2);
//            meshvec[0] = cmesh1;
//            meshvec[1] = cmesh2;
//            mphysics = MalhaCompMultifisica(meshvec, gmesh,hdivskeleton, ndiv);
//            
//            NDoFCond = mphysics->NEquations();
//            
//            //Resolver problema
//            TPZAnalysis analysis(mphysics);
//            
//            if(dim_problema==2){
//                
//                //                TPZSkylineStructMatrix skylstr(mphysics); //caso simetrico
//                //                skylstr.SetNumThreads(6);
//                //                analysis.SetStructuralMatrix(skylstr);
//                
//                int64_t neq = NDoFCond;
//                TPZVec<int64_t> skyline;
//                mphysics->Skyline(skyline);
//                TPZSkylMatrix<STATE> matsky(neq,skyline);
//                nNzeros = matsky.GetNelemts();
//                
//                TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(mphysics);
//                strmat.SetDecomposeType(ELDLt);
//                strmat.SetNumThreads(6);
//                analysis.SetStructuralMatrix(strmat);
//            }
//            else{
//                int64_t neq = NDoFCond;
//                TPZVec<int64_t> skyline;
//                mphysics->Skyline(skyline);
//                TPZSkylMatrix<STATE> matsky(neq,skyline);
//                nNzeros = matsky.GetNelemts();
//                
//                //TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(mphysics);
//                //strmat.SetDecomposeType(ELDLt);
//                TPZSymetricSpStructMatrix strmat(mphysics);
//                strmat.SetNumThreads(16);
//                analysis.SetStructuralMatrix(strmat);
//                
//                //                TPZSkylineStructMatrix skylstr(mphysics); //caso simetrico
//                //                analysis.SetStructuralMatrix(skylstr);
//            }
//            
//            TPZStepSolver<STATE> step;
//            step.SetDirect(ELDLt); //caso simetrico
//            //step.SetDirect(ELU);
//            analysis.SetSolver(step);
//            
//            
//#ifdef USING_BOOST
//            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
//#endif
//            analysis.Assemble();
//            //            std::stringstream sout;
//            //            analysis.StructMatrix()->Print("Matriz de RigidezBLABLA: ",sout,EMathematicaInput);F
//            
//#ifdef USING_BOOST
//            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
//#endif
//            analysis.Solve();
//            
//#ifdef USING_BOOST
//            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
//#endif
//            //REAL t1=0., t2=0., t3=0.;
//            
//            //            std::ofstream out("cmeshHib22.txt");
//            //            cmesh->Print(out);
//            
//            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
//            
//            if(ndiv!=5  && p!=5){
//                TPZManVector<std::string,10> scalnames(3), vecnames(2);
//                scalnames[0] = "Pressure";
//                scalnames[1] = "ExactPressure";
//                scalnames[2]="POrder";
//                vecnames[0]= "Flux";
//                vecnames[1]= "ExactFlux";
//                
//                
//                std::stringstream name;
//                name << "Solution_hdiv" <<ndiv<< ".vtk";
//                std::string paraviewfile(name.str());
//                analysis.DefineGraphMesh(dim_problema,scalnames,vecnames,paraviewfile);
//                analysis.PostProcess(0);
//                
//                //visualizar matriz no vtk
//                //                TPZFMatrix<REAL> vismat(100,100);
//                //                mphysics->ComputeFillIn(100,vismat);
//                //                VisualMatrixVTK(vismat,"matrixstruct.vtk");
//            }
//            
//            //            myerrorfile << ndiv <<  setw(13) << NDoF << setw(15)<< NDoFCond <<"    "<< (t2-t1) << "     "<< (t3-t2) <<"     "<<(t2-t1)+(t3-t2) << setw(18);
//            
//            TPZVec<STATE> ErroP;
//            TPZVec<STATE> ErroF;
//            ErrorH1(cmesh2, ErroP /*,myerrorfile*/);
//            ErrorHDiv(cmesh1, ErroF /*,myerrorfile*/);
//            
//            
//            REAL totalbanda = NDoFCond*NDoFCond;
//            REAL NumZeros = totalbanda - nNzeros;
//            REAL razao = NumZeros/totalbanda;
//            
//            
//            myerrorfile << ndiv <<  setw(13) << NDoF << setw(12) << NDoFCond << setw(13)<< NDoFCond*NDoFCond
//            << setw(15) << NumZeros << setw(12) << razao << "    " << (t2-t1) << "     " << (t3-t2) << "     "
//            << (t2-t1)+(t3-t2) << setw(12) << ErroP[1] << setw(15) << ErroF[1] <<std::endl;
//        }
//        
//        myerrorfile <<"\n-------------------------------------------------------------------------"<<std::endl;
//    }
//    
//    return 0;
//}
//
//
//////MUDANCAS DO PHIL
////int main(int argc, char *argv[])
////{
////#ifdef LOG4CXX
////    InitializePZLOG();
////#endif
////    
////    ///Refinamento
////    gRefDBase.InitializeUniformRefPattern(EOned);
////    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
////    
////    std::ofstream myerrorfile("Simulacao-Hdiv.txt",std::ios::app);
////    std::ofstream myerrorfile_fluxo("Simulacao-Hdiv-Fluxo.txt",std::ios::app);
////    myerrorfile<<"\nDADOS PARA O REFINAMENTO hp: Simulacao Hdiv"<<std::endl;
////    
////    
////    HDivPiola = 0;
////    bool hdivskeleton = false;
////    if(HDivPiola != 0)
////    {
////        hdivskeleton = false;
////    }
////    
////    
////    TPZGeoMesh *gmesh;
////    TPZCompMesh * cmesh1;
////    TPZCompMesh * cmesh2;
////    TPZCompMesh * mphysics;
////    int pini = 1;
////    for(int p = pini; p<pini+1; p++)
////    {
////        int pp = p;
////        if(HDivMaisMais){
////            pp = p+1;
////        }
////        
////        myerrorfile<<"\nORDEM p = "<<p <<"\n\n";
////        
////        myerrorfile << "ndiv" << setw(10) <<"NDoF"<< setw(12)<<"NDoFCond" << "     Entradas" <<"       NumZeros" <<
////        "       Razao" <<setw(19)<< "Assemble"<< setw(20)<<"Solve" << setw(20) <<"Ttotal" <<setw(12) <<"Error u" << setw(16)<<"Error gradU\n"<<std::endl;
////        
////        for(int ndiv=1; ndiv<5; ndiv++){
////            
//////            TPZGeoMesh *gmesh;
//////            TPZCompMesh * cmesh1;
//////            TPZCompMesh * cmesh2;
//////            TPZCompMesh * mphysics;
////            
////            if(hp_method){
////                if(dim_problema==2){
////                    gmesh = GMesh2D(fTriang);//malha geometrica
////                    UniformRefine(gmesh, flevel);
////                    RefiningNearCircunference(dim_problema, gmesh,ndiv,1);
////                    AjustarContorno(gmesh);
////                }
////                else{
////                    flevel = 2;
////                    gmesh = CreateOneCubo(flevel);
////                    //gmesh = CreateOneCuboWithTetraedrons(ndiv);
////                    RefiningNearCircunference(dim_problema, gmesh,ndiv,1);
////                    AjustarContorno(gmesh);
////                }
////            }else{
////                gmesh = CreateOneCubo(ndiv);
////            }
////            
////            
//////            {
//////                std::ofstream out2("gmesh.txt");
//////                gmesh->Print(out2);
//////                std::ofstream filemesh("MalhaGeometricaInicial.vtk");
//////                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
//////            }
////            
////            
////            int64_t NDoF=0, NDoFCond=0;
////            int64_t nNzeros=0;
////            
////            cmesh1= CMeshFlux(p, gmesh);
////            if(hp_method) {
////                Prefinamento(cmesh1, ndiv, p);
////            }
////            if(HDivMaisMais){
////                ChangeInternalConnectOrder(cmesh1);
////            }
////            
////            cmesh2 = CMeshPressure(pp, gmesh);
////            if(hp_method) {
////                Prefinamento(cmesh2, ndiv, pp);
////            }
////            {
////                std::ofstream out2("CMeshPressure.txt");
////                cmesh2->Print(out2);
////            }
////
////            
////            NDoF = cmesh1->NEquations() + cmesh2->NEquations();
////            
////            //malha multifisica
////            TPZManVector<TPZCompMesh *,2> meshvec(2);
////            meshvec[0] = cmesh1;
////            meshvec[1] = cmesh2;
////            mphysics = MalhaCompMultifisica(meshvec, gmesh,hdivskeleton, ndiv);
////            mphysics->Reference()->ResetReference();
////            mphysics->LoadReferences();
////            std::set<int64_t> cindex;
////            
//////            TPZCompEl *cel = mphysics->Element(0);
//////            TPZGeoEl *gel = cel->Reference();
//////            gel = gmesh->Element(274);
//////            TPZGeoElSide gelside(gel,gel->NSides()-1);
//////            TPZStack<TPZCompElSide> celstack;
//////            gelside.ConnectedCompElementList(celstack, false, false);
//////            for (int elst=0; elst<celstack.size(); elst++)
//////            {
//////                TPZCompEl *cel = celstack[elst].Element();
//////                int nc = cel->NConnects();
//////                for (int ic=0; ic<nc; ic++) {
//////                    cindex.insert(cel->ConnectIndex(ic));
//////                }
//////            }
//////            std::ostream_iterator< double > output( cout, " " );
//////            std::copy( cindex.begin(), cindex.end(), output );
////////            std::cout << "Connect indexes " << cindex << std::endl;
//////            std::cout << std::endl;
//////
////////            TPZCompMeshTools::GroupElements(mphysics);
////////            TPZCompMeshTools::CreatedCondensedElements(mphysics, true);
////////            mphysics->CleanUpUnconnectedNodes();
////////            mphysics->ExpandSolution();
////            
////            NDoFCond = mphysics->NEquations();
////            
////            //Resolver problema
////            TPZAnalysis analysis(mphysics);
////            
////            if(dim_problema==2){
////                
////                //                TPZSkylineStructMatrix skylstr(mphysics); //caso simetrico
////                //                skylstr.SetNumThreads(6);
////                //                analysis.SetStructuralMatrix(skylstr);
////                
////                int64_t neq = NDoFCond;
////                TPZVec<int64_t> skyline;
////                mphysics->Skyline(skyline);
////                TPZSkylMatrix<STATE> matsky(neq,skyline);
////                nNzeros = matsky.GetNelemts();
////                
////                TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(mphysics);
////                strmat.SetDecomposeType(ELDLt);
////                strmat.SetNumThreads(6);
////                analysis.SetStructuralMatrix(strmat);
////            }
////            else{
////                //Resolver problema
////                TPZAnalysis loc_analysis(mphysics,false);
////                int64_t neq = NDoFCond;
////                TPZVec<int64_t> skyline;
////                mphysics->Skyline(skyline);
////                TPZSkylMatrix<STATE> matsky(neq,skyline);
////                nNzeros = matsky.GetNelemts();
////                
//////                TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(mphysics);
////                
//////                strmat.SetDecomposeType(ELDLt);
//////                strmat.SetNumThreads(16);
//////                analysis.SetStructuralMatrix(strmat);
////                
////                TPZSkylineStructMatrix skylstr(mphysics); //caso simetrico
////                skylstr.SetNumThreads(8);
////                loc_analysis.SetStructuralMatrix(skylstr);
////                analysis.SetStructuralMatrix(skylstr);
////                
////                TPZSymetricSpStructMatrix strmat_p(mphysics);
////                strmat_p.SetNumThreads(16);
////                analysis.SetStructuralMatrix(strmat_p);
////
////                
////                TPZStepSolver<STATE> step_p;
////                step_p.SetDirect(ELDLt);
//////                TPZAutoPointer<TPZMatrix<STATE> > mat = loc_analysis.StructMatrix()->Create();
//////                TPZAutoPointer<TPZMatrix<STATE> > mat2 = mat->Clone();
//////                
//////                TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
//////                TPZStepSolver<STATE> *gmrs = new TPZStepSolver<STATE>(mat2);
//////                step->SetReferenceMatrix(mat2);
//////                step->SetDirect(ELDLt);
//////                gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
//////                TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
//////                TPZAutoPointer<TPZMatrixSolver<STATE> > autogmres = gmrs;
//////                
//////                loc_analysis.SetSolver(autogmres);
//////                loc_analysis.Assemble();
//////                loc_analysis.Solve();
////                TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
////
//////                for (std::set<int64_t>::iterator ic=cindex.begin(); ic!=cindex.end(); ic++) {
//////                    mphysics->ConnectVec()[*ic].Print(*mphysics,cout);
//////                }
////
//////                TPZCompMeshTools::UnGroupElements(mphysics);
//////                TPZCompMeshTools::UnCondensedElements(mphysics);
////                
////            }
//////            TPZCompMeshTools::GroupElements(mphysics);
//////            TPZCompMeshTools::CreatedCondensedElements(mphysics, true);
//////            mphysics->CleanUpUnconnectedNodes();
//////            mphysics->ExpandSolution();
////            
////            
////#ifdef USING_BOOST
////            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
////#endif
//////            TPZSkylineStructMatrix skylstr(mphysics); //caso simetrico
//////            skylstr.SetNumThreads(8);
//////            analysis.SetStructuralMatrix(skylstr);
//////
//////            analysis.Assemble();
//////            TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
//////            
////////            for (std::set<int64_t>::iterator ic=cindex.begin(); ic!=cindex.end(); ic++) {
////////                mphysics->ConnectVec()[*ic].Print(*mphysics,cout);
////////            }
//////
//////            analysis.LoadSolution(mphysics->Solution());
//////            //            std::stringstream sout;
//////            //            analysis.StructMatrix()->Print("Matriz de RigidezBLABLA: ",sout,EMathematicaInput);F
//////            
//////            TPZMatrixSolver<STATE> &solve = analysis.Solver();
//////            TPZMatrix<STATE> *mat = solve.Matrix().operator->();
//////            TPZFMatrix<STATE> locsol(analysis.Solution()),rhs;
//////            locsol.Resize(mat->Rows(), 1);
//////            mat->Multiply(locsol, rhs);
//////            rhs-= analysis.Rhs();
//////            
//////            std::cout << "residual norm " << Norm(rhs) << std::endl;
//////            
//////            {
//////                ofstream out("elementres.txt");
//////                analysis.PrintVectorByElement(out, rhs);
//////            }
////            
////#ifdef USING_BOOST
////            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
////#endif
////            analysis.Solve();
////            
////#ifdef USING_BOOST
////            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
////#endif
////            //REAL t1=0., t2=0., t3=0.;
////            {
////                std::ofstream out("mphysics.txt");
////                mphysics->Print(out);
////            }
////            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
////            
////            if(ndiv!=5  && p!=5){
////                TPZManVector<std::string,10> scalnames(3), vecnames(2);
////                scalnames[0] = "Pressure";
////                scalnames[1] = "ExactPressure";
////                scalnames[2]="POrder";
////                vecnames[0]= "Flux";
////                vecnames[1]= "ExactFlux";
////                
////                
////                std::stringstream name;
////                name << "Solution_hdiv" <<ndiv<< ".vtk";
////                std::string paraviewfile(name.str());
////                analysis.DefineGraphMesh(dim_problema,scalnames,vecnames,paraviewfile);
////                analysis.PostProcess(0);
////                
////                //visualizar matriz no vtk
//////                TPZFMatrix<REAL> vismat(100,100);
//////                mphysics->ComputeFillIn(100,vismat);
//////                VisualMatrixVTK(vismat,"matrixstruct.vtk");
////            }
////            
////            //            myerrorfile << ndiv <<  setw(13) << NDoF << setw(15)<< NDoFCond <<"    "<< (t2-t1) << "     "<< (t3-t2) <<"     "<<(t2-t1)+(t3-t2) << setw(18);
////            
////
////            {
////                std::ofstream sout("flux_cmesh.txt");
////                for (std::set<int64_t>::iterator ic=cindex.begin(); ic!=cindex.end(); ic++) {
////                    sout << "index " << *ic << " ";
////                    mphysics->ConnectVec()[*ic].Print(*mphysics,sout);
////                }
////
////                cmesh1->Print(sout);
////            }
////
////            std::cout << "Computing the error\n";
////            ofstream out("cmeshFlux.txt");
////            cmesh1->Print(out);
////
////            
////            TPZVec<STATE> ErroP;
////            TPZVec<STATE> ErroF;
////            ErrorH1(cmesh2, ErroP /*,myerrorfile*/);
////            ErrorHDiv(cmesh1, ErroF /*,myerrorfile*/);
////            
////            REAL totalbanda = NDoFCond*NDoFCond;
////            REAL NumZeros = totalbanda - nNzeros;
////            REAL razao = NumZeros/totalbanda;
////            
////            
////            myerrorfile << ndiv <<  setw(13) << NDoF << setw(12) << NDoFCond << setw(13)<< NDoFCond*NDoFCond
////            << setw(15) << NumZeros << setw(12) << razao << "    ";
////            
////#ifdef USING_BOOST
////            myerrorfile << (t2-t1) << "     " << (t3-t2) << "     " << (t2-t1)+(t3-t2);
////#endif
////            
////            myerrorfile << setw(12) << ErroP[1] << setw(15) << ErroF[1] <<std::endl;
////            
//////            myerrorfile_fluxo << cmesh2->Solution() << std::endl;
////            myerrorfile_fluxo << cmesh1->Solution() << std::endl;
////            
////            
////            //------------;
//////            cmesh1->CleanUp();
//////            cmesh2->CleanUp();
//////            delete cmesh1;
//////            delete cmesh2;
//////            delete gmesh;
////        }
////        
////        myerrorfile <<"\n-------------------------------------------------------------------------"<<std::endl;
////    }
////    
////    return 0;
////}
//
//
//TPZGeoMesh * CreateOneCuboWithTetraedrons(int ndiv)
//{
//    double dndiv = ndiv;
//    int nelem = (int) pow(2., dndiv);
//    
//    int tetraedra_2[6][4] = { {1,2,5,4}, {4,7,3,2}, {0,1,2,4}, {0,2,3,4}, {4,5,6,2},{4,6,7,2} };
//    
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    GenerateNodes(gmesh,nelem);
//    gmesh->SetDimension(3);
//    
//    for (int64_t i=0; i<nelem; i++) {
//        for (int64_t j=0; j<nelem; j++) {
//            for (int64_t k=0; k<nelem; k++) {
//                TPZManVector<int64_t,8> nodes(8,0);
//                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
//                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
//                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
//                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
//                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
//                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
//                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
//                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
//                
//                for (int el=0; el<6; el++)
//                {
//                    TPZManVector<int64_t,4> elnodes(4);
//                    int64_t index;
//                    for (int il=0; il<4; il++) {
//                        elnodes[il] = nodes[tetraedra_2[el][il]];
//                    }
//                    gmesh->CreateGeoElement(ETetraedro, elnodes, matId, index);
//                }
//            }
//        }
//    }
//    
//    gmesh->BuildConnectivity();
//    
//    // Boundary Conditions
//    const int numelements = gmesh->NElements();
//    //    const int bczMinus = -3, bczplus = -2, bcids = -1;
//    //    const int bczMinus = -1, bczplus = -1, bcids = -1;
//    
//    for(int el=0; el<numelements; el++)
//    {
//        TPZManVector <TPZGeoNode,4> Nodefinder(4);
//        TPZManVector <REAL,3> nodecoord(3);
//        TPZGeoEl *tetra = gmesh->ElementVec()[el];
//        TPZVec<int64_t> ncoordVec(0);
//        int64_t sizeOfVec = 0;
//        
//        // na face z = 0
//        for (int i = 0; i < 4; i++)
//        {
//            int64_t pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[2],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc0);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face y = 0
//        for (int i = 0; i < 4; i++)
//        {
//            int64_t pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[1],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc1);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face x = 1
//        for (int i = 0; i < 4; i++)
//        {
//            int64_t pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[0],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc2);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face y = 1
//        for (int i = 0; i < 4; i++)
//        {
//            int64_t pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[1],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc3);
//        }
//        
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face x = 0
//        for (int i = 0; i < 4; i++)
//        {
//            int64_t pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[0],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc4);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face z = 1
//        for (int i = 0; i < 4; i++)
//        {
//            int64_t pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[2],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc5);
//        }
//    }
//    
//    return gmesh;
//}
//
//TPZGeoMesh *CreateOneCubo(int ndiv)
//{
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    int nnodes = 8;
//    
//    gmesh->SetDimension(3);
//    gmesh->NodeVec().Resize(nnodes);
//    gmesh->SetMaxNodeId(8);
//    
//    TPZManVector<REAL,3> coord(3,0.);
//    int in = 0;
//    
//    //cubo [0,1]ˆ3
//    //c0
//    coord[0] = 0.0;
//    coord[1] = 0.0;
//    coord[2] = 0.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c1
//    coord[0] = 1.0;
//    coord[1] = 0.0;
//    coord[2] = 0.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c2
//    coord[0] = 1.0;
//    coord[1] = 1.0;
//    coord[2] = 0.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c3
//    coord[0] = 0.0;
//    coord[1] = 1.0;
//    coord[2] = 0.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    
//    //c4
//    coord[0] = 0.0;
//    coord[1] = 0.0;
//    coord[2] = 1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c5
//    coord[0] = 1.0;
//    coord[1] = 0.0;
//    coord[2] = 1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c6
//    coord[0] = 1.0;
//    coord[1] = 1.0;
//    coord[2] = 1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c7
//    coord[0] = 0.0;
//    coord[1] = 1.0;
//    coord[2] = 1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    gmesh->SetMaxElementId(in);
//    
//    int index = 0;
//    
//    TPZVec<int64_t> TopologyQuad(4);
//    
//    // bottom
//    TopologyQuad[0] = 0;
//    TopologyQuad[1] = 1;
//    TopologyQuad[2] = 2;
//    TopologyQuad[3] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc0,*gmesh);
//    index++;
//    
//    // Front
//    TopologyQuad[0] = 0;
//    TopologyQuad[1] = 1;
//    TopologyQuad[2] = 5;
//    TopologyQuad[3] = 4;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc1,*gmesh);
//    index++;
//    
//    // Rigth
//    TopologyQuad[0] = 1;
//    TopologyQuad[1] = 2;
//    TopologyQuad[2] = 6;
//    TopologyQuad[3] = 5;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc2,*gmesh);
//    index++;
//    
//    // Back
//    TopologyQuad[0] = 3;
//    TopologyQuad[1] = 2;
//    TopologyQuad[2] = 6;
//    TopologyQuad[3] = 7;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc3,*gmesh);
//    index++;
//    
//    // Left
//    TopologyQuad[0] = 0;
//    TopologyQuad[1] = 3;
//    TopologyQuad[2] = 7;
//    TopologyQuad[3] = 4;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc4,*gmesh);
//    index++;
//    
//    // Top
//    TopologyQuad[0] = 4;
//    TopologyQuad[1] = 5;
//    TopologyQuad[2] = 6;
//    TopologyQuad[3] = 7;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc5,*gmesh);
//    index++;
//    
//    TPZManVector<int64_t,8> TopolCubo(8,0);
//    TopolCubo[0] = 0;
//    TopolCubo[1] = 1;
//    TopolCubo[2] = 2;
//    TopolCubo[3] = 3;
//    TopolCubo[4] = 4;
//    TopolCubo[5] = 5;
//    TopolCubo[6] = 6;
//    TopolCubo[7] = 7;
//    
//    
//    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index, TopolCubo, matId, *gmesh);
//    index++;
//    
//    
//    gmesh->BuildConnectivity();
//    
//    /// gmesh para aqui
//    
//    TPZVec<TPZGeoEl *> sons;
//    for (int iref = 0; iref < ndiv; iref++) {
//        int nel = gmesh->NElements();
//        for (int iel = 0; iel < nel; iel++) {
//            TPZGeoEl *gel = gmesh->ElementVec()[iel];
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            gel->Divide(sons);
//        }
//    }
//    
//    
//    //    std::ofstream out("SingleCubeWithBcs.vtk");
//    //    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
//    
//    return gmesh;
//}
//
////TPZGeoMesh *CreateOneCubo(int ndiv)
////{
////    TPZGeoMesh *gmesh = new TPZGeoMesh;
////    int nnodes = 8;
////    
////    gmesh->SetDimension(3);
////    gmesh->NodeVec().Resize(nnodes);
////    gmesh->SetMaxNodeId(8);
////    
////    TPZManVector<REAL,3> coord(3,0.);
////    int in = 0;
////    
////    //cubo [0,1]ˆ3
////    //c0
////    coord[0] = 0.0;
////    coord[1] = 0.0;
////    coord[2] = 0.0;
////    gmesh->NodeVec()[in].SetCoord(coord);
////    gmesh->NodeVec()[in].SetNodeId(in);
////    in++;
////    //c1
////    coord[0] = 1.0;
////    coord[1] = 0.0;
////    coord[2] = 0.0;
////    gmesh->NodeVec()[in].SetCoord(coord);
////    gmesh->NodeVec()[in].SetNodeId(in);
////    in++;
////    //c2
////    coord[0] = 1.0;
////    coord[1] = 1.0;
////    coord[2] = 0.0;
////    gmesh->NodeVec()[in].SetCoord(coord);
////    gmesh->NodeVec()[in].SetNodeId(in);
////    in++;
////    //c3
////    coord[0] = 0.0;
////    coord[1] = 1.0;
////    coord[2] = 0.0;
////    gmesh->NodeVec()[in].SetCoord(coord);
////    gmesh->NodeVec()[in].SetNodeId(in);
////    in++;
////    
////    //c4
////    coord[0] = 0.0;
////    coord[1] = 0.0;
////    coord[2] = 1.0;
////    gmesh->NodeVec()[in].SetCoord(coord);
////    gmesh->NodeVec()[in].SetNodeId(in);
////    in++;
////    //c5
////    coord[0] = 1.0;
////    coord[1] = 0.0;
////    coord[2] = 1.0;
////    gmesh->NodeVec()[in].SetCoord(coord);
////    gmesh->NodeVec()[in].SetNodeId(in);
////    in++;
////    //c6
////    coord[0] = 1.0;
////    coord[1] = 1.0;
////    coord[2] = 1.0;
////    gmesh->NodeVec()[in].SetCoord(coord);
////    gmesh->NodeVec()[in].SetNodeId(in);
////    in++;
////    //c7
////    coord[0] = 0.0;
////    coord[1] = 1.0;
////    coord[2] = 1.0;
////    gmesh->NodeVec()[in].SetCoord(coord);
////    gmesh->NodeVec()[in].SetNodeId(in);
////    in++;
////    gmesh->SetMaxElementId(in);
////    
////    int index = 0;
////    
////    TPZVec<int64_t> TopologyQuad(4);
////    
////    // bottom
////    TopologyQuad[0] = 0;
////    TopologyQuad[1] = 1;
////    TopologyQuad[2] = 2;
////    TopologyQuad[3] = 3;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc0,*gmesh);
////    index++;
////    
////    // Front
////    TopologyQuad[0] = 0;
////    TopologyQuad[1] = 1;
////    TopologyQuad[2] = 5;
////    TopologyQuad[3] = 4;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc1,*gmesh);
////    index++;
////    
////    // Rigth
////    TopologyQuad[0] = 1;
////    TopologyQuad[1] = 2;
////    TopologyQuad[2] = 6;
////    TopologyQuad[3] = 5;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc2,*gmesh);
////    index++;
////    
////    // Back
////    TopologyQuad[0] = 3;
////    TopologyQuad[1] = 2;
////    TopologyQuad[2] = 6;
////    TopologyQuad[3] = 7;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc3,*gmesh);
////    index++;
////    
////    // Left
////    TopologyQuad[0] = 0;
////    TopologyQuad[1] = 3;
////    TopologyQuad[2] = 7;
////    TopologyQuad[3] = 4;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc4,*gmesh);
////    index++;
////    
////    // Top
////    TopologyQuad[0] = 4;
////    TopologyQuad[1] = 5;
////    TopologyQuad[2] = 6;
////    TopologyQuad[3] = 7;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,bc5,*gmesh);
////    index++;
////    
////    TPZManVector<int64_t,8> TopolCubo(8,0);
////    TopolCubo[0] = 0;
////    TopolCubo[1] = 1;
////    TopolCubo[2] = 2;
////    TopolCubo[3] = 3;
////    TopolCubo[4] = 4;
////    TopolCubo[5] = 5;
////    TopolCubo[6] = 6;
////    TopolCubo[7] = 7;
////    
////    
////    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index, TopolCubo, matId, *gmesh);
////    index++;
////    
////    
////    gmesh->BuildConnectivity();
////    
////    /// gmesh para aqui
////    
////    TPZVec<TPZGeoEl *> sons;
////    for (int iref = 0; iref < ndiv; iref++) {
////        int nel = gmesh->NElements();
////        for (int iel = 0; iel < nel; iel++) {
////            TPZGeoEl *gel = gmesh->ElementVec()[iel];
////            if (gel->HasSubElement()) {
////                continue;
////            }
////            gel->Divide(sons);
////        }
////    }
////    
////    
////    //    std::ofstream out("SingleCubeWithBcs.vtk");
////    //    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
////    
////    return gmesh;
////}
//
//
//void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem)
//{
//    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
//    for (int64_t i=0; i<=nelem; i++) {
//        for (int64_t j=0; j<=nelem; j++) {
//            for (int64_t k=0; k<=nelem; k++) {
//                TPZManVector<REAL,3> x(3);
//                x[0] = k*1./nelem;
//                x[1] = j*1./nelem;
//                x[2] = i*1./nelem;
//                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
//            }
//        }
//    }
//}
//
//bool MyDoubleComparer(REAL a, REAL b)
//{
//    if (IsZero(a-b)){
//        return true;
//    }
//    else{
//        return false;
//    }
//}
//
//
//TPZGeoMesh *GMesh2D(bool ftriang)
//{
//    
//    if(dim_problema!=2)
//    {
//        std::cout << "dimensao errada" << std::endl;
//        DebugStop();
//    }
//    
//    int Qnodes =  4;
//    
//    TPZGeoMesh * gmesh = new TPZGeoMesh;
//    gmesh->SetMaxNodeId(Qnodes-1);
//    gmesh->NodeVec().Resize(Qnodes);
//    TPZVec<TPZGeoNode> Node(Qnodes);
//    
//    gmesh->SetDimension(dim_problema);
//    
//    TPZVec <int64_t> TopolQuad(4);
//    TPZVec <int64_t> TopolTriang(3);
//    TPZVec <int64_t> TopolLine(2);
//    
//    //indice dos nos
//    int64_t id = 0;
//    
//    TPZManVector<REAL,3> coord(2,0.);
//    int in = 0;
//    //c0
//    coord[0] = 0.0;
//    coord[1] = 0.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c1
//    coord[0] =  1.0;
//    coord[1] = 0.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c2
//    coord[0] =  1.0;
//    coord[1] =  1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c3
//    coord[0] = 0.0;
//    coord[1] =  1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //indice dos elementos
//    id = 0;
//    
//    if(ftriang==true) // triangulo
//    {
//        TopolTriang[0] = 0;
//        TopolTriang[1] = 1;
//        TopolTriang[2] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
//        id++;
//        
//        TopolTriang[0] = 0;
//        TopolTriang[1] = 2;
//        TopolTriang[2] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
//        id++;
//        
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//        id++;
//        
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//        id++;
//        
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//        id++;
//        
//        TopolLine[0] = 3;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
//        id++;
//    }
//    else{
//        
//        TopolQuad[0] = 0;
//        TopolQuad[1] = 1;
//        TopolQuad[2] = 2;
//        TopolQuad[3] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//        id++;
//        
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//        id++;
//        
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//        id++;
//        
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//        id++;
//        
//        TopolLine[0] = 3;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
//    }
//    
//    gmesh->BuildConnectivity();
//    
//    return gmesh;
//}
//
//
//void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
//{
//    for(int D = 0; D < nDiv; D++)
//    {
//        int nels = gmesh->NElements();
//        for(int elem = 0; elem < nels; elem++)
//        {
//            TPZVec< TPZGeoEl * > filhos;
//            TPZGeoEl * gel = gmesh->ElementVec()[elem];
//            gel->Divide(filhos);
//        }
//    }
//    gmesh->BuildConnectivity();
//}
//
//
////refinamento uniforme em direcao ao no
//void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide){
//    
//    ///Refinamento
//    //    gRefDBase.InitializeUniformRefPattern(EOned);
//    //    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    
//    
//    for (int idivide = 0; idivide < divide; idivide++){
//        const int nels = gmesh->NElements();
//        TPZVec< TPZGeoEl * > allEls(nels);
//        for(int iel = 0; iel < nels; iel++){
//            allEls[iel] = gmesh->ElementVec()[iel];
//        }
//        
//        for(int iel = 0; iel < nels; iel++){
//            TPZGeoEl * gel = allEls[iel];
//            if(!gel) continue;
//            if(gel->HasSubElement()) continue;
//            int nnodes = gel->NNodes();
//            int found = -1;
//            for(int in = 0; in < nnodes; in++){
//                if(gel->NodePtr(in)->Id() == nodeAtOriginId){
//                    found = in;
//                    break;
//                }
//            }///for in
//            if(found == -1) continue;
//            
//            MElementType gelT = gel->Type();
//            TPZAutoPointer<TPZRefPattern> uniform = gRefDBase.GetUniformRefPattern(gelT);
//            if(!uniform){
//                DebugStop();
//            }
//            gel->SetRefPattern(uniform);
//            TPZVec<TPZGeoEl*> filhos;
//            gel->Divide(filhos);
//            
//        }///for iel
//    }//idivide
//    
//    gmesh->BuildConnectivity();
//}///void
//
//void PrefinamentoRibsHybridMesh(TPZCompMesh *cmesh){
//    const int nel = cmesh->NElements();
//    for(int iel = 0; iel < nel; iel++){
//        TPZCompEl * cel = cmesh->ElementVec()[iel];
//        if(!cel)continue;
//        TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
//        if(!face)continue;
//        TPZInterpolationSpace *TwoD = NULL;
//        TPZInterpolationSpace *Rib = NULL;
//        if(face->LeftElement()->Reference()->Dimension()==2){
//            TwoD = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
//            Rib = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
//        }
//        else{
//            Rib = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
//            TwoD = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
//        }
//        if(!Rib || !TwoD) DebugStop();
//        int porder2D = TwoD->MaxOrder();
//        int pOrder1D = Rib->MaxOrder();
//        if(pOrder1D<0) continue;//elemento de contorno
//        int neworder = pOrder1D > porder2D ? pOrder1D : porder2D;
//        Rib->PRefine(neworder);
//    }
//    
//}
//
//void SetPOrderRibsHybridMesh(TPZCompMesh *cmesh, int porder){
//    const int nel = cmesh->NElements();
//    for(int iel = 0; iel < nel; iel++){
//        TPZCompEl * cel = cmesh->ElementVec()[iel];
//        if(!cel)continue;
//        TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
//        if(!face)continue;
//        //TPZInterpolationSpace *TwoD = NULL;
//        TPZInterpolationSpace *Rib = NULL;
//        if(face->LeftElement()->Reference()->Dimension()==2){
//            //TwoD = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
//            Rib = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
//        }
//        else{
//            Rib = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
//            //TwoD = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
//        }
//        if(!Rib) DebugStop();
//        Rib->PRefine(porder);
//    }
//    
//}
//
//void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder){
//    if(ndiv<1) return;
//    
//    if(!rodarHdiv && !rodarH1)
//    {
//        cmesh->Reference()->ResetReference();
//    }
//    int nel = cmesh->NElements();
//    for(int iel = 0; iel < nel; iel++)
//    {
//        TPZCompEl *cel = cmesh->ElementVec()[iel];
//        if(!cel) continue;
//        
//        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
//        if(!sp) continue;
//        int level = sp->Reference()->Level();
//        TPZGeoEl * gel = sp->Reference();
//        if(gel->Dimension()==dim_problema)
//        {
//            if(dim_problema==2)
//            {
//                sp->PRefine(porder + (level-flevel) + (ndiv-1));
//            }else
//            {
//                //sp->PRefine(porder + (level-flevel)+ (ndiv-1));
//                sp->PRefine(porder + (level-flevel));
//            }
//        }
//    }
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    cmesh->ExpandSolution();
//}
//
//void ExactSolution(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
//    
//#ifdef SolutionPoly
//    PolyProblem(loc, u, du);
//    return;
//#endif
//    
//#ifdef SolutionShock
//    SolShockProblem(loc, u, du);
//    return;
//#endif
//    
//}
//void ForcingFunction(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &df){
//    
//#ifdef SolutionPoly
//    ForcingPolyProblem(pt, disp, df);
//    return;
//#endif
//
//#ifdef SolutionShock
//    ForcingShockProblem(pt, disp, df);
//    return;
//#endif
//    
//}
//
//void ForcingFunctionII(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
//    
//    TPZFMatrix<STATE> df;
//#ifdef SolutionPoly
//    ForcingPolyProblem(pt, disp, df);
//    return;
//#endif
//    
//#ifdef SolutionShock
//    ForcingShockProblem(pt, disp, df);
//    return;
//#endif
//    
//    
//}
//
//void SolShockProblem(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
//    
//    REAL x = loc[0];
//    REAL y = loc[1];
//    REAL z = loc[2];
//    
//    REAL x0 = 1.25, y0 = -0.25, z0 = -0.25;
//    REAL r0 = M_PI/3.;
//    REAL alpha = 0.;
//    
//    u.Resize(1, 0.);
//    du.Resize(3, 1);
//    du(0,0)=0., du(1,0)=0., du(2,0)=0.;
//    
//    REAL temp1;
//    REAL temp2;
//    REAL r;
//    REAL grad;
//    temp1 = 0.,temp2=0., r=0., grad=0.;
//    
//    if(dim_problema == 2)
//    {
//        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0);
//        r = sqrt(temp1);
//        alpha = 200.;
//        temp2 = (r - r0)*alpha;
//        u[0] = M_PI/2. - atan(temp2);
//        
//        temp2 =  r*(1./alpha + (r - r0)*(r - r0)*alpha);
//        temp1 = (x0-x);
//        grad=temp1/temp2;
//        if(rodarHdiv) grad *= -1.;
//        du(0,0) = grad;
//        
//        temp1 = (y0-y);
//        grad=temp1/temp2;
//        if(rodarHdiv) grad *= -1.;
//        du(1,0)=grad;
//    }
//    else
//    {
//        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
//        r = sqrt(temp1);
//        alpha = alpha_param;//5.;
//        temp2 = (r - r0)*alpha;
//        u[0] = M_PI/2. - atan(temp2);
//        
//        temp2 =  r*(1./alpha + (r - r0)*(r - r0)*alpha);
//        temp1 = (x0-x);
//        grad=temp1/temp2;
//        if(rodarHdiv) grad *= -1.;
//        du(0,0) = grad;
//        
//        temp1 = (y0-y);
//        grad=temp1/temp2;
//        if(rodarHdiv) grad *= -1.;
//        du(1,0) =grad;;
//        
//        temp1 = (z0-z);
//        grad=temp1/temp2;
//        if(rodarHdiv) grad *= -1.;
//        du(2,0) =grad;
//    }
//}
//
//void ForcingShockProblem(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &df)
//{
//    REAL x = pt[0];
//    REAL y = pt[1];
//    REAL z = pt[2];
//    
//    REAL x0 = 1.25, y0 = -0.25, z0 = -0.25;
//    REAL r0 = M_PI/3.;
//    REAL alpha = 0.;
//    
//    REAL temp1=0., temp2=0.,temp3=0.,r=0., sol=0.;
//    
//    if (dim_problema==2)
//    {
//        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0);
//        r = sqrt(temp1);
//        alpha = 200.;
//        
//        temp1 = (1./alpha) + (r0*r0 - r*r)*alpha;
//        temp2 = r*((r - r0)*(r - r0) + 1./(alpha*alpha));
//        temp3 = 1. + (r - r0)*(r - r0)*alpha*alpha;
//        
//        sol = temp1/(temp2*temp3);
//        
////        if(rodarH1 || rodarSIPGD){
////            sol *=-1.;
////        }
//        
//        disp[0] = sol;
//    }
//    else
//    {
//        temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
//        r = sqrt(temp1);
//        alpha = alpha_param;//5.;
//        
//        temp1 = (2./alpha) + 2.*r0*(r0-r)*alpha;
//        temp2 = r*((r - r0)*(r - r0) + 1./(alpha*alpha));
//        temp3 = 1. + (r - r0)*(r - r0)*alpha*alpha;
//        
//        sol = temp1/(temp2*temp3);
//       //if(rodarH1 || rodarSIPGD) sol *=-1.;
//        disp[0] = sol;
//    }
//    
//}
//
//// Polynomial
//void PolyProblem(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
//    
//    REAL x = pt[0];
//    REAL y = pt[1];
//    REAL z = pt[2];
//    
//    du.Resize(3, 1);
//    
//    u[0] = x*x*x + y*y*y + z*z*z;
//    REAL dudx, dudy, dudz;
//    
//    dudx = 3.0*x*x;
//    dudy = 3.0*y*y;
//    dudz = 3.0*z*z;
//    
//    du(0,0) = -dudx;
//    du(1,0) = -dudy;
//    du(2,0) = -dudz;
//    
//    return;
//}
//void ForcingPolyProblem(const TPZVec<REAL> &pt, TPZVec<STATE> &f, TPZFMatrix<STATE> &df){
//    
//    REAL x = pt[0];
//    REAL y = pt[1];
//    REAL z = pt[2];
//
//    
//    f[0] = -6.0*x - 6.0*y - 6.0*z;
//    return;
//}
//
//void ForcingShockProblem2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
//    
//    TPZFMatrix<STATE> df;
//#ifdef SolutionPoly
//    ForcingPolyProblem(pt, disp, df);
//    return;
//#endif
//    
//#ifdef SolutionShock
//    ForcingShockProblem(pt, disp, df);
//    return;
//#endif
//    
//
//}
//
//void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    TPZFMatrix<STATE> du(3,1);
//    
//#ifdef SolutionPoly
//    PolyProblem(loc,result,du);
//    return;
//#endif
//    
//#ifdef SolutionShock
//    SolShockProblem(loc,result,du);
//    return;
//#endif
//
//}
//
//void NeumannBC1(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    
//    REAL normal[3] = {0.,-1.,0.};
//    TPZManVector<STATE> u(1);
//    TPZFNMatrix<5,STATE> du(3,1);
//    
//#ifdef SolutionPoly
//    PolyProblem(loc,u,du);
//#endif
//    
//#ifdef SolutionShock
//    SolShockProblem(loc,u,du);
//#endif
//    
//    result.Resize(1);
//    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
//}
//
//void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    
//    REAL normal[3] = {1.,0.,0.};
//    TPZManVector<STATE> u(1);
//    TPZFNMatrix<5,STATE> du(3,1);
//    
//#ifdef SolutionPoly
//    PolyProblem(loc,u,du);
//#endif
//    
//#ifdef SolutionShock
//    SolShockProblem(loc,u,du);
//#endif
//    
//    result.Resize(1);
//    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
//}
//
//void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    REAL normal[3] = {0.,0.,1.};
//    TPZManVector<STATE> u(1);
//    TPZFNMatrix<5,STATE> du(3,1);
//    
//#ifdef SolutionPoly
//    PolyProblem(loc,u,du);
//#endif
//    
//#ifdef SolutionShock
//    SolShockProblem(loc,u,du);
//#endif
//    
//    result.Resize(1);
//    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
//}
//
//void NeumannAbaixo(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    REAL normal[3] = {0.,0.,-1.};
//    TPZManVector<STATE> u(1);
//    TPZFNMatrix<5,STATE> du(3,1);
//    
//    
//#ifdef SolutionPoly
//    PolyProblem(loc,u,du);
//#endif
//    
//#ifdef SolutionShock
//    SolShockProblem(loc,u,du);
//#endif
//    
//    result.Resize(1);
//    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
//}
//
//void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    TPZFNMatrix<10,REAL> fake(2,1);
//    result[0] = loc[0]*loc[0];
//}
//
////----------------------------------
//TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder, bool ismultiplierH1, int ndiv){
//    
//    //TPZCompEl::SetgOrder(porder);
//    TPZCompMesh *comp = new TPZCompMesh(&gmesh);
//    
//    comp->SetDimModel(dim_problema);
//    comp->SetDefaultOrder(porder);
//    
//    comp->ApproxSpace().CreateDisconnectedElements(true);
//    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
//    
//    
//    // Criar e inserir os materiais na malha
//    REAL beta = 6;
//    TPZMatDualHybridPoisson *mymaterial = new TPZMatDualHybridPoisson(matId,0.,beta);
//    TPZMaterial * automat(mymaterial);
//    comp->InsertMaterialObject(automat);
//    mymaterial->SetDimension(dim_problema);
//    
//    //Condicoes de contorno
//    TPZFMatrix<STATE> val1(2,2,0.),val2(1,1,0.);
//    
//    int int_order = 20;
//    if(ndiv>1){
//        int_order = 10;
//    }
//    
//    //vetor de carga: lada direita da equacao
//    TPZAutoPointer<TPZFunction<STATE> > forcefunction;
//    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingShockProblem);
//    dum->SetPolynomialOrder(int_order);
//    forcefunction = dum;
//    mymaterial->SetForcingFunction(forcefunction);
//    
//    //solucao exata
//    TPZAutoPointer<TPZFunction<STATE> > solexata;
//    solexata = new TPZDummyFunction<STATE>(SolShockProblem);
//    mymaterial->SetForcingFunctionExact(solexata);
//    
//    
//    //bc1
//    TPZMaterial *bnd = automat->CreateBC(automat, bc1, 1, val1, val2);
//    bnd->SetForcingFunction(NeumannBC1,int_order);
//    comp->InsertMaterialObject(bnd);
//    
//    //bc2
//    bnd = automat->CreateBC (automat, bc2, 1, val1, val2);
//    bnd->SetForcingFunction(NeumannBC2,int_order);
//    comp->InsertMaterialObject(bnd);
//    
//    //    //bc1
//    //    TPZMaterial *bnd = automat->CreateBC(automat, bc1, 0, val1, val2);
//    //    bnd->SetForcingFunction(Dirichlet);
//    //    comp->InsertMaterialObject(bnd);
//    //
//    //    //bc2
//    //    bnd = automat->CreateBC (automat, bc2, 0, val1, val2);
//    //    bnd->SetForcingFunction(Dirichlet);
//    //    comp->InsertMaterialObject(bnd);
//    
//    //bc3
//    bnd = automat->CreateBC (automat, bc3, 0, val1, val2);
//    bnd->SetForcingFunction(Dirichlet,int_order);
//    comp->InsertMaterialObject(bnd);
//    
//    //bc4
//    bnd = automat->CreateBC (automat, bc4, 0, val1, val2);
//    bnd->SetForcingFunction(Dirichlet,int_order);
//    comp->InsertMaterialObject(bnd);
//    
//    if(dim_problema==3)
//    {
//        //bc0
//        bnd = automat->CreateBC (automat, bc0, 0, val1, val2);
//        bnd->SetForcingFunction(Dirichlet,int_order);
//        comp->InsertMaterialObject(bnd);
//        
//        //bc5
//        bnd = automat->CreateBC (automat, bc5, 1, val1, val2);
//        bnd->SetForcingFunction(NeumannAcima,int_order);
//        comp->InsertMaterialObject(bnd);
//        
//        //bc5
//        //        bnd = automat->CreateBC (automat, bc5, 0, val1, val2);
//        //        bnd->SetForcingFunction(Dirichlet);
//        //        comp->InsertMaterialObject(bnd);
//    }
//    
//    
//    // Ajuste da estrutura de dados computacional
//    comp->ApproxSpace().CreateDisconnectedElements(true);
//    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
//    
//    std::set<int> matids;
//    matids.insert(matId);
//    comp->AutoBuild(matids);
//    
//    comp->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
//    matids.clear();
//    for (int i=1; i<=nbc; i++)
//    {
//        if(dim_problema==2){
//            matids.insert(-(i+1));
//        }else{
//            matids.insert(-i);
//        }
//    }
//    
//    comp->SetDefaultOrder(porder);
//    comp->AutoBuild(matids);
//    comp->SetDimModel(dim_problema);
//    
//    comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
//    comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
//    
//    comp->LoadReferences();
//    comp->ApproxSpace().CreateInterfaceElements(comp,true);
//    
//    //
//    //    matids.insert(1);
//    //    if(!fsolsuave){
//    //        for (int i=1; i<=nbc; i++) {
//    //            matids.insert(-i);
//    //        }
//    //    }else{
//    //        matids.insert(-1);
//    //    }
//    //
//    //    comp->ApproxSpace().Hybridize(*comp, matids, ismultiplierH1);
//    comp->AdjustBoundaryElements();
//    return comp;
//}
//
//
//void GroupElements(TPZCompMesh *cmesh, int dimproblema)
//{
//    cmesh->LoadReferences();
//    int nel = cmesh->NElements();
//    std::set<TPZCompEl *> celset;
//    for (int el=0; el<nel; el++) {
//        TPZCompEl *cel = cmesh->ElementVec()[el];
//        if (!cel) {
//            continue;
//        }
//        TPZGeoEl *gel = cel->Reference();
//        int dim = gel->Dimension();
//        if (dim ==dimproblema) {
//            celset.insert(cel);
//        }
//    }
//    std::set<int> elgroupindices;
//    
//    for (std::set<TPZCompEl *>::iterator it = celset.begin(); it != celset.end(); it++) {
//        
//        std::list<TPZCompEl *> group;
//        group.push_back(*it);
//        TPZCompEl *cel = *it;
//        TPZGeoEl *gel = cel->Reference();
//        int nsides = gel->NSides();
//        for (int is = 0; is<nsides; is++) {
//            if (gel->SideDimension(is) != dimproblema-1) {
//                continue;
//            }
//            TPZStack<TPZCompElSide> connected;
//            TPZCompElSide celside(cel,is);
//            celside.EqualLevelElementList(connected, false, false);
//            int neq = connected.NElements();
//            for (int eq=0; eq<neq; eq++) {
//                TPZCompElSide eqside = connected[eq];
//                TPZCompEl *celeq = eqside.Element();
//                TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(celeq);
//                if (!intface) {
//                    continue;
//                }
//                TPZCompEl *left = intface->LeftElement();
//                if (left == cel) {
//                    //put in the group
//                    group.push_back(intface);
//                }
//            }
//        }
//        
//        //#ifdef LOG4CXX
//        //        if (logger->isDebugEnabled())
//        //        {
//        //            std::stringstream sout;
//        //            for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
//        //                (*it)->Print(sout);
//        //            }
//        //            LOGPZ_DEBUG(logger, sout.str())
//        //        }
//        //#endif
//        
//        int64_t index;
//        TPZElementGroup *celgroup = new TPZElementGroup(*cmesh,index);
//        elgroupindices.insert(index);
//        for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
//            celgroup->AddElement(*it);
//        }
//    }
//    cmesh->ComputeNodElCon();
//    
//    for (std::set<int>::iterator it = elgroupindices.begin(); it!=elgroupindices.end(); it++) {
//        TPZCompEl *cel = cmesh->ElementVec()[*it];
//        TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel);
//    }
//    
//    cmesh->CleanUpUnconnectedNodes();
//    cmesh->ExpandSolution();
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        cmesh->Print(sout);
//        LOGPZ_DEBUG(logger, sout.str())
//    }
//#endif
//    
//    
//    //    cmesh->LoadReferences();
//    //    int nel = cmesh->NElements();
//    //    std::set<TPZCompEl *> celset;
//    //    for (int el=0; el<nel; el++) {
//    //        TPZCompEl *cel = cmesh->ElementVec()[el];
//    //        if (!cel) {
//    //            continue;
//    //        }
//    //        TPZGeoEl *gel = cel->Reference();
//    //        int dim = gel->Dimension();
//    //        if (dim ==2) {
//    //            celset.insert(cel);
//    //        }
//    //    }
//    //    std::set<int> elgroupindices;
//    //
//    //    for (std::set<TPZCompEl *>::iterator it = celset.begin(); it != celset.end(); it++) {
//    //
//    //        std::set<TPZCompEl *> group;
//    //        group.insert(*it);
//    //        TPZCompEl *cel = *it;
//    //        TPZGeoEl *gel = cel->Reference();
//    //        int nsides = gel->NSides();
//    //        for (int is = 0; is<nsides; is++) {
//    //            if (gel->SideDimension(is) != 1) {
//    //                continue;
//    //            }
//    //            TPZStack<TPZCompElSide> connected;
//    //            TPZCompElSide celside(cel,is);
//    //            celside.EqualLevelElementList(connected, false, false);
//    //
//    //            if(connected.size()==0){
//    //                celside.HigherLevelElementList(connected, 0, 0);
//    //            }
//    //
//    //            int neq = connected.NElements();
//    //            for (int eq=0; eq<neq; eq++) {
//    //                TPZCompElSide eqside = connected[eq];
//    //                TPZCompEl *celeq = eqside.Element();
//    //                TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(celeq);
//    //                if (!intface) {
//    //                    continue;
//    //                }
//    //                TPZCompEl *left = intface->LeftElement();
//    //                if (left == cel) {
//    //                    //put in the group
//    //                    group.insert(intface);
//    //                }
//    //            }
//    //        }
//    //#ifdef LOG4CXX
//    //        if (logger->isDebugEnabled())
//    //        {
//    //            std::stringstream sout;
//    //            for (std::set<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
//    //                (*it)->Print(sout);
//    //            }
//    //            LOGPZ_DEBUG(logger, sout.str())
//    //        }
//    //#endif
//    //        int64_t index;
//    //        TPZElementGroup *celgroup = new TPZElementGroup(*cmesh,index);
//    //        elgroupindices.insert(index);
//    //        for (std::set<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
//    //            celgroup->AddElement(*it);
//    //        }
//    //    }
//    //    cmesh->ComputeNodElCon();
//    //
//    //    for (std::set<int>::iterator it = elgroupindices.begin(); it!=elgroupindices.end(); it++) {
//    //        TPZCompEl *cel = cmesh->ElementVec()[*it];
//    //        TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel);
//    //    }
//    //
//    //    cmesh->CleanUpUnconnectedNodes();
//    //    cmesh->ExpandSolution();
//    //#ifdef LOG4CXX
//    //    if (logger->isDebugEnabled())
//    //    {
//    //        std::stringstream sout;
//    //        cmesh->Print(sout);
//    //        LOGPZ_DEBUG(logger, sout.str())
//    //    }
//    //#endif
//    
//}
//
//
////solucao suave:sin*sin
//void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp,  TPZFMatrix<STATE> &df){
//    double x = pt[0];
//    double y = pt[1];
//    disp.Resize(1);
//    disp[0]=0.;
//    df.Redim(1, 1);
//    disp[0]= 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
//}
//
//void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
//    
//    const REAL x = loc[0];
//    const REAL y = loc[1];
//    u.Resize(1, 0.);
//    du.Resize(3, 1.);
//    du(0,0)=du(1,0)=du(2,0)=0.;
//    
//    const REAL sol = sin(M_PI*x)*sin(M_PI*y);
//    u[0] = sol;
//    
//    du(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y);
//    du(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x);
//}
//
//void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
//    
//    int i;
//    bool isdefined = false;
//    
//    // Refinando no local desejado
//    int npoints = 1000;
//    TPZVec<REAL> point(3);
//    point[0] = 1.; point[1] = 0.0; point[2] = 0.0;
//    REAL r = 0.72;
//    TPZVec<TPZManVector<REAL> > Points(npoints);
//    //GetPointsOnCircunference(npoints,point,r,Points);
//    
//    if(ntyperefs==2) {
//        REAL radius = 0.19;
//        for(i=0;i<nref;i+=2) {
//            // To refine elements with center near to points than radius
//            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
//            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
//            if(nref < 5) radius *= 0.35;
//            else if(nref < 7) radius *= 0.2;
//            else radius *= 0.1;
//        }
//        if(i==nref) {
//            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
//        }
//    }
//    else {
//        REAL radius;
//        if(flevel==3) radius= 0.18;//0.22;
//        if(flevel==2) radius= 0.18;;
//        for(i=0;i<nref;i++) {
//            // To refine elements with center near to points than radius
//            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
//            if(nref < 6 && flevel==2) radius *= 0.7;//0.6
//            if(nref < 6 && flevel==3) radius *= 0.7;//0.6
//            //else radius *= 0.7;
//            //            else if(nref < 7) radius *= 0.3;
//            //            else radius *= 0.1;
//        }
//    }
//    // Constructing connectivities
//    //  gmesh->ResetConnectivities();
//    RegularizeMesh(gmesh);
//    gmesh->BuildConnectivity();
//}
//
//void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL,3> > &Points) {
//    Points.Resize(npoints);
//    TPZManVector<REAL> point(3,0.);
//    REAL angle = (2*M_PI)/npoints;
//    for(int i=0;i<npoints;i++) {
//        point[0] = center[0]+radius*cos(i*angle);
//        point[1] = center[1]+radius*sin(i*angle);
//        Points[i] = point;
//    }
//}
//
//void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
//    TPZManVector<REAL> centerpsi(3), center(3);
//    // Refinamento de elementos selecionados
//    TPZGeoEl *gel;
//    TPZVec<TPZGeoEl *> sub;
//    
//    int nelem = 0;
//    int ngelem=gmesh->NElements();
//    // na esquina inferior esquerda Nó = (0,-1,0)
//    while(nelem<ngelem) {
//        gel = gmesh->ElementVec()[nelem++];
//        if(/*gel->Dimension()!=dim ||*/ gel->HasSubElement()) continue;
//        gel->CenterPoint(gel->NSides()-1,centerpsi);
//        gel->X(centerpsi,center);
//        if(!isdefined) {
//            TPZVec<REAL> FirstNode(3,0.);
//            gel->CenterPoint(0,centerpsi);
//            gel->X(centerpsi,FirstNode);
//            REAL distancia = TPZGeoEl::Distance(center,FirstNode);
//            if(distancia > distance) distance = distancia;
//            isdefined = true;
//        }
//        REAL centerdist = TPZGeoEl::Distance(center,point);
//        if(fabs(r-centerdist) < distance) {
////            std::cout << "Dividing gel index " << gel->Index() << " nsubel " << gel->NSubElements()  <<std::endl;
//            gel->Divide(sub);
//        }
//    }
//}
//
//
//void RegularizeMesh(TPZGeoMesh *gmesh)
//{
//    bool changed = true;
//    while (changed)
//    {
//        changed = false;
//        int nel = gmesh->NElements();
//        for (int64_t el=0; el<nel; el++) {
//            TPZGeoEl *gel = gmesh->ElementVec()[el];
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            int dim = gel->Dimension();
//            if (dim != dim_problema) {
//                continue;
//            }
//            int nsides = gel->NSides();
//            int nrefined = 0;
//            int nsidedim = 0;
//            for (int is=0; is<nsides; is++) {
//                int sidedim = gel->SideDimension(is);
//                if (sidedim != dim-1) {
//                    continue;
//                }
//                nsidedim++;
//            }
//            for (int is=0; is<nsides; is++) {
//                int sidedim = gel->SideDimension(is);
//                if (sidedim != dim-1) {
//                    continue;
//                }
//                TPZGeoElSide thisside(gel,is);
//                TPZGeoElSide neighbour = thisside.Neighbour();
//                if (neighbour != thisside) {
//                    TPZStack<TPZGeoElSide> subelements;
//                    neighbour.GetSubElements2(subelements);
//                    int nsub = subelements.size();
//                    if (nsub > 0) {
//                        nrefined++;
//                    }
//                    for (int isub=0; isub<nsub; isub++) {
//                        TPZGeoElSide sub = subelements[isub];
//                        if (sub.Dimension() != dim-1) {
//                            continue;
//                        }
//                        if (sub.HasSubElement()) {
//                            TPZManVector<TPZGeoEl *> newsub;
//                            gel->Divide(newsub);
//                            changed = true;
//                            break;
//                        }
//                    }
//                }
//                if (gel->HasSubElement()) {
//                    break;
//                }
//            }
//            if (nrefined >= nsidedim-1) {
//                TPZManVector<TPZGeoEl *> newsub;
//                gel->Divide(newsub);
//                changed = true;
//            }
//        }
//    }
//    
//    
//    gmesh->BuildConnectivity();
//    int nel = gmesh->NElements();
//    for(int ie = 0; ie<nel; ie++){
//        TPZGeoEl *gel = gmesh->ElementVec()[ie];
//        if(!gel || gel->HasSubElement()) continue;
//        
//        int matid = gel->MaterialId();
//        if(matid>0) continue;
//        
//        int nside = gel->NSides();
//        TPZGeoElSide thisside(gel,nside-1);
//        TPZGeoElSide neighbour = thisside.Neighbour();
//        
//        int dimel = gel->Dimension();
//        int dimn = neighbour.Element()->Dimension();
//        if(dimn==dimel) continue;
//        
//        int nsubel = neighbour.Element()->NSubElements();
//        for(int is=0; is<nsubel; is++){
//            TPZGeoEl *neigel = neighbour.Element()->SubElement(is);
//            if(neigel->HasSubElement()) continue;
//            
//            int mylevel  = gel->Level();
//            int neiglevel = neigel->Level();
//            
//            if(mylevel!=neiglevel){
//                TPZVec< TPZGeoEl * > filhos;
//                gel->Divide(filhos);
//                break;
//            }
//        }
//    }
//}
//
//TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim, bool rodarSIPGD, int ndiv)
//{
//    /// criar materiais
//    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim_problema);
//    TPZMatLaplacian *material = new TPZMatLaplacian(matId,dim_problema);
//    
//    TPZMaterial * mat(material);
//    material->NStateVariables();
//    
//    if(rodarSIPGD)
//    {
//        material->SetSymmetric();
//        material->SetSolutionPenalty();
//
////        material->SetNonSymmetric();
////        material->SetNoPenalty();
//    }
//    
//    //solucao exata
//    TPZAutoPointer<TPZFunction<STATE> > solexata;
//    solexata = new TPZDummyFunction<STATE>(SolShockProblem);
//    material->SetForcingFunctionExact(solexata);
//    
//    int int_order = 20;
//    if(ndiv > 3 || dim_problema==3){
//        int_order = 10;
//    }
//
//    
//    //funcao do lado direito da equacao do problema
//    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingShockProblem2);
//    TPZAutoPointer<TPZFunction<STATE> > forcef;
//    dum->SetPolynomialOrder(int_order);
//    forcef = dum;
//    material->SetForcingFunction(forcef);
//    
//    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//    cmesh->SetDimModel(dim_problema);
//    cmesh->SetDefaultOrder(pOrder);
//    
//    cmesh->InsertMaterialObject(mat);
//    
//    ///Inserir condicao de contorno
//    TPZFMatrix<STATE> val1(2,2,0.), val2(1,1,0.);
//    val2(0,0) = 0.0;
//    
//    //bc1
//    TPZMaterial *bnd = mat->CreateBC(mat, bc1, 0, val1, val2);
//     bnd->SetForcingFunction(Dirichlet,int_order);
//    //bnd->SetForcingFunction(NeumannBC1,int_order);
//    cmesh->InsertMaterialObject(bnd);
//    
//    
//    //bc2
//    bnd = mat->CreateBC (mat, bc2, 1, val1, val2);
//    bnd->SetForcingFunction(NeumannBC2,int_order);
//    cmesh->InsertMaterialObject(bnd);
//    
//    
//    //    //bc1
//    //    TPZMaterial *bnd = mat->CreateBC(mat, bc1, 0, val1, val2);
//    //    bnd->SetForcingFunction(Dirichlet);
//    //    cmesh->InsertMaterialObject(bnd);
//    //
//    //    //bc2
//    //    bnd = mat->CreateBC (mat, bc2, 0, val1, val2);
//    //    bnd->SetForcingFunction(Dirichlet);
//    //    cmesh->InsertMaterialObject(bnd);
//    
//    
//    //bc3
//    bnd = mat->CreateBC (mat, bc3, 0, val1, val2);
//    bnd->SetForcingFunction(Dirichlet,int_order);
//    cmesh->InsertMaterialObject(bnd);
//    
//    //bc4
//    bnd = mat->CreateBC (mat, bc4, 0, val1, val2);
//    bnd->SetForcingFunction(Dirichlet,int_order);
//    cmesh->InsertMaterialObject(bnd);
//    
//    
//    if(dim_problema==3)
//    {
//        //bc0
//        bnd = mat->CreateBC (mat, bc0, 0, val1, val2);
//        bnd->SetForcingFunction(Dirichlet,int_order);
//        cmesh->InsertMaterialObject(bnd);
//        
//        //bc5
//        bnd = mat->CreateBC (mat, bc5, 1, val1, val2);
//        bnd->SetForcingFunction(NeumannAcima,int_order);
//        cmesh->InsertMaterialObject(bnd);
//    }
//    
//    if(rodarSIPGD){
//        cmesh->SetAllCreateFunctionsDiscontinuous();
//    }else{
//        cmesh->SetAllCreateFunctionsContinuous();
//    }
//    
//    //Ajuste da estrutura de dados computacional
//    cmesh->AutoBuild();
//    
//    //    cout<<"\nNumero total de Equacoes: "<<cmesh->NEquations()<<"\n";
//    //    //condensar
//    //    for (int64_t iel=0; iel<cmesh->NElements(); iel++) {
//    //        TPZCompEl *cel = cmesh->Element(iel);
//    //        if(!cel) continue;
//    //        TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel);
//    //    }
//    //    cout<<"Numero Equacoes apos condensar: "<<cmesh->NEquations()<<"\n";
//    //
//    //    cmesh->ExpandSolution();
//    //    cmesh->CleanUpUnconnectedNodes();
//    
//    if(rodarSIPGD)
//    {
//        int64_t nel = cmesh->NElements();
//        for(int64_t i=0; i<nel; i++){
//            TPZCompEl *cel = cmesh->ElementVec()[i];
//            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//            if(!celdisc) continue;
//            //            celdisc->SetConstC(1.);
//            //            celdisc->SetCenterPoint(0, 0.);
//            //            celdisc->SetCenterPoint(1, 0.);
//            //            celdisc->SetCenterPoint(2, 0.);
//            celdisc->SetFalseUseQsiEta();
//            if (dim_problema!=2) continue;
//            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//            {
//                if(fTriang==true) celdisc->SetTotalOrderShape();
//                else celdisc->SetTensorialShape();
//            }
//        }
//    }
//    
//    
//    
//    return cmesh;
//    
//}
//
//
//TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh, bool hdivskeleton, int ndiv){
//    
//    //Creating computational mesh for multiphysic elements
//    gmesh->ResetReference();
//    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
//    
//    //criando material
//    int dim =gmesh->Dimension();
//    
//    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(matId,dim);
//    if(hdivskeleton) material->UseSecondIntegrationByParts();
//    //TPZMixedPoisson *material = new TPZMixedPoisson(matId,dim);
//    
//    //incluindo os dados do problema
//    //    REAL coefk = 1.;
//    //    material->SetPermeability(coefk);
//    REAL coefvisc = 1.;
//    material->SetViscosity(coefvisc);
//    
//    //permeabilidade
//    TPZFMatrix<REAL> Ktensor(3,3,0.);
//    TPZFMatrix<REAL> InvK(3,3,0.);
//    Ktensor(0,0)=1.; Ktensor(1,1)=1.; Ktensor(2,2)=1.;
//    InvK=Ktensor;
//    material->SetPermeabilityTensor(Ktensor,InvK);
//    
//    int int_order = 20;
//    if(ndiv > 3 || dim_problema==3){
//        int_order = 10;
//    }
//    //solucao exata
//    TPZAutoPointer<TPZFunction<STATE> > solexata;
//    solexata = new TPZDummyFunction<STATE>(ExactSolution);
//    material->SetForcingFunctionExact(solexata);
//    
//    //funcao do lado direito da equacao do problema
//    TPZAutoPointer<TPZFunction<STATE> > force;
//    TPZDummyFunction<STATE> *dum;
//    dum = new TPZDummyFunction<STATE>(ForcingFunctionII);
//    dum->SetPolynomialOrder(int_order);
//    force = dum;
//    material->SetForcingFunction(force);
//    
//    //inserindo o material na malha computacional
//    TPZMaterial *mat(material);
//    mphysics->InsertMaterialObject(mat);
//    mphysics->SetDimModel(dim);
//    
//    ///Inserir condicao de contorno
//    TPZFMatrix<STATE> val1(2,2,0.), val2(1,1,0.);
//    val2(0,0) = 0.0;
//    
//    //bc1
//    TPZMaterial *bnd = mat->CreateBC(mat, bc1, 0, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > bc1exata;
//    TPZDummyFunction<STATE> *dum1;
//    dum1 = new TPZDummyFunction<STATE>(Dirichlet);//NeumannBC1
//    //    dum1 = new TPZDummyFunction<STATE>(Dirichlet);
//    dum1->SetPolynomialOrder(int_order);
//    bc1exata = dum1;
//    bnd->SetForcingFunction(bc1exata);
//    mphysics->InsertMaterialObject(bnd);
//    
//    //bc2
//    bnd = mat->CreateBC (mat, bc2, 1, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > bc2exata;
//    TPZDummyFunction<STATE> *dum2;
//    dum2 = new TPZDummyFunction<STATE>(NeumannBC2);
//    //dum2 = new TPZDummyFunction<STATE>(Dirichlet);
//    dum2->SetPolynomialOrder(int_order);
//    bc2exata = dum2;
//    bnd->SetForcingFunction(bc2exata);
//    mphysics->InsertMaterialObject(bnd);
//    
//    //bc3
//    bnd = mat->CreateBC (mat, bc3, 0, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > bc3exata;
//    TPZDummyFunction<STATE> *dum3;
//    dum3 = new TPZDummyFunction<STATE>(Dirichlet);
//    dum3->SetPolynomialOrder(int_order);
//    bc3exata = dum3;
//    bnd->SetForcingFunction(bc3exata);
//    mphysics->InsertMaterialObject(bnd);
//    
//    //bc4
//    bnd = mat->CreateBC (mat, bc4, 0, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > bc4exata;
//    TPZDummyFunction<STATE> *dum4;
//    dum4 = new TPZDummyFunction<STATE>(Dirichlet);
//    dum4->SetPolynomialOrder(int_order);
//    bc4exata = dum4;
//    bnd->SetForcingFunction(bc4exata);
//    mphysics->InsertMaterialObject(bnd);
//    
//    if(dim_problema==3)
//    {
//        //bc0
//        //        bnd = mat->CreateBC (mat, bc0, 0, val1, val2);
//        //        bnd->SetForcingFunction(Dirichlet);
//        //        mphysics->InsertMaterialObject(bnd);
//        
//        bnd = mat->CreateBC (mat, bc0, 0, val1, val2);
//        TPZAutoPointer<TPZFunction<STATE> > bc0exata;
//        TPZDummyFunction<STATE> *dum0;
//        dum0 = new TPZDummyFunction<STATE>(Dirichlet);
//        dum0->SetPolynomialOrder(int_order);
//        bc0exata = dum0;
//        bnd->SetForcingFunction(bc0exata);
//        mphysics->InsertMaterialObject(bnd);
//        
//        //bc5
//        //        bnd = mat->CreateBC (mat, bc5, 1, val1, val2);
//        //        bnd->SetForcingFunction(NeumannAcima);
//        //        mphysics->InsertMaterialObject(bnd);
//        
//        bnd = mat->CreateBC (mat, bc5, 1, val1, val2);
//        TPZAutoPointer<TPZFunction<STATE> > bc5exata;
//        TPZDummyFunction<STATE> *dum5;
//        dum5 = new TPZDummyFunction<STATE>(NeumannAcima);
//        //dum5 = new TPZDummyFunction<STATE>(Dirichlet);
//        dum5->SetPolynomialOrder(int_order);
//        bc5exata = dum5;
//        bnd->SetForcingFunction(bc5exata);
//        mphysics->InsertMaterialObject(bnd);
//        
//    }
//    
//    mphysics->SetAllCreateFunctionsMultiphysicElem();
//    
//    //Fazendo auto build
//    mphysics->AutoBuild();
//    mphysics->AdjustBoundaryElements();
//    mphysics->CleanUpUnconnectedNodes();
//    
//    if(hdivskeleton==false){
//        // Creating multiphysic elements into mphysics computational mesh
//        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
//        
//        
////        TPZCompMeshTools::GroupElements(mphysics);
//        TPZCompMeshTools::CreatedCondensedElements(mphysics, true);
//        
//        
////        //Condensacao Estatica
////        mphysics->Reference()->ResetReference();
////        mphysics->LoadReferences();
////        mphysics->SetDimModel(dim);
////        
////        // create condensed elements
////        // increase the NumElConnected of one pressure connects in order to prevent condensation
////        mphysics->ComputeNodElCon();
////        for (int64_t icel=0; icel < mphysics->NElements(); icel++) {
////            TPZCompEl  * cel = mphysics->Element(icel);
////            if(!cel) continue;
////            int nc = cel->NConnects();
////            for (int ic=0; ic<nc; ic++) {
////                TPZConnect &c = cel->Connect(ic);
////                if (c.LagrangeMultiplier() > 0) {
////                    c.IncrementElConnected();
////                    break;
////                }
////            }
////            new TPZCondensedCompEl(cel);
////        }
//    }
//    
//    //Creating multiphysic elements containing skeletal elements.
//    else
//    {
//        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//        mphysics->Reference()->ResetReference();
//        mphysics->LoadReferences();
//        
//        int64_t nel = mphysics->ElementVec().NElements();
//        
//        std::map<int64_t, int64_t> bctoel, eltowrap;
//        for (int64_t el=0; el<nel; el++) {
//            TPZCompEl *cel = mphysics->Element(el);
//            TPZGeoEl *gel = cel->Reference();
//            int matid = gel->MaterialId();
//            if (matid < 0) {
//                TPZGeoElSide gelside(gel,gel->NSides()-1);
//                TPZGeoElSide neighbour = gelside.Neighbour();
//                while (neighbour != gelside) {
//                    if (neighbour.Element()->Dimension() == dim && neighbour.Element()->Reference()) {
//                        // got you!!
//                        bctoel[el] = neighbour.Element()->Reference()->Index();
//                        break;
//                    }
//                    neighbour = neighbour.Neighbour();
//                }
//                if (neighbour == gelside) {
//                    DebugStop();
//                }
//            }
//        }
//        
//        TPZStack< TPZStack< TPZMultiphysicsElement *,7> > wrapEl;
//        for(int64_t el = 0; el < nel; el++)
//        {
//            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(mphysics->Element(el));
//            if(mfcel->Dimension()==dim) TPZBuildMultiphysicsMesh::AddWrap(mfcel, matId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
//        }
//        
//        for (int64_t el =0; el < wrapEl.size(); el++) {
//            TPZCompEl *cel = wrapEl[el][0];
//            int64_t index = cel->Index();
//            eltowrap[index] = el;
//        }
//        
//        meshvec[0]->CleanUpUnconnectedNodes();
//        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
//        
//        std::map<int64_t, int64_t>::iterator it;
//        for (it = bctoel.begin(); it != bctoel.end(); it++) {
//            int64_t bcindex = it->first;
//            int64_t elindex = it->second;
//            if (eltowrap.find(elindex) == eltowrap.end()) {
//                DebugStop();
//            }
//            int64_t wrapindex = eltowrap[elindex];
//            TPZCompEl *bcel = mphysics->Element(bcindex);
//            TPZMultiphysicsElement *bcmf = dynamic_cast<TPZMultiphysicsElement *>(bcel);
//            if (!bcmf) {
//                DebugStop();
//            }
//            wrapEl[wrapindex].Push(bcmf);
//            
//        }
//        
//        //------- Create and add group elements -------
//        int64_t index, nenvel;
//        nenvel = wrapEl.NElements();
//        TPZStack<TPZElementGroup *> elgroups;
//        for(int64_t ienv=0; ienv<nenvel; ienv++){
//            TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
//            elgroups.Push(elgr);
//            nel = wrapEl[ienv].NElements();
//            for(int jel=0; jel<nel; jel++){
//                elgr->AddElement(wrapEl[ienv][jel]);
//            }
//        }
//        
//        mphysics->ComputeNodElCon();
//        // create condensed elements
//        // increase the NumElConnected of one pressure connects in order to prevent condensation
//        for (int64_t ienv=0; ienv<nenvel; ienv++) {
//            TPZElementGroup *elgr = elgroups[ienv];
//            int nc = elgr->NConnects();
//            for (int ic=0; ic<nc; ic++) {
//                TPZConnect &c = elgr->Connect(ic);
//                if (c.LagrangeMultiplier() > 0) {
//                    c.IncrementElConnected();
//                    break;
//                }
//            }
//            TPZCondensedCompEl *condense = new TPZCondensedCompEl(elgr);
//        }
//    }
//    
//    mphysics->CleanUpUnconnectedNodes();
//    mphysics->ExpandSolution();
//    
//    return mphysics;
//}
//
//TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh)
//{
//    /// criar materiais
//    int dim = gmesh->Dimension();
//    TPZMatPoisson3d *material;
//    material = new TPZMatPoisson3d(matId,dim);
//    TPZMaterial * mat(material);
//    material->NStateVariables();
//    
//    
//    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//    cmesh->SetDimModel(dim);
//    
//    cmesh->SetAllCreateFunctionsHDiv();
//    
//    cmesh->InsertMaterialObject(mat);
//    
//    ///Criar condicoes de contorno
//    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//    TPZMaterial * BCond1 = material->CreateBC(mat, bc1, 0, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(mat, bc2, 0, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(mat, bc3, 0, val1, val2);
//    TPZMaterial * BCond4 = material->CreateBC(mat, bc4, 0, val1, val2);
//    
//    if(dim_problema==3)
//    {
//        //bc0
//        TPZMaterial * BCond0 = material->CreateBC(mat, bc0, 0, val1, val2);
//        cmesh->InsertMaterialObject(BCond0);
//        
//        //bc5
//        TPZMaterial * BCond5 = material->CreateBC(mat, bc5, 0, val1, val2);
//        cmesh->InsertMaterialObject(BCond5);
//    }
//    
//    
//    ///Inserir condicoes de contorno
//    cmesh->InsertMaterialObject(BCond1);
//    cmesh->InsertMaterialObject(BCond2);
//    cmesh->InsertMaterialObject(BCond3);
//    cmesh->InsertMaterialObject(BCond4);
//    
//    cmesh->SetDefaultOrder(pOrder);
//    cmesh->SetDimModel(dim);
//    
//    //    TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(matskeleton, dim-1, 1);
//    //    TPZMaterial * mat2(matskelet);
//    //    cmesh->InsertMaterialObject(mat2);
//    
//    cmesh->AutoBuild();//Ajuste da estrutura de dados computacional
//    
//    //#ifdef LOG4CXX
//    //	if(logdata->isDebugEnabled())
//    //	{
//    //        std::stringstream sout;
//    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
//    //        cmesh->Print(sout);
//    //        LOGPZ_DEBUG(logdata,sout.str())
//    //	}
//    //#endif
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    return cmesh;
//}
//
//TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh)
//{
//    /// criar materiais
//    int dim = gmesh->Dimension();
//    TPZMatPoisson3d *material;
//    material = new TPZMatPoisson3d(matId,dim);
//    material->NStateVariables();
//    
//    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//    cmesh->SetDimModel(dim);
//    TPZMaterial * mat(material);
//    cmesh->InsertMaterialObject(mat);
//    
//    
//    cmesh->SetDefaultOrder(pOrder);
//    cmesh->SetDimModel(dim);
//    
//    
//    cmesh->SetAllCreateFunctionsContinuous();
//    cmesh->ApproxSpace().CreateDisconnectedElements(true);
//    //cmesh->SetAllCreateFunctionsDiscontinuous();
//    
//    //Ajuste da estrutura de dados computacional
//    cmesh->AutoBuild();
//    
//    int ncon = cmesh->NConnects();
//    for(int i=0; i<ncon; i++)
//    {
//        TPZConnect &newnod = cmesh->ConnectVec()[i];
//        //newnod.SetPressure(true);
//        newnod.SetLagrangeMultiplier(1);
//    }
//    
//    int nel = cmesh->NElements();
//    for(int i=0; i<nel; i++){
//        TPZCompEl *cel = cmesh->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        if(!celdisc) continue;
//        celdisc->SetConstC(1.);
//        celdisc->SetCenterPoint(0, 0.);
//        celdisc->SetCenterPoint(1, 0.);
//        celdisc->SetCenterPoint(2, 0.);
//        celdisc->SetTrueUseQsiEta();
//        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//        {
//            if(fTriang==true) celdisc->SetTotalOrderShape();
//            else celdisc->SetTensorialShape();
//        }
//        
//    }
//    
//    cmesh->CleanUpUnconnectedNodes();
//    return cmesh;
//}
//
//void ErrorHDiv(TPZCompMesh *hdivmesh, TPZVec<STATE> &Error /*,std::ostream &out*/)
//{
//    int64_t nel = hdivmesh->NElements();
//    int dim = hdivmesh->Dimension();
//    TPZManVector<REAL,10> globerrors(10,0.);
//    TPZStack<REAL> vech;
//    
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = hdivmesh->ElementVec()[el];
//        if (!cel) {
//            continue;
//        }
//        
//        TPZGeoEl *gel = cel->Reference();
//        if (!gel || gel->Dimension() != dim) {
//            continue;
//        }
//        TPZManVector<REAL,10> elerror(10,0.);
//        cel->EvaluateError(ExactSolution, elerror, 0);
//        int nerr = elerror.size();
//        for (int i=0; i<nerr; i++) {
//            globerrors[i] += elerror[i]*elerror[i];
//        }
//        
//    }
//    
//    //out<< sqrt(globerrors[1])<<"\n";
//    //out << "Errors associated with HDiv space\n";
//    //out << "\nL2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
//    //out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<endl;
//    //out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<endl;
//    
//    Error.Resize(3, 0.);
//    Error[1] = sqrt(globerrors[1]);
//}
//
//void ErrorH1(TPZCompMesh *l2mesh, TPZVec<STATE> &Error /*,std::ostream &out*/)
//{
//    int64_t nel = l2mesh->NElements();
//    int dim = l2mesh->Dimension();
//    TPZManVector<REAL,10> globerrors(10,0.);
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = l2mesh->ElementVec()[el];
//        if (!cel) {
//            continue;
//        }
//        TPZGeoEl *gel = cel->Reference();
//        if (!gel || gel->Dimension() != dim) {
//            continue;
//        }
//        TPZManVector<REAL,10> elerror(10,0.);
//        elerror.Fill(0.);
//        cel->EvaluateError(ExactSolution, elerror, 0);
//        
//        int nerr = elerror.size();
//        globerrors.resize(nerr);
//        
//        for (int i=0; i<nerr; i++) {
//            globerrors[i] += elerror[i]*elerror[i];
//        }
//    }
//    
//    //out << sqrt(globerrors[1])<< setw(19);
//    //    out << "\n";
//    //    //out << "Errors associated with L2 or H1 space\n";
//    //    out << "\nH1 Norm = "    << sqrt(globerrors[0]) << endl;
//    //    out << "\nL2 Norm = "    << sqrt(globerrors[1]) << endl;
//    //    out << "\nSemi H1 Norm = "    << sqrt(globerrors[2]) << endl;
//    //    out << "\n=============================\n"<<endl;
//    
//    Error.Resize(3, 0.);
//    Error[1] = sqrt(globerrors[1]);
//}
//
//
//void AjustarContorno(TPZGeoMesh *gmesh)
//{
//    gmesh->BuildConnectivity();
//    int nel = gmesh->NElements();
//    for(int ie = 0; ie<nel; ie++){
//        TPZGeoEl *gel = gmesh->ElementVec()[ie];
//        if(!gel || gel->HasSubElement()) continue;
//        
//        int matid = gel->MaterialId();
//        if(matid>0) continue;
//        
//        int nside = gel->NSides();
//        TPZGeoElSide thisside(gel,nside-1);
//        TPZGeoElSide neighbour = thisside.Neighbour();
//        
//        int dimel = gel->Dimension();
//        int dimn = neighbour.Element()->Dimension();
//        if(dimn==dimel) continue;
//        
//        int nsubel = neighbour.Element()->NSubElements();
//        
//        for(int is=0; is < nsubel; is++){
//            TPZGeoEl *neigel = neighbour.Element()->SubElement(is);
//            if(neigel->HasSubElement()) continue;
//            
//            int mylevel  = gel->Level();
//            int neiglevel = neigel->Level();
//            
//            if(mylevel!=neiglevel){
//                TPZVec< TPZGeoEl * > filhos;
//                gel->Divide(filhos);
//                break;
//            }
//        }
//    }
//}
//
//void ChangeInternalConnectOrder(TPZCompMesh *mesh){
//    
//    int nEl= mesh-> NElements();
//    int dim = mesh->Dimension();
//    
//    for (int iel=0; iel<nEl; iel++) {
//        TPZCompEl *cel = mesh->ElementVec()[iel];
//        if (!cel) continue;
//        int ncon = cel->NConnects();
//        int corder = 0;
//        int nshape = 0;
//        int nshape2 = 0;
//        
//        if(cel->Dimension()== dim)
//        {
//            TPZConnect &conel = cel->Connect(ncon-1);
//            corder = conel.Order();
//            nshape = conel.NShape();
//            
//            int neworder = corder + 1;
//            int64_t cindex = cel->ConnectIndex(ncon-1);
//            conel.SetOrder(neworder,cindex);
//            
//            if(fTriang){
//                nshape2 = (corder + 2)*(corder + 2)-1;
//            }else{//Quadrilateral
//                nshape2 = 2*(corder + 1)*(corder + 2);
//            }
//            
//            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
//            intel->SetPreferredOrder(neworder);
//            nshape = intel->NConnectShapeF(ncon-1,neworder);
//            
//            if(dim==2){
//                if(nshape2!=nshape){
//                    DebugStop();
//                }
//            }
//            
//            conel.SetNShape(nshape);
//            mesh->Block().Set(conel.SequenceNumber(),nshape);
//        }
//    }
//    mesh->CleanUpUnconnectedNodes();
//    mesh->ExpandSolution();
//}
