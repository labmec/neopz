//
//  main_mixed.cpp
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
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
#include "tpzchangeel.h"
#include "tpzarc3d.h"
#include "tpzquadraticquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZGenGrid2D.h"

#include "pzlog.h"

#include "pzl2projection.h"
#include "pzmultiphysicselement.h"
#include "pzintel.h"
#include "TPZVTKGeoMesh.h"

#include "TPZMatDualHybridPoisson.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "tpzgeoelmapped.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZIntQuadQuarterPoint.h"
#include "tpzintpoints.h"
#include "pzquad.h"
#include "tpzhierarquicalgrid.h"
#include "TPZVecL2.h"

#include "Tools.h"

#include <iostream>
#include <math.h>
using namespace std;

//int const matId =1;
//int const bc1=-1;
//int const bc2=-2;
//int const bc3=-3;
//int const bc4=-4;
//int const bc5=-5;
//
//int const bc7=-7;
//int const bc8=-8;
//
//int const dirichlet =0;
//int const neumann = 1;

//TPZGeoMesh *GMesh();
//
////with hdiv
//TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
//TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
//TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson* &mymaterial, bool QuarterPointRule);
//TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder);
//TPZCompMesh *MalhaCompH1QP(TPZGeoMesh * gmesh,int ordem);
//TPZCompMesh *CMeshFluxL2(TPZGeoMesh *gmesh, int pOrder, int nodeAtOriginId);
//void ErroL2NoElemento(TPZCompMesh *hdivmesh, std::ostream &out,  int nodeAtOriginId);
//
//void TransferMatrixFromMeshes(TPZCompMesh *cmesh, TPZCompMesh *MFMesh, TPZAutoPointer< TPZMatrix<STATE> > matF,TPZAutoPointer< TPZMatrix<STATE> > matMP, int nodeAtOriginId);
//
////with hybrid method
//TPZGeoMesh * MalhaGeo(const int h);
//TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder,bool ismultiplierH1);
//void GroupElements(TPZCompMesh *cmesh);
//void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//
//
//void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide);
//void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
//void Prefinamento(TPZCompMesh * cmesh, int ndiv,int porder);
//
//void QuarterPointRef(TPZGeoMesh *gmesh, int nodeAtOriginId);
//
//void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
//void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//
//void SolExata3D(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
//void f_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &f);
//void Dirichlet_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannEsquerda_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannDireita_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannAcima_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//
//
//
//void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numThreads=0);
//void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
//void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out);
//void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out);
//void ErroHDivNoElemento(TPZCompMesh *hdivmesh, std::ostream &out,  int nodeAtOriginId);
//
//void ChangeSideConnectOrderConnects(TPZCompMesh *mesh, int reduction_value);
//void ChangeInternalConnectOrder(TPZCompMesh *mesh);
//
//void PrefinamentoRibsHybridMesh(TPZCompMesh *cmesh);
//
//void IntegrationRuleConvergence(bool intQuarterPoint);
//void NormMax(TPZFMatrix<STATE> A, REAL &val, int &indexi, int &indexj);
//void ChangeIntegrationRule(TPZCompMesh *cmesh, int porder,bool IsQPrule);
//void GlobalSubMatrix(TPZCompMesh *cmesh, TPZAutoPointer< TPZMatrix<STATE> > mat, int nodeAtOriginId, bool matInicial, std::ofstream &subMat);
//
//void DirectUniRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide);
//
//void ComputeCharacteristicHElSize(TPZGeoMesh * geometry, REAL & h_min, REAL & rho_min);
//TPZGeoMesh *GMeshTeste();
//
//void ParametricCircle(REAL &radius,REAL &theta, TPZManVector<REAL,3> & xcoor);
//TPZGeoMesh *MeshCircle( int ndiv);
//TPZGeoMesh *CurvedMesh();
//
//TPZGeoMesh *CurvedMesh2( int ndiv);
//TPZGeoMesh *VerticalExtrusion(TPZGeoMesh *malha2D );
//void ParametricfunctionZ(const TPZVec<REAL> &par, TPZVec<REAL> &X);

//#define SmoothSol

#ifdef PZ_LOG
static PZLogger logger("PiraHP.main");
#endif


#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#define print

bool runhdiv = true;//hibrido ou hdiv

//long ComputeAverageBandWidth(TPZCompMesh *cmesh){
//
//    TPZVec<int64_t> skyline;
//    cmesh->Skyline(skyline);
//    int nelemSkyline = 0;
//    long nel = skyline.NElements();
//    for(long i=0; i<nel; i++){
//        nelemSkyline += i-skyline[i]+1;
//        if(nelemSkyline < 0) std::cout << "negativo\n";
//    }
//    long averageband = nelemSkyline/skyline.NElements();
//    return averageband;
//}


int main(int argc, char *argv[])
{

    
    bool QuarterPoint = false; //geometria qp
    bool QuarterPointRule = false; //para regra de integracao qp
    
    bool CondenseQuarterPointRule = false;//para condensar os elementos usando a regra de integracao QP
    
    bool RefAfim=true;//refinamento afim
    
    bool HpRefine=false;//para p refinamento
    
    bool HDivMaisMais = true;//para RT+
    
    bool Is2DQ= false;
    int order_reduce = 0;
    
    int p = 1;
    int pq = p;
    int pp = p;
    // int order=0;
    if(HDivMaisMais){
        //pq +=1;
        pp +=1;
        order_reduce = 1;
    }
    
    ///Refinamento
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeRefPatterns(2);
    
    
    
    
    
    //string outputfile("Solution_mphysics");
    std::stringstream name;
    
    
    
    //std::ofstream myerrorfile("TesteT_4Uni4Direc.txt");
    std::ofstream myerrorfile("Malha3D.txt");
    
    for(int ndiv=2; ndiv<3; ndiv++)
    {
        int pq = p;
        int pp = p;
        
        
        TPZGeoMesh * gmesh = NULL;
        
        if (Is2DQ) {
            gmesh= CurvedMesh();//GMesh();
            
            int nodeAtOriginId=1;
            QuarterPointRef(gmesh, nodeAtOriginId);

           DirectUniRef(gmesh, nodeAtOriginId,ndiv);
            
            
        }else{
            TPZGeoMesh * geo_2D = CurvedMesh();
            gmesh = VerticalExtrusion(geo_2D);
        }
        UniformRefine(gmesh, ndiv);
        
        #ifdef print
                    {
                        std::ofstream malhaOut("Geometry.vtk");
                        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
                    }
        #endif
        
        
       
        //
        
        
        //        for(long el=0; el < gmesh->NElements(); el++)
        //        {
        //            TPZGeoEl *gel = gmesh->Element(el);
        //            gel->SetFather(-1);
        //
        //        }
        //
        //        QuarterPointRef(gmesh, 1);
        //
        //        DirectUniRef(gmesh, 1,ndiv);
        

        
        
//                {
//                    int ndirecRef=ndiv;
//                    set<int> MatId;
//                    MatId.insert(-6);
//                    for(int ref=0;ref<ndirecRef;ref++){
//                        TPZRefPatternTools::RefineDirectional(gmesh, MatId);
//                    }
//
//                }
        
        
        
        
        
        TPZCompMesh * cmesh1= CMeshFlux(gmesh, pq);
        TPZCompMesh * cmesh2= CMeshPressure(gmesh, pp);
        
        if(HpRefine){
            
            Prefinamento(cmesh1, ndiv,pq);
            Prefinamento(cmesh2, ndiv,pp);
            
        }
        
        if(HDivMaisMais)
        {
            ChangeInternalConnectOrder(cmesh1);
            
        }
        
        
        //        {
        //            std::ofstream out("cmeshPosPRef.txt");
        //            cmesh1->Print(out);
        //        }
        
        
        
        
        
        TPZVec<TPZCompMesh *> meshvec(2);
        meshvec[0] = cmesh1;
        meshvec[1] = cmesh2;
        TPZMixedPoisson * mymaterial;
        TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,mymaterial,QuarterPointRule);
        
        
        
        {
            long nEl=mphysics->NElements();
            for(int el=0; el< nEl;el++){

                TPZCompEl *compel=mphysics->Element(el);
                TPZGeoEl *gel=compel->Reference();
                if (gel->MaterialId()==-6) {
                    int ncon=compel->NConnects();
                    for(int con=0;con<ncon;con++){
                        TPZConnect &conec=compel->Connect(con);
                        conec.SetCondensed(true);

                    }
                }

            }


        }
        
        
        
        mphysics->InitializeBlock();
        
        
        
        TPZAnalysis anMP(mphysics);
        
        ResolverSistema(anMP, mphysics,3);
        long averageband_Depois = ComputeAverageBandWidth(mphysics);
        cout<<"averageband Depois de Reenumerar = " <<averageband_Depois<<std::endl;
        
        
        //pos-process
                TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                std::stringstream name;
                name << "HpQPTriang" <<ndiv<< ".vtk";
                std::string paraviewfile(name.str());
                PosProcessMultphysics(meshvec,  mphysics, anMP, paraviewfile);
        
        //Erro global
        
        myerrorfile << "\n-------------------------------------------------- \n";
        myerrorfile << "Ndiv = "<< ndiv << " p = " << p << "\n";
        myerrorfile << "DOF Total = " << cmesh1->NEquations() + cmesh2->NEquations()<< "\n";
        myerrorfile << "DOF Condensado = " << mphysics->NEquations() << "\n";
        TPZVec<REAL> erros(3);
        
        myerrorfile<<"\n\nErro da simulacao multifisica  para o fluxo";
        ErrorHDiv(cmesh1,myerrorfile);
        
        
        myerrorfile<<"\n\nErro da simulacao multifisica  para a pressao";
        ErrorL2(cmesh2,myerrorfile);
        
        cmesh1->CleanUp();
        cmesh2->CleanUp();
        delete cmesh1;
        delete cmesh2;
        delete gmesh;
        
        
        
        
        
        
        
        
        
        
    }
    
    return 0;
    
}



