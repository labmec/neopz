#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Common.h"
#include "pzgeoelbc.h"
#include "pzgengrid.h"
#include "pzcheckgeom.h"
#include "pzbndcond.h"
#include "TPZBuildSBFem.h"

#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

TPZAutoPointer<TPZGeoMesh> CreateGMesh(int nelx);

TPZCompMesh *BuildSBFem(TPZAutoPointer<TPZGeoMesh> gmesh, int nx, int porder);

void IntegrateDirect(TPZCompMesh *cmesh);

#ifdef _AUTODIFF
TElasticity2DAnalytic ElastExactLower;
TElasticity2DAnalytic ElastExactUpper;
#endif

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
#ifdef _AUTODIFF
    ElastExact.fProblemType = TElasticity2DAnalytic::ESquareRoot;
    ElastExact.fE = 10;
    ElastExact.fPoisson = 0.3;
    ElastExact.fPlaneStress = 0;
    ElastExactLower = ElastExact;
    ElastExactUpper = ElastExact;
    ElastExactLower.fProblemType = TElasticity2DAnalytic::ESquareRootLower;
    ElastExactUpper.fProblemType = TElasticity2DAnalytic::ESquareRootUpper;
#endif

    int maxnelxcount = 4;
    int maxporder = 2;
    int counter = 1;
    int nx = 2;
    TPZAutoPointer<TPZGeoMesh> gmesh = CreateGMesh(nx);
    if(0)
    {
        std::cout << "Plotting the geometric mesh\n";
        //                std::ofstream outg("GMesh3D.txt");
        //                gmesh->Print(outg);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }

    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
            for(int nelxcount = 0; nelxcount < maxnelxcount; nelxcount += 1)
            {
                TPZAutoPointer<TPZGeoMesh> locgmesh = new TPZGeoMesh(gmesh);
                {
                    TPZCheckGeom check(locgmesh.operator->());
                    check.UniformRefine(nelxcount);
                }
                if(1)
                {
                    std::cout << "Plotting the geometric mesh\n";
                    //                std::ofstream outg("GMesh3D.txt");
                    //                gmesh->Print(outg);
                    std::stringstream sout;
                    sout << "SBFem_Fem_Geometry." << counter << ".vtk";
                    std::ofstream out(sout.str());
                    TPZVTKGeoMesh vtk;
                    vtk.PrintGMeshVTK(locgmesh, out,true);
                }
                TPZCompMesh *SBFem = BuildSBFem(locgmesh, nx, POrder);
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    SBFem->Reference()->Print(sout);
                    SBFem->Print(sout);
                    LOGPZ_DEBUG(logger, sout.str())
                    std::ofstream out("CompMesh.vtk");
                    TPZVTKGeoMesh vtk;
                    vtk.PrintGMeshVTK(locgmesh, out,true);
                }
#endif
                
                std::cout << "nelx = " << nelxcount << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
                Analysis->SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                SolveSist(Analysis, SBFem);
                
                
                
                
                std::cout << "Post processing\n";
                //        ElasticAnalysis->Solution().Print("Solution");
                //        mphysics->Solution().Print("expandec");
#ifdef _AUTODIFF
                Analysis->SetExact(Elasticity_exact);

#endif
                //                ElasticAnalysis->SetExact(Singular_exact);
                
                TPZManVector<STATE> errors(3,0.);
                
                long neq = SBFem->Solution().Rows();
                
                if(0)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    vecnames.Push("Displacement");
                    scalnames.Push("SigmaX");
                    scalnames.Push("SigmaY");
                    scalnames.Push("TauXY");
                    scalnames.Push("EpsX");
                    scalnames.Push("EpsY");
                    scalnames.Push("EpsXY");
                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../EmbeddedSBFemElasticity2DSolution.vtk");
                    Analysis->PostProcess(3);
                }

                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                std::cout << "Plotting shape functions\n";
                if( nelxcount == 0)
                {
                    int numshape = 25;
                    if (numshape > SBFem->NEquations()) {
                        numshape = SBFem->NEquations();
                    }
                    TPZVec<long> eqindex(numshape);
                    for (int i=0; i<numshape; i++) {
                        eqindex[i] = i;
                    }
                    Analysis->ShowShape("SBFemSingular.vtk", eqindex);
                }

                
                std::cout << "Compute errors\n";
                
                Analysis->PostProcessError(errors);
                
//                VerifyShapeFunctionIntegrity(Analysis->Mesh());
                
//                IntegrateDirect(Analysis->Mesh());
                
                std::stringstream sout;
                sout << "../EmbeddedSBFem";
                sout << "Elastic2D.txt";
                
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "(* nx " << nx  << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
                TPZFMatrix<double> errmat(1,6);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                errmat(0,3) = 1./(10*pow(2.,nelxcount));
                errmat(0,4) = neq;
                errmat(0,5) = POrder;
                std::stringstream varname;
                varname << "Errmat[[" << nelxcount+1 << "]][[" << POrder << "]] = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                
                
                delete Analysis;
                delete SBFem;
                //                exit(-1);
            }
            //            exit(-1);

    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}






void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = gMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

void IntegrateDirect(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (elgr) {
            TPZStack<TPZCompEl *,5> elstack = elgr->GetElGroup();
            int nvol = elstack.size();
            TPZElementMatrix ekvol, efvol, ekgrp, efgrp;
            elgr->CalcStiff(ekgrp, efgrp);
            for (int iv=0; iv<nvol; iv++) {
                TPZCompEl *vcel = elstack[iv];
                TPZSBFemVolume *elvol = dynamic_cast<TPZSBFemVolume *>(vcel);
                TPZElementMatrix ek,ef;
                elvol->CalcStiff(ek, ef);
                if (iv==0) {
                    ekvol = ek;
                    efvol = ef;
                }
                else
                {
                    ekvol.fMat += ek.fMat;
                    efvol.fMat += ef.fMat;
                }
            }
            ekgrp.fMat.Print("EKGRP = ",std::cout,EMathematicaInput);
            ekvol.fMat.Print("EKVOL = ",std::cout,EMathematicaInput);
            break;
        }
    }

    
}

TPZAutoPointer<TPZGeoMesh> CreateGMesh(int nelx)
{
    // only even number of elements are allowed
    if(nelx%2 != 0)
    {
        DebugStop();
    }
    TPZManVector<REAL,4> x0(3,-1.),x1(3,1.);
    x0[0] = -1;
    x0[1] = -1;
    x1[0] = 1;
    x1[1] = 1;
    x0[2] = 0.;
    x1[2] = 0.;
    TPZManVector<int,4> nx(2,nelx);
    TPZGenGrid gengrid(nx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    
    //        OneQuad(gmesh);
    gengrid.Read(gmesh,Emat1);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc1);
    gengrid.SetBC(gmesh, 6, Ebc1);
    TPZManVector<REAL,3> start(3,0.), end(3,0.);
    start[0] = -1.;
    start[1] = 1.;
    end[0] = -1.;
    gengrid.SetBC(gmesh, start, end, Ebc2);
    start = end;
    end[1] = -1.;
    gengrid.SetBC(gmesh, start, end, Ebc3);
    long nnodes = gmesh->NNodes();
    gmesh->NodeVec().Resize(nnodes+nelx/2);
    REAL delx = 2./nelx;
    for (int in=0; in < nelx/2; in++) {
        TPZManVector<REAL,3> xco(3,0.);
        xco[0] = -1.+in*delx;
        gmesh->NodeVec()[nnodes+in].Initialize(xco, gmesh);
    }
    int minel = nelx*nelx/2;
    int maxel = nelx*(nelx+1)/2;
    for (long el = minel; el < maxel; el++) {
        long firstnode = el-nelx*nelx/2+nnodes;
        TPZGeoEl *gel = gmesh->Element(el);
        gel->SetNodeIndex(0, firstnode);
        if(firstnode+1 < nnodes+nelx/2)
        {
            gel->SetNodeIndex(1, firstnode+1);
        }
        if(el == nelx*nelx/2)
        {
            TPZGeoElSide gelside(gel,0);
            TPZGeoElSide neighbour(gelside.Neighbour());
            while (neighbour != gelside) {
                if (neighbour.Element()->MaterialId() == Ebc2) {
                    int sidenodelocindex = neighbour.SideNodeLocIndex(0);
                    neighbour.Element()->SetNodeIndex(sidenodelocindex, firstnode);
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    long index = (nelx-1)*(nelx/2)-1;
    TPZGeoEl *gel1 = gmesh->Element(index);
    TPZGeoElBC(gel1,7,ESkeleton);
    TPZGeoElBC(gel1,4,ESkeleton);
    TPZGeoEl *gel2 = gmesh->Element(index+1);
    TPZGeoElBC(gel2,4,ESkeleton);
    TPZGeoElBC(gel2,5,ESkeleton);
    TPZGeoEl *gel3 = gmesh->Element(index+nelx);
    TPZGeoElBC(gel3,6,ESkeleton);
    TPZGeoElBC(gel3,7,ESkeleton);
    TPZGeoEl *gel4 = gmesh->Element(index+nelx+1);
    TPZGeoElBC(gel4,5,ESkeleton);
    TPZGeoElBC(gel4,6,ESkeleton);
    gel1->RemoveConnectivities();
    delete gel1;
    gel2->RemoveConnectivities();
    delete gel2;
    gel3->RemoveConnectivities();
    delete gel3;
    gel4->RemoveConnectivities();
    delete gel4;
    return gmesh;
}

TPZCompMesh *BuildSBFem(TPZAutoPointer<TPZGeoMesh> gmesh, int nx, int porder)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    InsertMaterialObjects(cmesh, false, true);
    {
        TPZMaterial *mat = cmesh->FindMaterial(Ebc1);
        TPZBndCond *bndcond = dynamic_cast<TPZBndCond *>(mat);
        bndcond->SetType(0);
    }
    {
        TPZMaterial *mat = cmesh->FindMaterial(Ebc2);
        TPZBndCond *bndcond = dynamic_cast<TPZBndCond *>(mat);
        bndcond->SetType(0);
    }
    {
        TPZMaterial *mat = cmesh->FindMaterial(Ebc3);
        TPZBndCond *bndcond = dynamic_cast<TPZBndCond *>(mat);
        bndcond->SetType(0);
    }
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    std::map<int,int> matidtranslation;
    matidtranslation[ESkeleton] = Emat1;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
    TPZManVector<long,10> scalingcenters(1);
    scalingcenters[0] = ((nx+1)*(nx+1)-1)/2;
    long nel = gmesh->NElements();
    TPZManVector<long,10> elementgroup(nel,-1);
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel && gel->MaterialId() == ESkeleton) {
            elementgroup[el] = 0;
        }
    }
    build.SetPartitions(elementgroup, scalingcenters);
    build.BuildComputationalMeshFromSkeleton(*cmesh);
    return cmesh;
}

