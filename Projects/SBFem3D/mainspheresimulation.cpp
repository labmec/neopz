#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Common3D.h"
#include "TPZBuildSBFem.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzaxestools.h"

#include "pzgeoelbc.h"
#include "pzbndcond.h"
#include "pzelast3d.h"
#include "TPZVTKGeoMesh.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

void AddBoundaryElementsCook(TPZGeoMesh &gmesh);

void AddBoundaryElementsSphere(TPZGeoMesh &gmesh);

void SubstituteBoundaryConditionsSphere(TPZCompMesh &cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minrefskeleton = 0;
    int maxrefskeleton = 1;
    int minporder = 1;
    int maxporder = 2;
    int counter = 1;
    int numthreads = 0;
    ExactElast.fE = 1000;
    ExactElast.fPoisson = 0.33;
    for ( int POrder = minporder; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            
            std::string filename("../spheres_10_50_sbfemesh_32_8_1.txt");
            std::string vtkfilename;

            {
                int pos = filename.find(".txt");
                std::string truncate = filename;
                truncate.erase(pos);
                std::stringstream sout;
                sout << truncate << "_t" << numthreads << "_p" << POrder << "_href" << irefskeleton << ".vtk";
                vtkfilename = sout.str();
                
            }

            TPZManVector<long,1000> elpartitions;
            TPZVec<long> scalingcenterindices;
            TPZAutoPointer<TPZGeoMesh> gmesh =ReadUNSWSBGeoFile(filename, elpartitions, scalingcenterindices);
            AddBoundaryElementsSphere(gmesh);
            elpartitions.Resize(gmesh->NElements(), -1);
            
            
            std::map<int,int> matidtranslation;
            matidtranslation[ESkeleton] = Emat1;
            TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
            build.SetPartitions(elpartitions, scalingcenterindices);
            TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
            SBFem->SetDefaultOrder(POrder);
            bool scalarproblem = false;
            InsertMaterialObjects3D(SBFem, scalarproblem);
            {
                TPZElasticity3D *mat = dynamic_cast<TPZElasticity3D *>(SBFem->FindMaterial(Emat1));
                mat->SetMaterialDataHook(1000., 0.49999);
//                mat->SetMaterialDataHook(1000., 0.33);
            }
            SubstituteBoundaryConditionsSphere(*SBFem);
            build.BuildComputationalMeshFromSkeleton(*SBFem);
            
            long nelx = SBFem->NElements();
#ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                std::stringstream sout;
                SBFem->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            if(1)
            {
                std::ofstream outg("GMesh3D.txt");
                gmesh->Print(outg);
                std::ofstream out("Geometry3D.vtk");
                TPZVTKGeoMesh vtk;
                vtk.PrintGMeshVTK(gmesh, out,true);
            }

            std::cout << "nelx = " << nelx << std::endl;
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(Analysis, SBFem, numthreads);
            
            
            //                AnalyseSolution(SBFem);
            
            std::cout << "Post processing\n";
            
            TPZManVector<STATE> errors(3,0.);
            
            long neq = SBFem->Solution().Rows();
            
            
            if(1)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                vecnames.Push("State");
                scalnames.Push("StressX");
                scalnames.Push("StressY");
                scalnames.Push("StressZ");
                Analysis->DefineGraphMesh(3, scalnames, vecnames, vtkfilename);
                Analysis->PostProcess(1);
            }
            
            if(0)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            
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

long SBFemGroup(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *grp = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(grp) return el;
    }
    return -1;
}

void AddBoundaryElementsCook(TPZGeoMesh &gmesh)
{
    std::set<long> leftset, rightset;
    long nnodes = gmesh.NNodes();
    int dim = gmesh.Dimension();
    for (long in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
        if (abs(xco[0]) < 1.e-3) {
            leftset.insert(in);
        }
        if (abs(xco[0]-48.) < 1.e-3) {
            rightset.insert(in);
        }
    }
    long nelem = gmesh.NElements();
    for (long el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh.Element(el);
        if (gel->Dimension() != dim-1) {
            DebugStop();
        }
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            int nsidenodes = gel->NSideNodes(is);
            int nfoundleft = 0;
            int nfoundright = 0;
            for (int in=0; in<nsidenodes; in++) {
                long nodeindex = gel->SideNodeIndex(is, in);
                if (leftset.find(nodeindex) != leftset.end()) {
                    nfoundleft++;
                }
                if (rightset.find(nodeindex) != rightset.end()) {
                    nfoundright++;
                }
            }
            if (nfoundright == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc1);
            }
            else if (nfoundleft == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc2);
            }
            else
            {
                TPZGeoElSide gelside(gel,is);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if (neighbour == gelside) {
                    TPZGeoElBC(gelside,Ebc3);
                }
            }
        }
    }
}

void AddBoundaryElementsSphere(TPZGeoMesh &gmesh)
{
    std::set<long> xset, yset, zset, inner, outer;
    REAL inner_radius = 10.;
    REAL outer_radius = 50.;
    long nnodes = gmesh.NNodes();
    int dim = gmesh.Dimension();
    for (long in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
        if (abs(xco[0]-100.) < 1.e-3) {
            xset.insert(in);
        }
        if (abs(xco[1]-100.) < 1.e-3) {
            yset.insert(in);
        }
        if (abs(xco[2]-100.) < 1.e-3) {
            zset.insert(in);
        }
        for (int i=0; i<3; i++) {
            xco[i] -= 100.;
        }
        REAL radius = sqrt(xco[0]*xco[0]+xco[1]*xco[1]+xco[2]*xco[2]);
        
        if (abs(radius-inner_radius) < 1.e-1) {
            inner.insert(in);
        }
        if (abs(radius-outer_radius) < 1.e-1) {
            outer.insert(in);
        }
    }
    long nelem = gmesh.NElements();
    for (long el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh.Element(el);
        if (gel->Dimension() != dim-1) {
            DebugStop();
        }
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int nsidenodes = gel->NSideNodes(is);
            int nfoundx = 0;
            int nfoundy = 0;
            int nfoundz = 0;
            int nfoundinner = 0;
            int nfoundouter = 0;
            for (int in=0; in<nsidenodes; in++) {
                long nodeindex = gel->SideNodeIndex(is, in);
                if (xset.find(nodeindex) != xset.end()) {
                    nfoundx++;
                }
                if (yset.find(nodeindex) != yset.end()) {
                    nfoundy++;
                }
                if (zset.find(nodeindex) != zset.end()) {
                    nfoundz++;
                }
                if (inner.find(nodeindex) != inner.end()) {
                    nfoundinner++;
                }
                if (outer.find(nodeindex) != outer.end()) {
                    nfoundouter++;
                }
            }
            if (nfoundx == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc1);
            }
            else if (nfoundy == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc2);
            }
            else if (nfoundz == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc3);
            }
            else if (nfoundinner == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc4);
            }
            else if (nfoundouter == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,Ebc5);
            }
            else if(gelside == neighbour)
            {
                for (int in = 0; in < nsidenodes; in++) {
                    long index = gel->SideNodeIndex(is, in);
                    TPZManVector<REAL,3> xco(3);
                    gmesh.NodeVec()[index].GetCoordinates(xco);
                    REAL radius = sqrt(xco[0]*xco[0]+xco[1]*xco[1]+xco[2]*xco[2]);

                    std::cout << "in " << in << "xco " << xco << " radius " << radius <<  std::endl;
                }
                DebugStop();
            }
        }
    }
}

void SubstituteBoundaryConditionsSphere(TPZCompMesh &cmesh)
{
//    Ebc1 -> x
//    Ebc2 -> y
//    Ebc3 -> z
//    Ebc4 -> inner
//    Ebc5 -> outer
    
    {
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(cmesh.FindMaterial(Ebc1));
        bc->SetType(2);
        TPZFNMatrix<9,REAL> val1(3,3,0.), val2(3,1,0.);
        val1(0,0) = 1.e12;
        bc->Val1().Zero();
        bc->Val1() = val1;
        bc->Val2().Zero();
        TPZAutoPointer<TPZFunction<STATE> > zero;
        bc->TPZMaterial::SetForcingFunction(zero);
    }

    {
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(cmesh.FindMaterial(Ebc2));
        bc->SetType(2);
        TPZFNMatrix<9,REAL> val1(3,3,0.), val2(3,1,0.);
        val1(1,1) = 1.e12;
        bc->Val1() = val1;
        bc->Val2().Zero();
        TPZAutoPointer<TPZFunction<STATE> > zero;
        bc->TPZMaterial::SetForcingFunction(zero);
    }
    {
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(cmesh.FindMaterial(Ebc3));
        bc->SetType(2);
        TPZFNMatrix<9,REAL> val1(3,3,0.), val2(3,1,0.);
        val1(2,2) = 1.e12;
        bc->Val1() = val1;
        bc->Val2().Zero();
        TPZAutoPointer<TPZFunction<STATE> > zero;
        bc->TPZMaterial::SetForcingFunction(zero);
    }
    {
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(cmesh.FindMaterial(Ebc4));
        bc->SetType(4);
        TPZFNMatrix<4,REAL> val1(2,2,0.), val2(2,1,0.);
        val2(0,0) = 1;
        bc->Val1().Zero();
        bc->Val2().Zero();
        TPZAutoPointer<TPZFunction<STATE> > zero;
        bc->TPZMaterial::SetForcingFunction(zero);
    }
    {
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(cmesh.FindMaterial(Ebc5));
        bc->SetType(4);
        TPZFNMatrix<4,REAL> val1(2,2,0.), val2(2,1,0.);
        val2(0,0) = 1;
        bc->Val1().Zero();
        bc->Val2().Zero();
        TPZAutoPointer<TPZFunction<STATE> > zero;
        bc->TPZMaterial::SetForcingFunction(zero);
    }

}
