#include "Common3D.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "TPZMatElasticity2D.h"
#include "pzelast3d.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"

#include "TPZAcademicGeoMesh.h"
#include "pzgengrid.h"
#include "TPZBuildSBFem.h"

#include "TPZVTKGeoMesh.h"

#include "TPZSBFemElementGroup.h"
#include "pzinterpolationspace.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "pzgeoelbc.h"

TLaplaceExampleSmooth ExactLaplace;

TElasticity3DAnalytic ExactElast;
#include "boost/crc.hpp"

TPZVec<boost::crc_32_type::value_type> matglobcrc, eigveccrc, stiffcrc, matEcrc, matEInvcrc;
int gnumthreads = 0;

static void printvec(const std::string &name, TPZVec<boost::crc_32_type::value_type> &vec)
{
    std::ofstream out(name);
    long nel = vec.size();
    for (long el=0; el<nel; el++) {
        if(vec[el] != 0)
        {
            out << el << " " << vec[el] << std::endl;
        }
    }
}

void SolveSist(TPZAnalysis *an, TPZCompMesh *Cmesh, int numthreads)
{
    gnumthreads = numthreads;

    long nel = Cmesh->NElements();
    matglobcrc.Resize(nel, 0);
    eigveccrc.Resize(nel, 0);
    stiffcrc.Resize(nel, 0);
    matEcrc.Resize(nel, 0);
    matEInvcrc.Resize(nel, 0);
    std::stringstream matglob,eigvec,stiff,sol,matE,matEInv;
    matglob << "matglob_" << gnumthreads << "_" << nel << ".txt";
    eigvec << "eigvec_" << gnumthreads << "_" << nel << ".txt";
    stiff << "stiff_" << gnumthreads << "_" << nel << ".txt";
    sol << "sol_" << gnumthreads << "_" << nel << ".txt";
    matE << "matE_" << gnumthreads << "_" << nel << ".txt";
    matEInv << "matEInv_" << gnumthreads << "_" << nel << ".txt";
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
#ifdef USING_MKL
    TPZSkylineStructMatrix strmat(Cmesh);
//    TPZSymetricSpStructMatrix strmat(Cmesh);
#else
    TPZSkylineStructMatrix strmat(Cmesh);
#endif
    strmat.SetNumThreads(gnumthreads);
    an->SetStructuralMatrix(strmat);
    
    long neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an->SetSolver(step);
    
    
    try {
        an->Assemble();
    } catch (...) {
        printvec(matglob.str(), matglobcrc);
        printvec(eigvec.str(), eigveccrc);
        printvec(stiff.str(), stiffcrc);
        printvec(matE.str(), matEcrc);
        printvec(matEInv.str(), matEInvcrc);
        exit(-1);
    }

    printvec(matglob.str(), matglobcrc);
    printvec(eigvec.str(), eigveccrc);
    printvec(stiff.str(), stiffcrc);
    printvec(matE.str(), matEcrc);
    printvec(matEInv.str(), matEInvcrc);

#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    std::cout << "rhs norm " << Norm(an->Rhs()) << std::endl;
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
#endif
    
    {
        std::ofstream out(sol.str());
        an->Solution().Print("sol",out);
    }
}

void HarmonicNeumannLeft(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = -M_PI*exp(M_PI*x[0])*sin(M_PI*x[1]);
}

void HarmonicNeumannRight(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = M_PI*exp(M_PI*x[0])*sin(M_PI*x[1]);
}

void Harmonic_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    val[0] = exp(M_PI*xv[0])*sin(M_PI*xv[1]);
    deriv(0,0) = M_PI*val[0];
    deriv(1,0) = M_PI*exp(M_PI*xv[0])*cos(M_PI*xv[1]);
    
}


void InsertMaterialObjects3D(TPZCompMesh *cmesh, bool scalarproblem)
{
    // Plane strain assumption
    int planestress = 0;
    
    // Getting mesh dimension
    int dim = 2;
    int matId1 = Emat1;
    
    TPZMaterial *material;
    int nstate = 1;
    bool elasticity = false;
    if (scalarproblem == false) {
        elasticity = true;
    }
    if (elasticity)
    {
        TPZElasticity3D *matloc = new TPZElasticity3D(matId1);
        material = matloc;
        nstate = 3;
        //        REAL lamelambda = 1.0e9,lamemu = 0.5e3, fx= 0, fy = 0;
        matloc->SetMaterialDataHook(ExactElast.fE, ExactElast.fPoisson);
        matloc->SetForcingFunction(ExactElast.ForcingFunction());
        cmesh->InsertMaterialObject(matloc);
        TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
        {
            val1(0,0) = 0.01;
            val1(1,1) = 0.01;
            val1(2,2) = 0.01;
            TPZBndCond *BCond1 = material->CreateBC(material,Ebcpoint1,2, val1, val2);
            BCond1->TPZMaterial::SetForcingFunction(ExactElast.TensorFunction());
            cmesh->InsertMaterialObject(BCond1);
        }
        {
            val1(0,0) = 0.;
            val1(1,1) = 0.01;
            val1(2,2) = 0.01;
            TPZBndCond *BCond1 = material->CreateBC(material,Ebcpoint2,2, val1, val2);
            BCond1->TPZMaterial::SetForcingFunction(ExactElast.TensorFunction());
            cmesh->InsertMaterialObject(BCond1);
        }
        {
            val1(0,0) = 0.;
            val1(1,1) = 0.;
            val1(2,2) = 0.01;
            TPZBndCond *BCond1 = material->CreateBC(material,Ebcpoint3,2, val1, val2);
            BCond1->TPZMaterial::SetForcingFunction(ExactElast.TensorFunction());
            cmesh->InsertMaterialObject(BCond1);
        }
    }
    else
    {
        TPZMatLaplacian *matloc = new TPZMatLaplacian(matId1);
        matloc->SetForcingFunction(ExactLaplace.ForcingFunction());
        matloc->SetDimension(3);
        matloc->SetSymmetric();
        material = matloc;
        nstate = 1;
        cmesh->InsertMaterialObject(matloc);
        TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
        {
            val1(0,0) = 0.01;
            TPZBndCond *BCond1 = material->CreateBC(material,Ebcpoint1,2, val1, val2);
            BCond1->TPZMaterial::SetForcingFunction(ExactLaplace.TensorFunction());
            cmesh->InsertMaterialObject(BCond1);
        }
    }
    //material->SetBiotAlpha(Alpha);cade o metodo?
    
    
    TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond1;
    if(elasticity==0)
    {
        BCond1 = material->CreateBC(material,Ebc1,0, val1, val2);
        BCond1->SetForcingFunction(ExactLaplace.Exact());
    }
    else
    {
        BCond1 = material->CreateBC(material,Ebc1,1, val1, val2);
        BCond1->SetForcingFunction(ExactElast.TensorFunction());
    }
    
    {
        val1.Zero();
        val2.Zero();
        TPZMaterial *BCond2 = material->CreateBC(material, Ebc2, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }
    {
        val1.Zero();
        val2.Zero();
        TPZMaterial *BCond2 = material->CreateBC(material, Ebc3, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }
    {
        val1.Zero();
        val2.Zero();
        TPZMaterial *BCond2 = material->CreateBC(material, Ebc4, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }
    {
        val1.Zero();
        val2.Zero();
        TPZMaterial *BCond2 = material->CreateBC(material, Ebc5, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
    }

    
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BSkeleton = material->CreateBC(material,ESkeleton,1, val1, val2);
    
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BSkeleton);
    
}

TPZCompMesh *SetupSquareMesh3D(int nelx, int nrefskeleton, int porder, bool elasticityproblem)
{
    
    TPZAcademicGeoMesh acadgmesh(nelx,TPZAcademicGeoMesh::EPyramid);
    TPZManVector<int,6> bcids(6,-1);
    acadgmesh.SetBCIDVector(bcids);
    acadgmesh.SetMaterialId(EGroup);
    
    TPZAutoPointer<TPZGeoMesh> gmesh = acadgmesh.CreateGeoMesh();
    {
        long nnode = gmesh->NNodes();
        for (long n=0; n<nnode; n++) {
            TPZGeoNode *nptr = &gmesh->NodeVec()[n];
            TPZManVector<REAL,3> co(3);
            nptr->GetCoordinates(co);
            co[0] -= 1./2.;
            co[1] -= 1./2.;
            nptr->SetCoord(co);
        }
    }
    {
        TPZManVector<long,2> nodes(1,0);
        long index;
        gmesh->CreateGeoElement(EPoint, nodes, Ebcpoint1, index);
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[0].GetCoordinates(xco);
        std::cout << "xco " << xco << std::endl;
    }
    {
        TPZManVector<long,2> nodes(1,nelx);
        long index;
        gmesh->CreateGeoElement(EPoint, nodes, Ebcpoint2, index);
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[nelx].GetCoordinates(xco);
        std::cout << "xco " << xco << std::endl;
    }
    {
        TPZManVector<long,2> nodes(1,(nelx+1)*(nelx+1)-1);
        long index;
        gmesh->CreateGeoElement(EPoint, nodes, Ebcpoint3, index);
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[(nelx+1)*(nelx+1)-1].GetCoordinates(xco);
        std::cout << "xco " << xco << std::endl;
    }
    gmesh->BuildConnectivity();
    
    std::map<int,int> matmap;
    matmap[EGroup] = 1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    build.StandardConfiguration();
    build.DivideSkeleton(nrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 0;
    if (elasticityproblem) {
        problemtype = 0;
    }
    else
    {
        problemtype = 1;
    }
    InsertMaterialObjects3D(SBFem,problemtype);
    
    
    build.BuildComputationMesh(*SBFem);
    
    if(1)
    {
        std::ofstream outg("GMesh3D.txt");
        gmesh->Print(outg);
        std::ofstream out("Geometry3D.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
}

using namespace std;
/// Read a UNSWSBFem file
TPZGeoMesh *ReadUNSWSBGeoFile(const std::string &filename, TPZVec<long> &elpartition, TPZVec<long> &scalingcenterindices)
{
    std::ifstream file(filename);

    map<set<long> , long> midnode;
    string buf;
    getline(file,buf);
    if(!file) DebugStop();
    long nnodes, nvolumes;
    file >> nnodes >> nvolumes;
    elpartition.Resize(nvolumes*6, -1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    for (long in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        for (int i=0; i<3; i++) {
            file >> xco[i];
        }
        gmesh->NodeVec()[in].Initialize(xco, *gmesh);
    }
    long nothing;
    file >> nothing;
    for (long iv=0; iv<nvolumes; iv++) {
        int nfaces;
        file >> nfaces;
        for (int face = 0; face < nfaces; face++) {
            int elnnodes;
            file >> elnnodes;
            
            TPZManVector<long,10> nodes(elnnodes);
            for (int i=0; i<elnnodes; i++) {
                file >> nodes[i];
                nodes[i]--;
            }
            if (elnnodes == 3 || elnnodes == 4)
            {
                long index;
                MElementType eltype = ETriangle;
                if (elnnodes == 4) {
                    eltype = EQuadrilateral;
                }
                gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                elpartition[index] = iv;
            }
            else
            {
                set<long>  elnodes;
                TPZManVector<REAL,3> midxco(3,0.);
                for (int i=0; i<elnnodes; i++) {
                    elnodes.insert(nodes[i]);
                    TPZManVector<REAL,3> x(3);
                    gmesh->NodeVec()[nodes[i]].GetCoordinates(x);
//                    std::cout << "x " << x << endl;
                    for(int j=0; j<3; j++) midxco[j] += x[j]/elnnodes;
                }
                long midindex = -1;
                if (midnode.find(elnodes) == midnode.end()) {
                    midindex = gmesh->NodeVec().AllocateNewElement();
                    gmesh->NodeVec()[midindex].Initialize(midxco, *gmesh);
                    midnode[elnodes] = midindex;
                }
                else
                {
                    midindex = midnode[elnodes];
                }
                for (int triangle = 0; triangle <elnnodes; triangle++) {
                    TPZManVector<long,3> nodeindices(3);
                    for (int in=0; in<2; in++) {
                        nodeindices[in] = nodes[(triangle+in)%elnnodes];
                    }
                    nodeindices[2] = midindex;
                    long index;
                    gmesh->CreateGeoElement(ETriangle, nodeindices, ESkeleton, index);
                    elpartition[index] = iv;
                }
            }
        }
        if (elpartition.size() < gmesh->NElements()+100) {
            elpartition.Resize(elpartition.size()*2, -1);
        }
    }
    long nmidnodes = midnode.size();
    gmesh->NodeVec().Resize(nvolumes+nmidnodes+nnodes);
    scalingcenterindices.Resize(nvolumes, -1);
    for (long in=0; in<nvolumes; in++) {
        TPZManVector<REAL,3> xco(3);
        for (int i=0; i<3; i++) {
            file >> xco[i];
        }
        gmesh->NodeVec()[nnodes+nmidnodes+in].Initialize(xco, *gmesh);
        scalingcenterindices[in] = nnodes+nmidnodes+in;
    }
    {
        ofstream mirror("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, mirror);
    }
    elpartition.Resize(gmesh->NElements(), -1);
    gmesh->BuildConnectivity();
    return gmesh;
}



