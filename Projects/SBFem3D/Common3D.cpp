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
#include "tpzquadraticquad.h"
#include "tpzgeoblend.h"
#include "pzgeoelbc.h"

#ifdef _AUTODIFF
TLaplaceExample1 ExactLaplace;

TElasticity3DAnalytic ExactElast;
#endif

int gnumthreads = 0;

#ifdef USING_BOOST
#include "boost/crc.hpp"

TPZVec<boost::crc_32_type::value_type> matglobcrc, eigveccrc, stiffcrc, matEcrc, matEInvcrc;



static void printvec(const std::string &name, TPZVec<boost::crc_32_type::value_type> &vec)
{
    std::ofstream out(name);
    int64_t nel = vec.size();
    for (int64_t el=0; el<nel; el++) {
        if(vec[el] != 0)
        {
            out << el << " " << vec[el] << std::endl;
        }
    }
}

#endif

void SolveSist(TPZAnalysis *an, TPZCompMesh *Cmesh, int numthreads)
{
    gnumthreads = numthreads;

    int64_t nel = Cmesh->NElements();
#ifdef USING_BOOST
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
#endif
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
#ifdef USING_MKL
//    TPZSkylineStructMatrix strmat(Cmesh);
    TPZSymetricSpStructMatrix strmat(Cmesh);
#else
    TPZSkylineStructMatrix strmat(Cmesh);
#endif
    strmat.SetNumThreads(gnumthreads);
    an->SetStructuralMatrix(strmat);
    
    int64_t neq = Cmesh->NEquations();
    
    if(1 || neq > 20000)
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
#ifdef USING_BOOST
        printvec(matglob.str(), matglobcrc);
        printvec(eigvec.str(), eigveccrc);
        printvec(stiff.str(), stiffcrc);
        printvec(matE.str(), matEcrc);
        printvec(matEInv.str(), matEInvcrc);
#endif
        exit(-1);
    }

#ifdef USING_BOOST
    printvec(matglob.str(), matglobcrc);
    printvec(eigvec.str(), eigveccrc);
    printvec(stiff.str(), stiffcrc);
    printvec(matE.str(), matEcrc);
    printvec(matEInv.str(), matEInvcrc);

    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    std::cout << "rhs norm " << Norm(an->Rhs()) << std::endl;
    
    if(1 || neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;

    if(0)
    {
        std::ofstream out(sol.str());
        an->Solution().Print("sol",out);
    }
#endif

}


void InsertMaterialObjects3D(TPZCompMesh *cmesh, bool scalarproblem)
{
    
    // Getting mesh dimension
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
#ifdef _AUTODIFF
        matloc->SetMaterialDataHook(ExactElast.fE, ExactElast.fPoisson);
        matloc->SetForcingFunction(ExactElast.ForcingFunction());
#endif
        cmesh->InsertMaterialObject(matloc);
    }
    else
    {
        TPZMatLaplacian *matloc = new TPZMatLaplacian(matId1);
#ifdef _AUTODIFF
        matloc->SetForcingFunction(ExactLaplace.ForcingFunction());
#endif

        matloc->SetDimension(3);
        matloc->SetSymmetric();
        material = matloc;
        nstate = 1;
        cmesh->InsertMaterialObject(matloc);
        TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
        {
            val1(0,0) = 0.01;
            TPZBndCond *BCond1 = material->CreateBC(material,Ebcpoint1,2, val1, val2);
#ifdef _AUTODIFF
            BCond1->TPZMaterial::SetForcingFunction(ExactLaplace.TensorFunction());
#endif
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
#ifdef _AUTODIFF
        BCond1->SetForcingFunction(ExactLaplace.Exact());
#endif
    }
    else
    {
        BCond1 = material->CreateBC(material,Ebc1,1, val1, val2);
#ifdef _AUTODIFF
        BCond1->SetForcingFunction(ExactElast.TensorFunction());
#endif
    }
    
    if (elasticity) {
        {
            val1.Zero();
            val2.Zero();
            val1(0,0) = 0.01;
            val1(2,2) = 0.01;
            val1(1,1) = 0.01;
            TPZMaterial *BCond1 = material->CreateBC(material,Ebcpoint1 ,2, val1, val2);
#ifdef _AUTODIFF
            BCond1->SetForcingFunction(ExactElast.TensorFunction());
#endif
            cmesh->InsertMaterialObject(BCond1);
        }
        {
            val1.Zero();
            val2.Zero();
            val1(1,1) = 0.01;
            val1(2,2) = 0.01;
            TPZMaterial *BCond1 = material->CreateBC(material,Ebcpoint2 ,2, val1, val2);
#ifdef _AUTODIFF
            BCond1->SetForcingFunction(ExactElast.TensorFunction());
#endif
            cmesh->InsertMaterialObject(BCond1);
        }
        {
            val1.Zero();
            val2.Zero();
            val1(2,2) = 0.01;
            TPZMaterial *BCond1 = material->CreateBC(material,Ebcpoint3 ,2, val1, val2);
#ifdef _AUTODIFF
            BCond1->SetForcingFunction(ExactElast.TensorFunction());
#endif
            cmesh->InsertMaterialObject(BCond1);
        }
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
    
    TPZAcademicGeoMesh acadgmesh(nelx,TPZAcademicGeoMesh::EHexa);
    TPZManVector<int,6> bcids(6,-1);
    acadgmesh.SetBCIDVector(bcids);
    acadgmesh.SetMaterialId(EGroup);
    
    TPZAutoPointer<TPZGeoMesh> gmesh = acadgmesh.CreateGeoMesh();
    {
        int64_t nnode = gmesh->NNodes();
        for (int64_t n=0; n<nnode; n++) {
            TPZGeoNode *nptr = &gmesh->NodeVec()[n];
            TPZManVector<REAL,3> co(3);
            nptr->GetCoordinates(co);
            co[0] -= 1./2.;
            co[1] -= 1./2.;
            nptr->SetCoord(co);
        }
    }
    {
        TPZManVector<int64_t,2> nodes(1,0);
        int64_t index;
        gmesh->CreateGeoElement(EPoint, nodes, Ebcpoint1, index);
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[0].GetCoordinates(xco);
        std::cout << "xco " << xco << std::endl;
    }
    {
        TPZManVector<int64_t,2> nodes(1,nelx);
        int64_t index;
        gmesh->CreateGeoElement(EPoint, nodes, Ebcpoint2, index);
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[nelx].GetCoordinates(xco);
        std::cout << "xco " << xco << std::endl;
    }
    {
        TPZManVector<int64_t,2> nodes(1,(nelx+1)*(nelx+1)-1);
        int64_t index;
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
TPZGeoMesh *ReadUNSWSBGeoFile(const std::string &filename, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices)
{
    
    int maxvol = -1;
    
    std::ifstream file(filename);
    
    map<set<int64_t> , int64_t> midnode;
    string buf;
    getline(file,buf);
    if(!file) DebugStop();
    int64_t nnodes, nvolumes;
    file >> nnodes >> nvolumes;
    elpartition.Resize(nvolumes*6, -1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        for (int i=0; i<3; i++) {
            file >> xco[i];
        }
        gmesh->NodeVec()[in].Initialize(xco, *gmesh);
    }
#ifdef PZDEBUG
    std::set<int64_t> badvolumes;
#endif
    int64_t nothing;
    file >> nothing;
    for (int64_t iv=0; iv<nvolumes; iv++) {
#ifdef PZDEBUG
        map<set<int64_t>,int64_t> nodepairs;
#endif
        int nfaces;
        file >> nfaces;
        for (int face = 0; face < nfaces; face++) {
            int elnnodes;
            file >> elnnodes;
            
            TPZManVector<int64_t,10> nodes(elnnodes);
            for (int i=0; i<elnnodes; i++) {
                file >> nodes[i];
                nodes[i]--;
#ifdef PZDEBUG
                if (i>0) {
                    set<int64_t> edge;
                    edge.insert(nodes[i-1]);
                    edge.insert(nodes[i]);
                    nodepairs[edge]++;
                }
                if (i==elnnodes-1) {
                    set<int64_t> edge;
                    edge.insert(nodes[0]);
                    edge.insert(nodes[i]);
                    nodepairs[edge]++;
                }
#endif
            }
            
            // tototototo
            if (maxvol != -1 && iv >= maxvol) {
                continue;
            }
            if (elnnodes == 1)
            {
                int64_t index;
                MElementType eltype = EPoint;
                gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                elpartition[index] = iv;
                
            }
            else if (elnnodes == 2)
            {
                int64_t index;
                MElementType eltype = EOned;
                gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                elpartition[index] = iv;
                
            }
            else if (elnnodes == 3 || elnnodes == 4)
            {
                int64_t index;
                MElementType eltype = ETriangle;
                if (elnnodes == 4) {
                    eltype = EQuadrilateral;
                }
                gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                elpartition[index] = iv;
            }
            else if(elnnodes == 8)
            {
                int64_t index;
                new TPZGeoElRefPattern<pzgeom::TPZQuadraticQuad> (nodes, ESkeleton, *gmesh,  index);
                elpartition[index] = iv;

            }
            else if(elnnodes > 4)
            {
                set<int64_t>  elnodes;
                TPZManVector<REAL,3> midxco(3,0.);
                for (int i=0; i<elnnodes; i++) {
                    elnodes.insert(nodes[i]);
                    TPZManVector<REAL,3> x(3);
                    gmesh->NodeVec()[nodes[i]].GetCoordinates(x);
                    //                    std::cout << "x " << x << endl;
                    for(int j=0; j<3; j++) midxco[j] += x[j]/elnnodes;
                }
                int64_t midindex = -1;
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
                    TPZManVector<int64_t,3> nodeindices(3);
                    for (int in=0; in<2; in++) {
                        nodeindices[in] = nodes[(triangle+in)%elnnodes];
                    }
                    nodeindices[2] = midindex;
                    int64_t index;
                    gmesh->CreateGeoElement(ETriangle, nodeindices, ESkeleton, index);
                    elpartition[index] = iv;
                }
            }
            else
            {
                DebugStop();
            }
        }
#ifdef PZDEBUG
        bool suspicious = false;
        for (auto it = nodepairs.begin(); it != nodepairs.end(); it++) {
            if(it->second != 2) suspicious = true;
        }
        if (suspicious == true) {
            std::cout << "volume " << iv << " has no closure\n";
            badvolumes.insert(iv);
        }
#endif
        if (elpartition.size() < gmesh->NElements()+100) {
            elpartition.Resize(elpartition.size()*2, -1);
        }
    }
    // totototototo
    if (maxvol != -1) {
        nvolumes = maxvol;
    }
    int64_t nmidnodes = midnode.size();
    gmesh->NodeVec().Resize(nvolumes+nmidnodes+nnodes);
    scalingcenterindices.Resize(nvolumes, -1);
    for (int64_t in=0; in<nvolumes; in++) {
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
#ifdef PZDEBUG
    if (badvolumes.size()) {
        int64_t nel = gmesh->NElements();
        TPZManVector<REAL> elval(nel,0);
        for (int64_t el=0; el<nel; el++) {
            if (badvolumes.find(elpartition[el]) != badvolumes.end()) {
                elval[el] = 10.;
            }
        }
        {
            ofstream badel("gmesh_bad.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, badel, elval);
        }
    }
#endif
    elpartition.Resize(gmesh->NElements(), -1);
    std::cout << "Building element connectivity\n";
    gmesh->BuildConnectivity();
    return gmesh;
}



