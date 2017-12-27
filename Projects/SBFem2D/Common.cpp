#include "Common.h"

#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "TPZMatElasticity2D.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"

#include "pzgengrid.h"
#include "TPZBuildSBFem.h"

#include "TPZVTKGeoMesh.h"

#include "JSON.hpp"
#include "TPZSBFemElementGroup.h"
#include "pzinterpolationspace.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#ifdef _AUTODIFF
TElasticity2DAnalytic ElastExact;

TLaplaceExampleTimeDependent TimeLaplaceExact;

#endif

//TElasticity2DAnalytic::EDefState TElasticity2DAnalytic::fProblemType = TElasticity2DAnalytic::EStretchx;

void SolveSist(TPZAnalysis *an, TPZCompMesh *Cmesh)
{
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
    TPZSkylineStructMatrix strmat(Cmesh);
    //    TPZSymetricSpStructMatrix strmat(Cmesh);
    strmat.SetNumThreads(0);
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
    
    an->Assemble();
    
    //    std::ofstream andrade("../Andrade.mtx");
    //    andrade.precision(16);
    //    an->Solver().Matrix()->Print("Andrade",andrade,EMatrixMarket);
    //    std::cout << "Leaving Assemble\n";
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    
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
    
    
}

/// set the timestep of all SBFem Element groups
void SetSBFemTimestep(TPZCompMesh *CMesh, REAL delt)
{
    long nel = CMesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = CMesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!elgr) {
            continue;
        }
        if (delt > 0.) {
            elgr->SetComputeTimeDependent(delt);
        } else
        {
            elgr->SetComputeOnlyMassMatrix();
        }
    }
}
//    Compute a number of timesteps in parabolic analysis
void SolveParabolicProblem(TPZAnalysis *an, REAL delt, int nsteps, int numthreads)
{
    TPZCompMesh *Cmesh = an->Mesh();
    
    SetSBFemTimestep(Cmesh, delt);
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
#ifndef USING_MKL
    TPZSkylineStructMatrix strmat(Cmesh);
#else
    TPZSymetricSpStructMatrix strmat(Cmesh);
#endif
    strmat.SetNumThreads(numthreads);
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
    an->Assemble();
    
    //    std::ofstream andrade("../Andrade.mtx");
    //    andrade.precision(16);
    //    an->Solver().Matrix()->Print("Andrade",andrade,EMatrixMarket);
    //    std::cout << "Leaving Assemble\n";
    TPZAutoPointer<TPZMatrix<STATE> > stiff = an->Solver().Matrix();
    TPZFMatrix<STATE> rhs = an->Rhs();
    
    an->Solver().Matrix()->Zero();
    
    SetSBFemTimestep(Cmesh, 0.);
    
    an->Assemble();
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif

    TPZAutoPointer<TPZMatrix<STATE> > mass =  an->Solver().Matrix();
    
#ifdef USING_BOOST
    std::cout << "Time for assembly " << t2-t1 << std::endl;
#endif
    
    for (int istep = 0; istep < nsteps; istep++)
    {
        if(istep == 0 && neq > 20000)
        {
            std::cout << "Entering Solve neq = " << neq << "\n";
            std::cout.flush();
        }
        an->Rhs() = rhs;
        
        /**
         * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
         * @param x Is x on the above operation
         * @param y Is y on the above operation
         * @param z Is z on the above operation
         * @param alpha Is alpha on the above operation
         * @param beta Is beta on the above operation
         * @param opt Indicates if is Transpose or not
         */
//        virtual void MultAdd(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
//                             const TVar alpha=1., const TVar beta = 0., const int opt = 0) const;
        TPZFMatrix<STATE> rhstimestep;
        mass->MultAdd(an->Solution(), rhs, rhstimestep, 1./delt, 1.);
        an->Solve();

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


void InsertMaterialObjects(TPZCompMesh *cmesh, bool scalarproblem, bool applyexact)
{
    
    // Getting mesh dimension
    int dim = 2;
    int matId1 = Emat1;
    
    TPZMaterial *material;
    int nstate = 1;
    bool elasticity = false;
    if (!scalarproblem) {
        elasticity = true;
    }
    if (elasticity)
    {
        TPZMatElasticity2D *matloc = new TPZMatElasticity2D(matId1);
        material = matloc;
        nstate = 2;
        // Plane strain assumption
        //        REAL lamelambda = 1.0e9,lamemu = 0.5e3, fx= 0, fy = 0;
        REAL lamelambda = 0.,lamemu = 0.5e3, fx= 0, fy = 0;
        matloc->SetParameters(lamelambda,lamemu, fx, fy);
        TPZManVector<REAL,3> x(3,0.);
#ifdef _AUTODIFF
        // Setting up paremeters
		if (applyexact)
		{
			matloc->SetfPlaneProblem(ElastExact.fPlaneStress);
			matloc->SetElasticParameters(ElastExact.fE, ElastExact.fPoisson);
		}
#endif
        REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
        matloc->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
        
#ifdef _AUTODIFF
        if(applyexact)
        {
            matloc->SetForcingFunction(ElastExact.ForcingFunction());
        }
#endif
    }
    else
    {
        TPZMatLaplacian *matloc = new TPZMatLaplacian(matId1);
        matloc->SetDimension(2);
        matloc->SetSymmetric();
        material = matloc;
        nstate = 1;
    }
    //material->SetBiotAlpha(Alpha);cade o metodo?
    
    
    TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond1;
    if(!elasticity)
    {
        BCond1 = material->CreateBC(material,Ebc1,0, val1, val2);
    }
    else
    {
        BCond1 = material->CreateBC(material,Ebc1,1, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
            BCond1->SetForcingFunction(ElastExact.TensorFunction());
        }
#endif
    }
    
    val2(0,0) = 1.0*1000.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = 0;
    if(elasticity == 0)
    {
        BCond2 = material->CreateBC(material,Ebc2,1, val1, val2);
    }
    else
    {
        // mixed condition on the right side to desingularize the problem
        val1(0,0) = 1.;
        val1(1,1) = 1.;
        BCond2 = material->CreateBC(material,Ebc2,1, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
            BCond2->SetForcingFunction(ElastExact.TensorFunction());
        }
#endif
        val1.Zero();
    }
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond3;
    if(!elasticity)
    {
        BCond3 = material->CreateBC(material,Ebc3,0, val1, val2);
    }
    else
    {
        BCond3 = material->CreateBC(material,Ebc3,1, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
            BCond3->SetForcingFunction(ElastExact.TensorFunction());
        }
#endif
    }
    
    val2(0,0) = -1.0*1000.0;
    if(elasticity) val2(0,0) *=-1.;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond4 = material->CreateBC(material,Ebc4,1, val1, val2);
#ifdef _AUTODIFF
    if (applyexact) {
        BCond4->SetForcingFunction(ElastExact.TensorFunction());
    }

    if (elasticity && applyexact) {
        val1.Zero();
        val1(0,0) = 0.01;
        val1(1,1) = 0.01;
        TPZMaterial * BCond5 = material->CreateBC(material,EBCPoint1, 2, val1, val2);
        BCond5->SetForcingFunction(ElastExact.TensorFunction());
        val1(0,0) = 0.;
        TPZMaterial * BCond6 = material->CreateBC(material,EBCPoint2, 2, val1, val2);
        BCond6->SetForcingFunction(ElastExact.TensorFunction());
        cmesh->InsertMaterialObject(BCond5);
        cmesh->InsertMaterialObject(BCond6);
    }
#endif    
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BSkeleton = material->CreateBC(material,ESkeleton,1, val1, val2);
    
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BSkeleton);
    
}

TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder, bool scalarproblem, bool useexact)
{
    bool elasticityproblem = !scalarproblem;
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
    gengrid.Read(gmesh,EGroup);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);
    {
        TPZManVector<long,2> nodeindex(1);
        long index;
        nodeindex[0] = 0;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
        nodeindex[0] = nelx;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
        gmesh->BuildConnectivity();
    }
    
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
    InsertMaterialObjects(SBFem,!elasticityproblem, useexact);
    if(problemtype == 1)
    {
        TPZMaterial *BCond2 = SBFem->FindMaterial(Ebc2);
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(HarmonicNeumannRight);
        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
        BCond2->SetForcingFunction(autodummy);
    }
    if(problemtype == 1)
    {
        TPZMaterial *BCond4 = SBFem->FindMaterial(Ebc4);
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(HarmonicNeumannLeft);
        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
        BCond4->SetForcingFunction(autodummy);
    }
    
    
    build.BuildComputationMesh(*SBFem);
    
    if(1)
    {
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
}

TPZCompMesh *ReadJSonFile(const std::string &filename, int numrefskeleton, int pOrder)
{
    // read in json file
    std::ifstream myfile(filename);
    nlohmann::json json;
    myfile >> json;
    
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    // Part coor
    std::vector<std::vector<double>> coor = json["coor"]; // "coor" are in 2d vector
    int nnodesTotal = coor.size();
    gmesh->NodeVec().Resize(nnodesTotal);
    std::cout << "Coordinates ( " << nnodesTotal << " in total )" << std::endl;
    for (int i = 0; i < nnodesTotal; i++) {
        TPZManVector<REAL,3> co(3,0.);
        for (int j=0; j<2; j++) {
            co[j] = coor[i][j];
        }
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    
    
    
    
    
    // Part elem
    std::vector<nlohmann::json> elem = json["elem"]; // "elem" are in 1d vector with json object
    int nElem = elem.size();
    TPZVec<long> scalecenter(nElem,-1);
    for (long el=0; el<nElem; el++)
    {
        
        // sc idx
        int sc = elem[el]["sc"]; //sc idx
        scalecenter[el] = sc;
        // node idx
        std::vector<int> nodes = elem[el]["nodes"]; // nodes list
        int nnodes = nodes.size();
        TPZManVector<long> nodeindices(nnodes,-1);
        for (int k = 0; k < nnodes; k++)
        {
            nodeindices[k] = nodes[k];
        }
        // mat idx
        int mat = elem[el]["mat"];
        if (mat == Emat1) {
            mat = EGroup;
        }
        if (mat == Emat2) {
            mat = EGroup+1;
        }
        long index;
        switch(nnodes)
        {
            case 1:
                gmesh->CreateGeoElement(EPoint, nodeindices, mat, index);
                break;
            case 2:
                gmesh->CreateGeoElement(EOned, nodeindices, mat, index);
                break;
            case 3:
                gmesh->CreateGeoElement(ETriangle, nodeindices, mat, index);
                break;
            case 4:
                gmesh->CreateGeoElement(EQuadrilateral, nodeindices, mat, index);
                break;
            default:
                DebugStop();
                break;
        }
        gmesh->BuildConnectivity();
    }
    if(0)
    {
        std::map<int,TPZStack<int>> elementset;
        for (long el = 0; el<gmesh->NElements(); el++) {
            if (scalecenter[el] == -1) {
                continue;
            }
            elementset[scalecenter[el]].Push(el);
        }
        int materialindex = 100;
        for (std::map<int,TPZStack<int>>::iterator it = elementset.begin(); it != elementset.end(); it++) {
            long nel = it->second.NElements();
            long nodeindex = it->first;
            TPZManVector<REAL,3> xcenter(3);
            gmesh->NodeVec()[nodeindex].GetCoordinates(xcenter);
            for (long el=0; el<nel; el++) {
                long elindex = it->second[el];
                TPZGeoEl *gel = gmesh->Element(elindex);
                TPZManVector<REAL,3> xi(2),xco(3);
                gel->CenterPoint(gel->NSides()-1, xi);
                gel->X(xi,xco);
                long newnode = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newnode].Initialize(xco, gmesh);
                TPZManVector<long,3> cornerindexes(2);
                cornerindexes[0] = nodeindex;
                cornerindexes[1] = newnode;
                long index;
                gmesh->CreateGeoElement(EOned, cornerindexes, materialindex, index);
            }
            materialindex++;
        }
    }
    gmesh->BuildConnectivity();
    scalecenter.Resize(gmesh->NElements(), -1);
    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    matmap[EGroup+1] = Emat2;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    build.Configure(scalecenter);
    build.DivideSkeleton(numrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(pOrder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 0;
    
    bool applyexact = false;
    InsertMaterialObjects(SBFem,problemtype, applyexact);
    
    
    {
        TPZMaterial *mat = SBFem->FindMaterial(Emat1);
        TPZMaterial *mat2 = mat->NewMaterial();
        mat2->SetId(Emat2);
        TPZMatElasticity2D *matelas = dynamic_cast<TPZMatElasticity2D *>(mat2);
        matelas->SetElasticity(50, 0.);
        SBFem->InsertMaterialObject(mat2);
        TPZFNMatrix<4,STATE> val1(2,2,0.), val2(2,1,0.);
        // zero neumann at the bottom
        TPZMaterial *bnd = mat->CreateBC(mat, -1, 1, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        // zero neumann at the top
        bnd = mat->CreateBC(mat, -3, 1, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        val2(0,0) = 1.;
        bnd = mat->CreateBC(mat, -2, 1, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        val2.Zero();
        val1(1,1) = 1.;
        val1(0,0) = 1.;
        // remove rigid body modes
        bnd = mat->CreateBC(mat, -5, 2, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        val1(0,0) = 1.;
        val1(1,1) = 0.;
        bnd = mat->CreateBC(mat, -6, 2, val1, val2);
        SBFem->InsertMaterialObject(bnd);
        val1.Zero();
        val2.Zero();
        // traction to the left
        val2(0,0) = -1.;
        bnd = mat->CreateBC(mat, -4, 1, val1, val2);
        SBFem->InsertMaterialObject(bnd);
    }
    
    build.BuildComputationMesh(*SBFem);
    if(0)
    {
        long nel = SBFem->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = SBFem->Element(el);
            TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if (elgr) {
                TPZElementMatrix ek,ef;
                elgr->CalcStiff(ek, ef);
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (intel && intel->NConnects() ==3) {
                TPZGeoEl *ref = intel->Reference();
                TPZManVector<REAL,3> co(3),val(2,0.);
                ref->NodePtr(0)->GetCoordinates(co);
                val[0] = co[0]*0.01;
                long seqnum = intel->Connect(0).SequenceNumber();
                SBFem->Block().Put(seqnum, 0, 0, 0, val[0]);
                ref->NodePtr(1)->GetCoordinates(co);
                val[0] = co[0]*0.01;
                seqnum = intel->Connect(1).SequenceNumber();
                SBFem->Block().Put(seqnum, 0, 0, 0, val[0]);
            }
        }
        SBFem->LoadSolution(SBFem->Solution());
    }
    SBFem->LoadReferences();
    
    {
        std::ofstream out("JSonGeometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
        std::ofstream outg("JSonGeometry.txt");
        SBFem->Reference()->Print(outg);
        std::ofstream outc("JSonComp.txt");
        SBFem->Print(outc);
    }
    return SBFem;
}


TPZCompMesh *SetupOneArc(int numrefskeleton, int porder, REAL angle)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(4);
    TPZManVector<REAL,3> co(3,0.);
    gmesh->NodeVec()[0].Initialize(co, gmesh);
    co[0] = 1.;
    gmesh->NodeVec()[1].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(angle);
    co[1] = sin(angle);
    gmesh->NodeVec()[2].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(angle/2.);
    co[1] = sin(angle/2.);
    gmesh->NodeVec()[3].Initialize(co, gmesh);
    co.Fill(0.);
    TPZManVector<long,4> nodeindex(1,0);
    
    nodeindex[0] = 1;
    long elementid = 1;
    gmesh->CreateGeoElement(EPoint, nodeindex, Ebc1, elementid);
    
    nodeindex.Resize(3);
    // Definition of Arc coordenates
    // Create Geometrical Arc #1
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 3;
    elementid = 1;
    TPZGeoEl *arc = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (nodeindex, Ebc2, gmesh,elementid);
    
    nodeindex.Resize(4);
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 0;
    nodeindex[3] = 0;
    elementid = 2;
    TPZGeoEl *gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (nodeindex, EGroup, gmesh,elementid);
    
    gmesh->BuildConnectivity();
    
    //    gmesh->Print(std::cout);
    TPZManVector<REAL,3> xi(1),x(3);
    for (REAL s=-1.; s<=1.; s+= 1./10.) {
        xi[0] = s;
        arc->X(xi, x);
        std::cout << "xi " << xi << " x " << x << std::endl;
    }
    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    TPZManVector<long,5> elids(1,gblend->Index());
    build.AddPartition(elids, 0);
    
    build.DivideSkeleton(numrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
    bool applyexact = false;
    InsertMaterialObjects(SBFem,problemtype,applyexact);
    
    
    build.BuildComputationMesh(*SBFem);
    
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
    
}

void ElGroupEquations(TPZSBFemElementGroup *elgr, TPZVec<long> &equations)
{
    equations.Resize(0, 0);
    TPZCompMesh *cmesh = elgr->Mesh();
    int nc = elgr->NConnects();
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = elgr->Connect(ic);
        int blsize = c.NDof();
        long eqsize = equations.size();
        equations.Resize(eqsize+blsize, 0);
        long seqnum = c.SequenceNumber();
        for (int idf = 0; idf<blsize; idf++) {
            equations[eqsize+idf] = cmesh->Block().Position(seqnum)+idf;
        }
    }
}
/// Verify if the values of the shapefunctions corresponds to the value of ComputeSolution for all SBFemVolumeElements
void VerifyShapeFunctionIntegrity(TPZSBFemVolume *celv)
{
    TPZGeoEl *gel = celv->Reference();
    int dim = gel->Dimension();
    int nstate = celv->Connect(0).NState();
    TPZCompMesh *cmesh = celv->Mesh();
    int volside = gel->NSides()-1;
    TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cmesh->Element(celv->ElementGroupIndex()));
    TPZManVector<long> globeq;
    ElGroupEquations(elgr, globeq);
    TPZIntPoints *intpoints = gel->CreateSideIntegrationRule(volside, 3);
    cmesh->Solution().Zero();
    for (int ip=0; ip < intpoints->NPoints(); ip++) {
        TPZManVector<REAL,3> xi(gel->Dimension(),0.);
        TPZFNMatrix<32,REAL> phi,dphidxi;
        REAL weight;
        intpoints->Point(ip, xi, weight);
        celv->Shape(xi, phi, dphidxi);
        long neq = globeq.size();
        for (long eq=0; eq<neq; eq++) {
            long globindex = globeq[eq];
            cmesh->Solution().Zero();
            cmesh->Solution()(globindex,0) = 1.;
            cmesh->LoadSolution(cmesh->Solution());
            TPZSolVec sol;
            TPZGradSolVec dsol;
            TPZFNMatrix<9,REAL> axes(dim,3);
            celv->ComputeSolution(xi, sol, dsol, axes);
            REAL diffphi = 0., diffdphi = 0.;
            for (int istate = 0; istate < nstate; istate++) {
                diffphi += (sol[0][istate]-phi(eq*nstate+istate))*(sol[0][istate]-phi(eq*nstate+istate));
                for (int d=0; d<dim; d++) {
                    STATE diff = (dsol[0](d,istate)-dphidxi(d+istate*nstate,eq));
                    diffdphi += diff*diff;
                }
            }
            diffphi = sqrt(diffphi);
            diffdphi = sqrt(diffdphi);
            if (diffphi > 1.e-8 || diffdphi > 1.e-8) {
                std::cout << "Wrong shape function diffphi = " << diffphi << " diffdphi " << diffdphi << "\n";
            }
        }
    }
    delete intpoints;
}

/// Verify if the values of the shapefunctions corresponds to the value of ComputeSolution for all SBFemVolumeElements
void VerifyShapeFunctionIntegrity(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (elgr) {
            TPZStack<TPZCompEl *,5> elstack = elgr->GetElGroup();
            int nvol = elstack.size();
            for (int iv=0; iv<nvol; iv++) {
                TPZCompEl *vcel = elstack[iv];
                TPZSBFemVolume *elvol = dynamic_cast<TPZSBFemVolume *>(vcel);
                VerifyShapeFunctionIntegrity(elvol);
            }
        }
    }
}
