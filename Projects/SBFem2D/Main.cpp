#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Common.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "pzgengrid.h"
#include "TPZMatElasticity2D.h"
#include "TPZMatLaplacian.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "TPZLagrangeMultiplier.h"
#include "pzmat1dlin.h"

#include "pzanalysis.h"

#include "pzelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "pzmultiphysicselement.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

#include "TPZSSpStructMatrix.h"

#include "TPZBuildSBFem.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#include "tpzgeoelrefpattern.h"

#include "JSON.hpp"
void rect_mesh(int numEleVer = 5, double vert_domainsize = 50);

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include <cmath>
#include <set>
//#include <Accelerate/Accelerate.h>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif


/// Define a rotation tensor which will be applied to the boundary conditions and mesh coordinates
TPZFNMatrix<9,REAL> gRotate(3,3,0.);

void OneTriangle(TPZAutoPointer<TPZGeoMesh> gmesh)
{
    REAL coord[3][3] = {
        {-1,-1,0},
        {1,-1},
        {0,0}
    };
    gmesh->NodeVec().Resize(3);
    for (int i=0; i<3; i++) {
        TPZManVector<REAL> co(3);
        for(int c=0; c<3; c++) co[c] = coord[i][c];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    TPZManVector<long> indices(3);
    for(int i=0; i<3; i++) indices[i] = (i+0)%3;
    int matid = 1;
    long index;
    TPZGeoEl *gel = gmesh->CreateGeoElement(ETriangle, indices, matid, index);
    gmesh->BuildConnectivity();
    gel->CreateBCGeoEl(3, 2);
    gel->CreateBCGeoEl(4, 3);
    gel->CreateBCGeoEl(5, 4);
}

void OneQuad(TPZAutoPointer<TPZGeoMesh> gmesh)
{
    REAL coord[4][3] = {
        {-1,-1,0},
        {1,-1},
        {0,0},
        {0,0}
    };
    gmesh->NodeVec().Resize(4);
    for (int i=0; i<4; i++) {
        TPZManVector<REAL> co(3);
        for(int c=0; c<3; c++) co[c] = coord[i][c];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    TPZManVector<long> indices(4);
    for(int i=0; i<4; i++) indices[i] = (i+0)%4;
    int matid = 1;
    long index;
    TPZGeoEl *gel = gmesh->CreateGeoElement(EQuadrilateral, indices, matid, index);
    gmesh->BuildConnectivity();
    gel->CreateBCGeoEl(4, 2);
    gel->CreateBCGeoEl(5, 3);
    gel->CreateBCGeoEl(7, 4);
}

REAL mult[] = {10./45.,9./45.,8./45.,7./45.,6./45.,5./45.};

void SingularNeumann(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    REAL Lambda0 = 2./3.;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if(theta < 0.) theta += 2.*M_PI;
    val[0] = 0;
    for (int i=0; i<6; i++) {
        REAL Lambda = Lambda0*(i+1);
        val[0] += mult[i]*Lambda*pow(r,Lambda-1.)*cos(Lambda*theta);
    }
    std::cout << " x " << x << " theta " << theta << " val " << val[0] << std::endl;
}

void Singular_exact(const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL Lambda0 = 2./3.;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if (theta<0.) {
        theta += 2.*M_PI;
    }
    
    val[0] = 0;
    deriv.Zero();
    for (int i=0; i<6; i++) {
        REAL Lambda = Lambda0*(i+1);
        val[0] += mult[i]*pow(r,Lambda)*cos(Lambda*theta);
        deriv(0,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[0]*cos(Lambda*theta)+pow(r,Lambda-2)*(Lambda)*sin(Lambda*theta)*(x[1]));
        deriv(1,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[1]*cos(Lambda*theta)-pow(r,Lambda-2)*(Lambda)*sin(Lambda*theta)*(x[0]));
    }
    
    
}

void DirichletTestProblem(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    REAL theta = atan2(x[1],x[0]);
    if(theta < 0.) theta += 2.*M_PI;
    val[0] = 0;
    if(theta < M_PI/2.)
    {
        val[0] = 1.;
    }
    else if(theta < M_PI)
    {
        val[0] = 1.-(theta-M_PI/2.)/(M_PI/2);
    }
    else if(theta < 3.*M_PI/2.)
    {
        val[0] = 0.;
    }
    else
    {
        val[0] = (theta -3.*M_PI/2.)/(M_PI/2.);
    }
    std::cout << " x " << x << " theta " << theta*180/M_PI << " val " << val[0] << std::endl;
}



void Displacement(const TPZVec<REAL> &xv, TPZVec<STATE> &result)
{
    //1    result[0] = 100.+x[0]*x[0]+3.*x[0]*x[1]+4.*x[1]*x[1];
    //1    result[1] = 200.+5.*x[0]+6.*x[1]+2.*x[0]*x[0]+4.*x[0]*x[1]+6.*x[1]*x[1];
    // apply the rotation transposed
    REAL x = gRotate(0,0)*xv[0]+gRotate(1,0)*xv[1];
    REAL y = gRotate(0,1)*xv[0]+gRotate(1,1)*xv[1];
    //2    result[0] = exp(x-y)*(1.-x)*x*(1.-y)*y;
    //2    result[1] = sin(M_PI*x)*sin(M_PI*y);
    result[0] = gRotate(0,0)*x;
    result[1] = gRotate(1,0)*x;
    
}

void Force(const TPZVec<REAL> &xv, TPZVec<STATE> &result)
{
    REAL x = xv[0];
    REAL y = xv[1];
    //1    result[0] = 14.;
    //1    result[1] = 61./2.;
    REAL G = 0.5;
    REAL Lambda = 1.;
    
    
    //2    result[0]=exp(x-y)*x*(-1.+y)*(G*(4-4.*x+5*y+3*x*y)+(3+x)*y*Lambda)+M_PI*M_PI*(G+Lambda)*cos(M_PI*x)*cos(M_PI*y);
    //2    result[1] = -exp(x-y)*(-1.+x+x*x)*(1.+(-3.+y)*y)*(G+Lambda)-M_PI*M_PI*(3*G+Lambda)*sin(M_PI*x)*sin(M_PI*y);
    result[0] = 0;
    result[1] = 0;
}



void Exact_Sol(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL x = xv[0];
    REAL y = xv[1];
    val.resize(6);
    deriv.Resize(2, 6);
    
    REAL G = 0.5;
    REAL Lambda = 1.;
    
    deriv.Zero();
    val[0] = x;
    val[1] = 0;
    
    if(y<0.5)
    {
        G = 0.5;
        Lambda = 0.;
    }
    else
    {
        G = 5.;
        Lambda = 0.;
    }
    val[2] = 2.*G+Lambda;
    val[3] = 0.;
    val[4] = 0.;
    val[5] = Lambda;
    
    /* 1111111111111 */
    //    val[0] = 100.+x[0]*x[0]+3.*x[0]*x[1]+4.*x[1]*x[1];
    //    val[1] = 200.+5.*x[0]+6.*x[1]+2.*x[0]*x[0]+4.*x[0]*x[1]+6.*x[1]*x[1];
    //    // stresses
    //    val[2] = 6.+8.*x[0]+18.*x[1];
    //    val[3] = 2.5+3.5*x[0]+6.*x[1];
    //    val[4] = val[3];
    //    val[5] = 12.+10.*x[0]+27.*x[1];
    //
    //    deriv(0,0) = 2.*x[0]+3.*x[1];
    //    deriv(1,0) = 3.*x[0]+8.*x[1];
    //
    //    deriv(0,1) = 5.+4.*x[0]+4.*x[1];
    //    deriv(1,1) = 6.+4.*x[0]+12.*x[1];
    //
    //    deriv(0,2) = 8.;
    //    deriv(1,2) = 18.;
    //    deriv(0,3) = 3.5;
    //    deriv(1,3) = 6.;
    //    deriv(0,4) = deriv(0,3);
    //    deriv(1,4) = deriv(1,3);
    //    deriv(0,5) = 10.;
    //    deriv(1,5) = 27.;
    /* 1111111111111 */
    /* 2222222222222
     val[0] = exp(x-y)*(1.-x)*x*(1.-y)*y;
     val[1] = sin(M_PI*x)*sin(M_PI*y);
     
     val[2] = (2.*G+Lambda)*exp(x-y)*(-1.+x+x*x)*(-1.+y)*y+M_PI*cos(M_PI*y)*sin(M_PI*x);
     val[3] = G*(-exp(x-y)*(-1.+x)*x*(1.+(-3.+y)*y)+M_PI*cos(M_PI*x)*sin(M_PI*y));
     val[4] = val[3];
     val[5] = exp(x-y)*(-1.+x+x*x)*(-1.+y)*y*Lambda+(2.*G+Lambda)*M_PI*cos(M_PI*y)*sin(M_PI*x);
     
     deriv(0,0) = exp(x-y)*(-1.+x+x*x)*(-1.+y)*y;
     deriv(1,0) = -exp(x-y)*(-1+x)*x*(1.-3.*y+y*y);
     
     deriv(0,1) = M_PI*cos(M_PI*x)*sin(M_PI*y);
     deriv(1,1) = M_PI*cos(M_PI*y)*sin(M_PI*x);
     
     deriv(0,2) = exp(x-y)*x*(3.+x)*(-1.+y)*y*(2.*G+Lambda)+M_PI*M_PI*Lambda*cos(M_PI*x)*cos(M_PI*y);
     deriv(0,5) = exp(x-y)*x*(3.+x)*(-1.+y)*y*Lambda + M_PI*M_PI*(2.*G+Lambda)*cos(M_PI*x)*cos(M_PI*y);
     deriv(0,3) = -exp(x-y)*G*(-1.+x+x*x)*(1.+(-3.+y)*y)-G*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
     deriv(0,4) = deriv(0,3);
     
     deriv(1,2) = -exp(x-y)*(-1.+x+x*x)*(1.+(-3.+y)*y)*(2.*G+Lambda)-M_PI*M_PI*Lambda*sin(M_PI*x)*sin(M_PI*y);
     deriv(1,5) = -exp(x-y)*(-1.+x+x*x)*(1.+(-3.+y)*y)*Lambda-M_PI*M_PI*(2.*G+Lambda)*sin(M_PI*x)*sin(M_PI*y);
     deriv(1,3) = exp(x-y)*G*(-1.+x)*x*(-4.+y)*(-1.+y)+G*M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y);
     deriv(1,4) = deriv(1,3);
     22222222222222 */
}

TPZCompMesh *SetupRegularProblem(int nelx, int nrefskeleton, int porder)
{
    TPZManVector<REAL,4> x0(3,-1.),x1(3,1.);
    x0[0] = -1;
    x0[1] = 0;
    x1[0] = 1;
    x1[1] = 2;
    x0[2] = 0.;
    x1[2] = 0.;
    TPZManVector<int,4> nx(2,nelx);
    TPZGenGrid gengrid(nx,x0,x1);
    gengrid.SetElementType(ETriangle);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    
    //        OneQuad(gmesh);
    gengrid.Read(gmesh,EGroup);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);
    
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
    InsertMaterialObjects(SBFem,problemtype);
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
    
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
}

TPZCompMesh *SetupOneArc(int numrefskeleton, int porder)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(4);
    TPZManVector<REAL,3> co(3,0.);
    gmesh->NodeVec()[0].Initialize(co, gmesh);
    co[0] = 1.;
    gmesh->NodeVec()[1].Initialize(co, gmesh);
    co.Fill(0.);
    co[1] = -1.;
    gmesh->NodeVec()[2].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(M_PI_2);
    co[1] = sin(M_PI_2);
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
    
    gmesh->Print(std::cout);
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
    InsertMaterialObjects(SBFem,problemtype);
    
    {
        TPZMaterial *BCond2 = SBFem->FindMaterial(Ebc2);
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(SingularNeumann);
        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
        BCond2->SetForcingFunction(autodummy);
        TPZBndCond *BC1 = dynamic_cast<TPZBndCond *>(SBFem->FindMaterial(Ebc1));
        BC1->Val2()(0,0) = 1;
    }
    
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


void OutputFourtyFive(TPZCompMesh *cmesh, REAL radius);

TPZCompMesh *TestHeterogeneous(int numquadrant,TPZVec<REAL> &contrast, REAL radius, int numref, int porder)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    long nodind;
    TPZManVector<REAL,3> co(3,0.);
    nodind = gmesh->NodeVec().AllocateNewElement();
    long centernode = nodind;
    gmesh->NodeVec()[nodind].Initialize(co, gmesh);
    co[0] = radius;
    nodind = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[nodind].Initialize(co, gmesh);
    TPZManVector<long> scalecenters(2*numquadrant,-1);
    
    long lastnode = nodind;
    long firstnode = nodind;
    for (int quadrant=0; quadrant<numquadrant; quadrant++) {
        TPZManVector<long,4> nodes(3,0);
        nodes[0] = lastnode;
        REAL angle = M_PI*(quadrant+1)/2;
        co[0] = radius*cos(angle);
        co[1] = radius*sin(angle);
        if (quadrant == 3)
        {
            nodind = firstnode;
        }
        else
        {
            nodind = gmesh->NodeVec().AllocateNewElement();
        }
        gmesh->NodeVec()[nodind].Initialize(co, gmesh);
        nodes[1] = nodind;
        lastnode = nodind;
        angle = M_PI*(quadrant+1)/2-M_PI/4.;
        co[0] = radius*cos(angle);
        co[1] = radius*sin(angle);
        nodind = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[nodind].Initialize(co, gmesh);
        nodes[2] = nodind;
        long elementindex;
        TPZGeoEl *arc = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (nodes, Ebc1, gmesh,elementindex);
        nodes.resize(4);
        nodes[2] = centernode;
        nodes[3] = centernode;
        int matid = EGroup;
        matid = EGroup+1;
        TPZGeoEl *gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (nodes, matid, gmesh,elementindex);
        scalecenters[elementindex] = centernode;
        
    }
    
    gmesh->BuildConnectivity();
    
    gmesh->Print(std::cout);
    if(0)
    {
        long nel = gmesh->NElements();
        for (long el=0; el<nel; el++)
        {
            TPZGeoEl *gel = gmesh->Element(el);
            TPZGeoElRefPattern < pzgeom::TPZArc3D > *arc = dynamic_cast<TPZGeoElRefPattern < pzgeom::TPZArc3D > *>(gel);
            if (!arc) {
                continue;
            }
            std::cout << "Element index " << arc->Index() << std::endl;
            TPZManVector<REAL,3> xi(1),x(3);
            for (REAL s=-1.; s<=1.; s+= 1./10.) {
                xi[0] = s;
                arc->X(xi, x);
                std::cout << "xi " << xi << " x " << x << std::endl;
            }
        }
    }
    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    matmap[EGroup+1] = Emat2;
    matmap[EGroup+2] = Emat3;
    matmap[EGroup+3] = Emat4;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    build.Configure(scalecenters);
    build.DivideSkeleton(numref);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
    InsertMaterialObjects(SBFem,problemtype);
    
    {
        TPZMaterial *BCond1 = SBFem->FindMaterial(Ebc1);
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(DirichletTestProblem);
        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
        BCond1->SetForcingFunction(autodummy);
        
        TPZMatLaplacian *mat1 = dynamic_cast<TPZMatLaplacian *> (SBFem->FindMaterial(Emat1));
        TPZMatLaplacian *mat2 = dynamic_cast<TPZMatLaplacian *> (mat1->NewMaterial());
        TPZMatLaplacian *mat3 = dynamic_cast<TPZMatLaplacian *> (mat1->NewMaterial());
        TPZMatLaplacian *mat4 = dynamic_cast<TPZMatLaplacian *> (mat1->NewMaterial());
        STATE K,F;
        mat1->GetParameters(K, F);
        mat1->SetParameters(K*contrast[0], F);
        mat2->SetParameters(K*contrast[1], F);
        mat3->SetParameters(K*contrast[2], F);
        mat4->SetParameters(K*contrast[3], F);
        mat2->SetId(Emat2);
        mat3->SetId(Emat3);
        mat4->SetId(Emat4);
        SBFem->InsertMaterialObject(mat2);
        SBFem->InsertMaterialObject(mat3);
        SBFem->InsertMaterialObject(mat4);
    }
    
    build.BuildComputationMesh(*SBFem);
    
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
                TPZManVector<REAL,3> co(3),val(1);
                ref->NodePtr(0)->GetCoordinates(co);
                DirichletTestProblem(co, val);
                long seqnum = intel->Connect(0).SequenceNumber();
                SBFem->Block().Put(seqnum, 0, 0, 0, val[0]);
                ref->NodePtr(1)->GetCoordinates(co);
                DirichletTestProblem(co, val);
                seqnum = intel->Connect(1).SequenceNumber();
                SBFem->Block().Put(seqnum, 0, 0, 0, val[0]);
            }
        }
    }
    SBFem->LoadSolution(SBFem->Solution());
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    if(numref)
    {
        OutputFourtyFive(SBFem,radius);
    }
    return SBFem;
    
}

void OutputFourtyFive(TPZCompMesh *cmesh, REAL radius)
{
    TPZGeoMesh *gmesh = cmesh->Reference();
    cmesh->LoadReferences();
    long nelg = gmesh->NElements();
    TPZSBFemVolume *elfound = NULL;
    for (long el=0; el<nelg; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->NodeIndex(0) ==9) {
            TPZGeoElSide gelside(gel,0);
            TPZGeoElSide neighbour = gelside.Neighbour();
            TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(gel->Reference());
            if (vol) {
                elfound = vol;
                break;
            }
            while (neighbour != gelside)
            {
                TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(neighbour.Element()->Reference());
                if (vol) {
                    elfound = vol;
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if(elfound) break;
        }
    }
    if(!elfound)
    {
        DebugStop();
    }
    else
    {
        TPZFMatrix<std::complex<double>> phi,coef,eigvalmatrix;
        phi = elfound->Phi();
        coef = elfound->Coeficients();
        eigvalmatrix.Resize(1, phi.Cols());
        TPZManVector<std::complex<double> > eig = elfound->Eigenvalues();
        for (long i=0; i< eig.size(); i++) {
            eigvalmatrix(0,i) = eig[i];
        }
        std::ofstream out("Diagonal.nb");
        phi.Print("phi = ",out,EMathematicaInput);
        coef.Print("coef = ",out,EMathematicaInput);
        eigvalmatrix.Print("eig = ",out,EMathematicaInput);
        out << "a = Sum[Chop[phi[[1]][[i]] coef[[i]][[1]] ksi^(-eig[[1]][[i]])], {i,1, Length[phi[[1]]]}]\n"
        << "Plot[a, {ksi, 0, 1}, PlotRange -> All]\n"
        << "Plot[a, {ksi, 0, 0.00001}, PlotRange -> All]\n";
    }
}
int main(int argc, char *argv[])
{
    
    
    //tutorial();
    rect_mesh();
    
    
    
    gRotate.Identity();
    
    
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int maxnelx = 12;
    int numrefskeleton = 1;
    int maxporder = 2;
    int counter = 1;
    for(int nelx = 10; nelx < maxnelx; nelx *=2)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            for ( int POrder = 1; POrder < maxporder; POrder += 1)
            {
                
                TPZCompMesh *SBFem = SetupRegularProblem(nelx,irefskeleton,POrder);
                //                TPZCompMesh *SBFem = SetupOneArc(irefskeleton,POrder);
//                TPZCompMesh *SBFem = ReadJSonFile("rect.json",irefskeleton,POrder);
                int numquadrant = 4;
                REAL radius = 1.;
                TPZManVector<REAL> contrast(4,1.);
                contrast[0] = 100;
                contrast[2] = 100;
                //                    contrast[0] = 10;
                //                    TPZCompMesh *SBFem = TestHeterogeneous(numquadrant , contrast, radius, irefskeleton, POrder);
                TPZSBFemElementGroup *celgrp = 0;
                long nel = SBFem->NElements();
                for (long el=0; el<nel; el++) {
                    TPZSBFemElementGroup *cel = dynamic_cast<TPZSBFemElementGroup *>(SBFem->Element(el));
                    if(cel)
                    {
                        celgrp = cel;
                        break;
                    }
                }
                
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZAnalysis * ElasticAnalysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
                ElasticAnalysis->SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                SolveSist(ElasticAnalysis, SBFem);
                
                
                
                
                std::cout << "Post processing\n";
                //        ElasticAnalysis->Solution().Print("Solution");
                //        mphysics->Solution().Print("expandec");
                
                //                ElasticAnalysis->SetExact(Harmonic_exact);
                //                ElasticAnalysis->SetExact(Singular_exact);
                
                TPZManVector<STATE> errors(3,0.);
                
                long neq = SBFem->Solution().Rows();
                
                if(1)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    //                    scalnames.Push("State");
                    vecnames.Push("State");
                    scalnames.Push("SigmaX");
                    scalnames.Push("SigmaY");
//                    ElasticAnalysis->DefineGraphMesh(2, scalnames, vecnames, "../elasticityChecker.vtk");
                    ElasticAnalysis->DefineGraphMesh(2, scalnames, vecnames, "../elasticityRegular.vtk");
                    //                    ElasticAnalysis->DefineGraphMesh(2, scalnames, vecnames, "../RegularSolution.vtk");
                    //                    ElasticAnalysis->DefineGraphMesh(2, scalnames, vecnames, "../SingularSolution.vtk");
                    //                    ElasticAnalysis->DefineGraphMesh(2, scalnames, vecnames, "../Heterogeneous.vtk");
                    ElasticAnalysis->PostProcess(3);
                }
                
                if(1)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                //                ElasticAnalysis->PostProcessError(errors);
                
                
                
                std::stringstream sout;
                //                sout << "../Heterogeneous.txt";
                sout << "../CheckerboardDiagnostic.txt";
                //                sout << "../RegularSolution.txt";
                
                std::ofstream results(sout.str(),std::ios::app);
                //                results.precision(15);
                //                results << "nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << std::endl;
                //                TPZFMatrix<double> errmat(1,3);
                //                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                //                std::stringstream varname;
                //                varname << "Errmat_" << POrder << "_" << irefskeleton << " = (1/1000000)*";
                //                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                // for circular domain with contrast
                //                results << "nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << " contrast " << contrast << std::endl;
                
                if(0)
                {
                    std::multimap<REAL,REAL> eigmap;
                    TPZManVector<double> eigval = celgrp->EigenvaluesReal();
                    TPZFMatrix<double> coef = celgrp->CoeficientsReal();
                    for (int i=0; i<eigval.size(); i++) {
                        eigmap.insert(std::pair<double,double>(eigval[i],coef(i,0)));
                    }
                    for (std::multimap<double, double>::reverse_iterator it = eigmap.rbegin(); it!=eigmap.rend(); it++) {
                        results << it->first << "|" << it->second << " ";
                    }
                }
                //                results << std::endl;
                //                results << celgrp->EigenValues() << std::endl;
                
                std::cout << "Plotting shape functions\n";
                if(0 && irefskeleton == 0)
                {
                    int numshape = 25;
                    if (numshape > SBFem->NEquations()) {
                        numshape = SBFem->NEquations();
                    }
                    TPZVec<long> eqindex(numshape);
                    for (int i=0; i<numshape; i++) {
                        eqindex[i] = i;
                    }
                    ElasticAnalysis->ShowShape("Heterogeneous.vtk", eqindex);
                }
                
                delete ElasticAnalysis;
                delete SBFem;
                //                exit(-1);
            }
            //            exit(-1);
        }
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


