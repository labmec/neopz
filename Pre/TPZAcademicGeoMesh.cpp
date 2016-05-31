//
//  TPZAcademicGeoMesh.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/29/16.
//
//

#include "TPZAcademicGeoMesh.h"

#include "pzshapelinear.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

#include "pzcheckmesh.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzpoisson3d.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "TPZRefPatternTools.h"


#include <stdio.h>

#include "pzlog.h"

#include "pzgeoelbc.h"

/** Initialiazing file for Log4CXX for this project */
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.academicmesh"));
#endif



TPZAcademicGeoMesh::TPZAcademicGeoMesh() : fMeshType(EHexa), fBCNumbers(6,-1), fShouldDeform(false), fNumberElements(1)
{
    REAL coord[8][3] = {
        //        {0,0,0},
        //        {1,0,0},
        //        {1,1,0},
        //        {0,1,0},
        //        {0,0,1},
        //        {1,0,1},
        //        {1,1,1},
        //        {0,1,1}
        {-0.5,-0.5,-0.5},
        {2,-1,-1},
        {1.1,1.1,-0.1},
        {-1,2,-1},
        {-1,-1,2},
        {1.2,-0.2,1.2},
        {2,2,2},
        {-0.5,1.5,1.5}
    };
    fDeformed.NodeVec().Resize(8);
    TPZManVector<long,8> indices(8);
    for (int i=0; i<8; i++) {
        indices[i] = i;
        for (int c=0; c<3; c++) {
            fDeformed.NodeVec()[i].SetCoord(c, coord[i][c]);
        }
    }
    long index;
    fDeformed.CreateGeoElement(ECube, indices, 1, index);
}

/// Deform the geometric mesh according to the coordinates of fDeformed
void TPZAcademicGeoMesh::DeformGMesh(TPZGeoMesh &gmesh)
{
    long nnodes = gmesh.NodeVec().NElements();
    TPZManVector<REAL,3> xbefore(3),xafter(3);
    for (long nod=0; nod<nnodes; nod++) {
        gmesh.NodeVec()[nod].GetCoordinates(xbefore);
        for (int i=0; i<3; i++) {
            xbefore[i] = 2.*xbefore[i]-1.;
        }
        fDeformed.ElementVec()[0]->X(xbefore, xafter);
        gmesh.NodeVec()[nod].SetCoord(xafter);
    }
}

static int pyramid[2][5]=
{
    {0,1,2,3,4},
    {4,5,6,7,2}
};
static int tetraedra[2][4]=
{
    {1,2,5,4},
    {4,7,3,2}
};
static int tetraedra_2[6][4]=
{
    {1,2,5,4},
    {4,7,3,2},
    {0,1,2,4},
    {0,2,3,4},
    {4,5,6,2},
    {4,6,7,2}
};

static int prism[2][6] =
{
    {0,1,2,4,5,6},
    {0,2,3,4,6,7}
};


void TPZAcademicGeoMesh::GenerateNodes(TPZGeoMesh *gmesh)
{
    long nelem = fNumberElements;
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (long i=0; i<=nelem; i++) {
        for (long j=0; j<=nelem; j++) {
            for (long k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
}

TPZGeoMesh *TPZAcademicGeoMesh::PyramidalAndTetrahedralMesh()
{
    long nelem = fNumberElements;
    int MaterialId = fMaterialId;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh);
    
    for (long i=0; i<nelem; i++) {
        for (long j=0; j<nelem; j++) {
            for (long k=0; k<nelem; k++) {
                TPZManVector<long,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef LOG4CXX
                {
                    std::stringstream sout;
                    sout << "Pyramid and tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el=0; el<2; el++)
                {
                    TPZManVector<long,5> elnodes(5);
                    for (int il=0; il<5; il++) {
                        elnodes[il] = nodes[pyramid[el][il]];
                    }
                    long index;
                    gmesh->CreateGeoElement(EPiramide, elnodes, MaterialId, index);
                    elnodes.resize(4);
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    AddBoundaryElements(gmesh);
    if (fShouldDeform) {
        DeformGMesh(*gmesh);
    }
    return gmesh;
}

TPZGeoMesh *TPZAcademicGeoMesh::TetrahedralMesh()
{
    long nelem = fNumberElements;
    int MaterialId = fMaterialId;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh);
    
    for (long i=0; i<nelem; i++) {
        for (long j=0; j<nelem; j++) {
            for (long k=0; k<nelem; k++) {
                TPZManVector<long,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef LOG4CXX
                {
                    std::stringstream sout;
                    sout << "Tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el=0; el<6; el++)
                {
                    TPZManVector<long,4> elnodes(4);
                    long index;
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra_2[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    AddBoundaryElements(gmesh);
    if (fShouldDeform) {
        DeformGMesh(*gmesh);
    }
    return gmesh;
}


TPZGeoMesh *TPZAcademicGeoMesh::HexahedralMesh()
{
    long nelem = fNumberElements;
    int MaterialId = fMaterialId;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh);
    
    for (long i=0; i<nelem; i++) {
        for (long j=0; j<nelem; j++) {
            for (long k=0; k<nelem; k++) {
                TPZManVector<long,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Cube nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                long index;
                gmesh->CreateGeoElement(ECube, nodes, MaterialId, index);
            }
        }
    }
    gmesh->BuildConnectivity();
    AddBoundaryElements(gmesh);
    if (fShouldDeform) {
        DeformGMesh(*gmesh);
    }
    return gmesh;
}

/// verify if the faces without neighbour should be orthogonal to the main planes
void TPZAcademicGeoMesh::CheckConsistency(TPZGeoMesh *mesh)
{
    long nel = mesh->NElements();
    for(long el=0; el<nel; el++) {
        TPZGeoEl *gel = mesh->ElementVec()[el];
        int nsides = gel->NSides();
        for(int is=0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            if(gelside.Dimension() != 2) {
                continue;
            }
            if(gelside.Neighbour() != gelside) {
                continue;
            }
            TPZManVector<REAL,2> xi(2,0.);
            gelside.CenterPoint(xi);
            TPZFNMatrix<6,REAL> axes(2,3);
            TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
            REAL detjac;
            gelside.Jacobian(xi, jac, axes, detjac, jacinv);
            TPZManVector<REAL,3> x(3,0.);
            gelside.X(xi, x);
            TPZManVector<REAL,3> normal(3);
            normal[0] = fabs(axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1));
            normal[1] = fabs(-axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0));
            normal[2] = fabs(axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0));
            REAL tol = 1.e-6;
            REAL xmin = 1., xmax = 0.;
            int numtol = 0;
            for(int i=0; i<3; i++) {
                if(xmin > x[i]) xmin = x[i];
                if(xmax < x[i]) {
                    xmax = x[i];
                }
                if(normal[i] > tol) {
                    numtol++;
                }
            }
            if(numtol != 1) {
                DebugStop();
            }
            if(xmin > tol && xmax < 1.-tol) {
                DebugStop();
            }
        }
    }
}


int TPZAcademicGeoMesh::AddBoundaryElements(TPZGeoMesh *gmesh)
{
    long nel = gmesh->NElements();
    for(long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        int nsides = gel->NSides();
        for(int is=0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            if(gelside.Dimension() != 2) {
                continue;
            }
            if(gelside.Neighbour() != gelside) {
                continue;
            }
            TPZManVector<REAL,2> xi(2,0.);
            gelside.CenterPoint(xi);
            TPZFNMatrix<6,REAL> axes(2,3);
            TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
            REAL detjac;
            gelside.Jacobian(xi, jac, axes, detjac, jacinv);
            TPZManVector<REAL,3> x(3,0.);
            gelside.X(xi, x);
            TPZManVector<REAL,3> normal(3);
            normal[0] = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
            normal[1] = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
            normal[2] = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
            REAL tol = 1.e-6;
            REAL xmin = 1., xmax = 0.;
            int numfound = 0;
            int dir = -5;
            for (int i=0; i<3; i++) {
                if (fabs(normal[i] - 1.) < tol) {
                    dir = 1+i;
                    numfound++;
                }
                if (fabs(normal[i] + 1.) < tol) {
                    dir = -1-i;
                    numfound++;
                }
            }
            if (dir == -5) {
                DebugStop();
            }
            if(numfound != 1) {
                DebugStop();
            }
            switch (dir) {
                case -3:
                    TPZGeoElBC(gelside,fBCNumbers[0]);
                    break;
                case 1:
                    TPZGeoElBC(gelside,fBCNumbers[1]);
                    break;
                case 2:
                    TPZGeoElBC(gelside,fBCNumbers[2]);
                    break;
                case -1:
                    TPZGeoElBC(gelside,fBCNumbers[3]);
                    break;
                case -2:
                    TPZGeoElBC(gelside,fBCNumbers[4]);
                    break;
                case 3:
                    TPZGeoElBC(gelside,fBCNumbers[5]);
                    break;
                default:
                    DebugStop();
                    break;
            }
        }
    }
}


