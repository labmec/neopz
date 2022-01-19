/**
 * @file
 * @brief Contains the implementation of the TPZEllipse3D methods. 
 */

#include "tpzellipse3d.h"

#include "TPZGeoLinear.h"
#include "pzshapelinear.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pznoderep.h"
#include "pzgnode.h"
#include "pzreal.h"

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"

using namespace std;
using namespace pzgeom;
using namespace pztopology;
using namespace pzshape;

/// Tolerance value (Is zero)
const double tolerance = 1.E-9;

void TPZEllipse3D::SetAxes(TPZVec<REAL> Origin, TPZVec<REAL> SemiAxeX, TPZVec<REAL> SemiAxeY, TPZGeoMesh &gmesh)
{
#ifdef PZDEBUG
	if(SemiAxeX.size() != 3 || SemiAxeY.size() != 3 || Origin.size() != 3)
	{
		cout << "\nThe two semi-axes and origin of TPZEllipse3D must be in 3D!!!\n";
		DebugStop();
	}
	
	double dotXY = (SemiAxeX[0]*SemiAxeY[0] + SemiAxeX[1]*SemiAxeY[1] + SemiAxeX[2]*SemiAxeY[2]);
	if(fabs(dotXY) > tolerance)
	{
		cout << "\nThe two semi-axes of TPZEllipse3D must be orthogonal!!!\n";
		DebugStop();
	}
	
	double dotXX = sqrt(SemiAxeX[0]*SemiAxeX[0] + SemiAxeX[1]*SemiAxeX[1] + SemiAxeX[2]*SemiAxeX[2]);
	double dotYY = sqrt(SemiAxeY[0]*SemiAxeY[0] + SemiAxeY[1]*SemiAxeY[1] + SemiAxeY[2]*SemiAxeY[2]);
	if(dotXX < tolerance || dotYY < tolerance || fabs(dotXY) > tolerance)
	{
		cout << "Null semi-axe(s) on TPZEllipse3D!!! or semi-axes not orthogonal\n";
		DebugStop();
	}
#endif
	
	fOrigin = Origin;
    fsAxeX = sqrt(SemiAxeX[0]*SemiAxeX[0] + SemiAxeX[1]*SemiAxeX[1] + SemiAxeX[2]*SemiAxeX[2]);
    fsAxeY = sqrt(SemiAxeY[0]*SemiAxeY[0] + SemiAxeY[1]*SemiAxeY[1] + SemiAxeY[2]*SemiAxeY[2]);
    fAxes.Redim(3, 3);
	for(int i = 0; i < 3; i++)
	{
        fAxes(0,i) = SemiAxeX[i]/fsAxeX;
        fAxes(1,i) = SemiAxeY[i]/fsAxeY;
	}
    fAxes(2,0) = fAxes(0,1)*fAxes(1,2)-fAxes(0,2)*fAxes(1,1);
    fAxes(2,1) = -fAxes(0,0)*fAxes(1,2)+fAxes(0,2)*fAxes(1,0);
    fAxes(2,2) = fAxes(0,0)*fAxes(1,1)-fAxes(0,1)*fAxes(1,0);
	
#ifdef PZDEBUG
    TPZFNMatrix<9,REAL> ident(3,3);
    fAxes.MultAdd(fAxes, fAxes, ident,1.,0.,1);
    for (int i=0; i<3; i++) {
        ident(i,i) -= 1.;
    }
    REAL norm = Norm(ident);
    if (norm > 1.e-6) {
        fAxes.Print("Axes = ",std::cout,EMathematicaInput);
        DebugStop();
    }
#endif
    TPZManVector<REAL,3> xini(3),xfinal(3);
    gmesh.NodeVec()[fNodeIndexes[0]].GetCoordinates(xini);
    gmesh.NodeVec()[fNodeIndexes[1]].GetCoordinates(xfinal);
    fAngleIni = ComputeAngle(xini);
    fAngleFinal = ComputeAngle(xfinal);
    if (fAngleFinal < fAngleIni) {
        fAngleFinal += 2.*M_PI;
    }
    
}

/// Compute the angle of a coordinate as a function of the axes
double TPZEllipse3D::ComputeAngle(TPZVec<REAL> &co) const
{
    TPZFNMatrix<3,REAL> comat(3,1), localco(3,1);
    for(int i=0; i<3; i++) comat(i,0) = co[i]-fOrigin[i];
    fAxes.Multiply(comat, localco);
#ifdef PZDEBUG
    {
        if(fabs(localco[2]) > 1.e-6) DebugStop();
    }
#endif
    localco(0,0) /= fsAxeX;
    localco(1,0) /= fsAxeY;
#ifdef PZDEBUG
    {
        double r = localco(0,0)*localco(0,0)+localco(1,0)*localco(1,0);
        if (fabs(r-1.) > 1.e-6) {
            DebugStop();
        }
    }
#endif
    double angle = atan2(localco(1,0),localco(0,0));
    return angle;
}

template<class T>
void TPZEllipse3D::X(TPZFMatrix<REAL> &nodeCoord,TPZVec<T> &qsi,TPZVec<T> &x) const
{
#ifdef PZDEBUG
    {
        REAL angle = fAngleIni;
        TPZFNMatrix<3,REAL> xcoloc(3,1,0.),xcoglob(3,1,0.);
        xcoloc(0,0) = fsAxeX*cos(angle);
        xcoloc(1,0) = fsAxeY*sin(angle);
        int transp = 1;
        fAxes.MultAdd(xcoloc, xcoloc, xcoglob,1.,0.,transp);
        for (int i=0; i<3; i++) {
            x[i] = fOrigin[i]+xcoglob(i,0)-nodeCoord(i,0);
        }
        T norm = 0.;
        for (int i=0; i<3; i++) {
            norm += x[i]*x[i];
        }
        norm = sqrt(norm);
        if (norm > 1.e-6) {
            DebugStop();
        }
    }
    {
        REAL angle = fAngleFinal;
        TPZFNMatrix<3,REAL> xcoloc(3,1,0.),xcoglob(3,1,0.);
        xcoloc(0,0) = fsAxeX*cos(angle);
        xcoloc(1,0) = fsAxeY*sin(angle);
        int transp = 1;
        fAxes.MultAdd(xcoloc, xcoloc, xcoglob,1.,0.,transp);
        for (int i=0; i<3; i++) {
            x[i] = fOrigin[i]+xcoglob(i,0)-nodeCoord(i,1);
        }
        T norm = 0.;
        for (int i=0; i<3; i++) {
            norm += x[i]*x[i];
        }
        norm = sqrt(norm);
        if (norm > 1.e-6) {
            DebugStop();
        }
    }
#endif
    REAL delangle = fAngleFinal - fAngleIni;
    T angle = fAngleIni + (qsi[0]+1.)/2.*delangle;
    TPZFNMatrix<3,T> xcoloc(3,1,0.),xcoglob(3,1,0.);
    xcoloc(0,0) = fsAxeX*cos(angle);
    xcoloc(1,0) = fsAxeY*sin(angle);
    int transp = 1;
    T one(1.), zero(0.);
    TPZFMatrix<T> axloc(3,3);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            axloc(i,j) = fAxes.GetVal(i,j);
        }
    }
    axloc.MultAdd(xcoloc, xcoloc, xcoglob,one,zero,transp);
    for (int i=0; i<3; i++) {
        x[i] = fOrigin[i]+xcoglob(i,0);
    }
}

template<class T>
void TPZEllipse3D::GradX(TPZFMatrix<REAL> &cornerco, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
{
    REAL delangle = fAngleFinal - fAngleIni;
    T dangledqsi = 1./2.*delangle;
    T angle = fAngleIni + (par[0]+1.)/2.*delangle;
    TPZFNMatrix<3,T> dxcoloc(3,1,0.),dxcoglob(3,1,0.);
    dxcoloc(0,0) = -fsAxeX*sin(angle);
    dxcoloc(1,0) = fsAxeY*cos(angle);
    int transp = 1;
    TPZFMatrix<T> axloc(3,3);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            axloc(i,j) = fAxes.GetVal(i,j);
        }
    }
    axloc.MultAdd(dxcoloc, dxcoloc, dxcoglob,1.,0.,transp);
    for (int i=0; i<3; i++) {
        gradx(i,0) = dangledqsi*dxcoglob(i,0);
    }

}



void TPZEllipse3D::GetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix<REAL> &nodes)
{
	for(int i = 0; i < NNodes; i++)
	{
		const int64_t nodeindex = this->fNodeIndexes[i];
		TPZGeoNode * np = &(mesh.NodeVec()[nodeindex]);
		for(int j = 0; j < 3; j++)
		{
			nodes(j,i) = np->Coord(j);
		}
	}
}//void

void TPZEllipse3D::SetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix<REAL> &nodes)
{
	for(int i = 0; i < NNodes; i++)
	{
		const int64_t nodeindex = this->fNodeIndexes[i];
		TPZGeoNode * np = &(mesh.NodeVec()[nodeindex]);
		for(int j = 0; j < 3; j++)
		{
			np->SetCoord(j,nodes(j,i));
		}
	}
    TPZManVector<REAL,3> origin(fOrigin), Ax0(3), Ax1(3);
    for(int i=0; i<3; i++)
    {
        Ax0[i] = fAxes(0,i)*fsAxeX;
        Ax1[i] = fAxes(1,i)*fsAxeY;
    }
    SetAxes(origin, Ax0, Ax1, mesh);
}//void

#include <math.h>
void TPZEllipse3D::AdjustNodesCoordinates(TPZGeoMesh &mesh)
{
    TPZFNMatrix<6,REAL> nodes(3,2);
#ifdef PZDEBUG
    {
        REAL angle = fAngleIni;
        TPZFNMatrix<3,REAL> xcoloc(3,1,0.),xcoglob(3,1,0.);
        xcoloc(0,0) = fsAxeX*cos(angle);
        xcoloc(1,0) = fsAxeY*sin(angle);
        int transp = 1;
        fAxes.MultAdd(xcoloc, xcoloc, xcoglob,1.,0.,transp);
        for (int i=0; i<3; i++) {
            nodes(i,0) = fOrigin[i]+xcoglob(i,0);
        }
    }
    {
        REAL angle = fAngleFinal;
        TPZFNMatrix<3,REAL> xcoloc(3,1,0.),xcoglob(3,1,0.);
        xcoloc(0,0) = fsAxeX*cos(angle);
        xcoloc(1,0) = fsAxeY*sin(angle);
        int transp = 1;
        fAxes.MultAdd(xcoloc, xcoloc, xcoglob,1.,0.,transp);
        for (int i=0; i<3; i++) {
            nodes(i,1) = fOrigin[i]+xcoglob(i,0);
        }
    }
#endif
    SetNodesCoords(mesh, nodes);
}//void



#include "tpzgeoelmapped.h"
#include "tpzgeomid.h"


void TPZEllipse3D::ParametricDomainNodeCoord(int64_t node, TPZVec<REAL> &nodeCoord)
{
    if(node > this->NNodes)
    {
        DebugStop();
    }
    nodeCoord.Resize(Dimension, 0.);
    switch (node) {
        case (0):
        {
            nodeCoord[0] = -1.;
            break;
        }
        case (1):
        {
            nodeCoord[0] = 1.;
            break;
        }
        default:
        {
            DebugStop();
            break;
        }
    }
}

int TPZEllipse3D::ClassId() const{
    return Hash("TPZEllipse3D") ^ TPZNodeRep<2,pztopology::TPZLine>::ClassId() << 1;
}

void TPZEllipse3D::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    TPZFNMatrix<9,REAL> vander(3,3),rotation(3,3,0.), jac(3,3);
    // create a vandermonde matrix
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            vander(i,j) = pow(i+1.,j);
        }
    }
    vander.GramSchmidt(rotation, jac);
//    rotation.Print("rot = ",std::cout, EMathematicaInput);
    REAL ax1 = 3., ax2 = 1.;
    REAL angle1 = M_PI/5.5;
    REAL angle2 = M_PI*0.96;
    TPZFNMatrix<6,REAL> xloc(3,2,0.);
    xloc(0,0) = ax1*cos(angle1);
    xloc(1,0) = ax2*sin(angle1);
    xloc(0,1) = ax1*cos(angle2);
    xloc(1,1) = ax2*sin(angle2);
    TPZManVector<REAL,3> orig(3);
    for (int i=0; i<3; i++) {
        orig[i] = lowercorner[i]+3.;
        size[i] = 6;
    }
    TPZFNMatrix<6,REAL> xfinal(3,2);
    rotation.Multiply(xloc, xfinal);
    TPZManVector<REAL,3> axis1(3),axis2(3);
    for (int i=0; i<3; i++) {
        axis1[i] = rotation(i,0)*ax1;
        axis2[i] = rotation(i,1)*ax2;
        for (int j=0; j<2; j++) {
            xfinal(i,j) += orig[i];
        }
    }
    TPZManVector<int64_t,2> index(2);
    TPZManVector<REAL,3> xco(3);
    index[0] = gmesh.NodeVec().AllocateNewElement();
    for(int i=0; i<3; i++) xco[i] = xfinal(i,0);
    gmesh.NodeVec()[index[0]].Initialize(xco, gmesh);
    index[1] = gmesh.NodeVec().AllocateNewElement();
    for(int i=0; i<3; i++) xco[i] = xfinal(i,1);
    gmesh.NodeVec()[index[1]].Initialize(xco, gmesh);

    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> *gel = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D>(index,matid,gmesh);
    gel->Geom().SetAxes(orig, axis1, axis2, gmesh);
}


template class TPZRestoreClass< TPZGeoElRefPattern<TPZEllipse3D>>;
template class TPZGeoElRefLess<TPZEllipse3D>;
/*@orlandini : I REALLY dont know why is this here, so I have commented the following lines.
If it breaks something, I am sorry.*/
//template class pzgeom::TPZNodeRep<2,TPZEllipse3D>;
