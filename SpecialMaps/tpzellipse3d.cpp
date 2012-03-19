/**
 * @file
 * @brief Contains the implementation of the TPZEllipse3D methods. 
 */
/*
 *  Created by caju on 8/3/09.
 *  Copyright 2009 LabMeC. All rights reserved.
 *
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

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"

using namespace std;
using namespace pzgeom;
using namespace pztopology;
using namespace pzshape;

/// Tolerance value (Is zero)
const double tolerance = 1.E-10;
/// Pi value in radians
const double pi = 4.*atan(1.);

void TPZEllipse3D::SetAxes(TPZVec<REAL> Origin, TPZVec<REAL> SemiAxeX, TPZVec<REAL> SemiAxeY)
{
#if Debug
	if(SemiAxeX.size != 3 || SemiAxeY.size != 3 || Origin.size != 3)
	{
		cout << "\nThe two semi-axes and origin of TPZEllipse3D must be in 3D!!!\n";
		DebugStop();
	}
	
	double dotXY = sqrt(SemiAxeX[0]*SemiAxeY[0] + SemiAxeX[1]*SemiAxeY[1] + SemiAxeX[2]*SemiAxeY[2]);
	if(dotXY < tolerance)
	{
		cout << "\nThe two semi-axes of TPZEllipse3D must be orthogonal!!!\n";
		DebugStop();
	}
	
	double dotXX = sqrt(fSemiAxeX[0]*fSemiAxeX[0] + fSemiAxeX[1]*fSemiAxeX[1] + fSemiAxeX[2]*fSemiAxeX[2]);
	double dotYY = sqrt(fSemiAxeY[0]*fSemiAxeY[0] + fSemiAxeY[1]*fSemiAxeY[1] + fSemiAxeY[2]*fSemiAxeY[2]);
	if(dotXX < tolerance || dotYY < tolerance)
	{
		cout << "n\Null semi-axe(s) on TPZEllipse3D!!!\n";
		DebugStop();
	}
#endif
	
	fOrigin = Origin;
	fSemiAxeX.Resize(3);
	fSemiAxeY.Resize(3);
	for(int i = 0; i < 3; i++)
	{
		fSemiAxeX[i] = SemiAxeX[i];
		fSemiAxeY[i] = SemiAxeY[i];
	}
	
	fsAxeX = sqrt(fSemiAxeX[0]*fSemiAxeX[0] + fSemiAxeX[1]*fSemiAxeX[1] + fSemiAxeX[2]*fSemiAxeX[2]);
	fsAxeY = sqrt(fSemiAxeY[0]*fSemiAxeY[0] + fSemiAxeY[1]*fSemiAxeY[1] + fSemiAxeY[2]*fSemiAxeY[2]);
}

void TPZEllipse3D::X(TPZFMatrix<REAL> &nodeCoord,TPZVec<REAL> &qsi,TPZVec<REAL> &x) const
{
	TPZFNMatrix<6> ICnEllip(2,3,0.);
	TPZFNMatrix<6> IEllipCn(3,2,0.);
	TPZFNMatrix<3> viniCn(3,1,0.), vfinCn(3,1,0.);
	
	for(int i = 0; i < 3; i++)
	{
		//basis change matrix from R3 canonic base to R2 ellipse base
		ICnEllip(0,i) = fSemiAxeX[i]/fsAxeX;
		ICnEllip(1,i) = fSemiAxeY[i]/fsAxeY;
		
		//basis change matrix from R2 ellipse base to R3 canonic base
		IEllipCn(i,0) = fSemiAxeX[i]/fsAxeX;
		IEllipCn(i,1) = fSemiAxeY[i]/fsAxeY;
		
		viniCn(i,0) = nodeCoord(i,0) - fOrigin[i];
		vfinCn(i,0) = nodeCoord(i,1) - fOrigin[i];
	}
	
	// Basis change from canonic R3 base to R2 ellipse base
	//     ICnEllip.vCn = vEllip
	TPZFNMatrix<2> viniEllip(2,1,0.), vfinEllip(2,1,0.);
	ICnEllip.Multiply(viniCn,viniEllip);
	ICnEllip.Multiply(vfinCn,vfinEllip);
	
	//important verifications that avoid program crashes
	if(fabs(viniEllip(0,0)/fsAxeX) > 1.0)
	{
		cout << "\nInitial vector of TPZEllipse3D is out of range [-semiAxeX,+semiAxeX]!!!\n";
		cout << "See " << __PRETTY_FUNCTION__ << endl;
		DebugStop();
	}
	if(fabs(vfinEllip(0,0)/fsAxeX) > 1.0)
	{
		cout << "\nFinal vector of TPZEllipse3D is out of range [-semiAxeX,+semiAxeX]!!!\n";
		cout << "See " << __PRETTY_FUNCTION__ << endl;
		DebugStop();
	}
	
	//Func_varx = f(varx) = y
	double varx, Func_varx, vary;
	
	varx = viniEllip(0,0);
	Func_varx = fsAxeY * sqrt( 1.0 - (varx*varx) / (fsAxeX*fsAxeX) );
	vary = viniEllip(1,0);
	
	if(fabs(fabs(Func_varx) - fabs(vary)) > tolerance)
	{
		cout << "\nInitial node doesn't belong to an ellipse defined by the given axis on TPZEllipse3D!!!\n";
		DebugStop();
	}
	
	varx = vfinEllip(0,0);
	Func_varx = fsAxeY * sqrt( 1.0 - (varx*varx) / (fsAxeX*fsAxeX) );
	vary = vfinEllip(1,0);
	if(fabs(fabs(Func_varx) - fabs(vary)) > tolerance)
	{
		cout << "\nFinal node doesn't belong to an ellipse defined by the given axis on TPZEllipse3D!!!\n";
		cout << "See " << __PRETTY_FUNCTION__ << endl;
		DebugStop();
	}
	//end of verification
	
	//working in R2 ellipse base
	double QSI = qsi[0];
	double ang = Angle(QSI, viniEllip, vfinEllip);
	TPZFNMatrix<2> xEllip = EllipseR2equation(ang);//x in R2 ellipse base
	
	//working in R3 canonic base
	TPZFNMatrix<3> xCn(3,1,0.);//x in R3 canonic base
	//** Basis change from R2 ellipse base to canonic R3 base
	//     IEllipCn.vEllip = vCn
	IEllipCn.Multiply(xEllip, xCn);
	
	//The happy ending :)
	x.Resize(3);
	x[0] = xCn(0,0) + fOrigin[0];
	x[1] = xCn(1,0) + fOrigin[1];
	x[2] = xCn(2,0) + fOrigin[2];
}

void TPZEllipse3D::Jacobian(TPZFMatrix<REAL> &nodeCoord, TPZVec<REAL> &qsi, TPZFMatrix<REAL> &jac, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv) const
{
	TPZFNMatrix<6> ICnEllip(2,3,0.);
	TPZFNMatrix<6> IEllipCn(3,2,0.);
	TPZFNMatrix<3> viniCn(3,1,0.), vfinCn(3,1,0.);
	
	for(int i = 0; i < 3; i++)
	{
		//basis change matrix from R3 canonic base to R2 ellipse base
		ICnEllip(0,i) = fSemiAxeX[i]/fsAxeX;
		ICnEllip(1,i) = fSemiAxeY[i]/fsAxeY;
		
		//basis change matrix from R2 ellipse base to R3 canonic base
		IEllipCn(i,0) = fSemiAxeX[i]/fsAxeX;
		IEllipCn(i,1) = fSemiAxeY[i]/fsAxeY;
		
		viniCn(i,0) = nodeCoord(i,0) - fOrigin[i];
		vfinCn(i,0) = nodeCoord(i,1) - fOrigin[i];
	}
	
	//** Basis change from canonic R3 base to R2 ellipse base
	//     ICnEllip.vCn = vEllip
	TPZFNMatrix<2> viniEllip(2,1,0.), vfinEllip(2,1,0.);
	ICnEllip.Multiply(viniCn,viniEllip);
	ICnEllip.Multiply(vfinCn,vfinEllip);
	
	//important verifications that avoid program crashes
	if(fabs(viniEllip(0,0)/fsAxeX) > 1.0)
	{
		cout << "\nInitial vector of TPZEllipse3D is out of range [-semiAxeX,+semiAxeX]!!!\n";
		cout << "See " << __PRETTY_FUNCTION__ << endl;
		DebugStop();
	}
	if(fabs(vfinEllip(0,0)/fsAxeX) > 1.0)
	{
		cout << "\nFinal vector of TPZEllipse3D is out of range [-semiAxeX,+semiAxeX]!!!\n";
		cout << "See " << __PRETTY_FUNCTION__ << endl;
		DebugStop();
	}
	
	double varx, Func_varx, vary;
	
	varx = viniEllip(0,0);
	Func_varx = fsAxeY * sqrt( 1.0 - (varx*varx) / (fsAxeX*fsAxeX) );
	vary = viniEllip(1,0);
	
	if(fabs(fabs(Func_varx) - fabs(vary)) > tolerance)
	{
		cout << "\nInitial node doesn't belong to an ellipse defined by the given axis on TPZEllipse3D!!!\n";
		cout << "See " << __PRETTY_FUNCTION__ << endl;
		DebugStop();
	}
	
	varx = vfinEllip(0,0);
	Func_varx = fsAxeY * sqrt( 1.0 - (varx*varx) / (fsAxeX*fsAxeX) );
	vary = vfinEllip(1,0);
	if(fabs(fabs(Func_varx) - fabs(vary)) > tolerance)
	{
		cout << "\nFinal node doesn't belong to an ellipse defined by the given axis on TPZEllipse3D!!!\n";
		cout << "See " << __PRETTY_FUNCTION__ << endl;
		DebugStop();
	}
	//end of verification
	
	//working in R2 ellipse base
	double QSI = qsi[0];
	double ang = Angle(QSI, viniEllip, vfinEllip);
	TPZFNMatrix<2> jacEllip = DEllipseR2equationDang(ang);//jac in R2 ellipse base
	jacEllip(0,0) = jacEllip(0,0)*DAngleDqsi(viniEllip, vfinEllip);
	jacEllip(1,0) = jacEllip(1,0)*DAngleDqsi(viniEllip, vfinEllip);
	
	//working in R3 canonic base
	TPZFNMatrix<3> jacCn(3,1,0.);//jac in R3 canonic base
	//** Basis change from R2 ellipse base to canonic R3 base
	//     IEllipCn.vEllip = vCn
	IEllipCn.Multiply(jacEllip, jacCn);
	
	//The happy ending :)
	jac.Resize(1,1);
	jac(0,0) = sqrt( jacCn(0,0)*jacCn(0,0) + jacCn(1,0)*jacCn(1,0) + jacCn(2,0)*jacCn(2,0) );
	
	axes.Resize(1,3);
	for(int i = 0; i < 3; i++)
	{
		axes(0,i) = jacCn(i,0)/jac(0,0);
	}
	
	detjac = jac(0,0);
	
	jacinv.Resize(1,1);
	jacinv = 1./detjac;
}

void TPZEllipse3D::GetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix<REAL> &nodes)
{
	for(int i = 0; i < NNodes; i++)
	{
		const int nodeindex = this->fNodeIndexes[i];
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
		const int nodeindex = this->fNodeIndexes[i];
		TPZGeoNode * np = &(mesh.NodeVec()[nodeindex]);
		for(int j = 0; j < 3; j++)
		{
			np->SetCoord(j,nodes(j,i));
		}
	}
}//void

#include <math.h>
void TPZEllipse3D::AdjustNodesCoordinates(TPZGeoMesh &mesh)
{
	const int nnodes = NNodes;   //NNodes = 2
	TPZFNMatrix<3*nnodes> nodeCoord(3,nnodes);
	this->GetNodesCoords(mesh,nodeCoord);
	
	TPZFNMatrix<6> ICnEllip(2,3,0.);
	TPZFNMatrix<6> IEllipCn(3,2,0.);
	TPZFNMatrix<3> viniCn(3,1,0.), vfinCn(3,1,0.);
	
	for(int i = 0; i < 3; i++)
	{
		//basis change matrix from R3 canonic base to R2 ellipse base
		ICnEllip(0,i) = fSemiAxeX[i]/fsAxeX;
		ICnEllip(1,i) = fSemiAxeY[i]/fsAxeY;
		
		//basis change matrix from R2 ellipse base to R3 canonic base
		IEllipCn(i,0) = fSemiAxeX[i]/fsAxeX;
		IEllipCn(i,1) = fSemiAxeY[i]/fsAxeY;
		
		viniCn(i,0) = nodeCoord(i,0) - fOrigin[i];
		vfinCn(i,0) = nodeCoord(i,1) - fOrigin[i];
	}
	
	//** Basis change from canonic R3 base to R2 ellipse base, i.e.: ICnEllip.vCn = vEllip
	TPZFNMatrix<2> viniEllip(2,1,0.), vfinEllip(2,1,0.);
	ICnEllip.Multiply(viniCn,viniEllip);
	ICnEllip.Multiply(vfinCn,vfinEllip);
	
	double varx, vary, Y, val;
	int signal;
	
	//... initial node
	varx = ((fsAxeX*fsAxeY*viniEllip(0,0))/sqrt(fsAxeY*fsAxeY*viniEllip(0,0)*viniEllip(0,0) + fsAxeX*fsAxeX*viniEllip(1,0)*viniEllip(1,0)));
	//
	val = 1.0 - (varx*varx)/(fsAxeX*fsAxeX);
	if(val < 1.E-30) val = 0.;
	Y = viniEllip(1,0);
	(fabs(Y) < 1.E-30) ? signal = 0 : signal = Y/fabs(Y);
	//
	vary = signal*(fsAxeY * sqrt(val));
	viniEllip(0,0) = varx;
	viniEllip(1,0) = vary;
	//... end of initial node
	
	//... final node
	varx = ((fsAxeX*fsAxeY*vfinEllip(0,0))/sqrt(fsAxeY*fsAxeY*vfinEllip(0,0)*vfinEllip(0,0) + fsAxeX*fsAxeX*vfinEllip(1,0)*vfinEllip(1,0)));	
	//
	val = 1.0 - (varx*varx)/(fsAxeX*fsAxeX);
	if(val < 1.E-30) val = 0.;
	Y = vfinEllip(1,0);
	(fabs(Y) < 1.E-30) ? signal = 0 : signal = Y/fabs(Y);
	//
	vary = signal*(fsAxeY * sqrt(val));
	vfinEllip(0,0) = varx;
	vfinEllip(1,0) = vary;
	//... end of final node
	
	//** Basis change from R2 ellipse base to canonic R3 base, i.e.: IEllipCn.vEllip = vCn
	IEllipCn.Multiply(viniEllip, viniCn);
    IEllipCn.Multiply(vfinEllip, vfinCn);
	
	for(int coord = 0; coord < 3; coord++)
	{
		nodeCoord(coord,0) = viniCn(coord,0) + fOrigin[coord];
		nodeCoord(coord,1) = vfinCn(coord,0) + fOrigin[coord];
	}
	
	this->SetNodesCoords(mesh,nodeCoord);
}//void

TPZGeoEl *TPZEllipse3D::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
{
	if(side==2)
	{
		TPZManVector<int> nodes(3);
		nodes[0] = orig->SideNodeIndex(side,0); nodes[1] = orig->SideNodeIndex(side,1); nodes[2] = orig->SideNodeIndex(side,2); int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::ContainedSideLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::ContainedSideLocId(side,1)));
		TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	
	else if(side==0 || side==1)
	{
		TPZManVector<int> nodeindexes(1);
		nodeindexes[0] = orig->SideNodeIndex(side,0); int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	
	else PZError << "\nTPZGeoLinear::CreateBCGeoEl. Side = " << side << endl;
	return 0;
}

double TPZEllipse3D::Angle(double qsi, TPZFMatrix<REAL> &vini, TPZFMatrix<REAL> &vfin) const
{
	double angIni = atan2(vini(1,0)*fsAxeX, vini(0,0)*fsAxeY);
	double angFin = atan2(vfin(1,0)*fsAxeX, vfin(0,0)*fsAxeY);
	
	if(angFin < angIni)
	{
		angFin += 2.*pi;
	}
	
	//Linear Mapping
	return (  (1. - qsi)/2. * angIni  +  (1. + qsi)/2. * angFin  );
}

double TPZEllipse3D::DAngleDqsi(TPZFMatrix<REAL> &vini, TPZFMatrix<REAL> &vfin) const
{
	double angIni = atan2(vini(1,0)*fsAxeX, vini(0,0)*fsAxeY);
	double angFin = atan2(vfin(1,0)*fsAxeX, vfin(0,0)*fsAxeY);
	
	if(angFin < angIni)
	{
		angFin += 2.*pi;
	}
	
	//Linear Mapping
	return ( angFin/2. - angIni/2. );
}

TPZFMatrix<REAL> TPZEllipse3D::EllipseR2equation(double ang) const
{
	TPZFNMatrix<2> x(2,1,0.);
	x(0,0) = fsAxeX * cos(ang);
	x(1,0) = fsAxeY * sin(ang);
	
	return x;
}

TPZFMatrix<REAL> TPZEllipse3D::DEllipseR2equationDang(double ang) const
{
	TPZFNMatrix<2> x(2,1,0.);
	x(0,0) = - fsAxeX * sin(ang);
	x(1,0) =   fsAxeY * cos(ang);
	
	return x;
}


#include "tpzgeoelmapped.h"

/*
 * Creates a geometric element according to the type of the father element
 */

TPZGeoEl *TPZEllipse3D::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										 TPZVec<int>& nodeindexes,
										 int matid,
										 int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

/// Id for three dimensional ellipse element
#define TPZGEOELEMENTELLIPSE3DID 301
template<>
int TPZGeoElRefPattern<TPZEllipse3D>::ClassId() const
{
	return TPZGEOELEMENTELLIPSE3DID;
}

template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZEllipse3D>, TPZGEOELEMENTELLIPSE3DID>;


template class TPZGeoElRefLess<TPZEllipse3D>;
template class pzgeom::TPZNodeRep<2,TPZEllipse3D>;
