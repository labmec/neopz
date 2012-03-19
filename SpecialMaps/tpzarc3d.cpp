/**
 * @file
 * @brief Contains the implementation of the TPZArc3D methods. 
 */
#include "tpzarc3d.h"
#include "pzshapelinear.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pznoderep.h"
#include "pzgnode.h"
#include "pzreal.h"
#include "tpzgeoelmapped.h"

using namespace std;
using namespace pzgeom;
using namespace pztopology;
using namespace pzshape;


void TPZArc3D::ComputeAtributes(TPZFMatrix<REAL> &coord)
{
#ifdef DEBUG
	/** Cross[(mid-ini),(fin-ini)] */
	double CrossX, CrossY, CrossZ;
	try
	{
		CrossX = fabs(-(coord(1,1)*coord(2,0)) + coord(1,2)*coord(2,0) + coord(1,0)*coord(2,1) - coord(1,2)*coord(2,1) - coord(1,0)*coord(2,2) + coord(1,1)*coord(2,2));
		CrossY = fabs(coord(0,1)*coord(2,0) - coord(0,2)*coord(2,0) - coord(0,0)*coord(2,1) + coord(0,2)*coord(2,1) + coord(0,0)*coord(2,2) - coord(0,1)*coord(2,2));
		CrossZ = fabs(-(coord(0,1)*coord(1,0)) + coord(0,2)*coord(1,0) + coord(0,0)*coord(1,1) - coord(0,2)*coord(1,1) - coord(0,0)*coord(1,2) + coord(0,1)*coord(1,2));	
	}
	catch(...){
		std::cout << "Arc3D element with nodes coordinates non initialized!!!\n";
		DebugStop();
	}
	
	/** If Cross[(mid-ini),(fin-ini)] == 0, than the 3 given points are co-linear */
	if(CrossX <= 1.E-6 && CrossY <= 1.E-6 && CrossZ <= 1.E-6)
	{
		cout << "The 3 given points that define an TPZArc3D are co-linear!\n";
		cout << "Method aborted!";
		
		DebugStop();
	}
#endif
	
	/** fIBaseCn -> Basis Change Matrix: from Base(R2) to Canonic(R3) | fIBaseCn 1st column = BaseX | fIBaseCn 2nd column = BaseY */
	TPZFNMatrix<9> IBaseCnCP(3,3,0.), NotUsedHere(3,3,0.);
	fIBaseCn.Resize(3,3); fICnBase.Resize(3,3); fCenter3D.Resize(3);
	for(int i = 0; i < 3; i++)
	{
		IBaseCnCP(i,0) = coord(i,0) - coord(i,2);
		IBaseCnCP(i,1) = coord(i,1) - coord(i,2);
	}
	
	/** fIBaseCn 3rd column = BaseZ = Cross[BaseX,BaseY] |  */
	IBaseCnCP(0,2) = -coord(1,1)*coord(2,0) + coord(1,2)*coord(2,0) + coord(1,0)*coord(2,1) - coord(1,2)*coord(2,1) - coord(1,0)*coord(2,2) + coord(1,1)*coord(2,2);
	IBaseCnCP(1,2) =  coord(0,1)*coord(2,0) - coord(0,2)*coord(2,0) - coord(0,0)*coord(2,1) + coord(0,2)*coord(2,1) + coord(0,0)*coord(2,2) - coord(0,1)*coord(2,2);
	IBaseCnCP(2,2) = -coord(0,1)*coord(1,0) + coord(0,2)*coord(1,0) + coord(0,0)*coord(1,1) - coord(0,2)*coord(1,1) - coord(0,0)*coord(1,2) + coord(0,1)*coord(1,2);
	
	TPZFNMatrix<9> axest(IBaseCnCP.Rows(),IBaseCnCP.Cols());
	IBaseCnCP.GramSchmidt(axest,NotUsedHere);
	
	fIBaseCn = axest;
	axest.Transpose(&fICnBase);
	
	/** fICnBase -> Basis Change Matrix: from Canonic(R3) to Base(R2) | fICnBase(i,0) = Inverse[fIBaseCn] */
	//	fIBaseCn = fICnBase; fIBaseCn.Transpose();
	
	double Xa, Ya, Xb, Yb;
	ComputeR2Points(coord,Xa,Ya,Xb,Yb);
	
	/** Computing the Center Coordinates in R2 base */
	fXcenter = Xa/2.;
	fYcenter = (-Xa*Xb + Xb*Xb + Yb*Yb) / (2.*Yb);
	
	/** Computing the Center Coordinates in R3 base */
	TPZVec<REAL> Temp(3,0.);
	Temp[0] = fXcenter; Temp[1] = fYcenter; Temp[2] = 0.; double temp = 0.;
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++) temp += fIBaseCn(i,j)*Temp[j];
		fCenter3D[i] = temp + coord(i,2);
		temp = 0.;
	}
	
	/** Computing Radius */
	fRadius = sqrt( fXcenter*fXcenter + fYcenter*fYcenter );
	
	/** Computing (Center->First_Point) vector */
	finitialVector.Resize(2,0.);
	finitialVector[0] = Xa - fXcenter; finitialVector[1] = Ya - fYcenter;
}


/** This method compute the 3 given points with respect to R2 Basis */
void TPZArc3D::ComputeR2Points(TPZFMatrix<REAL> &coord, double &xa, double &ya, double &xb, double &yb)
{
	/** vector (ini - middle) written in new R2 base */
	TPZVec<REAL> Axe(3,0.), Temp(3,0.);
	for(int i = 0; i < 3; i++) Axe[i] = coord(i,0) - coord(i,2);
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++) Temp[i] += fICnBase(i,j)*Axe[j];
	}
	xa = Temp[0]; ya = Temp[1];
	Temp.Fill(0.);
	
	/** vector (final - middle) written in new R2 base */
	for(int i = 0; i < 3; i++) Axe[i] = coord(i,1) - coord(i,2);
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++) Temp[i] += fICnBase(i,j)*Axe[j];
	}
	xb = Temp[0]; yb = Temp[1];
	
	fAngle = ArcAngle(coord,xa, ya, xb, yb);
}

/** This method return the absolute angle with respect of the arc formed between (ini - center) and (fin - center), passing by midnode
 Note: (xm,ym) don't appear because this coordinates are always (0,0) - it's the origin of R2 basis */
double TPZArc3D::ArcAngle(TPZFMatrix<REAL> &coord, double xa, double ya, double xb, double yb) const
{
	REAL Xcenter, Ycenter, arcAngle, sinBig, cosBig;
	
	/** Computing the Center Coordinates in this R2 base M_PI */
	Xcenter = xa/2.;
	Ycenter = (-xa*xb + xb*xb + yb*yb) / (2.*yb);
	
	cosBig = (xb-Xcenter)*(xa-Xcenter)+(yb-Ycenter)*(ya-Ycenter);
	sinBig = (xa-Xcenter)*(yb-Ycenter)-(xb-Xcenter)*(ya-Ycenter);
	arcAngle = -atan2(sinBig, cosBig);
	return arcAngle;
	
	//old code (still here until above code be totally validated)
	REAL cos, Angle1, Angle2, Angle3;
	
	/** angle between (ini-center) 'n' (mid-center) */
	cos = (-((xa - Xcenter)*Xcenter) - (ya - Ycenter)*Ycenter) /
	(sqrt((xa - Xcenter)*(xa - Xcenter) + (ya - Ycenter)*(ya - Ycenter))*sqrt(Xcenter*Xcenter + Ycenter*Ycenter));
	Angle1 = acos(cos);
	
	/** angle between (fin-center) 'n' (mid-center) */
	cos = (-((xb - Xcenter)*Xcenter) - (yb - Ycenter)*Ycenter) /
	(sqrt((xb - Xcenter)*(xb - Xcenter) + (yb - Ycenter)*(yb - Ycenter))*sqrt(Xcenter*Xcenter + Ycenter*Ycenter));
	Angle2 = acos(cos);
	
	/** angle between (ini-center) 'n' (fin-center) */
	cos = ((xa - Xcenter)*(xb - Xcenter) + (ya - Ycenter)*(yb - Ycenter)) /
	(sqrt((xa - Xcenter)*(xa - Xcenter) + (ya - Ycenter)*(ya - Ycenter))*sqrt((xb - Xcenter)*(xb - Xcenter) + (yb - Ycenter)*(yb - Ycenter)));
	Angle3 = acos(cos);
	
	/** verification if midpoint is in smaller arc angle (<= pi) or in the bigger arc angle (>= pi)
	 Note: smaller and bigger arc angles reffers to the angle formed between
	 (ini-center) and (fin-center) vectors, where [smaller + bigger = 2PI] */
	if( fabs(Angle3/(Angle1 + Angle2)) >= 1. )
	{
		arcAngle = Angle3; /** Smaller Arc Angle = Angle3 */
	}
	else
	{
		arcAngle = (2.*M_PI - Angle3); /** Bigger Arc Angle = 2Pi - Angle3 */
	}
	
	return arcAngle;
}

/** Mapping -> result = f(NodesCoord,qsi) */
void TPZArc3D::X(TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
{
	/** Computing initialVector = (iniR2 - CenterR2) */
	TPZVec<REAL> MappedBASE2D(3,0.);
	
	TPZFMatrix<REAL> RotMatrix(2,2);
	double deflection = fAngle * fRadius * (loc[0] + 1.) / (2.*fRadius);
	RotMatrix(0,0) =  cos(deflection); RotMatrix(0,1) = sin(deflection);
	RotMatrix(1,0) = -sin(deflection); RotMatrix(1,1) = cos(deflection);
	
	/** MappedPoint_R2 = centerCoord + vectorRotated , where Vx = RotationMatrix . Va */
	double centerCoord, vectRotated = 0.;
	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < 2; j++) vectRotated += RotMatrix(i,j)*finitialVector[j];
		centerCoord = (1-i)*fXcenter + i*fYcenter;
		MappedBASE2D[i] = centerCoord + vectRotated;
		vectRotated = 0.;
	}
	
	/** Changing Basis of Obtained MappedPoint from R2 to R3 */
	MappedBASE2D[2] = 0.;
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			vectRotated += fIBaseCn.GetVal(i,j)*MappedBASE2D[j];
		}
		result[i] = vectRotated + nodes(i,2);
		if( fabs(result[i]) < 0.001) result[i] = 0.;
		vectRotated = 0.;
	}
}


void TPZArc3D::Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv) const
{
	jacobian.Resize(1,1); axes.Resize(1,3); jacinv.Resize(1,1);
	jacobian(0,0) = fAngle * fRadius/2.; jacinv(0,0) = 1./jacobian(0,0); detjac = jacobian(0,0);
	
	/** Computing Axes */
	TPZManVector< REAL > Vpc(3), Vpa(3), Vpb(3), Vt(3), OUT(3);
	
	TPZManVector< REAL > middle(1, 0.);
	X(coord,middle,OUT);
	
	/** Vector From MappedPoint to Ini */
	Vpa[0] = coord(0,0) - OUT[0]; Vpa[1] = coord(1,0) - OUT[1]; Vpa[2] = coord(2,0) - OUT[2];
	
	/** Vector From MappedPoint to Fin */
	Vpb[0] = coord(0,1) - OUT[0]; Vpb[1] = coord(1,1) - OUT[1]; Vpb[2] = coord(2,1) - OUT[2];
	
	X(coord,par,OUT);
	
	/** Vector From MappedPoint to Center */
	Vpc[0] = fCenter3D[0] - OUT[0]; Vpc[1] = fCenter3D[1] - OUT[1]; Vpc[2] = fCenter3D[2] - OUT[2];
	
	/** Tangent Vector From Point in the Arc */
	Vt[0] =  Vpa[1]*Vpb[0]*Vpc[1] - Vpa[0]*Vpb[1]*Vpc[1] + Vpa[2]*Vpb[0]*Vpc[2] - Vpa[0]*Vpb[2]*Vpc[2];
	Vt[1] = -Vpa[1]*Vpb[0]*Vpc[0] + Vpa[0]*Vpb[1]*Vpc[0] + Vpa[2]*Vpb[1]*Vpc[2] - Vpa[1]*Vpb[2]*Vpc[2];
	Vt[2] = -Vpa[2]*Vpb[0]*Vpc[0] + Vpa[0]*Vpb[2]*Vpc[0] - Vpa[2]*Vpb[1]*Vpc[1] + Vpa[1]*Vpb[2]*Vpc[1];
	
	double Vtnorm = 0.;
	for(int i = 0; i < 3; i++)
	{
		if( fabs(Vt[i]) < 1.E-12 ) Vt[i] = 0.;
		Vtnorm += Vt[i]*Vt[i];
	}
	if(Vtnorm < 0.) DebugStop();
	if(sqrt(Vtnorm) < 1e-16) DebugStop();
	for(int j = 0; j < 3; j++) axes(0,j) = Vt[j]/sqrt(Vtnorm);
}


TPZGeoEl *TPZArc3D::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
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


/**
 * Creates a geometric element according to the type of the father element
 */

TPZGeoEl *TPZArc3D::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									 TPZVec<int>& nodeindexes,
									 int matid,
									 int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

template class
TPZRestoreClass< TPZGeoElRefPattern<TPZArc3D>, TPZGEOELEMENTARC3DID>;

