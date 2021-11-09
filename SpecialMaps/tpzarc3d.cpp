/**
 * @file
 * @brief Contains the implementation of the TPZArc3D methods. 
 */
#include "tpzarc3d.h"
#include "pzshapelinear.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzvec_extras.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pznoderep.h"
#include "pzgnode.h"
#include "pzreal.h"
#include "tpzgeoelmapped.h"
#include <math.h>

using namespace std;
using namespace pzgeom;
using namespace pztopology;
using namespace pzshape;


void TPZArc3D::ComputeAtributes(TPZFMatrix<REAL> &coord)
{
#ifdef PZDEBUG
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
    
    
    
#ifdef PZDEBUG2
    
    TPZFNMatrix<9,REAL> nodes=coord;
    TPZManVector<REAL,1> loc(1,0.0);
    TPZManVector<REAL,3> result(3.);
    
    REAL resultlocminus1, resultloc0,resultlocplus1;
    
    loc[0] = -1.;
    X(nodes,loc,result);
    for (int i=0; i<3; i++) {
        result[i] -= nodes(i,0);
    }
    resultlocminus1 = Norm(result);
    
    //cout << result << endl;
    
    loc[0] = 1.;
    X(nodes,loc,result);
    for (int i=0; i<3; i++) {
        result[i] -= nodes(i,1);
    }
    resultlocplus1 = Norm(result);
    
    //cout << result << endl;
    loc[0] = 0.;
    X(nodes,loc,result);
    for (int i=0; i<3; i++) {
        result[i] -= nodes(i,2);
    }
    resultloc0 = Norm(result);
    //cout << result << endl;

    
    // testes
    /** If result* !=0 , o ponto mapeados e diferente da coordenada recebida  */
    if(resultlocminus1 >= 1.E-6 ) // loc = -1 ponto da esquerda do arco
    {
        cout << "Method aborted!";
        
        DebugStop();
    }
    if(resultlocplus1 >= 1.E-6) // loc = 1 ponto da direita do arco
    {
        cout << "Method aborted!";
        
        DebugStop();
    }
    if(resultloc0 >= 1.E-6) // loc = 0 ponto do centro do arco
    {

        cout << "Method aborted!";
        
        DebugStop();
    }
#endif
    
}


/** This method compute the 3 given points with respect to R2 Basis */
void TPZArc3D::ComputeR2Points(TPZFMatrix<REAL> &coord, double &xa, double &ya, double &xb, double &yb)
{
	/** vector (ini - middle) written in new R2 base */
	TPZManVector<REAL,3> Axe(3,0.), Temp(3,0.);
	int i, j;
	for(i = 0; i < 3; i++) Axe[i] = coord(i,0) - coord(i,2);
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++) Temp[i] += fICnBase(i,j)*Axe[j];
	}
	xa = Temp[0]; ya = Temp[1];
	Temp.Fill(0.);
	
	/** vector (final - middle) written in new R2 base */
	for(i = 0; i < 3; i++) Axe[i] = coord(i,1) - coord(i,2);
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++) Temp[i] += fICnBase(i,j)*Axe[j];
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
    
    REAL cosMid = (0.-Xcenter)*(xa-Xcenter)+(0.-Ycenter)*(0.-Ycenter);
    REAL sinMid = (xa-Xcenter)*(0.-Ycenter)-(0.-Xcenter)*(ya-Ycenter);
    REAL arcMid = -atan2(sinMid, cosMid);
    
    if (0. <= arcMid && arcMid <= arcAngle) {
        return arcAngle;
    }
    if (arcMid <= 0. && arcAngle <= arcMid)
    {
        return arcAngle;
    }
    if (arcAngle > 0) {
        arcAngle -= 2.*M_PI;
    }
    else
    {
        arcAngle += 2.*M_PI;
    }
#ifdef PZDEBUG
    /// the arc is always a positive number ??
    if(arcAngle < 0.) DebugStop();
#endif
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



void TPZArc3D::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    REAL coords[3][3] = {
        {1,0,0},
        {0,1,0},
        {M_SQRT1_2,M_SQRT1_2,0}
    };
    size[0] = 1.;
    size[1] = 1.;
    
    TPZManVector<int64_t,3> indexes(3);
    for (int i=0; i<3; i++) {
        TPZManVector<REAL,3> cods(3,0.);
        for (int j=0; j<3; j++) {
            cods[j] = lowercorner[j]+coords[i][j];
        }
        indexes[i] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[indexes[i]].Initialize(cods, gmesh);
    }
    TPZGeoElRefPattern<TPZArc3D> *gel = new TPZGeoElRefPattern<TPZArc3D>(indexes,matid,gmesh);
}


int TPZArc3D::ClassId() const{
    return Hash("TPZArc3D") ^ pzgeom::TPZNodeRep<3,pztopology::TPZLine>::ClassId() << 1;
}
template class
TPZRestoreClass< TPZGeoElRefPattern<TPZArc3D>>;
