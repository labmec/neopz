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

#include <array>

static std::array<REAL,3> Cross(const std::array<REAL,3>& a,
                           const std::array<REAL,3>& b)
{
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

#include <numeric>   // std::inner_product
#include <cmath>     // std::sqrt

REAL Norm(const std::array<REAL,3>& v)
{
    REAL sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    return std::sqrt(sum);
}

void Normalize(std::array<REAL,3>& v)
{
	REAL norm = Norm(v);
	for(auto &it : v) it /= norm;
}

void TPZArc3D::ComputeAtributes(TPZFMatrix<REAL> &coord)
{
#ifdef PZDEBUG
	/** verify is the points are not colinear */
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
	ComputeCenter(coord,fCenter3D);
	std::array<REAL,3> ax0;
	std::array<REAL,3> ax1;
	std::array<REAL,3> ax2;
	std::array<REAL,3> v1,v2;

	for(int i = 0; i < 3; i++)
	{
		ax0[i] = coord(i,0)-fCenter3D[i];
		ax1[i] = coord(i,1)-fCenter3D[i];
		v2[i] = coord(i,2)-fCenter3D[i];
	}
	v1 = ax1;
	fRadius = Norm(ax0);
	Normalize(ax0);
	ax2 = Cross(ax0,ax1);
	Normalize(ax2);
	/// ax2 should already have norm 1
	ax1 = Cross(ax2,ax0);
	for(int i=0; i<3; i++) {
		fRotationTensor(i,0) = ax0[i];
		fRotationTensor(i,1) = ax1[i];
		fRotationTensor(i,2) = ax2[i];
	}
	
	/** Computing Radius */
	fRadius = Norm(v1);
	
	/** Computing (Center->First_Point) vector */
	REAL cosalpha = std::inner_product(v1.begin(), v1.end(), ax0.begin(), 0.0);
    REAL sinalpha = std::inner_product(v1.begin(), v1.end(), ax1.begin(), 0.0);
	fAngle = atan2(sinalpha,cosalpha);

	REAL cosalpham = std::inner_product(v2.begin(), v2.end(), ax0.begin(), 0.0);
    REAL sinalpham = std::inner_product(v2.begin(), v2.end(), ax1.begin(), 0.0);
	REAL Anglem = atan2(sinalpham,cosalpham);
	if(!(fAngle*Anglem > 0 && Anglem/fAngle < 1.)) {
		if(fAngle < 0) fAngle += M_2_PI;
		else if(fAngle > 0) fAngle -= M_2_PI;
	}
    
    
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

void TPZArc3D::ComputeCenter(TPZFMatrix<REAL> &coord, TPZVec<REAL> &center) 
{
	std::array<REAL,3> cross;
	std::array<REAL,3> v0 = {coord(0,0),coord(1,0),coord(2,0)};
	std::array<REAL,3> v1 = {coord(0,1),coord(1,1),coord(2,1)};
	cross = Cross(v0,v1);
	TPZFNMatrix<9,REAL> Mat(3,3), F(3,1, 0.);
	for(int i=0; i<3; i++) {
		Mat(0,i) = coord(i,0)-coord(i,2);
		Mat(1,i) = coord(i,1)-coord(i,2);
		Mat(2,i) = cross[i];
	}
	for(int i=0; i<3; i++) {
		F(0,0) += (coord(i,0)+coord(i,2))*Mat(0,i)/2.;
		F(1,0) += (coord(i,1)+coord(i,2))*Mat(1,i)/2.;
		F(2,0) += coord(i,2)*cross[i];
	}
	Mat.Solve_LU(&F);
	for(int i = 0; i<3; i++) {
		center[i] = F(i,0);
	}

}


int TPZArc3D::ClassId() const{
    return Hash("TPZArc3D") ^ pzgeom::TPZNodeRep<3,pztopology::TPZLine>::ClassId() << 1;
}
template class
TPZRestoreClass< TPZGeoElRefPattern<TPZArc3D>>;
