/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticTetra methods. 
 */
#include "tpzquadratictetra.h"
#include "pzshapetetra.h"
#include "tpzgeoelmapped.h"

#include "pzlog.h"

#include "tpzgeomid.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.specialmaps.quadratictetra"));
#endif

using namespace pzgeom;
using namespace pztopology;
using namespace pzshape;
using namespace std;

TPZQuadraticTetra::~TPZQuadraticTetra()
{
}

void TPZQuadraticTetra::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{
	REAL qsi = pt[0], eta = pt[1], zeta = pt[2];
	
	phi(0,0)  =  (qsi + eta + zeta -1.) * (2.*qsi + 2.*eta + 2.*zeta - 1.);
	phi(1,0)  =  qsi  * (2.*qsi - 1.);
	phi(2,0)  =  eta  * (2.*eta - 1.);
	phi(3,0)  =  zeta * (2.*zeta - 1.);
	phi(4,0)  = -4.*qsi  * (qsi + eta + zeta -1.);
	phi(5,0)  =  4.*qsi  *  eta;
	phi(6,0)  = -4.*eta  * (qsi + eta + zeta -1.);
	phi(7,0)  = -4.*zeta * (qsi + eta + zeta -1.);
	phi(8,0)  =  4.*qsi  *  zeta;
	phi(9,0)  =  4.*eta  *  zeta;
	
	dphi(0,0) =  -3. + 4.*qsi + 4.*eta + 4.*zeta;
	dphi(1,0) =  -3. + 4.*qsi + 4.*eta + 4.*zeta;
	dphi(2,0) =  -3. + 4.*qsi + 4.*eta + 4.*zeta;
	dphi(0,1) =  -1. + 4.*qsi;
	dphi(1,1) =   0.;
	dphi(2,1) =   0.;
	dphi(0,2) =   0.;
	dphi(1,2) =  -1. + 4.*eta;
	dphi(2,2) =   0.;
	dphi(0,3) =   0.;
	dphi(1,3) =   0.;
	dphi(2,3) =  -1. + 4.*zeta;
	dphi(0,4) =  -4.*(-1. + 2.*qsi + eta + zeta);
	dphi(1,4) =  -4.*qsi;
	dphi(2,4) =  -4.*qsi;
	dphi(0,5) =   4.*eta;
	dphi(1,5) =   4.*qsi;
	dphi(2,5) =   0.;
	dphi(0,6) =  -4.*eta;
	dphi(1,6) =  -4.*(-1. + qsi + 2.*eta + zeta);
	dphi(2,6) =  -4.*eta;
	dphi(0,7) =  -4.*zeta;
	dphi(1,7) =  -4.*zeta;
	dphi(2,7) =  -4.*(-1. + qsi + eta + 2.*zeta);
	dphi(0,8) =   4.*zeta;
	dphi(1,8) =   0.;
	dphi(2,8) =   4.*qsi;
	dphi(0,9) =   0.;
	dphi(1,9) =   4.*zeta;
	dphi(2,9) =   4.*eta;
}

void TPZQuadraticTetra::X(TPZFMatrix<REAL> & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result)
{
	TPZFNMatrix<10> phi(10,1);
	TPZFNMatrix<30> dphi(3,10); Shape(loc,phi,dphi);
	for(int i = 0; i < 3; i++)
	{
		result[i] = 0.0;
		for(int j = 0; j < 10; j++) result[i] += phi(j,0)*coord(i,j);
	}
}

void TPZQuadraticTetra::Jacobian(TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv)
{
#ifdef DEBUG
	int nnodes = NNodes;
	if (nnodes != 10) { PZError << "TPZGeoTetrahedra.jacobian only implemented for\n10 nodes, NumberOfNodes = " << nnodes << "\n"; }
#endif
	
	jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
    jacobian.Zero(); axes.Zero();
    for(int d = 0; d < 3; d++) axes(d,d) = 1.;
    
	TPZFNMatrix<10> phi(10,1);
	TPZFNMatrix<30> dphi(3,10);
	Shape(param,phi,dphi);	
	for(int i = 0; i < 10; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			jacobian(j,0) += coord(j,i)*dphi(0,i);
			jacobian(j,1) += coord(j,i)*dphi(1,i);
			jacobian(j,2) += coord(j,i)*dphi(2,i);
		}
	}
	
	detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0);//- a02 a11 a20
	detjac += jacobian(0,1)*jacobian(1,2)*jacobian(2,0);//+ a01 a12 a20
	detjac += jacobian(0,2)*jacobian(1,0)*jacobian(2,1);//+ a02 a10 a21
	detjac -= jacobian(0,0)*jacobian(1,2)*jacobian(2,1);//- a00 a12 a21
	detjac -= jacobian(0,1)*jacobian(1,0)*jacobian(2,2);//- a01 a10 a22
	detjac += jacobian(0,0)*jacobian(1,1)*jacobian(2,2);//+ a00 a11 a22
    
    if(IsZero(detjac))
    {
#ifdef DEBUG
        std::stringstream sout;
        sout << "Singular Jacobian " << detjac;
        LOGPZ_ERROR(logger, sout.str())
#endif
        detjac = ZeroTolerance();
    }
	
	jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2))/detjac;//-a12 a21 + a11 a22
	jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))/detjac;// a02 a21 - a01 a22
	jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2))/detjac;//-a02 a11 + a01 a12
	jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))/detjac;// a12 a20 - a10 a22
	jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2))/detjac;//-a02 a20 + a00 a22
	jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))/detjac;// a02 a10 - a00 a12
	jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1))/detjac;//-a11 a20 + a10 a21
	jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))/detjac;// a01 a20 - a00 a21
	jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1))/detjac;//-a01 a10 + a00 a11
}


TPZGeoEl *TPZQuadraticTetra::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc)
{
	if(side < 0 || side > 14) { cout << "TPZGeoTetrahedra::CreateBCCompEl with bad side = " << side << "not implemented\n"; return 0; }
	if(side == 14) { cout << "TPZGeoTetrahedra::CreateBCCompEl with side = 14 not implemented\n"; return 0; }
	if(side < 4)
	{
		TPZManVector<long> nodeindexes(1);
		nodeindexes[0] = orig->NodeIndex(side);
		long index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	
	else if (side > 3 && side < 10)
	{// side = 4 a 9 : lados
		TPZManVector<long> nodes(2);
		nodes[0] = orig->SideNodeIndex(side,0);
		nodes[1] = orig->SideNodeIndex(side,1); 
		long index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::ContainedSideLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::ContainedSideLocId(side,1)));
		TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	
	else if (side > 9)
	{//side = 10 a 13 : faces
		TPZManVector<long> nodes(3); int in;
		
		for (in=0;in<3;in++) nodes[in] = orig->SideNodeIndex(side,in);
		long index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
		
		for (in=0;in<6;in++) TPZGeoElSide(gel,in).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::ContainedSideLocId(side,in)));
		TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	} 
	else PZError << "TPZGeoTetrahedra::CreateBCGeoEl. Side = " << side << endl;
	return 0;
}


/**
 * Creates a geometric element according to the type of the father element
 */

TPZGeoEl *TPZQuadraticTetra::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											  TPZVec<long>& nodeindexes,
											  int matid,
											  long& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

//void TPZQuadraticTetra::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
//{
//    if(node > this->NNodes)
//    {
//        DebugStop();
//    }
//    nodeCoord.Resize(Dimension, 0.);
//    switch (node) {
//        case (0):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 0.;
//            nodeCoord[2] = 0.;
//            break;
//        }
//        case (1):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 0.;
//            nodeCoord[2] = 0.;
//            break;
//        }
//        case (2):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 1.;
//            nodeCoord[2] = 0.;
//            break;
//        }
//        case (3):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 0.;
//            nodeCoord[2] = 1.;
//            break;
//        }
//        case (4):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.0;
//            nodeCoord[2] =  0.0;
//            break;
//        }
//        case (5):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  0.0;
//            break;
//        }
//        case (6):
//        {
//            nodeCoord[0] =  0.0;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  0.0;
//            break;
//        }
//        case (7):
//        {
//            nodeCoord[0] =  0.0;
//            nodeCoord[1] =  0.0;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        case (8):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.0;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        case (9):
//        {
//            nodeCoord[0] =  0.0;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        default:
//        {
//            DebugStop();
//            break;
//        }
//    }
//}

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"


///CreateGeoElement -> TPZQuadraticTetra

template<>
int TPZGeoElRefPattern<TPZQuadraticTetra>::ClassId() const {
	return TPZGEOELEMENTQUADRATICTETRAID;
}

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticTetra>, TPZGEOELEMENTQUADRATICTETRAID>;

template class TPZGeoElRefPattern<TPZQuadraticTetra>;
template class TPZGeoElRefLess<TPZQuadraticTetra>;
template class pzgeom::TPZNodeRep<10,TPZQuadraticTetra>;
