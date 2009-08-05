#include "tpzellipse.h"

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
using namespace pzshape;

void TPZEllipse::SetAxes(double a, double b)
{
    fA = a;
    fB = b;
}

/// Este método estabelece o mapeamento de um elemento linear contido num
/// hexaedro de dimensões (2.fA,2.fB,h) para a elipse inscrita no plano xy deste hexaedro
/// ou seja: elipse -> (semi-eixos fA e fB e extrudada no valor de h)
void TPZEllipse::X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result)
{
    //Mapeando o ponto do elemento mestre para a malha.
    TPZGeoLinear::X(nodes,loc,result);
	
    //Parametrizando para o elemento mestre do hexaedro (mantém coordenada z, i.e. result[2]).
    double  qsi = result[0]/fA;
    double  eta = result[1]/fB;
	
#if Debug
	if(fabs(qsi) > 1. || fabs(eta) > 1.)
	{
		cout << "\nThe linear element don't belong to the domain of the ellipse circumscribed hexahedron [(-fA,-fB,-h/2),(fA,fB,h/2)]\n";
		cout << "See TPZEllipse class!\n";
		exit(-1);
	}
#endif
	
    //Realizando o mapeamento no plano xy
    if(fabs(qsi) < 1.E-8 && fabs(eta) < 1.E-8)
    {
        result[0] = 0.;
        result[1] = 0.;
        return;
    }
    if( fabs(qsi) <= fabs(eta) )
    {
        result[0] = fA*fabs(eta)*qsi/sqrt(qsi*qsi + eta*eta);
        result[1] = fB*fabs(eta)*eta/sqrt(qsi*qsi + eta*eta);
        return;
    }
    if( fabs(qsi) > fabs(eta) )
    {
        result[0] = fA*fabs(qsi)*qsi/sqrt(qsi*qsi + eta*eta);
        result[1] = fB*fabs(qsi)*eta/sqrt(qsi*qsi + eta*eta);
        return;
    }
}

void TPZEllipse::Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv)
{
    TPZVec<REAL> result(3);
    TPZGeoLinear::X(coord,par,result);
    double  qsi = result[0]/fA;
    double  eta = result[1]/fB;
	
    TPZGeoLinear::Jacobian(coord, par, jacobian, axes, detjac, jacinv);
    TPZFMatrix JacLine(3,1);
	
    JacLine(0,0) = jacobian(0,0)*axes(0,0)/fA;
    JacLine(1,0) = jacobian(0,0)*axes(0,1)/fB;
    JacLine(2,0) = jacobian(0,0)*axes(0,2);
	
    TPZFMatrix JacEllipse(3,3);
	
    if(fabs(qsi) < 1.E-8 && fabs(eta) < 1.E-8)
    {
        JacEllipse(0,0) = fA;
        JacEllipse(0,1) = 0.;
        JacEllipse(0,2) = 0.;
		
        JacEllipse(1,0) = 0.;
        JacEllipse(1,1) = fB;
        JacEllipse(1,2) = 0.;
		
        JacEllipse(2,0) = 0.;
        JacEllipse(2,1) = 0.;
        JacEllipse(2,2) = 1.;
    }
    else if( fabs(qsi) <= fabs(eta) )
    {
        JacEllipse(0,0) =  fA*pow(eta*eta,1.5)/pow(eta*eta + qsi*qsi,1.5);
        JacEllipse(0,1) =  fA*eta*qsi*qsi*qsi/(fabs(eta)*pow(eta*eta + qsi*qsi,1.5));
        JacEllipse(0,2) =  0.;
		
        JacEllipse(1,0) = -fB*eta*fabs(eta)*qsi/pow(eta*eta + qsi*qsi,1.5);
        JacEllipse(1,1) =  fB*fabs(eta)*(eta*eta + 2.*qsi*qsi)/pow(eta*eta + qsi*qsi,1.5);
        JacEllipse(1,2) =  0.;
		
        JacEllipse(2,0) =  0.;
        JacEllipse(2,1) =  0.;
        JacEllipse(2,2) =  1.;
    }
    else if( fabs(qsi) > fabs(eta) )
    {
        JacEllipse(0,0) =  fA*fabs(qsi)*(2.*eta*eta + qsi*qsi)/pow(eta*eta + qsi*qsi,1.5);
        JacEllipse(0,1) = -fA*eta*qsi*fabs(qsi)/pow(eta*eta + qsi*qsi,1.5);
        JacEllipse(0,2) =  0.;
		
        JacEllipse(1,0) =  fB*eta*eta*eta*qsi/(fabs(qsi)*pow(eta*eta + qsi*qsi,1.5));
        JacEllipse(1,1) =  fB*pow(qsi*qsi,1.5)/pow(eta*eta + qsi*qsi,1.5);
        JacEllipse(1,2) =  0.;
		
        JacEllipse(2,0) =  0.;
        JacEllipse(2,1) =  0.;
        JacEllipse(2,2) =  1.;
    }
	
    TPZFMatrix JacTemp(3,1);
    JacEllipse.Multiply(JacLine,JacTemp);
    TPZFMatrix axest;
    JacTemp.GramSchmidt(axest,jacobian);
    axest.Transpose(&axes);
	
    detjac = jacobian(0,0);
    jacinv.Resize(1,1);
    jacinv(0,0) = 1./detjac;
}

void TPZEllipse::GetNodeCoord(TPZGeoMesh &mesh, TPZFMatrix &nodes){
	for(int i = 0; i < NNodes; i++){
		const int nodeindex = this->fNodeIndexes[i];
		TPZGeoNode * np = &(mesh.NodeVec()[nodeindex]);
		for(int j = 0; j < 3; j++){
			nodes(j,i) = np->Coord(j);
		}
	}
}///void


void TPZEllipse::AdjustNodeCoordinates(TPZGeoMesh &mesh){
	
	const int nnodes = NNodes;
	const int dimension = 1;
	TPZFNMatrix<3*nnodes> nodes(3,nnodes);
	this->GetNodeCoord(mesh,nodes);
	
	TPZManVector<REAL,3> qsi(dimension);
	TPZManVector<REAL,3> MappedX(3);
	
	for(int inode = 0; inode < nnodes; inode++){
		const int nodeindex = this->fNodeIndexes[inode];
		this->CenterPoint(inode, qsi);
		this->X(nodes,qsi,MappedX);
		for(int dim = 0; dim < 3; dim++){
			mesh.NodeVec()[nodeindex].SetCoord(dim,MappedX[dim]);
		}
		
	}///for inode
}///void

TPZGeoEl *TPZEllipse::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
{
	if(side==2)
	{
		TPZManVector<int> nodes(3);
		nodes[0] = orig->SideNodeIndex(side,0); nodes[1] = orig->SideNodeIndex(side,1); nodes[2] = orig->SideNodeIndex(side,2); int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,1)));
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

#include "tpzgeoelmapped.h"

///CreateGeoElement -> TPZEllipse
template< >
TPZGeoEl *TPZGeoElRefLess<TPZEllipse >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
	TPZGeoMesh &mesh = *(this->Mesh());
	if(!&mesh) return 0;
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

#define TPZGEOELEMENTELLIPSEID 30103
template<>
int TPZGeoElRefPattern<TPZEllipse>::ClassId() const
{
	return TPZGEOELEMENTELLIPSEID;
}
template class 
TPZRestoreClass< TPZGeoElRefPattern<TPZEllipse>, TPZGEOELEMENTELLIPSEID>;

template<>
TPZCompEl *(*TPZGeoElRefLess<TPZEllipse>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateLinearEl;

template class TPZGeoElRefLess<TPZEllipse>;
template class pzgeom::TPZNodeRep<2,TPZEllipse>;

