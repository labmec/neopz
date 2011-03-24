// TPZGeoTriangle.c: implementation of the TPZGeoTriangle class.
//
//////////////////////////////////////////////////////////////////////

#include "pzgeotriangle.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzshapetriang.h"
#include "pzgmesh.h"
//#include "pzgeoelrefless.h.h"
#include "pzlog.h"

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeotriangle"));
#endif

using namespace pzshape;
using namespace std;

namespace pzgeom {

const double tol = pzgeom_TPZNodeRep_tol;

void TPZGeoTriangle::Shape(TPZVec<REAL> &param,TPZFMatrix &phi,TPZFMatrix &dphi) {
	REAL qsi = param[0], eta = param[1];
	phi(0,0) = 1.-qsi-eta;
	phi(1,0) = qsi;
	phi(2,0) = eta;
	dphi(0,0) = dphi(1,0) = -1.;
	dphi(0,1) = dphi(1,2) =  1.;
	dphi(1,1) = dphi(0,2) =  0.;
}

void TPZGeoTriangle::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

        int spacedim = coord.Rows();
        jacobian.Resize(2,2); axes.Resize(2,3); jacinv.Resize(2,2);
	TPZFNMatrix<3> phi(3,1);
        TPZFNMatrix<6> dphi(2,3),axest(3,2);
	jacobian.Zero();
	Shape(param,phi,dphi);
        TPZFNMatrix<6> VecMatrix(3,2,0.);
        for(int i = 0; i < 3; i++) {
          for(int j = 0; j < spacedim; j++) {
              VecMatrix(j,0) += coord(j,i)*dphi(0,i);
              VecMatrix(j,1) += coord(j,i)*dphi(1,i);
          }
        }
        VecMatrix.GramSchmidt(axest,jacobian);
        axest.Transpose(&axes);
	detjac = jacobian(0,0)*jacobian(1,1)-jacobian(1,0)*jacobian(0,1);
    if(IsZero(detjac))
    {
        std::stringstream sout;
        sout << "Singular Jacobian " << detjac;
        LOGPZ_ERROR(logger, sout.str())
        detjac = ZeroTolerance();
    }
        if(detjac)
        {
          jacinv(0,0) =  jacobian(1,1)/detjac;
          jacinv(1,1) =  jacobian(0,0)/detjac;
          jacinv(0,1) = -jacobian(0,1)/detjac;
          jacinv(1,0) = -jacobian(1,0)/detjac;
        }
        else
        {
          jacinv.Zero();
        }
}

void TPZGeoTriangle::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){

	REAL spacephi[3],spacedphi[6];
	TPZFMatrix phi(3,1,spacephi,3);
	TPZFMatrix dphi(2,3,spacedphi,6);
	Shape(loc,phi,dphi);
        int space = coord.Rows();

	for(int i = 0; i < space; i++) {
               result[i] = 0.0;
               for(int j = 0; j < 3; j++) result[i] += phi(j,0)*coord(i,j);
	}
}

bool TPZGeoTriangle::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide) {

     double zero = 1.E-5;
     REAL qsi = InternalPar[0]; REAL eta = InternalPar[1];
     SidePar.Resize(1); JacToSide.Resize(1,2);
     if((qsi + eta - 1.) > 1.e-5 || qsi < -1.e-5 || eta < -1.e-5)
     {
         cout << "Point (qsi,eta) = (" << qsi << "," << eta << ") is out of TPZGeoTriangle Master Element Range!\n";
         cout << "See TPZGeoTriangle::MapToSide() method!\n";
		 DebugStop();
         exit(-1);
     }
	if(qsi < 0.) qsi = 0.;
	if(eta < 0.) eta = 0.;
	if(qsi+eta > 1.)
	{
		// O BUG ESTAVA AQUI!!!
		REAL qsieta = 1.-qsi-eta;
		qsi += qsieta/2.;
		eta += qsieta/2.;
	}
	bool regularmap = true;

     switch(side)
     {
          case 3:
               if(fabs(eta - 1.) < zero)
               {
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.;
				   regularmap = false;
               }
               else
               {
                    SidePar[0] = 2.*qsi/(1.-eta) - 1.;
                    JacToSide(0,0) = 2./(1.-eta); JacToSide(0,1) = 2.*qsi/((1.-eta)*(1.-eta));
               }
          break;

          case 4:
               if(qsi+eta < zero)
               {
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.;
				   regularmap = false;
               }
               else
               {
                    SidePar[0] = 1. - 2.*qsi/(qsi + eta);
                    JacToSide(0,0) = -2.*eta/((qsi+eta)*(qsi+eta)); JacToSide(0,1) = 2.*qsi/((qsi+eta)*(qsi+eta));
               }
          break;

          case 5:
               if(fabs(qsi - 1.) < zero)
               {
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.;
				   regularmap = false;
               }
               else
               {
                    SidePar[0] = 1. - 2.*eta/(1.-qsi);
                    JacToSide(0,0) = -2.*eta/((1.-qsi)*(1.-qsi)); JacToSide(0,1) = -2./(1.-qsi);
               }
          break;
     }
     if(side < 3 || side > 5)
     {
		 cout << "Cant compute MapToSide method in TPZGeoTriangle class!\nParameter (SIDE) must be 3, 4 or 5!\nMethod Aborted!\n"; 
		 DebugStop();
		 exit(-1);
     }
	return regularmap;
}

TPZGeoEl *TPZGeoTriangle::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
        if(side==6) {
		TPZManVector<int> nodes(3);
		int i;
		for (i=0;i<3;i++){
			nodes[i] = orig->SideNodeIndex(side,i);
		}
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
		int iside;
		for (iside = 0; iside <6; iside++){
			TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,iside)));
		}
		TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	else if(side>-1 && side<3) {
		TPZManVector<int> nodeindexes(1);
		nodeindexes[0] = orig->SideNodeIndex(side,0);
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	else if(side > 2 && side < 6) {
		TPZManVector<int> nodes(2);
		nodes[0] = orig->SideNodeIndex(side,0);
		nodes[1] = orig->SideNodeIndex(side,1);
		int index;
		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,0)));
		TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,1)));
		TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
		return gel;
	}
	else PZError << "TPZGeoTriangle::CreateBCGeoEl has no bc.\n";
	return 0;
}

void TPZGeoTriangle::FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint)
{
    ChangedPoint.Resize(OriginalPoint.NElements(),0.);
    ChangedPoint = OriginalPoint;

    switch(side)
    {
        case 3:
        {
            if(fabs(OriginalPoint[0]) <= tol && fabs(OriginalPoint[1]- 1.) <= tol)
            {
                ChangedPoint[0] = tol;
                ChangedPoint[1] = 1. - 2.*tol;
            }
            break;
        }

        case 4:
        {
            if(fabs(OriginalPoint[0]) <= tol && fabs(OriginalPoint[1]) <= tol)
            {
                ChangedPoint[0] = tol;
                ChangedPoint[1] = tol;
            }
            break;
        }

        case 5:
        {
            if(fabs(OriginalPoint[0] - 1.) <= tol && fabs(OriginalPoint[1]) <= tol)
            {
                ChangedPoint[0] = 1.-tol;
                ChangedPoint[1] = tol/2.;
            }
            break;
        }
    }
}

	/**
	 * Creates a geometric element according to the type of the father element
	 */
	TPZGeoEl *TPZGeoTriangle::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											TPZVec<int>& nodeindexes,
											int matid,
											int& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
	
};
