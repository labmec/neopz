/**
 * @file
 * @brief Contains the implementation of the TPZGeoTriangle methods. 
 */

#include "pzgeotriangle.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzshapetriang.h"
#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

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
	void TPZGeoTriangle::VecHdiv(TPZFMatrix & coord, TPZFMatrix & fNormalVec,TPZVec<int> &fVectorSide){
		if(coord.Rows()!=3)
		{
			cout<< "Erro na dimensão das linhas de coord"<< endl;
		}
		if(coord.Cols()!=3)
		{
			cout<< "Erro na dimensão das colunas de coord"<< endl;
		}
		TPZVec<REAL> p1(3), p2(3), p3(3),result(3);
		for(int j=0;j<3;j++)
		{
			p1[j]=coord(j,0);
			p2[j]=coord(j,1);
			p3[j]=coord(j,2);
		}
		fNormalVec.Resize(14, 3);
		fVectorSide.Resize(14);
		int count=0;
		
		//primeira face
		for(int j=0;j<3;j++)//v0
		{
			fNormalVec(0,j) = coord(j,0)- coord(j,2);
		}
		fVectorSide[count]=0;
		count++;
		for(int j=0;j<3;j++)//v1
		{
			fNormalVec(1,j) = coord(j,1)- coord(j,2);
		}
		fVectorSide[count]=1;
		count++;
		//v2
		ComputeNormal(p1,p2,p3,result);
		fNormalVec(2,0) = -result[0];
		fNormalVec(2,1) = -result[1];
		fNormalVec(2,2) = -result[2];
		fVectorSide[count]=3;
		count++;
		//segunda face
		for(int j=0;j<3;j++)//v3
		{
			fNormalVec(3,j) = coord(j,1)- coord(j,0);
		}
		fVectorSide[count]=1;
		count++;
		for(int j=0;j<3;j++)//v4
		{
			fNormalVec(4,j) = coord(j,2)- coord(j,0);
		}
		fVectorSide[count]=2;
		count++;
		//v5
		ComputeNormal(p2,p3,p1,result);
		fNormalVec(5,0) = -result[0];
		fNormalVec(5,1) = -result[1];
		fNormalVec(5,2) = -result[2];
		fVectorSide[count]=4;
		count++;
		//terceira face
		for(int j=0;j<3;j++)//v6
		{
			fNormalVec(6,j) = coord(j,2)- coord(j,1);
		}
		fVectorSide[count]=2;
		count++;
		for(int j=0;j<3;j++)//v7
		{
			fNormalVec(7,j) = coord(j,0)- coord(j,1);
		}
		fVectorSide[count]=0;
		count++;
		//v8
		ComputeNormal(p3,p1,p2,result);
		fNormalVec(8,0) = -result[0];
		fNormalVec(8,1) = -result[1];
		fNormalVec(8,2) = -result[2];
		fVectorSide[count]=5;
		count++;
		// internos tangentes
		for(int j=0;j<3;j++)//v9
		{
			fNormalVec(9,j) = coord(j,1)- coord(j,0);
		}  
		fVectorSide[count]=3;
		count++;
		for(int j=0;j<3;j++)//v10
		{
			fNormalVec(10,j) = coord(j,2)- coord(j,1);
		}	
		fVectorSide[count]=4;
		count++;
		
		for(int j=0;j<3;j++)//v11
		{
			fNormalVec(11,j) = coord(j,0)- coord(j,2);
		}	
		fVectorSide[count]=5;
		count++;
		//internos meio
		TPZVec<REAL> midle(3,0.);
		midle[0]=(1./3.)*(coord(0,2)+coord(0,0)+coord(0,1));
		midle[1]=(1./3.)*(coord(1,2)+coord(1,0)+coord(1,1));		
		midle[2]=(1./3.)*(coord(2,2)+coord(2,0)+coord(2,1));
		TPZFMatrix jacobian;
		TPZFMatrix axes;
		REAL detjac;
		TPZFMatrix jacinv;
		Jacobian(coord,midle,jacobian,axes,detjac,jacinv);
		fNormalVec(12,0)=axes(0,0);
		fNormalVec(12,1)=axes(0,1);
		fNormalVec(12,2)=axes(0,2);
		fNormalVec(13,0)=axes(1,0);
		fNormalVec(13,1)=axes(1,1);
		fNormalVec(13,2)=axes(1,2);
		fVectorSide[count]=6;
		fVectorSide[count+1]=6;
		//normalização
		for(int k=0;k<14;k++)
		{
			REAL temp=0.;
			temp=sqrt( fNormalVec(k,0)*fNormalVec(k,0) + fNormalVec(k,1)*fNormalVec(k,1) + fNormalVec(k,2)*fNormalVec(k,2));
			fNormalVec(k,0) *=1./temp;	
			fNormalVec(k,1) *=1./temp;	
		}
		// produto normal == 1
		for(int kk=0;kk<3;kk++)
		{
			REAL temp1=0.;
			REAL temp2=0.;
			temp1 =  fNormalVec(kk*3,0)*fNormalVec(kk*3+2,0) + fNormalVec(kk*3,1)*fNormalVec(kk*3+2,1);
			temp2 =  fNormalVec(kk*3+1,0)*fNormalVec(kk*3+2,0) + fNormalVec(kk*3+1,1)*fNormalVec(kk*3+2,1);
			fNormalVec(kk*3,0) *=1./temp1;	
			fNormalVec(kk*3,1) *=1./temp1;
			fNormalVec(kk*3+1,0) *=1./temp2;
			fNormalVec(kk*3+1,1) *=1./temp2;
		}	
#ifdef LOG4CXX
		{
			std::stringstream sout;
			fNormalVec.Print("fNormalVec", sout);		
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif 
		
	}
	
	void TPZGeoTriangle::VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result){
		if(v1.NElements()!=3||v2.NElements()!=3)
		{
			cout << " o tamanho do vetores eh diferente de 3"<< endl;
		}
		REAL x1=v1[0], y1=v1[1],z1=v1[2];
		REAL x2=v2[0], y2=v2[1],z2=v2[2];
		result.Resize(v1.NElements());
		result[0]=y1*z2-z1*y2;
		result[1]=z1*x2-x1*z2;	
		result[2]=x1*y2-y1*x2;	
	}
	
	void TPZGeoTriangle::ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result){
		TPZVec<REAL> v1(3);
		TPZVec<REAL> v2(3);
		TPZVec<REAL> normal(3);
		v1[0]=p1[0]-p2[0];
		v1[1]=p1[1]-p2[1];
		v1[2]=p1[2]-p2[2];
		v2[0]=p2[0]-p3[0];
		v2[1]=p2[1]-p3[1];
		v2[2]=p2[2]-p3[2];
		VectorialProduct(v1,v2,normal);
		VectorialProduct(v1,normal,result);	
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
		}
		if(qsi < 0.) qsi = 0.;
		if(eta < 0.) eta = 0.;
		if(qsi+eta > 1.)
		{
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
				TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::ContainedSideLocId(side,iside)));
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
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::ContainedSideLocId(side,0)));
			TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::ContainedSideLocId(side,1)));
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
