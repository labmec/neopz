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

#include "pzlog.h"

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeotriangle"));
#endif

using namespace pzshape;
using namespace std;

namespace pzgeom {
	
	const double tol = pzgeom_TPZNodeRep_tol;
	

	void TPZGeoTriangle::Shape(TPZVec<REAL> &param,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		REAL qsi = param[0], eta = param[1];
		phi(0,0) = 1.-qsi-eta;
		phi(1,0) = qsi;
		phi(2,0) = eta;
		dphi(0,0) = dphi(1,0) = -1.;
		dphi(0,1) = dphi(1,2) =  1.;
		dphi(1,1) = dphi(0,2) =  0.;
	}
	
	void TPZGeoTriangle::Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){
		
        int spacedim = coord.Rows();
        jacobian.Resize(2,2); axes.Resize(2,3); jacinv.Resize(2,2);
		TPZFNMatrix<3> phi(3,1);
        TPZFNMatrix<6> dphi(2,3),axest(3,2);
		jacobian.Zero();
		Shape(param,phi,dphi);
        TPZFNMatrix<6> VecMatrix(3,2,0.);
        for(int i = 0; i < 3; i++) {
			for(int j = 0; j < spacedim; j++) {
				VecMatrix(j,0) += coord.GetVal(j,i)*dphi(0,i);
				VecMatrix(j,1) += coord.GetVal(j,i)*dphi(1,i);
			}
        }
        VecMatrix.GramSchmidt(axest,jacobian);
        axest.Transpose(&axes);
		detjac = jacobian(0,0)*jacobian(1,1)-jacobian(1,0)*jacobian(0,1);
    REAL maxjac = 0.;
    for (int i=0; i<2; i++) {
      for (int j=0; j<2; j++) {
        maxjac = Max(maxjac,fabs(jacobian(i,j)));
      }
    }
        if(IsZero(maxjac) || IsZero(detjac/(maxjac*maxjac)))
		{
#ifdef PZDEBUG
			std::stringstream sout;
			sout << "Singular Jacobian " << detjac;
			LOGPZ_ERROR(logger, sout.str())
#endif
			detjac = ZeroTolerance();
		}
        
        jacinv(0,0) =  jacobian(1,1)/detjac;
        jacinv(1,1) =  jacobian(0,0)/detjac;
        jacinv(0,1) = -jacobian(0,1)/detjac;
        jacinv(1,0) = -jacobian(1,0)/detjac;
	}

	void TPZGeoTriangle::VecHdiv(TPZFMatrix<REAL> & coord, TPZFMatrix<REAL> & fNormalVec,TPZVec<int> &fVectorSide){
		if(coord.Rows()!=3)
		{
			cout<< "Erro na dimensao das linhas de coord"<< endl;
		}
		if(coord.Cols()!=3)
		{
			cout<< "Erro na dimensao das colunas de coord"<< endl;
		}
		TPZVec<REAL> p1(3), p2(3), p3(3),result(3);
		for(int j=0;j<3;j++)
		{
			p1[j]=coord.GetVal(j,0);
			p2[j]=coord.GetVal(j,1);
			p3[j]=coord.GetVal(j,2);
		}
		fNormalVec.Resize(14, 3);
		fVectorSide.Resize(14);
		int64_t count=0;
		
		//primeira face
		for(int j=0;j<3;j++)//v0
		{
			fNormalVec(0,j) = coord.GetVal(j,0)- coord.GetVal(j,2);
		}
		fVectorSide[count]=0;
		count++;
		for(int j=0;j<3;j++)//v1
		{
			fNormalVec(1,j) = coord.GetVal(j,1)- coord.GetVal(j,2);
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
			fNormalVec(3,j) = coord.GetVal(j,1)- coord.GetVal(j,0);
		}
		fVectorSide[count]=1;
		count++;
		for(int j=0;j<3;j++)//v4
		{
			fNormalVec(4,j) = coord.GetVal(j,2)- coord.GetVal(j,0);
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
			fNormalVec(6,j) = coord.GetVal(j,2)- coord.GetVal(j,1);
		}
		fVectorSide[count]=2;
		count++;
		for(int j=0;j<3;j++)//v7
		{
			fNormalVec(7,j) = coord.GetVal(j,0)- coord.GetVal(j,1);
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
			fNormalVec(9,j) = coord.GetVal(j,1)- coord.GetVal(j,0);
		}  
		fVectorSide[count]=3;
		count++;
		for(int j=0;j<3;j++)//v10
		{
			fNormalVec(10,j) = coord.GetVal(j,2)- coord.GetVal(j,1);
		}	
		fVectorSide[count]=4;
		count++;
		
		for(int j=0;j<3;j++)//v11
		{
			fNormalVec(11,j) = coord.GetVal(j,0)- coord.GetVal(j,2);
		}	
		fVectorSide[count]=5;
		count++;
		//internos meio
		TPZVec<REAL> midle(3,0.);
		midle[0]=(1./3.)*(coord.GetVal(0,2)+coord.GetVal(0,0)+coord.GetVal(0,1));
		midle[1]=(1./3.)*(coord.GetVal(1,2)+coord.GetVal(1,0)+coord.GetVal(1,1));		
		midle[2]=(1./3.)*(coord.GetVal(2,2)+coord.GetVal(2,0)+coord.GetVal(2,1));
		TPZFMatrix<REAL> jacobian;
		TPZFMatrix<REAL> axes;
		TPZFMatrix<REAL> jacinv;
        DebugStop();
		//Jacobian(coord,midle,jacobian,axes,detjac,jacinv);
		fNormalVec(12,0)=axes(0,0);
		fNormalVec(12,1)=axes(0,1);
		fNormalVec(12,2)=axes(0,2);
		fNormalVec(13,0)=axes(1,0);
		fNormalVec(13,1)=axes(1,1);
		fNormalVec(13,2)=axes(1,2);
		fVectorSide[count]=6;
		fVectorSide[count+1]=6;
		//normalizao
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
        if (logger->isDebugEnabled())
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
	
	TPZGeoEl *TPZGeoTriangle::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
        if(side==6) {
			TPZManVector<int64_t> nodes(3);
			int i;
			for (i=0;i<3;i++){
				nodes[i] = orig->SideNodeIndex(side,i);
			}
			int64_t index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
			int iside;
			for (iside = 0; iside <6; iside++){
				TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::ContainedSideLocId(side,iside)));
			}
			TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if(side>-1 && side<3) {
			TPZManVector<int64_t> nodeindexes(1);
			nodeindexes[0] = orig->SideNodeIndex(side,0);
			int64_t index;
			TPZGeoEl *gel = orig->CreateGeoElement(EPoint,nodeindexes,bc,index);
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if(side > 2 && side < 6) {
			TPZManVector<int64_t> nodes(2);
			nodes[0] = orig->SideNodeIndex(side,0);
			nodes[1] = orig->SideNodeIndex(side,1);
			int64_t index;
			TPZGeoEl *gel = orig->CreateGeoElement(EOned,nodes,bc,index);
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
	
	/** Creates a geometric element according to the type of the father element */
	TPZGeoEl *TPZGeoTriangle::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											   TPZVec<int64_t>& nodeindexes,
											   int matid,
											   int64_t& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
    
    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoTriangle::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<REAL,3> co(3),shift(3),scale(3);
        TPZManVector<int64_t,3> nodeindexes(3);
        for (int i=0; i<3; i++) {
            scale[i] = size[i]/3.;
            shift[i] = 1./2.+lowercorner[i];
        }
        
        for (int i=0; i<NCornerNodes; i++) {
            ParametricDomainNodeCoord(i, co);
            for (int j=0; j<co.size(); j++) {
                co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
            }
            nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
        }
        int64_t index;
        CreateGeoElement(gmesh, ETriangle, nodeindexes, matid, index);
    }
    
    int TPZGeoTriangle::ClassId() const{
        return Hash("TPZGeoTriangle") ^ TPZNodeRep<3, pztopology::TPZTriangle>::ClassId() << 1;
    }

    void TPZGeoTriangle::Read(TPZStream& buf, void* context) {
        TPZNodeRep<3, pztopology::TPZTriangle>::Read(buf, context);
    }

    void TPZGeoTriangle::Write(TPZStream& buf, int withclassid) const {
        TPZNodeRep<3, pztopology::TPZTriangle>::Write(buf, withclassid);
    }


};
