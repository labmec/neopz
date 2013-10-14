/**
 * @file
 * @brief Contains the implementation of the TPZGeoQuad methods. 
 */

#include "pzgeoquad.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeoquad"));
#endif

using namespace std;

namespace pzgeom {
	
	void TPZGeoQuad::Shape(TPZVec<REAL> &param,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
		REAL x=param[0], y=param[1];
		
		phi(0,0) = .25*(1.-x)*(1.-y);
		phi(1,0) = .25*(1.+x)*(1.-y);
		phi(2,0) = .25*(1.+x)*(1.+y);
		phi(3,0) = .25*(1.-x)*(1.+y);
		
		dphi(0,0) = .25*(y-1.);
		dphi(1,0) = .25*(x-1.);
		
		dphi(0,1) = .25*(1.-y);
		dphi(1,1) =-.25*(1.+x);
		
		dphi(0,2) = .25*(1.+y);
		dphi(1,2) = .25*(1.+x);
		
		dphi(0,3) =-.25*(1.+y);
		dphi(1,3) = .25*(1.-x);
	}
	//coord é uma matrix 3x4
	void TPZGeoQuad::VecHdiv(TPZFMatrix<REAL> & coord, TPZFMatrix<REAL> & fNormalVec,TPZVec<int> & fVectorSide){
		if(coord.Rows()!=3)
		{
			cout<< "Erro na dimensao das linhas de coord"<< endl;
		}
		if(coord.Cols()!=4)
		{
			cout<< "Erro na dimensao das colunas de coord"<< endl;
		}
		TPZVec<REAL> p1(3), p2(3), p3(3), p4(3),result(3);
		for(int j=0;j<3;j++)
		{
			p1[j]=coord.GetVal(j,0);
			p2[j]=coord.GetVal(j,1);
			p3[j]=coord.GetVal(j,2);
			p4[j]=coord.GetVal(j,3);
		}
		fNormalVec.Resize(18, 3);
		fVectorSide.Resize(18);
		long count=0;
		//primeira face
		for(int j=0;j<3;j++)//v0
		{
			fNormalVec(0,j) = coord.GetVal(j,0)- coord.GetVal(j,3);
			
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
		fVectorSide[count]=4;
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
			fNormalVec(4,j) = coord.GetVal(j,2)- coord.GetVal(j,3);
		}
		fVectorSide[count]=2;
		count++;
		//v5
		ComputeNormal(p2,p3,p4,result);
		fNormalVec(5,0) = -result[0];
		fNormalVec(5,1) = -result[1];
		fNormalVec(5,2) = -result[2];
		fVectorSide[count]=5;
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
			fNormalVec(7,j) = coord.GetVal(j,3)- coord.GetVal(j,0);
		}
		fVectorSide[count]=3;
		count++;
		//v8
		ComputeNormal(p3,p4,p1,result);
		fNormalVec(8,0) = -result[0];
		fNormalVec(8,1) = -result[1];
		fNormalVec(8,2) = -result[2];
		fVectorSide[count]=6;
		count++;
		//quarta face
		for(int j=0;j<3;j++)//v9
		{
			fNormalVec(9,j) = coord.GetVal(j,3)- coord.GetVal(j,2);
		}
		fVectorSide[count]=3;
		count++;
		for(int j=0;j<3;j++)//v10
		{
			fNormalVec(10,j) = coord.GetVal(j,0)- coord.GetVal(j,1);
		}
		fVectorSide[count]=0;
		count++;
		//v11
		ComputeNormal(p4,p1,p2,result);
		fNormalVec(11,0) = -result[0];
		fNormalVec(11,1) = -result[1];
		fNormalVec(11,2) = -result[2];
		fVectorSide[count]=7;
		count++;
		
		// internos tangentes
		for(int j=0;j<3;j++)//v12
		{
			fNormalVec(12,j) = coord.GetVal(j,1)- coord.GetVal(j,0);
		}  
		fVectorSide[count]=4;
		count++;
		for(int j=0;j<3;j++)//v13
		{
			fNormalVec(13,j) = coord.GetVal(j,2)- coord.GetVal(j,1);
		}	
		fVectorSide[count]=5;
		count++;
		for(int j=0;j<3;j++)//v14
		{
			fNormalVec(14,j) = coord.GetVal(j,3)- coord.GetVal(j,2);
		}	
		fVectorSide[count]=6;
		count++;
		for(int j=0;j<3;j++)//v15
		{
			fNormalVec(15,j) = coord.GetVal(j,0)- coord.GetVal(j,3);
		}	
		fVectorSide[count]=7;
		count++;
		//internos meio
		TPZVec<REAL> midle(3,0.);
		midle[0]=0.25*(coord.GetVal(0,2)+coord.GetVal(0,3)+coord.GetVal(0,0)+coord.GetVal(0,1));
		midle[1]=0.25*(coord.GetVal(1,2)+coord.GetVal(1,3)+coord.GetVal(1,0)+coord.GetVal(1,1));		
		midle[2]=0.25*(coord.GetVal(2,2)+coord.GetVal(2,3)+coord.GetVal(2,0)+coord.GetVal(2,1));
		TPZFMatrix<REAL> jacobian;
		TPZFMatrix<REAL> axes;
		REAL detjac;
		TPZFMatrix<REAL> jacinv;
		Jacobian(coord,midle,jacobian,axes,detjac,jacinv);
		fNormalVec(16,0)=axes(0,0);
		fNormalVec(16,1)=axes(0,1);
		fNormalVec(16,2)=axes(0,2);
		fNormalVec(17,0)=axes(1,0);
		fNormalVec(17,1)=axes(1,1);
		fNormalVec(17,2)=axes(1,2);
		fVectorSide[count]=8;
		fVectorSide[count+1]=8;
		
		
		//normalização
		for(int k=0;k<16;k++)
		{
			REAL temp=0.;
			temp = sqrt( fNormalVec(k,0)*fNormalVec(k,0) + fNormalVec(k,1)*fNormalVec(k,1) + fNormalVec(k,2)*fNormalVec(k,2));
			fNormalVec(k,0) *=1./temp;	
			fNormalVec(k,1) *=1./temp;	
		}
		// produto normal == 1
		for(int kk=0;kk<4;kk++)
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
			sout << " fVectorSide" << fVectorSide;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif 
		
	}
	void TPZGeoQuad::VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result){
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
	
	void TPZGeoQuad::ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result){
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
	
	
	void TPZGeoQuad::Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){
		
		jacobian.Resize(2,2); axes.Resize(2,3); jacinv.Resize(2,2);
		TPZFNMatrix<4> phi(4,1);
		TPZFNMatrix<8> dphi(2,4);
		TPZFNMatrix<6> axest(3,2);
		Shape(param,phi,dphi);
		jacobian.Zero();
        TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
		
		int spacedim = coord.Rows();
        
        for (int j=0; j<spacedim; j++) {
            minx[j] = coord.GetVal(j,0);
            maxx[j] = coord.GetVal(j,0);
        }
		TPZFNMatrix<6,REAL> VecMatrix(3,2,0.);
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < spacedim; j++) {
                minx[j] = minx[j] < coord.GetVal(j,i) ? minx[j]:coord.GetVal(j,i);
                maxx[j] = maxx[j] > coord.GetVal(j,i) ? maxx[j]:coord.GetVal(j,i);
				VecMatrix(j,0) += coord.GetVal(j,i)*dphi(0,i);
				VecMatrix(j,1) += coord.GetVal(j,i)*dphi(1,i);
			}
		}
        REAL delx = 0.;
        for (int j=0; j<spacedim; j++) {
            delx = delx > (maxx[j]-minx[j]) ? delx : (maxx[j]-minx[j]);
        }
        VecMatrix *= 1./delx;
		VecMatrix.GramSchmidt(axest,jacobian);
		axest.Transpose(&axes);
		detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
		
        if(IsZero(detjac))
		{
#ifdef DEBUG
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

        jacobian *= delx;
        jacinv *= 1./delx;
        detjac *= (delx*delx);
        
	}
	
	void TPZGeoQuad::X(const TPZFMatrix<REAL> & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
		REAL spacephi[4],spacedphi[8];
		int i,j;
		TPZFMatrix<REAL> phi(4,1,spacephi,4);
		TPZFMatrix<REAL> dphi(2,4,spacedphi,8);
		int space = coord.Rows();
		Shape(loc,phi,dphi);
		result.Fill(0.);
		for(i=0;i<space;i++) {
			for(j=0;j<4;j++)
				//result[i] += phi(j,0)*NodePtr(j)->coord.GetVal(i);
				result[i] += phi(j,0)*coord.GetVal(i,j);
		}
	}
	
	TPZGeoEl *TPZGeoQuad::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
		if(side==8) {//8
			TPZManVector<long> nodes(4);
			int i;
			for (i=0;i<4;i++) {
				nodes[i] = orig->SideNodeIndex(side,i);
			}
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EQuadrilateral,nodes,bc,index);
			//    TPZGeoElQ2d *gel = new TPZGeoElQ2d(nodes,bc,*orig->Mesh());
			int iside;
			for (iside = 0; iside <8; iside++){
				TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::ContainedSideLocId(side,iside)));
			}
			TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if(side>-1 && side<4) {//side = 0,1,2,3
			TPZManVector<long> nodeindexes(1);
			//    TPZGeoElPoint *gel;
			nodeindexes[0] = orig->SideNodeIndex(side,0);
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//    gel = new TPZGeoElPoint(nodeindexes,bc,*(orig->Mesh()));
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if(side>3 && side<8) {
			TPZManVector<long> nodes(2);
			nodes[0] = orig->SideNodeIndex(side,0);
			nodes[1] = orig->SideNodeIndex(side,1);
			//    TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::ContainedSideLocId(side,0)));
			TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::ContainedSideLocId(side,1)));
			TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else PZError << "TPZGeoQuad::CreateBCCompEl has no bc.\n";
		return 0;
	}
	
	/** Creates a geometric element according to the type of the father element */
	TPZGeoEl *TPZGeoQuad::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										   TPZVec<long>& nodeindexes,
										   int matid,
										   long& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
};
