/**
 * @file
 * @brief Contains declaration of TPZGeoElMapped class which implements a geometric element using its ancestral to compute ist jacobian.
 */

#ifndef TPZGEOELMAPPED_H
#define TPZGEOELMAPPED_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzaxestools.h"
#include "pzgmesh.h"

#include "pzlog.h"
#include <sstream>

/**
 * @ingroup geometry
 * @brief This class implements a geometric element which uses its ancestral to compute its jacobian. \ref geometry "Geometry"
 * @author Philippe R. B. Devloo
 */
/**
 * Its main intent is to make the division of specially mapped elements easier: \n 
 * if the coarse grid map is consistent, then so will all refined meshes
 */
template<class TBase>
class TPZGeoElMapped : public TBase {
public:
	typedef typename TBase::Geo Geo;
	TPZGeoElMapped() : TPZRegisterClassId(&TPZGeoElMapped::ClassId),
    TBase(), fCornerCo(Geo::Dimension,Geo::NNodes,0.)
	{
	}
	TPZGeoElMapped(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
	TPZRegisterClassId(&TPZGeoElMapped::ClassId),
    TBase(id,nodeindexes,matind,mesh), fCornerCo(Geo::Dimension,Geo::NNodes,0.)
	{
	}
	TPZGeoElMapped(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh) :
	TPZRegisterClassId(&TPZGeoElMapped::ClassId),
    TBase(nodeindices,matind,mesh), fCornerCo(Geo::Dimension,Geo::NNodes,0.)
	{
	}
	TPZGeoElMapped(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh,int64_t &index) :
	TPZRegisterClassId(&TPZGeoElMapped::ClassId),
    TBase(nodeindices,matind,mesh,index), fCornerCo(Geo::Dimension,Geo::NNodes,0.)
	{
	}
    
    TPZGeoElMapped(TPZGeoMesh &destmesh, const TPZGeoElMapped<TBase> &copy);
    
    TPZGeoElMapped(TPZGeoMesh &destmesh, const TPZGeoElMapped<TBase> &copy, std::map<int64_t,int64_t> &gl2lcNdIdx,
                   std::map<int64_t,int64_t> &gl2lcElIdx);
	
	~TPZGeoElMapped()
	{
	}
	
	public:
int ClassId() const override;

    
    virtual TPZGeoEl * Clone(TPZGeoMesh &DestMesh) const override;
    
	/** @} */
	
	virtual TPZGeoEl * ClonePatchEl(TPZGeoMesh &DestMesh,
									std::map<int64_t,int64_t> &gl2lcNdIdx,
									std::map<int64_t,int64_t> &gl2lcElIdx) const override;
	

    
    /** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override{
        TBase::Write(buf,withclassid);
        fCornerCo.Write(buf,0);
    }
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream& buf, void* context) override
    {
        TBase::Read(buf,context);
        fCornerCo.Read(buf,0);
    }
	
    virtual bool IsLinearMapping(int side) const override;
	
    /*
	virtual bool IsLinearMapping() const
	{
		TPZGeoEl *father = (TBase::fFatherIndex == -1) ? 0 : TBase::Mesh()->ElementVec()[TBase::fFatherIndex];
		if(father) return father->IsLinearMapping();
		else return TBase::IsLinearMapping();
	}
     */
	
	/** @brief Returns if is a TPZGeoElMapped< T > element */
	/** It is necessary due to the lack of dynamic cast for these elements */
	virtual bool IsGeoElMapped() const override {
		return true;
	}
    
    
	
	/** @brief Creates a geometric element according to the type of the father element */
	virtual TPZGeoEl *CreateGeoElement(MElementType type,
									   TPZVec<int64_t>& nodeindexes,
									   int matid,
									   int64_t& index) override;

	/** @brief Sets the father element index*/
	virtual void SetFatherIndex(int64_t fatherindex) override
	{
		TBase::SetFatherIndex(fatherindex);
		TPZGeoEl *father = TBase::Father();
		if(!father) return;
		TPZGeoEl *nextfather = 0;
		if(father) nextfather = father->Father();
		while(nextfather)
		{
			father = nextfather;
			nextfather = father->Father();
		}
		int in, nnodes = Geo::NNodes;
		TPZManVector<REAL,3> nodeX(3);
		TPZManVector<REAL,3> ptancestor(Geo::Dimension), aux(Geo::Dimension);
		REAL Tol;
		ZeroTolerance(Tol);
		if(Tol < 1.e-10) Tol = 1.e-10;
        
		//  Pira 18 maio 2009: nova implementação
        
        // @omar:: decreased for non linear mappings
        REAL epsilon = 0.1;
		for(in=0; in<nnodes; in++)
		{
			for(int id = 0; id < 3; id++){
				nodeX[id] = this->NodePtr(in)->Coord(id);
			}
			
			TPZGeoElSide gels(this,in);
			TPZGeoElSide nextfatherside = gels.Father2();
			TPZGeoElSide fatherside;
			while(nextfatherside.Element())
			{
				fatherside = nextfatherside;
				nextfatherside = nextfatherside.Father2();			
			}
			TPZTransform<> trzero(0);
			TPZTransform<> tr = this->BuildTransform2(in,father,trzero);
			TPZTransform<> trfather = fatherside.Element()->SideToSideTransform(fatherside.Side(),fatherside.Element()->NSides()-1);
			tr = trfather.Multiply(tr);
			TPZVec<REAL> zero(0);
			tr.Apply(zero, aux);
			//aux.Fill(0.);
			ptancestor.Fill(0.);
			father->ComputeXInverse(nodeX, aux, epsilon);
			int pointside = father->WhichSide(aux);
			TPZTransform<> project = father->Projection(pointside);
			project.Apply(aux,ptancestor);
			
#ifdef PZDEBUG
			{
				double ksidiff = 0.;
				for(int i = 0; i < ptancestor.NElements(); i++){
					ksidiff += (aux[i]-ptancestor[i])*(aux[i]-ptancestor[i]);
				}//i
				ksidiff = sqrt(ksidiff);
				if(ksidiff > 1e-8){
					std::cout.precision(12);
					std::cout << "\nError at " << __PRETTY_FUNCTION__ << __LINE__ << "\n";
					std::cout << "aux:\n";
					for(int i = 0; i < aux.NElements(); i++) std::cout << aux[i] << "\t";
					std::cout << "\nptancestor:\n";
					for(int i = 0; i < ptancestor.NElements(); i++) std::cout << ptancestor[i] << "\t";
					std::cout << "\n";
					DebugStop();
				}
				
				TPZManVector<REAL,3> fatherX(3);
				father->X(ptancestor,fatherX);
				double error = 0.;
				double diff;
				for(int i = 0; i < 3; i++){
					diff = fatherX[i]-nodeX[i];
					error += diff*diff;
				}
				error = sqrt(error);
				if(error > epsilon){
					std::cout << "\nError at " << __PRETTY_FUNCTION__ << __LINE__ << "\n";
					std::cout << "this->Index = " << this->Index() << "\n";
					std::cout << "aux:\n";
					for(int i = 0; i < aux.NElements(); i++) std::cout << aux[i] << "\t";
					std::cout << "\nptancestor:\n";
					for(int i = 0; i < ptancestor.NElements(); i++) std::cout << ptancestor[i] << "\t";
					std::cout << "\n";
					std::cout << "nodeX:\n";
					for(int i = 0; i < nodeX.NElements(); i++) std::cout << nodeX[i] << "\t";
					std::cout << "\nfatherX:\n";
					for(int i = 0; i < fatherX.NElements(); i++) std::cout << fatherX[i] << "\t";
					std::cout << "\n";
					DebugStop();
				}
				
				TPZManVector<REAL,3> fatherXaux(3);
				father->X(aux,fatherXaux);
				error = 0.;
				for(int i = 0; i < 3; i++){
					diff = fatherX[i]-fatherXaux[i];
					error += diff*diff;
				}

				error = sqrt(error);
				if(error > 1e-8){
#ifdef PZ_LOG
                    TPZLogger loggermapped("pz.mesh.geoelmapped");
                    if(loggermapped.isErrorEnabled()){
                        std::stringstream sout;
					
                        sout.precision(16);
                        sout << "\nError at " << __PRETTY_FUNCTION__ << __LINE__ << "\n";
                        sout << "this->Index = " << this->Index() << "\n";
                        sout << "Node number " << in << std::endl;
                        sout << "Node index " << this->NodeIndex(in)
                             << std::endl;
                        sout << "Father side " << this->FatherSide(in, 0)
                             << std::endl;
                        sout << "aux:\n";
                        for (int i = 0; i < aux.NElements(); i++)
                          sout << aux[i] << "\t";
                        sout << "\nptancestor:\n";
                        for (int i = 0; i < ptancestor.NElements(); i++)
                          sout << ptancestor[i] << "\t";
                        sout << "\n";
                        sout << "nodeX:\n";
                        for (int i = 0; i < nodeX.NElements(); i++)
                          sout << nodeX[i] << "\t";
                        sout << "\nfatherX:\n";
                        for (int i = 0; i < fatherX.NElements(); i++)
                          sout << fatherX[i] << "\t";
                        sout << "\nfatherXaux:\n";
                        for (int i = 0; i < fatherXaux.NElements(); i++)
                          sout << fatherXaux[i] << "\t";
                        std::cout << sout.str();
                        LOGPZ_ERROR(loggermapped, sout.str())
                    }
#endif
					DebugStop();
				}
				
			}
#endif
			for(int id=0; id<Geo::Dimension; id++)
			{
				fCornerCo(id,in) = ptancestor[id];
			}
			
		}//for
		
	}//method
	
    /** @brief Return the Jacobian matrix at the point*/
    virtual void GradXFad(TPZVec<Fad<REAL> > &qsi, TPZFMatrix<Fad<REAL> > &gradx) const  {
        DebugStop();
    }
//    /** @brief Return the Jacobian matrix at the point*/
//    virtual void GradX(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &gradx) const {
//
//            /// Creating Variables
//            TPZGeoEl *father = TBase::Father();
//            if(!father)
//            {
//                TBase::GradX(qsi,gradx);
//                return;
//            }
//        
//            TPZGeoEl *nextfather = 0;
//            if(father) nextfather = father->Father();
//            while(nextfather)
//            {
//                father = nextfather;
//                nextfather = father->Father();
//            }
//        
//            const int dim = Geo::Dimension;
//            const int father_dim = father->Dimension();
//            TPZManVector<REAL,3> ksibar(father_dim,0.0);
//        
//            TPZFNMatrix<9> gradxlocal(father_dim,dim);
//            Geo::X(fCornerCo,qsi,ksibar);
//            Geo::GradX(fCornerCo,qsi,gradxlocal);
//        
//            TPZFNMatrix<9> gradxfather;
//            father->GradX(ksibar, gradxfather);
//
//            // @brief Combining Variables
//            gradxfather.Multiply(gradxlocal, gradx);
//    }
    
    /** @brief Return the Jacobian matrix at the point*/
    void GradX(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &gradx) const  override {
        return TGradX(qsi, gradx);
    }
    void GradX(TPZVec<Fad<REAL>> &qsi, TPZFMatrix<Fad<REAL>> &gradx) const  override {
        return TGradX(qsi, gradx);
    }

//	/** @brief Returns the Jacobian matrix at the point (from son to father)*/
//	virtual void Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
//	{
//		/// Creating Variables
//		TPZGeoEl *father = TBase::Father();
//		if(!father) 
//		{
//			TBase::Jacobian(coordinate,jac,axes,detjac,jacinv);
//			return;
//		}
//		TPZGeoEl *nextfather = 0;
//		if(father) nextfather = father->Father();
//		while(nextfather)
//		{
//			father = nextfather;
//			nextfather = father->Father();
//		}
//		const int dim = Geo::Dimension;
//		TPZManVector<REAL,3> ksibar(father->Dimension());
//		TPZFNMatrix<dim*dim+1> jaclocal(dim,dim,0.),jacinvlocal(dim,dim,0.),jacfather(dim,dim,0.), jacinvfather(dim,dim,0.);
//		TPZFNMatrix<9> axeslocal(3,3,0.), axesfather(3,3,0.);
//		REAL detjaclocal, detjacfather;
//		
//        TPZFNMatrix<9> gradxlocal;
//        Geo::GradX(fCornerCo,coordinate,gradxlocal);
//		/// Processing Variables (isolated)
//		Geo::Jacobian(fCornerCo,coordinate,jaclocal,axeslocal,detjaclocal,jacinvlocal);    
//		Geo::X(fCornerCo,coordinate,ksibar);
//		father->Jacobian(ksibar,jacfather,axesfather,detjacfather,jacinvfather);
//		
//		/// @brief Combining Variables
//		TPZFNMatrix<9> aux(dim,dim);
//		
//		//jacinv
//		axeslocal.Resize(dim,dim); //reducing axes local to its correct dimension in this context
//		axeslocal.Multiply(jacinvfather,aux);
//		jacinvlocal.Multiply(aux,jacinv);
//		
//		//jac
//		axeslocal.Transpose();
//		axeslocal.Multiply(jaclocal,aux);
//		jacfather.Multiply(aux,jac);
//		
//		//detjac
//		detjac = detjaclocal*detjacfather;
//		
//		//axes
//		axes = axesfather;
//	}
	
	/** @brief Returns the coordinate in real space of the point coordinate in the master element space*/
    void X(TPZVec<REAL> &ksi,TPZVec<REAL> &result) const override {
	        return TX(ksi,result);
	}
    void X(TPZVec<Fad<REAL>> &ksi,TPZVec<Fad<REAL>> &result) const override {
        return TX(ksi,result);
    }

	virtual void Print(std::ostream & out = std::cout) override
	{
        TBase::Print(out);
		
        fCornerCo.Print("fCornerCo Print():",out);
	}
	
	/** @brief Avaliate the Jacobian 3D (3x3 size) by Expected Convergence Order*/
	virtual void JacobianConv(TPZVec< REAL > QsiEta,  TPZVec< REAL > &XYZ)
	{std::cout << "\n***USING THE JACOBIAN FOR 3D ELEMENTS METHOD***\n";
		TPZFMatrix<REAL> jacobian(3,3);
		TPZFMatrix<REAL> Axes(3,3);
		TPZFMatrix<REAL> InvJac(3,3);
		TPZVec< REAL > QsiEtaIni (3,1);
		QsiEtaIni[0] = QsiEta[0];
		QsiEtaIni[1] = QsiEta[1];
		QsiEtaIni[2] = QsiEta[2];
		const double deltaQsi  = 0.1;
		const double deltaEta  = 0.1;
		const double deltaZeta = 0.1;
		double alpha;
		
		std::cout << "\ninitial Qsi = " << QsiEtaIni[0] << " | initial Eta = " << QsiEtaIni[1] << " | initial Zeta = " << QsiEtaIni[2] <<"\n";
		std::cout << "deltaQsi = const = " << deltaQsi << " | deltaEta = const = " << deltaEta << " | deltaZeta = const = " << deltaZeta <<"\n\n";
		
		TPZVec< REAL > XYZaprox(3);
		TPZFMatrix<REAL> error(11,1,0.);
		int edge;
		double dX, dY, dZ;
		for(int i = 0; i <= 10; i++)
		{
			alpha = i/10.;
			X(QsiEtaIni,XYZ);//for aproximate compute
			//Jacobian(QsiEtaIni,jacobian,Axes,detJacobian,InvJac);
			
			dX = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta + jacobian.GetVal(0,2)*deltaZeta)*Axes(0,0) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta + jacobian.GetVal(1,2)*deltaZeta)*Axes(1,0) + alpha*( jacobian.GetVal(2,0)*deltaQsi + jacobian.GetVal(2,1)*deltaEta + jacobian.GetVal(2,2)*deltaZeta)*Axes(2,0);
			XYZaprox[0] = XYZ[0] + dX;
			
			dY = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta + jacobian.GetVal(0,2)*deltaZeta)*Axes(0,1) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta + jacobian.GetVal(1,2)*deltaZeta)*Axes(1,1) + alpha*( jacobian.GetVal(2,0)*deltaQsi + jacobian.GetVal(2,1)*deltaEta + jacobian.GetVal(2,2)*deltaZeta)*Axes(2,1);
			XYZaprox[1] = XYZ[1] + dY;
			
			dZ = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta + jacobian.GetVal(0,2)*deltaZeta)*Axes(0,2) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta + jacobian.GetVal(1,2)*deltaZeta)*Axes(1,2) + alpha*( jacobian.GetVal(2,0)*deltaQsi + jacobian.GetVal(2,1)*deltaEta + jacobian.GetVal(2,2)*deltaZeta)*Axes(2,2);
			XYZaprox[2] = XYZ[2] + dZ;
			
			QsiEta[0] = QsiEtaIni[0] + alpha*deltaQsi;
			QsiEta[1] = QsiEtaIni[1] + alpha*deltaEta;
			QsiEta[2] = QsiEtaIni[2] + alpha*deltaZeta;
			edge = i + 1;
			if((QsiEta[1]>1.00001-QsiEta[0]) || QsiEta[2]>1.00001-QsiEta[0]-QsiEta[1]) break;
			X(QsiEta,XYZ);//for real compute
			
			XYZ[0] -= XYZaprox[0];
			XYZ[1] -= XYZaprox[1];
			XYZ[2] -= XYZaprox[2];
			
			double XDiffNorm = sqrt(XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2]);
			error(int(i),0) = XDiffNorm;
		}
		
		std::cout << "ERROR Vector:\n";
		for(int j = 2; j < edge; j++)
		{
			std::cout << error(j,0) << "\n";
		}
		
		std::cout << "\nConvergence Order:\n";
		for(int j = 2; j < edge; j++)
		{
			std::cout << ( log(error(1,0)) - log(error(j,0)) )/( log(0.1)-log(j/10.) ) << "\n";
		}
		if(edge != 11) std::cout << "The direction track has touch the edge of the element.\nThe method has stopped!\n";
	}
	
private:
	
#ifdef WIN32
	TPZFNMatrix<24/*Geo::Dimension*Geo::NNodes*/> fCornerCo;
#else
	TPZFNMatrix<Geo::Dimension*Geo::NNodes> fCornerCo;
#endif
	
    /** @brief Compute the map of the point ksi to the ancestor ksibar and the gradient of the ancestor ksibar with respect to ksi */
    void KsiBar(TPZVec<REAL> &ksi, TPZVec<REAL> &ksibar, TPZFMatrix<REAL> &jac) const
    {
		const int dim = Geo::Dimension;
		TPZFNMatrix<Geo::NNodes> phi(Geo::NNodes,1,0.);
		TPZFNMatrix<dim*Geo::NNodes> dphi(dim,Geo::NNodes,0.);
		Geo::Shape(ksi,phi,dphi);
		jac.Redim(dim,dim);
		ksibar.Fill(0.);
		int in,id,jd;
		for(in=0; in<Geo::NNodes; in++)
		{
			for(id=0; id<dim; id++)
			{
				ksibar[id] += phi(in,0)*fCornerCo.GetVal(id,in);
				for(jd=0; jd<dim; jd++)
				{
					jac(id,jd) += dphi(jd,in)*fCornerCo.GetVal(id,in);
				}
			}
		}
    }
	
    /** @brief Compute the map of the point ksi to the ancestor ksibar */
    void KsiBar(TPZVec<REAL> &ksi, TPZVec<REAL> &ksibar) const{
        return TKsiBar(ksi,ksibar);
    }

    void KsiBar(TPZVec<Fad<REAL>> &ksi, TPZVec<Fad<REAL>> &ksibar) const{
        return TKsiBar(ksi,ksibar);
    }

	virtual TPZGeoEl *CreateBCGeoEl(int side, int bc) override {
		int ns = this->NSideNodes(side);
		TPZManVector<int64_t> nodeindices(ns);
		int in;
		for(in=0; in<ns; in++)
		{
			nodeindices[in] = this->SideNodeIndex(side,in);
		}
		int64_t index;
		
		TPZGeoMesh *mesh = this->Mesh();
		MElementType type = this->Type(side);
		
		TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
		TPZGeoElSide me(this,side);
		TPZGeoElSide newelside(newel,newel->NSides()-1);
		
		newelside.InsertConnectivity(me);
		newel->Initialize();
		
		return newel;
	}

protected:
    template<class T>
    void TGradX(TPZVec<T> &qsi, TPZFMatrix<T> &gradx) const {

        /// Creating Variables
        TPZGeoEl *father = TBase::Father();
        if(!father)
        {
            TBase::GradX(qsi,gradx);
            return;
        }
        TPZGeoEl *nextfather = father;
        while(nextfather)
        {
            father = nextfather;
            nextfather = father->Father();
        }

        TPZManVector<T,3> ksibar(father->Dimension());
        TPZFNMatrix<9,T> gradxlocal;
        Geo::GradX(fCornerCo,qsi,gradxlocal);
        Geo::X(fCornerCo,qsi,ksibar);
        TPZFNMatrix<9,T> gradxfather;
        father->GradX(ksibar, gradxfather);

        /// @brief Combining Variables
        gradxfather.Multiply(gradxlocal, gradx);
#ifdef PZ_LOG
        TPZLogger loggermapped("pz.mesh.geoelmapped");
        if(loggermapped.isDebugEnabled())
        {
            std::stringstream sout;
            gradxfather.Print("gradx father",sout);
            gradxlocal.Print("gradx local",sout);
            gradx.Print("gradx",sout);
            LOGPZ_DEBUG(loggermapped, sout.str())
        }
#endif
    }

    template<class T>
    void TX(TPZVec<T> &ksi,TPZVec<T> &result) const
    {
        TPZGeoEl *father = TBase::Father();

        if(!father)
        {
            return TBase::X(ksi,result);
        }

        else
        {
            TPZGeoEl *nextfather = 0;
            if(father) nextfather = father->Father();
            while(nextfather)
            {
                father = nextfather;
                nextfather = father->Father();
            }

            TPZManVector<T,3> ksibar(this->Dimension());
            KsiBar(ksi,ksibar);
            father->X(ksibar,result);
        }
    }

    /** @brief Compute the map of the point ksi to the ancestor ksibar */
    template<class T>
    void TKsiBar(TPZVec<T> &ksi, TPZVec<T> &ksibar) const
    {
        const int dim = Geo::Dimension;
        TPZFNMatrix<Geo::NNodes,T> phi(Geo::NNodes,1,0.);
        TPZFNMatrix<dim*dim+1,T> jac(dim,dim,0.);
        TPZFNMatrix<dim*Geo::NNodes+1,T> dphi(dim,Geo::NNodes,0.);
        Geo::TShape(ksi,phi,dphi);
        ksibar.Fill(0.);
        int in,id;
        for(in=0; in<Geo::NNodes; in++)
        {
            for(id=0; id<dim; id++)
            {
                ksibar[id] += phi(in,0)*fCornerCo.GetVal(id,in);
            }
        }
    }
};

template<class TBase>
inline bool TPZGeoElMapped<TBase>::IsLinearMapping(int side) const
{
    TPZGeoElSide fatherside = this->Father2(side);
    if (fatherside.Element()) {
        return fatherside.IsLinearMapping();
    }
    else
    {
        return TBase::IsLinearMapping(side);
    }
}

template<class TBase>
int TPZGeoElMapped<TBase>::ClassId() const{
    return Hash("TPZGeoElMapped") ^ TBase::ClassId() << 1;
}
#endif
