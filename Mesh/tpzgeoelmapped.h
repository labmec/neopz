/**
 * @file
 * @brief Contains declaration of TPZGeoElMapped class which implements a geometric element using its ancestral to compute ist jacobian.
 */
//
// C++ Interface: tpzgeoelmapped
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZGEOELMAPPED_H
#define TPZGEOELMAPPED_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzaxestools.h"
class TPZGeoMesh;

#include "pzlog.h"
#include <sstream>

#ifdef LOG4CXX
static LoggerPtr loggermapped(Logger::getLogger("pz.mesh.geoelmapped"));
#endif

/**
 @ingroup geometry
 @brief This class implements a geometric element which uses its ancestral to compute its jacobian. \ref geometry "Geometry"
 @author Philippe R. B. Devloo
 */
/**
 * Its main intent is to make the division of specially mapped elements easier: \n 
 * if the coarse grid map is consistent, then so will all refined meshes
 */
template<class TBase>
class TPZGeoElMapped : public TBase {
public:
	typedef typename TBase::Geo Geo;
	TPZGeoElMapped() : TBase(), fCornerCo(Geo::Dimension,Geo::NNodes,0.)
	{
	}
	TPZGeoElMapped(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
	TBase(id,nodeindexes,matind,mesh), fCornerCo(Geo::Dimension,Geo::NNodes,0.)
	{
	}
	TPZGeoElMapped(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh) :
	TBase(nodeindices,matind,mesh), fCornerCo(Geo::Dimension,Geo::NNodes,0.)
	{
	}
	TPZGeoElMapped(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index) :
	TBase(nodeindices,matind,mesh,index), fCornerCo(Geo::Dimension,Geo::NNodes,0.)
	{
	}
	
	~TPZGeoElMapped()
	{
	}
	
	virtual int ClassId() const;
    
    /** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid)
    {
        TBase::Write(buf,withclassid);
        fCornerCo.Write(buf,0);
    }
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context)
    {
        TBase::Read(buf,context);
        fCornerCo.Read(buf,0);
    }
	

	
	virtual bool IsLinearMapping() const
	{
		TPZGeoEl *father = (TBase::fFatherIndex == -1) ? 0 : TBase::Mesh()->ElementVec()[TBase::fFatherIndex];
		if(father) return father->IsLinearMapping();
		else return TBase::IsLinearMapping();
	}
	
	/** @brief Returns if is a TPZGeoElMapped< T > element
	 * It is necessary due to the lack of dynamic cast for
	 * these elements
	 */
	virtual bool IsGeoElMapped() const{
		return true;
	}
	
	/**
	 * @brief Creates a geometric element according to the type of the father element
	 */
	virtual TPZGeoEl *CreateGeoElement(MElementType type,
									   TPZVec<int>& nodeindexes,
									   int matid,
									   int& index);
	
	
	
	
	/** @brief Sets the father element index*/
	virtual void SetFather(int fatherindex)
	{
		TBase::SetFather(fatherindex);
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
		
		//  Pira 18 maio 2009: abaixo nao funciona para mapeamento nao-linear (nem bilinear)
		/*    for(in=0; in<nnodes; in++)
		 {
		 TPZTransform tr = Geo::SideToSideTransform(in,Geo::NSides-1);
		 TPZManVector<REAL,Geo::Dimension+1> ptin(0),ptout(Geo::Dimension,0.);
		 tr.Apply(ptin,ptout);
		 int nfs = father->NSides();
		 TPZGeoElSide thisside(this,Geo::NSides-1);
		 TPZGeoElSide ancestor(father,nfs-1);
		 TPZTransform trfat(Geo::Dimension);
		 TPZManVector<REAL,3> ptancestor(father->Dimension(),0.);
		 thisside.SideTransform3(ancestor,trfat);
		 trfat.Apply(ptout,ptancestor);
		 int id;
		 for(id=0; id<Geo::Dimension; id++)
		 {
		 fCornerCo(id,in) = ptancestor[id];
		 }
		 }
		 */
		
		//  Pira 18 maio 2009: nova implementação
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
			TPZTransform trzero(0);
			TPZTransform tr = this->BuildTransform2(in,father,trzero);
			TPZTransform trfather = fatherside.Element()->SideToSideTransform(fatherside.Side(),fatherside.Element()->NSides()-1);
			tr = trfather.Multiply(tr);
			TPZVec<REAL> zero(0);
			tr.Apply(zero, aux);
			//aux.Fill(0.);
			ptancestor.Fill(0.);
			father->ComputeXInverse(nodeX, aux,1e-10);
			int pointside = father->WhichSide(aux);
			TPZTransform project = father->Projection(pointside);
			project.Apply(aux,ptancestor);
			//      const int sideProjected = father->ProjectInParametricDomain(aux,ptancestor);
			
#ifdef DEBUG
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
				if(error > 1e-3){
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
				//	REAL size = father->CharacteristicSize();
				error = sqrt(error);
				if(error > 1e-8){
					std::stringstream sout;
					
					sout.precision(16);
					sout << "\nError at " << __PRETTY_FUNCTION__ << __LINE__ << "\n";
					sout << "this->Index = " << this->Index() << "\n";
					sout << "Node number " << in << std::endl;
					sout << "Node index " << this->NodeIndex(in) << std::endl;
					sout << "Father side " << this->FatherSide(in,0) << std::endl;
					sout << "aux:\n";
					for(int i = 0; i < aux.NElements(); i++) sout << aux[i] << "\t";
					sout << "\nptancestor:\n";
					for(int i = 0; i < ptancestor.NElements(); i++) sout << ptancestor[i] << "\t";
					sout << "\n";
					sout << "nodeX:\n";
					for(int i = 0; i < nodeX.NElements(); i++) sout << nodeX[i] << "\t";
					sout << "\nfatherX:\n";
					for(int i = 0; i < fatherX.NElements(); i++) sout << fatherX[i] << "\t";
					sout << "\nfatherXaux:\n";
					for(int i = 0; i < fatherXaux.NElements(); i++) sout << fatherXaux[i] << "\t";
					std::cout << sout.str();
#ifdef LOG4CXX
					father->Mesh()->Print(sout);
					LOGPZ_ERROR(loggermapped,sout.str())
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
	
	
	/** @brief Returns the Jacobian matrix at the point (from son to father)*/
	virtual void Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv)
	{
		/**
		 * @brief Creating Variables
		 */
		TPZGeoEl *father = TBase::Father();
		if(!father) 
		{
			TBase::Jacobian(coordinate,jac,axes,detjac,jacinv);
			return;
		}
		TPZGeoEl *nextfather = 0;
		if(father) nextfather = father->Father();
		while(nextfather)
		{
			father = nextfather;
			nextfather = father->Father();
		}
		const int dim = Geo::Dimension;
		TPZManVector<REAL,3> ksibar(father->Dimension());
		TPZFNMatrix<dim*dim+1> jaclocal(dim,dim,0.),jacinvlocal(dim,dim,0.),jacfather(dim,dim,0.), jacinvfather(dim,dim,0.);
		TPZFNMatrix<9> axeslocal(3,3,0.), axesfather(3,3,0.);
		REAL detjaclocal, detjacfather;
		
		/**
		 * @brief Processing Variables (isolated)
		 */
		Geo::Jacobian(fCornerCo,coordinate,jaclocal,axeslocal,detjaclocal,jacinvlocal);    
		Geo::X(fCornerCo,coordinate,ksibar);
		father->Jacobian(ksibar,jacfather,axesfather,detjacfather,jacinvfather);
		
		/**
		 * @brief Combining Variables
		 */
		TPZFNMatrix<9> aux(dim,dim);
		
		//jacinv
		axeslocal.Resize(dim,dim); //reducing axes local to its correct dimension in this context
		axeslocal.Multiply(jacinvfather,aux);
		jacinvlocal.Multiply(aux,jacinv);
		
		//jac
		axeslocal.Transpose();
		axeslocal.Multiply(jaclocal,aux);
		jacfather.Multiply(aux,jac);
		
		//detjac
		detjac = detjaclocal*detjacfather;
		
		//axes
		axes = axesfather;
		
	}
	
	/** @brief Returns the coordinate in real space of the point coordinate in the master element space*/
	virtual void X(TPZVec<REAL> &ksi,TPZVec<REAL> &result)
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
			
			TPZManVector<REAL,3> ksibar(this->Dimension());
			KsiBar(ksi,ksibar);
			father->X(ksibar,result);
		}
	}
	
	virtual void Print(std::ostream & out = std::cout)
	{
        TBase::Print(out);
		
        fCornerCo.Print("fCornerCo Print():",out);
	}
	
	/**Avaliate the Jacobian 2D (2x2 size) by Expected Convergence Order*/
	//   virtual void JacobianConv(TPZVec< REAL > QsiEta,  TPZVec< REAL > &XYZ)
	//   {std::cout << "***USING THE JACOBIAN FOR 2D ELEMENTS METHOD***\n\n";
	//           TPZFMatrix<REAL> jacobian(2,2);
	//           TPZFMatrix<REAL> Axes(3,3);
	//           REAL detJacobian;
	//           TPZFMatrix<REAL> InvJac(2,2);
	//           TPZVec< REAL > QsiEtaIni (2,1);
	//           QsiEtaIni[0] = QsiEta[0];
	//           QsiEtaIni[1] = QsiEta[1];
	//           const double deltaQsi = 0.1;
	//           const double deltaEta = 0.1;
	//           double alpha;
	// 
	//           std::cout << "\ninitial Qsi = " << QsiEtaIni[0] << " | initial Eta = " << QsiEtaIni[1] << "\n";
	//           std::cout << "deltaQsi = const = " << deltaQsi << " | deltaEta = const = " << deltaEta << "\n\n";
	// 
	//           TPZVec< REAL > XYZaprox(3);
	//           TPZFMatrix<REAL> error(11,1,0.);
	//           int edge;
	//           double dX, dY, dZ;
	//           for(int i = 0; i <= 10; i++)
	//           {
	//                 alpha = i/10.;
	//                 std::cout << "\nalpha = " << alpha << std::endl;
	//                 X(QsiEtaIni,XYZ);//for aproximate compute
	//                 std::cout << "f(x)                     : x = " << XYZ[0] << " | y = " << XYZ[1] << " | z = " << XYZ[2] << "\n";
	//                 Jacobian(QsiEtaIni,jacobian,Axes,detJacobian,InvJac);
	// 
	//                 dX = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,0) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,0);
	//                 XYZaprox[0] = XYZ[0] + dX;
	// 
	//                 dY = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,1) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,1);
	//                 XYZaprox[1] = XYZ[1] + dY;
	// 
	//                 dZ = alpha*( jacobian.GetVal(0,0)*deltaQsi + jacobian.GetVal(0,1)*deltaEta )*Axes(0,2) + alpha*( jacobian.GetVal(1,0)*deltaQsi + jacobian.GetVal(1,1)*deltaEta )*Axes(1,2);
	//                 XYZaprox[2] = XYZ[2] + dZ;
	// 
	//                 QsiEta[0] = QsiEtaIni[0] + alpha*deltaQsi;
	//                 QsiEta[1] = QsiEtaIni[1] + alpha*deltaEta;
	//                 edge = i;
	//                 if(QsiEta[1] > 1.00001 - QsiEta[0]) break;
	//                 X(QsiEta,XYZ);//for real compute
	// 
	//                 std::cout << "alpha.(axes.J).dx        : x = " << dX << " | y = " << dY << " | z = " << dZ << "\n";
	//                 std::cout << "f(x) + alpha.(axes.J).dx : x = " << XYZaprox[0] << " | y = " << XYZaprox[1] << " | z = " << XYZaprox[2] << "\n";
	//                 std::cout << "f(x + alpha.dx)          : x = " << XYZ[0] << " | y = " << XYZ[1] << " | z = " << XYZ[2] << "\n";
	// 
	//                 XYZ[0] -= XYZaprox[0];
	//                 XYZ[1] -= XYZaprox[1];
	//                 XYZ[2] -= XYZaprox[2];
	// 
	//                 double XDiffNorm = sqrt(XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2]);
	//                 error(int(i),0) = XDiffNorm;
	//           }
	// 
	//           std::cout << "\n|| [f(x + alpha.dx)] - [f(x) + alpha.(axes.J).dx] ||:\n\n"; error.Print();
	// 
	//           std::cout << "\nConvergence Order:\n";
	//           for(int j = 2; j < edge; j++)
	//           {
	//                 std::cout << ( log(error(1,0)) - log(error(j,0)) )/( log(0.1)-log(j/10.) ) << "\n";
	//           }
	//           if(edge != 10) std::cout << "The direction track has touch the edge of the element.\nThe method has stopped!\n";
	//   }
	
	/** @brief Avaliate the Jacobian 3D (3x3 size) by Expected Convergence Order*/
	virtual void JacobianConv(TPZVec< REAL > QsiEta,  TPZVec< REAL > &XYZ)
	{std::cout << "\n***USING THE JACOBIAN FOR 3D ELEMENTS METHOD***\n";
		TPZFMatrix<REAL> jacobian(3,3);
		TPZFMatrix<REAL> Axes(3,3);
		REAL detJacobian;
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
			Jacobian(QsiEtaIni,jacobian,Axes,detJacobian,InvJac);
			
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
	
    // Compute the map of the point ksi to the ancestor ksibar and the gradient of the ancestor ksibar with respect to ksi
    void KsiBar(TPZVec<REAL> &ksi, TPZVec<REAL> &ksibar, TPZFMatrix<REAL> &jac)
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
				ksibar[id] += phi(in,0)*fCornerCo(id,in);
				for(jd=0; jd<dim; jd++)
				{
					jac(id,jd) += dphi(jd,in)*fCornerCo(id,in);
				}
			}
		}
    }
	
    // Compute the map of the point ksi to the ancestor ksibar
    void KsiBar(TPZVec<REAL> &ksi, TPZVec<REAL> &ksibar)
    {
		const int dim = Geo::Dimension;
		TPZFNMatrix<Geo::NNodes> phi(Geo::NNodes,1,0.);
		TPZFNMatrix<dim*dim+1> jac(dim,dim,0.);
		TPZFNMatrix<dim*Geo::NNodes+1> dphi(dim,Geo::NNodes,0.);
		Geo::Shape(ksi,phi,dphi);
		ksibar.Fill(0.);
		int in,id;
		for(in=0; in<Geo::NNodes; in++)
		{
			for(id=0; id<dim; id++)
			{
				ksibar[id] += phi(in,0)*fCornerCo(id,in);
			}
		}
    }
	
	virtual TPZGeoEl *CreateBCGeoEl(int side, int bc){
		int ns = this->NSideNodes(side);
		TPZManVector<int> nodeindices(ns);
		int in;
		for(in=0; in<ns; in++)
		{
			nodeindices[in] = this->SideNodeIndex(side,in);
		}
		int index;
		
		TPZGeoMesh *mesh = this->Mesh();
		MElementType type = this->Type(side);
		
		TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
		TPZGeoElSide me(this,side);
		TPZGeoElSide newelside(newel,newel->NSides()-1);
		
		newelside.InsertConnectivity(me);
		newel->Initialize();
		
		return newel;
	}
	
};

/** @brief Creates geometric element of the specified type */
TPZGeoEl *CreateGeoElementMapped(TPZGeoMesh &mesh,
								 MElementType type,
								 TPZVec<int>& nodeindexes,
								 int matid,
								 int& index);

#endif
