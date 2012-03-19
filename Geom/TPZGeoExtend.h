/**
 * @file
 * @brief Contains the GPr class.
 */

#ifndef TPZGEOPRISMEXTENDH
#define TPZGEOPRISMEXTENDH

#include "pznoderep.h"

#include "pzvec.h"
#include "pzeltype.h"
#include "pzfmatrix.h"

#include <string>
#include "PrismExtend.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr loggernoderep2(Logger::getLogger("pz.geom.extend"));
#endif

template <class TVAR>
class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
	
	/** 
	 * @ingroup geometry
	 * @brief Geometric derivated class (template). \ref geometry "Geometry"
	 */
	template<class TFather, class Topology>
	class GPr : public TFather {
		
	public:
		
		typedef typename TFather::TMem FatMem;
		typedef typename TFather::Top FatTop;
		
		typedef Topology Top;
		
		enum {NNodes = 2*TFather::NNodes};
		
		int fNodeIndexes[TFather::NNodes];
		/** @brief Constructor with list of nodes */
		GPr(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TFather(nodeindexes, mesh)
		{
			int i;
			for(i=0; i<TFather::NNodes; i++) fNodeIndexes[i] = nodeindexes[TFather::NNodes+i];
		}
		
		/** @brief Empty constructor */
		GPr() : TFather()
		{
			int i;
			for(i=0; i<TFather::NNodes; i++) fNodeIndexes[i] = -1;
		}
		
		/** @brief Constructor with node map */
		GPr(const GPr<TFather,Topology> &cp,
			std::map<int,int> & gl2lcNdMap) : TFather(cp,gl2lcNdMap)
		{
			int i;
			for(i = 0; i < TFather::NNodes; i++)
			{
#ifdef DEBUG
				if (gl2lcNdMap.find(cp.fNodeIndexes[i+TFather::NNodes]) == gl2lcNdMap.end())
				{
					std::stringstream sout;
					sout << "ERROR in - " << __PRETTY_FUNCTION__
					<< " trying to clone a node " << i << " index " << cp.fNodeIndexes[i]
					<< " wich is not mapped";
					LOGPZ_ERROR(loggernoderep2,sout.str().c_str());
					fNodeIndexes[i] = -1;
					continue;
				}
#endif
				fNodeIndexes[i] = gl2lcNdMap [ cp.fNodeIndexes[i+TFather::NNodes] ];
			}
		}
		
		/** @brief Copy constructor */
		GPr(const GPr<TFather, Topology> &cp) : TFather(cp)
		{
			int i;
			for(i=0; i<TFather::NNodes; i++) fNodeIndexes[i] = cp.fNodeIndexes[i];
		}
		
		/** @brief Copy constructor */
		GPr(const GPr<TFather,Topology> &cp, TPZGeoMesh &) : TFather(cp)
		{
			int i;
			for(i=0; i<TFather::NNodes; i++) fNodeIndexes[i] = cp.fNodeIndexes[i];
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Pr:"+TFather::TypeName();} 
		
		static void X(TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);
		
		static void Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		static void Jacobian(TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
							 TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
		
		static void Jacobian(TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian);
		
		
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
		
		static void Diagnostic(TPZFMatrix<REAL> &coord);
		
		
	};
	
	template<class TFather, class Topology>
	inline void GPr<TFather, Topology>::Jacobian(TPZFMatrix<REAL> &coord,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian) 
	{
		
		int spacedim = coord.Rows();
		TPZFNMatrix<3*TFather::NNodes> coordlower(spacedim,TFather::NNodes),coordupper(spacedim,TFather::NNodes);
		TPZFNMatrix<100> jacobianlower,jacobianupper;
		TFather::Jacobian(coordupper,param,jacobianlower);
		TFather::Jacobian(coordupper,param,jacobianupper);
		REAL ksi = param[Top::Dimension-1];
		jacobian.Resize(spacedim,Top::Dimension);
		TPZManVector<REAL> xlower(spacedim),xupper(spacedim);
		TPZManVector<REAL> resultlower,resultupper;
		int i,j;
		for(i=0; i<spacedim; i++) for(j=0; j<TFather::NNodes; j++)
		{
			coordlower(i,j) = coord(i,j);
			coordupper(i,j) = coord(i,j+TFather::NNodes);
		}
		TFather::X(coordlower,param,xlower);
		TFather::X(coordupper,param,xupper);
		for(i=0; i<spacedim; i++)
		{
			for(j=0; j<Top::Dimension-1; j++)
			{
				jacobian(i,j) = jacobianlower(i,j)*(1.-ksi)/2.+jacobianupper(i,j)*(1.+ksi)/2.;
			}
			jacobian(i,j) = (xupper[i]-xlower[i])/2.;
		}
	}
	
	template<class TFather, class Topology>
	void GPr<TFather, Topology>::Jacobian(TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv)
	{
		TPZFNMatrix<16> axest,VecMatrix;
		Jacobian(coord,param,VecMatrix);
		VecMatrix.GramSchmidt(axest,jacobian);
		axest.Transpose(&axes);
		jacobian.DeterminantInverse(detjac,jacinv);
		
	}
	
	template<class TFather, class Topology>
	inline void GPr<TFather, Topology>::X(TPZFMatrix<REAL> &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result){
		
		
		int spacedim = coord.Rows();
		TPZFNMatrix<3*TFather::NNodes> lower(spacedim,TFather::NNodes),upper(spacedim,TFather::NNodes);
		TPZManVector<REAL> resultlower(spacedim),resultupper(spacedim);
		int i,j;
		for(i=0; i<spacedim; i++) for(j=0; j<TFather::NNodes; j++)
		{
			lower(i,j) = coord(i,j);
			upper(i,j) = coord(i,j+TFather::NNodes);
		}
		TFather::X(lower,loc,resultlower);
		TFather::X(upper,loc,resultupper);
		REAL ksi = loc[Top::Dimension-1];
		for(i=0; i<spacedim; i++)
		{
			result[i] = resultlower[i]*(1.-ksi)/2.+resultupper[i]*(1.+ksi)/2.;
		}
	}
	
	template<class TFather, class Topology>
	void GPr<TFather, Topology>::Diagnostic(TPZFMatrix<REAL> & coord)
	{
#ifdef LOG4CXX
		LoggerPtr logger(Logger::getLogger("pz.geom.pzgeoextend"));
		
		TPZIntPoints *integ = Top::CreateSideIntegrationRule(Top::NSides-1,3);
		int np = integ->NPoints();
		REAL w;
		TPZVec<REAL> pos(integ->Dimension());
		int ip;
		for(ip=0; ip<np; ip++)
		{
			integ->Point(ip,pos,w);
			TPZFMatrix<REAL> jacobian(coord.Rows(),Top::Dimension);
			Jacobian(coord,pos,jacobian);
			{
				std::stringstream sout;
				sout << "Testing Jacobian\n";
				sout << "position " <<  pos << " jacobian " << jacobian;
				LOGPZ_DEBUG(logger,sout.str());
			}
			TPZFMatrix<REAL> axes(Top::Dimension,coord.Rows());
			TPZFMatrix<REAL> jacinv(Top::Dimension,Top::Dimension);
			REAL detjac;
			Jacobian(coord,pos,jacobian,axes,detjac,jacinv);
			{
				std::stringstream sout;
				sout << "Testing Jacobian\n";
				sout << "position " <<  pos << " jacobian " << jacobian << 
				" axes " << axes << " detjac " << detjac << " jacinv " << jacinv;
				LOGPZ_DEBUG(logger,sout.str());
			}
			TPZVec<REAL> x(coord.Rows());
			X(coord,pos,x);
			{
				std::stringstream sout;
				sout << "Testing X\n";
				sout << "position " <<  pos << " x " << x;
				LOGPZ_DEBUG(logger,sout.str());
			}
		}

		delete integ;
		
#endif
		
	}
	
};

#endif
