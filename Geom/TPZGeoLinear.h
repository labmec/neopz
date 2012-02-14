/**
 * @file
 * @brief Contains the TPZGeoLinear class which implements the geometry of a one dimensional linear element.
 */
#ifndef TPZGEOLINEARH
#define TPZGEOLINEARH

#include "pznoderep.h"

#include "pzvec.h"
#include "pzeltype.h"
#include "pzfmatrix.h"
#include "tpzline.h"

#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @brief Implements the geometry of a one dimensional linear element. \ref geometry "Geometry"
	 */
	class TPZGeoLinear : public TPZNodeRep<2, pztopology::TPZLine> {
		
	public:
		enum {NNodes = 2};
		
		/** @brief Constructor with list of nodes */
		TPZGeoLinear(TPZVec<int> &nodeindexes) : TPZNodeRep<NNodes, pztopology::TPZLine>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoLinear() : TPZNodeRep<NNodes, pztopology::TPZLine>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoLinear(const TPZGeoLinear &cp,
					 std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoLinear(const TPZGeoLinear &cp) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoLinear(const TPZGeoLinear &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
		{
		}
		
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Linear";} 
		
		/* brief compute the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
        /* @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
		static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);
		
		/** @brief Returns the projection of a given point from "NSide - 1" side to "side". */
		static bool MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide);
		
		static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);
		
		static void Jacobian(TPZFMatrix &nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
							 TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);
		
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
		
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid,
										  int& index);
		
		
	};
	
	// VERSAO ORIGINAL
	//   TPZFNMatrix<9> phi(NNodes, pztopology::TPZLine,1);
	//   TPZFNMatrix<18> dphi(1,NNodes, pztopology::TPZLine);
	//   Shape(param,phi,dphi);
	
	//   int ic;
	//   TPZManVector<REAL,3> v1(3,0.);
	//   REAL mod1 = 0.;
	
	//   for(int i=0; i < NNodes, pztopology::TPZLine; i++) {
	//     for(ic = 0; ic < 3; ic++) {
	//       v1[ic] += coord(ic,i)*dphi(0,i);
	//     }
	//   }
	
	//   for(ic=0; ic<3; ic++) {
	//     mod1 += v1[ic]*v1[ic];
	//   }
	//   mod1 = sqrt(mod1);
	//   jacobian(0,0) = mod1;
	//   detjac = mod1;
	//   jacinv(0,0) = 1./detjac;
	
	
	//  axes.Zero();
	//   for(ic=0; ic<3; ic++) {
	//     axes(0,ic) = v1[ic]/mod1;
	//   }
	
	inline void TPZGeoLinear::X(TPZFMatrix &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result){
		
		int ic;
		REAL xi = loc[0];
		int nrow = coord.Rows();
		for(ic=0; ic<nrow; ic++) result[ic] = coord(ic,0)*(1.-xi)*0.5+coord(ic,1)*(1.+xi)*0.5;
		
		//   TPZFNMatrix<9> phi(NNodes, pztopology::TPZLine,1);
		//   TPZFNMatrix<18> dphi(1,NNodes, pztopology::TPZLine);
		//   Shape(loc,phi,dphi);
		//   int in,ic;
		//   for(in=0; in<3; in++) result[in] = 0.;
		//   for(in = 0; in < NNodes, pztopology::TPZLine; in++) {
		//     for(ic=0; ic<3 ; ic++) {
		//       result[ic] += coord(ic,in)*phi(in,0);
		//     }
		//   }
		
	}
	
	inline bool TPZGeoLinear::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide) {
		TPZTransform Transf = pztopology::TPZLine::SideToSideTransform(TPZGeoLinear::NSides - 1, side);
		SidePar.Resize(SideDimension(side));
		Transf.Apply(InternalPar,SidePar);
		
		JacToSide = Transf.Mult();
		return true;
	}
	
};

#endif
