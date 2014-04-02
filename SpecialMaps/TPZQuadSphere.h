//
//  TPZQuadSphere.h
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#ifndef __PZ__TPZQuadSphere__
#define __PZ__TPZQuadSphere__

#include <iostream>
#include "pzgeoquad.h"

namespace pzgeom {
	
	class TPZQuadSphere : public TPZGeoQuad
	{
	private:
		REAL fR;
		TPZVec<REAL> fxc;
		
	public:
		
		/** @brief Constructor with list of nodes */
		TPZQuadSphere(TPZVec<long> &nodeindexes) : TPZGeoQuad(nodeindexes), fR(0.), fxc(3,0.)
		{
		}
		
		/** @brief Empty constructor */
		TPZQuadSphere() : TPZGeoQuad(), fR(0.), fxc(3,0.)
		{
		}
		
		/** @brief Constructor with node map */
		TPZQuadSphere(const TPZQuadSphere &cp,
									std::map<long,long> & gl2lcNdMap) : TPZGeoQuad(cp,gl2lcNdMap), fR(0.), fxc(3,0.)
		{
		}
		
		/** @brief Copy constructor */
		TPZQuadSphere(const TPZQuadSphere &cp) : TPZGeoQuad(cp), fR(cp.fR), fxc(cp.fxc)
		{
		}
		
		/** @brief Copy constructor */
		TPZQuadSphere(const TPZQuadSphere &cp, TPZGeoMesh &) : TPZGeoQuad(cp), fR(cp.fR), fxc(cp.fxc)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "QuadSphere";}
		
		void SetData(const REAL &R, const TPZVec<REAL> &xc)
		{
			fR = R;
			fxc = xc;
#ifdef DEBUG
			if (R <= 0) {
				PZError << "R must be positive!\n";
				DebugStop();
			}
			if (xc.NElements() != 3) {
				PZError << "XCenter must have 3 coordinates!\n";
				DebugStop();
			}
#endif
		}
		
		/* @brief Computes the coordinate of a point given in parameter space */
		void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
		{
			TPZFNMatrix<3*NNodes> coord(3,NNodes);
			CornerCoordinates(gel, coord);
			X(coord,loc,result);
		}
		
		/* @brief Computes the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
		{
			TPZFNMatrix<3*NNodes> coord(3,NNodes);
			CornerCoordinates(gel, coord);
			Jacobian(coord, param, jacobian, axes, detjac, jacinv);
		}
		
		/** @brief Computes the jacobian*/
		void Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;
		
		void X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const;
		
		
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
		
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
																			TPZVec<long>& nodeindexes,
																			int matid,
																			long& index);
		
		void Read(TPZStream &buf,void *context)
		{
			std::cout << __PRETTY_FUNCTION__ << "PLEASE IMPLEMENT ME!!!\n";
			DebugStop();
		}
		
		void Write(TPZStream &buf)
		{
			std::cout << __PRETTY_FUNCTION__ << "PLEASE IMPLEMENT ME!!!\n";
			DebugStop();
		}
		
		
	};
	
	
}

#endif /* defined(__PZ__TPZQuadSphere__) */
