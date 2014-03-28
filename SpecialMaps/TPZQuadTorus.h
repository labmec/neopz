//
//  TPZQuadTorus.h
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#ifndef __PZ__TPZQuadTorus__
#define __PZ__TPZQuadTorus__

#include <iostream>
#include "pzgeoquad.h"

namespace pzgeom {
    
    class TPZQuadTorus : public TPZGeoQuad
    {
        TPZFNMatrix<8,REAL> fPhiTheta;
    public:
        
        /** @brief Constructor with list of nodes */
		TPZQuadTorus(TPZVec<long> &nodeindexes) : TPZGeoQuad(nodeindexes), fPhiTheta(4,2,0.)
		{
		}
		
		/** @brief Empty constructor */
		TPZQuadTorus() : TPZGeoQuad(), fPhiTheta(4,2,0.)
		{
		}
		
		/** @brief Constructor with node map */
		TPZQuadTorus(const TPZQuadTorus &cp,
				   std::map<long,long> & gl2lcNdMap) : TPZGeoQuad(cp,gl2lcNdMap), fPhiTheta(cp.fPhiTheta)
		{
		}
		
		/** @brief Copy constructor */
		TPZQuadTorus(const TPZQuadTorus &cp) : TPZGeoQuad(cp), fPhiTheta(cp.fPhiTheta)
		{
		}
		
		/** @brief Copy constructor */
		TPZQuadTorus(const TPZQuadTorus &cp, TPZGeoMesh &) : TPZGeoQuad(cp), fPhiTheta(cp.fPhiTheta)
		{
		}
        
        void SetData(const TPZFMatrix<REAL> &phitheta)
        {
#ifdef DEBUG
            if (phitheta.Rows() != 4 || phitheta.Cols() != 2) {
                DebugStop();
            }
#endif
            fPhiTheta = phitheta;
        }

		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "TorusQuad";}
		
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
            std::cout << __PRETTY_FUNCTION__ << "PLEASE IMPLEMENT ME!!!\n";
            DebugStop();
            TPZGeoQuad::Jacobian(gel, param, jacobian , axes, detjac, jacinv);
        }
        
		void X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            std::cout << __PRETTY_FUNCTION__ << "PLEASE IMPLEMENT ME!!!\n";
            DebugStop();
            TPZGeoQuad::X(nodes,loc,result);
        }
		
		
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
		
        void Read(TPZStream &buf,void *context)
        {
            pzgeom::TPZGeoQuad::Read(buf,0);
        }
        
        void Write(TPZStream &buf)
        {
            pzgeom::TPZGeoQuad::Write(buf);
		}

		
	};

    
}

#endif /* defined(__PZ__TPZQuadTorus__) */
