//
//  TPZTriangleSphere.h
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#ifndef __PZ__TPZTriangleSphere__
#define __PZ__TPZTriangleSphere__

#include <iostream>
#include "pzgeotriangle.h"

namespace pzgeom {
    
    class TPZTriangleSphere : public TPZGeoTriangle
    {
        //int fNumWaves;
       
        TPZVec<REAL> fXc;
		 REAL fR;

    public:
        
        /** @brief Constructor with list of nodes */
		TPZTriangleSphere(TPZVec<long> &nodeindexes) : TPZGeoTriangle(nodeindexes), fXc(0), fR(0)
		{
		}
		
		/** @brief Empty constructor */
		TPZTriangleSphere() : TPZGeoTriangle(), fXc(0), fR(0)
		{
		}
		
		/** @brief Constructor with node map */
		TPZTriangleSphere(const TPZTriangleSphere &cp,
				   std::map<long,long> & gl2lcNdMap) : TPZGeoTriangle(cp,gl2lcNdMap), fXc(cp.fXc), fR(cp.fR)
		{
		}
		
		/** @brief Copy constructor */
		TPZTriangleSphere(const TPZTriangleSphere &cp) : TPZGeoTriangle(cp), fXc(cp.fXc), fR(cp.fR)
		{
		}
		
		/** @brief Copy constructor */
		TPZTriangleSphere(const TPZTriangleSphere &cp, TPZGeoMesh &) : TPZGeoTriangle(cp), fXc(cp.fXc), fR(cp.fR)
		{
		}
        
        void SetData(TPZVec<REAL> &Xc, const REAL R)
        {
#ifdef DEBUG
            if (Xc.size() != 3 || R == 0.0 ) {
                DebugStop();
            }
#endif
            fXc = Xc;
            fR = R;
        }

		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Wavy";}
		
		/* @brief Computes the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
        /* @brief Computes the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;
        
		void X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const;
		
		
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
		
        void Read(TPZStream &buf,void *context)
        {
            pzgeom::TPZGeoTriangle::Read(buf,0);
            buf.Read(&fR,1);
        }
        
        void Write(TPZStream &buf)
        {
            pzgeom::TPZGeoTriangle::Write(buf);
            buf.Write(&fR,1);
		}

		
	};

    
}

#endif /* defined(__PZ__TPZTriangleSphere__) */
