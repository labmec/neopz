//
//  TPZTriangleTorus.h
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#ifndef __PZ__TPZTriangleTorus__
#define __PZ__TPZTriangleTorus__

#include <iostream>
#include "pzgeotriangle.h"

namespace pzgeom {
    
    class TPZTriangleTorus : public TPZGeoTriangle
    {
        int fNumWaves;
        TPZManVector<REAL,3> fWaveDir;

    public:
        
        /** @brief Constructor with list of nodes */
		TPZTriangleTorus(TPZVec<long> &nodeindexes) : TPZGeoTriangle(nodeindexes), fNumWaves(0), fWaveDir()
		{
		}
		
		/** @brief Empty constructor */
		TPZTriangleTorus() : TPZGeoTriangle(), fNumWaves(0), fWaveDir()
		{
		}
		
		/** @brief Constructor with node map */
		TPZTriangleTorus(const TPZTriangleTorus &cp,
				   std::map<long,long> & gl2lcNdMap) : TPZGeoTriangle(cp,gl2lcNdMap), fNumWaves(cp.fNumWaves), fWaveDir(cp.fWaveDir)
		{
		}
		
		/** @brief Copy constructor */
		TPZTriangleTorus(const TPZTriangleTorus &cp) : TPZGeoTriangle(cp), fNumWaves(cp.fNumWaves), fWaveDir(cp.fWaveDir)
		{
		}
		
		/** @brief Copy constructor */
		TPZTriangleTorus(const TPZTriangleTorus &cp, TPZGeoMesh &) : TPZGeoTriangle(cp), fNumWaves(cp.fNumWaves), fWaveDir(cp.fWaveDir)
		{
		}
        
        void SetData(TPZVec<REAL> &wavedir, int numwaves)
        {
#ifdef DEBUG
            if (wavedir.size() != 3) {
                DebugStop();
            }
#endif
            fWaveDir = wavedir;
            fNumWaves = numwaves;
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
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            std::cout << __PRETTY_FUNCTION__ << "PLEASE IMPLEMENT ME!!!\n";
            DebugStop();
            TPZGeoTriangle::Jacobian(gel, param, jacobian , axes, detjac, jacinv);
        }
        
		void X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            std::cout << __PRETTY_FUNCTION__ << "PLEASE IMPLEMENT ME!!!\n";
            DebugStop();
            TPZGeoTriangle::X(nodes,loc,result);
        }
		
		
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
		
        void Read(TPZStream &buf,void *context)
        {
            pzgeom::TPZGeoTriangle::Read(buf,0);
            buf.Read(&fNumWaves,1);
            TPZSaveable::ReadObjects<3>(buf, fWaveDir);
        }
        
        void Write(TPZStream &buf)
        {
            pzgeom::TPZGeoTriangle::Write(buf);
            buf.Write(&fNumWaves,1);
            TPZSaveable::WriteObjects(buf, fWaveDir);
		}

		
	};

    
}

#endif /* defined(__PZ__TPZTriangleTorus__) */
