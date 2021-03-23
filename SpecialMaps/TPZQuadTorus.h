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
	private:
		
		REAL fR = 0.;
		REAL fr = 0.;
        
        TPZManVector<REAL> fOrigin;
		
        TPZFNMatrix<12,REAL> fPhiTheta;
    public:
        
        public:
int ClassId() const override;

        static bool IsLinearMapping(int side)
        {
            return false;
        }
        /** @brief Constructor with list of nodes */
		TPZQuadTorus(TPZVec<int64_t> &nodeindexes) : TPZGeoQuad(nodeindexes), fOrigin(3,0.), fPhiTheta(3,4,0.)
		{
		}
		
		/** @brief Empty constructor */
		TPZQuadTorus() : TPZGeoQuad(), fOrigin(3,0.), fPhiTheta(3,4,0.)
		{
		}
		
		/** @brief Constructor with node map */
		TPZQuadTorus(const TPZQuadTorus &cp,
				   std::map<int64_t,int64_t> & gl2lcNdMap) : TPZGeoQuad(cp,gl2lcNdMap), fR(cp.fR), fr(cp.fr), fOrigin(cp.fOrigin), fPhiTheta(cp.fPhiTheta)
		{
		}
		
		/** @brief Copy constructor */
		TPZQuadTorus(const TPZQuadTorus &cp) : TPZGeoQuad(cp), fR(cp.fR), fr(cp.fr), fOrigin(cp.fOrigin), fPhiTheta(cp.fPhiTheta)
		{
		}
		
		/** @brief Copy constructor */
		TPZQuadTorus(const TPZQuadTorus &cp, TPZGeoMesh &) : TPZGeoQuad(cp), fR(cp.fR), fr(cp.fr), fOrigin(cp.fOrigin), fPhiTheta(cp.fPhiTheta)
		{
		}
        
        TPZQuadTorus &operator=(const TPZQuadTorus &cp)
        {
            TPZGeoQuad::operator=(cp);
            fR = cp.fR;
            fr = cp.fr;
            fOrigin = cp.fOrigin;
            fPhiTheta = cp.fPhiTheta;
            return *this;
        }
        
        void SetDataPhiTheta(const TPZFMatrix<REAL> &phitheta)
        {
#ifdef PZDEBUG
            if (phitheta.Rows() != 2 || phitheta.Cols() != 4) {
                DebugStop();
            }
#endif
            fPhiTheta = phitheta;
        }
		
        void SetDataRadius(const REAL &R, const REAL &r)
        {
#ifdef PZDEBUG
            if (R < r) 
			{
                DebugStop();
            }
#endif
            fR = R;
            fr = r;			
        }
        
        void SetOrigin(TPZVec<REAL> &origin)
        {
            fOrigin = origin;
        }

		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "QuadTorus";}
		
		/* @brief Computes the coordinate of a point given in parameter space */
        
        template<class T>
        void GradX(TPZFMatrix<REAL> &cornerco, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            
            
            TPZFNMatrix<6,T> DxDphi(3,2,0.), gradphi(2,2);
            TPZManVector<T,3> ft(3,0.);
            TPZGeoQuad::X(fPhiTheta,par,ft);
            TPZGeoQuad::GradX(fPhiTheta, par, gradphi);
            
            DxDphi(0,0) = -cos(ft[1]) * sin(ft[0]);
            DxDphi(0,1) = -(3. + cos(ft[0])) * sin(ft[1]);
            DxDphi(1,0) = -sin(ft[1]) * sin(ft[0]);
            DxDphi(1,1) = cos(ft[1]) * (3. + cos(ft[0]));
            DxDphi(2,0) = cos(ft[0]);
            DxDphi(2,1) = 0.;
            DxDphi.Multiply(gradphi, gradx);
            
            
        }
        /* @brief Computes the jacobian of the map between the master element and deformed element */
//        void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;

        template<class T>
		void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZManVector<T,2> resloc(2);
            TPZGeoQuad::X(this->fPhiTheta,loc,resloc);
            
            result[0] = (fR + fr*cos(resloc[0]))*cos(resloc[1]);
            result[1] = (fR + fr*cos(resloc[0]))*sin(resloc[1]);
            result[2] = fr*sin(resloc[0]);
            
        }

		
		// static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

		// /** @brief Creates a geometric element according to the type of the father element */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);
		
        void Read(TPZStream& buf, void* context) override
        {
            pzgeom::TPZGeoQuad::Read(buf,0);
        }
        
        void Write(TPZStream &buf, int withclassid) const override
        {
            pzgeom::TPZGeoQuad::Write(buf, withclassid);
		}

        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        

	};
	

    
}

#endif /* defined(__PZ__TPZQuadTorus__) */
