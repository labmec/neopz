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
		
		REAL fR;
		REAL fr;
		
        TPZFNMatrix<12,REAL> fPhiTheta;
    public:
        
        public:
virtual int ClassId() const;

        
        /** @brief Constructor with list of nodes */
		TPZQuadTorus(TPZVec<long> &nodeindexes) : TPZGeoQuad(nodeindexes), fPhiTheta(3,4,0.)
		{
		}
		
		/** @brief Empty constructor */
		TPZQuadTorus() : TPZGeoQuad(), fPhiTheta(3,4,0.)
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
        
        void SetDataPhiTheta(const TPZFMatrix<REAL> &phitheta)
        {
#ifdef PZDEBUG
            if (phitheta.Rows() != 3 || phitheta.Cols() != 4) {
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

		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "TorusQuad";}
		
		/* @brief Computes the coordinate of a point given in parameter space */
        template<class T>
        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
        
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            TPZFNMatrix<9,T> GradPhi(3,3,0.);
            TPZGeoQuad::GradX(gel, par, GradPhi);
            TPZFNMatrix<6,T> DxDphi(3,3,0.);
            TPZManVector<T,3> ft(3,0.);
            TPZGeoQuad::X(fPhiTheta,par,ft);
            DxDphi(0,0) = -cos(ft[1]) * sin(ft[0]);
            DxDphi(0,1) = -(3. + cos(ft[0])) * sin(ft[1]);
            DxDphi(1,0) = -sin(ft[1]) * sin(ft[0]);
            DxDphi(1,1) = cos(ft[1]) * (3. + cos(ft[0]));
            DxDphi(2,0) = cos(ft[0]);
            DxDphi(2,1) = 0.;
            DxDphi.Multiply(GradPhi, gradx);
            
            TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
            
            int spacedim = fPhiTheta.Rows();
            
            for (int j=0; j<spacedim; j++) {
                minx[j] = fPhiTheta.GetVal(j,0);
                maxx[j] = fPhiTheta.GetVal(j,0);
            }
            
            for(int i = 0; i < 4; i++) {
                for(int j = 0; j < spacedim; j++) {
                    minx[j] = minx[j] < fPhiTheta.GetVal(j,i) ? minx[j]:fPhiTheta.GetVal(j,i);
                    maxx[j] = maxx[j] > fPhiTheta.GetVal(j,i) ? maxx[j]:fPhiTheta.GetVal(j,i);
                }
            }
            REAL delx = 0.;
            for (int j=0; j<spacedim; j++) {
                delx = delx > (maxx[j]-minx[j]) ? delx : (maxx[j]-minx[j]);
            }
            gradx *= 1./delx;

        }
		
        /* @brief Computes the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;

        template<class T>
		void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            
            TPZGeoQuad::X(this->fPhiTheta,loc,result);
            TPZVec <T> toro(3,0.0);
            
            toro[0] = (fR + fr*cos(result[0]))*cos(result[1]);
            toro[1] = (fR + fr*cos(result[0]))*sin(result[1]);
            toro[2] = fr*sin(result[0]);		
            result=toro;
            
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
        
        virtual void Write(TPZStream &buf) const
        {
            pzgeom::TPZGeoQuad::Write(buf);
		}

		
	};
	

    
}

#endif /* defined(__PZ__TPZQuadTorus__) */
