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
	
    template<class GeomQuad = pzgeom::TPZGeoQuad>
	class TPZQuadSphere : public GeomQuad
	{
	private:
		REAL fR;
		TPZVec<REAL> fxc;
		
	public:
		
		/** @brief Constructor with list of nodes */
		TPZQuadSphere(TPZVec<long> &nodeindexes) : GeomQuad(nodeindexes), fR(0.), fxc(3,0.)
		{
		}
		
		/** @brief Empty constructor */
		TPZQuadSphere() : GeomQuad(), fR(0.), fxc(3,0.)
		{
		}
		
		/** @brief Constructor with node map */
		TPZQuadSphere(const TPZQuadSphere &cp,
									std::map<long,long> & gl2lcNdMap) : GeomQuad(cp,gl2lcNdMap), fR(0.), fxc(3,0.)
		{
		}
		
		/** @brief Copy constructor */
		TPZQuadSphere(const TPZQuadSphere &cp) : GeomQuad(cp), fR(cp.fR), fxc(cp.fxc)
		{
		}
		
		/** @brief Copy constructor */
		TPZQuadSphere(const TPZQuadSphere &cp, TPZGeoMesh &) : GeomQuad(cp), fR(cp.fR), fxc(cp.fxc)
		{
		}
        
		/** @brief declare geometry as blended element */
        bool IsGeoBlendEl() const;

		
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
            TPZManVector<REAL,3> xqsi(3,0.); // will store (x,y,z) from (qsi,eta)
            TPZManVector<REAL,3> xqsiLxc(3,0.); // will store (x,y,z)-xc
            GeomQuad::X(gel,loc,xqsi); // gives the map (qsi,eta) to (x,y,z)
            
            REAL norm = 0.;
            for (int i = 0; i < 3; i++) { // Does xqsi-xc and calculates its norm
                xqsiLxc[i] = xqsi[i] - fxc[i];
                norm += xqsiLxc[i] * xqsiLxc[i];
            }
            norm = sqrt(norm);
            
            for (int i = 0; i < 3; i++) {
                result[i] = fxc[i] + xqsiLxc[i] * fR / norm;
            }
                
            
		}
		
		/* @brief Computes the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
		{
           
            
            // will first do dxdqsi and then d(phi*v)/dx, finaly d(phi*v)/dx * dxdqsi, where phi = 1/norm(xqsi - xc) and v = (xqsi - xc) * fR
            
            TPZManVector<REAL,3> xqsi(3,0.); // will store (x,y,z) from (qsi,eta)
            TPZManVector<REAL,3> xqsiLxc(3,0.); // will store (x,y,z)-xc
            GeomQuad::X(gel,param,xqsi); // gives the map (qsi,eta) to (x,y,z)
            REAL norm = 0.;
            for (int i = 0; i < 3; i++) { // Does xqsi-xc and calculates its norm
                xqsiLxc[i] = xqsi[i] - fxc[i];
                norm += xqsiLxc[i] * xqsiLxc[i];
            }
            norm = sqrt(norm);
            
            TPZFNMatrix<6,REAL> dxdqsi(3,2,0.); // But it is a (3,2) matrix. It is set (3,3) because of the later products
            GeomQuad::Jacobian(gel, param, jacobian, axes, detjac, jacinv); // first calculate the derivative dxdqsi (note a lot of dummies in the parameters)
            TPZFMatrix<REAL> axest;
            axes.Transpose(&axest);
            axest.Multiply(jacobian, dxdqsi);
            jacobian.Zero();
            
            
            TPZFNMatrix<3,REAL> gradphi(3,1,0.), v(3,1,0.); // here phi = 1/norm(xqsi - xc) and v = (xqsi - xc) * fR
            TPZFNMatrix<9,REAL> gradv(3,3,0.); // here v = (xqsi - xc) * fR
            REAL phi = 1./norm;
            
            for (int i = 0; i < 3; i++) {
                v(i,0) = xqsiLxc[i] * fR;
                gradv(i,i) = fR;
                gradphi(i,0) = - (1. / (norm*norm*norm) ) * xqsiLxc[i];
            }
            
            TPZFNMatrix <9,REAL> DphivDx(3,3,0.); // will store d(phi*v)/dx
            DphivDx = TensorProd(v,gradphi) + phi*gradv;
            
            TPZFNMatrix <6,REAL> GradX(3,2,0.); // stores d(phi*v)/dx * dxdqsi
            DphivDx.Multiply(dxdqsi, GradX);
            
            /*
             // Amplifying
             TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
             
             int spacedim = coord.Rows();
             
             for (int j=0; j<spacedim; j++) {
             minx[j] = coord.GetVal(j,0);
             maxx[j] = coord.GetVal(j,0);
             }
             TPZFNMatrix<6,REAL> VecMatrix(3,2,0.);
             for(int i = 0; i < 4; i++) {
             for(int j = 0; j < spacedim; j++) {
             minx[j] = minx[j] < coord.GetVal(j,i) ? minx[j]:coord.GetVal(j,i);
             maxx[j] = maxx[j] > coord.GetVal(j,i) ? maxx[j]:coord.GetVal(j,i);
             }
             }
             REAL delx = 0.;
             for (int j=0; j<spacedim; j++) {
             delx = delx > (maxx[j]-minx[j]) ? delx : (maxx[j]-minx[j]);
             }
             
             GradX *= 1./delx;
             */
            
            //GradX.Print("GradX");
            
            GradX.GramSchmidt(axest,jacobian);
            axest.Transpose(&axes);
            detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
            
            if(IsZero(detjac))
            {
                detjac = ZeroTolerance();
            }
            
            jacinv(0,0) =  jacobian(1,1)/detjac;
            jacinv(1,1) =  jacobian(0,0)/detjac;
            jacinv(0,1) = -jacobian(0,1)/detjac;
            jacinv(1,0) = -jacobian(1,0)/detjac;
            
            /*
             // Left without the amplification for learning purposes
             jacobian *= delx;
             jacinv *= 1./delx;
             detjac *= (delx*delx);
             */
            
		}
        
        static TPZFMatrix<REAL> TensorProd(TPZFMatrix<REAL> &vec1, TPZFMatrix<REAL> &vec2)
        {
            TPZFNMatrix<9,REAL> res(3,3,0.);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    res(i,j) = vec1(i,0) * vec2(j,0);
                }
            }
            return res;
        }
		
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
            // ???
//            pzgeom::TPZGeoQuad::Read(buf, 0);
            buf.Read(&fR,1);
            
		}
		
		void Write(TPZStream &buf)
		{
			std::cout << __PRETTY_FUNCTION__ << "PLEASE IMPLEMENT ME!!!\n";
			DebugStop();
            // ???
//            pzgeom::TPZGeoQuad::Write(buf);
            buf.Write(&fR,1);
		}
		
        bool IsLinearMapping() const
        {
            return false;
        }
		
	};
	
	
}

#endif /* defined(__PZ__TPZQuadSphere__) */
