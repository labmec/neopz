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
        REAL fR = 0;
        REAL fr = 0;
        
        TPZManVector<REAL> fOrigin;
        
        TPZFNMatrix<12,REAL> fPhiTheta;

    public:

        public:
int ClassId() const override;

        
        /** @brief Constructor with list of nodes */
		TPZTriangleTorus(TPZVec<int64_t> &nodeindexes) : TPZGeoTriangle(nodeindexes), fR(0), fr(0), fOrigin(3,0), fPhiTheta(2,3,0.)
		{
		}
		
		/** @brief Empty constructor */
		TPZTriangleTorus() : TPZGeoTriangle(), fR(0), fr(0), fOrigin(3,0), fPhiTheta(2,3,0.)
		{
		}
		
		/** @brief Constructor with node map */
		TPZTriangleTorus(const TPZTriangleTorus &cp,
				   std::map<int64_t,int64_t> & gl2lcNdMap) : TPZGeoTriangle(cp,gl2lcNdMap), fR(cp.fR), fr(cp.fr), fOrigin(cp.fOrigin), fPhiTheta(cp.fPhiTheta)
		{
		}
		
		/** @brief Copy constructor */
		TPZTriangleTorus(const TPZTriangleTorus &cp) : TPZGeoTriangle(cp), fR(cp.fR), fr(cp.fr), fOrigin(cp.fOrigin), fPhiTheta(cp.fPhiTheta)
		{
		}
		
		/** @brief Copy constructor */
		TPZTriangleTorus(const TPZTriangleTorus &cp, TPZGeoMesh &) : TPZGeoTriangle(cp), fR(cp.fR), fr(cp.fr), fOrigin(cp.fOrigin), fPhiTheta(cp.fPhiTheta)
		{
		}
        
        TPZTriangleTorus &operator=(const TPZTriangleTorus &cp)
        {
            TPZGeoTriangle::operator=(cp);
            fR = cp.fR;
            fr = cp.fr;
            fOrigin = cp.fOrigin;
            fPhiTheta = cp.fPhiTheta;
            return *this;
        }
        
        void SetDataPhiTheta(const TPZFMatrix<REAL> &phitheta)
        {
#ifdef PZDEBUG
            if (phitheta.Rows() != 2 || phitheta.Cols() != 3) {
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
		static std::string TypeName() { return "TriangleTorus";}
		
		/* @brief Computes the coordinate of a point given in parameter space */
        template<class T>
        void X(TPZFMatrix<REAL> &coord,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZGeoTriangle::X(this->fPhiTheta,loc,result);
            TPZVec <T> toro(3,0.0);
            
            toro[0] = (fR + fr*cos(result[0]))*cos(result[1]);
            toro[1] = (fR + fr*cos(result[0]))*sin(result[1]);
            toro[2] = fr*sin(result[0]);		
            result=toro;
        }
        
        template<class T>
        void GradX(TPZFMatrix<REAL> &cornerco, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            TPZFNMatrix<6,T> DxDphi(3,2,0.), gradphi(2,2);
            TPZManVector<T,3> ft(3,0.);
            TPZGeoTriangle::X(fPhiTheta,par,ft);
            TPZGeoTriangle::GradX(fPhiTheta, par, gradphi);
            
            DxDphi(0,0) = -cos(ft[1]) * sin(ft[0]);
            DxDphi(0,1) = -(3. + cos(ft[0])) * sin(ft[1]);
            DxDphi(1,0) = -sin(ft[1]) * sin(ft[0]);
            DxDphi(1,1) = cos(ft[1]) * (3. + cos(ft[0]));
            DxDphi(2,0) = cos(ft[0]);
            DxDphi(2,1) = 0.;
            DxDphi.Multiply(gradphi, gradx);

            
        }
		
        /* @brief Computes the jacobian of the map between the master element and deformed element */
//        void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
//        {
//            std::cout << __PRETTY_FUNCTION__ << "PLEASE IMPLEMENT ME!!!\n";
//            DebugStop();
//            //TPZGeoTriangle::Jacobian(gel, param, jacobian , axes, detjac, jacinv);
//        }
        
        template<class T>
		void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            
            TPZManVector<T,2> resloc(2);
            TPZGeoTriangle::X(this->fPhiTheta,loc,resloc);
            
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
            pzgeom::TPZGeoTriangle::Read(buf,0);
            buf.Read(&fR);
            buf.Read(&fr);
            fPhiTheta.Read(buf,0);
        }
        
        void Write(TPZStream &buf, int withclassid) const override
        {
            pzgeom::TPZGeoTriangle::Write(buf, withclassid);
            buf.Write(&fR);
            buf.Write(&fr);
            fPhiTheta.Write(buf, 0);
		}

        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        

	};

    
}

#endif /* defined(__PZ__TPZTriangleTorus__) */
