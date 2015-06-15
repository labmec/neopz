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
    
    template <class GeomTriang = pzgeom::TPZGeoTriangle>
    class TPZTriangleSphere : public GeomTriang
    {
        //int fNumWaves;
       
        TPZVec<REAL> fXc;
		 REAL fR;

    public:
        
        /** @brief Constructor with list of nodes */
		TPZTriangleSphere(TPZVec<long> &nodeindexes) : GeomTriang(nodeindexes), fXc(0.), fR(0.)
		{
		}
		
		/** @brief Empty constructor */
		TPZTriangleSphere() : GeomTriang(), fXc(0.), fR(0.)
		{
		}
		
		/** @brief Constructor with node map */
		TPZTriangleSphere(const TPZTriangleSphere &cp,
				   std::map<long,long> & gl2lcNdMap) : GeomTriang(cp,gl2lcNdMap), fXc(cp.fXc), fR(cp.fR)
		{
		}
		
		/** @brief Copy constructor */
		TPZTriangleSphere(const TPZTriangleSphere &cp) : GeomTriang(cp), fXc(cp.fXc), fR(cp.fR)
		{
		}
		
		/** @brief Copy constructor */
		TPZTriangleSphere(const TPZTriangleSphere &cp, TPZGeoMesh &) : GeomTriang(cp), fXc(cp.fXc), fR(cp.fR)
		{
		}
        
        void SetData(const REAL R, TPZVec<REAL> &Xc)
        {
#ifdef DEBUG
            if (Xc.size() != 3 || R == 0.0 ) {
                DebugStop();
            }
#endif
            fXc = Xc;
            fR = R;
        }
        
        /** @brief declare geometry as blended element */
        bool IsGeoBlendEl() const;

		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "TPZTriangleSphere";}
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }
		
		/* @brief Computes the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            X(coord,loc,result);
            
            
            TPZManVector<REAL,3> xqsi(3,0.); // will store (x,y,z) from (qsi,eta)
            TPZManVector<REAL,3> XtminusXc(3,0.0); // will store (x,y,z)-xc
            GeomTriang::X(gel,loc,xqsi);
            TPZManVector<REAL,3> Xc = fXc;
            
            REAL NormValue = 0.0;
            
            for (int i = 0; i < 3; i++) {
                XtminusXc[i] = xqsi[i] - Xc[i];
                NormValue += XtminusXc[i]*XtminusXc[i];
            }
            
//            XtminusXc[0,0)= xqsi[0] - Xc[0];
//            XtminusXc(1,0)= xqsi[1] - Xc[1];
//            XtminusXc(2,0)= xqsi[2] - Xc[2];
            
            NormValue = sqrt(NormValue);
            
            TPZManVector<REAL,3> Xsphere(3,0.0);
            for(int i=0; i<3; i++)
            {
                result[i] = (fR/NormValue)*(XtminusXc[i])+ Xc[i];
            }
//            Xsphere[0] = (fR/NormValue)*(XtminusXc(0,0)+Xc[0]);
//            Xsphere[1] = (fR/NormValue)*(XtminusXc(1,0)+Xc[1]);
//            Xsphere[2] = (fR/NormValue)*(XtminusXc(2,0)+Xc[2]);		
//            result=Xsphere;
            
            
        }
		
        /* @brief Computes the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZManVector<REAL,3> XTriangle(3,0.0);
            TPZFNMatrix<9,REAL> GradOneoverNorm(3,1,0.0);
            TPZFNMatrix<9,REAL> TensorXtGradX(3,2,0.0);
            TPZFNMatrix<9,REAL> GradXtScaled(3,3,0.0);
            
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
            
            GeomTriang::X(gel,param,XTriangle);
            TPZFMatrix<REAL> XtminusXc(3,1,0.0);
            TPZVec<REAL> Xc= fXc;
            
            for (int i = 0; i < XTriangle.size(); i++)
            { XtminusXc(i,0)= XTriangle[i] - Xc[i]; }
            REAL NormValue = Norm(XtminusXc);
            
            GeomTriang::Jacobian(gel, param, jacobian , axes, detjac, jacinv);
            TPZFNMatrix<6> axest(3,2);
            axes.Transpose(&axest);
            
            TPZFNMatrix<9,REAL> GradXt;
            axest.Multiply(jacobian, GradXt);
            
            GradOneoverNorm(0,0) = XtminusXc(0,0)*GradXt(0,0)+XtminusXc(1,0)*GradXt(1,0)+XtminusXc(2,0)*GradXt(2,0);
            GradOneoverNorm(1,0) = XtminusXc(0,0)*GradXt(0,1)+XtminusXc(1,0)*GradXt(1,1)+XtminusXc(2,0)*GradXt(2,1);
            GradOneoverNorm = -(fR/(NormValue*NormValue*NormValue))*GradOneoverNorm;
            
            TensorXtGradX(0,0)= XtminusXc(0,0)*GradOneoverNorm(0,0);
            TensorXtGradX(0,1)= XtminusXc(0,0)*GradOneoverNorm(1,0);
            TensorXtGradX(1,0)= XtminusXc(1,0)*GradOneoverNorm(0,0);
            TensorXtGradX(1,1)= XtminusXc(1,0)*GradOneoverNorm(1,0);
            TensorXtGradX(2,0)= XtminusXc(2,0)*GradOneoverNorm(0,0);
            TensorXtGradX(2,1)= XtminusXc(2,0)*GradOneoverNorm(1,0);
            
            TPZFNMatrix<9,REAL> VecMatrix;
            VecMatrix=((fR/(NormValue))*GradXt)+TensorXtGradX;
            
            
            
            //-----
			/*
			TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
			
			int spacedim = coord.Rows();
			
			for (int j=0; j<spacedim; j++) { minx[j] = coord.GetVal(j,0); maxx[j] = coord.GetVal(j,0);}
			
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < spacedim; j++) {
					minx[j] = minx[j] < coord.GetVal(j,i) ? minx[j]:coord.GetVal(j,i);
					maxx[j] = maxx[j] > coord.GetVal(j,i) ? maxx[j]:coord.GetVal(j,i);
				}
			}
            
			REAL delx = 0.;
			for (int j=0; j<spacedim; j++) {
				delx = delx > (maxx[j]-minx[j]) ? delx : (maxx[j]-minx[j]);
			}
            
			VecMatrix *= 1./delx;
            */
            //-----
            
            
            
            
			VecMatrix.GramSchmidt(axest,jacobian);
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
        
		//void X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const;
		
		
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
