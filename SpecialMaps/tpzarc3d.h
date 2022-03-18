/**
 * @file
 * @brief Contains the TPZArc3D class which implements three dimensional arc.
 */

#ifndef TPZARC3D_H
#define TPZARC3D_H

#include "pzgeoel.h"
#include "pznoderep.h"
#include "tpzline.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeomid.h"

#include <iostream>

namespace pzgeom
{
	
	/** 
	 * @author Paulo Cesar de Alvarenga Lucci (Caju)
	 * @ingroup geometry
	 * @brief Implements three dimensional arc. \ref geometry "Geometry"
	 * @since 2007
	 */
	class TPZArc3D : public pzgeom::TPZNodeRep<3,pztopology::TPZLine> {
		
	public:
        typedef pztopology::TPZLine Top;
		/** @brief Number of nodes (connects) */
		enum {NNodes = 3};

		/** @brief Copy constructor with map of nodes */
		TPZArc3D(const TPZArc3D &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : 
            TPZRegisterClassId(&TPZArc3D::ClassId),
            pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap){
			this->fICnBase = cp.fICnBase;
			this->fIBaseCn = cp.fIBaseCn;
			this->fCenter3D = cp.fCenter3D;
                this->finitialVector = cp.finitialVector;
                
			this->fRadius = cp.fRadius;
                this->fAngle = cp.fAngle;
                this->fXcenter = cp.fXcenter;
                this->fYcenter = cp.fYcenter;
		}
        /*
        TPZFNMatrix<9> fICnBase, fIBaseCn;
        TPZManVector< REAL,3 > fCenter3D, finitialVector;
#ifdef REALpzfpcounter
        double fAngle, fRadius, fXcenter, fYcenter;
#else
        REAL fAngle, fRadius, fXcenter, fYcenter;
#endif
*/
		/** @brief Default constructor */
		TPZArc3D() : TPZRegisterClassId(&TPZArc3D::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(),fICnBase(3,3),fIBaseCn(3,3), fCenter3D(3,0.), finitialVector(3,0.),
            fAngle(0.), fRadius(0.), fXcenter(0.), fYcenter(0.) {
		}
		/** @brief Copy constructor */
		TPZArc3D(const TPZArc3D &cp) : TPZRegisterClassId(&TPZArc3D::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp){
			this->fICnBase = cp.fICnBase;
			this->fIBaseCn = cp.fIBaseCn;
			this->fCenter3D = cp.fCenter3D;
            this->finitialVector = cp.finitialVector;
			this->fRadius = cp.fRadius;
            this->fAngle = cp.fAngle;
            this->fXcenter = cp.fXcenter;
            this->fYcenter = cp.fYcenter;
		}

        /** @brief Copy constructor */
        TPZArc3D& operator=(const TPZArc3D &cp)
        {
            pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>::operator=(cp);
            this->fICnBase = cp.fICnBase;
            this->fIBaseCn = cp.fIBaseCn;
            this->fCenter3D = cp.fCenter3D;
            this->finitialVector = cp.finitialVector;
            this->fRadius = cp.fRadius;
            this->fAngle = cp.fAngle;
            this->fXcenter = cp.fXcenter;
            this->fYcenter = cp.fYcenter;
            return *this;
        }

        /** @brief Another copy constructor */
		TPZArc3D(const TPZArc3D &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZArc3D::ClassId),
        pzgeom::TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
			this->fICnBase  = cp.fICnBase;
			this->fIBaseCn  = cp.fIBaseCn;
			this->fCenter3D = cp.fCenter3D;
            this->finitialVector = cp.finitialVector;
			this->fRadius   = cp.fRadius;
            this->fAngle = cp.fAngle;
            this->fXcenter = cp.fXcenter;
            this->fYcenter = cp.fYcenter;
		}
		
		TPZArc3D(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZArc3D::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes), fICnBase(3,3), fIBaseCn(3,3), fCenter3D(3,0.), finitialVector(3,0.),
        fAngle(0.), fRadius(0.), fXcenter(0.), fYcenter(0.) {
			int64_t nnod = nodeindexes.NElements();
			if(nnod != 3)
			{
				std::cout << "Arc geometry created with " << nnod << " nodes, bailing out\n";
				DebugStop();
			}
		}
		
        TPZArc3D(TPZFMatrix<REAL> &coord) : TPZRegisterClassId(&TPZArc3D::ClassId), pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(), fICnBase(3,3), fIBaseCn(3,3), fCenter3D(3,0.), finitialVector(3,0.),
            fAngle(0.), fRadius(0.), fXcenter(0.), fYcenter(0.){
			ComputeAtributes(coord);
		}
		
		/** @brief Initialize the internal data structure of the arc using the coordinates of the nodes */
		void Initialize(TPZGeoEl *refel)
		{
			int nnod = 3;
			TPZFNMatrix<9,REAL> coord(3,nnod);
			int nod, co;
			for(nod=0; nod<3; nod++)
			{
				for(co=0; co<3; co++)
				{
					coord(co,nod) = refel->NodePtr(nod)->Coord(co);
				}
			}
			ComputeAtributes(coord);
		}

        template<class T>
		void X(TPZFMatrix<REAL> &coord,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            /** Computing initialVector = (iniR2 - CenterR2) */
            TPZManVector<T,3> MappedBASE2D(3,0.);
            
            TPZFNMatrix<4,T> RotMatrix(2,2);
            T deflection = fAngle * fRadius * (loc[0] + 1.) / (2.*fRadius);
            RotMatrix(0,0) =  cos(deflection); RotMatrix(0,1) = sin(deflection);
            RotMatrix(1,0) = -sin(deflection); RotMatrix(1,1) = cos(deflection);
            
            /** MappedPoint_R2 = centerCoord + vectorRotated , where Vx = RotationMatrix . Va */
            T centerCoord, vectRotated = 0.;
            for(int i = 0; i < 2; i++)
            {
                vectRotated = 0.;
                for(int j = 0; j < 2; j++) vectRotated += RotMatrix(i,j)*finitialVector[j];
                centerCoord = (1-i)*fXcenter + i*fYcenter;
                MappedBASE2D[i] = centerCoord + vectRotated;
            }
            
            /** Changing Basis of Obtained MappedPoint from R2 to R3 */
            MappedBASE2D[2] = 0.;
            for(int i = 0; i < 3; i++)
            {
                vectRotated = 0.;
                for(int j = 0; j < 3; j++)
                {
                    vectRotated += fIBaseCn.GetVal(i,j)*MappedBASE2D[j];
                }
                result[i] = vectRotated + coord(i,2);
            }

        }
//        void Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv) const;
		
//        template<class T>
//        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            X(coord,loc,result);
//        }
        
        template<class T>
        void GradX(TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            
            /** Computing Axes */
            TPZManVector< T > Vpc(3,0.), Vpa(3,0.), Vpb(3,0.), Vt(3,0.), OUTv(3,0.);
            
            TPZManVector< T > middle(1, 0.);
            X(coord,middle,OUTv);
            
            /** Vector From MappedPoint to Ini */
            Vpa[0] = coord(0,0) - OUTv[0]; Vpa[1] = coord(1,0) - OUTv[1]; Vpa[2] = coord(2,0) - OUTv[2];
            
            /** Vector From MappedPoint to Fin */
            Vpb[0] = coord(0,1) - OUTv[0]; Vpb[1] = coord(1,1) - OUTv[1]; Vpb[2] = coord(2,1) - OUTv[2];
            
            X(coord,par,OUTv);
            
            /** Vector From MappedPoint to Center */
            Vpc[0] = fCenter3D[0] - OUTv[0]; Vpc[1] = fCenter3D[1] - OUTv[1]; Vpc[2] = fCenter3D[2] - OUTv[2];
            
            /** Tangent Vector From Point in the Arc */
            gradx(0) =  Vpa[1]*Vpb[0]*Vpc[1] - Vpa[0]*Vpb[1]*Vpc[1] + Vpa[2]*Vpb[0]*Vpc[2] - Vpa[0]*Vpb[2]*Vpc[2];
            gradx(1) = -Vpa[1]*Vpb[0]*Vpc[0] + Vpa[0]*Vpb[1]*Vpc[0] + Vpa[2]*Vpb[1]*Vpc[2] - Vpa[1]*Vpb[2]*Vpc[2];
            gradx(2) = -Vpa[2]*Vpb[0]*Vpc[0] + Vpa[0]*Vpb[2]*Vpc[0] - Vpa[2]*Vpb[1]*Vpc[1] + Vpa[1]*Vpb[2]*Vpc[1];
            
            T Vtnorm = 0.;
            for(int i = 0; i < 3; i++)
            {
                if( fabs(Vt[i]) < 1.E-12 ) Vt[i] = 0.;
                Vtnorm += gradx(i)*gradx(i);
            }
            if(Vtnorm < 0.) DebugStop();
//            if(sqrt(Vtnorm) < 1e-16) DebugStop();
            T scale = fAngle*fRadius/(2.*sqrt(Vtnorm));
            for(int j = 0; j < 3; j++) gradx(j) = gradx(j)*scale;

        }
		
//        void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
//        }

		
        static std::string TypeName() { return "Arc3D";}
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }
        
		// static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);
		void Print(std::ostream &out) const
		{
			pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Print(out);
			out << "fCenter3D " << fCenter3D << " finitialVector " << finitialVector << std::endl;
			out << "fAngle " << fAngle << " fRadius " << fRadius << " fXcenter " << fXcenter << " fYcenter " << fYcenter << std::endl;
			fICnBase.Print("fICnBase", out);
			fIBaseCn.Print("fIBaseCn", out);
		}
        
        void Read(TPZStream& buf, void* context) override
        {
            pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Read(buf, context);
            fICnBase.Read(buf,0);
            fIBaseCn.Read(buf,0);
            buf.Read<3>( fCenter3D);
            buf.Read<3>(finitialVector);
            buf.Read(&fAngle,1);
            buf.Read(&fRadius,1);
            buf.Read(&fXcenter,1);
            buf.Read(&fYcenter,1);
        }
        
        void Write(TPZStream &buf, int withclassid) const override
        {
            pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Write(buf, withclassid);
            fICnBase.Write(buf,0);
            fIBaseCn.Write(buf,0);
            buf.Write( fCenter3D);
            buf.Write(finitialVector);
            buf.Write(&fAngle,1);
            buf.Write(&fRadius,1);
            buf.Write(&fXcenter,1);
            buf.Write(&fYcenter,1);
		}
        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
    public:
int ClassId() const override;
    public:
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        

	protected:
		
		void ComputeAtributes(TPZFMatrix<REAL> &coord);
		void ComputeR2Points(TPZFMatrix<REAL> &coord, double &xa, double &ya, double &xb, double &yb);
		double ArcAngle(TPZFMatrix<REAL> &coord, double xa, double ya, double xb, double yb) const;
		
		/**
		 * @name Atributes
		 * @{
		 */
		 
		TPZFNMatrix<9> fICnBase, fIBaseCn;
		TPZManVector< REAL,3 > fCenter3D, finitialVector;
#ifdef REALpzfpcounter
		double fAngle, fRadius, fXcenter, fYcenter;
#else
		REAL fAngle, fRadius, fXcenter, fYcenter;
#endif
		/** @} */
		
	};
	
};

#endif
