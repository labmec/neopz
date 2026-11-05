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
     * @author Philippe Devloo (refactor in 2025)
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
			this->fRotationTensor = cp.fRotationTensor;
			this->fCenter3D = cp.fCenter3D;
                
			this->fRadius = cp.fRadius;
            this->fAngle = cp.fAngle;
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
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(),fRotationTensor(3,3), fCenter3D(3,0.), 
            fAngle(0.), fRadius(0.) {
		}
		/** @brief Copy constructor */
		TPZArc3D(const TPZArc3D &cp) : TPZRegisterClassId(&TPZArc3D::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp){
			this->fRotationTensor = cp.fRotationTensor;
			this->fCenter3D = cp.fCenter3D;
			this->fRadius = cp.fRadius;
            this->fAngle = cp.fAngle;
		}

        /** @brief Copy constructor */
        TPZArc3D& operator=(const TPZArc3D &cp)
        {
            pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>::operator=(cp);
            this->fRotationTensor = cp.fRotationTensor;
            this->fCenter3D = cp.fCenter3D;
            this->fRadius = cp.fRadius;
            this->fAngle = cp.fAngle;
            return *this;
        }

        /** @brief Another copy constructor */
		TPZArc3D(const TPZArc3D &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZArc3D::ClassId),
        pzgeom::TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
			this->fRotationTensor  = cp.fRotationTensor;
			this->fCenter3D = cp.fCenter3D;
			this->fRadius   = cp.fRadius;
            this->fAngle = cp.fAngle;
		}
		
		TPZArc3D(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZArc3D::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes), fRotationTensor(3,3), fCenter3D(3,0.), 
        fAngle(0.), fRadius(0.){
			int64_t nnod = nodeindexes.NElements();
			if(nnod != 3)
			{
				std::cout << "Arc geometry created with " << nnod << " nodes, bailing out\n";
				DebugStop();
			}
		}

        TPZArc3D(TPZFMatrix<REAL> &coord) : TPZRegisterClassId(&TPZArc3D::ClassId), pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(), fRotationTensor(3,3), fCenter3D(3,0.),
            fAngle(0.), fRadius(0.) {
			ComputeAtributes(coord);
		}
		
		/** @brief Initialize the internal data structure of the arc using the coordinates of the nodes */
		void Initialize(TPZGeoEl *refel)
		{
            if(fNodeIndexes[2] == fNodeIndexes[0]) {
                return;
            }
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

        void SetAttributes(TPZGeoMesh *gmesh, const TPZVec<REAL> &center, const TPZFMatrix<REAL> &rotationmatrix, REAL radius, REAL openingangle) {
            fCenter3D = center;
            fRotationTensor = rotationmatrix;
            fRadius = radius;
            fAngle = openingangle;
#ifdef PZDEBUG
            if(fRadius <= 0.) {
                std::cout << "TPZArc3D::SetAttributes radius <= 0 "<< fRadius << std::endl;
                DebugStop();
            }
            TPZManVector<REAL,3> xnode(3,0.), xcomp(3,0.);
            TPZVec<REAL> loc(1);
            for(int i=-1; i<=1; i+=2) {
                loc[0] = i;
                TPZFMatrix<REAL> coord(3,3,0.);
                X(coord,loc, xcomp);
                gmesh->NodeVec()[fNodeIndexes[(i+1)/2]].GetCoordinates(xnode);
                REAL dist = 0.;
                for(int j=0; j<3; j++) {
                    dist += (xnode[j]-xcomp[j])*(xnode[j]-xcomp[j]);
                }
                dist = sqrt(dist);
                if(dist > 1.e-6) {
                    std::cout << "TPZArc3D::SetAttributes node coordinates do not match "<< dist << std::endl;
                    DebugStop();
                }
            }
#endif
        }

        template<class T>
		void X(TPZFMatrix<REAL> &coord,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            T angle = (loc[0]+1.)*fAngle/2.;
            for(int i=0; i<3; i++) {
                T vecdir = cos(angle)*fRotationTensor(i,0)+sin(angle)*fRotationTensor(i,1);
                result[i] = fCenter3D[i]+fRadius*vecdir;
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
            T dangle = fAngle/2.;
            T angle = (par[0]+1.)*fAngle/2.;
            for(int i=0; i<3; i++) {
                T vecdir = -sin(angle)*fRotationTensor(i,0)+cos(angle)*fRotationTensor(i,1);
                gradx(i,0) = fRadius*dangle*vecdir;
            }


        }
		

		
        static std::string TypeName() { return "Arc3D";}
        
        static bool IsLinearMapping(int side)
        {
            return false;
        }
        
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);
		void Print(std::ostream &out) const
		{
			pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Print(out);
			out << "fCenter3D " << fCenter3D <<  std::endl;
			out << "fAngle " << fAngle << " fRadius " << fRadius << std::endl;
			fRotationTensor.Print("fRotationTensor", out);
		}
        
        void Read(TPZStream& buf, void* context) override
        {
            pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Read(buf, context);
            fRotationTensor.Read(buf,0);
            buf.Read<3>( fCenter3D);
            buf.Read(&fAngle,1);
            buf.Read(&fRadius,1);
        }
        
        void Write(TPZStream &buf, int withclassid) const override
        {
            pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Write(buf, withclassid);
            fRotationTensor.Write(buf,0);
            buf.Write( fCenter3D);
            buf.Write(&fAngle,1);
            buf.Write(&fRadius,1);
		}
        
        //virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
    public:
int ClassId() const override;
    public:
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        

	protected:
		/// @brief Compute the configuration of the arc based on 3 points
		/// @param coord is composed of 3 points defined column wise
        /// the first and second column contain the endpoints of the arc
        /// the third point is an intermediate point through which the arc will pass
		void ComputeAtributes(TPZFMatrix<REAL> &coord);
		
        static void ComputeCenter(TPZFMatrix<REAL> &coord, TPZVec<REAL> &center);
		/**
		 * @name Atributes
		 * @{
		 */
		/// @brief Rotation matrix. 
        /// The first column points from the center to the first point
        /// The second column is orthogonal to the first and to the normal of the arc plane
        /// the third column is the normal to the arc plane
		TPZFNMatrix<9> fRotationTensor;
		/// @brief The center of the arc in R3 coordinates
		TPZManVector< REAL,3 > fCenter3D;
#ifdef REALpzfpcounter
		double fAngle, fRadius;
#else
		/// @brief The opening angle of the arc
        REAL fAngle;
        /// @brief The radius of the arc
		REAL fRadius;
#endif
		/** @} */
		
	};
	
};

#endif
