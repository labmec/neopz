/**
 * @file
 * @brief Contains the TPZEllipse3D class which defines a linear geometric element which maps a line segment to an ellipse.
 */

#ifndef TPZELLIPSE3D_H
#define TPZELLIPSE3D_H

#include "pzgeoel.h"
#include "pznoderep.h"
#include "tpzline.h"

#include <iostream>

class TPZGeoMesh;

namespace pzgeom
{
    
	/** 
	 * @brief Defines a linear geometric element which maps a line segment to an ellipse. \ref geometry "Geometry"
	 * @ingroup geometry
	 */
	class TPZEllipse3D : public pzgeom::TPZNodeRep<2,pztopology::TPZLine> {
		
	public:
        typedef pztopology::TPZLine Top;
		/** @brief Number of nodes (connects) */
		enum {NNodes = 2};
		/** @brief It is not linear mapping */
                
                public:
int ClassId() const override;

        
        virtual void ParametricDomainNodeCoord(int64_t node, TPZVec<REAL> &nodeCoord);
        
		/** @brief Constructor */
		TPZEllipse3D(const TPZEllipse3D &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZEllipse3D::ClassId),pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap),
            fsAxeX(cp.fsAxeX), fsAxeY(cp.fsAxeY), fAxes(cp.fAxes), fOrigin(cp.fOrigin), fAngleIni(cp.fAngleIni), fAngleFinal(cp.fAngleFinal)
        {
		}
		/** @brief Default constructor */
		TPZEllipse3D() : TPZRegisterClassId(&TPZEllipse3D::ClassId),pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(),
        fsAxeX(0.), fsAxeY(0.), fAxes(), fOrigin(),fAngleIni(0.), fAngleFinal(0.)
        {
		}
		/** @brief Destructor */
		virtual ~TPZEllipse3D()
		{
		}
		/** @brief Copy constructor */
		TPZEllipse3D(const TPZEllipse3D &cp) : TPZRegisterClassId(&TPZEllipse3D::ClassId),pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp),
        fsAxeX(cp.fsAxeX), fsAxeY(cp.fsAxeY), fAxes(cp.fAxes), fOrigin(cp.fOrigin), fAngleIni(cp.fAngleIni), fAngleFinal(cp.fAngleFinal)
        {
		}
		/** @brief Copy constructor */		
		TPZEllipse3D(const TPZEllipse3D &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZEllipse3D::ClassId),pzgeom::TPZNodeRep<NNodes, pztopology::TPZLine>(cp),
        fsAxeX(cp.fsAxeX), fsAxeY(cp.fsAxeY), fAxes(cp.fAxes), fOrigin(cp.fOrigin), fAngleIni(cp.fAngleIni), fAngleFinal(cp.fAngleFinal)
        {
		}
		/** @brief Constructor with node indexes given */
		TPZEllipse3D(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZEllipse3D::ClassId),pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes),
        fsAxeX(0.), fsAxeY(0.), fAxes(), fOrigin(),fAngleIni(0.), fAngleFinal(0.)
        {
		}
        
        TPZEllipse3D &operator=(const TPZEllipse3D &cp)
        {
            pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>::operator=(cp);
            fsAxeX = cp.fsAxeX;
            fsAxeY = cp.fsAxeY;
            fAxes = cp.fAxes;
            fOrigin = cp.fOrigin;
            fAngleIni = cp.fAngleIni;
            fAngleFinal = cp.fAngleFinal;
            return *this;
        }
		
		/**
		 * @brief Adjust node coordinates in case of does not belong to the
		 * ellipse arc defined by the given origin and semiAxes
		 */
		virtual void AdjustNodesCoordinates(TPZGeoMesh &mesh);
		
		/**
		 * @brief Origin defines the translation of ellipse while semi-axes defines the rotation of ellipse, i.:
		 * @param Origin origin coordinates
		 * @param SemiAxeX is the vector that defines the direction of X axes of ellipse, whose value is its norm.
		 * @param SemiAxeY is the vector that defines the direction of Y axes of ellipse, whose value is its norm.
		 */
		void SetAxes(TPZVec<REAL> Origin, TPZVec<REAL> SemiAxeX, TPZVec<REAL> SemiAxeY, TPZGeoMesh &gmesh);
		
		/* brief compute the coordinate of a point given in parameter space */
//        template<class T>
//        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            X(coord,loc,result);
//        }
//
        template<class T>
        void GradX(TPZFMatrix<REAL> &cornerco, TPZVec<T> &par, TPZFMatrix<T> &gradx) const;
    
        
		/**
		 * Este metodo estabelece o mapeamento de um elemento linear de dominio [-1,1]
		 * para o arco de elipse definido pelos nohs inicial e final, os quais devem pertencer ao lugar geometrico da elipse,
		 * ========> SEGUINDO A REGRA DA MAO DIREITA (sentido anti-horario em relacao ao sistema local definido pelos vetores dos semi-eixos X e Y).
		 */
        template<class T>
		void X(TPZFMatrix<REAL> &nodeCoord,TPZVec<T> &qsi,TPZVec<T> &x) const;
		
		
        static bool IsLinearMapping(int side)
        {
            return false;
        }
        
		static std::string TypeName() { return "Ellipse3d";}
		// static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);
		
		/** @brief Creates a geometric element according to the type of the father element */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);
		
		void GetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix<REAL> &nodes);
        
		void SetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix<REAL> &nodes);
		
		double sAxeX()
		{
			return this->fsAxeX;
		}
		double sAxeY()
		{
			return this->fsAxeY;
		}
		TPZVec<REAL> Origin()
		{
			return this->fOrigin;
		}
		
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        

	private:
		
		
        /// Compute the angle of a coordinate as a function of the axes
        double ComputeAngle(TPZVec<REAL> &co) const;
        
		
		double fsAxeX;//norma de fSemiAxeX
		
		double fsAxeY;//norma de fSemiAxeY
        
        /// Rotation matrix where the axes correspond to rows of the matrix
        // for any point on the 3D ellips with coordinate co:
        // fAxes*(co-fOrigin) = coordinates in the local system (the third coordinate should be zero)
        TPZFNMatrix<9,REAL> fAxes;
		
		TPZManVector<REAL,3> fOrigin;
        
        double fAngleIni = 0.;
        double fAngleFinal = 0.;
        
	};
	
};

#endif
