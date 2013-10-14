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
		/** @brief Number of nodes (connects) */
		enum {NNodes = 2};
		/** @brief It is not linear mapping */
        
        virtual void ParametricDomainNodeCoord(long node, TPZVec<REAL> &nodeCoord);
        
		bool IsLinearMapping() const
		{
			return false;
		}
		/** @brief Constructor */
		TPZEllipse3D(const TPZEllipse3D &cp,std::map<long,long> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap){
		}
		/** @brief Default constructor */
		TPZEllipse3D() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(){
		}
		/** @brief Destructor */
		virtual ~TPZEllipse3D()
		{
		}
		/** @brief Copy constructor */
		TPZEllipse3D(const TPZEllipse3D &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp){
		}
		/** @brief Copy constructor */		
		TPZEllipse3D(const TPZEllipse3D &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
		}
		/** @brief Constructor with node indexes given */
		TPZEllipse3D(TPZVec<long> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes){
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
		void SetAxes(TPZVec<REAL> Origin, TPZVec<REAL> SemiAxeX, TPZVec<REAL> SemiAxeY);
		
		/* brief compute the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
        /* @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
		/**
		 * Este metodo estabelece o mapeamento de um elemento linear de dominio [-1,1]
		 * para o arco de elipse definido pelos nohs inicial e final, os quais devem pertencer ao lugar geometrico da elipse,
		 * ========> SEGUINDO A REGRA DA MAO DIREITA (sentido anti-horario em relacao ao sistema local definido pelos vetores dos semi-eixos X e Y).
		 */
		void X(TPZFMatrix<REAL> &nodeCoord,TPZVec<REAL> &qsi,TPZVec<REAL> &x) const;
		
		/** @brief Computes the jacobian matrix to X mapping */
		void Jacobian(TPZFMatrix<REAL> &nodeCoord, TPZVec<REAL> &qsi, TPZFMatrix<REAL> &jac, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv) const;
		
		static std::string TypeName() { return "Linear";}
		static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);
		
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
		
		void GetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix<REAL> &nodes);
		void SetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix<REAL> &nodes);
		
		TPZVec<REAL> SemiAxeX()
		{
			return this->fSemiAxeX;
		}
		double sAxeX()
		{
			return this->fsAxeX;
		}
		TPZVec<REAL> SemiAxeY()
		{
			return this->fSemiAxeY;
		}
		double sAxeY()
		{
			return this->fsAxeY;
		}
		TPZVec<REAL> Origin()
		{
			return this->fOrigin;
		}
		
	private:
		
		/**
		 * Este metodo retorna o angulo[ang(vini),ang(vfin)] que corresponde a um qsi[-1,1]
		 * @brief Returns the angle corresponding to qsi
		 */
		double Angle(double qsi, TPZFMatrix<REAL> &vini, TPZFMatrix<REAL> &vfin) const;
		
		/** @brief Returns the derivate of the angle with respect to qsi */
		double DAngleDqsi(TPZFMatrix<REAL> &vini, TPZFMatrix<REAL> &vfin) const;
		
		/** @brief Compute de equation of the elipse as function of the angle */
		TPZFMatrix<REAL> EllipseR2equation(double ang) const;
		
		/**
		 * Derivada da equacao em R2 da elipse em funcao com respeito ao angulo
		 */
		TPZFMatrix<REAL> DEllipseR2equationDang(double ang) const;
		
		/**
		 * fSemiAxeX @param - Half of Ellipse Axe in X direction
		 * fSemiAxeY @param - Half of Ellipse Axe in Y direction
		 */
		TPZVec<REAL> fSemiAxeX;
		double fsAxeX;//norma de fSemiAxeX
		
		TPZVec<REAL> fSemiAxeY;
		double fsAxeY;//norma de fSemiAxeY
		
		TPZVec<REAL> fOrigin;
	};
	
};

#endif
