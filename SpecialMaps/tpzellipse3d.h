/**
 * @file
 * @brief Contains the TPZEllipse3D class which defines a linear geometric element which maps a line segment to an ellipse.
 */
/*
 *  Created by caju on 03/aug/09.
 *  Copyright 2009 LabMeC. All rights reserved.
 *
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
		
		enum {NNodes = 2};
		
		bool IsLinearMapping() const
		{
			return false;
		}
		
		TPZEllipse3D(const TPZEllipse3D &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap){
		}
		
		TPZEllipse3D() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(){
		}
		
		virtual ~TPZEllipse3D()
		{
			
		}
		
		TPZEllipse3D(const TPZEllipse3D &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp){
		}
		
		TPZEllipse3D(const TPZEllipse3D &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
		}
		
		TPZEllipse3D(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes){
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
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) const
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
		void X(TPZFMatrix &nodeCoord,TPZVec<REAL> &qsi,TPZVec<REAL> &x) const;
		
		/**
		 * Este metodo estabelece a matriz jacobiana do mapeamento X (acima declarado)
		 * @brief Computes the jacobian matrix to X mapping
		 */
		void Jacobian(TPZFMatrix &nodeCoord, TPZVec<REAL> &qsi, TPZFMatrix &jac, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv) const;
		
		static std::string TypeName() { return "Linear";}
		static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);
		
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid,
										  int& index);
		
		void GetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix &nodes);
		void SetNodesCoords(TPZGeoMesh &mesh, TPZFMatrix &nodes);
		
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
		
	private://metodos utilizados apenas pelos metodos publicos desta classe!
		
		/**
		 * Este metodo retorna o angulo[ang(vini),ang(vfin)] que corresponde a um qsi[-1,1]
		 * @brief Returns the angle corresponding to qsi
		 */
		double Angle(double qsi, TPZFMatrix &vini, TPZFMatrix &vfin) const;
		
		/**
		 * Este metodo retorna a derivada do angulo[ang(vini),ang(vfin)] com respeito a qsi[-1,1]
		 * @brief Returns the derivate of the angle with respect to qsi
		 */
		double DAngleDqsi(TPZFMatrix &vini, TPZFMatrix &vfin) const;
		
		/**
		 * Equacao em R2 da elipse em funcao do angulo[ 0 , 2pi]
		 * @brief Compute de equation of the elipse as function of the angle
		 */
		TPZFMatrix EllipseR2equation(double ang) const;
		
		/**
		 * Derivada da equacao em R2 da elipse em funcao com respeito ao angulo
		 */
		TPZFMatrix DEllipseR2equationDang(double ang) const;
		
	private:
		
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
