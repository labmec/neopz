/*
 *  tpzellipse3d.h
 *  NeoPZ
 *
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

using namespace std;
using namespace pzgeom;
using namespace pztopology;

class TPZEllipse3D : public TPZNodeRep<2,TPZLine> {

public:

	enum {NNodes = 2};

    bool IsLinearMapping() const
	{
		return false;
	}

    TPZEllipse3D(const TPZEllipse3D &cp,std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap){
    }

    TPZEllipse3D() : TPZNodeRep<NNodes,pztopology::TPZLine>(){
    }

    virtual ~TPZEllipse3D()
    {

    }

    TPZEllipse3D(const TPZEllipse3D &cp) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp){
    }

    TPZEllipse3D(const TPZEllipse3D &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
    }

	TPZEllipse3D(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes){
	}

    /**
	 * Adjust node coordinates in case of does not belong to the
	 * ellipse arc defined by the given origin and semiAxes
     */
	virtual void AdjustNodesCoordinates(TPZGeoMesh &mesh);

    /**
     * Origin defines the translation of ellipse while semi-axes defines the rotation of ellipse, i.:
     * - SemiAxeX is the vector that defines the direction of X axes of ellipse, whose value is its norm.
	 * - SemiAxeY is the vector that defines the direction of Y axes of ellipse, whose value is its norm.
	 */
	void SetAxes(TPZVec<REAL> Origin, TPZVec<REAL> SemiAxeX, TPZVec<REAL> SemiAxeY);

	 /**
      * Este metodo estabelece o mapeamento de um elemento linear de dominio [-1,1]
      * para o arco de elipse definido pelos nohs inicial e final, os quais devem pertencer ao lugar geometrico da elipse,
      * ========> SEGUINDO A REGRA DA MAO DIREITA (sentido anti-horario em relacao ao sistema local definido pelos vetores dos semi-eixos X e Y).
      */
    void X(TPZFMatrix &nodeCoord,TPZVec<REAL> &qsi,TPZVec<REAL> &x);

     /**
      * Este metodo estabelece a matriz jacobiana do mapeamento X (acima declarado)
      */
    void Jacobian(TPZFMatrix &nodeCoord, TPZVec<REAL> &qsi, TPZFMatrix &jac, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);

    static std::string TypeName() { return "Linear";}
    static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);

	/**
	 * Creates a geometric element according to the type of the father element
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
    */
    double Angle(double qsi, TPZFMatrix &vini, TPZFMatrix &vfin);

   /**
    * Este metodo retorna a derivada do angulo[ang(vini),ang(vfin)] com respeito a qsi[-1,1]
    */
    double DAngleDqsi(TPZFMatrix &vini, TPZFMatrix &vfin);

   /**
    * Equacao em R2 da elipse em funcao do angulo[ 0 , 2pi]
    */
    TPZFMatrix EllipseR2equation(double ang);

    /**
     * Derivada da equacao em R2 da elipse em funcao com respeito ao angulo
     */
    TPZFMatrix DEllipseR2equationDang(double ang);

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

#endif
