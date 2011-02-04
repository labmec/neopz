#ifndef TPZGEOBLEND_H
#define TPZGEOBLEND_H

#include "pzgeotriangle.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

template <class TGeo>
class TPZGeoBlend : public TGeo {

public:

     bool IsLinearMapping() const { 
          return false; 
     }
     bool IsGeoBlendEl() const { 
          return true; 
     }

     TPZGeoBlend(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TGeo(nodeindexes,mesh) {
     }

     TPZGeoBlend() : TGeo() {
     }

     TPZGeoBlend(const TPZGeoBlend &cp,std::map<int,int> & gl2lcNdMap) : TGeo(cp,gl2lcNdMap) {
     }

     TPZGeoBlend(const TPZGeoBlend &cp) : TGeo(cp) {
     }

     TPZGeoBlend(const TPZGeoBlend &cp, TPZGeoMesh &) : TGeo(cp) {
     }

     void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform &trans);

     TPZGeoElSide Neighbour(int side) {
          return fNeighbours[side-TGeo::NNodes];
     }

     TPZTransform TransfBetweenNeigh(int side) {
          return fTrans[side - TGeo::NNodes];
     }

     /**
      * returns the type name of the element
      */
    static std::string TypeName() { return TGeo::TypeName();} 

    void X(TPZFMatrix & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);

    void Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

    /**Print all relevant data of the element to cout*/
    void Print(std::ostream & out = std::cout);

    /**
    / Este método percorre todos os lados do elemento fornecido (el) e, para cada lado,
    / percorre todos os vizinhos até retornar a si próprio. Se neste processo encontrar algum elemento que
    / seja quadrático e não seja do tipo TPZGeoBlend, recorre ao método SetNeighbourInfo() para estabelecer que este
    / elemento encontrado será seu vizinho pelo respectivo lado.
   */
    void Initialize(TPZGeoEl *refel);

    void Initialize(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh);

   /**
    * Method which creates a geometric boundary condition 
    * element based on the current geometric element, 
    * a side and a boundary condition number
    */
  TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);

  TPZGeoEl *CreateBCGeoBlendEl(TPZGeoEl *orig,int side,int bc);

  TPZGeoEl *CreateGeoBlend(TPZGeoMesh &mesh, MElementType type, TPZVec<int>& nodeindexes, int matid, int& index);

	
public:
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<int>& nodeindexes,
									  int matid,
									  int& index);
	
	

protected:

	/// Project the InternalPar parameter to the parameter of the neighbour along side. Return true if the map is nonsingular
    bool MapToNeighSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &NeighPar, TPZFMatrix &JacNeighSide);
    TPZGeoElSide fNeighbours[1+TGeo::NSides - TGeo::NNodes];
    TPZTransform fTrans[1+TGeo::NSides - TGeo::NNodes];
};

#endif
