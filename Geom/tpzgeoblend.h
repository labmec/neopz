/**
 * @file
 * @brief Contains the TPZGeoBlend class which implements a blending map from curved boundaries to the interior of the element.
 */

#ifndef TPZGEOBLEND_H
#define TPZGEOBLEND_H

#include "pzgeotriangle.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeoel.h"
#include "pznoderep.h"
#include "pzgeoelside.h"

#include <iostream>

namespace pzgeom 
{
	
	/**
	 * @ingroup geometry
	 * @brief Implements a blending map from curved boundaries to the interior of the element. \ref geometry "Geometry"
	 * @author Paulo Cesar de Alvarenga Lucci
	 * @since 2007
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
		
		/** @brief Constructor with list of nodes */
		TPZGeoBlend(TPZVec<long> &nodeindexes) : TGeo(nodeindexes) {
		}
		/** @brief Empty constructor */
		TPZGeoBlend() : TGeo() {
		}
		
		/** @brief Constructor with node map */
		TPZGeoBlend(const TPZGeoBlend &cp,std::map<long,long> & gl2lcNdMap) : TGeo(cp,gl2lcNdMap) {
		}
		/** @brief Copy constructor */
		TPZGeoBlend(const TPZGeoBlend &cp) : TGeo(cp) {
		}
		/** @brief Copy constructor */
		TPZGeoBlend(const TPZGeoBlend &cp, TPZGeoMesh &) : TGeo(cp) {
		}
        
        void Read(TPZStream &buf,void *context)
        {
            TGeo::Read(buf,context);
            for (int is=0; is < 1+TGeo::NSides - TGeo::NNodes; is++) {
                fNeighbours[is].Read(buf);
            }
            for (int is=0; is < 1+TGeo::NSides - TGeo::NNodes; is++) {
                fTrans[is].Read(buf);
            }
        }
        
        void Write(TPZStream &buf)
        {
            TGeo::Write(buf);
            for (int is=0; is < 1+TGeo::NSides - TGeo::NNodes; is++) {
                fNeighbours[is].Write(buf);
            }
            for (int is=0; is < 1+TGeo::NSides - TGeo::NNodes; is++) {
                fTrans[is].Write(buf);
            }
		}

		void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform &trans);
		
		TPZGeoElSide Neighbour(int side,TPZGeoMesh *gmesh) const {
            if (side < TGeo::NNodes) {
                DebugStop();
            }
			return TPZGeoElSide(fNeighbours[side-TGeo::NNodes],gmesh);
		}
		
		TPZTransform TransfBetweenNeigh(int side) const {
			return fTrans[side - TGeo::NNodes];
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return TGeo::TypeName();} 
		/** @brief Get the coordinates of the point at geometric elements from coordinates of the parametric point at the master element */
		void X(const TPZGeoEl &gel, TPZVec<REAL>& par, TPZVec<REAL> &result) const;
		/** @brief Computes the Jacobian for parametric point at master element */
		void Jacobian(const TPZGeoEl &gel, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;
		
		/** @brief Print all relevant data of the element to cout*/
		void Print(std::ostream & out = std::cout) const;
		
		/**
		 * @brief Initialize the element checking connectivity on all sides.
		 *
		 / Este método percorre todos os lados do elemento fornecido (el) e, para cada lado,
		 / percorre todos os vizinhos até retornar a si próprio. Se neste processo encontrar algum elemento que
		 / seja quadrático e não seja do tipo TPZGeoBlend, recorre ao método SetNeighbourInfo() para estabelecer que este
		 / elemento encontrado será seu vizinho pelo respectivo lado.
		 */
		void Initialize(TPZGeoEl *refel);
		
		//void Initialize(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh);
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);
		
		TPZGeoEl *CreateBCGeoBlendEl(TPZGeoEl *orig,int side,int bc);
		
		TPZGeoEl *CreateGeoBlend(TPZGeoMesh &mesh, MElementType type, TPZVec<long>& nodeindexes, int matid, long& index);
		
        
		
	public:
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
		
	protected:
		
		/// Project the InternalPar parameter to the parameter of the neighbour along side. Return true if the map is nonsingular
		bool MapToNeighSide(int side, int sidedim, TPZVec<REAL> &InternalPar, TPZVec<REAL> &NeighPar, TPZFMatrix<REAL> &JacNeighSide) const;
		/** @brief Vector of indexes of the neighbours */
		TPZGeoElSideIndex fNeighbours[1+TGeo::NSides - TGeo::NNodes];
		TPZTransform fTrans[1+TGeo::NSides - TGeo::NNodes];
	};
	
};
#endif
