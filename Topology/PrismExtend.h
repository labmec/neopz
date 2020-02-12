/**
 * @file
 * @brief Contains the Pr class(template) which defines the Prismatic extension of a topology.
 */

#ifndef TTPZGEOMPRISMEXTENDHPP
#define TTPZGEOMPRISMEXTENDHPP

#pragma once

#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "tpzprinteg.h"
#include "tpzpoint.h"
#include "pzaxestools.h"

namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @brief Defines the Prismatic extension of a topology. \ref topology "Topology"
	 */
	template<class TFather>
	class Pr :
	public TFather
	{
	public:
		
		/** @brief Enumerate for topological characteristics */
		enum {NCornerNodes = 2*TFather::NCornerNodes, NSides = 3*TFather::NSides, Dimension = TFather::Dimension+1};
		
		/** @brief Default constructor */
		Pr();
		
		/** @brief Default destructor */
		virtual ~Pr();
		
		/** @name About sides of the topological element
		 * @{ */
		
		/**
		 * @brief Returns the dimension of the side
		 */
		static int SideDimension(int side);
		
		/** @brief Get all sides with lower dimension on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
		/** @brief Get all sides with lower dimension but equal to DimTarget on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);
		
		/**
		 * @brief Returns all sides whose closure contains side
		 * @param side smaller dimension side
		 * @param high vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);
		/**
		 * @brief Returns the number of nodes (not connectivities) associated with a side
		 */
		static int NSideNodes(int side);
		/**
		 * @brief Returns the local node number of the node "node" along side "side"
		 */
		static int SideNodeLocId(int side, int node);
		/**
		 * @brief Returns the number of connects of the element (3)
		 */
		static int NumSides();
		
		/**
		 * @brief Returns the number of nodes (not connectivities) associated with a side
		 */
		static int NContainedSides(int side);
		/**
		 * @brief Returns the local connect number of the connect "c" along side "side"
		 */
		static int ContainedSideLocId(int side, int c);
		
		/** @} */
		
		/**
		 * @brief Returns the barycentric coordinates in the master element space of the original element
		 */
		static void CenterPoint(int side, TPZVec<REAL> &center);
		
		/** @name About Transformations
		 * @{ */
		
		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom ot
		 * the side sideto
		 * @param sidefrom side where the point resides
		 * @param sideto side whose closure contains sidefrom
		 */
		static TPZTransform<> SideToSideTransform(int sidefrom, int sideto);
		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side side from which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform<> object
		 * @see the class TPZTransform
		 */
		static TPZTransform<> TransformSideToElement(int side);
		/**
		 * @brief Returns the transformation which transform a point from the interior of the element to the side
		 * @param side side to which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform<> object
		 * @see the class TPZTransform
		 */
		static TPZTransform<> TransformElementToSide(int side);
		
		/**
		 * @brief Computes the linear map from an internal point to the parameter space of the side
		 * returns the jacobian of the transformation
		 */
		static void MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);
		
		/** @} */
		
		/** @name Methods related over numeric integration
		 * @{ */
		
		/**
		 * @brief Create an integration rule over side
		 * @param side Side to create integration rule
		 * @param order Order of the integration rule to be created
		 */
		static TPZIntPoints *CreateSideIntegrationRule(int side, int order);
		/** @brief Typedef to numerical integration rule */
		typedef typename TFather::IntruleType fatintruletype;
		/** @brief Typedef to graphical element type */
		typedef TPZPrInteg<fatintruletype> IntruleType;
		
		/** @} */
		
		/**
		 * @brief Returns the type of the element as specified in file pzeltype.h
		 */
		static std::string StrType();
		
		/**
		 * @brief Returns the type of the element as specified in file pzeltype.h
		 */
		static std::string StrType(int side);
		
		/** @brief Volume of the master element (measure) */
		static REAL RefElVolume();
		
		/**
		 * @brief Uses log4cxx to print the results of all methods
		 */
		static void Diagnostic();
	};
	
}

#endif
