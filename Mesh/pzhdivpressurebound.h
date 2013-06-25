//
//  pzhdivpressurebound.h
//  PZ
//
//  Created by Agnaldo Farias on 25/06/13.
//
//

/**
 * @file
 * @brief Contains declaration of TPZCompElHDivPressureBound class which implements a generic computational element (HDiv-Pressure scope variant).
 */


#ifndef __PZ__pzhdivpressurebound__
#define __PZ__pzhdivpressurebound__

#include <iostream>

#include "pzelchdivbound2.h"

/** \addtogroup CompElement */
/** @{ */
/**
 * @brief Implements a generic computational element to HDiv-Pressure scope. \ref CompElement "Computational Element"
 * @author Agnaldo Farias
 * @since Jun 25, 2013.
 */
/**
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivPressureBound : public TPZCompElHDivBound2<TSHAPE> {

    /** @brief Defines the interpolation order for pressure variable*/
	int fPressureOrder;
    
public:
    
    TPZCompElHDivPressureBound(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);
    
    /**
	 * @brief Constructor used to generate patch mesh...\n
	 * Generates a map of connect index from global mesh to clone mesh
	 */
	TPZCompElHDivPressureBound(TPZCompMesh &mesh, const TPZCompElHDivBound2<TSHAPE> &copy, std::map<int,int> & gl2lcConMap, std::map<int,int> & gl2lcElMap);
    
    
    TPZCompElHDivPressureBound(TPZCompMesh &mesh, const TPZCompElHDivPressureBound<TSHAPE> &copy);
    
    /** @brief Default constructor */
	TPZCompElHDivPressureBound();
    
	/** @brief Default destructor */
	virtual ~TPZCompElHDivPressureBound();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZCompElHDivPressureBound<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> & gl2lcConMap,std::map<int,int>&gl2lcElMap) const
	{
		return new TPZCompElHDivPressureBound<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}

    
    
    virtual int NConnects() const;

};

#endif /* defined(__PZ__pzhdivpressurebound__) */
