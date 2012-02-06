/**
 * @file
 * @brief Contains the functions to create different computational elements (one- two- three-dimensional).
 */
#ifndef CREATECONTINUOUSHPP
#define CREATECONTINUOUSHPP
/*
 *  pzcreatecontinuous.h
 *  NeoPZ
 *
 *  Created by Philippe Devloo on 20/11/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */

class TPZCompEl;
class TPZCompMesh;
class TPZGeoEl;
class TPZCompEl;
class TPZCompMesh;
#include <set>
#include "pzvec.h"
//#include "pzcmesh.h"

typedef TPZCompEl *(*TCreateFunction)(TPZGeoEl *el,TPZCompMesh &mesh,int &index);
/*
 * @brief Administer the creation of approximation spaces
 */
class TPZCreateApproximationSpace
{
    /**
     * @brief Function pointer which determines what type of computational element will be created
     */
    TPZCompEl *(*fp[8])(TPZGeoEl *el,TPZCompMesh &mesh,int &index);

public:
    
    TPZCreateApproximationSpace()
    {
        SetAllCreateFunctionsContinuous();
    }
    
    /*
     * @brief Create discontinuous approximation spaces
     */
    void SetAllCreateFunctionsDiscontinuous();
    /*
     * @brief Create continuous approximation spaces
     */
	void SetAllCreateFunctionsContinuous();
    /*
     * @brief Create a discontinuous approximation space with referred elements
     */
	void SetAllCreateFunctionsDiscontinuousReferred();
    /*
     * @brief Create a continuous approximation space with referred elements
     */
	void SetAllCreateFunctionsContinuousReferred();
    /*
     * @brief Create an approximation space with HDiv elements
     */
	void SetAllCreateFunctionsHDiv();
    /*
     * @brief Create approximation spaces corresponding to the space defined by cel
     */
	void SetAllCreateFunctions(TPZCompEl &cel, TPZCompMesh *mesh);
    /*
     * @brief Create an approximation space based on multiphysics elements
     */
	void SetAllCreateFunctionsMultiphysicElem();
    /*
     * @brief Create an approximation space with continous elements with memory
     */
    void SetAllCreateFunctionsContinuousWithMem();
    
    /**
     * @brief Set custom function pointers
     */
    void SetCreateFunctions(TPZVec<TCreateFunction> &createfuncs);
    
    /*
     * @brief Create a computational element using the function pointer for the topology
     */
    TPZCompEl *CreateCompEl(TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
    
	/**
	 * @brief Creates the computational elements, and the degree of freedom nodes
	 */ 
	/** Only element of material id in the set<int> will be created */
	static void BuildMesh(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs);
	
	/** @brief Creates the computational elements, and the degree of freedom nodes */
	static void AutoBuild(TPZCompMesh &cmesh);
    
	/** @brief Creates the computational elements, and the degree of freedom nodes */
	/**
	 * Elements created may be TPZInterpolatedElement or TPZCompElDisc. \n
	 * indices contains the type of the element. Element type are given by the enumerate MCreationType.
	 */
	static void AutoBuildContDisc(const TPZVec<TPZGeoEl*> &continuous, const TPZVec<TPZGeoEl*> &discontinuous);
    
    /**
     * @brief Encapsulate the elements in condensed computational elements
     */
    static void CondenseLocalEquations(TPZCompMesh &cmesh);
    
    /**
     * @brief Undo the encapsulate elements
     */
    static void UndoCondenseLocalEquations(TPZCompMesh &cmesh);
    
    /**
     * @brief transform in low order Raviar Tomas
     */
    static void MakeRaviartTomas(TPZCompMesh &cmesh);
    
    /**
     * @brief transform in low order Raviar Tomas
     */
    static void UndoMakeRaviartTomas(TPZCompMesh &cmesh);
    
    
};

/** @brief Creates computational point element */
TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational linear element */
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational quadrilateral element */
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational triangular element */
TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational cube element */
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational prismal element */
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational pyramidal element */
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational tetrahedral element */
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif