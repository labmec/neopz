//
//  TRMPhaseInterfaceTransport.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPhaseInterfaceTransport.h"

/**
 * Empty Constructor
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport() : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>()
{
    
}

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh.
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport(int matid) : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>(matid)
{
    
}


/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport(const TRMPhaseInterfaceTransport &mat) : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>(mat)
{
    
}

/**
 * Destructor
 */
TRMPhaseInterfaceTransport::~TRMPhaseInterfaceTransport()
{
    DebugStop();
}

/** Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered
 * as necessary. Each derived class may optimize performance by selecting
 * only the necessary data.
 * @since April 10, 2007
 */
void TRMPhaseInterfaceTransport::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    DebugStop();
}

void TRMPhaseInterfaceTransport::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    DebugStop();
}



/** print out the data associated with the material */
void TRMPhaseInterfaceTransport::Print(std::ostream &out)
{
    DebugStop();
}


/** returns the variable index associated with the name */
int TRMPhaseInterfaceTransport::VariableIndex(const std::string &name)
{
    DebugStop();
}


/** returns the number of variables associated with the variable
 indexed by var.  var is obtained by calling VariableIndex */
int TRMPhaseInterfaceTransport::NSolutionVariables(int var)
{
    DebugStop();
}


/** Computes the divergence over the parametric space */
void TRMPhaseInterfaceTransport::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU)
{
    DebugStop();
}


/** returns the solution associated with the var index based on
 * the finite element approximation */
void TRMPhaseInterfaceTransport::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout)
{
    DebugStop();
}





// Contribute Methods being used

/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}


/**
 * It computes a contribution to the load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}



/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    DebugStop();
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
    DebugStop();
}


/**
 * Unique identifier for serialization purposes
 */
int TRMPhaseInterfaceTransport::ClassId() const
{
    DebugStop();
}


/**
 * Save the element data to a stream
 */
void TRMPhaseInterfaceTransport::Write(TPZStream &buf, int withclassid)
{
    DebugStop();
}


/**
 * Read the element data from a stream
 */
void TRMPhaseInterfaceTransport::Read(TPZStream &buf, void *context)
{
    DebugStop();
}



/// Copy the n+1 data to the n data
void TRMPhaseInterfaceTransport::UpdateMemory()
{
    DebugStop();
}

