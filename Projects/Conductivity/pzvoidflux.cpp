//
//  pzvoidflux.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/3/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzvoidflux.h"
#include "pzbndcond.h"

/** Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered
 * as necessary. Each derived class may optimize performance by selecting
 * only the necessary data.
 * @since April 10, 2007
 */
void TPZVoidFlux::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(true);
    data.fNeedsNeighborSol = false;
    data.fNeedsNeighborCenter = false;
    data.fNeedsNormal = false;

}

/** Fill material data parameter with necessary requirements for the
 * ContributeInterface method. Here, in base class, all requirements are considered
 * as necessary. Each derived class may optimize performance by selecting
 * only the necessary data.
 * @since April 10, 2007
 */
void TPZVoidFlux::FillDataRequirementsInterface(TPZMaterialData &data)
{
    data.fNeedsSol = false;
    data.fNeedsNeighborSol = false;
    data.fNeedsHSize = true;
    data.fNeedsNeighborCenter = true;
    data.fNeedsNormal = false;
    if(fLinearContext == false){
        data.fNeedsNeighborSol = true;
    }
}

/// returns the number of state variables associated with the material
int TPZVoidFlux::NStateVariables() 
{
    return 1;
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZVoidFlux::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef)
{
    
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
void TPZVoidFlux::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc)
{
    std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
    DebugStop();
}

/**
 * It computes a contribution to the residual vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the residual vector
 * @since April 16, 2007
 */
void TPZVoidFlux::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef)
{
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
void TPZVoidFlux::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ef, TPZBndCond &bc)
{
    std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
    DebugStop();    
}




/// computes a contribution to stiffness matrix and load vector at one integration point
/**
 * @param data [in] all data needed to compute the stiffness matrix
 * @param weight [in] weight of the integration point
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZVoidFlux::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef)
{
    int nleft = data.phil.Rows();
    int nright = data.phir.Rows();
    if (nleft != 1 || nright != 1) {
        DebugStop();
    }
    int il,jl,ir,jr;
    for (il=0; il<nleft; il++) {
        for (jl=0; jl<nleft; jl++) {
            ek(il,jl) += weight*fConductivity*data.phil(il,0)*data.phil(jl,0);
        }
    }
    for (il=0; il<nleft; il++) {
        for (jr=0; jr<nright; jr++) {
            ek(il,nleft+jr) -= weight*fConductivity*data.phil(il,0)*data.phir(jr,0);
        }
    }
    for (ir=0; ir<nright; ir++) {
        for (jl=0; jl<nleft; jl++) {
            ek(ir+nleft,jl) -= weight*fConductivity*data.phir(ir,0)*data.phil(jl,0);
        }
    }
    for (ir=0; ir<nright; ir++) {
        for (jr=0; jr<nright; jr++) {
            ek(ir+nleft,jr+nleft) += weight*fConductivity*data.phir(ir,0)*data.phir(jr,0);
        }
    }
}

/**
 * It computes a contribution to residual vector at one integration point
 * @param data [in]
 * @param weight [in]
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZVoidFlux::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef)
{
    std::cout << __PRETTY_FUNCTION__ << " not implemented yet\n";
    DebugStop();
}

/**
 * It computes a contribution to stiffness matrix and load vector at one BC integration point
 * @param data [in]
 * @param weight [in]
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition object
 * @since April 16, 2007
 */
void TPZVoidFlux::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc)
{
    int nleft = data.phil.Rows();
    if (nleft != 1) {
        DebugStop();
    }
    switch (bc.Type()) {
            /// Dirichlet pressure boundary
        case 0:
        {
            int il,jl;
            for (il=0; il<nleft; il++) {
                for (jl=0; jl<nleft; jl++) {
                    ek(il,jl) += weight*fConductivity*data.phil(il,0)*data.phil(jl,0);
                }
                ef(il,0) += weight*fConductivity*bc.Val2()(0,0);
            }
        }
            break;
        case 1:
        {
            int il;
            for (il=0; il<nleft; il++) {
                ef(il,0) += weight*bc.Val2()(0,0);
            }
        }
        default:
            std::cout << __PRETTY_FUNCTION__ << "boundary condition of unknown type " << bc.Type() << std::endl;
            DebugStop();
            break;
    }
}

/**
 * It computes a contribution to residual vector at one BC integration point
 * @param data [in]
 * @param weight [in]
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition object
 * @since April 16, 2007
 */
void TPZVoidFlux::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef,TPZBndCond &bc)
{
    std::cout << __PRETTY_FUNCTION__ << " not implemented yet\n";
    DebugStop();
}

/**
 * Dicontinuous galerkin materials implement contribution of discontinuous elements and interfaces.
 * Interfaces may be conservative or not conservative. It is important to agglomeration techniques
 * when using multigrid pre-conditioner. Conservative interfaces into agglomerate elements do not
 * need to be computed. However non-conservative interfaces must be computed in all multigrid levels.
 * Default is non-conservative, because of the computation of a conservative interface into an agglomerate
 * does not ruin the solution.
 * @since Feb 05, 2004
 */
int TPZVoidFlux::IsInterfaceConservative()
{
    return 1;
}


/** print out the data associated with the material*/
void TPZVoidFlux::Print(std::ostream &out)
{
    TPZDiscontinuousGalerkin::Print(out);
    out << "Material conductivity " << fConductivity << std::endl;
}

/**returns the variable index associated with the name*/
int TPZVoidFlux::VariableIndex(const std::string &name)
{
    return TPZDiscontinuousGalerkin::VariableIndex(name);
}

/** returns the number of variables associated with the variable
 indexed by var.  var is obtained by calling VariableIndex*/
int TPZVoidFlux::NSolutionVariables(int var)
{
    return TPZDiscontinuousGalerkin::NSolutionVariables(var);
}

/**returns the solution associated with the var index based on
 * the finite element approximation*/
void TPZVoidFlux::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
    return TPZDiscontinuousGalerkin::Solution(data, var , Solout);
}

/**
 * Unique identifier for serialization purposes
 */
int TPZVoidFlux::ClassId() const
{
    return TPZMATERIALVOIDFLUX;
}

template class TPZRestoreClass<TPZVoidFlux,TPZMATERIALVOIDFLUX>;

/**
 * Save the element data to a stream
 */
void TPZVoidFlux::Write(TPZStream &buf, int withclassid)
{
    TPZDiscontinuousGalerkin::Write(buf,withclassid);
    buf.Write(&fConductivity);
}

/**
 * Read the element data from a stream
 */
void TPZVoidFlux::Read(TPZStream &buf, void *context)
{
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fConductivity);
}

/// create another material of the same type
TPZAutoPointer<TPZMaterial> TPZVoidFlux::NewMaterial()
{
    return new TPZVoidFlux(*this);
}

/// Read data of the material from a istream (file data)
void TPZVoidFlux::SetData(std::istream &data)
{
    TPZDiscontinuousGalerkin::SetData(data);
    data >> fConductivity;
}


