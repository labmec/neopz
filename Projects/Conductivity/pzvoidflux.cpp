//
//  pzvoidflux.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/3/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzvoidflux.h"
#include "pzbndcond.h"

REAL fluxfunction(REAL ratio)
{
    REAL expval = log(ratio);
    REAL logvalue = -3.04-2.21*expval-0.51*expval*expval-0.0597*expval*expval*expval-0.002638*expval*expval*expval*expval;
    REAL result = exp(logvalue);
    return result;
}

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
    data.fNeedsNormal = true;
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
void TPZVoidFlux::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef)
{
    int nleft = dataleft.phi.Rows();
    int nright = dataright.phi.Rows();
    if (nleft != 1 || nright != 1) {
        std::cout << __PRETTY_FUNCTION__ << " works only for constant shape functions\n";
        DebugStop();
    }
    TPZManVector<REAL,3> delleft(3), delright(3);
    for (int i=0; i<3; i++) {
        delleft[i] = data.x[i]-dataleft.XCenter[i];
        delright[i] = data.x[i] - dataright.XCenter[i];
    }
    REAL distleft(0.),distright(0.);
    for (int i=0; i<3; i++) {
        distleft += delleft[i]*data.normal[i];
        distright += delright[i]*data.normal[i];
    }
    REAL sumdist = fabs(distleft)+fabs(distright);
    REAL ratio = fBridgeSize/sumdist;
    REAL mult = fluxfunction(ratio);
    int il,jl,ir,jr;
    for (il=0; il<nleft; il++) {
        for (jl=0; jl<nleft; jl++) {
            ek(il,jl) += weight*fConductivity*mult*dataleft.phi(il,0)*dataleft.phi(jl,0);
        }
    }
    for (il=0; il<nleft; il++) {
        for (jr=0; jr<nright; jr++) {
            ek(il,nleft+jr) -= weight*fConductivity*mult*dataleft.phi(il,0)*dataright.phi(jr,0);
        }
    }
    for (ir=0; ir<nright; ir++) {
        for (jl=0; jl<nleft; jl++) {
            ek(ir+nleft,jl) -= weight*fConductivity*mult*dataright.phi(ir,0)*dataleft.phi(jl,0);
        }
    }
    for (ir=0; ir<nright; ir++) {
        for (jr=0; jr<nright; jr++) {
            ek(ir+nleft,jr+nleft) += weight*fConductivity*mult*dataright.phi(ir,0)*dataright.phi(jr,0);
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
void TPZVoidFlux::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix &ef)
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
void TPZVoidFlux::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc)
{
    int nleft = dataleft.phi.Rows();
    if (nleft != 1) {
        DebugStop();
    }
    TPZManVector<REAL,3> delleft(3);
    for (int i=0; i<3; i++) {
        delleft[i] = data.x[i]-dataleft.XCenter[i];
    }
    REAL distleft(0.);
    for (int i=0; i<3; i++) {
        distleft += delleft[i]*data.normal[i];
    }
    REAL sumdist = fabs(distleft);
    REAL ratio = fBridgeSize/sumdist/2.;
    REAL mult = fluxfunction(ratio);

    switch (bc.Type()) {
            /// Dirichlet pressure boundary
        case 0:
        {
            int il,jl;
            for (il=0; il<nleft; il++) {
                for (jl=0; jl<nleft; jl++) {
                    ek(il,jl) += weight*fConductivity*mult*dataleft.phi(il,0)*dataleft.phi(jl,0);
                }
                ef(il,0) += weight*fConductivity*mult*bc.Val2()(0,0);
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
void TPZVoidFlux::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ef,TPZBndCond &bc)
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
void TPZVoidFlux::SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<REAL> &Solout)
{
    return TPZDiscontinuousGalerkin::SolutionDisc(data, dataleft, dataright, var , Solout);
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
    buf.Write(&fBridgeSize);
}

/**
 * Read the element data from a stream
 */
void TPZVoidFlux::Read(TPZStream &buf, void *context)
{
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fConductivity);
    buf.Read(&fBridgeSize);
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
    data >> fBridgeSize;
}


