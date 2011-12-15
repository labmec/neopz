//
//  pzcondensedcompel.cpp
//  PZ
//
//  Created by Philippe Devloo on 12/9/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzcondensedcompel.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcondensedcompel"));
#endif



TPZCondensedCompEl::~TPZCondensedCompEl()
{
    fMesh->ElementVec()[fIndex] = fReferenceCompEl;
    delete fReferenceCompEl;
}

/**
 * @brief Set the index i to node inode
 * @param inode node to set index
 * @param index index to be seted
 */
void TPZCondensedCompEl::SetConnectIndex(int inode, int index)
{
    LOGPZ_ERROR(logger,"SetConnectIndex should never be called")
    DebugStop();
}

/**
 * @brief Method for creating a copy of the element in a patch mesh
 * @param mesh Patch clone mesh
 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
 * @param gl2lcElMap map the computational elements
 */
/**
 * Otherwise of the previous clone function, this method don't
 * copy entire mesh. Therefore it needs to map the connect index
 * from the both meshes - original and patch
 */
TPZCompEl *TPZCondensedCompEl::ClonePatchEl(TPZCompMesh &mesh,
                                std::map<int,int> & gl2lcConMap,
                                std::map<int,int> & gl2lcElMap) const
{
    TPZCompEl *cel = fReferenceCompEl->ClonePatchEl(mesh,gl2lcConMap,gl2lcElMap);
    TPZCondensedCompEl *result = new TPZCondensedCompEl(cel);
    result->fIndexes = fIndexes;
    result->fCondensed = fCondensed;
    return result;
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param sol finite element solution
 * @param dsol solution derivatives
 * @param axes axes associated with the derivative of the solution
 */
void TPZCondensedCompEl::ComputeSolution(TPZVec<REAL> &qsi,
                             TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes)
{
    fReferenceCompEl->ComputeSolution(qsi,sol,dsol,axes);
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
 * This method will function for both volumetric and interface elements
 * @param qsi master element coordinate of the interface element
 * @param normal vector
 * @param leftsol finite element solution
 * @param dleftsol solution derivatives
 * @param leftaxes axes associated with the left solution
 * @param rightsol finite element solution
 * @param drightsol solution derivatives
 * @param rightaxes axes associated with the right solution
 */
void TPZCondensedCompEl::ComputeSolution(TPZVec<REAL> &qsi,
                             TPZVec<REAL> &normal,
                             TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                             TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes)
{
    fReferenceCompEl->ComputeSolution(qsi,normal,leftsol,dleftsol,leftaxes,rightsol,drightsol,rightaxes);
}

/**
 * @brief Computes solution and its derivatives in local coordinate qsi
 * @param qsi master element coordinate
 * @param phi matrix containing shape functions compute in qsi point
 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
 * @param axes [in] axes indicating the direction of the derivatives
 * @param sol finite element solution
 * @param dsol solution derivatives
 */
void TPZCondensedCompEl::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                             const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol)
{
    fReferenceCompEl->ComputeSolution(qsi,phi,dphix,axes,sol,dsol);
}

void TPZCondensedCompEl::Resequence()
{
    TPZStack<int> condensed;
    TPZStack<int> notcondensed;
    int nint=0,next=0;
    int ncon = NConnects();
    for (int i=0; i<ncon ; ++i) {
        TPZConnect &c = Connect(i);
        if (c.IsCondensed()) {
            if(c.HasDependency())
            {
                std::cout << "Not Implemented yet\n";
                DebugStop();
            }
            condensed.Push(i);
            nint += c.NDof();
        }
        else
        {
            notcondensed.Push(i);
            next += c.NDof();
        }
    }
    int ncond = condensed.size();
    for (int i=0; i<ncond; ++i) {
        fIndexes[i] = condensed[i];
    }
    for (int i=0; i<notcondensed.size(); ++i) {
        fIndexes[i+ncond] = notcondensed[i];
    }
    fCondensed.Redim(nint+next,nint);
}

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZCondensedCompEl::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    fReferenceCompEl->CalcStiff(ek,ef);
    ek.PermuteGather(fIndexes);
    ef.PermuteGather(fIndexes);
    int dim = ek.fMat.Rows();
    for (int i=0; i<dim ; ++i) {
        for (int j=0; j<dim ; ++j) {
            fCondensed(i,j) = ek.fMat(i,j);
        }
    }
//    fCondensed = ek.fMat;
    fCondensed.SetF(ef.fMat);
}


/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
void TPZCondensedCompEl::CalcResidual(TPZElementMatrix &ef)
{
    fReferenceCompEl->CalcResidual(ef);
    ef.PermuteGather(fIndexes);
    fCondensed.SetF(ef.fMat);
}


