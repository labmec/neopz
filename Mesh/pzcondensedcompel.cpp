//
//  pzcondensedcompel.cpp
//  PZ
//
//  Created by Philippe Devloo on 12/9/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzcondensedcompel.h"
#include "pzlog.h"
#include "pzstepsolver.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcondensedcompel"));
#endif

TPZCondensedCompEl::TPZCondensedCompEl(TPZCompEl *ref)
{
    if(!ref)
    {
        DebugStop();
    }
    fReferenceCompEl = ref;
    fMesh = ref->Mesh();
    SetReference(ref->Reference()->Index());
    SetIndex(ref->Index());
    fMesh->ElementVec()[fIndex] = this;
    int ncon = ref->NConnects();
    fIndexes.resize(ncon);
    for (int ic=0; ic< ncon ; ic++) {
        fIndexes[ic] = ic;
    }
    Resequence();
}


TPZCondensedCompEl::~TPZCondensedCompEl()
{
    if (fMesh->ElementVec()[fIndex] == this) {
        fMesh->ElementVec()[fIndex] = fReferenceCompEl;
        delete fReferenceCompEl;
    }
}

/**
 * @brief create a copy of the condensed computational element in the other mesh
 */
TPZCondensedCompEl::TPZCondensedCompEl(const TPZCondensedCompEl &copy, TPZCompMesh &mesh)
{
    TPZCompEl *ref = fReferenceCompEl->Clone(mesh);
    if(!ref)
    {
        DebugStop();
    }
    fReferenceCompEl = ref;
    fMesh = ref->Mesh();
    SetReference(ref->Reference()->Index());
    SetIndex(ref->Index());
    fMesh->ElementVec()[fIndex] = this;
    int ncon = ref->NConnects();
    fIndexes.resize(ncon);
    for (int i=0; i<ncon ; ++i) {
        fIndexes[i] = i;
    }
}


/**
 * @brief unwrap the condensed element from the computational element and delete the condensed element
 */
void TPZCondensedCompEl::Unwrap()
{
    fMesh->ElementVec()[fIndex] = fReferenceCompEl;
    int ncon = NConnects();
    for (int ic=0; ic<ncon ; ic++) {
        Connect(ic).SetCondensed(false);
    }
    delete this;
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
                             TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes)
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
                             TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
                             TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes)
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
void TPZCondensedCompEl::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                             const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol)
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
        if(c.NElConnected() == 1)
        {
            c.SetCondensed(true);
        }
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
    TPZAutoPointer<TPZMatrix<REAL> > k00 = new TPZFMatrix<REAL>(nint,nint,0.);
    TPZStepSolver<REAL> *step = new TPZStepSolver<REAL>(k00);
    step->SetDirect(ECholesky);
    fCondensed.SetSolver(step);
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
    fCondensed.Zero();
    int dim = ek.fMat.Rows();
    for (int i=0; i<dim ; ++i) {
        for (int j=0; j<dim ; ++j) {
            fCondensed(i,j) = ek.fMat(i,j);
        }
    }
//    fCondensed = ek.fMat;
    fCondensed.SetF(ef.fMat);
    const TPZFMatrix<REAL> &k11 = fCondensed.K11Red();
    const TPZFMatrix<REAL> &f1 = fCondensed.F1Red();
    int dim0 = dim-k11.Rows();
    for (int i=dim0; i<dim; i++) {
        ef.fMat(i,0) = f1.GetVal(i-dim0,0);
        for (int j=dim0; j<dim; j++) {
            ek.fMat(i,j) = k11.GetVal(i-dim0,j-dim0);
        }
    }
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
    const TPZFMatrix<REAL> &f1 = fCondensed.F1Red();
    int dim1 = f1.Rows();
    int dim = ef.fMat.Rows();
    int dim0 = dim-dim1;
    for (int i= dim0; i<dim; i++) {
        ef.fMat(i,0) = f1.GetVal(i-dim0,0);
    }
}

/**
 * @brief Prints element data
 * @param out Indicates the device where the data will be printed
 */
void TPZCondensedCompEl::Print(std::ostream &out) const
{
    out << "Output for a condensed element\n";
    TPZCompEl::Print(out);
    out << "Internal index resequencing " << fIndexes << std::endl;
    fCondensed.Print("Condensed matrix",out);
}


/**
 * @brief Loads the solution within the internal data structure of the element
 */ 
/** Is used to initialize the solution of connect objects with dependency
 * Is also used to load the solution within SuperElements
 */
void TPZCondensedCompEl::LoadSolution()
{
    // initialize the solution of the constrained connects
    TPZCompEl::LoadSolution();
    // compute the solution of the internal equations
    int dim0=0, dim1=0;
    int nc = NConnects(),nc0 = 0, nc1 = 0;
    int ic;
    for (ic=0; ic<nc ; ic++) {
        TPZConnect &con = Connect(ic);
        int sz = con.NShape()*con.NState();
        if (con.IsCondensed()) {
#ifdef DEBUG
            if (dim1) {
                DebugStop();
            }
#endif
            dim0 += sz;
            nc0++;
        }
        else
        {
            dim1 += sz;
            nc1++;
        }
    }
    TPZBlock<REAL> &bl = Mesh()->Block();
    int count = 0;
    TPZFMatrix<REAL> u1(dim1,1,0.);
    TPZFMatrix<REAL> elsol(dim0+dim1,1,0.);
    for (ic=nc0; ic<nc ; ic++) {
        TPZConnect &c = Connect(ic);
        int seqnum = c.SequenceNumber();
        int blsize = bl.Size(seqnum);
        for (int ibl=0; ibl<blsize; ibl++) {
            u1(count++,0) = bl(seqnum,ibl,0,0);
        }
    }
    fCondensed.UGlobal(u1, elsol);
    count = 0;
    for (ic=0; ic<nc0 ; ic++) {
        TPZConnect &c = Connect(ic);
        int seqnum = c.SequenceNumber();
        int blsize = bl.Size(seqnum);
        for (int ibl=0; ibl<blsize; ibl++) {
            bl(seqnum,ibl,0,0) = elsol(count++,0);
        }
    }
}
