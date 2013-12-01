//
//  pzreducedspace.cpp
//  PZ
//
//  Created by Philippe Devloo on 7/30/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "pzreducedspace.h"
#include "tpzcompmeshreferred.h"
#include "pzmaterial.h"
#include "pzelmat.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZInterpolationSpace"));
#endif

/** @brief Default constructor */
TPZReducedSpace::TPZReducedSpace() : TPZInterpolationSpace()
{
    
}

/** @brief Default destructor */
TPZReducedSpace::~TPZReducedSpace()
{
    
}

/** @brief Puts a copy of the element in the referred mesh */
TPZReducedSpace::TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy) : TPZInterpolationSpace(mesh,copy)
{
    
}

/** @brief Puts a copy of the element in the patch mesh */
TPZReducedSpace::TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy, std::map<long,long> &gl2lcElMap) : TPZInterpolationSpace(mesh,copy,gl2lcElMap)
{
    
}

/** @brief Copy of the element in the new mesh whit alocated index */
TPZReducedSpace::TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy, long &index) : TPZInterpolationSpace(mesh,copy,index)
{
    
}

/**
 * @brief Create a computational element within mesh
 * @param mesh mesh wher will be created the element
 * @param gel geometrical element to insert
 * @param index new elemen index
 */
/** Inserts the element within the data structure of the mesh */
TPZReducedSpace::TPZReducedSpace(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) : TPZInterpolationSpace(mesh,gel,index)
{
    
}

/** @brief It returns the shapes number of the element */
int TPZReducedSpace::NShapeF() const
{
    TPZConnect &c = Connect(0);
    return c.NShape();
}

/** @brief Returns the number of shapefunctions associated with a connect*/
int TPZReducedSpace::NConnectShapeF(int inod) const
{
#ifdef DEBUG
    if (inod != 0) {
        DebugStop();
    }
#endif
    TPZConnect &c = Connect(0);
    return c.NShape();
}

/** @brief Returns the max order of interpolation. */
int TPZReducedSpace::MaxOrder()
{
    TPZInterpolationSpace *intel = ReferredIntel();
    return intel->MaxOrder();
}

/** 
 * @brief Computes the shape function set at the point x. 
 * @param qsi point in master element coordinates
 * @param phi vector of values of shapefunctions, dimension (numshape,1)
 * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
 */
/**
 * This method uses the order of interpolation
 * of the element along the sides to compute the number of shapefunctions
 */
void TPZReducedSpace::Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi)
{
    DebugStop();
}

/** 
 * @brief Computes the shape function set at the point x. 
 * @param qsi point in master element coordinates
 * @param phi vector of values of shapefunctions, dimension (numshape,1)
 * @param dphix matrix of derivatives of shapefunctions, dimension (dim,numshape)
 */
/**
 * This method uses the order of interpolation
 * of the element along the sides to compute the number of shapefunctions
 */
void TPZReducedSpace::ShapeX(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphix, TPZFMatrix<REAL> &axes)
{
#ifdef STATE_COMPLEX
    DebugStop();
#else
    TPZInterpolationSpace *intel = ReferredIntel();
    TPZSolVec sol;
    TPZGradSolVec dsol;
    intel->ComputeSolution(qsi, sol, dsol, axes);
    int nsol = sol.size();
    int nstate = sol[0].size();
    int dim = axes.Rows();
    phi.Resize(nstate, nsol);
    dphix.Resize(nstate*dim, nsol);
    for (int isol =0; isol<nsol; isol++) {
        for (int istate=0; istate<nstate; istate++) {
            phi(istate,isol) = sol[isol][istate];
            for (int id=0; id<dim; id++) {
                dphix(id+istate*dim,isol) = dsol[isol](id,istate);
            }
        }
    }
#endif
}


void TPZReducedSpace::ShapeX(TPZVec<REAL> &qsi,TPZMaterialData &data)
{
#ifdef STATE_COMPLEX
    DebugStop();
#else
    TPZInterpolationSpace *intel = ReferredIntel();
    
    intel->ComputeSolution(qsi, data.sol, data.dsol, data.axes);
    long nsol = data.sol.size();
    int nstate = data.sol[0].size();
    int dim = data.axes.Rows();
    data.phi.Resize(nstate, nsol);
    data.dphix.Resize(nstate*dim, nsol);
    for (long isol =0; isol<nsol; isol++) {
        for (int istate=0; istate<nstate; istate++) {
            data.phi(istate,isol) = data.sol[isol][istate];
            for (int id=0; id<dim; id++) {
                data.dphix(id+istate*dim,isol) = data.dsol[isol](id,istate);
            }
        }
    }
#endif
}


void TPZReducedSpace::ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data){
    ShapeX(qsi,data);
}

void TPZReducedSpace::ComputeSolution(TPZVec<REAL> &qsi,TPZMaterialData &data){
    
    ComputeSolution(qsi, data.phi, data.dphix, data.axes,data.sol,data.dsol);
}

/**
 * @brief Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */
void TPZReducedSpace::InitMaterialData(TPZMaterialData &data)
{
    data.fShapeType = TPZMaterialData::EVecShape;
    TPZMaterial *mat = Material();
    mat->FillDataRequirements(data);
}

/** @brief Compute and fill data with requested attributes */
void TPZReducedSpace::ComputeRequiredData(TPZMaterialData &data,
                                 TPZVec<REAL> &qsi)
{
    data.intGlobPtIndex = -1;
    ShapeX(qsi, data.phi, data.dphix, data.axes);
    
    if (data.fNeedsSol) {
        ComputeSolution(qsi, data.phi, data.dphix, data.axes, data.sol, data.dsol);
    }
    if (data.fNeedsHSize){
		data.HSize = 2.*this->InnerRadius();
	}//fNeedHSize
    
    if(data.fNeedsNormal){
        this->ComputeNormal(data);
    }
    data.x.Resize(3., 0.);
    Reference()->X(qsi, data.x);
    
    int dim = Reference()->Dimension();
    data.jacobian.Resize(dim,dim);
    data.jacinv.Resize(dim,dim);
    Reference()->Jacobian(qsi, data.jacobian, data.axes, data.detjac, data.jacinv);
}


/** @brief Initialize element matrix in which is computed CalcStiff */
void TPZReducedSpace::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    TPZMaterial *mat = this->Material();
	const int numdof = 1;
	const int ncon = this->NConnects();
#ifdef DEBUG
    if (ncon != 1) {
        DebugStop();
    }
#endif
	const int nshape = this->NShapeF();
	const int numeq = nshape*numdof;
    const int numloadcases = mat->NumLoadCases();
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = numdof;
	ef.fNumStateVars = numdof;
	int i;
	for(i=0; i<ncon; i++){
        unsigned int nshape = NConnectShapeF(i);
#ifdef DEBUG
        TPZConnect &c = Connect(i);
        if(c.NShape() != nshape || c.NState() != numdof)
        {
            DebugStop();
        }
#endif
		ek.fBlock.Set(i,nshape*numdof);
		ef.fBlock.Set(i,nshape*numdof);
	}
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
    
}

/** @brief Initialize element matrix in which is computed in CalcResidual */
void TPZReducedSpace::InitializeElementMatrix(TPZElementMatrix &ef)
{
    TPZMaterial *mat = this->Material();
	const int numdof = 1;
	const int ncon = this->NConnects();
#ifdef DEBUG
    if (ncon != 1) {
        DebugStop();
    }
#endif
	const int nshape = this->NShapeF();
	const int numeq = nshape*numdof;
    const int numloadcases = mat->NumLoadCases();
	ef.fMat.Redim(numeq,numloadcases);
	ef.fBlock.SetNBlocks(ncon);
	ef.fNumStateVars = numdof;
	int i;
	for(i=0; i<ncon; i++){
        unsigned int nshape = NConnectShapeF(i);
#ifdef DEBUG
        TPZConnect &c = Connect(i);
        if(c.NShape() != nshape || c.NState() != numdof)
        {
            DebugStop();
        }
#endif
		ef.fBlock.Set(i,nshape*numdof);
	}
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
	}    
}

/** @brief Save the element data to a stream */
void TPZReducedSpace::Write(TPZStream &buf, int withclassid)
{
    TPZInterpolationSpace::Write(buf, withclassid);
}

/** @brief Read the element data from a stream */
void TPZReducedSpace::Read(TPZStream &buf, void *context)
{
    TPZInterpolationSpace::Read(buf, context);
}

TPZInterpolationSpace *TPZReducedSpace::ReferredIntel() const
{
    TPZCompMesh *cmesh = Mesh();
    TPZCompMeshReferred *cmeshref = dynamic_cast<TPZCompMeshReferred *>(cmesh);
    
#ifdef DEBUG
    if (!cmeshref) {
        DebugStop();
    }
#endif
    
    TPZCompEl *cel = cmeshref->ReferredEl(Index());
    
#ifdef DEBUG
    if (!cel) {
        DebugStop();
    }
#endif
    
    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);

#ifdef DEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    return intel;
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
void TPZReducedSpace::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                                const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol)
{
    const int dim = this->Reference()->Dimension();
    const int numdof = this->Material()->NStateVariables();
    
#ifdef DEBUG
    const int ncon = this->NConnects();
    if (ncon != 1) {
        DebugStop();
    }
#endif
    
    TPZFMatrix<STATE> &MeshSol = Mesh()->Solution();
    long numbersol = MeshSol.Cols();
    sol.Resize(numbersol);
    dsol.Resize(numbersol);
	
    for (long is=0 ; is<numbersol; is++) {
        sol[is].Resize(numdof);
        sol[is].Fill(0.);
        dsol[is].Redim(dim, numdof);
        dsol[is].Zero();
    }
	
    TPZBlock<STATE> &block = Mesh()->Block();
    long d;
    TPZConnect *df = &this->Connect(0);
    long dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    long pos = block.Position(dfseq);
    for(int jn=0; jn<dfvar; jn++) {
#ifdef DEBUG
        {
            if(phi.Rows() != numdof)
            {
                DebugStop();
            }
            if (phi.Cols() != dfvar) {
                DebugStop();
            }
            if (dphix.Rows() != dim*numdof) {
                DebugStop();
            }
            if (dphix.Cols() != dfvar) {
                DebugStop();
            }
            
        }
#endif
        for (long is=0; is<numbersol; is++) {
            for(d=0; d<numdof; d++){
                sol[is][d%numdof] += (STATE)phi(d,jn)*MeshSol(pos+jn,is);
            }
            for(d=0; d<dim*numdof; d++){
                dsol[is](d%dim,d/dim) += (STATE)dphix(d,jn)*MeshSol(pos+jn,is);
            }
        }
    }
}

static TPZCompEl * CreateReducedElement(TPZGeoEl *gel,TPZCompMesh &mesh,long &index)
{
    return new TPZReducedSpace(mesh,gel,index);
}

void TPZReducedSpace::SetAllCreateFunctionsReducedSpace(TPZCompMesh *cmesh)
{
    TPZManVector<TCreateFunction,10> functions(8);
    functions[EPoint] = CreateReducedElement;
	functions[EOned] = CreateReducedElement;
	functions[EQuadrilateral] = CreateReducedElement;
	functions[ETriangle] = CreateReducedElement;
	functions[EPrisma] = CreateReducedElement;
	functions[ETetraedro] = CreateReducedElement;
	functions[EPiramide] = CreateReducedElement;
	functions[ECube] = CreateReducedElement;
    cmesh->ApproxSpace().SetCreateFunctions(functions);
}

TPZCompEl* TPZReducedSpace::Clone(TPZCompMesh &mesh) const{
    return new TPZReducedSpace(mesh, *this);
}

TPZCompEl * TPZReducedSpace::ClonePatchEl (TPZCompMesh &mesh, std::map< long, long > &gl2lcConMap, std::map< long, long > &gl2lcElMap) const {
    return new TPZReducedSpace(mesh,*this,gl2lcElMap);
}

#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "tpzgraphelt2dmapped.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel.h"


void TPZReducedSpace::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
	TPZGeoEl *ref = Reference();
	int mat = Material()->Id();
	int nsides = ref->NSides();
	
	if(dimension == 2 && mat > 0){
		if(nsides == 9){
			new TPZGraphElQ2dd(this,&grmesh);
			return;
		}
		if(nsides == 7){
			new TPZGraphElT2dMapped(this,&grmesh);
			return;
		}
	}//2d
	
	if(dimension == 3 && mat > 0){
		if(nsides == 27){
			new TPZGraphElQ3dd(this,&grmesh);
			return;
		}//cube
		if(nsides == 21){
			new TPZGraphElPrismMapped(this,&grmesh);
			return;
		}//prism
		if(nsides == 15){
			new TPZGraphElT3d(this,&grmesh);
			return;
		}//tetra
		if(nsides == 19){
			new TPZGraphElPyramidMapped(this,&grmesh);
			return;
		}//pyram
	}//3d
	
	if(dimension == 1 && mat > 0){
		new TPZGraphEl1dd(this,&grmesh);
	}//1d
}
