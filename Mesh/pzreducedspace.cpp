//
//  pzreducedspace.cpp
//  PZ
//
//  Created by Philippe Devloo on 7/30/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "pzcmesh.h"
#include "pzreducedspace.h"
#include "pzmultiphysicselement.h"
#include "TPZMaterial.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatLoadCases.h"
#include "pzelmat.h"
#include "pzlog.h"
#include "pzerror.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZInterpolationSpace");
#endif

/** @brief Default constructor */
TPZReducedSpace::TPZReducedSpace() : TPZRegisterClassId(&TPZReducedSpace::ClassId),
TPZInterpolationSpace()
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" should be reimplemented without TPZCompMeshReferred\n";
    PZError<<"Aborting...";
}

/** @brief Default destructor */
TPZReducedSpace::~TPZReducedSpace()
{
    
}


/** @brief Puts a copy of the element in the patch mesh */
TPZReducedSpace::TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy, std::map<int64_t,int64_t> &gl2lcElMap) : TPZRegisterClassId(&TPZReducedSpace::ClassId),
TPZInterpolationSpace(mesh,copy,gl2lcElMap)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" should be reimplemented without TPZCompMeshReferred\n";
    PZError<<"Aborting...";
    DebugStop();
}

/** @brief Copy of the element in the new mesh whit alocated index */
TPZReducedSpace::TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy) : TPZRegisterClassId(&TPZReducedSpace::ClassId),
TPZInterpolationSpace(mesh,copy)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" should be reimplemented without TPZCompMeshReferred\n";
    PZError<<"Aborting...";
    DebugStop();
}

/**
 * @brief Create a computational element within mesh
 * @param mesh mesh wher will be created the element
 * @param gel geometrical element to insert
 * @param index new elemen index
 */
/** Inserts the element within the data structure of the mesh */
TPZReducedSpace::TPZReducedSpace(TPZCompMesh &mesh, TPZGeoEl *gel) : TPZRegisterClassId(&TPZReducedSpace::ClassId),
TPZInterpolationSpace(mesh,gel)
{
    std::cout << "Creating reduced with dim " << gel->Dimension() << '\n';
}

/** @brief It returns the shapes number of the element */
int TPZReducedSpace::NShapeF() const
{
    TPZConnect &c = Connect(0);
    return c.NShape();
}

/** @brief Returns the number of shapefunctions associated with a connect*/
int TPZReducedSpace::NConnectShapeF(int inod, int order) const
{
#ifdef PZDEBUG
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
    TPZInterpolationSpace *intel = ReferredIntel();
    TPZMaterialDataT<STATE> inteldata;
    inteldata.fNeedsSol=true;
    constexpr bool hasPhi{false};
    intel->ComputeSolution(qsi,inteldata,hasPhi);
    TPZSolVec<STATE> &sol = inteldata.sol;
    TPZGradSolVec<STATE> &dsol = inteldata.dsol;
    int nsol = sol.size();
    int nstate = sol[0].size();
    int dim = axes.Rows();
    phi.Resize(nsol,nstate);
    dphix.Resize(nstate*dim,nsol);
    for (int isol =0; isol<nsol; isol++) {
        for (int istate=0; istate<nstate; istate++) {
            phi(isol,istate) = sol[isol][istate];
            for (int id=0; id<dim; id++) {
                dphix(id+istate*dim,isol) = dsol[isol](id,istate);
            }
        }
    }
}


void TPZReducedSpace::ShapeX(TPZVec<REAL> &qsi,TPZMaterialDataT<STATE> &data)
{
    TPZInterpolationSpace *intel = ReferredIntel();

    {
        TPZMaterialDataT<STATE> inteldata;
        inteldata.fNeedsSol=true;
        constexpr bool hasPhi{false};
        intel->ComputeSolution(qsi,inteldata,hasPhi);
        data.sol = std::move(inteldata.sol);
        data.dsol = std::move(inteldata.dsol);
    }
    int64_t nsol = data.sol.size();
    int nstate = data.sol[0].size();
    int dim = data.axes.Rows();
    data.phi.Resize(nsol,nstate);
    data.dphix.Resize(nstate*dim,nsol);
    for (int64_t isol =0; isol<nsol; isol++) {
        for (int istate=0; istate<nstate; istate++) {
            data.phi(isol,istate) = data.sol[isol][istate];
            for (int id=0; id<dim; id++) {
                data.dphix(id+istate*dim,isol) = data.dsol[isol](id,istate);
            }
        }
    }
}


void TPZReducedSpace::ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data){
    auto *tmp = dynamic_cast<TPZMaterialDataT<STATE>*>(&data);
    if(!tmp){
        PZError<<__PRETTY_FUNCTION__;
        PZError<< "is not available to complex types yet.\nAborting...\n";
        DebugStop();
    }
    ShapeX(qsi,*tmp);
}

/**
 * @brief Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */
void TPZReducedSpace::InitMaterialData(TPZMaterialData &data)
{
    data.fShapeType = TPZMaterialData::EVecShape;
    auto *mat =
        dynamic_cast<TPZMatSingleSpace*>(this->Material());
    mat->FillDataRequirements(data);
    TPZConnect &c = Connect(0);
    int nshape = c.NShape();
    int nstate = c.NState();
    int dim = Reference()->Dimension();
    data.phi.Resize(nshape, nstate);
    data.dphix.Resize(dim*nstate,nshape);
    data.jacobian.Resize(dim,dim);
    data.jacinv.Resize(dim,dim);
    data.axes.Resize(dim, 3);
}

/** @brief Compute and fill data with requested attributes */
void TPZReducedSpace::ComputeRequiredData(TPZMaterialDataT<STATE> &data,
                                 TPZVec<REAL> &qsi)
{
    data.intGlobPtIndex = -1;
    int dim = Reference()->Dimension();
    Reference()->Jacobian(qsi, data.jacobian, data.axes, data.detjac, data.jacinv);
    ShapeX(qsi, data.phi, data.dphix, data.axes);

    if (data.fNeedsSol) {
        ReallyComputeSolution(data);
    }
    if (data.fNeedsHSize){
		data.HSize = 2.*this->InnerRadius();
	}//fNeedHSize
    
    if(data.fNeedsNormal){
        this->ComputeNormal(data);
    }
    data.x.Resize(3., 0.);
    Reference()->X(qsi, data.x);
    
}


/** @brief Initialize element matrix in which is computed CalcStiff */
void TPZReducedSpace::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    TPZMaterial *mat = this->Material();
	const int ncon = this->NConnects();
#ifdef PZDEBUG
    if (ncon != 1) {
        DebugStop();
    }
#endif
	const int nshape = this->NShapeF();
	const int numeq = nshape;
    const int numloadcases = [mat](){
        if (auto *tmp = dynamic_cast<TPZMatLoadCasesBase*>(mat); tmp){
            return tmp->NumLoadCases();
        }else{
            return 1;
        }
    }();
	ek.Matrix().Redim(numeq,numeq);
	ef.Matrix().Redim(numeq,numloadcases);
    auto &ekBlock = ek.Block();
    auto &efBlock = ef.Block();
	ekBlock.SetNBlocks(ncon);
	efBlock.SetNBlocks(ncon);
    ek.fMesh = Mesh();
    ef.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
    ef.fType = TPZElementMatrix::EF;
	int i;
	for(i=0; i<ncon; i++){
        unsigned int nshape = Connect(i).NShape();
		ekBlock.Set(i,nshape);
		efBlock.Set(i,nshape);
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
#ifdef PZDEBUG
    if (ncon != 1) {
        DebugStop();
    }
#endif
	const int nshape = this->NShapeF();
	const int numeq = nshape*numdof;
    const int numloadcases = [mat](){
        if (auto *tmp = dynamic_cast<TPZMatLoadCasesBase*>(mat); tmp){
            return tmp->NumLoadCases();
        }else{
            return 1;
        }
    }();
	ef.Matrix().Redim(numeq,numloadcases);
	ef.Block().SetNBlocks(ncon);
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
	int i;
	for(i=0; i<ncon; i++){
        unsigned int nshape = Connect(i).NShape();
		ef.Block().Set(i,nshape*numdof);
	}
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
	}    
}

/** @brief Save the element data to a stream */
void TPZReducedSpace::Write(TPZStream &buf, int withclassid) const
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
    return fReferred;
//     TPZCompMeshReferred *cmeshref = dynamic_cast<TPZCompMeshReferred *>(cmesh);
    
// #ifdef PZDEBUG
//     if (!cmeshref) {
//         DebugStop();
//     }
// #endif
    
//     TPZCompEl *cel = cmeshref->ReferredEl(Index());
    
// #ifdef PZDEBUG
//     if (!cel) {
//         DebugStop();
//     }
// #endif
    
//     TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
    
    
//     TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement *>(cel);

// #ifdef PZDEBUG
//     if (!intel && !mf_cel) {
//         DebugStop();
//     }
// #endif
    
//     if (intel) {
//         return intel;
//     }
    
//     TPZInterpolationSpace * intel_mf = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(0)); //@omar:: garbage solution
//     if(intel_mf){
//         return intel_mf;
//     }
    
//     return intel;
}

void TPZReducedSpace::ReallyComputeSolution(TPZMaterialDataT<STATE>& data)
{
    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &dphix = data.dphix;
    const TPZFMatrix<REAL> &axes = data.axes;
    TPZSolVec<STATE> &sol = data.sol;
    TPZGradSolVec<STATE> &dsol = data.dsol;
    const int dim = axes.Rows();//this->Reference()->Dimension();
    const int nstate = this->Material()->NStateVariables();
    
#ifdef PZDEBUG
    const int ncon = this->NConnects();
    if (ncon != 1) {
        DebugStop();
    }
#endif
    
//    TPZInterpolationSpace *intel = ReferredIntel();
//    TPZSolVec sol_t;
//    TPZGradSolVec dsol_t;
//    TPZFMatrix<REAL> axes_t = axes;
//    intel->ComputeSolution(qsi, sol_t, dsol_t, axes_t);
    
    int nsol = phi.Rows();//sol_t.size();
//    int numdof = sol_t[0].size();
//    int dim = axes_t.Rows();
    
    TPZFMatrix<STATE> &MeshSol = Mesh()->Solution();
    int64_t numbersol = MeshSol.Cols();
    int64_t numberdof = MeshSol.Rows();
    sol.Resize(numbersol);
    dsol.Resize(numbersol);
	
    for (int64_t is=0 ; is<numbersol; is++) {
        sol[is].Resize(nstate);
        sol[is].Fill(0.);
        dsol[is].Redim(dim, nstate);
    }
    
    TPZBlock &block = Mesh()->Block();
    TPZConnect *df = &this->Connect(0);
    int64_t dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int64_t pos = block.Position(dfseq);
    
#ifdef PZDEBUG
    {
        if(nsol * nstate != dfvar)
        {
            DebugStop();
        }
        
    }
    
    
#endif
    
    for(int ib=0; ib < nsol; ib++) {

        for (int64_t is=0; is<numbersol; is++) {
            
            for(int64_t iv = 0; iv < nstate; iv++){
                sol[is][iv%nstate] += (STATE)phi.GetVal(ib,iv)*MeshSol(pos+ib*nstate+iv,is);
                
                for(int64_t id = 0; id < dim; id++){
                    dsol[is](id,iv%nstate) += (STATE)dphix.GetVal(id+iv*dim,ib)*MeshSol(pos+ib*nstate+iv,is);
                }
            }
        }
    }
    
}

static TPZCompEl * CreateReducedElement(TPZGeoEl *gel,TPZCompMesh &mesh)
{
    return new TPZReducedSpace(mesh,gel);
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

TPZCompEl * TPZReducedSpace::ClonePatchEl (TPZCompMesh &mesh, std::map< int64_t, int64_t > &gl2lcConMap, std::map< int64_t, int64_t > &gl2lcElMap) const {
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
#include "pzgraphmesh.h"


void TPZReducedSpace::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
	TPZGeoEl *ref = Reference();
	int matid = Material()->Id();
	int nsides = ref->NSides();
    bool to_postpro = grmesh.Material_Is_PostProcessed(matid);
	if(dimension == 2 && to_postpro){
		if(nsides == 9){
			new TPZGraphElQ2dd(this,&grmesh);
			return;
		}
		if(nsides == 7){
			new TPZGraphElT2dMapped(this,&grmesh);
			return;
		}
	}//2d
	
	if(dimension == 3 && to_postpro){
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
	
	if(dimension == 1 && to_postpro){
		new TPZGraphEl1dd(this,&grmesh);
	}//1d
}

int TPZReducedSpace::ClassId() const{
    return Hash("TPZReducedSpace") ^ TPZInterpolationSpace::ClassId() << 1;
}
