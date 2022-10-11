#include "TPZCompElUnitaryLagrange.h"
#include "pzcompel.h"
#include "pzcmesh.h"
#include "TPZElementMatrixT.h"

#include "pzmultiphysicselement.h"

TPZCompElUnitaryLagrange::TPZCompElUnitaryLagrange(TPZCompMesh &mesh, TPZGeoEl *reference, TPZCompElSide &wrapSide, TPZCompElSide &lagrangeSide) :
    TPZCompElDisc(mesh, reference) {
    this->fConnectIndexes.resize(2);

    //Initializes the connect indexes structure - set the first connect as the same as Lagrange element and the second as the wrap element
    fConnectIndexes[0] = lagrangeSide.Element()-> ConnectIndex(0);
    fConnectIndexes[1] = wrapSide.Element()-> ConnectIndex(0);
    //Resize the dependency matrix and set as identity. Its size is always 1x1 since the semi-hybridization corresponds to the constant flux,
    //i.e., a connect with only one function.
    fDepMatrix.Resize(1,1);
    fDepMatrix(0,0) = 1.;

    //Gets the side orient from the wrap's volumetric element neighbour
    TPZGeoElSide volside = wrapSide.Reference().operator--();
    TPZCompElSide cvolside = volside.Reference();
    TPZCompEl *celvol = cvolside.Element();

    //Sets the side orient
    TPZMultiphysicsElement *celmulti = dynamic_cast<TPZMultiphysicsElement *> (celvol); 
    TPZInterpolationSpace *celv = dynamic_cast<TPZInterpolationSpace *> (celmulti->Element(0)); 
    fSideOrient = celv->GetSideOrient(cvolside.Side());

    //Gets the correct dependency matrix if the wrap element has dependency
    TPZConnect &c = wrapSide.Element()->Connect(0);
    if (c.HasDependency()){
        auto *ptr = c.FirstDepend();
        int count = 0;
        while(ptr){
            // std::cout << "Dependency matrix = " << ptr->fDepMatrix << std::endl;
            if (count == 1) fDepMatrix = ptr->fDepMatrix;
            count++;
            ptr = ptr->fNext;
        }
    }
}

int TPZCompElUnitaryLagrange::NConnects() const {
	return 2;
}

int TPZCompElUnitaryLagrange::NSideConnects(int side) const{
//     // o erro provavelmente esta aqui
// 	if(TSHAPE::SideDimension(side)<= this->Dimension()-2) return 0;
// 	if(TSHAPE::SideDimension(side)== this->Dimension()-1) return 2;
// 	if(TSHAPE::SideDimension(side)== this->Dimension()) {
//         int ncon = 1;
//         return ncon;
//     }
// #ifdef PZ_LOG
// 	{
// 		std::stringstream sout;
// 		sout << __PRETTY_FUNCTION__ << "Side: " << side <<"unhandled case ";
// 		LOGPZ_ERROR(logger,sout.str())
// 	}
// #endif
	return -1;

}

int TPZCompElUnitaryLagrange::NConnectShapeF(int connect, int order)const
{
    return 1;
}


int64_t TPZCompElUnitaryLagrange::ConnectIndex(int con) const{
    if (con < 0 || con >=2) DebugStop();
	return this->fConnectIndexes[con];
}

int TPZCompElUnitaryLagrange::SideConnectLocId(int node,int side) const {
    // if (TSHAPE::Dimension == 2){
    //     return 2*(side-TSHAPE::NCornerNodes);
    // } else if (TSHAPE::Dimension == 3){
    //     return 2*(side-(TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1));
    // } else {
    //     DebugStop();
    // }
    
    return -1;
}


void TPZCompElUnitaryLagrange::SetConnectIndex(int i, int64_t connectindex)
{
	this->fConnectIndexes[i] = connectindex;
}

template<class TVar>
void TPZCompElUnitaryLagrange::CalcStiffInternal(TPZElementMatrixT<TVar> &ek,TPZElementMatrixT<TVar> &ef){
    InitializeElementMatrix(ek, ef);

    //Compute the stiffness matrix (constant Lagrange multiplier), scaled by the dependency matrix.
    //In this case it corresponds to a scalar relating the areas of the wrap and large elements
    int fSize = ek.fMat.Rows();
    int nvar = fSize/2;
    for (int i = 0; i < nvar; i++){
        ek.fMat(i,nvar + i) = TVar(fSideOrient) / fDepMatrix(0,0);
        ek.fMat(nvar + i,i) = TVar(fSideOrient) / fDepMatrix(0,0);
    }

    ef.fMat.Zero();
}


void TPZCompElUnitaryLagrange::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef){
	int numdof = 1;
    TPZMaterial *mat = this->Material();
    const int numloadcases = [mat](){
        if (auto *tmp = dynamic_cast<TPZMatLoadCasesBase*>(mat); tmp){
            return tmp->NumLoadCases();
        }else{
            return 1;
        }
    }();
	const int ncon = this->NConnects();
    
    ek.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
    
	ek.Block().SetNBlocks(ncon);
	ef.Block().SetNBlocks(ncon);
	int i;
    int numeq=0;
	for(i=0; i<ncon; i++){
        TPZConnect &c = Connect(i);
        int nshape = c.NShape();
        int nstate = c.NState();
        
		ek.Block().Set(i,nshape*nstate);
		ef.Block().Set(i,nshape*nstate);
        numeq += nshape*nstate;
	}
	ek.Matrix().Redim(numeq,numeq);
	ef.Matrix().Redim(numeq,numloadcases);
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
}//void