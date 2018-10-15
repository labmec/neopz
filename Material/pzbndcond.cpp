/**
 * \file
 * @brief Contains implementations of the TPZBndCond methods.
 */

#include "pzbndcond.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.bndcond"));
#endif

int TPZBndCond::TPZ_BCDefine::ClassId() const {
    return Hash("TPZBndCond::TPZ_BCDefine");
}

void TPZBndCond::TPZ_BCDefine::Read(TPZStream& buf, void* context) {
    fBCVal2.Read(buf, context);
    fForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fForcingFunctionExact = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTimeDependentForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTimedependentFunctionExact = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fBCForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTimedependentBCForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
}

void TPZBndCond::TPZ_BCDefine::Write(TPZStream& buf, int withclassid) const {
    fBCVal2.Write(buf, withclassid);
    TPZPersistenceManager::WritePointer(fForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fForcingFunctionExact.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fTimeDependentForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fTimedependentFunctionExact.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fBCForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fTimedependentBCForcingFunction.operator ->(), &buf);
}

void TPZBndCond::Clone(std::map<int, TPZMaterial * > &matvec) {
	int matid = Id();
	
	TPZMaterial * refmaterial = Material();
	TPZMaterial * newrefmaterial = NULL;
	int refmatid = 0;
	if(refmaterial) {
		refmaterial->Clone(matvec);
		refmatid = refmaterial->Id();
		newrefmaterial = matvec[refmatid];
	}
	std::map<int, TPZMaterial * >::iterator matit;
	matit = matvec.find(matid);
	if(matit == matvec.end())
	{
		TPZMaterial * newmat = (new TPZBndCond(*this, newrefmaterial));
		matvec[matid] = newmat;
	}
}

void TPZBndCond::InterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZSolVec &rightu,TPZSolVec &jump){
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial);
	
	if(!mat) return;
	if(fForcingFunction) {
		TPZManVector<STATE> result(fBCVal2.Rows());
		fForcingFunction->Execute(x,result);
		int i;
		for(i=0; i<fBCVal2.Rows(); i++) {
			fBCVal2(i,0) = result[i];
		}
	}
	
	if(leftu.NElements() == 0) {
		mat->BCInterfaceJump(x, rightu, *this, jump);
		return;
	}
	
	if(rightu.NElements() == 0) {
		mat->BCInterfaceJump(x, leftu, *this, jump);
		return;
	}
	
	PZError << __PRETTY_FUNCTION__ << " - Huge problem. Both leftu and rightu contain elements. Wich one is the actual element neighbour to the Dirichlet boundary ?" << std::endl;
	DebugStop();
	
}//InterfaceJump

int TPZBndCond::ClassId() const{
    return Hash("TPZBndCond") ^ TPZDiscontinuousGalerkin::ClassId() << 1;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZBndCond>;
#endif

void TPZBndCond::Write(TPZStream &buf, int withclassid) const {
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(fBCs);
    buf.Write(&fType);
    fBCVal1.Write(buf, withclassid);
    fBCVal2.Write(buf, withclassid);
    TPZPersistenceManager::WritePointer(fMaterial, &buf);
}

void TPZBndCond::Read(TPZStream &buf, void *context){
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(fBCs);
    buf.Read(&fType);
    fBCVal1.Read(buf, context);
    fBCVal2.Read(buf, context);
    fMaterial = dynamic_cast<TPZMaterial *>(TPZPersistenceManager::GetInstance(&buf));
    if (!fMaterial) {
        std::cout << " reading a boundary condition without material object!!\n";
#ifdef LOG4CXX
        LOGPZ_FATAL(logger, "reading a boundary condition without material object!!");
#endif
    }
}

void TPZBndCond::ContributeInterfaceErrors( TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
										   REAL weight,
										   TPZVec<STATE> &nkL,
										   TPZVec<STATE> &nkR,
										   int &errorid){
	
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial);
	if(!mat) return;
	
    UpdateBCValues(data.x);
	
	if(dataleft.sol.NElements() < dataright.sol.NElements()){
		//		data.InvertLeftRightData();
        for(int i=0; i<3; i++) data.normal[i] *= -1.;
		mat->ContributeInterfaceBCErrors(data,dataright,weight,nkR,*this, errorid);
	}
	else {
		mat->ContributeInterfaceBCErrors(data,dataleft,weight,nkL,*this, errorid);
	}
	
}

void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    this->fMaterial->ContributeBC(data, weight, ek, ef, *this);
}

//----
void TPZBndCond::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    UpdateBCValues(datavec[0].x);
	this->fMaterial->ContributeBC(datavec,weight,ek,ef,*this);
}
//----

void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
	UpdateBCValues(data.x);
	this->fMaterial->ContributeBC(data,weight,ef,*this);
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
	if(!mat) DebugStop();
    UpdateBCValues(data.x);

	
	if(dataleft.phi.Rows() == 0){//it meanst right data has been filled
		//left data should be filled instead of right data
        // shouldn't we invert the normal vector?
        for (int i=0; i<3; i++) data.normal[i] *= -1.;
        mat->ContributeBCInterface(data,dataright,weight,ek,ef,*this);
	}
    else
    {
        mat->ContributeBCInterface(data,dataleft,weight,ek,ef,*this);
    }
}//void

void TPZBndCond::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
    if(!mat) DebugStop();
	UpdateBCValues(data.x);
	
    int nel = dataleft.size();
    if(!nel) DebugStop();
    
	if(dataleft[0].phi.Rows() == 0){//it meanst right data has been filled
		//left data should be filled instead of right data
        // shouldn't we invert the normal vector?
        for (int i=0; i<3; i++) data.normal[i] *= -1.;
        
        mat->ContributeBCInterface(data,dataright,weight,ek,ef,*this);
	}
    else
    {
        mat->ContributeBCInterface(data,dataleft,weight,ek,ef,*this);
    }

    
}

void TPZBndCond::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef){
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
    UpdateBCValues(data.x);
	if(dataleft.phi.Rows() == 0){//it meanst right data has been filled
		//left data should be filled instead of right data
        for(int i=0; i<3; i++) data.normal[i] *= -1.;
        mat->ContributeBCInterface(data,dataright,weight,ef,*this);
		//		data.InvertLeftRightData();
	} else
    {
        mat->ContributeBCInterface(data,dataleft,weight,ef,*this);
    }
}

void TPZBndCond::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	DebugStop();//nothing to be done here
}


void TPZBndCond::UpdateBCValues(TPZManVector<REAL,3> & x){
	if(fForcingFunction){
		TPZManVector<STATE,3> result(fBCVal2.Rows(),0.);
        TPZFNMatrix<9,STATE> gradu(Dimension(),fBCVal2.Rows());
        // use gradu to update the Neumann boundary condition
		fForcingFunction->Execute(x,result,gradu);
        
#ifdef PZDEBUG
        if (fBCVal2.Rows() != result.size()) {
            DebugStop();
        }
#endif
		int i;
		for(i=0; i<fBCVal2.Rows(); i++) {
			fBCVal2(i,0) = result[i];
		}
	}
	
	if( this->fValFunction ){
		TPZManVector<STATE> result(this->fBCVal2.Rows(),0.);
		this->fValFunction(x, this->fBCVal1, result, this->fType );
		int i;
		for(i = 0; i < this->fBCVal2.Rows(); i++) {
			this->fBCVal2(i,0) = result[i];
		}
	}//if
    
    int nbc = fBCs.size();
    for (int ibc=0; ibc<nbc; ibc++) {
        if (fBCs[ibc].fForcingFunction) {
            TPZManVector<STATE> result(fBCVal2.Rows(),0.);
            fBCs[ibc].fForcingFunction->Execute(x,result);
            int i;
            for(i=0; i<fBCVal2.Rows(); i++) {
                fBCs[ibc].fBCVal2(i,0) = result[i];
            }

        }
    }
}


void TPZBndCond::UpdateBCValues(TPZVec<TPZMaterialData> &datavec){
	
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}


void TPZBndCond::FillDataRequirements(TPZMaterialData &data){
	if(!fMaterial)
	{
		PZError << "\nUnable to call TPZBndCond::fMaterial::FillDataRequirements - fMaterial pointer is null!\n";
		return;
	}
	fMaterial->FillBoundaryConditionDataRequirement(fType,data);
	if(fLinearContext == false || fType == 50){
		data.fNeedsSol = true;
	}
}

void TPZBndCond::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
	if(!fMaterial)
	{
		PZError << "\nUnable to call TPZBndCond::fMaterial::FillDataRequirements - fMaterial pointer is null!\n";
		return;
	}
	fMaterial->FillBoundaryConditionDataRequirement(fType,datavec);
    int nref = datavec.size();
    
	if(fLinearContext == false){
        for(int iref=0; iref<nref; iref++){
            datavec[iref].fNeedsSol = true;
        }
	}
}
