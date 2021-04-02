/**
 * \file
 * @brief Contains implementations of the TPZBndCond methods.
 */

#include "pzbndcond.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.bndcond");
#endif

int TPZBndCond::TPZ_BCDefine::ClassId() const {
    return Hash("TPZBndCond::TPZ_BCDefine");
}

void TPZBndCond::TPZ_BCDefine::Read(TPZStream& buf, void* context) {
    fBCVal2.Read(buf, context);
    fForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fExactSol = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
}

void TPZBndCond::TPZ_BCDefine::Write(TPZStream& buf, int withclassid) const {
    fBCVal2.Write(buf, withclassid);
    TPZPersistenceManager::WritePointer(fForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fExactSol.operator ->(), &buf);
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
	TPZMaterial *mat = dynamic_cast<TPZMaterial *>(this->fMaterial);
	
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
    return Hash("TPZBndCond") ^ TPZMaterial::ClassId() << 1;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZBndCond>;
#endif

void TPZBndCond::Write(TPZStream &buf, int withclassid) const {
    TPZMaterial::Write(buf, withclassid);
    buf.Write(fBCs);
    buf.Write(&fType);
    fBCVal1.Write(buf, withclassid);
    fBCVal2.Write(buf, withclassid);
    TPZPersistenceManager::WritePointer(fMaterial, &buf);
}

void TPZBndCond::Read(TPZStream &buf, void *context){
    TPZMaterial::Read(buf, context);
    buf.Read(fBCs);
    buf.Read(&fType);
    fBCVal1.Read(buf, context);
    fBCVal2.Read(buf, context);
    fMaterial = dynamic_cast<TPZMaterial *>(TPZPersistenceManager::GetInstance(&buf));
    if (!fMaterial) {
        std::cout << " reading a boundary condition without material object!!\n";
#ifdef PZ_LOG
        LOGPZ_FATAL(logger, "reading a boundary condition without material object!!");
#endif
    }
}

void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    this->fMaterial->ContributeBC(data, weight, ek, ef, *this);
}

//void TPZBndCond::ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc){
//    DebugStop();
//}
//----

void TPZBndCond::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

    
	this->fMaterial->ContributeBC(datavec,weight,ek,ef,*this);
}
//----

void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
	this->fMaterial->ContributeBC(data,weight,ef,*this);
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	TPZMaterial *mat = dynamic_cast<TPZMaterial *>(fMaterial);
	if(!mat) DebugStop();

	
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

void TPZBndCond::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    TPZMaterial *mat = dynamic_cast<TPZMaterial *>(fMaterial);
    if(!mat) DebugStop();
	
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
	TPZMaterial *mat = dynamic_cast<TPZMaterial *>(fMaterial);
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

void TPZBndCond::FillDataRequirementsInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavec_left, std::map<int, TPZMaterialData> &datavec_right) {
    data.fNeedsNormal = true;
    int nref_left = datavec_left.size();
    for(auto &it : datavec_left)
    {
        it.second.SetAllRequirements(false);
        it.second.fNeedsSol = true;
        it.second.fNeedsNormal = true;
    }
    int nref_right = datavec_right.size();
    for(auto &it : datavec_right)
    {
        it.second.SetAllRequirements(false);
        it.second.fNeedsSol = true;
        it.second.fNeedsNormal = true;
    }
}
