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

int TPZBndCond::ClassId() const
{
	return TPZBNDCONDID;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZBndCond,TPZBNDCONDID>;
#endif

void TPZBndCond::Write(TPZStream &buf, int withclassid)
{
	TPZMaterial::Write(buf, withclassid);
	buf.Write(&fType, 1);
	fBCVal1.Write(buf, 0);
	fBCVal2.Write(buf, 0);
	int MatId =fMaterial->Id();
	buf.Write(&MatId, 1);
}

void TPZBndCond::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf, context);
	buf.Read(&fType, 1);
	fBCVal1.Read(buf, 0);
	fBCVal2.Read(buf, 0);
	int MatId;
	buf.Read(&MatId,1);
	TPZCompMesh * pCM = (TPZCompMesh * )/*dynamic_cast<TPZCompMesh *>*/(context);
	fMaterial = pCM->FindMaterial(MatId);
	if(!fMaterial)
	{
		std::cout << " reading a boundary condition without material object!!\n";
#ifdef LOG4CXX
		LOGPZ_FATAL(logger,"reading a boundary condition without material object!!");
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
	
    this->UpdateBCValues( data );
	
	if(dataleft.sol.NElements() < dataright.sol.NElements()){
		//		data.InvertLeftRightData();
        for(int i=0; i<3; i++) data.normal[i] *= -1.;
		mat->ContributeInterfaceBCErrors(data,dataright,weight,nkR,*this, errorid);
	}
	else {
		mat->ContributeInterfaceBCErrors(data,dataleft,weight,nkL,*this, errorid);
	}
	
}


void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
  TPZBndCond copy(*this);
	copy.UpdateBCValues(data);
	int numbersol = data.sol.size();
	//clone meshes required analysis
	int typetmp = copy.fType;
	if (copy.fType == 50){
		int i;
#ifdef DEBUG2
		{
			std::stringstream sout;
			sout << __PRETTY_FUNCTION__ << data.sol << " " << data.x;
			LOGPZ_DEBUG(logger,sout.str().c_str());
		}
#endif
		for (i = 0; i <data.sol.NElements(); i++){
      for (int is=0; is<numbersol; is++) {
        copy.fBCVal2(i,0) = (STATE)gBigNumber*data.sol[is][i];
      }
			copy.fBCVal1(i,i) = gBigNumber;
		}
		copy.fType = 2;
	}
	
	this->fMaterial->ContributeBC(data,weight,ek,ef,copy);
	copy.fType = typetmp;
}

//----
void TPZBndCond::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	
  TPZBndCond copy(*this);
	int typetmp = copy.fType;
	if (fType == 50) {
				int i;
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
		{
			for(int iref=0; iref < datavec.size(); iref++){
				std::stringstream sout;
				sout << __PRETTY_FUNCTION__ << datavec[iref].sol << " " << datavec[iref].x;
				LOGPZ_DEBUG(logger,sout.str().c_str());
			}
		}
#endif
		for (i = 0; i <datavec[0].sol[0].NElements(); i++){
					copy.fBCVal2(i,0) = ((STATE)gBigNumber)*datavec[0].sol[0][i];
					copy.fBCVal1(i,i) = ((STATE)gBigNumber);
        }
        copy.fType = 2;
	}
	this->fMaterial->ContributeBC(datavec,weight,ek,ef,copy);
	copy.fType = typetmp;
}
//----

void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
    TPZBndCond copy(*this);
	copy.UpdateBCValues(data);
	int numbersol = data.sol.size();
	//clone meshes required analysis
	int typetmp = copy.fType;
	if (copy.fType == 50){
		int i;
#ifdef DEBUG2
		{
			std::stringstream sout;
			sout << __PRETTY_FUNCTION__ << data.sol << " " << data.x;
			LOGPZ_DEBUG(logger,sout.str().c_str());
		}
#endif
		for (i = 0; i <data.sol.NElements(); i++){
            for (int is=0; is<numbersol; is++) {
                copy.fBCVal2(i,0) = (STATE)gBigNumber*data.sol[is][i];
            }
			copy.fBCVal1(i,i) = gBigNumber;
		}
		copy.fType = 2;
	}
	
	this->fMaterial->ContributeBC(data,weight,ef,copy);
	copy.fType = typetmp;
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
	if(!mat) DebugStop();// return;
	this->UpdateBCValues(data);
	
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
	if(!mat) DebugStop();// return;
	//this->UpdateBCValues(data);
	
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
	if(!mat) DebugStop();//return;
	this->UpdateBCValues(data);
	
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


void TPZBndCond::UpdateBCValues(TPZMaterialData &data){
	if(fForcingFunction){
		TPZManVector<STATE> result(fBCVal2.Rows(),0.);
		fForcingFunction->Execute(data.x,result);
		int i;
		for(i=0; i<fBCVal2.Rows(); i++) {
			fBCVal2(i,0) = result[i];
		}
	}
	
	if( this->fValFunction ){
		TPZManVector<STATE> result(this->fBCVal2.Rows(),0.);
		this->fValFunction( data.x, this->fBCVal1, result, this->fType );
		int i;
		for(i = 0; i < this->fBCVal2.Rows(); i++) {
			this->fBCVal2(i,0) = result[i];
		}
	}//if
    
    int nbc = fBCs.size();
    for (int ibc=0; ibc<nbc; ibc++) {
        if (fBCs[ibc].fForcingFunction) {
            TPZManVector<STATE> result(fBCVal2.Rows(),0.);
            fBCs[ibc].fForcingFunction->Execute(data.x,result);
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
