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

void TPZBndCond::Clone(std::map<int, TPZAutoPointer<TPZMaterial> > &matvec) {
	int matid = Id();
	
	TPZAutoPointer<TPZMaterial> refmaterial = Material();
	TPZAutoPointer<TPZMaterial> newrefmaterial;
	int refmatid = 0;
	if(refmaterial) {
		refmaterial->Clone(matvec);
		refmatid = refmaterial->Id();
		newrefmaterial = matvec[refmatid];
	}
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	matit = matvec.find(matid);
	if(matit == matvec.end())
	{
		TPZAutoPointer<TPZMaterial> newmat = TPZAutoPointer<TPZMaterial>(new TPZBndCond(*this, newrefmaterial));
		matvec[matid] = newmat;
	}
}

//#ifdef _AUTODIFF

/*
 void TPZBndCond::ContributeEnergy(TPZVec<REAL> &x,
 TPZVec<FADFADREAL> &sol, TPZVec<FADFADREAL> &dsol,
 FADFADREAL &U, REAL weight)
 {
 
 int typetmp = fType;
 if (fType == 50){
 int i;
 for (i=0;i<sol.NElements();i++){
 fBCVal2(i,0) = gBigNumber*sol[i].val().val();
 fBCVal1(i,i) = gBigNumber;
 }
 fType = 2;
 }
 fMaterial->ContributeBCEnergy(x,sol,U,weight,*this);
 fType = typetmp;
 
 
 }
 */
//#endif

void TPZBndCond::InterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZSolVec &rightu,TPZSolVec &jump){
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial.operator ->());
	
	if(!mat) return;
	if(fForcingFunction) {
		TPZManVector<REAL> result(fBCVal2.Rows());
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
										   TPZVec<REAL> &nkL,
										   TPZVec<REAL> &nkR,
										   int &errorid){
	
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial.operator ->());
	if(!mat) return;
	
    this->UpdataBCValues( data );
	
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
	this->UpdataBCValues(data);
	
	//clone meshes required analysis
	/*  {
	 std::stringstream sout;
	 sout << __PRETTY_FUNCTION__ << "bc type " <<  fType << " x " << data.x;
	 LOGPZ_DEBUG(logger,sout.str().c_str());
	 }*/
    int numbersol = data.sol.size();
	int typetmp = fType;
	if (fType == 50){
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
                fBCVal2(i,0) = gBigNumber*data.sol[is][i];
                fBCVal1(i,i) = gBigNumber;
            }
		}
		fType = 2;
	}
	
	this->fMaterial->ContributeBC(data,weight,ek,ef,*this);
	fType = typetmp;
}

//----
void TPZBndCond::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	
	//this->UpdataBCValues(datavec);
	
	int typetmp = fType;
	if (fType == 50) {
		//		int i;
#ifdef DEBUG2
		{
			for(int iref=0; iref < datavec.size(); iref++){
				std::stringstream sout;
				sout << __PRETTY_FUNCTION__ << datavec[iref].sol << " " << datavec[iref].x;
				LOGPZ_DEBUG(logger,sout.str().c_str());
			}
		}
#endif
		//for (i = 0; i <data.sol.NElements(); i++){
		//			fBCVal2(i,0) = gBigNumber*data.sol[i];
		//			fBCVal1(i,i) = gBigNumber;
		//		}
		//		fType = 2;
	}
	
	this->fMaterial->ContributeBC(datavec,weight,ek,ef,*this);
	fType = typetmp;
}
//----

void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
	this->UpdataBCValues(data);
	int numbersol = data.sol.size();
	//clone meshes required analysis
	int typetmp = fType;
	if (fType == 50){
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
                fBCVal2(i,0) = gBigNumber*data.sol[is][i];
            }
			fBCVal1(i,i) = gBigNumber;
		}
		fType = 2;
	}
	
	this->fMaterial->ContributeBC(data,weight,ef,*this);
	fType = typetmp;
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	DebugStop();//nothing to be done here
}

void TPZBndCond::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial.operator ->());
	if(!mat) DebugStop();// return;
	this->UpdataBCValues(data);
	
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

void TPZBndCond::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef){
	TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial.operator ->());
	if(!mat) DebugStop();//return;
	this->UpdataBCValues(data);
	
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


void TPZBndCond::UpdataBCValues(TPZMaterialData &data){
	if(fForcingFunction){
		TPZManVector<REAL> result(fBCVal2.Rows(),0.);
		fForcingFunction->Execute(data.x,result);
		int i;
		for(i=0; i<fBCVal2.Rows(); i++) {
			fBCVal2(i,0) = result[i];
		}
	}
	
	if( this->fValFunction ){
		TPZManVector<REAL> result(this->fBCVal2.Rows(),0.);
		this->fValFunction( data.x, this->fBCVal1, result, this->fType );
		int i;
		for(i = 0; i < this->fBCVal2.Rows(); i++) {
			this->fBCVal2(i,0) = result[i];
		}
	}//if
}

void TPZBndCond::UpdataBCValues(TPZVec<TPZMaterialData> &datavec){
	
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
	if(fLinearContext == false){
		data.fNeedsSol = true;
	}
}
