/**
 * \file
 * @brief Contains implementations of the TPZCoupledTransportDarcyBC methods.
 */

#include "pzcoupledtransportdarcyBC.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"

TPZCoupledTransportDarcyBC::TPZCoupledTransportDarcyBC(TPZCoupledTransportDarcy *material, int id) : TPZBndCond(){
	this->SetId(id);
	this->fMaterial = material;
	this->fMaterials[0] = NULL;
	this->fMaterials[1] = NULL;
}

TPZCoupledTransportDarcyBC::~TPZCoupledTransportDarcyBC(){}

void TPZCoupledTransportDarcyBC::Contribute(TPZMaterialData &data,
                                            REAL weight,
                                            TPZFMatrix &ek,
                                            TPZFMatrix &ef){
	
	
	TPZBndCond * bc = this->GetCurrentMaterial();
	if (!bc) return;
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	this->UpdateConvectionDir(data.dsol[0]);
	bc->Contribute(data, weight, ek, ef);
}


void TPZCoupledTransportDarcyBC::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                                     REAL weight,
                                                     TPZFMatrix &ef) {
	
    int numbersol = dataleft.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	
	TPZBndCond * bc = this->GetCurrentMaterial();
	if (!bc) return;
	this->UpdateConvectionDirInterface(dataleft.dsol[0], dataright.dsol[0], dataleft.phi, dataright.phi);
	bc->ContributeInterface(data,dataleft, dataright, weight, ef);
}

void TPZCoupledTransportDarcyBC::
ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
					REAL weight,
					TPZFMatrix &ek,
					TPZFMatrix &ef){
	
    int numbersol = dataleft.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	// TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	// TPZFMatrix &phi = data.phi;
	TPZFMatrix &phiL = dataleft.phi;
	TPZFMatrix &phiR = dataright.phi;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	TPZFMatrix &dsolL=dataleft.dsol[0];
	TPZFMatrix &dsolR=dataright.dsol[0];
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	
	TPZBndCond * bc = this->GetCurrentMaterial();
	if (!bc) return;
	this->UpdateConvectionDirInterface(dsolL, dsolR, phiL, phiR);
	bc->ContributeInterface(data, dataleft, dataright, weight, ek, ef);
}

void TPZCoupledTransportDarcyBC::UpdateConvectionDir(TPZFMatrix &dsol){
	TPZCoupledTransportDarcy * mat = dynamic_cast<TPZCoupledTransportDarcy*>(this->Material().operator ->());
	if (!mat){
		PZError << __PRETTY_FUNCTION__ << " FATAL ERROR" << std::endl;
		exit(-1);
	}
	mat->UpdateConvectionDir(dsol);
}

void TPZCoupledTransportDarcyBC::UpdateConvectionDirInterface(TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                                              TPZFMatrix &phiL, TPZFMatrix &phiR){
	TPZCoupledTransportDarcy * mat = dynamic_cast<TPZCoupledTransportDarcy*>(this->Material().operator ->());
	if (!mat){
		PZError << __PRETTY_FUNCTION__ << " FATAL ERROR" << std::endl;
		exit(-1);
	}
	if (phiL.Rows()) mat->UpdateConvectionDir(dsolL);
	if (phiR.Rows()) mat->UpdateConvectionDir(dsolR);
}
