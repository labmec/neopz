/**
 * \file
 * @brief Contains the TPZBCTension class which implements a tension boundary condition.
 */

// $Id: pzbctension.h,v 1.13 2009-09-01 19:44:46 phil Exp $

#ifndef BCTENSIONHPP
#define BCTENSIONHPP

#include "pzbndcond.h"
#include "TPZPlacaOrthotropic.h"
#include "TPZMulticamadaOrtho.h"


template <class T, int N>
class TPZManVector;
class TPZInterpolatedElement;
class TPZMulticamadaOrthotropic;


/**
 * @ingroup material
 * @brief Class which implements a tension boundary condition, where the tensor is computed from a finite element analysis
 */
class TPZBCTension : public TPZBndCond {
	
	TPZMulticamadaOrthotropic *fMultCam;
	int fCamada;
	REAL fSign;
	
private:
	
	//  TPZInterpolatedElement *fReference;
	
	public :
    
    ~TPZBCTension(){}
	
	TPZBCTension(TPZAutoPointer<TPZMaterial> &material,int id,int type,TPZFMatrix<REAL> &val1,TPZFMatrix<REAL> &val2, REAL sign, TPZMulticamadaOrthotropic *mult, int camada);
	
	virtual int NFluxes(){ return Material()->NFluxes(); }
	
	int NStateVariables() { return Material()->NStateVariables(); }
	
	void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef) {

		TPZFMatrix<REAL> &dphi = data.dphix;
		// TPZFMatrix<REAL> &dphiL = data.dphixl;
		// TPZFMatrix<REAL> &dphiR = data.dphixr;
		TPZFMatrix<REAL> &phi = data.phi;
		// TPZFMatrix<REAL> &phiL = data.phil;
		// TPZFMatrix<REAL> &phiR = data.phir;
		// TPZManVector<REAL,3> &normal = data.normal;
		TPZManVector<REAL,3> &x = data.x;
		//int POrder=data.p;
		//int LeftPOrder=data.leftp;
		//int RightPOrder=data.rightp;
		TPZVec<REAL> &sol=data.sol[0];
		// TPZVec<REAL> &solL=data.soll;
		// TPZVec<REAL> &solR=data.solr;
		TPZFMatrix<REAL> &dsol=data.dsol[0];
		// TPZFMatrix<REAL> &dsolL=data.dsoll;
		// TPZFMatrix<REAL> &dsolR=data.dsolr;
		//REAL faceSize=data.HSize;
		TPZFMatrix<REAL> &jacinv = data.jacinv;
		TPZFMatrix<REAL> &axes = data.axes;
		
		
		
		int typekeep = fType;
		if(fType == 4) {
			TPZManVector<REAL,3> normal(3);
			normal[0] = fSign*(axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1));
			normal[1] = fSign*(axes(0,2)*axes(1,0)-axes(0,0)*axes(1,2));
			normal[2] = fSign*(axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0));
			TPZFNMatrix<9> tensor(3,3);
			fMultCam->Tensor(x,fCamada,tensor);
			int i,j;
			Val2().Zero();
			for(i=0; i<3; i++) {
				for(j=0; j<3; j++) {
					Val2()(i,0) += tensor(i,j)*normal[j];
				}
			}
			fType = 1;
			Material()->ContributeBC(data,weight,ek,ef,*this);
			fType = typekeep;
		} else {
			TPZMaterialData data;
			data.x = x;
			data.jacinv = jacinv;
			data.sol[0] = sol;
			data.dsol[0] = dsol;
			data.axes = axes;
			data.phi = phi;
			data.dphix = dphi;
			TPZBndCond::Contribute(data,weight,ek,ef);
		}
		
	}
	void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<REAL> &ef)
	{
		TPZBndCond::Contribute(data,weight,ef);
	}

};

#endif


