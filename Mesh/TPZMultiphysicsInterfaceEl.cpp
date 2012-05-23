/**
 * @file
 * @brief Contains the implementation of the Multiphysic interface methods
 * @author Agnaldo
 * @since 10/26/11.
 */

#include "TPZMultiphysicsInterfaceEl.h"
#include "pzelmat.h"
#include "pzinterpolationspace.h"
#include "pzmaterial.h"


TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement() : TPZCompEl(),fLeftElSideVec(0), fRightElSideVec(0)
{
}

TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) : TPZCompEl(mesh, ref, index),fLeftElSideVec(0), fRightElSideVec(0)
{
}

TPZMultiphysicsInterfaceElement::~TPZMultiphysicsInterfaceElement(){
}

void TPZMultiphysicsInterfaceElement::ComputeSideTransform(TPZManVector<TPZCompElSide> &Neighbor, TPZManVector<TPZTransform> &transf)
{
    TPZGeoEl *gel = Reference();
    int side = gel->NSides()-1;
    TPZGeoElSide thisside(gel,side);
    int numneigh = Neighbor.size();
    for (int in=0; in<numneigh; in++) {
        TPZGeoElSide gelside = Neighbor[in].Reference();
        if(! thisside.NeighbourExists(gelside))
        {
            DebugStop();
        }
        TPZTransform tr(thisside.Dimension());
        thisside.SideTransform3(gelside, tr);
        transf[in] = tr;
    }
}//ComputeSideTransform

/**
 * Add elements to the list of left and right elements
 */
void TPZMultiphysicsInterfaceElement::AddLeftRightEement(TPZCompElSide &leftel, TPZCompElSide &rightel, int index)
{
    int sz = fLeftElSideVec.size();
    if (index >= sz) {
        fLeftElSideVec.Resize(index+1);
        fRightElSideVec.Resize(index+1);
    }
    fLeftElSideVec[index] = leftel;
    fRightElSideVec[index] = rightel;
}

void TPZMultiphysicsInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
	TPZMaterial * material = Material();
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
		ek.Reset();
		ef.Reset();
		return;
	}
	
	InitializeElementMatrix(ek,ef);
	
	if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
	/*
	TPZVec<TPZMaterialData> datavec;
	const int nref = fElementVec.size(); 
	datavec.resize(nref);
	InitMaterialData(datavec);
	
	TPZManVector<TPZTransform> trvec;
	AffineTransform(trvec);
	
	int dim = Dimension();
	TPZAutoPointer<TPZIntPoints> intrule;
	
	TPZManVector<REAL,3> intpoint(dim,0.), intpointtemp(dim,0.);
	REAL weight = 0.;
	
	TPZVec<int> ordervec;
	ordervec.resize(nref);
	for (int iref=0;  iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
		datavec[iref].p = msp->MaxOrder();
		ordervec[iref] = datavec[iref].p; 
	}
	int order = material->IntegrationRuleOrder(ordervec);
	
	TPZGeoEl *ref = this->Reference();
	intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
	
	TPZManVector<int,3> intorder(dim,order);
	intrule->SetOrder(intorder);	
	int intrulepoints = intrule->NPoints();
	
	TPZFMatrix<REAL> jac, axe, jacInv;
	REAL detJac; 
	for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
	{		
		intrule->Point(int_ind,intpointtemp,weight);
		ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
		weight *= fabs(detJac);
		for (int iref=0; iref<fElementVec.size(); iref++)
		{			
			TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref]);
			trvec[iref].Apply(intpointtemp, intpoint);
			
			msp->ComputeShape(intpoint, datavec[iref].x, datavec[iref].jacobian, datavec[iref].axes, 
							  datavec[iref].detjac, datavec[iref].jacinv, datavec[iref].phi, datavec[iref].dphix);
			datavec[iref].intPtIndex = int_ind;
			msp->ComputeRequiredData(datavec[iref], intpoint);
		}
		material->Contribute(datavec,weight,ek.fMat,ef.fMat);
	}//loop over integratin points
	*/
	
}//CalcStiff

void TPZMultiphysicsInterfaceElement::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
	const int ncon = this->NConnects();
	int numeq = 0;
	int ic;
	
	for(ic=0; ic<ncon; ic++)
	{
		numeq += Connect(ic).NDof(*Mesh());
	}
	
	int nref = this->fLeftElSideVec.size();
	int nstate;
	//nstate=1;
    int numloadcases = 1;
	for (int iref=0; iref<nref; iref++) {
		
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fLeftElSideVec[iref].Element());
        TPZMaterial *mat = msp->Material();
		nstate += mat->NStateVariables();
        numloadcases = mat->NumLoadCases();        
		msp  = dynamic_cast <TPZInterpolationSpace *>(fRightElSideVec[iref].Element());
        mat = msp->Material();
		nstate += mat->NStateVariables();
        numloadcases = mat->NumLoadCases();        
	}
	
	const int numstate = nstate;
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = numstate;
	ef.fNumStateVars = numstate;
	
	int i;
	for(i=0; i<ncon; i++){
        int ndof = Connect(i).NDof(*Mesh());
#ifdef DEBUG
        TPZConnect &c = Connect(i);
        if (c.NShape()*c.NState() != ndof) {
            DebugStop();
        }
#endif
		ek.fBlock.Set(i,ndof);
		ef.fBlock.Set(i,ndof);
	}
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ek.fConnect)[i] = ConnectIndex(i);
		(ef.fConnect)[i] = ConnectIndex(i);
	}
	
}//void



