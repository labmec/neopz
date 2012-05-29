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
#include "pzmultiphysicselement.h"
#include "tpzintpoints.h"


TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement() : TPZCompEl(),fLeftElSide(0), fRightElSide(0)
{
}

TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int &index,
                                                                    TPZCompElSide leftside, TPZCompElSide rightside) : TPZCompEl(mesh, ref, index),fLeftElSide(leftside), fRightElSide(rightside)
{
}

/** @brief create a copy of the given element */
TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy) : TPZCompEl(mesh,copy)
{
    TPZCompElSide left = copy.Left();
    TPZCompElSide right = copy.Right();
    if (!left || !right ) {
        DebugStop();
    }
    int leftside = left.Side();
    int rightside = right.Side();
    int leftindex = left.Element()->Index();
    int rightindex = right.Element()->Index();
    TPZCompEl *leftel = mesh.ElementVec()[leftindex];
    TPZCompEl *rightel = mesh.ElementVec()[rightindex];
    if (!leftel || !rightel) {
        DebugStop();
    }
    fLeftElSide = TPZCompElSide(leftel,leftside);
    fRightElSide = TPZCompElSide(rightel,rightside);
}

/** @brief create a copy of the given element using index mapping */
TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy, std::map<int,int> & gl2lcConMap,
                                std::map<int,int> & gl2lcElMap)
{
    TPZCompElSide left = copy.Left();
    TPZCompElSide right = copy.Right();
    if (!left || !right ) {
        DebugStop();
    }
    int leftside = left.Side();
    int rightside = right.Side();
    int leftindex = left.Element()->Index();
    int rightindex = right.Element()->Index();
    if (gl2lcElMap.find(leftindex) == gl2lcElMap.end() || gl2lcElMap.find(rightindex) == gl2lcElMap.end()) {
        DebugStop();
    }
    TPZCompEl *leftel = mesh.ElementVec()[gl2lcElMap[leftindex]];
    TPZCompEl *rightel = mesh.ElementVec()[gl2lcElMap[rightindex]];
    if (!leftel || !rightel) {
        DebugStop();
    }
    fLeftElSide = TPZCompElSide(leftel,leftside);
    fRightElSide = TPZCompElSide(rightel,rightside);
    
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
void TPZMultiphysicsInterfaceElement::SetLeftRightEement(TPZCompElSide &leftel, TPZCompElSide &rightel)
{
    fLeftElSide = leftel;
    fRightElSide = rightel;
}

/** @brief Returns the number of nodes of the element */
int TPZMultiphysicsInterfaceElement::NConnects() const
{
    return fLeftElSide.Element()->NConnects() + fRightElSide.Element()->NConnects();
}

/**
 * @brief Returns the index of the ith connectivity of the element
 * @param i connectivity index who want knows
 */
int TPZMultiphysicsInterfaceElement::ConnectIndex(int i) const
{
#ifdef DEBUG
    if (i < 0) {
        DebugStop();
    }
#endif
    int nleft = fLeftElSide.Element()->NConnects();
    if (i < nleft) {
        return fLeftElSide.Element()->ConnectIndex(i);
    }
    int nright = fRightElSide.Element()->NConnects();
    if (i < nleft+nright) {
        return fRightElSide.Element()->ConnectIndex(i-nleft);
    }
    DebugStop();
    return -1;
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
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
#ifdef DEBUG
    if (!leftel || !rightel) {
        DebugStop();
    }
#endif
    TPZManVector<TPZMaterialData,6> datavecleft,datavecright;
    TPZMaterialData data;
    InitMaterialData(datavecleft, leftel);
    InitMaterialData(datavecright, rightel);
    InitMaterialData(data);
    int intleftorder = leftel->IntegrationOrder();
    int intrightorder = rightel->IntegrationOrder();
    int integrationorder = MAX(intleftorder, intrightorder);
    TPZGeoEl *gel = Reference();
    int dimension = gel->Dimension();
    int thisside = gel->NSides()-1;
    TPZFNMatrix<9,REAL> jac(dimension,dimension),axes(dimension,3), jacInv(dimension,dimension);
    
    TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(thisside, integrationorder);
    TPZManVector<REAL,3> Point(dimension), leftPoint(leftel->Dimension()), rightPoint(rightel->Dimension());
    TPZGeoElSide neighleft(fLeftElSide.Reference()), neighright(fRightElSide.Reference());
    TPZTransform trleft(dimension),trright(dimension);
    TPZGeoElSide gelside(this->Reference(),thisside);
    gelside.SideTransform3(neighleft, trleft);
    gelside.SideTransform3(neighright, trright);
    
    int nintpoints = intrule->NPoints();
    for (int ip =0; ip<nintpoints; ip++) {
        REAL weight;
        intrule->Point(ip, Point, weight);
        REAL detjac;
        gel->Jacobian(Point, jac, axes, detjac, jacInv);
    }
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
        TPZConnect &c = Connect(ic);
		numeq += c.NShape()*c.NState();
	}
	
    TPZMultiphysicsElement *mfcel_left = dynamic_cast<TPZMultiphysicsElement *>(fLeftElSide.Element());
    TPZMultiphysicsElement *mfcel_right = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    if (! mfcel_left || !mfcel_right) {
        DebugStop();
    }
	//nstate=1;
    int numloadcases = 1;
    TPZMultiphysicsElement *msp  = dynamic_cast <TPZMultiphysicsElement *>(fLeftElSide.Element());
    if (!msp) {
        DebugStop();
    }
    TPZMaterial *mat = msp->Material();
    int nstate = mat->NStateVariables();
    numloadcases = mat->NumLoadCases();        
	
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = nstate;
	ef.fNumStateVars = nstate;
	
	int i;
	for(i=0; i<ncon; i++)
    {
        TPZConnect &c = Connect(i);
        int ndof = Connect(i).NShape()*c.NState();
#ifdef DEBUG
        if (c.NDof(*Mesh()) != ndof) {
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


/** @brief Initialize the material data for the neighbouring element */
void TPZMultiphysicsInterfaceElement::InitMaterialData(TPZVec<TPZMaterialData> &data, TPZMultiphysicsElement *mfcel)
{
    
}

/** @brief initialize the material data for the geometric data */
void TPZMultiphysicsInterfaceElement::InitMaterialData(TPZMaterialData &data)
{
    
}


