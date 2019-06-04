/**
 * @file
 * @brief Contains the implementation of the Multiphysic interface methods
 * @author Agnaldo
 * @since 10/26/11.
 */

#include "TPZMultiphysicsInterfaceEl.h"
#include "pzelmat.h"
#include "pzinterpolationspace.h"
#include "TPZMaterial.h"
#include "pzmultiphysicselement.h"
#include "tpzintpoints.h"
#include "pzdiscgal.h"
#include "pzmultiphysicscompel.h"
#include "pzgeoel.h"

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
#include "TPZMultiphysicsCompMesh.h"


TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement() : TPZRegisterClassId(&TPZMultiphysicsInterfaceElement::ClassId),
TPZCompEl(),fLeftElSide(0), fRightElSide(0)
{
}

TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index,
                                                                    TPZCompElSide leftside, TPZCompElSide rightside) : 
TPZRegisterClassId(&TPZMultiphysicsInterfaceElement::ClassId),TPZCompEl(mesh, ref, index)
{
	
	ref->SetReference(this);
#ifdef PZDEBUG
    TPZMaterial *mat = mesh.FindMaterial(ref->MaterialId());
    if (!mat) {
        DebugStop();
    }
#endif
	ref->IncrementNumInterfaces();
    
    TPZMultiphysicsElement * mp_left = dynamic_cast<TPZMultiphysicsElement * >(leftside.Element());
    TPZMultiphysicsElement * mp_right = dynamic_cast<TPZMultiphysicsElement * >(rightside.Element());
    
#ifdef PZDEBUG
    if (!mp_left || !mp_right) {
        DebugStop();
    }
#endif
    
    int left_n_meshes = mp_left->NMeshes();
    int right_n_meshes = mp_right->NMeshes();
    
    fLeftElIndices.Resize(left_n_meshes);
    fRightElIndices.Resize(right_n_meshes);
    
    for (int iref = 0; iref < left_n_meshes; iref++) {
        fLeftElIndices[iref] = iref;
    }

    for (int iref = 0; iref < right_n_meshes; iref++) {
        fRightElIndices[iref] = iref;
    }
    
	this->SetLeftRightElement(leftside, rightside);
	this->IncrementElConnected();
    this->CreateIntegrationRule();
}

void TPZMultiphysicsInterfaceElement::IncrementElConnected(){
	const int ncon = this->NConnects();
	for(int i = 0; i < ncon; i++){
		int64_t index = this->ConnectIndex(i);
		fMesh->ConnectVec()[index].IncrementElConnected();
	}
}

/** @brief create a copy of the given element */
TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy) : 
TPZRegisterClassId(&TPZMultiphysicsInterfaceElement::ClassId),TPZCompEl(mesh,copy)
{
    TPZCompElSide left = copy.Left();
    TPZCompElSide right = copy.Right();
    if (!left || !right ) {
        DebugStop();
    }
    int leftside = left.Side();
    int rightside = right.Side();
    int64_t leftindex = left.Element()->Index();
    int64_t rightindex = right.Element()->Index();
    TPZCompEl *leftel = mesh.ElementVec()[leftindex];
    TPZCompEl *rightel = mesh.ElementVec()[rightindex];
    if (!leftel || !rightel) {
        DebugStop();
    }
    fLeftElSide = TPZCompElSide(leftel,leftside);
    fRightElSide = TPZCompElSide(rightel,rightside);
    SetLeftRightElement(fLeftElSide, fRightElSide);
}

/** @brief create a copy of the given element using index mapping */
TPZMultiphysicsInterfaceElement::TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy, std::map<int64_t,int64_t> & gl2lcConMap,
                                                                 std::map<int64_t,int64_t> & gl2lcElMap) : 
TPZRegisterClassId(&TPZMultiphysicsInterfaceElement::ClassId),TPZCompEl(mesh,copy,gl2lcElMap)
{
    /// constructor not implemented right
    DebugStop();
    TPZCompElSide left = copy.Left();
    TPZCompElSide right = copy.Right();
    if (!left || !right ) {
        DebugStop();
    }
    int leftside = left.Side();
    int rightside = right.Side();
    int64_t leftindex = left.Element()->Index();
    int64_t rightindex = right.Element()->Index();
    if (gl2lcElMap.find(leftindex) == gl2lcElMap.end() || gl2lcElMap.find(rightindex) == gl2lcElMap.end()) {
        DebugStop();
    }
    TPZCompEl *leftel = mesh.ElementVec()[gl2lcElMap[leftindex]];
    TPZCompEl *rightel = mesh.ElementVec()[gl2lcElMap[rightindex]];
    if (!leftel || !rightel) {
        DebugStop();
    }
    SetLeftRightElement(TPZCompElSide(leftel,leftside), TPZCompElSide(rightel,rightside));
    
}




TPZMultiphysicsInterfaceElement::~TPZMultiphysicsInterfaceElement(){
    if (Reference()) {
        Reference()->ResetReference();
    }
}

void TPZMultiphysicsInterfaceElement::ComputeSideTransform(TPZManVector<TPZCompElSide> &Neighbor, TPZManVector<TPZTransform<> > &transf)
{
    TPZGeoEl *gel = Reference();
    int side = gel->NSides()-1;
    TPZGeoElSide thisside(gel,side);
    int64_t numneigh = Neighbor.size();
    for (int64_t in=0; in<numneigh; in++) {
        TPZGeoElSide gelside = Neighbor[in].Reference();
        if(! thisside.NeighbourExists(gelside))
        {
            DebugStop();
        }
        TPZTransform<> tr(thisside.Dimension());
        thisside.SideTransform3(gelside, tr);
        transf[in] = tr;
    }
}//ComputeSideTransform

/**
 * Add elements to the list of left and right elements
 */
void TPZMultiphysicsInterfaceElement::SetLeftRightElement(const TPZCompElSide &leftel, const TPZCompElSide &rightel)
{
    fLeftElSide = leftel;
    fRightElSide = rightel;
    this->SetLeftRightElementIndices(fLeftElIndices, fRightElIndices);
}

/**
 * Set indices to the list of left and right elements
 */
void TPZMultiphysicsInterfaceElement::SetLeftRightElementIndices(const TPZVec<int64_t> &leftindices, const TPZVec<int64_t> &rightindices)
{
    if(! fLeftElSide || ! fRightElSide){
        DebugStop();
    };
    
    fLeftElIndices=leftindices;
    fRightElIndices=rightindices;
    TPZMultiphysicsElement *LeftEl = dynamic_cast<TPZMultiphysicsElement*>(fLeftElSide.Element());
    TPZMultiphysicsElement *RightEl = dynamic_cast<TPZMultiphysicsElement*>(fRightElSide.Element());

    int64_t nleftmeshes = LeftEl->NMeshes();
    int64_t nrightmeshes = RightEl->NMeshes();

    int64_t nleftindices = leftindices.size();
    int64_t nrightindices = rightindices.size();

    //Number of connects in each element
    TPZManVector<int64_t,5> nclvec(nleftmeshes,0);
    TPZManVector<int64_t,5> ncrvec(nrightmeshes,0);
    TPZManVector<int64_t,5> first_left_c_index(nleftmeshes + 1,0);
    TPZManVector<int64_t,5> first_right_c_index(nrightmeshes + 1,0);

    int64_t ncl=0, ncr=0;
    
    //left side
    for (int64_t iref = 0; iref<nleftmeshes; iref++) {
        TPZCompEl *Left = LeftEl->Element(iref);
        if(Left){
            nclvec[iref] = Left->NConnects();
            bool is_active = LeftEl->IsActiveApproxSpaces(iref);
            if (is_active) {
                first_left_c_index[iref+1] = first_left_c_index[iref] + nclvec[iref];
            }
        }
    }
    
    //right side
    for (int64_t iref = 0; iref<nrightmeshes; iref++) {
        TPZCompEl *Right = RightEl->Element(iref);
        if(Right){
            ncrvec[iref] = Right->NConnects();
            bool is_active = RightEl->IsActiveApproxSpaces(iref);
            if (is_active) {
                first_right_c_index[iref+1] = first_right_c_index[iref] + ncrvec[iref];
            }
        }
    }
    
    //left side
    for(int64_t i = 0; i < nleftindices; i++){
        int iref = leftindices[i];
        TPZCompEl *Left = LeftEl->Element(iref);
        if(Left){
            if (LeftEl->IsActiveApproxSpaces(iref)) {
                ncl += nclvec[iref];
            }
        }
    }
    
    //right side
    for(int64_t i = 0; i < nrightindices; i++){
        int iref = rightindices[i];
        TPZCompEl *Right = RightEl->Element(iref);
        if(Right){
            if (RightEl->IsActiveApproxSpaces(iref)) {
                ncr += ncrvec[iref];
            }
        }
    }
    
    /// Considering active elements
    fConnectIndexes.Resize(ncl+ncr);
    int64_t count = 0;
    for(int64_t i = 0; i < nleftindices; i++){
        int iref = leftindices[i];
        bool is_active = LeftEl->IsActiveApproxSpaces(iref);
        if (is_active) {
            for (int ic=0; ic < nclvec[iref]; ic++) {
                fConnectIndexes[count++] = LeftEl->ConnectIndex(first_left_c_index[iref]+ic);
            }
        }
    }
    
    for(int64_t i = 0; i < nrightindices; i++){
        int iref = rightindices[i];
        bool is_active = RightEl->IsActiveApproxSpaces(iref);
        if (is_active) {
            for (int ic=0; ic < ncrvec[iref]; ic++) {
                fConnectIndexes[count++] = RightEl->ConnectIndex(first_right_c_index[iref]+ic);
            }
        }
    }
    
    if (count != fConnectIndexes.size() ) {
        DebugStop();
    }
}



/**
 * Get left and right elements
 */
void TPZMultiphysicsInterfaceElement::GetLeftRightElement(TPZCompElSide &leftel, TPZCompElSide &rightel)
{
	leftel = fLeftElSide;
	rightel = fRightElSide;
}

/** @brief Returns the number of nodes of the element */
int TPZMultiphysicsInterfaceElement::NConnects() const
{
    return fConnectIndexes.size();
}

/**
 * @brief Returns the index of the ith connectivity of the element
 * @param i connectivity index who want knows
 */
int64_t TPZMultiphysicsInterfaceElement::ConnectIndex(int i) const
{

#ifdef PZDEBUG
    if (i < 0 || i >= fConnectIndexes.size()) {
        DebugStop();
    }
#endif
    return fConnectIndexes[i];
}


#include "pzmultiphysicscompel.h"
void TPZMultiphysicsInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    
	TPZDiscontinuousGalerkin  * material = dynamic_cast<TPZDiscontinuousGalerkin *> (this->Material());
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
    TPZGeoEl *leftgel = leftel->Reference();
    TPZGeoEl *rightgel = rightel->Reference();
#ifdef PZDEBUG
    if (!leftel || !rightel) {
        DebugStop();
    }
#endif
       
    TPZManVector<TPZMaterialData,6> datavecleft,datavecright;
    TPZMaterialData data;
    InitMaterialData(data, datavecleft, datavecright);
    
    TPZManVector<TPZTransform<REAL>,6> leftcomptr, rightcomptr;
    leftel->AffineTransform(leftcomptr);
    rightel->AffineTransform(rightcomptr);

    int nmesh =datavecleft.size();
    for(int id = 0; id<nmesh; id++){
        datavecleft[id].p = 0;
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(leftel->Element(id));
        if (msp)
        {
            datavecleft[id].p =msp->MaxOrder();
        }
    }

    TPZManVector<int> intleftorder;
    leftel->PolynomialOrder(intleftorder);
    TPZManVector<int> intrightorder;
    rightel->PolynomialOrder(intrightorder);
    int integrationorder = material->GetIntegrationOrder(intleftorder, intrightorder);
    TPZGeoEl *gel = Reference();
    int dimension = gel->Dimension();
    int thisside = gel->NSides()-1;
    TPZFNMatrix<9,REAL> jac(dimension,dimension),axes(dimension,3), jacInv(dimension,dimension);
    
    TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(thisside, integrationorder);
    TPZManVector<REAL,3> Point(dimension), leftPoint(leftel->Dimension()), rightPoint(rightel->Dimension());
    TPZGeoElSide neighleft(fLeftElSide.Reference()), neighright(fRightElSide.Reference());
    TPZTransform<> trleft(dimension),trright(dimension);
    TPZGeoElSide gelside(this->Reference(),thisside);
    // compute the transformation between neighbours
    gelside.SideTransform3(neighleft, trleft);
    gelside.SideTransform3(neighright, trright);
    
    TPZTransform<> leftloctr = leftgel->SideToSideTransform(neighleft.Side(), leftgel->NSides()-1);
    TPZTransform<> rightloctr = rightgel->SideToSideTransform(neighright.Side(), rightgel->NSides()-1);
    // transform from the element to the interior of the neighbours
    trleft = leftloctr.Multiply(trleft);
    trright = rightloctr.Multiply(trright);
    
    
    int nintpoints = intrule->NPoints();
    for (int ip =0; ip<nintpoints; ip++) {
        REAL weight;
        data.intLocPtIndex = ip;
        intrule->Point(ip, Point, weight);
        ComputeRequiredData(data, Point);
        weight *= fabs(data.detjac);
        trleft.Apply(Point, leftPoint);
        leftel->ComputeRequiredData(leftPoint, leftcomptr, datavecleft, fLeftElIndices);
        trright.Apply(Point, rightPoint);
        rightel->ComputeRequiredData(rightPoint, rightcomptr, datavecright, fRightElIndices);
        
        data.x = datavecleft[0].x;
        material->ContributeInterface(data, datavecleft, datavecright, weight, ek.fMat, ef.fMat);
    }	
	
}//CalcStiff

void TPZMultiphysicsInterfaceElement::CalcStiff(TPZElementMatrix &ef)
{
    TPZDiscontinuousGalerkin  * material = dynamic_cast<TPZDiscontinuousGalerkin *> (this->Material());
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        ef.Reset();
        return;
    }
    
    InitializeElementMatrix(ef);
    
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    TPZGeoEl *leftgel = leftel->Reference();
    TPZGeoEl *rightgel = rightel->Reference();
#ifdef PZDEBUG
    if (!leftel || !rightel) {
        DebugStop();
    }
#endif
    
    TPZManVector<TPZMaterialData,6> datavecleft,datavecright;
    TPZMaterialData data;
    InitMaterialData(data, datavecleft, datavecright);
    
    TPZManVector<TPZTransform<> > leftcomptr, rightcomptr;
    leftel->AffineTransform(leftcomptr);
    rightel->AffineTransform(rightcomptr);
    
    InitMaterialData(data);
    int nmesh =datavecleft.size();
    for(int id = 0; id<nmesh; id++){
        datavecleft[id].fNeedsNormal=true;
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(leftel->Element(id));
        datavecleft[id].p =msp->MaxOrder();
    }
    data.fNeedsHSize=true;
    
    TPZManVector<int> intleftorder;
    leftel->PolynomialOrder(intleftorder);
    TPZManVector<int> intrightorder;
    rightel->PolynomialOrder(intrightorder);
    int integrationorder = material->GetIntegrationOrder(intleftorder, intrightorder);
    TPZGeoEl *gel = Reference();
    int dimension = gel->Dimension();
    int thisside = gel->NSides()-1;
    TPZFNMatrix<9,REAL> jac(dimension,dimension),axes(dimension,3), jacInv(dimension,dimension);
    
    TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(thisside, integrationorder);
    TPZManVector<REAL,3> Point(dimension), leftPoint(leftel->Dimension()), rightPoint(rightel->Dimension());
    TPZGeoElSide neighleft(fLeftElSide.Reference()), neighright(fRightElSide.Reference());
    TPZTransform<> trleft(dimension),trright(dimension);
    TPZGeoElSide gelside(this->Reference(),thisside);
    // compute the transformation between neighbours
    gelside.SideTransform3(neighleft, trleft);
    gelside.SideTransform3(neighright, trright);
    
    TPZTransform<> leftloctr = leftgel->SideToSideTransform(neighleft.Side(), leftgel->NSides()-1);
    TPZTransform<> rightloctr = rightgel->SideToSideTransform(neighright.Side(), rightgel->NSides()-1);
    // transform from the element to the interior of the neighbours
    trleft = leftloctr.Multiply(trleft);
    trright = rightloctr.Multiply(trright);
    
    
    int nintpoints = intrule->NPoints();
    for (int ip =0; ip<nintpoints; ip++) {
        REAL weight;
        data.intLocPtIndex = ip;
        intrule->Point(ip, Point, weight);
        ComputeRequiredData(data, Point);
        weight *= fabs(data.detjac);
        trleft.Apply(Point, leftPoint);
        leftel->ComputeRequiredData(leftPoint, leftcomptr, datavecleft, fLeftElIndices);
        trright.Apply(Point, rightPoint);
        rightel->ComputeRequiredData(rightPoint, rightcomptr, datavecright, fRightElIndices);
        
        data.x = datavecleft[0].x;
        material->ContributeInterface(data, datavecleft, datavecright, weight, ef.fMat);
    }	
    
}//CalcStiff

const TPZIntPoints & TPZMultiphysicsInterfaceElement::GetIntegrationRule() const
{
    if (!fIntegrationRule) {
        DebugStop();
    }
    return *fIntegrationRule;
}


int TPZMultiphysicsInterfaceElement::ComputeIntegrationOrder() const {

    TPZMultiphysicsElement *left = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *right = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    if (!left || !right) return -1;
    
    int left_order = left->ComputeIntegrationOrder();
    int right_order = right->ComputeIntegrationOrder();
    int int_order = MAX(left_order, right_order);
    
    return int_order;
}


void TPZMultiphysicsInterfaceElement::CreateIntegrationRule()
{
    if (fIntegrationRule) {
        delete fIntegrationRule;
    }
    
    TPZDiscontinuousGalerkin  * material = dynamic_cast<TPZDiscontinuousGalerkin *> (this->Material());
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        DebugStop();
    }
    
    
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    
#ifdef PZDEBUG
    if (!leftel || !rightel) {
        DebugStop();
    }
#endif
    
    TPZManVector<int> intleftorder;
    leftel->PolynomialOrder(intleftorder);
    TPZManVector<int> intrightorder;
    rightel->PolynomialOrder(intrightorder);
    int integrationorder = material->GetIntegrationOrder(intleftorder, intrightorder);
    TPZGeoEl *gel = Reference();
    int thisside = gel->NSides()-1;
    
    fIntegrationRule = gel->CreateSideIntegrationRule(thisside, integrationorder);
}

void TPZMultiphysicsInterfaceElement::ComputeRequiredData(TPZVec<REAL> &intpointtemp, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialData> &datavec)
{
    DebugStop();
}//ComputeRequiredData

void TPZMultiphysicsInterfaceElement::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    ek.fMesh = Mesh();
    ef.fMesh = ek.fMesh;
    ek.fType = TPZElementMatrix::EK;
    ef.fType = TPZElementMatrix::EF;
	const int ncon = this->NConnects();
	int64_t numeq = 0;
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
	
	int i;
	for(i=0; i<ncon; i++)
    {
        TPZConnect &c = Connect(i);
        int ndof = Connect(i).NShape()*c.NState();
#ifdef PZDEBUG
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

void TPZMultiphysicsInterfaceElement::InitializeElementMatrix(TPZElementMatrix &ef)
{

    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
    const int ncon = this->NConnects();
    int64_t numeq = 0;
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
    
    ef.fMat.Redim(numeq,numloadcases);
    ef.fBlock.SetNBlocks(ncon);
    
    int i;
    for(i=0; i<ncon; i++)
    {
        TPZConnect &c = Connect(i);
        int ndof = Connect(i).NShape()*c.NState();
#ifdef PZDEBUG
        if (c.NDof(*Mesh()) != ndof) {
            DebugStop();
        }
#endif
        ef.fBlock.Set(i,ndof);
    }
    ef.fConnect.Resize(ncon);
    for(i=0; i<ncon; i++){
        (ef.fConnect)[i] = ConnectIndex(i);
    }
    
}//void

void TPZMultiphysicsInterfaceElement::ComputeCenterNormal(TPZVec<REAL> &normal) const{
    
    TPZGeoEl *gel = Reference();
    int dim = gel->Dimension();
    int nsides = gel->NSides();
    TPZManVector<REAL> center(dim);
    gel->CenterPoint(nsides-1 , center);
    TPZGeoElSide gelside(gel,nsides-1);
    gelside.Normal(center, fLeftElSide.Element()->Reference(), fRightElSide.Element()->Reference(), normal);
}

void TPZMultiphysicsInterfaceElement::Print(std::ostream &out) const {
	
	TPZCompEl::Print(out);
	out << "\nInterface element : \n";
	if(!LeftElement() || !LeftElement()->Reference()) out << "\tNULL LeftElement - this is inconsistent\n";
	else{
                out << "\tLeft Computational Index: " << LeftElement()->Index() << std::endl;
		out << "\tLeft Geometric Index: " << LeftElement()->Reference()->Index() << std::endl;
		out << "\tLeft Geometric Id: " << LeftElement()->Reference()->Id() << std::endl;
		out << "\tElement Dimension " << LeftElement()->Reference()->Dimension() << std::endl;
	}
	
	if(!RightElement() || !RightElement()->Reference()) out << "\tNULL RightElement - this is inconsistent";
	else{
                out << "\tRight Computational Index: " << RightElement()->Index() << std::endl;
		out << "\tRight Geometric Index: " << RightElement()->Reference()->Index() << std::endl;
		out << "\tRight Geometric Id: " << RightElement()->Reference()->Id() << std::endl;
		out << "\tElement Dimension " << RightElement()->Reference()->Dimension() << std::endl;
	}
	
	out << "\tMaterial id : " << Reference()->MaterialId() << std::endl;
	
    TPZMaterial *mat = Material();
    TPZMaterialData data;
    mat->FillDataRequirements(data);
    if (mat && data.fNeedsNormal)
    {
        TPZVec<REAL> center_normal;
        ComputeCenterNormal(center_normal);
        out << "\tNormal vector (at center point): ";
        out << "(" << center_normal[0] << "," << center_normal[1] << "," << center_normal[2] << ")\n";
        
    }
}

///** @brief Initialize the material data for the neighbouring element */
//void TPZMultiphysicsInterfaceElement::InitMaterialData(TPZVec<TPZMaterialData> &data, TPZMultiphysicsElement *mfcel,TPZVec<int64_t> *indices)
//{
//    data.resize(mfcel->NMeshes());
//    mfcel->InitMaterialData(data,indices);
//}

/** @brief Initialize the material data structures */
void TPZMultiphysicsInterfaceElement::InitMaterialData(TPZMaterialData &center_data, TPZVec<TPZMaterialData> &data_left, TPZVec<TPZMaterialData> &data_right){

    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
  //  TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(Mesh());
    int n_meshes = leftel->NMeshes();
    data_left.resize(n_meshes);
    data_right.resize(n_meshes);
    
    TPZVec<int64_t> *leftindices(0), *rightindices(0);
    if (fLeftElIndices.size()) {
        leftindices = &fLeftElIndices;
    }
    if (fRightElIndices.size()) {
        rightindices = &fRightElIndices;
    }
    
    leftel->InitMaterialData(data_left,leftindices);
    rightel->InitMaterialData(data_right,rightindices);
    
    TPZMaterial * mat = this->Material();
    mat->FillDataRequirementsInterface(center_data, data_left, data_right);
    
    
}

/** @brief initialize the material data for the geometric data */
void TPZMultiphysicsInterfaceElement::InitMaterialData(TPZMaterialData &data)
{
    TPZGeoEl *gel = Reference();
    int dim = gel->Dimension();
    int nsides = gel->NSides();
    TPZManVector<REAL> center(dim);
    gel->CenterPoint(nsides-1 , center);
    TPZGeoElSide gelside(gel,nsides-1);
    TPZMaterial *mat = Material();
    if (mat) {
        mat->FillDataRequirements(data);
    }
    if (data.fNeedsNormal)
    {
        gelside.Normal(center, fLeftElSide.Element()->Reference(), fRightElSide.Element()->Reference(), data.normal);
    }
    data.axes.Redim(dim,3);
    data.jacobian.Redim(dim,dim);
	data.jacinv.Redim(dim,dim);
	data.x.Resize(3);
}

/** @brief Compute the data needed to compute the stiffness matrix at the integration point */
void TPZMultiphysicsInterfaceElement::ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &point)
{
    data.intGlobPtIndex = -1;
    TPZGeoEl *gel = Reference();
    TPZGeoElSide gelside(gel,gel->NSides()-1);
    gel->Jacobian(point, data.jacobian, data.axes, data.detjac, data.jacinv);
    
    TPZMaterial *mat = Material();
    if (mat) {
        mat->FillDataRequirementsInterface(data);
    }

    if (data.fNeedsNormal)
    {
        
        if (gelside.Dimension() == gelside.Element()->Dimension()-1) {
            gelside.Normal(point, data.normal);
        }else{
            gelside.Normal(point, fLeftElSide.Element()->Reference(), fRightElSide.Element()->Reference(), data.normal);
        }
    }
    
    if (data.fNeedsHSize){
		const int dim = this->Dimension();
		REAL faceSize;
		if (dim == 0){//it means I am a point
//            DebugStop();
            faceSize = 0.;
		}
		else{
			faceSize = 2.*this->Reference()->ElementRadius();//Igor Mozolevski's suggestion. It works well for elements with small aspect ratio
		}
		data.HSize = faceSize;
	}

}

void TPZMultiphysicsInterfaceElement::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)
{
	TPZGeoEl *ref = Reference();
	if (ref->Dimension() != dimension) {
		return;
	}
	
	TPZMaterial * material = Material();
    
    TPZManVector<std::string,4> scalarnames, vecnames;
    scalarnames = grmesh.ScalarNames();
    vecnames = grmesh.VecNames();
    for (int i=0; i<scalarnames.size(); i++) {
        if (material->VariableIndex(scalarnames[i]) == -1) {
            return;
        }
    }
    for (int i=0; i<vecnames.size(); i++) {
        if (material->VariableIndex(vecnames[i]) == -1) {
            return;
        }
    }
	int matid = material->Id();
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

void TPZMultiphysicsInterfaceElement::Solution(TPZVec<REAL> &qsi, int var,TPZVec<STATE> &sol)
{
	
	if(var >= 100) {
		TPZCompEl::Solution(qsi,var,sol);
		return;
	}
	
	TPZMaterial * material = this->Material();
	if(!material){
		sol.Resize(0);
		return;
	}
	
	if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
	
	TPZCompElSide LeftSide;
	TPZCompElSide RightSide;
	this->GetLeftRightElement(LeftSide, RightSide);
	
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (LeftSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(RightSide.Element());
		
    TPZManVector<TPZMaterialData,6> datavecleft,datavecright;
    TPZMaterialData data;
    InitMaterialData(data, datavecleft, datavecright);
	
    TPZManVector<TPZTransform<> > leftcomptr, rightcomptr;
    leftel->AffineTransform(leftcomptr);
    rightel->AffineTransform(rightcomptr);
    InitMaterialData(data);	
	TPZTransform<> lefttr;
	TPZTransform<> righttr;	
	
	//		Integration points in left and right elements: making transformations to neighbour elements
	this->ComputeSideTransform(LeftSide, lefttr);
	this->ComputeSideTransform(RightSide, righttr);	
	
	TPZVec<REAL> myqsi;
	myqsi.resize(qsi.size());	
	lefttr.Apply(qsi, myqsi);
	lefttr.Apply(qsi, myqsi);
	
	leftel->ComputeRequiredData(myqsi, leftcomptr, datavecleft,fLeftElIndices);
	rightel->ComputeRequiredData(myqsi, rightcomptr, datavecright,fRightElIndices);
		
	material->Solution(data,datavecleft,datavecright,var, sol,LeftSide.Element(),RightSide.Element());
}

void TPZMultiphysicsInterfaceElement::ComputeSideTransform(TPZCompElSide &Neighbor, TPZTransform<> &transf){
	TPZGeoEl * neighel = Neighbor.Element()->Reference();
	const int dim = this->Dimension();
	TPZTransform<> LocalTransf(dim);
	TPZGeoElSide thisgeoside(this->Reference(), this->Reference()->NSides()-1);
	TPZGeoElSide neighgeoside(neighel, Neighbor.Side());
	thisgeoside.SideTransform3(neighgeoside, LocalTransf);
	
	TPZGeoElSide highdim(neighel, neighel->NSides()-1);
	transf = neighgeoside.SideToSideTransform(highdim).Multiply(LocalTransf);
}//ComputeSideTransform


int TPZMultiphysicsInterfaceElement::ClassId() const{
    return Hash("TPZMultiphysicsInterfaceElement") ^ TPZCompEl::ClassId() << 1;
}
