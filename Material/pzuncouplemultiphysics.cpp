/**
 * @file
 * @brief Contains implementations of the TPZUncoupledMultiPhysics methods.
 */

#include "pzuncoupledmultiphysics.h"
#include "pzmaterialdata.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzadmchunk.h"
#include "tpzintpoints.h"

#include "pzlog.h"

#ifdef LOG4CXX
#ifdef PZDEBUG
#define DEBUG2
#endif
static LoggerPtr logger(Logger::getLogger("pz.uncoupledmultiphysics"));
#endif


using namespace std;


TPZUncoupledMultiPhysics::TPZUncoupledMultiPhysics() {
}

TPZUncoupledMultiPhysics::TPZUncoupledMultiPhysics(int id) : TPZMaterial(id) {
}

TPZUncoupledMultiPhysics::~TPZUncoupledMultiPhysics()
{
}


TPZUncoupledMultiPhysics::TPZUncoupledMultiPhysics(const TPZUncoupledMultiPhysics &material) : TPZMaterial(material)
{
    fReferredMaterials = material.fReferredMaterials;
}


void TPZUncoupledMultiPhysics::FillDataRequirements(TPZMaterialData &data)
{
    fReferredMaterials[0]->FillDataRequirements(data);
}

void TPZUncoupledMultiPhysics::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++ )
	{
        fReferredMaterials[i]->FillDataRequirements(datavec[i]);
	}
	
}

void TPZUncoupledMultiPhysics::Print(std::ostream & out) {
	out << std::endl << "Material Id = " << fId << std::endl;
    for (int i=0; i<fReferredMaterials.size(); i++) {
        fReferredMaterials[i]->Print();
    }
}

int TPZUncoupledMultiPhysics::VariableIndex(const std::string &name) {
    for (int i=0; i<fReferredMaterials.size(); i++) {
        int varindex = fReferredMaterials[i]->VariableIndex(name);
        if (varindex != -1) {
            varindex = i*10000+varindex;
            return varindex;
        }
    }
    return -1;
}

int TPZUncoupledMultiPhysics::NSolutionVariables(int index) {
    int matid = index/10000;
    int varindex = index%10000;
    return fReferredMaterials[matid]->NSolutionVariables(varindex);
}

void TPZUncoupledMultiPhysics::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    int matindex = var/10000;
    int varindex = var%10000;
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	fReferredMaterials[matindex]->Solution(data, varindex, Solout);
}

void TPZUncoupledMultiPhysics::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	if (datavec.size()==1) {
		this->Solution(datavec[0], var, Solout);
	}
	else {
		this->Solution(datavec, var, Solout);
	}
}


TPZBndCond *TPZUncoupledMultiPhysics::CreateBC(TPZMaterial* reference, int id, int typ, TPZFMatrix<STATE> &val1, TPZFMatrix<STATE> &val2) {
    DebugStop();
	return NULL;
}

void TPZUncoupledMultiPhysics::SetData(std::istream &data) {
	PZError << "TPZUncoupledMultiPhysics::SetData is called.\n";
	data >> fId;
}

TPZUncoupledMultiPhysics * TPZUncoupledMultiPhysics::NewMaterial() {
	PZError << "TPZUncoupledMultiPhysics::NewMaterial is called.\n";
	return 0;
}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZUncoupledMultiPhysics::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition material
 * @since October 07, 2011
 */
void TPZUncoupledMultiPhysics::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
}

void TPZUncoupledMultiPhysics::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->Contribute(data, weight, fakeek, ef);
}

void TPZUncoupledMultiPhysics::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeBC(data, weight, fakeek, ef, bc);
}

void TPZUncoupledMultiPhysics::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	int nref=datavec.size();
	if (nref== 1) {
		this->Contribute(datavec[0], weight, ek,ef);
	}
}

void TPZUncoupledMultiPhysics::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, 
							   TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	int nref=datavec.size();
	if (nref== 1) {
		this->ContributeBC(datavec[0], weight, ek,ef,bc);
	}
}

void TPZUncoupledMultiPhysics::Clone(std::map<int, TPZMaterial * >&matvec) {
	int matid = Id();
	std::map<int, TPZMaterial * >::iterator matit;
	matit = matvec.find(matid);
	if(matit != matvec.end()) return;
	TPZUncoupledMultiPhysics * newmat = NewMaterial();
	newmat->SetForcingFunction(TPZMaterial::fForcingFunction);
	matvec[matid] = newmat;
}

/** Get the order of the integration rule necessary to integrate an
 * element with polinomial order p */
int TPZUncoupledMultiPhysics::IntegrationRuleOrder(int elPMaxOrder) const
{   
    return fReferredMaterials[0]->IntegrationRuleOrder(elPMaxOrder);
    
}

int TPZUncoupledMultiPhysics::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const
{
    int order = 0;
    for (int i=0; i<fReferredMaterials.size(); i++) {
        int reforder = fReferredMaterials[i]->IntegrationRuleOrder(elPMaxOrder[i]);
        order = order < reforder? reforder : order;
    }
    return order;
}

int TPZUncoupledMultiPhysics::ClassId() const{
    return Hash("TPZUncoupledMultiPhysics") ^ TPZMaterial::ClassId() << 1;
}

/* Saves the element data to a stream */
void TPZUncoupledMultiPhysics::Write(TPZStream &buf, int withclassid) const
{
	TPZMaterial::Write(buf,withclassid);
}

/* Reads the element data from a stream */
void TPZUncoupledMultiPhysics::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
}

