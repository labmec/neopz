/**
 * \file
 * @brief Contains implementations of the TPZMaterial methods.
 */

#include "pzmaterial.h"
#include "pzmaterialid.h"
#include "pzmaterialdata.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzadmchunk.h"
#include "tpzintpoints.h"

#include "pzlog.h"

#ifdef LOG4CXX
#ifdef DEBUG
#define DEBUG2
#endif
static LoggerPtr logger(Logger::getLogger("pz.material"));
#endif

TPZVec< void(*) ( TPZVec<REAL> &, TPZVec<STATE>& ) > GFORCINGVEC;

using namespace std;
REAL TPZMaterial::gBigNumber = 1.e12;


TPZMaterial::TPZMaterial(){
	this->fId = -666;
	this->fForcingFunction = NULL;
	this->fLinearContext = true;
}

TPZMaterial::TPZMaterial(int id) {
	this->SetId(id);
	fForcingFunction = 0;
	this->fLinearContext = true;
}

TPZMaterial::~TPZMaterial()
{
}


TPZMaterial::TPZMaterial(const TPZMaterial &material) {
	fId = material.fId;
	fForcingFunction = material.fForcingFunction;
	fLinearContext = material.fLinearContext;
}

void TPZMaterial::SetLinearContext(bool IsLinear){
	fLinearContext = IsLinear;
}

void TPZMaterial::FillDataRequirements(TPZMaterialData &data)
{
	data.SetAllRequirements(true);
	data.fNeedsNeighborSol = false;
	data.fNeedsNeighborCenter = false;
	data.fNeedsNormal = false;
}

void TPZMaterial::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++ )
	{
		datavec[i].SetAllRequirements(true);
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = false;
	}
	
}

void TPZMaterial::Print(std::ostream & out) {
	out << std::endl << "Material Id = " << fId << std::endl;
}

int TPZMaterial::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"state")) return 0;
	if(!strcmp(name.c_str(),"POrder")) return 99;
	if(!strcmp(name.c_str(),"Error")) return 100;
	if(!strcmp(name.c_str(),"TrueError")) return 101;
	if(!strcmp(name.c_str(),"EffectivityIndex")) return 102;
	
	if(!strcmp(name.c_str(),"L2Error")) return 103;
	if(!strcmp(name.c_str(),"SemiH1Error")) return 104;
	if(!strcmp(name.c_str(),"H1Error")) return 105;
	
	if(!strcmp(name.c_str(),"L2ErrorPerArea")) return 106;
	if(!strcmp(name.c_str(),"SemiH1ErrorPerArea")) return 107;
	if(!strcmp(name.c_str(),"H1ErrorPerArea")) return 108;
	if(!strcmp(name.c_str(),"dudxErrorPerArea")) return 109;
	if(!strcmp(name.c_str(),"dudyErrorPerArea")) return 110;
	if(!strcmp(name.c_str(),"ContDisc")) return 111;
	
	
	std::cout << __PRETTY_FUNCTION__ << " Variable " << name << " not found\n";
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Variable " << name << " not found";
		LOGPZ_ERROR(logger,sout.str())
	}
#endif
	DebugStop();
	return -1;
}

int TPZMaterial::NSolutionVariables(int index) {
	if(index == 0) return NStateVariables();
	if(index == 99) return 1;
	if(index == 100) return 1;
	if(index == 101) return 1;
	if(index == 102) return 1;
	if (index == 103) return 1;
	if (index == 104) return 1;
	if (index == 105) return 1;
	if (index == 106) return 1;
	if (index == 107) return 1;
	if (index == 108) return 1;
	if (index == 109) return 1;
	if (index == 110) return 1;
	if (index == 111) return 1;
	PZError << "TPZMaterial::NSolutionVariables called index = " << index << "\n";
    DebugStop();
	return 0;
}

void TPZMaterial::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout){
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
}

void TPZMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
	if (datavec.size()==1) {
		this->Solution(datavec[0], var, Solout);
	}
	else {
		this->Solution(datavec, var, Solout);
	}
}

void TPZMaterial::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &/*DSol*/,TPZFMatrix<REAL> &/*axes*/,int var,
						   TPZVec<REAL> &Solout){
	DebugStop();
	//	if(var == 0) Solout = Sol;
	//	else if(var == 99 || var == 100 || var == 101 || var == 102) {
	//  	PZError << "TPZMaterial var = "<< var << " the element should treat this case\n";
	//		Solout[0] = Sol[0]; // = 0.;
	//	} else Solout.Resize(0);
}

TPZBndCond *TPZMaterial::CreateBC(TPZAutoPointer<TPZMaterial> &reference, int id, int typ, TPZFMatrix<STATE> &val1, TPZFMatrix<STATE> &val2) {
	return new TPZBndCond(reference,id,typ,val1,val2);
}

void TPZMaterial::SetData(std::istream &data) {
	PZError << "TPZMaterial::SetData is called.\n";
	data >> fId;
}

TPZAutoPointer<TPZMaterial> TPZMaterial::NewMaterial() {
	PZError << "TPZMaterial::NewMaterial is called.\n";
	return 0;
}

void TPZMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->Contribute(data, weight, fakeek, ef);
}

void TPZMaterial::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeBC(data, weight, fakeek, ef, bc);
}

void TPZMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	int nref=datavec.size();
	if (nref== 1) {
		this->Contribute(datavec[0], weight, ek,ef);
	}
}

void TPZMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, 
							   TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	int nref=datavec.size();
	if (nref== 1) {
		this->ContributeBC(datavec[0], weight, ek,ef,bc);
	}
}

void TPZMaterial::Clone(std::map<int, TPZAutoPointer<TPZMaterial> >&matvec) {
	int matid = Id();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	matit = matvec.find(matid);
	if(matit != matvec.end()) return;
	TPZAutoPointer<TPZMaterial> newmat = NewMaterial();
	newmat->SetForcingFunction(TPZMaterial::fForcingFunction);
	matvec[matid] = newmat;
}

/* 
 * Get the order of the integration rule necessary to integrate an
 * element with polinomial order p
 */
int TPZMaterial::IntegrationRuleOrder(int elPMaxOrder) const
{
    return 2*elPMaxOrder;
}


int TPZMaterial::IntegrationRuleOrder(TPZVec<int> elPMaxOrder) const
{
	int pmax = 0;
	for (int ip=0;  ip<elPMaxOrder.size(); ip++) 
	{
		if(elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];  
	}
	
	return  2*pmax;
}

/*
 void TPZMaterial::SetIntegrationRule(TPZAutoPointer<TPZIntPoints> rule,
 int elPMaxOrder,
 int elDimension){
 TPZManVector<int,3> p2(elDimension,2*elPMaxOrder);
 rule->SetOrder(p2);
 if(this->HasForcingFunction()) {
 TPZManVector<int,3> order(elDimension,rule->GetMaxOrder());
 rule->SetOrder(order);
 }
 }
 */

//#ifdef _AUTODIFF

//void TPZMaterial::ContributeEnergy(TPZVec<REAL> &x,
//	TPZVec<FADFADREAL> &sol, TPZVec<FADFADREAL> &dsol,
//	FADFADREAL &U, REAL weight)
//{
//	PZError << "\nEnergy Contribution not implemented\n";
//}

//void TPZMaterial::ContributeBCEnergy(TPZVec<REAL> & x,
//	TPZVec<FADFADREAL> & sol, FADFADREAL &U,
//	REAL weight, TPZBndCond &bc)
//{
//	PZError << "\nBC Energy Contribution not implemented\n";
//}

//#endif

int TPZMaterial::ClassId() const
{
	return TPZMATERIALID;
}

/*
 Saves the element data to a stream
 */
void TPZMaterial::Write(TPZStream &buf, int withclassid)
{
    if(ClassId() == TPZMATERIALID)
    {
        DebugStop();
    }
	TPZSaveable::Write(buf,withclassid);
	buf.Write(&fId,1);
	buf.Write(&gBigNumber,1);
    fForcingFunction->Write(buf, 1);
    /*
	 int forcingIdx = -1;
	 if (fForcingFunction)
	 {
	 for (forcingIdx=0;forcingIdx<GFORCINGVEC.NElements();forcingIdx++)
	 {
	 if (GFORCINGVEC[ forcingIdx ] == fForcingFunction) break;
	 }
	 if ( forcingIdx == GFORCINGVEC.NElements() ) forcingIdx = -1;
	 }
	 #ifdef DEBUG2
	 {
	 std::stringstream sout;
	 sout << __PRETTY_FUNCTION__ << " writing forcing function index " << forcingIdx;
	 LOGPZ_DEBUG( logger,sout.str().c_str() );
	 }
	 #endif
	 buf.Write( &forcingIdx,1 );
     */
}

/*
 Reads the element data from a stream
 */
void TPZMaterial::Read(TPZStream &buf, void *context)
{
	TPZSaveable::Read(buf,context);
	buf.Read(&fId,1);
	buf.Read(&gBigNumber,1);
    TPZSaveable *sav = TPZSaveable::Restore(buf, context);
    {
        TPZFunction *func = dynamic_cast<TPZFunction *>(sav);
        if(!func) 
        {
            DebugStop();
        }
        fForcingFunction = func;
    }
    /*
	 int forcingIdx = -1;
	 buf.Read( &forcingIdx,1 );
	 #ifdef DEBUG2
	 {
	 std::stringstream sout;
	 sout << " Read forcing function index " << forcingIdx;
	 LOGPZ_DEBUG( logger,sout.str().c_str() );
	 }
	 #endif
	 
	 if ( forcingIdx > -1 && forcingIdx < GFORCINGVEC.NElements() )
	 {
	 fForcingFunction = GFORCINGVEC[ forcingIdx ] ;
	 #ifdef DEBUG2
	 {
	 std::stringstream sout;
	 sout << " Seting forcing function index " << forcingIdx;
	 LOGPZ_DEBUG( logger,sout.str().c_str() );
	 }
	 #endif
	 }
	 */
}

