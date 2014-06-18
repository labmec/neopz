/**
 * @file
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

TPZVec< void(*) (const TPZVec<REAL> &, TPZVec<STATE>& ) > GFORCINGVEC;

using namespace std;
REAL TPZMaterial::gBigNumber = 1.e12;


TPZMaterial::TPZMaterial() : fNumLoadCases(1), fPostProcIndex(0) {
	this->fId = -666;
	this->fForcingFunction = NULL;
	this->fLinearContext = true;
}

TPZMaterial::TPZMaterial(int id) : fId(id), fNumLoadCases(1), fPostProcIndex(0) {
	this->SetId(id);
	fForcingFunction = 0;
	this->fLinearContext = true;
}

TPZMaterial::~TPZMaterial()
{
}


TPZMaterial::TPZMaterial(const TPZMaterial &material) {
	fId = material.fId;
    fNumLoadCases = material.fNumLoadCases;
    fPostProcIndex = material.fPostProcIndex;
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
    out << __PRETTY_FUNCTION__ << std::endl;
	out << std::endl << "Material Id = " << fId << std::endl;
    out << "Linear context " << fLinearContext << std::endl;
    out << "Num loadcases " << fNumLoadCases << std::endl;
    out << "Big number " << gBigNumber << std::endl;
    
    if (!fForcingFunction) {
        out << "Has no forcing function\n";
    }
    else {
        out << "Forcing function\n";
        fForcingFunction->Print(out);
    }
    if (!fForcingFunctionExact) {
        out << "Has no exact forcing function\n";
    }
    else {
        out << "Forcing function exact\n";
        fForcingFunctionExact->Print(out);
    }
    if (!fTimeDependentForcingFunction) {
        out << "Has no time dependent forcing function\n";
    }
    else {
        out << "Time dependent forcing function\n";
        fTimeDependentForcingFunction->Print(out);
    }
    if (!fTimedependentFunctionExact) {
        out << "No time dependent forcing function exact\n";
    }
    else {
        out << "Time dependent forcing function exact\n";
        fTimedependentFunctionExact->Print(out);
    }
    
}

int TPZMaterial::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"state")) return 0;
	if(!strcmp(name.c_str(),"State")) return 0;
	if(!strcmp(name.c_str(),"Solution")) return 0;
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
#ifdef STATE_COMPLEX
	if(index == 0) return NStateVariables()*2;    
#else
	if(index == 0) return NStateVariables();
#endif
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

void TPZMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
}

void TPZMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
    this->Solution(datavec, var, Solout);
}

void TPZMaterial::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout)
{
		this->Solution(data,dataleftvec,datarightvec, var, Solout);
}

void TPZMaterial::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl *left, TPZCompEl *ritgh)
{
	this->Solution(data,dataleftvec,datarightvec, var, Solout, left, ritgh);
}

void TPZMaterial::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &/*DSol*/,TPZFMatrix<REAL> &/*axes*/,int var,
						   TPZVec<STATE> &Solout){
#ifdef STATE_COMPLEX
    if(var == 0) 
    {
        Solout[0] = Sol[0].real();
        Solout[1] = Sol[0].imag();
    }
    else if(var == 99 || var == 100 || var == 101 || var == 102) {
        PZError << "TPZMaterial var = "<< var << " the element should treat this case\n";
        Solout[0] = Sol[0].real(); // = 0.;
    } 
	else 
    {
        Solout.Resize(0);
    }
#else
    if(var == 0) Solout = Sol;
    else if(var == 99 || var == 100 || var == 101 || var == 102) {
    PZError << "TPZMaterial var = "<< var << " the element should treat this case\n";
        Solout[0] = Sol[0]; // = 0.;
    } else Solout.Resize(0);
#endif
}

TPZBndCond *TPZMaterial::CreateBC(TPZMaterial * reference, int id, int typ, TPZFMatrix<STATE> &val1, TPZFMatrix<STATE> &val2) {
	return new TPZBndCond(reference,id,typ,val1,val2);
}

void TPZMaterial::SetData(std::istream &data) {
	PZError << "TPZMaterial::SetData is called.\n";
	data >> fId;
}

TPZMaterial * TPZMaterial::NewMaterial() {
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

void TPZMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
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

void TPZMaterial::Clone(std::map<int, TPZMaterial * >&matvec) {
	int matid = Id();
	std::map<int, TPZMaterial * >::iterator matit;
	matit = matvec.find(matid);
	if(matit != matvec.end()) return;
	TPZMaterial * newmat = NewMaterial();
	newmat->SetForcingFunction(TPZMaterial::fForcingFunction);
	matvec[matid] = newmat;
}

/** Get the order of the integration rule necessary to integrate an
 * element with polinomial order p */
int TPZMaterial::IntegrationRuleOrder(int elPMaxOrder) const
{   
//    if(fForcingFunction){return 20;}
    return 2*elPMaxOrder;
    
}

int TPZMaterial::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const
{
    int order = 0;
    if(fForcingFunction){
        order = fForcingFunction->PolynomialOrder();
    }
    
	int pmax = 0;
	for (int ip=0;  ip<elPMaxOrder.size(); ip++) 
	{
		if(elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];  
	}
    if (pmax < order) {
        pmax = order;
    }
	
	return  2*pmax;
}

int TPZMaterial::ClassId() const
{
	return TPZMATERIALID;
}

/* Saves the element data to a stream */
void TPZMaterial::Write(TPZStream &buf, int withclassid)
{
    if(ClassId() == TPZMATERIALID)
    {
        DebugStop();
    }
	TPZSaveable::Write(buf,withclassid);
	buf.Write(&fId,1);
	buf.Write(&gBigNumber,1);
    if(fForcingFunction)
    {
        fForcingFunction->Write(buf, 1);
    }
    else {
        int minone = -1;
        buf.Write(&minone);
    }
    buf.Write(&fNumLoadCases);
    int linearcontext = fLinearContext;
    buf.Write(&linearcontext);
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

/* Reads the element data from a stream */
void TPZMaterial::Read(TPZStream &buf, void *context)
{
	TPZSaveable::Read(buf,context);
	buf.Read(&fId,1);
	buf.Read(&gBigNumber,1);
    TPZSaveable *sav = TPZSaveable::Restore(buf, context);
    if(sav)
    {
        TPZFunction<STATE> *func = dynamic_cast<TPZFunction<STATE> *>(sav);
        if(!func) 
        {
            DebugStop();
        }
        fForcingFunction = func;
    }
    buf.Read(&fNumLoadCases);
    int linearcontext;
    buf.Read(&linearcontext);
    if (linearcontext) {
        fLinearContext = true;
    }
    else {
        fLinearContext = false;
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

