/**
 * @file
 * @brief Contains implementations of the TPZMaterial methods.
 */

#include "TPZMaterial.h"
#include "pzmaterialdata.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzadmchunk.h"
#include "tpzintpoints.h"

#include "pzlog.h"
#include "TPZPersistenceManager.h"

#ifdef LOG4CXX
#ifdef PZDEBUG
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
    this->fForcingFunctionExact = NULL;
    this->fTimeDependentForcingFunction = NULL;
    this->fTimedependentFunctionExact = NULL;
    this->fBCForcingFunction = NULL;
    this->fTimedependentBCForcingFunction = NULL;
	this->fLinearContext = true;
    this->fBCForcingFunction = NULL;
    
}

TPZMaterial::TPZMaterial(int id) : fId(id), fNumLoadCases(1), fPostProcIndex(0) {
	this->SetId(id);
    this->fForcingFunction = NULL;
    this->fForcingFunctionExact = NULL;
    this->fTimeDependentForcingFunction = NULL;
    this->fTimedependentFunctionExact = NULL;
    this->fBCForcingFunction = NULL;
    this->fTimedependentBCForcingFunction = NULL;
    this->fLinearContext = true;
    this->fBCForcingFunction = NULL;

}

TPZMaterial::~TPZMaterial()
{
    this->fId = -999;
}


TPZMaterial::TPZMaterial(const TPZMaterial &material) {
	fId = material.fId;
    fForcingFunction = material.fForcingFunction;
    fForcingFunctionExact = material.fForcingFunctionExact;
    fTimeDependentForcingFunction = material.fTimeDependentForcingFunction;
    fTimedependentFunctionExact = material.fTimedependentFunctionExact;
    fBCForcingFunction = material.fBCForcingFunction;
    fTimedependentBCForcingFunction = material.fTimedependentBCForcingFunction;
    fLinearContext = material.fLinearContext;
    fNumLoadCases = material.fNumLoadCases;
    fPostProcIndex = material.fPostProcIndex;
}

TPZMaterial &TPZMaterial::operator=(const TPZMaterial &material)
{
    fId = material.fId;
    fForcingFunction = material.fForcingFunction;
    fForcingFunctionExact = material.fForcingFunctionExact;
    fTimeDependentForcingFunction = material.fTimeDependentForcingFunction;
    fTimedependentFunctionExact = material.fTimedependentFunctionExact;
    fBCForcingFunction = material.fBCForcingFunction;
    fTimedependentBCForcingFunction = material.fTimedependentBCForcingFunction;
    fLinearContext = material.fLinearContext;
    fNumLoadCases = material.fNumLoadCases;
    fPostProcIndex = material.fPostProcIndex;
    return *this;
}

void TPZMaterial::GetSolutionDimensions(uint64_t &u_len,
                                       uint64_t &du_row,
                                       uint64_t &du_col)
{
  PZError << __PRETTY_FUNCTION__ << std::endl;
  PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
  DebugStop();
  u_len = -1; du_row = -1; du_col = -1;
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
	out << "Material Id = " << fId << std::endl;
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
    if(!strcmp(name.c_str(),"GradState")) return 1;
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
    if(!strcmp(name.c_str(),"MaterialId")) return 98;
	
	
//	std::cout << __PRETTY_FUNCTION__ << " Variable " << name << " not found\n";
	
#ifdef LOG4CXX2
	{
		std::stringstream sout;
		sout << "Variable " << name << " not found";
		LOGPZ_ERROR(logger,sout.str())
	}
#endif
	return -1;
}

int TPZMaterial::NSolutionVariables(int index) {
#ifdef STATE_COMPLEX
	if(index == 0) return NStateVariables()*2;    
#else
	if(index == 0) return NStateVariables();
#endif
	if(index == 99 || index == 98) return 1;
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
    int nvec = datavec.size();
    int numdata = 0;
    int dataindex = -1;
    for (int iv=0; iv<nvec; iv++) {
        if(datavec[iv].fShapeType != TPZMaterialData::EEmpty)
        {
            numdata++;
            dataindex = iv;
        }
    }
    if (numdata == 1) {
        Solution(datavec[dataindex], var, Solout);
        return;
    }
    DebugStop();
}

void TPZMaterial::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout)
{
		this->Solution(data,dataleftvec,datarightvec, var, Solout);
}

void TPZMaterial::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl *left, TPZCompEl *ritgh)
{
	this->Solution(data,dataleftvec,datarightvec, var, Solout, left, ritgh);
}

void TPZMaterial::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
						   TPZVec<STATE> &Solout){
    if(var == 98){
        Solout[0] = this->Id();
        return;
    }
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
    } else if(var == 1)
    {
        Solout.resize(Sol.size()*3);
        int64_t nsol = Sol.size();
        Solout.Fill(0.);
        int64_t dim = axes.Rows();
        for (int64_t is=0; is<nsol; is++) {
            for (int64_t d=0; d<dim; d++) {
                for (int64_t jco=0; jco<3; jco++) {
                    Solout[jco+3*is] += axes(d,jco)*DSol(d,is);
                }
            }
        }
    } else
    {
        DebugStop();
        Solout.Resize(0);
    }
#endif
}

TPZBndCond *TPZMaterial::CreateBC(TPZMaterial * reference, int id, int typ, const TPZFMatrix<STATE> &val1, const TPZFMatrix<STATE> &val2) {
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
    int ndif = 0;
    int onemat = 0;
    for (int ir = 0; ir < nref; ir++) {
        int nphis=datavec[ir].phi.Rows();
        if (datavec[ir].phi.Rows()) {
            onemat = ir;
            ndif++;
        }
    }
	if (ndif == 1) {
		this->Contribute(datavec[onemat], weight, ek,ef);
	}
    else
    {
        DebugStop();
    }
}

void TPZMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
							   TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	int nref=datavec.size();
	if (nref== 1) {
		this->ContributeBC(datavec[0], weight, ek,ef,bc);
	}
    else
    {
        DebugStop();
    }
}

void TPZMaterial::Clone(std::map<int, TPZMaterial * >&matvec) {
	int matid = Id();
	std::map<int, TPZMaterial * >::iterator matit;
	matit = matvec.find(matid);
	if(matit != matvec.end()) return;
	TPZMaterial * newmat = NewMaterial();
	matvec[matid] = newmat;
}

/** Get the order of the integration rule necessary to integrate an
 * element with polinomial order p */
int TPZMaterial::IntegrationRuleOrder(int elPMaxOrder) const
{
    int order = 0;
    if(fForcingFunction){
        order = fForcingFunction->PolynomialOrder();
    }
    
    int pmax = elPMaxOrder;
    int integrationorder = 2*pmax;
    if (pmax < order) {
        integrationorder = order+pmax;
    }
    return  integrationorder;
}

int TPZMaterial::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const {
    int order = 0;
    if (fForcingFunction) {
        order = fForcingFunction->PolynomialOrder();
    }

    int pmax = 0;
    for (int ip = 0; ip < elPMaxOrder.size(); ip++) {
        if (elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];
    }
    int integrationorder = 2 * pmax;
    if (pmax < order) {
        integrationorder = order + pmax;
    }

    return integrationorder;
}

int TPZMaterial::ClassId() const{
    return Hash("TPZMaterial");
}

/* Saves the element data to a stream */
void TPZMaterial::Write(TPZStream &buf, int withclassid) const {
    buf.Write(&fId, 1);
    buf.Write(&gBigNumber, 1);
    TPZPersistenceManager::WritePointer(fForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fForcingFunctionExact.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fTimeDependentForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fTimedependentFunctionExact.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fBCForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fTimedependentBCForcingFunction.operator ->(), &buf);
    buf.Write(fLinearContext);
    buf.Write(&fNumLoadCases);
    buf.Write(&fPostProcIndex);
}

/* Reads the element data from a stream */
void TPZMaterial::Read(TPZStream &buf, void *context) {
    buf.Read(&fId, 1);
    buf.Read(&gBigNumber, 1);
    fForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fForcingFunctionExact = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTimeDependentForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTimedependentFunctionExact = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fBCForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTimedependentBCForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    buf.Read(fLinearContext);
    buf.Read(&fNumLoadCases);
    buf.Read(&fPostProcIndex);
}

