/**
 * @file
 */

#include "pzpostprocmat.h"
//#include "poroelastoplasticid.h"
#include "pzbndcond.h"

#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr postprocLogger(Logger::getLogger("material.pzPostProcMat"));
#endif


TPZPostProcMat::TPZPostProcMat() :TPZRegisterClassId(&TPZPostProcMat::ClassId),TPZDiscontinuousGalerkin()
{
	fVars.Resize(0);	
	fDimension = -1;
}

TPZPostProcMat::TPZPostProcMat(int64_t id) : TPZRegisterClassId(&TPZPostProcMat::ClassId),TPZDiscontinuousGalerkin(id)
{
	fVars.Resize(0);	
	fDimension = -1;
}

TPZPostProcMat::TPZPostProcMat(const TPZPostProcMat &mat) : TPZRegisterClassId(&TPZPostProcMat::ClassId),TPZDiscontinuousGalerkin(mat), fVars(mat.fVars), fDimension(mat.fDimension)
{
}

TPZPostProcMat::~TPZPostProcMat()
{
#ifdef PZDEBUG
    std::cout << "TPZPostProcMat:: Material Id = " << Id() << ", it is being deleted. " << std::endl;
#endif

}

void TPZPostProcMat::Print(std::ostream &out)
{
	out << this->Name();
	out << "\n Base material Data:\n";
	TPZDiscontinuousGalerkin::Print(out);
	out << "Dimension " << fDimension << std::endl;
	int64_t nVars = fVars.NElements();
	out << "\n Post Process Variables\n";
	for(int64_t i = 0; i < nVars; i++)
	{
		out << fVars[i].fName << " of size " << fVars[i].fNumEq << " and index " << fVars[i].fIndex << std::endl;
	}
}

int TPZPostProcMat::VariableIndex(const std::string &name)
{
	int64_t i, nVars = fVars.NElements();
	
	i = 0;

	while(i < nVars && strcmp(fVars[i].fName.c_str(), name.c_str()))i++;
	
	if(i >= nVars)
	{
		PZError << "TPZPostProcMat::Variable " << name << " not found\n";
		return -1; // variable not found
	}
	
	return fVars[i].fIndex;
}

int TPZPostProcMat::NSolutionVariables(int var)
{
	int64_t i, nVars = fVars.NElements();
	
	i = 0;

	while(i < nVars && var != fVars[i].fIndex)i++;
	
	if(i >= nVars)return -1; // variable not found
	
	return fVars[i].fNumEq;
}

int TPZPostProcMat::NStateVariables()
{
	int64_t i, nVars = fVars.NElements(), size = 0;
	for(i = 0; i < nVars; i++)size += fVars[i].fNumEq;
	return size;
}

void TPZPostProcMat::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
	
#ifdef LOG4CXX_keep
  {
    std::stringstream sout;
    sout << ">>> TPZPostProcMat::Solution() *** called for variable index = " << var;
    LOGPZ_DEBUG(postprocLogger,sout.str().c_str());
  }
#endif
	
	int64_t i, nVars = fVars.NElements(), offset = 0;
	
	i = 0;

	while(i < nVars && var != fVars[i].fIndex)
	{
		offset += fVars[i].fNumEq;
		i++;
	}
	
	if(i >= nVars)return; // variable not found
	
	int64_t numeq = fVars[i].fNumEq;
	
	Solout.Resize(numeq);	

	for(i = 0; i < numeq; i++)Solout[i] = data.sol[0][offset+i];
}

void TPZPostProcMat::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	
  TPZFMatrix<REAL> &phi = data.phi;
  TPZVec<STATE> &sol = data.sol[0];
  int nstate = NStateVariables();
			
  int64_t nshape = phi.Rows();
  int64_t i, j, i_var;

  for(i = 0; i < nshape; i++)
	 for(j = 0; j < nshape; j++)
		ek(i,j) += phi(i,0) * phi(j,0);
	
  for(i = 0; i < nstate; i++)
	{
		int64_t eqOffset = i*nshape;
		for(i_var = 0; i_var < nshape; i_var++)
			ef(eqOffset+i_var,0) += (STATE)phi(i_var,0) * sol[i];
	}
}

void TPZPostProcMat::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcMat::Contribute(ef) should never be called\n";
    return;
}

void TPZPostProcMat::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcMat::ContributeBC() should never be called\n";
	return;
}

void TPZPostProcMat::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef, TPZFMatrix<STATE> &ek){
  // do nothing
}

void TPZPostProcMat::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef, TPZFMatrix<STATE> &ek,TPZBndCond &bc){
  // do nothing
}

int TPZPostProcMat::ClassId() const{
    return Hash("TPZPostProcMat") ^ TPZDiscontinuousGalerkin::ClassId() << 1;
}

std::string TPZPostProcMat::Name() {
    return "TPZPostProcMat";
}

void TPZPostProcMat::Write(TPZStream &buf, int withclassid) const
{
	
	TPZDiscontinuousGalerkin::Write(buf, withclassid);
	buf.Write(fVars);
	buf.Write(&fDimension);
}

void TPZPostProcMat::Read(TPZStream &buf, void *context)
{
    	TPZDiscontinuousGalerkin::Read(buf, context);
	
	buf.Read(fVars);
	buf.Read(&fDimension);
}

void TPZPostProcMat::FillDataRequirements(TPZMaterialData &data){
  	
	TPZMaterial::FillDataRequirements(data);
	data.SetAllRequirements(false);
	data.fNeedsSol = true;	
}

void TPZPostProcMat::GetPostProcessVarIndexList(TPZVec<int> & varIndexList)
{
	int64_t i, n = fVars.NElements();
	varIndexList.Resize(n);
	
	for(i = 0; i < n; i++)varIndexList[i] = fVars[i].fIndex;
}

void TPZPostProcMat::SetPostProcessVarIndexList(TPZVec<std::string> & varIndexNames, TPZMaterial * pRefMat)
{
	if(!pRefMat)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcMat::SetPostProcessVarIndexList() without valid reference material to post process\n";
		return;
	}
	
	int64_t i, n = varIndexNames.NElements(), k = 0;
	int varindex;
	fVars.Resize(n);
	
	for(i = 0; i < n; i++)
	{
		varindex = pRefMat->VariableIndex(varIndexNames[i]);
		if(varindex >= 0)
		{
			fVars[k].fIndex = varindex;
			fVars[k].fNumEq = pRefMat->NSolutionVariables(varindex);
			fVars[k].fName  = varIndexNames[i];
            
			k++;
		}
	}
    
	
	fVars.Resize(k);
	fDimension = pRefMat->Dimension();
	
}

int TPZPostProcMat::Dimension() const
{
	return fDimension;
}

int TPZPostProcVar::ClassId() const {
    return Hash("TPZPostProcVar");
}

void TPZPostProcVar::Read(TPZStream& buf, void* context) {
    buf.Read(&fIndex);
    buf.Read(&fName);
    buf.Read(&fNumEq);    
}

void TPZPostProcVar::Write(TPZStream& buf, int withclassid) const {
    buf.Write(&fIndex);
    buf.Write(&fName);
    buf.Write(&fNumEq);    
}

template class TPZRestoreClass<TPZPostProcVar>;
template class TPZRestoreClass<TPZPostProcMat>;
