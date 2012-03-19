//$Id: pzpostprocmat.cpp,v 1.7 2010-11-23 18:57:31 diogo Exp $

#include "pzpostprocmat.h"
//#include "poroelastoplasticid.h"
#include "pzbndcond.h"

#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr postprocLogger(Logger::getLogger("material.pzPostProcMat"));
#endif


TPZPostProcMat::TPZPostProcMat() : /*TPZMaterial*/TPZDiscontinuousGalerkin()
{
	fVars.Resize(0);	
	fDimension = -1;
}

TPZPostProcMat::TPZPostProcMat(int id) : /*TPZMaterial*/TPZDiscontinuousGalerkin(id)
{
	fVars.Resize(0);	
	fDimension = -1;
}

TPZPostProcMat::TPZPostProcMat(const TPZPostProcMat &mat) : /*TPZMaterial*/TPZDiscontinuousGalerkin(mat), fVars(mat.fVars), fDimension(mat.fDimension)
{
}

TPZPostProcMat::~TPZPostProcMat()
{
	PZError << "\nMaterial " << Id() << " killed\n";

}

void TPZPostProcMat::Print(std::ostream &out)
{
	out << this->Name();
	out << "\n Base material Data:\n";
	TPZDiscontinuousGalerkin::Print(out);
	out << "Dimension " << fDimension << std::endl;
	int nVars = fVars.NElements();
	out << "\n Post Process Variables\n";
	for(int i = 0; i < nVars; i++)
	{
		out << fVars[i].fName << " of size " << fVars[i].fNumEq << " and index " << fVars[i].fIndex << std::endl;
	}
}

int TPZPostProcMat::VariableIndex(const std::string &name)
{
	int i, nVars = fVars.NElements();
	
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
	int i, nVars = fVars.NElements();
	
	i = 0;

	while(i < nVars && var != fVars[i].fIndex)i++;
	
	if(i >= nVars)return -1; // variable not found
	
	return fVars[i].fNumEq;
}

int TPZPostProcMat::NStateVariables()
{
	int i, nVars = fVars.NElements(), size = 0;
	for(i = 0; i < nVars; i++)size += fVars[i].fNumEq;
	return size;
}

void TPZPostProcMat::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
	
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << ">>> TPZPostProcMat::Solution() *** called for variable index = " << var;
    LOGPZ_DEBUG(postprocLogger,sout.str().c_str());
  }
#endif
	
	int i, nVars = fVars.NElements(), offset = 0;
	
	i = 0;

	while(i < nVars && var != fVars[i].fIndex)
	{
		offset += fVars[i].fNumEq;
		i++;
	}
	
	if(i >= nVars)return; // variable not found
	
	int numeq = fVars[i].fNumEq;
	
	Solout.Resize(numeq);	

	for(i = 0; i < numeq; i++)Solout[i] = data.sol[0][offset+i];
}

void TPZPostProcMat::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{

//#ifdef LOG4CXX
//  {
//    std::stringstream sout;
//    sout << ">>> TPZPostProcMat::Contribute ***";
//    LOGPZ_DEBUG(postprocLogger,sout.str().c_str());
//  }
//#endif
	
  TPZFMatrix<REAL> &phi = data.phi;
  TPZVec<REAL> &sol = data.sol[0];
  int nstate = NStateVariables();
			
  int nshape = phi.Rows();
	
  TPZFMatrix<REAL> L2(nshape,nshape,0.);
  int i, j, i_var;

  for(i = 0; i < nshape; i++)
	 for(j = 0; j < nshape; j++)
		ek(i,j) += phi(i,0) * phi(j,0);
	
  for(i = 0; i < nstate; i++)
	{
		int eqOffset = i*nshape;
		for(i_var = 0; i_var < nshape; i_var++)
			ef(eqOffset+i_var,0) += phi(i_var,0) * sol[i];
	}
}

void TPZPostProcMat::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef)
{
	PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcMat::Contribute(ef) should never be called\n";
	return;
}

void TPZPostProcMat::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc)
{
	PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcMat::ContributeBC() should never be called\n";
	return;
}

void TPZPostProcMat::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ef, TPZFMatrix<REAL> &ek){
  // do nothing
}

void TPZPostProcMat::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<REAL> &ef, TPZFMatrix<REAL> &ek,TPZBndCond &bc){
  // do nothing
}

int TPZPostProcMat::ClassId() const
{
	return TPZPOSTPROCMAT_ID;
}

std::string TPZPostProcMat::Name()
{
	return "TPZPostProcMat"; 
}

void TPZPostProcMat::Write(TPZStream &buf, int withclassid)
{
	TPZSaveable::Write(buf, withclassid);
	
	TPZMaterial::Write(buf, 0);
	
	int i, nVars = fVars.NElements();
	
	buf. Write(&nVars, 1);
	
	for(i = 0; i < nVars; i++)
	{
		buf. Write(&fVars[i].fIndex, 1);
		buf. Write(&fVars[i].fNumEq, 1);
		int size = strlen(fVars[i].fName.c_str());
		if(size > 255) size = 255;
		buf. Write(fVars[i].fName.c_str(), size);
	}
}

void TPZPostProcMat::Read(TPZStream &buf, void *context)
{
    TPZSaveable::Read(buf, context);
	
	TPZMaterial::Read(buf, context);
	
	int i, nVars;
	
	buf.Read(&nVars, 1);
	
	for(i = 0; i < nVars; i++)
	{
		buf. Read(&fVars[i].fIndex, 1);
		buf. Read(&fVars[i].fNumEq, 1);
		int size;
		buf. Read(&size, 1);
		char name[256];
		buf. Read(name, size);
		fVars[i].fName = name;
	}
}

void TPZPostProcMat::FillDataRequirements(TPZMaterialData &data){
  	
	TPZMaterial::FillDataRequirements(data);
	data.SetAllRequirements(false);
	data.fNeedsSol = true;	
}

void TPZPostProcMat::GetPostProcessVarIndexList(TPZVec<int> & varIndexList)
{
	int i, n = fVars.NElements();
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
	
	int i, n = varIndexNames.NElements(), k = 0;
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

int TPZPostProcMat::Dimension()
{
	return fDimension;
}
