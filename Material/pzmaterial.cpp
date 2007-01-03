//METHODS DEFINITION FOR CLASS TPZMaterial

#include "pzmaterial.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzadmchunk.h"

#include "pzlog.h"

#ifdef LOG4CXX
  #ifdef DEBUG
    #define DEBUG2
  #endif
  static LoggerPtr logger(Logger::getLogger("pz.material"));
#endif

using namespace std;
REAL TPZMaterial::gBigNumber = 1.e12;

TPZMaterial::TPZMaterial(int id) {
   fId = id;
   fForcingFunction = 0;
}

TPZMaterial::TPZMaterial(const TPZMaterial &material) {
   fId = material.fId;
   fForcingFunction = material.fForcingFunction;
}

void TPZMaterial::Print(std::ostream & out) {
  out << std::endl << "Material Id = " << fId << std::endl;
}

int TPZMaterial::VariableIndex(char *name) {
   if(!strcmp(name,"state")) return 0;
   if(!strcmp(name,"POrder")) return 99;
   if(!strcmp(name,"Error")) return 100;
   if(!strcmp(name,"TrueError")) return 101;
   if(!strcmp(name,"EffectivityIndex")) return 102;


   return -1;
}

int TPZMaterial::NSolutionVariables(int index) {
   if(index == 0) return NStateVariables();
   if(index == 99) return 1;
   if(index == 100) return 1;
   if(index == 101) return 1;
   if(index == 102) return 1;
   PZError << "TPZMaterial::NSolutionVariables called index = " << index << "\n";
   return 0;
}

void TPZMaterial::Solution(TPZVec<REAL> &Sol,TPZFMatrix &/*DSol*/,TPZFMatrix &/*axes*/,int var,
			   TPZVec<REAL> &Solout){
   if(var == 0) Solout = Sol;
   else if(var == 99 || var == 100 || var == 101 || var == 102) {
      //  	PZError << "TPZMaterial var = "<< var << " the element should treat this case\n";
      Solout[0] = Sol[0]; // = 0.;
   } else Solout.Resize(0);
}

TPZBndCond *TPZMaterial::CreateBC(TPZAutoPointer<TPZMaterial> &reference, int id, int typ, TPZFMatrix &val1, TPZFMatrix &val2) {
   return new TPZBndCond(reference,id,typ,val1,val2);
}

void TPZMaterial::SetData(std::istream &data) {
   PZError << "TPZMaterial::SetData is called.\n";
   data >> fId;
}

TPZMaterial *TPZMaterial::NewMaterial() {
   PZError << "TPZMaterial::NewMaterial is called.\n";
   return 0;
}
void TPZMaterial::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv, TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ef){
   TPZFMatrix ek(ef.Rows(),ef.Rows(),0.);
   Contribute(x,jacinv,sol,dsol,weight,axes,phi,dphi,ek,ef);
}

void TPZMaterial::Clone(std::map<int, TPZAutoPointer<TPZMaterial> >&matvec) {
   int matid = Id();
   std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
   matit = matvec.find(matid);
   if(matit != matvec.end()) return;
   TPZAutoPointer<TPZMaterial> newmat = TPZAutoPointer<TPZMaterial>(NewMaterial());
   newmat->SetForcingFunction(TPZMaterial::fForcingFunction);
   matvec[matid] = newmat;
}

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

#define TPZMATERIALID 300
int TPZMaterial::ClassId() const
{
  return TPZMATERIALID;
}

  /**
  Save the element data to a stream
  */
void TPZMaterial::Write(TPZStream &buf, int withclassid)
{
  TPZSaveable::Write(buf,withclassid);
  buf.Write(&fId,1);
}
  
  /**
  Read the element data from a stream
  */
void TPZMaterial::Read(TPZStream &buf, void *context)
{
  TPZSaveable::Read(buf,context);
  buf.Read(&fId,1);
}
