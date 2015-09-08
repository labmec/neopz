/**
 * @file
 * @brief Contains implementations of the TPZMaterial methods.
 */

#include "TPZVecL2.h"
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



TPZVecL2::TPZVecL2() : TPZMaterial() {
}

TPZVecL2::TPZVecL2(int id) : TPZMaterial(id) {
}

TPZVecL2::~TPZVecL2()
{
}


TPZVecL2::TPZVecL2(const TPZVecL2 &material) : TPZMaterial(material) {
}



void TPZVecL2::Print(std::ostream & out) {
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMaterial::Print(out);
    
}

int TPZVecL2::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"state")) return 0;
	if(!strcmp(name.c_str(),"State")) return 0;
	if(!strcmp(name.c_str(),"Solution")) return 0;
	
    return TPZMaterial::VariableIndex(name );
}

int TPZVecL2::NSolutionVariables(int index) {
#ifdef STATE_COMPLEX
	if(index == 0) return NStateVariables()*2;    
#else
	if(index == 0) return 3;
#endif
    return TPZMaterial::NSolutionVariables(index);
}

void TPZVecL2::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
}

void TPZVecL2::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	DebugStop();
}

void TPZVecL2::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout)
{
		this->Solution(data,dataleftvec,datarightvec, var, Solout);
}

void TPZVecL2::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl *left, TPZCompEl *right)
{
	this->Solution(data,dataleftvec,datarightvec, var, Solout, left, right);
}

void TPZVecL2::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
						   TPZVec<STATE> &Solout){
    if(var == 0) Solout = Sol;
    else
    {
        TPZMaterial::Solution(Sol, DSol,axes,var,Solout);
    }
}


TPZMaterial * TPZVecL2::NewMaterial() {
    return new TPZVecL2(*this);
}

void TPZVecL2::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->Contribute(data, weight, fakeek, ef);
}

void TPZVecL2::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeBC(data, weight, fakeek, ef, bc);
}

void TPZVecL2::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	int nref=datavec.size();
	if (nref== 1) {
		this->Contribute(datavec[0], weight, ek,ef);
	}
}




int TPZVecL2::ClassId() const
{
	return TPZVECL2ID;
}

/* Saves the element data to a stream */
void TPZVecL2::Write(TPZStream &buf, int withclassid)
{
	TPZMaterial::Write(buf,withclassid);
}

/* Reads the element data from a stream */
void TPZVecL2::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
 * @param datavec [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 */
void TPZVecL2::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    TPZManVector<STATE,3> force(3);
    force.Fill(0.);
    if(fForcingFunction) {
        fForcingFunction->Execute(data.x,force);
    }
    
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = data.phi;
    
    int phrq;
    phrq = data.fVecShapeIndex.NElements();
    
    //Calculate the matrix contribution for flux. Matrix A
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = data.fVecShapeIndex[iq].first;
        int ishapeind = data.fVecShapeIndex[iq].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = data.fNormalVec(id,ivecind);
        }
        STATE ff = 0.;
        for (int i=0; i<3; i++) {
            ff += ivec(i,0)*force[i];
        }
        
        ef(iq,0) += weight*ff*phiQ(ishapeind,0);
        
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = data.fVecShapeIndex[jq].first;
            int jshapeind = data.fVecShapeIndex[jq].second;
            
            for(int id=0; id<3; id++){
                jvec(id,0) = data.fNormalVec(id,jvecind);
            }
            
            //jvecZ.Print("mat1 = ");
            REAL prod1 = ivec(0,0)*jvec(0,0) + ivec(1,0)*jvec(1,0) + ivec(2,0)*jvec(2,0);
            ek(iq,jq) += weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
           
        }
    }
}

void TPZVecL2::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){

    
    TPZFMatrix<REAL>  &phiQ = data.phi;
    int phrq = phiQ.Rows();
    
    //AQUIFRANREAL v2;
  STATE v2;
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> res(3);
        bc.ForcingFunction()->Execute(data.x,res);
        v2 = res[0];
    }else
    {
        v2 = bc.Val2()(0,0);
    }
    
    switch (bc.Type()) {
        case 0 :		// Dirichlet condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2*phiQ(iq,0)*weight;
            }
            break;
            
        case 1 :			// Neumann condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                ef(iq,0)+= gBigNumber*v2*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                }
            }
            break;
            
        case 2 :			// mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
                ef(iq,0) += v2*phiQ(iq,0)*weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq,jq) += weight*bc.Val1()(0,0)*phiQ(iq,0)*phiQ(jq,0);
                }
            }
            
            break;
    }

}

