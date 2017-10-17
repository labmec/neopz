/**
 * @file
 * @brief Contains implementations of the TPZMatLaplacian methods.
 */

#include "TPZMatLaplacianLagrange.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"
#include "pzaxestools.h"

#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.TPZMatLaplacian"));
#endif


using namespace std;

TPZMatLaplacianLagrange::TPZMatLaplacianLagrange(int nummat, int dim) : TPZMatLaplacian(nummat,dim)
{
    fIsDPGPhil = false;
    fUseMDP = false;
}

TPZMatLaplacianLagrange::TPZMatLaplacianLagrange() : TPZMatLaplacian()
{
    fIsDPGPhil = false;
    fUseMDP = false;
}


TPZMatLaplacianLagrange & TPZMatLaplacianLagrange::operator=(const TPZMatLaplacianLagrange &copy)
{
	TPZMatLaplacian::operator = (copy);
    fIsDPGPhil = copy.fIsDPGPhil;
    fUseMDP = copy.fUseMDP;
	return *this;
}

TPZMatLaplacianLagrange::TPZMatLaplacianLagrange(const TPZMatLaplacianLagrange &copy) : TPZMatLaplacian(copy)
{
    this->operator=(copy);
}

TPZMatLaplacianLagrange::~TPZMatLaplacianLagrange() {
}


void TPZMatLaplacianLagrange::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
    TPZMatLaplacian::Print(out);
}

void TPZMatLaplacianLagrange::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
	if (fIsDPGPhil){
        ContributeDPGPhil(datavec, weight, ek, ef);
        return;
    }
    
    if (fUseMDP){
        ContributeMDP(datavec, weight, ek, ef);
        return;
    }
    
    //contribuicoes referente a malha fina: 
    TPZFMatrix<REAL>  &phif = datavec[0].phi;
    TPZFMatrix<REAL> &dphif = datavec[0].dphix;
    TPZVec<REAL>  &x = datavec[0].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr_f = phif.Rows();
    
    STATE XfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
        XfLoc = res[0];
    }
    
    //Produto interno (e,y) = gard(e).grad(y) + e*y. Obter espaco teste otimo
    for(int in = 0; in < phr_f; in++)
    {
        int kd;
        ef(in, 0) += (STATE)weight*XfLoc*(STATE)phif(in,0);
        for(int jn = 0; jn < phr_f; jn++)
        {
            ek(in,jn) += (STATE)weight*((STATE)(phif(in,0)*phif(jn,0)));
            
            for(kd=0; kd<fDim; kd++)
            {
                ek(in,jn) += (STATE)weight*((STATE)(dphif(kd,in)*dphif(kd,jn)));
            }
        }
    }
    
    
    for (int in=0; in<phr_f; in++)
    {
        ek(phr_f,in) -= (STATE)phif(in,0)*weight;
        ek(in,phr_f) -= (STATE)phif(in,0)*weight;
    }
    
    //dessingularizando a matriz
    ek(phr_f+1,phr_f) += weight;
    ek(phr_f,phr_f+1) += weight;
    //ek(phr,phr) -= 1.*(STATE)weight;


    //Equacao de Poisson: espaco coarce
    TPZFMatrix<REAL>  &phic = datavec[3].phi;
    TPZFMatrix<REAL> &dphic = datavec[3].dphix;
    int phr_c = phic.Rows();

    for( int in = 0; in < phr_f; in++ ) {
        int kd;
        for(int jn = 0; jn < phr_c; jn++)
        {
            for(kd=0; kd<fDim; kd++)
            {
                //B
                ek(in, phr_f + 2 + jn) += (STATE)weight*(fK*(STATE)(dphif(kd,in)*dphic(kd,jn)));
                //B^T
                ek(phr_f + 2 + jn, in) += (STATE)weight*(fK*(STATE)(dphic(kd,jn)*dphif(kd,in)));
            }
        }
    }
    
#ifdef PZDEBUG
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
#endif
    
}

void TPZMatLaplacianLagrange::ContributeDPGPhil(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    //contribuicoes referente a malha fina:
    TPZFMatrix<REAL>  &phif = datavec[0].phi;
    TPZFMatrix<REAL> &dphif = datavec[0].dphix;
    TPZVec<REAL>  &x = datavec[0].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr_f = phif.Rows();
    
    STATE XfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
        XfLoc = res[0];
    }
    
    //Obter espaco teste otimo: obter e in V
    //Operador original b(e,y) = gard(e).grad(y), y in V
    for(int in = 0; in < phr_f; in++)
    {
        int kd;
        ef(in, 0) += (STATE)weight*XfLoc*(STATE)phif(in,0);
        for(int jn = 0; jn < phr_f; jn++)
        {
            for(kd=0; kd<fDim; kd++)
            {
                ek(in,jn) += (STATE)weight*(fK*(STATE)(dphif(kd,in)*dphif(kd,jn)));
            }
        }
    }
    
    
    for (int in=0; in<phr_f; in++)
    {
        ek(phr_f,in) -= (STATE)phif(in,0)*weight;
        ek(in,phr_f) -= (STATE)phif(in,0)*weight;
    }
    
    //dessingularizando a matriz
    ek(phr_f+1,phr_f) += weight;
    ek(phr_f,phr_f+1) += weight;
    //ek(phr,phr) -= 1.*(STATE)weight;
    
    //Espaco coarce: obter (u,q) in Xh
    TPZFMatrix<REAL>  &phic = datavec[3].phi;
    TPZFMatrix<REAL> &dphic = datavec[3].dphix;
    int phr_c = phic.Rows();
    
    //matrix B: ((u,q),y) = gard(u).grad(y) + u*y + q[y], y in V.
    for(int in = 0; in < phr_f; in++) {
        int kd;
        for(int jn = 0; jn < phr_c; jn++)
        {
            //B
            ek(in, phr_f + 2 + jn) += (STATE)weight*((STATE)(phif(in,0)*phic(jn,0)));
            
            //B^T
            //ek(phr_f + 2 + jn, in) += (STATE)weight*((STATE)(phif(in,0)*phic(jn,0)));
            
            for(kd=0; kd<fDim; kd++)
            {
                //B
                ek(in, phr_f + 2 + jn) += (STATE)weight*((STATE)(dphif(kd,in)*dphic(kd,jn)));
                //B^T
                //ek(phr_f + 2 + jn, in) += (STATE)weight*((STATE)(dphic(kd,jn)*dphif(kd,in)));
            }
        }
    }
    
    //B^T: ((z,r),e) = gard(z).grad(e) + z*e + r[e], for all (z,r) in Xh
    for(int in = 0; in < phr_c; in++) {
        int kd;
        for(int jn = 0; jn < phr_f; jn++)
        {
            ek(phr_f + 2 + in, jn) += (STATE)weight*((STATE)(phic(in,0)*phif(jn,0)));
            
            for(kd=0; kd<fDim; kd++)
            {
                ek(phr_f + 2 + in, jn) += (STATE)weight*((STATE)(dphic(kd,in)*dphif(kd,jn)));
            }
        }
    }

#ifdef PZDEBUG
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
#endif

}

void TPZMatLaplacianLagrange::ContributeMDP(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    //contribuicoes referente a malha fina:
    TPZFMatrix<REAL>  &phif = datavec[0].phi;
    TPZFMatrix<REAL> &dphif = datavec[0].dphix;
    TPZVec<REAL>  &x = datavec[0].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr_f = phif.Rows();
    
    STATE XfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
        XfLoc = res[0];
    }
    
    //Obter espaco teste otimo: obter e in V
    //Operador original a(ur,vr) = gard(ur).grad(vr), vr in Vr
    for(int in = 0; in < phr_f; in++)
    {
        int kd;
        ef(in, 0) += (STATE)weight*XfLoc*(STATE)phif(in,0);
        for(int jn = 0; jn < phr_f; jn++)
        {
            for(kd=0; kd<fDim; kd++)
            {
                ek(in,jn) += (STATE)weight*(fK*(STATE)(dphif(kd,in)*dphif(kd,jn)));
            }
        }
    }
    
    
    for (int in=0; in<phr_f; in++)
    {
        ek(phr_f,in) -= (STATE)phif(in,0)*weight;
        ek(in,phr_f) -= (STATE)phif(in,0)*weight;
    }
    
    //dessingularizando a matriz
    ek(phr_f+1,phr_f) += weight;
    ek(phr_f,phr_f+1) += weight;
    //ek(phr,phr) -= 1.*(STATE)weight;
    
    //Produto interno: (u,v) = grad(u)grad(v) + uv
    //Espaco coarce: obter u in Xh
    //(ur,v) - (u,v) = 0, v in Xh.
    TPZFMatrix<REAL>  &phic = datavec[3].phi;
    TPZFMatrix<REAL> &dphic = datavec[3].dphix;
    int phr_c = phic.Rows();
    
    //(ur,v)
    for( int in = 0; in < phr_c; in++ ) {
        int kd;
        for(int jn = 0; jn < phr_f; jn++)
        {
            //(ur,v)
            ek(phr_f + 2 + in, jn) += (STATE)weight*((STATE)(phic(in,0)*phif(jn,0)));
            
            for(kd=0; kd<fDim; kd++)
            {
                //(grad(ur),grad(v))
                ek(phr_f + 2 + in, jn) += (STATE)weight*((STATE)(dphic(kd,in)*dphif(kd,jn)));
            }
        }
    }
    
    //-(u,v)
    for( int in = 0; in < phr_c; in++ ) {
        int kd;
         //ef(phr_f + 2 +in, 0) += (STATE)weight*XfLoc*(STATE)phic(in,0);
        for(int jn = 0; jn < phr_c; jn++)
        {
            //-(u,v)
            ek(phr_f + 2 + in, phr_f + 2 + jn) += (-1.)*(STATE)weight*((STATE)(phic(in,0)*phic(jn,0)));
            
            for(kd=0; kd<fDim; kd++)
            {
                //-(grad(u),grad(v))
                ek(phr_f + 2 + in, phr_f + 2 + jn) += (-1.)*(STATE)weight*((STATE)(dphic(kd,in)*dphic(kd,jn)));
            }
        }
    }

    
#ifdef PZDEBUG
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
#endif
}


void TPZMatLaplacianLagrange::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	
	TPZFMatrix<REAL>  &phif = datavec[0].phi;
    TPZFMatrix<REAL>  &phic = datavec[3].phi;
    //	TPZFMatrix<REAL> &axes = data.axes;
    
	int phrf = phif.Rows();
    int phrc = phic.Rows();
	short in,jn;
	STATE v2[1];
	v2[0] = bc.Val2()(0,0);

	if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
		TPZManVector<STATE> res(1);
		bc.ForcingFunction()->Execute(datavec[3].x,res);       // dphi(i,j) = dphi_j/dxi
		v2[0] = res[0];
	}
    
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
			for(in = 0 ; in < phrc; in++) {
				ef(phrf+in,0) += (STATE)(gBigNumber* phic(in,0)*weight)*v2[0];
				for (jn = 0 ; jn < phrc; jn++) {
					ek(phrf+in,phrf+jn) += gBigNumber * phic(in,0)*phic(jn,0)*weight;
				}
			}
			break;
            
        case 10 :			// NeumannDirichlet: Neumann fine mesh and Dirochlet coarse mesh
            //Neumann
            for(in = 0 ; in < phif.Rows(); in++) {
				ef(in,0) += v2[0]*(STATE)(phif(in,0)*weight);
			}
            
            //Dirichlet
			for(in = 0 ; in < phrc; in++) {
				ef(phrf+in,0) += (STATE)(gBigNumber* phic(in,0)*weight)*v2[0];
				for (jn = 0 ; jn < phrc; jn++) {
					ek(phrf+in,phrf+jn) += gBigNumber * phic(in,0)*phic(jn,0)*weight;
				}
			}
			break;
            
            
		case 1 :			// Neumann condition
			for(in = 0 ; in < phif.Rows(); in++) {
                //std::cout<<"nao implementado\n";
                //DebugStop();
				ef(in,0) += v2[0]*(STATE)(phif(in,0)*weight);
			}
			break;
		case 2 :		// mixed condition
			for(in = 0 ; in < phif.Rows(); in++) {
                std::cout<<"nao implementado\n";
                DebugStop();
				//ef(in, 0) += v2[0]*(STATE)(phif(in, 0)*weight);
				//for (jn = 0 ; jn < phif.Rows(); jn++) {
					//ek(in,jn) += bc.Val1()(0,0)*(STATE)(phif(in,0)*phif(jn,0)*weight);     // peso de contorno => integral de contorno
				//}
			}
			break;
	}

	if (this->IsSymetric()) {//only 1.e-3 because of bignumbers.
		if ( !ek.VerifySymmetry( 1.e-3 ) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
}


int TPZMatLaplacianLagrange::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str()))        return  1;
	if(!strcmp("Derivative",name.c_str()))      return  2;
	if(!strcmp("ExactSolution",name.c_str()))   return  3;
    if(!strcmp("ErrorEstimatorDPG",name.c_str())) return 4;
	return TPZMaterial::VariableIndex(name);
}

int TPZMatLaplacianLagrange::NSolutionVariables(int var){
	if(var==1 || var==3 || var==4) return 1;
	if(var == 2) return fDim;
    
	return TPZMaterial::NSolutionVariables(var);
}



void TPZMatLaplacianLagrange::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    if (var == 1) {
        Solout[0] = datavec[3].sol[0][0];
        return;
    }
    
    if (var == 4) {
        Solout[0] = datavec[0].sol[0][0];
        return;
    }

    
    TPZVec<STATE> solExata(1);
	TPZFMatrix<STATE> flux(2,1);
    //Exact soluion
	if(var == 3){
		fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}//var6

}

void TPZMatLaplacianLagrange::Write(TPZStream &buf, int withclassid) const{
	TPZMatLaplacian::Write(buf, withclassid);
}

void TPZMatLaplacianLagrange::Read(TPZStream &buf, void *context){
	TPZMatLaplacian::Read(buf, context);
}

template class TPZRestoreClass<TPZMatLaplacianLagrange>;
