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

TPZMatLaplacianLagrange::TPZMatLaplacianLagrange(int nummat, int dim, bool finermesh) : TPZMatLaplacian(nummat,dim)
{
    fIsFinerMesh = finermesh;
}

TPZMatLaplacianLagrange::TPZMatLaplacianLagrange() : TPZMatLaplacian()
{
}


TPZMatLaplacianLagrange & TPZMatLaplacianLagrange::operator=(const TPZMatLaplacianLagrange &copy)
{
	TPZMatLaplacian::operator = (copy);
    fIsFinerMesh = copy.fIsFinerMesh;
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

void TPZMatLaplacianLagrange::Contribute(TPZVec<TPZMaterialData> &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
	
    if(fIsFinerMesh==false)
    {
        ContributeCoarse(data, weight,ek,ef);
        return;
    }
    TPZFMatrix<REAL>  &phi = data[0].phi;
    TPZFMatrix<REAL> &dphi = data[0].dphix;
    TPZVec<REAL>  &x = data[0].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE XfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        TPZFMatrix<STATE> dres(Dimension(),1);
        fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        XfLoc = res[0];
    }
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * XfLoc*(STATE)phi(in,0);
        for(int jn = 0; jn < phr; jn++)
        {
            for(kd=0; kd<fDim; kd++)
            {
                ek(in,jn) += (STATE)weight*(fK*(STATE)(dphi(kd,in)*dphi(kd,jn)));
            }
        }
    }
    for (int in=0; in<phr; in++)
    {
        ek(phr,in) -= (STATE)phi(in,0)*weight;
        ek(in,phr) -= (STATE)phi(in,0)*weight;
    }
    ek(phr+1,phr) += weight;
    ek(phr,phr+1) += weight;
    //ek(phr,phr) -= 1.*(STATE)weight;
    
#ifdef DEBUG
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
#endif
    
}

void TPZMatLaplacianLagrange::ContributeCoarse(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
	
    TPZFMatrix<REAL>  &phi = datavec[3].phi;
    TPZFMatrix<REAL> &dphi = datavec[3].dphix;
    int phr = phi.Rows();
    
    
    //inner product
    for( int in = 0; in < phr; in++ ) {
        int kd;
        for(int jn = 0; jn < phr; jn++)
        {
            ek(in,jn) += (STATE)phi(in,0)*phi(jn,0)*weight;
            
            for(kd=0; kd<fDim; kd++)
            {
                ek(in,jn) += (STATE)weight*(STATE)(dphi(kd,in)*dphi(kd,jn));
            }
        }
    }
}


void TPZMatLaplacianLagrange::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    if (var == 0) {
        Solout[0] = datavec[0].sol[0][0];//+datavec[2].sol[0][0];
    }
    else
    {
        TPZMaterial::Solution(datavec, var, Solout);
    }
}

void TPZMatLaplacianLagrange::Write(TPZStream &buf, int withclassid){
	TPZMatLaplacian::Write(buf, withclassid);
}

void TPZMatLaplacianLagrange::Read(TPZStream &buf, void *context){
	TPZMatLaplacian::Read(buf, context);
}

template class TPZRestoreClass<TPZMatLaplacianLagrange,TPZMatLaplacianLagrangeID>;
