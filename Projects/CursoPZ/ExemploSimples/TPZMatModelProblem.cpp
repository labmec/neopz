#include "TPZMatModelProblem.h"
#include "pzbndcond.h"

TPZMatModelProblem::TPZMatModelProblem(int id) : TPZMaterial(id)
{
    
}

/** @brief Default constructor */
TPZMatModelProblem::TPZMatModelProblem() : TPZMaterial()
{
    
}


TPZMatModelProblem::TPZMatModelProblem(const TPZMatModelProblem &mat) : TPZMaterial(mat)
{
    
}

TPZMatModelProblem::~TPZMatModelProblem()
{
    
}

void TPZMatModelProblem::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL> &x = data.x;
    int nshape=phi.Rows();
    
    for(int i = 0 ; i<nshape ; i++)
    {
        const STATE rhs = x[0]*phi(i,0);
        ef(i,0) += rhs*weight;
        for(int j=0;j<nshape;j++)
        {
            const STATE stiff = dphi(0,i)*dphi(0,j)+phi(i,0)*phi(j,0);
            ek(i,j) += stiff*weight;
        }
    }
    
    
}

void TPZMatModelProblem::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    TPZFMatrix<REAL> &phi = data.phi;
    int nshape=phi.Rows();
    
    REAL BIG = TPZMaterial::gBigNumber;//sera posto na matriz K
    STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
    STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
    
    switch (bc.Type()) {
        case 0: // Dirichlet
            for(int i = 0 ; i<nshape ; i++)
            {
                const STATE rhs = phi(i,0)*BIG*v2;
                ef(i,0) += rhs*weight;
                for(int j=0;j<nshape;j++)
                {
                    const STATE stiff = phi(i,0)*phi(j,0);
                    ek(i,j) += stiff*weight*BIG;
                }
            }
            break;
            
        case 1: // Neumann
            for(int i = 0 ; i<nshape ; i++)
            {
                const STATE rhs = phi(i,0)*v2;
                ef(i,0) += rhs*weight;
            }
            break;
            
        case 2: // Mista
            for(int i = 0 ; i<nshape ; i++)
            {
                const STATE rhs = phi(i,0)*v2;
                ef(i,0) += rhs*weight;
                for(int j=0;j<nshape;j++)
                {
                    const STATE stiff = v1*phi(i,0)*phi(j,0);
                    ek(i,j) += stiff*weight;
                }
            }
            break;
            
        default:
            DebugStop();
            break;
    }
    
}

int TPZMatModelProblem::ClassId() const{
    return Hash("TPZMatModelProblem") ^ TPZMaterial::ClassId() << 1;
}

template class
TPZRestoreClass<TPZMatModelProblem>;
