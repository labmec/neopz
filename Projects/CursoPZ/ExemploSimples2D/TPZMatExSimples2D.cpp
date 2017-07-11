#include "TPZMatExSimples2D.h"
#include "pzbndcond.h"

TPZMatExSimples2D::TPZMatExSimples2D(int id) : TPZMaterial(id)
{
    
}

/** @brief Default constructor */
TPZMatExSimples2D::TPZMatExSimples2D() : TPZMaterial()
{
    
}


TPZMatExSimples2D::TPZMatExSimples2D(const TPZMatExSimples2D &mat) : TPZMaterial(mat)
{
    
}

TPZMatExSimples2D::~TPZMatExSimples2D()
{
    
}

void TPZMatExSimples2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    int nshape = phi.Rows();
    
    for(int i = 0 ; i<nshape ; i++)
    {
        const STATE rhs = 0;
        ef(i,0) += rhs*weight;
        for(int j=0;j<nshape;j++)
        {
            const STATE stiff = dphi(0,i)*dphi(0,j)+dphi(1,i)*dphi(1,j);
            ek(i,j) += stiff*weight;
        }
    }
    
    
}

void TPZMatExSimples2D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
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













