#include "Elasticity/TPZHybridElasticity3D.h"
#include "pzaxestools.h"
#include "pzvec_extras.h"

void TPZHybridElasticity3D::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    if(data.size() > 4)
    {
        std::cout << "Please implement me!\n";
        DebugStop();
    }
    TPZElasticity3D::Contribute(data[1],weight,ek,ef);
    if(data.size() > 2) {
        int nphi = (int) data[1].phi.Rows();
        int64_t nek = ek.Rows();
        TPZManVector<REAL,3> delx(3,0.);
        delx = data[2].x - data[2].XCenter;
        TPZFNMatrix<18,REAL> rb(3, 6,0.);
        TPZFNMatrix<120,REAL> phivec(3, nek,0.);
        for(int i=0; i<3; i++) rb(i,i) = 1.;
        rb(1,3) = -delx[2]; rb(2,3) = delx[1];
        rb(0,4) = delx[2]; rb(1,4) = -delx[0];
        rb(0,5) = -delx[1]; rb(1,5) = delx[0];
        for(int ip=0; ip<nphi; ip++) {
            for(int i=0; i<3; i++) {
                phivec(i,3*ip+i) = data[1].phi(ip,0);
            }
        }
        for(int i=0; i<3; i++) {
            for(int j=0; j<6; j++) {
                phivec(i,3*nphi+6+j) = -rb(i,j);
            }
        }
        ek.AddContribution(3*nphi, 0, rb, 1, phivec, 0, weight);
        ek.AddContribution(0, 3*nphi, phivec, 1, rb, 0, weight);
    }
}

void TPZHybridElasticity3D::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                                   TPZBndCondT<STATE> &bc) {
    
    TPZFMatrix<REAL>  &phi_u = datavec[1].phi;
    TPZFMatrix<REAL>  &phi_flux = datavec[0].phi;
    //    TPZFMatrix<REAL> &axes = data.axes;
    int64_t phr_primal = phi_u.Rows();
    int64_t phr_hybrid = phi_flux.Rows();
    bool primal = true;   /// weather pressure or flux is hybridized
    TPZManVector<REAL,3> x(3), normal;
    if(phr_hybrid)
    {
        primal = false;
        x = datavec[0].x;
        normal = datavec[0].normal;
    }
    else
    {
        x = datavec[1].x;
        normal = datavec[1].normal;
    }
    short in,jn;
    const int dim = 3;
    TPZManVector<REAL,3> v2 = bc.Val2();

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE, 4> u(dim);
        TPZFNMatrix<4, STATE> gradu(dim, dim), deform_voight(6,1), sigma_voight(6,1), sigma(3,3);
        bc.ForcingFunctionBC()(x, u, gradu);
        ComputeStressTensor(sigma, gradu);
        if(bc.Type() == 0) {
            v2 = u;
        } else if(bc.Type() == 1) {
            v2[0] = 0.;
            v2[1] = 0.;
            v2[2] = 0.;
            for (int i = 0; i < dim; i++) {
                v2[0] += sigma(0,i) * normal[i];
                v2[1] += sigma(1,i) * normal[i];
                v2[2] += sigma(2,i) * normal[i];
            }
        } else if(bc.Type() == 2) {
           DebugStop(); // I have to think how to implement this
        }
    }

    if(primal)
    {
        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr_primal; in++) {
                    for(int idf = 0; idf<2; idf++) {
                        ef(2*in+idf,0) += (STATE)(fBigNumber* phi_u(in,0) * weight) * v2[idf];
                        for (jn = 0 ; jn < phr_primal; jn++) {
                            ek(2*in+idf,2*jn+idf) += fBigNumber * phi_u(in,0) * phi_u(jn,0) * weight;
                        }
                    }
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phr_primal; in++) {
                    for(int idf = 0; idf<2; idf++)
                        ef(2*in+idf,0) += v2[idf] * (STATE)(phi_u(in,0) * weight);
                }
                break;
            default:
                DebugStop();
        }
    } else
    {
        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr_hybrid; in++) {
                    for(int idf = 0; idf<2; idf++) {
                        ef(2*in+idf,0) += v2[idf] * (STATE)(phi_flux(in,0) * weight);
                    }
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phr_hybrid; in++) {
                    for(int idf = 0; idf<2; idf++) {
                        ef(2*in+idf,0) -= (STATE)(fBigNumber* phi_flux(in,0) * weight) * v2[idf];
                        for (jn = 0 ; jn < phr_hybrid; jn++) {
                            ek(2*in+idf,2*jn+idf) += fBigNumber * phi_flux(in,0) * phi_flux(jn,0) * weight;
                        }
                    }
                }
                break;
        }
    }
}

void TPZHybridElasticity3D::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity3D::FillDataRequirements(data[1]);
}

void TPZHybridElasticity3D::FillBoundaryConditionDataRequirements(int type,
                                                           TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int64_t nsp = datavec.size();
    for(int i=0; i<nsp; i++) {
        auto &data = datavec[i];
        data.SetAllRequirements(false);
        if (type == 3 || type == 1) {
            data.fNeedsNormal = true;
        }
        if (HasForcingFunction()) {
            data.fNeedsNormal = true;
        }
    }
}




void TPZHybridElasticity3D::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data,
                               int var, TPZVec<STATE> &Solout)
{
    TPZElasticity3D::Solution(data[1], var, Solout);
}


TPZMaterial* TPZHybridElasticity3D::NewMaterial() const
{
    return new TPZHybridElasticity3D(*this);
}


int TPZHybridElasticity3D::ClassId() const
{
    return Hash("TPZHybridElasticity3D") ^
        TBase::ClassId() << 1;
}

void TPZHybridElasticity3D::Read(TPZStream &buf, void *context)
{
    TPZElasticity3D::Read(buf,context);
}
    
void TPZHybridElasticity3D::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity3D::Write(buf,withclassid);
}


void TPZHybridElasticity3D::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                             TPZVec<REAL> &values) {
    TPZElasticity3D::Errors(data[1],values);
}

