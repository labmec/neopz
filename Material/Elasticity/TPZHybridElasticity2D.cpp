//
//  TPZHybridElasticity2D.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 22/10/25.
//

#include "TPZHybridElasticity2D.h"

#include "pzaxestools.h"
#include <array>
#include <iostream>
#include <iomanip>


// Plane stress D matrix
static void D_plane_stress(double E, double nu, TPZFMatrix<STATE> &D) {
    const double denom = 1.0 - nu*nu;      // 1 - nu^2
    const double C = E / denom;
    D(0,0) = C * 1.0;
    D(0,1) = C * nu;
    D(1,0) = C * nu;
    D(1,1) = C * 1.0;
    D(2,2) = C * (1.0 - nu) / 2.0;
    D(0,2) = D(1,2) = D(2,0) = D(2,1) = 0.0;
}

// Plane strain D matrix
static void D_plane_strain(double E, double nu, TPZFMatrix<STATE> &D) {
    const double denom = (1.0 + nu) * (1.0 - 2.0*nu); // (1+nu)(1-2nu)
    const double C = E / denom;
    D(0,0) = C * (1.0 - nu);
    D(0,1) = C * nu;
    D(1,0) = C * nu;
    D(1,1) = C * (1.0 - nu);
    D(2,2) = C * (1.0 - 2.0*nu) / 2.0;
    D(0,2) = D(1,2) = D(2,0) = D(2,1) = 0.0;
}

TPZHybridElasticity2D::TPZHybridElasticity2D(int id, STATE E, STATE nu, STATE fx, STATE fy, int planestress) :
                               TPZMatCombinedSpacesT<STATE>(), TPZElasticity2D(id,E,nu,fx,fy,planestress) {}

TPZHybridElasticity2D::TPZHybridElasticity2D(int id) :
                                        TPZElasticity2D(id) {}





int TPZHybridElasticity2D::ClassId() const {
    return Hash("TPZHybridElasticity2D") ^ TPZElasticity2D::ClassId() << 1;
}

TPZMaterial *TPZHybridElasticity2D::NewMaterial() const {
    return new TPZHybridElasticity2D(*this);
}

void TPZHybridElasticity2D::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    TPZElasticity2D::Print(out);
}

/** @name Contribute */
/** @{ */
/**
 * @brief It computes a contribution to the stiffness matrix
 * and load vector at one integration point.
 * @param[in] datavec stores all input data
 * @param[in] weight is the weight of the integration rule
 * @param[out] ek is the element matrix
 * @param[out] ef is the rhs vector
 */
void TPZHybridElasticity2D::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                        REAL weight,TPZFMatrix<STATE> &ek,
                        TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    int64_t phr = dphi.Cols();
    TPZFMatrix<REAL> dphixyz(Dimension(),phr);
    TPZAxesTools<REAL>::Axes2XYZ(dphi,dphixyz,datavec[1].axes);
    
    
    TPZFNMatrix<30,REAL> BMat(3, phr*2, 0.), DBMat(3, phr, 0.);
    TPZFNMatrix<9,REAL> D(3, 3, 0.);
    ComputeD(datavec[0].x, D);
    for (int64_t i = 0; i<phr; i++) {
        BMat(0,2*i) = dphixyz(0,i);
        BMat(2,2*i) = dphixyz(1,i);
        BMat(1,2*i+1) = dphixyz(1,i);
        BMat(2,2*i+1) = dphixyz(0,i);
    }
    D.Multiply(BMat, DBMat);
    ek.AddContribution(0, 0, BMat, 1, DBMat, 0, weight);

    if(ek.Rows() == 2*phr) return;
    if(ek.Rows() != 2*phr+6) DebugStop();

    TPZFMatrix<REAL> &phi = datavec[1].phi;
    TPZVec<REAL> &xcenter = datavec[2].XCenter;
    TPZVec<REAL> &x = datavec[0].x;
    if(datavec.size() >=4){
        for (int in =0; in < phr; in++) {
            ek(2*phr,2*in) += weight*phi(in,0);//lambda*phi
            ek(2*phr+1,2*in+1) += weight*phi(in,0);
            ek(2*phr+2,2*in) += -weight*(x[1]-xcenter[1])*phi(in,0);
            ek(2*phr+2,2*in+1) += weight*(x[0]-xcenter[0])*phi(in,0);
            ek(2*in,2*phr) += weight*phi(in,0);//lambda*phi
            ek(2*in+1,2*phr+1) += weight*phi(in,0);
            ek(2*in,2*phr+2) += -weight*(x[1]-xcenter[1])*phi(in,0);
            ek(2*in+1,2*phr+2) += weight*(x[0]-xcenter[0])*phi(in,0);
        }
        TPZFNMatrix<9,REAL> rigcouple(3, 3,0.);
        rigcouple(0,0) = weight;
        rigcouple(1,1) = weight;
        rigcouple(0,2) = -(x[1]-xcenter[1])*weight;
        rigcouple(1,2) = (x[0]-xcenter[0])*weight;
        rigcouple(2,0) = rigcouple(0,2);
        rigcouple(2,1) = rigcouple(1,2);
        rigcouple(2,2) = ((x[0]-xcenter[0])*(x[0]-xcenter[0])+(x[1]-xcenter[1])*(x[1]-xcenter[1]))*weight;
        for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
            ek(2*phr+i,2*phr+3+j) += -rigcouple(i,j);
            ek(2*phr+3+i,2*phr+j) += -rigcouple(j,i);
        }
    }
#ifdef PZDEBUG
    if(fForcingFunction && fAnisotropicForcingFunction) {
        PZError << "TPZElasticity2D::Contribute error both forcing functions are set\n";
        DebugStop();
    }
#endif
    TPZManVector<STATE,3> floc(ff);
    if(fAnisotropicForcingFunction) {
        TPZManVector<STATE,3> res(3,0.);
        fAnisotropicForcingFunction(datavec[0].x, fConstitutiveLaw, res);
        floc[0] = res[0];
        floc[1] = res[1];
    }
    else if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
        TPZManVector<STATE,3> res(3,0.);
        fForcingFunction(datavec[0].x,res);
        floc[0] = res[0];
        floc[1] = res[1];
        floc[2] = res[2];
    }
    {
        int64_t efc=ef.Cols();
        for( int64_t in = 0; in < phr; in++ ) {
            for (int col = 0; col < efc; col++)
            {
                ef(2*in, col) += weight * (floc[0]*phi(in,0) - dphixyz(0,in)*fPreStressXX - dphixyz(1,in)*fPreStressXY);  // direcao x
                ef(2*in+1, col) += weight * (floc[1]*phi(in,0) - dphixyz(0,in)*fPreStressXY - dphixyz(1,in)*fPreStressYY);// direcao y <<<----
            }
        }

    }

    //equacoes de restricao de pressao media
    if(datavec.size() >4) {
        DebugStop();
    }
}

void TPZHybridElasticity2D::ContributeBC(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &axes = data.axes;
    int64_t phr = phi.Rows();

    TPZManVector<REAL,2> v2 = bc.Val2();
    int dim = 2;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE, 4> rhs_val(2);
        TPZFNMatrix<4, STATE> mat_val(dim, 2), deform_voight(3,1), sigma_voight(3,1), sigma(2,2);
        bc.ForcingFunctionBC()(data.x, rhs_val, mat_val);
        deform_voight(0,0) = mat_val(0,0);
        deform_voight(1,0) = mat_val(1,1);
        deform_voight(2,0) = mat_val(0,1)+mat_val(1,0);
        TPZFNMatrix<9,REAL> D(3, 3);
        if(fPlaneStress) D_plane_stress(fE_def, fnu_def, D);
        else D_plane_strain(fE_def, fnu_def, D);
        D.Multiply(deform_voight, sigma_voight);
        sigma(0,0) = sigma_voight(0,0);
        sigma(1,1) = sigma_voight(1,0);
        sigma(0,1) = sigma(1,0) = sigma_voight(2,0);
        TPZManVector<REAL,3> normal(3,0.);
        for (int i = 0; i < dim; i++) {
            normal[i] = data.normal[i];
        }
        if(bc.Type() == 0) {
            v2 = rhs_val;
        } else if(bc.Type() == 1) {
            v2[0] = 0.;
            v2[1] = 0.;
            for (int i = 0; i < dim; i++) {
                v2[0] += sigma(0,i) * normal[i];
                v2[1] += sigma(1,i) * normal[i];
            }
        } else if(bc.Type() == 2) {
           DebugStop(); // I have to think how to implement this
        }
    }

    switch (bc.Type()) {
        case 1 : // Neumann condition
            for (int in = 0; in < phr; in++) {
                for(int idf=0; idf<2; idf++) {
                    ef(2*in+idf, 0) += (STATE) (TPZMaterial::fBigNumber * phi(in, 0) * weight) * v2[idf];
                    for (int jn = 0; jn < phr; jn++) {
                        int jdf = idf;
                        ek(2*in+idf, 2*jn+jdf) += TPZMaterial::fBigNumber * phi(in, 0) * phi(jn, 0) * weight;
                    }
                }
            }
            break;
        case 0 : // Dirichlet condition
            for (int in = 0; in < phi.Rows(); in++) {
                for(int idf = 0; idf<2; idf++) {
                    ef(2*in+idf, 0) += v2[idf] * (STATE) (phi(in, 0) * weight);
                }
            }
            break;
        default:
            PZError << __PRETTY_FUNCTION__
                    << "\nBoundary condition type not implemented. Please use one of the following:\n"
                    << "\t 0: Dirichlet\n"
                    << "\t 1: Neumann\n"
                    << "\t 2: Robin\n";
            DebugStop();
    }
}


void TPZHybridElasticity2D::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc)
{
    TPZFMatrix<REAL>  &phi_u = datavec[1].phi;
    TPZFMatrix<REAL>  &phi_flux = datavec[0].phi;
    //    TPZFMatrix<REAL> &axes = data.axes;
    int phr_primal = phi_u.Rows();
    int phr_hybrid = phi_flux.Rows();
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
    int dim = 2;
    TPZManVector<REAL,2> v2 = bc.Val2();

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE, 4> rhs_val(2);
        TPZFNMatrix<4, STATE> mat_val(dim, 2), deform_voight(3,1), sigma_voight(3,1), sigma(2,2);
        bc.ForcingFunctionBC()(x, rhs_val, mat_val);
        deform_voight(0,0) = mat_val(0,0);
        deform_voight(1,0) = mat_val(1,1);
        deform_voight(2,0) = mat_val(0,1)+mat_val(1,0);
        TPZFNMatrix<9,REAL> D(3, 3);
        if(fPlaneStress) D_plane_stress(fE_def, fnu_def, D);
        else D_plane_strain(fE_def, fnu_def, D);
        D.Multiply(deform_voight, sigma_voight);
        sigma(0,0) = sigma_voight(0,0);
        sigma(1,1) = sigma_voight(1,0);
        sigma(0,1) = sigma(1,0) = sigma_voight(2,0);
        if(bc.Type() == 0) {
            v2 = rhs_val;
        } else if(bc.Type() == 1) {
            v2[0] = 0.;
            v2[1] = 0.;
            for (int i = 0; i < dim; i++) {
                v2[0] += sigma(0,i) * normal[i];
                v2[1] += sigma(1,i) * normal[i];
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

void TPZHybridElasticity2D::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    if(!fExactSol) return;

    TPZElasticity2D::Errors(data[1],errors);
    return;
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE> u_exact(1);
    TPZFNMatrix<9,STATE> du_exact;


    if(this->fExactSol){

        this->fExactSol(data[1].x,u_exact,du_exact);
    }

    REAL pressure = data[1].sol[0][0];

    // errors[0] norm L2 || u ||_l2

    errors[0] = (pressure-u_exact[0])*(pressure-u_exact[0]);//exact error pressure

    // errors[1] Semi norm H1 || grad u ||_l2

    TPZManVector<STATE,3> sol(1),dsol(3,0.);

    TPZFMatrix<REAL> &dsolaxes = data[1].dsol[0];
    TPZFNMatrix<9,REAL> flux(3,0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, flux, data[1].axes);
    errors[1] = 0.;
    int dim = 2;
    for(int id=0; id<dim; id++) {
        REAL diff = fabs(flux(id,0) - du_exact(id,0));
        errors[1]  += diff*diff;
    }

    // error[2] H1 norm

    errors[2] = errors[0] +errors[1];

    // error[3] Energy norm || u ||_e = a(u,u)= int_K K gradu.gradu dx

    STATE KPerm = 1.;


    TPZFNMatrix<9,REAL> gradpressure(dim,1),Kgradu(dim,1);
    for (int i=0; i<dim; i++) {
        gradpressure(i,0) = du_exact(i,0);
        Kgradu(i,0) = gradpressure(0)*KPerm;
    }

    REAL energy = 0.;
    for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++) {
            double cperm =0.;
            if(i==j)
                cperm = KPerm;
            energy += cperm*fabs(flux(j,0) - du_exact(j,0))*fabs(flux(i,0) - du_exact(i,0));
        }
    }

    errors[3] = energy;
}

void TPZHybridElasticity2D::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {

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

/**@}*/
