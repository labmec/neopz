//
//  TPZMatLaplacianHybrid.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#include "TPZMatLaplacianHybrid.h"
#include "pzbndcond.h"
#include "pzaxestools.h"

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid(int matid, int dim)
: TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian(matid,dim)
{
    
}

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid() :
TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian()
{
    
}

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid(const TPZMatLaplacian &copy) :
TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian(copy)
{
    
}

TPZMatLaplacianHybrid::~TPZMatLaplacianHybrid()
{
    
}

TPZMatLaplacianHybrid &TPZMatLaplacianHybrid::operator=(const TPZMatLaplacianHybrid &copy)
{
    TPZMatLaplacian::operator=(copy);
    return *this;
}

TPZMaterial *TPZMatLaplacianHybrid::NewMaterial()
{
    return new TPZMatLaplacianHybrid(*this);
}

int TPZMatLaplacianHybrid::ClassId() const
{
    return Hash("TPZMatLaplacianHybrid") ^ TPZMatLaplacian::ClassId() << 1;
}


void TPZMatLaplacianHybrid::Write(TPZStream &buf, int withclassid) const
{
    TPZMatLaplacian::Write(buf,withclassid);
}

void TPZMatLaplacianHybrid::Read(TPZStream &buf, void *context)
{
    TPZMatLaplacian::Read(buf,context);
}

int TPZMatLaplacianHybrid::VariableIndex(const std::string &name)
{

    if(name == "Pressure") return 44;
    if(name == "PressureExact") return 45;

    if(name == "Flux") return 10;
    if(name == "ExactFlux") return 13;

    if(name == "ExactFluxShiftedOrigin") return 23;

    return -1;
}

int TPZMatLaplacianHybrid::NSolutionVariables(int var){
    if(var == 44 || var==45) return 1;
    if(var == 10 || var==13 || var == 23) return fDim;
    
    else{
        
        return TPZMatLaplacian::NSolutionVariables(var);
        return 0;
        
    }
}


void TPZMatLaplacianHybrid::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
     datavec[1] L2 mesh (phi's)
     datavec[0] Hdiv mesh,
     datavec[2] Interface Mesh
     datavec[3] Interface Mesh
     
     Implement the matrix
     |Sk Ck^T |  = |f1|
     |Ck  0   |    |f2|
     Sk = int_K K graduk.gradv dx = int_K K gradphi_i.gradphi_j dx
     CK = int_partialK lambda_k*uk dx = int_K phi_i dx
     f1 = int_K f*v dx = int_K f*phi_j dx
     ck = int_partialK phi_i*mu_j dx
     f2 = int_partialK g*mu_j dx
     
     **/
    
    TPZFMatrix<REAL>  &phi = datavec[1].phi;
    TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    TPZVec<REAL>  &x = datavec[1].x;

    int phr = phi.Rows();
    
    STATE fXfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        //fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fForcingFunction->Execute(x,res);
        fXfLoc = res[0];
    }
    
    STATE KPerm = fK;
    if (fPermeabilityFunction) {
        TPZFNMatrix<9,STATE> perm, invperm;
        TPZManVector<STATE,3> func;
        TPZFNMatrix<18,STATE> dfunc(6,3,0.);
        fPermeabilityFunction->Execute(x, func, dfunc);
        KPerm = dfunc(0,0);
    }
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        
        //matrix Sk
        for( int jn = 0; jn < phr; jn++ ) {
            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight*(KPerm*(STATE)(dphi(kd,in)*dphi(kd,jn)));
            }
        }
    }
    for (int in =0; in < phr; in++) {
        ek(phr,in) += weight*phi(in,0);//lambda*phi
        ek(in,phr) += weight*phi(in,0);
    }
    //equacoes de restricao de pressao media
    ek(phr,phr+1) -= weight;
    ek(phr+1,phr) -= weight;
    
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry(1.e-10) ) std::cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << std::endl;
    }
}

void TPZMatLaplacianHybrid::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL>  &phi = datavec[1].phi;
    TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    TPZVec<REAL>  &x = datavec[1].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE fXfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        //fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fForcingFunction->Execute(x,res);
        fXfLoc = res[0];
    }
    
    STATE KPerm = fK;
    if (fPermeabilityFunction) {
        TPZFNMatrix<9,STATE> perm, invperm;
        TPZManVector<STATE,3> func;
        TPZFNMatrix<18,STATE> dfunc(6,3,0.);
        fPermeabilityFunction->Execute(x, func, dfunc);
        KPerm = dfunc(0,0);
    }
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        for(kd=0; kd<fDim; kd++) {
            ef(in,0) -= (STATE)weight*(fK*(STATE)(dphi(kd,in)*datavec[1].dsol[0](kd,0)));
        }
    }
    ef(phr,0) += weight*(-datavec[1].sol[0][0]+ datavec[3].sol[0][0]);
    ef(phr+1,0) += weight*datavec[2].sol[0][0];

}

void TPZMatLaplacianHybrid::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<REAL>  &phi_u = datavec[1].phi;
    TPZFMatrix<REAL>  &phi_flux = datavec[0].phi;
    //    TPZFMatrix<REAL> &axes = data.axes;
    int phr_primal = phi_u.Rows();
    int phr_hybrid = phi_flux.Rows();
    bool primal = true;
    TPZManVector<REAL,3> x(3);
    if(phr_hybrid)
    {
        primal = false;
        x = datavec[0].x;
    }
    else
    {
        x = datavec[1].x;
    }
    short in,jn;
    STATE v2[1];
    v2[0] = bc.Val2()(0,0);
    
    if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
        TPZManVector<STATE> res(1);
        TPZFNMatrix<3,STATE> dres(3,1);
        bc.ForcingFunction()->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        v2[0] = res[0];
    }
    
    if(primal)
    {
        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in,0) += (STATE)(gBigNumber* phi_u(in,0) * weight) * v2[0];
                    for (jn = 0 ; jn < phr_primal; jn++) {
                        ek(in,jn) += gBigNumber * phi_u(in,0) * phi_u(jn,0) * weight;
                    }
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in,0) += v2[0] * (STATE)(phi_u(in,0) * weight);
                }
                break;
            case 2 :        // mixed condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in, 0) += v2[0] * (STATE)(phi_u(in, 0) * weight);
                    for (jn = 0 ; jn < phi_u.Rows(); jn++) {
                        ek(in,jn) += bc.Val1()(0,0) * (STATE)(phi_u(in,0) * phi_u(jn,0) * weight);     // peso de contorno => integral de contorno
                    }
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
                    ef(in,0) += v2[0] * (STATE)(phi_flux(in,0) * weight);
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phr_hybrid; in++) {
                    ef(in,0) += (STATE)(gBigNumber* phi_flux(in,0) * weight) * v2[0];
                    for (jn = 0 ; jn < phr_hybrid; jn++) {
                        ek(in,jn) += gBigNumber * phi_flux(in,0) * phi_flux(jn,0) * weight;
                    }
                }
                break;
            case 2 :        // mixed condition
                DebugStop();
                for(in = 0 ; in < phr_hybrid; in++) {
                    ef(in, 0) += v2[0] * (STATE)(phi_flux(in, 0) * weight);
                    for (jn = 0 ; jn < phi_flux.Rows(); jn++) {
                        ek(in,jn) += 1./bc.Val1()(0,0) * (STATE)(phi_flux(in,0) * phi_flux(jn,0) * weight);     // peso de contorno => integral de contorno
                    }
                }
                break;
        }
    }
}

void TPZMatLaplacianHybrid::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    /**
     datavec[1] L2 mesh
     datavec[0] Hdiv mesh,
     datavec[2] Interface Mesh
     datavec[3] Interface Mesh
     **/
    if(var == 0)
    {
        TPZMatLaplacian::Solution(datavec[1],var,Solout);
        return;
    }
    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    
    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> grad(fDim,1,0.), fluxinv(fDim,1),gradu(fDim,1,0);//no TPZAnalytic solution grad é 3x1
    
    if(fForcingFunctionExact)
    {
        this->fForcingFunctionExact->Execute(datavec[1].x, pressexact,grad);
        
        for(int i = 1; i<fDim ; i++){
            
            gradu(i,0)=grad(i,0);
        }
        
        
    }
    
    PermTensor.Multiply(gradu, fluxinv);
    

    switch (var)
    {
  

        case 44://PressureFem
            Solout[0] = datavec[1].sol[0][0];
            break;
        case 45://Pressure Exact
            Solout[0] = pressexact[0];
            break;
        case 10:
        case 13:
        case 23:
            TPZMatLaplacian::Solution(datavec[1],var,Solout);
            break;
        default:
            DebugStop();
    }
}


void TPZMatLaplacianHybrid::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    /**
     datavec[1] L2 mesh (phi's)
     datavec[0] Hdiv mesh,
     datavec[2] Interface Mesh
     datavec[3] Interface Mesh

     error[0] = L2 norm
     error[1] = semi H1 norm
     error[2] = H1 norm
     error[3] = energy norm
     
 
     **/
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    

    
    if(this->fForcingFunctionExact){
        
        this->fForcingFunctionExact->Execute(data[1].x,u_exact,du_exact);
    }
    

    REAL pressure = data[1].sol[0][0];
    

    
    // errors[0] norm L2 || u ||_l2

    errors[0] = (pressure-u_exact[0])*(pressure-u_exact[0]);//exact error pressure

    // errors[1] Semi norm H1 || grad u ||_l2
    
    TPZManVector<STATE,3> sol(1),dsol(3,0.);
    
    TPZFMatrix<REAL> &dsolaxes = data[1].dsol[0];
    TPZFNMatrix<9,REAL> flux(3,0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, flux, data[1].axes);
    
    for(int id=0; id<fDim; id++) {
        REAL diff = fabs(flux(id,0) - du_exact(id,0));
        errors[1]  += diff*diff;
    }
    
    // error[2] H1 norm
    
    errors[2] = errors[0] +errors[1];
    
    // error[3] Energy norm || u ||_e = a(u,u)= int_K K gradu.gradu dx
 
        TPZFNMatrix<9,REAL> PermTensor = fTensorK;
        TPZFNMatrix<9,REAL> gradpressure(fDim,1),Kgradu(fDim,1);
        for (int i=0; i<fDim; i++) {
            gradpressure(i,0) = du_exact(i,0);
        }
        PermTensor.Multiply(gradpressure,Kgradu);


    
    REAL energy = 0.;
    for (int i=0; i<fDim; i++) {
        for (int j=0; j<fDim; j++) {
            energy += PermTensor(i,j)*fabs(flux(j,0) - du_exact(j,0))*fabs(flux(i,0) - du_exact(i,0));
        }
    }
    
    errors[3] = energy;


    
}
