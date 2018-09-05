//
//  TPZMatElasticity2D.cpp
//  PZ
//
//  Created by Omar on 10/27/14.
//
//


#include <iostream>
#include <string>
#include "TPZMatElasticity2D.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif


TPZMatElasticity2D::TPZMatElasticity2D()
: TPZRegisterClassId(&TPZMatElasticity2D::ClassId), TPZMaterial()
{
    m_E = 0.;
    m_nu = 0.;
    m_lambda = 0.;
    m_mu = 0.;
    m_f.resize(2);
    m_f[0]=0.;
    m_f[1]=0.;
    m_plane_stress = 1;
    m_s0_xx = 0.0;
    m_s0_xy = 0.0;
    m_s0_yy = 0.0;
    m_s0_zz = 0.0;
    
}

TPZMatElasticity2D::TPZMatElasticity2D(int matid)
: TPZRegisterClassId(&TPZMatElasticity2D::ClassId), TPZMaterial(matid)
{
    m_E = 0.;
    m_nu = 0.;
    m_lambda = 0.;
    m_mu = 0.;
    m_f.resize(2);
    m_f[0]=0.;
    m_f[1]=0.;
    m_plane_stress = 1;
    m_s0_xx = 0.0;
    m_s0_xy = 0.0;
    m_s0_yy = 0.0;
    m_s0_zz = 0.0;
}

TPZMatElasticity2D::TPZMatElasticity2D(int matid, REAL E, REAL nu, REAL fx, REAL fy, int plainstress)
: TPZRegisterClassId(&TPZMatElasticity2D::ClassId), TPZMaterial(matid)
{
    m_E = E;
    m_nu = nu;
    m_lambda = (E*nu)/((1+nu)*(1-2*nu));
    m_mu = E/(2*(1+nu));
    m_f.resize(2);
    m_f[0]=fx;
    m_f[1]=fy;
    m_plane_stress = plainstress;
    m_s0_xx = 0.0;
    m_s0_xy = 0.0;
    m_s0_yy = 0.0;
    m_s0_zz = 0.0;
}

TPZMatElasticity2D::~TPZMatElasticity2D()
{
}


TPZMatElasticity2D::TPZMatElasticity2D(const TPZMatElasticity2D &copy)
: TPZRegisterClassId(&TPZMatElasticity2D::ClassId),  TPZMaterial(copy)
{
    m_E = copy.m_E;
    m_nu = copy.m_nu;
    m_lambda = copy.m_lambda;
    m_mu = copy.m_mu;
    m_f.resize(copy.m_f.size());
    for (int i = 0; i < copy.m_f.size(); i++) {
        m_f[i] = copy.m_f[i];
    }
    m_plane_stress = copy.m_plane_stress;
    m_s0_xx = copy.m_s0_xx;
    m_s0_xy = copy.m_s0_xy;
    m_s0_yy = copy.m_s0_yy;
    m_s0_zz = copy.m_s0_zz;
}

TPZMatElasticity2D & TPZMatElasticity2D::operator=(const TPZMatElasticity2D &copy)
{
	TPZMaterial::operator = (copy);
    m_E = copy.m_E;
    m_nu = copy.m_nu;
    m_lambda = copy.m_lambda;
    m_mu = copy.m_mu;
    m_s0_xx = copy.m_s0_xx;
    m_s0_xy = copy.m_s0_xy;
    m_s0_yy = copy.m_s0_yy;
    m_s0_zz = copy.m_s0_zz;
    m_f.resize(copy.m_f.size());
    for (int i = 0; i < copy.m_f.size(); i++) {
        m_f[i] = copy.m_f[i];
    }
    m_plane_stress = copy.m_plane_stress;
    return *this;
}

int TPZMatElasticity2D::NStateVariables() {
    return 2;
}


void TPZMatElasticity2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef) {
    
   
    if (data.fShapeType == TPZMaterialData::EVecShape) {
        ContributeVec(data,weight,ek,ef);
        return;
    }
    // Getting weight functions
    TPZFMatrix<REAL>  &phi_u     =  data.phi;
    TPZFMatrix<REAL> &dphi_u     =  data.dphix;
    int n_phi_u = phi_u.Rows();
    int first_u  = 0;
    
    TPZFNMatrix<40,REAL> grad_phi_u(3,n_phi_u);
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, grad_phi_u, data.axes);
    TPZFMatrix<STATE> dsol_u    = data.dsol[0];
    
    TPZManVector<STATE,3> grad_p(m_f);
    
    if(this->HasForcingFunction())
    {
        fForcingFunction->Execute(data.x,grad_p);
    }
    
    REAL dvdx,dvdy,dudx,dudy;
    REAL duxdx,duxdy,duydx,duydy;
    
    //  Gradient for ux
    duxdx = dsol_u(0,0)*data.axes(0,0)+dsol_u(1,0)*data.axes(1,0); // dux/dx
    duxdy = dsol_u(0,0)*data.axes(0,1)+dsol_u(1,0)*data.axes(1,1); // dux/dy
    
    //  Gradient for uy
    duydx = dsol_u(0,1)*data.axes(0,0)+dsol_u(1,1)*data.axes(1,0); // duy/dx
    duydy = dsol_u(0,1)*data.axes(0,1)+dsol_u(1,1)*data.axes(1,1); // duy/dy
    
    for(int iu = 0; iu < n_phi_u; iu++ )
    {
        dvdx = grad_phi_u(0,iu);
        dvdy = grad_phi_u(1,iu);
        
        ef(2*iu + first_u)     +=    weight * (grad_p[0] * phi_u(iu, 0) - (dvdx*m_s0_xx + dvdy*m_s0_xy) );    // x direction
        ef(2*iu+1 + first_u)   +=    weight * (grad_p[1] * phi_u(iu, 0) - (dvdx*m_s0_xy + dvdy*m_s0_yy) );    // y direction
        
        if (m_plane_stress == 1)
        {
            /* Plain stress state */
            ef(2*iu + first_u)           += weight*((4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu))*dvdx*duxdx      + (2*m_mu)*dvdy*duxdy);
            
            ef(2*iu + first_u)           += weight*((2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu))*dvdx*duydy         + (2*m_mu)*dvdy*duydx);
            
            ef(2*iu+1 + first_u)         += weight*((2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu))*dvdy*duxdx         + (2*m_mu)*dvdx*duxdy);
            
            ef(2*iu+1 + first_u)         += weight*((4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu))*dvdy*duydy     + (2*m_mu)*dvdx*duydx);
        }
        else
        {
            /* Plain Strain State */
            ef(2*iu + first_u)           += weight*  ((m_lambda + 2*m_mu)*dvdx*duxdx  + (m_mu)*dvdy*(duxdy));
            
            ef(2*iu + first_u)           += weight*  (m_lambda*dvdx*duydy            + (m_mu)*dvdy*(duydx));
            
            ef(2*iu+1 + first_u)         += weight*  (m_lambda*dvdy*duxdx            + (m_mu)*dvdx*(duxdy));
            
            ef(2*iu+1 + first_u)         += weight*  ((m_lambda + 2*m_mu)*dvdy*duydy  + (m_mu)*dvdx*(duydx));
        }
        
        for(int ju = 0; ju < n_phi_u; ju++)
        {
            
            dudx = grad_phi_u(0,ju);
            dudy = grad_phi_u(1,ju);
            
            if (this->m_plane_stress == 1)
            {
                /* Plain stress state */
                ek(2*iu + first_u, 2*ju + first_u)	     += weight*((4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu))*dvdx*dudx		+ (m_mu)*dvdy*dudy);
                
                ek(2*iu + first_u, 2*ju+1 + first_u)       += weight*((2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu))*dvdx*dudy			+ (m_mu)*dvdy*dudx);
                
                ek(2*iu+1 + first_u, 2*ju + first_u)       += weight*((2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu))*dvdy*dudx			+ (m_mu)*dvdx*dudy);
                
                ek(2*iu+1 + first_u, 2*ju+1 + first_u)     += weight*((4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu))*dvdy*dudy		+ (m_mu)*dvdx*dudx);
            }
            else
            {
                /* Plain Strain State */
                ek(2*iu + first_u,2*ju + first_u)         += weight*	((m_lambda + 2*m_mu)*dvdx*dudx	+ (m_mu)*dvdy*dudy);
                
                ek(2*iu + first_u,2*ju+1 + first_u)       += weight*	(m_lambda*dvdx*dudy			+ (m_mu)*dvdy*dudx);
                
                ek(2*iu+1 + first_u,2*ju + first_u)       += weight*	(m_lambda*dvdy*dudx			+ (m_mu)*dvdx*dudy);
                
                ek(2*iu+1 + first_u,2*ju+1 + first_u)     += weight*	((m_lambda + 2*m_mu)*dvdy*dudy	+ (m_mu)*dvdx*dudx);
                
            }
        }
    }

}

void TPZMatElasticity2D::ContributeVec(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if (data.fShapeType != TPZMaterialData::EVecShape) {
        DebugStop();
    }
    
    // Getting weight functions
    TPZFMatrix<REAL> &dphiU     =  data.dphix;
    int phrU = dphiU.Cols();
    
    TPZFNMatrix<200,REAL> dudaxes(2,dphiU.Cols()), dvdaxes(2,dphiU.Cols()), dudx(3,phrU), dvdx(3,phrU);
    for (int i=0; i<2; i++) {
        for (int j=0; j<phrU; j++) {
            dudaxes(i,j) = dphiU(i,j);
            dvdaxes(i,j) = dphiU(2+i,j);
        }
    }
    TPZAxesTools<REAL>::Axes2XYZ(dudaxes, dudx, data.axes);
    TPZAxesTools<REAL>::Axes2XYZ(dvdaxes, dvdx, data.axes);


    for(int iu = 0; iu < phrU; iu++ )
    {
        TPZFNMatrix<4,REAL> gradv(2,2);
        gradv(0,0) = dudx(0,iu);
        gradv(0,1) = dudx(1,iu);
        gradv(1,0) = dvdx(0,iu);
        gradv(1,1) = dvdx(1,iu);
        
        for(int ju = 0; ju < phrU; ju++)
        {
            TPZFNMatrix<4,REAL> gradu(2,2);
            gradu(0,0) = dudx(0,ju);
            gradu(0,1) = dudx(1,ju);
            gradu(1,0) = dvdx(0,ju);
            gradu(1,1) = dvdx(1,ju);
            
            if (this->m_plane_stress == 1)
            {
                /* Plain stress state
                 \sigma_x = E/(1-\nu\nu) (\epsilon_x + \nu \epsilon_y)
                 \sigma_x = \frac{4\mu(\lambda+\mu)}{\lambda+2\mu)}\epsilon_x + \frac{2\mu\lambda}{\lambda+2\mu} \epsilon_y
                 \sigma_y = E/(1-\nu\nu) (\epsilon_y + \nu \epsilon_x)
                 \sigma_y = \frac{4\mu(\lambda+\mu)}{\lambda+2\mu)}\epsilon_y + \frac{2\mu\lambda}{\lambda+2\mu} \epsilon_x
                 \tau_{xy} = \frac{E}{1+\nu} \epsilon_{xy}
                 \tau_{xy} = \frac{1}{2\mu} \epsilon_{xy}
                 */
                TPZFNMatrix<4,REAL> sigma_u(2,2);
                sigma_u(0,0) = m_E/(1-m_nu*m_nu) *(gradu(0,0)+m_nu*gradu(1,1));
                sigma_u(1,1) = m_E/(1-m_nu*m_nu) *(gradu(1,1)+m_nu*gradu(0,0));
                sigma_u(0,1) = m_E/(2.*(1+m_nu))*(gradu(0,1)+gradu(1,0));
                sigma_u(1,0) = sigma_u(0,1);
                ek(iu, ju)         += weight*(sigma_u(0,0)*gradv(0,0)+sigma_u(1,1)*gradv(1,1)+
                                              sigma_u(1,0)*gradv(1,0)+sigma_u(0,1)*gradv(0,1));
            }
            else
            {
                /* Plain Strain State */
                
                /*
                 \sigma_x = \frac{E}{(1+\nu)(1-2\nu)}((1-\nu)\epsilon_x + \nu\epsilon_y)
                 \sigma_x = (\lambda+2\mu)\epsilon_x + \lambda\epsilon_y
                 \sigma_y = \frac{E}{(1+\nu)(1-2\nu)}((1-\nu)\epsilon_y + \nu\epsilon_x)
                 \sigma_y = (\lambda+2\mu)\epsilon_y + \lambda\epsilon_x
                 \tau_{xy} = \frac{E}{1+\nu} \epsilon_{xy}
                 \tau_{xy} = \frac{1}{2\mu} \epsilon_{xy}
                 */
                TPZFNMatrix<4,REAL> sigma_u(2,2);
                sigma_u(0,0) = (m_lambda+2.*m_mu)*gradu(0,0)+m_lambda*gradu(1,1);
                sigma_u(1,1) = (m_lambda+2.*m_mu)*gradu(1,1)+m_lambda*gradu(0,0);
                sigma_u(0,1) = m_mu*(gradu(0,1)+gradu(1,0));
                sigma_u(1,0) = sigma_u(0,1);
                STATE energy = (sigma_u(0,0)*gradv(0,0)+sigma_u(1,1)*gradv(1,1)+
                                sigma_u(1,0)*gradv(1,0)+sigma_u(0,1)*gradv(0,1));
                ek(iu, ju)  += weight*energy;
                
            }
        }
    }
    //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    this->ContributeVec(data,weight,ef);

}
void TPZMatElasticity2D::ContributeVec(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    // Getting weight functions
    TPZFMatrix<REAL>  &phi_u =  data.phi;
    int n_phi_u = phi_u.Rows()/2;
    TPZFMatrix<REAL> &dphiU     =  data.dphix;

    TPZFNMatrix<200,REAL> dudaxes(2,dphiU.Cols()), dvdaxes(2,dphiU.Cols()), dudx(2,n_phi_u), dvdx(2,n_phi_u);
    for (int i=0; i<2; i++) {
        for (int j=0; j<n_phi_u; j++) {
            dudaxes(i,j) = dphiU(i,j);
            dvdaxes(i,j) = dphiU(2+i,j);
        }
    }
    TPZAxesTools<REAL>::Axes2XYZ(dudaxes, dudx, data.axes);
    TPZAxesTools<REAL>::Axes2XYZ(dvdaxes, dvdx, data.axes);
    
    TPZManVector<STATE,3> sol_u =data.sol[0];
    TPZFNMatrix<4,STATE> dsol_xy(2,2), dsol_u = data.dsol[0];
    
    TPZAxesTools<STATE>::Axes2XYZ(dsol_u, dsol_xy, data.axes);
    
    TPZManVector<STATE,3> p(m_f);
//    TPZFNMatrix<4,STATE> grad_p(2,2,0.0);
    
    if(this->HasForcingFunction())
    {
        fForcingFunction->Execute(data.x,p);
    }
    
    for(int iu = 0; iu < n_phi_u; iu++ )
    {
        TPZFNMatrix<4,REAL> gradv(2,2);
        gradv(0,0) = dudx(0,iu);
        gradv(0,1) = dudx(1,iu);
        gradv(1,0) = dvdx(0,iu);
        gradv(1,1) = dvdx(1,iu);

        //          Vector Force right hand term
        ef(iu)     +=    weight*(p[0]*phi_u(2*iu, 0) + p[1]*phi_u(2*iu+1,0)
                                 -m_s0_xx*gradv(0,0)-m_s0_yy*gradv(1,1)-m_s0_xy*(gradv(0,1)+gradv(1,0)));
        
        if (m_plane_stress == 1)
        {
            /* Plain stress state */
            TPZFNMatrix<4,REAL> sigma_u(2,2);
            sigma_u(0,0) = m_E/(1-m_nu*m_nu) *(dsol_xy(0,0)+m_nu*dsol_xy(1,1));
            sigma_u(1,1) = m_E/(1-m_nu*m_nu) *(dsol_xy(1,1)+m_nu*dsol_xy(0,0));
            sigma_u(0,1) = m_E/(2.*(1+m_nu))*(dsol_xy(0,1)+dsol_xy(1,0));
            sigma_u(1,0) = sigma_u(0,1);

            ef(iu) -= weight*(gradv(0,0)*sigma_u(0,0)+gradv(1,1)*sigma_u(1,1)+gradv(1,0)*sigma_u(1,0)+gradv(0,1)*sigma_u(0,1));
        }
        else
        {
            /* Plain Strain State */
            TPZFNMatrix<4,REAL> sigma_u(2,2);
            sigma_u(0,0) = (m_lambda+2.*m_mu)*dsol_xy(0,0)+m_lambda*dsol_xy(1,1);
            sigma_u(1,1) = (m_lambda+2.*m_mu)*dsol_xy(1,1)+m_lambda*dsol_xy(0,0);
            sigma_u(0,1) = m_mu*(dsol_xy(0,1)+dsol_xy(1,0));
            sigma_u(1,0) = sigma_u(0,1);
            ef(iu) -= weight*(gradv(0,0)*sigma_u(0,0)+gradv(1,1)*sigma_u(1,1)+gradv(1,0)*sigma_u(1,0)+gradv(0,1)*sigma_u(0,1));

        }
    }
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    

}

void TPZMatElasticity2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) {
    
    if (data.fShapeType == TPZMaterialData::EVecShape) {
        ContributeVec(data, weight, ef);
        return;
    }
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phi_u =  data.phi;
    int n_phi_u = phi_u.Rows();
    int first_u  = 0;
    
    TPZFNMatrix<40,REAL> grad_phi_u(3,n_phi_u);
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, grad_phi_u, data.axes);
    TPZFMatrix<STATE> dsol_u    = data.dsol[0];
    
    TPZManVector<STATE,3> grad_p(m_f);
    
    if(this->HasForcingFunction())
    {
        fForcingFunction->Execute(data.x,grad_p);
    }
    
    REAL dvdx,dvdy;
    REAL duxdx,duxdy,duydx,duydy;
    
    //  Gradient for ux
    duxdx = dsol_u(0,0)*data.axes(0,0)+dsol_u(1,0)*data.axes(1,0); // dux/dx
    duxdy = dsol_u(0,0)*data.axes(0,1)+dsol_u(1,0)*data.axes(1,1); // dux/dy
    
    //  Gradient for uy
    duydx = dsol_u(0,1)*data.axes(0,0)+dsol_u(1,1)*data.axes(1,0); // duy/dx
    duydy = dsol_u(0,1)*data.axes(0,1)+dsol_u(1,1)*data.axes(1,1); // duy/dy

    
    for(int iu = 0; iu < n_phi_u; iu++ )
    {
        dvdx = grad_phi_u(0,iu);
        dvdy = grad_phi_u(1,iu);
        
         ef(2*iu + first_u)     +=    weight * (grad_p[0] * phi_u(iu, 0) - (dvdx*m_s0_xx + dvdy*m_s0_xy) );    // x direction
         ef(2*iu+1 + first_u)   +=    weight * (grad_p[1] * phi_u(iu, 0) - (dvdx*m_s0_xy + dvdy*m_s0_yy) );    // y direction
        
        if (m_plane_stress == 1)
        {
            /* Plain stress state */
            ef(2*iu + first_u)           += weight*((4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu))*dvdx*duxdx      + (2*m_mu)*dvdy*duxdy);
            
            ef(2*iu + first_u)           += weight*((2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu))*dvdx*duydy         + (2*m_mu)*dvdy*duydx);
            
            ef(2*iu+1 + first_u)         += weight*((2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu))*dvdy*duxdx         + (2*m_mu)*dvdx*duxdy);
            
            ef(2*iu+1 + first_u)         += weight*((4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu))*dvdy*duydy     + (2*m_mu)*dvdx*duydx);
        }
        else
        {
            /* Plain Strain State */
            ef(2*iu + first_u)           += weight*  ((m_lambda + 2*m_mu)*dvdx*duxdx  + (m_mu)*dvdy*(duxdy));
            
            ef(2*iu + first_u)           += weight*  (m_lambda*dvdx*duydy            + (m_mu)*dvdy*(duydx));
            
            ef(2*iu+1 + first_u)         += weight*  (m_lambda*dvdy*duxdx            + (m_mu)*dvdx*(duxdy));
            
            ef(2*iu+1 + first_u)         += weight*  ((m_lambda + 2*m_mu)*dvdy*duydy  + (m_mu)*dvdx*(duydx));
        }
    }
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    
}

/// compute the stress tensor as a function of the solution gradient
void TPZMatElasticity2D::ComputeSigma(const TPZFMatrix<STATE> &dudx, TPZFMatrix<STATE> &sigma)
{
#ifdef PZDEBUG
    if (dudx.Rows() < 2 || dudx.Cols() != 2 || sigma.Rows() != 2 || sigma.Cols() != 2) {
        DebugStop();
    }
#endif

    if (m_plane_stress == 1)
    {
        sigma(0,0) = (4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu))*dudx.g(0,0)+(2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu))*dudx.g(1,1);
        sigma(0,1) = (m_mu)*(dudx.g(1,0)+dudx.g(0,1));
        sigma(1,0) = sigma(0,1);
        sigma(1,1) = (2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu))*dudx.g(0,0)+(4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu))*dudx.g(1,1);
        /* Plain stress state */
    }
    else
    {
        sigma(0,0) = (m_lambda + 2*m_mu)*dudx.g(0,0)+m_lambda*dudx.g(1,1);
        sigma(1,0) = m_mu*(dudx.g(1,0)+dudx.g(0,1));
        sigma(0,1) = sigma(1,0);
        sigma(1,1) = (m_lambda + 2*m_mu)*dudx.g(1,1)+m_lambda*dudx.g(0,0);
        /* Plain Strain State */
    }

}



void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{

    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<STATE,3> sol_u = data.sol[0];
    TPZFMatrix<STATE> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    REAL uy = sol_u[1];
    
    TPZFNMatrix<4,STATE> val1loc(bc.Val1()),val2loc(bc.Val2());
    
    if (bc.HasForcingFunction()) {
        TPZManVector<STATE,2> val2vec(2);
        
        bc.ForcingFunction()->Execute(data.x, val2vec, val1loc);
        val2loc(0,0) = val2vec[0];
        val2loc(1,0) = val2vec[1];
        // we assume type 2 is displacement value is weakly imposed
        if(bc.Type() == 2)
        {
            val1loc = bc.Val1();
            for (int i=0; i<2; i++) {
                val2loc(i,0) = 0.;
                for (int j=0; j<2; j++) {
                    val2loc(i,0) += val1loc(i,j)*val2vec[j];
                }
            }
        }
        if(bc.Type() == 1)
        {
            for (int i=0; i<2; i++) {
                val2loc(i,0) = 0.;
                for (int j=0; j<2; j++) {
                    val2loc(i,0) += val1loc(i,j)*data.normal[j];
                }
            }
        }
    }
    int phru = phiu.Rows();
    short in,jn;
    STATE v2[3];
    TPZFMatrix<STATE> &v1 = val1loc;
    v2[0] = val2loc(0,0);	//	Ux displacement or Tnx
    v2[1] = val2loc(1,0);	//	Uy displacement or Tny
    
    //	Here each digit represent an individual boundary condition corresponding to each state variable.
    //	0 means Dirichlet condition on x-y
    //	1 means Neumann condition
    //	7 means Dirichlet condition on x
    //	8 means Dirichlet condition on y
    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    switch (bc.Type())
    {
        case 0 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)      += BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)       += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
        case 1 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0*v2[1]*phiu(in,0)*weight;		//	Tny
            }
            break;
        }
        case 2 :
        {
            //	Mixed condition for each state variable no used here
            //	Elasticity Equation
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += v1(i,j)*data.sol[0][j];
            }
            
            for(in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += weight * (v2[0]-res(0,0)) * phiu(in,0);
                ef(2*in+1,0) += weight * (v2[1]-res(1,0)) * phiu(in,0);
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf< this->Dimension(); jdf++)
                    {
                        ek(2*in+idf,2*jn+jdf) += v1(idf,jdf)*phiu(in,0)*phiu(jn,0)*weight;
                        //      Not Complete with val2? HERE! PHIL!!!!
                        //      DebugStop();
                    }
                }
            }
            
            break;
        }
        case 3 :
        {
            //	Null Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)      += BIGNUMBER*( v2[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= BIGNUMBER*( v2[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)       += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
        case 4 :
        {
            //	Stress Field as Neumann condition for each state variable
            //	Elasticity Equation
            
            for(in = 0; in < this->Dimension(); in ++){
                v2[in] = ( v1(in,0) * data.normal[0] + v1(in,1) * data.normal[1]);
            }
            
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;        //	Tnx
                ef(2*in+1,0)	+= -1.0*v2[1]*phiu(in,0)*weight;		//	Tny
            }
            
            break;
        }
        case 5 :
            //	Normal Pressure condition Pressure value Should be inserted in v2[0]
            //	Elasticity Equation
            {
                TPZFNMatrix<2,STATE> res(2,1,0.);
                for(int i=0; i<2; i++) for(int j=0; j<2; j++)
                {
                    res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
                }
                for(int in = 0 ; in < phru; in++)
                {
                    ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
                    ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
                    for(int jn=0; jn< phru; jn++)
                    {
                        for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf < this->Dimension(); jdf++)
                        {
                            ek(2*in+idf,2*jn+jdf) += v1(idf,jdf)*data.normal[idf]*data.normal[jdf]*phiu(in,0)*phiu(jn,0)*weight;
                            //      Not Complete with val2? HERE! PHIL!!!!
                            //      DebugStop();
                        }
                    }
                }
            }
            break;
        case 6 :
            //	Normal Pressure condition Pressure value Should be inserted in v2[0]
            //	Elasticity Equation
            {
                TPZFNMatrix<2,STATE> res(2,1,0.);
                for(int i=0; i<2; i++) for(int j=0; j<2; j++)
                {
                    res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
                }
                for(int in = 0 ; in < phru; in++)
                {
                    ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
                    ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
                    for(int jn=0; jn< phru; jn++)
                    {
                        for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf < this->Dimension(); jdf++)
                        {
                            ek(2*in+idf,2*jn+jdf) += v1(idf,jdf)*data.normal[idf]*data.normal[jdf]*phiu(in,0)*phiu(jn,0)*weight;
                            //      Not Complete
                            //      DebugStop();
                        }
                    }
                }
            }
            break;
        case 7 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)		+= BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                }
            }
            
            break;
        }
        case 8 :
        {
            //	Dirichlet condition for uy
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in+1,0)	+= BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
        default:
        {
            PZError << "TPZMatElasticity2D::ContributeBC error - Wrong boundary condition type" << std::endl;
            DebugStop();
        }
            break;
    }

}

void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<STATE,3> sol_u = data.sol[0];
    TPZFMatrix<STATE> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    REAL uy = sol_u[1];
    
    int phru = phiu.Rows();
    short in;
    STATE v2[3]; TPZFMatrix<STATE> &v1 = bc.Val1();
    v2[0] = bc.Val2()(0,0);	//	Ux displacement or Tnx
    v2[1] = bc.Val2()(1,0);	//	Uy displacement or Tny
    
    //	Here each digit represent an individual boundary condition corresponding to each state variable.
    //	0 means Dirichlet condition on x-y
    //	1 means Neumann condition
    //	7 means Dirichlet condition on x
    //	8 means Dirichlet condition on y
    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    switch (bc.Type())
    {
        case 0 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)      += BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;	// y displacement Value

            }
            
            break;
        }
        case 1 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0*v2[1]*phiu(in,0)*weight;		//	Tny
            }
            break;
        }
        case 2 :
        {
            //	Mixed condition for each state variable no used here
            //	Elasticity Equation
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += bc.Val1()(i,j)*data.sol[0][j];
            }
            
            for(in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += weight * (v2[0]-res(0,0)) * phiu(in,0);
                ef(2*in+1,0) += weight * (v2[1]-res(1,0)) * phiu(in,0);
                
            }
            
            break;
        }
        case 3 :
        {
            //	Null Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)      += BIGNUMBER*(0.0 - v2[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= BIGNUMBER*(0.0 - v2[1])*phiu(in,0)*weight;	// y displacement Value

            }
            
            break;
        }
        case 4 :
        {
            //	Stress Field as Neumann condition for each state variable
            //	Elasticity Equation
            
            for(in = 0; in < this->Dimension(); in ++){ v2[in] = ( v1(in,0) * data.normal[0] + v1(in,1) * data.normal[1]);}
            
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;        //	Tnx
                ef(2*in+1,0)	+= -1.0*v2[1]*phiu(in,0)*weight;		//	Tny
            }
            
            break;
        }
        case 5 :
            //	Normal Pressure condition Pressure value Should be inserted in v2[0]
            //	Elasticity Equation
        {
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
            }
            for(int in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
                ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
            }
        }
            break;
        case 6 :
            //	Normal Pressure condition Pressure value Should be inserted in v2[0]
            //	Elasticity Equation
        {
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
            }
            for(int in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
                ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
            }
        }
            break;
        case 7 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)		+= BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;	// X displacement Value
            }
            
            break;
        }
        case 8 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in+1,0)	+= BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;	// y displacement Value
            }
            
            break;
        }
        default:
        {
            PZError << "TPZMatElasticity2D::ContributeBC error - Wrong boundary condition type" << std::endl;
            DebugStop();
        }
            break;
    }
    
}


void TPZMatElasticity2D::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNeighborSol = true;
    data.fNeedsNeighborCenter = false;
    data.fNeedsNormal = true;
}

void TPZMatElasticity2D::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}


void TPZMatElasticity2D::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Plane Problem (m_plane_stress = 0, for Plane Strain conditions) " << m_plane_stress << std::endl;
    out << "Properties for elasticity: \n";
    out << "\t Young modulus   = "											<< m_E		<< std::endl;
    out << "\t Poisson Ratio   = "											<< m_nu		<< std::endl;
    out << "\t First Lamé Parameter   = "									<< m_lambda	<< std::endl;
    out << "\t Second Lamé Parameter   = "									<< m_mu		<< std::endl;
    out << "\t Body force vector B {X-direction, Y-direction}   = "			<< m_f[0] << ' ' << m_f[1]   << std::endl;
    out << "\t m_s0_xx   = "			<< m_s0_xx << std::endl;
    out << "\t m_s0_xy   = "			<< m_s0_xy << std::endl;
    out << "\t m_s0_yy   = "			<< m_s0_yy << std::endl;
    out << "\t m_s0_zz   = "			<< m_s0_zz << std::endl;
    out << "Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
    
}

/** Returns the variable index associated with the name */
int TPZMatElasticity2D::VariableIndex(const std::string &name)
{
    //	Elasticity Variables
    if(!strcmp("Displacement",name.c_str()))				return	1;
    if(!strcmp("SolidPressure",name.c_str()))				return	2;
    if(!strcmp("SigmaX",name.c_str()))						return	3;
    if(!strcmp("SigmaY",name.c_str()))						return	4;
    if(!strcmp("SigmaZ",name.c_str()))						return	5;
    if(!strcmp("TauXY",name.c_str()))						return	6;
    if(!strcmp("EpsX",name.c_str()))                        return    7;
    if(!strcmp("EpsY",name.c_str()))                        return    8;
    if(!strcmp("EpsZ",name.c_str()))                        return    9;
    if(!strcmp("EpsXY",name.c_str()))                        return    10;
//    PZError << "TPZMatElasticity2D::VariableIndex Error\n";
    
    return TPZMaterial::VariableIndex(name);
}

/**
 * Save the element data to a stream
 */
void TPZMatElasticity2D::Write(TPZStream &buf, int withclassid) const
{
    TPZMaterial::Write(buf,withclassid);
    buf.Write(&m_E);
    buf.Write(&m_nu);
    buf.Write(&m_lambda);
    buf.Write(&m_mu);
    buf.Write( m_f);
    buf.Write(&m_s0_xx);
    buf.Write(&m_s0_xy);
    buf.Write(&m_s0_yy);
    buf.Write(&m_s0_zz);
    buf.Write(&m_plane_stress);
    
}

/**
 * Read the element data from a stream
 */
void TPZMatElasticity2D::Read(TPZStream &buf, void *context)
{
    TPZMaterial::Read(buf,context);
    buf.Read(&m_E);
    buf.Read(&m_nu);
    buf.Read(&m_lambda);
    buf.Read(&m_mu);
    buf.Read( m_f);
    buf.Read(&m_s0_xx);
    buf.Read(&m_s0_xy);
    buf.Read(&m_s0_yy);
    buf.Read(&m_s0_zz);
    buf.Read(&m_plane_stress);
    
}

int TPZMatElasticity2D::NSolutionVariables(int var){
    if(var == 1)	return 3;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    if(var == 6)	return 1;
    if(var == 7)    return 1;
    if(var == 8)    return 1;
    if(var == 9)    return 1;
    if(var == 10)    return 1;

    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZMatElasticity2D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize(this->NSolutionVariables(var));
    
    TPZManVector<STATE,3> SolU, SolP;
    TPZFNMatrix <6,STATE> DSolU, DSolP;
    TPZFNMatrix <9> axesU, axesP;
    
    TPZVec<REAL> ptx(3);
    TPZVec<STATE> solExata(3);
    TPZFMatrix<STATE> flux(5,1);
    
    if (data.sol.size() != 1) {
        DebugStop();
    }
    
    SolU	=	data.sol[0];
    DSolU	=	data.dsol[0];
    axesU	=	data.axes;
    
    
    //	Displacements
    if(var == 1 || var == 0){
        Solout[0] = SolU[0];
        Solout[1] = SolU[1];
        if(var==1) Solout[2] = 0.0;
        return;
    }
    
    
    REAL epsx;
    REAL epsy;
    REAL epsxy;
    REAL SigX;
    REAL SigY;
    REAL SigZ;
    REAL Tau, DSolxy[2][2];
    REAL divu;
    
    DSolxy[0][0] = DSolU(0,0)*axesU(0,0)+DSolU(1,0)*axesU(1,0); // dUx/dx
    DSolxy[1][0] = DSolU(0,0)*axesU(0,1)+DSolU(1,0)*axesU(1,1); // dUx/dy
    
    DSolxy[0][1] = DSolU(0,1)*axesU(0,0)+DSolU(1,1)*axesU(1,0); // dUy/dx
    DSolxy[1][1] = DSolU(0,1)*axesU(0,1)+DSolU(1,1)*axesU(1,1); // dUy/dy
    
    divu = DSolxy[0][0]+DSolxy[1][1]+0.0;	
    
    epsx = DSolxy[0][0];// du/dx
    epsy = DSolxy[1][1];// dv/dy
    epsxy = 0.5*(DSolxy[1][0]+DSolxy[0][1]);
    REAL C11 = 4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu);
    REAL C22 = 2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu);
    
    if (this->m_plane_stress)
    {
        SigX = C11*epsx+C22*epsy;
        SigY = C11*epsy+C22*epsx;
        SigZ = 0.0;
        Tau = 2.0*m_mu*epsxy;
    }
    else
    {
        SigX = ((m_lambda + 2*m_mu)*(epsx) + (m_lambda)*epsy);
        SigY = ((m_lambda + 2*m_mu)*(epsy) + (m_lambda)*epsx);
        SigZ = m_lambda*divu;
        Tau = 2.0*m_mu*epsxy;
    }
    
    
    //	Hydrostatic stress
    if(var == 2) 
    {
        Solout[0] = SigX+SigY+SigZ;
        return;
    }
    
    //	Effective Stress x-direction
    if(var == 3) {
        Solout[0] = SigX + m_s0_xx;
        return;
    }
    
    //	Effective Stress y-direction	
    if(var == 4) {
        Solout[0] = SigY + m_s0_yy;
        return;
    }
    
    //	Effective Stress y-direction
    if(var == 5) {
        Solout[0] = SigZ + m_s0_zz;
        return;
    }
    
    //	Shear Stress	
    if(var == 6) {
        Solout[0] = Tau + m_s0_xy;
        return;
    }
    
    // epsx
    if (var == 7) {
        Solout[0] = epsx;
    }
    
    // epsy
    if (var == 8) {
        Solout[0] = epsy;
    }
    
    // epsz
    if (var == 9) {
        if (m_plane_stress) {
            Solout[0] = -m_nu*(epsx+epsy);
        }
        else
        {
            Solout[0] = 0.;
        }
    }
    
    // epsxy
    if (var == 10) {
        Solout[0] = epsxy;
    }
}

int TPZMatElasticity2D::ClassId() const{
    return Hash("TPZMatElasticity2D") ^ TPZMaterial::ClassId() << 1;
}

/**
 * @brief Computes the error due to the difference between the interpolated flux \n
 * and the flux computed based on the derivative of the solution
 */
void TPZMatElasticity2D::Errors(TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
                    TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
                    TPZVec<STATE> &uexact, TPZFMatrix<STATE> &duexact,
                    TPZVec<REAL> &val)
{
    TPZFNMatrix<9,STATE> dudx(2,2), stress(2,2), stressexact(2,2);
    TPZAxesTools<STATE>::Axes2XYZ(dsol, dudx, axes);
    ComputeSigma(dudx, stress);
    ComputeSigma(duexact, stressexact);
    REAL L2 = 0.;
    L2 = (sol[0]-uexact[0])*(sol[0]-uexact[0])+(sol[1]-uexact[1])*(sol[1]-uexact[1]);
    REAL H1 = 0.;
    REAL energy = 0.;
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++) {
            H1 += (dudx(i,j)-duexact(i,j))*(dudx(i,j)-duexact(i,j));
            energy += (stress(i,j)-stressexact(i,j))*(dudx(i,j)-duexact(i,j));
        }
    }
    val[0] = energy;
    val[1] = L2;
    val[2] = H1;
}
