//
//  TPZMonoPhaseWell.cpp
//  PZ
//
//  Created by omar duran on 25/05/2015.
//
//

#include "TPZMonoPhaseWell.h"
#include "pzmatrix.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"

/** @brief Constructor */
TPZMonoPhaseWell::TPZMonoPhaseWell(int id): TPZMaterial(id)
{
    
    /** @brief Well cross sectional area */
    fAp=0.0508*0.0508*M_PI;
    
    /** @brief Well diameter */
    fd=0.0508;
    
    /** @brief Pipe roughness */
    fepsilon=0.00006;
    
    /** @brief Pipe inclination */
    ftheta=90.0*M_PI/180.0;
    
    /** @brief atribute for derivative computations */
    fdelta=1.0e-10;
    
    /** @brief time step */
    fdt=1.0;
    
    /** @brief Gravity constant */
    fg = 9.81;
    
    fNextStep = false;
    
}

/** @brief Destructor */
TPZMonoPhaseWell::~TPZMonoPhaseWell()
{
    /** @brief Well cross sectional area */
    fAp=0.0508*0.0508*M_PI;
    
    /** @brief Well diameter */
    fd=0.0508;
    
    /** @brief Pipe roughness */
    fepsilon=0.00006;
    
    /** @brief Pipe inclination */
    ftheta=90.0*M_PI/180.0;
    
    /** @brief atribute for derivative computations */
    fdelta=1.0e-10;
    
    /** @brief time step */
    fdt=1.0;
    
    /** @brief Gravity constant */
    fg = 9.81;
    
}


/** @brief Copy constructor */
TPZMonoPhaseWell::TPZMonoPhaseWell(TPZMonoPhaseWell & copy){
    
    fAp=copy.fAp;
    fd=copy.fd;
    fepsilon=copy.fepsilon;
    ftheta=copy.ftheta;
    fdelta=copy.fdelta;
    fdt=copy.fdt;
    fg=copy.fg;
    
}


/**
 * @name Contribute methods (weak formulation)
 * @{
 */

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
 * @param datavec [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 */
void TPZMonoPhaseWell::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int wblock = 0;
    int pblock = 1;
    int c0 = 0;
    
    TPZFMatrix<REAL> wphis = datavec[wblock].phi;
    TPZFMatrix<REAL> pphis = datavec[pblock].phi;
    TPZFMatrix<REAL> dwphis = datavec[wblock].dphix;
    TPZFMatrix<REAL> dpphis = datavec[pblock].dphix;
    
    // n time step
    TPZManVector<STATE,3> wa0 = datavec[wblock].sol[c0];
    TPZManVector<STATE,3> pa0 = datavec[pblock].sol[c0];
    TPZFMatrix<STATE> dw0 = datavec[wblock].dsol[c0];
    TPZFMatrix<STATE> dp0 = datavec[pblock].dsol[c0];
    
//    REAL dwds0 = dw0(0,0)*datavec[wblock].axes(0,0)+dw0(0,0)*datavec[wblock].axes(0,1)+dw0(0,0)*datavec[wblock].axes(0,2);
    
    // Computing properties at (w,p) conditions
    REAL w0 = wa0[0];
    REAL p0 = pa0[0];
    
    REAL f0, dfdw0, dfdp0;
    REAL rho0, drhodw0, drhodp0;
    REAL mu0, dmudw0, dmudp0;
    this->friction(f0, p0, w0, dfdw0, dfdp0);
    this->Rho(rho0, p0, drhodw0, drhodp0);
    this->Mu(mu0, p0, dmudw0, dmudp0);
    
//    REAL ek0 = ((w0*w0)/(rho0*rho0)) * drhodp0;
//    REAL ekn0 = (2.0/rho0);
    
    
//    // n+1 time step
//    TPZManVector<REAL,3> wa1 = datavec[wblock].sol[c1];
//    TPZManVector<REAL,3> pa1 = datavec[pblock].sol[c1];
//    TPZFMatrix<REAL> dw1 = datavec[wblock].dsol[c1];
//    TPZFMatrix<REAL> dp1 = datavec[pblock].dsol[c1];
//    //REAL dwds1 = dw1(0,0)*datavec[wblock].axes(0,0)+dw1(1,0)*datavec[wblock].axes(1,0)+dw1(2,0)*datavec[wblock].axes(2,0);
//    
//    // Computing properties at (w,p) conditions
//    REAL w1 = wa1[0];
//    REAL p1 = pa1[0];
//    
//    REAL f1, dfdw1, dfdp1;
//    REAL rho1, drhodw1, drhodp1;
//    this->friction(f1, p1, w1, dfdw1, dfdp1);
//    this->Rho(rho1, p1, drhodw1, drhodp1);
    
    
    int nwphis = wphis.Rows();
    int npphis = pphis.Rows();
    
    if (fNextStep) {
        
        // Conservation of linear momentum like HagenPoiseille
        for(int iw = 0; iw < nwphis; iw++ )
        {
            REAL dvdsi = dwphis(0,iw)*datavec[wblock].axes(0,0)+dwphis(0,iw)*datavec[wblock].axes(0,1)+dwphis(0,iw)*datavec[wblock].axes(0,2);
            
            for(int jw = 0; jw < nwphis; jw++ )
            {
                ek(iw,jw) += weight * ( 0.0/(fdt) * wphis(jw,0) * wphis(iw,0)
                                        + wphis(jw,0) * ( (64 * mu0) / (8.0*fd)) * (M_PI*fd/fAp) * wphis(iw,0));
            }
            
            for(int jp = 0; jp < npphis; jp++ )
            {
                ek(iw,jp + nwphis) += weight * ( - dvdsi * pphis(jp,0) + drhodp0 * fg * sin(ftheta) * pphis(jp,0) );
            }
            
        }
        
        // Conservation of mass
        for(int ip = 0; ip < npphis; ip++ )
        {
            
            for(int jw = 0; jw < nwphis; jw++ )
            {
                REAL dvdsj = dwphis(0,jw)*datavec[wblock].axes(0,0)+dwphis(0,jw)*datavec[wblock].axes(0,1)+dwphis(0,jw)*datavec[wblock].axes(0,2);
                
                ek(ip + nwphis,jw) += weight * ( -1.0 * dvdsj * pphis(ip,0) );
            }
            
            for(int jp = 0; jp < npphis; jp++ )
            {
                ek(ip + nwphis, jp + nwphis) += weight * ( ((-1.0/fdt) * (drhodp0 * pphis(jp,0) /* - rho0 */ ) )* pphis(ip,0) );
            }
        }
        
    }
    
//    // Conservation of linear momentum like HagenPoiseille
//    for(int iw = 0; iw < nwphis; iw++ )
//    {
//        ef(iw,0) += weight * (-1.0/(fdt) * (/* w1 */ w0) * wphis(iw,0));
//    }
//    
//    // Conservation of mass
//    for(int ip = 0; ip < npphis; ip++ )
//    {
//        ef(ip + nwphis,0) += weight * ( ((-1.0/fdt) * (/* rho1 */ - rho0 ) + dwds0 )* pphis(ip,0) );
//    }
    
    this->Contribute(datavec, weight, ef);
}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
 * @param datavec [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 */
void TPZMonoPhaseWell::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    int wblock = 0;
    int pblock = 1;
    int c0 = 0;
    
    
    TPZFMatrix<REAL> wphis = datavec[wblock].phi;
    TPZFMatrix<REAL> pphis = datavec[pblock].phi;
    TPZFMatrix<REAL> dwphis = datavec[wblock].dphix;
    TPZFMatrix<REAL> dpphis = datavec[pblock].dphix;
    
    
    // n time step
    TPZManVector<STATE,3> wa0 = datavec[wblock].sol[c0];
    TPZManVector<STATE,3> pa0 = datavec[pblock].sol[c0];
    TPZFMatrix<STATE> dw0 = datavec[wblock].dsol[c0];
    TPZFMatrix<STATE> dp0 = datavec[pblock].dsol[c0];
    REAL dwds0 = dw0(0,0)*datavec[wblock].axes(0,0)+dw0(0,0)*datavec[wblock].axes(0,1)+dw0(0,0)*datavec[wblock].axes(0,2);
    
    // Computing properties at (w,p) conditions
    REAL w0 = wa0[0];
    REAL p0 = pa0[0];
    
    REAL f0, dfdw0, dfdp0;
    REAL rho0, drhodw0, drhodp0;
    REAL mu0, dmudw0, dmudp0;
    this->friction(f0, p0, w0, dfdw0, dfdp0);
    this->Rho(rho0, p0, drhodw0, drhodp0);
    this->Mu(mu0, p0, dmudw0, dmudp0);
    
//    REAL ek0 = ((w0*w0)/(rho0*rho0)) * drhodp0;
//    REAL ekn0 = (2.0/rho0);
    

//    // n+1 time step
//    TPZManVector<REAL,3> wa1 = datavec[wblock].sol[c1];
//    TPZManVector<REAL,3> pa1 = datavec[pblock].sol[c1];
//    TPZFMatrix<REAL> dw1 = datavec[wblock].dsol[c1];
//    TPZFMatrix<REAL> dp1 = datavec[pblock].dsol[c1];
//    REAL dwds1 = dw1(0,0)*datavec[wblock].axes(0,0)+dw1(1,0)*datavec[wblock].axes(1,0)+dw1(2,0)*datavec[wblock].axes(2,0);
//    
//    // Computing properties at (w,p) conditions
//    REAL w1 = wa1[0];
//    REAL p1 = pa1[0];
//    
//    REAL f1, dfdw1, dfdp1;
//    REAL rho1, drhodw1, drhodp1;
//    this->friction(f1, p1, w1, dfdw1, dfdp1);
//    this->Rho(rho1, p1, drhodw1, drhodp1);
    
    
    int nwphis = wphis.Rows();
    int npphis = pphis.Rows();
    
    if (fNextStep) {
        
        // Conservation of linear momentum like HagenPoiseille
        for(int iw = 0; iw < nwphis; iw++ )
        {
            REAL dvdsi = dwphis(0,iw)*datavec[wblock].axes(0,0)+dwphis(0,iw)*datavec[wblock].axes(0,1)+dwphis(0,iw)*datavec[wblock].axes(0,2);
            
            ef(iw,0) += weight * ( (1.0/(fdt) * (w0) * wphis(iw,0)
                                    + (w0 * ( (64 * mu0) / (8.0*fd)) * (M_PI*fd/fAp)) + rho0 * fg * sin(ftheta)) * wphis(iw,0)
                                    - dvdsi * p0 );
        }
        
        // Conservation of mass
        for(int ip = 0; ip < npphis; ip++ )
        {
            ef(ip + nwphis,0) += weight * ( ((-1.0/fdt) * (rho0 /* - rho0 */ ) - dwds0 )* pphis(ip,0) );
        }
        
        
        return;
    }

    // Conservation of linear momentum like HagenPoiseille
    for(int iw = 0; iw < nwphis; iw++ )
    {
        ef(iw,0) += weight * (-1.0/(fdt) * (/* w1 */ w0) * wphis(iw,0));
    }

    TPZManVector<STATE,1> fvalue(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[pblock].x,fvalue);
    }
    
    // Conservation of mass
    for(int ip = 0; ip < npphis; ip++ )
    {
        ef(ip + nwphis,0) += weight * ( ( (-1.0/fdt) * (/* rho1 */ - rho0 ) - dwds0 + fvalue[0] )* pphis(ip,0) );
    }



}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
 * to multiphysics simulation.
 * @param datavec [in]  stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition material
 * @since October 18, 2011
 */
void TPZMonoPhaseWell::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    if (!fNextStep) {
        
        return;
    }
    
    int wblock = 0;
    int pblock = 1;
    int c0 = 0;
    int c1 = 0;
    
    TPZFMatrix<REAL> wphis = datavec[wblock].phi;
    TPZFMatrix<REAL> pphis = datavec[pblock].phi;
    TPZFMatrix<REAL> dwphis = datavec[wblock].dphix;
    TPZFMatrix<REAL> dpphis = datavec[pblock].dphix;
    
    // n time step
    TPZManVector<STATE,3> wa0 = datavec[wblock].sol[c0];
    TPZManVector<STATE,3> pa0 = datavec[pblock].sol[c0];
    TPZFMatrix<STATE> dw0 = datavec[wblock].dsol[c0];
    TPZFMatrix<STATE> dp0 = datavec[pblock].dsol[c0];
//    REAL dwds0 = dw0(0,0)*datavec[wblock].axes(0,0)+dw0(1,0)*datavec[wblock].axes(1,0)+dw0(2,0)*datavec[wblock].axes(2,0);
    
    // Computing properties at (w,p) conditions
    REAL w0 = wa0[0];
    REAL p0 = pa0[0];
    
    REAL f0, dfdw0, dfdp0;
    REAL rho0, drhodw0, drhodp0;
    this->friction(f0, p0, w0, dfdw0, dfdp0);
    this->Rho(rho0, p0, drhodw0, drhodp0);
    
//    REAL ek0 = ((w0*w0)/(rho0*rho0)) * drhodp0;
//    REAL ekn0 = (2.0/rho0);
    
    // n+1 time step
    TPZManVector<STATE,3> wa1 = datavec[wblock].sol[c1];
    TPZManVector<STATE,3> pa1 = datavec[pblock].sol[c1];
    TPZFMatrix<STATE> dw1 = datavec[wblock].dsol[c1];
    TPZFMatrix<STATE> dp1 = datavec[pblock].dsol[c1];
    //REAL dwds1 = dw1(0,0)*datavec[wblock].axes(0,0)+dw1(1,0)*datavec[wblock].axes(1,0)+dw1(2,0)*datavec[wblock].axes(2,0);
    
    // Computing properties at (w,p) conditions
    REAL w1 = wa1[0];
    REAL p1 = pa1[0];
    
    REAL f1, dfdw1, dfdp1;
    REAL rho1, drhodw1, drhodp1;
    this->friction(f1, p1, w1, dfdw1, dfdp1);
    this->Rho(rho1, p1, drhodw1, drhodp1);
    
    
    int nwphis = wphis.Rows();
    int npphis = pphis.Rows();
    
    REAL value = bc.Val2()(0,0);

    switch (bc.Type())
    {
        case 0 :            // Dirichlet condition
        {
            for(int iw=0; iw<nwphis; iw++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iw,0) += weight * value * wphis(iw,0);
            }
        }
        break;
            
        case 1 :            // Neumann condition
        {
            
            for(int iw=0; iw<nwphis; iw++)
            {
                ef(iw,0) += gBigNumber * (w0 - value) * wphis(iw,0);
                for (int jw=0; jw < nwphis; jw++) {
                    
                    ek(iw,jw) += weight * gBigNumber * wphis(iw,0) * wphis(jw,0);
                }
            }
            
            for(int iw=0; iw<nwphis; iw++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iw,0) += weight * p0 * wphis(iw,0);
                for(int jp=0; jp< npphis; jp++)
                {
                    ek(iw, jp + nwphis) += weight * pphis(jp,0) * wphis(iw,0);
                }
            }
            
        }
        break;
    }

    
}

/** @} */

/**
 * @brief Fill material data parameter with necessary requirements for the
 * @since April 10, 2007
 */
/**
 * Contribute method. Here, in base class, all requirements are considered as necessary.
 * Each derived class may optimize performance by selecting only the necessary data.
 */
void TPZMonoPhaseWell::FillDataRequirements(TPZMaterialData &data){
    return;
}

/**
 * @brief Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered as necessary.
 * Each derived class may optimize performance by selecting only the necessary data.
 */
void TPZMonoPhaseWell::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++ )
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
    }
    
}

/**
 * @brief Fill material data parameter with necessary requirements for the
 * ContributeBC method. Here, in base class, all requirements are considered as necessary.
 * Each derived class may optimize performance by selecting only the necessary data.
 */
void TPZMonoPhaseWell::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsNeighborSol = true;
    }
}




/** @brief Print out the data associated with the material */
void TPZMonoPhaseWell::Print(std::ostream &out){
    
}

/** @brief Returns the variable index associated with the name */
int TPZMonoPhaseWell::VariableIndex(const std::string &name){
    if(!strcmp("Momentum",name.c_str()))        return  1;
    if(!strcmp("Pressure",name.c_str()))        return  2;
    if(!strcmp("Velocity",name.c_str()))        return  3;
    if(!strcmp("Density",name.c_str()))        return  4;
    
    return TPZMaterial::VariableIndex(name);
    
}

int TPZMonoPhaseWell::NSolutionVariables(int var){
    if(var == 1) return 1;
    if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPZMonoPhaseWell::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    int wblock = 0;
    int pblock = 1;
    int c0 = 0;
    
    Solout.Resize(this->NSolutionVariables(var));
    
    
    // at each time step
    TPZManVector<STATE,3> wa0 = datavec[wblock].sol[c0];
    TPZManVector<STATE,3> pa0 = datavec[pblock].sol[c0];
//    TPZFMatrix<REAL> dw0 = datavec[wblock].dsol[c0];
//    TPZFMatrix<REAL> dp0 = datavec[pblock].dsol[c0];
//    REAL dwds0 = dw0(0,0)*datavec[wblock].axes(0,0)+dw0(0,0)*datavec[wblock].axes(0,1)+dw0(0,0)*datavec[wblock].axes(0,2);
    
    // Computing properties at (w,p) conditions
    REAL w0 = wa0[0];
    REAL p0 = pa0[0];
    
    REAL f0, dfdw0, dfdp0;
    REAL rho0, drhodw0, drhodp0;
    REAL mu0, dmudw0, dmudp0;
    this->friction(f0, p0, w0, dfdw0, dfdp0);
    this->Rho(rho0, p0, drhodw0, drhodp0);
    this->Mu(mu0, p0, dmudw0, dmudp0);
    
    
    switch (var) {
        case 1:
        {
            Solout[0] = w0;
        }
            break;
        case 2:
        {
            Solout[0] = p0;
        }
            break;
        case 3:
        {
            Solout[0] = w0/rho0;
        }
            break;
        case 4:
        {
            Solout[0] = rho0;
        }
            break;
        default:
            break;
    }
    
    
}




/** @brief friction factor using Haaland friction factor because it is a explicit expression of f */
void TPZMonoPhaseWell::friction(REAL &f, REAL P, REAL w, REAL &dfdw, REAL &dfdP){
    
    REAL mu, dmudw, dmudP;
    REAL Re,Rem,Rep,dRedw;
    this->Mu(mu,P,dmudw,dmudP);
    
    Re = w * fd / mu ;
    
//    if (Re >= 2000.0) {
//        std::cout << " Reynolds number is turbulent =  "  << Re << std::endl;
//    }
//    
//    if (fabs(Re) <= 1.0e-12) {
//        std::cout << " Reynolds number is zero =  "  << Re << std::endl;
//        f = 0.0;
//        dfdw = 0.0;
//        return;
//    }
    
    if (Re <= 2000.0) {
    // Using Laminar friction factor
        Rem= (w - fdelta) * fd / mu;
        Rep= (w + fdelta) * fd / mu;
        dRedw = (Rep - Rem)/(2.0*fdelta);
        f = 64.0/Re;
        dfdw = -(64.0*dRedw)/(Re*Re);
    }else
    {
    // Using turbulent Haaland friction factor
        Rem= (w - fdelta) * fd / mu;
        Rep= (w + fdelta) * fd / mu;
        dRedw = (Rep - Rem)/(2.0*fdelta);
        f = 64.0/Re;
        dfdw = -(64.0*dRedw)/(Re*Re);
    }

    return;
}


/** @brief Fluid density  */
void TPZMonoPhaseWell::Rho(REAL &rho, REAL P, REAL &drhodw, REAL &drhodP){
    
    REAL c      = 10.0e-10;
    REAL rhoref = 1000.0;
    REAL Pref = 1.0e6;
    rho = rhoref * (1.0 + c*(P-Pref));
    drhodw = 0.0;
    drhodP = rhoref*c;
    
}

/** @brief Fluid density  */
void TPZMonoPhaseWell::Mu(REAL &mu, REAL P, REAL &dmudw, REAL &dmudP){
    mu = 0.001;
    dmudw = 0.0;
    dmudP = 0.0;
}

/** @{
 * @name Save and Load methods
 */

/** @brief Unique identifier for serialization purposes */
int TPZMonoPhaseWell::ClassId() const{
    return Hash("TPZMonoPhaseWell") ^ TPZMaterial::ClassId() << 1;
}

/** @brief Saves the element data to a stream */
void TPZMonoPhaseWell::Write(TPZStream &buf, int withclassid) const{
    
}

/** @brief Reads the element data from a stream */
void TPZMonoPhaseWell::Read(TPZStream &buf, void *context){
    
}

/** @} */
