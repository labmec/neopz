/**
 * @file
 * @brief Contains implementations of the TPZElasticityMaterial methods.
 */

#include "TPZElasticity2DHybrid.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logdata("pz.material.elasticity.data");
#endif

#include <fstream>
using namespace std;

TPZElasticity2DHybrid::TPZElasticity2DHybrid() : TPZRegisterClassId(&TPZElasticity2DHybrid::ClassId),
TPZElasticityMaterial(0) {
}

TPZElasticity2DHybrid::TPZElasticity2DHybrid(int id) : TPZRegisterClassId(&TPZElasticity2DHybrid::ClassId),
TPZElasticityMaterial(id) {
}

TPZElasticity2DHybrid::TPZElasticity2DHybrid(int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress)
: TPZRegisterClassId(&TPZElasticity2DHybrid::ClassId), TPZElasticityMaterial(num,E,nu,fx,fy,plainstress) {
}

TPZElasticity2DHybrid::~TPZElasticity2DHybrid() {
}

void TPZElasticity2DHybrid::Contribute(TPZVec<TPZMaterialData> &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    TPZMaterialData::MShapeFunctionType shapetype = data[0].fShapeType;
    if(shapetype==data[0].EVecShape){
        DebugStop();
        return;
    }
    TPZFMatrix<REAL> &dphi = data[0].dphix;
    TPZFMatrix<REAL> &phi = data[0].phi;
    TPZFMatrix<REAL> &axes=data[0].axes;
    
    TPZFMatrix<REAL> &phirigidbodymode = data[1].phi;
    
    int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
    phc = phi.Cols();
    phr = phi.Rows();
    int phr_rigid = phirigidbodymode.Rows();
    dphc = dphi.Cols();
    dphr = dphi.Rows();
    efr = ef.Rows();
    efc = ef.Cols();
    ekr = ek.Rows();
    ekc = ek.Cols();
    if(phc != 1 || dphr != 2 || phr != dphc || 2*phr+2*phr_rigid != ekr){
        PZError << "\nTPZElasticityMaterial.contr, inconsistent input data : \n" <<
        "phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
        " phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
        dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
        << ek.Cols() <<
        "\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
        << ef.Cols() << "\n";
        DebugStop();
        return;
        //		PZError.show();
    }
    TPZManVector<STATE,3> locforce(ff);
    if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
        TPZManVector<STATE,3> res(3);
        fForcingFunction->Execute(data[0].x,res);
        locforce = res;
    }
    
    REAL E(fE_def), nu(fnu_def);
    
    if (fElasticity) {
        TPZManVector<STATE,2> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity->Execute(data[0].x, result, Dres);
        E = result[0];
        nu = result[1];
    }
    
    REAL Eover1MinNu2 = E/(1-nu*nu);
    REAL Eover21PlusNu = E/(2.*(1+nu));
    

    TPZFNMatrix<4,STATE> du(2,2);
    /*
     * Plane strain materials values
     */
    REAL nu1 = 1. - nu;//(1-nu)
    REAL nu2 = (1.-2.*nu)/2.;
    REAL F = E/((1.+nu)*(1.-2.*nu));
    STATE epsx, epsy,epsxy,epsz = 0.;
    TPZFNMatrix<9,STATE> DSolxy(2,2);
    // dudx - dudy
    DSolxy(0,0) = data[0].dsol[0](0,0)*axes(0,0)+data[0].dsol[0](1,0)*axes(1,0);
    DSolxy(1,0) = data[0].dsol[0](0,0)*axes(0,1)+data[0].dsol[0](1,0)*axes(1,1);
    // dvdx - dvdy
    DSolxy(0,1) = data[0].dsol[0](0,1)*axes(0,0)+data[0].dsol[0](1,1)*axes(1,0);
    DSolxy(1,1) = data[0].dsol[0](0,1)*axes(0,1)+data[0].dsol[0](1,1)*axes(1,1);
    epsx = DSolxy(0,0);// du/dx
    epsy = DSolxy(1,1);// dv/dy
    epsxy = 0.5*(DSolxy(1,0)+DSolxy(0,1));
    REAL lambda = GetLambda(E,nu);
    REAL mu = GetMU(E,nu);
    if (fPlaneStress) {
        epsz = -lambda*(epsx+epsy)/(lambda+2.*mu);
    }
    else
    {
        epsz = 0.;
    }
    //    epsz = data[1].sol[0][0];
    STATE SigX = lambda*(epsx+epsy+epsz)+2.*mu*epsx + fPreStressXX;
    STATE SigY = lambda*(epsx+epsy+epsz)+2.*mu*epsy + fPreStressYY;
    
    STATE SigZ = lambda*(epsx+epsy+epsz)+2.*mu*epsz;
    STATE TauXY = 2*mu*epsxy+fPreStressXY;
    
    
    for( int in = 0; in < phr; in++ ) {
        du(0,0) = dphi(0,in)*axes(0,0)+dphi(1,in)*axes(1,0);//dvx
        du(1,0) = dphi(0,in)*axes(0,1)+dphi(1,in)*axes(1,1);//dvy
        
        
        for (int col = 0; col < efc; col++)
        {
            ef(2*in,   col) += weight * (locforce[0]*phi(in,0) - du(0,0)*(SigX) - du(1,0)*(TauXY));  // direcao x
            ef(2*in+1, col) += weight * (locforce[1]*phi(in,0) - du(0,0)*(TauXY) - du(1,0)*(SigY));  // direcao y <<<----
        }
        for( int jn = 0; jn < phr; jn++ ) {
            du(0,1) = dphi(0,jn)*axes(0,0)+dphi(1,jn)*axes(1,0);//dux
            du(1,1) = dphi(0,jn)*axes(0,1)+dphi(1,jn)*axes(1,1);//duy
            
            
            if (fPlaneStress != 1){
                /* Plane Strain State */
                ek(2*in,2*jn) += weight * (
                                           nu1 * du(0,0)*du(0,1)+ nu2 * du(1,0)*du(1,1)
                                           ) * F;
                
                ek(2*in,2*jn+1) += weight * (
                                             nu*du(0,0)*du(1,1)+ nu2*du(1,0)*du(0,1)
                                             ) * F;
                
                ek(2*in+1,2*jn) += weight * (
                                             nu*du(1,0)*du(0,1)+ nu2*du(0,0)*du(1,1)
                                             ) * F;
                
                ek(2*in+1,2*jn+1) += weight * (
                                               nu1*du(1,0)*du(1,1)+ nu2*du(0,0)*du(0,1)
                                               ) * F;
            }
            else{
                DebugStop();
                /* Plain stress state */
                ek(2*in,2*jn) += weight * (
                                           Eover1MinNu2 * du(0,0)*du(0,1)+ Eover21PlusNu * du(1,0)*du(1,1)
                                           );
                
                ek(2*in,2*jn+1) += weight * (
                                             Eover1MinNu2*nu*du(0,0)*du(1,1)+ Eover21PlusNu*du(1,0)*du(0,1)
                                             );
                
                ek(2*in+1,2*jn) += weight * (
                                             Eover1MinNu2*nu*du(1,0)*du(0,1)+ Eover21PlusNu*du(0,0)*du(1,1)
                                             );
                
                ek(2*in+1,2*jn+1) += weight * (
                                               Eover1MinNu2*du(1,0)*du(1,1)+ Eover21PlusNu*du(0,0)*du(0,1)
                                               );
            }
        }
    }
    //constantes 1
    for (int in=0; in<phr; in++) {
        ek(2*phr,2*in) -= (STATE)phi(in,0)*weight;
        ek(2*phr+1,2*in+1) -= (STATE)phi(in,0)*weight;
        ek(2*in,2*phr) -= phi(in,0)*weight;
        ek(2*in+1,2*phr+1) -= phi(in,0)*weight;
        // u*(x-x0)-v*(y-y0)
        ek(2*phr+2,2*in) -= (phi(in,0)*phirigidbodymode(1))*weight;
        ek(2*phr+2,2*in+1) -= (-phi(in,0)*phirigidbodymode(0))*weight;
        ek(2*in,2*phr+2) -= (phi(in,0)*phirigidbodymode(1))*weight;
        ek(2*in+1,2*phr+2) -= (-phi(in,0)*phirigidbodymode(0))*weight;
    }
    
    //constante 2
    ek(2*phr+3,2*phr) += weight;
    ek(2*phr,2*phr+3) += weight;
    ek(2*phr+4,2*phr+1) += weight;
    ek(2*phr+1,2*phr+4) += weight;
    ek(2*phr+5,2*phr+2) += weight;
    ek(2*phr+2,2*phr+5) += weight;
    //ek(phr,phr) -= 1.*(STATE)weight;
    
}

//void TPZElasticity2DHybrid::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
//    TPZElasticityMaterial::Solution(datavec[0],var,Solout);
//}



void TPZElasticity2DHybrid::ContributeBC(TPZMaterialData &data,REAL weight,
										 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
    
    
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype==data.EVecShape){
        ContributeVecShapeBC(data,weight,ek, ef,bc);
        return;
    }
    
    TPZManVector<STATE,3> val2(2,0.);
    val2[0] = bc.Val2()(0,0);
    val2[1] = bc.Val2()(1,0);
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE,3> result(3,0.);
        TPZFNMatrix<9,STATE> resval1(2,2,0.);
        bc.ForcingFunction()->Execute(data.x, result, resval1);
        val2[0] = result[0];
        val2[1] = result[1];
    }
	TPZFMatrix<REAL> &phi = data.phi;
     int dim = Dimension();

	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    
	int phr = phi.Rows();
	short in,jn;
    
    if (ef.Cols() != bc.NumLoadCases()) {
        DebugStop();
    }
	
//		In general when the problem is  needed to stablish any convention for ContributeBC implementations

//     REAL v2[2];
// 	v2[0] = bc.Val2()(0,0);
// 	v2[1] = bc.Val2()(1,0);
    int nstate = NStateVariables();

    TPZFMatrix<STATE> &v1 = bc.Val1();


    switch (bc.Type()) {
        case 1 :			// Neumann condition
        {
            for(in = 0 ; in < phr; in++) {
                for (int il = 0; il<NumLoadCases(); il++)
                {
                    REAL v2[2];
                    v2[0] = val2[0];
                    v2[1] = val2[0];
                    ef(2*in,il)   += BIGNUMBER * v2[0] * phi(in,0) * weight;        // forced v2 displacement
                    ef(2*in+1,il) += BIGNUMBER * v2[1] * phi(in,0) * weight;        // forced v2 displacement
                }
                for (jn = 0 ; jn < phi.Rows(); jn++)
                {
                    ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
                }
            }
        }
            break;
            
        case 0 :		// Dirichlet condition
        {
            for (in = 0; in < phr; in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    ef(2*in,il) += -val2[0] * phi(in,0) * weight;        // force in x direction
                    ef(2*in+1,il) +=  -val2[1] * phi(in,0) * weight;      // force in y direction
                }
            }
        }
            break;
            
        case 2 :		// Mixed Condition
        {
            DebugStop();
            for(in = 0 ; in < phi.Rows(); in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(2*in,il) += v2(0,0) * phi(in,0) * weight;        // force in x direction
                    ef(2*in+1,il) += v2(1,0) * phi(in,0) * weight;      // forced in y direction
                }
#ifdef PZDEBUG
                if(bc.Val1()(0,1) != 0. || bc.Val1()(1,0) != 0.) DebugStop();
#endif
                for (jn = 0 ; jn < phi.Rows(); jn++) {
                    ek(2*in,2*jn) += 1./bc.Val1()(0,0) * phi(in,0) * phi(jn,0) * weight;         // peso de contorno => integral de contorno
                    ek(2*in+1,2*jn) += bc.Val1()(1,0) * phi(in,0) * phi(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += 1./bc.Val1()(1,1) * phi(in,0) * phi(jn,0) * weight;
                    ek(2*in,2*jn+1) += bc.Val1()(0,1) * phi(in,0) * phi(jn,0) * weight;
                }
            }   // este caso pode reproduzir o caso 0 quando o deslocamento
            
            break;
            
        }
    }
}

//void TPZElasticity2DHybrid::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
//                               TPZFMatrix<STATE> &ef, TPZBndCond &bc){
//    bc.Contribute(datavec[0],weight,ek,ef);
//    ContributeBC(datavec[0],weight,ek,ef,bc);
//}

TPZElasticity2DHybrid::TPZElasticity2DHybrid(const TPZElasticity2DHybrid &copy) :
TPZRegisterClassId(&TPZElasticity2DHybrid::ClassId), TPZElasticityMaterial(copy)
{
}


int TPZElasticity2DHybrid::ClassId() const{
    return Hash("TPZElasticity2DHybrid") ^ TPZElasticityMaterial::ClassId() << 1;
}

template class TPZRestoreClass<TPZElasticity2DHybrid>;

void TPZElasticity2DHybrid::Read(TPZStream &buf, void *context)
{
	TPZElasticityMaterial::Read(buf,context);
	
}

void TPZElasticity2DHybrid::Write(TPZStream &buf, int withclassid) const
{
	TPZElasticityMaterial::Write(buf,withclassid);
	
}

