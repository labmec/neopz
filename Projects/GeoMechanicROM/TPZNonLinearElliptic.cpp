//
//  TPZNonLinearElliptic.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZNonLinearElliptic.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.permeabilityc"));
#endif


TPZNonLinearElliptic::TPZNonLinearElliptic():TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin>() {

    fDim = 2;
    fnstate = 1;
    fmu_0 = 1.0;
    fmu_1 = 1.0;
    fmu_2 = 1.0;
    
}

TPZNonLinearElliptic::TPZNonLinearElliptic(int matid, int dim): TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> (matid)   {

    fDim = dim;
    fnstate = 1;    
    fmu_0 = 1.0;
    fmu_1 = 1.0;
    fmu_2 = 1.0;
    
}

TPZNonLinearElliptic::~TPZNonLinearElliptic(){
}


int TPZNonLinearElliptic::NStateVariables() {
    return fnstate;
}


void TPZNonLinearElliptic::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef){
    
    int u_b = 0;
    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    
    TPZFMatrix<REAL>    &dphiu   =   datavec[u_b].dphix;
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    
    
    TPZFNMatrix<300,REAL> Grad_phiu_xy;
    TPZFNMatrix<9,REAL> axes_u_T, Gradu_xy;
    
    axes_u.Transpose(&axes_u_T);
    axes_u.Multiply(dphiu,Grad_phiu_xy,1/* Transpose axes_u */);
    axes_u.Multiply(du,Gradu_xy,1/* Transpose axes_u */);
    
    int nphi_u = phiu.Rows();
    int first_u = 0;
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        return;
    }
    
    TPZManVector<STATE,1> f(1,0.0);
    
    if(HasfTimedependentForcingFunction())
    {
        TPZFMatrix<double> gradf;
        REAL time = 0.;
        TimeDependentForcingFunction()->Execute(datavec[u_b].x, time, f, gradf);
    }
    
    for (int iu = 0; iu < nphi_u; iu++) {
        
        REAL dot = 0.0;
        for (int i = 0; i < fDim; i++) {
            dot += Gradu_xy(i,0)*Grad_phiu_xy(i,iu);
        }
        
//        ef(iu + first_u, 0)   += weight * ( fmu_0 * dot + (fmu_1 * ((exp(fmu_2*u[0])-1.0)/fmu_2) - f[0]) * phiu(iu, 0) );
//        ef(iu + first_u, 0)   += weight * ( fmu_0 * dot + (fmu_1 * ((fmu_2*u[0]*u[0]-1.0)/fmu_2) - f[0]) * phiu(iu, 0) );
        ef(iu + first_u, 0)   += weight * ( fmu_0 * dot + (0.0*fmu_1 * u[0] - f[0]/fmu_2) * phiu(iu, 0) );// Linear case
        
        for (int ju = 0; ju < nphi_u; ju++) {
            
            REAL dot = 0.0;
            for (int i = 0; i < fDim; i++) {
                dot += Grad_phiu_xy(i,ju)*Grad_phiu_xy(i,iu);
            }
            
            ek(iu + first_u, ju + first_u)   += weight * (fmu_0 * dot + fmu_1 * exp(fmu_2*u[0]) * phiu(ju, 0) * phiu(iu, 0));
//            ek(iu + first_u, ju + first_u)   += weight * (fmu_0 * dot +  2.0 * fmu_1 * u[0] * phiu(ju, 0) * phiu(iu, 0));
            ek(iu + first_u, ju + first_u)   += weight * (fmu_0 * dot + 0.0*fmu_1 * phiu(ju, 0) * phiu(iu, 0));// Linear case
            
        }
        
    }
    
}

void TPZNonLinearElliptic::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    
}

void TPZNonLinearElliptic::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int u_b = 0;
    
    TPZFMatrix<REAL>  &phiu = datavec[u_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    
    int phru = phiu.Rows();
    short in,jn;
    REAL v[3];
    v[0] = bc.Val2()(0,0);	//	u or flux
    
    REAL time = this->SimulationData()->t();
    REAL Value = bc.Val2()(0,0);
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,3> f(1);
        TPZFMatrix<REAL> gradf;
        bc.TimedependentBCForcingFunction()->Execute(datavec[u_b].x, time, f, gradf);
        v[0] = f[0];	//	u or flux
    }
    else{
        Value = bc.Val2()(0,0);
    }
    
    // Dirichlet in Pressure
    switch (bc.Type())
    {
        case 0 :
        {
            for(in = 0 ; in < phru; in++)
            {
                ef(in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// u
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ek(in,jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
                }
            }
            
            break;
        }
            
        case 1 :
        {
            for(in = 0 ; in <phru; in++)
            {
                ef(in,0)		+= -1.0*v[0]*phiu(in,0)*weight;		//	flux
            }
            
            break;
        }
            
        default:
        {
            DebugStop();
        }
            break;
    }
    
    
}

void TPZNonLinearElliptic::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)

{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNeighborSol = true;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = true;
    }
}

void TPZNonLinearElliptic::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;
    }
}

void TPZNonLinearElliptic::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Properties for TPZNonLinearElliptic: \n";
    out << "\t Parameter 1   = "									<< fmu_1		<< std::endl;
    out << "\t Parameter 2   = "									<< fmu_2		<< std::endl;
    out << "Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
    
}

/** Returns the variable index associated with the name */
int TPZNonLinearElliptic::VariableIndex(const std::string &name)
{
    //	Elasticity Variables
    if(!strcmp("u",name.c_str()))				return	1;
    if(!strcmp("sigma",name.c_str()))             return	2;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZNonLinearElliptic::NSolutionVariables(int var){
    if(var == 1)	return 1;
    if(var == 2)	return fDim;
    
    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZNonLinearElliptic::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    int u_b = 0;

    // Getting the space functions
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_u(1,3,0.0);
    
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    
    //	u
    if(var == 1){
        Solout[0] = u[0];
        return;
    }

    //	sigma
    if(var == 2) {
        Solout[0] = Grad_u(0,0);
        Solout[1] = Grad_u(0,1);
        return;
    }
    
}
