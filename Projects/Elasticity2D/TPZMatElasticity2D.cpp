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


TPZMatElasticity2D::TPZMatElasticity2D():TPZMaterial(), ff(0), fnu(0.), fPlaneStress(0) {
    fE = 0.;
    fmatId = 0;
    ff.resize(2);
    ff[0]=0.;
    ff[1]=0.;
    fPlaneStress = 1.;
    
}

TPZMatElasticity2D::TPZMatElasticity2D(int matid):TPZMaterial(matid), ff(0), fnu(0.),fPlaneStress(0) {
    fE = 0.;
    ff.resize(2);
    ff[0]=0.;
    ff[1]=0.;
    fPlaneStress = 1;
    fmatId = matid;
}

TPZMatElasticity2D::~TPZMatElasticity2D(){
}


int TPZMatElasticity2D::NStateVariables() {
    return 2;
}


void TPZMatElasticity2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef) {
    
   
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiU     =  data.phi;
    TPZFMatrix<REAL> &dphiU     =  data.dphix;
    int phrU = phiU.Rows();
    
    int FirstU  = 0;

    TPZManVector<REAL,3> sol_u =    data.sol[0];
    
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    
    REAL LambdaL, MuL;
    
    // Functions computed at point x_{k} for each integration point
    LambdaL     = flambda;
    MuL         = fmu;
    
        //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
        //  Contribution of domain integrals for Jacobian matrix
        //  Elasticity Block (Equation for elasticity )
        //	Elastic equation
        //	Linear strain operator
        //	Ke Matrix
        TPZFMatrix<REAL>	du(2,2);
        for(int iu = 0; iu < phrU; iu++ )
        {
            //	Derivative for Vx
            du(0,0) = dphiU(0,iu)*data.axes(0,0)+dphiU(1,iu)*data.axes(1,0);
            //	Derivative for Vy
            du(1,0) = dphiU(0,iu)*data.axes(0,1)+dphiU(1,iu)*data.axes(1,1);
            
            for(int ju = 0; ju < phrU; ju++)
            {
                //	Derivative for Ux
                du(0,1) = dphiU(0,ju)*data.axes(0,0)+dphiU(1,ju)*data.axes(1,0);
                //	Derivative for Uy
                du(1,1) = dphiU(0,ju)*data.axes(0,1)+dphiU(1,ju)*data.axes(1,1);
                
                if (this->fPlaneStress == 1)
                {
                    /* Plain stress state */
                    ek(2*iu + FirstU, 2*ju + FirstU)	     += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(0,0)*du(0,1)		+ (2*MuL)*du(1,0)*du(1,1));
                    
                    ek(2*iu + FirstU, 2*ju+1 + FirstU)       += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(0,0)*du(1,1)			+ (2*MuL)*du(1,0)*du(0,1));
                    
                    ek(2*iu+1 + FirstU, 2*ju + FirstU)       += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(1,0)*du(0,1)			+ (2*MuL)*du(0,0)*du(1,1));
                    
                    ek(2*iu+1 + FirstU, 2*ju+1 + FirstU)     += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(1,0)*du(1,1)		+ (2*MuL)*du(0,0)*du(0,1));
                }
                else
                {
                    /* Plain Strain State */
                    ek(2*iu + FirstU,2*ju + FirstU)         += weight*	((LambdaL + 2*MuL)*du(0,0)*du(0,1)	+ (MuL)*du(1,0)*du(1,1));
                    
                    ek(2*iu + FirstU,2*ju+1 + FirstU)       += weight*	(LambdaL*du(0,0)*du(1,1)			+ (MuL)*du(1,0)*du(0,1));
                    
                    ek(2*iu+1 + FirstU,2*ju + FirstU)       += weight*	(LambdaL*du(1,0)*du(0,1)			+ (MuL)*du(0,0)*du(1,1));
                    
                    ek(2*iu+1 + FirstU,2*ju+1 + FirstU)     += weight*	((LambdaL + 2*MuL)*du(1,0)*du(1,1)	+ (MuL)*du(0,0)*du(0,1));
                    
                }
            }
        }
        //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    this->Contribute(data,weight,ef);
}

void TPZMatElasticity2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) {
    
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiU =  data.phi;
    
    TPZFMatrix<REAL> &dphiU = data.dphix;
    int phrU = phiU.Rows();
    int FirstU  = 0;
    
    TPZManVector<REAL,3> sol_u =data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    REAL LambdaL, MuL;
    
    LambdaL = flambda;
    MuL     = fmu;
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  Contribution of domain integrals for Residual Vector
    //  Elastic equation
    //  Linear strain operator
    //  Ke Matrix
    TPZFMatrix<REAL>    du(2,2);
    TPZFMatrix<REAL>    dux(2,2);
    TPZFMatrix<REAL>    duy(2,2);
    // Required check out of this implementation
    //  Derivative for Ux
    dux(0,1) = dsol_u(0,0)*data.axes(0,0)+dsol_u(1,0)*data.axes(1,0); // dUx/dx
    dux(1,1) = dsol_u(0,0)*data.axes(0,1)+dsol_u(1,0)*data.axes(1,1); // dUx/dy
    
    //  Derivative for Uy
    duy(0,1) = dsol_u(0,1)*data.axes(0,0)+dsol_u(1,1)*data.axes(1,0); // dUy/dx
    duy(1,1) = dsol_u(0,1)*data.axes(0,1)+dsol_u(1,1)*data.axes(1,1); // dUy/dy
    
    for(int iu = 0; iu < phrU; iu++ )
    {
        //  Derivative for Vx
        du(0,0) = dphiU(0,iu)*data.axes(0,0)+dphiU(1,iu)*data.axes(1,0);
        //  Derivative for Vy
        du(1,0) = dphiU(0,iu)*data.axes(0,1)+dphiU(1,iu)*data.axes(1,1);
        
//          Vector Force right hand term
             ef(2*iu + FirstU)     +=    weight*ff[0]*phiU(iu, 0);
             ef(2*iu+1 + FirstU)   +=    weight*ff[1]*phiU(iu, 0);
        
        if (fPlaneStress == 1)
        {
            /* Plain stress state */
            ef(2*iu + FirstU)           += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(0,0)*dux(0,1)      + (2*MuL)*du(1,0)*dux(1,1));
            
            ef(2*iu + FirstU)           += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(0,0)*duy(1,1)         + (2*MuL)*du(1,0)*duy(0,1));
            
            ef(2*iu+1 + FirstU)         += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(1,0)*dux(0,1)         + (2*MuL)*du(0,0)*dux(1,1));
            
            ef(2*iu+1 + FirstU)         += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(1,0)*duy(1,1)     + (2*MuL)*du(0,0)*duy(0,1));
        }
        else
        {
            /* Plain Strain State */
            ef(2*iu + FirstU)           += weight*  ((LambdaL + 2*MuL)*du(0,0)*dux(0,1)  + (MuL)*du(1,0)*(dux(1,1)));
            
            ef(2*iu + FirstU)           += weight*  (LambdaL*du(0,0)*duy(1,1)            + (MuL)*du(1,0)*(duy(0,1)));
            
            ef(2*iu+1 + FirstU)         += weight*  (LambdaL*du(1,0)*dux(0,1)            + (MuL)*du(0,0)*(dux(1,1)));
            
            ef(2*iu+1 + FirstU)         += weight*  ((LambdaL + 2*MuL)*du(1,0)*duy(1,1)  + (MuL)*du(0,0)*(duy(0,1)));
        }
    }
    
    REAL Biotalpha=1.0;
    REAL Pressure = 1.0;
    
    //  Matrix Qc
    //  Coupling matrix for Elastic equation
    for(int in = 0; in < phrU; in++ )
    {
        du(0,0) = dphiU(0,in)*data.axes(0,0)+dphiU(1,in)*data.axes(1,0);
        du(1,0) = dphiU(0,in)*data.axes(0,1)+dphiU(1,in)*data.axes(1,1);
        
        ef(2*in + FirstU)    += (-1.0)*Biotalpha*weight*(Pressure*du(0,0));
        ef(2*in+1 + FirstU)  += (-1.0)*Biotalpha*weight*(Pressure*du(1,0));
    }
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    
}

void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    if (bc.Val2().Rows() != 2)
    {
        std::cout << " Error:: This material need boundary conditions for ux, uy and p (pore pressure).\n";
        std::cout << " give me one matrix with this form Val2(3,1).\n";
        DebugStop();
    }
    
    if (bc.Val1().Rows() != 2)
    {
        std::cout << " Error:: This material need boundary conditions for ux, uy and p (pore pressure).\n";
        DebugStop();
    }
    
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<REAL,3> sol_u = data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    REAL uy = sol_u[1];
    
    int phru = phiu.Rows();
    short in,jn;
    STATE v2[3];
    v2[0] = bc.Val2()(0,0);	//	Ux displacement or Tnx
    v2[1] = bc.Val2()(1,0);	//	Uy displacement or Tny
    
    //	Here each digit represent an individual boundary condition corresponding to each state variable.
    //	0 means Dirichlet condition on x-y
    //	1 means Neumann condition
    //	2 means Dirichlet condition on x
    //	3 means Dirichlet condition on y
    
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
                ef(2*in,0)		+= BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
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
                ef(2*in,0)		+= v2[0]*phiu(in,0)*weight;		//	Tnx
                ef(2*in+1,0)	+= v2[1]*phiu(in,0)*weight;		//	Tny
            }
            break;
        }
        case 2 :
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
        case 3 :
        {
            //	Dirichlet condition for each state variable
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
            std::cout << "Bc not impemented!" << std::endl;
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
}


void TPZMatElasticity2D::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Plane Problem (fPlaneStress = 0, for Plane Strain conditions) " << fPlaneStress << std::endl;
    out << "Properties for Poroelasticity: \n";
    out << "\t Young modulus   = "											<< fE		<< std::endl;
    out << "\t Poisson Ratio   = "											<< fnu		<< std::endl;
    out << "\t First Lamé Parameter   = "									<< flambda	<< std::endl;
    out << "\t Second Lamé Parameter   = "									<< fmu		<< std::endl;
    out << "\t Body force vector B {X-direction, Y-direction}   = "			<< ff[0] << ' ' << ff[1]   << std::endl;
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
    
    
//    if(!strcmp("EStrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainVol;
//    if(!strcmp("EStrainXX",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXX;
//    if(!strcmp("EStrainYY",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYY;
//    if(!strcmp("EStrainZZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainZZ;
//    if(!strcmp("EStrainXY",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXY;
//    if(!strcmp("EStrainXZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXZ;
//    if(!strcmp("EStrainYZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYZ;
//    if(!strcmp("EElStrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainVol;
//    if(!strcmp("EElStrainXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXX;
//    if(!strcmp("EElStrainYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYY;
//    if(!strcmp("EElStrainZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainZZ;
//    if(!strcmp("EElStrainXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXY;
//    if(!strcmp("EElStrainXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXZ;
//    if(!strcmp("EElStrainYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYZ;
//    if(!strcmp("EPlStrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainVol;
//    if(!strcmp("EPlStrainXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXX;
//    if(!strcmp("EPlStrainYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYY;
//    if(!strcmp("EPlStrainZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainZZ;
//    if(!strcmp("EPlStrainXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXY;
//    if(!strcmp("EPlStrainXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXZ;
//    if(!strcmp("EPlStrainYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYZ;
//    if(!strcmp("EPlStrainSqJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2;
//    if(!strcmp("EPlStrainSqJ2El",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2El;
//    if(!strcmp("EPlAlpha",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlAlpha;
//    if(!strcmp("EDisplacementX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementX;
//    if(!strcmp("EDisplacementY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementY;
//    if(!strcmp("EDisplacementZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementZ;
//    if(!strcmp("EDisplacementTotal",	name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementTotal;
//    if(!strcmp("ETotStressI1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressI1;
//    if(!strcmp("ETotStressJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressJ2;
//    if(!strcmp("ETotStressXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXX;
//    if(!strcmp("ETotStressYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYY;
//    if(!strcmp("ETotStressZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressZZ;
//    if(!strcmp("ETotStressXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXY;
//    if(!strcmp("ETotStressXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXZ;
//    if(!strcmp("ETotStressYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYZ;
//    if(!strcmp("ETotStress1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress1;
//    if(!strcmp("ETotStress2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress2;
//    if(!strcmp("ETotStress3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress3;
//    if(!strcmp("EEffStressI1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressI1;
//    if(!strcmp("EEffStressJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressJ2;
//    if(!strcmp("EEffStressXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXX;
//    if(!strcmp("EEffStressYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYY;
//    if(!strcmp("EEffStressZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressZZ;
//    if(!strcmp("EEffStressXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXY;
//    if(!strcmp("EEffStressXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXZ;
//    if(!strcmp("EEffStressYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYZ;
//    if(!strcmp("EEffStress1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress1;
//    if(!strcmp("EEffStress2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress2;
//    if(!strcmp("EEffStress3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress3;
//    if(!strcmp("EYieldSurface1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface1;
//    if(!strcmp("EYieldSurface2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface2;
//    if(!strcmp("EYieldSurface3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface3;
//    if(!strcmp("EPOrder",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPOrder;
//    if(!strcmp("ENSteps",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ENSteps;
//    if(!strcmp("EPorePressure",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPorePressure;
//    if(!strcmp("EMatPorosity",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatPorosity;
//    if(!strcmp("EMatE",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatE;
//    if(!strcmp("EMatPoisson",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatPoisson;
    
    //return TPZMatWithMem<TMEM>::VariableIndex(name);
    PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
    return -1;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZMatElasticity2D::NSolutionVariables(int var){
    if(var == 1)	return 3;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainVol)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXX	)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYY)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainZZ)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXY)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXZ)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYZ)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainVol)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXX)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainZZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainVol)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXX)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainZZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2El)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlAlpha)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementX)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementTotal)	 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressI1)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressJ2)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXX)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressZZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress1)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress2)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress3)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressI1)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressJ2)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXX)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressZZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXY)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYZ)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress1)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress2)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress3)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface1)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface2)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface3)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPOrder)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ENSteps)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPorePressure)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatPorosity)		 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatE)			 return -1;
//    if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatPoisson)		 return -1;
    
    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZMatElasticity2D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    TPZManVector<STATE,3> SolU, SolP;
    TPZFNMatrix <6,STATE> DSolU, DSolP;
    TPZFNMatrix <9> axesU, axesP;
    
    TPZVec<REAL> ptx(3);
    TPZVec<STATE> solExata(3);
    TPZFMatrix<STATE> flux(5,1);
    
    if (data.sol.size() != 1) {
        DebugStop();
    }
    
    SolU=data.sol[0];
    DSolU=data.dsol[0];
    axesU=data.axes;
    
    
    //	Displacements
    if(var == 1){
        Solout[0] = SolU[0];
        Solout[1] = SolU[1];
        Solout[2] = 0.0;
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
    REAL C11 = 4*(fmu)*(flambda+fmu)/(flambda+2*fmu);
    REAL C22 = 2*(fmu)*(flambda)/(flambda+2*fmu);
    
    if (this->fPlaneStress)
    {
        SigX = C11*epsx+C22*epsy;
        SigY = C11*epsy+C22*epsx;
        SigZ = 0.0;
        Tau = 2.0*fmu*epsxy;
    }
    else
    {
        SigX = ((flambda + 2*fmu)*(epsx) + (flambda)*epsy);
        SigY = ((flambda + 2*fmu)*(epsy) + (flambda)*epsx);
        SigZ = flambda*divu;
        Tau = 2.0*fmu*epsxy;		
    }
    
    REAL R = sqrt( ((SigX - SigY)/2)*((SigX - SigY)/2) + Tau*Tau);
    REAL C = (SigX + SigY)/2;
    REAL Sigma1 = C + R;
    REAL Sigma2 = C - R;
    
    // Sigma1 is the largest stress regadless of sign
    if (abs(Sigma2) > abs(Sigma1)) 
    {
        REAL tmp = Sigma1;
        Sigma1 = Sigma2;
        Sigma2 = tmp;
    }
    
    //	Hydrostatic stress
    if(var == 2) 
    {
        Solout[0] = SigX+SigY+SigZ;
        return;
    }
    
    //	Effective Stress x-direction
    if(var == 3) {
        Solout[0] = SigX;
        return;
    }
    
    //	Effective Stress y-direction	
    if(var == 4) {
        Solout[0] = SigY;
        return;
    }
    
    //	Effective Stress y-direction
    if(var == 5) {
        Solout[0] = SigZ;
        return;
    }
    
    //	Shear Stress	
    if(var == 6) {
        Solout[0] = Tau;
        return;
    }
    
}


