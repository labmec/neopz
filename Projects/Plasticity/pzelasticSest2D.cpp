//
//  pzelasticSest2D.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/23/14.
//
//

#include "pzelasticSest2D.h"

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D():TPZMatElasticity2D()
{
    fBiotAlpha = 0.;
    fZDeformation = 0.;
}

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(int id, REAL E, REAL nu, REAL fx, REAL fy, int plainstress):TPZMatElasticity2D(id,E,nu,fx,fy,plainstress)
{
    fBiotAlpha = 0.;
    fZDeformation = 0.;
}

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(int id):TPZMatElasticity2D(id)
{
    fBiotAlpha = 0.;
    fZDeformation = 0.;
}

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(const TPZElasticityMaterialSest2D &copy) : TPZMatElasticity2D(copy)
{
    fBiotAlpha = copy.fBiotAlpha;
    fZDeformation = copy.fZDeformation;
}

TPZElasticityMaterialSest2D::~TPZElasticityMaterialSest2D()
{
    
}

TPZElasticityMaterialSest2D & TPZElasticityMaterialSest2D::operator=(const TPZElasticityMaterialSest2D &copy)
{
    TPZMatElasticity2D::operator=(copy);
    fBiotAlpha = copy.fBiotAlpha;
    fZDeformation = copy.fZDeformation;
    return *this;
}

int TPZElasticityMaterialSest2D::VariableIndex(const std::string &name)
{
    if(!strcmp("StrainVol",		name.c_str()))  return TPZElasticityMaterialSest2D::EStrainVol;
    if(!strcmp("StrainXX",		name.c_str()))  return TPZElasticityMaterialSest2D::EStrainXX;
    if(!strcmp("StrainYY",		name.c_str()))  return TPZElasticityMaterialSest2D::EStrainYY;
    if(!strcmp("StrainZZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EStrainZZ;
    if(!strcmp("StrainXY",		name.c_str()))  return TPZElasticityMaterialSest2D::EStrainXY;
    if(!strcmp("StrainXZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EStrainXZ;
    if(!strcmp("StrainYZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EStrainYZ;
    if(!strcmp("ElStrainVol",		name.c_str()))  return TPZElasticityMaterialSest2D::EElStrainVol;
    if(!strcmp("ElStrainXX",		name.c_str()))  return TPZElasticityMaterialSest2D::EElStrainXX;
    if(!strcmp("ElStrainYY",		name.c_str()))  return TPZElasticityMaterialSest2D::EElStrainYY;
    if(!strcmp("ElStrainZZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EElStrainZZ;
    if(!strcmp("ElStrainXY",		name.c_str()))  return TPZElasticityMaterialSest2D::EElStrainXY;
    if(!strcmp("ElStrainXZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EElStrainXZ;
    if(!strcmp("ElStrainYZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EElStrainYZ;
    if(!strcmp("PlStrainVol",		name.c_str()))  return TPZElasticityMaterialSest2D::EPlStrainVol;
    if(!strcmp("PlStrainXX",		name.c_str()))  return TPZElasticityMaterialSest2D::EPlStrainXX;
    if(!strcmp("PlStrainYY",		name.c_str()))  return TPZElasticityMaterialSest2D::EPlStrainYY;
    if(!strcmp("PlStrainZZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EPlStrainZZ;
    if(!strcmp("PlStrainXY",		name.c_str()))  return TPZElasticityMaterialSest2D::EPlStrainXY;
    if(!strcmp("PlStrainXZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EPlStrainXZ;
    if(!strcmp("PlStrainYZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EPlStrainYZ;
    if(!strcmp("PlStrainSqJ2",		name.c_str()))  return TPZElasticityMaterialSest2D::EPlStrainSqJ2;
    if(!strcmp("PlStrainSqJ2El",		name.c_str()))  return 100;//return TPZElasticityMaterialSest2D::EPlStrainSqJ2El;
    if(!strcmp("PlAlpha",			name.c_str()))  return TPZElasticityMaterialSest2D::EPlAlpha;
    if(!strcmp("DisplacementX",		name.c_str()))  return TPZElasticityMaterialSest2D::EDisplacementX;
    if(!strcmp("DisplacementY",		name.c_str()))  return TPZElasticityMaterialSest2D::EDisplacementY;
    if(!strcmp("DisplacementZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EDisplacementZ;
    if(!strcmp("DisplacementTotal",	name.c_str()))  return TPZElasticityMaterialSest2D::EDisplacementTotal;
    if(!strcmp("TotStressI1",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStressI1;
    if(!strcmp("TotStressJ2",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStressJ2;
    if(!strcmp("TotStressXX",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStressXX;
    if(!strcmp("TotStressYY",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStressYY;
    if(!strcmp("TotStressZZ",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStressZZ;
    if(!strcmp("TotStressXY",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStressXY;
    if(!strcmp("TotStressXZ",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStressXZ;
    if(!strcmp("TotStressYZ",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStressYZ;
    if(!strcmp("TotStress1",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStress1;
    if(!strcmp("TotStress2",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStress2;
    if(!strcmp("TotStress3",		name.c_str()))  return TPZElasticityMaterialSest2D::ETotStress3;
    if(!strcmp("EffStressI1",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStressI1;
    if(!strcmp("EffStressJ2",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStressJ2;
    if(!strcmp("EffStressXX",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStressXX;
    if(!strcmp("EffStressYY",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStressYY;
    if(!strcmp("EffStressZZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStressZZ;
    if(!strcmp("EffStressXY",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStressXY;
    if(!strcmp("EffStressXZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStressXZ;
    if(!strcmp("EffStressYZ",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStressYZ;
    if(!strcmp("EffStress1",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStress1;
    if(!strcmp("EffStress2",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStress2;
    if(!strcmp("EffStress3",		name.c_str()))  return TPZElasticityMaterialSest2D::EEffStress3;
    if(!strcmp("YieldSurface1",		name.c_str()))  return TPZElasticityMaterialSest2D::EYieldSurface1;
    if(!strcmp("YieldSurface2",		name.c_str()))  return TPZElasticityMaterialSest2D::EYieldSurface2;
    if(!strcmp("YieldSurface3",		name.c_str()))  return TPZElasticityMaterialSest2D::EYieldSurface3;
    if(!strcmp("POrder",			name.c_str()))  return TPZElasticityMaterialSest2D::EPOrder;
    if(!strcmp("NSteps",			name.c_str()))  return TPZElasticityMaterialSest2D::ENSteps;
    if(!strcmp("PorePressure",		name.c_str()))  return TPZElasticityMaterialSest2D::EPorePressure;
    if(!strcmp("MatPorosity",		name.c_str()))  return TPZElasticityMaterialSest2D::EMatPorosity;
    if(!strcmp("MatE",			name.c_str()))  return TPZElasticityMaterialSest2D::EMatE;
    if(!strcmp("MatPoisson",		name.c_str()))  return TPZElasticityMaterialSest2D::EMatPoisson;
    
    //return TPZMatWithMem<TMEM>::VariableIndex(name);
    PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
    DebugStop();
    return -1;
}

int TPZElasticityMaterialSest2D::NSolutionVariables(int var)
{
    if(var == TPZElasticityMaterialSest2D::EStrainVol)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EStrainXX)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EStrainYY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EStrainZZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EStrainXY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EStrainXZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EStrainYZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EElStrainVol)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EElStrainXX)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EElStrainYY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EElStrainZZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EElStrainXY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EElStrainXZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EElStrainYZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainVol)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainXX)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainYY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainZZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainXY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainXZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainYZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainSqJ2)	 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlStrainSqJ2El)	 return 1;
    if(var == TPZElasticityMaterialSest2D::EPlAlpha)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EDisplacementX)	 return 1;
    if(var == TPZElasticityMaterialSest2D::EDisplacementY)	 return 1;
    if(var == TPZElasticityMaterialSest2D::EDisplacementZ)	 return 1;
    if(var == TPZElasticityMaterialSest2D::EDisplacementTotal)	 return 2;
    if(var == TPZElasticityMaterialSest2D::ETotStressI1)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStressJ2)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStressXX)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStressYY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStressZZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStressXY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStressXZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStressYZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStress1)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStress2)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ETotStress3)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStressI1)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStressJ2)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStressXX)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStressYY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStressZZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStressXY)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStressXZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStressYZ)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStress1)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStress2)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EEffStress3)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EYieldSurface1)	 return 1;
    if(var == TPZElasticityMaterialSest2D::EYieldSurface2)	 return 1;
    if(var == TPZElasticityMaterialSest2D::EYieldSurface3)	 return 1; // Should never be called
    if(var == TPZElasticityMaterialSest2D::EPOrder)		 return 1;
    if(var == TPZElasticityMaterialSest2D::ENSteps)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EPorePressure)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EMatPorosity)		 return 1;
    if(var == TPZElasticityMaterialSest2D::EMatE)			 return 1;
    if(var == TPZElasticityMaterialSest2D::EMatPoisson)		 return 1;
    if(var == 100) return 1;
    return TPZElasticityMaterialSest2D::NSolutionVariables(var);
}


/// Compute post processed values
void TPZElasticityMaterialSest2D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    TPZManVector<STATE,3> SolU, SolP(1,0.0);
    TPZFNMatrix <6,STATE> DSolU, DSolP(3,1,0.0);
    TPZFNMatrix <9> axesU, axesP;
    SolU	=	data.sol[0];
    DSolU	=	data.dsol[0];
    axesU	=	data.axes;
    
    if(this->HasForcingFunction())
    {
        fForcingFunction->Execute(data.x,SolP,DSolP);
    }
    
    REAL epsx, epsy, epsxy, epsz, SigX, SigY, SigZ, Tau, DSolxy[2][2], divu;
    DSolxy[0][0] = DSolU(0,0)*axesU(0,0)+DSolU(1,0)*axesU(1,0); // dUx/dx
    DSolxy[1][0] = DSolU(0,0)*axesU(0,1)+DSolU(1,0)*axesU(1,1); // dUx/dy
    DSolxy[0][1] = DSolU(0,1)*axesU(0,0)+DSolU(1,1)*axesU(1,0); // dUy/dx
    DSolxy[1][1] = DSolU(0,1)*axesU(0,1)+DSolU(1,1)*axesU(1,1); // dUy/dy
    
    epsx = DSolxy[0][0];	// du/dx
    epsy = DSolxy[1][1];	// dv/dy
    epsxy = 0.5*(DSolxy[1][0]+DSolxy[0][1]);
    REAL C11 = 4*(fmu)*(flambda+fmu)/(flambda+2*fmu);
    REAL C22 = 2*(fmu)*(flambda)/(flambda+2*fmu);
    
    if (this->fPlaneStress)
    {
        epsz = -(flambda/(flambda+2*fmu))*(epsx+epsy);
        divu = epsx+epsy+epsz;
        SigX = C11*epsx+C22*epsy;
        SigY = C11*epsy+C22*epsx;
        SigZ = 0.0;
        Tau = 2.0*fmu*epsxy;
    }
    else
    {
        epsz = -(flambda/(flambda+2*fmu))*(epsx+epsy);
        divu = epsx+epsy+epsz;
        SigX = ((flambda + 2*fmu)*(epsx) + (flambda)*epsy);
        SigY = ((flambda + 2*fmu)*(epsy) + (flambda)*epsx);
        SigZ = flambda*divu;
        Tau = 2.0*fmu*epsxy;
    }
    
    TPZTensor<REAL> totalStress;
    totalStress.XX()	=	SigX;
    totalStress.YY()	=	SigY;
    totalStress.ZZ()	=	SigZ;
    totalStress.XY()	=	Tau;
    totalStress.YZ()	=	0.0;
    
    TPZTensor<REAL> Sigma;
    Sigma.XX()	=	SigX - fBiotAlpha * SolP[0];
    Sigma.YY()	=	SigY - fBiotAlpha * SolP[0];
    Sigma.ZZ()	=	SigZ - fBiotAlpha * SolP[0];
    Sigma.XY()	=	Tau;
    Sigma.YZ()	=	0.0;
    
    REAL R = sqrt( ((SigX - SigY)/2)*((SigX - SigY)/2) + Tau*Tau);
    REAL C = (SigX + SigY)/2;
    REAL Sigma1 = C + R;
    REAL Sigma2 = C - R;
    
    switch (var) {
            // Total Strain
        case EStrainVol:{
            Solout[0] = divu;}
            break;
        case EStrainXX:{
            Solout[0] = epsx;}
            break;
        case EStrainYY:{
            Solout[0] = epsy;}
            break;
        case EStrainZZ:{
            Solout[0] = epsz;}
            break;
        case EStrainXY:{
            Solout[0] = epsxy;}
            break;
        case EStrainXZ:{
            Solout[0] = 0.0;}
            break;
        case EStrainYZ:{
            Solout[0] = 0.0;}
            break;
            // Elastic Strain
        case EElStrainVol:{
            Solout[0] = divu;}
            break;
        case EElStrainXX:{
            Solout[0] = epsx;}
            break;
        case EElStrainYY:{
            Solout[0] = epsy;}
            break;
        case EElStrainZZ:{
            Solout[0] = epsz;}
            break;
        case EElStrainXY:{
            Solout[0] = epsxy;}
            break;
        case EElStrainXZ:{
            Solout[0] = 0.0;}
            break;
        case EElStrainYZ:{
            Solout[0] = 0.0;}
            break;
            // Plastic Strain
        case EPlStrainVol:{
            Solout[0] = 0.0;}
            break;
        case EPlStrainXX:{
            Solout[0] = 0.0;}
            break;
        case EPlStrainYY:{
            Solout[0] = 0.0;}
            break;
        case EPlStrainZZ:{
            Solout[0] = 0.0;}
            break;
        case EPlStrainXY:{
            Solout[0] = 0.0;}
            break;
        case EPlStrainXZ:{
            Solout[0] = 0.0;}
            break;
        case EPlStrainYZ:{
            Solout[0] = 0.0;}
            break;
            // SqJ2 and alpha
        case EPlStrainSqJ2:{
            Solout[0] = 0.0;}
            break;
        case EPlStrainSqJ2El:{
            DebugStop();}
            break;
        case EPlAlpha:{
            Solout[0] = 0.0;}
            break;
            // Displacement
        case EDisplacementX:{
            Solout[0] = SolU[0];}
            break;
        case EDisplacementY:{
            Solout[0] = SolU[1];}
            break;
        case EDisplacementZ:{
            Solout[0] = 0.;} // Always usnig plane strain
            break;
        case EDisplacementTotal:{
            for (int i = 0; i < 2; i++){
                Solout[i] = SolU[i];}
        }
            break;
            // Total Stress
        case ETotStressI1:{
            Solout[0] = totalStress.I1();}
            break;
        case ETotStressJ2:{
            Solout[0] = totalStress.J2();}
            break;
        case ETotStressXX:{
            Solout[0] = totalStress.XX();}
            break;
        case ETotStressYY:{
            Solout[0] = totalStress.YY();}
            break;
        case ETotStressZZ:{
            Solout[0] = totalStress.ZZ();}
            break;
        case ETotStressXY:{
            Solout[0] = totalStress.XY();}
            break;
        case ETotStressXZ:{
            Solout[0] = totalStress.XZ();}
            break;
        case ETotStressYZ:{
            Solout[0] = totalStress.YZ();}
            break;
        case ETotStress1:
        {
            TPZTensor<STATE> eigenval;
            totalStress.EigenValue(eigenval);
            Solout[0] = eigenval.XX();
        }
            break;
        case ETotStress2:
        {
            TPZTensor<STATE> eigenval;
            totalStress.EigenValue(eigenval);
            Solout[0] = eigenval.YY();
        }
            break;
        case ETotStress3:
        {
            TPZTensor<STATE> eigenval;
            totalStress.EigenValue(eigenval);
            Solout[0] = eigenval.ZZ();
        }
            break;
            // Effective stress
        case EEffStressI1:{
            Solout[0] = Sigma.I1();}
            break;
        case EEffStressJ2:{
            Solout[0] = Sigma.J2();}
            break;
        case EEffStressXX:{
            Solout[0] = Sigma.XX();}
            break;
        case EEffStressYY:{
            Solout[0] = Sigma.YY();}
            break;
        case EEffStressZZ:{
            Solout[0] = Sigma.ZZ();}
            break;
        case EEffStressXY:{
            Solout[0] = Sigma.XY();}
            break;
        case EEffStressXZ:{
            Solout[0] = Sigma.XZ();}
            break;
        case EEffStressYZ:{
            Solout[0] = Sigma.YZ();}
            break;
        case EEffStress1:
        {
            TPZTensor<STATE> eigenval;
            Sigma.EigenValue(eigenval);
            Solout[0] = eigenval.XX();
        }
            break;
        case EEffStress2:
        {
            TPZTensor<STATE> eigenval;
            Sigma.EigenValue(eigenval);
            Solout[0] = eigenval.YY();
        }
            break;
        case EEffStress3:
        {
            TPZTensor<STATE> eigenval;
            Sigma.EigenValue(eigenval);
            Solout[0] = eigenval.ZZ();
        }
            break;
            // Yield Surface
        case EYieldSurface1:
        {
            Solout[0] = 0.0;
        }
            break;
        case EYieldSurface2:
        {
            Solout[0] = 0.0;
        }
            break;
        case EYieldSurface3:
        {
            Solout[0] = 0.0;
        }
            break;
            // Simulation
        case EPOrder:{
            Solout[0] = data.p;}
            break;
        case ENSteps:{
            Solout[0] = 0;}
            break;
            // Pore pressure
        case EPorePressure:{
            Solout[0] = SolP[0];}
            break;
            // Material
        case EMatPorosity:{
            Solout[0] = 6378;} // AQUINATHAN raio da terra}
            break;
        case EMatE:{
            Solout[0] = this->fE;}
            break;
        case EMatPoisson:{
            Solout[0] = this->fnu;}
            break;
        default:{
            DebugStop();}
            break;
    }
    
}


void TPZElasticityMaterialSest2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    //  Contribution of domain integrals for Jacobian matrix
    //  Elasticity Block (Equation for elasticity )
    //	Elastic equation
    //	Linear strain operator
    //	Ke Matrix
    TPZMatElasticity2D::Contribute(data, weight, ek,ef);
    //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiU =  data.phi;
    TPZFMatrix<REAL> &dphiU = data.dphix;
    int phrU = phiU.Rows();
    int FirstU  = 0;
    
    TPZManVector<REAL,3> sol_u =data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    REAL LambdaL, MuL, Pressure;
    
    LambdaL = flambda;
    MuL     = fmu;
    
    TPZVec<STATE> P(2,0.0);
    TPZFMatrix<STATE> GradP(2,1,0.0);
    
    if(this->HasForcingFunction())
    {
        fForcingFunction->Execute(data.x,P,GradP);
        Pressure = P[0];
    }
    
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
    
    //  Matrix Qc
    //  Coupling matrix for Elastic equation
    for(int in = 0; in < phrU; in++ )
    {
        du(0,0) = dphiU(0,in)*data.axes(0,0)+dphiU(1,in)*data.axes(1,0);
        du(1,0) = dphiU(0,in)*data.axes(0,1)+dphiU(1,in)*data.axes(1,1);
        
        ef(2*in + FirstU)    += (-1.0)* fBiotAlpha * weight * (Pressure*du(0,0));
        ef(2*in+1 + FirstU)  += (-1.0)* fBiotAlpha * weight * (Pressure*du(1,0));
    }
    
}
void TPZElasticityMaterialSest2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    
    TPZMatElasticity2D::Contribute(data,weight,ef);
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiU =  data.phi;
    TPZFMatrix<REAL> &dphiU = data.dphix;
    int phrU = phiU.Rows();
    int FirstU  = 0;
    
    TPZManVector<REAL,3> sol_u =data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    REAL LambdaL, MuL, Pressure;
    
    LambdaL = flambda;
    MuL     = fmu;
    
    TPZVec<STATE> P(2,0.0);
    TPZFMatrix<STATE> GradP(2,1,0.0);
    
    if(this->HasForcingFunction())
    {
        fForcingFunction->Execute(data.x,P,GradP); //AQUIPHIL o que o omar recebe da TPBRBiotForce
        Pressure = P[0];
    }
    
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
    
    //  Matrix Qc
    //  Coupling matrix for Elastic equation
    for(int in = 0; in < phrU; in++ )
    {
        du(0,0) = dphiU(0,in)*data.axes(0,0)+dphiU(1,in)*data.axes(1,0);
        du(1,0) = dphiU(0,in)*data.axes(0,1)+dphiU(1,in)*data.axes(1,1);
        
        ef(2*in + FirstU)    += (-1.0)* fBiotAlpha * weight * (Pressure*du(0,0));
        ef(2*in+1 + FirstU)  += (-1.0)* fBiotAlpha * weight * (Pressure*du(1,0));
    }
    
}


void TPZElasticityMaterialSest2D::Write(TPZStream &buf, int withclassid)
{
    TPZMatElasticity2D::Write(buf,withclassid);
    buf.Write(&fBiotAlpha);
    buf.Write(&fZDeformation);
}

void TPZElasticityMaterialSest2D::Read(TPZStream &buf, void *context)
{
    TPZMatElasticity2D::Read(buf,context);
    buf.Read(&fBiotAlpha);
    buf.Read(&fZDeformation);
}

void TPZElasticityMaterialSest2D::Print(std::ostream &out)
{
    out << "TPZElasticityMaterialSest2D :" << std::endl;
    out << "\t fBiotAlpha   = "			<< fBiotAlpha << std::endl;
    out << "\t fZDeformation   = "      << fZDeformation<< std::endl;
    out << "Father properties :";
    TPZMatElasticity2D::Print(out);
    out << "\n";
}
