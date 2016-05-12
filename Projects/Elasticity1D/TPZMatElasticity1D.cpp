//
//  TPZMatElasticity1D.cpp
//  PZ
//
//  Created by Nathalia on 4/15/16.
//
//

#include "TPZMatElasticity1D.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity1D"));
#endif





TPZMatElasticity1D::TPZMatElasticity1D():TPZMaterial()
{
    flambda = 0.;
    fmu = 0.;
    fb.resize(2);
    fb[0]=0.;
    fb[1]=0.;
    fLineStrain=0.;

}


TPZMatElasticity1D::TPZMatElasticity1D(int matid, REAL Lambda, REAL Mu, REAL b, int LineStrain):TPZMaterial(matid)
{
    flambda = 0.;
    fmu = 0.;
    fb.resize(2);
    fb[0]=0.;
    fb[1]=0.;
    fLineStrain=0.;
}


TPZMatElasticity1D::TPZMatElasticity1D(int matid):TPZMaterial(matid)
{
    flambda = 0.;
    fmu = 0.;
    fb.resize(2);
    fb[0]=0.;
    fb[1]=0.;
    fLineStrain=0.;

}



TPZMatElasticity1D & TPZMatElasticity1D::operator=(const TPZMatElasticity1D &copy)
{
    TPZMaterial::operator = (copy);
    flambda = copy.flambda;
    fmu = copy.fmu;
    fb.resize(copy.fb.size());
    for (int i = 0; i < copy.fb.size(); i++) {
        fb[i] = copy.fb[i];
    }
    fLineStrain = copy.fLineStrain;
    return *this;

}



TPZMatElasticity1D::~TPZMatElasticity1D()
{
}


TPZMatElasticity1D::TPZMatElasticity1D(const TPZMatElasticity1D & cp):TPZMaterial(cp)
{
    flambda = cp.flambda;
    fmu = cp.fmu;
    fb.resize(cp.fb.size());
    for (int i = 0; i < cp.fb.size(); i++) {
        fb[i] = cp.fb[i];
    }
    fLineStrain = cp.fLineStrain;

}


int TPZMatElasticity1D::NStateVariables() // retorna 1, pois tenho apenas uma variavel de estado
{
    return 1;
}




void TPZMatElasticity1D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

    
    //Getting weight functions
    TPZFNMatrix<220,REAL> & phi_u  = data.phi;
    TPZFNMatrix<660,REAL> & dphi_u_axes = data.dphix;
    TPZFNMatrix<9,REAL> axes = data.axes;
    
    
    int nphi = phi_u.Rows();
    
    int FirstU  = 0;
    
    REAL LambdaL, MuL;
    LambdaL     = flambda;
    MuL         = fmu;
    
    // Functions computed at point x_{k} for each integration point

    // Calcula "du" multiplicando as funcoes pelas coordenadas (axes), para o plano desejado (nesse caso x)
    
    TPZFMatrix<REAL>	Gradvi(1,1),Gradvj(1,1);
    for(int iu = 0; iu < nphi; iu++ )
    {
        //	Derivative for Vx, somente em x pois eh 1D (em y e z serao = 0)
        Gradvi(0,0) = dphi_u_axes(0,iu)*data.axes(0,0)+dphi_u_axes(0,iu)*data.axes(0,1)+dphi_u_axes(0,iu)*data.axes(0,2);
        
        for(int ju = 0; ju < nphi; ju++)
        {
            //	Derivative for Ux, somente em x pois eh 1D (em y e z serao = 0)
            Gradvj(0,0) = dphi_u_axes(0,ju)*data.axes(0,0)+dphi_u_axes(0,ju)*data.axes(0,1)+dphi_u_axes(0,ju)*data.axes(0,2);

        if (this->fLineStrain == 1)
        {
             /* Line Strain State */
            ek(1*iu + FirstU,1*ju + FirstU)         += weight*(LambdaL+2.0*MuL)*Gradvi(0,0)*Gradvj(0,0);
            
        }
    }
}
//  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    this->Contribute(data,weight,ef);
}


void TPZMatElasticity1D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){

  
    //Getting weight functions
    TPZFNMatrix<220,REAL> & phi_u  = data.phi;
    TPZFNMatrix<660,REAL> & dphi_u_axes = data.dphix;
    TPZFNMatrix<9,REAL> axes = data.axes;
    TPZFNMatrix<15,REAL> & du_daxes = data.dsol[0];
    
    int nphi = phi_u.Rows();
    
    TPZManVector<REAL,3> sol_u =data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];

    int FirstU  = 0;
    
    REAL LambdaL, MuL;
    LambdaL     = flambda;
    MuL         = fmu;

    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  Contribution of domain integrals for Residual Vector
    //  Elastic equation
    //  Linear strain operator
    //  Ke Matrix
    
    TPZFMatrix<REAL>    du(2,2);
    TPZFMatrix<REAL>    dux(2,2);

    //  Derivative for Ux
    TPZFMatrix<REAL>	Gradvi(1,1),Gradu(1,1);
    
    //	Derivative for Ux, somente em x pois eh 1D (em y e z serao = 0)
    Gradu(0,0) = du_daxes(0,0)*data.axes(0,0)+du_daxes(0,0)*data.axes(0,1)+du_daxes(0,0)*data.axes(0,2);
    
    for(int iu = 0; iu < nphi; iu++ )
    {
        //	Derivative for Vx, somente em x pois eh 1D (em y e z serao = 0)
        Gradvi(0,0) = dphi_u_axes(0,iu)*data.axes(0,0)+dphi_u_axes(0,iu)*data.axes(0,1)+dphi_u_axes(0,iu)*data.axes(0,2);

        //  Vector Force right hand term
        ef(1*iu + FirstU)     +=    -1.0*weight*fb[0]*data.x[0]*phi_u(iu, 0);    // direcao x

        if (this->fLineStrain == 1)
        {
            /* Line Strain State */
            
            ef(1*iu + FirstU)           += weight*  ((LambdaL + 2.0*MuL)*Gradvi(0,0)*Gradu(0,0));
            

        }
    }

}



void TPZMatElasticity1D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){

    
    TPZFNMatrix<220,REAL> & phi_u  = data.phi;
    TPZManVector<REAL,3> sol_u =data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    
    int nphi = phi_u.Rows();
    short in,jn;
    STATE v2[3]; // ja criado para 3 direcoes
    v2[0] = bc.Val2()(0,0);	//	Ux displacement or Tnx

    //	Here each digit represent an individual boundary condition corresponding to each state variable.
    //	0 means Dirichlet condition on x
    //	1 means Neumann condition
    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    
    switch (bc.Type())
    {
        case 0 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            
            for(in = 0 ; in < nphi; in++)
            {
                //	Contribution for load Vector
                ef(1*in,0)      += BIGNUMBER*(ux - v2[0])*phi_u(in,0)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < nphi; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(1*in,1*jn)       += BIGNUMBER*phi_u(in,0)*phi_u(jn,0)*weight;	// X displacement
                }
            }
            
            break;
        }

        case 1 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <nphi; in++)
            {
                //	Normal Tension Components on neumann boundary
                ef(1*in,0)      += -1.0*v2[0]*phi_u(in,0)*weight;		//	Tnx
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


void TPZMatElasticity1D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){

    
    TPZFNMatrix<220,REAL> & phi_u  = data.phi;
    TPZManVector<REAL,3> sol_u =data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    
    int nphi = phi_u.Rows();
    short in;
    STATE v2[3]; // ja criado para 3 direcoes
    v2[0] = bc.Val2()(0,0);	//	Ux displacement or Tnx
    
    //	Here each digit represent an individual boundary condition corresponding to each state variable.
    //	0 means Dirichlet condition on x
    //	1 means Neumann condition

    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    
    switch (bc.Type())
    {
        case 0 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            
            for(in = 0 ; in < nphi; in++)
            {
                //	Contribution for load Vector
                ef(1*in,0)      += BIGNUMBER*(ux - v2[0])*phi_u(in,0)*weight;	// X displacement Value
                
            }
            
            break;
        }
            
            
        case 1 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <nphi; in++)
            {
                //	Normal Tension Components on neumann boundary
                ef(1*in,0)      += -1.0*v2[0]*phi_u(in,0)*weight;		//	Tnx
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


//****************************************//

void TPZMatElasticity1D::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNeighborSol = true;
    data.fNeedsNeighborCenter = false;
    data.fNeedsNormal = true;
}

void TPZMatElasticity1D::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}


void TPZMatElasticity1D::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "1D Problem (for Line Strain conditions) " << fLineStrain << std::endl;
    out << "Properties for elasticity: \n";
    out << "\t First Lamé Parameter   = "									<< flambda	<< std::endl;
    out << "\t Second Lamé Parameter   = "									<< fmu		<< std::endl;
    out << "\t Body force vector B {X-direction}   = "                      << fb[0] << ' ' << std::endl;
    out << "Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
    
}

//*****************************************//

/** Returns the variable index associated with the name */
int TPZMatElasticity1D::VariableIndex(const std::string &name)
{
    //	Elasticity Variables
    if(!strcmp("Displacement",name.c_str()))				return	1;
    if(!strcmp("SigmaX",name.c_str()))						return	3;
    PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
    return -1;
    
    return TPZMaterial::VariableIndex(name);
}


/**
 * Save the element data to a stream
 */
void TPZMatElasticity1D::Write(TPZStream &buf, int withclassid)
{
    TPZMaterial::Write(buf,withclassid);
    buf.Write(&flambda);
    buf.Write(&fmu);
    TPZSaveable::WriteObjects(buf, fb);
    
}


/**
 * Read the element data from a stream
 */
void TPZMatElasticity1D::Read(TPZStream &buf, void *context)
{
    TPZMaterial::Read(buf,context);
    buf.Read(&flambda);
    buf.Read(&fmu);
    TPZSaveable::ReadObjects(buf, fb);
    
}



//************************************//

int TPZMatElasticity1D::NSolutionVariables(int var){
    if(var == 1)	return 3;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    if(var == 6)	return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}


//******************************************//

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZMatElasticity1D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
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
    if(var == 1){
        Solout[0] = SolU[0];
        return;
    }
    
    
    REAL epsx;
    REAL SigX;
    REAL DSolxy[1][1];
    REAL divu;
    
    DSolxy[0][0] = DSolU(0,0)*axesU(0,0); //+DSolU(1,0)*axesU(1,0); // dUx/dx

    divu = DSolxy[0][0];
    
    epsx = DSolxy[0][0];// du/dx

    
    if (this->fLineStrain)
    {
       SigX = ((flambda + 2*fmu)*(epsx));

    }
    
    
    //	Hydrostatic stress
    if(var == 2)
    {
        Solout[0] = SigX;
        return;
    }
    
    //	Effective Stress x-direction
    if(var == 3) {
        Solout[0] = SigX;
        return;
    }

    
}





//
//#endif /* TPZMatElasticity1D_hpp */
