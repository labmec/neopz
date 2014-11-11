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


TPZMatElasticity2D::TPZMatElasticity2D():TPZMaterial()
{
    fE = 0.;
    fnu = 0.;
    flambda = 0.;
    fmu = 0.;
    ff.resize(3);
    ff[0]=0.;
    ff[1]=0.;
    ff[2]=0.;
    fPlaneStress = 1.;
    fPreStressXX = 0.0;
    fPreStressXY = 0.0;
    fPreStressYY = 0.0;
    fPreStressZZ = 0.0;
    
}

TPZMatElasticity2D::TPZMatElasticity2D(int matid):TPZMaterial(matid)
{
    fE = 0.;
    fnu = 0.;
    flambda = 0.;
    fmu = 0.;
    ff.resize(3);
    ff[0]=0.;
    ff[1]=0.;
    ff[2]=0.;
    fPlaneStress = 1.;
    fPreStressXX = 0.0;
    fPreStressXY = 0.0;
    fPreStressYY = 0.0;
    fPreStressZZ = 0.0;
}

TPZMatElasticity2D::TPZMatElasticity2D(int matid, REAL E, REAL nu, REAL fx, REAL fy, int plainstress):TPZMaterial(matid)
{
    fE = E;
    fnu = nu;
    flambda = (E*nu)/((1+nu)*(1-2*nu));
    fmu = E/(2*(1+nu));
    ff.resize(3);
    ff[0]=fx;
    ff[1]=fy;
    ff[2]=0.0;
    fPlaneStress = plainstress;
    fPreStressXX = 0.0;
    fPreStressXY = 0.0;
    fPreStressYY = 0.0;
    fPreStressZZ = 0.0;
}

TPZMatElasticity2D::~TPZMatElasticity2D()
{
}


TPZMatElasticity2D::TPZMatElasticity2D(const TPZMatElasticity2D &copy) : TPZMaterial(copy)
{
    fE = copy.fE;
    fnu = copy.fnu;
    flambda = copy.flambda;
    fmu = copy.fmu;
    ff.resize(copy.ff.size());
    for (int i = 0; i < copy.ff.size(); i++) {
        ff[i] = copy.ff[i];
    }
    fPlaneStress = copy.fPlaneStress;
    fPreStressXX = copy.fPreStressXX;
    fPreStressXY = copy.fPreStressXY;
    fPreStressYY = copy.fPreStressYY;
    fPreStressZZ = copy.fPreStressZZ;
}

TPZMatElasticity2D & TPZMatElasticity2D::operator=(const TPZMatElasticity2D &copy)
{
	TPZMaterial::operator = (copy);
    fE = copy.fE;
    fnu = copy.fnu;
    flambda = copy.flambda;
    fmu = copy.fmu;
    fPreStressXX = copy.fPreStressXX;
    fPreStressXY = copy.fPreStressXY;
    fPreStressYY = copy.fPreStressYY;
    fPreStressZZ = copy.fPreStressZZ;
    ff.resize(copy.ff.size());
    for (int i = 0; i < copy.ff.size(); i++) {
        ff[i] = copy.ff[i];
    }
    fPlaneStress = copy.fPlaneStress;
    return *this;
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
    REAL LambdaL, MuL, Pressure;
    
    LambdaL = flambda;
    MuL     = fmu;
    
    TPZVec<STATE> P(1,0.0);
    TPZFMatrix<STATE> GradP(2,1,0.0);
    
    if(this->HasffBCForcingFunction())
    {
	fForcingFunction->Execute(data.x,P,GradP);
	Pressure = P[0];
    }
    
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
             ef(2*iu + FirstU)     +=    weight*ff[0]*phiU(iu, 0)- (du(0,0)*fPreStressXX + du(1,0)*fPreStressXY);    // direcao x
             ef(2*iu+1 + FirstU)   +=    weight*ff[1]*phiU(iu, 0)- (du(0,0)*fPreStressXY + du(1,0)*fPreStressYY);    // direcao y
        
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
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    
}

void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    /*
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
     */
    
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<REAL,3> sol_u = data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    REAL uy = sol_u[1];
    
    int phru = phiu.Rows();
    short in,jn;
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
                res(i,0) += bc.Val1()(i,j)*data.sol[0][j];
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
                ef(2*in,0)      += BIGNUMBER*(0.0 - v2[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= BIGNUMBER*(0.0 - v2[1])*phiu(in,0)*weight;	// y displacement Value
                
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
                    for(int jn=0; jn< phru; jn++)
                    {
                        for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf < this->Dimension(); jdf++)
                        {
                            ek(2*in+idf,2*jn+jdf) += bc.Val1()(idf,jdf)*data.normal[idf]*data.normal[jdf]*phiu(in,0)*phiu(jn,0)*weight;
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
                            ek(2*in+idf,2*jn+jdf) += bc.Val1()(idf,jdf)*data.normal[idf]*data.normal[jdf]*phiu(in,0)*phiu(jn,0)*weight;
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
            PZError << "TPZMatElasticity2D::ContributeBC error - Wrong boundary condition type" << std::endl;
            DebugStop();
        }
            break;
    }

}

void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    /*
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
     */
    
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<REAL,3> sol_u = data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    REAL uy = sol_u[1];
    
    int phru = phiu.Rows();
    short in,jn;
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
    out << "Plane Problem (fPlaneStress = 0, for Plane Strain conditions) " << fPlaneStress << std::endl;
    out << "Properties for elasticity: \n";
    out << "\t Young modulus   = "											<< fE		<< std::endl;
    out << "\t Poisson Ratio   = "											<< fnu		<< std::endl;
    out << "\t First Lamé Parameter   = "									<< flambda	<< std::endl;
    out << "\t Second Lamé Parameter   = "									<< fmu		<< std::endl;
    out << "\t Body force vector B {X-direction, Y-direction}   = "			<< ff[0] << ' ' << ff[1]   << std::endl;
    out << "\t fPreStressXX   = "			<< fPreStressXX << std::endl;
    out << "\t fPreStressXY   = "			<< fPreStressXY << std::endl;
    out << "\t fPreStressYY   = "			<< fPreStressYY << std::endl;
    out << "\t fPreStressZZ   = "			<< fPreStressZZ << std::endl;
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
    PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
    return -1;
    
    return TPZMaterial::VariableIndex(name);
}

/**
 * Save the element data to a stream
 */
void TPZMatElasticity2D::Write(TPZStream &buf, int withclassid)
{
    TPZMaterial::Write(buf,withclassid);
    buf.Write(&fE);
    buf.Write(&fnu);
    buf.Write(&flambda);
    buf.Write(&fmu);
    TPZSaveable::WriteObjects(buf, ff);
    buf.Write(&fPreStressXX);
    buf.Write(&fPreStressXY);
    buf.Write(&fPreStressYY);
    buf.Write(&fPreStressZZ);
    buf.Write(&fPlaneStress);
    
}

/**
 * Read the element data from a stream
 */
void TPZMatElasticity2D::Read(TPZStream &buf, void *context)
{
    TPZMaterial::Read(buf,context);
    buf.Read(&fE);
    buf.Read(&fnu);
    buf.Read(&flambda);
    buf.Read(&fmu);
    TPZSaveable::ReadObjects(buf, ff);
    buf.Read(&fPreStressXX);
    buf.Read(&fPreStressXY);
    buf.Read(&fPreStressYY);
    buf.Read(&fPreStressZZ);
    buf.Read(&fPlaneStress);
    
}

int TPZMatElasticity2D::NSolutionVariables(int var){
    if(var == 1)	return 3;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    
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


