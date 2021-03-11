/*
 *  pzporoelastic2d.cpp
 *  PZ
 *
 *  Created by Agnaldo on 11/28/11.
 *	Modified and Improved by Omar Duran on 11/28/11.
 *  Copyright 2012 L@bMeC. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include "pzporoelastic2d.h"
#include "pzelasmat.h" 
#include "pzbndcond.h"
#include "pzaxestools.h"

// Numerical Approximation of Reservoir Fault Stability with Linear Poroelasticty
// Duran O.
// Defining State Variables Formulation Displacement/Pressure u/p
const int StateVarUx = 0;
const int StateVarUy = 1;
const int StateVarPressure = 2;

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic.data"));
#endif

TPZPoroElastic2d::EState TPZPoroElastic2d::gState = ECurrentState;

TPZPoroElastic2d::TPZPoroElastic2d():TPZMaterial(), ff(0), fnu(0.), falpha(0.), fk(0.), fvisc(0.), fPlaneStress(0) {
	fE = 0.;
	fDim = 2;
	fmatId = 0;
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
	fPlaneStress = 1.;
	ftheta = 1.0;
	
}

TPZPoroElastic2d::TPZPoroElastic2d(int matid, int dim):TPZMaterial(matid), ff(0), fnu(0.), falpha(0.), fk(0.), fvisc(0.),fPlaneStress(0) {
	fE = 0.;
	fDim = dim;
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
	fPlaneStress = 1;
	fmatId = matid;
	ftheta = 1.0;	
}

TPZPoroElastic2d::~TPZPoroElastic2d(){
}

int TPZPoroElastic2d::ClassId() const{
    return Hash("TPZPoroElastic2d") ^ TPZMaterial::ClassId() << 1;
}




void TPZPoroElastic2d::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef) {
	
	// The finite element formulation at element level is implemented in this method
	// Numerical Approximation of Reservoir Fault Stability with Linear Poroelasticty
	// Duran O.	Devloo P.
	
	int nref =  datavec.size();
	
	//	for 3D implementation remove this
	if (nref != 2 ) {
		std::cout << " Error.!! datavec size must be 2, this is a two dimension problem, right!\n";
		DebugStop();
	}
	
	//	Matrix definitions
	//	The total nodes for displacements is equal phiu.Rows()
	//	The total nodes for displacements is equal phip.Rows()	
	
	//	Setting the size of first block of first problem. 
	//	Linear Elastic problem
	int phcu, phru, dphcu, dphru;
	
	TPZFMatrix<REAL>	&phiu	=	datavec[0].phi;
	TPZFMatrix<REAL>	&dphiu	=	datavec[0].dphix;
	TPZFMatrix<REAL>	&axes	=	datavec[0].axes;
	phru	=	phiu.Rows();
	phcu	=	phiu.Cols();
	dphcu	=	dphiu.Cols();
	dphru	=	dphiu.Rows();
	
	// Data validation
	if(phcu != 1 || dphru != 2 || phru != dphcu) 
	{
		PZError 
		<<	"\n Inconsistent Data for Elasticity system : \n" 
		<<	"	phi.Cols() = "		<<	phiu.Cols() 
		<<	"	dphi.Cols() = "	<<	dphiu.Cols() 
		<<	"	phi.Rows = " 		<<	phiu.Rows() 
		<<	"	dphi.Rows = " 	<<	dphiu.Rows() 
		<<	"	\n";
		return;
	}
	
	//	Setting the size of first block of first problem. 
	//	Diffusion problem
	
	TPZFMatrix<REAL>	&phip	=	datavec[1].phi;
	TPZFMatrix<REAL>	&dphip	=	datavec[1].dphix;
	TPZFMatrix<REAL>	du(2,2);
	int phrp = phip.Rows();	
	int efr, efc, ekr, ekc;  
	efr	= ef.Rows();
	efc	= ef.Cols();
	ekr	= ek.Rows();
	ekc	= ek.Rows();
	
	// Data validation
	if(ekr != (2*phru + phrp) || ekc != (2*phru + phrp) || efr != (2*phru + phrp) || efc != 1)
	{
		PZError 
		<<	"	\n Inconsistent input data : \n" 
		<<	"	\nek.Rows() = "	<< ek.Rows() 
		<<	"	ek.Cols() = " << ek.Cols() 
		<<	"	\nef.Rows() = " << ef.Rows() 
		<<	"	ef.Cols() = " << ef.Cols() 
		<<	"	\n";
		return;
	}
	
	//	Contribution for Stiffness Matrix
	if(gState == ECurrentState)
	{	 

		TPZVec <double> StateVariable(3,0.0);		
	
		//	Elastic equation
		//	Linear strain operator
		//	Ke Matrix
		for(int in = 0; in < phru; in++ )
		{
			//	Derivative calculations for Ux
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			//	Derivative calculations for Uy			
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			//	Fu Vector Force right hand term
			ef(2*in, 0)		+= weight*ff[0]*phiu(in, 0) + weight*StateVariable[StateVarUx]; 
			ef(2*in+1, 0)	+= weight*ff[1]*phiu(in, 0) + weight*StateVariable[StateVarUy];
			
			for(int jn = 0; jn < phru; jn++)
			{
				//	Derivative calculations for Vx				
				du(0,1) = dphiu(0,jn)*axes(0,0)+dphiu(1,jn)*axes(1,0);
				//	Derivative calculations for Vy				
				du(1,1) = dphiu(0,jn)*axes(0,1)+dphiu(1,jn)*axes(1,1);
				
				if (fPlaneStress == 1)
				{
					/* Plain stress state */ 
					ek(2*in,2*jn)		+= weight*((4*(fmu)*(flambda+fmu)/(flambda+2*fmu))*du(0,0)*du(0,1)		+ (2*fmu)*du(1,0)*du(1,1));
					
					ek(2*in,2*jn+1)		+= weight*((2*(fmu)*(flambda)/(flambda+2*fmu))*du(0,0)*du(1,1)			+ (2*fmu)*du(1,0)*du(0,1));
					
					ek(2*in+1,2*jn)		+= weight*((2*(fmu)*(flambda)/(flambda+2*fmu))*du(1,0)*du(0,1)			+ (2*fmu)*du(0,0)*du(1,1));
					
					ek(2*in+1,2*jn+1)	+= weight*((4*(fmu)*(flambda+fmu)/(flambda+2*fmu))*du(1,0)*du(1,1)		+ (2*fmu)*du(0,0)*du(0,1));					
				}
				else
				{
					/* Plain Strain State */  
					ek(2*in,2*jn)		+= weight*	((flambda + 2*fmu)*du(0,0)*du(0,1)	+ (fmu)*du(1,0)*du(1,1));
					
					ek(2*in,2*jn+1)		+= weight*	(flambda*du(0,0)*du(1,1)			+ (fmu)*du(1,0)*du(0,1));
					
					ek(2*in+1,2*jn)		+= weight*	(flambda*du(1,0)*du(0,1)			+ (fmu)*du(0,0)*du(1,1));
					
					ek(2*in+1,2*jn+1)	+= weight*	((flambda + 2*fmu)*du(1,0)*du(1,1)	+ (fmu)*du(0,0)*du(0,1));					
					
				}
			}
		}		

		//	Matrix Qc
		//	Coupling matrix 
		for(int in = 0; in < phru; in++ )
		{
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			for(int jn = 0; jn < phrp; jn++)
			{
				ek(2*in,2*phru+jn) += (-1.)*falpha*weight*(phip(jn,0)*du(0,0));		
				ek(2*in+1,2*phru+jn) += (-1.)*falpha*weight*(phip(jn,0)*du(1,0));							
			}
		}
		
		//	Matrix QcˆT		
		//	Coupling matrix transpose
		for(int in = 0; in < phru; in++ )
		{
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			for(int jn = 0; jn < phrp; jn++)
			{
				ek(2*phru+jn,2*in) += (-1.)*falpha*weight*(phip(jn,0)*du(0,0));		
				ek(2*phru+jn,2*in+1) += (-1.)*falpha*weight*(phip(jn,0)*du(1,0));
				
			}
		}
		
		
		//	Diffusion Equation
		//	Compresibility and Permeability  matrix
		const REAL DeltaT = fTimeStep;
		for(int in = 0; in < phrp; in++)
		{
			//	Fp Vector Mass right hand term			
			ef(in+2*phru, 0) += weight*DeltaT*StateVariable[StateVarPressure];			
			for(int jn = 0; jn < phrp; jn++)
			{	
				//	S Matrix				
				ek(in+2*phru, jn+2*phru) += (-1.0)*weight*fSe*phip(in,0)*phip(jn,0);	

				//	H Matrix				
				for(int kd=0; kd<fDim; kd++) 
				{
					ek(in+2*phru, jn+2*phru) +=(-1.0)*weight*(fk/fvisc)*DeltaT*ftheta*dphip(kd,in)*dphip(kd,jn);
				}
			}
		}
	}
	
	
	//	Mass Matrix Calculations
	if(gState == ELastState)
	{				
		TPZFMatrix<STATE> ekk(ek.Rows(),ek.Rows(),0.0);
		const REAL DeltaT = fTimeStep;	
		
		//	Matrix QcˆT		
		//	Coupling matrix transpose
		for(int in = 0; in < phru; in++ )
		{
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			for(int jn = 0; jn < phrp; jn++)
			{
				ek(2*phru+jn,2*in) += (-1.)*falpha*weight*(phip(jn,0)*du(0,0));		
				ek(2*phru+jn,2*in+1) += (-1.)*falpha*weight*(phip(jn,0)*du(1,0));
				
			}
		}
		
		//	Diffusion Equation
		//	Compresibility and Permeability  matrix
		for(int in = 0; in < phrp; in++)
		{	
			for(int jn = 0; jn < phrp; jn++)
			{	
				//	S Matrix
				ek(in+2*phru, jn+2*phru) += (-1.0)*weight*fSe*phip(in,0)*phip(jn,0);
				
				//	H Matrix				
				for(int kd=0; kd<fDim; kd++) 
				{
					ek(in+2*phru, jn+2*phru) +=(1.0)*weight*(fk/fvisc)*DeltaT*(1-ftheta)*dphip(kd,in)*dphip(kd,jn);
				}
			}
		}
	}
	
}

//	Here is implemented the contribution for boundary conditions
void TPZPoroElastic2d::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) 
{
	
	//	This state calculate the contribution for the Mass Matrix
	if(gState == ELastState)
	{
		return;
	}
			
	int nref =  datavec.size();
	if (nref != 2) {
		std::cout << " Error:: datavec size must to be equal to 2 \n";
		DebugStop();
	}
	if (bc.Val2().Rows() != 3) 
	{
		std::cout << " Error:: This material need boundary conditions for ux, uy and p (pore pressure).\n";
		std::cout << " give me one matrix with this form Val2(3,1).\n";
		DebugStop();
	}
	
	if (bc.Val1().Rows() != 3) 
	{
		std::cout << " Error:: This material need boundary conditions for ux, uy and p (pore pressure).\n";
		DebugStop();
	}
	

	TPZFMatrix<REAL>  &phiu = datavec[0].phi;
	TPZFMatrix<REAL>  &phip = datavec[1].phi;
	
	int phru = phiu.Rows();
	int phrp = phip.Rows();
	short in,jn;
	STATE v2[3];
	v2[0] = bc.Val2()(0,0);	//	Ux displacement
	v2[1] = bc.Val2()(1,0);	//	Uy displacement
	v2[2] = bc.Val2()(2,0);	//	Pressure
		
//	Here each digit represent an individual boundary condition corresponding to each state variable.
//	0 means Dirichlet condition
//	1 means Neumman conditino
//	2 means Mixed condition
//	Hence 000 means Dirichlet for all state variables or 010 means Dirichlet for Ux,P and Neumman for Uy, and so forth.

	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
	const REAL DeltT = fTimeStep;	
	switch (bc.Type()) 
	{
		case 0 :
		{
			//	Dirichlet condition for each state variable
			//	Elasticity Equation
			for(in = 0 ; in < phru; in++) 
			{
				//	Contribution for load Vector
				ef(2*in,0)		+= BIGNUMBER*v2[0]*phiu(in,0)*weight;	// X displacement Value      
				ef(2*in+1,0)	+= BIGNUMBER*v2[1]*phiu(in,0)*weight;	// y displacement Value 
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					//	Contribution for Stiffness Matrix
					ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
					ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
				}
			}
						
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Contribution for load Vector				
				ef(in+2*phru,0)		+= BIGNUMBER*v2[2]*phip(in,0)*weight;	// P Pressure Value
				
				for (jn = 0 ; jn < phrp; jn++) 
				{
					//	Contribution for Stiffness Matrix					
					ek(in+2*phru,jn+2*phru)		+= BIGNUMBER*phip(in,0)*phip(jn,0)*weight;	// P Pressure
				}
			}
			break;
		}	
		case 111 :
		{
			//	Neumann condition for each state variable
			//	Elasticity Equation
			for(in = 0 ; in <phru; in++) 
			{           
				//	Normal Tension Components on neumman boundary
				ef(2*in,0)		+= v2[0]*phiu(in,0)*weight;		//	Tnx
				ef(2*in+1,0)	+= v2[1]*phiu(in,0)*weight;		//	Tny
			}

			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Normal Flux on neumman boundary				
				ef(in+2*phru,0)	+= v2[2]*DeltT*phip(in,0)*weight;	// Qnormal
			}
			break;
		}	
		case 222 : 
		{
			//	Mixed condition for each state 			
			//	Elasticity Equation			
			for(in = 0 ; in < phru; in++) 
			{
				ef(2*in, 0)				+= v2[0] * phiu(in, 0) * weight;
				ef(2*in+1, 0)			+= v2[1] * phiu(in, 0) * weight;
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					ek(2*in,2*jn)		+= bc.Val1()(0,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn)		+= bc.Val1()(1,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn+1)	+= bc.Val1()(1,1)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in,2*jn+1)		+= bc.Val1()(0,1)*phiu(in,0)*
					phiu(jn,0)*weight;
				}
			} 
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				ef(in+2*phru, 0)			+= v2[2] * phip(in, 0)*weight;
				for (jn = 0 ; jn < phrp; jn++) 
				{
					ek(in+2*phru,jn+2*phru)	+= bc.Val1()(2,0)*phip(in,0)*
					phip(jn,0)*weight;
				}
			}	
			break;
		}
			
		case 1:
		{		
			//	Dirichlet(Ux,Uy)	Neumman(P)			
			//	Elasticity Equation
			for(in = 0 ; in < phru; in++) 
			{
				//	Contribution for load Vector
				ef(2*in,0)		+= BIGNUMBER*v2[0]*phiu(in,0)*weight;	// X displacement Value      
				ef(2*in+1,0)	+= BIGNUMBER*v2[1]*phiu(in,0)*weight;	// y displacement Value 
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					//	Contribution for Stiffness Matrix
					ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
					ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
				}
			}
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Normal Flux on neumman boundary				
				ef(in+2*phru,0)	+= v2[2]*DeltT*phip(in,0)*weight;	// Qnormal
			}
			break;
		}	
		case 11:
		{		
			//	Dirichlet(Ux)		Neumman(Uy,P)			
			//	Elasticity Equation
			for(in = 0 ; in < phru; in++) 
			{
				//	Contribution for load Vector
				ef(2*in,0)		+= BIGNUMBER*v2[0]*phiu(in,0)*weight;	// X displacement Value  
				//	Normal Tension Components on neumman boundary
				ef(2*in+1,0)	+= v2[1]*phiu(in,0)*weight;				//	Tny				

				for (jn = 0 ; jn < phru; jn++) 
				{
					//	Contribution for Stiffness Matrix
					ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
				}
			}				
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Normal Flux on neumman boundary				
				ef(in+2*phru,0)	+= v2[2]*DeltT*phip(in,0)*weight;	// Qnormal
			}
			break;
		}			
			
		case 101:
		{		
			//	Neumman(Ux)	Dirichlet(Uy)	Neumman(P)		
			//	Elasticity Equation
			for(in = 0 ; in < phru; in++) 
			{
				//	Contribution for load Vector      
				//	Normal Tension Components on neumman boundary
				ef(2*in,0)		+= v2[0]*phiu(in,0)*weight;				//	Tnx							
				ef(2*in+1,0)	+= BIGNUMBER*v2[1]*phiu(in,0)*weight;	// y displacement Value
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					//	Contribution for Stiffness Matrix
					ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
				}
			}
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Normal Flux on neumman boundary				
				ef(in+2*phru,0)	+= v2[2]*DeltT*phip(in,0)*weight;	// Qnormal
			}
			break;
		}
			
		case 110:
		{		
			//	Neumman(Ux,Uy)	Dirichlet(P)			
			//	Elasticity Equation
			for(in = 0 ; in <phru; in++) 
			{           
				//	Normal Tension Components on neumman boundary
				ef(2*in,0)		+= v2[0]*phiu(in,0)*weight;		//	Tnx
				ef(2*in+1,0)	+= v2[1]*phiu(in,0)*weight;		//	Tny
			}
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Contribution for load Vector				
				ef(in+2*phru,0)		+= BIGNUMBER*v2[2]*phip(in,0)*weight;	// P Pressure Value
				
				for (jn = 0 ; jn < phrp; jn++) 
				{
					//	Contribution for Stiffness Matrix					
					ek(in+2*phru,jn+2*phru)		+= BIGNUMBER*phip(in,0)*phip(jn,0)*weight;	// P Pressure
				}
			}
			break;
		}
			
		case 100:
		{		
			//	Neumman(Ux)		Dirichlet(Uy,P)			
			//	Elasticity Equation
			for(in = 0 ; in < phru; in++) 
			{
				//	Contribution for load Vector
				//	Normal Tension Components on neumman boundary
				ef(2*in,0)		+= v2[0]*phiu(in,0)*weight;		//	Tnx				
				ef(2*in+1,0)	+= BIGNUMBER*v2[1]*phiu(in,0)*weight;	// y displacement Value 
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					//	Contribution for Stiffness Matrix
					ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
				}
			}			

			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Contribution for load Vector				
				ef(in+2*phru,0)		+= BIGNUMBER*v2[2]*phip(in,0)*weight;	// P Pressure Value
				
				for (jn = 0 ; jn < phrp; jn++) 
				{
					//	Contribution for Stiffness Matrix					
					ek(in+2*phru,jn+2*phru)		+= BIGNUMBER*phip(in,0)*phip(jn,0)*weight;	// P Pressure
				}
			}
			break;
		}
			
		case 10:
		{		
			//	Dirichlet(Ux)		Neumman(Uy,P)			
			//	Elasticity Equation
			for(in = 0 ; in < phru; in++) 
			{
				//	Contribution for load Vector
				ef(2*in,0)		+= BIGNUMBER*v2[0]*phiu(in,0)*weight;	// X displacement Value  
				//	Normal Tension Components on neumman boundary
				ef(2*in+1,0)	+= v2[1]*phiu(in,0)*weight;				//	Tny				
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					//	Contribution for Stiffness Matrix
					ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
				}
			}				
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Contribution for load Vector				
				ef(in+2*phru,0)		+= BIGNUMBER*v2[2]*phip(in,0)*weight;	// P Pressure Value
				
				for (jn = 0 ; jn < phrp; jn++) 
				{
					//	Contribution for Stiffness Matrix					
					ek(in+2*phru,jn+2*phru)		+= BIGNUMBER*phip(in,0)*phip(jn,0)*weight;	// P Pressure
				}
			}
			break;
		}			
					
		case 220 : 
		{
			//	Mixed condition for each state 			
			//	Elasticity Equation			
			for(in = 0 ; in < phru; in++) 
			{
				ef(2*in, 0)				+= v2[0] * phiu(in, 0) * weight;
				ef(2*in+1, 0)			+= v2[1] * phiu(in, 0) * weight;
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					ek(2*in,2*jn)		+= bc.Val1()(0,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn)		+= bc.Val1()(1,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn+1)	+= bc.Val1()(1,1)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in,2*jn+1)		+= bc.Val1()(0,1)*phiu(in,0)*
					phiu(jn,0)*weight;
				}
			} 
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Contribution for load Vector				
				ef(in+2*phru,0)		+= BIGNUMBER*v2[2]*phip(in,0)*weight;	// P Pressure Value
				
				for (jn = 0 ; jn < phrp; jn++) 
				{
					//	Contribution for Stiffness Matrix					
					ek(in+2*phru,jn+2*phru)		+= BIGNUMBER*phip(in,0)*phip(jn,0)*weight;	// P Pressure
				}
			}
			break;
		}				
		case 221 : 
		{
			//	Mixed condition for each state 			
			//	Elasticity Equation			
			for(in = 0 ; in < phru; in++) 
			{
				ef(2*in, 0)				+= v2[0] * phiu(in, 0) * weight;
				ef(2*in+1, 0)			+= v2[1] * phiu(in, 0) * weight;
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					ek(2*in,2*jn)		+= bc.Val1()(0,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn)		+= bc.Val1()(1,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn+1)	+= bc.Val1()(1,1)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in,2*jn+1)		+= bc.Val1()(0,1)*phiu(in,0)*
					phiu(jn,0)*weight;
				}
			} 
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Normal Flux on neumman boundary				
				ef(in+2*phru,0)	+= v2[2]*DeltT*phip(in,0)*weight;	// Qnormal
			}
			break;
		}
		case 3220 : 
		{
			//	Mixed condition for each state 			
			//	Elasticity Equation			
			for(in = 0 ; in < phru; in++) 
			{
//				ef(2*in, 0)				+= v2[0] * phiu(in, 0) * weight * (n[0]);
//				ef(2*in+1, 0)			+= v2[1] * phiu(in, 0) * weight * (n[0]);
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					ek(2*in,2*jn)		+= bc.Val1()(0,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn)		+= bc.Val1()(1,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn+1)	+= bc.Val1()(1,1)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in,2*jn+1)		+= bc.Val1()(0,1)*phiu(in,0)*
					phiu(jn,0)*weight;
				}
			} 
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Contribution for load Vector				
				ef(in+2*phru,0)		+= BIGNUMBER*v2[2]*phip(in,0)*weight;	// P Pressure Value
				
				for (jn = 0 ; jn < phrp; jn++) 
				{
					//	Contribution for Stiffness Matrix					
					ek(in+2*phru,jn+2*phru)		+= BIGNUMBER*phip(in,0)*phip(jn,0)*weight;	// P Pressure
				}
			}
			break;
		}				
		case 3221 : 
		{
			//	Mixed condition for each state 			
			//	Elasticity Equation			
			for(in = 0 ; in < phru; in++) 
			{
		//		ef(2*in, 0)				+= v2[0] * phiu(in, 0) * weight * (n[0]);
		//		ef(2*in+1, 0)			+= v2[1] * phiu(in, 0) * weight * (n[1]);
				
				for (jn = 0 ; jn < phru; jn++) 
				{
//					ek(2*in,2*jn)		+= bc.Val1()(0,0)*phiu(in,0)*
//					phiu(jn,0)*weight * (n[0]);
//					
//					ek(2*in+1,2*jn)		+= bc.Val1()(1,0)*phiu(in,0)*
//					phiu(jn,0)*weight * (n[0]);
					
//					ek(2*in+1,2*jn+1)	+= bc.Val1()(1,1)*phiu(in,0)*
//					phiu(jn,0)*weight * (n[1]);
//					
//					ek(2*in,2*jn+1)		+= bc.Val1()(0,1)*phiu(in,0)*
//					phiu(jn,0)*weight * (n[1]);
				}
			} 
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Normal Flux on neumman boundary				
				ef(in+2*phru,0)	+= v2[2]*DeltT*phip(in,0)*weight;	// Qnormal
			}
			break;
		}
		case 3021 : 
		{
			//	Mixed condition for each state 			
			//	Elasticity Equation			
			for(in = 0 ; in < phru; in++) 
			{
//				ef(2*in+1, 0)			+= v2[1] * phiu(in, 0) * weight * (n[0]);
				
				for (jn = 0 ; jn < phru; jn++) 
				{
//					ek(2*in,2*jn)		+= bc.Val1()(0,0)*phiu(in,0)*
//					phiu(jn,0)*weight;
//					
//					ek(2*in+1,2*jn)		+= bc.Val1()(1,0)*phiu(in,0)*
//					phiu(jn,0)*weight;
					
//					ek(2*in+1,2*jn+1)	+= bc.Val1()(1,1)*phiu(in,0)*
//					phiu(jn,0)*weight;
//					
//					ek(2*in,2*jn+1)		+= bc.Val1()(0,1)*phiu(in,0)*
//					phiu(jn,0)*weight;
				}
			} 
			
			//	Dirichlet(Ux,Uy)	Neumman(P)			
			//	Elasticity Equation
			for(in = 0 ; in < phru; in++) 
			{
				//	Contribution for load Vector
				ef(2*in,0)		+= BIGNUMBER*v2[0]*phiu(in,0)*weight;	// X displacement Value       
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					//	Contribution for Stiffness Matrix
					ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
				}
			}			
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Contribution for load Vector				
				ef(in+2*phru,0)		+= BIGNUMBER*v2[2]*phip(in,0)*weight;	// P Pressure Value
				
				for (jn = 0 ; jn < phrp; jn++) 
				{
					//	Contribution for Stiffness Matrix					
					ek(in+2*phru,jn+2*phru)		+= BIGNUMBER*phip(in,0)*phip(jn,0)*weight;	// P Pressure
				}
			}
			break;
		}				
		case 3201 : 
		{
			//	Mixed condition for each state 			
			//	Elasticity Equation			
			for(in = 0 ; in < phru; in++) 
			{
//				ef(2*in, 0)				+= v2[0] * phiu(in, 0) * weight * (n[0]);
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					ek(2*in,2*jn)		+= bc.Val1()(0,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
					ek(2*in+1,2*jn)		+= bc.Val1()(1,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					
				}
			} 
			
			//	Dirichlet(Ux,Uy)	Neumman(P)			
			//	Elasticity Equation
			for(in = 0 ; in < phru; in++) 
			{
				//	Contribution for load Vector     
				ef(2*in+1,0)	+= BIGNUMBER*v2[1]*phiu(in,0)*weight;	// y displacement Value 
				
				for (jn = 0 ; jn < phru; jn++) 
				{
					//	Contribution for Stiffness Matrix
					ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
				}
			}			
			
			//	Diffusion Equation 
			for(in = 0 ; in < phrp; in++) 
			{
				//	Normal Flux on neumman boundary				
				ef(in+2*phru,0)	+= v2[2]*DeltT*phip(in,0)*weight;	// Qnormal
			}
			break;
		}							
		case 5://Normal load -> Pressure
		{
//			for(in = 0 ; in < phip.Rows(); in++)
//			{
//				ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight * (data.normal[0]);
//				ef(nstate*in+1,0) += v2[0] * phi(in,0) * weight * (data.normal[1]);
//				ef(nstate*in+2,0) += v2[0] * phi(in,0) * weight * (data.normal[2]);
//			}
		}			
			break;			
	}
	
}

void TPZPoroElastic2d::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)

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


void TPZPoroElastic2d::Print(std::ostream &out) 
	{
	out << "Material Name : " << Name() << "\n";
	out << "Plane Problem (fPlaneStress = 0, for Plane Strain conditions) " << fPlaneStress << std::endl;		
	out << "Properties for Poroelasticity: \n";
	out << "\t Young modulus   = "											<< fE		<< std::endl;
	out << "\t Undrained Young modulus   = "								<< fEu		<< std::endl;
	out << "\t Poisson Ratio   = "											<< fnu		<< std::endl;
	out << "\t Undarined Poisson Ratio   = "								<< fnuu		<< std::endl;		
	out << "\t First Lamé Parameter   = "									<< flambda	<< std::endl;
	out << "\t Second Lamé Parameter   = "									<< fmu		<< std::endl;		
	out << "\t Undrained First Lamé Parameter   = "							<< flambdau	<< std::endl;
	out << "\t Biot coefficient   = "										<< falpha	<< std::endl;		
	out << "\t Body force vector B {X-direction, Y-direction}   = "			<< ff[0] << ' ' << ff[1]   << std::endl;		
	out << "Properties for Diffusion: \n";		
	out << "\t Permeability   = "											<< fk		<< std::endl;		
	out << "\t Fluid Viscosity   = "										<< fvisc	<< std::endl;
	out << "\t Constrained specific storage at constant strain Se = "		<< fSe		<< std::endl;		
	out << "Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
		
}

/** Returns the variable index associated with the name */
int TPZPoroElastic2d::VariableIndex(const std::string &name)
{
	//	Elasticity Variables
	if(!strcmp("Displacement",name.c_str()))				return	1;
	if(!strcmp("SolidPressure",name.c_str()))				return	2;
	if(!strcmp("SigmaX",name.c_str()))						return	3;
	if(!strcmp("SigmaY",name.c_str()))						return	4;
	if(!strcmp("TauXY",name.c_str()))						return	5;
	if(!strcmp("TotalSigmaX",name.c_str()))					return	22;
	if(!strcmp("TotalSigmaY",name.c_str()))					return	23;	
	if(!strcmp("DisplacementX",name.c_str()))				return	8;
	if(!strcmp("DisplacementY",name.c_str()))				return	9;
	
	//	Diffusion Variables
	if(!strcmp("FluidPressure",name.c_str()))				return	6;
	if(!strcmp("FluxVector",name.c_str()))					return	7;
	
	//	Exact Solution Variables
	if(!strcmp("EDisplacementX",name.c_str()))				return	11;
	if(!strcmp("EDisplacementY",name.c_str()))				return	12;
	if(!strcmp("ESIGX",name.c_str()))						return	15;
	if(!strcmp("ESIGY",name.c_str()))						return	16;
	if(!strcmp("ETAUXY",name.c_str()))						return	17;
	if(!strcmp("EPressure",name.c_str()))					return	10;
	if(!strcmp("EFluxVector",name.c_str()))					return	18;

	//	Total mass content Defined in: Fluid mass sources and point forces in linear elastic diffusive solids (John Rudnicki 1986)
	if(!strcmp("TotalMassContent",name.c_str()))			return	13;	
		
	// Stress Arching Ratio
	if(!strcmp("SArchX",name.c_str()))						return	19;
	if(!strcmp("SArchY",name.c_str()))						return	20;
	if(!strcmp("SArchXY",name.c_str()))						return	21;	
	
	// Reservoir Stress Path
	if(!strcmp("GammaX",name.c_str()))						return	24;
	if(!strcmp("GammaY",name.c_str()))						return	25;
	if(!strcmp("GammaXY",name.c_str()))						return	26;
	
	// Exact Reservoir Stress Path
	if(!strcmp("EGammaX",name.c_str()))						return	27;
	if(!strcmp("EGammaY",name.c_str()))						return	28;
	if(!strcmp("EGammaXY",name.c_str()))					return	29;
	
	// Principal Space
	if(!strcmp("Sigma1",name.c_str()))						return	30;
	if(!strcmp("Sigma2",name.c_str()))						return	31;
	if(!strcmp("Theta1",name.c_str()))						return	32;
	if(!strcmp("Theta2",name.c_str()))						return	33;	
	if(!strcmp("I1",name.c_str()))							return	34;
	if(!strcmp("J2",name.c_str()))							return	35;		
		
	return TPZMaterial::VariableIndex(name);
}

int TPZPoroElastic2d::NSolutionVariables(int var){
	if(var == 1)	return 3;
	if(var == 2)	return 1;
	if(var == 3)	return 1;
	if(var == 4)	return 1;
	if(var == 5)	return 1;
	if(var == 6)	return 1;
	if(var == 7)	return 3;
	if(var == 8)	return 1;
	if(var == 9)	return 1;
	if(var == 10)	return 1;
	if(var == 11)	return 1;
	if(var == 12)	return 1;
	if(var == 13)	return 1;
	if(var == 14)	return 3;
	if(var == 15)	return 1;	
	if(var == 16)	return 1;	
	if(var == 17)	return 1;
	if(var == 18)	return 3;
	if(var == 19)	return 1;	
	if(var == 20)	return 1;	
	if(var == 21)	return 1;
	if(var == 22)	return 1;	
	if(var == 23)	return 1;
	if(var == 24)	return 1;
	if(var == 25)	return 1;	
	if(var == 26)	return 1;
	if(var == 29)	return 1;
	if(var == 28)	return 1;	
	if(var == 27)	return 1;
	if(var == 30)	return 1;	
	if(var == 31)	return 1;
	if(var == 32)	return 1;
	if(var == 33)	return 1;	
	if(var == 34)	return 1;
	if(var == 35)	return 1;	
	return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZPoroElastic2d::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZManVector<STATE,3> SolU, SolP;
	TPZFNMatrix <6,STATE> DSolU, DSolP;
	TPZFNMatrix <9> axesU, axesP;
	
	TPZVec<REAL> ptx(3);
    TPZVec<STATE> solExata(3);
	TPZFMatrix<STATE> flux(5,1);
    
    if (datavec[0].sol.size() != 1) {
        DebugStop();
    }
	
	SolU=datavec[0].sol[0];
	DSolU=datavec[0].dsol[0];
	axesU=datavec[0].axes;
	SolP=datavec[1].sol[0];
	DSolP=datavec[1].dsol[0];
	axesP=datavec[1].axes;
	
	
	//	Fill data for Non-permeable region without Diffusive domain
	if (SolP.size() == 0) {
		SolP.Resize(1, 0.0);
		
		DSolP.Resize(2,3);		
		datavec[1].dsol[0].Resize(2,3);
		DSolP=datavec[1].dsol[0];			
		
		datavec[1].x.Resize(3,0.0);
		datavec[1].x = datavec[0].x;
		
		axesP.Resize(3,3);
		datavec[1].axes.Resize(3,3);
		axesP=datavec[1].axes;
	}
	
	//	Displacements
	if(var == 1){
		Solout[0] = SolU[0];
		Solout[1] = SolU[1];
		Solout[2] = 0.0;
		return;
	}
	
	//	x-displacement
	if(var == 8){
		Solout[0] = SolU[0];
		return;
	}
	
	//	y-displacement
	if(var == 9){
		Solout[0] = SolU[1];
		return;
	}
	
	//	Pore pressure excess
	if(var == 6) {
		Solout[0] = SolP[0];
		return;
	}
	
	//	Darcy's velocity
	if (var == 7)
	{ 
		int id;
		TPZManVector<STATE> dsolp(2,0);
		dsolp[0] = datavec[1].dsol[0](0,0)*datavec[1].axes(0,0)+datavec[1].dsol[0](1,0)*datavec[1].axes(1,0);
		dsolp[1] = datavec[1].dsol[0](0,0)*datavec[1].axes(0,1)+datavec[1].dsol[0](1,0)*datavec[1].axes(1,1);			
		for(id=0 ; id<fDim; id++) 
		{
			Solout[id] = -1. * this->fK * dsolp[id];
		}
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
	REAL TotalMass;
						
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
		SigX = ((flambda + 2*fmu)*(epsx) + (flambda)*epsy - falpha*SolP[0]);
		SigY = ((flambda + 2*fmu)*(epsy) + (flambda)*epsx - falpha*SolP[0]);	
		SigZ = flambda*divu - falpha*SolP[0];
		Tau = 2.0*fmu*epsxy;		
	}
	
	REAL R = sqrt( ((SigX - SigY)/2)*((SigX - SigY)/2) + Tau*Tau);
	REAL C = (SigX + SigY)/2;
	REAL Sigma1 = C + R;
	REAL Sigma2 = C - R;
	
	// Sigma1 is the largest stress regadless of sign
	if (fabs(Sigma2) > fabs(Sigma1))
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

	//	Shear Stress	
	if(var == 5) {
		Solout[0] = Tau;
		return;
	}
		
	//	Total Stress x-direction	
	if(var ==22) {
		Solout[0] = SigX + falpha*SolP[0];
		return;
	}
	
	//	Total Stress y-direction	
	if(var == 23) {
		Solout[0] = SigY + falpha*SolP[0];
		return;
	}	
	
	// Total mass fluid content 
	if(var == 13) {
		TotalMass = falpha*divu + fSe*SolP[0];
		Solout[0] = TotalMass;
		return;
	}
	
	//	Exact Displacements
	if(var == 14){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		Solout[0] = solExata[1];
		Solout[1] = solExata[2];
		Solout[2] = 0.;
		return;
	}
		
	//	Exact x-displacement 
	if(var == 11){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		Solout[0] = solExata[1];
		return;
	}
	
	//	Exact y-displacement
	if(var == 12){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		Solout[0] = solExata[2];
		return;
	}	

	//	Exact solution for pressure excess
	if(var == 10){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}	
	
	//	Exact Darcy's velocity
	if(var == 18){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		Solout[0] = flux(3,0);
		Solout[1] = flux(4,0);
		Solout[2] = 0.0;		
		return;
	}		
	
	//	Exact Effective Stress x-direction
	if(var == 15){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		Solout[0] = flux(0,0);
		return;
	}
	
	//	Exact Effective Stress y-direction
	if(var == 16){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		Solout[0] = flux(1,0);
		return;
	}
	
	//	Exact Shear Stress	
	if(var == 17){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		Solout[0] = flux(2,0);
		return;
	}
	
	//	Stress Arching ratio X
	if(var == 19){
		Solout[0] = (SigX+falpha*SolP[0])/((2*fmu)/(flambda+2*fmu));
		return;
	}
	
	//	Stress Arching ratio Y
	if(var == 20){
		Solout[0] = (SigY+falpha*SolP[0])/((2*fmu)/(flambda+2*fmu));
		return;
	}
	
	//	Stress Arching ratio XY
	if(var == 21){
		Solout[0] = (Tau+falpha*SolP[0])/((2*fmu)/(flambda+2*fmu));
		return;
	}	
	
	
	//	Reservoir Stress Path Gamma ratio X
	if(var == 24){
		if (SolP[0] <= 1.0e-10 ) 
		{
			Solout[0] = 0;
		}else
		{
			Solout[0] = (SigX/SolP[0]);		
		}			
		return;
	}
	
	//	Reservoir Stress Path Gamma Y
	if(var == 25){
		if (SolP[0] <= 1.0e-10 ) 
		{
			Solout[0] = 0;
		}else
		{
			Solout[0] = (SigY/SolP[0]);		
		}			
		return;
	}
	
	//	Reservoir Stress Path Gamma XY
	if(var == 26){
		if (SolP[0] <= 1.0e-10 ) 
		{
			Solout[0] = 0;
		}else
		{
			Solout[0] = (Tau/SolP[0]);		
		}			
		return;
	}
	
	//	Exact Reservoir Stress Path Gamma ratio X
	if(var == 27){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);	
		if (solExata[0] <= 1.0e-10 ) 
		{
			Solout[0] = 0;
		}else
		{
			Solout[0] = (flux(0,0)/solExata[0]);			
		}
		return;
	}
	
	//	Exact Reservoir Stress Path Gamma Y
	if(var == 28){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		if (solExata[0] <= 1.0e-10 ) 
		{
			Solout[0] = 0;
		}else
		{
			Solout[0] = (flux(1,0)/solExata[0]);		
		}		
		return;
	}
	
	//	Exact Reservoir Stress Path Gamma XY
	if(var == 29){
		fTimedependentFunctionExact->Execute(datavec[1].x, fTimeValue, solExata,flux);
		if (solExata[0] <= 1.0e-10 ) 
		{
			Solout[0] = 0;
		}else
		{
			Solout[0] = (flux(2,0)/solExata[0]);		
		}	
		return;
	}	
	
	//	Principal Stress
	if(var == 30){
		Solout[0]=Sigma1;
		return;
	}
	
	//	Secondary Stress
	if(var == 31){
		Solout[0]=Sigma2;
		return;
	}	
	
	//	Fist Invariant
	if(var == 34){
		REAL I1 = SigX+SigY;
		Solout[0]=I1;
		return;
	}
	
	//	Desviatoric Second Invariant
	if(var == 35){
		REAL J2 = ((SigX + SigY)*(SigX + SigY) - (3*(-SigX*SigX - SigY*SigY + (SigX + SigY)*(SigX + SigY) - 2*Tau*Tau))/2.)/2.;
		Solout[0]=J2;
		return;
	}
	
}

void TPZPoroElastic2d::ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZPoroElastic2d::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	DebugStop();
}

