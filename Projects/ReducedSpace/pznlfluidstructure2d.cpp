//
//  pznlfluidstructure2d.cpp
//  PZ
//
//  Created by Agnaldo Farias on 9/17/12.
//
//

#include "pznlfluidstructure2d.h"


#include <iostream>
#include <string>

#include "pzelasmat.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzintel.h"

#include "pzlog.h"
#ifdef LOG4CXX
static PZLogger logdata("pz.material.elastpressure");
#endif


TPZNLFluidStructure2d::EState TPZNLFluidStructure2d::gState = ECurrentState;

TPZNLFluidStructure2d::TPZNLFluidStructure2d() : TPZMaterial()
{
	fmatId = 0;
	fE = 0.;
	fPoiss = 0.;
	fVisc = 0.;
	fDim = 2;
	fPlaneStress = 1.;
	
	ff.resize(fDim);
	ff[0] = 0.;
	ff[1] = 0.;
}

TPZNLFluidStructure2d::TPZNLFluidStructure2d(int matid, int dim, REAL young, REAL poiss, REAL visc): TPZMaterial(matid)
{
	fmatId = matid;
	fE = young;
	fPoiss = poiss;
	fVisc = visc;
	fDim = dim;
	fPlaneStress = 1.;
	
	ff.resize(fDim);
	ff[0] = 0.;
	ff[1] = 0.;
}

TPZNLFluidStructure2d::~TPZNLFluidStructure2d()
{
}

void TPZNLFluidStructure2d::Print(std::ostream &out)
{
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "\t E   = " << fE << std::endl;
	out << "\t nu   = " << fPoiss << std::endl;
	out << "\t Forcing function F   = " << ff[0] << ' ' << ff[1]   << std::endl;
	out << "altura da fratura fHf "<< globFractInputData.Hf() << std::endl;
	out << "Viscosidade do fluido fvisc "<< fVisc << std::endl;
	out << "Carter fCl " << globFractInputData.Cl() << std::endl;
	out << "Pressao estatica fPe " << globFractInputData.Pe() << std::endl;
	out << "Pressao de referencia (Carter) fPref " << globFractInputData.Pref() << std::endl;
	out << "Spurt loss fvsp " << globFractInputData.vsp() << std::endl;
	out << "2D problem " << fPlaneStress << std::endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}


void TPZNLFluidStructure2d::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
	if(gState == ELastState)
	{
		return; 
	}
	
	int nref =  datavec.size();
	if (nref != 2) {
		std::cout << " Error.!! the size of datavec is equal to two\n";
    std::cout << " datavec[0]->elasticity and datavec[1]->pressure\n";
		DebugStop();
	}
	
	TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
	if(shapetype!=datavec[0].EVecShape && datavec[0].phi.Cols()!=0)
	{
		std::cout << " The space to elasticity problem must be reduced space.\n";
		DebugStop();
	}
	
	//Calculate the matrix contribution for elastic problem.
	TPZFMatrix<REAL> &dphi_u = datavec[0].dphix;
	TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	TPZFMatrix<REAL> &axes=datavec[0].axes;
	
	TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	TPZFMatrix<REAL> &dsol_u = datavec[0].dsol[0];
	
	int phcu = phi_u.Cols();
	int efcu = ef.Cols();
	
	if(fForcingFunction)
	{// phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE> res(3);
		fForcingFunction->Execute(datavec[0].x,res);
		ff[0] = res[0];
		ff[1] = res[1];
		ff[2] = res[2];
	}
	
	TPZFNMatrix<4,STATE> dphix_i(2,1),dphiy_i(2,1), dphix_j(2,1),dphiy_j(2,1), dsolx_j(2,1), dsoly_j(2,1);
	/*
	 * Plain strain materials values
	 */
	REAL young = this->fE;
	REAL poisson = this->fPoiss;
	
	REAL fEover1MinNu2 = young/(1-poisson*poisson);  ///4G(lamb+G)/(lamb+2G)
	REAL fEover21PlusNu = young/(2.*(1.+poisson));/*fE/(2.*(1+fnu));*/ ///2G=2mi
	REAL nu1 = 1 - poisson;//(1-nu)
	REAL nu2 = (1.-2.*poisson)/2.;
	REAL F = young/((1.+poisson)*(1.-2.*poisson));
	
	for(int in = 0; in < phcu; in++)
	{
		dphix_i(0,0) = dphi_u(0,in)*axes(0,0)+dphi_u(1,in)*axes(1,0);//dphix/dx
		dphix_i(1,0) = dphi_u(0,in)*axes(0,1)+dphi_u(1,in)*axes(1,1);//dphix/dy
		dphiy_i(0,0) = dphi_u(2,in)*axes(0,0)+dphi_u(3,in)*axes(1,0);//dphiy/dx
		dphiy_i(1,0) = dphi_u(2,in)*axes(0,1)+dphi_u(3,in)*axes(1,1);//dphiy/dy
		
		//Residuo
		for (int col = 0; col < efcu; col++)
		{
			//termo f*u
			ef(in,col) += weight * (ff[0]*phi_u(0, in) + ff[1]*phi_u(1, in));
			
			//termo prestress*gradu
			ef(in,col) += (-1.) * weight * ( dphix_i(0,0) * globFractInputData.PreStressXX() +
																			dphix_i(1,0) * globFractInputData.PreStressXY() +
																			dphiy_i(0,0) * globFractInputData.PreStressXY() +
																			dphiy_i(1,0) * globFractInputData.PreStressYY() );
			
			//termos da matriz
			dsolx_j(0,0) = dsol_u(0,0)*axes(0,0) + dsol_u(1,0)*axes(1,0);
			dsolx_j(1,0) = dsol_u(0,0)*axes(0,1) + dsol_u(1,0)*axes(1,1);
			dsoly_j(0,0) = dsol_u(0,1)*axes(0,0) + dsol_u(1,1)*axes(1,0);
			dsoly_j(1,0) = dsol_u(0,1)*axes(0,1) + dsol_u(1,1)*axes(1,1);
			
			
			if (fPlaneStress != 1)
			{/* Plain Strain State */
				ef(in,col) += (-1.) * weight * (nu1*dphix_i(0,0)*dsolx_j(0,0) + nu2*dphix_i(1,0)*dsolx_j(1,0) +
																				
																				poisson*dphix_i(0,0)*dsoly_j(1,0) + nu2*dphix_i(1,0)*dsoly_j(0,0) +
																				
																				poisson*dphiy_i(1,0)*dsolx_j(0,0) + nu2*dphiy_i(0,0)*dsolx_j(1,0) +
																				
																				nu1*dphiy_i(1,0)*dsoly_j(1,0) + nu2*dphiy_i(0,0)*dsoly_j(0,0)) * F;
			}
			else
			{/* Plain stress state */
				ef(in,col) += (-1.) * weight * (fEover1MinNu2*dphix_i(0,0)*dsolx_j(0,0) + fEover21PlusNu*dphix_i(1,0)*dsolx_j(1,0) +
																				
																				fEover1MinNu2*poisson*dphix_i(0,0)*dsoly_j(1,0) + fEover21PlusNu*dphix_i(1,0)*dsoly_j(0,0) +
																				
																				fEover1MinNu2*poisson*dphiy_i(1,0)*dsolx_j(0,0) + fEover21PlusNu*dphiy_i(0,0)*dsolx_j(1,0) +
																				
																				fEover1MinNu2*dphiy_i(1,0)*dsoly_j(1,0) + fEover21PlusNu*dphiy_i(0,0)*dsoly_j(0,0));
			}
		}//fim para o residuo
		
		//Matriz Tangente (Jacobiana)
		for( int jn = 0; jn < phcu; jn++ )
		{
			dphix_j(0,0) = dphi_u(0,jn)*axes(0,0) + dphi_u(1,jn)*axes(1,0);
			dphix_j(1,0) = dphi_u(0,jn)*axes(0,1) + dphi_u(1,jn)*axes(1,1);
			dphiy_j(0,0) = dphi_u(2,jn)*axes(0,0) + dphi_u(3,jn)*axes(1,0);
			dphiy_j(1,0) = dphi_u(2,jn)*axes(0,1) + dphi_u(3,jn)*axes(1,1);
			
			
			if (fPlaneStress != 1)
			{/* Plain Strain State */
				ek(in,jn) += weight*(nu1*dphix_i(0,0)*dphix_j(0,0) + nu2*dphix_i(1,0)*dphix_j(1,0) +
														 
														 poisson*dphix_i(0,0)*dphiy_j(1,0) + nu2*dphix_i(1,0)*dphiy_j(0,0) +
														 
														 poisson*dphiy_i(1,0)*dphix_j(0,0) + nu2*dphiy_i(0,0)*dphix_j(1,0) +
														 
														 nu1*dphiy_i(1,0)*dphiy_j(1,0) + nu2*dphiy_i(0,0)*dphiy_j(0,0))*F;
			}
			else
			{/* Plain stress state */
				ek(in,jn) += weight*(fEover1MinNu2*dphix_i(0,0)*dphix_j(0,0) + fEover21PlusNu*dphix_i(1,0)*dphix_j(1,0) +
														 
														 fEover1MinNu2*dphix_i(0,0)*dphiy_j(1,0) + fEover21PlusNu*dphix_i(1,0)*dphiy_j(0,0) +
														 
														 fEover1MinNu2*dphiy_i(1,0)*dphix_j(0,0) + fEover21PlusNu*dphiy_i(0,0)*dphix_j(1,0) +
														 
														 fEover1MinNu2*dphiy_i(1,0)*dphiy_j(1,0) + fEover21PlusNu*dphiy_i(0,0)*dphiy_j(0,0));
			}
		}//fim para a tangente
	}
	
	// Calculate the matrix contribution for pressure
	ContributePressure(datavec, weight, ek, ef);
}

void TPZNLFluidStructure2d::ContributePressure(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
	if(!datavec[1].phi) return;
	
	TPZFMatrix<REAL>  & phi_p = datavec[1].phi;
	TPZFMatrix<REAL>  & dphi_p = datavec[1].dphix;
	TPZManVector<REAL,3> sol_p = datavec[1].sol[0];
	TPZFMatrix<REAL> & dsol_p = datavec[1].dsol[0];
	
	TPZFMatrix<REAL> & phi_u = datavec[0].phi;
	
	int nsolu = datavec[0].sol.NElements();
	REAL w = 0.;
	
	for(int s = 0; s < nsolu; s++)
	{
		TPZManVector<REAL,3> sol_u = datavec[0].sol[s];
		REAL uy = sol_u[1];
		w += 2.*uy;
	}
	
	int phipRows = phi_p.Rows();
	int phiuCols = phi_u.Cols();
	
	REAL visc = this->fVisc;
	REAL deltaT = globFractInputData.actDeltaT();
	
	if(gState == ECurrentState) //current state (n+1): Matrix stiffnes
	{
		REAL actQl = globFractInputData.QlFVl(datavec[1].gelElId, sol_p[0]);
		REAL actdQldp = globFractInputData.dQlFVl(datavec[1].gelElId, sol_p[0]);
		
		for(int in = 0; in < phipRows; in++)
		{
			//----Residuo----
			//termo (wˆ3/(12*mi))*gradP * gradVp
			ef(phiuCols+in,0) += (-1.) * weight * (w*w*w/(12.*visc)) * dsol_p(0,0) * dphi_p(0,in);
			
			//termo w/deltaT * Vp
			ef(phiuCols+in,0) += (-1.) * weight * w/deltaT * phi_p(in,0);
			
			//termo 2Ql * Vp
			ef(phiuCols+in,0) += (-1.) * weight * (2.*actQl) * phi_p(in,0);
			
			
			//------Matriz tangente-----
			for(int jn = 0; jn < phiuCols; jn++)
			{
				//termo D[ (wˆ3/(12*mi))*gradP * gradVp , w ]
				ek(phiuCols+in, jn) += (+1.) * weight * ( 3.*w*w/(12.*visc) * (2.*phi_u(1,jn)) ) * dsol_p(0,0) * dphi_p(0,in);
				
				//termo D[ w/deltaT * Vp , w ]
				ek(phiuCols+in, jn) += (+1.) * weight * ( 2./deltaT * phi_u(1,jn) ) * phi_p(in,0);
			}
			for(int jn = 0; jn < phipRows; jn++)
			{
				//termo D[ (wˆ3/(12*mi))*gradP * gradVp , p ]
				ek(phiuCols+in, phiuCols+jn) += (+1.) * weight * (w*w*w/(12.*visc)) * dphi_p(0,in) * dphi_p(0,jn);
				
				//termo D[ 2Ql * Vp , p]
				ek(phiuCols+in, phiuCols+jn) += (+1.) * weight * (2.*actdQldp) * phi_p(in,0) * phi_p(jn,0);
			}
		}
	}
	
	//Last state (n): Matrix mass
	if(gState == ELastState)
	{
		for(int phip = 0; phip < phipRows; phip++)
		{            
			//termo w/deltaT * Vp
			ef(phiuCols+phip,0) += (+1.) * weight * w/deltaT * phi_p(phip,0);
		}
	}
}

void TPZNLFluidStructure2d::ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
{
	if(gState == ELastState)
	{
		return;
	}
	
	TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	
	TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	const REAL big  = TPZMaterial::gBigNumber;
	int phc = phi_u.Cols();
	
	for(int in = 0 ; in < phc; in++)
	{
		for(int il = 0; il < fNumLoadCases; il++)
		{
			//termo big*u*v do vetor de carga
			TPZFNMatrix<3,STATE> v2 = bc.Val2(il);
			ef(in,il) += big * ( v2(0,il)*phi_u(0,in) + v2(1,il)*phi_u(1,in) ) * weight;
			
			//termo big*u*v da matriz
			ef(in,il) += (-1.)*big * ( phi_u(0,in)*sol_u[0] + phi_u(1,in)*sol_u[1] ) * weight;
		}
		
		for (int jn = 0; jn < phc; jn++)
		{
			ek(in,jn) += big * ( phi_u(0,in)*phi_u(0,jn) + phi_u(1,in)*phi_u(1,jn) ) * weight;
		}
	}
}

void TPZNLFluidStructure2d::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef,TPZBndCond &bc)
{
	if(gState == ELastState)
	{
		return;
	}
	
	//    REAL auxvar, G;
	//    auxvar = 0.817*(1.-fPoiss)*globFractInputData.Hf();
	//    G = fE/(2.*(1. + fPoiss));
	REAL factor = 0.;//G/auxvar;
	
	TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	
	TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	
	REAL sol_un = sol_u[0]*datavec[0].normal[0] + sol_u[1]*datavec[0].normal[1];
	int nc_u = phi_u.Cols();
	
	TPZFMatrix<REAL> &phi_p = datavec[1].phi;
	TPZManVector<REAL,3> sol_p = datavec[1].sol[0];
	int phrp = phi_p.Rows();
	
	
	for (int in = 0; in < nc_u; in++)
	{
		//--- residuo ----
		for (int il = 0; il <fNumLoadCases; il++)
		{
			ef(in,il) += weight * factor * ( phi_u(0,in)*sol_un*datavec[0].normal[0] + phi_u(1,in)*sol_un*datavec[0].normal[1] )
			- weight * ( phi_u(0,in)*sol_p[0]*datavec[0].normal[0] + phi_u(1,in)*sol_p[0]*datavec[0].normal[1] );
		}
		
		//----- matriz tangente -----
		//termo k(phi_ix*phi_jun*nx + phi_iy*phi_jun*ny)
		for (int jn = 0; jn < nc_u; jn++)
		{
			REAL phi_jun = phi_u(0,jn)*datavec[0].normal[0] + phi_u(1,jn)*datavec[0].normal[1];
			ek(in,jn) += (-1.)*weight*factor*(phi_u(0,in)*phi_jun*datavec[0].normal[0] + phi_u(1,in)*phi_jun*datavec[0].normal[1]);
		}
		
		for (int jp = 0; jp < phrp; jp++)
		{
			ek(in,jp+nc_u) += weight*(phi_u(0,in)*phi_p(jp,0)*datavec[0].normal[0] + phi_u(1,in)*phi_p(jp,0)*datavec[0].normal[1]);
		}
	}
}

void TPZNLFluidStructure2d::ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc)
{
	if(gState == ELastState)
	{
		return;
	}
	
	TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	
	TPZManVector<REAL,3> sol_u =datavec[0].sol[0];
	int phc = phi_u.Cols();
	
	for(int in = 0 ; in < phc; in++)
	{
		for (int il = 0; il <fNumLoadCases; il++)
		{
			TPZFNMatrix<3,STATE> v2 = bc.Val2(il);
			ef(in,il)+= weight * ( v2(0,il)*phi_u(0,in) + v2(1,il)*phi_u(1,in) );
			
			//termo da matriz
			ef(in,il) += (-1.)*( bc.Val1()(0,0)*phi_u(0,in)*sol_u[0]*weight
													
													+ bc.Val1()(1,0)*phi_u(1,in)*sol_u[0]*weight
													
													+ bc.Val1()(0,1)*phi_u(0,in)*sol_u[1]*weight
													
													+ bc.Val1()(1,1)*phi_u(1,in)*sol_u[1]*weight );
		}
		for (int jn = 0; jn <phc; jn++)
		{
			
			ek(in,jn) += bc.Val1()(0,0)*phi_u(0,in)*phi_u(0,jn)*weight
			
			+ bc.Val1()(1,0)*phi_u(1,in)*phi_u(0,jn)*weight
			
			+ bc.Val1()(0,1)*phi_u(0,in)*phi_u(1,jn)*weight
			
			+ bc.Val1()(1,1)*phi_u(1,in)*phi_u(1,jn)*weight;
		}
	}
}


void TPZNLFluidStructure2d::ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc)
{
	if(gState == ELastState)
	{
		return;
	}
	
	TPZFMatrix<REAL> & phi_u = datavec[0].phi;
	int c_u = phi_u.Cols();
	
	TPZFMatrix<REAL> & phi_p = datavec[1].phi;
	int  phrp = phi_p.Rows();
	
	REAL Qinj = bc.Val2()(0,0);
	for(int in = 0; in < phrp; in++)
	{
		ef(in+c_u,0) += (-1.) * weight * Qinj * phi_p(in,0);
	}
}


void TPZNLFluidStructure2d::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc)
{
	TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
	if(shapetype!=datavec[0].EVecShape && datavec[0].phi.Cols()!= 0)
	{
		std::cout << " The space to elasticity problem must be reduced space.\n";
		DebugStop();
	}
	
	switch (bc.Type())
	{
		case 1:// Neumann condition in  both elasticity and pressure equations
		{
			ApplyNeumann_U(datavec, weight,ek,ef,bc);
			ContributePressure(datavec, weight, ek, ef);
			break;
		}
		case 2:// Dirichlet condition only on the elasticity equation
		{
			// Calculate the matrix contribution for pressure
			ApplyDirichlet_U(datavec, weight, ek,ef,bc);
			break;
		}
		case 3:// Mixed condition only on the elasticity equation
		{
			ApplyMixed_U(datavec, weight, ek,ef,bc);
			break;
		}
		case 4:// Neumann condition only on the pressure equation pressure
		{
			ApplyNeumann_P(datavec, weight, ek,ef,bc);
			break;
		}
		default:
		{
			DebugStop();
		}
	}
}

int TPZNLFluidStructure2d::VariableIndex(const std::string &name){
	if(!strcmp("Pressure",name.c_str()))        return 1;
	if(!strcmp("MinusKGradP",name.c_str()))     return 2;
	if(!strcmp("DisplacementX",name.c_str()))   return 3;
	if(!strcmp("DisplacementY",name.c_str()))   return 4;
	if(!strcmp("SigmaX",name.c_str()))          return 5;
	if(!strcmp("SigmaY",name.c_str()))          return 6;
	if(!strcmp("Displacement",name.c_str()))    return 7;
	if(!strcmp("W",name.c_str()))               return 8;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZNLFluidStructure2d::NSolutionVariables(int var){
	if(var == 1) return 1;
	if(var == 2) return 1;
	if(var == 3) return 1;
	if(var == 4) return 1;
	if(var == 5) return 1;
	if(var == 6) return 1;
	if(var == 7) return fDim;
	if(var == 8) return 1;
	
	return TPZMaterial::NSolutionVariables(var);
}


void TPZNLFluidStructure2d::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
	
	Solout.Resize(this->NSolutionVariables(var));
	
	TPZVec<REAL> SolP, SolU;
	TPZFMatrix<> DSolP, DSolU;
	TPZFMatrix<> axesP, axesU;
	SolP=datavec[1].sol[0];
	DSolP=datavec[1].dsol[0];
	axesP=datavec[1].axes;
	
	SolU=datavec[0].sol[0];
	DSolU=datavec[0].dsol[0];
	axesU=datavec[0].axes;
	
	REAL young, poisson, sigmaConf, visc, Hf;
	young = fE;
	poisson = fPoiss;
	sigmaConf = globFractInputData.SigmaConf();
	visc = fVisc;
	Hf = globFractInputData.Hf();
	
	if(var == 1)
	{
		//Aqui nao funciona pois nao vem datavec[1] preenchido
		//para elemento de contorno (que eh o nosso caso)!
		DebugStop();
		
		if(!datavec[1].phi)
		{
			return;
		}
		Solout[0] = SolP[0];
		return;
	}
	else if(var == 2)
	{
		//Aqui nao funciona pois nao vem datavec[1] preenchido
		//para elemento de contorno (que eh o nosso caso)!
		DebugStop();
		
		if(!datavec[1].phi)
		{
			return;
		}
		//        REAL G, un;
		//        G = young/(2.*(1.+poisson));
		//        un = 0.817*(1-poisson)*(SolP[0]-sigmaConf)*Hf/G;
		REAL factor = 0.;//(un*un*un)/(12.*visc);
		
		int id;
		TPZFNMatrix<9,REAL> dsoldx;
		TPZAxesTools<REAL>::Axes2XYZ(DSolP, dsoldx, axesP);
		for(id=0 ; id<1; id++)
		{
			Solout[id] = -1.*factor*dsoldx(id,0);
		}
		return;
	}
	else if(var == 3)//function (state variable ux)
	{
		Solout[0] = SolU[0];
		return;
	}
	else if(var == 4)//function (state variable uy)
	{
		Solout[0] = SolU[1];
		return;
	}
	else if(var == 5 || var == 6)
	{
		REAL DSolxy[2][2];
		REAL SigX, SigY, epsx, epsy;
		
		DSolxy[0][0] = DSolU(0,0)*axesU(0,0)+DSolU(1,0)*axesU(1,0);
		DSolxy[1][1] = DSolU(0,1)*axesU(0,1)+DSolU(1,1)*axesU(1,1);
		
		epsx = DSolU(0,0);// du/dx
		epsy = DSolU(1,1);// dv/dy
		REAL Gmodule = young/(1-poisson*poisson);
		
		if (this->fPlaneStress)
		{
			SigX = Gmodule*(epsx+poisson*epsy) + globFractInputData.PreStressXX();
			SigY = Gmodule*(poisson*epsx+epsy) + globFractInputData.PreStressYY();
		}
		else
		{
			SigX = young/((1.-2.*poisson)*(1.+poisson))*((1.-poisson)*epsx+poisson*epsy) + globFractInputData.PreStressXX();
			SigY = young/((1.-2.*poisson)*(1.+poisson))*(poisson*epsx+(1.-poisson)*epsy) + globFractInputData.PreStressYY();
		}
		if(var == 5)
		{
			Solout[0] = SigX;
			return;
		}
		else if(var == 6)
		{
			Solout[0] = SigY;
			return;
		}
	}
	else if(var == 7)
	{
		Solout[0] = SolU[0];
		Solout[1] = SolU[1];
		return;
	}
	else if(var == 8)
	{
		Solout[0] = 2.*SolU[1];
	}
}

void TPZNLFluidStructure2d::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
		datavec[i].SetAllRequirements(false);
		datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = true;
	}
}

void TPZNLFluidStructure2d::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
		datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNormal = true;
	}
}





