//
//  TPZPlasticFrac2D.cpp
//  PZ
//
//  Created by Agnaldo Farias on 9/17/12.
//
//

#include "TPZPlasticFrac2D.h"


#include <iostream>
#include <string>

#include "pzelasmat.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzintel.h"
#include "TPZElasticResponse.h"

#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
#endif




template<class T,class TMEM>
TPZPlasticFrac2D<T,TMEM>::TPZPlasticFrac2D() : TPZMatElastoPlastic2D<T,TMEM>()
{
	fmatId = 0;
	fE = 0.;
	fPoiss = 0.;
	fVisc = 0.;
	fDim = 2;
	fPlaneStress = 1;
	
	ff.resize(fDim);
	ff[0] = 0.;
	ff[1] = 0.;
	this->SetCurrentState();
  fSetRunPlasticity = false;
  flastElastFunction = new TPZLastElastFunction;
}

template<class T,class TMEM>
TPZPlasticFrac2D<T,TMEM>::TPZPlasticFrac2D(int matid, int dim, REAL young, REAL poiss, REAL visc) : TPZMatElastoPlastic2D<T,TMEM>(matid,1)
{
	fmatId = matid;
	fE = young;
	fPoiss = poiss;
	fVisc = visc;
	fDim = dim;
	fPlaneStress = 1; // it has to be one because of matelastoplastic initializer!!!!!
	
	ff.resize(2);
	ff[0] = 0.;
	ff[1] = 0.;
	this->SetCurrentState();
  fSetRunPlasticity = false;
  flastElastFunction = new TPZLastElastFunction;
}

template<class T,class TMEM>
TPZPlasticFrac2D<T,TMEM>::~TPZPlasticFrac2D()
{
}

template<class T,class TMEM>
int TPZPlasticFrac2D<T,TMEM>::NStateVariables()
{
	return 1;
}

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::Print(std::ostream &out)
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
template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::SetUpdateToUseFullU(bool update)
{
  fUpdateToUseFullDiplacement = update;
}

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
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
  
  //TPZFNMatrix<40> ekP(ek), efP(ef);
  if (fSetRunPlasticity)
  {
    ContributePlastic(datavec[0],weight,ek,ef);
   	ContributePressure(datavec, weight, ek, ef);
    return;
  }

	//Calculate the matrix contribution for elastic problem.
	TPZFMatrix<REAL> &dphi_u = datavec[0].dphix;
	TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	TPZFMatrix<REAL> &axes=datavec[0].axes;
	
	TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	TPZFMatrix<REAL> &dsol_u = datavec[0].dsol[0];
	
	int phcu = phi_u.Cols();
	int efcu = ef.Cols();
	
	if(this->fForcingFunction)
	{// phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE> res(3);
		this->fForcingFunction->Execute(datavec[0].x,res);
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
														 
														 fEover1MinNu2*poisson*dphix_i(0,0)*dphiy_j(1,0) + fEover21PlusNu*dphix_i(1,0)*dphiy_j(0,0) +
														 
														 fEover1MinNu2*poisson*dphiy_i(1,0)*dphix_j(0,0) + fEover21PlusNu*dphiy_i(0,0)*dphix_j(1,0) +
														 
														 fEover1MinNu2*dphiy_i(1,0)*dphiy_j(1,0) + fEover21PlusNu*dphiy_i(0,0)*dphiy_j(0,0));
			}
		}//fim para a tangente
	}
  

  // Old test to verify if matrixes are equal in case of non plastification
  /*
  REAL zeroK = 0.;
  REAL zeroF = 0.;
  for (int i = 0; i < phcu; i++) {
    zeroF += fabs(efP(i,0)-ef(i,0));
    for (int j = 0; j < phcu; j++) {
      zeroK += fabs(ekP(i,j)-ek(i,j));
    }
  }
  if (zeroF > 1.e-8 || zeroK > 1.e-8) {
    //DebugStop();
  }
   */
  
  
#ifdef LOG4CXX
  if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::Contribute ***";
		sout << " Resultant rhs vector:\n" << ef;
		sout << " Resultant stiff vector:\n" << ek;
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
	
	// Calculate the matrix contribution for pressure
	ContributePressure(datavec, weight, ek, ef);
}

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::ContributePlastic(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
	TPZFMatrix<REAL> &dphi = data.dphix, dphiXY;
	TPZFMatrix<REAL> &phi  = data.phi;
	//TPZFMatrix<REAL> &axes = data.axes, axesT;
	dphiXY = dphi;
	//axes.Transpose(&axesT);
	//axesT.Multiply(dphi,dphiXY);
	
	const int phc = phi.Cols();
	
	TPZFNMatrix<4>  Deriv(2,2);
	TPZFNMatrix<9> Dep(3,3);
	TPZFNMatrix<3>  DeltaStrain(3,1);
	TPZFNMatrix<3>  Stress(3,1);
  int ptindex = data.intGlobPtIndex;
  
  
  //	feclearexcept(FE_ALL_EXCEPT);
  //	int res = fetestexcept(FE_ALL_EXCEPT);
  //	if(res)
  //	{
  //		std::cout << " \n " << __PRETTY_FUNCTION__ <<"\n NAN DETECTED \n";
  //		DebugStop();
  //	}
  //
  
  if (fUpdateToUseFullDiplacement){ // SO I CAN RUN PLASTICITY WITH u AND NOT DeltaU
    TPZMatWithMem<TMEM>::MemItem(ptindex).fPlasticState.fEpsT.Zero();
    int solsize = data.sol[0].size();
		for(int i=0; i<solsize; i++)
    {
      TPZMatWithMem<TMEM>::MemItem(ptindex).fDisplacement[i] = 0.;
    }
    return;
  }
  
  if (TPZMatWithMem<TMEM>::fUpdateMem && data.sol.size() > 1)
  {
    // Loop over the solutions if update memory is true
    TPZSolVec locsol(data.sol);
    TPZGradSolVec locdsol(data.dsol);
    int numsol = locsol.size();
    
    for (int is=0; is<numsol; is++)
    {
      data.sol[0] = locsol[is];
      data.dsol[0] = locdsol[is];
      
      this->ComputeDeltaStrainVector(data, DeltaStrain);
      this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    }
  }
  else
  {
    this->ComputeDeltaStrainVector(data, DeltaStrain);
    this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
  }
#ifdef MACOS
  feclearexcept(FE_ALL_EXCEPT);
  if(fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT	)) {
    std::cout << "division by zero reported\n";
    DebugStop();
  }
#endif
	

  
#ifdef LOG4CXX
  if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
		sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
		sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
		sout << "\ndata.axes = " << data.axes;
		sout << "\nDep " <<endl;
		sout << Dep(0,0) << "\t" << Dep(0,1) << "\t" << Dep(0,2) <<"\n";
		sout << Dep(1,0) << "\t" << Dep(1,1) << "\t" << Dep(1,2) <<"\n";
		sout << Dep(2,0) << "\t" << Dep(2,1) << "\t" << Dep(2,2) <<"\n";
		
		sout << "\nStress " <<endl;
		sout << Stress(0,0) << "\t" << Stress(1,0) << "\t" << Stress(2,0) <<"\n";
		
		sout << "\nDELTA STRAIN " <<endl;
		sout << DeltaStrain(0,0) << "\t" << DeltaStrain(1,0) << "\t" << DeltaStrain(2,0) <<"\n";
		sout << "data.phi" << data.phi;
		
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
  /*
   //NAN detector
   res = fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW );
   if(res)
   {
   std::cout << " \n " << __PRETTY_FUNCTION__ <<"\n NAN DETECTED \n";
   DebugStop();
   }
   */
  ptindex = 0;
	//int nstate = NStateVariables();
	REAL val;/*,val1,val2,val3,val4*/;
	
	int in;
	for(in = 0; in < phc; in++)
	{
		
		val  = this->fRhoB * this->fForce[0] * phi(0,in);
		val -= Stress(0,0) * dphiXY(0,in); //dphixdx
		val -= Stress(2,0) * dphiXY(1,in); //dphixdy
		ef(in,0) += weight * val;
    
		val  = this->fRhoB * this->fForce[1] * phi(1,in);
		val -= Stress(2,0) * dphiXY(2,in); //dphiydx
		val -= Stress(1,0) * dphiXY(3,in); //dphiydy
		ef(in,0) += weight * val;
		
		for( int jn = 0; jn < phc; jn++)
		{
      
      /* EXPLANATION
       phi and dphi are ordered in the following way for each reduced space i
       phi(0,i) = phix
       phi(1,i) = phiy
       
       dphiXY(0,i) = dphixdx
       dphiXY(1,i) = dphixdy
       dphiXY(2,i) = dphiydx
       dphiXY(3,i) = dphiydy
       */
			
			
			val  = 2. * Dep(0,0) * dphiXY(0,in)*dphiXY(0,jn);//dphixdxI*dphixdxJ
			val +=      Dep(0,2) * dphiXY(0,in)*dphiXY(1,jn);//dphixdxI*dphixdyJ
			val += 2. * Dep(2,0) * dphiXY(1,in)*dphiXY(0,jn);//dphixdyI*dphixdxJ
			val +=      Dep(2,2) * dphiXY(1,in)*dphiXY(1,jn);//dphixdyI*dphixdyJ
			val *= 0.5;
			ek(in,jn) += weight * val;
			
			val  =      Dep(0,2) * dphiXY(0,in)*dphiXY(2,jn);//dphixdxI*dphiydxJ
			val += 2. * Dep(0,1) * dphiXY(0,in)*dphiXY(3,jn);//dphixdxI*dphiydyJ
			val +=      Dep(2,2) * dphiXY(1,in)*dphiXY(2,jn);//dphixdyI*dphiydxJ
			val += 2. * Dep(2,1) * dphiXY(1,in)*dphiXY(3,jn);//dphixdyI*dphiydyJ
			val *= 0.5;
			ek(in,jn) += weight * val;
      
			val  = 2. * Dep(2,0) * dphiXY(2,in)*dphiXY(0,jn);//dphiydxI*dphixdxJ
			val +=      Dep(2,2) * dphiXY(2,in)*dphiXY(1,jn);//dphiydxI*dphixdyJ
			val += 2. * Dep(1,0) * dphiXY(3,in)*dphiXY(0,jn);//dphiydyI*dphixdxJ
			val	+=      Dep(1,2) * dphiXY(3,in)*dphiXY(1,jn);//dphiydyI*dphixdyJ
			val *= 0.5;
			ek(in,jn) += weight * val;
      
			val  =      Dep(2,2) * dphiXY(2,in)*dphiXY(2,jn);//dphiydxI*dphiydxJ
			val += 2. * Dep(2,1) * dphiXY(2,in)*dphiXY(3,jn);//dphiydxI*dphiydyJ
			val +=      Dep(1,2) * dphiXY(3,in)*dphiXY(2,jn);//dphiydyI*dphiydxJ
			val += 2. * Dep(1,1) * dphiXY(3,in)*dphiXY(3,jn);//dphiydyI*dphiydyJ
			val *= 0.5;
			ek(in,jn) += weight * val;
			
		}
	}
  
	
#ifdef LOG4CXX
  if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "<<< TPZPlasticFrac2D<T,TMEM>::ContributePlastic ***";
		sout << " Resultant rhs vector:\n" << ef;
		sout << " Resultant stiff vector:\n" << ek;
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif

}


template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::ContributePressure(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
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
    const STATE wconst = 0.01; // AQUINATHAN
    
		REAL actQl = globFractInputData.QlFVl(datavec[1].gelElId, sol_p[0]);
		REAL actdQldp = globFractInputData.dQlFVl(datavec[1].gelElId, sol_p[0]);
		
		for(int in = 0; in < phipRows; in++)
		{
			//----Residuo----
			//termo (wˆ3/(12*mi))*gradP * gradVp
      //ef(phiuCols+in,0) += (-1.) * weight * (wconst*wconst*wconst/(12.*visc)) * dsol_p(0,0) * dphi_p(0,in); w constante
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
				//ek(phiuCols+in, phiuCols+jn) += (+1.) * weight * (wconst*wconst*wconst/(12.*visc)) * dphi_p(0,in) * dphi_p(0,jn); w constante
        
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
      if (this->flastElastFunction->flastElastCMesh) {
        REAL uy = 0.;
        this->flastElastFunction->Execute(datavec[1].x,uy);
        w = 2.*uy;
      }
			ef(phiuCols+phip,0) += (+1.) * weight * w/deltaT * phi_p(phip,0);
		}
	}
	
#ifdef LOG4CXX
	if (logger->isDebugEnabled()) {
		std::stringstream str;
		str << "\n------- Contribute da Pressure -------" << std::endl;
		str << "GeoElId = " << datavec[1].gelElId << std::endl; 
		ek.Print("ek",str);
		ef.Print("ef",str);
		LOGPZ_DEBUG(logger,str.str())		
	}
#endif
}

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
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
		for(int il = 0; il < this->fNumLoadCases; il++)
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

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef,TPZBndCond &bc)
{
	if(gState == ELastState)
	{
		return;
	}
	
  REAL auxvar, G, factor;
  auxvar = 0.817*(1.-fPoiss)*globFractInputData.Hf();
  G = fE/(2.*(1. + fPoiss));
  bool EnglandGreen = true;
  if (EnglandGreen) {
    factor = G/auxvar;
  }
  else{
    factor = 0.;//G/auxvar;
  }

	
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
		for (int il = 0; il < this->fNumLoadCases; il++)
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

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc)
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
		for (int il = 0; il < this->fNumLoadCases; il++)
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

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc)
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

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc)
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

template<class T,class TMEM>
int TPZPlasticFrac2D<T,TMEM>::VariableIndex(const std::string &name){
  if(fSetRunPlasticity){ // Soh acontece
    return TPZMatElastoPlastic<T,TMEM>::VariableIndex(name);
    return TPZMaterial::VariableIndex(name);
  }
  
	if(!strcmp("Pressure",name.c_str()))        return 1;
	if(!strcmp("MinusKGradP",name.c_str()))     return 2;
	if(!strcmp("DisplacementX",name.c_str()))   return 3;
	if(!strcmp("DisplacementY",name.c_str()))   return 4;
	if(!strcmp("SigmaX",name.c_str()))          return 5;
	if(!strcmp("SigmaY",name.c_str()))          return 6;
	if(!strcmp("Displacement",name.c_str()))    return 7;
	if(!strcmp("W",name.c_str()))               return 8;
	if(!strcmp("Displacement",name.c_str()))    return 9;
  if(!strcmp("YieldSurface",name.c_str()))    return 10;
	
	return TPZMaterial::VariableIndex(name);
}

template<class T,class TMEM>
int TPZPlasticFrac2D<T,TMEM>::NSolutionVariables(int var){
  if(fSetRunPlasticity){ // Soh acontece
    return TPZMatElastoPlastic<T,TMEM>::NSolutionVariables(var);
    return TPZMaterial::NSolutionVariables(var);
  }
  
	if(var == 1) return 1;
	if(var == 2) return 1;
	if(var == 3) return 1;
	if(var == 4) return 1;
	if(var == 5) return 1;
	if(var == 6) return 1;
	if(var == 7) return fDim;
	if(var == 8) return 1;
 	if(var == 9) return 2;
  if(var == 10) return T::fNYields::NYield;
	
	return TPZMaterial::NSolutionVariables(var);
}

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){

	Solout.Resize(this->NSolutionVariables(var));
  
  if(fSetRunPlasticity && var != 7){ // Soh acontece
    TPZMatElastoPlastic2D<T,TMEM>::Solution(datavec[0],var,Solout);
    return;
  }
	
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
  else if(var == 9)
  {
    Solout[0] = SolU[0];
    Solout[1] = SolU[1];
  }
  else if(var == 10)
  {
    /*
    TPZTensor<REAL> & EpsT = TPZMatWithMem<TMEM>::fMemory[intPt].fPlasticState.fEpsT;
    TPZTensor<STATE> epsElastic(EpsT);
    epsElastic-=TPZMatWithMem<TMEM>::fMemory[intPt].fPlasticState.fEpsP;
    plasticloc.Phi(epsElastic,Solout);
     */
    Solout[0] = -6378.;
  }
}

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
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

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
		datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNormal = true;
	}
}

template<class T,class TMEM>
void TPZPlasticFrac2D<T,TMEM>::SetRunPlasticity(bool IsPlasticity)
{
  fSetRunPlasticity = IsPlasticity;
}


#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZSandlerDimaggio.h"

//template class TPZPlasticFrac2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZElastoPlasticMem>;
template class TPZPlasticFrac2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZPlasticFrac2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>;
//template class TPZPlasticFrac2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>;