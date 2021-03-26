/**
 * @file
 * @brief Contains implementations of the TPZNLElasticityMaterial methods.
 */

#include "pznlelasmat.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.elasticity.data"));
#endif

#include <fstream>
using namespace std;

TPZNLElasticityMaterial::TPZNLElasticityMaterial() : TPZMaterial(0) {
	fE	= -1.;  // Young modulus
	fnu	= -1.;   // poisson coefficient
	ff[0]	= 0.; // X component of the body force
	ff[1]	= 0.; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
	fEover1MinNu2 = -1.;  //G = E/2(1-nu);
	fEover21PlusNu = -1.;//E/(1-nu)
  
	
	//Added by Cesar 2001/03/16
	fPreStressXX = 0.;  //Prestress in the x direction
	fPreStressYY = 0.;  //Prestress in the y direction
	fPreStressXY = 0.;  //Prestress in the z direction
	fPreStressZZ = 0.;  //Prestress in the z direction
	fPlaneStress = -1;
  
  // Added by Philippe 2012
  fPostProcIndex = 0;
}

TPZNLElasticityMaterial::TPZNLElasticityMaterial(int id) : TPZMaterial(id) {
	fE	= -1.;  // Young modulus
	fnu	= -1.;   // poisson coefficient
	ff[0]	= 0.; // X component of the body force
	ff[1]	= 0.; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
	fEover1MinNu2 = -1.;  //G = E/2(1-nu);
	fEover21PlusNu = -1.;//E/(1-nu)
  
	
	//Added by Cesar 2001/03/16
	fPreStressXX = 0.;  //Prestress in the x direction
	fPreStressYY = 0.;  //Prestress in the y direction
	fPreStressXY = 0.;  //Prestress in the z direction
	fPreStressZZ = 0.;  //Prestress in the z direction
	fPlaneStress = -1;
  
  // Added by Philippe 2012
  fPostProcIndex = 0;
}

TPZNLElasticityMaterial::TPZNLElasticityMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress) : TPZMaterial(num) {
	
	fE	= E;  // Young modulus
	fnu	= nu;   // poisson coefficient
	ff[0]	= fx; // X component of the body force
	ff[1]	= fy; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
	fEover1MinNu2 = E/(1-fnu*fnu);  //G = E/2(1-nu);
	fEover21PlusNu = E/(2.*(1+fnu));//E/(1-nu)
	
	//Added by Cesar 2001/03/16
	fPreStressXX = 0.;  //Prestress in the x direction
	fPreStressYY = 0.;  //Prestress in the y direction
	fPreStressXY = 0.;  //Prestress in the z direction
	fPreStressZZ = 0.;  //Prestress in the z direction
	fPlaneStress = plainstress;
  // Added by Philippe 2012
  fPostProcIndex = 0;
}

TPZNLElasticityMaterial::~TPZNLElasticityMaterial() {
}



void TPZNLElasticityMaterial::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "\tE   = " << fE   << endl;
	out << "\tnu   = " << fnu   << endl;
	out << "\tF   = " << ff[0] << ' ' << ff[1]   << endl;
	out << "\t PreStress: \n"
	<< "Sigma xx = \t" << fPreStressXX << "\t"
	<< "Sigma yy = \t" << fPreStressYY << "\t"
	<< "Sigma xy = \t" << fPreStressXY << "Sigma zz = \t" << fPreStressZZ << endl;
}

//Added by Cesar 2001/03/16
void TPZNLElasticityMaterial::SetPreStress(REAL Sigxx, REAL Sigyy, REAL Sigxy, REAL Sigzz){
	fPreStressXX = Sigxx;
	fPreStressYY = Sigyy;
	fPreStressXY = Sigxy;
  fPreStressZZ = Sigzz;
}


void TPZNLElasticityMaterial::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
  
  TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
  if(shapetype==data.EVecShape){
    ContributeVecShape(data,weight,ek, ef);
    return;
  }
  
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &axes=data.axes;
	
  TPZVec<REAL> &sol = data.sol[0];
  TPZFMatrix<REAL> &dsol = data.dsol[0];
  
	int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphi.Cols();
	dphr = dphi.Rows();
	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
	if(phc != 1 || dphr != 2 || phr != dphc ||
	   ekr != phr*2 || ekc != phr*2 ||
	   efr != phr*2 ){
		PZError << "\nTPZNLElasticityMaterial.contr, inconsistent input data : \n" <<
		"phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
    << ek.Cols() <<
		"\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
    << ef.Cols() << "\n";
		return;
		//		PZError.show();
	}
	if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE,3> res(3);
		fForcingFunction->Execute(data.x,res);
		ff[0] = res[0];
		ff[1] = res[1];
		ff[2] = res[2];
	}
	
	TPZFNMatrix<2,STATE> dphix(2,1), dphiy(2,1), du(2,1),dv(2,1), dphixj(2,1), dphiyj(2,1);
	/*
	 * Plain strain materials values
	 */
	REAL nu1 = 1. - fnu;//(1-nu)
	REAL nu2 = (1.-2.*fnu)/2.;
	REAL F = fE/((1.+fnu)*(1.-2.*fnu));
	
	for( int in = 0; in < phr; in++ ) {
		dphix(0,0) = dphi(0,in)*axes(0,0)+dphi(1,in)*axes(1,0);//dvx
		dphix(1,0) = dphi(0,in)*axes(0,1)+dphi(1,in)*axes(1,1);//dvy
    dphiy = dphix;
    
		du(0,0) = dsol(0,0)*axes(0,0)+dsol(1,0)*axes(1,0);
		du(1,0) = dsol(0,0)*axes(0,1)+dsol(1,0)*axes(1,1);
		dv(0,0) = dsol(0,1)*axes(0,0)+dsol(1,1)*axes(1,0);
		dv(1,0) = dsol(0,1)*axes(0,1)+dsol(1,1)*axes(1,1);
    
    for (int col = 0; col < efc; col++)
    {
      ef(2*in, col) += (-1.) * weight * (ff[0]*phi(in,0) - dphix(0,0)*fPreStressXX - dphix(1,0)*fPreStressXY);  // direcao x
      ef(2*in+1, col) += (-1.) * weight * (ff[1]*phi(in,0) - dphiy(0,0)*fPreStressXY - dphiy(1,0)*fPreStressYY);// direcao y <<<----
      
      if (fPlaneStress != 1) {
        DebugStop();
      }
      else {
        ef(2*in, col) += (-1.) * weight * ( (fEover1MinNu2 * du(0,0) + fEover1MinNu2 * fnu * dv(1,0) ) * dphix(0,0)
                                           + fEover21PlusNu * (du(1,0)+dv(0,0)) * dphix(1,0) );  // direcao x
        
        ef(2*in+1, col) += (-1.) * weight * ( (fEover1MinNu2 * fnu * du(0,0) + fEover1MinNu2 * dv(1,0) ) * dphiy(1,0)
                                             + fEover21PlusNu * (du(1,0)+dv(0,0)) * dphiy(0,0) );// direcao y <<<----
      }
    }
		
		
		for( int jn = 0; jn < phr; jn++ ) {
      dphixj(0,0) = dphi(0,jn)*axes(0,0)+dphi(1,jn)*axes(1,0);
      dphixj(1,0) = dphi(0,jn)*axes(0,1)+dphi(1,jn)*axes(1,1);
      dphiyj = dphixj;
			
			
			if (fPlaneStress != 1){
				/* Plain Strain State */
				ek(2*in,2*jn) += weight * (
                                   nu1 * dphix(0,0)*dphix(0,1)+ nu2 * dphix(1,0)*dphix(1,1)
                                   ) * F;
				
				ek(2*in,2*jn+1) += weight * (
                                     fnu*dphix(0,0)*dphix(1,1)+ nu2*dphix(1,0)*dphix(0,1)
                                     ) * F;
				
				ek(2*in+1,2*jn) += weight * (
                                     fnu*dphix(1,0)*dphix(0,1)+ nu2*dphix(0,0)*dphix(1,1)
                                     ) * F;
				
				ek(2*in+1,2*jn+1) += weight * (
                                       nu1*dphix(1,0)*dphix(1,1)+ nu2*dphix(0,0)*dphix(0,1)
                                       ) * F;
			}
			else{
				/* Plain stress state */
				ek(2*in,2*jn) += weight * (
                                   fEover1MinNu2 * dphix(0,0)*dphixj(0,0)+ fEover21PlusNu * dphix(1,0)*dphixj(1,0)
                                   );
				
				ek(2*in,2*jn+1) += weight * (
                                     fEover1MinNu2*fnu*dphix(0,0)*dphiyj(1,0)+ fEover21PlusNu*dphix(1,0)*dphiyj(0,0)
                                     );
				
				ek(2*in+1,2*jn) += weight * (
                                     fEover1MinNu2*fnu*dphiy(1,0)*dphixj(0,0)+ fEover21PlusNu*dphiy(0,0)*dphixj(1,0)
                                     );
				
				ek(2*in+1,2*jn+1) += weight * (
                                       fEover1MinNu2*dphiy(1,0)*dphiyj(1,0)+ fEover21PlusNu*dphiy(0,0)*dphiyj(0,0)
                                       );
			}
		}
	}
	
  //#ifdef LOG4CXX
  //	if(logdata->isDebugEnabled())
  //	{
  //		std::stringstream sout;
  //		ek.Print("ek_elastmat = ",sout,EMathematicaInput);
  //		ef.Print("ef_elastmat = ",sout,EMathematicaInput);
  //		LOGPZ_DEBUG(logdata,sout.str())
  //	}
  //#endif
	
}



void TPZNLElasticityMaterial::FillDataRequirements(TPZMaterialData &data)
{
  data.fNeedsSol = true;
  data.fNeedsNormal = false;
  
}

void TPZNLElasticityMaterial::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
  data.fNeedsSol = true;
  data.fNeedsNormal = true;
  if (type == 4 || type == 5 || type == 6) {
    data.fNeedsNormal = true;
  }
}

void TPZNLElasticityMaterial::ContributeVecShape(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
  TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &axes=data.axes;
	
	int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphi.Cols();
	dphr = dphi.Rows();
	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
	
	if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE> res(3);
		fForcingFunction->Execute(data.x,res);
		ff[0] = res[0];
		ff[1] = res[1];
		ff[2] = res[2];
	}
	
	TPZFNMatrix<4,STATE> dphix_i(2,1),dphiy_i(2,1), dphix_j(2,1), dphiy_j(2,1);
	/*
	 * Plain strain materials values
	 */
	REAL nu1 = 1 - fnu;//(1-nu)
	REAL nu2 = (1-2*fnu)/2;
	REAL F = fE/((1+fnu)*(1-2*fnu));
  
	for( int in = 0; in < phc; in++ )
  {
		dphix_i(0,0) = dphi(0,in)*axes(0,0)+dphi(1,in)*axes(1,0);
		dphix_i(1,0) = dphi(0,in)*axes(0,1)+dphi(1,in)*axes(1,1);
		dphiy_i(0,0) = dphi(2,in)*axes(0,0)+dphi(3,in)*axes(1,0);
		dphiy_i(1,0) = dphi(2,in)*axes(0,1)+dphi(3,in)*axes(1,1);
		
    for (int col = 0; col < efc; col++)
    {
      ef(in,col) += weight*(   ff[0] * phi(0, in)- dphix_i(0,0)*fPreStressXX - dphix_i(1,0)*fPreStressXY
                            + ff[1] * phi(1, in)- dphiy_i(0,0)*fPreStressYY - dphiy_i(1,0)*fPreStressXY);
    }
		for( int jn = 0; jn < phc; jn++ ) {
      
      dphix_j(0,0) = dphi(0,jn)*axes(0,0)+dphi(1,jn)*axes(1,0);
      dphix_j(1,0) = dphi(0,jn)*axes(0,1)+dphi(1,jn)*axes(1,1);
      dphiy_j(0,0) = dphi(2,jn)*axes(0,0)+dphi(3,jn)*axes(1,0);
      dphiy_j(1,0) = dphi(2,jn)*axes(0,1)+dphi(3,jn)*axes(1,1);
			
			
			if (fPlaneStress != 1){
				/* Plain Strain State */
				ek(in,jn) += weight*(nu1*dphix_i(0,0)*dphix_j(0,0) + nu2*dphix_i(1,0)*dphix_j(1,0) +
                             
                             fnu*dphix_i(0,0)*dphiy_j(1,0) + nu2*dphix_i(1,0)*dphiy_j(0,0) +
                             
                             fnu*dphiy_i(1,0)*dphix_j(0,0) + nu2*dphiy_i(0,0)*dphix_j(1,0) +
                             
                             nu1*dphiy_i(1,0)*dphiy_j(1,0) + nu2*dphiy_i(0,0)*dphiy_j(0,0))*F;
			}
			else{
				/* Plain stress state */
        
        ek(in,jn) += weight*(fEover1MinNu2*dphix_i(0,0)*dphix_j(0,0) + fEover21PlusNu*dphix_i(1,0)*dphix_j(1,0) +
                             
                             fEover1MinNu2*dphix_i(0,0)*dphiy_j(1,0) + fEover21PlusNu*dphix_i(1,0)*dphiy_j(0,0) +
                             
                             fEover1MinNu2*dphiy_i(1,0)*dphix_j(0,0) + fEover21PlusNu*dphiy_i(0,0)*dphix_j(1,0) +
                             
                             fEover1MinNu2*dphiy_i(1,0)*dphiy_j(1,0) + fEover21PlusNu*dphiy_i(0,0)*dphiy_j(0,0));
      }
		}
	}
}


void TPZNLElasticityMaterial::ContributeBC(TPZMaterialData &data,REAL weight,
                                           TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
  
  
  TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
  if(shapetype==data.EVecShape){
    ContributeVecShapeBC(data,weight,ek, ef,bc);
    return;
  }
  
	TPZFMatrix<REAL> &phi = data.phi;
  int dim = Dimension();
  
	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
  
	int phr = phi.Rows();
	short in,jn;
	REAL v2[2];
	v2[0] = bc.Val2()(0,0);
	v2[1] = bc.Val2()(1,0);
	
  //		In general when the problem is  needed to stablish any convention for ContributeBC implementations
  
  //     REAL v2[2];
  // 	v2[0] = bc.Val2()(0,0);
  // 	v2[1] = bc.Val2()(1,0);
  
	TPZFMatrix<STATE> &v1 = bc.Val1();
  int nstate = NStateVariables();
  
	
  REAL dif0 = (data.sol[0][0]-v2[0]);
  REAL dif1 = (data.sol[0][1]-v2[1]);
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
		{
			for(in = 0 ; in < phr; in++) {
				ef(2*in,0)   += - BIGNUMBER * dif0 * phi(in,0) * weight;        // forced v2 displacement
				ef(2*in+1,0) += - BIGNUMBER * dif1 * phi(in,0) * weight;        // forced v2 displacement
				for (jn = 0 ; jn < phi.Rows(); jn++)
        {
					ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
					ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
				}
			}
		}
			break;
			
		case 1 :		// Neumann condition
		{
      for (in = 0; in < phr; in++)
      {
        for (int il = 0; il <fNumLoadCases; il++)
        {
          TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
          ef(2*in,il) += (-1.) * (- v2(0,0)) * phi(in,0) * weight;        // force in x direction
          ef(2*in+1,il) += (-1.) * (- v2(1,0)) * phi(in,0) * weight;      // force in y direction
        }
      }
		}
			break;
			
		case 2 :		// Mixed Condition
		{
			for(in = 0 ; in < phi.Rows(); in++)
      {
        for (int il = 0; il <fNumLoadCases; il++)
        {
          TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
          ef(2*in,il) += v2(0,0) * phi(in,0) * weight;        // force in x direction
          ef(2*in+1,il) += v2(1,0) * phi(in,0) * weight;      // forced in y direction
        }
				
				for (jn = 0 ; jn < phi.Rows(); jn++) {
					ek(2*in,2*jn) += bc.Val1()(0,0) * phi(in,0) *
					phi(jn,0) * weight;         // peso de contorno => integral de contorno
					ek(2*in+1,2*jn) += bc.Val1()(1,0) * phi(in,0) *
					phi(jn,0) * weight;
					ek(2*in+1,2*jn+1) += bc.Val1()(1,1) * phi(in,0) *
					phi(jn,0) * weight;
					ek(2*in,2*jn+1) += bc.Val1()(0,1) * phi(in,0) *
					phi(jn,0) * weight;
				}
			}   // este caso pode reproduzir o caso 0 quando o deslocamento
      break;
      
    case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
		{
			for(in = 0 ; in < phr; in++) {
				ef(2*in,0)   += - BIGNUMBER * dif0 * phi(in,0) * weight;        // forced v2 displacement
				//ef(2*in+1,0) += - BIGNUMBER * dif1 * phi(in,0) * weight;        // forced v2 displacement
				for (jn = 0 ; jn < phi.Rows(); jn++)
				{
					ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
					//ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
				}
			}
		}
		break;
      
      
    case 4: // stressField Neumann condition
      for(in = 0; in < dim; in ++)
      {
        v2[in] =  ( v1(in,0) * data.normal[0] +
                   v1(in,1) * data.normal[1]);
      }
      // The normal vector points towards the neighbour. The negative sign is there to
      // reflect the outward normal vector.
      for(in = 0 ; in < phi.Rows(); in++) {
        ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
        ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
        //	cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
        //	cout << "val2:  " << v2[0]  << endl;
      }
      break;
      
		case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
      {
        TPZFNMatrix<2,STATE> res(2,1,0.);
        for(in = 0 ; in < phi.Rows(); in++)
        {
          ef(nstate*in+0,0) += (v2[0]*data.normal[0]) * phi(in,0) * weight ;
          ef(nstate*in+1,0) += (v2[0]*data.normal[1]) * phi(in,0) * weight ;
          for(jn=0; jn<phi.Rows(); jn++)
          {
            for(int idf=0; idf<2; idf++) for(int jdf=0; jdf<2; jdf++)
            {
              ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*data.normal[idf]*data.normal[jdf]*phi(in,0)*phi(jn,0)*weight;
              //BUG FALTA COLOCAR VAL2
              //                        DebugStop();
            }
          }
          
        }
      }
			break;
			
		case 6://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
      {
        TPZFNMatrix<2,STATE> res(2,1,0.);
        for(in = 0 ; in < phi.Rows(); in++)
        {
          ef(nstate*in+0,0) += (v2[0]*data.normal[0]) * phi(in,0) * weight ;
          ef(nstate*in+1,0) += (v2[0]*data.normal[1]) * phi(in,0) * weight ;
          for(jn=0; jn<phi.Rows(); jn++)
          {
            for(int idf=0; idf<2; idf++) for(int jdf=0; jdf<2; jdf++)
            {
              ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
              //BUG FALTA COLOCAR VAL2
              //                        DebugStop();
            }
          }
          
        }
        
      }
      break;
      
      
      
      
    }      // �nulo introduzindo o BIGNUMBER pelos valores da condi�o
  } // 1 Val1 : a leitura �00 01 10 11
}


void TPZNLElasticityMaterial::ContributeVecShapeBC(TPZMaterialData &data,REAL weight,
                                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
  
  TPZFMatrix<REAL> &phi = data.phi;
  
	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
  
	int phc = phi.Cols();
	short in,jn;
	
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
			for(in = 0 ; in < phc; in++) {
        for (int il = 0; il <fNumLoadCases; il++)
        {
          TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
          
          ef(in,il) += weight*BIGNUMBER*(v2(0,il)*phi(0,in) + v2(1,il) * phi(1,in));
        }
				for (jn = 0 ; jn < phc; jn++) {
          
          ek(in,jn) += weight*BIGNUMBER*(phi(0,in)*phi(0,jn) + phi(1,in)*phi(1,jn));
				}
			}
			break;
			
		case 1 :			// Neumann condition
      for (in = 0; in < phc; in++)
      {
        for (int il = 0; il <fNumLoadCases; il++)
        {
          TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
          ef(in,il)+= weight*(v2(0,il)*phi(0,in) + v2(1,il)*phi(1,in));
        }
      }
			break;
			
		case 2 :		// condicao mista
			for(in = 0 ; in < phc; in++)
      {
        for (int il = 0; il <fNumLoadCases; il++)
        {
          TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
          ef(in,il) += weight * (v2(0,il)*phi(0,in) + v2(1,il)*phi(1,in));
        }
				
				for (jn = 0; jn <phc; jn++) {
          
          ek(in,jn) += bc.Val1()(0,0)*phi(0,in)*phi(0,jn)*weight
          
          + bc.Val1()(1,0)*phi(1,in)*phi(0,jn)*weight
          
          + bc.Val1()(0,1)*phi(0,in)*phi(1,jn)*weight
          
          + bc.Val1()(1,1)*phi(1,in)*phi(1,jn)*weight;
				}
			}// este caso pode reproduzir o caso 0 quando o deslocamento
	}      //  eh nulo introduzindo o BIGNUMBER pelos valores da condicao
}


/** Returns the variable index associated with the name. */
int TPZNLElasticityMaterial::VariableIndex(const std::string &name){
  
  
  /*
   if(!strcmp("Displacement",             name.c_str()))  return TPZNLElasticityMaterial::EDisplacement;
   if(!strcmp("DisplacementX",            name.c_str()))  return TPZNLElasticityMaterial::EDisplacementX;
   if(!strcmp("DisplacementY",            name.c_str()))  return TPZNLElasticityMaterial::EDisplacementY;
   if(!strcmp("DisplacementZ",            name.c_str()))  return TPZNLElasticityMaterial::EDisplacementZ;
   if(!strcmp("NormalStress",             name.c_str()))  return TPZNLElasticityMaterial::ENormalStress;
   if(!strcmp("ShearStress",              name.c_str()))  return TPZNLElasticityMaterial::EShearStress;
   if(!strcmp("NormalStrain",             name.c_str()))  return TPZNLElasticityMaterial::ENormalStrain;
   if(!strcmp("ShearStrain",              name.c_str()))  return TPZNLElasticityMaterial::EShearStrain;
   if(!strcmp("PrincipalStress",          name.c_str()))  return TPZNLElasticityMaterial::EPrincipalStress;
   if(!strcmp("Stress1",                  name.c_str()))  return TPZNLElasticityMaterial::EStress1;
   if(!strcmp("PrincipalStrain",          name.c_str()))  return TPZNLElasticityMaterial::EPrincipalStrain;
   if(!strcmp("Strain1",                  name.c_str()))  return TPZNLElasticityMaterial::EStrain1;
   if(!strcmp("PrincipalStressDirection1",name.c_str()))  return TPZNLElasticityMaterial::EPrincipalStressDirection1;
   if(!strcmp("PrincipalStressDirection2",name.c_str()))  return TPZNLElasticityMaterial::EPrincipalStressDirection2;
   if(!strcmp("PrincipalStressDirection3",name.c_str()))  return TPZNLElasticityMaterial::EPrincipalStressDirection3;
   if(!strcmp("I1Stress",                 name.c_str()))  return TPZNLElasticityMaterial::EI1Stress;
   if(!strcmp("J2Stress",                 name.c_str()))  return TPZNLElasticityMaterial::EJ2Stress;
   if(!strcmp("I1J2Stress",               name.c_str()))  return TPZNLElasticityMaterial::EI1J2Stress;
   if(!strcmp("DirStress",                name.c_str()))  return TPZNLElasticityMaterial::EDirStress;
   if(!strcmp("DirStrain",                name.c_str()))  return TPZNLElasticityMaterial::EDirStrain;
   if(!strcmp("VolElasticStrain",         name.c_str()))  return TPZNLElasticityMaterial::EVolElasticStrain;
   if(!strcmp("VolPlasticStrain",         name.c_str()))  return TPZNLElasticityMaterial::EVolPlasticStrain;
   if(!strcmp("VolTotalStrain",           name.c_str()))  return TPZNLElasticityMaterial::EVolTotalStrain;
   if(!strcmp("VolTEPStrain",             name.c_str()))  return TPZNLElasticityMaterial::EVolTEPStrain;
   if(!strcmp("Alpha",                    name.c_str()))  return TPZNLElasticityMaterial::EAlpha;
   if(!strcmp("PlasticSteps",             name.c_str()))  return TPZNLElasticityMaterial::EPlasticSteps;
   if(!strcmp("YieldSurface",             name.c_str()))  return TPZNLElasticityMaterial::EYield;
   if(!strcmp("TotalPlasticStrain",     name.c_str()))  return TPZNLElasticityMaterial::ENormalPlasticStrain;
   if(!strcmp("EMisesStress",     name.c_str()))  return TPZNLElasticityMaterial::EMisesStress;
   PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
   return -1;
   */
  
	if(!strcmp("displacement",name.c_str()))     return 9;
	if(!strcmp("Displacement",name.c_str()))     return 9;
	if(!strcmp("DisplacementMem",name.c_str()))     return 9;
	if(!strcmp("Pressure",name.c_str()))         return 1;
	if(!strcmp("MaxStress",name.c_str()))        return 2;
	if(!strcmp("PrincipalStress1",name.c_str())) return 3;
	if(!strcmp("PrincipalStress2",name.c_str())) return 4;
	if(!strcmp("SigmaX",name.c_str()))           return 5;
	if(!strcmp("SigmaY",name.c_str()))           return 6;
	if(!strcmp("TauXY",name.c_str()))            return 8;//Cedric
	if(!strcmp("Strain",name.c_str()))           return 11;//Philippe
	if(!strcmp("SigmaZ",name.c_str()))           return 12;//Philippe
  
	if(!strcmp("sig_x",name.c_str()))            return 5;
	if(!strcmp("sig_y",name.c_str()))            return 6;
	if(!strcmp("tau_xy",name.c_str()))           return 8;//Cedric
	if(!strcmp("Displacement6",name.c_str()))    return 7;
	if(!strcmp("Stress",name.c_str()))           return 10;
	if(!strcmp("Flux",name.c_str()))           return 10;
  if(!strcmp("J2",name.c_str()))           return 20;
  if(!strcmp("I1",name.c_str()))           return 21;
  if(!strcmp("J2Stress",name.c_str()))           return 20;
  if(!strcmp("I1Stress",name.c_str()))           return 21;
  if(!strcmp("Alpha",name.c_str()))        return 22;
  if(!strcmp("PlasticSqJ2",name.c_str()))        return 22;
  if(!strcmp("PlasticSqJ2El",name.c_str()))        return 22;
  if(!strcmp("YieldSurface",name.c_str()))        return 27;
  if(!strcmp("NormalStress",name.c_str()))        return 23;
  if(!strcmp("ShearStress",name.c_str()))        return 24;
  if(!strcmp("NormalStrain",name.c_str()))        return 25;
  if(!strcmp("ShearStrain",name.c_str()))        return 26;
  
  
  
	//   cout << "TPZNLElasticityMaterial::VariableIndex Error\n";
	return TPZMaterial::VariableIndex(name);
}

/** Returns the number of variables associated with the variable indexed by var. */
int TPZNLElasticityMaterial::NSolutionVariables(int var){
  
	switch(var) {
		case 0:
			return 2;
		case 1:
		case 2:
			return 1;
		case 3:
		case 4:
			return 2;
		case 5:
		case 6:
		case 8:
			return 1;
		case 7:
			return 6;
		case 9:
			return 3;
		case 10 : //Stress Tensor
			return 3;
    case 11 : //Strain Tensor
      return 3;
      // SigZ
    case 12:
      return 1;
    case 20:
      return 1;
    case 21:
      return 1;
    case 22:
      return 1;
    case 23:
    case 24:
    case 25:
    case 26:
    case 27:
      return 3;
		default:
			return TPZMaterial::NSolutionVariables(var);
	}
}

void TPZNLElasticityMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
  int numbersol = data.dsol.size();
  int ipos = 0;
  if (fPostProcIndex < numbersol) {
    ipos = fPostProcIndex;
  }
  
  TPZVec<STATE> &Sol = data.sol[ipos];
  TPZFMatrix<STATE> &DSol = data.dsol[ipos];
  TPZFMatrix<REAL> &axes = data.axes;
	
	REAL epsx;
	REAL epsy;
	REAL epsxy;
  REAL epsz = 0.;
	REAL SigX;
	REAL SigY;
  REAL SigZ;
	REAL TauXY,aux,Sig1,Sig2,angle,DSolxy[2][2];
  
  // dudx - dudy
	DSolxy[0][0] = DSol(0,0)*axes(0,0)+DSol(1,0)*axes(1,0);
	DSolxy[1][0] = DSol(0,0)*axes(0,1)+DSol(1,0)*axes(1,1);
	// dvdx - dvdy
	DSolxy[0][1] = DSol(0,1)*axes(0,0)+DSol(1,1)*axes(1,0);
	DSolxy[1][1] = DSol(0,1)*axes(0,1)+DSol(1,1)*axes(1,1);
  
  epsx = DSolxy[0][0];// du/dx
  epsy = DSolxy[1][1];// dv/dy
  epsxy = 0.5*(DSolxy[1][0]+DSolxy[0][1]);
  
  REAL lambda = GetLambda();
  REAL mu = GetMU();
  if (this->fPlaneStress) {
    REAL lambda = GetLambda();
    REAL mu = GetMU();
    epsz = -lambda*(epsx+epsy)/(lambda+2.*mu);
  }
  else {
    epsz = 0.;
  }
  TauXY = 2*mu*epsxy+fPreStressXY;
#ifdef PZDEBUG
  REAL TauXY2 = fE*epsxy/(1.+fnu)+fPreStressXY;
#ifdef REALfloat
  if (fabs(TauXY-TauXY2) > 1.e-10) {
    DebugStop();
  }
#else
  if (fabs(TauXY-TauXY2) > 1.e-6) {
    DebugStop();
  }
#endif
#endif
  if (this->fPlaneStress){
    SigX = fEover1MinNu2*(epsx+fnu*epsy)+fPreStressXX;
    SigY = fEover1MinNu2*(fnu*epsx+epsy)+fPreStressYY;
    SigZ = fPreStressZZ;
  }
  else
  {
    SigX = fE/((1.-2.*fnu)*(1.+fnu))*((1.-fnu)*epsx+fnu*epsy)+fPreStressXX;
    SigY = fE/((1.-2.*fnu)*(1.+fnu))*(fnu*epsx+(1.-fnu)*epsy)+fPreStressYY;
    SigZ = fPreStressZZ+lambda*(epsx+epsy);
  }
  
	switch(var) {
		case 0:
			//numvar = 2;
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			break;
		case 7:
			//numvar = 6;
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			Solout[2] = 0.;
			Solout[3] = 0.;
			Solout[4] = 0.;
			Solout[5] = 0.;
			break;
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 8:
		case 10:
			
			//numvar = 1;
			Solout[0] = SigX+SigY+SigZ;
      // Pressure variable
			if(var == 1) {
				Solout[0] = SigX+SigY+SigZ;
				return;
			}
      // TauXY variable
			if(var == 8) {
				Solout[0] = TauXY;
				return;
			}
			if(var ==5) {
				Solout[0] = SigX;
				return;
			}
			if(var == 6) {
				Solout[0] = SigY;
				return;
			}
			aux = sqrt(0.25*(SigX-SigY)*(SigX-SigY)
                 +(TauXY)*(TauXY));
			// Philippe 13/5/99
			//         if(abs(Tau) < 1.e-10 && abs(SigY-SigX) < 1.e-10) angle = 0.;
			if(fabs(TauXY) < 1.e-10 && fabs(SigY-SigX) < 1.e-10) angle = 0.;
			else angle = atan2(2*TauXY,SigY-SigX)/2.;
			Sig1 = 0.5*(SigX+SigY)+aux;
			Sig2 = 0.5*(SigX+SigY)-aux;
			if(var == 3 ){
				//numvar = 2;
				Solout[0] = Sig1*cos(angle);
				Solout[1] = Sig1*sin(angle);
				return;
			}
			if(var == 4 ) {
				//numvar = 2;
				Solout[0] = -Sig2*sin(angle);
				Solout[1] = Sig2*cos(angle);
				return;
			}
			if(var == 2) {
				REAL sigmax;
				sigmax = (fabs(Sig1) < fabs(Sig2))? fabs(Sig2) : fabs(Sig1);
				Solout[0] = sigmax;
				return;
			}
			if (var ==10)
			{
				Solout[0] = SigX;
				Solout[1] = SigY;
				Solout[2] = TauXY;
				return;
			}
			cout << "Very critical error TPZNLElasticityMaterial::Solution\n";
			exit(-1);
			//         Solout[0] /= 0.;
			break;
		case 9:
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			Solout[2] = 0.;
			break;
    case 11:
      Solout[0] = epsx;
			Solout[1] = epsy;
			Solout[2] = epsxy;
      break;
    case 12:
      Solout[0] = SigZ;
      break;
      
    case 20:
    {
      
      REAL J2 = (pow(SigX + SigY,2) - (3*(-pow(SigX,2) - pow(SigY,2) + pow(SigX + SigY,2) - 2*pow(TauXY,2)))/2.)/2.;
      
      Solout[0]=J2;
      break;
    }
    case 21:
    {
      REAL I1 = SigX+SigY;
      Solout[0]=I1;
      break;
    }
    case 22:
      Solout[0] = 0.;
      break;
    case 23:
      // normal stress
      Solout[0] = SigX;
      Solout[1] = SigY;
      Solout[2] = SigZ;
      break;
    case 24:
      // shear stress
      Solout[0] = TauXY;
      Solout[1] = 0.;
      Solout[2] = 0.;
      break;
    case 25:
      Solout[0] = epsx;
      Solout[1] = epsy;
      Solout[2] = epsz;
      break;
    case 26:
      Solout[0] = epsxy;
      Solout[1] = 0.;
      Solout[2] = 0.;
      break;
    case 27:
      Solout[0] = 0.;
      Solout[1] = 0.;
      Solout[2] = 0.;
      break;
		default:
			TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
			break;
	}
}

void TPZNLElasticityMaterial::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
                                     TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, 
                                     TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
	values[0] = 0.;
	TPZVec<REAL> sigma(3,0.),sigma_exact(3,0.);
	REAL sigx,sigy,sigxy,gamma;
	TPZFMatrix<STATE> du(dudx.Rows(),dudx.Cols());
	du(0,0) = dudx(0,0)*axes(0,0)+dudx(1,0)*axes(1,0);
	du(1,0) = dudx(0,0)*axes(0,1)+dudx(1,0)*axes(1,1);
	du(0,1) = dudx(0,1)*axes(0,0)+dudx(1,1)*axes(1,0);
	du(1,1) = dudx(0,1)*axes(0,1)+dudx(1,1)*axes(1,1);
	
	//tens�s aproximadas : uma forma
	gamma = du(1,0)+du(0,1);
	sigma[0] = fPreStressXX+fEover1MinNu2*(du(0,0)+fnu*du(1,1));
	sigma[1] = fPreStressYY+fEover1MinNu2*(fnu*du(0,0)+du(1,1));
	sigma[2] = fPreStressXY+fE*0.5/(1.+fnu)*gamma;
	
	//exata
	gamma = du_exact(1,0)+du_exact(0,1);
	sigma_exact[0] = fEover1MinNu2*(du_exact(0,0)+fnu*du_exact(1,1));
	sigma_exact[1] = fEover1MinNu2*(fnu*du_exact(0,0)+du_exact(1,1));
	sigma_exact[2] = fE*0.5/(1.+fnu)*gamma;
	sigx  = (sigma[0] - sigma_exact[0]);
	sigy  = (sigma[1] - sigma_exact[1]);
	sigxy = (sigma[2] - sigma_exact[2]);
  
	//values[0] = calculo do erro estimado em norma Energia
	values[0] = fE*(sigx*sigx + sigy*sigy + 2*fnu*sigx*sigy)/(1-fnu*fnu);
	values[0] = (values[0] + .5*fE*sigxy*sigxy/(1+fnu));
	
	//values[1] : erro em norma L2 em tens�s
	//values[1] = sigx*sigx + sigy*sigy + sigxy*sigxy;
	
	//values[1] : erro em norma L2 em deslocamentos
	values[1] = pow((REAL)fabs(u[0] - u_exact[0]),(REAL)2.0)+pow((REAL)fabs(u[1] - u_exact[1]),(REAL)2.0);
	
	//values[2] : erro estimado na norma H1
  REAL SemiH1 =0.;
  for(int i = 0; i < 2; i++) for(int j = 0; j < 2; j++) SemiH1 += (du(i,j) - du_exact(i,j)) * (du(i,j) - du_exact(i,j));
	values[2] = values[1] + SemiH1;
}


TPZNLElasticityMaterial::TPZNLElasticityMaterial(const TPZNLElasticityMaterial &copy) :
TPZMaterial(copy),
fE(copy.fE),
fnu(copy.fnu),
fEover21PlusNu(copy.fEover21PlusNu),
fEover1MinNu2(copy.fEover1MinNu2),
fPreStressXX(copy.fPreStressXX),
fPreStressYY(copy.fPreStressYY),
fPreStressXY(copy.fPreStressXY),
fPreStressZZ(copy.fPreStressZZ)
{
	ff[0]=copy.ff[0];
	ff[1]=copy.ff[1];
	ff[2]=copy.ff[2];
	fPlaneStress = copy.fPlaneStress;
  // Added by Philippe 2012
  fPostProcIndex = copy.fPostProcIndex;
  
}


int TPZNLElasticityMaterial::ClassId() const{
    return Hash("TPZNLElasticityMaterial") ^ TPZMaterial::ClassId() << 1;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZNLElasticityMaterial>;
#endif

void TPZNLElasticityMaterial::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
	buf.Read(&fE,1);
	buf.Read(&fnu,1);
	buf.Read(&fEover21PlusNu,1);
	buf.Read(&fEover1MinNu2,1);
	buf.Read(&fPreStressXX,1);
	buf.Read(&fPreStressYY,1);
	buf.Read(&fPreStressXY,1);
	buf.Read(&fPreStressZZ,1);
	
	buf.Read(ff,3);
	buf.Read(&fPlaneStress,1);
  buf.Read(&fPostProcIndex);
	
}

void TPZNLElasticityMaterial::Write(TPZStream &buf, int withclassid) const
{
	TPZMaterial::Write(buf,withclassid);
	buf.Write(&fE,1);
	buf.Write(&fnu,1);
	buf.Write(&fEover21PlusNu,1);
	buf.Write(&fEover1MinNu2,1);
	buf.Write(&fPreStressXX,1);
	buf.Write(&fPreStressYY,1);
	buf.Write(&fPreStressXY,1);
	buf.Write(&fPreStressZZ,1);
	
	buf.Write(ff,3);
	buf.Write(&fPlaneStress,1);
  buf.Write(&fPostProcIndex);
	
}

