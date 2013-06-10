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
static LoggerPtr logdata(Logger::getLogger("pz.material.elastpressure"));
#endif

TPZNLFluidStructure2d::EState TPZNLFluidStructure2d::gState = ECurrentState;

TPZNLFluidStructure2d::TPZNLFluidStructure2d() : TPZDiscontinuousGalerkin(){
    
    fmatId = 0;
    fDim = 2;
    fPlaneStress = 1.;
    
    ff.resize(fDim);
	ff[0] = 0.;
	ff[1] = 0.;
}

TPZNLFluidStructure2d::TPZNLFluidStructure2d(int matid, int dim): TPZDiscontinuousGalerkin(matid){
    
    fmatId = matid;
    fDim = dim;
    fPlaneStress = 1.;
    
    ff.resize(fDim);
	ff[0] = 0.;
	ff[1] = 0.;
}

TPZNLFluidStructure2d::~TPZNLFluidStructure2d(){
}


int TPZNLFluidStructure2d::NStateVariables() {
	return 1;
}

void TPZNLFluidStructure2d::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "\t E   = " << fInputData->E() << std::endl;
	out << "\t nu   = " << fInputData->Poisson()   << std::endl;
	out << "\t Forcing function F   = " << ff[0] << ' ' << ff[1]   << std::endl;
    out << "altura da fratura fHf "<< fInputData->Hf() << std::endl;
	out << "Viscosidade do fluido fvisc "<< fInputData->Visc() << std::endl;
    out << "Carter fCl " << fInputData->Cl() << std::endl;
    out << "Pressao estatica fPe " << fInputData->Pe() << std::endl;
    out << "Pressao de referencia (Carter) fPref " << fInputData->Pref() << std::endl;
    out << "Spurt loss fvsp " << fInputData->vsp() << std::endl;
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
	if (nref != 2 ) {
		std::cout << " Error.!! the size of datavec is equal to two\n";
        std::cout << " datavec[0]->elasticity and datavec[1]->pressure\n";
		DebugStop();
	}
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
    if(shapetype!=datavec[0].EVecShape && datavec[0].phi.Cols()!=0){
        
        std::cout << " The space to elasticity problem must be reduced space.\n";
		DebugStop();
    }
    
    //Calculate the matrix contribution for elastic problem.
    TPZFMatrix<REAL> &dphi_u = datavec[0].dphix;
	TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	TPZFMatrix<REAL> &axes=datavec[0].axes;
    
    TPZManVector<REAL,3> sol_u=datavec[0].sol[0];
	TPZFMatrix<REAL> &dsol_u=datavec[0].dsol[0];
    
	int phcu,efcu;
	phcu = phi_u.Cols();
	efcu = ef.Cols();
	
	if(fForcingFunction) {// phi(in, 0) :  node in associated forcing function
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
    REAL young = fInputData->E();
    REAL poisson = fInputData->Poisson();

    REAL fEover1MinNu2 = young/(1-poisson*poisson);  ///4G(lamb+G)/(lamb+2G)
    REAL fEover21PlusNu = 2.*young/(2.*(1.+poisson));/*fE/(2.*(1+fnu));*/ ///2G=2mi
	REAL nu1 = 1 - poisson;//(1-nu)
	REAL nu2 = (1.-2.*poisson)/2.;
	REAL F = young/((1.+poisson)*(1.-2.*poisson));
    
	for(int in = 0; in < phcu; in++)
    {
		dphix_i(0,0) = dphi_u(0,in)*axes(0,0)+dphi_u(1,in)*axes(1,0);
		dphix_i(1,0) = dphi_u(0,in)*axes(0,1)+dphi_u(1,in)*axes(1,1);
		dphiy_i(0,0) = dphi_u(2,in)*axes(0,0)+dphi_u(3,in)*axes(1,0);
		dphiy_i(1,0) = dphi_u(2,in)*axes(0,1)+dphi_u(3,in)*axes(1,1);
		
        //Residuo
        for (int col = 0; col < efcu; col++)
        {
            //termo f*u
            ef(in,col) += weight*(ff[0]*phi_u(0, in) + ff[1]*phi_u(1, in));
            
            //termos da matriz
            dsolx_j(0,0) = dsol_u(0,0)*axes(0,0)+dsol_u(1,0)*axes(1,0);
            dsolx_j(1,0) = dsol_u(0,0)*axes(0,1)+dsol_u(1,0)*axes(1,1);
            dsoly_j(0,0) = dsol_u(0,1)*axes(0,0)+dsol_u(1,1)*axes(1,0);
            dsoly_j(1,0) = dsol_u(0,1)*axes(0,1)+dsol_u(1,1)*axes(1,1);
            
			
			if (fPlaneStress != 1){/* Plain Strain State */
                ef(in,col) += (-1.)*weight*(nu1*dphix_i(0,0)*dsolx_j(0,0) + nu2*dphix_i(1,0)*dsolx_j(1,0) +
                                     
                                     poisson*dphix_i(0,0)*dsoly_j(1,0) + nu2*dphix_i(1,0)*dsoly_j(0,0) +
                                     
                                     poisson*dphiy_i(1,0)*dsolx_j(0,0) + nu2*dphiy_i(0,0)*dsolx_j(1,0) +
                                     
                                     nu1*dphiy_i(1,0)*dsoly_j(1,0) + nu2*dphiy_i(0,0)*dsoly_j(0,0))*F;
			}
			else{/* Plain stress state */
                ef(in,col) += (-1.)*weight*(fEover1MinNu2*dphix_i(0,0)*dsolx_j(0,0) + fEover21PlusNu*dphix_i(1,0)*dsolx_j(1,0) +
                                     
                                     fEover1MinNu2*dphix_i(0,0)*dsoly_j(1,0) + fEover21PlusNu*dphix_i(1,0)*dsoly_j(0,0) +
                                     
                                     fEover1MinNu2*dphiy_i(1,0)*dsolx_j(0,0) + fEover21PlusNu*dphiy_i(0,0)*dsolx_j(1,0) +
                                     
                                     fEover1MinNu2*dphiy_i(1,0)*dsoly_j(1,0) + fEover21PlusNu*dphiy_i(0,0)*dsoly_j(0,0));
            }
		}//fim para o residuo
        
        //Matriz Tangente (Jacobiana)
        for( int jn = 0; jn < phcu; jn++ ) {
            
            dphix_j(0,0) = dphi_u(0,jn)*axes(0,0)+dphi_u(1,jn)*axes(1,0);
            dphix_j(1,0) = dphi_u(0,jn)*axes(0,1)+dphi_u(1,jn)*axes(1,1);
            dphiy_j(0,0) = dphi_u(2,jn)*axes(0,0)+dphi_u(3,jn)*axes(1,0);
            dphiy_j(1,0) = dphi_u(2,jn)*axes(0,1)+dphi_u(3,jn)*axes(1,1);
			
			
			if (fPlaneStress != 1){/* Plain Strain State */
                ek(in,jn) += weight*(nu1*dphix_i(0,0)*dphix_j(0,0) + nu2*dphix_i(1,0)*dphix_j(1,0) +
                                     
                                     poisson*dphix_i(0,0)*dphiy_j(1,0) + nu2*dphix_i(1,0)*dphiy_j(0,0) +
                                     
                                     poisson*dphiy_i(1,0)*dphix_j(0,0) + nu2*dphiy_i(0,0)*dphix_j(1,0) +
                                     
                                     nu1*dphiy_i(1,0)*dphiy_j(1,0) + nu2*dphiy_i(0,0)*dphiy_j(0,0))*F;
			}
			else{/* Plain stress state */
                ek(in,jn) += weight*(fEover1MinNu2*dphix_i(0,0)*dphix_j(0,0) + fEover21PlusNu*dphix_i(1,0)*dphix_j(1,0) +
                                     
                                     fEover1MinNu2*dphix_i(0,0)*dphiy_j(1,0) + fEover21PlusNu*dphix_i(1,0)*dphiy_j(0,0) +
                                     
                                     fEover1MinNu2*dphiy_i(1,0)*dphix_j(0,0) + fEover21PlusNu*dphiy_i(0,0)*dphix_j(1,0) +
                                     
                                     fEover1MinNu2*dphiy_i(1,0)*dphiy_j(1,0) + fEover21PlusNu*dphiy_i(0,0)*dphiy_j(0,0));
            }
		}//fim para a tangente

	}
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//		std::stringstream sout;
//		ek.Print("ek_reduced = ",sout,EMathematicaInput);
//		ef.Print("ef_reduced = ",sout,EMathematicaInput);
//		LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
    
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
    TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
    
    int phrp = phi_p.Rows();
    int phiuCols = phi_u.Cols();
    if(phiuCols != 1)
    {
        DebugStop();//Nao fiz!!!
    }
    int nPhiU = 1;
    
    REAL uy = sol_u[1];
    REAL w = 2.*uy;
    
    REAL visc = fInputData->Visc();
    REAL deltaT = fInputData->deltaT();
    
	if(gState == ECurrentState) //current state (n+1): Matrix stiffnes
    {
        REAL actQl = this->QlFVl(datavec[1].gelElId, sol_p[0]);
        REAL actdQldp = this->dQlFVl(datavec[1].gelElId, sol_p[0]);
        
        for(int in = 0; in < phrp; in++)
        {
            //----Residuo----
            //termo (wˆ3/(12*mi))*gradP * gradVp
            ef(nPhiU+in,0) += (-1.) * weight * (w*w*w/(12.*visc)) * dsol_p(0,0) * dphi_p(0,in);
            
            //termo w/deltaT * Vp
            ef(nPhiU+in,0) += (-1.) * weight * w/deltaT * phi_p(in,0);
            
            //termo 2Ql * Vp
            ef(nPhiU+in,0) += (-1.) * weight * (2.*actQl) * phi_p(in,0);

            
            //------Matriz tangente-----
            for(int jn = 0; jn < nPhiU; jn++)
            {
                //termo D[ (wˆ3/(12*mi))*gradP * gradVp , w ]
                ek(nPhiU+in, jn) += (+1.) * weight * ( 3.*w*w/(12.*visc) * (2.*phi_u(1,jn)) ) * dsol_p(0,0) * dphi_p(0,in);

                //termo D[ w/deltaT * Vp , w ]
                ek(nPhiU+in, jn) += (+1.) * weight * ( 2./deltaT * phi_u(1,jn) ) * phi_p(in,0);
            }
            
            for(int jn = 0; jn < phrp; jn++)
            {
                //termo D[ (wˆ3/(12*mi))*gradP * gradVp , p ]
                ek(nPhiU+in, nPhiU+jn) += (+1.) * weight * (w*w*w/(12.*visc)) * dphi_p(0,in) * dphi_p(0,jn);
                
                //termo D[ 2Ql * Vp , p]
                ek(nPhiU+in, nPhiU+jn) += (+1.) * weight * (2.*actdQldp) * phi_p(in,0) * phi_p(jn,0);
            }
        }
    }
    
    //Last state (n): Matrix mass
	if(gState == ELastState)
    {
        for(int in = 0; in < phrp; in++)
        {            
            //termo w/deltaT * Vp
            ef(nPhiU+in,0) += (+1.) * weight * w/deltaT * phi_p(in,0);
        }
    }
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//		std::stringstream sout;
//		ek.Print("ek_reduced = ",sout,EMathematicaInput);
//		ef.Print("ef_reduced = ",sout,EMathematicaInput);
//		LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
}

void TPZNLFluidStructure2d::ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc){
    
    //  if(gState == ELastState) return;
    
    TPZFMatrix<REAL> &phi_u = datavec[0].phi;
    TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	const REAL big  = TPZMaterial::gBigNumber;
	int phc = phi_u.Cols();
    
    for(int in = 0 ; in < phc; in++) {
        for (int il = 0; il <fNumLoadCases; il++)
        {
            //termo big*u*v do vetor de carga
            TPZFNMatrix<3,STATE> v2 = bc.Val2(il);
            ef(in,il) += big*(v2(0,il)*phi_u(0,in) + v2(1,il)*phi_u(1,in))*weight;
            
            //termo big*u*v da matriz
            ef(in,il) += (-1.)*big*(phi_u(0,in)*sol_u[0] + phi_u(1,in)*sol_u[1])*weight;
        }
        
        for (int jn = 0; jn < phc; jn++) {
            
            ek(in,jn) += big*(phi_u(0,in)*phi_u(0,jn) + phi_u(1,in)*phi_u(1,jn))*weight;
        }
    }
}

void TPZNLFluidStructure2d::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef,TPZBndCond &bc){
    
    if(gState == ELastState) return;
    REAL auxvar = 0.817*(1-fInputData->Poisson())*fInputData->Hf();
    REAL factor = 0.;//fG/auxvar;
    
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
            ef(in,il) += weight * factor * (phi_u(0,in)*sol_un*datavec[0].normal[0] + phi_u(1,in)*sol_un*datavec[0].normal[1])
                         - weight * (phi_u(0,in)*sol_p[0]*datavec[0].normal[0] + phi_u(1,in)*sol_p[0]*datavec[0].normal[1]);
        }
        
        //----- matriz tangente -----
        //termo k(phi_ix*phi_jun*nx + phi_iy*phi_jun*ny)
        for (int jn = 0; jn <nc_u; jn++)
        {
            REAL phi_jun = phi_u(0,jn)*datavec[0].normal[0] + phi_u(1,jn)*datavec[0].normal[1];
            ek(in,jn) += (-1.)*weight*factor*(phi_u(0,in)*phi_jun*datavec[0].normal[0] + phi_u(1,in)*phi_jun*datavec[0].normal[1]);
        }
        
        for (int jp = 0; jp <phrp; jp++)
        {
            ek(in,jp+nc_u) += weight*(phi_u(0,in)*phi_p(jp,0)*datavec[0].normal[0] + phi_u(1,in)*phi_p(jp,0)*datavec[0].normal[1]);
        }
    }
        
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//		std::stringstream sout;
//		ek.Print("ek_reduced = ",sout,EMathematicaInput);
//		ef.Print("ef_reduced = ",sout,EMathematicaInput);
//		LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
    
}

void TPZNLFluidStructure2d::ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc){
    
    if(gState == ELastState) return;
    
    TPZFMatrix<REAL> &phi_u = datavec[0].phi;
    TPZManVector<REAL,3> sol_u =datavec[0].sol[0];
	int phc = phi_u.Cols();
    
    for(int in = 0 ; in < phc; in++)
    {
        for (int il = 0; il <fNumLoadCases; il++)
        {
            TPZFNMatrix<3,STATE> v2 = bc.Val2(il);
            ef(in,il)+= weight*(v2(0,il)*phi_u(0,in) + v2(1,il)*phi_u(1,in));
            
            //termo da matriz
            ef(in,il) += (-1.)*(bc.Val1()(0,0)*phi_u(0,in)*sol_u[0]*weight
            
            + bc.Val1()(1,0)*phi_u(1,in)*sol_u[0]*weight
            
            + bc.Val1()(0,1)*phi_u(0,in)*sol_u[1]*weight
            
            + bc.Val1()(1,1)*phi_u(1,in)*sol_u[1]*weight);
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

void TPZNLFluidStructure2d::ApplyDirichlet_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc){
    
    if(gState == ELastState) return;
    
    TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	int c_u = phi_u.Cols();
    
    TPZFMatrix<REAL> &phi_p = datavec[1].phi;
    TPZManVector<REAL,3> sol_p = datavec[1].sol[0];
    int phrp = phi_p.Rows();
    
    for(int in = 0; in<phrp; in++)
    {
        //termo do vetor de carga
        ef(in+c_u,0)+=gBigNumber*bc.Val2()(2,0)*phi_p(in,0)*weight;
        
        //termo da matriz
        ef(in+c_u,0)+=(-1.)*gBigNumber*sol_p[0]*phi_p(in,0)*weight;
        
        for(int jn = 0; jn<phrp; jn++)
        {
            ek(in+c_u,jn+c_u)+=gBigNumber*phi_p(in,0)*phi_p(jn,0)*weight;
        }
    }
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//		std::stringstream sout;
//		phi_p.Print("phi = ",sout,EMathematicaInput);
//        bc.Val2().Print("val2 = ",sout);
//        sout << datavec[1].sol[0][0];
//		LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
}

void TPZNLFluidStructure2d::ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc){
    
    if(gState == ELastState) return;
    
    TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	int c_u = phi_u.Cols();
    
    TPZFMatrix<REAL> &phi_p = datavec[1].phi;
    int  phrp = phi_p.Rows();

    REAL Qinj = bc.Val2()(2,0);
    for(int in = 0; in < phrp; in++)
    {
        ef(in+c_u,0) += (-1.) * weight * Qinj * phi_p(in,0);
    }
}

void TPZNLFluidStructure2d::ApplyMixed_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc){
    
    if(gState == ELastState) return;
    
    TPZFMatrix<REAL> &phi_u = datavec[0].phi;
	int c_u = phi_u.Cols();
    
    TPZFMatrix<REAL> &phi_p = datavec[1].phi;
    TPZManVector<REAL,3> sol_p = datavec[1].sol[0];
    int phrp = phi_p.Rows();
    
    for(int in=0; in<phrp; in++)
    {
        //termo do vetor de carga
        ef(in+c_u,0)+=bc.Val2()(2,0)*phi_p(in,0)*weight;
        
        //termo da matriz
        ef(in+c_u,0) += (-1.)*bc.Val1()(2,0)*phi_p(in,0)*sol_p[0]*weight;
        
        for(int jn = 0; jn<phrp; jn++)
        {
            ek(in+c_u,jn+c_u) += bc.Val1()(2,0)*phi_p(in,0)*phi_p(jn,0)*weight;
        }
    }
}

void TPZNLFluidStructure2d::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
    if(shapetype!=datavec[0].EVecShape && datavec[0].phi.Cols()!=0)
    {
        std::cout << " The space to elasticity problem must be reduced space.\n";
		DebugStop();
    }
	
	switch (bc.Type())
    {
		case 0 :			// Dirichlet condition in  both elasticity and pressure equations
			
            ApplyDirichlet_U(datavec, weight, ek,ef,bc);
            ApplyDirichlet_P(datavec, weight, ek,ef,bc);
            break;
			
		case 1 :			// Neumann condition in  both elasticity and pressure equations
            ApplyNeumann_U(datavec, weight,ek,ef,bc);
            ContributePressure(datavec, weight, ek, ef);
            break;
			
		case 2 :            // Mixed condition in  both elasticity and pressure equations
            ApplyMixed_U(datavec, weight, ek,ef,bc);
            ApplyMixed_P(datavec, weight, ek,ef,bc);
            break;
            
        case 10:        // Neumann condition only on the elasticity equation
        {
            // Calculate the matrix contribution for pressure
            ApplyNeumann_U(datavec, weight,ek, ef,bc);
            break;
        }
        case 11:        // Dirichlet condition only on the elasticity equation
        {
            // Calculate the matrix contribution for pressure
             ApplyDirichlet_U(datavec, weight, ek,ef,bc);
            break;
        }
            
        case 20:        // Mixed condition only on the elasticity equation
            ApplyMixed_U(datavec, weight, ek,ef,bc);
            break;
            
        case 21:        // Neumann condition only on the pressure equation pressure
            ApplyNeumann_P(datavec, weight, ek,ef,bc);
            break;
            
        case 22:        // Dirichlet condition only on the pressure equation pressure
            ApplyDirichlet_P(datavec, weight, ek,ef,bc);
            break;
    }
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//		std::stringstream sout;
//		ek.Print("ek_reduced = ",sout,EMathematicaInput);
//		ef.Print("ef_reduced = ",sout,EMathematicaInput);
//		LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
    
}

int TPZNLFluidStructure2d::VariableIndex(const std::string &name){
    if(!strcmp("Pressure",name.c_str()))        return  1;
    if(!strcmp("MinusKGradP",name.c_str()))     return  2;
    if(!strcmp("DisplacementX",name.c_str()))  return 3;
	if(!strcmp("DisplacementY",name.c_str()))  return 4;
    if(!strcmp("SigmaX",name.c_str()))        return  5;
	if(!strcmp("SigmaY",name.c_str()))        return  6;
    
    if(!strcmp("Displacement",name.c_str()))        return  7;
    
    if(!strcmp("W",name.c_str()))        return  8;
    
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
    
    REAL young = fInputData->E();
    REAL poisson = fInputData->Poisson();
    REAL sigmaConf = fInputData->SigmaConf();
    REAL Hf = fInputData->Hf();
    
    if(var == 1)
    {
        if(!datavec[1].phi) return;
		Solout[0] = SolP[0];
		return;
	}//var1
    
    if (var == 2)
    {
        REAL fG = 1.;//AQUICAJU!!!
        if(!datavec[1].phi) return;
        REAL un = 0.817*(1-poisson)*(SolP[0]-sigmaConf)*Hf/fG;
        REAL factor = 0.;//(un*un*un)/(12.*fvisc);
        
		int id;
		TPZFNMatrix<9,REAL> dsoldx;
		TPZAxesTools<REAL>::Axes2XYZ(DSolP, dsoldx, axesP);
		for(id=0 ; id<1; id++){
			Solout[id] = -1.*factor*dsoldx(id,0);
		}
		return;
	}//var2
    
    //function (state variable ux)
	if(var == 3)
    {
		Solout[0] = SolU[0];
		return;
	}//var3
	
	//function (state variable uy)
	if(var == 4)
    {
		Solout[0] = SolU[1];
		return;
	}//var4
    
    if(var == 5 || var == 6)
    {
        REAL DSolxy[2][2];
        REAL SigX, SigY, epsx, epsy;
        
        DSolxy[0][0] = DSolU(0,0)*axesU(0,0)+DSolU(1,0)*axesU(1,0);
        DSolxy[1][1] = DSolU(0,1)*axesU(0,1)+DSolU(1,1)*axesU(1,1);
    //    epsx = DSolxy[0][0];// du/dx
    //	epsy = DSolxy[1][1];// dv/dy
        
        epsx = DSolU(0,0);// du/dx
        epsy = DSolU(1,1);// dv/dy
        REAL Gmodule = young/(1-poisson*poisson);
        
        
        if (this->fPlaneStress){
            SigX = Gmodule*(epsx+poisson*epsy);
            SigY = Gmodule*(poisson*epsx+epsy);
        }
        else{
            SigX = young/((1.-2.*poisson)*(1.+poisson))*((1.-poisson)*epsx+poisson*epsy);
            SigY = young/((1.-2.*poisson)*(1.+poisson))*(poisson*epsx+(1.-poisson)*epsy);
        }
        
        if(var == 5) {
            Solout[0] = SigX;
            return;
        }//var5
        
        if(var == 6) {
            Solout[0] = SigY;
            return;
        }//var6
    }
    if(var == 7){
		Solout[0] = SolU[0];
		Solout[1] = SolU[1];
		return;
	}//var1
    
    if(var == 8)
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

void TPZNLFluidStructure2d::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
        datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNormal = true;
	}
}



////////////////////////////////////////////////////////////////// Leakoff

REAL TPZNLFluidStructure2d::VlFtau(REAL pfrac, REAL tau)
{
    REAL Cl = fInputData->Cl();
    REAL sigmaConf = fInputData->SigmaConf();
    REAL Pe = fInputData->Pe();
    REAL Pref = fInputData->Pref();
    REAL vsp = fInputData->vsp();
    
    REAL Clcorr = Cl * sqrt((pfrac + sigmaConf - Pe)/Pref);
    REAL Vl = 2. * Clcorr * sqrt(tau) + vsp;
    
    return Vl;
}

REAL TPZNLFluidStructure2d::FictitiousTime(REAL VlAcum, REAL pfrac)
{
    REAL Cl = fInputData->Cl();
    REAL sigmaConf = fInputData->SigmaConf();
    REAL Pe = fInputData->Pe();
    REAL Pref = fInputData->Pref();
    REAL vsp = fInputData->vsp();
    
    REAL tStar = 0.;
    if(VlAcum > vsp)
    {
        REAL Clcorr = Cl * sqrt((pfrac + sigmaConf - Pe)/Pref);
        tStar = (VlAcum - vsp)*(VlAcum - vsp)/( (2. * Clcorr) * (2. * Clcorr) );
    }
    
    return tStar;
}

REAL TPZNLFluidStructure2d::QlFVl(int gelId, REAL pfrac)
{
    std::map<int,REAL>::iterator it = fGelId_vl.find(gelId);
    if(it == fGelId_vl.end())
    {
        fGelId_vl[gelId] = 0.;//Nao coloque vsp! Eh ZERO mesmo!
        it = fGelId_vl.find(gelId);
    }
    REAL VlAcum = it->second;
    
    REAL deltaT = fInputData->deltaT();
    
    REAL tStar = FictitiousTime(VlAcum, pfrac);
    REAL Vlnext = VlFtau(pfrac, tStar + deltaT);
    REAL Ql = (Vlnext - VlAcum)/deltaT;

    return Ql;
}

REAL TPZNLFluidStructure2d::dQlFVl(int gelId, REAL pfrac)
{
    std::map<int,REAL>::iterator it = fGelId_vl.find(gelId);
    if(it == fGelId_vl.end())
    {
        fGelId_vl[gelId] = 0.;
        it = fGelId_vl.find(gelId);
    }
    REAL VlAcum = it->second;
    
    REAL deltaPfrac = fabs(pfrac/10000.);
    if(deltaPfrac < 1.E-10)
    {
        deltaPfrac = 1.E-10;
    }
    else if(deltaPfrac > 1.E-3)
    {
        deltaPfrac = 1.E-3;
    }
    
    REAL deltaT = fInputData->deltaT();
    /////////////////////////////////////////////////Ql maior
    REAL pfracUP = pfrac + deltaPfrac;
    REAL tStar1 = FictitiousTime(VlAcum, pfracUP);
    REAL Vlnext1 = VlFtau(pfracUP, tStar1 + deltaT);
    REAL Ql1 = (Vlnext1 - VlAcum )/deltaT;
    //...
    
    /////////////////////////////////////////////////Ql menor
    REAL pfracDOWN = pfrac - deltaPfrac;
    REAL tStar0 = FictitiousTime(VlAcum, pfracDOWN);
    REAL Vlnext0 = VlFtau(pfracDOWN, tStar0 + deltaT);
    REAL Ql0 = (Vlnext0 - VlAcum)/deltaT;
    //...
    
    REAL dQldpfrac = (Ql1-Ql0)/(2.*deltaPfrac);

    return dQldpfrac;
}

std::ofstream outVl("vl.txt");
void TPZNLFluidStructure2d::UpdateLeakoff(TPZCompMesh * cmesh)
{
    std::map<int,REAL>::iterator it;
    
    int outVlCount = 0;
    for(int i = 0;  i < cmesh->ElementVec().NElements(); i++)
    {
        ///////////////////////
        TPZCompEl * cel = cmesh->ElementVec()[i];
        
        #ifdef DEBUG
        if(!cel)
        {
            DebugStop();
        }
        #endif
        
        TPZGeoEl * gel = cel->Reference();
        
        if(gel->Dimension() != 1)
        {
            continue;
        }
        
        TPZInterpolatedElement * sp = dynamic_cast <TPZInterpolatedElement*> (cel);
        if(!sp)
        {
            continue;
        }
        
        #ifdef DEBUG
        if(!gel || gel->MaterialId() != globPressureMatId)
        {
            DebugStop();
        }
        #endif
        
        TPZVec<REAL> qsi(1,0.);
        cel->Reference()->CenterPoint(cel->Reference()->NSides()-1, qsi);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        
        sp->ComputeShape(qsi, data);
        sp->ComputeSolution(qsi, data);

        REAL pfrac = data.sol[0][0];
        ///////////////////////
        
        it = fGelId_vl.find(gel->Id());
        
        #ifdef DEBUG
        if(it == fGelId_vl.end())
        {
            DebugStop();
        }
        #endif
        
        REAL deltaT = fInputData->deltaT();
        
        REAL VlAcum = it->second;
        REAL tStar = FictitiousTime(VlAcum, pfrac);
        REAL Vlnext = VlFtau(pfrac, tStar + deltaT);
        
        it->second = Vlnext;
        if(outVlCount == 0)
        {
            outVl << VlAcum << "\n";
        }
        outVlCount++;
    }
    
    #ifdef DEBUG
    if(outVlCount != fGelId_vl.size())
    {
        DebugStop();
    }
    #endif
}
