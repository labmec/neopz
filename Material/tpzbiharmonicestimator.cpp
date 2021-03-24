/**
 * @file
 * @brief Contains implementations of the TPZBiharmonicEstimator methods.
 */

#include "tpzbiharmonicestimator.h"
#include "pzbiharmonic.h"
#include "pzbndcond.h"

#include <cmath>
using namespace std;

TPZBiharmonicEstimator::TPZBiharmonicEstimator(int nummat, STATE f):
TPZRegisterClassId(&TPZBiharmonicEstimator::ClassId), TPZBiharmonic( nummat,  f)
{
	this->SetExactSolutions(NULL, NULL);
}

TPZBiharmonicEstimator::~TPZBiharmonicEstimator()
{
}

void TPZBiharmonicEstimator::SetExactSolutions(
											   void (*fp)(TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv),
											   void (*fd)(TPZVec<REAL> &locdual,TPZVec<STATE> &valdual,TPZFMatrix<STATE> &derivdual)){
	this->fPrimalExactSol = fp;
	this->fDualExactSol = fd;
}

//static STATE Pi = 4.*atan(1.);
// std::ofstream CoutPontual ( "solPontual.txt" );
void TPZBiharmonicEstimator::ContributeErrorsDual(TPZMaterialData &data,
												  REAL weight,
												  TPZVec<REAL> &nk){
	TPZVec<STATE> u_exactp(1);
	TPZFMatrix<STATE> du_exactp(9,1);
	TPZVec<STATE> u_exactd(1);
	TPZFMatrix<STATE> du_exactd(9,1);
	this->fPrimalExactSol(data.x,u_exactp,du_exactp);
	this->fDualExactSol(data.x,u_exactd,du_exactd);
	
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZManVector<REAL,3> x = data.x;
	TPZVec<STATE> sol  = data.sol[0];
	TPZFMatrix<STATE>   dsol = data.dsol[0];
	//   REAL Laplac = dsol(2,0);
	STATE LaplacLaplac = dsol(8,0);
	//   REAL DuDx = dsol(0,0);
	//   REAL DuDy = dsol(1,0);
	
    if (this->fForcingFunction)
	{
		TPZManVector<STATE, 1> result(1);
		this->fForcingFunction->Execute(x, result);
		this->fXf = result[0];
	}
	
	
	// std::cout.precision(16);
	// if(fabs(x[0]-0.5)<0.001 && fabs(x[1]-0.5)<0.001){
	// CoutPontual << "x: "<< x[0]<<","<<x[1]<<std::endl; 
	// CoutPontual << "sol: "<< sol[0]<<std::endl;
	// }
	STATE f = this->fXf;
	
	STATE Z  = sol[2];
	STATE Z1 = sol[3];
	nk[0] += weight *( f -  LaplacLaplac )*(Z1-Z);
	nk[1] += weight *( f -  LaplacLaplac )*(u_exactd[0]-Z);//segundo a teoria deve dar o erro exato do funcional(a menos do erro na integracao numerica)!
}

void TPZBiharmonicEstimator::ContributeInterfaceErrorsDual(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
														   REAL weight,
														   TPZVec<STATE> &nkL, 
														   TPZVec<STATE> &nkR ){
	TPZManVector<REAL,3> normal = data.normal;
	TPZManVector<REAL,3> x = data.x;
	int LeftPOrder = dataleft.p;
	int RightPOrder = dataright.p;
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZVec<STATE> solL=dataleft.sol[0];
	TPZVec<STATE> solR=dataright.sol[0];
	TPZFMatrix<STATE> dsolL=dataleft.dsol[0];
	TPZFMatrix<STATE> dsolR=dataright.dsol[0];
	STATE faceSize=data.HSize;
	
	STATE alpha=gSigmaA*(pow((STATE)LeftPOrder,(STATE)gL_alpha)+pow((STATE)RightPOrder,(STATE)gL_alpha)) /(2.*pow(faceSize,(STATE)gM_alpha));
	STATE betta=gSigmaB*(pow((STATE)LeftPOrder,(STATE)gL_betta)+pow((STATE)RightPOrder,(STATE)gL_betta)) /(2.*pow(faceSize,(STATE)gM_betta));
	
	// solucoes exatas
	TPZVec<STATE> u_exactp(1);
	TPZFMatrix<STATE> du_exactp(9,1);
	TPZVec<STATE> u_exactd(1);
	TPZFMatrix<STATE> du_exactd(9,1);
	this->fPrimalExactSol(data.x,u_exactp,du_exactp);
	this->fDualExactSol(data.x,u_exactd,du_exactd);
	
	
	// u -> left e right
	STATE uL = solL[0];
	STATE uR = solR[0];
	STATE graduL[2];
	graduL[0] = dsolL(0,0);
	graduL[1] = dsolL(1,0);
	STATE graduR[2];
	graduR[0] = dsolR(0,0);
	graduR[1] = dsolR(1,0);
	STATE LaplacianoUL =dsolL(2,0);
	STATE LaplacianoUR =dsolR(2,0);
	STATE  GRADLaplacianoUL[2];
	GRADLaplacianoUL[0] =dsolL(3,0);
	GRADLaplacianoUL[1] =dsolL(4,0);
	STATE  GRADLaplacianoUR[2];
	GRADLaplacianoUR[0] =dsolR(3,0);
	GRADLaplacianoUR[1] =dsolR(4,0);
	
	// dualhp -> left e right
	STATE dualhpL = solL[2];
	STATE dualhpR = solR[2];
	STATE graddualhpL[2];
	graddualhpL[0] = dsolL(0,2);
	graddualhpL[1] = dsolL(1,2);
	STATE graddualhpR[2];
	graddualhpR[0] = dsolR(0,2);
	graddualhpR[1] = dsolR(1,2);
	STATE LaplacianodualhpL = dsolL(2,2);
	STATE LaplacianodualhpR = dsolR(2,2);
	STATE  GRADLaplacianodualhpL[2];
	GRADLaplacianodualhpL[0] = dsolL(3,2);
	GRADLaplacianodualhpL[1] = dsolL(4,2);
	STATE  GRADLaplacianodualhpR[2];
	GRADLaplacianodualhpR[0] = dsolR(3,2);
	GRADLaplacianodualhpR[1] = dsolR(4,2);
	
	// dualhpPlus -> left e right
	STATE dualhpPlusL = solL[3];
	STATE dualhpPlusR = solR[3];
	STATE graddualhpPlusL[2];
	graddualhpPlusL[0] = dsolL(0,3);
	graddualhpPlusL[1] = dsolL(1,3);
	STATE graddualhpPlusR[2];
	graddualhpPlusR[0] = dsolR(0,3);
	graddualhpPlusR[1] = dsolR(1,3);
	STATE LaplacianodualhpPlusL = dsolL(2,3);
	STATE LaplacianodualhpPlusR = dsolR(2,3);
	STATE  GRADLaplacianodualhpPlusL[2];
	GRADLaplacianodualhpPlusL[0] = dsolL(3,3);
	GRADLaplacianodualhpPlusL[1] = dsolL(4,3);
	STATE  GRADLaplacianodualhpPlusR[2];
	GRADLaplacianodualhpPlusR[0] = dsolR(3,3);
	GRADLaplacianodualhpPlusR[1] = dsolR(4,3);

	// sao seis termos com left e right
	
	// primeiro termo (s/ penalizacao)
	nkL[0]+=-1.0*weight*
	(LaplacianoUL-LaplacianoUR)*
	0.5*(normal[0]*(graddualhpPlusL[0]-graddualhpL[0])+normal[1]*(graddualhpPlusL[1]-graddualhpL[1]));
	nkR[0]+=-1.0*weight*
	(LaplacianoUL-LaplacianoUR)*
	(/*-*/0.5)*(normal[0]*(graddualhpPlusR[0]-graddualhpR[0])+normal[1]*(graddualhpPlusR[1]-graddualhpR[1]));
	// segundo termo (s/ penalizacao)
	nkL[0]+=weight*((GRADLaplacianoUL[0]-GRADLaplacianoUR[0])*normal[0]+(GRADLaplacianoUL[1]-GRADLaplacianoUR[1])*normal[1])*0.5*(dualhpPlusL-dualhpL);
	nkR[0]+=/*-*/weight*((GRADLaplacianoUL[0]-GRADLaplacianoUR[0])*normal[0]+(GRADLaplacianoUL[1]-GRADLaplacianoUR[1])*normal[1])*0.5*(dualhpPlusR-dualhpR);
	// terceiro termo (s/ penalizacao)
	nkL[0]+=-weight*(uL-uR)*0.5*(normal[0]*(GRADLaplacianodualhpPlusL[0]-GRADLaplacianodualhpL[0])+normal[1]*(GRADLaplacianodualhpPlusL[1]-GRADLaplacianodualhpL[1]));
	nkR[0]+=-weight*(uL-uR)*0.5*(/*-*/normal[0]*(GRADLaplacianodualhpPlusR[0]-GRADLaplacianodualhpR[0])/*-*/+normal[1]*(GRADLaplacianodualhpPlusR[1]-GRADLaplacianodualhpR[1]));
	// quarto termo ( s/ penalizacao)
	nkL[0]+=weight*((graduL[0]-graduR[0])*normal[0]+(graduL[1]-graduR[1])*normal[1])*0.5*(LaplacianodualhpPlusL-LaplacianodualhpL);
	nkR[0]+=/*-*/weight*((graduL[0]-graduR[0])*normal[0]+(graduL[1]-graduR[1])*normal[1])*0.5*(LaplacianodualhpPlusR-LaplacianodualhpR);
	
	// quinto termo
	nkL[0]+=-weight*alpha*(uL-uR)*(dualhpPlusL-dualhpL);
	nkR[0]+=-weight*alpha*(uL-uR)*(-1.)*(dualhpPlusR-dualhpR);
	//sexto termo
	nkL[0]+=-weight*betta*((graduL[0]-graduR[0])*normal[0]+(graduL[1]-graduR[1])*normal[1]) * (normal[0]*(graddualhpPlusL[0]-graddualhpL[0])+normal[1]*(graddualhpPlusL[1]-graddualhpL[1]));
	nkR[0]+=-weight*betta*((graduL[0]-graduR[0])*normal[0]+(graduL[1]-graduR[1])*normal[1])*
	(-1.)*(normal[0]*(graddualhpPlusR[0]-graddualhpR[0])+normal[1]*(graddualhpPlusR[1]-graddualhpR[1]));
	
	// com sol's  exatas
	
	// primeiro termo
	nkL[1]+=-weight*(LaplacianoUL-LaplacianoUR)*0.5*(normal[0]*(du_exactd[0]-graddualhpL[0])+normal[1]*(du_exactd[1]-graddualhpL[1]));
	nkR[1]+=-weight*(LaplacianoUL-LaplacianoUR)*0.5*(normal[0]*(du_exactd[0]-graddualhpR[0])+normal[1]*(du_exactd[1]-graddualhpR[1]));
	// segundo termo
	nkL[1]+=weight*((GRADLaplacianoUL[0]-GRADLaplacianoUR[0])*normal[0]+(GRADLaplacianoUL[1]-GRADLaplacianoUR[1])*normal[1])*0.5*(u_exactd[0]-dualhpL);
	nkR[1]+=weight*((GRADLaplacianoUL[0]-GRADLaplacianoUR[0])*normal[0]+(GRADLaplacianoUL[1]-GRADLaplacianoUR[1])*normal[1])*0.5*(u_exactd[0]-dualhpR);
	// terceiro termo
	nkL[1]+=-weight*(uL-uR)*0.5*(normal[0]*(du_exactd[3]-GRADLaplacianodualhpL[0])+normal[1]*(du_exactd[4]-GRADLaplacianodualhpL[1]));
	nkR[1]+=-weight*(uL-uR)*0.5*(normal[0]*(du_exactd[3]-GRADLaplacianodualhpR[0])+normal[1]*(du_exactd[4]-GRADLaplacianodualhpR[1]));
	// quarto termo
	nkL[1]+=weight*((graduL[0]-graduR[0])*normal[0]+(graduL[1]-graduR[1])*normal[1])*0.5*(du_exactd[2]-LaplacianodualhpL);
	nkR[1]+=weight*((graduL[0]-graduR[0])*normal[0]+(graduL[1]-graduR[1])*normal[1])*0.5*(du_exactd[2]-LaplacianodualhpR);
	// quinto termo
	nkL[1]+=-weight*alpha*(uL-uR)*(u_exactd[0]-dualhpL);
	nkR[1]+=-weight*alpha*(uL-uR)*(-1.)*(u_exactd[0]-dualhpR);
	//sexto termo
	nkL[1]+=-weight* /*betaL*/ ((graduL[0]-graduR[0])*normal[0]+(graduL[1]-graduR[1])*normal[1]) * (normal[0]*(du_exactd[0]-graddualhpL[0])+normal[1]*(du_exactd[1]-graddualhpL[1]));
	nkR[1]+=-weight* /*betaR*/ ((graduL[0]-graduR[0])*normal[0]+(graduL[1]-graduR[1])*normal[1]) * (-1.)*(normal[0]*(du_exactd[0]-graddualhpR[0])+normal[1]*(du_exactd[1]-graddualhpR[1]));//acertar betas!!!!!!!!
	
}

void TPZBiharmonicEstimator::ContributeInterfaceBCErrorsDual(TPZMaterialData &data, TPZMaterialData &dataleft,
															 REAL weight,
															 TPZVec<STATE> &nk,
															 TPZBndCond &bc){
	TPZFMatrix<REAL> dphi = data.dphix;
	TPZFMatrix<REAL> dphiL = dataleft.dphix;
	TPZFMatrix<REAL> phi = data.phi;
	TPZFMatrix<REAL> phiL = dataleft.phi;
	TPZManVector<REAL,3> normal = data.normal;
	TPZManVector<REAL,3> x = data.x;
	int LeftPOrder=dataleft.p;
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZVec<STATE> solL=dataleft.sol[0];
	TPZFMatrix<STATE> dsolL=dataleft.dsol[0];
	STATE faceSize=data.HSize;
	
	TPZVec<STATE> u_exactp(1);
	TPZFMatrix<STATE> du_exactp(9,1);
	TPZVec<STATE> u_exactd(1);
	TPZFMatrix<STATE> du_exactd(9,1);
	this->fDualExactSol(data.x,u_exactd,du_exactd);
	this->fPrimalExactSol(data.x,u_exactp,du_exactp);
	
	STATE alpha = gSigmaA*pow((STATE)LeftPOrder, (STATE)gL_alpha) /  pow(faceSize, (STATE)gM_alpha);
	STATE betta = gSigmaB*pow((STATE)LeftPOrder, (STATE)gL_betta) /  pow(faceSize, (STATE)gM_betta);
	
	STATE u = dataleft.sol[0][0];
	STATE gradu[2];
	gradu[0] = dataleft.dsol[0](0,0);
	gradu[1] = dataleft.dsol[0](1,0);
	
	STATE dualhp = dataleft.sol[0][2];//pq a dual hp está na terceira coluna
	STATE GRADdualhp[2];//
	STATE GRADLaplacdualhp[2];
	GRADdualhp[0] = dataleft.dsol[0](0,2);
	GRADdualhp[1] = dataleft.dsol[0](1,2);
	STATE Laplacdualhp = dataleft.dsol[0](2,2);
	GRADLaplacdualhp[0] = dataleft.dsol[0](3,2);
	GRADLaplacdualhp[1] = dataleft.dsol[0](4,2);
	
	STATE dualhpPlus = dataleft.sol[0][3];// 
	STATE GRADdualhpPlus[2];//
	STATE GRADLaplacdualhpPlus[2];//
	GRADdualhpPlus[0]= dataleft.dsol[0](0,3);//
	GRADdualhpPlus[1]= dataleft.dsol[0](1,3);//
	STATE LaplacdualhpPlus = dataleft.dsol[0](2,3);
	GRADLaplacdualhpPlus[0] = dataleft.dsol[0](3,3);
	GRADLaplacdualhpPlus[1] = dataleft.dsol[0](4,3);
	
	
	STATE gradZ[2];
	gradZ[0] = GRADdualhp[0]; gradZ[1] = GRADdualhp[1];
	STATE gradZPlus[2];
	gradZPlus[0] = GRADdualhpPlus[0]; gradZPlus[1] = GRADdualhpPlus[1];
	
	STATE ValBC_D;
	STATE ValBC_N;
	ValBC_D = bc.Val2()(0,0);// Dirichlet  
	ValBC_N = bc.Val2()(1,0);// Neumann
	
	STATE ResD = ValBC_D - u;
	STATE ResN = ValBC_N - (gradu[0]*normal[0]+gradu[1]*normal[1]);
	
	nk[0] += weight * ResD * ( (GRADLaplacdualhpPlus[0]*normal[0]+GRADLaplacdualhpPlus[1]*normal[1]) - (GRADLaplacdualhp[0]*normal[0]+GRADLaplacdualhp[1]*normal[1]) );
	nk[0] -= weight*ResN*(LaplacdualhpPlus-Laplacdualhp);
	nk[0] += weight*alpha*ResD*(dualhpPlus-dualhp);
	nk[0] += weight*betta*ResN*( (GRADdualhpPlus[0]*normal[0]+GRADdualhpPlus[1]*normal[1]) - (GRADdualhp[0]*normal[0]+GRADdualhp[1]*normal[1]) );
	
	
	nk[1] += weight * ResD * ( (du_exactd[3]*normal[0]+du_exactd[4]*normal[1]) - (GRADLaplacdualhp[0]*normal[0]+GRADLaplacdualhp[1]*normal[1]) );
	nk[1] -= weight*ResN*(du_exactd[2]-Laplacdualhp);
	nk[1] += weight*alpha*ResD*(u_exactd[0]-dualhp);
	nk[1] += weight*betta*ResN*( (du_exactd[0]*normal[0]+du_exactd[1]*normal[1]) - (GRADdualhp[0]*normal[0]+GRADdualhp[1]*normal[1]) );
}

/** @brief Initializing variable for computing L2 error */
STATE L2ErrorPrimal = 0.;
/** @brief Initializing variable for computing L2 error to dual form. */
STATE L2ErrorDual = 0.;

void TPZBiharmonicEstimator::ContributeErrorsSimple(TPZMaterialData &data,
													REAL weight,
													TPZVec<REAL> &nk){
	int phis = data.phi.Rows();
	int dphir = data.dphix.Rows();
	OrderSolution(data);
	TPZFNMatrix<100,REAL> axes(3,3,0.),jacinv(2,2,0.);
    TPZFNMatrix<200,STATE> ek(data.phi.Rows(),data.phi.Rows(),0.),ef(data.phi.Rows(),1,0.);
	this->Contribute(data,weight,ek,ef);
	data.phi.Resize(phis,1);
	data.dphix.Resize(dphir,phis);
	//       std:: cout<< "x=("<< data.x[0]<<","<< data.x[1]<<")"<< "  "<< "xsim=(" << 1.-data.x[0]<<","<< 1.-data.x[1]<< ")";
	//       ek.Print("Teste", std::cout, EFormatted); 
	
	nk.Resize(10);
	
	// L2ErrorPrimal+= ((dphi(1,3)-dphi(1,1))*(dphi(1,3)-dphi(1,1))+(dphi(0,3)-dphi(0,1))*(dphi(0,3)-dphi(0,1))+(phi(3,0)-phi(1,0))*(phi(3,0)-phi(1,0)))*weight;
	// L2ErrorDual  += ((dphi(1,2)-dphi(1,0))*(dphi(1,2)-dphi(1,0))+(dphi(0,2)-dphi(0,0))*(dphi(0,2)-dphi(0,0))+(phi(2,0)-phi(0,0))*(phi(2,0)-phi(0,0)))*weight;
	
	nk[0] += ek(0,0);//norma da energia do erro aprox. dual 
	nk[3] += ek(1,1);// norma da energia do erro aprox. primal
	nk[6] += ek(0,1);// estimador de erro (p e p+2)
	nk[1] += ek(2,2);// norma da energia do erro dual exato
	nk[4] += ek(3,3);// norma da energia do erro primal exato 
	nk[7] += ek(2,3);// estimador de erro (p e exato)
	nk[2] += ek(4,4);// norma da energia do erro dual exato da sol. enriquecida
	nk[5] += ek(5,5);// norma da energia do erro primal exato da sol. enriquecida
	nk[8] += ek(4,5);// estimador de erro (p+2 e exato)
	nk[9] += 0.; // indice de efetividade
	// std:: cout<< "  "<< nk[0]<<"  "<< nk[1]<<"  "<< nk[2]<< "  "<< nk[3]<< "  " <<nk[4]<< "  "<< nk[5]<< std::endl;
}



void TPZBiharmonicEstimator::ContributeInterfaceErrorsSimple(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
															 REAL weight,
															 TPZVec<STATE> &nkL, 
															 TPZVec<STATE> &nkR){
	//   return;
	OrderSolutionLeft(data,dataleft);
	OrderSolutionRight(data,dataright);
	int sz = dataleft.phi.Rows()+dataright.phi.Rows();
	int ss = dataleft.phi.Rows(); //inclui isso pq a matriz eh szXsz mas os blocos são sz/2 X sz/2
	TPZFNMatrix<100,REAL> axesleft(3,3,0.),axesright;
    TPZFNMatrix<100,STATE> ef(sz,1,0.),ek(sz,sz,0.);
	this->ContributeInterface(data,dataleft,dataright,weight,ek,ef);
	
	//       std:: cout<< "x=("<< data.x[0]<<","<< data.x[1]<<")"; //"  "<< "xsim=(" << 1.-data.x[0]<<","<< 1.-data.x[1]<< ")";
	//       ek.Print("Teste", std::cout, EFormatted); 
	nkL.Resize(10);
	nkR.Resize(10);
	
	nkL[0] += 0.5*(ek(0,0)+ek(ss,0)+ek(0,ss)+ek(ss,ss));   // norma da energia do erro aprox. dual-Left
	nkL[3] += 0.5*(ek(1,1)+ek(ss+1,1)+ek(1,1+ss)+ek(1+ss,1+ss)); // norma da energia do erro aprox. primal-Left
	nkL[6] += 0.5*(ek(0,1)+ek(ss,ss+1)+ek(1,ss)+ek(ss+1,0));   // estimador (p e p+2)-Left  // O erro está no segundo termo
	nkL[1] += 0.5*(ek(2,2)+ek(2+ss,2+ss)+ek(2+ss,2)+ek(2,2+ss)); // norma da energia do erro exato da solucao dual-Left
	nkL[4] += 0.5*(ek(3,3)+ek(3+ss,3+ss)+ek(3+ss,3)+ek(3,3+ss)); // norma da energia do erro exato da solucao primal-Left
	nkL[7] += 0.5*(ek(2,3)+ek(ss+2,ss+3)+ek(3,2+ss)+ek(3+ss,2)); // estimador(p e exata)-Left
	nkL[2] += 0.5*(ek(4,4)+ek(4+ss,4+ss)+ek(4+ss,4)+ek(4,4+ss));// norma da energia do erro dual-left exato da sol. enriquecida
	nkL[5] += 0.5*(ek(5,5)+ek(ss+5,ss+5)+ek(5+ss,5)+ek(5,5+ss));// norma da energia do erro primal-left exato da sol. enriquecida
	nkL[8] += 0.5*(ek(4,5)+ek(ss+4,ss+5)+ek(5,4+ss)+ek(5+ss,4)); // estimador(p+2 e exata)-Left
	nkL[9] += 0.;      // indíce de efetividade
	
	nkR = nkL;
	
	//       nkR[0] += ek(ss,ss)    +0.5*(ek(0,ss)+ek(ss,0));      // norma da energia do erro aprox. dual-Right
	//       nkR[3] += ek(1+ss,1+ss)+0.5*(ek(1,ss+1)+ek(ss+1,1));// norma da energia do erro aprox. dual-Right
	//       nkR[6] += ek(ss,ss+1)+(ek(ss+1,0) +ek(1,ss))/2.;    // estimador(p e p+2)-Right  // O erro está no segundo termo
	//       nkR[1] += ek(ss+2,ss+2)+0.5*(ek(2,2+ss)+ek(2+ss,2));// norma da energia do erro exato da sol. dual Right
	//       nkR[4] += ek(ss+3,ss+3)+0.5*(ek(3,3+ss)+ek(3+ss,3));// norma da energia do erro exato da sol. primal Right
	//       nkR[7] += ek(ss+2,ss+3)+(ek(3,2+ss)+ek(3+ss,2))/2.;// estimador (p e exato) - Right
	//       nkR[2] += ek(ss+4,ss+4)+0.5*(ek(4,4+ss)+ek(4+ss,4));// norma da energia do erro dual-right exato da sol. enriquecida
	//       nkR[5] += ek(ss+5,ss+5)+0.5*(ek(5,5+ss)+ek(5+ss,5));// norma da energia do erro primal-right exato da sol. enriquecida
	//       nkR[8] += ek(ss+4,ss+5)+(ek(5,4+ss)+ek(5+ss,4))/2.;// estimador (p+2 e exato) - Right
	//       nkR[9] += 0.;    // indice de efetividade
}




void TPZBiharmonicEstimator::ContributeInterfaceBCErrorsSimple( TPZMaterialData &data, TPZMaterialData &dataleft,
															   REAL weight,
															   TPZVec<STATE> &nk,
															   TPZBndCond &bc){
	// return;
	OrderSolutionLeft(data,dataleft);
//	OrderSolutionRight(data,dataright);
	TPZFNMatrix<100,STATE> ek(dataleft.phi.Rows(),dataleft.phi.Rows(),0.),ef(dataleft.phi.Rows(),1,0.);
	
	this->ContributeBCInterface(data,dataleft,weight,ek,ef,bc);
	//      std:: cout<< "x=("<< data.x[0]<<","<< data.x[1]<<")";
	//       ek.Print("Teste", std::cout, EFormatted); 
	nk.Resize(10);
	
	nk[0] += ek(0,0);//norma da energia do erro aprox. dual 
	nk[1] += ek(1,1);// norma da energia do erro aprox.primal 
	nk[2] += ek(0,1);// estimador de erro (p e p+2)
	nk[3] += ek(2,2);// norma da energia do erro dual exato 
	nk[4] += ek(3,3);// norma da energia do erro primal exato 
	nk[5] += ek(2,3);// estimador de erro (p e exato)
	nk[6] += ek(4,4);// norma da energia do erro dual exato enriquecido
	nk[7] += ek(5,5);// norma da energia do erro primal exato enriquecido
	nk[8] += ek(4,5);// estimador de erro (p+2 e exato)
	nk[9] += 0.; // indice de efetividade
	// std:: cout<< "  "<< nk[0]<<"  "<< nk[1]<<"  "<< nk[2]<< "  "<< nk[3]<< "  " <<nk[4]<< "  "<< nk[5]<< std::endl;
}



void TPZBiharmonicEstimator::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u, TPZFMatrix<STATE> &dudx,
									TPZFMatrix<REAL> &axes,
									TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,
									TPZVec<REAL> &values) {
	TPZVec<STATE> sol(1), dsol(9,0.);
	
	Solution(u,dudx,axes,1,sol);
	Solution(u,dudx,axes,2,dsol);
	TPZVec<STATE> tmp_1(1,0.);
	Psi(x, tmp_1);
	values[0]  = sol[0]*tmp_1[0];// parte da implentação do funcional, funcional com núcleo
	values[1]  = (sol[0] - u_exact[0])*(sol[0] - u_exact[0]); // error em norma L2 
	values[2] = 0.;                                          // erro em semi norma H1
	for(int id=0; id<2; id++) {
		values[2]  += (dsol[id] - du_exact(id,0))*(dsol[id] - du_exact(id,0));// erro em semi norma H1
	}//(-grad(u).n)v
	values[3] = values[1]+values[2];  // erro em norma H1
}


void TPZBiharmonicEstimator::Psi(TPZVec<REAL> &x, TPZVec<STATE> &pisci) {
	pisci[0]/*=8*exp(10*(2 - 2*x[0] + pow(x[0],2) - 2*x[1] + pow(x[1],2)))*
			 (1 - 46*x[1] + 1129*pow(x[1],2) - 120000*pow(x[0],7)*pow(-1 + x[1],2)*pow(x[1],2) + 20000*pow(x[0],8)*pow(-1 + x[1],2)*pow(x[1],2) - 
			 2526*pow(x[1],3) + 2043*pow(x[1],4) - 800*pow(x[1],5) + 200*pow(x[1],6) + 
			 x[0]*(-46 + 2116*x[1] - 31836*pow(x[1],2) + 76000*pow(x[1],3) - 73880*pow(x[1],4) + 36800*pow(x[1],5) - 9200*pow(x[1],6)) - 
			 800*pow(x[0],5)*(1 - 46*x[1] + 996*pow(x[1],2) - 2260*pow(x[1],3) + 1910*pow(x[1],4) - 800*pow(x[1],5) + 200*pow(x[1],6)) + 
			 200*pow(x[0],6)*(1 - 46*x[1] + 1986*pow(x[1],2) - 4240*pow(x[1],3) + 2900*pow(x[1],4) - 800*pow(x[1],5) + 200*pow(x[1],6)) + 
			 pow(x[0],2)*(1129 - 31836*x[1] + 322916*pow(x[1],2) - 838160*pow(x[1],3) + 1045960*pow(x[1],4) - 796800*pow(x[1],5) + 397200*pow(x[1],6) - 
			 120000*pow(x[1],7) + 20000*pow(x[1],8)) - 2*pow(x[0],3)*
			 (1263 - 38000*x[1] + 419080*pow(x[1],2) - 1066400*pow(x[1],3) + 1264600*pow(x[1],4) - 904000*pow(x[1],5) + 424000*pow(x[1],6) - 
			 120000*pow(x[1],7) + 20000*pow(x[1],8)) + pow(x[0],4)*
			 (2043 - 73880*x[1] + 1045960*pow(x[1],2) - 2529200*pow(x[1],3) + 2604400*pow(x[1],4) - 1528000*pow(x[1],5) + 580000*pow(x[1],6) - 
			 120000*pow(x[1],7) + 20000*pow(x[1],8)));*/
	
	// Experimento 1
	=pow(2.,8)*( 24.*(1.-x[0])*(1.-x[0])*x[0]*x[0] 
				+ 24.*(1.-x[1])*(1.-x[1])*x[1]*x[1]
				+ 2*(2*(1-x[0])*(1-x[0])-8*(1-x[0])*x[0]+2*x[0]*x[0])*(1-x[1])*(1-x[1])
				- 8*(2*(1-x[0])*(1-x[0])-8*(1-x[0])*x[0]+2*x[0]*x[0])*(1-x[1])* x[1] 
				+ 2* (2* (1 - x[0])* (1 - x[0]) - 8* (1 - x[0])* x[0] + 2* x[0]* x[0])* x[1]* x[1]
				+ (2*(1-x[0])*(1-x[0])-8*(1-x[0])*x[0]+2*x[0]*x[0])*(2* (1 - x[1])* (1 - x[1]) 
																	 - 8* (1 - x[1])* x[1] + 2* x[1]* x[1]));//validacao
}// função núcleo do funcional

void TPZBiharmonicEstimator::OrderSolution(TPZMaterialData &data)
{
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	data.phi.Resize(6,1);
	data.dphix.Resize(data.dsol[0].Rows(),6);
	data.phi(0,0) = data.sol[0][3]- data.sol[0][2];//erro aprox. dual
	data.phi(1,0) = data.sol[0][1]- data.sol[0][0];//erro aprox. primal 
	data.phi(2,0) = 0.;//erro exato dual p
	data.phi(3,0) = 0.;//erro exato primal p 
	data.phi(4,0) = 0.;//erro exato dual p enriquecido 
	data.phi(5,0) = 0.;//erro exato primal p enriquecido 
	
	int nd = data.dphix.Rows();
	int id;
	for(id=0; id<nd; id++)
	{
		data.dphix(id,0) = data.dsol[0](id,3)-data.dsol[0](id,2);// derivada do erro aprox dual 
		data.dphix(id,1) = data.dsol[0](id,1)-data.dsol[0](id,0);// derivada do erro aprox. primal
		data.dphix(id,2) = 0. ;// derivada do erro exato dual
		data.dphix(id,3) = 0. ;// derivada do erro exato primal
		data.dphix(id,4) = 0. ;// derivada do erro exato dual enriquecido
		data.dphix(id,5) = 0. ;// derivada do erro exato primal enriquecido
	}
	TPZVec<STATE> u_exactp(1);
	TPZFMatrix<STATE> du_exactp(nd,1);
	TPZVec<STATE> u_exactd(1);
	TPZFMatrix<STATE> du_exactd(nd,1);
	// erro exato primal 
	if (this->fPrimalExactSol )
	{ this->fPrimalExactSol(data.x,u_exactp,du_exactp);
		data.phi(3,0)=u_exactp[0]-data.sol[0][0];
		for(id=0; id<nd; id++)
		{ data.dphix(id,3) = du_exactp(id,0)- data.dsol[0](id,0);// derivada do erro primal exato 
			data.dphix(id,5) = du_exactp(id,0)- data.dsol[0](id,1);// derivada do erro primal exato enriquecido
		}
	}
	
	//erro exato dual
	if (this->fDualExactSol )
	{ this->fDualExactSol(data.x,u_exactd,du_exactd);
		data.phi(2,0)=u_exactd[0]-data.sol[0][2];
		// std::cout<< "dual exato= "<< u_exactd[0]<< std::endl;
		// std::cout<< "dual aprx p2= "<< data.sol[2]<< std::endl;
		// std::cout<< "dual aprx q3= "<< data.sol[3]<< std::endl;
		for(id=0; id<nd; id++)
		{ data.dphix(id,2) = du_exactd(id,0)- data.dsol[0](id,2);// derivada do erro dual exato 
			data.dphix(id,4) = du_exactd(id,0)- data.dsol[0](id,3);// derivada do erro dual exato enriquecido 
		}
	}
}

void TPZBiharmonicEstimator::OrderSolutionLeft(TPZMaterialData &data, TPZMaterialData &dataleft)
{
    int numbersol = dataleft.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	dataleft.phi.Resize(6,1);
	dataleft.dphix.Resize(dataleft.dsol[0].Rows(),6);
	dataleft.phi(0,0) = dataleft.sol[0][3]-dataleft.sol[0][2];//erro aprox. dual
	dataleft.phi(1,0) = dataleft.sol[0][1]-dataleft.sol[0][0];//erro aprox. primal 
	dataleft.phi(2,0) = 0.;//erro exato dual   
	dataleft.phi(3,0) = 0.;//erro exato primal  
	dataleft.phi(4,0) = 0.;//erro exato dual enriquecido
	dataleft.phi(5,0) = 0.;//erro exato primal enriquecido
	
	int nd = dataleft.dphix.Rows();
	int id;
	for(id=0; id<nd; id++)
	{
		dataleft.dphix(id,0) = dataleft.dsol[0](id,3)-dataleft.dsol[0](id,2);// derivada do erro aprox dual 
		dataleft.dphix(id,1) = dataleft.dsol[0](id,1)-dataleft.dsol[0](id,0);// derivada do erro aprox. primal
		dataleft.dphix(id,2) = 0. ;// derivada do erro exato dual
		dataleft.dphix(id,3) = 0. ;// derivada do erro exato primal
		dataleft.dphix(id,4) = 0. ;// derivada do erro exato dual enriquecido
		dataleft.dphix(id,5) = 0. ;// derivada do erro exato primal enriquecido 
	}
	
	TPZVec<STATE> u_exactp(1);
	TPZFMatrix<STATE> du_exactp(nd,1);
	TPZVec<STATE> u_exactd(1);
	TPZFMatrix<STATE> du_exactd(nd,1);
	
	// erro exato primal 
	if (this->fPrimalExactSol )
	{
		this->fPrimalExactSol(data.x,u_exactp,du_exactp);
		dataleft.phi(3,0)=u_exactp[0]-dataleft.sol[0][0];
		dataleft.phi(5,0)=u_exactp[0]-dataleft.sol[0][1];
		for(id=0; id<nd; id++)
		{ dataleft.dphix(id,3) = du_exactp(id,0)- dataleft.dsol[0](id,0);// derivada do erro primal exato
			dataleft.dphix(id,5) = du_exactp(id,0)- dataleft.dsol[0](id,1);// derivada do erro primal exato enriquecido 
		}
	}
	
	//erro exato dual
	if (this->fDualExactSol )
	{
		this->fDualExactSol(data.x,u_exactd,du_exactd);
		dataleft.phi(2,0)=u_exactd[0]-dataleft.sol[0][2];
		dataleft.phi(4,0)=u_exactd[0]-dataleft.sol[0][3];
		for(id=0; id<nd; id++)
		{ dataleft.dphix(id,2) = du_exactd(id,0)- dataleft.dsol[0](id,2);// derivada do erro dual exato 
			dataleft.dphix(id,4) = du_exactd(id,0)- dataleft.dsol[0](id,3);// derivada do erro dual exato enriquecido
		}
	}
}

void TPZBiharmonicEstimator::OrderSolutionRight(TPZMaterialData &data, TPZMaterialData &dataright)
{
    int numbersol = dataright.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	dataright.phi.Resize(6,1);
	dataright.dphix.Resize(dataright.dsol[0].Rows(),6);
	dataright.phi(0,0) = dataright.sol[0][3]-dataright.sol[0][2];//erro aprox. dual
	dataright.phi(1,0) = dataright.sol[0][1]-dataright.sol[0][0];//erro aprox. primal 
	dataright.phi(2,0) = 0.;//erro exato dual   
	dataright.phi(3,0) = 0.;//erro exato primal  
	dataright.phi(4,0) = 0.;//erro exato dual enriquecido   
	dataright.phi(5,0) = 0.;//erro exato primal  enriquecido
	
	int nd = dataright.dphix.Rows();
	int id;
	for(id=0; id<nd; id++)
	{
		dataright.dphix(id,0) = dataright.dsol[0](id,3)-dataright.dsol[0](id,2);// derivada do erro aprox dual 
		dataright.dphix(id,1) = dataright.dsol[0](id,1)-dataright.dsol[0](id,0);// derivada do erro aprox. primal
		dataright.dphix(id,2) = 0. ;// derivada do erro exato dual
		dataright.dphix(id,3) = 0. ;// derivada do erro exato primal
		dataright.dphix(id,4) = 0. ;// derivada do erro exato dual enriquecido
		dataright.dphix(id,5) = 0. ;// derivada do erro exato primal enriquecido
		
	}
	
	TPZVec<STATE> u_exactp(1);
	TPZFMatrix<STATE> du_exactp(8,1);
	TPZVec<STATE> u_exactd(1);
	TPZFMatrix<STATE> du_exactd(8,1);
	
	// erro exato primal 
	if (this->fPrimalExactSol )
	{ this->fPrimalExactSol(data.x,u_exactp,du_exactp);
		dataright.phi(3,0)=u_exactp[0]-dataright.sol[0][0];
		dataright.phi(5,0)=u_exactp[0]-dataright.sol[0][1];
		for(id=0; id<nd; id++)
		{ dataright.dphix(id,3) = du_exactp(id,0)- dataright.dsol[0](id,0);// derivada do erro primal exato 
			dataright.dphix(id,5) = du_exactp(id,0)- dataright.dsol[0](id,1);// derivada do erro primal exato enriquecido
		}
	}
	
	//erro exato dual
	if (this->fDualExactSol )
	{ this->fDualExactSol(data.x,u_exactd,du_exactd);
		dataright.phi(2,0)=u_exactd[0]-dataright.sol[0][2];
		dataright.phi(4,0)=u_exactd[0]-dataright.sol[0][3];
		for(id=0; id<nd; id++)
		{ dataright.dphix(id,2) = du_exactd(id,0)- dataright.dsol[0](id,2);// derivada do erro dual exato 
			dataright.dphix(id,4) = du_exactd(id,0)- dataright.dsol[0](id,3);// derivada do erro dual exato enriquecido
		}
	}
	
}

int TPZBiharmonicEstimator::ClassId() const{
    return Hash("TPZBiharmonicEstimator") ^ TPZBiharmonic::ClassId() << 1;
}
