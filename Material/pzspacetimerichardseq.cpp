/**
 * \file
 * @brief Contains implementations of the TPZSpaceTimeRichardsEq methods.
 */
//$Id: pzspacetimerichardseq.cpp,v 1.6 2011-05-11 02:24:19 phil Exp $

#include "pzspacetimerichardseq.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzaxestools.h"

#include <cmath>
using namespace std;

/** @brief Inicializing local variable TCoeff */
double TCoeff = 1./60.;
/** @brief Inicializing local variable LCoeff */
double LCoeff = 1000.;
/** @brief Inicializing loval variable deltaDerivada */
double deltaDerivada = 1.e-3;

TPZSpaceTimeRichardsEq::TPZSpaceTimeRichardsEq(): TPZMaterial()
{
}


TPZSpaceTimeRichardsEq::TPZSpaceTimeRichardsEq(int id): TPZMaterial(id)
{
}


TPZSpaceTimeRichardsEq::TPZSpaceTimeRichardsEq(int matid, REAL Alpha, REAL N, REAL ThetaS, REAL ThetaR, REAL Ks) 
: TPZMaterial(matid){
	this->Set(Alpha, N, ThetaS, ThetaR, Ks);
}


TPZSpaceTimeRichardsEq::~TPZSpaceTimeRichardsEq()
{
}

void TPZSpaceTimeRichardsEq::Set(REAL Alpha, REAL N, REAL ThetaS, REAL ThetaR, REAL Ks){
	this->fAlpha = Alpha;
	this->fN = N;
	this->fThetaS = ThetaS;
	this->fThetaR = ThetaR;
	this->fKs = Ks;
}

int TPZSpaceTimeRichardsEq::Dimension(){
	return 2;
}

int TPZSpaceTimeRichardsEq::NStateVariables(){
	return 1;
}

void TPZSpaceTimeRichardsEq::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
	TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &dphi = data.dphix;
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	const REAL sol = data.sol[0][0];
	
	const REAL BetaBarT = 0*LCoeff*data.detjac/2.; //beta=(0,1)
	
	TPZFNMatrix<2> dsol(2,1,0.);
	TPZAxesTools<REAL>::Axes2XYZ(data.dsol[0], dsol, data.axes);
	
	const int phr = phi.Rows();
	int i, j;
	
	for(i = 0; i < phr; i++){
		const REAL BetaBarGradV = BetaBarT*dphi(1,i);
		ef(i,0) += -1.*weight *( -1.*sol*dphi(1,i) +1.*dsol(1,0)*BetaBarGradV + /*(K/C)*/(dsol(0,0))*dphi(0,i));
		for(j = 0; j < phr; j++){
			ek(i,j) += weight * ( -1.*phi(j,0)*dphi(1,i)+dphi(1,j)*BetaBarGradV + dphi(0,i)*dphi(0,j) );
		}
	}//for i
	
	// ek.Identity();ek*=LCoeff;
	
}//Contribute

void TPZSpaceTimeRichardsEq::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
	
	const REAL v2 = bc.Val2()(0,0);
	TPZFMatrix<REAL> &phi = data.phi;
	const int phr = phi.Rows();
	int in, jn;
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	
	switch (bc.Type()){
			
			// Dirichlet condition
		case 0 : {
			for(in = 0 ; in < phr; in++) {
				ef(in,0) += weight * ( gBigNumber * phi(in,0) * (v2 - data.sol[0][0]) );
				for (jn = 0 ; jn < phr; jn++) {
					ek(in,jn) +=  gBigNumber * phi(in,0) * phi(jn,0) * weight;
				}
			}
			break;
		}
			
			// Neumann condition
		case 1:{
			// please implement me
		}
			break;
			
			// outflow condition
		case 3 : { 
			
			const REAL sol = data.sol[0][0];
			//       const REAL C = this->C_Coef(sol);
			//       REAL ConvDir[2] = {0., C};
			REAL ConvDir[2] = {0., 1.}; 
			REAL normal[2];
			normal[0] = data.axes(0,1);
			normal[1] = -1.*data.axes(0,0);
			
			REAL ConvNormal = ConvDir[0]*normal[0] + ConvDir[1]*normal[1];
			if(ConvNormal > 0.) {
				for(int il = 0; il < phr; il++) {
					for(int jl = 0; jl < phr; jl++) {
						ek(il,jl) += weight * ConvNormal * phi(il)*phi(jl);
					}
					ef(il,0) += -1. * weight * ConvNormal * phi(il) * sol;
				}
			}
			else{
				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
			}  
		}
			break;
			
		default:{
			std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
		}
	}//switch
	
}//ContributeBC

REAL TPZSpaceTimeRichardsEq::C_Coef(REAL sol){
	
	sol = sol/LCoeff;
	
	//filter
	//   if(sol > -0.3){
	//     sol = -0.3;
	//   }
	//   if(sol < -12.){
	//     sol = -12.;
	//   }
	
	REAL n = this->fN;
	REAL m = 1.-(1./n);
	REAL TR = this->fThetaR;
	REAL TS = this->fThetaS;
	REAL a = this->fAlpha;
	REAL result = (m*n*pow(pow(a,(REAL)2.)*pow(sol,(REAL)2.),n/2.)*pow(1./(1. + pow(pow(a,(REAL)2.)*pow(sol,(REAL)2.),n/2.)),1. + m)*(TR - TS))/sol;
	
	//filter
	//   if(sol > -0.3 || sol < -12.) return sol*sol;
	
	return result/LCoeff;
}


REAL TPZSpaceTimeRichardsEq::Se(REAL sol){
	REAL n = this->fN;
	REAL m = 1.-(1./n);
	REAL result = pow((REAL)(1./(1.+pow(fabs(this->fAlpha*sol),n))),m);
	return result;
}

REAL TPZSpaceTimeRichardsEq::Theta(REAL Se){
	REAL result = Se*(fThetaS-fThetaR)+fThetaR;
	return result;
}

/** @brief Return sign of the real A. If A is closed to zero up to tolerance tol returns zero (in absolute value) */
int Sign(double A, double tol){
	if(fabs(A) < tol) return 0;
	if(A > 0.) return +1;
	return -1;
}

REAL TPZSpaceTimeRichardsEq::DCDsol(REAL sol){
	REAL sol1 = sol*(1.-deltaDerivada);
	REAL antes =   this->C_Coef(sol1);
	REAL sol2 = sol * (1.+deltaDerivada);
	REAL depois = this->C_Coef(sol2);
	REAL result = (depois-antes)/(sol2-sol1);
	return result;
}

void TPZSpaceTimeRichardsEq::AnalysisOfParameters(REAL sol0, REAL solL, char* filename){
	int np = 100;
	TPZFMatrix<REAL> C(np,2), K(np,2), dCdSol(np,2), dKdsol(np,2);
	double delta = (solL - sol0)/(np-1);
	double sol;
	for(int i = 0; i < np; i++){
		sol = sol0+i*delta;
		C(i,0) = sol;
		C(i,1) = this->C_Coef(sol);
		K(i,0) = sol;
		K(i,1) = this->K_Coef(sol);
		dCdSol(i,0) = sol;
		dCdSol(i,1) = this->DCDsol(sol);
		dKdsol(i,0) = sol;
		dKdsol(i,1) = this->DKDsol(sol);
	}
	
	std::ofstream file(filename);
	C.Print("C=",file,EMathematicaInput);
	K.Print("K=",file,EMathematicaInput);
	dCdSol.Print("dCdSol=",file,EMathematicaInput);
	dKdsol.Print("dKdsol=",file,EMathematicaInput);
	
	K.Print("K=",file);
	
}


REAL TPZSpaceTimeRichardsEq::K_Coef(REAL sol){
	
	//   return -sol * 1e-12;
	sol = sol / LCoeff;
	
	//filter
	//   if(sol > -0.3){
	//     sol = -0.3;
	//   }
	//   if(sol < -12.){
	//     sol = -12.;
	//   }
	
	REAL n = this->fN;
	REAL m = 1.-(1./n);
	REAL Se = this->Se(sol);
	REAL Ks = this->fKs;
	REAL result = Ks*sqrt(Se)*pow(((REAL)1.)-pow(((REAL)1.)-pow(Se,((REAL)1.)/m),m),(REAL)2.);
	
	return result*LCoeff/TCoeff;
}


REAL TPZSpaceTimeRichardsEq::DKDsol(REAL sol){
	
    REAL sol1 = sol*(1.-deltaDerivada);
    REAL antes = this->K_Coef(sol1);
    REAL sol2 = sol * (1.+deltaDerivada);
    REAL depois = this->K_Coef(sol2);
    REAL resposta = (depois-antes)/(sol2-sol1);
    return resposta;
	
	sol = sol / LCoeff;
	REAL n = this->fN;
	REAL m = 1.-(1./n);  
	REAL Se = this->Se(sol);
	REAL DkDSe = 0.5 * this->fKs * (1.-pow(((REAL)1.)-pow(Se,((REAL)1.)/m),m))*(4.*pow(Se,((REAL)-0.5)+(((REAL)1.)/m))*pow(((REAL)1.)-pow(Se,((REAL)1.)/m),m-((REAL)1.))-(-1.+pow(((REAL)1.)-pow(Se,((REAL)1.)/m),m))/sqrt(Se));
	REAL DSeDsol = -(1./sol)*m*n*pow(fabs(this->fAlpha*sol),n)*pow(1./(1.+pow(fabs(this->fAlpha*sol),n)),m+1.);
	REAL dkdsol = DkDSe * DSeDsol;
	
	/*  if(fabs(dkdsol - resposta) > 1e-3){
	 std::cout << "peguei erro 1 no DKDsol\n";
	 }
	 if(Sign(dkdsol,1e-8) != Sign(resposta,1e-8)){
	 std::cout << "peguei erro 2 no DKDsol\n";
	 }*/
	
	//   if(fabs(dkdsol) < 1.){
	//     dkdsol = 1.;
	//   }
	
	return dkdsol*LCoeff/(TCoeff*LCoeff);
}

/*
 #include <iostream>
 #include "pzcmesh.h"
 #include "pzgmesh.h"
 #include "tpzautopointer.h"
 #include "pzanalysis.h"
 #include "pzfstrmatrix.h"
 #include "TPZParFrontStructMatrix.h"
 #include "TPZFrontNonSym.h"
 #include "pzskylstrmatrix.h"
 #include "pzbstrmatrix.h"
 #include "pzstepsolver.h"
 #include "pzcompel.h"
 #include "pzinterpolationspace.h"
 #include "pznonlinanalysis.h"
 
 using namespace std;
 
 const int ndiv = 6-2;
 const double LTotal = 1.*LCoeff;
 const double TTotal = 86400.*TCoeff;
 
 void TPZSpaceTimeRichardsEq::DirichletT0(TPZVec<REAL> &x, TPZVec<REAL> &f) {
 const double DeltaX = LTotal/pow(2.,ndiv);
 f[0] = -10.*LCoeff;
 if(x[0] > (LTotal-DeltaX)){
 const double xr = x[0] - (LTotal-DeltaX);
 double result = -10.*LCoeff + ((-0.75*LCoeff +10.*LCoeff)/DeltaX)*xr;
 f[0] = result;
 }
 }
 
 int TPZSpaceTimeRichardsEq::main(){
 TPZCompMesh * cmesh = TPZSpaceTimeRichardsEq::CreateMesh(LTotal, TTotal, 1, ndiv);
 
 TPZNonLinearAnalysis an(cmesh,cout);
 TPZBandStructMatrix matrix(cmesh);
 an.SetStructuralMatrix(matrix);
 TPZStepSolver step;
 step.SetDirect(ELDLt);
 an.SetSolver(step);
 const int neq = an.Mesh()->NEquations();
 std::cout << "Numero de equacoes = " << neq << std::endl;
 cout << "Banda = " << cmesh->BandWidth() << "\n";
 cout.flush();
 
 TPZMaterial::gBigNumber = 1.e12;
 for(int i = 0; i < an.Solution().Rows(); i++) an.Solution()(i,0) = -10.*(1.+1e-6*i)*LCoeff;
 an.LoadSolution(an.Solution());
 an.IterativeProcess(cout, 1e-10, 50, true, true);
 
 //  int nnodes = cmesh->Reference()->NodeVec().NElements();
 //  for(int i = 0; i < nnodes; i++){
 //    cmesh->Reference()->NodeVec()[i].SetCoord(0, cmesh->Reference()->NodeVec()[i].Coord(0) * 86400.);
 //  }
 
 TPZVec<std::string> scalnames(1);
 scalnames[0] = "state";
 TPZVec<std::string> vecnames(0);
 std::string plotfile = "richards.dx";
 an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
 an.PostProcess(1);  
 
 TPZGeoMesh *gmesh = cmesh->Reference();
 delete cmesh;
 delete gmesh;
 
 return EXIT_SUCCESS;  
 
 }//main
 
 TPZCompMesh * TPZSpaceTimeRichardsEq::CreateMesh(REAL L, REAL Time, int p, int ndiv){
 
 const REAL Alpha = 3.35;
 const REAL N = 2.;
 const REAL ThetaS = 0.368;
 const REAL ThetaR = 0.102;
 const REAL Ks = 9.22e-5;
 
 const int nnodes = 4;
 REAL co[nnodes][2] = {{0.,0.},{L,0.},{L,Time},{0.,Time}};
 const int nelem = 1; 
 int indices[nelem][4] = {{0,1,2,3}};
 const int nelbc = 4;
 int contorno[nelbc][3] = {{0,1,-1},{1,2,-2},{2,3,-3},{3,0,-4}};
 
 TPZGeoMesh *gmesh = new TPZGeoMesh();
 for(int nod=0; nod<nnodes; nod++) {
 int nodind = gmesh->NodeVec().AllocateNewElement();
 TPZManVector<REAL,2> coord(2);
 coord[0] = co[nod][0];
 coord[1] = co[nod][1];
 gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
 }
 
 const int matid = 1;
 for(int el=0; el<nelem; el++) {
 TPZManVector<int,4> nodind(4);
 for(int nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
 int index;
 gmesh->CreateGeoElement(EQuadrilateral,nodind,matid,index);
 }
 
 for(int el=0; el<nelbc; el++) {
 TPZManVector<int,2> nodind(2);
 for(int nod=0; nod<2; nod++) nodind[nod]=contorno[el][nod];
 int bcmatid = contorno[el][2];
 int index;
 gmesh->CreateGeoElement(EOned,nodind,bcmatid,index);
 }
 
 gmesh->BuildConnectivity();
 
 for(int idiv = 0; idiv < ndiv; idiv++){
 TPZVec<TPZGeoEl*> filhos;
 int n = gmesh->NElements();
 for(int i = 0; i < n; i++){
 TPZGeoEl * gel = gmesh->ElementVec()[i];
 if (!gel) continue;
 if(gel->HasSubElement()) continue;
 gel->Divide(filhos);
 }
 }
 
 //  TPZCompEl::SetgOrder(p); 
 TPZCompMesh * cmesh = new TPZCompMesh(gmesh);  
 cmesh->SetDefaultOrder(p);
 cmesh->SetDimModel(2);
 
 TPZAutoPointer<TPZMaterial> mat = new TPZSpaceTimeRichardsEq(matid,Alpha,N,ThetaS,ThetaR,Ks);
 
 int nstate = 1;
 TPZFMatrix<REAL> val1(nstate,nstate,0.),val2(nstate,1,0.);
 
 val2(0,0) = -10.*LCoeff;
 TPZAutoPointer<TPZMaterial> bcT0 ( mat->CreateBC(mat,-1, 0,val1,val2) );
 bcT0->SetForcingFunction(DirichletT0);
 val2(0,0) = -0.75*LCoeff;
 TPZAutoPointer<TPZMaterial> bcXL = mat->CreateBC(mat,-2, 0,val1,val2);
 //no value needed for outflow bc
 TPZAutoPointer<TPZMaterial> bcOutFlow = mat->CreateBC(mat,-3, 3,val1,val2);
 val2(0,0) = -10.*LCoeff;
 TPZAutoPointer<TPZMaterial> bcX0 = mat->CreateBC(mat,-4, 0,val1,val2);
 
 cmesh->InsertMaterialObject(mat);
 cmesh->InsertMaterialObject(bcT0);
 cmesh->InsertMaterialObject(bcXL);
 cmesh->InsertMaterialObject(bcOutFlow);
 cmesh->InsertMaterialObject(bcX0);
 
 //  TPZCompEl::SetgOrder(p);
 cmesh->SetDefaultOrder(p);
 cmesh->AutoBuild();
 
 return cmesh;
 
 }
 */