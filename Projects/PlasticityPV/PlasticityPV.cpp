/**
 * @file
 * @brief Tests for sub structuration
 * @author Philippe Devloo
 * @since 2006
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzmat1dlin.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZElasticResponse.h"
#include "tpzautopointer.h"
#include "TPZTimer.h" 

#include "TPZYCMohrCoulombPV.h"
#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"


#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
#endif

REAL NormVec(TPZManVector<REAL,3> &vec1);

REAL mypow(REAL a, int n);

void TesteFAD();

void TaylorCheck();

void TaylorCheck2();

void TaylorCheck3(); // Tomara que o ultimo!!

void CurvaFig12PlasticPV();

void CurvaFig12Diogo();

void DepPlasticPV();

void DepPlasticPVMC();

TPZFNMatrix<9> FromEgToMat(TPZManVector<REAL,3> egva, TPZManVector<TPZTensor<REAL>, 3 > &Eigenvec);

TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat);

int main()
{
	InitializePZLOG();
  CurvaFig12PlasticPV();
//	const REAL Phi = M_PI/9., Psi = M_PI/9., c = 9.35;
//	TPZElasticResponse ER;
//	ER.SetUp(1000, 0.25);
//	
//	TPZYCMohrCoulombPV *MohrCoulombPV = new TPZYCMohrCoulombPV(Phi,Psi,c,ER);
//	
//	TPZManVector< TFad<3,REAL>, 3> sigtrialfad(3),sigprojfad(3);
//	for (int i = 0 ; i < 3 ; i++){
//		sigtrialfad[i].fastAccessDx(i) = 1;
//	}
//	sigtrialfad[0].val() = 151;
//	sigtrialfad[1].val() = 160.;
//	sigtrialfad[2].val() = 149.;
//	
//	TPZYCMohrCoulombPV::TComputeSequence toto;
//	toto.fGamma.resize(2);
//	toto.fGamma[0] = 0.;
//	toto.fGamma[1] = 0.;
//	
//	REAL epsbarfake = 0.;
//	MohrCoulombPV->ReturnMapApex<TFad<3,REAL> >(sigtrialfad, sigprojfad, toto, epsbarfake);
//	TPZFNMatrix<9> tang;
//	MohrCoulombPV->ComputeApexTangent(tang, epsbarfake);
//	tang.Print("tang:");
//	std::cout << "sigprojfad = " << sigprojfad << std::endl;
	
	TPZTensor<REAL> teste;
	teste.XX() = 2.;
	teste.YY() = 3.;
	teste.ZZ() = 2.;
	teste.YZ() = 0.;
	
	TPZTensor<REAL>::TPZDecomposed decomp;
	teste.EigenSystemJacobi(decomp);
	
	
	
	DepPlasticPVMC();
//	TaylorCheck();
//	TaylorCheck2();
//	TaylorCheck3();	
	return 0;
}

void DepPlasticPV()
{

	const REAL A = 0.25, B = 0.67, C = 0.18, D = 0.67, K = 66.67, G = 40., W = 0.066, R = 2.5, Phi = 0., N = 0., Psi = 1., kappa_0 = 0.0;
	TPZSandlerExtended materialmodel(A, B, C, D, K, G, W, R, Phi, N, Psi,kappa_0);
	TPZTensor<REAL> epsT,Sigma;
	TPZElasticResponse ER;
	ER.SetUp(100., 0.25);
	STATE kproj,kprev,epv=0.;
	kproj=0.;
	kprev=0.13;
	materialmodel.Firstk(epv,kprev);
	TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> stepPV(kprev);
	stepPV.fYC = materialmodel;
	stepPV.fER = ER;
	
	TPZTensor<REAL> eps, deps, sigma;
	eps.XX() = -0.001;
	eps.YY() = -0.002;
	eps.ZZ() = -0.003;
//	eps.XY() = -23. * 0.00001;
//	eps.XZ() = -24. * 0.001;
//	eps.YZ() = -65. * 0.001;
	deps.XX() = -0.00001;
	deps.YY() = -0.0000;
	deps.ZZ() = -0.0000;
//	deps.XY() = 0.0001;
//    deps.XZ() = 0.63 * 0.00002;
//	deps.YZ() = 0.64686 * 0.00001;
	TPZFNMatrix<36> Dep(6,6,0.);
	TPZManVector<REAL,5> conv;
	stepPV.TaylorCheck(eps, deps, kprev,conv);
	//stepPV.ApplyStrainComputeDep(eps, sigma, Dep); 	
	Dep.Print("Dep:");
}

void DepPlasticPVMC()
{

	const REAL Phi = M_PI/9., Psi = M_PI/9., c = 9.35;
	TPZElasticResponse ER;
	ER.SetUp(100., 0.25);
	TPZYCMohrCoulombPV materialmodel(Phi,Psi,c,ER);
	TPZTensor<REAL> epsT,Sigma;

	TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> stepPV(0.);
	stepPV.fYC = materialmodel;
	stepPV.fER = ER;
	
	TPZTensor<REAL> eps, deps;
	REAL kprev = 0.;
	eps.XX() = -0.01;
	eps.YY() = -0.01;
	eps.ZZ() = -0.8;
	//eps.XY() = 0.01;
	//eps.XZ() = -24. * 0.001;
	//eps.YZ() = -65. * 0.001;
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\n-------------------------- BEGINING OF TAYLORCHECK TEST USING PRINCIPAL VALUES --------------------------" << endl;
		
		TPZManVector<REAL,5> conv;
		deps.XX() = eps.XX();
		deps *= 1/20;
		stepPV.TaylorCheck(eps, deps, kprev, conv);
		str << "\nXX:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;
		
		deps.XX() = 0.;
		deps.YY() = eps.YY();
		deps *= 1/20.;
		stepPV.TaylorCheck(eps, deps, kprev, conv);
		str << "\nYY:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;

		deps.YY() = 0.;
		deps.ZZ() = eps.ZZ();
		deps *= 1/20.;
		stepPV.TaylorCheck(eps, deps, kprev, conv);
		str << "\nZZ:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;

		deps.ZZ() = 0.;
		deps.XY() = -0.001;
		stepPV.TaylorCheck(eps, deps, kprev, conv);
		str << "\nXY:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;

		deps.XY() = 0.;
		deps.XZ() = -0.001;
		stepPV.TaylorCheck(eps, deps, kprev, conv);
		str << "\nXZ:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;

		deps.XZ() = 0.;
		deps.YZ() = -0.001;
		stepPV.TaylorCheck(eps, deps, kprev, conv);
		str << "\nYZ:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;
		
		str << "\n-------------------------- END OF TAYLORCHECK TEST USING PRINCIPAL VALUES --------------------------" << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	

//	TPZTensor<REAL> eps, sigma;
//	eps.XX() = 1.e-7;
//	eps.YY() = 1.e-7;
//	eps.ZZ() = 1.e-7;
//	//	eps.XY() = -23. * 0.001;
//	//	eps.XZ() = -24. * 0.001;
//	//	eps.YZ() = -65. * 0.001;
//	TPZFNMatrix<36> Dep(6,6,0.);
//	stepPV.ApplyStrainComputeDep(eps, sigma, Dep);
//	Dep.Print("Dep:");

}

void CurvaFig12PlasticPV()
{
	const REAL A = 0.25, B = 0.67, C = 0.18, D = 0.67, K = 66.67, G = 40., W = 0.066, R = 2.5, Phi = 0., N = 0., Psi = 1., kappa_0 = 0.13;
	TPZSandlerExtended materialmodel(A, B, C, D, K, G, W, R, Phi, N, Psi, kappa_0);
	TPZTensor<REAL> epsT,Sigma;
	TPZElasticResponse ER;
	ER.SetUp(100., 0.25);
	STATE kproj,kprev,epv=0.;
	kproj=0.;
	kprev=0.13;
	materialmodel.Firstk(epv,kprev);
	
  REAL cohesion = A - C;
  REAL PhiMC = B*C;
  
  TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> stepPV(kprev);
	stepPV.fYC = materialmodel;
  stepPV.fER = ER;
  
  TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> MCPV;
  MCPV.fYC.SetUp(PhiMC, PhiMC, cohesion, ER);
  MCPV.fER = ER;
  
	REAL deltaEps= -0.0000135;
	std::list<std::pair<REAL,REAL> > loadcicle;
 	std::list<std::pair<REAL,REAL> > unloadcicle;
	for (int i = 0; i < 500; i++) {
		if (i == 250) {
			deltaEps *= -1.;
		}
		//epsT.YY() = iniEpsAxi*i;
		//epsT.ZZ() = iniEpsAxi*i;
		epsT.Print(std::cout);
		MCPV.ApplyStrainComputeSigma(epsT,Sigma);
    if (deltaEps < 0){
  		loadcicle.push_back(std::pair<REAL,REAL>(-epsT.XX(),-Sigma.XX()));
    }
    else{
  		unloadcicle.push_back(std::pair<REAL,REAL>(-epsT.XX(),-Sigma.XX()));
    }

		//loadcicle[-epsT.XX()] = -Sigma.XX();
		epsT.XX() += deltaEps;
	}
	
	ofstream out("LoadCicleMC.nb");
	std::list<std::pair<REAL,REAL> >::iterator it = loadcicle.begin();
	out << "loadcicle=" << "{";
	out << "{" << it->first << "," << it->second << "}";
	it++;
	for (; it != loadcicle.end(); it++) {
		out << ",{" << it->first << "," << it->second << "}";
	}
	out << "};" << endl;
  
  out << "unloadcicle=" << "{";
  it = unloadcicle.begin();
	out << "{" << it->first << "," << it->second << "}";
	it++;
	for (; it != unloadcicle.end(); it++) {
		out << ",{" << it->first << "," << it->second << "}";
	}
	out << "};" << endl;

	out << "loadcicleplot=ListPlot[loadcicle,Joined->True,PlotStyle->Thick];" << endl;
 	out << "unloadcicleplot=ListPlot[unloadcicle,Joined->True,PlotStyle->{Red,Thick}];" << endl;
  out << "Show[loadcicleplot,unloadcicleplot]" << endl;
}

void CurvaFig12Diogo()
{
	TPZManVector<STATE,3> epsPnext(3),epsT(3),deleps(3),epssol(3),deltaepsP(3),sigproj(3),sigtrial(3),deltasigma(3);
	TPZSandlerExtended materialmodel(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1,0.0);
	
	ofstream outfile("FIGURA_12x.nb");
	outfile << "VecFig12={";
	
	deleps[0]= -0.00135;
	
	for (int k=0; k<3;k++) {
		deltaepsP[k]=0.;
		epsT[0]=0.;
		epsPnext[0]=0.;
		epssol[0]=0.;
	}
	
	STATE kproj,kprev,epv=0.;
	TPZFMatrix<STATE> GradSigma;
	
	kproj=0.;
	kprev=0.13;
	materialmodel.Firstk(epv,kprev);
	for(int i=0;i<65;i++)
	{
		for (int k=0; k<3;k++) {
			epssol[k]=epsT[k]-epsPnext[k];
		}
		
		materialmodel.ApplyStrainComputeElasticStress(epssol, sigtrial);
		materialmodel.ProjectSigmaDep(sigtrial,kprev,sigproj,kproj,GradSigma);
		
		outfile <<"{" << -epsT[0]<< "," << -sigproj[0] << "}";
		if (i!=64) {
			outfile << ",";
		}
		
		if(i==50)
		{
			deleps[0]*=-1;
		}
		
		for (int k=0; k<3;k++) {
			deltasigma[k]=sigtrial[k]-sigproj[k];
		}
		
		materialmodel.ApplyStressComputeElasticStrain(deltasigma, deltaepsP);
		
		for (int k=0; k<3;k++) {
			epsPnext[k]+=deltaepsP[k];
			epsT[k]+=deleps[k];
		}
		kprev=kproj;
		
	}
	outfile << "};\nListPlot[VecFig12,Joined->True]";
}

void TaylorCheck3() // Tomara que o ultimo!!
{
	const REAL Phi = M_PI/9., Psi = M_PI/9., c = 9.35;
	TPZElasticResponse ER;
	ER.SetUp(1000, 0.25);
	TPZYCMohrCoulombPV *MCPV = new TPZYCMohrCoulombPV(Phi,Psi,c,ER);
	
	TPZTensor <REAL> eps, eps1, eps2, deps, sigtrtemp;
	eps.XX() = -50. * 0.01;
	eps.YY() = -40. * 0.01;
	eps.ZZ() = 45 * 0.01;
	eps.XY() = -23. * 0.001;
	eps.XZ() = -24. * 0.001;
	eps.YZ() = -25. * 0.001;
	deps.XX() = 3.86274 * 0.01;
	deps.YY() = 1.13727 * 0.01;
	deps.ZZ() = 1. * 0.01;
	deps.XY() = 0.62686 * 0.01;
	deps.XZ() = 0.63686 * 0.0001;
	deps.YZ() = 0.64686 * 0.0001;
	TPZTensor<REAL>::TPZDecomposed DecompEps;
	eps.EigenSystemJacobi(DecompEps);
	TPZFNMatrix <9> epsMat(3,3,0.),depsMat(3,3,0.), sigtrMat(3,3,0.);
	epsMat = eps;
	depsMat = deps;
	TPZVec<REAL> eigenvaluestemp(3);
	TPZFMatrix <REAL> EigenvectorsTrue(3,3,0.);
	int64_t numiter = 1000;
	REAL tol = 1.e-6;

	TPZManVector <TFad<3,REAL>,3> sigtr(3), sigpr(3);
	TPZManVector <REAL,3> sigtr1(3), sigtr2(3), sigpr1(3), sigpr2(3);
	
	ER.Compute(eps, sigtrtemp);
	sigtrMat = sigtrtemp;
	sigtrMat.SolveEigensystemJacobi(numiter, tol, eigenvaluestemp, EigenvectorsTrue);
	sigtrMat.Print("sigtrmat");
	std::cout << "eigenvalues = " << eigenvaluestemp << endl;
	std::cout << "eigenvalues modo 2 = " << DecompEps.fEigenvalues << endl;
	EigenvectorsTrue.Print("asd");
	DecompEps.fEigentensors[0].Print(cout);
	DecompEps.fEigentensors[1].Print(cout);
	DecompEps.fEigentensors[2].Print(cout);
	
	//sigtrtemp.Print(std::cout);
	TPZTensor<REAL>::TPZDecomposed Decomp;
	sigtrtemp.EigenSystemJacobi(Decomp);
	TPZManVector<TPZTensor<REAL>, 3 > &Eigenvec = Decomp.fEigentensors;
	TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
	TPZManVector<TPZFMatrix<REAL>,3> epsegveFromProj(3);
	for (int i = 0; i < 3; i++) 
	{
		EigenvecMat[i] = Eigenvec[i];
		epsegveFromProj[i].Resize(3,1);
		for (int j = 0; j < 3; j++) {
			epsegveFromProj[i](j,0) = EigenvecMat[i](j,0);
		}
	}
	for (int i = 0; i < 3; i++) {
		REAL normvec = 0.;
		normvec = NormVecOfMat(epsegveFromProj[i]);
		for (int j = 0; j < 3; j++) {
			epsegveFromProj[i](j,0) /= normvec;
		}
	}
	epsegveFromProj[0].Print("vec1");
	epsegveFromProj[1].Print("vec2");
	epsegveFromProj[2].Print("vec3");
	
	for (int i = 0 ; i < 3; i++) {
		sigtr[i].fastAccessDx(i) = 1;
		sigtr[i].val() = Decomp.fEigenvalues[i];
	}
	
	TPZYCMohrCoulombPV::TComputeSequence toto;
	toto.fGamma.resize(1);
	toto.fGamma[0] = 0.;
	
	//Teste Right Edge
	
	TPZYCMohrCoulombPV::TComputeSequence toto2;
	toto2.fGamma.resize(2);
	toto2.fGamma[0] = 0.;
	toto2.fGamma[1] = 0.;
	cout << "sigtr = " << sigtr << endl;
	REAL epsbarfake = 0.;
	MCPV->ReturnMapPlane<TFad <3,REAL> >(sigtr,sigpr,toto2, epsbarfake);
	cout << "sigpr = " << sigpr << endl;

	//MCPV->ReturnMapPlane<TFad <3,REAL> >(sigtr,sigpr,toto,ER);
	cout << "sigpr" << sigpr << endl;
//	std::cout << "sigtr = " << sigtr << std::endl;
//	std::cout << "sigpr = " << sigpr << std::endl;
	TPZManVector<REAL,3> sigprmanvec(3,0.);	
	for (int i = 0; i < 3; i++) sigprmanvec[i] = sigpr[i].val();
	
	TPZFNMatrix<9> FADMat(3,3,0.);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			FADMat(i,j) = sigpr[i].dx(j);
		}
	}
	
	// Criando matriz tangente
	TPZFNMatrix<36> dSigDe(6,6,0.);
	
	//Montando a matriz tangente
	int kival[] = {0,0,0,1,1,2};
	int kjval[] = {0,1,2,1,2,2};
	REAL G = ER.G();
	REAL lambda = ER.Lambda();
	
	for (int k = 0; k < 6; k++) 
	{
		int ki, kj;
		ki = kival[k];
		kj = kjval[k];
		for (int i = 0; i < 3; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				for (int l = 0; l < 6; l++) 
				{
					REAL temp = 2 * G * EigenvecMat[j](kj,ki); // * EigenvecMat[j](j,ki);
					
					if (ki == kj) 
					{
						temp += lambda;
					}
					else {
						temp *= 2.;
					}

					temp *= FADMat(i,j);
					dSigDe(l,k) += temp * Eigenvec[i][l];
				}///l
			}///j
		}///i		
	}///k
//	dSigDe.Print("dSigDe");
	
	REAL deigensig = 0., deigeneps = 0.;
	
	TPZManVector<TPZFMatrix<REAL>,3> epsegve(3,0.);
	for (int i = 0; i < 3; i++) {
		epsegve[i].Redim(3,1);
		for (int j = 0; j < 3; j++) {
			epsegve[i](j,0) = EigenvectorsTrue(i,j);  
		}		
	}
/*	EigenvectorsTrue.Print("aqq");
	epsegve[0].Print("asd");
	epsegve[1].Print("asd");
	epsegve[2].Print("asd");*/
	
	TPZFNMatrix<6,REAL> factorMat(3,3,0.);
	for (int i = 0; i < 2; i++) {
		for (int j = i+1; j<3 ; j++) {
			deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
			deigensig = sigprmanvec[i] - sigprmanvec[j];
				
			TPZFNMatrix<9,REAL> tempMat(3,1,0.);
			depsMat.Multiply(epsegveFromProj[i], tempMat);
//			tempMat.Print("tempMat");
//			epsegve[i].Print("epsegve i");
//			epsegve[j].Print("epsegve j");
			REAL deij = InnerVecOfMat(tempMat,epsegveFromProj[j]);
			
			REAL factor = deigensig * deij / deigeneps;
			tempMat.Redim(3, 3);
//			tempMat.Print("tempMat");
                        TPZManVector<REAL, 3> epsegveFromProji(3), epsegveFromProjj(3);
                        for (unsigned int k = 0; k < 3; ++k) {
                            epsegveFromProji[k] = epsegveFromProj[i](k,0);
                            epsegveFromProjj[k] = epsegveFromProj[j](k,0);
                        }

			tempMat = ProdT(epsegveFromProji,epsegveFromProjj) + ProdT(epsegveFromProjj,epsegveFromProji);
			factorMat += tempMat * factor;
//			tempMat.Print("temMat");
//			factorMat.Print("factorMat");
		}
	}
	TPZFNMatrix <6> factorVoight = FromMatToVoight(factorMat);
//	factorVoight.Print("qwe");
	
	REAL scale = 0.01;
	REAL alphatable[] = {0.1,0.2,0.3,0.4,0.5,0.6};
	for (int i = 0; i < 6; i++) {
		alphatable[i] *= scale;
	}
	for (int ia = 0 ; ia < 5; ia++) {
		REAL alpha1 = alphatable[ia];
		REAL alpha2 = alphatable[ia+1];
		eps1.Scale(0.);
		eps2.Scale(0.);
		eps1 = eps;
		eps2 = eps;
		eps1.Add(deps, alpha1);
		eps2.Add(deps, alpha2);
		
		TPZTensor<REAL>::TPZDecomposed Decomp1, Decomp2;
		ER.Compute(eps1, sigtrtemp);
		sigtrtemp.EigenSystemJacobi(Decomp1);
		for (int i = 0 ; i < 3; i++) {
			sigtr1[i] = Decomp1.fEigenvalues[i];
		}
		
		ER.Compute(eps2, sigtrtemp);
		sigtrtemp.EigenSystemJacobi(Decomp2);
		for (int i = 0 ; i < 3; i++) {
			sigtr2[i] = Decomp2.fEigenvalues[i];
		}		
		
		toto.fGamma[0] = 0.;
		MCPV->ReturnMapPlane<REAL>(sigtr1,sigpr1,toto,epsbarfake);
		toto.fGamma[0] = 0.;
		MCPV->ReturnMapPlane<REAL>(sigtr2,sigpr2,toto,epsbarfake);
		
		TPZFNMatrix <6> deps1(6,1,0.),deps2(6,1,0.);
		TPZFNMatrix <9> depsMat(3,3,0.);
		depsMat = deps;
		deps1 = FromMatToVoight(depsMat);
		deps2 = FromMatToVoight(depsMat);
//		deps1.Print("deps1");
//		deps2.Print("deps2");

		TPZFNMatrix <6> tanmult1(6,1,0.), tanmult2(6,1,0.);
		dSigDe.Multiply(deps1, tanmult1);
		dSigDe.Multiply(deps2, tanmult2);
		//tanmult1.Print("tanmult1");
		//tanmult2.Print("tanmult2");
		
		tanmult1 += factorVoight;
		tanmult2 += factorVoight;
		
		for (int i = 0 ; i < 6; i++) {
			tanmult1(i,0) *= alpha1;
			tanmult2(i,0) *= alpha2;
		}
		
//		tanmult1.Print("tamult1");
//		tanmult2.Print("tanmul2");
			
		TPZFNMatrix <6> sigprMat(6,1,0.),sigpr1Mat(6,1,0.),sigpr2Mat(6,1,0.);
		sigprMat = FromMatToVoight(FromEgToMat(sigprmanvec, Decomp.fEigentensors));
//		sigprMat.Print("sigprojmat");
//		cout << "sigpr = " << sigprmanvec << endl;				
		
		sigpr1Mat = FromMatToVoight(FromEgToMat(sigpr1, Decomp1.fEigentensors));
//		sigpr1Mat.Print("sigprojmat1");
//		cout << "sigpr1 = " << sigpr1 << endl;				
		
		sigpr2Mat = FromMatToVoight(FromEgToMat(sigpr2, Decomp2.fEigentensors));
//		sigpr2Mat.Print("sigprojmat2");
//		cout << "sigpr2 = " << sigpr2 << endl;	
		

		
		TPZFNMatrix<6> error1(6,1,0.), error2(6,1,0.);
		for (int i = 0 ; i < 6; i++) {
			error1(i,0) = sigpr1Mat(i,0) - sigprMat(i,0) - tanmult1(i,0);
			error2(i,0) = sigpr2Mat(i,0) - sigprMat(i,0) - tanmult2(i,0);
		}
//		error1.Print("error1");
//		error2.Print("error2");
		REAL norm1, norm2;
		norm1 = NormVecOfMat(error1);
		norm2 = NormVecOfMat(error2);
		//norm1 = fabs(error1(0,0));
		//norm2 = fabs(error2(0,0));
		REAL n;
		n = ( log(norm1) - log(norm2) ) / ( log(alpha1) - log(alpha2) );
		cout << "n = " << n << endl;
		
	}
	
}

void TaylorCheck2()
{
	const REAL Phi = M_PI/9., Psi = M_PI/9., c = 9.35;
	TPZElasticResponse ER;
	ER.SetUp(1000, 0.25);
	TPZYCMohrCoulombPV *MCPV = new TPZYCMohrCoulombPV(Phi,Psi,c,ER);
	
	TPZTensor <REAL> eps, eps1, eps2, deps, sigtrtemp;
	eps.XX() = -50. * 0.01;
	eps.YY() = -40. * 0.01;
	eps.ZZ() = 45 * 0.01;
	deps.XX() = -10. * 0.0001;
	deps.YY() = -11. * 0.0001;
	deps.ZZ() = -12. * 0.0001;
	
	TPZManVector <TFad<3,REAL>,3> sigtr(3), sigpr(3);
	TPZManVector <REAL,3> sigtr1(3), sigtr2(3), sigpr1(3), sigpr2(3);
	
	ER.Compute(eps, sigtrtemp);
	sigtrtemp.Print(std::cout);
	TPZTensor<REAL>::TPZDecomposed Decomp;
	sigtrtemp.EigenSystemJacobi(Decomp);
	TPZManVector<TPZTensor<REAL>, 3 > &Eigentensors = Decomp.fEigentensors;
	TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
	for (int i = 0; i < 3; i++) 
	{
		EigenvecMat[i] = Eigentensors[i];
	}
	
	for (int i = 0 ; i < 3; i++) {
		sigtr[i].fastAccessDx(i) = 1;
		sigtr[i].val() = Decomp.fEigenvalues[i];
	}
	
	TPZYCMohrCoulombPV::TComputeSequence toto;
	toto.fGamma.resize(1);
	toto.fGamma[0] = 0.;
	
	REAL epsbarfake = 0.;
	MCPV->ReturnMapPlane<TFad <3,REAL> >(sigtr,sigpr,toto, epsbarfake);
	std::cout << "sigtr = " << sigtr << std::endl;
	std::cout << "sigpr = " << sigpr << std::endl;
	TPZManVector<REAL,3> sigprmanvec(3,0.);	
	for (int i = 0; i < 3; i++) sigprmanvec[i] = sigpr[i].val();
	
	TPZFNMatrix<9> FADMat(3,3,0.);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			FADMat(i,j) = sigpr[i].dx(j);
		}
	}
	
	// Criando matriz tangente
	TPZFNMatrix<36> dSigDe(6,6,0.);
	
	//Montando a matriz tangente
	int kival[] = {0,0,0,1,1,2};
	int kjval[] = {0,1,2,1,2,2};
	REAL G = ER.G();
	REAL lambda = ER.Lambda();
	
	for (int k = 0; k < 6; k++) 
	{
		int ki, kj;
		ki = kival[k];
		kj = kjval[k];
		for (int i = 0; i < 3; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				for (int l = 0; l < 6; l++) 
				{
					REAL temp = 2 * G * EigenvecMat[j](ki,kj);// * EigenvecMat[j](j,ki);
					
					if (ki == kj) 
					{
						temp += lambda;
					}
					temp *= FADMat(i,j);
					dSigDe(l,k) += temp * Eigentensors[i][l];
				}///l
			}///j
		}///i		
	}///k
	
	//Realizando a correcao pelos giros rigidos
	
	
	
	REAL alphatable[] = {0.1,0.2,0.3,0.4,0.5,0.6};
	for (int ia = 0 ; ia < 5; ia++) {
		REAL alpha1 = alphatable[ia];
		REAL alpha2 = alphatable[ia+1];
		eps1.Scale(0.);
		eps2.Scale(0.);
		eps1 = eps;
		eps2 = eps;
		eps1.Add(deps, alpha1);
		eps2.Add(deps, alpha2);

		TPZTensor<REAL>::TPZDecomposed Decomp1, Decomp2;
		ER.Compute(eps1, sigtrtemp);
		sigtrtemp.EigenSystemJacobi(Decomp1);
		for (int i = 0 ; i < 3; i++) {
			sigtr1[i] = Decomp1.fEigenvalues[i];
		}
		
		ER.Compute(eps2, sigtrtemp);
		sigtrtemp.EigenSystemJacobi(Decomp2);
		for (int i = 0 ; i < 3; i++) {
			sigtr2[i] = Decomp2.fEigenvalues[i];
		}		
		
		toto.fGamma[0] = 0.;
		MCPV->ReturnMapPlane<REAL>(sigtr1,sigpr1,toto,epsbarfake);
		toto.fGamma[0] = 0.;
		MCPV->ReturnMapPlane<REAL>(sigtr2,sigpr2,toto,epsbarfake);
		
		TPZFNMatrix <6> deps1(6,1,0.),deps2(6,1,0.);
		TPZFNMatrix <9> depsMat(3,3,0.);
		depsMat = deps;
		deps1 = FromMatToVoight(depsMat);
		deps2 = FromMatToVoight(depsMat);
		for (int i = 0 ; i < 6; i++) {
			deps1(i,0) *= alpha1;
			deps2(i,0) *= alpha2;
		}
		TPZFNMatrix <6> tanmult1(6,1,0.), tanmult2(6,1,0.);
		dSigDe.Multiply(deps1, tanmult1);
		dSigDe.Multiply(deps2, tanmult2);
		//tanmult1.Print("tanmult1");
		//tanmult2.Print("tanmult2");
		
		TPZFNMatrix <6> sigprMat(6,1,0.),sigpr1Mat(6,1,0.),sigpr2Mat(6,1,0.);
		sigprMat = FromMatToVoight(FromEgToMat(sigprmanvec, Decomp.fEigentensors));
		//sigprMat.Print("sigprojmat");
		//cout << "sigpr = " << sigprmanvec << endl;				
		
		sigpr1Mat = FromMatToVoight(FromEgToMat(sigpr1, Decomp1.fEigentensors));
		//sigpr1Mat.Print("sigprojmat1");
		//cout << "sigpr1 = " << sigpr1 << endl;				

		sigpr2Mat = FromMatToVoight(FromEgToMat(sigpr2, Decomp2.fEigentensors));
		//sigpr2Mat.Print("sigprojmat2");
		//cout << "sigpr2 = " << sigpr2 << endl;	
		
		TPZFNMatrix<6> error1(6,1,0.), error2(6,1,0.);
		for (int i = 0 ; i < 6; i++) {
			error1(i,0) = sigpr1Mat(i,0) - sigprMat(i,0) - tanmult1(i,0);
			error2(i,0) = sigpr2Mat(i,0) - sigprMat(i,0) - tanmult2(i,0);
		}
		REAL norm1, norm2;
		norm1 = NormVecOfMat(error1);
		norm2 = NormVecOfMat(error2);
		REAL n;
		n = log(norm1/norm2 ) / log (alpha1 / alpha2);
		cout << "n = " << n << endl;
		
	}
}

void TaylorCheck()
{
	const REAL Phi = M_PI/9., Psi = M_PI/9., c = 9.35;
	TPZElasticResponse ER;
	ER.SetUp(1000, 0.25);
	TPZYCMohrCoulombPV *MohrCoulombPV = new TPZYCMohrCoulombPV(Phi,Psi,c,ER);
	
	TPZManVector< TFad<3,REAL>, 3> sigtrialfad(3), sigtrialfad1(3),sigtrialfad2(3), sigprojfad(3), sigprojfad1(3), sigprojfad2(3);
	for (int i = 0 ; i < 3 ; i++){
		sigtrialfad[i].fastAccessDx(i) = 1;
		sigtrialfad1[i].fastAccessDx(i) = 1;
		sigtrialfad2[i].fastAccessDx(i) = 1;
	}
	sigtrialfad[0].val() = -138.888889;
	sigtrialfad[1].val() = -138.888889;
	sigtrialfad[2].val() = -555.555556;
	
	TPZManVector <REAL, 3> dsigtrial(3,0.);
	dsigtrial[0] = -10.;
	dsigtrial[1] = -11.;
	dsigtrial[2] = -12.;
	REAL alphatable[] = {0.1,0.2,0.3,0.4,0.5,0.6};

	TPZYCMohrCoulombPV::TComputeSequence toto;
	toto.fGamma.resize(1);
	toto.fGamma[0] = 0.;

	for (int ia = 0; ia < 5; ia++) {
		// Ponto 1;
		REAL alpha1 = alphatable[0];
		REAL alpha2 = alphatable[ia+1];
		for (int i = 0; i < 3; i++) {
			sigtrialfad1[i].val() = sigtrialfad[i].val() + alpha1*dsigtrial[i];
			sigtrialfad2[i].val() = sigtrialfad[i].val() + alpha2*dsigtrial[i];
		}
		
		REAL epsbarfake = 0.;
		MohrCoulombPV->ReturnMapPlane<TFad<3,REAL> >(sigtrialfad, sigprojfad, toto,epsbarfake);
		MohrCoulombPV->ReturnMapPlane<TFad<3,REAL> >(sigtrialfad1, sigprojfad1, toto,epsbarfake);
		MohrCoulombPV->ReturnMapPlane<TFad<3,REAL> >(sigtrialfad2, sigprojfad2, toto,epsbarfake);//
//		std::cout << "sigtrialfad:\n" << sigtrialfad << "\nsigprojfad:\n" << sigprojfad << std::endl;
//		std::cout << "sigtrialfad1:\n" << sigtrialfad1 << "\nsigprojfad1:\n" << sigprojfad1 << std::endl;
//		std::cout << "sigtrialfad2:\n" << sigtrialfad2 << "\nsigprojfad2:\n" << sigprojfad2 << std::endl;
		
		TPZFNMatrix<9> FadProj(3,3,0.);
		for (int i = 0 ; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				FadProj(i,j) = sigprojfad[i].dx(j);
			}
		}	
		
		TPZFNMatrix<3> dsigMat(3,1,0.), tanmult1(3,1,0.),tanmult2(3,1,0.);
		for (int i = 0; i < 3; i++) {
			dsigMat(i,0) = dsigtrial[i] * alpha1;
		}
		
		FadProj.Multiply(dsigMat, tanmult1);
		
		for (int i = 0; i < 3; i++) {
			dsigMat(i,0) = dsigtrial[i] * alpha2;
		}	
		FadProj.Multiply(dsigMat, tanmult2);
		
		//tanmult1.Print("TanMult1:");
		
		TPZManVector<REAL,3> errorvec1(3,0.),errorvec2(3,0.);
		for (int i = 0; i < 3; i++) {
			errorvec1[i] = sigprojfad1[i].val() - sigprojfad[i].val() - tanmult1(i,0);
			errorvec2[i] = sigprojfad2[i].val() - sigprojfad[i].val() - tanmult2(i,0);
		}
		REAL norm1, norm2;
		norm1 = NormVec(errorvec1);
		norm2 = NormVec(errorvec2);
		
//		REAL n = ( log(norm1) - log(norm2) ) / ( log(alpha1) - log(alpha2) ) ;
		REAL n = ( log(norm1/norm2) ) / ( log(alpha1/alpha2) ) ;
		std::cout << "TaylorCheck parameter = " << n << std::endl;
	}
	
}

REAL NormVec(TPZManVector<REAL,3> &vec)
{
	REAL norm = 0.;
	for (int i = 0; i < vec.NElements(); i++) {
		norm += vec[i] * vec[i];
	}
	norm = sqrt(norm);
	return norm;
}

REAL mypow(REAL a, int n)
{
	if (n == 0) return 1.;
	return (a * mypow(a,n-1));
}

void TesteFAD()
{
	// Aprendendo a usar FAD
	TFad<3,REAL> g, a, b,c;
	
	std::cout << a << std::endl;
	
	std::cout << b << std::endl;
	
	std::cout << c << std::endl;
	
	std::cout << g << std::endl;
	
	a.fastAccessDx(0) = 1;
	b.fastAccessDx(1) = 1;
	c.fastAccessDx(2) = 1;
	
	std::cout << a << std::endl;
	
	std::cout << b << std::endl;
	
	std::cout << c << std::endl;
	
	std::cout << g << std::endl;
	
	a.val() = 4.;
	b.val() = 2.;
	c.val() = 3.;
	
	g = a*b*b*sin(c);
	
	std::cout << a << std::endl;
	
	std::cout << b << std::endl;
	
	std::cout << c << std::endl;
	
	std::cout << g << std::endl;
	
}

TPZFNMatrix<9> FromEgToMat(TPZManVector<REAL,3> egva, TPZManVector<TPZTensor<REAL>, 3 > &Eigenvec)
{
	TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
	for (int i = 0; i < 3; i++) 
	{
		EigenvecMat[i] = Eigenvec[i];
	}
	
	TPZFNMatrix <9> mat(3,3,0.);
	for (int k = 0; k < 3; k++) {
		for (int i = 0 ; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				mat(i,j) += egva[k] * EigenvecMat[k](i,j);
			}
		}		
	}
	return mat;	
}



