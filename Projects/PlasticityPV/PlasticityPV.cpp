/**
 * @file
 * @brief Tests for sub structuration
 * @author Philippe Devloo
 * @since 2006
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
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

#include "TPZYCMohrCoulombPV.h"
#include "pzsandlerextPV.h"
#include "TPZPlasticStepPV.h"


#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
#endif

REAL NormVec(TPZManVector<REAL,3> &vec1);

REAL NormVecOfMat(TPZFNMatrix <9> mat);

REAL mypow(REAL a, int n);

void TesteFAD();

REAL InnerVecOfMat(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2);

void TaylorCheck();

void TaylorCheck2();

void TaylorCheck3(); // Tomara que o ultimo!!

void TestePlasticStepPV();

TPZFNMatrix<9> FromEgToMat(TPZManVector<REAL,3> egva, TPZManVector<TPZTensor<REAL>, 3 > &Eigenvec);

TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat);

TPZFMatrix<REAL> ProdT(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2);

int main()
{
	
	TestePlasticStepPV();
	
	
	//TaylorCheck3();
	
	// TESTE DO MOHRCOULOMBPV
	TPZYCMohrCoulombPV *MohrCoulombPV = new TPZYCMohrCoulombPV;
	
	// Gerando tensor elastico
	TPZElasticResponse ER;
	ER.SetUp(1000, 0.2);
	TPZTensor<REAL> epstotalT, sigtrialT, sigtrial;
	epstotalT.XX() = -50 * 0.01;
	epstotalT.YY() = -40 * 0.01;
	ER.Compute(epstotalT, sigtrialT);
	sigtrialT.Print(std::cout);
	TPZTensor<REAL>::TPZDecomposed SigDec;

	sigtrialT.EigenSystem(SigDec);
	SigDec.Print(cout);	
	
	TPZManVector<REAL, 3> sigtrialvec(3,0.), sigprojectvec(3,0.);
	TPZYCMohrCoulombPV::TComputeSequence toto;
	toto.fGamma.resize(1);
	
		
	sigtrialvec = SigDec.fEigenvalues;
	TPZManVector<TPZTensor<REAL>, 3 > &Eigenvec = SigDec.fEigenvectors;
	TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
	
	for (int i = 0; i < 3; i++) 
	{
		EigenvecMat[i] = Eigenvec[i];
		EigenvecMat[i].Print("asd");
		std::cout << Eigenvec[i][0] << std::endl;
		std::cout << Eigenvec[i][1] << std::endl;
		std::cout << Eigenvec[i][2] << std::endl;
		std::cout << Eigenvec[i][3] << std::endl;
		std::cout << Eigenvec[i][4] << std::endl;
		std::cout << Eigenvec[i][5] << std::endl;
	}
	//sigtrialvec[0] = 500.;
	
	TPZManVector<TFad<3,REAL>,3 > sigtrialfad(3), sigprojfad(3);
	for (int i = 0; i < 3; i++) 
	{
		sigtrialfad[i].val() = sigtrialvec[i];
		sigtrialfad[i].fastAccessDx(i) = 1;
	}
	
	
	MohrCoulombPV->ReturnMapPlane<REAL>(sigtrialvec, sigprojectvec, toto, ER);
	
	toto.fGamma[0] = 0.;
	MohrCoulombPV->ReturnMapPlane<TFad<3,REAL> >(sigtrialfad, sigprojfad, toto, ER);
	std::cout << "sigtrialfad:\n" << sigtrialfad << "\nsigprojfad:\n" << sigprojfad << std::endl;
	
	//
	
	//TaylorCheck();
	
	TPZFNMatrix<9> FadProj(3,3,0.);
	for (int i = 0; i < 3; i++) 
	{
		for (int j = 0; j < 3; j++) 
		{
			FadProj(i,j) = sigprojfad[i].dx(j);
		}
	}
	FadProj.Print("FadProj");
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
					REAL temp = 2 * G * EigenvecMat[j](j,kj) * EigenvecMat[j](j,ki);
				
					if (ki == kj) 
					{
						temp += lambda;
					}
					temp *= FadProj(i,j);
					dSigDe(l,k) += temp * Eigenvec[i][l];
				}///l
			}///j
		}///i		
	}///k
	
	
	
	//REAL critplat = MohrCoulombPV->PhiPlane<REAL>(sigprojectvec);
	//std::cout << "CritPlat = " << critplat << std::endl;
	
	//  critplat = MohrCoulombPV->PhiPlane<REAL>(sigmaproject);
	//	std::cout << "CritPlat = " << critplat << std::endl;
	
	
	// Testando Left Edge
	/*
	toto.fGamma.Resize(2,0.);
	toto.fGamma.Fill(0.);
	std::cout << "toto.fGamma = " << toto.fGamma << std::endl;
	sigtrialvec[0] = 500.;
	sigtrialvec[1] = 0.;
	sigtrialvec[2] = 0.;	
	sigprojectvec.Fill(0.);
	std::cout << "sigprojectvec = " << sigprojectvec << std::endl;
	MohrCoulombPV->ReturnMapLeftEdge(sigtrialvec, sigprojectvec, toto);
	MohrCoulombPV->ReturnMapRightEdge(sigtrialvec, sigprojectvec, toto);
	*/
	//	template<class T>
	//	bool ReturnMapLeftEdge(const typename TPZTensor<T>::TPZDecomposed &sigma_trial, typename TPZTensor<T>::TPZDecomposed &sigma_projected,
	//												 TComputeSequence &memory)
	
	// Será que está na superfície?
	//    T PhiPlane(typename TPZTensor<T>::TPZDecomposed &sigma) const
	
	
	
	//	template<class T>
	//	bool ReturnMapPlane(const typename TPZTensor<T>::TPZDecomposed &sigma_trial, typename TPZTensor<T>::TPZDecomposed &sigma_projected, 
	//											TComputeSequence &memory)
	
	
	
	// FIM TESTE DO MOHRCOULOMBPV
	
}

void TestePlasticStepPV()
{
	TPZSandlerExtended materialmodel(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1);
	TPZTensor<REAL> epsT,Sigma;
	TPZElasticResponse ER;
	ER.SetUp(100., 0.25);
	REAL Alpha = 0.;
	TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> stepPV(Alpha);
	stepPV.fYC = materialmodel;
	stepPV.fER = ER;
	
	REAL iniEpsAxi = -0.0001;
	std::map<REAL,REAL> loadcicle;
	for (int i = 0; i < 1000; i++) {
		epsT.XX() = iniEpsAxi*i;
		stepPV.ApplyStrainComputeSigma(epsT,Sigma);
		loadcicle[-epsT.XX()] = -Sigma.XX();
	}
	
	ofstream out("LoadCicle.nb");
	std::map<REAL,REAL>::iterator it = loadcicle.begin();
	out << "loadcicle=" << "{";
	out << "{" << it->first << "," << it->second << "}";
	it++;
	for (; it != loadcicle.end(); it++) {
		out << ",{" << it->first << "," << it->second << "}";
	}
	out << "};" << endl;
	out << "ListPlot[loadcicle,Joined->True,PlotStyle->Thick]" << endl;
}

void TaylorCheck3() // Tomara que o ultimo!!
{
	TPZYCMohrCoulombPV *MCPV = new TPZYCMohrCoulombPV;
	TPZElasticResponse ER;
	ER.SetUp(1000, 0.25);
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
	eps.EigenSystem(DecompEps);
	TPZFNMatrix <9> epsMat(3,3,0.),depsMat(3,3,0.), sigtrMat(3,3,0.);
	epsMat = eps;
	depsMat = deps;
	TPZVec<REAL> eigenvaluestemp(3);
	TPZFMatrix <REAL> EigenvectorsTrue(3,3,0.);
	long numiter = 1000;
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
	DecompEps.fEigenvectors[0].Print(cout);
	DecompEps.fEigenvectors[1].Print(cout);
	DecompEps.fEigenvectors[2].Print(cout);
	
	//sigtrtemp.Print(std::cout);
	TPZTensor<REAL>::TPZDecomposed Decomp;
	sigtrtemp.EigenSystem(Decomp);
	TPZManVector<TPZTensor<REAL>, 3 > &Eigenvec = Decomp.fEigenvectors;
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
	MCPV->ReturnMapRightEdge<TFad <3,REAL> >(sigtr,sigpr,toto2,ER);
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
			tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
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
		sigtrtemp.EigenSystem(Decomp1);
		for (int i = 0 ; i < 3; i++) {
			sigtr1[i] = Decomp1.fEigenvalues[i];
		}
		
		ER.Compute(eps2, sigtrtemp);
		sigtrtemp.EigenSystem(Decomp2);
		for (int i = 0 ; i < 3; i++) {
			sigtr2[i] = Decomp2.fEigenvalues[i];
		}		
		
		toto.fGamma[0] = 0.;
		MCPV->ReturnMapPlane<REAL>(sigtr1,sigpr1,toto,ER);
		toto.fGamma[0] = 0.;
		MCPV->ReturnMapPlane<REAL>(sigtr2,sigpr2,toto,ER);
		
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
		sigprMat = FromMatToVoight(FromEgToMat(sigprmanvec, Decomp.fEigenvectors));
//		sigprMat.Print("sigprojmat");
//		cout << "sigpr = " << sigprmanvec << endl;				
		
		sigpr1Mat = FromMatToVoight(FromEgToMat(sigpr1, Decomp1.fEigenvectors));
//		sigpr1Mat.Print("sigprojmat1");
//		cout << "sigpr1 = " << sigpr1 << endl;				
		
		sigpr2Mat = FromMatToVoight(FromEgToMat(sigpr2, Decomp2.fEigenvectors));
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
	TPZYCMohrCoulombPV *MCPV = new TPZYCMohrCoulombPV;
	TPZElasticResponse ER;
	ER.SetUp(1000, 0.25);
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
	sigtrtemp.EigenSystem(Decomp);
	TPZManVector<TPZTensor<REAL>, 3 > &Eigenvec = Decomp.fEigenvectors;
	TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
	for (int i = 0; i < 3; i++) 
	{
		EigenvecMat[i] = Eigenvec[i];
	}
	
	for (int i = 0 ; i < 3; i++) {
		sigtr[i].fastAccessDx(i) = 1;
		sigtr[i].val() = Decomp.fEigenvalues[i];
	}
	
	TPZYCMohrCoulombPV::TComputeSequence toto;
	toto.fGamma.resize(1);
	toto.fGamma[0] = 0.;
	
	MCPV->ReturnMapPlane<TFad <3,REAL> >(sigtr,sigpr,toto,ER);
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
					dSigDe(l,k) += temp * Eigenvec[i][l];
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
		sigtrtemp.EigenSystem(Decomp1);
		for (int i = 0 ; i < 3; i++) {
			sigtr1[i] = Decomp1.fEigenvalues[i];
		}
		
		ER.Compute(eps2, sigtrtemp);
		sigtrtemp.EigenSystem(Decomp2);
		for (int i = 0 ; i < 3; i++) {
			sigtr2[i] = Decomp2.fEigenvalues[i];
		}		
		
		toto.fGamma[0] = 0.;
		MCPV->ReturnMapPlane<REAL>(sigtr1,sigpr1,toto,ER);
		toto.fGamma[0] = 0.;
		MCPV->ReturnMapPlane<REAL>(sigtr2,sigpr2,toto,ER);
		
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
		sigprMat = FromMatToVoight(FromEgToMat(sigprmanvec, Decomp.fEigenvectors));
		//sigprMat.Print("sigprojmat");
		//cout << "sigpr = " << sigprmanvec << endl;				
		
		sigpr1Mat = FromMatToVoight(FromEgToMat(sigpr1, Decomp1.fEigenvectors));
		//sigpr1Mat.Print("sigprojmat1");
		//cout << "sigpr1 = " << sigpr1 << endl;				

		sigpr2Mat = FromMatToVoight(FromEgToMat(sigpr2, Decomp2.fEigenvectors));
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
	TPZYCMohrCoulombPV *MohrCoulombPV = new TPZYCMohrCoulombPV;
	
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

	TPZElasticResponse ER;
	ER.SetUp(1000, 0.2);
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
		
		MohrCoulombPV->ReturnMapPlane<TFad<3,REAL> >(sigtrialfad, sigprojfad, toto, ER);
		MohrCoulombPV->ReturnMapPlane<TFad<3,REAL> >(sigtrialfad1, sigprojfad1, toto, ER);
		MohrCoulombPV->ReturnMapPlane<TFad<3,REAL> >(sigtrialfad2, sigprojfad2, toto, ER);//
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

REAL NormVecOfMat(TPZFNMatrix <9> mat)
{
	REAL norm = 0.;
	for (int i = 0; i < mat.Rows(); i++) {
		norm += mat(i,0) * mat(i,0);
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

TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat)
{
	TPZFNMatrix <6> voi(6,1,0.);
	int k = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i ; j < 3; j++) {
			voi(k++,0) = mat(i,j);
		}
	}
	return voi;	
}

REAL InnerVecOfMat(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
{
	REAL dot = 0.;
	for (int i = 0; i < m1.Rows(); i++) {
		dot += m1(i,0) * m2(i,0);
	}
	return dot;
}

TPZFMatrix<REAL> ProdT(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
{
	TPZFMatrix<REAL> mat(3,3,0.);
	for (int i = 0; i < 3; i++) {
		for (int j = 0 ; j < 3; j++) {
			mat(i,j) = m1(i,0) * m2(j,0);
		}
	}
	return mat;
}