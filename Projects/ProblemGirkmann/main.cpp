#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iostream>
#include <string>
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "tools.h"
#include <fstream>
#include <sstream>
using namespace std;

//ofstream outfile("Resultados.txt");
int main(int argc, char *argv[])
{
	ofstream outfile("Resultados.txt");
	outfile.precision(16);
	
	gRefDBase.InitializeAllUniformRefPatterns();
	gRefDBase.InitializeRefPatterns();
	
	std::string logs("log4cxx.girkmann");
	InitializePZLOG(logs);
	
	
	REAL rc = 15.;
	REAL h = 0.06;
	REAL pi = atan(1.)*4.;
	REAL alpha = 2.*pi/9.;
	REAL a = 0.6;
	REAL b = 0.5;
	REAL pesoEspecifico = 32.69*1000.; //N/m3
	REAL E = 20.59*10.e9;//N/m^2
	REAL poisson = 0.;
	bool anelleve = true;
	tools * Shell = new tools(rc, h, alpha, a, b, pesoEspecifico,E,poisson,anelleve);
	
//	int ndiv = 30;
//	//int nh = 1;
//	int nDiv_DirectRefR = 3;
//	int nDiv_DirectRefL = 4;
//	int nDiv_DirectRefp = 0;
//	REAL sim = 1.;
//	REAL pen = 0.;//100.*E;
//	int p = 3;	
//	for(int nh = 0; nh<= 4; nh++) 
//	{
//		REAL Mreal =-37.45 /*-41.12*/;
//		REAL Qreal =942.5 /*943.8*/;
//		
//		//------------- RESOLUCAO MALHA CONTINUA ---------------------------------------------------------------------------------------
//		//TPZGeoMesh * gmesh = Shell->MalhaGeoGen(ndiv, nDiv_DirectRefR, nDiv_DirectRefL,nDiv_DirectRefp, false, 4);
//		//		//ofstream arg("gmeshContin.txt");
//		//		//gmesh->Print(arg);
//		//		
//		//		TPZCompMesh * cmesh = Shell->MalhaCompGen(gmesh, p);
//		//		TPZAnalysis an(cmesh);
//		//		//ofstream arg4("cmeshContin.txt");
//		//		//cmesh->Print(arg4);
//		//		
//		//		Shell->SolveSist(an, cmesh);
//		//		
//		//		TPZVec<REAL> QeM(9,0.);
//		//		QeM = Shell->CalcCortMomento(cmesh);
//		//		delete cmesh;
//		//		cmesh = 0;
//		//		delete gmesh;
//		//		gmesh = 0;
//		//		
//		//		
//		outfile <<"=====Calculo com, nh = " << nh << ", Ndiv = "<< ndiv <<" e p = " << p <<endl;/*"", RefDir = "<< nDiv_DirectRefL<<" e RefDirp = "<< nDiv_DirectRefp<<" ====" <<endl;*/
//		//		outfile<<"RESULTADOS SEM INTERFACE"<<endl;
//		//		outfile <<" MomentoR = "<<QeM[0]<<"  Erro Relativo ==> " << fabs((Mreal - QeM[0])/Mreal)*100.<<" %"<<endl;
//		//		outfile <<" MomentoL = "<<QeM[2]<<"  Erro Relativo ==> " << fabs((Mreal - QeM[2])/Mreal)*100.<<" %"<<endl;
//		//		outfile <<" MomentoMedia = "<<QeM[4]<<"  Erro Relativo ==> " << fabs((Mreal - QeM[4])/Mreal)*100.<<" %"<<endl;
//		//	      outfile <<endl;
//		//		outfile <<" CortanteR = "<<QeM[1]<<"  Erro Relativo ==> " << fabs((Qreal - QeM[1])/Qreal)*100.<<" %"<<endl;
//		//		outfile <<" CortanteL = "<<QeM[3]<<"  Erro Relativo ==> " << fabs((Qreal - QeM[3])/Qreal)*100.<<" %"<<endl;
//		//		outfile <<" CortanteMedia= "<<QeM[5]<<"  Erro Relativo ==> " << fabs((Qreal - QeM[5])/Qreal)*100.<<" %"<<endl;
//		//		outfile <<endl;
//		//		outfile <<" T1zR = "<<QeM[6]<<endl;
//		//		outfile <<" T1zL = "<<QeM[7]<<endl;
//		//		outfile <<" T1z = "<<QeM[8]<<endl;
//		//		outfile <<endl;
//		
//		//cmesh->Print();
//		
//		//------------- RESOLUCAO MALHA DESCONTINUA ---------------------------------------------------------------------------------------
//		TPZGeoMesh * gmesh2 = Shell->MalhaGeoGen(ndiv, nDiv_DirectRefR, nDiv_DirectRefL, nDiv_DirectRefp, true, 4);
//		//ofstream arg3("gmeshdiscont.txt");
//		//gmesh2->Print(arg3);
//		
//		Shell->RefinamentoUniforme(*gmesh2, nh);
//		
//		ofstream arg3("gmeshdiscont.txt");
//		gmesh2->Print(arg3);
//		
//		TPZCompMesh * cmesh_changed = Shell->MalhaCompMeshWithInterface(gmesh2, p,sim, pen);
//		ofstream arg2("cmeshdiscont.txt");
//		cmesh_changed->Print(arg2);
//		
//		int nEq = cmesh_changed->NEquations();
//		cout << "\n Numero de Equacoes = " << nEq<<endl;
//		
//		TPZAnalysis an_ch(cmesh_changed);
//		Shell->SolveSist(an_ch, cmesh_changed);
//		
//		
//		//cmesh_changed->Print();
//		
//		//Shell->PrintInterface(cmesh_changed);
//		
//		TPZVec<REAL> QeM2;
//		QeM2 = Shell->CalcCortMomento(cmesh_changed);
//		
//		delete cmesh_changed;
//		cmesh_changed = 0;
//		delete gmesh2;
//		gmesh2 = 0;
//		
//		outfile<<"RESULTADOS COM INTERFACE"<<endl;
//		outfile <<" MomentoR = "<<QeM2[0]<<"  Erro Relativo ==> " << fabs((Mreal - QeM2[0])/Mreal)*100.<<" %"<<endl;
//		outfile <<" MomentoL = "<<QeM2[2]<<"  Erro Relativo ==> " << fabs((Mreal - QeM2[2])/Mreal)*100.<<" %"<<endl;
//		outfile <<" MomentoMedia = "<<QeM2[4]<<"  Erro Relativo ==> " << fabs((Mreal - QeM2[4])/Mreal)*100.<<" %"<<endl;
//		outfile <<endl;
//		outfile <<" CortanteR = "<<QeM2[1]<<"  Erro Relativo ==> " << fabs((Qreal - QeM2[1])/Qreal)*100.<<" %"<<endl;
//		outfile <<" CortanteL = "<<QeM2[3]<<"  Erro Relativo ==> " << fabs((Qreal - QeM2[3])/Qreal)*100.<<" %"<<endl;
//		outfile <<" CortanteMedia= "<<QeM2[5]<<"  Erro Relativo ==> " << fabs((Qreal - QeM2[5])/Qreal)*100.<<" %"<<endl;
//		outfile <<endl;
//		outfile <<" T1zR = "<<QeM2[6]<<endl;
//		outfile <<" T1zL = "<<QeM2[7]<<endl;
//		outfile <<" T1z = "<<QeM2[8]<<endl;
//		outfile <<endl;
//	}
	
	
	//------------Validacao da formulacao fraca--------------------------------------
	//int num = 1;
	//	REAL E = 10;
	//	REAL nu = 5;
	//	REAL fx = 1.;
	//	REAL fy = 1.;
	//	REAL simet =-1.;
	//	REAL penal = 1.;
	//	TPZElasticityAxiMaterial ElastAxiMat(num, E, nu, fx, fy, simet, penal);
	//	
	//	TPZMaterialData MatData;	
	//	MatData.phil.Resize(1,1), MatData.phir.Resize(1,1);
	//	MatData.dphixl.Resize(2,1);
	//	MatData.dphixr.Resize(2,1);
	//	MatData.axesleft.Redim(3,3), MatData.axesright.Redim(3,3);
	//	MatData.normal.Resize(3);
	//	
	//	MatData.phil(0,0)=2.;  MatData.phir(0,0)=3.;
	//	MatData.dphixl(0,0)=-1.; MatData.dphixl(1,0)=2.;
	//	MatData.dphixr(0,0)=-1.; MatData.dphixr(1,0)=2.;
	//	MatData.axesleft(0,0)=1.; MatData.axesleft(1,1)=1.;MatData.axesleft(2,2)=1.;
	//	MatData.axesright(0,0)=1.; MatData.axesright(1,1)=1.;MatData.axesright(2,2)=1.;
	//	MatData.normal[0]=1.;MatData.normal[1]=1.;MatData.normal[2]=0.;
	//	
	//	
	//	vector<REAL> Orig(3);		Orig[0]  = 0.;	Orig[1]  = 0.;	Orig[2]  = 0.;
	//	vector<REAL> AxisZ(3);		AxisZ[0] = 0.;	AxisZ[1] = 1.;	AxisZ[2] = 0.;
	//	vector<REAL> AxisR(3);		AxisR[0] = 1.;	AxisR[1] = 0.;	AxisR[2] = 0.;
	//	ElastAxiMat.SetOrigin(Orig, AxisZ, AxisR);
	//	
	//	REAL peso = 1.;
	//	TPZFMatrix KIJ, FJ;
	//	KIJ.Resize(4,4);
	//	
	//	ElastAxiMat.ContributeInterface(MatData, peso,KIJ, FJ);
	//	KIJ.Print();	
	
    return EXIT_SUCCESS;
}
