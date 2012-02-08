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

#include "TPZVTKGeoMesh.h"

using namespace std;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.girkmannproblem"));
#endif

ofstream outfile("Resultados.txt");
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	gRefDBase.InitializeAllUniformRefPatterns();
	gRefDBase.InitializeRefPatterns();
	//std::ifstream inFileRefPat("refpatternsAgnaldo.txt");
	//gRefDBase.ReadRefPatternDBase(inFileRefPat);
	
	outfile.precision(16);
	
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
	
	int ndiv = 240;
	int nDiv_DirectRefCasca = 1;
	int nDiv_DirectRefAnel = 4;
	int nDiv_DirectRefPonto = 0;
	REAL sim = 1.;/*-1.;*/
	REAL pen = 10.;/*10.e10;*/
		
	for(int nh = 1; nh<= 1; nh++) 
	{
		for (int p = 2; p < 9; p++)
		{
			REAL Mreal =-37.45 /*-41.12*/;
			REAL Qreal =942.5 /*943.8*/;

			//------------- RESOLUCAO MALHA CONTINUA ---------------------------------------------------------------------------------------
			/*TPZGeoMesh * gmesh = Shell->MalhaGeoGen(ndiv, nDiv_DirectRefCasca, nDiv_DirectRefAnel,nDiv_DirectRefPonto, false, 4);
			ofstream arg("gmeshContin.txt");
			gmesh->Print(arg);
		
			ofstream file("malhageometrica.vtk");
			TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
		
				
			TPZCompMesh * cmesh = Shell->MalhaCompGen(gmesh, p);
			TPZAnalysis an(cmesh);
			ofstream arg4("cmeshContin.txt");
			cmesh->Print(arg4);
				
			Shell->SolveSist(an, cmesh,sim);
				
			TPZVec<REAL> QeM(9,0.);
			QeM = Shell->CalcCortMomento(cmesh);
		
			//cmesh->Print();
			delete cmesh;
			cmesh = 0;
			delete gmesh;
			gmesh = 0;
		*/	
			outfile <<"\n Termo de simetria, sim = " << sim <<endl;	
			outfile <<"=====Calculo com, nh = " << nh << ", Ndiv = "<< ndiv <<" e p = " << p <<endl;
			outfile<<" RefDirAnel = "<< nDiv_DirectRefAnel<<"   RefDirCasca = "<< nDiv_DirectRefCasca<<" e   RefDirPonto = "<< nDiv_DirectRefPonto<<endl<<endl;
		/*	outfile<<"RESULTADOS SEM INTERFACE"<<endl;
			outfile <<" MomentoR = "<<QeM[0]<<"  Erro Relativo ==> " << fabs((Mreal - QeM[0])/Mreal)*100.<<" %"<<endl;
			outfile <<" MomentoL = "<<QeM[2]<<"  Erro Relativo ==> " << fabs((Mreal - QeM[2])/Mreal)*100.<<" %"<<endl;
			outfile <<" MomentoMedia = "<<QeM[4]<<"  Erro Relativo ==> " << fabs((Mreal - QeM[4])/Mreal)*100.<<" %"<<endl;
			outfile <<endl;
			outfile <<" CortanteR = "<<QeM[1]<<"  Erro Relativo ==> " << fabs((Qreal - QeM[1])/Qreal)*100.<<" %"<<endl;
			outfile <<" CortanteL = "<<QeM[3]<<"  Erro Relativo ==> " << fabs((Qreal - QeM[3])/Qreal)*100.<<" %"<<endl;
			outfile <<" CortanteMedia= "<<QeM[5]<<"  Erro Relativo ==> " << fabs((Qreal - QeM[5])/Qreal)*100.<<" %"<<endl;
			outfile <<endl;
			outfile <<" T1zR = "<<QeM[6]<<endl;
			outfile <<" T1zL = "<<QeM[7]<<endl;
			outfile <<" T1z = "<<QeM[8]<<endl;
			outfile <<endl;
		
		*/
		
			//------------- RESOLUCAO MALHA DESCONTINUA ---------------------------------------------------------------------------------------
			TPZGeoMesh * gmesh2 = Shell->MalhaGeoGen(ndiv, nDiv_DirectRefCasca, nDiv_DirectRefAnel, nDiv_DirectRefPonto, true, 4);
							
			Shell->RefinamentoUniforme(*gmesh2, nh);
		
			ofstream arg3("gmeshdiscont.txt");
			gmesh2->Print(arg3);
		
			TPZCompMesh * cmesh_changed = Shell->MalhaCompMeshWithInterface(gmesh2, p,sim, pen);
			ofstream arg2("cmeshdiscont.txt");
			cmesh_changed->Print(arg2);
		
			int nEq = cmesh_changed->NEquations();
			cout << "\n Numero de Equacoes = " << nEq<<endl;
		
			TPZAnalysis an_ch(cmesh_changed);
			Shell->SolveSist(an_ch, cmesh_changed,sim);
			
			
			//cmesh_changed->Print();
		
			Shell->PrintInterface(cmesh_changed);
			
			TPZVec<REAL> QeM2(9,0.);
			QeM2 = Shell->CalcCortMomento(cmesh_changed);
					
			delete cmesh_changed;
			cmesh_changed = 0;
			delete gmesh2;
			gmesh2 = 0;
		
			outfile<<"RESULTADOS COM INTERFACE"<<endl;
			outfile <<" MomentoR = "<<QeM2[0]<<"  Erro Relativo ==> " << fabs((Mreal - QeM2[0])/Mreal)*100.<<" %"<<endl;
			outfile <<" MomentoL = "<<QeM2[2]<<"  Erro Relativo ==> " << fabs((Mreal - QeM2[2])/Mreal)*100.<<" %"<<endl;
			outfile <<" MomentoMedia = "<<QeM2[4]<<"  Erro Relativo ==> " << fabs((Mreal - QeM2[4])/Mreal)*100.<<" %"<<endl;
			outfile <<endl;
			outfile <<" CortanteR = "<<QeM2[1]<<"  Erro Relativo ==> " << fabs((Qreal - QeM2[1])/Qreal)*100.<<" %"<<endl;
			outfile <<" CortanteL = "<<QeM2[3]<<"  Erro Relativo ==> " << fabs((Qreal - QeM2[3])/Qreal)*100.<<" %"<<endl;
			outfile <<" CortanteMedia= "<<QeM2[5]<<"  Erro Relativo ==> " << fabs((Qreal - QeM2[5])/Qreal)*100.<<" %"<<endl;
			outfile <<endl;
			outfile <<" T1zR = "<<QeM2[6]<<endl;
			outfile <<" T1zL = "<<QeM2[7]<<endl;
			outfile <<" T1z = "<<QeM2[8]<<endl;
			outfile <<endl;
	
		}
	}
	
		
	std::ofstream outFileRefPat("refpatternsAgnaldo.txt");
	gRefDBase.WriteRefPatternDBase(outFileRefPat);
	
    return EXIT_SUCCESS;
}
