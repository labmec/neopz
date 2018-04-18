 /***************************************************************************
 *   Copyright (C) 2009 by joao *
 *   joao@joao-laptop   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"


#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "tpzsparseblockdiagonalstructmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzbndmat.h" 
#include "pzfstrmatrix.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "tpzautopointer.h"
#include "math.h"
#include "pzlog.h"
#include <sstream>
#include <fstream>
#include "tpzcopysolve.h"

#include "TPZTimer.h"
#include "TPZRefPattern.h"
#include "pzgeoquad.h"
#include "pzespmat.h"
#include "pzgeoelbc.h"



// para o projeto LNCC1 com esse main consideramos o problema de poisson em duas dimensões. Embora mais simples esse projeto ajuda a entender o funcionamento básico do FEM e as funcionalidades do PZ para tal

using namespace std;
//#undef LOG4CXX
#ifdef LOG4CXX
static LoggerPtr loggermain(Logger::getLogger("pz.lncc1"));
#endif



/*----------DECLARAÇÃO DOS METODOS-----------*/
TPZGeoMesh * MalhaGeo ( const int h );
TPZCompMesh *CreateCompMesh ( TPZGeoMesh &gmesh, int porder );
void Solve ( TPZAnalysis &an );
void SolveIterative(TPZAnalysis &an);
void PosProcess(TPZAnalysis& an_p1,ofstream& AE, int h, int plncc1);
void PrintGeoMeshVTKWithSomeData(TPZGeoMesh *gmesh,char *filename);

//vetor de carga e condições de contorno
//static REAL pi=4.*atan(1.);
void Forcing1( const TPZVec<REAL> &x, TPZVec<REAL> &disp );
void CC1 ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC2 ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC3 ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC4 ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void ExactSolution( const TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix<REAL> &deriv);

int main() 
{
	InitializePZLOG("/Users/Joao/WORK/PROJETOS/NeoPZ/Projects/LNCC1/log4cxx.cfg");
	int h=0;
	int p=1;
	
	std::ofstream AE ( "results.txt" );
	AE.precision(16);
	
	for(h=1;h<2;h++)
	{
        std::cout<< "h= "<< h<< std::endl;
	TPZGeoMesh *gmesh = MalhaGeo(h);

    TPZCompMesh *cmesh = CreateCompMesh(*gmesh,p);
	
	TPZAnalysis an(cmesh);
		
	Solve( an );

	PosProcess(an, AE, h, p);
	}
	
	return EXIT_SUCCESS;
} 


TPZGeoMesh * MalhaGeo/*Quadrilateros*/( const int h ){
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}}; 
	TPZGeoEl *elvec[1];//[2];
	int nnode = 4;
	int nelem = 1;//2;
	int nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}
	
	int el;
	
	//int nodindAll[2][3]={{0,1,2},{2,3,0}};
	int nodindAll[1][4]={{0,1,2,3}};
	for ( el=0; el<nelem; el++ )
	{
		TPZVec<int> nodind(4);//(3);
		nodind[0]=nodindAll[el][0];
		nodind[1]=nodindAll[el][1];
		nodind[2]=nodindAll[el][2];
		nodind[3]=nodindAll[el][3];		
		int index;
		//elvec[el] = gmesh->CreateGeoElement (ETriangle,nodind,1,index );
		elvec[el] = gmesh->CreateGeoElement (EQuadrilateral,nodind,1,index );		
	}
//	int nodindAll[1][4]={{0,1,2,3}};
//	for ( el=0; el<nelem; el++ )
//	{
//		TPZVec<int> nodind(4);
//		nodind[0]=nodindAll[el][0];
//		nodind[1]=nodindAll[el][1];
//		nodind[2]=nodindAll[el][2];
//		nodind[3]=nodindAll[el][3];
//		int index;
//		elvec[el] = gmesh->CreateGeoElement (EQuadrilateral,nodind,1,index );
//	}
	
	
	gmesh->BuildConnectivity();
	
	TPZGeoElBC gbc1( elvec[0],4,-1);
	TPZGeoElBC gbc2( elvec[0],5,-2);
	TPZGeoElBC gbc3( elvec[0],6,-3);
	TPZGeoElBC gbc4( elvec[0],7,-4);
	
//	TPZGeoElBC gbc1 ( elvec[0],3,-1);
//	TPZGeoElBC gbc2 ( elvec[0],4,-2);
//	TPZGeoElBC gbc3 ( elvec[1],3,-3);
//	TPZGeoElBC gbc4 ( elvec[1],4,-4);
	
//	TPZGeoElBC gbc5 ( elvec[0],0,-5);// um ponto (a origem)
	

	
	///Refinamento uniforme
	for ( int ref = 0; ref < h; ref++ )
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ )
		{
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
		}//for i
	}//ref
    
    ofstream arcg ( "gmesh.txt" );
	gmesh->Print ( arcg );
    
	return gmesh;
}

TPZGeoMesh * MalhaGeo2Triang( const int h ){
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}}; 
	TPZGeoEl *elvec[2];
	int nnode = 4;
	int nelem = 2;
	int nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}
	
	int el;
	
	int nodindAll[2][3]={{0,1,2},{2,3,0}};
	//int nodindAll[1][4]={{0,1,2,3}};
	for ( el=0; el<nelem; el++ )
	{
		TPZVec<int> nodind(3);
		nodind[0]=nodindAll[el][0];
		nodind[1]=nodindAll[el][1];
		nodind[2]=nodindAll[el][2];
//		nodind[3]=nodindAll[el][3];		
		int index;
		elvec[el] = gmesh->CreateGeoElement (ETriangle,nodind,1,index );
		//elvec[el] = gmesh->CreateGeoElement (EQuadrilateral,nodind,1,index );		
	}
	
	gmesh->BuildConnectivity();
	
//	TPZGeoElBC gbc1 ( elvec[0],4,-1);
//	TPZGeoElBC gbc2 ( elvec[0],5,-2);
//	TPZGeoElBC gbc3 ( elvec[0],6,-3);
//	TPZGeoElBC gbc4 ( elvec[0],7,-4);
	
	TPZGeoElBC gbc1 ( elvec[0],3,-1);
	TPZGeoElBC gbc2 ( elvec[0],4,-2);
	TPZGeoElBC gbc3 ( elvec[1],3,-3);
	TPZGeoElBC gbc4 ( elvec[1],4,-4);
	
	//	TPZGeoElBC gbc5 ( elvec[0],0,-5);// um ponto (a origem)
	
	
	///Refinamento uniforme
	for ( int ref = 0; ref < h; ref++ )
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ )
		{
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			gel->Divide ( filhos );
		}//for i
	}//ref
	
	ofstream arcg ( "gmesh.txt" );
	gmesh->Print ( arcg );

	return gmesh;
}

//malha Cross triangulos
TPZGeoMesh * MalhaGeoCross( const int h ){
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	REAL co[5][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.},{.5,.5}}; 
	int nnode = 5;
	int nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}
	
	int nelem = 4;
	TPZGeoEl *elvec[nelem];
	int el;
	
	int nodindAll[4][3]={{0,1,4},{1,2,4},{2,3,4},{3,0,4}};
	for ( el=0; el<nelem; el++ )
	{
		TPZVec<int> nodind(3);
		nodind[0]=nodindAll[el][0];
		nodind[1]=nodindAll[el][1];
		nodind[2]=nodindAll[el][2];	
		int index;
		elvec[el] = gmesh->CreateGeoElement (ETriangle,nodind,1,index, 0 );
	}
	
	gmesh->BuildConnectivity();
	
	//	TPZGeoElBC gbc1 ( elvec[0],4,-1);
	//	TPZGeoElBC gbc2 ( elvec[0],5,-2);
	//	TPZGeoElBC gbc3 ( elvec[0],6,-3);
	//	TPZGeoElBC gbc4 ( elvec[0],7,-4);
	
	TPZGeoElBC gbc1 ( elvec[0],3,-1);
	TPZGeoElBC gbc2 ( elvec[1],3,-2);
	TPZGeoElBC gbc3 ( elvec[2],3,-3);
	TPZGeoElBC gbc4 ( elvec[3],3,-4);
	
	//	TPZGeoElBC gbc5 ( elvec[0],0,-5);// um ponto (a origem)
	
	
	///Refinamento uniforme
	for ( int ref = 0; ref < h; ref++ )
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ )
		{
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			gel->Divide ( filhos );
		}//for i
	}//ref
	
	
	
	ofstream arcg ( "gmesh.txt" );
	gmesh->Print ( arcg );
	
	return gmesh;
}
 
TPZCompMesh *CreateCompMesh ( TPZGeoMesh &gmesh, int porder ){
	TPZCompEl::SetgOrder ( porder );
	TPZCompMesh *result = new TPZCompMesh( &gmesh );
	result->SetDimModel ( 2 );
	
	result->SetAllCreateFunctionsDiscontinuous();
//	result->SetAllCreateFunctionsContinuous();//default, não precisaria ser setado novamente
	
	TPZMatPoisson3d *material ;
	material = new TPZMatPoisson3d ( 1,2 );
	
	TPZVec<REAL> convdir(3);//direção da convecção
	convdir[0]=0.0;//sqrt(2.)/2.;
	convdir[1]=0.0;//sqrt(2.)/2.;
	REAL diff= 1.;
	REAL conv=0.;
	material-> SetParameters(diff, conv, convdir);
    TPZAutoPointer<TPZFunction<STATE> > ExactSol = new TPZDummyFunction<STATE> (ExactSolution);
	material->SetForcingFunctionExact(ExactSol);
    //material->SetNonSymmetric();
	//material->SetSymmetric();
	TPZMaterial *mat ( material );


	
	TPZAutoPointer<TPZFunction<STATE> > BC1 = new TPZDummyFunction<STATE> (CC1);
	TPZAutoPointer<TPZFunction<STATE> > BC2 = new TPZDummyFunction<STATE> (CC2);
	TPZAutoPointer<TPZFunction<STATE> > BC3 = new TPZDummyFunction<STATE> (CC3);
	TPZAutoPointer<TPZFunction<STATE> > BC4 = new TPZDummyFunction<STATE> (CC4);		

	TPZAutoPointer<TPZFunction<STATE> > LoadVector = new TPZDummyFunction<STATE> (Forcing1);
	material->SetForcingFunction ( LoadVector );
	result->InsertMaterialObject ( mat );
	
	TPZFMatrix<REAL> val1 ( 1,1,0. ), val2 ( 1,1,0. );// 0 é Dirichlet, 1 é Neumann, 2 é Robin(implementada apenas no Contínuo)

	TPZMaterial *bnd1 = material->CreateBC ( mat,-1,1, val1, val2 );
	TPZMaterial *bnd2 = material->CreateBC ( mat,-2,0, val1, val2 );
	TPZMaterial *bnd3 = material->CreateBC ( mat,-3,1, val1, val2 );
	TPZMaterial *bnd4 = material->CreateBC ( mat,-4,0, val1, val2 );
	
	bnd1->SetForcingFunction ( BC1 );
	bnd2->SetForcingFunction ( BC2 );
	bnd3->SetForcingFunction ( BC3 );
	bnd4->SetForcingFunction ( BC4 );
	
	result->InsertMaterialObject ( bnd1 ); 
	result->InsertMaterialObject ( bnd2 );
	result->InsertMaterialObject ( bnd3 );
	result->InsertMaterialObject ( bnd4 );
	
	result->AutoBuild();

	
	//
	ofstream arc ( "cmesh.txt" );
	arc << "NEquation = "<< result->NEquations() << endl;
	result->Print ( arc );
	ofstream arcg ( "gmesh.txt" );
	result->Reference()->Print ( arcg );
	//
	result->SetName("CMesh1");
	//
	return result;
}

void SolveIterative( TPZAnalysis &an ){
	TPZCompMesh *malha = an.Mesh();
	
	//TPZBandStructMatrix mat(malha);
	//TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);
	//TPZParFrontStructMatrix mat(malha);
	//TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
	//TPZFStructMatrix mat( malha );
	TPZSpStructMatrix mat( malha );
	TPZStepSolver<REAL> solv;
		
	TPZFMatrix<REAL> fakeRhs(malha->NEquations(),1);
	//TPZSkylineStructMatrix PrecondMatrix(malha);
	//TPZParSkylineStructMatrix PrecondMatrix(malha, 1);
	TPZFrontStructMatrix<TPZFrontNonSym<REAL> >PrecondMatrix(malha);
//	PrecondMatrix.SetNumberOfThreads(4);
	//TPZBlockDiagonalStructMatrix PrecondMatrix(malha);
	
	TPZAutoPointer<TPZGuiInterface> guiInterface1 = new TPZGuiInterface();//
	TPZStepSolver<REAL> precond(PrecondMatrix.CreateAssemble(fakeRhs, guiInterface1));

	precond.SetDirect(ELU);//ELU //ECholesky // ELDLt
	
 	solv.SetGMRES( 100000, 50, precond, 1.e-15, 0 );
	//solv.SetCG(10000, precond, 1.e-16, 0);	
	//solv.SetJacobi(1000, 1.e-15, 0)
	//solv.SetBiCGStab();
	//solv.SetSOR();
	//solv.SetSSOR();
	
	cout << "ELU " << endl;
	an.SetSolver ( solv );
	an.SetStructuralMatrix ( mat );
	cout << endl;
	an.Solution().Redim ( 0,0 );
	cout << "Assemble " << endl;
	an.Assemble();
	an.Solve();
	
	cout << endl;
	cout << "No equacoes = " << malha->NEquations() << endl;
	cout << "Banda = " << malha->BandWidth() << endl;
}

void Solve ( TPZAnalysis &an ){
	TPZCompMesh *malha = an.Mesh();

	TPZBandStructMatrix mat(malha);
	//TPZSkylineStructMatrix mat(malha);// requer decomposição simétrica, não pode ser LU!
	//TPZBlockDiagonalStructMatrix mat(malha);//ok!
	//TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
	//TPZFStructMatrix mat( malha );// ok! matriz estrutural cheia
	//TPZSpStructMatrix mat( malha );//matriz estrutural esparsa (???? NÃO FUNCIONOU !!!!!!!!!!)
	TPZStepSolver<REAL> solv;
	solv.SetDirect (  ELU );//ECholesky);// ELU , ELDLt , 

	
//	cout << "ELDLt " << endl;
	an.SetSolver ( solv );
	an.SetStructuralMatrix ( mat );
	cout << endl;
	an.Solution().Redim ( 0,0 );
	cout << "Assemble " << endl;
	an.Assemble();
//	std::ofstream fileout("rigidez.txt");
//	an.Solver().Matrix()->Print("Rigidez", fileout, EMathematicaInput);
	
	an.Solve();
	cout << endl;
	cout << "No equacoes = " << malha->NEquations() << endl;
}

void PosProcess(TPZAnalysis& an ,ofstream& AE,  int h, int p){
	
	AE<<"------------------------------------"<<endl;
	AE << "hxh, h="<< h<< ", p="<< p<< endl;
	AE << "DOF = " <<  an.Mesh()->NEquations()    << endl;
	AE<<"------------------------------------"<<endl;
	
	/// VETOR SOLUCAO
	/// 1. Coeficientes da Solucao Numerica 
	//         ofstream filep1("Solution_primal1.out");
	//         TPZFMatrix<REAL> toprintp1 = an_p1.Solution();
	//         toprintp1.Print("solution", filep1);
	///GRAFICOS
	///2. Saida para dx ou vtk
	TPZVec<std::string> scalar_names(2);
	TPZVec<std::string> vec_names(0);
	
	scalar_names[0] = "Solution";
	scalar_names[1] = "ExactSolution";	

	char buf[256] ;
	sprintf(buf,"iteracao_p%d_h%d.vtk",p,h); 
	int dim=2;
	an.DefineGraphMesh(dim, scalar_names, vec_names,buf);
	int ExtraResolution = 0;
	an.PostProcess(ExtraResolution);// o parametro é a resolução que será acrescida na visualização.
	
	///ERROS
	an.SetExact(ExactSolution);
	TPZVec<REAL> pos;
	an.PostProcess(pos,AE); // Calcula os erros.(pzanalysis.cpp)
	
}

void PrintGeoMeshVTKWithSomeData(TPZGeoMesh *gmesh,char *filename) {
	int i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		if(gmesh->ElementVec()[i])
			DataElement[i] = (gmesh->ElementVec()[i])->Id();
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}

void Forcing1( const TPZVec<REAL> &x, TPZVec<REAL> &disp )
{
	disp[0]=1.;
}
void CC1 ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-1
{
	f[0] = 0.;
}
void CC2 ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-2
{
	f[0] = 0.;
}
void CC3 ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-3
{
	f[0] = 0.;
}
void CC4 ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-4
{
	f[0] = 1.;
}
void ExactSolution( const TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix<REAL> &deriv)
{
	u[0] = 1.-x[0];
}
	

