/***************************************************************************
 *   Copyright (C) 2006 by Philippe Devloo   *
 *   phil@fec.unicamp.br   *
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
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrprecond.h"
#include "pzdohrstructmatrix.h"
#include "pzstepsolver.h"
#include "pzcompel.h"

#include "pzelast3d.h"
#include "pzbndcond.h"

#include "tpzdohrassembly.h"

#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "tpzpairstructmatrix.h"
#include "pzviscoelastic.h"
#include "TPZTimer.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeotetrahedra.h"

#include "tpzarc3d.h"
//

#include "pzvtkmesh.h"

#include "pzlog.h"

#include <fstream>
#include <string>

#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
#endif
#include "tpzdohrmatrix.h"

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh);
void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh);
TPZGeoMesh *MalhaPredio();
TPZGeoMesh *MalhaCubo();
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);

using namespace std;


int main(int argc, char *argv[])
{
	/* Quando se está usando o tal log4cxx */
	//	InitializePZLOG("log4cxx.cfg");
	InitializePZLOG();
	
	TPZTimer total;
	total.start();
	std::cout << "COMECA O TEMPO"<< std::endl;
	
	int dim = 2;
	int maxlevel = 5;
	int sublevel = 3;
	int plevel = 1;
	TPZPairStructMatrix::gNumThreads = 2;
	int numthreads = 2;
	//	tempo.fNumthreads = numthreads;					// alimenta timeTemp com o numero de threads
	TPZGeoMesh *gmesh = 0;
	{
		TPZCompEl::SetgOrder(plevel);
		
		TPZAutoPointer<TPZCompMesh> cmesh;
		
		if(0)
		{
			TPZGenSubStruct sub(dim,maxlevel,sublevel);
			cmesh = sub.GenerateMesh();
			cmesh->SetDimModel(dim);
			gmesh = cmesh->Reference();
		}
		else 
		{
			dim = 3;
			gmesh = MalhaPredio();
			//gmesh = MalhaCubo();
			cmesh = new TPZCompMesh(gmesh);
			cmesh->SetDimModel(3);
			InsertViscoElasticity(cmesh);
			cmesh->AutoBuild();
		}
		
		std::cout << "Numero de equacoes " << cmesh->NEquations() << std::endl;
		
		int numthread_assemble = 16;
		int numthread_decompose = 16;
		TPZAutoPointer<TPZCompMesh> cmeshauto(cmesh);
		TPZDohrStructMatrix dohrstruct(cmeshauto,numthread_assemble,numthread_decompose);
		
		dohrstruct.IdentifyExternalConnectIndexes();
		
		std::cout << "Substructuring the mesh\n";
		//		TPZfTime timetosub; // init of timer
		dohrstruct.SubStructure(8);
		//		tempo.ft0sub = timetosub.ReturnTimeDouble();  // end of timer
		//		std::cout << tempo.ft0sub << std::endl;
		
		//		sub.SubStructure();
#ifdef LOG4CXX2
		{
			std::stringstream str;
			cmesh->Print(str);
			LOGPZ_DEBUG(logger,str.str());
		}
#endif
		
		
		dohrstruct.SetNumThreads(numthreads);
		
		TPZAutoPointer<TPZGuiInterface> gui;
		TPZFMatrix rhs(cmesh->NEquations(),1,0.);
		TPZAutoPointer<TPZMatrix> dohr = dohrstruct.CreateAssemble(rhs, gui);
		TPZAutoPointer<TPZMatrix> precond = dohrstruct.Preconditioner();
		
		
		TPZFMatrix diag(dohr->Rows(),1,5.), produto(dohr->Rows(),1);
		std::cout << "Numero de equacoes " << dohr->Rows() << std::endl;
		//		tempo.fNumEqCoarse = dohr->Rows();											// alimenta timeTemp com o numero de equacoes coarse
		dohr->Multiply(diag,produto);
		
		TPZDohrMatrix<TPZDohrSubstructCondense> *dohrptr = dynamic_cast<TPZDohrMatrix<TPZDohrSubstructCondense> *> (dohr.operator->());
		if (!dohrptr) {
			DebugStop();
		}
		dohrptr->AdjustResidual(produto);
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			produto.Print("O valor do produto", sout );
			diag.Print("O valor da diagonal",sout);
			LOGPZ_DEBUG(loggerconverge,sout.str())
		}
#endif
		diag.Zero();
		TPZStepSolver pre(precond);
		pre.SetMultiply();
		TPZStepSolver cg(dohr); 
		//  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const REAL tol,const int FromCurrent);
		
		cg.SetCG(500,pre,1.e-8,0);
		
		
		//		TPZfTime timetosolve; // init of timer
		cg.Solve(rhs,diag);
		//		tempo.ft6iter = timetosolve.ReturnTimeDouble(); // end of timer
		//		cout << "Total: " << tempo.ft6iter << std::endl;
		
		//		cout << "Tempos para multiplicacao: " << tempo.fMultiply << std::endl;
		//		cout << "Tempos para precondicionamento: " << tempo.fPreCond << std::endl;
		
		string FileName;
		FileName = "Times_in_Line.txt";
		ofstream OutputFile;
		
		//		bool shouldprint = tempo.NeedsHeader(FileName);			// verify the need of a header
		//		OutputFile.open(FileName.c_str(), ios::app);					// creates the file
		//		if (shouldprint == true) tempo.PrintHeader(OutputFile);		// prints the header if It is the first time the program is executed
		
		//		tempo.PrintLine(OutputFile);		// print all the information in one line
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			diag.Print("Resultado do processo antes do ajuste",sout);
			LOGPZ_INFO(loggerconverge,sout.str())
		}
#endif
		
		dohrptr->AddInternalSolution(diag);
		
		typedef std::list<TPZAutoPointer<TPZDohrSubstructCondense> > subtype;
		const subtype &sublist = dohrptr->SubStructures(); 
		subtype::const_iterator it = sublist.begin();
		int subcount=0;
		while (it != sublist.end()) 
		{
			TPZFMatrix subext,subu;
			dohrptr->fAssembly->Extract(subcount,diag,subext);
			(*it)->UGlobal(subext,subu);
			TPZCompMesh *submesh = SubMesh(cmeshauto, subcount);
			submesh->LoadSolution(subu);
			
			//Atualizando memoria do material
			std::map<int ,TPZAutoPointer<TPZMaterial> > materialmap(submesh->MaterialVec());
			std::map<int ,TPZAutoPointer<TPZMaterial> >::iterator itmat;
			for (itmat = materialmap.begin(); itmat != materialmap.end() ; itmat++) 
			{
				TPZAutoPointer<TPZMaterial> mat = itmat->second;
				TPZViscoelastic *vmat = dynamic_cast< TPZViscoelastic *> (mat.operator->());
				if(vmat)
				{
					vmat->SetUpdateMem();
				}
			}			
			subcount++;
			it++;
		}
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			diag.Print("Resultado do processo iterativo",sout);
			LOGPZ_INFO(loggerconverge,sout.str())
	   	}
#endif		
		
		TPZAutoPointer<TPZMaterial> mat = cmeshauto->FindMaterial(1);
		int nstate = mat->NStateVariables();
		int nscal = 0, nvec = 0;
		if(nstate ==1) 
		{
            nscal = 1;
		}
		else
		{
			nvec = 1;
		}
		TPZManVector<std::string> scalnames(nscal),vecnames(nvec);
		if(nscal == 1)
		{
			scalnames[0]="state";            
		}
		else
		{
			vecnames[0] = "state";
		}
		std::string postprocessname("dohrmann_visco.vtk");
		TPZVTKGraphMesh vtkmesh(cmesh.operator->(),dim,mat,scalnames,vecnames);
		vtkmesh.SetFileName(postprocessname);
		vtkmesh.SetResolution(1);
		int numcases = 1;
		
		
		// Iteracoes de tempo
		int istep = 0, nsteps = 2;
		vtkmesh.DrawMesh(numcases);
		vtkmesh.DrawSolution(istep, 1.);
		std::cout << "To seguindo!!!" << std::endl;
		for (istep = 1 ; istep < nsteps ; istep++)
		{
			TPZAutoPointer<TPZGuiInterface> guifake;
			dohrstruct.Assemble( rhs, guifake);
			cg.Solve(rhs,diag);			
			vtkmesh.DrawMesh(numcases);
			vtkmesh.DrawSolution(istep, 1.);			
		}
		
		
	}
	
	total.stop();
	std::cout << "TEMPO = " << total.seconds() << std::endl;
    
	
	
	delete gmesh;
	return EXIT_SUCCESS;
}


int main2(int argc, char *argv[])
{
	/* Quando se está usando o tal log4cxx */
	InitializePZLOG("log4cxx.cfg");
	
	/*
	 TPZFMatrix teste(2,2);
	 TPZFMatrix parte;
	 teste(0,0)=1;
	 teste(0,1)=2;
	 teste(1,0)=3;
	 teste(1,1)=4;
	 teste.GetSub(0,0,2,1,parte);
	 cout << parte << endl;
	 cout << "Hello, world!" << endl;
	 */
	
	
	/**
	 TPZDohrSubstruct meuobjeto;
	 TPZDohrMatrix *matriz = new TPZDohrMatrix();
	 TPZDohrPrecond *precond = new TPZDohrPrecond();
	 TPZStepSolver dohrprecond(precond);
	 dohrprecond.SetMultiply();
	 TPZStepSolver cg(matriz);
	 cg.SetCG(10,dohrprecond,1.e-7,1);
	 TPZFMatrix rhs,result;
	 cg.Solve(rhs,result);*/
	//meuobjeto.
	/*  int dim = 2;
	 TPZGenSubStruct sub(dim,6,3);*/
	//	int dim = 2;
	//	int maxlevel = 6;
	//	int sublevel = 3;
	//	int plevel = 3;
	int dim = 2;
	int maxlevel = 4;
	int sublevel = 1;
	int plevel = 2;
	TPZGenSubStruct sub(dim,maxlevel,sublevel);
	int nk = 8;
	//	int ik;
	//	for(ik=1; ik<nk; ik++)
	//	{
	//		sub.fK[ik] = 1.;//50.*ik;
	//	}
	//sub.fMatDist = TPZGenSubStruct::RandomMat;
	
	TPZCompEl::SetgOrder(plevel);
	
	sub.GenerateMesh();
	
	/*
	 TPZAutoPointer<TPZDohrAssembly> dohrassembly = new TPZDohrAssembly;
	 TPZDohrMatrix<TPZDohrSubstruct> *dohrptr = new TPZDohrMatrix<TPZDohrSubstruct>(dohrassembly);
	 TPZAutoPointer<TPZMatrix> dohr(dohrptr);
	 sub.InitializeDohr(dohr,dohrassembly);
	 
	 
	 // loop over the substructures
	 // This is a lengthy process which should run on the remote processor
	 //	void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstruct> substruct,  TPZDohrAssembly &dohrassembly);
	 
	 
	 dohrptr->Initialize();
	 #ifdef LOG4CXX
	 {
	 std::stringstream sout;
	 dohrptr->Print("DohrMatrix without condensation", sout);
	 LOGPZ_DEBUG(logger,sout.str())
	 }
	 #endif
	 TPZDohrPrecond<TPZDohrSubstruct> *precondptr = new TPZDohrPrecond<TPZDohrSubstruct>(*dohrptr,dohrassembly);
	 precondptr->Initialize();
	 TPZAutoPointer<TPZMatrix> precond(precondptr);
	 
	 */
	
	TPZAutoPointer<TPZDohrAssembly> dohrassembly2 = new TPZDohrAssembly;
	TPZDohrMatrix<TPZDohrSubstructCondense> *dohrptr2 = new TPZDohrMatrix<TPZDohrSubstructCondense>(dohrassembly2);
	dohrptr2->SetNumThreads(4);
	TPZAutoPointer<TPZMatrix> dohr2(dohrptr2);
	sub.InitializeDohrCondense(dohr2,dohrassembly2);
	dohrptr2->Initialize();
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		dohr2->Print("The dohr matrix condensed",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	
	
	
#ifdef LOG4CXX
	std::stringstream sout;
	sout << "Three dimensional substructures, maxlevel " << maxlevel << " level of substructures " << sublevel << std::endl;
	sout << "Number of substructures " << dohrptr2->SubStructures().size() << std::endl;
	sout << "Interpolation order " << plevel;
	LOGPZ_DEBUG(loggerconverge,sout.str());
#endif
	
	
	TPZDohrPrecond<TPZDohrSubstructCondense> *precondptr2 = new TPZDohrPrecond<TPZDohrSubstructCondense>(*dohrptr2,dohrassembly2);
	precondptr2->Initialize();
	
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Printing after creating the preconditioner\n";
		dohrptr2->Print("After creating the preconditioner",sout);
		LOGPZ_DEBUG(loggerconverge,sout.str());
	}
#endif
	
	
	
	TPZAutoPointer<TPZMatrix> precond2(precondptr2);
	
	
	TPZFMatrix diag(dohr2->Rows(),1,5.), produto(dohr2->Rows(),1), produto2(dohr2->Rows(),1);
	precondptr2->Multiply(diag,produto);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		produto.Print("O valor do produto", sout );
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif

	precondptr2->Multiply(diag,produto2);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		produto2.Print("O valor do produto2", sout );
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif
	//#define TOTAL
#ifdef TOTAL
	{
		int dim=dohr->Rows();
		
		TPZFMatrix teste1(dim,dim);
		teste1.Zero();
		TPZFMatrix teste2(dim,dim);
		teste2.Zero();
		
		int i,j;
		TPZFMatrix col(dim,1); //column of the identity matrix
		TPZFMatrix resul(dohr->Rows(),1);
		for (i=0;i<dim;i++) {
			col.Zero();
			col(i,0)=1;
			precondptr->MultAdd(col,col,resul,1,0,0,1);
			for (j=0;j<dim;j++) {
				teste1(i,j) = resul(j,0);
				teste2(i,j) = resul(j,0);
			}
		}
		teste1.Transpose();
		teste1 -= teste2;
		std::cout << "Norma da diferenca das matrizes " << Norm(teste1) << std::endl;
	}
#endif
	std::cout << "Numero de equacoes " << dohr2->Rows() << std::endl;
	//  produto.Print("The value of the product is");
#ifndef MAKEINTERNAL
	diag(0,0) = 0.;
#endif
	dohr2->Multiply(diag,produto);
	dohrptr2->AdjustResidual(produto);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		produto.Print("O valor do produto", sout );
		diag.Print("O valor da diagonal",sout);
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif
	diag.Zero();
	TPZStepSolver pre(precond2);
	pre.SetMultiply();
	TPZStepSolver cg(dohr2);
	//  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const REAL tol,const int FromCurrent);
	
	cg.SetCG(500,pre,1.e-8,0);
	cg.Solve(produto,diag);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		diag.Print("Resultado do processo antes do ajuste",sout);
		LOGPZ_INFO(loggerconverge,sout.str())
	}
#endif
	
	dohrptr2->AddInternalSolution(diag);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		diag.Print("Resultado do processo iterativo",sout);
		LOGPZ_INFO(loggerconverge,sout.str())
	}
#endif
	//diag.Print("Resultado do solve");
	/* Solve
	 TPZFMatrix *teste = new TPZFMatrix(2,2);
	 (*teste)(0,0)=1;
	 (*teste)(0,1)=2;
	 (*teste)(1,0)=3;
	 (*teste)(1,1)=4;
	 TPZStepSolver coef;
	 coef.SetMatrix(teste);
	 coef.SetDirect(ELU);
	 TPZFMatrix resul(2,2);
	 resul(0,0)=2;
	 resul(0,1)=3;
	 resul(1,0)=4;
	 resul(1,1)=5;
	 TPZFMatrix res(2,2);
	 coef.Solve(resul,res);
	 cout << res << endl;*/
	
	
	return EXIT_SUCCESS;
}

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1;
	REAL E = 1.e6;
	REAL poisson = 0.3;
	TPZManVector<REAL> force(3,0.);
	force[1] = 20.;
	TPZElasticity3D *elast = new TPZElasticity3D(nummat,E,poisson,force);
	TPZAutoPointer<TPZMaterial> elastauto(elast);
	TPZFMatrix val1(3,3,0.),val2(3,1,0.);
	TPZBndCond *bc = elast->CreateBC(elastauto, -1, 0, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto(bc);
	mesh->InsertMaterialObject(elastauto);
	mesh->InsertMaterialObject(bcauto);
}

void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1;
	REAL Ela = 1.e6;
	REAL poisson = 0.2;
	TPZManVector<REAL> force(3,0.);
	force[1] = 20.;
	REAL lambdaV = 0, muV = 0, alphaT = 0;
	lambdaV = 111.3636;
	muV = 455.4545;
	alphaT = 0.1;	
	TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat,Ela,poisson,lambdaV,muV,alphaT,force);
	TPZAutoPointer<TPZMaterial> viscoelastauto(viscoelast);
	TPZFMatrix val1(3,3,0.),val2(3,1,0.);
	TPZBndCond *bc = viscoelast->CreateBC(viscoelastauto, -1, 0, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto(bc);
	mesh->InsertMaterialObject(viscoelastauto);
	mesh->InsertMaterialObject(bcauto);	
}

TPZGeoMesh *MalhaPredio()
{
	//int nBCs = 1;
	int numnodes=-1;
	int numelements=-1;
	
	string FileName;
	FileName = "../8andares02.txt";
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		ifstream read (FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements") countelements = false;
			if(countelements) numelements++;
		}
	}
	
	TPZGeoMesh * gMesh = new TPZGeoMesh;
	
	gMesh -> NodeVec().Resize(numnodes);
	
	TPZVec <int> TopolTetra(4);
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int nodeId = 0, elementId = 0, matElId = 1;
	
	ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	for(in=0; in<numnodes; in++)
	{ 
		read >> nodeId;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		Node[nodeId-1].SetNodeId(nodeId);
		Node[nodeId-1].SetCoord(0,nodecoordX);
		Node[nodeId-1].SetCoord(1,nodecoordY);
		Node[nodeId-1].SetCoord(2,nodecoordZ);
		gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
		
		
	}
	
	{
		
		read.close();
		read.open(FileName.c_str());
		
		
		
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		
		int el;
		int matBCid = -1;
		//std::set<int> ncoordz; //jeitoCaju
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> TopolTetra[0]; //node 1
			read >> TopolTetra[1]; //node 2
			read >> TopolTetra[2]; //node 3
			read >> TopolTetra[3]; //node 4
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
			TopolTetra[0]--;
			TopolTetra[1]--;
			TopolTetra[2]--;
			TopolTetra[3]--;
			
			int index = el;
			
			TPZGeoEl * tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		}
		
		gMesh->BuildConnectivity();
		// Colocando as condicoes de contorno

		for(el=0; el<numelements; el++)
		{
			TPZManVector <TPZGeoNode,4> Nodefinder(4);
			TPZManVector <REAL,3> nodecoord(3);
			TPZGeoEl *tetra = gMesh->ElementVec()[el];
			// na face z = 0
			TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				int pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[2] == 0.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = pos;
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,matBCid);		
			}
		}
	}
	
	
	
	// identificando as superficies que terao cond de contorno. Coord z dos 3 nos = 0
	//	for (int el = 0; el < numnodes-1; el++) 
	//	{
	//		Nodefind[el] = gMesh->NodeVec()[el];
	//
	//	}
	//	Nodefind.Print(std::cout);
	//	std::cout.flush();
	
	//TPZGeoElBC(TPZGeoEl *el,int side,int matid, TPZGeoMesh &mesh);
	//TPZGeoElBC(TPZGeoElSide &elside,int matid, TPZGeoMesh &mesh);
	
	ofstream arg("malhaPZ.txt");
	gMesh->Print(arg);
	ofstream predio("GeoPredio.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true); 
	
	return gMesh;
	
}


TPZGeoMesh *MalhaCubo()
{
	//int nBCs = 1;
	int numnodes=-1;
	int numelements=-1;
	
	string FileName;
	FileName = "../cube1.txt";
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		ifstream read (FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements") countelements = false;
			if(countelements) numelements++;
		}
	}
	
	TPZGeoMesh * gMesh = new TPZGeoMesh;
	
	gMesh -> NodeVec().Resize(numnodes);
	
	TPZManVector <int> TopolTetra(4);
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int nodeId = 0, elementId = 0, matElId = 1;
	
	ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	for(in=0; in<numnodes; in++)
	{ 
		read >> nodeId;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		Node[nodeId-1].SetNodeId(nodeId);
		Node[nodeId-1].SetCoord(0,nodecoordX);
		Node[nodeId-1].SetCoord(1,nodecoordY);
		Node[nodeId-1].SetCoord(2,nodecoordZ);
		gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
		
		//(Node, idbcnode, *gMesh);	
		
	}
	
	{
		read.close();
		read.open(FileName.c_str());
		
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		
		int el;
		int neumann1 = -4, neumann2 = -5, dirp2 = -6;
        int index = 0;
		//std::set<int> ncoordz; //jeitoCaju
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> TopolTetra[0]; //node 1
			read >> TopolTetra[1]; //node 2
			read >> TopolTetra[2]; //node 3
			read >> TopolTetra[3]; //node 4
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
			TopolTetra[0]--;
			TopolTetra[1]--;
			TopolTetra[2]--;
			TopolTetra[3]--;
			
			int index = el;
			
			//TPZGeoEl * tetra = gMesh->CreateGeoElement(ETetraedro, TopolTetra, matElId, index);
			TPZGeoEl * tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
			
		}
		
		gMesh->BuildConnectivity();
		
		// Colocando as condicoes de contorno
		for(el=0; el<numelements; el++)
		{
			
			
			TPZManVector <TPZGeoNode,4> Nodefinder(4);
			TPZManVector <REAL,3> nodecoord(3);
			TPZGeoEl *tetra = gMesh->ElementVec()[el];
			
			// na face x = 1
			TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				int pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[0] == 1.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = pos;
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,neumann1);		
			}
			
			// Na face x = -1
			ncoordzVec.Resize(0);
			sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				int pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];
				
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[0] == -1.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = pos;
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,neumann2);		
			}
			
		}
		
		TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
		yz[0] = 1.;
		z[2] = -1;
		int bcidxyz = -1, bcidyz = -2, bcidz = -3;
		SetPointBC(gMesh, xyz, bcidxyz);
		SetPointBC(gMesh, yz, bcidyz);
		SetPointBC(gMesh, z, bcidz);
		
		//gMesh->BuildConnectivity();
	}
	
	ofstream arg("malhaPZ.txt");
	gMesh->Print(arg);
	
	std::ofstream out("Cube.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gMesh, out, true);
	
	return gMesh;
}

/// Generate a boundary geometric element at the indicated node
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
{
	// look for an element/corner node whose distance is close to start
	TPZGeoNode *gn1 = gr->FindNode(x);
	int iel;
	int nelem = gr->ElementVec().NElements();
	TPZGeoEl *gel;
	for (iel = 0; iel<nelem; iel++) {
		gel = gr->ElementVec()[iel];
		if(!gel) continue;
		int nc = gel->NCornerNodes();
		int c;
		for (c=0; c<nc; c++) {
			TPZGeoNode *gn = gel->NodePtr(c);
			if (gn == gn1) {
				break;
			}
		}
		if (c<nc) {
			TPZGeoElBC(gel, c, bc);
			return;
		}
	}
}
