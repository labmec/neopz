/*
 *  substructeigen.cpp
 *  SubstructEigen
 *
 *  Created by Philippe Devloo on 11/15/08.
 *  Copyright 2008 UNICAMP. All rights reserved.
 *
 */

#include <iostream>
#include <sstream>
#include "substructeigen.h"
#include "pzfmatrix.h"
#include "pzgengrid.h"
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzmat2dlin.h"
#include "pzbndcond.h"
#include "pzsubcmesh.h"
#include "pzelmat.h"
#include "pzinterpolationspace.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzanalysis.h"
#include "pzgraphmesh.h"

#include "pzstepsolver.h"
#include "GenerateSquare.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.main"));
static LoggerPtr loggerres(Logger::getLogger("pz.main.result"));
#endif

//#define EXTCOMP

#include "tconfig.h"
#include "tlistconfig.h"

void CriarMateriaisSteklov(TPZCompMesh &mesh, TConfig &cf);

void CriarMateriaisFlux(TPZCompMesh &mesh, TConfig &cf);

TPZFMatrix ProjecaoL2(TPZCompMesh &mesh, int subelindex, TConfig &cf);

void TransferirElementos(TPZCompMesh &mesh, int materialid, int &elindex);

void ResequenceAsSubmesh(TPZSubCompMesh *submesh);

void CondensedEquations(TPZCompMesh &comp, int materialid , TPZVec<int> &origin);

void RefinarMalha(TPZCompMesh &comp,int materialid, int nref, int porder);

void RefineSquare(TPZCompMesh &cmesh, int nref, int porderstart);

void PlotEigenvectors(TPZSubCompMesh *element, TListConfig &listconfig, TConfig &cf, TPZFMatrix &solution , TPZVec<int> &condense);

void ExecuteConfig(TListConfig &listconfig, TConfig &cf, TPZGeoMesh &init);

int main()
{
	// hrefborder pborder href porder
	int allconf[][4] =
	{
		{0,1,0,1},
		{0,2,0,2},
		{0,3,0,3},
		{0,4,0,4},
		{0,5,0,5},
		{0,3,6,4},
		{0,3,7,4},
		{0,3,8,4}
	};
	int numconf = 5;

#ifdef LOG4CXX
	InitializePZLOG();
#endif
	TPZGeoMesh geo;
	geo.InitializeRefPatterns();

	TListConfig listconfig("CurvedQuadrilateralTEST",TListConfig::CurvedQuadrilateral);
	int ic;
	for(ic=0; ic<numconf; ic++)
	{
		TConfig &conf = listconfig.GenerateConfig(allconf[ic]);
		std::cout << "Executing " << conf.AsList() << std::endl;
		ExecuteConfig(listconfig,conf,geo);
	}

#ifdef LOG4CXX
	{
		LOGPZ_INFO(loggerres,listconfig.PrintFinal());
	}
#endif

	return 0;
}

void ExecuteConfig(TListConfig &listconfig, TConfig &cf, TPZGeoMesh &cp)
{
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "(*\nnref = "<< cf.nref << "\nporder = " << cf.porder << "\nnrefborder = "<< cf.nrefborder << "\npborder = " << cf.pborder << "*)";
		LOGPZ_INFO(loggerres,sout.str());
	}
#endif

	TPZGeoMesh *geo = listconfig.GenerateMesh(cp,cf);

#ifdef LOG4CXX
	{
	std::stringstream sout;
	geo->Print(sout);
	LOGPZ_DEBUG(logger,sout.str());
	}
#endif

	TPZCompMesh *comp = new TPZCompMesh(geo);
	CriarMateriaisFlux(*comp,cf);
	comp->SetDefaultOrder(cf.pborder);
	comp->AutoBuild();
	int nel = comp->NElements();
	int iel;
	for(iel=0; iel<nel; iel++)
	{
		TPZCompEl *cel = comp->ElementVec()[iel];
		if(cel && cel->Reference()->Dimension() == 2)
		{
			TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
			intel->PRefine(cf.porder);
		}
	}
	comp->ExpandSolution();
	// NOW WE NEED TO ADJUST THE INTERPOLATION ORDER OF THE BOUNDARY ELEMENTS
	if(listconfig.fConfigType == TListConfig::SquareAdaptive)
	{
#ifdef LOG4CXX
		{
			std::stringstream sout;
			comp->Print(sout);
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif
		RefineSquare(*comp, cf.nref, cf.porder);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			comp->Reference()->Print(sout);
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif
	}

#ifdef LOG4CXX
	{
	std::stringstream sout;
	comp->Print(sout);
	LOGPZ_DEBUG(logger,sout.str());
	}
#endif

	if(listconfig.fConfigType != TListConfig::SquareAdaptive)
	{

		RefinarMalha(*comp,cf.materialborder1,cf.nrefborder,cf.pborder);
		RefinarMalha(*comp,cf.materialborder2,cf.nrefborder,cf.pborder);
		RefinarMalha(*comp,cf.materialcondense,cf.nrefborder,cf.pborder);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Apos refinamento da fronteira\n";
			geo->Print(sout);
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Apos refinamento da fronteira\n";
			comp->Print(sout);
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif
		RefinarMalha(*comp,cf.materialid,cf.nref,cf.porder);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Apos refinamento\n";
			comp->Print(sout);
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif
	}
	int subelindex = 0;
	TransferirElementos(*comp, cf.materialid, subelindex);

#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Apos transferencia\n";
		comp->Print(sout);
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	// criar um objeto analysis para calcular uma solucao
	TPZAnalysis analysis(comp);
	TPZSubCompMesh *subc = dynamic_cast<TPZSubCompMesh *> (comp->ElementVec()[subelindex]);
	ResequenceAsSubmesh(subc);
	TPZAutoPointer<TPZStructMatrix> str = new TPZFStructMatrix(comp);
	analysis.SetStructuralMatrix(str);
	TPZStepSolver step;
	step.SetDirect(ECholesky);
	analysis.SetSolver(step);

	// reordenar as equacoes para por fazer a enumeracao global corresponder a equacao do elemento subestruturado
	// identificar as equacoes correspondente a condicao dirichlet
	TPZManVector<int> origin;
	CondensedEquations(*comp,cf.materialcondense,origin);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << cf.condense() << " = {";
		int i;
		for(i=0; i<origin.NElements()-1; i++)
		{
			sout << origin[i] << ",";
		}
		if(i<origin.NElements())
		{
			sout << origin[i];
		}
		sout << "};";
		LOGPZ_INFO(loggerres,sout.str())
	}
#endif
	// run the problem, create the plotfile
	analysis.Run();
	TPZVec<std::string> scalnames(1),vecnames(0);
	scalnames[0] = "state";
	analysis.DefineGraphMesh(2,scalnames,vecnames,cf.SolutionFileName(listconfig.fName));
	analysis.PostProcess(4);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << cf.rhs() << " = ";
		analysis.Rhs().Print(" ",sout,EMathematicaInput);
		sout << cf.rhs() << " = Flatten[" << cf.rhs() << "];\n";
		sout << "(*\nnref = "<< cf.nref << "\nporder = " << cf.porder << "\nnrefborder = "<< cf.nrefborder << "\npborder = " << cf.pborder << "*)\n";
		//		sout << "Eig" << pborder << nrefborder << porder << nref << " = Eigenvalues[ek" << pborder << nrefborder << porder << nref << "]\n";
		LOGPZ_INFO(loggerres,sout.str());
	}
#endif

	// print the stiffness matrix corresponding to the condensed problem
	TPZElementMatrix ek(comp,TPZElementMatrix::EK),ef(comp,TPZElementMatrix::EF);
	comp->ElementVec()[subelindex]->CalcStiff(ek,ef);
	PlotEigenvectors(subc, listconfig, cf, analysis.Solution() , origin);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << cf.ek() << " = ";
		ek.fMat.Print(" ",sout,EMathematicaInput);
		sout << "(*\nnref = "<< cf.nref << "\nporder = " << cf.porder << "\nnrefborder = "<< cf.nrefborder << "\npborder = " << cf.pborder << "*)\n";
//		sout << "Eig" << pborder << nrefborder << porder << nref << " = Eigenvalues[ek" << pborder << nrefborder << porder << nref << "]\n";
		LOGPZ_INFO(loggerres,sout.str());
	}
#endif
	TPZCompEl *substruct = comp->ElementVec()[subelindex];
	comp->ElementVec()[subelindex] = 0;
	CriarMateriaisSteklov(*comp, cf);
	TPZFStructMatrix fullstruct(comp);
	TPZFMatrix rhs;
	TPZMatrix *Mass = fullstruct.CreateAssemble(rhs);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << cf.em() << " = ";
		Mass->Print(" ",sout,EMathematicaInput);
		LOGPZ_INFO(loggerres,sout.str());
	}
#endif
	comp->ElementVec()[subelindex] = substruct;
	delete Mass;

#ifdef EXTCOMP
	TPZFMatrix L2 = ProjecaoL2(*comp, subelindex, cf);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << cf.emEXT() << " = ";
		L2.Print(" ",sout,EMathematicaInput);
		sout << cf.Eigenvalues();
		LOGPZ_INFO(loggerres,sout.str());
	}
#endif
#endif
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "(*\nnref = "<< cf.nref << "\nporder = " << cf.porder << "\nnrefborder = "<< cf.nrefborder << "\npborder = " << cf.pborder << "*)\n";
		sout << cf.Eigenvalues();
		LOGPZ_INFO(loggerres,sout.str());
	}
#endif
	geo->ResetReference();
	delete comp;
		delete geo;
}


void RefineSquare(TPZCompMesh &cmesh, int nref, int porderstart)
{
//	TPZGeoMesh *gmesh = cmesh.Reference();
	int iref;
	for(iref = 0; iref<nref+1; iref++)
	{
		// find quadrilateral elements which contain node 1
		int nel = cmesh.NElements();
		std::list<TPZCompEl *> right;
		int iel;
		for(iel=0; iel<nel; iel++)
		{
			TPZCompEl *cel = cmesh.ElementVec()[iel];
			if(!cel) continue;
			TPZGeoEl *gel = cel->Reference();
			if(!gel) continue;
			if(gel->NodeIndex(0) == 1 || gel->NodeIndex(1) == 1)
			{
				right.push_back(cel);
			}
		}
		std::list<TPZCompEl *>::iterator it;
		for(it = right.begin(); it!= right.end(); it++)
		{
			TPZCompEl *cel = *it;
			TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> (cel);
			if(intel)
			{
				if(iref < nref)
				{
					intel->PRefine(porderstart+iref);
					int index = cel->Index();
					TPZManVector<int> subel;
					cel->Divide(index,subel);
				}
				else
				{
					intel->PRefine(1);
				}
			}
		}
	}
	cmesh.ExpandSolution();
}

TPZGeoMesh *LShape(TPZGeoMesh &model, TConfig &cf)
{
	TPZVec<REAL> x0(2,0.),x1(2,1.);
	TPZVec<int> nx(2,2);
	TPZGenGrid gen(nx,x0,x1,1,0.);
	TPZGeoMesh *geo = new TPZGeoMesh(model);
	gen.Read(*geo);
	gen.SetBC(geo,0,cf.materialborder1);
	gen.SetBC(geo,1,cf.materialborder2);
	gen.SetBC(geo,2,cf.materialborder2);
	gen.SetBC(geo,3,cf.materialborder1);
	TPZGeoEl *gel =  geo->ElementVec()[0];
	geo->DeleteElement(gel, 0);
	geo->DeleteElement(geo->ElementVec()[4], 4);
	geo->DeleteElement(geo->ElementVec()[10], 10);
/*
	TPZVec<REAL> one(3,0.), two(3,0.);
	one[0] = 0.5;
	two[0] = 0.5;
	two[1] = 0.5;
	gen.SetBC(geo,two,one,cf.materialcondense);
 */
	return geo;

}


void CriarMateriaisSteklov(TPZCompMesh &mesh, TConfig &cf)
{
	TPZMat2dLin *mat2dlin = new TPZMat2dLin(cf.materialid);
	TPZFMatrix xk(1,1,1.),xc(1,1,0.),xf(1,1,0.);
	mat2dlin->SetMaterial(xk, xc, xf);
	TPZAutoPointer<TPZMaterial> mat(mat2dlin);
	mesh.InsertMaterialObject(mat);
	TPZFMatrix val1(1,1,1.),val2(1,1,0.);
	TPZBndCond *bc = new TPZBndCond(mat,cf.materialborder1,2,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
	bc = new TPZBndCond(mat,cf.materialborder2,2,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
	bc = new TPZBndCond(mat,cf.materialborder3,2,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
	bc = new TPZBndCond(mat,cf.materialcorner,0,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
	bc = new TPZBndCond(mat,cf.materialcondense,1,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
}

void CriarMateriaisFlux(TPZCompMesh &mesh, TConfig &cf)
{
	TPZMat2dLin *mat2dlin = new TPZMat2dLin(cf.materialid);
	TPZFMatrix xk(1,1,1.),xc(1,1,0.),xf(1,1,0.);
	mat2dlin->SetMaterial(xk, xc, xf);
	TPZAutoPointer<TPZMaterial> mat(mat2dlin);
	mesh.InsertMaterialObject(mat);
	TPZFMatrix val1(1,1,1.),val2(1,1,0.);
	TPZBndCond *bc = new TPZBndCond(mat,cf.materialborder1,1,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
	val2(0,0) = 1.;
	bc = new TPZBndCond(mat,cf.materialborder2,1,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
	val1(0,0) = 1.;
	bc = new TPZBndCond(mat,cf.materialborder3,2,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
	val2(0,0) = 0.;
	bc = new TPZBndCond(mat,cf.materialcorner,1,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
	bc = new TPZBndCond(mat,cf.materialcondense,0,val1,val2);
	{
		TPZAutoPointer<TPZMaterial> bcauto(bc);
		mesh.InsertMaterialObject(bcauto);
	}
}

TPZFMatrix ProjecaoL2(TPZCompMesh &mesh, int subelindex, TConfig &cf)
{
	TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(mesh.ElementVec()[subelindex]);
	TPZMat2dLin *mat2dlin = dynamic_cast<TPZMat2dLin *>(sub->FindMaterial(cf.materialid).operator->());
	TPZFMatrix xk(1,1,0.),xc(1,1,1.),xf(1,1,0.);
	mat2dlin->SetMaterial(xk, xc, xf);
	TPZSkylineStructMatrix skylstr(sub);
	TPZFMatrix rhs;
	TPZMatrix *skyl = skylstr.CreateAssemble(rhs);
	int mesheq = mesh.NEquations();
	int condeseq = ((TPZCompMesh *) sub)->NEquations();
	TPZFMatrix sol(mesheq,1,0.);
	TPZFMatrix condens(condeseq ,mesheq);
	int ieq;
	for(ieq=0; ieq<mesheq; ieq++)
	{
		sol.Zero();
		sol(ieq,0) = 1.;
		mesh.LoadSolution(sol);
		TPZFMatrix subsol;
		subsol=((TPZCompMesh *)sub)->Solution();
		int j;
		for(j=0; j<condeseq; j++)
		{
			condens(j,ieq) = subsol(j,0);
		}
	}
	TPZFMatrix condensTR = condens;
	condensTR.Transpose();
	TPZFMatrix A,B;
	skyl->Multiply(condens, A);
	condensTR.Multiply(A, B);
	return B;

}


void TransferirElementos(TPZCompMesh &mesh, int materialid,int  &elindex)
{
	TPZSubCompMesh *submesh = new TPZSubCompMesh(mesh,elindex);
	int i, nelem;
	nelem = mesh.NElements();
	for(i=0; i<nelem; i++)
	{
		TPZCompEl *cel = mesh.ElementVec()[i];
		TPZAutoPointer<TPZMaterial> mat = cel->Material();
		if(mat && mat->Id() == materialid)
		{
			submesh->TransferElement(&mesh,i);
		}
	}
	submesh->MakeAllInternal();
	mesh.CleanUpUnconnectedNodes();
	ResequenceAsSubmesh(submesh);
}

void ResequenceAsSubmesh(TPZSubCompMesh *submesh)
{
	int nindices = submesh->NConnects();
	TPZManVector<int> permute(nindices);
	int ic;
	bool good=true;
	for(ic=0; ic<nindices; ic++)
	{
		int seqnum = submesh->Connect(ic).SequenceNumber();
		if(seqnum < nindices)
		{
			permute[seqnum] = ic;
		}
		else
		{
			good=false;
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "After transferring elements permutation " << permute;
		if(!good)
		{
			submesh->Mesh()->Print(sout);
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	submesh->Mesh()->Permute(permute);
}

void RefinarMalha(TPZCompMesh &comp,int materialid, int nref, int porder)
{
	int iref;
	for(iref = 0; iref<nref; iref++)
	{
		std::set<int> elref;
		int i, nel = comp.NElements();
		for(i=0; i<nel; i++)
		{
			TPZCompEl *cel = comp.ElementVec()[i];
			if(!cel) continue;
			TPZGeoEl *gel = cel->Reference();
			if(!gel) continue;
			if(gel->MaterialId() != materialid) continue;
			elref.insert(i);
		}
		std::set<int>::iterator it;
		for(it=elref.begin(); it!=elref.end(); it++)
		{
			TPZVec<int> subel;
			TPZCompEl *cel = comp.ElementVec()[*it];
			TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> (cel);
			if(intel) intel->PRefine(porder);
			cel->Divide(*it,subel);
		}
	}
	comp.ExpandSolution();
}

void CondensedEquations(TPZCompMesh &comp, int materialid , TPZVec<int> &origin)
{
	int nel = comp.ElementVec().NElements();
	int iel;
	std::set<int> exclude;
	for(iel=0; iel<nel; iel++)
	{
		TPZCompEl *cel = comp.ElementVec()[iel];
		if(!cel) continue;
		TPZGeoEl *gel = cel->Reference();
		if(!gel) continue;
		int matid = gel->MaterialId();
		if(matid == materialid)
		{
			int ncon;
			ncon = cel->NConnects();
			int ic;
			for(ic=0; ic<ncon; ic++)
			{
				TPZConnect &no = cel->Connect(ic);
				int seq = no.SequenceNumber();
				int firsteq = comp.Block().Position(seq);
				int size = comp.Block().Size(seq);
				int eq;
				for(eq=0; eq<size; eq++)
				{
					exclude.insert(firsteq+eq);
				}
			}
		}
	}
	int excludesize = exclude.size();
	origin.Resize(comp.NEquations()-excludesize);
	int count = 0,i;
	for(i=0; i<origin.NElements(); i++)
	{
		while(exclude.count(count)) count++;
		origin[i] = count++;
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Condensed equations " << origin;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

void PlotEigenvectors(TPZSubCompMesh *subc, TListConfig &listconfig, TConfig &cf, TPZFMatrix &solution , TPZVec<int> &condense)
{
	std::string eigenvectorfile = "Eigenvectors/" + cf.ExportFileName(listconfig.fName);
	std::ifstream sin(eigenvectorfile.c_str());
	if(sin)
	{
		TPZAnalysis *an = subc->GetAnalysis();
		int nvec = condense.NElements();
		int neq = solution.Rows();
		// Posprocessamento
		TPZVec<std::string> scalnames(1);
		scalnames[0] = "state";    //nome das vari�eis que se quer p�-processar
		TPZVec<std::string> vecnames(0);
		an->DefineGraphMesh(2, scalnames, vecnames, cf.EigenvectorFileName(listconfig.fName));
		int ivec,maxivec(nvec);
		if(maxivec > 40) maxivec = 40;
		for(ivec=0; ivec<maxivec; ivec++)
		{
			TPZVec<REAL> eigenvec(nvec,0.);
			TPZFMatrix eigenvecExpand(neq,1,0.);
			int i;
			for(i=0; i<nvec; i++)
			{
				sin >> eigenvec[i];
			}
			if(!sin)
			{
				return;
			}
			for(i=0; i<nvec; i++)
			{
				eigenvecExpand(condense[i],0.) = eigenvec[i];
			}
#ifdef LOG4CXX
			{
				std::stringstream sout;
				solution.Print("Finite element solution",sout);
				eigenvecExpand.Print("Eigenvector ",sout);
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			subc->Mesh()->LoadSolution(eigenvecExpand);
			//Define-se para a an�ise as vari�eis a p�-processar
			//Executa os c�culos para gera�o dos resultados de p�-processamento
			an->GraphMesh(2)->SetResolution(3);
			an->GraphMesh(2)->DrawSolution(ivec,ivec*1.);
		}
	}
}
