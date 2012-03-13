/*
 *  untitled.cpp
 *  PZ
 *
 *  Created by Denise de Siqueira on 6/9/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif



#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"

#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzquadrilateral.h"
#include "pzgnode.h"
#include "tpzarc3d.h"
#include "TPZGeoElement.h"
#include "tpzcompmeshreferred.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"
#include "pzinterpolationspace.h"
#include "TPZInterfaceEl.h"
#include "pzfstrmatrix.h"
#include "pzelasAXImat.h" 
#include "TPZRefLinear.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzgeopoint.h"
#include "pzpoisson3d.h"
#include "pzmaterialcoupling.h"
#include "pzelchdiv.h"
#include "pzsubcmesh.h"
#include "TPZRefPatternDataBase.h"
#include "TPZFrontStructMatrix.h"

#include "pzfunction.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("Acoplamento.main"));

#endif

using namespace std;

const int quadmat1=  1; // Parte inferior do quadrado
const int quadmat2 =  2; // Parte superior do quadrado

const int matInterface = 4;
const int quadmat3=3;// Material de interface
//contorno em omega1
const int mat1BC1 = -1; //Contorno 
const int mat1BC2= -2;
const int mat1BC3= -3;

//contorno em omega2

const int mat2BC1 = -4; //Contorno 
const int mat2BC2= -5;
const int mat2BC3= -6;



const int dirichlet = 0;
const int neumann = 1;
const int mista = 2;
const REAL Pi=4.*atan(1.);




void SolveLU ( TPZAnalysis &an );

void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file)
{
		file.clear();
		int nelements = gmesh->NElements();
		
		std::stringstream node, connectivity, type;
		
		//Header
		file << "# vtk DataFile Version 3.0" << std::endl;
		file << "TPZGeoMesh VTK Visualization" << std::endl;
		file << "ASCII" << std::endl << std::endl;
		
		file << "DATASET UNSTRUCTURED_GRID" << std::endl;
		file << "POINTS ";
		
		int actualNode = -1, size = 0, nVALIDelements = 0;
		
		for(int el = 0; el < nelements; el++)
		{		
				if(gmesh->ElementVec()[el]->Type() == EPoint)//Exclude Lines and Arc3D
				{
						continue;
				}
				if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
				{
						continue;
				}
				if(gmesh->ElementVec()[el]->HasSubElement())
				{
						continue;
				}
				
				int elNnodes = gmesh->ElementVec()[el]->NNodes();
				size += (1+elNnodes);
				connectivity << elNnodes;
				
				for(int t = 0; t < elNnodes; t++)
				{
						for(int c = 0; c < 3; c++)
						{
								double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
								node << coord << " ";
						}			
						node << std::endl;
						
						actualNode++;
						connectivity << " " << actualNode;
				}
				connectivity << std::endl;
				
				int elType = -1;
				switch (gmesh->ElementVec()[el]->Type())
				{
						case (ETriangle):
						{
								elType = 5;
								break;				
						}
						case (EQuadrilateral ):
						{
								elType = 9;
								break;				
						}
						case (ETetraedro):
						{
								elType = 10;
								break;				
						}
						case (EPiramide):
						{
								elType = 14;
								break;				
						}
						case (EPrisma):
						{
								elType = 13;
								break;				
						}
						case (ECube):
						{
								elType = 12;
								break;				
						}
						default:
						{
								//ElementType NOT Found!!!
								DebugStop();
								break;	
						}
				}

				type << elType << std::endl;
				nVALIDelements++;
		}
		node << std::endl;
		actualNode++;
		file << actualNode << " float" << std::endl << node.str();
		
		file << "CELLS " << nVALIDelements << " ";
		
		file << size << std::endl;
		file << connectivity.str() << std::endl;
		
		file << "CELL_TYPES " << nVALIDelements << std::endl;
		file << type.str();
		
		file.close();
}

void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file)
{
    refp->PrintVTK(file);
}

/**
 *Criar malha geometrica
 *h: Numero de refinamento uniforme
 */
TPZGeoMesh * MalhaGeoGen(int h, REAL &anglo);

TPZCompMesh *MalhaCompGen(TPZGeoMesh * gMesh, int p);
TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
void PrintInterface(TPZCompMesh *malha);
void SaddlePermute(TPZCompMesh * cmesh);
TPZGeoMesh * MalhaGeoT(const int h);
/**
 *Matrix de rotacao
 */
TPZFMatrix MatrixR(REAL ang);

void Forcing1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp) {
	/*	double x = pt[0];
		double y = pt[1];
		disp[0]= -2.*(1.-x*x) -2.*(1.-y*y);
		*/
		double x = pt[0];
		double y = pt[1];
		disp[0]= -2.*pow(Pi,2.)*sin(Pi*x)*sin(Pi*y);
		return;
}

void SolExata(TPZVec<REAL> &pt, TPZVec<REAL> &p, TPZFMatrix &flux ) {
		double x = pt[0];
		double y = pt[1];
		TPZVec<REAL> disp;
    p[0]= sin(Pi*x)*sin(Pi*y);
		flux(0,0)= -Pi*cos(Pi*x)*sin(Pi*y);
		flux(1,0)= - Pi*cos(Pi*y)*sin(Pi*x);
		flux(2,0)=2.*pow(Pi,2.)*sin(Pi*x)*sin(Pi*y);//coloco o divergente aq para testar
		
		return;		
		
		
}
void CC1(TPZVec<REAL> &pt, TPZVec<REAL> &func){//neumann
		func[0]=0;
		
		
}

int main()
{
		
		std::ofstream erro("CoupleSemEnrriqP2h5.txt");
//		REAL PI=4.*atan(1);
	//REAL a = PI/2.;
		
		for(int jorder=2;jorder<3;jorder++){
				erro<<"order "<<jorder<<std::endl;
		for(int nref=5;nref<6;nref++){
				
				erro<< "num ref "<<nref<<std::endl;
	// TPZGeoMesh *gmesh=MalhaGeoGen(nref, a);
				TPZGeoMesh *gmesh=MalhaGeoT(nref);
	//	std::ofstream file("MalhaAcoplado.vtk");
		//PrintGMeshVTK( gmesh, file);
		
		TPZCompMesh *comp=MalhaCompGen(gmesh, jorder);
		
			//PrintInterface(comp);
		TPZAnalysis analysis(comp);//resolve o problema
		//aplicar permutacao para os elementos de Hdiv
		
		TPZAdmChunkVector<TPZCompEl *> elvec = comp->ElementVec();
		
		int iel;
		int nel = elvec.NElements();
		for (iel=0; iel<nel; iel++) {
				TPZCompEl *cel = comp->ElementVec()[iel];
				if(!cel) continue;
				TPZGeoEl *gel = cel->Reference();
				if(gel->MaterialId() ==quadmat1)
				{
						SaddlePermute(comp);
				}
		}
		
		SolveLU(analysis);
		//std::ofstream CoupleSol("CoupleSolution.txt");
		//analysis.Solution().Print("CoupleSolution",CoupleSol);
		//std::ofstream CoupleWithCon("SolutionWithConects.txt");
		//analysis.Print( "SolutionWithConects" ,  CoupleWithCon);
		
		//calculo do erro
		
				TPZVec<REAL> calcErro;
				analysis.SetExact(*SolExata);
		    analysis.PostProcess(calcErro,erro);
		
		//4. visualizacao grafica usando vtk
		//		TPZVec<std::string> scalnames(1), vecnames(1);
				
		 
		// scalnames[0] = "Pressure";
				//scalnames[1] = "PressureOmega2";
		//		scalnames[1] = "ExactPressure";
			//	scalnames[2] = "Solution";
		 		 
		// vecnames[0] = "Flux";
				//vecnames[1] = "Derivate";
		//		vecnames[1] = "ExactFlux";
		 		 		 
		 
		// std::string plotfile("SolCouplig.vtk");
		// const int dim = 2;
		/* int div = 2;
				
				char buf[256] ;
				sprintf(buf,"SolCoupleSemEnrriq_jorder%d_nref%d.vtk",jorder,nref); 
				analysis.DefineGraphMesh(dim,scalnames,vecnames,buf);
		// analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
		 analysis.PostProcess(div);
		 
		 */
		}
		
		}
				
		return 0;
}

TPZGeoMesh *MalhaGeoGen(int h, REAL &anglo)
{
		int Qnodes = 9;
		
		TPZGeoMesh * gMesh = new TPZGeoMesh;
		gMesh->SetMaxNodeId(Qnodes-1);
		
			
		gMesh->NodeVec().Resize(Qnodes);
		
		TPZVec <int> TopolQuad(4);
		TPZVec <int> TopolLine(2);
		TPZVec <int> TopolPoint(1);
		
		TPZVec<TPZGeoNode> Node(Qnodes);
		
		//indice dos nos
		TPZFMatrix mA(2,6,0.), mR(2,2), mRA(2,9);
	
	/*	
		mA.Put(0,0,0.); mA.Put(1,0,0.);//(0,0)      
		mA.Put(0,1,0.5); mA.Put(1,1,0.);//(0.5,0)
		mA.Put(0,2,1.); mA.Put(1,2,0.);//(1,0)
		//mA.Put(0,3,0.); mA.Put(1,3,0.5);//(0,0.5)
		//mA.Put(0,4,0.5); mA.Put(1,4,0.5);//(0.5,0.5)
	//	mA.Put(0,5,1.); mA.Put(1,5,0.5);//(1,0.5)
		mA.Put(0,3,0.); mA.Put(1,3,1.);//(0,1)
		mA.Put(0,4,0.5); mA.Put(1,4,1);//(0.5,1)
		mA.Put(0,5,1.); mA.Put(1,5,1.);//(1,1)
		
		mR = MatrixR(anglo);
		mR.Multiply(mA, mRA);
		
		int id;
		id = 0;
		for (int j=0; j<6;j++) {
				
				Node[id].SetNodeId(id);
				Node[id].SetCoord(0 , mRA.GetVal(0, j)+0.);//coord r
				Node[id].SetCoord(1 , mRA.GetVal(1, j));//coord z
				gMesh->NodeVec()[id] = Node[id];
				id++;
		}
		id--;
		*/
		//Criar ns
		const int nnode = 6;//AQUI
		//const int nelem = 2;
		const int dim = 2;//AQUI
		
		REAL co[nnode][dim] = {{0.,0.},{0.5,0.},{1.,0.},{1.,1.},{0.5,1.},{0.,1.}};//{{-1.,-1.},{0.,-1.},{1.,-1.},{1.,1.},{0.,1.},{-1.,1.}};////AQUI
		int indices[1][nnode];//como serao enumerados os nos
		
		
		for(int i = 0; i < nnode; i++)
		{  
				indices[0][i] = i;
		}
		
		
		int nod;
		TPZVec<REAL> coord(dim);
		for(nod=0; nod<nnode; nod++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				
				for(int d = 0; d < dim; d++)
				{
						coord[d] = co[nod][d];
				}
				gMesh->NodeVec()[nod].Initialize(nodind,coord,*gMesh);
		}
		
		///--------------------------
		
		//indice dos elementos
			
	//	int index;
		
		TPZVec<int> nodind1(4);
		nodind1[0] = 0;
		nodind1[1] = 1;
		nodind1[2] = 4;
		nodind1[3] = 5;
		
		int id=0;
		
		new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,nodind1, quadmat1,*gMesh);
		id ++;
		
		TPZVec<int> nodind2(4);
		nodind2[0] = 1;
		nodind2[1] = 2;
		nodind2[2] = 3;
		nodind2[3] = 4;
		new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,nodind2, quadmat2,*gMesh);
		id ++;
		
		//-----------------------------------elemento de interface---------------------------
		TPZVec<int> nodind3(2);
		
		nodind3[0]=1;
		nodind3[1]=4;
		
		
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,nodind3,matInterface,*gMesh);
		id++;
		

		
		gMesh->AddInterfaceMaterial(quadmat1, quadmat2, quadmat3);
		gMesh->AddInterfaceMaterial(quadmat2, quadmat1, quadmat3);
		
		
		
		//--------------------------------elementos de contorno----------------------------------------	
				//omega1
		
				TopolLine[0] = 0;
				TopolLine[1] = 1;
				new TPZGeoElRefPattern<pzgeom:: TPZGeoLinear > (id,TopolLine,mat1BC1,*gMesh);//omega 1
				id++;
				
				TopolLine[0] = 4;
				TopolLine[1] = 5;
				new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat1BC2,*gMesh);
				id++;
		
				TopolLine[0] = 5;
				TopolLine[1] = 0;
				new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat1BC3,*gMesh);
				id++;
		
				// em omega 2
				TopolLine[0] = 1;
				TopolLine[1] = 2;
				new TPZGeoElRefPattern<pzgeom:: TPZGeoLinear > (id,TopolLine,mat2BC1,*gMesh);
				id++;
				
				
		
				TopolLine[0] = 2;
				TopolLine[1] = 3;
				new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat2BC2,*gMesh);
				id++;
				
				TopolLine[0] = 3;
				TopolLine[1] = 4;
				new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat2BC3,*gMesh);
				id++;
							
				gMesh->BuildConnectivity();
		
		//	Refinamento uniforme
		for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
				TPZVec<TPZGeoEl *> filhos;
				int n = gMesh->NElements();
				for(int i = 0; i < n; i++){
						TPZGeoEl * gel = gMesh->ElementVec()[i];
																			
						if(!gel->HasSubElement())
						{
								gel->Divide(filhos);
						}		
					
				}
				}
				
			//	ofstream arg("gmesh.txt");
			//	gMesh->Print(arg);
		
		
		return gMesh;	
}

TPZFMatrix MatrixR(REAL ang)
{
		TPZFMatrix r(2,2,0.);
		r(0,0) = cos(ang); r(0,1) = -sin(ang);
		r(1,0) = sin(ang);	r(1,1) =  cos(ang);
		
		return r;
}
TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
		if(!gel->Reference() && gel->NumInterfaces() == 0)
				return new TPZInterfaceElement(mesh,gel,index);
		
#ifdef LOG4CXX
		{
				std::stringstream sout;
				sout<<"elemento de interface "<<std::endl;
				LOGPZ_DEBUG(logger, sout.str().c_str());
		}
#endif 		 		
		
		return NULL;
}
TPZCompMesh *MalhaCompGen(TPZGeoMesh * gMesh, int porder)
{
		
		// Criar os tipos de materiais 
		TPZMatPoisson3d *mat1 = new TPZMatPoisson3d(quadmat1,2);//material para omega1
		TPZMatPoisson3d *mat2 = new TPZMatPoisson3d(quadmat2,2);//material para omega2
		TPZMaterialCoupling *matCouple = new TPZMaterialCoupling(quadmat3,1);//material de interface
		
		
		
		TPZAutoPointer<TPZMaterial> automat(mat1);
		TPZAutoPointer<TPZMaterial> automat2(mat2);
		
		TPZAutoPointer<TPZMaterial> automat3(matCouple);
		
		//relacao entre malha computacional e geometrica
		
		
		TPZCompEl::SetgOrder(porder);
		TPZCompMesh *comp = new TPZCompMesh(gMesh);
		
		int orderhdiv=comp->GetDefaultOrder();
		
			
		
		
		//inserindo os materiais na malha computacional	
		comp->InsertMaterialObject(automat);
		comp->InsertMaterialObject(automat2);
		comp->InsertMaterialObject(automat3);
		
		
		
		
		//setando forcing function para os dois materiais- Omega1 e 2
        TPZAutoPointer<TPZFunction> force = new TPZDummyFunction(Forcing1);
	  mat1->SetForcingFunction(force);
		mat2->SetForcingFunction(force);
		mat1->SetForcingFunctionExact(SolExata);
		mat2->SetForcingFunctionExact(SolExata);
		
		TPZFMatrix k(1,1,0.),f(1,1,0.);
		//k(0,0)=BIG;
		//	k(1,1)=BIG;
		
		//condicoes em omega1
		TPZMaterial *bnd =  automat->CreateBC(automat, mat1BC1, dirichlet,k,f);
		TPZMaterial *bnd1 =  automat->CreateBC(automat, mat1BC2, dirichlet,k,f);
		TPZMaterial *bnd2 =  automat->CreateBC(automat, mat1BC3, dirichlet,k,f);
		
			//condicoes em omega2
		TPZMaterial *bnd3 =  automat2->CreateBC(automat2, mat2BC1, dirichlet,k,f);
		//TPZFMatrix k2(1,1,0.),f2(1,1,1.);
		TPZMaterial *bnd4 =  automat2->CreateBC(automat2, mat2BC2, dirichlet,k,f);
		TPZMaterial *bnd5 =  automat2->CreateBC(automat2, mat2BC3, dirichlet,k,f);
		
		
		comp->InsertMaterialObject(bnd);
		comp->InsertMaterialObject(bnd1);
		comp->InsertMaterialObject(bnd2);
		comp->InsertMaterialObject(bnd3);
		comp->InsertMaterialObject(bnd4);
		comp->InsertMaterialObject(bnd5);
	//	comp->InsertMaterialObject(bnd6);
		
		//AQUI: criar espacos de aproximacao para Omega1
		
		comp->SetDefaultOrder(orderhdiv);
		comp->SetAllCreateFunctionsHDiv();
		
		//AQUI: AutoBuild(mat1)
		set<int> SETmat1;
		SETmat1.insert(quadmat1);
		SETmat1.insert(mat1BC1);
		SETmat1.insert(mat1BC2);
		SETmat1.insert(mat1BC3);
		comp->AutoBuild(SETmat1);
		
		comp->AdjustBoundaryElements();
		comp->CleanUpUnconnectedNodes();
		comp->SetName("Malha Para HDiv");
		
		ofstream arg1("CmeshHdiv.txt");
		comp->Print(arg1);
		gMesh->ResetReference();
		//setar ordem em omega 2
		comp->SetDefaultOrder(porder);
		//AQUI: criar espacos de aproximacao para Omega2
		
		comp->SetAllCreateFunctionsContinuous();
		//comp->SetAllCreateFunctionsDiscontinuous();
		//AQUI: AutoBuild(mat2)
		set<int> SETmat2;
		SETmat2.insert(quadmat2);
		SETmat2.insert(mat2BC1);
		SETmat2.insert(mat2BC2);
		SETmat2.insert(mat2BC3);
		SETmat2.insert(matInterface);//coloquei a interface aqui
		comp->AutoBuild(SETmat2);
		comp->AdjustBoundaryElements();
		comp->CleanUpUnconnectedNodes();
		
		comp->SetName("Malha Para H1");
		ofstream arg("CmeshPosH1.txt");
		comp->Print(arg);
		comp->LoadReferences(); 
		
		
		//AQUI: AutoBuild(matinterface)
	//	set<int> SETmat3;
	//	SETmat3.insert(matInterface);
	//	comp->AutoBuild(SETmat3);
	//	comp->AdjustBoundaryElements();
	//	comp->CleanUpUnconnectedNodes();
		
		//AQUI: Criar elemento de interface
		for(int el = 0; el < comp->ElementVec().NElements(); el++)
		{
				TPZCompEl * compEl = comp->ElementVec()[el];
				if(!compEl) continue;
				int matId = compEl->Reference()->MaterialId();
				if(matId== matInterface/*(matId== quadmat1)||(matId== quadmat2)||(matId== matInterface)*/)
				{
						compEl->Reference()->ResetReference();
				}
				
		}
		
		for(int el = 0; el < comp->ElementVec().NElements(); el++)
		{
				TPZCompEl * compEl = comp->ElementVec()[el];
				if(!compEl) continue;
				int index = compEl ->Index();
				if((compEl) && (compEl->Dimension() == 2) && (compEl->Reference()->MaterialId() == quadmat2))
				{
						TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(comp->ElementVec()[index]);
						int nsides = compEl->Reference()->NSides();
						TPZStack<TPZCompElSide> neighbours;
						
						int isini;//inicial 
						if (nsides==7) isini=3;
						else isini=4;
						
						for (int is=isini; is<nsides-1; is++) {
								neighbours.Resize(0);
								TPZCompElSide celside(compEl,is);
								celside.EqualLevelElementList(neighbours, 1, 1);
								int nneig = neighbours.NElements();
								if (nneig>1) DebugStop();
								if (nneig==1) {
										int idneig = neighbours[0].Reference().Element()->MaterialId();
										if (idneig==quadmat1) {
												/*TPZInterfaceElement * FaceElem =*/ InterpEl->CreateInterface(is,true);
												
										}
								}
						} 
				}
		}
		
		ofstream arg3("CmeshPosInterface.txt");
		comp->Print(arg3);
		
		return comp;
		
}
void SolveLU ( TPZAnalysis &an ){
		// Com matriz mal condicionada a solução é poluída com resíduos,
		// tanto com LU quanto Choleski. Isso resulta em não simetrias.
		TPZCompMesh *malha = an.Mesh();
	//	TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
		TPZFStructMatrix mat( malha );
		//	TPZSpStructMatrix mat( malha );
		TPZStepSolver solv;
		
		solv.SetDirect ( ELU );
		//		solv.SetDirect(ECholesky);
		
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
		// cout << "Banda = " << malha->BandWidth() << endl;
}

void PrintInterface(TPZCompMesh *malha)
{
		int el;
		const int nelem = malha->NElements();
		
		TPZInterfaceElement *face;
		
		///Loop sobre elementos computacionais 
		for(el=0; el < nelem; el++) 
				{
						TPZCompEl * Cel = malha->ElementVec()[el];
						if(!Cel) continue;
						int matId = Cel->Material()->Id();
						
						face = dynamic_cast<TPZInterfaceElement *> (Cel);
						if(face){
								int index = Cel->Index();
								cout<<endl;
								cout << "Ok ... O elemento de Index "<< index <<" e matId = "<< matId<<" eh de  interface"<<endl;
								Cel->Reference()->Print();
								cout <<endl;
								}
					}if(!face) cout << "nao tem interface"<< endl;
}

void SaddlePermute(TPZCompMesh * cmesh){
		
#ifdef LOG4CXX
		{
				std::stringstream sout;
				sout<< "Implementando permutacao para problemas de ponto de sela"<< std::endl;
				LOGPZ_DEBUG(logger, sout.str().c_str());
		}
#endif
		TPZVec<int> permute;
		int numinternalconnects = cmesh->NIndependentConnects();
  	permute.Resize(numinternalconnects,0);
		
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (cmesh);
		if(submesh)
		{
				int nexternal = submesh->NConnects();
				numinternalconnects -= nexternal;
		}
		//	else {
		//		DebugStop();
		//	}
		
		int jperm=0;
		int nel=cmesh->ElementVec().NElements();
		for (int jel=0; jel<nel; jel++) {
				
				for (int ip=0; ip<permute.NElements(); ip++) {
						permute[ip]=ip;
				}
				
				TPZCompEl *elvec=cmesh->ElementVec()[jel];
				
				
				//	int idtroca=0;
				int eqmax=0;
				if(!elvec)continue;
				int ncon=elvec->NConnects();
								
				
				//	if(ncon==1) continue;
				int eqpress=elvec->Connect(ncon-1).SequenceNumber();
				for (int icon=0; icon< ncon-1; icon++) {
						TPZConnect &coel=elvec->Connect(icon);
						int eqflux=coel.SequenceNumber();
						if (eqflux >= numinternalconnects) {
								continue;
						}
						eqmax = max(eqmax,eqflux);
				}
				
				
				if(eqpress<eqmax) {
						
						permute[eqpress]=eqmax;
						
				}
				
				
				for ( jperm = eqpress+1; jperm<=eqmax; jperm++) {
						permute[jperm]=jperm-1;
						
				}
				/*
				 #ifdef LOG4CXX
				 {
				 std::stringstream sout;
				 sout << "vetor SaddlePermute  do elemento - "<<jel<< " - " <<permute;
				 LOGPZ_DEBUG(logger, sout.str().c_str());
				 }
				 #endif
				 */
				cmesh->Permute(permute);
				
		}		
		
		
}
TPZGeoMesh * MalhaGeoT(const int h){//malha triangulo
		
		TPZGeoMesh *gmesh = new TPZGeoMesh();
		
		//Criar ns
		
		const int nnode = 6;//AQUI
		const int dim = 2;//AQUI
		
		REAL co[nnode][dim] = {{0.,0.},{0.5,0.},{1.,0.},{1.,1.},{0.5,1.},{0.,1.}};//{{-1.,-1.},{0.,-1.},{1.,-1.},{1.,1.},{0.,1.},{-1.,1.}};////AQUI
		int indices[1][nnode];//como serao enumerados os nos
		
		
		for(int i = 0; i < nnode; i++)
		{  
				indices[0][i] = i;
		}
				
		
		int nod;
		TPZVec<REAL> coord(dim);
		for(nod=0; nod<nnode; nod++) {
				int nodind = gmesh->NodeVec().AllocateNewElement();
				
				for(int d = 0; d < dim; d++)
				{
						coord[d] = co[nod][d];
				}
				gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
		}
		//Criacao de elementos
		
		
		//indice dos elementos
		
		//int index;
		
		TPZVec<int> nodind1(3);
		nodind1[0] = 0;
		nodind1[1] = 1;
		nodind1[2] = 5;
				
		int id=0;
		
		new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,nodind1, quadmat1,*gmesh);
		id ++;
		
		TPZVec<int> nodind2(3);
		nodind2[0] = 4;
		nodind2[1] = 5;
		nodind2[2] = 1;
				new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,nodind2, quadmat1,*gmesh);
		id ++;
		
		TPZVec<int> nodind3(3);
		nodind3[0] = 1;
		nodind3[1] = 2;
		nodind3[2] = 4;
		new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,nodind3, quadmat2,*gmesh);
		id ++;
		
		TPZVec<int> nodind4(3);
		nodind4[0] = 3;
		nodind4[1] = 4;
		nodind4[2] = 2;
		new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,nodind4, quadmat2,*gmesh);
		id ++;
		
		//-----------------------------------elemento de interface---------------------------
		TPZVec<int> nodind5(2);
		
		nodind5[0]=1;
		nodind5[1]=4;
		
		
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,nodind5,matInterface,*gmesh);
		id++;
		
		
		
		gmesh->AddInterfaceMaterial(quadmat1, quadmat2, quadmat3);
		gmesh->AddInterfaceMaterial(quadmat2, quadmat1, quadmat3);
		
		
		
		//--------------------------------elementos de contorno----------------------------------------	
		//omega1
		TPZVec<int> TopolLine(2);
		TopolLine[0] = 0;
		TopolLine[1] = 1;
		new TPZGeoElRefPattern<pzgeom:: TPZGeoLinear > (id,TopolLine,mat1BC1,*gmesh);//omega 1
		id++;
		
		TopolLine[0] = 4;
		TopolLine[1] = 5;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat1BC2,*gmesh);
		id++;
		
		TopolLine[0] = 5;
		TopolLine[1] = 0;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat1BC3,*gmesh);
		id++;
		
		// em omega 2
		TopolLine[0] = 1;
		TopolLine[1] = 2;
		new TPZGeoElRefPattern<pzgeom:: TPZGeoLinear > (id,TopolLine,mat2BC1,*gmesh);
		id++;
		
		
		
		TopolLine[0] = 2;
		TopolLine[1] = 3;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat2BC2,*gmesh);
		id++;
		
		TopolLine[0] = 3;
		TopolLine[1] = 4;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat2BC3,*gmesh);
		id++;

		
		gmesh->BuildConnectivity();
		
				
	
		const std::string nameref;
		
		TPZAutoPointer<TPZRefPattern> ref;
		
		
		//	Refinamento uniforme
		for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
				TPZVec<TPZGeoEl *> filhos;
				int n = gmesh->NElements();
				for(int i = 0; i < n; i++){
						TPZGeoEl * gel = gmesh->ElementVec()[i];
					
						if(!gel->HasSubElement() )
						{
								gel->Divide(filhos);
						}	
				}}
		
		
#ifdef LOG4CXX
		{
				std::stringstream sout;
				gmesh->Print(sout);
				LOGPZ_DEBUG(logger, sout.str().c_str());
		}
#endif 		 		 
		return gmesh;
}
