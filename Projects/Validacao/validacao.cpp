/*
 * @file
 * @author Denise de Siqueira
 * @since 6/9/11.
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzanalysiserror.h"
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "pzmatrix.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "pzinterpolationspace.h"
#include "pzsubcmesh.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzbndcond.h"
#include "TPZMaterial.h"
#include "tpzquadrilateral.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include <algorithm>

#include "TPZRefPattern.h"


#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("HdivTestes.main"));

#endif

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
const REAL Pi=4.*atan(1.);
void SolveLU ( TPZAnalysis &an );
void NormEnergia(TPZCompMesh *malha, std::ofstream &erro);
REAL NormEnergiaH1(TPZCompMesh *malha, std::ofstream &erro);
TPZGeoMesh * MalhaGeoT(const int h);
TPZGeoMesh * MalhaGeoQ(const int h);
TPZGeoMesh * MalhaGeo(const int h);
TPZGeoMesh * MalhaGeoQ2(const int h);
TPZCompMeshReferred *CreateCompMesh2d(TPZGeoMesh &gmesh,int porder);
int64_t SubStructure(TPZCompMesh *cmesh, int materialid);
void SaddlePermute(TPZCompMesh * cmesh);
void Forcing1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
	double x = pt[0];
	double y = pt[1];
	disp[0]= 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);//2.*pow(Pi,2.)*cos(Pi*y)*sin(Pi*x);//(1.)*8.;//-2.*exp(x)*(1. + 4.*x + pow(x,2.))*(-1. + pow(y,2.));//(exp(x)*(-3. + pow(y,2.) + x*(-4. + x + (4. + x)*pow(y,2.))));//2.*(1.-x*x) +2.*(1.-y*y); //
	return;
}
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux ) {
	double x = pt[0];
	double y = pt[1];
	TPZVec<REAL> disp;
    p[0]= sin(Pi*x)*sin(Pi*y);//-sin(Pi*x)*cos(Pi*y);//(1-x*x)*(1-y*y)*exp(x);//(1-x*x)*(1-y*y);//Solucao
	flux(0,0)= (-1.)*Pi*cos(Pi*x)*sin(Pi*y);//Pi*cos(Pi*x)*cos(Pi*y);//2.*exp(x)*x*(1. - pow(y,2.)) - exp(x)*(1. - pow(x,2.))*(1. - pow(y,2.));//2*x*(1-y*y);//
	flux(1,0)=  (-1.)*Pi*cos(Pi*y)*sin(Pi*x);//Pi*sin(Pi*y)*sin(Pi*x);//2.*exp(x)*(1. - pow(x,2.))*y;//2*(1-x*x)*y; dy
	flux(2,0)=2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);//-2.*pow(Pi,2.)*sin(Pi*x)*cos(Pi*y);//coloco o divergetne aq para testar
	
	
}
void CC1(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
	//double x=pt[0];
	//double y=pt[1];
	f[0] = 0.;//2*(1-x*x);// 
	
}
void CC2(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
	//double x=pt[0];
	double y=pt[1];
	f[0] = Pi*cos(Pi*y);//0.;//2*(1-x*x);// 
	
}
void CC3(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
	//double x=pt[0];
	//double y=pt[1];
	f[0]=0.;//2.*exp(x)*(1. - pow(x,2.));	//0.;//	
}
void CC4(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
	//double x=pt[0];
	double y=pt[1];
	f[0]=-Pi*cos(Pi*y);//2.*exp(x)*(1. - pow(x,2.));	//0.;//	
}
TPZGeoMesh * MalhaGeo2(const int h);
void SolGraf(TPZCompMesh *malha, std::ofstream &GraficoSol);
TPZCompMesh *CompMeshPAdap(TPZGeoMesh &gmesh,int porder);

/** Resolver o problema do tipo 
 * -Laplac(u) = 0
 * du/dn = lambda u em todo contorno
 */

void PrintMesh(TPZCompMesh *cmesh)
{
    int nel = cmesh->NElements();
    int iel;
    for(iel=0; iel<nel; iel++){
        
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout<<"\n Elemento computacional " << iel <<std::endl;
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif

        //TPZInterpolationSpace *intel;
        TPZInterpolatedElement *intel;
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        TPZGeoEl *geo = cel->Reference();
        intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(intel)
        {
//            TPZMaterialData data;
//            intel->InitMaterialData(data);
//            TPZElementMatrix ek(cmesh,TPZElementMatrix::EK),ef(cmesh,TPZElementMatrix::EF);
//            intel->InitializeElementMatrix(ek, ef);
//            intel->CalcStiff(ek,ef);
            
            //int nshapel = intel->NShapeF();
            //int dim = intel->Reference()->Dimension();
            
            int ncon = intel->NConnects();
            for (int i=0; i<ncon; i++)
            {
                int mlado = i + /*pzgeom::TPZGeoQuad::NSides*/ geo->NSides()- ncon;
                int ordinterp = intel->EffectiveSideOrder(mlado);
                
                intel->SetSideOrder(mlado, ordinterp+i);
#ifdef LOG4CXX
                int indexcon = intel->ConnectIndex(i);
                int preforder = intel->PreferredSideOrder(mlado);
                int nshape = intel->NConnectShapeF(i,preforder);
                int nsidescon = intel->NSideConnects(mlado);
                int newsideorde = intel->EffectiveSideOrder(mlado);
                {
                    std::stringstream sout;
                    sout<<"\n Connect " << indexcon << " Connect index = " << indexcon << " numero shapes = " <<nshape<<std::endl;
                    sout<< "\n lado associado ao connect "<< indexcon << " = "<< mlado << " ordem do polinomio associado ao lado = " << preforder<<" ordem de interpolacao do lado = " << ordinterp <<" new ordem de interpolacao do lado = " << newsideorde << " numero lados conectados = " <<nsidescon<<std::endl;
                    LOGPZ_DEBUG(logger,sout.str())
                }
#endif

            }
            
            TPZVec<int> ord;
            intel->GetInterpolationOrder(ord);
         
            
#ifdef LOG4CXX
			//            TPZManVector<REAL> xi(dim,0.);
			//            TPZFNMatrix<100> phi(nshape,1,0.),dphi(2,nshape,0.);
			//            intel->Shape(xi,phi,dphi);
            {
                std::stringstream sout;
                sout<<"\n num elementos  " << ord.size();
                sout<<"\n ordem de interpolacao  "<<ord<<std::endl;
                sout<<"\n =======================================\n"<<std::endl;
               // ord.Print(sout);
//                phi.Print("Shape functions ",sout);
//                dphi.Print("Derivative shape functions ",sout);
//                ek.Print(sout);
//                ef.Print(sout);
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
        }
    }
}

using namespace std;
int main()
{
	
#ifdef LOG4CXX
	{
		InitializePZLOG();
		std::stringstream sout;
		sout<< "Validacoes do codigo Hdiv"<<endl;
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif
	std::ofstream erro("Caulotaxa.txt");
	//std::ofstream GraficoSol("SolGraf.txt");
	//	std::ofstream CalcSolExata("CalSolExata.txt");
	TPZVec<REAL> calcErro;
	for (int h=1; h<3; h++) {
		
		erro<<"refinamento "<<h <<std::endl;
		//	erro<< " Flux exato " << "\t "<<" Flux aprox "<<std::endl;//"P exata " << " \t " <<"P aprox " << "\t " << " Flux exato " << "\t "<<" Flux aprox "<<std::endl;
		for(int porder=6;porder<10;porder++){
			erro<<std::endl;
			erro<< "\n Ordem polinomial "<<porder<<std::endl;
			//1. Criacao da malha geom. e computacional
			TPZGeoMesh *gmesh = MalhaGeoT(h);
    std::ofstream file("MalhaGeom.vtk");
		PrintGMeshVTK( gmesh, file);
				TPZCompMeshReferred *cmesh=CreateCompMesh2d(*gmesh,porder);
				
				int nDofs=cmesh->NEquations();
				erro<< "\n NDofs "<<nDofs<<std::endl;
			
			//TPZCompMesh *cmesh = CompMeshPAdap(*gmesh2,porder);
			
			//cmesh->LoadReferences();//mapeia para a malha geometrica lo
			
			//TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
			
	
            
			//2. Resolve o problema
			TPZAnalysis analysis(cmesh);
			//SaddlePermute(cmesh);
			SolveLU ( analysis );
		//	ofstream file("Solutout");
    //        analysis.Solution().Print("solution", file);
			//Resolver o sistema linear
			 //TPZFStructMatrix str(cmesh);
			// analysis.SetStructuralMatrix(str);
			// TPZStepSolver step;
			// step.SetDirect(ELU);
			// analysis.SetSolver(step);
			// analysis.Run();
			 
			//Pos processamento
//			std::ofstream SolPoissonHdiv("Solucao.txt");
//			analysis.Solution().Print("Solucao",SolPoissonHdiv);
//			std::ofstream SolP("teste.txt");
//			analysis.Print( "SolTeste" ,  SolP);
			/*----
			 TPZVec< REAL > p(1),pto(3);
			 TPZVec<REAL> fluxo(3);
			 TPZManVector<REAL,100> fluxoAprox(3);
			 
			 
			 int nNodes=gmesh2->NodeVec().NElements();//NNodes();
			 TPZVec<REAL> VecSol(nNodes);
			 for (int jnode=0; jnode< nNodes; jnode++) {
			 
			 pto[0]=	gmesh2->NodeVec()[jnode].Coord(0);
			 pto[1]=	gmesh2->NodeVec()[jnode].Coord(1);
			 //	pto[2]=	gmesh2->NodeVec()[jnode].Coord(2);
			 SolExata(pto, p, fluxo);
			 VecSol[jnode]=p[0];
			 CalcSolExata<<"{ "<<pto<<","<<p[0]<< "},";
			 //	CalcSolExata<<p[0]<< ",";
			 }
			 
			 */				
			
			
			
			
			//3. Calcula o erro usando norma energia
			
			//	NormEnergiaH1(cmesh,erro);
			//		NormEnergia(cmesh,erro);
			
			/*
			 int nelem=cmesh->ElementVec().NElements();
			 for (int jel=0; jel<nelem; jel++) {
			 #ifdef LOG4CXX
			 {
			 std::stringstream sout;
			 sout<< "El 0"<<jel<< std::endl;
			 LOGPZ_DEBUG(logger,sout.str())
			 }
			 #endif
			 
			 TPZCompEl *elvec=cmesh->ElementVec()[jel];
			 if(!elvec)continue;
			 int ncon=elvec->NConnects();
			 if(ncon==1) continue;
			 for (int icon=0; icon< ncon; icon++) {
			 TPZConnect &coel=elvec->Connect(icon);
			 int eqflux=coel.SequenceNumber();
			 #ifdef LOG4CXX
			 {
			 std::stringstream sout;
			 sout<< " connect "<<icon<< " seq number "<< eqflux<<std::endl;
			 LOGPZ_DEBUG(logger,sout.str())
			 }
			 #endif
			 
			 }
			 
			 }
			 */			
			
			
			//-----
			//TPZVec<REAL> calcErro;
			analysis.SetExact(*SolExata);
		  analysis.PostProcess(calcErro,erro);
			
//			//4. visualizacao grafica usando vtk
//			 TPZVec<std::string> scalnames(2), vecnames(2);
//			 
//			// scalnames[0] = "Divergence";
//			 scalnames[0] = "Pressure";
//			 scalnames[1] = "ExactPressure";
//			 //	scalnames[2] = "ExactDiv";
//			 
//			 
//			 vecnames[0] = "ExactFlux";
//			 vecnames[1] = "Flux";
//			 //scalnames[2] = "Divergence";
//			 
//			 
//			 //vecnames[0] = "Derivate";
//			 
//			 
//			 std::string plotfile("GraficoSolution.vtk");
//			 const int dim = 2;
//			 int div = 2;
//			 analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//			 analysis.PostProcess(div);
//			 
			
		}}
	
	
	return 0;
}



TPZCompMesh *CompMeshPAdap(TPZGeoMesh &gmesh,int porder){
		
		//int NgeoEl=gmesh.NElements();
//		for(int iel=0;iel<NgeoEl;iel++){
//				TPZGeoEl *geoEl=gmesh.ElementVec()[iel];
//				int ElId=geoEl->MaterialId();
//				if (ElId==40 || ElId==42) {
//						TPZCompEl::SetgOrder(porder);
//				}
//				if (ElId==41 || ElId==43){
//						TPZCompEl::SetgOrder(porder+1);
//				}
//
//		
//		}

		
		TPZCompEl::SetgOrder(porder);
		TPZCompMesh *comp = new TPZCompMesh(&gmesh);
		
		// Criar e inserir os materiais na malha
		TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
		TPZMaterial * automat(mat);
		comp->InsertMaterialObject(automat);
		
	
    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1, 5);
		mat->SetForcingFunction(force1);
		TPZAutoPointer<TPZFunction<STATE> > exata1 = new TPZDummyFunction<STATE>(SolExata, 5);
		mat->SetForcingFunctionExact(exata1);
		
		comp->SetAllCreateFunctionsHDivPressure();
		
		// Ajuste da estrutura de dados computacional
		comp->AutoBuild();
		
		
		///Criar condicoes de contorno
		
		TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
		TPZFMatrix<STATE> val11(1,1,0.), val22(1,1,0.);
		TPZMaterial *bnd = automat->CreateBC (automat,-1,0,val1,val2);
		TPZMaterial *bnd2 = automat->CreateBC (automat,-2,0,val1,val2);
		TPZMaterial *bnd3 = automat->CreateBC (automat,-3,0,val1,val2);
		TPZMaterial *bnd4 = automat->CreateBC (automat,-4,0,val1,val2);
		TPZMaterial *bnd5 = automat->CreateBC (automat,-5,0,val1,val2);
		TPZMaterial *bnd6 = automat->CreateBC (automat,-6,0,val1,val2);
		TPZMaterial *bnd7 = automat->CreateBC (automat,-7,0,val1,val2);
		TPZMaterial *bnd8 = automat->CreateBC (automat,-8,0,val1,val2);
		
		///Inserir condicoes de contorno
		comp->InsertMaterialObject(bnd);
		comp->InsertMaterialObject(bnd2);
		comp->InsertMaterialObject(bnd3);
		comp->InsertMaterialObject(bnd4);	
		comp->InsertMaterialObject(bnd5);
		comp->InsertMaterialObject(bnd6);
		comp->InsertMaterialObject(bnd7);
		comp->InsertMaterialObject(bnd8);
		
  	//espaco de aproximacao
		comp->SetAllCreateFunctionsHDivPressure();
		
		// Ajuste da estrutura de dados computacional
		comp->AutoBuild();
		comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
		comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
		comp->SetName("Malha Computacional Original");
		
#ifdef LOG4CXX
		{
				std::stringstream sout;
				comp->Print(sout);
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		
		
    return comp;
		
}
TPZCompMeshReferred *CreateCompMesh2d(TPZGeoMesh &gmesh,int porder){
	TPZCompEl::SetgOrder(porder);
	TPZCompMeshReferred *comp = new TPZCompMeshReferred(&gmesh);
    int dim=2;
    comp->SetDimModel(dim);
	// Criar e inserir os materiais na malha
	TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,dim);
	TPZMaterial * automat(mat);
	comp->InsertMaterialObject(automat);
	
    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1, 5);
	mat->SetForcingFunction(force1);
	TPZAutoPointer<TPZFunction<STATE> > exata1 = new TPZDummyFunction<STATE>(SolExata, 5);
	mat->SetForcingFunctionExact(exata1);
	///Inserir condicoes de contorno
	
	TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
	TPZFMatrix<STATE> val11(1,1,0.), val22(1,1,0.);
	TPZMaterial *bnd = automat->CreateBC (automat,-1,0,val1,val2);//1
	TPZMaterial *bnd2 = automat->CreateBC (automat,-2,0,val1,val2);
	TPZMaterial *bnd3 = automat->CreateBC (automat,-3,0,val1,val2);//1
	TPZMaterial *bnd4 = automat->CreateBC (automat,-4,0,val1,val2);
	
//	TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(CC1, 5);
//    //TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(CC2, 5);
//	bnd->SetForcingFunction(fCC1);
//	bnd2->SetForcingFunction(fCC1);
//	bnd3->SetForcingFunction(fCC1);
//	bnd4->SetForcingFunction(fCC1);
	
    comp->SetAllCreateFunctionsHDivPressure();
	comp->InsertMaterialObject(bnd);
	comp->InsertMaterialObject(bnd2);
	comp->InsertMaterialObject(bnd3);
	comp->InsertMaterialObject(bnd4);	
	//espaco de aproximacao
	//comp->SetAllCreateFunctionsHDiv();
		
	//comp->SetAllCreateFunctionsContinuous();
	
	// Ajuste da estrutura de dados computacional
	comp->AutoBuild();
	comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
	comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
	comp->SetName("Malha Computacional Original");
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		comp->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	
	
    return comp;
	
}


TPZGeoMesh * MalhaGeoT(const int h){//malha triangulo
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();

    //Criar ns
	const int nnode = 4;//AQUI
//	const int nelem = 2;
	TPZGeoEl *elvec[2]; //nelem	
	const int dim = 2;//AQUI
    gmesh->SetDimension(dim);

	REAL co[nnode][dim] ={{0.,0.},{1.,0.},{1.,1.},{0.,1.}};//{{-1.,0.},{1.,0.},{1.,1.},{-1.,1.}};// {{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};//
	int indices[2][nnode];//como serao enumerados os nos
	
	//el 1
	indices[0][0] = 0;
	indices[0][1] = 1;
	indices[0][2] = 3;
	//el2
	indices[1][0] = 2;
	indices[1][1] = 3;
	indices[1][2] = 1;
	
	
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
	
	
	TPZVec<int64_t> nodind1(3);
	TPZVec<int64_t> nodind2(3);
	for(int i=0; i<3; i++){
		nodind1[i] = indices[0][i];
		nodind2[i] = indices[1][i];
	}
	
	int64_t index;
	elvec[0] = gmesh->CreateGeoElement(ETriangle,nodind1,1,index); //AQUI
	elvec[1] = gmesh->CreateGeoElement(ETriangle,nodind2,1,index); //AQUI
	
	
	gmesh->BuildConnectivity();
	
	//	como usar os padroes de refinamento do caju..nao posso rodar ainda pq tenho o Pz antigo..nao tenho este gRefDBase.
	//{
	//			TPZAutoPointer< TPZRefPattern > refExemplo = gRefDBase->FinfRefPattern("a_Quad_Rib_Side_4_5_6");
	//				if(refExemplo)
	//			{
	//						TPZGeoEl * elExemplo = gmesh->Elementvec()[10];
	//						elExemplo->SetRefPattern(refExemplo);
	//						elEmemplo->Divide(filhos);
	//				}
	
	//	}
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec[0],3,-1);// condicao de fronteira tipo -1: 
	TPZGeoElBC gbc2(elvec[0],5,-2);// condicao de fronteira tipo -2: 
	
	TPZGeoElBC gbc3(elvec[1],3,-3);// condicao de fronteira tipo -3: 
	TPZGeoElBC gbc4(elvec[1],5,-4);// condicao de fronteira tipo -4: 
	
	const std::string nameref;
	
	TPZAutoPointer<TPZRefPattern> ref;
	
	
	//	Refinamento uniforme
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
		}}
		//refinamento diferencialvel
		{
				
				TPZVec<TPZGeoEl *> filhos;
				int n = gmesh->NElements();
				
				
				for(int i = 0; i < n; i++){	
						TPZGeoEl * gel = gmesh->ElementVec()[i];
						if(!gel->HasSubElement() && gel->Dimension()==2 && i%2==0)
								
						{
								gel->Divide(filhos);
						}		
				}
		}
		
		//refinamento 1D--irei refinar tambem os elementos 1D
		
		{
				
				TPZVec<TPZGeoEl *> filhos;
				int n = gmesh->NElements();
				
				
				for(int i = 0; i < n; i++){	
						TPZGeoEl * gel = gmesh->ElementVec()[i];
						if (gel->Dimension()!=1) {
								continue;
						}
						TPZGeoElSide Elside=gel->Neighbour(2);
						TPZGeoEl *NeighEl=Elside.Element();
						if (NeighEl->HasSubElement()) {
								gel->Divide(filhos);
						}
						
						
						
				}
		}
		
		
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif 		 		 
	return gmesh;
}
TPZGeoMesh * MalhaGeo/*QUADRILATEROS*/ ( const int h )
{
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	REAL co[4][2] = {{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};//{{0.,0.},{1.,0.},{1.,1.},{0.,1.}};
	int indices[1][4] = {{0,1,2,3}};
	
	int nnode = 4;
	const int64_t nelem = 1;
	TPZGeoEl *elvec[1]; //nelem
	int nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}
	
	int64_t el;
	for ( el=0; el<nelem; el++ )
	{
		TPZVec<int64_t> nodind ( 4 );
		for ( nod=0; nod<4; nod++ ) nodind[nod]=indices[el][nod];
		int64_t index;
		elvec[el] = gmesh->CreateGeoElement ( EQuadrilateral,nodind,1,index );
	}
	
	gmesh->BuildConnectivity();
	
	TPZGeoElBC gbc1 ( elvec[0],4,-1 );// condicao de fronteira tipo -1: (x,y=0)
	TPZGeoElBC gbc2 ( elvec[0],5,-2 );// condicao de fronteira tipo -2: (x=1,y)
	TPZGeoElBC gbc3 ( elvec[0],6,-3 );// condicao de fronteira tipo -3: (x,y=1)
	TPZGeoElBC gbc4 ( elvec[0],7,-4 );// condicao de fronteira tipo -4: (x=0,y)
	
	///Refinamento uniforme
	for ( int ref = 0; ref < h; ref++ )
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int64_t n = gmesh->NElements();
		for ( int64_t i = 0; i < n; i++ )
		{
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			//if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}	
		}//for i
	}//ref
		//refinamento diferencialvel
//		{
//				
//				TPZVec<TPZGeoEl *> filhos;
//				int n = gmesh->NElements();
//				
//				
//				for(int i = 0; i < n; i++){	
//						TPZGeoEl * gel = gmesh->ElementVec()[i];
//						if(!gel->HasSubElement() && gel->Dimension()==2 && i%2==0)
//								
//						{
//								gel->Divide(filhos);
//						}		
//				}
//		}
//		
//		//refinamento 1D--irei refinar tambem os elementos 1D
//		
//		{
//				
//				TPZVec<TPZGeoEl *> filhos;
//				int n = gmesh->NElements();
//				
//				
//				for(int i = 0; i < n; i++){	
//						TPZGeoEl * gel = gmesh->ElementVec()[i];
//						if (gel->Dimension()!=1) {
//								continue;
//						}
//						TPZGeoElSide Elside=gel->Neighbour(2);
//						TPZGeoEl *NeighEl=Elside.Element();
//						if (NeighEl->HasSubElement()) {
//								gel->Divide(filhos);
//						}
//						
//						
//						
//				}
//		}
//		
		
	return gmesh;
}

TPZGeoMesh * MalhaGeoQ(const int h){//malha quadrilatera
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	//Criar nos
	const int64_t nnode = 4;//AQUI
	const int dim = 2;//AQUI
	
	REAL co[nnode][dim] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};//{{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};//
	int64_t indices[1][nnode];//={0,1,2,3};//como serao enumerados os nos
	
	
	
	for(int64_t i = 0; i < nnode; i++){
		indices[0][i] = i;
	}
	
	
	int64_t nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		int64_t nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos
	
	
	TPZVec<int64_t> nodind(4);
	for(int i=0; i<4; i++){
		nodind[i] = indices[0][i];
	}
	
	int64_t index;
	TPZGeoEl *elvec = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index); //AQUI
	
	//	gmesh->BuildConnectivity();
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec,4,-1);// condicao de fronteira tipo -1: (x,y=0)
	TPZGeoElBC gbc2(elvec,5,-2);// condicao de fronteira tipo -2: (x=1,y)
	TPZGeoElBC gbc3(elvec,6,-3);// condicao de fronteira tipo -3: (x,y=1)
	TPZGeoElBC gbc4(elvec,7,-4);// condicao de fronteira tipo -4: (x=0,y)
	gmesh->BuildConnectivity();
	
	const std::string nameref;
	
	TPZAutoPointer<TPZRefPattern> ref;
	//gmesh->RefPatternList(EQuadrilateral);
	
	//	Refinamento uniforme
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int64_t n = gmesh->NElements();
		for(int64_t i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			
			
		}
		
		
	}
		
		//refinamento diferencialvel
		{
				
				TPZVec<TPZGeoEl *> filhos;
				int64_t n = gmesh->NElements();
				
				
				for(int64_t i = 0; i < n; i++){	
						TPZGeoEl * gel = gmesh->ElementVec()[i];
						if(!gel->HasSubElement() && gel->Dimension()==2 && i%2==0)
								
						{
								gel->Divide(filhos);
						}		
				}
		}
		
		//refinamento 1D--irei refinar tambem os elementos 1D
		
		{
				
				TPZVec<TPZGeoEl *> filhos;
				int64_t n = gmesh->NElements();
				
				
				for(int64_t i = 0; i < n; i++){	
						TPZGeoEl * gel = gmesh->ElementVec()[i];
						if (gel->Dimension()!=1) {
								continue;
						}
						TPZGeoElSide Elside=gel->Neighbour(2);
						TPZGeoEl *NeighEl=Elside.Element();
						if (NeighEl->HasSubElement()) {
								gel->Divide(filhos);
						}
						
						
						
				}
		}
		
	
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif 		 		 
	return gmesh;
} 
TPZGeoMesh * MalhaGeo2(const int h){//malha quadrilatera com 2 elementos
	TPZGeoMesh *gmesh = new TPZGeoMesh();
		int nelem=4;
	TPZGeoEl *elvec[4]; //nelem
	//Criar ns
	const int nnode = 9;//AQUI
	const int dim = 2;//AQUI
	
		REAL co[nnode][dim] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.},{0.5,0},{1.,0.5},{0.5,1.},{0.,0.5},{0.5,0.5}};//{{-1.,-1},{1.,-1},{1.,1.},{-1.,1.},{0.,-1.},{0.,1.}};
	
	
		int nodindAll[4][nnode]={{0,4,8,7},{4,1,5,8},{8,5,2,6},{7,8,6,3}};//como serao enumerados os nos
	
	
	int nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		int64_t nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos
		int matId=10;
		int64_t index;
	for ( int el=0; el<nelem; el++ )
	{
		TPZVec<int64_t> nodind(4);
		nodind[0]=nodindAll[el][0];
		nodind[1]=nodindAll[el][1];
		nodind[2]=nodindAll[el][2];
		nodind[3]=nodindAll[el][3];
		
		elvec[el] = gmesh->CreateGeoElement (EQuadrilateral,nodind,matId,index );
			
		//	matId++;
			index++;

			
	}
		
		
		

	
	
	gmesh->BuildConnectivity();
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec[0],4,-1);
  TPZGeoElBC gbc2(elvec[0],7,-2);
	TPZGeoElBC gbc3(elvec[1],4,-3);
	TPZGeoElBC gbc4(elvec[1],5,-4);
	TPZGeoElBC gbc5(elvec[2],5,-5);
	TPZGeoElBC gbc6(elvec[2],6,-6);
	TPZGeoElBC gbc7(elvec[3],6,-7);
  TPZGeoElBC gbc8(elvec[3],7,-8);
	
	
	const std::string nameref;
	
	TPZAutoPointer<TPZRefPattern> ref;
	//gmesh->RefPatternList(ETriangle);
	
	//Refinamento uniforme
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			
			
		}
		
		
	}
	
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif 		 		 
	return gmesh;
	
}
void SolGraf(TPZCompMesh *malha, std::ofstream &GraficoSol){
	const int nelem = malha->NElements();
	
	//TPZFMatrix<REAL> sol(4,1,0.);
	TPZManVector<REAL,4> solP(1);
	TPZManVector<REAL,4> solF(3,0.);
	
	TPZFMatrix<REAL> axes(3,1,0.);
	TPZFMatrix<REAL> phi;
	TPZFMatrix<REAL> dphix;
	
	//TPZManVector<REAL> sol(1);
	TPZSolVec sol;
	sol[0].Resize(1);
	//TPZFMatrix<REAL> dsol(3,1,0.);
	TPZGradSolVec dsol;
	dsol[0].Redim(3,1);
	
	///Percorrer todos elementos 
	for(int el=0; el < nelem; el++){
		TPZCompEl * Cel = malha->ElementVec()[el];
		
		if(!Cel) continue;
		
		
		TPZGeoEl *gel = Cel->Reference();
		if(	malha->ElementVec()[el]->Dimension()== 2) {
			TPZIntPoints *intr = gel-> CreateSideIntegrationRule(gel->NSides() - 1,2);//side and order 
			
			for(int in=0; in < intr->NPoints(); in++)
			{
				REAL peso;
				TPZVec<REAL> pto(2),pto2(1);
				intr->Point(in, pto, peso);
				
				TPZFMatrix<REAL> jac(2,2),invjac(2,2),axes(3,3);
				REAL jacdet;
				gel->Jacobian(pto,jac,axes,jacdet,invjac);
				
				//	malha->LoadSolution(sol);
				
				
				TPZManVector< REAL,3 > xco(3);
                TPZManVector<STATE> p(1);
				TPZFMatrix<STATE> fluxo(3,0);
				TPZManVector<STATE,4> solF(3,0.);
				
				gel->X(pto,xco);
				SolExata(xco, p, fluxo);
				Cel->ComputeSolution(xco,sol,dsol,axes);
				GraficoSol<<"{ "<<xco[0]<< " , "<< xco[1]<< ","<<sol[0]<<"},"<<std::endl;
				
				
			}
			delete intr;
			
		}
		else continue;
		
		
	}
	
	
}
/*void NormEnergia(TPZCompMesh *malha, std::ofstream &erro){
 const int nelem = malha->NElements();
 REAL NormaP=0.,NormaF=0.,auxNormaP=0.,auxNormaF=0.;
 REAL NormaPexact=0.,NormaFexact=0.,auxNormaPexact=0.,auxNormaFexact=0.;
 
 //TPZFMatrix<REAL> sol(4,1,0.);
 TPZManVector<REAL,4> solP(1);
 
 TPZManVector<REAL,4> solF(3,0.);
 TPZManVector<REAL,3> fatErroF(3,0.);
 TPZFNMatrix<660> dsol(3,1,0.);
 TPZFMatrix<REAL> axes(3,1,0.);
 TPZFMatrix<REAL> phi;
 TPZFMatrix<REAL> dphix;
 ///Percorrer todos elementos 
 for(int el=0; el < nelem; el++){
 TPZCompEl * Cel = malha->ElementVec()[el];
 
 if(!Cel) continue;
 
 
 TPZGeoEl *gel = Cel->Reference();
 if(	malha->ElementVec()[el]->Dimension()== 2) {
 TPZIntPoints *intr = gel-> CreateSideIntegrationRule(gel->NSides() - 1,5);//side and order 
 
 for(int in=0; in < intr->NPoints(); in++)
 {
 REAL peso;
 TPZVec<REAL> pto(2),pto2(1);
 intr->Point(in, pto, peso);
 
 TPZFMatrix<REAL> jac(2,2),invjac(2,2),axes(3,3);
 REAL jacdet;
 gel->Jacobian(pto,jac,axes,jacdet,invjac);
 
 
 
 TPZVec< REAL > xco(3), p(1);
 TPZVec<REAL> fluxo(3);
 TPZManVector<REAL,4> solF(3,0.);
 
 gel->X(pto,xco);
 SolExata(xco, p, fluxo);
 
 Cel->Solution(xco, 10, solF); //quero primeiro a pressao, fluxo eh 10
 //Cel->Solution(xco, 11, solP);F
 auxNormaF=(fluxo[0])*(fluxo[0])+(fluxo[1])*(fluxo[1])+(fluxo[2])*(fluxo[2]);//(solF[0])*(solF[0])+(solF[1])*(solF[1])+(solF[2])*(solF[2]);
 auxNormaFexact=p[0]*p[0];//(fluxo[0])*(fluxo[0])+(fluxo[1])*(fluxo[1])+(fluxo[2])*(fluxo[2]);
 
 
 #ifdef LOG4CXX
 {
 std::stringstream sout;
 //		sout<<"pto "<<xco<< " pressao "<< p<< " pressao aprox "<< solP<<" erro "<< erroP<<std::endl;
 sout<<" fluxo "<< fluxo<< " fluxo aprox "<< solF<<" erro "<< fatErroF <<std::endl;
 LOGPZ_DEBUG(logger, sout.str().c_str());
 }
 #endif 	
 NormaF+=auxNormaF*peso*fabs(jacdet);
 NormaFexact+=auxNormaFexact*peso*fabs(jacdet);
 
 }
 delete intr;
 
 }
 else continue;
 
 
 }
 erro<<sqrt(NormaFexact)<< " \t  " <<sqrt(NormaF) << endl;
 //erro<< sqrt(intPexata)<< " \t  " <<sqrt(ePressure) << " \t  " << sqrt(intFexata)<<" \t " << sqrt(fluxoHdiv)<< endl;
 //		return sqrt(ePressure);	
 
 
 }
 
 REAL NormEnergiaH1(TPZCompMesh *malha, std::ofstream &erro){
 int el;
 const int nelem = malha->NElements();
 REAL valor = 0.,valor2=0.,valor4=0.,valor3=0.;
 REAL ProdInt4=0.;
 ///Percorrer todos elementos da malha computacional  
 for(el=0; el < nelem; el++) 
 {
 TPZCompEl * Cel = malha->ElementVec()[el];
 if(!Cel) continue;
 
 TPZGeoEl *gel = Cel->Reference();
 if(	malha->ElementVec()[el]->Dimension()== 2) {
 
 TPZIntPoints *intr = gel-> CreateSideIntegrationRule(gel->NSides() - 1,10);//side and order 
 
 for(int in=0; in < intr->NPoints(); in++)
 {
 REAL peso;
 TPZVec<REAL> pto(2);
 intr->Point(in, pto, peso);
 
 TPZFMatrix<REAL> jac(2,2),invjac(2,2),axes(3,3);
 REAL jacdet;
 gel->Jacobian(pto,jac,axes,jacdet,invjac);
 TPZManVector<REAL> sol(1);
 TPZFMatrix<REAL> dsol(3,1,0.);
 Cel->ComputeSolution(pto, sol, dsol, axes);
 
 TPZVec< REAL > xco(3), p(1);
 TPZVec<REAL> fluxo(3);
 gel->X(pto,xco);
 SolExata(xco, p, fluxo);
 //erro<< "pto "<< xco<< " fluxo "<<fluxo<<std::endl;
 //erro<< "pto "<< xco<< " fluxoAprox "<<dsol<<std::endl;
 REAL ProdInt = (p[0]-sol[0])*(p[0]-sol[0]);
 ProdInt4 =(fluxo[0]-dsol(0,0))*(fluxo[0]-dsol(0,0))+(fluxo[1]-dsol(1,0))*(fluxo[1]-dsol(1,0));
 //	erro<< (fluxo[0]-dsol(0,0))<<" ---- "<<(fluxo[1]-dsol(1,0))<<std::endl;
 //REAL ProdInt = p[0]*p[0];
 //REAL ProdInt2 = sol[0]*sol[0];
 //REAL ProdInt3 = dsol(0,0)*dsol(0,0)+dsol(1,0)*dsol(1,0);
 //ProdInt4 =fluxo[0]*fluxo[0]+fluxo[1]*fluxo[1]+fluxo[2]*fluxo[2];
 valor += ProdInt*peso*fabs(jacdet);//norma p exata
 //	valor2 += ProdInt2*peso*fabs(jacdet);//norma p aprox
 //	valor3 += ProdInt3*peso*fabs(jacdet);//norma flux aprox
 valor4 += ProdInt4*peso*fabs(jacdet);//norma flux exato
 }
 
 delete intr;
 }
 else continue;
 
 }
 erro<< sqrt(valor)<< " \t  " << sqrt(valor4)<< endl;
 // erro<< sqrt(valor)<< " \t  " <<sqrt(valor2) << " \t  " << sqrt(valor4)<< "\t "<<sqrt(valor3)<<endl;
 //	erro<< sqrt(valor)<<endl;
 return sqrt(valor);
 
 }*/

void SaddlePermute(TPZCompMesh * cmesh){
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout<< "Implementando permutacao para problemas de ponto de sela"<< std::endl;
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif
	TPZVec<int64_t> permute;
	int64_t numinternalconnects = cmesh->NIndependentConnects();
  	permute.Resize(numinternalconnects,0);
	
	TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (cmesh);
	if(submesh)
	{
		int64_t nexternal = submesh->NConnects();
		numinternalconnects -= nexternal;
	}
	//	else {
	//		DebugStop();
	//	}
	
	int64_t jperm=0;
	int64_t nel=cmesh->ElementVec().NElements();
	for (int64_t jel=0; jel<nel; jel++) {
		
		for (int64_t ip=0; ip<permute.NElements(); ip++) {
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
void SolveLU ( TPZAnalysis &an ){
	// Com matriz mal condicionada a solução é poluída com resíduos,
	// tanto com LU quanto Choleski. Isso resulta em não simetrias.
	TPZCompMesh *malha = an.Mesh();
	//TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
	TPZFStructMatrix mat( malha );
	//	TPZSpStructMatrix mat( malha );
	TPZStepSolver<STATE> solv;
	
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
TPZGeoMesh * MalhaGeoQ2(const int h){//malha triangulo
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	//Criar ns
	const int nnode = 9;//AQUI
	const int nelem = 4;
	TPZGeoEl *elvec[4];	//nelem
	const int dim = 2;//AQUI
	
	REAL co[nnode][dim] = {{-1.,-1},{1.,-1},{1.,1.},{-1.,1.},{0.,-1.},{1.,0.},{0.,1.},{-1.,0.},{0.,0.}};
	int indices[nelem][nnode];//como serao enumerados os nos
	
	//el 1
	indices[0][0] = 0;
	indices[0][1] = 4;
	indices[0][2] = 8;
	indices[0][3] = 7;
	//el2
	indices[1][0] = 4;
	indices[1][1] = 1;
	indices[1][2] = 5;
	indices[1][3] = 8;
	//el3
	indices[2][0] = 8;
	indices[2][1] = 5;
	indices[2][2] = 2;
	indices[2][3] = 6;
	//el4
	indices[3][0] = 7;
	indices[3][1] = 8;
	indices[3][2] = 6;
	indices[3][3] = 3;
	
	
	int nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		int64_t nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos
	
	//	TPZVec<int> nodind1(4);
	//	TPZVec<int> nodind2(4);
	TPZVec<int64_t> nodind1(4);
	for (int iel=0; iel<nelem; iel++) {
		
		for(int i=0; i<4; i++){
			nodind1[i] = indices[iel][i];
			
		}
		int64_t index;
		elvec[iel] = gmesh->CreateGeoElement(EQuadrilateral,nodind1,1,index);
	}
	
	//	int index;
	//	elvec[0] = gmesh->CreateGeoElement(ETriangle,nodind1,1,index); //AQUI
	//	elvec[1] = gmesh->CreateGeoElement(ETriangle,nodind2,1,index); //AQUI
	
	
	gmesh->BuildConnectivity();
	
	//	como usar os padroes de refinamento do caju..nao posso rodar ainda pq tenho o Pz antigo..nao tenho este gRefDBase.
	//{
	//			TPZAutoPointer< TPZRefPattern > refExemplo = gRefDBase->FinfRefPattern("a_Quad_Rib_Side_4_5_6");
	//				if(refExemplo)
	//			{
	//						TPZGeoEl * elExemplo = gmesh->Elementvec()[10];
	//						elExemplo->SetRefPattern(refExemplo);
	//						elEmemplo->Divide(filhos);
	//				}
	
	//	}
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec[0],4,-1); 
	TPZGeoElBC gbc2(elvec[0],7,-2);
	
	TPZGeoElBC gbc3(elvec[1],4,-3); 
	TPZGeoElBC gbc4(elvec[1],5,-4); 
	
	TPZGeoElBC gbc5(elvec[2],5,-5); 
	TPZGeoElBC gbc6(elvec[2],6,-6);
	
	TPZGeoElBC gbc7(elvec[3],6,-7); 
	TPZGeoElBC gbc8(elvec[3],7,-8);
	
	
	
	const std::string nameref;
	
	TPZAutoPointer<TPZRefPattern> ref;
	
	
	/*	Refinamento uniforme
	 for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
	 TPZVec<TPZGeoEl *> filhos;
	 int n = gmesh->NElements();
	 for(int i = 0; i < n; i++){
	 TPZGeoEl * gel = gmesh->ElementVec()[i];
	 if(!gel->HasSubElement())
	 {
	 gel->Divide(filhos);
	 }		
	 }}
	 */
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif 		 		 
	return gmesh;
}
int64_t SubStructure(TPZCompMesh *cmesh, int materialid)
{
	int64_t index;
	TPZSubCompMesh *submesh = new TPZSubCompMesh(*cmesh,index);//alocacao de memoria...o constructor do tpzsubcompmesh é inicializado com o parametro index que sera o numero de elementos computacionais da malha
	
	int64_t nelem = cmesh->NElements();
	int64_t iel;
	for(iel = 0; iel<nelem; iel++)
	{
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if(!cel || cel == submesh) continue;
		TPZMaterial * celmat = cel->Material();
		if(!celmat) continue;
		int matid = celmat->Id();
		if(matid == materialid)
		{
			submesh->TransferElement(cmesh, iel);//cada elemento de matId=1 eh transferido para a submalha..cada subelemento tera todas as caracterisitcas de um elemento computacional(?)
		}
	}
	cmesh->ComputeNodElCon();//verifica os no connectados
	submesh->MakeAllInternal();//faz as conexoes da malha com os nos internos
	cmesh->CleanUpUnconnectedNodes();//deleta os nos q nao tem elementos conectados
	
	//	submesh->SetAnalysisSkyline(numThreads4Assemble, guiInterface);
	// submesh->SetAnalysis();
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Mesh after substructuring\n";
		cmesh->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	return index;
	
}
