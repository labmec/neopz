/*
 *  Tools.cpp
 *  PZ
 *
 *  Created by Denise de Siqueira on 9/5/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "Tools.h"
#include "pzpoisson3d.h"
#include"pzgeoelbc.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"




#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("HpAdaptivity.main"));

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
const REAL Pi=4.*atan(1.);
void Forcing1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    double x = pt[0];
    double y = pt[1];
    //disp[0]=-1/(4.*pow(pow(x,2) + pow(y,2),0.75));
    disp[0]= 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);//2.*pow(Pi,2.)*cos(Pi*y)*sin(Pi*x);//(1.)*8.;//-2.*exp(x)*(1. + 4.*x + pow(x,2.))*(-1. + pow(y,2.));//(exp(x)*(-3. + pow(y,2.) + x*(-4. + x + (4. + x)*pow(y,2.))));//2.*(1.-x*x) +2.*(1.-y*y); //
    return;
}
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux ) {
    double x = pt[0];
    double y = pt[1];
    p[0]= sin(Pi*x)*sin(Pi*y);//(1.)*pow(pow(x,2) + pow(y,2),0.25);-sin(Pi*x)*cos(Pi*y);//(1-x*x)*(1-y*y)*exp(x);//(1-x*x)*(1-y*y);//Solucao
    flux(0,0)= (-1.)*Pi*cos(Pi*x)*sin(Pi*y);//(1.)*x/(2.*pow(pow(x,2) + pow(y,2),0.75));//Pi*cos(Pi*x)*cos(Pi*y);//2.*exp(x)*x*(1. - pow(y,2.)) - exp(x)*(1. - pow(x,2.))*(1. - pow(y,2.));//2*x*(1-y*y);//
    flux(1,0)= (-1.)*Pi*cos(Pi*y)*sin(Pi*x);//(1.)*y/(2.*pow(pow(x,2) + pow(y,2),0.75)); //Pi*sin(Pi*y)*sin(Pi*x);//2.*exp(x)*(1. - pow(x,2.))*y;//2*(1-x*x)*y; dy
    flux(2,0)= 2*Pi*Pi*sin(Pi*x)*sin(Pi*y);//pow(pow(x,2) + pow(y,2),0.25);//-2.*pow(Pi,2.)*sin(Pi*x)*cos(Pi*y);//coloco o divergetne aq para testar
    
    
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

TPZCompMesh *CompMeshPAdap(TPZGeoMesh &gmesh,int porder,bool prefine){
		

		TPZCompMesh *comp = new TPZCompMesh(&gmesh);
		
		
		comp->SetDefaultOrder(porder);
		// Criar e inserir os materiais na malha
		TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
		TPZMaterial * automat(mat);
		comp->InsertMaterialObject(automat);
		
		
    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
		mat->SetForcingFunction(force1);
		TPZAutoPointer<TPZFunction<STATE> > exata1 = new TPZDummyFunction<STATE>(SolExata);
		mat->SetForcingFunctionExact(exata1);
			
		///Criar condicoes de contorno
		
		TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
		TPZFMatrix<STATE> val11(1,1,0.), val22(1,1,0.);
		TPZMaterial *bnd = automat->CreateBC (automat,-1,0,val1,val2);
		TPZMaterial *bnd2 = automat->CreateBC (automat,-2,0,val1,val2);
		TPZMaterial *bnd3 = automat->CreateBC (automat,-3,0,val1,val2);
		TPZMaterial *bnd4 = automat->CreateBC (automat,-4,0,val1,val2);
		
			comp->SetAllCreateFunctionsHDiv();
	//	comp->SetAllCreateFunctionsContinuous();
		
		///Inserir condicoes de contorno
		comp->InsertMaterialObject(bnd);
		comp->InsertMaterialObject(bnd2);
		comp->InsertMaterialObject(bnd3);
		comp->InsertMaterialObject(bnd4);	
    
		
		// Ajuste da estrutura de dados computacional
		comp->AutoBuild();
		comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
		comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
		
		comp->SetName("Malha Computacional com ordem inicializada");	
#ifdef LOG4CXX
		{
				std::stringstream sout;
				comp->Print(sout);
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	
		
		int nel = comp->NElements();
    int iel;
    for(iel=0; iel<nel; iel++){
				
				TPZInterpolatedElement *intel;
				TPZCompEl *cel = comp->ElementVec()[iel];

				if (prefine) {
						
				intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                    if(intel)
                    {
                        int ncon = intel->NConnects();
                        for (int i=0; i<ncon; i++)
                        {
                                int indexcon = intel->ConnectIndex(i);
                                if (indexcon%2==0) {
                                        intel->PRefine(porder+1);
                                        //comp->ExpandSolution();
                                }
                                
                        
                        }
                    }
				}
		}
		  	
		comp->ExpandSolution();
		
		comp->SetName("Malha Computacional com ordens diferentes");	
#ifdef LOG4CXX
		{
				std::stringstream sout;
				comp->Print(sout);
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		
		
    return comp;
		
}



TPZGeoMesh * MalhaGeoTetraedro(const int h,bool hrefine){//malha triangulo
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    //Criar ns
    const int nnode = 4;//AQUI
    const int nelem = 1;
    TPZGeoEl *elvec[nelem];
    const int dim = 3;//AQUI
    
    REAL co[nnode][dim] ={{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
    int indices[3][nnode];//como serao enumerados os nos
    
    //el 1
    indices[0][0] = 0;
    indices[0][1] = 1;
    indices[0][2] = 2;
    indices[0][3] = 3;
    
    
    int nod;
    TPZVec<REAL> coord(dim);
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        nodind = nnode-(nodind+3)%nnode;
        
        for(int d = 0; d < dim; d++)
        {
            coord[d] = co[nod][d];
        }
        gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
    }
    //Criacao de elementos
    
    
		TPZVec<long> nodind1(4);
		TPZVec<long> nodind2(4);
		for(int i=0; i<4; i++){
				nodind1[i] = indices[0][i];
        
    }
    
		long index;
		elvec[0] = gmesh->CreateGeoElement(ETetraedro,nodind1,1,index); //AQUI
    
    
    
    gmesh->BuildConnectivity();
    
    
    const std::string nameref;
    
    TPZAutoPointer<TPZRefPattern> ref;
    
    
		//	Refinamento uniforme
		for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
				TPZVec<TPZGeoEl *> filhos;
				long n = gmesh->NElements();
				for(long i = 0; i < n; i++){
						TPZGeoEl * gel = gmesh->ElementVec()[i];
						if(!gel->HasSubElement())
						{
								gel->Divide(filhos);
						}		
				}}
    
    if (hrefine) {
        //refinamento diferencialvel
        {

            TPZVec<TPZGeoEl *> filhos;
            long n = gmesh->NElements();
                       
            for(long i = 0; i < n; i++){	
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
            long n = gmesh->NElements();
            for(long i = 0; i < n; i++){
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



void SolveLU ( TPZAnalysis &an ){

		TPZCompMesh *malha = an.Mesh();
		//TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
		TPZFStructMatrix mat( malha );
		//	TPZSpStructMatrix mat( malha );
		TPZStepSolver<STATE> solv;
		
		solv.SetDirect ( ELU );
		//		solv.SetDirect(ECholesky);
		
		std::cout << "ELU " << std::endl;
		an.SetSolver ( solv );
		an.SetStructuralMatrix ( mat );
		std::cout << std::endl;
		an.Solution().Redim ( 0,0 );
		std::cout << "Assemble " << std::endl;
		an.Assemble();
		an.Solve();
		std::cout << std::endl;
		std::cout << "No equacoes = " << malha->NEquations() << std::endl;
		// cout << "Banda = " << malha->BandWidth() << endl;
}

void ChangeP(TPZCompMesh * cmesh, TPZCompEl * cel, int newP)
{
		//cmesh->RemoveCompElfromPorderContainer(cel);
		//alterar P
		TPZInterpolationSpace * sp = dynamic_cast<TPZInterpolatedElement*>(cel);
		if(sp)
				{
						sp->PRefine(newP);
						}
		else 
				{
						//FUDEU!!!
						DebugStop();
						}
		//cmesh->AddCompElfromPorderContainer(cel);
}

