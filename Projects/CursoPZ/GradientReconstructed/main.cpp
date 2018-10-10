/**
 * @file
 * @brief Tutorial program showing a method to reconstruction gradient for a solution precomputed
 */

#include <iostream>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzreferredcompel.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"
#include "pzpoisson3d.h"
#include "pzmat1dlin.h"

#include "pzgradient.h"
#include "pzl2projection.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "TPZGenSpecialGrid.h"

#include <iostream>
#include <math.h>
using namespace std;

/*
 *Projeto para validar a reconstucao do gradiente
 *Equacao: du2dx2 + du2dy2 = 0 em (0,1)x(0,1)
 *Solucao: u(x,y) = a*x + b*y
 */

REAL const coef_a = 2.;
REAL const coef_b = 4.;
TPZVec<REAL> gradU(2,0.);
TPZVec<REAL> normal_plano(2,0.);

//const int NN = 300;

int const matId =1;
//int const matIdL2Proj = 2;
int const bc0=-1; //em y=0
int const bc1=-2; //em x=1
int const bc2=-3; //em y=1
int const bc3=-4; //em x=0

int const dirichlet =0;
//int const neumann = 1;
//int const mixed = 2;

TPZFMatrix<REAL> MatrixR(REAL ang);
TPZGeoMesh *GMesh(int triang_elements, REAL angle, REAL origX, REAL origY, int nh);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int dim,int pOrder,int discont);

void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Forcingbc1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Forcingbc2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Forcingbc3(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
void PosProcessamento(TPZAnalysis &an,TPZCompMesh *Cmesh, std::string plotfile,TPZFMatrix<REAL> &gradients);

//void GradientReconstructionByGreenFormula(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var=0);

/*
 *@brief Method to reconstruction gradient by Least Squares
 *@param cel in:  Computational element
 *@param center out:  Center point of the element
 *@param solalfa out: Value of the approximate solution (uh) at the center point
 *@param grad out: Value of the gradient reconstructed
 *@param var in: Index of the Solution Variable
 */
void GradReconstructionByLeastSquares(TPZCompEl *cel,TPZFMatrix<REAL> &Grad,int var);
void GradReconstructionByLeastSquares_Self(TPZCompEl *cel,TPZFMatrix<REAL> &Grad,int var);

/*
 *@brief Method to replace the solution by finite element method from the reconstruction of the gradient.
 *Using the L2 projection
 *@param cmesh in:  Computational mesh
 *@param var in: Index of the Solution Variable
 *@param matid_l2proj in: Id of the l2 projection material
 */
void ProjectionGradientReconstructedInFESpace(TPZCompMesh *cmesh,int var, int matid_l2proj);

void SaidaMathGradiente(TPZFMatrix<REAL> gradients);

//Trocar todos os elementos do cmesh apontando para o material TPZL2ProjectionFromGradient
void ChangeMaterialIdIntoCompElement(TPZCompEl *cel, int oldmatid, int newmatid);

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=false, const int matidtodivided=1);

// MAIN FUNCTION
int main(int argc, char *argv[]) {
	
#ifdef LOG4CXX
    InitializePZLOG();
#endif
	
	// Initializing uniform refinements for quadrilaterals and triangles
	gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	
    gradU[0]=coef_a;
    gradU[1]=coef_b;
    normal_plano[0]=coef_a;
    normal_plano[1]=coef_b;
    
    REAL anglo = 0;   // M_PI/4.;
    REAL x0 = -1.;
    REAL y0= 2.;
	
    time_t starttime=0, endtime=0;
    int dim=2;
    int p = 3;
    for(int discont=0;discont<2;discont++) {
		for(int nrefs=0;nrefs<2;nrefs++) {
			for(int type=0;type<2;type++) {
				cout << "\n\nGRADIENT RECONSTRUCTION:";
				cout << "Case: continuous: " << (!discont) << "  nrefs: " << nrefs;
				if(type) cout << " Triangular element. ";
				else cout << " Quadrilateral element. ";
				cout << endl;
				
				ctime(&starttime);
				// geometric mesh (initial)
				TPZGeoMesh *gmesh = GMesh(type,anglo,x0,y0,nrefs);
				UniformRefinement(nrefs+1,gmesh,2);
				
				ofstream arg1("gmesh_inicial.txt");
				gmesh->Print(arg1);
				
				// First computational mesh
				TPZCompMesh * cmesh= CMesh(gmesh,dim,p,discont);
				ofstream arg2("cmesh_inicial.vtk");
				TPZVTKGeoMesh::PrintGMeshVTK(gmesh,arg2);
				
				// Solving
				TPZAnalysis an(cmesh);
				mySolve(an,cmesh);
				
				ctime(&endtime);
				cout << "\n Time solving: " << endtime-starttime << endl;
				
				// Pos processing
				TPZStack<std::string> scalarnames, vecnames;
				std::string filename = "Poisson2D";
				char pp[16];
				sprintf(pp,"%d_%d_%d.vtk",discont,nrefs,type);
				filename += pp;
				scalarnames.Push("POrder");
				scalarnames.Push("Solution");
				scalarnames.Push("Pressure");
				
				vecnames.Push("Derivative");
				vecnames.Push("Flux");
				an.DefineGraphMesh(dim,scalarnames,vecnames,filename);
				an.PostProcess(0,dim);				
				
				// Computing approximation of gradient
				TPZFMatrix<REAL> gradients(dim,1);
				for(int i=0;i<cmesh->NElements();i++) {
					TPZCompEl *cel = cmesh->ElementVec()[i];
					if(!cel || cel->Dimension()!=dim) continue;
					gradients.Zero();
					//					if(!discont) 
					GradReconstructionByLeastSquares(cel,gradients,0);
					//					else 
					//						GradReconstructionByLeastSquares_Self(cel,gradients,0);
					cout << "Element " << i << " Grad = (";
					for(int j=0;j<dim;j++) {
						cout << gradients[j];
						if(j<dim-1) cout << ",";
					}
					cout << ")\n";
				}
				delete cmesh;
				delete gmesh;
			}
		}
	}
    return 0;
}

// COMPUTING RECONSTRUCTION OF GRADIENTS TO APPROXIMATED SOLUTION
void GradReconstructionByLeastSquares(TPZCompEl *cel,TPZFMatrix<REAL> &Grad,int var) {
	// Nada sera realizado para elementos con dimension diferente de la dimension del problema
	if(!cel) return;
	int dim = cel->Mesh()->Dimension();
	if(cel->Dimension()!=dim) return;
	TPZStack<TPZCompElSide> neighs;
	int nneighs;
	int nstates = cel->Material()->NSolutionVariables(var);
	int k, side, nsides = cel->Reference()->NSides()-1;   // Desconsiderando o proprio elemento (lado)
	
	TPZManVector<REAL> centerpsi(3,0.0);
	TPZManVector<REAL,3> center(3,0.0), centerbeta(3,0.0);
	TPZManVector<STATE> solalfa(nstates,0.0), solbeta(nstates,0.0);
	TPZFMatrix<REAL> A(dim,dim);    // Linear System matrix
	TPZFMatrix<REAL> B(dim,1,0.);   // Linear System vector
	
	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	
	// Limpiando las matrizes
	A.Zero();
	// Encontramos el centro del elemento corriente cel
	TPZGeoEl* gelalfa = cel->Reference();
	gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
	center.Fill(0.);
	gelalfa->X(centerpsi,center);
	cel->Solution(centerpsi,var,solalfa);
	neighs.Resize(0);
	// Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
	for(side = 0; side < nsides; side++) {
		TPZCompElSide celside(cel,side);
		celside.ConnectedElementList(neighs,0,0);
	}
	nneighs = neighs.NElements();
	if(!nneighs) return;
	
	// If exist neighboard clean duplicated elements and store only index of neighboard
	TPZStack<int> realneighs;
	int kk,jj;
	int id;   // = neighs[0].Element()->Index();
	for(kk=0;kk<nneighs;kk++) {
		id=neighs[kk].Element()->Index();
		for(jj=0;jj<realneighs.NElements();jj++) {
			if(id == realneighs[jj])
				break;
		}
		if(jj==realneighs.NElements() && cel->Mesh()->ElementVec()[id]->Dimension() == dim)
			realneighs.Push(id);
	}
	nneighs = realneighs.NElements();
	// si no hay vecinos continuamos con el siguiente elemento
	if(!nneighs) return;
	
	// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente			
	// Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
	// y el valor de la solucion en su centro solbeta
	DeltaH.Redim(nneighs,dim);
	DeltaHTranspose.Redim(dim,nneighs);
	DifSol.Redim(nneighs,1);
	// Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
	for(int ineighs=0;ineighs<nneighs;ineighs++) {
		TPZGeoEl * gelbeta = cel->Mesh()->ElementVec()[realneighs[ineighs]]->Reference();
		if(!gelbeta)
			DebugStop();
		centerpsi.Fill(0.0);
		centerbeta.Fill(0.0);
		gelbeta->CenterPoint(gelbeta->NSides()-1,centerpsi);
		gelbeta->X(centerpsi,centerbeta);
		gelbeta->Reference()->Solution(centerpsi,var,solbeta);
		
		for(k=0;k<dim;k++)
			DeltaH(ineighs,k) = centerbeta[k] - center[k];
		DifSol(ineighs,0) = solbeta[var] - solalfa[var];
	}
	// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u) 
	DeltaH.Transpose(&DeltaHTranspose);
	Grad = DeltaHTranspose*DifSol;
	A = DeltaHTranspose*DeltaH;
	A.SolveDirect(Grad,ELU);
}
void GradReconstructionByLeastSquares_Self(TPZCompEl *cel,TPZFMatrix<REAL> &Grad,int var) {
	// Nada sera realizado para elementos con dimension diferente de la dimension del problema
	if(!cel) return;
	int dim = cel->Mesh()->Dimension();
	if(cel->Dimension()!=dim) return;
	int nstates = cel->Material()->NSolutionVariables(var);
	int side, nsides = cel->Reference()->NSides()-1;   // Desconsiderando o proprio elemento (lado)
	
	TPZManVector<REAL> centerpsi(3,0.0);
	TPZManVector<REAL,3> center(3,0.0), centerbeta(3,0.0);
	TPZManVector<STATE> solalfa(nstates,0.0), solbeta(nstates,0.0);
	TPZFMatrix<REAL> A(dim,dim);    // Linear System matrix
	TPZFMatrix<REAL> B(dim,1,0.);   // Linear System vector
	
	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	
	// Limpiando las matrizes
	A.Zero();
	// Encontramos el centro del elemento corriente cel
	TPZGeoEl* gelalfa = cel->Reference();
	gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
	center.Fill(0.);
	gelalfa->X(centerpsi,center);
	cel->Solution(centerpsi,var,solalfa);
	
	// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente			
	// Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
	// y el valor de la solucion en su centro solbeta
	DeltaH.Redim(nsides,dim);
	DeltaHTranspose.Redim(dim,nsides);
	DifSol.Redim(nsides,1);
	// Para cada lado calculamos los deltaH (desde el centro del elemento al centro del lado de dimension menor a él
	// y el valor de la solucion en su centro solbeta
	DeltaH.Redim(nsides,dim);
	DeltaHTranspose.Redim(dim,nsides);
	DifSol.Redim(nsides,1);
	// Procuramos todos los puntos medios de cada lado del elemento y calculamos baseados en los valores de la solucion sobre ellos
	for(side = 0; side < nsides; side++) {
		centerpsi.Fill(0.0);
		centerbeta.Fill(0.0);
		cel->Reference()->CenterPoint(side,centerpsi);
		cel->Reference()->X(centerpsi,centerbeta);
		cel->Solution(centerpsi,var,solbeta);
		for(int k=0;k<dim;k++)
			DeltaH(side,k) = centerbeta[k] - center[k];
		DifSol(side,0) = solbeta[var] - solalfa[var];
		
	}
	// Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
	// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u) 
	DeltaH.Transpose(&DeltaHTranspose);
	Grad = DeltaHTranspose*DifSol;
	A = DeltaHTranspose*DeltaH;
	A.Print("A");
	Grad.Print("Grad");
	A.SolveDirect(Grad,ELU);
}


// CONSTRUCTION OF MESHES - GEOMETRIC AND COMPUTATIONAL
TPZGeoMesh *GMesh(int triang_elements, REAL angle, REAL origX, REAL origY, int nh){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
    TPZFMatrix<REAL> mA(2,4,0.), mR(2,2), mRA(2,4);
    mA.Put(0,0,origX); mA.Put(1,0,origY);
    mA.Put(0,1,origX+1.); mA.Put(1,1,origY);
    mA.Put(0,2,origX+1.); mA.Put(1,2,origY+1.);
    mA.Put(0,3,origX); mA.Put(1,3,origY+1.);
    
    mR = MatrixR(angle);
	mR.Multiply(mA, mRA);
    
    int64_t id;
	id = 0;
	for (int j=0; j<4;j++) {
		
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0, mRA.GetVal(0, j));//coord x
		Node[id].SetCoord(1, mRA.GetVal(1, j));//coord y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	//indice dos elementos
	id = 0;
    
    if(triang_elements==1)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
    for ( int ref = 0; ref < nh; ref++ ) {
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ) {
			TPZGeoEl * gel = gmesh->ElementVec()[i];
            gel->Divide (filhos);
		}
	}
	
	return gmesh;
}
TPZCompMesh *CMesh(TPZGeoMesh *gmesh,int dim, int pOrder,int discont)
{
    /// criar materiais
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
    
    REAL diff=1.;
    REAL conv =0.;
    TPZVec<REAL> convdir;
    convdir.Resize(dim,0.);
    convdir[0]=0.;
    material-> SetParameters(diff, conv, convdir);
    
    REAL ff=0.;
    material->SetInternalFlux(ff);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > fCC0 = new TPZDummyFunction<STATE>(Forcingbc0, 5);
    TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(Forcingbc2, 5);
    TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(Forcingbc1, 5);
    TPZAutoPointer<TPZFunction<STATE> > fCC3 = new TPZDummyFunction<STATE>(Forcingbc3, 5);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    
    
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    BCond0->SetForcingFunction(fCC0);
    BCond2->SetForcingFunction(fCC2);
    BCond1->SetForcingFunction(fCC1);
    BCond3->SetForcingFunction(fCC3);
    
	if(discont) cmesh->SetAllCreateFunctionsDiscontinuous();
	else cmesh->SetAllCreateFunctionsContinuous();
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->ExpandSolution();
	cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;
}

void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= coef_a*x + coef_b*y;
}

void Forcingbc1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= coef_a*x + coef_b*y;
}

void Forcingbc2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= coef_a*x + coef_b*y;
}

void Forcingbc3(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    disp[0]= coef_a*x + coef_b*y;
}

#include "pzbstrmatrix.h"
void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	TPZBandStructMatrix full(Cmesh); //caso nao-simetrico
	//TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	//step.SetDirect(ELDLt); //caso simetrico
	step.SetDirect(ELU);//caso nao simetrico
	an.SetSolver(step);
	an.Run();
}


void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial, const int matidtodivided) {
	TPZManVector<TPZGeoEl*> filhos;
    for(int D=0; D<nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {    
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
			if(!gel || gel->HasSubElement())
				continue;
			if(dim > 0 && gel->Dimension() != dim) continue;
			if(!allmaterial){
				if(gel->MaterialId() == matidtodivided){
					gel->Divide(filhos);
				}
			}
			else{
				gel->Divide(filhos);
			}
        }
    }
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}

TPZFMatrix<REAL> MatrixR(REAL ang)
{
	TPZFMatrix<REAL> r(2,2,0.);
	r(0,0) = cos(ang); r(0,1) = -sin(ang);
	r(1,0) = sin(ang);	r(1,1) =  cos(ang);
    
	return r;
}

