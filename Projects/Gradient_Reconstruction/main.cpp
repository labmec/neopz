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

#include <iostream>
#include <math.h>
using namespace std;

/**
 * Projeto para validar a reconstucao do gradiente
 * Equacao: du2dx2 + du2dy2 = 0 em (0,1)x(0,1)
 * Solucao: u(x,y) = a*x + b*y
 */

REAL const coef_a = 2.;
REAL const coef_b = 4.;
TPZVec<REAL> gradU(2,0.);
TPZVec<REAL> normal_plano(2,0.);

const int NN = 300;

int const matId =1;
int const bc0=-1; //em y=0
int const bc1=-2; //em x=1
int const bc2=-3; //em y=1
int const bc3=-4; //em x=0

int const dirichlet =0;
int const neumann = 1;
int const mixed = 2;

TPZFMatrix<REAL> MatrixR(REAL ang);
TPZGeoMesh *GMesh(int triang_elements, REAL angle, REAL origX, REAL origY, int nh);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);

void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc2(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc3(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);

void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
void PosProcessamento(TPZAnalysis &an, std::string plotfile,TPZFMatrix<REAL> &gradients);

void GradientReconstructionByGreenFormula(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var=0);
void GradientReconstructionByLeastSquares(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var=0,bool continuous = false);

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
    
    REAL anglo = M_PI/4.;
    REAL x0 = -1.;
    REAL y0= 5.;
    
    int p = 4;
	bool continuous = false;
	//primeira malha
	do {
		for(int nrefs=1;nrefs<4;nrefs++) {
			for(int type=0;type<2;type++) {
				cout << "\nCase: continuous: " << continuous << "  nrefs: " << nrefs << "   type: " << type << endl;
				// geometric mesh (initial)
				TPZGeoMesh * gmesh = GMesh(type,anglo,x0,y0,nrefs);
				ofstream arg1("gmesh_inicial.txt");
				gmesh->Print(arg1);
				
				// First computational mesh
				TPZCompMesh * cmesh= CMesh(gmesh,p);
				ofstream arg2("cmesh_inicial.txt");
				cmesh->Print(arg2);
				
				// Solving
				TPZAnalysis an(cmesh);
				mySolve(an, cmesh);
				
				// Computing approximation of gradient
				/** 
				 * @brief Method to reconstruct a gradient after run Solve of the analysis
				 * @param cmesh Computational mesh with solution */
				TPZFMatrix<REAL> gradients;
				GradientReconstructionByLeastSquares(gradients,cmesh,0,0,continuous);
				gradients.Print();

				// Print gradient reconstructed
				string plotfile("GradientAndSolution.vtk");
				PosProcessamento(an,plotfile,gradients);
			}
		}
		continuous = !continuous;
	} while(continuous);
	
    return 0;
}

/** Reconstrucción del gradiente utilizando la linearizacion (Taylor) de la solución para los centros de todos los elementos vecinos */
/** Formula: u(xbi,ybi,zbi) = u(xa,ya,za) + a*(xbi-xa) + b*(ybi-ya) + c*(zbi-za)  ->  donde Grad(u) ~= (a,b,c) */
/** (xa,ya,za) es el centro del elemento donde queremos aproximar o gradiente de u */
/** (xbi,ybi,zbi) son los centros de los elementos vecinos al elemento corriente por alguno de sus lados, e enumerados por i */
void GradientReconstructionByLeastSquares(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var,bool continuous) {
	int i, nstates=0;
	TPZCompEl *cel;
	int dim = cmesh->Dimension();
	for(i=0;i<cmesh->NElements();i++) {
		cel = cmesh->ElementVec()[i];
		if(cel && cel->Dimension() == dim) {
			nstates = cel->Material()->NSolutionVariables(var);
			break;
		}
	}
	
	// Redimensionando a matriz dos gradientes
	int nelem = cmesh->NElements();
    gradients.Redim(nelem,2*dim);
	
	int k, side;
	int counter = 0;
	
	TPZStack<TPZCompElSide> neighs;
	int nneighs;
	
	TPZManVector<REAL> normal(3,0.0);
	TPZManVector<REAL> centerpsi(3,0.0);
	TPZManVector<REAL> center(3,0.0), centerbeta(3,0.0);
	TPZManVector<REAL> solalfa(nstates,0.0), solbeta(nstates,0.0);
	
	TPZFMatrix<REAL> A(dim,dim);    // Linear System matrix
	TPZFMatrix<REAL> B(dim,1,0.);   // Linear System vector
	
	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> DeltaH(nneighs,dim,0.);
	TPZFMatrix<REAL> DeltaHTranspose(dim,nneighs,0.);
	TPZFMatrix<REAL> DifSol(nneighs,1,0.);
	REAL Grad;
	
	// Calculando el gradiente por elemento computacional
	for(i=0;i<nelem;i++) {
		cel = cmesh->ElementVec()[i];
		// Nada sera realizado para elementos con dimension diferente de la dimension del problema
		if(!cel || cel->Dimension()!=dim) continue;
		
		// Limpiando las matrizes
		A.Zero(); B.Zero();
		// Encontramos el centro del elemento corriente cel
		TPZGeoEl* gelalfa = cel->Reference();
		gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
		center.Fill(0.);
		gelalfa->X(centerpsi,center);
		cel->Solution(centerpsi,var,solalfa);
		
		// PREFERENCIAL PARA CASOS DE CALCULO CON FUNCIONES DISCONTINUAS - Pues utiliza los valores de la solución en los elementos vecinos
		if(!continuous) {
			neighs.Resize(0);
			// Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
			for(side = cel->Reference()->NCornerNodes(); side < cel->NConnects(); side++) {
				TPZCompElSide celside(cel,side);
				celside.ConnectedElementList(neighs,1,0);
			}
			nneighs = neighs.NElements();
			// si no hay vecinos continuamos con el siguiente elemento
			if(!nneighs) continue;
			// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente			
			// Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
			// y el valor de la solucion en su centro solbeta
			DeltaH.Redim(nneighs,dim);
			DeltaHTranspose.Redim(dim,nneighs);
			DifSol.Redim(nneighs,1);
			// Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
			for(int ineighs=0;ineighs<nneighs;ineighs++) {
				TPZGeoEl* gelbeta = neighs[ineighs].Element()->Reference();
				if(!gelbeta)
					DebugStop();
				centerpsi.Fill(0.0);
				centerbeta.Fill(0.0);
				gelbeta->CenterPoint(gelbeta->NSides()-1,centerpsi);
				gelbeta->X(centerpsi,centerbeta);
				gelbeta->Reference()->Solution(centerpsi,var,solbeta);
				for(k=0;k<dim;k++)
					DeltaH(ineighs,k) = centerbeta[k] - center[k];
				DifSol(ineighs,0) = solbeta[n_var] - solalfa[n_var];
			}
		}
		else {
			int nsides = cel->NConnects()-1;
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
				for(k=0;k<dim;k++)
					DeltaH(side,k) = centerbeta[k] - center[k];
				DifSol(side,0) = solbeta[n_var] - solalfa[n_var];
				
			}
		}
		// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u) 
		DeltaH.Transpose(&DeltaHTranspose);
		B = DeltaHTranspose*DifSol;
		A = DeltaHTranspose*DeltaH;
		A.SolveDirect(B,ELU);
		
		// Normalizando el vector gradiente
		Grad = 0.0;
		for(k=0;k<dim;k++)
			Grad += (B(k,0)*B(k,0));
		// Almacenando los gradientes encontrados
		for(k=0;k<dim;k++) {
			if(!IsZero(B(k))) {
				gradients(counter,k) = B(k,0)/sqrt(Grad);
			}
			gradients(counter,dim+k) = center[k];
		}
		counter++;
	}
	// Redimensionando la matriz de los gradientes
	gradients.Resize(counter,2*dim);
}

void GradientReconstructionByGreenFormula(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var) {
	int i, nstates;
	TPZCompEl *cel;
	for(i=0;i<cmesh->NElements();i++) {
		cel = cmesh->ElementVec()[i];
		if(cel->Dimension() == cmesh->Dimension())
			nstates = cel->Material()->NSolutionVariables(var);
	}
    
    gradients.Redim(cmesh->NElements(),cmesh->Dimension());
    
	REAL integral;
	int side;
	TPZManVector<REAL> normal(3,0.0);
	TPZManVector<REAL> center(3,0.0);
	TPZManVector<REAL> solalfa(nstates,0.), solbeta(nstates,0.);
	TPZManVector<REAL> centerside(3,0.);
	
	REAL measure, sidemeasure;
	for(i=0;i<cmesh->NElements();i++) {
		cel = cmesh->ElementVec()[i];
		
		// WRONG
		
		
		if(!cel || cel->Dimension()!=cmesh->Dimension()) continue;
		TPZStack<TPZCompElSide> neighs;
		measure = cel->VolumeOfEl();
		for(side = cel->Reference()->NCornerNodes(); side < cel->NConnects()-1; side++) {
			neighs.Resize(0);
			TPZCompElSide celside(cel,side);
			celside.ConnectedElementList(neighs,1,0);
			if(!neighs.NElements())
				DebugStop();
			// Aqui deveria ser feita a integracao sobre o lado
			TPZManVector<REAL> centersidepsi(3,0.);
			celside.Reference().CenterPoint(centersidepsi);
			celside.Reference().X(centersidepsi,centerside);
			cel->Solution(centerside,var,solalfa);
			sidemeasure = celside.Reference().Area();
			// Calculo da normal no ponto do meio do lado, exterior ao elemento
			celside.Reference().Normal(centersidepsi,cel->Reference(),neighs[0].Reference().Element(),normal);
			integral = sidemeasure*solalfa[n_var]/measure;
			// Incrementando no gradients para o elemento alfa o valor por vizinho
			for(int j=0;j<cel->Reference()->Dimension();j++){
				gradients(i,j) += integral*normal[j];
			}
		}
	}
}

void PosProcessamento(TPZAnalysis &an, std::string plotfile,TPZFMatrix<REAL> &gradients){
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
    vecnames[0] = "Derivate";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	// Carregar gradients como soluçao na malla computacional
	// Imprimir solucao novamente - Acertar gradients dependente dos graus de liberdade no cmesh
	// cmesh->LoadSolution(gradients);
	//an.PostProcess(div,dim);
	
    //	std::ofstream out("malha.txt");
	//	an.Print("nothing",out);
}

TPZFMatrix<REAL> MatrixR(REAL ang)
{
	TPZFMatrix<REAL> r(2,2,0.);
	r(0,0) = cos(ang); r(0,1) = -sin(ang);
	r(1,0) = sin(ang);	r(1,1) =  cos(ang);
    
	return r;
}

TPZGeoMesh *GMesh(int triang_elements, REAL angle, REAL origX, REAL origY, int nh){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
    TPZVec <int> TopolTriang(3);
	TPZVec <int> TopolLine(2);
	
	//indice dos nos
    TPZFMatrix<REAL> mA(2,4,0.), mR(2,2), mRA(2,4);
    mA.Put(0,0,origX); mA.Put(1,0,origY);
    mA.Put(0,1,origX+1.); mA.Put(1,1,origY);
    mA.Put(0,2,origX+1.); mA.Put(1,2,origY+1.);
    mA.Put(0,3,origX); mA.Put(1,3,origY+2.);
    
    mR = MatrixR(angle);
	mR.Multiply(mA, mRA);
    
    int id;
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
    
    for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			//if (gel->Dimension() == 2) gel->Divide (filhos);
            gel->Divide (filhos);
		}//for i
	}//ref
	
	return gmesh;
}


TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
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
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > fCC0 = new TPZDummyFunction<STATE>(Forcingbc0);
    TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(Forcingbc2);
    TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(Forcingbc1);
    TPZAutoPointer<TPZFunction<STATE> > fCC3 = new TPZDummyFunction<STATE>(Forcingbc3);
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    
	
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	
	//    val2(0,0)=coef_a;
	//    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,neumann, val1, val2);
	//    val2(0,0)=-coef_a;
	//    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,neumann, val1, val2);
    
    BCond0->SetForcingFunction(fCC0);
    BCond2->SetForcingFunction(fCC2);
    BCond1->SetForcingFunction(fCC1);
//    BCond3->SetForcingFunction(fCC3);
    
	cmesh->SetAllCreateFunctionsContinuous();
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;
}


void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= coef_a*x + coef_b*y;
}

void Forcingbc1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= coef_a*x + coef_b*y;
}

void Forcingbc2(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= coef_a*x + coef_b*y;
}

void Forcingbc3(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
    double x = pt[0];
    double y = pt[1];
    disp[0]= coef_a*x + coef_b*y;
}


void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

