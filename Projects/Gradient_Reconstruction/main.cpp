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
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <math.h>
using namespace std;

#ifdef USING_BOOST
#include <boost/math/special_functions/erf.hpp> //Required for erfc function on windows
#endif

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
int const matIdL2Proj = 2;
int const matInterf = 3;
int const bc0=-1; //em y=0
int const bc1=-2; //em x=1
int const bc2=-3; //em y=1
int const bc3=-4; //em x=0

int const dirichlet =0;
int const neumann = 1;
int const mixed = 2;


void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);
void RefiningNearLine(int dim,TPZGeoMesh *gmesh,int nref);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance);



TPZFMatrix<REAL> MatrixR(REAL ang);
TPZGeoMesh *GMesh(int triang_elements, REAL angle, REAL origX, REAL origY, int nh);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, bool isdiscontinuous);

TPZGeoMesh *GMesh2();
TPZCompMesh *CMesh2(TPZGeoMesh *gmesh, int pOrder,bool isdiscontinuous);

void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc2(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc3(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);

void SolucaoExata(const TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix<REAL> &deriv);

void ForcingF(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);

void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
void PosProcessamento(TPZAnalysis &an,TPZCompMesh *Cmesh, std::string plotfile);

void GradientReconstructionByGreenFormula(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var=0);

/*
 *@brief Method to reconstruction gradient by Least Squares
 *@param cel in:  Computational element
 *@param center out:  Center point of the element
 *@param solalfa out: Value of the approximate solution (uh) at the center point
 *@param grad out: Value of the gradient reconstructed
 *@param var in: Index of the Solution Variable
 */
void GradReconstructionByLeastSquares(TPZCompEl *cel, TPZManVector<REAL,3> &center, TPZManVector<REAL> &solalfa, TPZFMatrix<REAL> &Grad, int var);
void GradReconstByLeastSquaresOnlyContEl(TPZCompEl *cel,TPZManVector<REAL,3> &center, TPZManVector<REAL> &solalfa,  TPZFMatrix<REAL> &grad, int var);
//void GradientReconstructionByLeastSquares(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var=0,bool continuous = false);
void PosProcessGradientReconstruction(TPZCompMesh *cmesh,int var,TPZFMatrix<REAL> &datagradients);

/*
 *@brief Method to replace the solution by finite element method from the reconstruction of the gradient.
 *Using the L2 projection
 *@param cmesh in:  Computational mesh
 *@param var in: Index of the Solution Variable
 *@param matid_l2proj in: Id of the l2 projection material
 */
void ProjectionGradientReconstructedInFESpace(TPZCompMesh *cmesh,int var, int matid_l2proj);

void AssembleGlobalMatrix(TPZCompEl *el, TPZElementMatrix &ek, TPZElementMatrix &ef,TPZMatrix<STATE> & stiffmatrix, TPZFMatrix<STATE> &rhs);

void SaidaMathGradiente(TPZFMatrix<REAL> gradients);

//Trocar todos os elementos do cmesh apontando para o material TPZL2ProjectionFromGradient
void ChangeMaterialIdIntoCompElement(TPZCompEl *cel, int oldmatid, int newmatid);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
	
	// Initializing uniform refinements for quadrilaterals and triangles
    //gRefDBase.InitializeRefPatterns();
	gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	
    gradU[0]=coef_a;
    gradU[1]=coef_b;
    normal_plano[0]=coef_a;
    normal_plano[1]=coef_b;
    
    int p = 1;
    char saida[256];
    ofstream erro("erro.txt");
	int dim = 2;
    int MaxRefs = 2;
    for(int nrefs=1;nrefs<MaxRefs;nrefs++)
    {
         erro << "\n\nRESOLVENDO COM FEM: p = " << p << "  E h = "<< nrefs << "\n";
        
        // geometric mesh (initial)
        //TPZGeoMesh * gmesh = GMesh(2,anglo,x0,y0,nrefs);
        TPZGeoMesh * gmesh = GMesh2();
		// Refining near the points belong a circunference with radio r - maxime distance radius
		RefiningNearLine(dim,gmesh,nrefs);
        
		if(nrefs == MaxRefs-1) {
			sprintf(saida,"initialgmesh_%d.vtk",nrefs);
			PrintGeoMeshVTKWithDimensionAsData(gmesh,saida);
		}
        
        // First computational mesh
        //TPZCompMesh * cmesh= CMesh(gmesh,p,true);
        TPZCompMesh * cmesh= CMesh2(gmesh,p,true);
        
        // Solving
        TPZAnalysis an(cmesh);
        mySolve(an,cmesh);
        //TPZFMatrix<REAL> oldsolution = cmesh->Solution();
        
        //2. Plotar as normas do erros
        an.SetExact(SolucaoExata);
        TPZVec<REAL> posproc;
        cout<<"ERROS FEM\n";
        //an.PostProcess(posproc,cout);
        an.PostProcessError(posproc,erro);
        cout<<endl;

        
        // Print solution by FEM
		char saidaVTK[256];
		sprintf(saidaVTK,"FEMSolution_%d.vtk",nrefs);
        string plotfile(saidaVTK);
        PosProcessamento(an,cmesh,plotfile);
        
        erro << "\n\nRESOLVENDO COM RECONST. GRADIENT: p = " << p << "  E h = "<< nrefs << "\n";
       
        //l2 projection of the gradient into finite element space  
        ProjectionGradientReconstructedInFESpace(cmesh,1, matIdL2Proj);

        an.LoadSolution(cmesh->Solution());
        //2. Plotar as normas do erros
        an.SetExact(SolucaoExata);
        TPZVec<REAL> posproc2;
        cout<<"ERROS L2\n";
        an.PostProcessError(posproc2,erro);
        cout<<endl;

		sprintf(saidaVTK,"L2PROJSolution_%d.vtk",nrefs);
        string plotfile2(saidaVTK);
        PosProcessamento(an,cmesh,plotfile2);

		//delete cmesh;
		//delete gmesh;
        cmesh->CleanUp();
        delete gmesh;

        
    }
	erro.close();
    return 0;
}

void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename) {
	int i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		if(gmesh->ElementVec()[i])
			DataElement[i] = (gmesh->ElementVec()[i])->Dimension();
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}
void RefiningNearLine(int dim,TPZGeoMesh *gmesh,int nref) {
	
	int i;
	
	// Refinando no local desejado
	TPZVec<REAL> point(3);
	point[0] = 1.; point[1] = point[2] = 0.0;
	REAL r = 0.0;
	
	REAL radius = 0.75;
	for(i=0;i<nref;i++) {
		// To refine elements with center near to points than radius
		RefineGeoElements(dim,gmesh,point,r,radius);
		radius *= 0.75;
	}
	// Constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();

	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		REAL centerdist = abs(center[0] - point[0]);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
}


void ChangeMaterialIdIntoCompElement(TPZCompEl *cel, int oldmatid, int newmatid) {
    
    // Changes material Id only elements with required id (matid)
    if(cel->Material()->Id() != oldmatid) return;

    //mudar o material id
    TPZGeoEl *gel;
    gel = cel->Reference();
    gel->SetMaterialId(newmatid);
}

ofstream outfile1("SaidaGradiente.nb");
void SaidaMathGradiente(TPZFMatrix<REAL> gradients){
    
    outfile1<<" data = {";
    for(int i = 0;  i< gradients.Rows(); i++)
    {
        outfile1<<"{" << gradients(i,1)<< ", " << gradients(i,2)<<"}";
        if(i!= gradients.Rows()-1) outfile1<<", ";
        if(i== gradients.Rows()-1) outfile1<<"};"<<std::endl;
    }
    outfile1<<"ListPlot[data, Joined -> True, PlotRange -> All, Frame -> True]"<<endl;
}



/** Reconstrucción del gradiente utilizando la linearizacion (Taylor) de la solución para los centros de todos los elementos vecinos */
/** Formula: u(xbi,ybi,zbi) = u(xa,ya,za) + a*(xbi-xa) + b*(ybi-ya) + c*(zbi-za)  ->  donde Grad(u) ~= (a,b,c) */
/** (xa,ya,za) es el centro del elemento donde queremos aproximar o gradiente de u */
/** (xbi,ybi,zbi) son los centros de los elementos vecinos al elemento corriente por alguno de sus lados, e enumerados por i */

void GradientReconstructionByLeastSquares(TPZCompEl *cel,TPZManVector<REAL,3> &center, TPZManVector<REAL> &solalfa,  TPZFMatrix<REAL> &grad, int var) {
    
    int dim;
    dim = cel->Mesh()->Dimension();
    
    // Nada sera realizado para elementos com dimensao diferente da dimensao do problema
    if(!cel || cel->Dimension()!=dim) DebugStop();
    
	int nstates=0;
    nstates = cel->Material()->NSolutionVariables(var);
    
    int k, side;
    
	TPZStack<TPZCompElSide> neighs;
	int nneighs;
	
    center.Resize(3, 0.);
    solalfa.Resize(nstates,0.0);
	TPZManVector<REAL> centerpsi(3,0.0), centerbeta(3,0.0), solbeta(nstates,0.0);;
	
	TPZFMatrix<REAL> A(dim,dim);  // Linear System matrix
	grad.Redim(dim,1);   
	
	// matrizes para aplicar o metodo dos minimos quadrados
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	
    // Encontrando o centro do elemento atual (cel)
    TPZGeoEl* gelalfa = cel->Reference();
    gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
    center.Fill(0.);
    gelalfa->X(centerpsi,center);
    cel->Solution(centerpsi,var,solalfa);
    
    neighs.Resize(0);
    
    // Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
    for(side = 0; side <cel->Reference()->NSides()-1; side++)
    {
        TPZCompElSide celside(cel,side);
        celside.ConnectedElementList(neighs,0,0);
    }
    
    // si no hay vecinos continuamos con el siguiente elemento
    nneighs = neighs.NElements();
    if(!nneighs) DebugStop();
    
    std::set<TPZCompEl *> neighscel;
    for(int i =0; i<nneighs; i++)
    {
        TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(neighs[i].Element());
        if(!InterpEl || InterpEl->Dimension()!=dim) continue;
        neighscel.insert(neighs[i].Element());
    }
    
    nneighs=neighscel.size();
    
    // si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente
    // Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
    // y el valor de la solucion en su centro solbeta
    DeltaH.Redim(nneighs,dim);
    DeltaHTranspose.Redim(dim,nneighs);
    DifSol.Redim(nneighs,1);
    
    // Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
    int ineighs=-1;
    int counter=0;
    std::set<TPZCompEl *>::iterator it;
    for(it=neighscel.begin(); it!=neighscel.end(); ++it)
    {
        //(*it)->Print();
        ineighs++;
        TPZGeoEl * gelbeta = (*it)->Reference();
        
        if(!gelbeta) DebugStop();
        
        centerpsi.Fill(0.0);
        centerbeta.Fill(0.0);
        gelbeta->CenterPoint(gelbeta->NSides()-1,centerpsi);
        gelbeta->X(centerpsi,centerbeta);
        gelbeta->Reference()->Solution(centerpsi,var,solbeta);
        
        for(k=0;k<dim;k++)
        {
            DeltaH(ineighs,k) = centerbeta[k] - center[k];
        }
        DifSol(ineighs,0) = solbeta[nstates-1] - solalfa[nstates-1];
        
        counter ++;
    }
     
    // Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    grad = DeltaHTranspose*DifSol;
    A = DeltaHTranspose*DeltaH;
    if(counter > 0){
        A.SolveDirect(grad,ELU);
    }
}

void GradReconstByLeastSquaresOnlyContEl(TPZCompEl *cel,TPZManVector<REAL,3> &center, TPZManVector<REAL> &solalfa,  TPZFMatrix<REAL> &grad, int var) {
    
    int dim;
    dim = cel->Mesh()->Dimension();
    
    
    // Nada sera realizado para elementos com dimensao diferente da dimensao do problema
    if(!cel || cel->Dimension()!=dim) DebugStop();
    
	int nstates=0;
    nstates = cel->Material()->NSolutionVariables(var);
    
    int k, side;
    
	TPZStack<TPZCompElSide> neighs;
	int nneighs;
	
    center.Resize(3, 0.);
    solalfa.Resize(nstates,0.0);
	TPZManVector<REAL> centerpsi(3,0.0), centerbeta(3,0.0), solbeta(nstates,0.0);;
	
	TPZFMatrix<REAL> A(dim,dim);    // Linear System matrix
	grad.Redim(dim,1);   // Linear System vector
	
	// matrizes para aplicar o metodo dos minimos quadrados
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	
    // Encontrando o centro do elemento atual (cel)
    TPZGeoEl* gelalfa = cel->Reference();
    gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
    center.Fill(0.);
    gelalfa->X(centerpsi,center);
    cel->Solution(centerpsi,var,solalfa);
    
    neighs.Resize(0);
    
    // Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
    for(side = 0; side < cel->Reference()->NSides()-1; side++) {
        TPZCompElSide celside(cel,side);
        celside.ConnectedElementList(neighs,0,0);
    }
    nneighs = neighs.NElements();
    
    // si no hay vecinos continuamos con el siguiente elemento
    if(!nneighs) DebugStop();
    
    // si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente
    // Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
    // y el valor de la solucion en su centro solbeta
    DeltaH.Redim(nneighs,dim);
    DeltaHTranspose.Redim(dim,nneighs);
    DifSol.Redim(nneighs,1);
    
    // Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
    int nsides = cel->Reference()->NSides()-1;
    
    // Para cada lado calculamos los deltaH (desde el centro del elemento al centro del lado de dimension menor a él
    // y el valor de la solucion en su centro solbeta
    DeltaH.Redim(nsides,dim);
    DeltaHTranspose.Redim(dim,nsides);
    DifSol.Redim(nsides,1);
    
    // Procuramos todos los puntos medios de cada lado del elemento y calculamos baseados en los valores de la solucion sobre ellos
    for(side = 0; side < nsides; side++)
    {
        centerpsi.Fill(0.0);
        centerbeta.Fill(0.0);
        cel->Reference()->CenterPoint(side,centerpsi);
        cel->Reference()->X(centerpsi,centerbeta);
        cel->Solution(centerpsi,var,solbeta);
        for(k=0;k<dim;k++)
        {
            DeltaH(side,k) = centerbeta[k] - center[k];
        }
        DifSol(side,0) = solbeta[nstates-1] - solalfa[nstates-1];
    }
    
    // Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
    DeltaH.Transpose(&DeltaHTranspose);
    grad = DeltaHTranspose*DifSol;
    A = DeltaHTranspose*DeltaH;
    A.SolveDirect(grad,ELU);
}


void ProjectionGradientReconstructedInFESpace(TPZCompMesh *cmesh,int var, int matid_l2proj){
    
    // Redimensionando a matriz dos dados da reconstruca de gradientes
    int dim  = cmesh->Dimension();
    int nelem = cmesh->NElements();
    
    TPZManVector<REAL,3> center;
    TPZManVector<REAL> solalfa;
    TPZFMatrix<REAL> Grad;
    
    //criar ponteiro para TPZFunction
    TPZGradient *pGrad = new TPZGradient;
    TPZAutoPointer<TPZFunction<STATE> > fp(pGrad);
    
    //Criar matrix de rigidez e vetor de carga
    int numloadcases;
	unsigned int im;
    for(im=0; im<cmesh->MaterialVec().size(); im++){
        if(!cmesh->MaterialVec()[im]) continue;
        numloadcases = cmesh->MaterialVec()[im]->NumLoadCases();
        break;
    }
    
    int neq = cmesh->NEquations();
    TPZFMatrix<STATE> rhs;
    rhs.Redim(neq,numloadcases);
    
    //TPZBandStructMatrix stmatrix(cmesh);
    TPZSkylineStructMatrix stmatrix(cmesh);
    TPZMatrix<STATE> *stiffmatrix = stmatrix.Create();
    
    int matid;
    
    for(int i=0; i<nelem; i++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZElementMatrix ek(cel->Mesh(), TPZElementMatrix::EK);
        TPZElementMatrix ef(cel->Mesh(), TPZElementMatrix::EF);
        
        // Nada sera realizado para elementos con dimension diferente de la dimension del problema
        if(cel->Dimension()!=dim) continue;
            
        matid = cel->Material()->Id();
        
        //gradient reconstruction
        GradientReconstructionByLeastSquares(cel, center, solalfa, Grad,var);
        
        //set data of the gradient reconstructed
        pGrad->SetData(center, Grad, solalfa[0]);
        
        //change material id current to material id of the L2ProjectionMaterial 
        ChangeMaterialIdIntoCompElement(cel, matid, matid_l2proj);
        
        //set forcing function of l2 projection material
        TPZMaterial *mat = cel->Material();
        mat->SetForcingFunction(fp);
        
        //load the matrix ek and vector ef of the element
        cel->CalcStiff(ek,ef);
        
        //assemble pos l2 projection
        AssembleGlobalMatrix(cel, ek, ef, *stiffmatrix, rhs);
        
        //Return for original material and current solution of the mesh
        ChangeMaterialIdIntoCompElement(cel, matid_l2proj, matid);
    }
    
        
//#ifdef LOG4CXX
//    if(logdata->isDebugEnabled())
//    {
//        std::stringstream sout;
//        stiffmatrix->Print("Matriz de Rigidez: ",sout,EMathematicaInput);
//        rhs.Print("Right Handside", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str())
//    }
//#endif
    
    //Solve linear system and transfer the solution to computational mesh
    TPZStepSolver<REAL> step;
    //step.SetDirect(ELU);
    step.SetDirect(ELDLt);
    step.SetMatrix(stiffmatrix);
    TPZFMatrix<REAL> result;
    step.Solve(rhs, result);
    cmesh->LoadSolution(result);
    
    //stiffmatrix->SolveDirect(rhs,ELU);
    //cmesh->LoadSolution(rhs);

}

void AssembleGlobalMatrix(TPZCompEl *el, TPZElementMatrix &ek, TPZElementMatrix &ef,TPZMatrix<STATE> & stiffmatrix, TPZFMatrix<STATE> &rhs){
    
//    
//#ifdef CHECKCONSISTENCY
//    extern TPZCheckConsistency stiffconsist("ElementStiff");
//    stiffconsist.SetOverWrite(true);
//    bool result;
//    result = stiffconsist.CheckObject(ek.fMat);
//    if(!result)
//    {
//        globalresult = false;
//        std::stringstream sout;
//        sout << "element " << iel << " computed differently";
//        LOGPZ_ERROR(loggerCheck,sout.str())
//    }
//    
//#endif
    
    
    if(!el->HasDependency()) {
        ek.ComputeDestinationIndices();
        
        stiffmatrix.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                
#ifdef LOG4CXX
        if(logdata->isDebugEnabled())
        {
            std::stringstream sout;
            ek.Print(sout);
            ef.Print(sout);
            LOGPZ_DEBUG(logdata,sout.str())
        }
#endif
    } else {
        // the element has dependent nodes
        ek.ApplyConstraints();
        ef.ApplyConstraints();
        ek.ComputeDestinationIndices();
        
        stiffmatrix.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        
//#ifdef LOG4CXX
//        if(logdata->isDebugEnabled())
//        {
//            std::stringstream sout;
//            ek.Print(sout);
//            ef.Print(sout);
//            LOGPZ_DEBUG(logdata,sout.str())
//        }
//#endif
    }
}

void PosProcessGradientReconstruction(TPZCompMesh *cmesh,int var,TPZFMatrix<REAL> &datagradients){
    
    // Redimensionando a matriz dos dados da reconstruca de gradientes
    int dim  = cmesh->Dimension();
    int nelem = cmesh->NElements();
    datagradients.Redim(nelem,3*dim);
    
    TPZManVector<REAL,3> center;
    TPZManVector<REAL> solalfa;
    TPZFMatrix<REAL> Grad;
    
	int i, k;
	int counter = 0;

   
    TPZCompEl *cel;
    // Calculando el gradiente por elemento computacional
    for(i=0;i<nelem;i++) {
        
        cel = cmesh->ElementVec()[i];
        
        // Nada sera realizado para elementos con dimension diferente de la dimension del problema
        if(cel->Dimension()!=dim) continue;
        
        GradientReconstructionByLeastSquares(cel, center, solalfa, Grad,var);
		
        //data of the vector gradiente
        for(k=0;k<dim;k++){
            if(!k) datagradients(counter,0) = cel->Index();//Id do elemento
            datagradients(counter,dim+k) = center[k];//centro do elemento
            if(!IsZero(Grad(k,0))) datagradients(counter,2*dim+k) = Grad(k,0);//valor do gradiente
        }
        
		counter++;
    	}
    // Redimensionando la matriz de los gradientes
    datagradients.Resize(counter,3*dim);    
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

void PosProcessamento(TPZAnalysis &an, TPZCompMesh *Cmesh, std::string plotfile){
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
	scalnames[0] = "Solution";
    scalnames[1] = "ExactPressure";
    vecnames[0] = "Derivative";
    vecnames[1] = "ExactFlux";
    
	const int dim = Cmesh->Dimension();
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
}

TPZFMatrix<REAL> MatrixR(REAL ang)
{
	TPZFMatrix<REAL> r(2,2,0.);
	r(0,0) = cos(ang); r(0,1) = -sin(ang);
	r(1,0) = sin(ang);	r(1,1) =  cos(ang);
    
	return r;
}

#include "pzgeoelbc.h"

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
    mA.Put(0,1,origX+2.); mA.Put(1,1,origY);
    mA.Put(0,2,origX+2.); mA.Put(1,2,origY+2.);
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
    
    for ( int ref = 0; ref < nh; ref++ ) {
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ) {
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			if (gel->Dimension() == 2) gel->Divide (filhos);
            gel->Divide (filhos);
		}//for i
	}//ref
	
	return gmesh;
}

TPZGeoMesh *GMesh2(){
    
    int Qnodes = 6;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolLine(2);
    TPZVec <int> TopolPoint(1);
	
	//indice dos nos
    int id;
	id = 0;
    REAL dx = 1.;
	for (int i=0; i<Qnodes/2;i++) {
		
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0, i*dx);//coord x
		Node[id].SetCoord(1, 0.);//coord y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	for (int i=0; i<Qnodes/2;i++) {
		
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0, 2. - i*dx);//coord x
		Node[id].SetCoord(1, 2.);//coord y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	//indice dos elementos
	id = 0;
    
    for(int i = 0; i<0.5*Qnodes-1; i++){
        TopolQuad[0] = i;
        TopolQuad[1] = i+1;
        TopolQuad[2] = (Qnodes-2)-i;
        TopolQuad[3] = (Qnodes-1)-i;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
    }
	gmesh->BuildConnectivity();
    
	// Boundary elements
	TPZGeoElBC gbc10(gmesh->ElementVec()[0],4,bc0);
	TPZGeoElBC gbc11(gmesh->ElementVec()[0],6,bc2);
	TPZGeoElBC gbc12(gmesh->ElementVec()[0],7,bc3);
	TPZGeoElBC gbc13(gmesh->ElementVec()[1],4,bc0);
	TPZGeoElBC gbc14(gmesh->ElementVec()[1],5,bc1);
	TPZGeoElBC gbc15(gmesh->ElementVec()[1],6,bc2);

	gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
	
	return gmesh;
}
    /*------------------ fazer refinamento direcional ------------------------------------
	 TopolLine[0] = 1;
	 TopolLine[1] = 4;
	 new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
	 id++;
	 
	 //    TopolPoint[0] = 2;
	 //    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,matId+1,*gmesh);
	 //    id++;
	 //    
	 //    TopolPoint[0] = 7;
	 //    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,matId+1,*gmesh);
	 //    id++;
	 
	 
	 for(int i = 0; i<Qnodes/2-1; i++){
	 TopolLine[0] = i;
	 TopolLine[1] = i+1;
	 new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	 id++;
	 }
	 
	 TopolLine[0] = (Qnodes/2)-1;
	 TopolLine[1] = Qnodes/2;
	 new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
	 id++;
	 
	 for(int i = 0; i<Qnodes/2-1; i++){
	 TopolLine[0] = (Qnodes/2) + i;
	 TopolLine[1] = (Qnodes/2+1) + i;
	 new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
	 id++;
	 }
	 
	 TopolLine[0] = Qnodes-1;
	 TopolLine[1] = 0;
	 new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
	 
	 gmesh->BuildConnectivity();

	 set<int> SETmatRef;
    for(int j = 0; j < nrefdir; j++)
    {
        int nel = gmesh->NElements();
        for (int iref = 0; iref < nel; iref++)
        {
            TPZVec<TPZGeoEl*> filhos;
            TPZGeoEl * gelref = gmesh->ElementVec()[iref];
            if(!gelref) continue;
            SETmatRef.insert(matId + 1);
            //TPZRefPatternTools::RefineDirectional(gelref, SETmatRef);
            TPZRefPatternTools::RefineUniformIfNeighMat(gelref, SETmatRef);
        }		
    }
    
    //refinamento uniforme
    for ( int ref = 0; ref < nh; ref++ ) {
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ) {
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			if (gel->Dimension() == 2) gel->Divide (filhos);
            gel->Divide (filhos);
		}//for i
	}//ref
	
	return gmesh;
}*/


TPZCompMesh *CMesh2(TPZGeoMesh *gmesh, int pOrder,bool isdiscontinuous)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    
    REAL diff=1.;
    REAL conv =0.;
    TPZVec<REAL> convdir;
    convdir.Resize(dim,0.);
    convdir[0]=0.;
    
    material-> SetParameters(diff, conv, convdir);
    material->SetNoPenalty();
    material->SetNonSymmetric();
    
    TPZAutoPointer<TPZFunction<STATE> > myforce = new TPZDummyFunction<STATE>(ForcingF);
    material->SetForcingFunction(myforce);
    
    TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(SolucaoExata);
	material->SetForcingFunctionExact(solExata);
    
    TPZMaterial * mat(material);
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    //TPZAutoPointer<TPZFunction<STATE> > fCC0 = new TPZDummyFunction<STATE>(Forcingbc);
    //TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(Forcingbc);
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    
    val2(0,0) =  0.0886226925452758;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    
    val2(0,0) =  -0.0886226925452758;
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    //BCond0->SetForcingFunction(fCC0);
    //BCond2->SetForcingFunction(fCC2);
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,material->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
    
    if (isdiscontinuous==true)
    {
        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
        for(int el = 0; el < cmesh->ElementVec().NElements(); el++)
        {
            TPZCompEl * compEl = cmesh->ElementVec()[el];
            if(!compEl) continue;
            int index = compEl ->Index();
            if(compEl->Dimension() == cmesh->Dimension())
            {
                TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[index]);
                if(!InterpEl) continue;
                InterpEl->CreateInterfaces(false);
            }
        }
    }
    
    else {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    }
	
    return cmesh;
}

#include "pzbuildmultiphysicsmesh.h"
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder,bool isdiscontinuous)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    
    REAL diff=1.;
    REAL conv =0.;
    TPZVec<REAL> convdir;
    convdir.Resize(dim,0.);
    convdir[0]=0.;
    
    material-> SetParameters(diff, conv, convdir);
    material->SetNoPenalty();
    material->SetNonSymmetric();
   // material->SetPenaltyConstant(1.);
    
    REAL ff=0.;
    material->SetInternalFlux(ff);
    
    TPZMaterial * mat(material);
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    TPZAutoPointer<TPZFunction<STATE> > fCC0 = new TPZDummyFunction<STATE>(Forcingbc0);
    TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(Forcingbc1);
    TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(Forcingbc2);
    TPZAutoPointer<TPZFunction<STATE> > fCC3 = new TPZDummyFunction<STATE>(Forcingbc3);
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    val2(0,0)=0.;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    val2(0,0)=2.;
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    val2(0,0)=2.;
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
//    
    BCond0->SetForcingFunction(fCC0);
    BCond1->SetForcingFunction(fCC1);
    BCond2->SetForcingFunction(fCC2);
    BCond3->SetForcingFunction(fCC3);
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,material->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
    
    if (isdiscontinuous==true)
    {
        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
        for(int el = 0; el < cmesh->ElementVec().NElements(); el++)
        {
            TPZCompEl * compEl = cmesh->ElementVec()[el];
            if(!compEl) continue;
            int index = compEl ->Index();
            if(compEl->Dimension() == cmesh->Dimension())
            {
                TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[index]);
                if(!InterpEl) continue;
                InterpEl->CreateInterfaces(false);
            }
        }
    }
    
    else {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    }
	
    return cmesh;
}

void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]=  coef_a*x + coef_b*y;
    
}

void Forcingbc(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
    
    double x = pt[0];
#ifdef USING_BOOST
    disp[0]=0.05*sqrt(M_PI)*(boost::math::erf(10.*(x-1.)));
#else
    disp[0]=0.05*sqrt(M_PI)*erf(10.*(x-1.));
#endif
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

void ForcingF(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
    
    double x = pt[0];
    double aux1 = -100.*(1 - 2.*x + x*x);
    double aux2 = exp(aux1);
    disp[0]=-(200.*x - 200.)*aux2;
}

void SolucaoExata(const TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix<REAL> &deriv){
    
    deriv(0,0)=0.;
    deriv(1,0)=0.;
    sol[0]=0.;
    
    double x = pt[0];
#ifdef USING_BOOST
    sol[0] = 0.05*sqrt(M_PI)*(boost::math::erf(10.*(x-1.)));
#else
    sol[0] = 0.05*sqrt(M_PI)*erf(10.*(x-1.));
#endif
    REAL val = -100.*(x-1.)*(x-1.);
    deriv(0,0)=exp(val);
}

#include "TPZSkylineNSymStructMatrix.h"
void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	TPZBandStructMatrix full(Cmesh); //caso nao-simetrico
	//TPZSkylineStructMatrix full(Cmesh); //caso simetrico
    //TPZSkylineNSymStructMatrix full(Cmesh);//caso nao-simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	//step.SetDirect(ELDLt); //caso simetrico
	step.SetDirect(ELU);//caso nao simetrico
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

