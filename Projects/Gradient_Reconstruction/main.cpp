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

#include "TPZGenGrid2D.h"
#include "pzgradient.h"
#include "pzl2projection.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "pzgeoelbc.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "TPZReadGIDGrid.h"
#include "TPZExtendGridDimension.h"


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

// Alfa -> Coefficient of the arctang argument
REAL ALFA = 50.;

REAL ValueK = 1000;

void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &points,REAL r,REAL &distance,bool &isdefined);

void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);
void RefiningNearLine(int dim,TPZGeoMesh *gmesh,int nref);


TPZFMatrix<REAL> MatrixR(REAL ang);
TPZGeoMesh *GMesh(int triang_elements, REAL angle, REAL origX, REAL origY, int nh);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, bool isdiscontinuous);

TPZGeoMesh *GMesh2();
TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int pOrder,bool isdiscontinuous);

void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Forcingbc1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Forcingbc2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Forcingbc3(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Forcingbc(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

TPZGeoMesh *CreateGeoMesh(MElementType typeel);
TPZGeoMesh *CreateGeoMesh(std::string &nome);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);

REAL MeanCell(TPZCompEl *cel,int IntOrder);
REAL ExactSolution(int dim,const TPZVec<REAL> &pt);
void GradExactSolution(int dim,const TPZVec<REAL> &x, TPZVec<REAL> &dsol);

void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

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
void GradientReconstructionByLeastSquares(TPZCompEl *cel,TPZManVector<REAL,3> &center,TPZVec<REAL> &Grad);

// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given
void PrintDataMeshVTK(TPZCompMesh *cmesh, char *filename,TPZFMatrix<REAL> &elData);
void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

void GradReconstByLeastSquaresOnlyContEl(TPZCompEl *cel,TPZManVector<REAL,3> &center, TPZManVector<REAL> &solalfa,  TPZFMatrix<REAL> &grad, int var);
void PosProcessGradientReconstruction(TPZCompMesh *cmesh,TPZFMatrix<REAL> &datagradients);

void AssembleGlobalMatrix(TPZCompEl *el, TPZElementMatrix &ek, TPZElementMatrix &ef,TPZMatrix<STATE> & stiffmatrix, TPZFMatrix<STATE> &rhs);

void SaidaMathGradiente(TPZFMatrix<REAL> gradients,int ref,int typeel,ofstream &out,bool uniform = true);

//Trocar todos os elementos do cmesh apontando para o material TPZL2ProjectionFromGradient
void ChangeMaterialIdIntoCompElement(TPZCompEl *cel, int oldmatid, int newmatid);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

int MaxRefs = 6;
int InitRefs = 5;

int main(int argc, char *argv[]) {

	
	// Initializing uniform refinements for quadrilaterals and triangles
    //gRefDBase.InitializeRefPatterns();
	gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);

    int p = 1;
    char saida[256];
	int dim = 0;
	MElementType typeel = EOned;
	int nrefs;
	
	// One dimensional case - Uniform refinement on regular mesh
	// Output file in Mathematica format
	sprintf(saida,"Grad2MathUnifMeshRUnif1D.nb");
	ofstream outfilemath(saida);
	TPZGeoMesh *gmesh;
	TPZCompMesh *cmesh;
	for(nrefs=InitRefs;nrefs<MaxRefs;nrefs++)
	{
		// geometric mesh (initial)
		gmesh = CreateGeoMesh(typeel);
		// Refining near the points belong a circunference with radio r - maxime distance radius
		UniformRefine(gmesh,nrefs);
		
		// computational mesh
		cmesh= CreateCMesh(gmesh,p,true);
		
		// Computing gradient reconstructed
		TPZFMatrix<REAL> gradients;
		PosProcessGradientReconstruction(cmesh,gradients);
		
		// Printing to VTK
		sprintf(saida,"Grad_UnifMeshRUnif_H%d_E%d.vtk",nrefs,typeel);
		PrintDataMeshVTK(cmesh,saida,gradients);
		// Printing to Mathematica
		SaidaMathGradiente(gradients,nrefs,typeel,outfilemath);
		
		cmesh->CleanUp();
		delete cmesh;
		delete gmesh;
	}
	outfilemath.close();
	// One dimensional case - refinement localized on regular mesh
	// Output file in Mathematica format
	sprintf(saida,"Grad2MathUnifMeshRNoUnif1D.nb");
	outfilemath.open(saida);
	for(nrefs=InitRefs;nrefs<MaxRefs;nrefs++)
	{
		// geometric mesh (initial)
		gmesh = CreateGeoMesh(typeel);
		// Refining near the points belong a circunference with radio r - maxime distance radius
		RefiningNearCircunference(dim,gmesh,nrefs,1);
		
		// computational mesh
		cmesh= CreateCMesh(gmesh,p,true);
		
		// Computing gradient reconstructed
		TPZFMatrix<REAL> gradients;
		PosProcessGradientReconstruction(cmesh,gradients);
		
		// Printing to VTK
		sprintf(saida,"Grad_UnifMeshRUnif_H%d_E%d.vtk",nrefs,typeel);
		PrintDataMeshVTK(cmesh,saida,gradients);
		// Printing to Mathematica
		SaidaMathGradiente(gradients,nrefs,typeel,outfilemath);
		
		cmesh->CleanUp();
		delete cmesh;
		delete gmesh;
	}
	outfilemath.close();
	
	// Utilizando para um refinamento localizado sobre uma malha uniforme
	// Output file in Mathematica format
	sprintf(saida,"Grad2MathUnifMeshRUnif2D.nb");
	outfilemath.open(saida);
	// Bi-dimensional case: triangles and quadrilaterals, on regular mesh and uniform refinement
	while(1) {
		if(typeel==ETriangle) typeel = EQuadrilateral;
		if(typeel==EOned) typeel = ETriangle;
		for(nrefs=InitRefs;nrefs<MaxRefs;nrefs++)
		{
			// geometric mesh (initial)
			gmesh = CreateGeoMesh(typeel);
			// Refining near the points belong a circunference with radio r - maxime distance radius
			UniformRefine(gmesh,nrefs);
			// computational mesh
			cmesh= CreateCMesh(gmesh,p,true);
			
			// Computing gradient reconstructed
			TPZFMatrix<REAL> gradients;
			PosProcessGradientReconstruction(cmesh,gradients);
			
			// Printing to VTK
			sprintf(saida,"Grad_UnifMeshRUnif_H%d_E%d.vtk",nrefs,typeel);
			PrintDataMeshVTK(cmesh,saida,gradients);
			// Printing to Mathematica
			SaidaMathGradiente(gradients,nrefs,typeel,outfilemath);
			cmesh->CleanUp();
			delete cmesh;
			delete gmesh;
		}
		if(typeel==EQuadrilateral) break;
	}
	outfilemath.close();

	// Utilizando para um refinamento localizado sobre uma malha uniforme
	// Output file in Mathematica format
	sprintf(saida,"Grad2MathUnifMeshRNoUnif2D.nb");
	outfilemath.open(saida);
	while(1) {
		if(typeel==ETriangle) typeel = EQuadrilateral;
		if(typeel==EOned) typeel = ETriangle;
		for(nrefs=InitRefs;nrefs<MaxRefs;nrefs++)
		{
			// geometric mesh (initial)
			gmesh = CreateGeoMesh(typeel);
			// Refining near the points belong a circunference with radio r - maxime distance radius
			RefiningNearCircunference(dim,gmesh,nrefs,1);
			// First computational mesh
			cmesh= CreateCMesh(gmesh,p,true);
			
			// Computing gradient reconstructed
			TPZFMatrix<REAL> gradients;
			PosProcessGradientReconstruction(cmesh,gradients);
			// Printing to VTK
			sprintf(saida,"Grad_UnifMeshRNoUnif_H%d_E%d.vtk",nrefs,typeel);
			PrintDataMeshVTK(cmesh,saida,gradients);
			// Printing to Mathematica
			SaidaMathGradiente(gradients,nrefs,typeel,outfilemath,false);
			
			cmesh->CleanUp();
			delete cmesh;
			delete gmesh;
		}
		if(typeel==EQuadrilateral) break;
	}
	outfilemath.close();

	// Utilizando para um refinamento uniforme sobre uma malha gerada pelo GID
	// Output file in Mathematica format
	sprintf(saida,"Grad2Math_GIDRUnif2D.nb");
	outfilemath.open(saida);
	
	while(1) {
		if(typeel==ETriangle) typeel = EQuadrilateral;
		if(typeel==EOned) typeel = ETriangle;
		for(nrefs=InitRefs;nrefs<MaxRefs;nrefs++)
		{
			// geometric mesh (initial)
			std::string nombre = PZSOURCEDIR;
			if(!typeel)
				nombre += "/Projects/Gradient_Reconstruction/RegionQuadrada.dump";
			else 
				nombre += "/Projects/Gradient_Reconstruction/RegionQuadradaT.dump";
			gmesh = CreateGeoMesh(nombre);
			if(!gmesh) break;
			// Uniform refinements
			UniformRefine(gmesh,nrefs);
			
			// First computational mesh
			cmesh= CreateCMesh(gmesh,p,true);
			
			// Computing gradient reconstructed
			TPZFMatrix<REAL> gradients;
			PosProcessGradientReconstruction(cmesh,gradients);
			
			// Printing to VTK
			sprintf(saida,"Grad_GIDRUnif_H%d_E%d.vtk",nrefs,typeel);
			PrintDataMeshVTK(cmesh,saida,gradients);
			// Printing to Mathematica
			SaidaMathGradiente(gradients,nrefs,typeel,outfilemath);
			
			cmesh->CleanUp();
			delete cmesh;
			delete gmesh;
			
		}
		if(typeel==EQuadrilateral) break;
	}
	outfilemath.close();
	
	// Utilizando para um refinamento localizado sobre uma malha gerada pelo GID
	// Output file in Mathematica format
	sprintf(saida,"Grad2Math_GIDRNoUnif2D.nb");
	outfilemath.open(saida);
	
	while(1) {
		if(typeel==ETriangle) typeel = EQuadrilateral;
		if(typeel==EOned) typeel = ETriangle;
		for(nrefs=InitRefs;nrefs<MaxRefs;nrefs++)
		{
			// geometric mesh (initial)
			std::string nombre = PZSOURCEDIR;
			if(!typeel)
				nombre += "/Projects/Gradient_Reconstruction/RegionQuadrada.dump";
			else 
				nombre += "/Projects/Gradient_Reconstruction/RegionQuadradaT.dump";
			gmesh = CreateGeoMesh(nombre);
			if(!gmesh) break;
			// Refining near the points belong a circunference with radio r - maxime distance radius
			RefiningNearCircunference(dim,gmesh,nrefs,1);    
			// First computational mesh
			cmesh= CreateCMesh(gmesh,p,true);
			
			// Computing gradient reconstructed
			TPZFMatrix<REAL> gradients;
			PosProcessGradientReconstruction(cmesh,gradients);
			
			// Printing to VTK
			sprintf(saida,"Grad_GIDRNoUnif_H%d_E%d.vtk",nrefs,typeel);
			PrintDataMeshVTK(cmesh,saida,gradients);
			// Printing to Mathematica
			SaidaMathGradiente(gradients,nrefs,typeel,outfilemath,false);
			
			cmesh->CleanUp();
			delete cmesh;
			delete gmesh;
			
		}
		if(typeel==EQuadrilateral) break;
	}
	outfilemath.close();

	// Three dimensional case
	sprintf(saida,"Grad2Math_3D_RUnif.nb");
	outfilemath.open(saida);
	while(1) {
		typeel = ECube;
		for(nrefs=InitRefs;nrefs<MaxRefs;nrefs++)
		{
			// geometric mesh (initial)
			gmesh = CreateGeoMesh(typeel);
			// Refining near the points belong a circunference with radio r - maxime distance radius
			UniformRefine(gmesh,nrefs);
			
			// computational mesh
			cmesh= CreateCMesh(gmesh,p,true);
			
			// Computing gradient reconstructed
			TPZFMatrix<REAL> gradients;
			PosProcessGradientReconstruction(cmesh,gradients);
			
			// Printing to VTK
			sprintf(saida,"Grad3D_UnifMeshRUnif_H%d_E%d.vtk",nrefs,typeel);
			PrintDataMeshVTK(cmesh,saida,gradients);
			// Printing to Mathematica
			SaidaMathGradiente(gradients,nrefs,typeel,outfilemath);
			
			cmesh->CleanUp();
			delete cmesh;
			delete gmesh;
			
		}
		if(typeel==ECube) break;
	}
	outfilemath.close();
	
	// Three dimensional case
	sprintf(saida,"Grad2Math_3D_RNoUnif.nb");
	outfilemath.open(saida);
	while(1) {
		typeel = ECube;
		for(nrefs=InitRefs;nrefs<MaxRefs;nrefs++)
		{
			// geometric mesh (initial)
			gmesh = CreateGeoMesh(typeel);
			// Refining near the points belong a circunference with radio r - maxime distance radius
			RefiningNearCircunference(dim,gmesh,nrefs,1);
			// computational mesh
			cmesh= CreateCMesh(gmesh,p,true);
			
			// Computing gradient reconstructed
			TPZFMatrix<REAL> gradients;
			PosProcessGradientReconstruction(cmesh,gradients);
			
			// Printing to VTK
			sprintf(saida,"Grad3D_UnifMeshRUnif_H%d_E%d.vtk",nrefs,typeel);
			PrintDataMeshVTK(cmesh,saida,gradients);
			// Printing to Mathematica
			SaidaMathGradiente(gradients,nrefs,typeel,outfilemath);
			
			cmesh->CleanUp();
			delete cmesh;
			delete gmesh;
			
		}
		if(typeel==ECube) break;
	}
	outfilemath.close();
	
	return 0;
}

REAL ExactSolution(int dim,const TPZVec<REAL> &x) {
	REAL result = 0.;
	if(dim == 1) {
		TPZVec<REAL> C0(1,-0.25);
		
		REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]));
		result = atan(ALFA * ( R0 - 0.75) );
	}
	else if(dim==2) {
		TPZVec<REAL> C0(2,-0.25);
		
		REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]));
		result = atan(ALFA * ( R0 - 0.75) );
	}
	else if(dim == 3) {
		TPZVec<REAL> C0(3,-0.25);
		
		REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
		result = atan(ALFA * ( R0 - 0.75) );
	}
	else {
		REAL F = 2*sqrt(ValueK);
		REAL arc = F*((0.25*0.25) - (x[0] - 0.5)*(x[0] - 0.5) - (x[1] - 0.5)*(x[1] - 0.5));
		REAL prodx = x[0]*(x[0]-1.);
		REAL prody = x[1]*(x[1]-1.);
		REAL prod = prodx*prody;
		result = (8*prod*(1+(2./M_PI)*(atan(arc))));
	}
	return result;
}
void GradExactSolution(int dim,const TPZVec<REAL> &x, TPZVec<REAL> &dsol) {
	if(dim == 1) {
		TPZVec<REAL> C0(1,-0.25);
		REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]));
		REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(2.))*(R0-sqrt(2.)));
		if(IsZero(den))
			DebugStop();
		dsol[0] = (ALFA*(x[0]-C0[0]))/den;
	}
	else if(dim==2) {
		TPZVec<REAL> C0(2,-0.25);
		REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]));
		REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(2.))*(R0-sqrt(2.)));
		if(IsZero(den))
			DebugStop();
		dsol[0] = (ALFA*(x[0]-C0[0]))/den;
		dsol[1] = (ALFA*(x[1]-C0[1]))/den;
	}
	else if(dim==3) {
		TPZVec<REAL> C0(3,-0.25);
		REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
		REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(2.))*(R0-sqrt(2.)));
		if(IsZero(den))
			DebugStop();
		dsol[0] = (ALFA*(x[0]-C0[0]))/den;
		dsol[1] = (ALFA*(x[1]-C0[1]))/den;
		dsol[2] = (ALFA*(x[2]-C0[2]))/den;
	}
	else {
		REAL F = 2*sqrt(ValueK);
		REAL arc = F*((0.25*0.25) - (x[0] - 0.5)*(x[0] - 0.5) - (x[1] - 0.5)*(x[1] - 0.5));
		REAL prodx = x[0]*(x[0]-1.);
		REAL prody = x[1]*(x[1]-1.);
		REAL prod = prodx*prody;
		REAL temp = prody*(2*x[0]-1.)*(M_PI + 2*atan(arc));
		REAL frac = 2*prod*F*(1.-2*x[0]);
		frac = frac/(1+arc*arc);
		dsol[0] = (8./M_PI)*(temp + frac);
		temp = prodx*(2*x[1]-1.)*(M_PI + 2*atan(arc));
		frac = 2*prod*F*(1.-2*x[1]);
		frac = frac/(1+arc*arc);
		dsol[1] = (8./ M_PI)*(temp + frac);
	}
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {    
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
	// Re-constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}

void SaidaMathGradiente(TPZFMatrix<REAL> gradients,int nref,int typeel,ofstream &outfile,bool uniform) {
	int i, j;
    int dim = (gradients.Cols()-2)/2;
	TPZManVector<REAL> x(3), xunit(3), xorth(3);
	
	TPZFMatrix<REAL> Grads(gradients.Rows(),7);

	REAL r, Grad, temp;
	TPZVec<REAL> center(dim);
	TPZVec<REAL> GradExact(dim,0.0);
	REAL error = 0.0;
    for(i=0;i<gradients.Rows();i++)
    {
		r = temp = 0.0;
		if(dim==2) {
			for(j=0;j<dim;j++) {
				center[j] = gradients(i,j);
				x[j] = gradients(i,j)-0.5;
				r += x[j]*x[j];
				temp += gradients(i,dim+j)*gradients(i,dim+j);
			}
		}
		else if(dim==3) {
			for(j=0;j<dim;j++) {
				center[j] = gradients(i,j);
				x[j] = gradients(i,j)+0.25;
				r += x[j]*x[j];
				temp += gradients(i,dim+j)*gradients(i,dim+j);
			}
		}
		// Calculando o gradiente exato para cada centro
		GradExactSolution(dim,center,GradExact);
		Grad = sqrt(temp);
		temp = r;
		r = sqrt(temp);
		for(j=0;j<dim;j++) {
			if(IsZero(r)) 
				continue;
			xunit[j] = x[j]/r;
		}
		// Vetor ortogonal
		if(dim==2) {
			xorth[0] = -xunit[1];
			xorth[1] = xunit[0];
		}
		else if(dim==3) {
			xorth = xunit;
		}
        Grads(i,0) = r;
		Grads(i,1) = Grad;
		Grads(i,2) = 0.0;
		for(j=0;j<dim;j++) 
			Grads(i,2) += xunit[j]*gradients(i,dim+j);
		Grads(i,3) = 0.0;
		for(j=0;j<dim;j++) 
			Grads(i,3) = xorth[j]*gradients(i,dim+j);
		temp = 0.0;
		for(j=0;j<dim;j++) 
			temp += (GradExact[j]*GradExact[j]);
		Grads(i,4) = sqrt(temp);
		Grads(i,5) = 0.0;
		for(j=0;j<dim;j++) 
			Grads(i,5) += xunit[j]*GradExact[j];
		Grads(i,6) = 0.0;
		for(j=0;j<dim;j++) 
			Grads(i,6) += xorth[j]*GradExact[j];
		error += ((Grads(i,4)-Grad)*(Grads(i,4)-Grad)*gradients(i,2*dim+1));
    }
	
	// Printing in mathematica format
	char name[256];
	int countbyrow = 0;

	sprintf(name,"GradNormH%dE%d",nref,typeel);
	outfile << "(*H-Refinement = " << nref << "  Type of element = " << typeel << " *)" << std::endl;
    outfile << name << " = {";
    for(i=0;i<gradients.Rows();i++)
    {
        outfile << "{" << Grads(i,0) << ", " << Grads(i,1) << ", " << Grads(i,2) << ", " << Grads(i,3) << ", " << Grads(i,4) << ", " << Grads(i,5) << ", " << Grads(i,6) << "}";
		countbyrow++;
        if(i != gradients.Rows()-1) {
			if(countbyrow==3) {
				outfile << ", " << endl;
				countbyrow = 0;
			}
			else
				outfile << ", ";
		}
        if(i == gradients.Rows()-1) outfile << "};" << std::endl;
    }
	// Creating lists of the results Norm(GradR), GradRU.UnitVec, GradRU.OrthVec, Norm(Grad), GradU.UnitVec and GradU.OrthVec
    outfile << "NormGradRU" << nref << typeel << " = Table[{"<< name << "[[i,1]],"<< name <<"[[i,2]]},{i,1,Length["<<name<<"]}];" << endl;
    outfile << "GradRUPointV" << nref << typeel << " = Table[{"<< name << "[[i,1]],"<< name <<"[[i,3]]},{i,1,Length["<<name<<"]}];" << endl;
    outfile << "GradRUPointVOrth" << nref << typeel << " = Table[{"<< name << "[[i,1]],"<< name <<"[[i,4]]},{i,1,Length["<<name<<"]}];" << endl;
	outfile << "Error" << nref << "E" << typeel << " = " << sqrt(error) << ";" << endl;
	
	if(nref==MaxRefs-1) {
		// List plot of the Norm of gradients incremented with norm gradient exact
		outfile << "NormGradUExact" << typeel << " = Table[{"<< name << "[[i,1]],"<< name <<"[[i,5]]},{i,1,Length["<<name<<"]}];" << endl;
		outfile << "ListPlot[{";
		for(i=InitRefs-1;i<nref;i++) {
			outfile << "NormGradRU" << i+1 << typeel;
			if(i!=nref-1) outfile << ",";
			else outfile << ",NormGradUExact" << typeel << "}";
		}
		outfile << ",DataRange-> {";
		for(i=InitRefs-1;i<nref;i++) {
			outfile << i+2 << ",";
			if(i==nref-1) outfile << 4*(i+2) << "}";
		}
		outfile << ",PlotRange->All,Frame->True,AxesLabel->{\"r\",\"||Grad||\"}]" << endl << endl;
		// List plot of the scalar product gradR with unitV incremented with scalar product gradient exact with unitV
		outfile << "GradUPointVExact" << typeel << " = Table[{"<< name << "[[i,1]],"<< name <<"[[i,6]]},{i,1,Length["<<name<<"]}];" << endl;
		outfile << "ListPlot[{";
		for(i=InitRefs-1;i<nref;i++) {
			outfile << "GradRUPointV" << i+1 << typeel;
			if(i!=nref-1) outfile << ",";
			else outfile << ",GradUPointVExact" << typeel << "}";
		}
		outfile << ",DataRange-> {";
		for(i=InitRefs-1;i<nref;i++) {
			outfile << i+2 << ",";
			if(i==nref-1) outfile << 4*(i+2) << "}";
		}
		outfile << ",PlotRange->All,Frame->True,AxesLabel->{\"r\",\"GradU.Vunit\"}]" << endl << endl;
		// List plot of the scalar product gradR with V orthogonal incremented with scalar product gradient exact with V orthogonal
		outfile << "GradUPointVOrthExact" << typeel << " = Table[{"<< name << "[[i,1]],"<< name <<"[[i,7]]},{i,1,Length["<<name<<"]}];" << endl;
		outfile << "ListPlot[{";
		for(i=InitRefs-1;i<nref;i++) {
			outfile << "GradRUPointVOrth" << i+1 << typeel;
			if(i!=nref-1) outfile << ",";
			else outfile << ",GradUPointVOrthExact" << typeel << "}";
		}
		outfile << ",DataRange-> {";
		for(i=InitRefs-1;i<nref;i++) {
			outfile << i+2 << ",";
			if(i==nref-1) outfile << 4*(i+2) << "}";
		}
		outfile << ",PlotRange->All,Frame->True,AxesLabel->{\"r\",\"GradU.Vorth\"}]" << endl << endl;
		// Vector of the errors on number of elements
		outfile << "NElementsE" << typeel << " = {";
		for(int k=InitRefs-1; k < nref; k++) {
			outfile << "Length[NormGradRU" << k+1 << typeel;
			if(k!=nref-1) outfile << "],";
			else outfile << "]}";
		}
		outfile << endl;
		outfile << "ErrorsE" << typeel << " = {";
		for(int k =InitRefs-1; k < nref; k++) {
			outfile << "Error" << k+1 << "E" << typeel;
			if(k!=nref-1) outfile << ",";
			else outfile << "}";
		}
		outfile << endl;
		if(uniform) {
			outfile << "ListPlot[Table[{1./Power[NElementsE" << typeel << "[[i]],1./" << dim << "],ErrorsE" << typeel << "[[i]]},{i,1,Length[NElementsE" << typeel;
			outfile << "]}],Joined->True,AxesLabel->{\"h\",\"Error\"},PlotMarkers->Automatic,AxesOrigin ->{0.0,0.0}]" << endl << endl;
		}
		else {
			outfile << "ListPlot[Table[{Log[NElementsE" << typeel << "[[i]]],Log[ErrorsE" << typeel << "[[i]]]},{i,1,Length[NElementsE" << typeel;
			outfile << "]}],Joined->True,AxesLabel->{\"h\",\"Error\"},PlotMarkers->Automatic,AxesOrigin ->{0.0,0.0}]" << endl << endl;
		}
	}
	outfile << endl;
}

TPZGeoMesh *CreateGeoMesh(std::string &archivo) {
	
	// Ejemplo uni-dimensional para la generacion de una malla para un reservatorio 
	TPZReadGIDGrid grid;
	TPZGeoMesh *meshgrid = grid.GeometricGIDMesh(archivo);
	if(!meshgrid->NElements())
		return 0;
	
	return meshgrid;
}

TPZGeoMesh *CreateGeoMesh(MElementType typeel) {
	TPZGeoMesh* gmesh = new TPZGeoMesh;
	REAL MaxX = 1.;
	TPZManVector<REAL> x0(3,0.), x1(3,MaxX);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
	
	switch(typeel) {
		case EOned:
		{
			x1[1] = x1[2] = 0.;
			int Qnodes = 2;
			
			gmesh->SetMaxNodeId(Qnodes-1);
			gmesh->NodeVec().Resize(Qnodes);
			TPZVec<TPZGeoNode> Node(Qnodes);
			
			TPZVec <int64_t> TopolLine(2);
			TPZVec <int64_t> TopolPoint(1);
			
			int64_t id = 0;
			for (int j=0; j<2;j++) {
				Node[id].SetNodeId(id);
				if(!j) Node[id].SetCoord(x0);//coord x
				else Node[id].SetCoord(x1);
				gmesh->NodeVec()[id] = Node[id];
				id++;
			}
			
			//indice dos elementos
			id = 0;
			
			TopolLine[0] = 0;
			TopolLine[1] = 1;
			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
			id++;
			
			TopolPoint[0] = 0;
			new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc0,*gmesh);
			id++;
			TopolPoint[0] = 1;
			new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc1,*gmesh);
			gmesh->BuildConnectivity();
			
			for ( int ref = 0; ref < 2; ref++ ) {
				TPZVec<TPZGeoEl *> filhos;
				int n = gmesh->NElements();
				for ( int i = 0; i < n; i++ ) {
					TPZGeoEl *gel = gmesh->ElementVec() [i];
					if (gel->Dimension() == 1) gel->Divide (filhos);
					gel->Divide (filhos);
				}//for i
			}//ref
		}
			break;
		case ETriangle:
		{
			x1[2] = 0.;
			TPZManVector<int> nx(2,1);   // subdivisions in X and in Y. 
			TPZGenGrid2D gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1
			gen.SetElementType(MMeshType::ETriangular);       // typeel = 0 means rectangular elements, typeel = 1 means triangular elements
			gen.Read(gmesh,matId);             // generating grid in gmesh
			gmesh->BuildConnectivity();
			TPZGeoElBC gbc10(gmesh->ElementVec()[0],3,bc0);
			TPZGeoElBC gbc11(gmesh->ElementVec()[0],4,bc1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[0],5,bc0);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case EQuadrilateral:
		{
			x1[2] = 0.;
			TPZManVector<int> nx(2,1);   // subdivisions in X and in Y. 
			TPZGenGrid2D gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1
			if(typeel==ETriangle) gen.SetElementType(MMeshType::ETriangular);       // typeel = 0 means rectangular elements, typeel = 1 means triangular elements
			else gen.SetElementType(EQuadrilateral);       // typeel = 0 means rectangular elements, typeel = 1 means triangular elements
			gen.Read(gmesh,matId);             // generating grid in gmesh
			gmesh->BuildConnectivity();
			TPZGeoElBC gbc10(gmesh->ElementVec()[0],4,bc0);
			TPZGeoElBC gbc11(gmesh->ElementVec()[0],5,bc1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[0],6,bc0);
			TPZGeoElBC gbc13(gmesh->ElementVec()[0],7,bc0);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case ECube:
		{
			TPZGeoMesh *gmeshtemp = NULL;
			x1[2] = 0.;
			TPZManVector<int> nx(2,1);   // subdivisions in X and in Y. 
			TPZGenGrid2D gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1
			gen.SetElementType(EQuadrilateral);       // typeel = 0 means rectangular elements, typeel = 1 means triangular elements
			gen.Read(gmeshtemp,matId);             // generating grid in gmesh
			
			REAL InitialH = MaxX;
			TPZExtendGridDimension gmeshextend(gmeshtemp,InitialH);
			gmesh = gmeshextend.ExtendedMesh(1,bc0,bc0);
			gmesh->BuildConnectivity();
			TPZGeoElBC gbc21(gmesh->ElementVec()[0],21,bc1);
			TPZGeoElBC gbc22(gmesh->ElementVec()[0],22,bc1);
			TPZGeoElBC gbc23(gmesh->ElementVec()[0],23,bc1);
			TPZGeoElBC gbc24(gmesh->ElementVec()[0],24,bc1);
			
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		default:
			return 0;
	}
	
	return gmesh;
}

void GradientReconstructionByLeastSquares(TPZCompEl *cel,TPZManVector<REAL,3> &center,TPZVec<REAL> &Grad) {
    TPZFMatrix<REAL> grad;
    int dim;
    dim = cel->Mesh()->Dimension();
    
    // Nada sera realizado para elementos com dimensao diferente da dimensao do problema
    if(!cel || cel->Dimension()!=dim) DebugStop();
	REAL solalfa;
	REAL solbeta;
    
    int k, side;
	// Integration order to calculate cell mean of solution
	int intOrder = 2;
    
	TPZStack<TPZCompElSide> neighs;
	int nneighs;
	
    center.Resize(3, 0.);
	TPZManVector<REAL> centerpsi(3,0.0), centerbeta(3,0.0);
	
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
	
	solalfa = MeanCell(cel,intOrder);
    
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
    int64_t counter=0;
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
        solbeta = MeanCell(gelbeta->Reference(),intOrder);
        
        for(k=0;k<dim;k++)
        {
            DeltaH(ineighs,k) = centerbeta[k] - center[k];
        }
        DifSol(ineighs,0) = solbeta - solalfa;
        counter ++;
    }
    
    // Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    grad = DeltaHTranspose*DifSol;
    A = DeltaHTranspose*DeltaH;
    if(counter > 0)
        A.SolveDirect(grad,ELU);
	
	// Return gradient vector
	Grad.Resize(dim);
	for(int j=0;j<dim;j++)
		Grad[j] = grad(j,0);
}
void PosProcessGradientReconstruction(TPZCompMesh *cmesh,TPZFMatrix<REAL> &datagradients){
    
    // Redimensionando a matriz dos dados da reconstruca de gradientes
    int dim  = cmesh->Dimension();
    int64_t nelem = cmesh->NElements();
    datagradients.Redim(nelem,2*dim+2);
	
    TPZManVector<REAL,3> center;
    TPZManVector<REAL> Grad(dim);
    
	int64_t i, k;
	int64_t counter = 0;
	
	
    TPZCompEl *cel;
    // Calculando el gradiente por elemento computacional
    for(i=0;i<nelem;i++) {
        
        cel = cmesh->ElementVec()[i];
        
        // Nada sera realizado para elementos con dimension diferente de la dimension del problema
        if(!cel || cel->Dimension()!=dim) continue;
		center.Fill(0.0);
		Grad.Fill(0.0);
        
        GradientReconstructionByLeastSquares(cel, center, Grad);
		
        //data of the vector gradiente
        for(k=0;k<dim;k++){
            datagradients(counter,k) = center[k];//centro do elemento
            datagradients(counter,dim+k) = Grad[k];//valor do gradiente
        }
		// Increment a last column to store volume of the finite element
		datagradients(counter,2*dim) = cel->Index();
		datagradients(counter,2*dim+1) = cel->VolumeOfEl();
        
		counter++;
	}
    // Redimensionando la matriz de los gradientes
	k = datagradients.Cols();
    datagradients.Resize(counter,k);
}

// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given
void PrintDataMeshVTK(TPZCompMesh *cmesh, char *filename,TPZFMatrix<REAL> &elData)
{
	std::ofstream file(filename);
#ifdef PZDEBUG
	if(!file.is_open())
		DebugStop();
#endif
	
	int dim = cmesh->Dimension();
	TPZGeoMesh *gmesh = cmesh->Reference();
	int64_t nelements = elData.Rows();
	
	std::stringstream connectivity, type, cellval1, cellval2, cellval3;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int64_t t, c, el;
	int64_t actualNode = -1, size = 0, nVALIDelements = 0;
	int64_t counternodes = gmesh->NNodes();
	TPZGeoEl *gel;
	TPZVec<REAL> centerpsi(3);
	TPZManVector<REAL> center(3);
	TPZManVector<REAL> gradient(3);
	int64_t counter = 0;
	
	for(el = 0; el < nelements; el++)
	{
		gel = cmesh->ElementVec()[elData(el,2*dim)]->Reference();
		if(!gel || gel->Reference()->Dimension()!=dim)
			continue;
		
		MElementType elt = gel->Type();
		int elNnodes = MElementType_NNodes(elt);
        
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(t = 0; t < elNnodes; t++)
		{
			actualNode = gel->NodeIndex(t);
			if(actualNode < 0) 
				DebugStop();
			
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
		type << elType << std::endl;
		REAL gradN, Norm, tempN = 0.0, temp = 0.0;
		TPZVec<REAL> v(3);
		for(c=0;c<dim;c++) {
			v[c] = elData(counter,c);
			tempN += v[c]*v[c];
			gradient[c] = elData(counter,dim+c);
			temp += gradient[c]*gradient[c];
		}
		Norm = sqrt(tempN);
		gradN = sqrt(temp);
		for(c=0;c<dim;c++) {
			if(IsZero(Norm)) v[c] = 0.;
			else v[c] /= Norm;
		}
		TPZVec<REAL> vort(3);
		vort[0] = -v[1];
		vort[1] = v[0];
		vort[2] = 0.0;

		cellval1 << gradN << std::endl;
		if(dim == 2) {
			cellval2 << vort[0]*gradient[0] + vort[1]*gradient[1] << std::endl;
			cellval3 << v[0]*gradient[0] + v[1]*gradient[1] << std::endl;
		}
		else {
			cellval2 << elData(counter,1) << std::endl;
			cellval3 << elData(counter,0) << std::endl;
		}			
		counter++;
		nVALIDelements++;
	}
	
	// Printing all nodes of the mesh
	file << counternodes << " float" << std::endl;
	for(t=0;t<counternodes;t++) {
		TPZGeoNode *node = &(gmesh->NodeVec()[t]);
		for(c = 0; c < 3; c++) {
			double coord = node->Coord(c);
			file << coord << " ";
		}			
		file << std::endl;
	}
	
	file << std::endl << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str() << std::endl;
	
	file << "CELL_DATA" << " " << nVALIDelements << std::endl;
	int nsubstructures = 3;
	file << "FIELD FieldData "<< nsubstructures << std::endl;
	
	file << "Substruct1 1 " << nVALIDelements << " float" << std::endl;
	
	file << cellval1.str();
		
	file << "Substruct2 1 " << nVALIDelements << " float" << std::endl;
	
	file << cellval2.str();
	
	file << "Substruct3 1 " << nVALIDelements << " float" << std::endl;
	
	file << cellval3.str();

	file.close();
}

void PrintCompMeshVTKWithGradientAsData(TPZCompMesh *cmesh,char *filename,TPZFMatrix<REAL> &elData) {
	int64_t i, size = cmesh->NElements();
	int64_t counter = 0;
	for(i=0;i<size;i++) {
		TPZCompEl *cel = cmesh->ElementVec()[i];
		if(!cel || cel->Dimension()!=cmesh->Dimension())
			continue;
		counter++;
	}
	if(counter != elData.Rows() || elData.Cols() != cmesh->Dimension())
		DebugStop();
	// Printing geometric mesh to visualization in Paraview
	PrintDataMeshVTK(cmesh, filename, elData);
}

REAL MeanCell(TPZCompEl *cel,int IntOrder) {
	TPZIntPoints *pointIntRule = ((TPZInterpolatedElement*)cel)->Reference()->CreateSideIntegrationRule((cel->Reference()->NSides())-1,IntOrder);
	int it, npoints = pointIntRule->NPoints();
	int dim = cel->Mesh()->Dimension();
	REAL integral = 0.0;
	TPZManVector<REAL> point(3,0.);
	TPZManVector<REAL> xpoint(3,0.);
	REAL weight;
	for (it=0;it<npoints;it++){
		pointIntRule->Point(it,point,weight);
		weight /= cel->Reference()->RefElVolume();
		cel->Reference()->X(point,xpoint);
		integral += weight * ExactSolution(dim,xpoint);
	}
	//REAL area = cel->Reference()->Volume();
	return integral;
}

void RefiningNearLine(int dim,TPZGeoMesh *gmesh,int nref) {
	
	int i;
	
	// Refinando no local desejado
	TPZManVector<REAL> point(3);
	point[0] = .5; point[1] = point[2] = 0.0;
	REAL r = 0.0;
	
	REAL radius = 0.9;
	for(i=0;i<nref;i++) {
		// To refine elements with center near to points than radius
		RefineGeoElements(dim,gmesh,point,r,radius);
		radius *= 0.8;
	}
	// Constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
	
	int i;
	
	// Refinando no local desejado
	TPZVec<REAL> point(3);
	point[0] = point[1] = point[2] = -0.25;
	if(dim==1) point[1] = point[2] = 0.0;
	else if(dim==2) point[2] = 0.0;
	REAL r = sqrt(2.);
	
	if(ntyperefs==2) {
		REAL radius = 0.2;
		for(i=0;i<nref;i+=2) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius);
			RefineGeoElements(dim,gmesh,point,r,radius);
			if(nref < 5) radius *= 0.35;
			else if(nref < 7) radius *= 0.2;
			else radius *= 0.1;
		}
		if(i==nref) {
			RefineGeoElements(dim,gmesh,point,r,radius);
		}
	}
	else {
		REAL radius = 0.4;
		for(i=0;i<nref+1;i++) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius);
			radius *= 0.6;
		}
	}
	// Constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			REAL distancia = TPZGeoEl::Distance(center,FirstNode);
			if(distancia > distance) distance = distancia;
			isdefined = true;
		}
		REAL centerdist = TPZGeoEl::Distance(center,point);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
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
		REAL centerdist = TPZGeoEl::Distance(center,point);
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

void GradReconstByLeastSquaresOnlyContEl(TPZCompEl *cel,TPZManVector<REAL,3> &center, TPZManVector<STATE> &solalfa,  TPZFMatrix<REAL> &grad, int var) {
    
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
	TPZManVector<REAL> centerpsi(3,0.0), centerbeta(3,0.0);
    TPZManVector<STATE> solbeta(nstates,0.0);;
	
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

void AssembleGlobalMatrix(TPZCompEl *el, TPZElementMatrix &ek, TPZElementMatrix &ef,TPZMatrix<STATE> & stiffmatrix, TPZFMatrix<STATE> &rhs){
    
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
    }
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
	TPZManVector<STATE> solalfa(nstates,0.), solbeta(nstates,0.);
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
    mA.Put(0,1,origX+2.); mA.Put(1,1,origY);
    mA.Put(0,2,origX+2.); mA.Put(1,2,origY+2.);
    mA.Put(0,3,origX); mA.Put(1,3,origY+2.);
    
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
	
	TPZVec <int64_t> TopolQuad(4);
	TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
	
	//indice dos nos
    int64_t id;
	id = 0;
    REAL dx = .5;
	for (int i=0; i<Qnodes/2;i++) {
		
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0, i*dx);//coord x
		Node[id].SetCoord(1, 0.);//coord y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	for (int i=0; i<Qnodes/2;i++) {
		
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0, 1. - i*dx);//coord x
		Node[id].SetCoord(1, 1.);//coord y
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

TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int pOrder,bool isdiscontinuous)
{
    /// criar materiais
	int64_t ngelem = gmesh->NElements();
	TPZGeoEl *gel;
	int dim = 0;
	for(int64_t j=0;j<ngelem;j++) {
		gel = gmesh->ElementVec()[j];
		if(gel->MaterialId() > 0 && dim < gel->Dimension())
			dim = gel->Dimension();
	}
	TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    
    REAL diff=1.;
    REAL conv =0.;
    TPZVec<REAL> convdir;
    convdir.Resize(dim,0.);
    convdir[0]=0.;
    
    material-> SetParameters(diff, conv, convdir);
    material->SetNoPenalty();
    material->SetNonSymmetric();
    
    TPZAutoPointer<TPZFunction<STATE> > myforce = new TPZDummyFunction<STATE>(ForcingF,5);
    material->SetForcingFunction(myforce);
    
    TPZMaterial * mat(material);
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    
    val2(0,0) =  0.0886226925452758;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    
    val2(0,0) =  -0.0886226925452758;
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
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
    
    TPZAutoPointer<TPZFunction<STATE> > fCC0 = new TPZDummyFunction<STATE>(Forcingbc0,5);
    TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(Forcingbc1,5);
    TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(Forcingbc2,5);
    TPZAutoPointer<TPZFunction<STATE> > fCC3 = new TPZDummyFunction<STATE>(Forcingbc3,5);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
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

void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]=  coef_a*x + coef_b*y;
    
}

void Forcingbc(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
#ifdef USING_BOOST
    disp[0]=0.05*sqrt(M_PI)*boost::math::erf(10.*(x-1.));
#else
//#include "specialFunctions.h"
    disp[0]=0.05*sqrt(M_PI)*erf(10.*(x-1.));
#endif
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

void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    double aux1 = -100.*(1 - 2.*x + x*x);
    double aux2 = exp(aux1);
    disp[0]=-(200.*x - 200.)*aux2;
}

void SolucaoExata(const TPZVec<REAL> &pt, TPZVec<STATE> &sol) {
	REAL x = pt[0];
    sol[0] = -100.*(x-1.)*(x-1.);
}

#include "TPZSkylineNSymStructMatrix.h"
void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	TPZBandStructMatrix full(Cmesh); //caso nao-simetrico
	//TPZSkylineStructMatrix full(Cmesh); //caso simetrico
    //TPZSkylineNSymStructMatrix full(Cmesh);//caso nao-simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	//step.SetDirect(ELDLt); //caso simetrico
	step.SetDirect(ELU);//caso nao simetrico
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

