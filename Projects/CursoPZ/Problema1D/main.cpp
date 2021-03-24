
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzcheckgeom.h"

#include "pzfmatrix.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZMaterial.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "MiViga1D.h"

#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

std::string Archivo = PZSOURCEDIR;

TPZGeoMesh *CreateGeoMesh(std::string &nome);
TPZGeoMesh *CreateQuadrilateralMesh();
TPZGeoMesh *CreateQuadrilateralMesh2();
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);
TPZCompMesh *CreateMesh2D(TPZGeoMesh *gmesh);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void GradientReconstruction(TPZCompMesh *cmesh,TPZFMatrix<REAL> &gradients);

int main() {
    int counter = 3;
	int resolution = 8;

    TPZGeoMesh *gmesh;
    if(counter==1)
        gmesh = CreateQuadrilateralMesh();
    else
        gmesh = CreateQuadrilateralMesh2();
	int p = 2;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMesh2D(gmesh);
    
	if(counter!=1) {
        TPZInterpolatedElement *el = (TPZInterpolatedElement*)(cmesh->ElementVec()[0]);
        TPZVec<int64_t> subels;
        el->Divide(el->Index(),subels);
        el = (TPZInterpolatedElement*)(cmesh->ElementVec()[subels[3]]);
        el->Divide(el->Index(),subels);
        cmesh->AutoBuild();
        
        cmesh->AdjustBoundaryElements();
        cmesh->ExpandSolution();
        cmesh->CleanUpUnconnectedNodes();
    }

    TPZAnalysis an(cmesh);
    TPZSkylineStructMatrix strskyl(cmesh);
    an.SetStructuralMatrix(strskyl);
	TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
	direct->SetDirect(ECholesky);
	an.SetSolver(*direct);
	delete direct;
	direct = 0;
   // an.Run();

    int nrows = cmesh->Solution().Rows();
    int cols = cmesh->Solution().Cols();
    if(counter==1) {
        for(int i=0;i<nrows;i++) {
            for(int j=0;j<cols;j++) {
                cmesh->Solution().PutVal(i,j,0.);
            }
        }
        cmesh->Solution().PutVal(5,0,.25);
    }
	// Post processing
    TPZStack<std::string> scalarnames, vecnames;
    scalarnames.Push("Solution");
    scalarnames.Push("POrder");
    char namef[256];
    if(counter==1) {
        sprintf(namef,"Solution2D_%d.vtk",counter);
        an.DefineGraphMesh(2,scalarnames,vecnames,namef);
        
        an.PostProcess(resolution,2);
    }
    if(counter==2) {
        for(int k=0;k<nrows;k++) {
            for(int i=0;i<nrows;i++) {
                for(int j=0;j<cols;j++) {
                    cmesh->Solution().PutVal(i,j,0.);
                }
            }
			if(k==33)
				cmesh->Solution().PutVal(33,0,1.);
			else if(k==17) {
				cmesh->Solution().PutVal(17,0,1.);
				cmesh->Solution().PutVal(37,0,0.5);
			}
			else
	            cmesh->Solution().PutVal(k,0,.25);
            sprintf(namef,"Solution2D_%d_%d.vtk",counter,k);
            an.DefineGraphMesh(2,scalarnames,vecnames,namef);
            an.PostProcess(resolution,2);
		}
    }
    if(counter==3) {
        for(int i=0;i<nrows;i++) {
            for(int j=0;j<cols;j++) {
                cmesh->Solution().PutVal(i,j,0.);
            }
        }
        for(int k=0;k<nrows;k++) {
			if(k==2)
                cmesh->Solution().PutVal(2,0,0.1875-0.125);
            else if(k==33)
				cmesh->Solution().PutVal(33,0,0.1875);
			else if(k==17) {
				cmesh->Solution().PutVal(17,0,0.25);
				cmesh->Solution().PutVal(37,0,0.125);
			}
			else if(k==32)
	            cmesh->Solution().PutVal(k,0,0.109375-0.09375);
			else if(k==21)
	            cmesh->Solution().PutVal(k,0,0.234375-(0.09375+0.125));
		}
        sprintf(namef,"Solution2D_%d.vtk",counter);
        an.DefineGraphMesh(2,scalarnames,vecnames,namef);
        an.PostProcess(resolution,2);
    }
    
    return 0;
}
// uni-dimensional problem for elasticity
int main1D() {
    
	Archivo += "/Projects/CursoPZ/Problema1D/";
	Archivo += "Viga1D.dump";

	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMesh(Archivo);
    gmesh->Print();
	UniformRefine(gmesh,1);

	// Creating computational mesh (approximation space and materials)
	int p = 3;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMesh(gmesh);
    cmesh->Print();
	// Solving linear equations
	// Initial steps
	TPZAnalysis an (cmesh);
	TPZSkylineStructMatrix strskyl(cmesh);
	an.SetStructuralMatrix(strskyl);
	// Solver (is your choose)
	TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
	direct->SetDirect(ECholesky);
	an.SetSolver(*direct);
	delete direct;
	direct = 0;
/*
	// Caso no simetrico
//	TPZFStructMatrix full(cmesh);
	TPZBandStructMatrix full(cmesh);
	an.SetStructuralMatrix(full);
	an.Solution().Zero();
	TPZStepSolver<STATE> step;
	step.SetDirect(ELU);
	an.SetSolver(step);
	*/
	an.Run();
	
    
    TPZFMatrix<REAL> gradients;
	// Redimension of gradients storage
	GradientReconstruction(cmesh,gradients);
	gradients.Print("gradients");
    
	// Post processing
	TPZManVector<std::string> scalarnames(0), vecnames(2);
	vecnames[0] = "desplazamiento";
	vecnames[1] = "giro";
	an.DefineGraphMesh(1,scalarnames,vecnames,"MiSolucion1D.vtk");

	an.PostProcess(2);
    return 0;
}

TPZGeoMesh *CreateQuadrilateralMesh() {
    REAL co[4][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.}
    };
    int64_t indices[1][4] = {{0,1,2,3}};
    
    const int nelem = 1;
    int nnode = 4;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int64_t nod;
    for(nod=0; nod<nnode; nod++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int64_t el;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;
}
TPZGeoMesh *CreateQuadrilateralMesh2() {
    REAL co[4][3] = {
        {1.,0.,0.},
        {2.,0.,0.},
        {2.,1.,0.},
        {1.,1.,0.}
    };
    int64_t indices[1][4] = {{0,1,2,3}};
    
    const int nelem = 1;
    int nnode = 4;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int64_t nod;
    for(nod=0; nod<nnode; nod++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int64_t el;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;
}
TPZCompMesh *CreateMesh2D(TPZGeoMesh *gmesh) {
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
//    cmesh->SetAllCreateFunctionsDiscontinuous();
	// Creating Poisson material
    int dim = 2;
    int MaterialId = 1;
	TPZMaterial *mat = new TPZMatPoisson3d(MaterialId,dim);
	TPZVec<REAL> convd(3,0.);
	((TPZMatPoisson3d *)mat)->SetParameters(1.,0.,convd);
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
//	cmesh->SetAllCreateFunctionsContinuous();
    
	// Boundary conditions
	// Dirichlet
	TPZFMatrix<STATE> val1(dim,dim,0.),val2(dim,1,0.);
	TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
	cmesh->InsertMaterialObject(bnd);
	
	cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}


void GradientReconstruction(TPZCompMesh *cmesh,TPZFMatrix<REAL> &gradients) {
    
    gradients.Redim(cmesh->NElements(),1);
    
	TPZCompEl *cel;
	int i = 0, side, nneighs;
	TPZVec<REAL> normal(3,0.0);
	TPZVec<REAL> center(3,0.);
	TPZVec<STATE> solalfa(3,0.), solbeta(3,0.);
	
	REAL measure, sidemeasure;
	for(i=0;i<cmesh->NElements();i++) {
		cel = cmesh->ElementVec()[i];
		if(!cel /*|| cel->Dimension()!=cmesh->Dimension()*/) continue;
		TPZStack<TPZCompElSide> neighs;
		for(side = 0; side < cel->Reference()->NSides(); side++) {
			neighs.Resize(0);
			TPZGeoElSide gelside(cel->Reference(),side);
			if(gelside.Dimension() != cel->Dimension() - 1) continue;
			gelside.EqualorHigherCompElementList2(neighs,1,0);
            
			if(!neighs.NElements()) continue;
			measure = cel->VolumeOfEl();
			// for testing, we are using the solution on centroid
			TPZVec<REAL> centeralfa(3,0.);
			TPZGeoElSide gelalfa(cel->Reference(),cel->Reference()->NSides() -1);
			gelalfa.CenterPoint(centeralfa);
			cel->Solution(centeralfa,1,solalfa);
			for(nneighs=0;nneighs<neighs.NElements();nneighs++) {
                
				if(neighs[nneighs].Element() == cel || neighs[nneighs].Element()->Dimension() != cel->Dimension()) continue;
				sidemeasure = 1.;   //neighs[nneighs].Reference().Area();   //gelside->Area();
				neighs[nneighs].Reference().CenterPoint(center);
				neighs[nneighs].Reference().Normal(center,gelside.Element(),neighs[nneighs].Reference().Element(),normal);
                
				// for testing, we are using the solution on centroid
				TPZGeoElSide gelbeta(neighs[nneighs].Reference().Element(),neighs[nneighs].Reference().Element()->NSides() -1);
				gelbeta.CenterPoint(centeralfa);
				gelbeta.Element()->Reference()->Solution(centeralfa,1,solbeta);
				
				// Incrementando no gradients para o elemento alfa o valor por vizinho
                for(int j=0;j<cel->Reference()->Dimension();j++){
                    gradients(i,j) += (sidemeasure * normal[j] *(solbeta[1] - solalfa[1]));
                }
			}
		}
        for(int j=0;j<cel->Reference()->Dimension();j++){
            if(!IsZero(measure))
                gradients(i,j) /= (2*measure);
        }
	}
}

//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh(std::string &archivo) {

	// Ejemplo uni-dimensional para la generacion de una malla para un reservatorio 
	TPZReadGIDGrid grid;
	TPZGeoMesh *meshgrid = grid.GeometricGIDMesh(archivo);

	int nelem = meshgrid->NElements();
	if(!nelem)
		return 0;

	// Buscando 

	// Criando los elementos puntuales en los extremos
	TPZGeoEl *gel = meshgrid->ElementVec()[0];
	TPZGeoElBC gbc1(gel,0,3);
	gel = meshgrid->ElementVec()[nelem-1];
	TPZGeoElBC gbc2(gel,1,2);

	// Re construindo la conectividad en la malha geometrica
	meshgrid->BuildConnectivity();

	return meshgrid;
}

//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());

    // Creating elasticity material
    TPZMaterial * mat = new TPZMiViga1D(1,2000.,0.3,0.4);
    cmesh->InsertMaterialObject(mat);
	((TPZMiViga1D *)mat)->SetLoad(0.);

	// Creating four boundary condition
    TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);
	val2(0,0) = 0.;
	val2(1,0) = 1.0;
	TPZMaterial *bcRight;
	TPZMaterial *bcLeft;
	// Condicion Dirichlet desplazamiento y giro en el extremo inicial de la viga (v,phi) = (0,0)
    bcLeft = mat->CreateBC(mat,2,0,val1,val2);
    cmesh->InsertMaterialObject(bcLeft);
	// Condicion de Neumann colocando un peso en el extremo
	val2(0,0) = 1.;
	val2(1,0) = 0.2;
    bcRight = mat->CreateBC(mat,3,1,val1,val2);
	cmesh->InsertMaterialObject(bcRight);

	// Inserting boundary conditions into computational mesh
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
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
//	gmesh->BuildConnectivity();
}
