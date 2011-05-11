/*
 * HDiv.cpp
 *
 *  Created on: Jun 29, 2009
 *      Author: phil
 */

#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
//#include "pzerror.h"
#include "pzgeoel.h"
//#include "pzcosys.h"
#include "pzmatrix.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "pzpoisson3d.h"
#include "pzcheckgeom.h"
#include "pzbndcond.h"
#include "pzelmat.h"
#include "pzsubcmesh.h"
#include "pzgengrid.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("main.hdiv"));
#endif
TPZCompMesh *CreateMesh2d(int nelrow);

// print the element stiffness matrices
void PrintMesh(TPZCompMesh *cmesh);

// substructure the mesh
// put all elements with materialid in the submesh
int SubStructure(TPZCompMesh *cmesh, int materialid);

void ValFunction(TPZVec<REAL> &loc, TPZFMatrix &Val1, TPZVec<REAL> &Val2, int &BCType);

int main()
{
	InitializePZLOG();
	TPZCompEl::SetgOrder(2);
	int nelrow = 5;
	TPZCompMesh *cmesh = CreateMesh2d(nelrow);
	PrintMesh(cmesh);
	int submeshindex = -1;
//	if(nelrow > 1)
    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
	if(1)
	{
		submeshindex = SubStructure(cmesh,1);
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (cmesh->ElementVec()[submeshindex]);
		submesh->SetNumberRigidBodyModes(1);
		cmesh->ExpandSolution();
#ifdef LOG4CXX
        {
            std::stringstream sout;
            cmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        int numthreads = 2;
		submesh->SetAnalysisFrontal(numthreads,gui);
																  
	}
	cmesh->LoadReferences();
	PrintMesh(cmesh);
	TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
	TPZAnalysis analysis(cmesh);
	TPZFStructMatrix str(cmesh);

	// zero the submesh index so that it won't be assembled
	if(submeshindex >= 0) cmesh->ElementVec()[submeshindex] = 0;
	
	TPZFMatrix rhs;
	// this will be the matrix with only the boundary condition
	TPZAutoPointer<TPZMatrix> A = str.CreateAssemble(rhs,gui);
	
	// keep only the subcmesh
	cmesh->ElementVec() = elvec;
	int nel = elvec.NElements();
	int iel;
	for(iel=0; iel<nel; iel++)
	{
		if (iel != submeshindex) {
			cmesh->ElementVec()[iel]=0;
		}
	}
#ifdef LOG4CXX
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	TPZAutoPointer<TPZMatrix> B = str.CreateAssemble(rhs,gui);
	{
		std::ofstream eig("steklovhdivQ2h5.nb");
		A->Print("Acondense = ",eig,EMathematicaInput);
		B->Print("Bcondense = ",eig,EMathematicaInput);
	}	
	// restore the original state
	cmesh->ElementVec() = elvec;
	
/*	analysis.SetStructuralMatrix(str);
	TPZStepSolver step;
	step.SetDirect(ELU);
	analysis.SetSolver(step);
	analysis.Run();
 */
	{
        std::ifstream in("Eigvec.txt");
		int nel = analysis.Solution().Rows();
		int i;
		for (i=0; i<nel; i++) {
			in >> analysis.Solution()(i,0);
		}
		cmesh->LoadSolution(analysis.Solution());
	}
	analysis.Solution().Print("solution",std::cout);
	
	return 0;
}

//*************************************
//************Option 0*****************
//*******Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh2d(int nelrow){

  //malha quadrada de nr x nc
  const	int numrel = nelrow;
  const	int numcel = nelrow;
	TPZManVector<int,2> nx(2,nelrow);
	TPZManVector<REAL,2> x0(2,0.),x1(2,1.);
//	TPZGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numl = 1, REAL rot = 0.5);
	TPZGenGrid gen(nx,x0,x1,1,0.);
	// criar um objeto tipo malha geometrica
	TPZGeoMesh *geomesh = new TPZGeoMesh();
	gen.Read(*geomesh);
	
	gen.SetBC(geomesh, 0, -2);
	gen.SetBC(geomesh, 1, -2);
	gen.SetBC(geomesh, 2, -2);
	gen.SetBC(geomesh, 3, -2);

#ifdef LOG4CXX
  {
	  TPZFNMatrix<100> normals;
	  TPZManVector<int> indexes;
	  geomesh->ElementVec()[0]->ComputeNormals(normals,indexes);
	  std::stringstream sout;
	  normals.Print("Normal Vectors",sout);
	  sout << std::endl << "Side indexes " << indexes << std::endl;
	  geomesh->ElementVec()[0]->Print(sout);
	  LOGPZ_INFO(logger,sout.str())
  }
#endif
  //Divisão dos elementos
  TPZVec<TPZGeoEl *> sub;

  geomesh->BuildConnectivity();

#ifdef LOG4CXX
	{
		std::stringstream sout;
		geomesh->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
  // Criação da malha computacional
  TPZCompMesh *comp = new TPZCompMesh(geomesh);

  // Criar e inserir os materiais na malha
  TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
  TPZAutoPointer<TPZMaterial> automat(mat);
  //mat->problema = 2;
  comp->InsertMaterialObject(automat);


  // Condições de contorno
  // Dirichlet
	TPZFMatrix val1(1,1,0.),val2(1,1,15.);
	TPZMaterial *bnd = automat->CreateBC (automat,-1,0,val1,val2);
	comp->InsertMaterialObject(bnd);

	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-2,2,val1,val2);
	TPZBndCond *bndcond = dynamic_cast<TPZBndCond *> (bnd);
	bndcond->SetValFunction(ValFunction);
	comp->InsertMaterialObject(bnd);

	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-3,2,val1,val2);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-4,2,val1,val2);
	comp->InsertMaterialObject(bnd);
	
	comp->SetAllCreateFunctionsHDiv();
  // Ajuste da estrutura de dados computacional
  comp->AutoBuild();


  comp->AdjustBoundaryElements();
  comp->CleanUpUnconnectedNodes();
//    comp->Print(cout);
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

// print the element stiffness matrices
void PrintMesh(TPZCompMesh *cmesh)
{
	int nel = cmesh->NElements();
	int iel;
	for(iel=0; iel<nel; iel++)
	{
		TPZInterpolationSpace *intel;
		intel = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[iel]);
		if(intel)
		{
			TPZMaterialData data;
			intel->InitMaterialData(data);
			TPZElementMatrix ek(cmesh,TPZElementMatrix::EK),ef(cmesh,TPZElementMatrix::EF);
			intel->InitializeElementMatrix(ek, ef);
			intel->CalcStiff(ek,ef);
			int nshape = intel->NShapeF();
			int dim = intel->Reference()->Dimension();
			TPZManVector<REAL> xi(dim,0.);
			TPZFNMatrix<100> phi(nshape,1,0.),dphi(2,nshape,0.);
			intel->Shape(xi,phi,dphi);
#ifdef LOG4CXX
			{
				std::stringstream sout;
				phi.Print("Shape functions ",sout);
				dphi.Print("Derivative shape functions ",sout);
				ek.Print(sout);
				ef.Print(sout);
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
		}
	}
}

// substructure the mesh
// put all elements with materialid in the submesh
int SubStructure(TPZCompMesh *cmesh, int materialid)
{
	int index;
	TPZSubCompMesh *submesh = new TPZSubCompMesh(*cmesh,index);
	
	int nelem = cmesh->NElements();
	int iel;
	for(iel = 0; iel<nelem; iel++)
	{
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if(!cel || cel == submesh) continue;
		TPZAutoPointer<TPZMaterial> celmat = cel->Material();
		if(!celmat) continue;
		int matid = celmat->Id();
		if(matid == materialid)
		{
			submesh->TransferElement(cmesh, iel);
		}
	}
	submesh->MakeAllInternal();
	cmesh->CleanUpUnconnectedNodes();
	
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

void ValFunction(TPZVec<REAL> &loc, TPZFMatrix &Val1, TPZVec<REAL> &Val2, int &BCType)
{
	BCType = 2;
	Val1.Redim(1, 1);
	Val1(0,0) = 1.;
	Val2[0] = loc[0];
}
