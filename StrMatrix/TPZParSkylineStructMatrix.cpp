/**
 * @file
 * @brief Contains the implementation of the TPZParSkylineStructMatrix methods. 
 */
#include "TPZFrontStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzskylmatpar.h"
#include "pzvec.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzmat2dlin.h"
#include "pzbndcond.h"

#include "pzanalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"

#include "pzdxmesh.h"
#include <fstream>
using namespace std;

#include "pzelmat.h"
#include "TPZFrontStructMatrix.h"


TPZParSkylineStructMatrix::TPZParSkylineStructMatrix(const TPZParSkylineStructMatrix &cp) : TPZSkylineStructMatrix(cp){
	fNumThreads = cp.fNumThreads;
}

TPZParSkylineStructMatrix::TPZParSkylineStructMatrix(TPZCompMesh *mesh, int numthreads) : TPZSkylineStructMatrix(mesh)
{
	fNumThreads = numthreads;
}

TPZStructMatrix * TPZParSkylineStructMatrix::Clone(){
    return new TPZParSkylineStructMatrix(*this);
}

TPZMatrix<REAL> * TPZParSkylineStructMatrix::Create(){
    int neq = fMesh->NEquations();
    TPZVec<int> skyline;
	if (fOnlyInternal) {
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (fMesh);
		if (submesh) {
			submesh->SkylineInternal(skyline);
		}
		else {
			fMesh->Skyline(skyline);
		}
	}
	else {
		fMesh->Skyline(skyline);
	}
    if(HasRange())
    {
		neq = fMaxEq-fMinEq;
		FilterSkyline(skyline);
    }
	else {
		// This statement is needed for compatibility with TPZSubCompMesh The number of equations of the stiffness matrix corresponds to "only" the internal nodes
		neq = skyline.NElements();
	}
    return new TPZSkylParMatrix<REAL>(neq,skyline,fNumThreads);
}
TPZMatrix<REAL> * TPZParSkylineStructMatrix::CreateAssemble(TPZFMatrix<REAL> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix<REAL> *mat = Create();
	rhs.Redim(mat->Rows(),1);
	Assemble(*mat,rhs,guiInterface);
	return mat;
}
int TPZParSkylineStructMatrix::main() {
	
	int refine=5;
	int order=6;
	
	TPZGeoMesh gmesh;
	TPZCompMesh cmesh(&gmesh);
	double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},
		{0.,1.,0.}};
	
	int i,j;
	TPZVec<REAL> coord(3,0.);
	for(i=0; i<4; i++) {
		// initializar as coordenadas do no em um vetor
		for (j=0; j<3; j++) coord[j] = coordstore[i][j];
		
		// identificar um espa�o no vetor onde podemos armazenar
		// este vetor
		
		// initializar os dados do n�
		gmesh.NodeVec ()[i].Initialize (i,coord,gmesh);
	}
	int el;
	//TPZGeoEl *gel;
	for(el=0; el<1; el++) {
		
		// initializar os indices dos n�s
		TPZVec<int> indices(4);
		for(i=0; i<4; i++) indices[i] = i;
		// O proprio construtor vai inserir o elemento na malha
		int index;
		/*gel = */gmesh.CreateGeoElement(EQuadrilateral, indices,1,index);
	}
	gmesh.BuildConnectivity ();
	
	TPZVec<TPZGeoEl *> subel;
	//gel->Divide(subel);
	
	
	/*
	 cout << "Refinement ";
	 cin >> refine;
	 cout << endl;
	 */
	UniformRefine(refine,gmesh);
	
	
	TPZMat2dLin *meumat = new TPZMat2dLin(1);
	TPZFMatrix<REAL> xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
	meumat->SetMaterial (xk,xc,xf);
	TPZAutoPointer<TPZMaterial> meumatptr = meumat;
	cmesh.InsertMaterialObject(meumatptr);
	
	TPZFMatrix<REAL> val1(1,1,0.),val2(1,1,0.);
	TPZAutoPointer<TPZMaterial> bnd = meumat->CreateBC (meumatptr,-4,0,val1,val2);
	cmesh.InsertMaterialObject(bnd);
	
	
	/*
	 cout << "Interpolation order ";
	 cin >> order;
	 cout << endl;
	 */
	//	TPZCompEl::gOrder = order;
	cmesh.SetDefaultOrder(order);
	
	cmesh.AutoBuild();
	//	cmesh.AdjustBoundaryElements();
	cmesh.InitializeBlock();
	
	ofstream output("outputSkyPar.dat");
	//	ofstream output2("outputNon.dat");
	cmesh.Print(output);
	TPZAnalysis an(&cmesh,output);
	//	TPZAnalysis an2(&cmesh,output);
	
	TPZVec<int> numelconnected(cmesh.NEquations(),0);
	int ic;
	//cout << "N�mero de Equa��es -> " << cmesh.NEquations() << endl;
	//cout.flush();
	
	ofstream out("cmeshBlock_out.txt");
	//	cmesh.Print(out);
	//	cmesh.Block().Print("Block",out);
	for(ic=0; ic<cmesh.ConnectVec().NElements(); ic++) {
		TPZConnect &cn = cmesh.ConnectVec()[ic];
		if(cn.HasDependency() || cn.IsCondensed()) continue;
		int seqn = cn.SequenceNumber();
		if(seqn < 0) continue;
		int firsteq = cmesh.Block().Position(seqn);
		int lasteq = firsteq+cmesh.Block().Size(seqn);
		int ind;
		int temp = cmesh.ConnectVec()[ic].NElConnected();
		for(ind=firsteq;ind<lasteq;ind++) {
			numelconnected[ind] = temp;//cmesh.ConnectVec()[ic].NElConnected();
		}
	}
	//	//cout << "nequations " << numelconnected.NElements();
	//	for(ic=0;ic<numelconnected.NElements(); ic++) //cout << numelconnected[ic] <<' ';
	//	//cout << endl;
	//	//cout.flush();
	
	//	TPZFrontMatrix<TPZFileEqnStorage, TPZFrontNonSym> *mat = new TPZFrontMatrix<TPZFileEqnStorage, TPZFrontNonSym>(cmesh.NEquations());
	//TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym> *mat = new TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym>(cmesh.NEquations());
	//TPZFrontMatrix<TPZStackEqnStorage> *mat = new TPZFrontMatrix<TPZStackEqnStorage>(cmesh.NEquations());
	
	const int numthreads = 2;
	TPZParSkylineStructMatrix mat(&cmesh,numthreads);
	
	//   TPZFStructMatrix mat2(&cmesh);
	//  mat->SetNumElConnected(numelconnected);
	//mat = CreateAssemble();
	/*int threads = 0;
	 cout << "Number of Threads  ";
	 cin >> threads;
	 cout << endl;*/
	
	//mat.SetNumberOfThreads(threads);
	//     mat.SetNumberOfThreads(1);
	
	an.SetStructuralMatrix(mat);
	//	an2.SetStructuralMatrix(mat2);
	
	TPZStepSolver<REAL> sol;
	//	sol.SetDirect(ELU);
	sol.SetDirect(ECholesky);
	//	TPZStepSolver sol2;
	//	sol2.SetDirect(ECholesky);
	//	sol.SetDirect(ELU);
	
	
	an.SetSolver(sol);
	//     an2.SetSolver(sol2);
	//	mat->SetNumElConnected(numelconnected);
	//	mat->SetFileName("longhin.bin");
	//	an.Solver().SetDirect(ELU);
	//	mat->FinishWriting();
	//  mat->SetFileName('r',"longhin.bin");
	//	//cout << "******************************************************************************************************AQUI 1" << endl;
	an.Run(output);
	an.Print("solution of Skyline solver", output);
	//	//cout << "******************************************************************************************************AQUI 2" << endl;
	//	an2.Run(output2);
	//	an2.Print("solution of frontal solver", output2);
	/*
	 TPZVec<char *> scalnames(1);
	 scalnames[0] = "state";
	 
	 TPZVec<char *> vecnames(0);
	 
	 TPZDXGraphMesh graph(&cmesh,2,meumat,vecnames,scalnames);
	 ofstream *dxout = new ofstream("poisson.dx");
	 graph.SetOutFile(*dxout);
	 graph.SetResolution(0);
	 
	 //an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
	 //an.Print("FEM SOLUTION ",output);
	 //an.PostProcess(1);
	 int istep = 0,numstep=1;
	 
	 graph.DrawMesh(numstep+1);
	 graph.DrawSolution(0,0);
	 
	 TPZAnalysis an2(&cmesh,output);
	 TPZFMatrix<REAL> *full = new TPZFMatrix(cmesh.NEquations(),cmesh.NEquations(),0.);
	 an2.SetMatrix(full);
	 an2.Solver().SetDirect(ELU);
	 an2.Run(output);
	 an2.Print("solution of full matrix", output);
	 
	 //	full->Print("full decomposed matrix");
	 */
	output.flush();
	cout.flush();
	return 0;
	
}

