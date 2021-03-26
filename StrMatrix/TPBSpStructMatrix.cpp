/**
 * @file
 * @brief Contains the implementation of the TPBSpStructMatrix methods. 
 */

#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzstrmatrix.h"

#include "pzgeoelbc.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

#include "pzanalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"

#include "pzdxmesh.h"
#include <fstream>
#include "pzelmat.h"
#include "pzysmp.h"

#include "pzbndcond.h"


using namespace std;

#ifndef STATE_COMPLEX
#include "pzmat2dlin.h"

int TPBSpStructMatrix::main() {
	
	int refine=5;
	int order=5;
	
	TPZGeoMesh gmesh;
	TPZCompMesh cmesh(&gmesh);
	double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},
		{0.,1.,0.}};
	
	int i,j;
	TPZVec<REAL> coord(3,0.);
	for(i=0; i<4; i++) {
		// initializar as coordenadas do no em um vetor
		for (j=0; j<3; j++) coord[j] = coordstore[i][j];
		
		// identificar um espa� no vetor onde podemos armazenar
		// este vetor
		gmesh.NodeVec ().AllocateNewElement ();
		
		// initializar os dados do n�
		gmesh.NodeVec ()[i].Initialize (i,coord,gmesh);
	}
	int el;
	TPZGeoEl *gel;
	for(el=0; el<1; el++) {
		
		// initializar os indices dos n�
		TPZVec<int64_t> indices(4);
		for(i=0; i<4; i++) indices[i] = i;
		// O proprio construtor vai inserir o elemento na malha
		int64_t index;
		gel = gmesh.CreateGeoElement(EQuadrilateral, indices,1,index);
	}
	gmesh.BuildConnectivity ();
	
	TPZVec<TPZGeoEl *> subel;
	//gel->Divide(subel);
	
	
	
	cout << "Refinement ";
	cin >> refine;
	cout << endl;
	
    DebugStop();
//	UniformRefine(refine,gmesh);
	
	
	TPZGeoElBC gelbc(gel,4,-4);
	TPZMat2dLin *meumat = new TPZMat2dLin(1);
	TPZFMatrix<STATE> xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
	meumat->SetMaterial (xk,xc,xf);
	TPZMaterial * meumatptr(meumat);
	cmesh.InsertMaterialObject(meumatptr);
	
	TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
	TPZMaterial * bnd = meumat->CreateBC (meumatptr,-4,0,val1,val2);
	
	cmesh.InsertMaterialObject(bnd);
	
	
	
	cout << "Interpolation order ";
	cin >> order;
	cout << endl;
	
	//     TPZCompEl::gOrder = order;
	cmesh.SetDefaultOrder(order);
	
	cmesh.AutoBuild();
	//	cmesh.AdjustBoundaryElements();
	cmesh.InitializeBlock();
	
	ofstream output("outputPar.dat");
	//	ofstream output2("outputNon.dat");
	//cmesh.Print(output);
	TPZAnalysis an(&cmesh,true,output);
	//	TPZAnalysis an2(&cmesh,output);
	
	TPZVec<int> numelconnected(cmesh.NEquations(),0);
	int64_t ic;
	//cout << "Nmero de Equa�es -> " << cmesh.NEquations() << endl;
	//cout.flush();
	
	ofstream out("cmeshBlock_out.txt");
	//	cmesh.Print(out);
	//	cmesh.Block().Print("Block",out);
	for(ic=0; ic<cmesh.ConnectVec().NElements(); ic++) {
		TPZConnect &cn = cmesh.ConnectVec()[ic];
		if(cn.HasDependency()) continue;
		int64_t seqn = cn.SequenceNumber();
		if(seqn < 0) continue;
		int64_t firsteq = cmesh.Block().Position(seqn);
		int64_t lasteq = firsteq+cmesh.Block().Size(seqn);
		int64_t ind;
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
	
	//TPZParFrontStructMatrix<TPZFrontSym> mat(&cmesh);
	TPBSpStructMatrix mat(&cmesh);
	
	//   TPZFStructMatrix mat2(&cmesh);
	//  mat->SetNumElConnected(numelconnected);
	//mat = CreateAssemble();
	
	//int threads=3;
	//cout << "Number of Threads  ";
	//cin >> threads;
	//cout << endl;
	
	//mat.SetNumberOfThreads(threads);
	//mat.SetNumberOfThreads(1);
	
	an.SetStructuralMatrix(mat);
	//	an2.SetStructuralMatrix(mat2);
	
	TPZStepSolver<STATE> sol;
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
	//an.Print("solution of frontal solver", output);
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
	 TPZFMatrix<STATE> *full = new TPZFMatrix(cmesh.NEquations(),cmesh.NEquations(),0.);
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

#endif

TPZStructMatrix * TPBSpStructMatrix::Clone(){
    return new TPBSpStructMatrix(*this);
}
TPZMatrix<STATE> * TPBSpStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    int64_t neq = fMesh->NEquations();
    if(fMesh->FatherMesh()) {
		cout << "TPZSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
		return new TPZFYsmpMatrix<STATE>(0,0);
    }
    TPZMatrix<STATE> *stiff = Create();//new TPZFYsmpMatrix(neq,neq);
    rhs.Redim(neq,1);
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    Assemble(*stiff,rhs, guiInterface);
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    return stiff;
}
TPZMatrix<STATE> * TPBSpStructMatrix::Create(){
    //checked
    
    int64_t neq = fEquationFilter.NActiveEquations();
    TPZFYsmpMatrix<STATE> * mat = new TPZFYsmpMatrix<STATE>(neq,neq);
	
    /**Rearange elements order*/
	//    TPZVec<int> elorder(fMesh->NEquations(),0);
	
	
    /**
     *Longhin implementation
	 */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
	//    int nnodes = 0;
    fMesh->ComputeElGraph(elgraph,elgraphindex);
    /**Creates a element graph*/
    TPZRenumbering metis;
    metis.SetElementsNodes(elgraphindex.NElements() -1 ,fMesh->NIndependentConnects());
    metis.SetElementGraph(elgraph,elgraphindex);
	
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    metis.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    int64_t i;
    int64_t nblock = nodegraphindex.NElements()-1;
    int64_t totalvar = 0;
    int64_t totaleq = 0;
    for(i=0;i<nblock;i++){
		int64_t iblsize = fMesh->Block().Size(i);
		int64_t iblpos = fMesh->Block().Position(i);
        int64_t numactive = fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        if (numactive != iblsize) {
            DebugStop();
        }
		totaleq += iblsize;
		int64_t icfirst = nodegraphindex[i];
		int64_t iclast = nodegraphindex[i+1];
		int64_t j;
		//longhin
		totalvar+=iblsize*iblsize;
		for(j=icfirst;j<iclast;j++) {
			int64_t col = nodegraph[j];
			int64_t colsize = fMesh->Block().Size(col);
			int64_t colpos = fMesh->Block().Position(col);
            int64_t numactive = fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
			totalvar += iblsize*colsize;
		}
    }
	
    int64_t ieq = 0;
    int64_t pos = 0;
	
    nblock=fMesh->NIndependentConnects();
	
    int64_t * Eq = new int64_t[totaleq+1];
    int64_t * EqCol = new int64_t[totalvar/2];
    STATE * EqValue = new STATE[totalvar/2];
    for(i=0;i<nblock;i++){
		int64_t iblsize = fMesh->Block().Size(i);
		int64_t iblpos = fMesh->Block().Position(i);
        int64_t numactive = fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
		int64_t ibleq;
		for(ibleq=0; ibleq<iblsize; ibleq++) {
			Eq[ieq] = pos;
			if(ieq%2) {
				ieq++;
				continue;
			}
			int64_t colsize = fMesh->Block().Size(i);
			int64_t colpos = fMesh->Block().Position(i);
			int64_t jbleq;
			for(jbleq=0; jbleq<colsize; jbleq++) {
				/**It can also be implemented using half the size of both columns and data vectors*/
				EqCol[pos] = -1;//colpos;
				EqValue[pos] = 0.;
				colpos++;
				pos++;
			}
			
			int64_t icfirst = nodegraphindex[i];
			int64_t iclast = nodegraphindex[i+1];
			int64_t j;
			for(j=icfirst;j<iclast;j++) {
				int64_t col = nodegraph[j];
				colsize = fMesh->Block().Size(col);
				colpos = fMesh->Block().Position(col);
                int64_t numactive = fEquationFilter.NumActive(colpos, colpos+colsize);
                if (!numactive) {
                    continue;
                }
                for(jbleq=0; jbleq<colsize; jbleq++) {
					EqCol[pos] = -1;//colpos;
					EqValue[pos] = 0.;
					colpos++;
					pos++;
				}
			}
			ieq++;
		}
    }
    Eq[ieq] = pos;
	/*    for(i=0;i<totalvar;i++){
	 if(i<totaleq+1){
	 cout << i <<  " " << Eq[i] << " "<< EqCol[i] << " " << EqValue[i] << endl;
	 }else{
	 cout << i <<  " " << " "<< EqCol[i] << " " << EqValue[i] << endl;
	 }
	 }
	 */
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}
TPBSpStructMatrix::TPBSpStructMatrix(TPZCompMesh *mesh) : TPZRegisterClassId(&TPBSpStructMatrix::ClassId),
TPZSpStructMatrix(mesh)
{
}

int TPBSpStructMatrix::ClassId() const{
    return Hash("TPBSpStructMatrix") ^ TPZSpStructMatrix::ClassId() << 1;
}