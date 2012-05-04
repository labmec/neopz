/**
 * @file
 * @brief Contains the implementation of the TPZSpStructMatrix methods. 
 */

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
#include "pzmetis.h"
#include "pzbndcond.h"
#include "TPZTimer.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.StrMatrix"));
#endif

using namespace std;

TPZStructMatrix * TPZSpStructMatrix::Clone(){
    return new TPZSpStructMatrix(*this);
}
TPZMatrix<STATE> * TPZSpStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,
                                              TPZAutoPointer<TPZGuiInterface> guiInterface){
	
    LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble starting");
	
    int neq = fMesh->NEquations();
    if(fMesh->FatherMesh()) {
		cout << "TPZSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
		return new TPZFYsmpMatrix<STATE>(0,0);
    }
    TPZMatrix<STATE> *stiff = Create();//new TPZFYsmpMatrix(neq,neq);
	//    TPZFYsmpMatrix *mat = dynamic_cast<TPZFYsmpMatrix *> (stiff);
    rhs.Redim(neq,1);
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    TPZTimer before("Assembly of a sparse matrix");
    before.start();
    LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble calling Assemble()");
	Assemble(*stiff,rhs,guiInterface);
    before.stop();
    std::cout << __PRETTY_FUNCTION__ << " " << before << std::endl;
    //    mat->ComputeDiagonal();
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble exiting");
    return stiff;
}
TPZMatrix<STATE> * TPZSpStructMatrix::Create(){
    int neq = fMesh->NEquations();
    if(HasRange())
    {
		neq = fMaxEq-fMinEq;
    }
    else
    {
		fMinEq = 0;
		fMaxEq = neq;
    }
	/*    if(fMesh->FatherMesh()) {
	 TPZSubCompMesh *smesh = (TPZSubCompMesh *) fMesh;
	 neq = smesh->NumInternalEquations();
	 }*/
    TPZFYsmpMatrix<STATE> * mat = new TPZFYsmpMatrix<STATE>(neq,neq);
	
    /**
     *Longhin implementation
	 */
    TPZStack<int> elgraph;
    TPZVec<int> elgraphindex;
    //    int nnodes = 0;
    fMesh->ComputeElGraph(elgraph,elgraphindex);
    /**Creates a element graph*/
    TPZMetis metis(elgraphindex.NElements() -1 ,fMesh->NIndependentConnects());
    metis.SetElementGraph(elgraph,elgraphindex);
	
    TPZVec<int> nodegraph;
    TPZVec<int> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    metis.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    int i;
    int nblock = nodegraphindex.NElements()-1;
    int totalvar = 0;
    int totaleq = 0;
    for(i=0;i<nblock;i++){
		int iblsize = fMesh->Block().Size(i);
		int iblpos = fMesh->Block().Position(i);
		if(iblpos < fMinEq || iblpos >= fMaxEq) continue;
		totaleq += iblsize;
		int icfirst = nodegraphindex[i];
		int iclast = nodegraphindex[i+1];
		int j;
		//longhin
		totalvar+=iblsize*iblsize;
		for(j=icfirst;j<iclast;j++) {
			int col = nodegraph[j];
			int colsize = fMesh->Block().Size(col);
			int colpos = fMesh->Block().Position(col);
			if(colpos < fMinEq || colpos >= fMaxEq) continue;
			totalvar += iblsize*colsize;
		}
    }
	
    int ieq = 0;
    int pos = 0;
	
    nblock=fMesh->NIndependentConnects();
	
    int * Eq = new int[totaleq+1];
    int * EqCol = new int[totalvar];
    STATE * EqValue = new STATE [totalvar];
    for(i=0;i<nblock;i++){
		int iblsize = fMesh->Block().Size(i);
		int iblpos = fMesh->Block().Position(i);
		if(iblpos < fMinEq || iblpos >= fMaxEq) continue;
		int ibleq;
		for(ibleq=0; ibleq<iblsize; ibleq++) {
			Eq[ieq] = pos;
			int colsize,colpos,jbleq;
			int diagonalinsert = 0;
			int icfirst = nodegraphindex[i];
			int iclast = nodegraphindex[i+1];
			int j;
			for(j=icfirst;j<iclast;j++) {
				int col = nodegraph[j];
				if(!diagonalinsert && col > i)
				{
					diagonalinsert = 1;
					int colsize = fMesh->Block().Size(i);
					int colpos = fMesh->Block().Position(i);
					int jbleq;
					for(jbleq=0; jbleq<colsize; jbleq++) {
						//             if(colpos+jbleq == ieq) continue;
						EqCol[pos] = colpos+jbleq-fMinEq;
						EqValue[pos] = 0.;
						//            colpos++;
						pos++;
					}
				}
				colsize = fMesh->Block().Size(col);
				colpos = fMesh->Block().Position(col);
				if(colpos < fMinEq || colpos >= fMaxEq) continue;
				for(jbleq=0; jbleq<colsize; jbleq++) {
					EqCol[pos] = colpos-fMinEq;
					EqValue[pos] = 0.;
					colpos++;
					pos++;
				}
			}
			if(!diagonalinsert)
			{
				diagonalinsert = 1;
				int colsize = fMesh->Block().Size(i);
				int colpos = fMesh->Block().Position(i);
				int jbleq;
				for(jbleq=0; jbleq<colsize; jbleq++) {
					//             if(colpos+jbleq == ieq) continue;
					EqCol[pos] = colpos+jbleq-fMinEq;
					EqValue[pos] = 0.;
					//            colpos++;
					pos++;
				}
			}
			ieq++;
		}
    }
    Eq[ieq] = pos;
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}

TPZSpStructMatrix::TPZSpStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{}

#ifndef STATE_COMPLEX
#include "pzmat2dlin.h"

int TPZSpStructMatrix::main() {
	
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
		
		// initializar os dados do n�       gmesh.NodeVec ()[nodeindex].Initialize (i,coord,gmesh);
	}
	int el;
	TPZGeoEl *gel;
	for(el=0; el<1; el++) {
		
		// initializar os indices dos nos
		TPZVec<int> indices(4);
		for(i=0; i<4; i++) indices[i] = i;
		// O proprio construtor vai inserir o elemento na malha
		//       gel = new TPZGeoElQ2d(el,indices,1,gmesh);
		int index;
		gel = gmesh.CreateGeoElement(EQuadrilateral,indices,1,index);
	}
	gmesh.BuildConnectivity ();
	
	TPZVec<TPZGeoEl *> subel;
	
	cout << "Refinement ";
	cin >> refine;
	cout << endl;
	
	UniformRefine(refine,gmesh);
	
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
	
	//TPZCompEl::gOrder = order;
	cmesh.SetDefaultOrder(order);
	
	cmesh.AutoBuild();
	//	cmesh.AdjustBoundaryElements();
	cmesh.InitializeBlock();
	
	ofstream output("outputPar.dat");
	TPZAnalysis an(&cmesh,output);
	
	TPZVec<int> numelconnected(cmesh.NEquations(),0);
	TPZSpStructMatrix mat(&cmesh);
	
	an.SetStructuralMatrix(mat);
	
	TPZStepSolver<STATE> sol;
	sol.SetJacobi(100,1.e-5,0);
	
	
	an.SetSolver(sol);
	an.Run(output);
	output.flush();
	cout.flush();
	return 0;
	
}
#endif