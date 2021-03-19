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
#include "pzbndcond.h"
#include "TPZTimer.h"

#include "pzlog.h"
#ifdef LOG4CXX
static PZLogger logger("pz.StrMatrix");
#endif

using namespace std;

TPZStructMatrix * TPZSpStructMatrix::Clone(){
    return new TPZSpStructMatrix(*this);
}
TPZMatrix<STATE> * TPZSpStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,
                                              TPZAutoPointer<TPZGuiInterface> guiInterface){

#ifdef LOG4CXX
    if(logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble starting")
    }
#endif
	
    int64_t neq = fMesh->NEquations();
    if(fMesh->FatherMesh()) {
		cout << "TPZSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
		return new TPZFYsmpMatrix<STATE>(0,0);
    }
    TPZMatrix<STATE> *stiff = Create();//new TPZFYsmpMatrix(neq,neq);
    TPZFYsmpMatrix<STATE> *mat = dynamic_cast<TPZFYsmpMatrix<STATE> *> (stiff);
    rhs.Redim(neq,1);
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    TPZTimer before("Assembly of a sparse matrix");
    before.start();
#ifdef LOG4CXX
    if (logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble calling Assemble()");
    }
#endif
    Assemble(*stiff,rhs,guiInterface);
    before.stop();
    std::cout << __PRETTY_FUNCTION__ << " " << before << std::endl;
//    mat->ComputeDiagonal();
    //    mat->ComputeDiagonal();
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
#ifdef LOG4CXX
    if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble exiting");
#endif
    return stiff;
}
TPZMatrix<STATE> * TPZSpStructMatrix::Create(){
    int64_t neq = fEquationFilter.NActiveEquations();
	/*    if(fMesh->FatherMesh()) {
	 TPZSubCompMesh *smesh = (TPZSubCompMesh *) fMesh;
	 neq = smesh->NumInternalEquations();
	 }*/
    TPZFYsmpMatrix<STATE> * mat = new TPZFYsmpMatrix<STATE>(neq,neq);
	
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
    
#ifdef LOG4CXX2
    if(logger.isDebugEnabled()){
        std::stringstream sout;
        sout << "Node graph \n";
        metis.TPZRenumbering::Print(nodegraph, nodegraphindex);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    /**vector sizes*/
    int64_t nblock = nodegraphindex.NElements()-1;
    // number of values in the sparse matrix
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(int64_t i=0;i<nblock;i++){
		int64_t iblsize = fMesh->Block().Size(i);
		int64_t iblpos = fMesh->Block().Position(i);
        int64_t numactive = fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
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
	
#ifdef LOG4CXX
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Number of equations " << totaleq << " number of nonzero s " << totalvar;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
	
    nblock=fMesh->NIndependentConnects();
	
    TPZManVector<int64_t,400> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<STATE> EqValue(totalvar);
    for(int64_t i=0;i<nblock;i++){
		int64_t iblsize = fMesh->Block().Size(i);
		int64_t iblpos = fMesh->Block().Position(i);
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t ij=0; ij<iblsize; ij++) {
            rowdestindices[ij] = iblpos+ij;
        }
        fEquationFilter.Filter(rowdestindices);

		int64_t ibleq;
        // working equation by equation
        // rowdestindices contains the equation number of each element in the block number "i"
		for(ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            int rowind = rowdestindices[ibleq];
//            if (rowind != pos) {
//                DebugStop();
//            }
			Eq[ieq] = pos;
			int64_t colsize,colpos,jbleq;
			int64_t diagonalinsert = 0;
			int64_t icfirst = nodegraphindex[i];
			int64_t iclast = nodegraphindex[i+1];
			int64_t j;
			for(j=icfirst;j<iclast;j++)
            {
                // col is the block linked to block "i"
				int64_t col = nodegraph[j];
                
                // force the diagonal block to be inserted
                // the nodegraph does not contain the pointer to itself
				if(!diagonalinsert && col > i)
				{
					diagonalinsert = 1;
					int64_t colsize = fMesh->Block().Size(i);
					int64_t colpos = fMesh->Block().Position(i);
                    TPZManVector<int64_t> destindices(colsize);
                    for (int64_t i=0; i<colsize; i++) {
                        destindices[i] = colpos+i;
                    }
                    fEquationFilter.Filter(destindices);
					int64_t jbleq;
					for(jbleq=0; jbleq<destindices.size(); jbleq++) {
						//             if(colpos+jbleq == ieq) continue;
						EqCol[pos] = destindices[jbleq];
						EqValue[pos] = 0.;
						//            colpos++;
                        // pos is the position within EqCol or EqVal where we will assemble
						pos++;
					}
				}
				colsize = fMesh->Block().Size(col);
				colpos = fMesh->Block().Position(col);
                // optimization statement : if all equations in the range are inactive -> continue
                if (fEquationFilter.NumActive(colpos, colpos+colsize) == 0) {
                    continue;
                }
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
				for(jbleq=0; jbleq<destindices.size(); jbleq++) {
					EqCol[pos] = destindices[jbleq];
					EqValue[pos] = 0.;
					colpos++;
					pos++;
				}
			}
            // all elements are below (last block certainly)
			if(!diagonalinsert)
			{
				diagonalinsert = 1;
				int64_t colsize = fMesh->Block().Size(i);
				int64_t colpos = fMesh->Block().Position(i);
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
				int64_t jbleq;
				for(jbleq=0; jbleq<destindices.size(); jbleq++) {
					//             if(colpos+jbleq == ieq) continue;
					EqCol[pos] = destindices[jbleq];
					EqValue[pos] = 0.;
					//            colpos++;
					pos++;
				}
			}
			ieq++;
		}
    }
    Eq[ieq] = pos;
    if(pos != totalvar)
    {
        DebugStop();
    }
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}

TPZSpStructMatrix::TPZSpStructMatrix() : TPZStructMatrix(){
}

TPZSpStructMatrix::TPZSpStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{}

int TPZSpStructMatrix::ClassId() const{
    return Hash("TPZSpStructMatrix") ^ TPZStructMatrix::ClassId() << 1;
}

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
		TPZVec<int64_t> indices(4);
		for(i=0; i<4; i++) indices[i] = i;
		// O proprio construtor vai inserir o elemento na malha
		//       gel = new TPZGeoElQ2d(el,indices,1,gmesh);
		int64_t index;
		gel = gmesh.CreateGeoElement(EQuadrilateral,indices,1,index);
	}
	gmesh.BuildConnectivity ();
	
	TPZVec<TPZGeoEl *> subel;
	
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
	
	//TPZCompEl::gOrder = order;
	cmesh.SetDefaultOrder(order);
	
	cmesh.AutoBuild();
	//	cmesh.AdjustBoundaryElements();
	cmesh.InitializeBlock();
	
	ofstream output("outputPar.dat");
	TPZAnalysis an(&cmesh,true,output);
	
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
