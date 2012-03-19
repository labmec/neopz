/**
 * @file
 * @brief Contains the implementation of the TPZFrontStructMatrix methods. 
 */

#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontStructMatrix.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzconnect.h"
#include "pzadmchunk.h"

#include "pzsmfrontalanal.h"

#include "pzmat2dlin.h"

#include "pzanalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "TPZFrontMatrix.h"

#include "pzdxmesh.h"
#include <fstream>

using namespace std;

#include "pzelmat.h"
#include "pzbndcond.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.frontstructmatrix"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
#endif


template <class front>
void TPZFrontStructMatrix<front>::GetNumElConnected(TPZVec <int> &numelconnected){
	int ic;
	
	fMesh->ComputeNodElCon();
	
	for(ic=0; ic<fMesh->ConnectVec().NElements(); ic++) {
		TPZConnect &cn = fMesh->ConnectVec()[ic];
		if(cn.HasDependency() || cn.IsCondensed()) continue;
		int seqn = cn.SequenceNumber();
		if(seqn < 0) continue;
		int firsteq = fMesh->Block().Position(seqn);
		int lasteq = firsteq+fMesh->Block().Size(seqn);
		if(firsteq < fMinEq || firsteq >= fMaxEq) continue;
		int ind;
		for(ind=firsteq;ind<lasteq;ind++) numelconnected[ind-fMinEq] = fMesh->ConnectVec()[ic].NElConnected();
	}
}

template<class front>
TPZFrontStructMatrix<front>::TPZFrontStructMatrix(TPZCompMesh *mesh): TPZStructMatrix(mesh) {
	f_quiet = 0;
}


template<class front>
TPZFrontStructMatrix<front>::~TPZFrontStructMatrix(){}



template<class front>
TPZMatrix<REAL> * TPZFrontStructMatrix<front>::Create(){
	
    /* TPZVec <int> numelconnected(fMesh->NEquations(),0);
	 // TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym> *mat = new TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym>(fMesh->NEquations());
	 
     GetNumElConnected(numelconnected);
     //mat->SetNumElConnected(numelconnected);
     return mat;
     */
	return 0;
}

template<class front>
TPZStructMatrix * TPZFrontStructMatrix<front>::Clone(){
	
	TPZFrontStructMatrix<front>* result =  new TPZFrontStructMatrix<front>(*this);
	return result;
}
template<class front>
void TPZFrontStructMatrix<front>::OrderElement()//TPZVec<int> &elorder)
{
	int numelconnected = 0;
	int nconnect = fMesh->ConnectVec().NElements();
	int ic;
	//firstelconnect contains the first element index in the elconnect vector
	TPZVec<int> firstelconnect(nconnect+1);
	firstelconnect[0] = 0;
	for(ic=0; ic<nconnect; ic++) {
		numelconnected += fMesh->ConnectVec()[ic].NElConnected();
		firstelconnect[ic+1] = firstelconnect[ic]+fMesh->ConnectVec()[ic].NElConnected();
	}
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout<<"numelconnected " << numelconnected << endl;
		sout<< "firstelconnect "<< firstelconnect;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//cout << "numelconnected " << numelconnected << endl;
	//cout << "firstelconnect ";
	//  for(ic=0; ic<nconnect; ic++) cout << firstelconnect[ic] << ' ';
  	TPZVec<int> elconnect(numelconnected,-1);
  	int el;
  	TPZCompEl *cel;
  	for(el=0; el<fMesh->ElementVec().NElements(); el++) {
  		cel = fMesh->ElementVec()[el];
  		if(!cel) continue;
  		TPZStack<int> connectlist;
  		cel->BuildConnectList(connectlist);
  		int nc = connectlist.NElements();
  		int ic;
  		for(ic=0; ic<nc; ic++) {
  			int cindex = connectlist[ic];
  			elconnect[firstelconnect[cindex]] = el;
  			firstelconnect[cindex]++;
  		}
  	}
	//  for(ic=0; ic<numelconnected; ic++) cout << elconnect[ic] << endl;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout<< "elconnect "<< elconnect;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  	firstelconnect[0] = 0;
  	for(ic=0; ic<nconnect; ic++) {
  		firstelconnect[ic+1] = firstelconnect[ic]+fMesh->ConnectVec()[ic].NElConnected();
  	}
	//cout << "elconnect\n";
	//  int no;
	for(int no=0; no< fMesh->ConnectVec().NElements(); no++) {
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout<< "no numero " << no << ' ' << " seq num " << fMesh->ConnectVec()[no].SequenceNumber() << ' ';
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		//cout << "no numero " << no << ' ' << " seq num " << fMesh->ConnectVec()[no].SequenceNumber() << ' ';
		//       for(ic=firstelconnect[no]; ic<firstelconnect[no+1];ic++) cout << elconnect[ic] << ' ';
		//cout << endl;
	}
	
	
  	fElementOrder.Resize(fMesh->ElementVec().NElements(),-1);
  	fElementOrder.Fill(-1);
  	TPZVec<int> nodeorder(fMesh->ConnectVec().NElements(),-1);
  	firstelconnect[0] = 0;
  	for(ic=0; ic<nconnect; ic++) {
  		int seqnum = fMesh->ConnectVec()[ic].SequenceNumber();
  		if(seqnum >= 0) nodeorder[seqnum] = ic;
  	}
	//  cout << "nodeorder ";
	/*for(ic=0; ic<fMesh->ConnectVec().NElements(); ic++) cout << nodeorder[ic] << ' ';
	 cout << endl;
	 cout.flush();*/
  	int seq;
  	int elsequence = 0;
  	TPZVec<int> elorderinv(fMesh->ElementVec().NElements(),-1);
  	for(seq=0; seq<nconnect; seq++) {
  		ic = nodeorder[seq];
  		if(ic == -1) continue;
  		int firstind = firstelconnect[ic];
  		int lastind = firstelconnect[ic+1];
  		int ind;
  		for(ind=firstind; ind<lastind; ind++) {
  			el = elconnect[ind];
			if(el == -1) {
				continue;
			}
  			if(elorderinv[el]==-1) elorderinv[el] = elsequence++;
  		}
  	}
	//  cout << "elorderinv ";
	//  for(seq=0;seq<fMesh->ElementVec().NElements();seq++) cout << elorderinv[seq] << ' ';
	//  cout << endl;
  	elsequence = 0;
  	for(seq=0;seq<fMesh->ElementVec().NElements();seq++) {
  		if(elorderinv[seq] == -1) continue;
  		fElementOrder[elorderinv[seq]] = seq;
  	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout<< "elorderinv " << elorderinv<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//  cout << "elorder" << endl;
	//  for(ic=0; ic<fMesh->ElementVec().NElements(); ic++) cout << elorder[ic] << endl;
	
}

template<class front>
TPZMatrix<REAL> * TPZFrontStructMatrix<front>::CreateAssemble(TPZFMatrix<REAL> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
	
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
	TPZManVector <int> numelconnected(neq,0);
	//TPZFrontMatrix<TPZStackEqnStorage, front> *mat = new TPZFrontMatrix<TPZStackEqnStorage, front>(neq);//(fMesh->NEquations());
	
	TPZFrontMatrix<TPZFileEqnStorage, front> *mat = new TPZFrontMatrix<TPZFileEqnStorage, front>(neq);
	
	// if the frontal matrix is applied to a submesh, we assume there may be rigid body modes
	TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (fMesh);	
	if (subcmesh) {
		int nrigid = subcmesh->NumberRigidBodyModes();
		if (nrigid > 0) {
			mat->GetFront().SetNumRigidBodyModes(nrigid);
			
		}
	}
	
	GetNumElConnected(numelconnected);
	mat->SetNumElConnected(numelconnected);
	
	OrderElement();
	
	Assemble(*mat,rhs,guiInterface);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
        mat->FinishWriting();
        mat->ReOpen();
		//mat->Print("Frontal matrix", sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	mat->FinishWriting();
	mat->ReOpen();
	return mat;
}

template<class front>
void TPZFrontStructMatrix<front>::AssembleNew(TPZMatrix<REAL> & stiffness, TPZFMatrix<REAL> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	
	int iel;
	int numel = 0, nelem = fMesh->NElements();
	TPZElementMatrix ek(fMesh,TPZElementMatrix::EK),ef(fMesh,TPZElementMatrix::EF);
	TPZManVector<int> destinationindex(0);
	TPZManVector<int> sourceindex(0);
	
	TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
	
	
	/**Rearange elements order*/
	TPZVec<int> elorder(fMesh->NEquations(),0);
	
	OrderElement();
	
	
	for(iel=0; iel < nelem; iel++) {
		
		if(guiInterface) if(guiInterface->AmIKilled()){
			break;
		}
		
		if(fElementOrder[iel] < 0) continue;
		TPZCompEl *el = elementvec[fElementOrder[iel]];
		if(!el) continue;
		//		int dim = el->NumNodes();
		
		//Builds elements stiffness matrix
		el->CalcStiff(ek,ef);
		ek.ComputeDestinationIndices();
		FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,fMinEq,fMaxEq);
		//ek.fMat->Print(out);
		//ef.fMat->Print();
		if(!f_quiet)
		{
			if(!(numel%20)) cout << endl << numel;
			//    if(!(numel%20)) cout << endl;
			cout << '*';
			cout.flush();
		}
		numel++;
		
		if(!el->HasDependency()) {
			stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
			rhs.AddFel(ef.fMat,ek.fSourceIndex, ek.fDestinationIndex);
		}
		else {
			//ek.Print(*fMesh,cout);
			ek.ApplyConstraints();
			ef.ApplyConstraints();
			stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
			rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
			/*
			 if(ek.fConstrMat->Decompose_LU() != -1) {
			 el->ApplyConstraints(ek,ef);
			 ek.Print(*this,check);
			 check.flush();
			 }
			 */
		}
		
	}//fim for iel
	if(!f_quiet)
	{
		cout << endl;
	}
}


template<class front>
void TPZFrontStructMatrix<front>::Assemble(TPZMatrix<REAL> & stiffness, TPZFMatrix<REAL> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
	
	int iel;
	int numel = 0, nelem = fMesh->NElements();
	TPZElementMatrix ek(fMesh,TPZElementMatrix::EK),ef(fMesh,TPZElementMatrix::EF);
	
	TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
	
	
	/**Rearange elements order*/
	TPZVec<int> elorder(fMesh->NEquations(),0);
	
	OrderElement();
	
	for(iel=0; iel < nelem; iel++) {
		
		if(fElementOrder[iel] < 0) continue;
		TPZCompEl *el = elementvec[fElementOrder[iel]];
		if(!el) continue;
		TPZAutoPointer<TPZMaterial> mat = el->Material();
		if(mat)
		{
			int matid = mat->Id();
			if(this->ShouldCompute(matid) == false)
			{
				continue;
			}
		}
		else
		{
			TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
			if(!submesh)
			{
				continue;
			}
		}
		
		//		int dim = el->NumNodes();
		
		//Builds elements stiffness matrix
		el->CalcStiff(ek,ef);
		
		
		std::cout<< " assemblando elemento frontal " << iel <<std::endl;
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			ek.fMat.Print("Element stiffness Frontal",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		AssembleElement(el, ek, ef, stiffness, rhs);
#ifdef LOG4CXX
		if(loggerel->isDebugEnabled())
		{
			std::stringstream sout;
			ek.fMat.Print("Element stiffness depois de assemblada ",sout);
			LOGPZ_DEBUG(loggerel,sout.str())
		}
#endif
		
		if(!f_quiet)
		{
			cout << '*';
			if(!(numel%20)) {
				cout << " " << (100*iel/nelem) << "% Elements assembled " << endl;
				cout.flush();
			}
		}
		numel++;
		
	}//fim for iel
	
}

//Verificar declaracao dos parametros !!!!!
template<class front>
void TPZFrontStructMatrix<front>::AssembleElement(TPZCompEl * el, TPZElementMatrix & ek, TPZElementMatrix & ef, TPZMatrix<REAL> & stiffness, TPZFMatrix<REAL> & rhs){
	
	
	if(!el->HasDependency()) {
		//ek.fMat->Print("stiff has no constraint",test);
		//ef.fMat->Print("rhs has no constraint",test);
		//test.flush();
		ek.ComputeDestinationIndices();
		this->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,fMinEq,fMaxEq);
#ifdef LOG4CXX
		{
            std::stringstream sout;
            sout << "Element index " << el->Index() << " Unconstrained destination index " << ek.fDestinationIndex;
            LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		//ek.Print(*fMesh,cout);
		stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
		rhs.AddFel(ef.fMat,ek.fSourceIndex, ek.fDestinationIndex);                 //  ??????????? Erro
	}
	else
	{
        ek.ApplyConstraints();
        ef.ApplyConstraints();
        ek.ComputeDestinationIndices();
        FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,fMinEq,fMaxEq);
#ifdef LOG4CXX
        {
			std::stringstream sout;
			sout << "Element index " << el->Index() << " Constrained destination index " << ek.fDestinationIndex;
			LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
	}
}



/*!
 \fn TPZFrontStructMatrix::SetQuiet(int quiet)
 */
template<class front>
void TPZFrontStructMatrix<front>::SetQuiet(int quiet)
{
	this->f_quiet = quiet;
}

template<class front>
int TPZFrontStructMatrix<front>::main() {
	int refine = 5;
	int order = 2;
	
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
		/*gel = */gmesh.CreateGeoElement(EQuadrilateral, indices, 1, index);
	}
	gmesh.BuildConnectivity ();
	
	TPZVec<TPZGeoEl *> subel;
	//gel->Divide(subel);
	
	
	
	cout << "Refinement ";
	cin >> refine;
	cout << endl;
	UniformRefine(refine,gmesh);
	
	
	
	TPZMat2dLin *meumat = new TPZMat2dLin(1);
	TPZFMatrix<REAL> xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
	meumat->SetMaterial (xk,xc,xf);
	TPZAutoPointer<TPZMaterial> meumatptr(meumat);
	cmesh.InsertMaterialObject(meumatptr);
	
	TPZFMatrix<REAL> val1(1,1,0.),val2(1,1,0.);
	TPZAutoPointer<TPZMaterial> bnd = meumat->CreateBC (meumatptr,-4,0,val1,val2);
	cmesh.InsertMaterialObject(bnd);
	
	
	
	cout << "Interpolation order ";
	cin >> order;
	cout << endl;
	
	//	TPZCompEl::gOrder = order;
	cmesh.SetDefaultOrder(order);
	
	cmesh.AutoBuild();
	//	cmesh.AdjustBoundaryElements();
	cmesh.InitializeBlock();
	
	ofstream output("outputNon.dat");
	//	ofstream output2("outputNon.dat");
	cmesh.Print(output);
	TPZAnalysis an(&cmesh,output);
	//	TPZAnalysis an2(&cmesh,output);
	
	TPZVec<int> numelconnected(cmesh.NEquations(),0);
	int ic;
	cout << "N�mero de Equa��es -> " << cmesh.NEquations() << endl;
	cout.flush();
	
	ofstream out("cmeshBlock_out.txt");
	//	cmesh.Print(out);
	//	cmesh.Block().Print("Block",out);
	for(ic=0; ic<cmesh.ConnectVec().NElements(); ic++) {
		TPZConnect &cn = cmesh.ConnectVec()[ic];
		if(cn.HasDependency()) continue;
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
	//	 << "nequations " << numelconnected.NElements();
	//	for(ic=0;ic<numelconnected.NElements(); ic++) cout << numelconnected[ic] <<' ';
	//	cout << endl;
	//	cout.flush();
	
	//	TPZFrontMatrix<TPZFileEqnStorage, TPZFrontNonSym> *mat = new TPZFrontMatrix<TPZFileEqnStorage, TPZFrontNonSym>(cmesh.NEquations());
	//TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym> *mat = new TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym>(cmesh.NEquations());
	//TPZFrontMatrix<TPZStackEqnStorage> *mat = new TPZFrontMatrix<TPZStackEqnStorage>(cmesh.NEquations());
	
	TPZFrontStructMatrix<TPZFrontSym> mat(&cmesh);
	
	//   TPZFStructMatrix mat2(&cmesh);
	//  mat->SetNumElConnected(numelconnected);
	//mat = CreateAssemble();
	
	
	an.SetStructuralMatrix(mat);
	//	an2.SetStructuralMatrix(mat2);
	
	TPZStepSolver<REAL> sol;
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
	//	cout << "******************************************************************************************************AQUI 1" << endl;
	an.Run(output);
	an.Print("solution of frontal solver", output);
	//	cout << "******************************************************************************************************AQUI 2" << endl;
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

/**
 * Resequence the connects according to the element order
 **/
template<class front>
void TPZFrontStructMatrix<front>::AdjustSequenceNumbering()
{
	int nconnect = this->fMesh->ConnectVec().NElements();
	TPZManVector<int> permute(nconnect);
	fMesh->ComputeNodElCon();
	int i;
	for(i=0; i<nconnect; i++)
	{
		permute[i] = i;
	}
	TPZCompEl *cel;
	int nelem = fElementOrder.NElements();
	int el;
	int connectcount = 0;
	for(i=0; i<nelem; i++)
	{
		el = fElementOrder[i];
		cel = fMesh->ElementVec()[el];
		if(!cel) continue;
		std::set<int> indepconnects, depconnects;
		cel->BuildConnectList(indepconnects,depconnects);
		std::set<int>::iterator it;
		for(it=indepconnects.begin(); it != indepconnects.end(); it++)
		{
			TPZConnect &nod = fMesh->ConnectVec()[*it];
			int nelcon = nod.NElConnected()-1;
			int seqnum = nod.SequenceNumber();
			if(nelcon == 0) permute[seqnum]= connectcount++;
			nod.DecrementElConnected();
		}
	}
	for(i=0; i<nconnect; i++)
	{
		TPZConnect &nod = fMesh->ConnectVec()[i];
		if(nod.SequenceNumber() < 0 || nod.NElConnected() <= 0) continue;
		if(permute[nod.SequenceNumber()] < connectcount)
		{
			std::cout << __PRETTY_FUNCTION__ << " very fishy\n";
			DebugStop();
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " permutation " << permute;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fMesh->Permute(permute);
}

class TPZFrontSym;
class TPZFrontNonSym;

template class TPZFrontStructMatrix<TPZFrontSym>;
template class TPZFrontStructMatrix<TPZFrontNonSym>;
