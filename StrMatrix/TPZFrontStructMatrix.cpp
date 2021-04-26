/**
 * @file
 * @brief Contains the implementation of the TPZFrontStructMatrix methods. 
 */

#include "TPZFrontStructMatrix.h"
#include "TPZFrontMatrix.h"
// #include "TPZFrontSym.h"
// #include "TPZFrontNonSym.h"
#include "pzsubcmesh.h"
#include "TPZElementMatrixT.h"
#include "TPZMaterial.h"
#include "TPZStackEqnStorage.h"

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.frontstructmatrix");
static TPZLogger loggerel("pz.strmatrix.element");
#endif

using namespace std;

template <class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::GetNumElConnected(TPZVec <int> &numelconnected){
	int64_t ic;
	
	this->fMesh->ComputeNodElCon();
	
	for(ic=0; ic<this->fMesh->ConnectVec().NElements(); ic++) {
		TPZConnect &cn = this->fMesh->ConnectVec()[ic];
		if(cn.HasDependency() || cn.IsCondensed()) continue;
		int64_t seqn = cn.SequenceNumber();
		if(seqn < 0) continue;
		int64_t firsteq = this->fMesh->Block().Position(seqn);
		int64_t lasteq = firsteq+this->fMesh->Block().Size(seqn);
        int64_t numactive = this->fEquationFilter.NumActive(firsteq, lasteq);
		if(!numactive) continue;
        if (numactive != lasteq-firsteq) {
            DebugStop();
        }
        TPZManVector<int64_t> firstind(numactive),destindex(numactive);
        for (int64_t i = 0; i<numactive; i++) {
            firstind[i] = i;
            destindex[i] = firsteq+i;
        }
        this->fEquationFilter.Filter(firstind, destindex);
		for(int64_t ind=0;ind<destindex.size();ind++) 
			numelconnected[destindex[ind] ] = this->fMesh->ConnectVec()[ic].NElConnected();
	}
}



template<class TFront, class TVar, class TPar>
TPZMatrix<TVar> *TPZFrontStructMatrix<TFront,TVar,TPar>::Create(){
	
    /* TPZVec <int> numelconnected(fMesh->NEquations(),0);
	 // TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym> *mat = new TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym>(fMesh->NEquations());
	 
     GetNumElConnected(numelconnected);
     //mat->SetNumElConnected(numelconnected);
     return mat;
     */
	return nullptr;
}

template<class TFront, class TVar, class TPar>
TPZStructMatrix * TPZFrontStructMatrix<TFront,TVar,TPar>::Clone(){
	
	auto* result =  new TPZFrontStructMatrix<TFront,TVar,TPar>(*this);
	return result;
}
template<class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::OrderElement()//TPZVec<int> &elorder)
{
	int64_t numelconnected = 0;
	int64_t nconnect = this->fMesh->ConnectVec().NElements();
	int64_t ic;
	//firstelconnect contains the first element index in the elconnect vector
	TPZVec<int64_t> firstelconnect(nconnect+1);
	firstelconnect[0] = 0;
	for(ic=0; ic<nconnect; ic++) {
		numelconnected += this->fMesh->ConnectVec()[ic].NElConnected();
		firstelconnect[ic+1] = firstelconnect[ic]+this->fMesh->ConnectVec()[ic].NElConnected();
	}
#ifdef PZDEBUG
    TPZVec<int64_t> firstel_copy(firstelconnect);
#endif
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
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
  	TPZVec<int64_t> elconnect(numelconnected,-1);
  	int64_t el;
  	TPZCompEl *cel;
  	for(el=0; el<this->fMesh->ElementVec().NElements(); el++) {
  		cel = this->fMesh->ElementVec()[el];
  		if(!cel) continue;
  		TPZStack<int64_t> connectlist;
  		cel->BuildConnectList(connectlist);
  		int64_t nc = connectlist.NElements();
  		int64_t ic;
  		for(ic=0; ic<nc; ic++) {
  			int64_t cindex = connectlist[ic];
#ifdef PZDEBUG
            if(firstelconnect[cindex] >= firstel_copy[cindex+1])
            {
                std::cout << "firstelconnect " << firstelconnect << std::endl;
                std::cout << "firstel_copy " << firstel_copy << std::endl;
                DebugStop();
            }
#endif
  			elconnect[firstelconnect[cindex]] = el;
  			firstelconnect[cindex]++;
  		}
  	}
	//  for(ic=0; ic<numelconnected; ic++) cout << elconnect[ic] << endl;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout<< "elconnect "<< elconnect;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  	firstelconnect[0] = 0;
  	for(ic=0; ic<nconnect; ic++) {
  		firstelconnect[ic+1] = firstelconnect[ic]+this->fMesh->ConnectVec()[ic].NElConnected();
  	}
	//cout << "elconnect\n";
	//  int no;
	for(int64_t no=0; no< this->fMesh->ConnectVec().NElements(); no++) {
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout<< "Node index " << no << ' ' << " seq num " << this->fMesh->ConnectVec()[no].SequenceNumber() << ' ';
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		//cout << "no numero " << no << ' ' << " seq num " << this->fMesh->ConnectVec()[no].SequenceNumber() << ' ';
		//       for(ic=firstelconnect[no]; ic<firstelconnect[no+1];ic++) cout << elconnect[ic] << ' ';
		//cout << endl;
	}
	
	
  	fElementOrder.Resize(this->fMesh->ElementVec().NElements(),-1);
  	fElementOrder.Fill(-1);
  	TPZVec<int> nodeorder(this->fMesh->ConnectVec().NElements(),-1);
  	firstelconnect[0] = 0;
  	for(ic=0; ic<nconnect; ic++) {
  		int64_t seqnum = this->fMesh->ConnectVec()[ic].SequenceNumber();
  		if(seqnum >= 0) nodeorder[seqnum] = ic;
  	}
	//  cout << "nodeorder ";
	/*for(ic=0; ic<this->fMesh->ConnectVec().NElements(); ic++) cout << nodeorder[ic] << ' ';
	 cout << endl;
	 cout.flush();*/
  	int64_t seq;
  	int64_t elsequence = 0;
  	TPZVec<int64_t> elorderinv(this->fMesh->ElementVec().NElements(),-1);
  	for(seq=0; seq<nconnect; seq++) {
  		ic = nodeorder[seq];
  		if(ic == -1) continue;
  		int64_t firstind = firstelconnect[ic];
  		int64_t lastind = firstelconnect[ic+1];
  		int64_t ind;
  		for(ind=firstind; ind<lastind; ind++) {
  			el = elconnect[ind];
			if(el == -1) {
				continue;
			}
  			if(elorderinv[el]==-1) elorderinv[el] = elsequence++;
  		}
  	}
	//  cout << "elorderinv ";
	//  for(seq=0;seq<this->fMesh->ElementVec().NElements();seq++) cout << elorderinv[seq] << ' ';
	//  cout << endl;
  	elsequence = 0;
  	for(seq=0;seq<this->fMesh->ElementVec().NElements();seq++) {
  		if(elorderinv[seq] == -1) continue;
  		fElementOrder[elorderinv[seq]] = seq;
  	}
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
        sout << "element order " << fElementOrder << std::endl;
		sout<< "elorderinv " << elorderinv << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//  cout << "elorder" << endl;
	//  for(ic=0; ic<this->fMesh->ElementVec().NElements(); ic++) cout << elorder[ic] << endl;
	
}

template<class TFront, class TVar, class TPar>
TPZMatrix<TVar> * TPZFrontStructMatrix<TFront,TVar,TPar>::CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPar::InitCreateAssemble();
    int64_t neq = this->fEquationFilter.NActiveEquations();
	TPZManVector <int> numelconnected(neq,0);
	TPZFrontMatrix<TVar,TPZStackEqnStorage<TVar>, TFront> *mat = new TPZFrontMatrix<TVar,TPZStackEqnStorage<TVar>, TFront>(neq);//(this->fMesh->NEquations());
	
//	TPZFrontMatrix<TVar,TPZFileEqnStorage<TVar>, TFront> *mat = new TPZFrontMatrix<TVar,TPZFileEqnStorage<TVar>, TFront>(neq);
    mat->GetFront().SetDecomposeType(fDecomposeType);
	// if the frontal matrix is applied to a submesh, we assume there may be rigid body modes
	TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (this->fMesh);	
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
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
        mat->FinishWriting();
        mat->ReOpen();
//		mat->Print("Frontal matrix", sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	mat->FinishWriting();
	mat->ReOpen();
	return mat;
}

template<class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::AssembleNew(TPZMatrix<TVar> & stiffness, TPZFMatrix<TVar> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	
	int64_t iel;
	int64_t numel = 0, nelem = this->fMesh->NElements();
	TPZElementMatrixT<TVar> ek(this->fMesh,TPZElementMatrix::EK),ef(this->fMesh,TPZElementMatrix::EF);
	TPZManVector<int64_t> destinationindex(0);
	TPZManVector<int64_t> sourceindex(0);
	
	TPZAdmChunkVector<TPZCompEl *> &elementvec = this->fMesh->ElementVec();
	
	
	/**Rearange elements order*/
	TPZVec<int> elorder(this->fMesh->NEquations(),0);
	
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
		this->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
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
			//ek.Print(*this->fMesh,cout);
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


template<class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::Assemble(TPZBaseMatrix & stiff_base, TPZBaseMatrix & rhs_base, TPZAutoPointer<TPZGuiInterface> guiInterface){
    if(!dynamic_cast<TPZMatrix<TVar>*>(&stiff_base) ||
       dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" incompatible types. Aborting...\n";
        DebugStop();
    }
	auto& stiffness = dynamic_cast<TPZMatrix<TVar>&>(stiff_base);
    auto& rhs = dynamic_cast<TPZFMatrix<TVar>&>(rhs_base);
	int64_t iel;
	int64_t numel = 0, nelem = this->fMesh->NElements();
	TPZElementMatrixT<TVar> ek(this->fMesh,TPZElementMatrix::EK),ef(this->fMesh,TPZElementMatrix::EF);
	
	TPZAdmChunkVector<TPZCompEl *> &elementvec = this->fMesh->ElementVec();
	
	
	/**Rearange elements order*/
	TPZVec<int> elorder(this->fMesh->NEquations(),0);
	
	OrderElement();
	
	for(iel=0; iel < nelem; iel++) {
		
        int64_t elindex = fElementOrder[iel];
		if(elindex < 0) continue;
		TPZCompEl *el = elementvec[elindex];
		if(!el) continue;
		TPZMaterial * mat = el->Material();
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
		}
		
		//		int dim = el->NumNodes();
		
		//Builds elements stiffness matrix
		el->CalcStiff(ek,ef);
		
		
		std::cout<< " assemblando elemento frontal " << iel <<std::endl;
		
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
		{
			std::stringstream sout;
			ek.fMat.Print("Element stiffness Frontal",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		AssembleElement(el, ek, ef, stiffness, rhs);
#ifdef PZ_LOG
		if(loggerel.isDebugEnabled())
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
template<class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::AssembleElement(TPZCompEl * el, TPZElementMatrix & ekb, TPZElementMatrix & efb, TPZMatrix<TVar> & stiffness, TPZFMatrix<TVar> & rhs){
	//TODOCOMPLEX
    auto &ek =
		dynamic_cast<TPZElementMatrixT<TVar>&>(ekb);
	auto &ef =
		dynamic_cast<TPZElementMatrixT<TVar>&>(efb);
	
	if(!el->HasDependency()) {
		//ek.fMat->Print("stiff has no constraint",test);
		//ef.fMat->Print("rhs has no constraint",test);
		//test.flush();
		ek.ComputeDestinationIndices();
		this->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
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
        this->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
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
template<class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::SetQuiet(int quiet)
{
	this->f_quiet = quiet;
}

/**
 * Resequence the connects according to the element order
 **/
template<class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::AdjustSequenceNumbering()
{
	int64_t nconnect = this->fMesh->ConnectVec().NElements();
	TPZManVector<int64_t> permute(nconnect);
	this->fMesh->ComputeNodElCon();
	int64_t i;
	for(i=0; i<nconnect; i++)
	{
		permute[i] = i;
	}
	TPZCompEl *cel;
	int64_t nelem = fElementOrder.NElements();
	int64_t el;
	int64_t connectcount = 0;
	for(i=0; i<nelem; i++)
	{
		el = fElementOrder[i];
        if(el<0) continue;
		cel = this->fMesh->ElementVec()[el];
		if(!cel) continue;
		std::set<int64_t> indepconnects, depconnects;
		cel->BuildConnectList(indepconnects,depconnects);
		std::set<int64_t>::iterator it;
		for(it=indepconnects.begin(); it != indepconnects.end(); it++)
		{
			TPZConnect &nod = this->fMesh->ConnectVec()[*it];
			int nelcon = nod.NElConnected()-1;
			int64_t seqnum = nod.SequenceNumber();
			if(nelcon == 0) permute[seqnum]= connectcount++;
			nod.DecrementElConnected();
		}
	}
	for(i=0; i<nconnect; i++)
	{
		TPZConnect &nod = this->fMesh->ConnectVec()[i];
		if(nod.SequenceNumber() < 0 || nod.NElConnected() <= 0) continue;
		if(permute[nod.SequenceNumber()] < connectcount)
		{
			std::cout << __PRETTY_FUNCTION__ << " very fishy\n";
			DebugStop();
		}
	}
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " permutation " << permute;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	this->fMesh->Permute(permute);
}

template<class TFront, class TVar, class TPar>
int TPZFrontStructMatrix<TFront,TVar,TPar>::ClassId() const{
    return Hash("TPZFrontStructMatrix") ^
        (TPZStructMatrixT<TVar>::ClassId()<< 1 ) ^
        TFront().ClassId() << 2 ^
        TPar::ClassId() << 3;
}

template<class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}
template<class TFront, class TVar, class TPar>
void TPZFrontStructMatrix<TFront,TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZFrontStructMatrix<TPZFrontSym<STATE>,STATE,TPZStructMatrixOR<STATE>>;
template class TPZFrontStructMatrix<TPZFrontNonSym<STATE>,STATE,TPZStructMatrixOR<STATE>>;
template class TPZFrontStructMatrix<TPZFrontSym<STATE>,STATE,TPZStructMatrixOT<STATE>>;
template class TPZFrontStructMatrix<TPZFrontNonSym<STATE>,STATE,TPZStructMatrixOT<STATE>>;
template class TPZFrontStructMatrix<TPZFrontSym<STATE>,STATE,TPZStructMatrixTBBFlow<STATE>>;
template class TPZFrontStructMatrix<TPZFrontNonSym<STATE>,STATE,TPZStructMatrixTBBFlow<STATE>>;


template class TPZFrontStructMatrix<TPZFrontSym<CSTATE>,CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZFrontStructMatrix<TPZFrontNonSym<CSTATE>,CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZFrontStructMatrix<TPZFrontSym<CSTATE>,CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZFrontStructMatrix<TPZFrontNonSym<CSTATE>,CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZFrontStructMatrix<TPZFrontSym<CSTATE>,CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;
template class TPZFrontStructMatrix<TPZFrontNonSym<CSTATE>,CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;
