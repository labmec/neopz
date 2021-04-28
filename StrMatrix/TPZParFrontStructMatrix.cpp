/**
 * @file
 * @brief Contains the implementation of the TPZParFrontStructMatrix methods. 
 */

#include "TPZParFrontStructMatrix.h"
#include "TPZFrontMatrix.h"

#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"
#include "TPZParFrontMatrix.h"
#include "pzsubcmesh.h"
#include "TPZMaterial.h"
#include "TPZElementMatrixT.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZFileEqnStorage.h"


#include <thread>             // std::thread
#include <mutex>              // std::mutex, std::unique_lock
#include <condition_variable> // std::condition_variable


#ifdef PZ_LOG
#include "pzlog.h"
static TPZLogger logger("pz.strmatrix.frontstructmatrix");

#endif
/// Semaphore which controls threads assembling elements
std::mutex mutex_element_assemble;
/// Semaphore which controls thread assembling global matrices
std::mutex mutex_global_assemble;
/// Semaphore which controls condensed assembling
std::condition_variable condassemble;
/// Semaphore
std::condition_variable stackfull;

template<class TFront, class TVar, class TPar>
TPZParFrontStructMatrix<TFront,TVar,TPar>::TPZParFrontStructMatrix(TPZCompMesh *mesh): TPZFrontStructMatrix<TFront,TVar,TPar>(mesh)
{
	fMaxStackSize = 500;
	TPZStructMatrix::SetNumThreads(3);
}

template<class TFront, class TVar, class TPar>
TPZParFrontStructMatrix<TFront,TVar,TPar>::~TPZParFrontStructMatrix()
{

}

template<class TFront, class TVar, class TPar>
TPZStructMatrix * TPZParFrontStructMatrix<TFront,TVar,TPar>::Clone(){
	TPZParFrontStructMatrix<TFront,TVar,TPar> * mat = new TPZParFrontStructMatrix<TFront,TVar,TPar>(*this);
	return mat;
	;
	
}

template<class TFront, class TVar, class TPar>
void *TPZParFrontStructMatrix<TFront,TVar,TPar>::ElementAssemble(void *t){
    
    
	TPZParFrontStructMatrix<TFront,TVar,TPar> *parfront = (TPZParFrontStructMatrix<TFront,TVar,TPar> *) t;
	
	TPZAdmChunkVector<TPZCompEl *> &elementvec = parfront->fMesh->ElementVec();
	
	
	while(parfront->fCurrentElement < parfront->fNElements) {
		
		/**
		 *Lock mutex and search for an avilable element
		 *A global variable to be updated whenever a element is processed
		 */
		
		//Lock a mutex and get an element number
        std::unique_lock<std::mutex> assemble_lock(mutex_element_assemble);
//        assemble_lock.lock();

		//Stack is full and process must wait here!
		if(parfront->felnum.NElements()==parfront->fMaxStackSize)
        {
			/*          cout << "    Stack full" << endl;
			 cout << "    Waiting" << endl;
			 cout.flush();*/
			//cout << "Mutex unlocked on Condwait" << endl;
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
			{
				std::stringstream sout;
				sout << "Entering cond_wait because of stack overflow ";
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
            stackfull.wait(assemble_lock);
			//cout << "Mutex LOCKED leaving Condwait" << endl;
			
        }
		
		//cout << "Locking mutex_element_assemble" << endl;
		//cout.flush();
		int64_t local_element = parfront->fCurrentElement;
		if(local_element==parfront->fNElements)
        {
            assemble_lock.unlock();
            return 0;
        }
		/*          cout << "All element matrices assembled" << endl;
		 return 0;
		 }
		 */
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout << "Computing element " << local_element;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		parfront->fCurrentElement++;
		//Unlock mutex and / or send a broadcast
		//cout << "Computing Element " << parfront->fCurrentElement << endl;
		//cout << "Unlocking mutex_element_assemble" << endl;
		//cout.flush();
        assemble_lock.unlock();

		
		if(parfront->fElementOrder[local_element] < 0) continue;
		TPZCompEl *el = elementvec[parfront->fElementOrder[local_element]];
		if(!el) continue;
		//		int dim = el->NumNodes();
		
		//Builds elements stiffness matrix
		auto *ek =
            new TPZElementMatrixT<TVar>(parfront->fMesh,TPZElementMatrix::EK);
		auto *ef =
            new TPZElementMatrixT<TVar>(parfront->fMesh,TPZElementMatrix::EF);
		
		el->CalcStiff(*ek, *ef);
		//Locks a mutex and adds element contribution to frontmatrix
		//if mutex is locked go to condwait waiting for an specific condvariable
		// este mutex deve ser outro mutex -> mutexassemble
        {
            std::lock_guard<std::mutex> global_ass(mutex_global_assemble);
            //cout << "Locking mutex_global_assemble" << endl;
            //cout << "Pushing variables to the stack" << endl;
            //cout.flush();
            
            // colocar ek, ef e o element_local no stack
            parfront->felnum.Push(local_element);
            parfront->fekstack.Push(ek);
            parfront->fefstack.Push(ef);
            
    #ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Pushing element " << local_element << " on the stack, stack sizes " << parfront->felnum.NElements() << " " << parfront->fekstack.NElements() << " " << parfront->fefstack.NElements();
                LOGPZ_DEBUG(logger,sout.str())
            }
    #endif
            // Outro thread procura se stack contem o proximo elemento
            // qdo nao encontra entra em condwait de acordo com condassemble
            // isso ocorre num outro processo
            // Uma vez que uma nova ek foi adicionada ao stack
            // chame broadcast para acordar o thread que faz assemblagem global
            
            //cout << "Unlocking mutex_global_assemble" << endl;
            //cout << "Broadcasting condassemble" << endl;
            //cout.flush();
            
            /*     if(!(parfront->fCurrentElement%20)){
             cout << endl << "Computing " << parfront->fCurrentElement << " on thread " << std::this_thread::get_id() << endl;
             cout << " " << (100*parfront->fCurrentElement/parfront->fNElements) << "% Elements computed" << "     " << (100*parfront->fCurrentAssembled/parfront->fNElements) << "% Elements assembled" << endl;
             }
             cout << '*';
             cout.flush();
             */
            //Alterado cond_broadcast para cond_signal
            //invertendo a sequencia das chamadas
            condassemble.notify_all();
        }
		// o thread de assemblagem utiliza mutexassemble
		// e feito em outro thread     AssembleElement(el, ek, ef, stiffness, rhs);
		
		if(parfront->fGuiInterface) if(parfront->fGuiInterface->AmIKilled()){
			break;
		}
	}//fim for iel
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Falling through";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	std::cout << __PRETTY_FUNCTION__ << " Falling through \n";
    
	return NULL;
	
}

#include "arglib.h"

clarg::argInt num_threads("-ntdec", "Number of threads to decompose in TPZParFrontStructMatrix.", 6);

template<class TFront, class TVar, class TPar>
void *TPZParFrontStructMatrix<TFront,TVar,TPar>::GlobalAssemble(void *t){
    
	TPZParFrontStructMatrix<TFront,TVar,TPar> *parfront = (TPZParFrontStructMatrix<TFront,TVar,TPar> *) t;
	TPZAdmChunkVector<TPZCompEl *> &elementvec = parfront->fMesh->ElementVec();
	while(parfront->fCurrentAssembled < parfront->fNElements) {
		
#ifndef USING_ATLAS
        const int unit = (parfront->fNElements/200 == 0 ? 1 : parfront->fNElements/200);
        if(!(parfront->fCurrentAssembled%unit))std::cout << "*";
		std::cout.flush();
		if(!(parfront->fCurrentAssembled%(20*unit)) && parfront->fCurrentAssembled)
        {
			if(parfront->fCurrentElement!=parfront->fNElements){
				std::cout << " Element " << parfront->fCurrentElement << " " << (100*parfront->fCurrentElement/parfront->fNElements) << "% Elements computed " << (100*parfront->fCurrentAssembled/parfront->fNElements) << "% Elements assembled " << std::endl;
				std::cout.flush();
			}else{
				std::cout << " " << (100*parfront->fCurrentAssembled/parfront->fNElements) << "% Elements assembled " << std::endl;
				std::cout.flush();
			}
		}
#endif
		/**
		 *Lock mutex and search for an available element
		 *A global variable to be updated whenever a element is processed
		 */
		
		//Lock a mutex and get an element number
		int64_t local_element = parfront->fCurrentAssembled;
		parfront->fCurrentAssembled++;
		
		if(parfront->fElementOrder[local_element] < 0) continue;
		TPZCompEl *el = elementvec[parfront->fElementOrder[local_element]];
		if(!el) continue;
		TPZMaterial * mat = el->Material();
		if(mat)
		{
			int matid = mat->Id();
			if(parfront->ShouldCompute(matid) == false)
			{
				continue;
			}
		}
		else 
		{
			TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
            TPZElementGroup *elgrp = dynamic_cast<TPZElementGroup *>(el);
            TPZCondensedCompEl *condel = dynamic_cast<TPZCondensedCompEl *>(el);
			if(!submesh && ! elgrp && ! condel)
			{
				continue;
			}
		}
		
		//Searches for next element
		int64_t i=0;
		int64_t aux = -1;
		TPZElementMatrix *ekaux = 0, *efaux = 0;
        {
            std::unique_lock<std::mutex> global_lock(mutex_global_assemble);
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Acquired mutex_global_assemble";
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            while(aux != local_element){
                while(i < parfront->felnum.NElements()) {
                    if(parfront->felnum[i] == local_element){
                        TPZElementMatrix *ektemp, *eftemp;
                        aux = parfront->felnum[i];
                        ekaux = parfront->fekstack[i];
                        efaux = parfront->fefstack[i];
                        int64_t itemp = parfront->felnum.Pop();
                        ektemp = parfront->fekstack.Pop();
                        eftemp = parfront->fefstack.Pop();
                        if(parfront->felnum.NElements()<parfront->fMaxStackSize){
                            stackfull.notify_all();
                        }
                        if(i < parfront->felnum.NElements()) {
                            
                            parfront->felnum[i] = itemp;
                            parfront->fekstack[i]=ektemp;
                            parfront->fefstack[i]=eftemp;
                        }
                        break;
                    }
                    i++;
                }
                if(aux!=local_element){
                    i=0;
    #ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Waiting on condassemble";
                        LOGPZ_DEBUG(logger,sout.str())
                    }
    #endif
                    condassemble.wait(global_lock);
                }
            }
    #ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Unlocking mutex_global_assemble";
                LOGPZ_DEBUG(logger,sout.str())
            }
    #endif
            // global_lock is being destroyed...
        }
		parfront->AssembleElement(el, *ekaux, *efaux, *parfront->fStiffness, *parfront->fRhs);
		if(parfront->fCurrentAssembled == parfront->fNElements)
		{
#ifdef STACKSTORAGE
            TPZParFrontMatrix<TVar,TPZStackEqnStorage<TVar>, TFront> *mat = dynamic_cast< TPZParFrontMatrix<TVar, TPZStackEqnStorage<TVar>, front>* > (parfront->fStiffness);
#else
            TPZParFrontMatrix<TVar, TPZFileEqnStorage<TVar>, TFront> *mat = dynamic_cast<TPZParFrontMatrix<TVar, TPZFileEqnStorage<TVar>, TFront>* > (parfront->fStiffness);
            
#endif
			mat->FinishWriting();
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
			{
				std::stringstream sout;
				sout << "fFinishedComputing set to 1";
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
		}
		
		delete ekaux;
		delete efaux;
		
		if(parfront->fGuiInterface) if(parfront->fGuiInterface->AmIKilled()){
#ifdef STACKSTORAGE
			TPZParFrontMatrix<TVar, TPZStackEqnStorage<TVar>, front> *mat = dynamic_cast<TPZParFrontMatrix<TVar, TPZStackEqnStorage<TVar>, front>* > (parfront->fStiffness);
#else
            TPZParFrontMatrix<TVar, TPZFileEqnStorage<TVar>, TFront> *mat = dynamic_cast<TPZParFrontMatrix<TVar, TPZFileEqnStorage<TVar>, TFront>* > (parfront->fStiffness);
            
#endif
			mat->FinishWriting();
			break;
		}
		
		
	}//fim for iel
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Terminating assemble thread";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	std::cout << "Matrix assembled\n";
	std::cout.flush();

	return 0;
}

template<class TFront, class TVar, class TPar>
void TPZParFrontStructMatrix<TFront,TVar,TPar>::Assemble(TPZMatrix<TVar> & matref, TPZFMatrix<TVar> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
	this->fGuiInterface = guiInterface;
	
#ifdef STACKSTORAGE
	TPZParFrontMatrix<TVar, TPZStackEqnStorage<TVar>, TFront> *mat = dynamic_cast<TPZParFrontMatrix<TVar, TPZStackEqnStorage<TVar>, TFront> *>(&matref);
#else
    TPZParFrontMatrix<TVar, TPZFileEqnStorage<TVar>, TFront> *mat = dynamic_cast<TPZParFrontMatrix<TVar, TPZFileEqnStorage<TVar>, TFront>* > (&matref);
    
#endif
	if(!mat)
	{
		std::cout << __PRETTY_FUNCTION__ << " we are in serious trouble : wrong type of matrix"<< std::endl;
		DebugStop();
	}
	
	int nthreads;
	//cout << "Number of Threads " << endl;
	//cin >> nthreads;
	//fNThreads = nthreads;
    if (this->fNumThreads < 3) {
        this->fNumThreads = 3;
    }
	std::cout << "Number of Threads " << this->fNumThreads << std::endl;
	nthreads = this->fNumThreads;
	std::cout.flush();
	//int nthreads = fNThreads+1;
	
    std::vector<std::thread> allthreads(nthreads);
	int i;
	
	TPZVec <int> numelconnected(this->fMesh->NEquations(),0);
	//TPZFrontMatrix<TPZStackEqnStorage, TFront> *mat = new TPZFrontMatrix<TPZStackEqnStorage, TFront>(fMesh->NEquations());
	
	//TPZFrontMatrix<TPZFileEqnStorage, TFront> *mat = new TPZFrontMatrix<TPZFileEqnStorage, TFront>(fMesh->NEquations());
	
	
	
	
	//TPZParFrontMatrix<TPZStackEqnStorage, TFront> *mat = new TPZParFrontMatrix<TPZStackEqnStorage, TFront>(this->fMesh->NEquations());
	fNElements = this->fMesh->NElements();
	
	this->OrderElement();
	
//	this->AdjustSequenceNumbering();
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        this->fMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	
	this->GetNumElConnected(numelconnected);
	
	mat->SetNumElConnected(numelconnected);
	
    if (num_threads.was_set())
        mat->GetFront().ProductTensorMTInitData(num_threads.get_value());
    else
        mat->GetFront().ProductTensorMTInitData( nthreads - 2 ); // Here it initializes the multthread decomposition (comment to deactivate. Remember to coment the Finish also)
	
	fStiffness = mat;
	fRhs = &rhs;
	fCurrentElement = 0;
	fCurrentAssembled = 0;
    
	/*
	 *Triger 'n' threads passing Assemble as argument
	 */
	// try{
    allthreads[nthreads-1] = std::thread(this->GlobalAssemble,this);
//                                        this->GlobalAssemble, this, __FUNCTION__);
//	if(!res[nthreads-1]){
//#ifdef VC
//		cout << "GlobalAssemble Thread created Successfuly "<< allthreads[nthreads-1].x << endl;
//#else
//		std::cout << "GlobalAssemble Thread created Successfuly "<< allthreads[nthreads-1] << std::endl;
//#endif
//		std::cout.flush();
//	}else{
//#ifdef VC
//		std::cout << "GlobalAssemble Thread Fail "<< allthreads[nthreads-1].x << std::endl;
//#else
//		std::cout << "GlobalAssemble Thread Fail "<< allthreads[nthreads-1] << std::endl;
//#endif
//		std::cout.flush();
//		//          DebugStop();
//	}
    allthreads[nthreads-2] = std::thread(mat->WriteFile,mat);
//                                        mat->WriteFile, mat, __FUNCTION__);
//	if(!res[nthreads-2]){
//#ifdef VC
//		cout << "WriteFile Thread created Successfuly "<< allthreads[nthreads-2].x << endl;
//#else
//		std::cout << "WriteFile Thread created Successfuly "<< allthreads[nthreads-2] << std::endl;
//#endif
//		std::cout.flush();
//	}else{
//#ifdef VC
//		std::cout << "WriteFile Thread Fail "<< allthreads[nthreads-2].x << std::endl;
//#else
//		std::cout << "WriteFile Thread Fail "<< allthreads[nthreads-2] << std::endl;
//#endif
//		std::cout.flush();
//		//          DebugStop();
//	}
	
	for(i=0;i<nthreads-2;i++){
        allthreads[i] = std::thread(this->ElementAssemble,this);
//		if(!res[i]){
//#ifdef VC
//			cout << "ElementAssemble Thread "<< i+1 <<  " created Successfuly "<< allthreads[i].x << endl;
//#else
//			std::cout << "ElementAssemble Thread "<< i+1 <<  " created Successfuly "<< allthreads[i] << std::endl;
//#endif
//			std::cout.flush();
//		}else{
//			std::cout << "Error " << res[i] << "\t";
//#ifdef VC
//			std::cout << "ElementAssemble Thread "<< i+1 << " Fail " << allthreads[i].x << std::endl;
//#else
//			std::cout << "ElementAssemble Thread "<< i+1 << " Fail " << allthreads[i] << std::endl;
//#endif
//			std::cout.flush();
//		}
	}
	for(i=0;i<nthreads;i++) {
        allthreads[i].join();
	}
	
	mat->GetFront().ProductTensorMTFinish(); // Here it ends the multthread decomposition (comment to deactivate. Remember to coment the initialization also)
	fStiffness = 0;
	fRhs = 0;
    
}

template<class TFront, class TVar, class TPar>
TPZMatrix<TVar> * TPZParFrontStructMatrix<TFront,TVar,TPar>::CreateAssemble(TPZFMatrix<TVar> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPar::InitCreateAssemble();
	
	//TPZFrontMatrix<TPZStackEqnStorage, TFront> *mat = new TPZFrontMatrix<TPZStackEqnStorage, TFront>(fMesh->NEquations());
	
	//TPZFrontMatrix<TPZFileEqnStorage, TFront> *mat = new TPZFrontMatrix<TPZFileEqnStorage, TFront>(fMesh->NEquations());
	int64_t neq = this->fEquationFilter.NActiveEquations();
	
	//
#ifdef STACKSTORAGE
	TPZParFrontMatrix<TVar, TPZStackEqnStorage<TVar>, TFront> *mat = new TPZParFrontMatrix<TVar, TPZStackEqnStorage<TVar>, TFront>(neq);
#else
    TPZParFrontMatrix<TVar, TPZFileEqnStorage<TVar>, TFront> *mat = new TPZParFrontMatrix<TVar, TPZFileEqnStorage<TVar>, TFront>(neq);
#endif
    if (this->fDecomposeType != ENoDecompose)
    {
        mat->GetFront().SetDecomposeType(this->fDecomposeType);
    }
	rhs.Redim(neq,1);
	
	Assemble(*mat,rhs,guiInterface);
	return mat;
	
}

template<class TFront, class TVar, class TPar>
int TPZParFrontStructMatrix<TFront,TVar,TPar>::ClassId() const{
    return Hash("TPZParFrontStructMatrix") ^
        TPZFrontStructMatrix<TFront,TVar,TPar>::ClassId()<<1;
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZParFrontStructMatrix<TPZFrontSym<STATE>,STATE,TPZStructMatrixOR<STATE>>;
template class TPZParFrontStructMatrix<TPZFrontNonSym<STATE>,STATE,TPZStructMatrixOR<STATE>>;
template class TPZParFrontStructMatrix<TPZFrontSym<STATE>,STATE,TPZStructMatrixOT<STATE>>;
template class TPZParFrontStructMatrix<TPZFrontNonSym<STATE>,STATE,TPZStructMatrixOT<STATE>>;
template class TPZParFrontStructMatrix<TPZFrontSym<STATE>,STATE,TPZStructMatrixTBBFlow<STATE>>;
template class TPZParFrontStructMatrix<TPZFrontNonSym<STATE>,STATE,TPZStructMatrixTBBFlow<STATE>>;

template class TPZParFrontStructMatrix<TPZFrontSym<CSTATE>,CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZParFrontStructMatrix<TPZFrontNonSym<CSTATE>,CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZParFrontStructMatrix<TPZFrontSym<CSTATE>,CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZParFrontStructMatrix<TPZFrontNonSym<CSTATE>,CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZParFrontStructMatrix<TPZFrontSym<CSTATE>,CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;
template class TPZParFrontStructMatrix<TPZFrontNonSym<CSTATE>,CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;



