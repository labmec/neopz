/**
 * @file
 * @brief Contains the implementation of the TPZParFrontStructMatrix methods. 
 */

#include "TPZFrontMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"

#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"

#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzbndcond.h"

#include "pzanalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"

#include "TPZParFrontMatrix.h"

#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

#include "pzdxmesh.h"
#include <fstream>

#include "pzelmat.h"

#include "TPZFileEqnStorage.h"
#include "pzlog.h"


#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("pz.strmatrix.frontstructmatrix"));

#endif
/// Semaphore which controls threads assembling elements
pthread_mutex_t mutex_element_assemble = PTHREAD_MUTEX_INITIALIZER;
/// Semaphore which controls thread assembling global matrices
pthread_mutex_t mutex_global_assemble = PTHREAD_MUTEX_INITIALIZER;

/// Semaphore which controls condensed assembling
pthread_cond_t condassemble = PTHREAD_COND_INITIALIZER;
/// Semaphore
pthread_cond_t stackfull = PTHREAD_COND_INITIALIZER;


template<class front>
TPZParFrontStructMatrix<front>::TPZParFrontStructMatrix(TPZCompMesh *mesh): TPZFrontStructMatrix<front>(mesh)
{
	fMaxStackSize = 500;
	TPZStructMatrix::SetNumThreads(3);
}

template<class front>
TPZParFrontStructMatrix<front>::TPZParFrontStructMatrix(const TPZParFrontStructMatrix &copy): TPZFrontStructMatrix<front>(copy), fMaxStackSize(copy.fMaxStackSize)
{
}

template<class front>
TPZParFrontStructMatrix<front>::~TPZParFrontStructMatrix()
{

}

template<class front>
TPZStructMatrix * TPZParFrontStructMatrix<front>::Clone(){
	TPZParFrontStructMatrix<front> * mat = new TPZParFrontStructMatrix<front>(*this);
	return mat;
	;
	
}

template<class front>
void *TPZParFrontStructMatrix<front>::ElementAssemble(void *t){
    
    
	TPZParFrontStructMatrix<front> *parfront = (TPZParFrontStructMatrix<front> *) t;
	
	TPZAdmChunkVector<TPZCompEl *> &elementvec = parfront->fMesh->ElementVec();
	
	
	while(parfront->fCurrentElement < parfront->fNElements) {
		
		/**
		 *Lock mutex and search for an avilable element
		 *A global variable to be updated whenever a element is processed
		 */
		
		//Lock a mutex and get an element number
		
		PZ_PTHREAD_MUTEX_LOCK(&mutex_element_assemble, "TPZParFrontStructMatrix<front>::ElementAssemble()");
		
		//Stack is full and process must wait here!
		if(parfront->felnum.NElements()==parfront->fMaxStackSize)
        {
			/*          cout << "    Stack full" << endl;
			 cout << "    Waiting" << endl;
			 cout.flush();*/
			//cout << "Mutex unlocked on Condwait" << endl;
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
			{
				std::stringstream sout;
				sout << "Entering cond_wait because of stack overflow ";
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			PZ_PTHREAD_COND_WAIT(&stackfull,&mutex_element_assemble,"TPZParFrontStructMatrix<front>::ElementAssemble()");
			//cout << "Mutex LOCKED leaving Condwait" << endl;
			
        }
		
		//cout << "Locking mutex_element_assemble" << endl;
		//cout.flush();
		int64_t local_element = parfront->fCurrentElement;
		if(local_element==parfront->fNElements)
        {
            PZ_PTHREAD_MUTEX_UNLOCK(&mutex_element_assemble, "TPZParFrontStructMatrix<front>::ElementAssemble()");
            return 0;
        }
		/*          cout << "All element matrices assembled" << endl;
		 return 0;
		 }
		 */
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
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
		PZ_PTHREAD_MUTEX_UNLOCK(&mutex_element_assemble, "TPZParFrontStructMatrix<front>::ElementAssemble()");
        
		
		if(parfront->fElementOrder[local_element] < 0) continue;
		TPZCompEl *el = elementvec[parfront->fElementOrder[local_element]];
		if(!el) continue;
		//		int dim = el->NumNodes();
		
		//Builds elements stiffness matrix
		TPZElementMatrix *ek = new TPZElementMatrix(parfront->fMesh,TPZElementMatrix::EK);
		TPZElementMatrix *ef = new TPZElementMatrix(parfront->fMesh,TPZElementMatrix::EF);
		
		el->CalcStiff(*ek, *ef);
		//Locks a mutex and adds element contribution to frontmatrix
		//if mutex is locked go to condwait waiting for an specific condvariable
		// este mutex deve ser outro mutex -> mutexassemble
		
		PZ_PTHREAD_MUTEX_LOCK(&mutex_global_assemble, "TPZParFrontStructMatrix<front>::ElementGlobalAssemble()");
		//cout << "Locking mutex_global_assemble" << endl;
		//cout << "Pushing variables to the stack" << endl;
		//cout.flush();
		
		// colocar ek, ef e o element_local no stack
		parfront->felnum.Push(local_element);
		parfront->fekstack.Push(ek);
		parfront->fefstack.Push(ef);
		
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
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
		 cout << endl << "Computing " << parfront->fCurrentElement << " on thread " << pthread_self() << endl;
		 cout << " " << (100*parfront->fCurrentElement/parfront->fNElements) << "% Elements computed" << "     " << (100*parfront->fCurrentAssembled/parfront->fNElements) << "% Elements assembled" << endl;
		 }
		 cout << '*';
		 cout.flush();
		 */
		//Alterado cond_broadcast para cond_signal
		//invertendo a sequencia das chamadas
		PZ_PTHREAD_COND_BROADCAST(&condassemble,"TPZParFrontStructMatrix<front>::CondAssemble()");
		PZ_PTHREAD_MUTEX_UNLOCK(&mutex_global_assemble,"TPZParFrontStructMatrix<front>::ElementGlobalAssemble()");
        
		// o thread de assemblagem utiliza mutexassemble
		// e feito em outro thread     AssembleElement(el, ek, ef, stiffness, rhs);
		
		if(parfront->fGuiInterface) if(parfront->fGuiInterface->AmIKilled()){
			break;
		}
	}//fim for iel
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
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

template<class front>
void *TPZParFrontStructMatrix<front>::GlobalAssemble(void *t){
    
	TPZParFrontStructMatrix<front> *parfront = (TPZParFrontStructMatrix<front> *) t;
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
		PZ_PTHREAD_MUTEX_LOCK(&mutex_global_assemble,"TPZParFrontStructMatrix<front>::GlobalAssemble()");
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
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
						PZ_PTHREAD_COND_BROADCAST(&stackfull,"TPZParFrontStructMatrix<front>::GlobalAssemble()");
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
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
				{
					std::stringstream sout;
					sout << "Waiting on condassemble";
					LOGPZ_DEBUG(logger,sout.str())
				}
#endif
				PZ_PTHREAD_COND_WAIT(&condassemble, &mutex_global_assemble, "TPZParFrontStructMatrix<front>::GlobalAssemble()");
			}
		}
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
		{
			std::stringstream sout;
			sout << "Unlocking mutex_global_assemble";
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		PZ_PTHREAD_MUTEX_UNLOCK(&mutex_global_assemble,"TPZParFrontStructMatrix<front>::GlobalAssemble()");
		parfront->AssembleElement(el, *ekaux, *efaux, *parfront->fStiffness, *parfront->fRhs);
		if(parfront->fCurrentAssembled == parfront->fNElements)
		{
#ifdef STACKSTORAGE
            TPZParFrontMatrix<STATE,TPZStackEqnStorage<STATE>, front> *mat = dynamic_cast< TPZParFrontMatrix<STATE, TPZStackEqnStorage<STATE>, front>* > (parfront->fStiffness);
#else
            TPZParFrontMatrix<STATE, TPZFileEqnStorage<STATE>, front> *mat = dynamic_cast<TPZParFrontMatrix<STATE, TPZFileEqnStorage<STATE>, front>* > (parfront->fStiffness);
            
#endif
			mat->FinishWriting();
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
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
			TPZParFrontMatrix<STATE, TPZStackEqnStorage<STATE>, front> *mat = dynamic_cast<TPZParFrontMatrix<STATE, TPZStackEqnStorage<STATE>, front>* > (parfront->fStiffness);
#else
            TPZParFrontMatrix<STATE, TPZFileEqnStorage<STATE>, front> *mat = dynamic_cast<TPZParFrontMatrix<STATE, TPZFileEqnStorage<STATE>, front>* > (parfront->fStiffness);
            
#endif
			mat->FinishWriting();
			break;
		}
		
		
	}//fim for iel
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
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

template<class front>
void TPZParFrontStructMatrix<front>::Assemble(TPZMatrix<STATE> & matref, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
	this->fGuiInterface = guiInterface;
	
#ifdef STACKSTORAGE
	TPZParFrontMatrix<STATE, TPZStackEqnStorage<STATE>, front> *mat = dynamic_cast<TPZParFrontMatrix<STATE, TPZStackEqnStorage<STATE>, front> *>(&matref);
#else
    TPZParFrontMatrix<STATE, TPZFileEqnStorage<STATE>, front> *mat = dynamic_cast<TPZParFrontMatrix<STATE, TPZFileEqnStorage<STATE>, front>* > (&matref);
    
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
	
	pthread_t *allthreads = new pthread_t[nthreads];
	int *res = new int[nthreads];
	int i;
	
	TPZVec <int> numelconnected(this->fMesh->NEquations(),0);
	//TPZFrontMatrix<TPZStackEqnStorage, front> *mat = new TPZFrontMatrix<TPZStackEqnStorage, front>(fMesh->NEquations());
	
	//TPZFrontMatrix<TPZFileEqnStorage, front> *mat = new TPZFrontMatrix<TPZFileEqnStorage, front>(fMesh->NEquations());
	
	
	
	
	//TPZParFrontMatrix<TPZStackEqnStorage, front> *mat = new TPZParFrontMatrix<TPZStackEqnStorage, front>(this->fMesh->NEquations());
	fNElements = this->fMesh->NElements();
	
	this->OrderElement();
	
//	this->AdjustSequenceNumbering();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
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
	//pthread_create(&allthreads[fNThreads-1],NULL,this->GlobalAssemble, this);
	// try{
	res[nthreads-1] = PZ_PTHREAD_CREATE(&allthreads[nthreads-1], NULL,
                                        this->GlobalAssemble, this, __FUNCTION__);
	if(!res[nthreads-1]){
#ifdef VC
		cout << "GlobalAssemble Thread created Successfuly "<< allthreads[nthreads-1].x << endl;
#else
		std::cout << "GlobalAssemble Thread created Successfuly "<< allthreads[nthreads-1] << std::endl;
#endif
		std::cout.flush();
	}else{
#ifdef VC
		std::cout << "GlobalAssemble Thread Fail "<< allthreads[nthreads-1].x << std::endl;
#else
		std::cout << "GlobalAssemble Thread Fail "<< allthreads[nthreads-1] << std::endl;
#endif
		std::cout.flush();
		//          DebugStop();
	}
	res[nthreads-2] = PZ_PTHREAD_CREATE(&allthreads[nthreads-2], NULL, 
                                        mat->WriteFile, mat, __FUNCTION__);
	if(!res[nthreads-2]){
#ifdef VC
		cout << "WriteFile Thread created Successfuly "<< allthreads[nthreads-2].x << endl;
#else
		std::cout << "WriteFile Thread created Successfuly "<< allthreads[nthreads-2] << std::endl;
#endif
		std::cout.flush();
	}else{
#ifdef VC
		std::cout << "WriteFile Thread Fail "<< allthreads[nthreads-2].x << std::endl;
#else
		std::cout << "WriteFile Thread Fail "<< allthreads[nthreads-2] << std::endl;
#endif
		std::cout.flush();
		//          DebugStop();
	}
	
	for(i=0;i<nthreads-2;i++){
        res[i] = PZ_PTHREAD_CREATE(&allthreads[i], NULL, 
                                   this->ElementAssemble, this, __FUNCTION__);
		if(!res[i]){
#ifdef VC
			cout << "ElementAssemble Thread "<< i+1 <<  " created Successfuly "<< allthreads[i].x << endl;
#else
			std::cout << "ElementAssemble Thread "<< i+1 <<  " created Successfuly "<< allthreads[i] << std::endl;
#endif
			std::cout.flush();
		}else{
			std::cout << "Error " << res[i] << "\t";
#ifdef VC
			std::cout << "ElementAssemble Thread "<< i+1 << " Fail " << allthreads[i].x << std::endl;
#else
			std::cout << "ElementAssemble Thread "<< i+1 << " Fail " << allthreads[i] << std::endl;
#endif
			std::cout.flush();
		}
	}
	for(i=0;i<nthreads;i++) {
        PZ_PTHREAD_JOIN(allthreads[i], NULL, __FUNCTION__);
	}
	
	delete[] allthreads;// fThreadUsed, fDec;
	delete[] res;
	mat->GetFront().ProductTensorMTFinish(); // Here it ends the multthread decomposition (comment to deactivate. Remember to coment the initialization also)
	fStiffness = 0;
	fRhs = 0;
    
}


#ifndef STATE_COMPLEX
#include "pzmat2dlin.h"

template<class front>
int TPZParFrontStructMatrix<front>::main() {
	
	int refine=1;
	int order=1;
	
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
		
		// initializar os dados do n�		gmesh.NodeVec ()[i].Initialize (i,coord,gmesh);
	}
	int el;
	//TPZGeoEl *gel;
	for(el=0; el<1; el++) {
		
		// initializar os indices dos n�
		TPZVec<int64_t> indices(4);
		for(i=0; i<4; i++) indices[i] = i;
		// O proprio construtor vai inserir o elemento na malha
		int64_t index;
		/*gel = */gmesh.CreateGeoElement(EQuadrilateral,indices,1,index);
	}
	gmesh.BuildConnectivity ();
	
	TPZVec<TPZGeoEl *> subel;
	//gel->Divide(subel);
	
	
	
	std::cout << "Refinement ";
	std::cin >> refine;
	std::cout << std::endl;
	DebugStop();
	//UniformRefine(refine,gmesh);
	
	
	TPZMat2dLin *mat2d = new TPZMat2dLin(1);
	TPZFMatrix<STATE> xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
	mat2d->SetMaterial (xk,xc,xf);
	auto * meumat = (TPZMaterial*) mat2d;
	cmesh.InsertMaterialObject(meumat);
	
	TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
	auto * bnd = (TPZMaterial*) meumat->CreateBC (meumat,-4,0,val1,val2);
	cmesh.InsertMaterialObject(bnd);
	
	
	
	std::cout << "Interpolation order ";
	std::cin >> order;
	std::cout << std::endl;
	
	//	TPZCompEl::gOrder = order;
	cmesh.SetDefaultOrder(order);
	
	cmesh.AutoBuild();
	//	cmesh.AdjustBoundaryElements();
	cmesh.InitializeBlock();
	
	std::ofstream output("outputPar.dat");
	//	ofstream output2("outputNon.dat");
	cmesh.Print(output);
	TPZAnalysis an(&cmesh,true,output);
	//	TPZAnalysis an2(&cmesh,output);
	
	TPZVec<int> numelconnected(cmesh.NEquations(),0);
	int64_t ic;
	//cout << "Nmero de Equa�es -> " << cmesh.NEquations() << endl;
	//cout.flush();
	
	std::ofstream out("cmeshBlock_out.txt");
	//	cmesh.Print(out);
	//	cmesh.Block().Print("Block",out);
	for(ic=0; ic<cmesh.ConnectVec().NElements(); ic++) {
		TPZConnect &cn = cmesh.ConnectVec()[ic];
		if(cn.HasDependency() || cn.IsCondensed()) continue;
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
	
	TPZParFrontStructMatrix<TPZFrontSym<STATE> > mat(&cmesh);
	
	//   TPZFStructMatrix mat2(&cmesh);
	//  mat->SetNumElConnected(numelconnected);
	//mat = CreateAssemble();
	int threads;
	std::cout << "Number of Threads  ";
	std::cin >> threads;
	std::cout << std::endl;
	
	mat.SetNumThreads(threads);
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
	an.Print("solution of frontal solver", output);
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
	std::cout.flush();
	return 0;
	
}
#endif

template<class front>
TPZMatrix<STATE> * TPZParFrontStructMatrix<front>::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
	
	//TPZFrontMatrix<TPZStackEqnStorage, front> *mat = new TPZFrontMatrix<TPZStackEqnStorage, front>(fMesh->NEquations());
	
	//TPZFrontMatrix<TPZFileEqnStorage, front> *mat = new TPZFrontMatrix<TPZFileEqnStorage, front>(fMesh->NEquations());
	int64_t neq = this->fEquationFilter.NActiveEquations();
	
	//
#ifdef STACKSTORAGE
	TPZParFrontMatrix<STATE, TPZStackEqnStorage<STATE>, front> *mat = new TPZParFrontMatrix<STATE, TPZStackEqnStorage<STATE>, front>(neq);
#else
    TPZParFrontMatrix<STATE, TPZFileEqnStorage<STATE>, front> *mat = new TPZParFrontMatrix<STATE, TPZFileEqnStorage<STATE>, front>(neq);
#endif
    if (this->fDecomposeType != ENoDecompose)
    {
        mat->GetFront().SetDecomposeType(this->fDecomposeType);
    }
	rhs.Redim(neq,1);
	
	Assemble(*mat,rhs,guiInterface);
	return mat;
	
}


template<class TVar>
class TPZFrontSym;
template<class TVar>
class TPZFrontNonSym;

template class TPZParFrontStructMatrix<TPZFrontSym<STATE> >;
template class TPZParFrontStructMatrix<TPZFrontNonSym<STATE> >;

