/**
 * @file
 * @brief Contains the implementation of the TPZStackEqnStorage methods.
 */

#include "TPZStackEqnStorage.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzreal.h"
#include <stdio.h>

using namespace std;

template<class TVar>
void TPZStackEqnStorage<TVar>::SetBlockSize(){}

template<class TVar>
void TPZStackEqnStorage<TVar>::ReOpen(){}

template<class TVar>
void TPZStackEqnStorage<TVar>::Print(const char *name, std::ostream& out) const {
    int64_t i, loop_limit;
	loop_limit=fEqnStack.NElements();
    out <<  "Number of entries on EqnStack  "<< fEqnStack.NElements() << endl;
    for(i=0;i<loop_limit;i++) fEqnStack[i].Print(name, out);
}
template<class TVar>
void TPZStackEqnStorage<TVar>::Reset()
{
	fEqnStack.Resize(0);
}
template<class TVar>
void TPZStackEqnStorage<TVar>::Backward(TPZFMatrix<TVar> &f, DecomposeType dec) const
{
	int64_t i, stack_size;
	stack_size=fEqnStack.NElements();
	for(i=stack_size-1;i>=0;i--){
		fEqnStack[i].EqnBackward(f, dec);
	}
	
}
template<class TVar>
void TPZStackEqnStorage<TVar>::Forward(TPZFMatrix<TVar> &f, DecomposeType dec) const
{
	int64_t i, stack_size;
	stack_size=fEqnStack.NElements();
	for(i=0;i<stack_size;i++){
		fEqnStack[i].EqnForward(f, dec);
	}
	
}

template<class TVar>
void TPZStackEqnStorage<TVar>::AddEqnArray(TPZEqnArray<TVar> *EqnArray)
{
	fEqnStack.Push(*EqnArray);
}

template<class TVar>
TPZStackEqnStorage<TVar>::TPZStackEqnStorage() : TPZRegisterClassId(&TPZStackEqnStorage<TVar>::ClassId)
{
}

template<class TVar>
void TPZStackEqnStorage<TVar>::Zero()
{
	fEqnStack.Resize(0);
}

template<class TVar>
TPZStackEqnStorage<TVar>::~TPZStackEqnStorage()
{
}

template<class TVar>
void TPZStackEqnStorage<TVar>::main()
{
}

template<class TVar>
TPZStackEqnStorage<TVar>::TPZStackEqnStorage(char option, const char *name) : TPZRegisterClassId(&TPZStackEqnStorage<TVar>::ClassId)
{
	
}

template<class TVar>
void TPZStackEqnStorage<TVar>::OpenGeneric(char option, const char * name){}
template<class TVar>
void TPZStackEqnStorage<TVar>::ReadBlockPositions() {}
template<class TVar>
void TPZStackEqnStorage<TVar>::FinishWriting(){}
template<class TVar>
std::string TPZStackEqnStorage<TVar>::GetStorage() {return "Stack Storage";}

template class TPZStackEqnStorage<float>;
template class TPZStackEqnStorage<double>;
template class TPZStackEqnStorage<long double>;

template class TPZStackEqnStorage<std::complex<float> >;
template class TPZStackEqnStorage<std::complex<double> >;
template class TPZStackEqnStorage<std::complex<long double> >;
