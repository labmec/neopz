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

void TPZStackEqnStorage::SetBlockSize(){}

void TPZStackEqnStorage::ReOpen(){}

void TPZStackEqnStorage::Print(const char *name, std::ostream& out) const {
    int i, loop_limit;
	loop_limit=fEqnStack.NElements();
    out <<  "Number of entries on EqnStack  "<< fEqnStack.NElements() << endl;
    for(i=0;i<loop_limit;i++) fEqnStack[i].Print(name, out);
}
void TPZStackEqnStorage::Reset()
{
	fEqnStack.Resize(0);
}
void TPZStackEqnStorage::Backward(TPZFMatrix<REAL> &f, DecomposeType dec) const
{
	int i, stack_size;
	stack_size=fEqnStack.NElements();
	for(i=stack_size-1;i>=0;i--){
		fEqnStack[i].EqnBackward(f, dec);
	}
	
}
void TPZStackEqnStorage::Forward(TPZFMatrix<REAL> &f, DecomposeType dec) const
{
	int i, stack_size;
	stack_size=fEqnStack.NElements();
	for(i=0;i<stack_size;i++){
		fEqnStack[i].EqnForward(f, dec);
	}
	
}

void TPZStackEqnStorage::AddEqnArray(TPZEqnArray *EqnArray)
{
	fEqnStack.Push(*EqnArray);
}

TPZStackEqnStorage::TPZStackEqnStorage()
{
}

void TPZStackEqnStorage::Zero()
{
	fEqnStack.Resize(0);
}

TPZStackEqnStorage::~TPZStackEqnStorage()
{
}

void TPZStackEqnStorage::main()
{
}

TPZStackEqnStorage::TPZStackEqnStorage(char option, const char *name)
{
	
}
//void TPZStackEqnStorage::SetFileName(const char *name){}
void TPZStackEqnStorage::OpenGeneric(char option, const char * name){}
void TPZStackEqnStorage::ReadBlockPositions() {}
void TPZStackEqnStorage::FinishWriting(){}
std::string TPZStackEqnStorage::GetStorage() {return "Stack Storage";}

