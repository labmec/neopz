/**
 * @file
 * @brief Contains the implementation of the TPZEqnArray methods.
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include "tpzeqnarray.h"
#include "pzreal.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void TPZEqnArray::SetNonSymmetric(){
    fSymmetric=EIsNonSymmetric;
}
int TPZEqnArray::IsSymmetric(){
    return fSymmetric;
}
void TPZEqnArray::SetSymmetric(){
    fSymmetric = EIsSymmetric;
}
TPZEqnArray::~TPZEqnArray(){
	
}
TPZEqnArray::TPZEqnArray() : fEqStart(), fEqNumber(), fEqValues(), fIndex() {
	fEqStart.Push(0);
	fNumEq=0;
	fLastTerm=0;
	fSymmetric=EIsUndefined;
	
}
void TPZEqnArray::Print(const char *name, std::ostream& out)
{
	if(name) out << name << endl;
	int i,j;
	for(i=0;i<fNumEq;i++){
		int aux_limit;
		
		if(i==fNumEq-1){
			aux_limit=fEqValues.NElements();
		}else aux_limit=fEqStart[i+1];
		for(j=fEqStart[i]; j<aux_limit ; j++) {
			out << "col = " << fIndex[j] << "   ";
		}
		out << endl;
		for(j=fEqStart[i]; j<aux_limit ; j++) {
			out << fEqValues[j] << "  ";
		}
		out << endl;
		out << endl;
	}
}


void TPZEqnArray::Write(char * outputfile){}
void TPZEqnArray::Read(char * inputfile){}
void TPZEqnArray::EndEquation(){
	fEqStart.Push(fLastTerm);
}
/*void TPZEqnArray::AddTerm(int col, REAL val){
 fIndex.Push(col);
 fEqValues.Push(val);
 fLastTerm++;
 } 
 */
void TPZEqnArray::BeginEquation(int eq){
	fEqNumber.Push(eq);
	fNumEq++;
}
void TPZEqnArray::Reset(){
	fEqStart.Resize(0);
	fEqStart.Push(0);
	fEqValues.Resize(0);
	fIndex.Resize(0);
	fNumEq=0;
	fLastTerm=0;
	fSymmetric=EIsUndefined;
	
}
void TPZEqnArray::EqnBackward(TPZFMatrix & U, DecomposeType dec){
	int n;
	if(IsSymmetric()==EIsSymmetric){
    	for(n=fNumEq-1; n>=0; n--) {
    		int index = fEqStart[n]; 
    		
			//		if(fEqStart.NElements()==n+1){
			//			fEqStart.Expand(n+2);
			//		}
    		int last_term = fEqStart[n + 1];
    		REAL acum = 0.;
    		int i;
    		for(i = index + 1; i < last_term; i++) acum -= U(fIndex[i],0) * fEqValues[i];
			
            /** Case LU - LDLT - Cholesky*/
            U(fIndex[index],0) += acum;
			
            /** Case Cholesky */
            if(dec == ECholesky){
               	U(fIndex[index],0) /= fEqValues[index];
    		}
    	}
	}else if(IsSymmetric()==EIsNonSymmetric){
        for(n=fNumEq-2; n>=0; n-=2) {
            int index = fEqStart[n]; 
            
			//		if(fEqStart.NElements()==n+1){
			//			fEqStart.Expand(n+2);
			//		}
            int last_term = fEqStart[n + 1];
            REAL acum = 0.;
            int i;
            for(i = index + 1; i < last_term; i++) acum -= U(fIndex[i],0) * fEqValues[i];
			
            /** Case LU - LDLT - Cholesky*/
            U(fIndex[index],0) += acum;
			
            /** Case Cholesky */
            if(dec == ECholesky || dec == ELU){
                U(fIndex[index],0) /= fEqValues[index];
            }
        }
    }
}

void TPZEqnArray::EqnForward(TPZFMatrix & F, DecomposeType dec){
	int j;
	if(IsSymmetric()==EIsSymmetric){
		
		for(j=0; j<fNumEq; j++) {
			int index = fEqStart[j];
			if(fEqValues[index] == 0.){
				cout << "Diagonal Value = 0 >> Aborting on Index = " << index << " Equation = " << j << endl;
				DebugStop();
			}          
			int last_term;
			if(j==fNumEq-1){
				last_term=fEqValues.NElements();
			}else last_term = fEqStart[j+1];
			//if(fEqStart.NElements()
			
			
			/** cholesky e lu */
			if(dec==ECholesky || dec==ELU) {
				F(fIndex[index],0) /= fEqValues[index];
			}
			/** ldlt cholesky lu */
			REAL udiag = F(fIndex[index],0);
			
			
			int i;
			//+1 ou +2
			for(i = index + 1; i < last_term; i++) F(fIndex[i],0) -= udiag * fEqValues[i];
			/** finalizacao para ldlt */
			if(dec == ELDLt) F(fIndex[index],0) /= fEqValues[index];
		}
	}else if(IsSymmetric()==EIsNonSymmetric){
		for(j=1; j<fNumEq; j+=2) {
			int index = fEqStart[j];
			
			if(fEqValues[index] == 0.){
				cout << "Diagonal Value = 0 >> Aborting on Index = " << index << " Equation = " << j << endl;
				DebugStop();
			}          
			
			int last_term;
			if(j==fNumEq-1){
				last_term=fEqValues.NElements();
			}else last_term = fEqStart[j+1];
			//if(fEqStart.NElements()
			
			
			/** cholesky e lu */
			if(dec==ECholesky || dec==ELU) {
				F(fIndex[index],0) /= fEqValues[index];
			}
			/** ldlt cholesky lu */
			REAL udiag = F(fIndex[index],0);
			
			
			int i;
			for(i = index + 1; i < last_term; i++) F(fIndex[i],0) -= udiag * fEqValues[i];
			/** finalizacao para ldlt */
			if(dec == ELDLt) F(fIndex[index],0) /= fEqValues[index];
		}
	}
}

void TPZEqnArray::Write(FILE * outputfile){
	/** Number of equations */
	fwrite(&fNumEq,sizeof(int),1,outputfile);
	//cout << ftell(outputfile) << endl;
	/** Last term added*/
	fwrite(&fLastTerm,sizeof(int),1,outputfile);
	//cout << ftell(outputfile) << endl;
	/** TPZStack fEqStart data */
	int aux = fEqStart.NElements();
	if(aux == 0) 
	{
		std::cout << __PRETTY_FUNCTION__ << __LINE__ << std::endl;
	}
	fwrite(&aux,sizeof(int),1,outputfile);
	fwrite(&fEqStart[0],sizeof(int), aux ,outputfile);
	
	/** TPZStack fEqNumber data */
	aux = fEqNumber.NElements();
	if(aux == 0) 
	{
		std::cout << __PRETTY_FUNCTION__ << __LINE__ << std::endl;
	}
	fwrite(&aux,sizeof(int),1,outputfile);
	fwrite(&fEqNumber[0],sizeof(int), aux ,outputfile);
	
	/** TPZStack fIndex data */
	aux = fIndex.NElements();
	fwrite(&aux,sizeof(int),1,outputfile);
	fwrite(&fIndex[0],sizeof(int), aux ,outputfile);
	
	/** TPZStack fEqValues data */
	aux = fEqValues.NElements();
	fwrite(&aux,sizeof(int),1,outputfile);
	fwrite(&fEqValues[0],sizeof(REAL), aux ,outputfile);
}
void TPZEqnArray::Read(FILE * inputfile){
	/** Number of equations */
	fread(&fNumEq,sizeof(int),1,inputfile);
	//cout << ftell(inputfile) << endl; 
	/** Last term added*/
	fread(&fLastTerm,sizeof(int),1,inputfile);
	//cout << ftell(inputfile) << endl;
	int aux=0;
	/** TPZStack fEqStart data */
	fread(&aux,sizeof(int),1,inputfile);
	fEqStart.Resize(aux);
	fread(&fEqStart[0],sizeof(int), aux ,inputfile);
	
	/** TPZStack fEqNumber data */
	fread(&aux,sizeof(int),1,inputfile);
	fEqNumber.Resize(aux);
	fread(&fEqNumber[0],sizeof(int), aux ,inputfile);
	
	/** TPZStack fIndex data */
	fread(&aux,sizeof(int),1,inputfile);
	fIndex.Resize(aux);
	fread(&fIndex[0],sizeof(int), aux ,inputfile);
	
	/** TPZStack fEqValues data */
	fread(&aux,sizeof(int),1,inputfile);
	fEqValues.Resize(aux);
	fread(&fEqValues[0],sizeof(REAL), aux ,inputfile);
	
	this->fSymmetric = EIsSymmetric;
	if(fNumEq && fNumEq%2==0 && fEqNumber[0]==fEqNumber[1]) fSymmetric = EIsNonSymmetric;
}


void TPZEqnArray::main()
{
	char  filename[20];
	
	cout << "Entre o nome do Arquivo\n";
	cin >> filename;
	ofstream output(filename,ios::app);
	
	TPZFMatrix MatrixA(10,10);
	int i, j; 
	for(i=0;i<10;i++) {
		for(j=i;j<10;j++) {
			int random = rand();
			double rnd = (random*10.)/RAND_MAX;
			MatrixA(i,j)=rnd;
			MatrixA(j,i)=MatrixA(i,j);
			if(i==j) MatrixA(i,j)=6000.;
		}
	}
	
	
	MatrixA.Print("Teste 1");
	MatrixA.Decompose_Cholesky();
	MatrixA.Print("Decomposta");
	
	TPZEqnArray Test;
	Test.Reset();
	for(i=0;i<10;i++) {
		Test.BeginEquation(i);
		for(j=i;j<10;j++) {
			Test.AddTerm(j, MatrixA(i,j));
		}
		Test.EndEquation(); 
	}
	
	TPZFMatrix rhs(10,1), rhs2(10,1);
	//Inicializar rhs:
	rhs.Zero();
	rhs(0,0) = 1.;
	rhs2=rhs;
	
	Test.Print("Result", cout);
	
	MatrixA.Subst_Forward(&rhs);
	DecomposeType decType = ECholesky;
	Test.EqnForward(rhs2, decType);
	
	rhs.Print("FMatrix Decomposition"); 
	rhs2.Print("FrontalMatrix Decomposition"); 
	
}
