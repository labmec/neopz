/**
 * @file
 * @brief Contains the implementation of the TPZMTAssemble methods. 
 */
//$Id: TPZMTAssemble.cpp,v 1.4 2011-05-11 03:11:17 phil Exp $

#include "TPZMTAssemble.h"
#include "pzstrmatrix.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzintel.h"

#include "pzgnode.h"
#include "TPZTimer.h"

#include "pzmaterial.h"

using namespace std;

std::vector< std::pair< TPZElementMatrix *, SMTAssembleResidual*> > TPZMTAssemble::gComputedEF;

TPZMTAssemble::TPZMTAssemble()
{
}

TPZMTAssemble::~TPZMTAssemble()
{
}

void * TPZMTAssemble::ExecuteAssembleResidualMT(void * ExtData){
	
	SMTAssembleResidual * data = static_cast<SMTAssembleResidual*>(ExtData);
	
	TPZCompEl *el = data->compel;
	if(!el) return NULL;
	
	if(data->MaterialIds){
		TPZAutoPointer<TPZMaterial> mat = el->Material();
		if (!mat) return NULL;
		int matid = mat->Id();
		if (data->MaterialIds->find(matid) == data->MaterialIds->end()) return NULL;
	}//if
	
	TPZCompMesh * mesh = el->Mesh();
	TPZElementMatrix * ef = new TPZElementMatrix(mesh, TPZElementMatrix::EF);
	
	el->CalcResidual(*ef);
	
	pair<TPZElementMatrix*, SMTAssembleResidual*> mypair(ef,data);
	TPZMTAssemble::gComputedEF[data->index] = mypair;
	
	return NULL;
}

void TPZMTAssemble::AssembleMT(TPZFMatrix<REAL> & rhs, TPZCompMesh &mesh, int mineq, int maxeq, std::set<int> *MaterialIds){
	int iel;
	const int nelem = mesh.NElements();
	
	const int nthreads = 3;
	pthread_t allthreads[nthreads];
	
	TPZMTAssemble::gComputedEF.resize(nthreads);
	
	for(iel=0; iel < nelem;) {
		
		for(int ithread = 0; ithread < nthreads; ithread++){
			allthreads[ithread] = NULL;
			std::pair< TPZElementMatrix *, SMTAssembleResidual * > nullPair(NULL,NULL);
			TPZMTAssemble::gComputedEF[ithread] = nullPair;
		}
		
		for(int ithread = 0; ithread < nthreads && iel < nelem; ithread++, iel++){
			TPZCompEl *el = mesh.ElementVec()[iel];
			
			if(!el){
				ithread--;
				continue;
			}
			
			if(MaterialIds){
				TPZAutoPointer<TPZMaterial> mat = el->Material();
				if (!mat){
					ithread--;
					continue;
				}
				int matid = mat->Id();
				if (MaterialIds->find(matid) == MaterialIds->end()){
					ithread--;
					continue;
				} 
			}//if
			
			SMTAssembleResidual * data = new SMTAssembleResidual(el, &rhs, mineq, maxeq, MaterialIds, ithread);
			pthread_create(&allthreads[ithread],NULL,ExecuteAssembleResidualMT, data);
			
		}//threads
		
		for(int i=0;i<nthreads;i++){
			if(!allthreads[i]) continue;
			pthread_join(allthreads[i], NULL);
		}
		
		
		TPZMTAssemble::ContributeEFs();
		
	}//fim for iel
	
}
// ofstream effilemt("MTEF.txt");
void TPZMTAssemble::ContributeEFs(){
	
	for(int i = 0; i < (int)TPZMTAssemble::gComputedEF.size(); i++){
		
		TPZElementMatrix * ef = TPZMTAssemble::gComputedEF[i].first;
		SMTAssembleResidual * data = TPZMTAssemble::gComputedEF[i].second;
		
		if(!ef || !data) continue;
		
		TPZCompEl * el = data->compel;
		const int mineq = data->mineq;
		const int maxeq = data->maxeq;
		
		if(!el->HasDependency()) {
			ef->ComputeDestinationIndices();
			if(mineq != -1 & maxeq != -1)
			{
				TPZStructMatrix::FilterEquations(ef->fSourceIndex,ef->fDestinationIndex,mineq,maxeq);
			}
			data->rhs->AddFel(ef->fMat, ef->fSourceIndex, ef->fDestinationIndex);
		} else {
			// the element has dependent nodes
			ef->ApplyConstraints();
			ef->ComputeDestinationIndices();
			if(mineq != -1 & maxeq != -1)
			{
				TPZStructMatrix::FilterEquations(ef->fSourceIndex,ef->fDestinationIndex,mineq,maxeq);
			}
			data->rhs->AddFel(ef->fConstrMat,ef->fSourceIndex,ef->fDestinationIndex);
		}
		
		delete ef;
		delete data;
		
	}//for
	
}//void

