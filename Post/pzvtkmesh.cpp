/**
 * @file
 * @brief Contains the implementation of the TPZVTKGraphMesh methods. 
 */

#include <sstream>
#include "pzvtkmesh.h"
#include "pzcmesh.h"
#include "TPZMaterial.h"
#include "pzgraphnode.h"
#include "pzgraphel.h"
#include "TPZSimpleTimer.h"

using namespace std;


TPZVTKGraphMesh::TPZVTKGraphMesh(TPZCompMesh *cmesh, int dimension, TPZVTKGraphMesh *graph) :
TPZGraphMesh(cmesh, dimension, graph->fMaterialIds, graph->ScalarNames(),graph->VecNames(), graph->TensorNames()) {
	fNumCases = graph->fNumCases;
	fNumSteps = graph->fNumSteps;
	fStyle = EVTKStyle;
}

TPZVTKGraphMesh::TPZVTKGraphMesh(TPZCompMesh *cmesh, int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames): TPZGraphMesh(cmesh, dimension, matids,scalnames,vecnames,tensnames) {
    fNumCases = 0;
    fNumSteps = 0;
    fStyle = EVTKStyle;
    fVecNames = vecnames;
    fScalarNames = scalnames;
    fTensorNames = tensnames;
}

void TPZVTKGraphMesh::DrawMesh(int numcases) {
	
	fNumCases = numcases;
	fNumSteps = 0;
}

void TPZVTKGraphMesh::DrawSolution(int step, REAL time){

  std::set<int> matids = MaterialIds(); /// partial solution
  if(matids.size() == 0) {
    cout << __PRETTY_FUNCTION__
         <<"\nno material found\n";
    return;
  }
  
  //choose a random material for checking variables' names
  TPZMaterial * matp = fCompMesh->FindMaterial(*(matids.begin()));
  if (!matp) DebugStop();
    
	if(fOutFile.is_open()){
    fOutFile.close();
  }
  
	{
		std::string fn =
      fFileName.substr(0,fFileName.size()-4) + ".scal_vec." + std::to_string(step) + ".vtk";
		fOutFile.open(fn.c_str());
	}
  
	(fOutFile) << "# vtk DataFile Version 3.0" << endl;
	(fOutFile) << "Arquivo gerado por PZ" << endl;
	(fOutFile) << "ASCII\n" << endl;
	(fOutFile) << "DATASET UNSTRUCTURED_GRID" << endl;

  {
    TPZSimpleTimer timer("DrawNodes");
    DrawNodes();
  }

  {
    TPZSimpleTimer timer("DrawConnectivity");
    DrawConnectivity(ECube);
  }
  
	const int numscal = fScalarNames.NElements();
	const int numvec = fVecNames.NElements();
  const int numtens = fTensorNames.NElements();
  
	if(numscal || numvec || numtens){(fOutFile) << "POINT_DATA " << NPoints() << endl;}
	if(numscal){
    TPZSimpleTimer timer("scal");
    TPZManVector<int> scalind(0);
    scalind.Fill(-1,0,numscal);
    scalind.Resize(numscal);
    for(auto n=0; n<numscal; n++) {
      scalind[n] = matp->VariableIndex( fScalarNames[n]);
      if (scalind[n] == -1) {
        std::cout << fScalarNames[n] << " not recognized as post processing name\n";
      }
    }
    for(auto n=0; n<numscal; n++){
      if (scalind[n] != -1) {
        TPZSimpleTimer timer("scal_each");
        (fOutFile) << "SCALARS " << fScalarNames[n] << " float" << endl << "LOOKUP_TABLE default\n";
        const int64_t nnod = fNodeMap.NElements();
        for(auto i=0;i<nnod;i++) {
          TPZGraphNode *node = &fNodeMap[i];
          if(node) node->DrawSolution(scalind[n], EVTKStyle);
        }
        (fOutFile) << std::endl;
      }
    }
  }
	if(numvec){
    TPZSimpleTimer timer("vec");
    TPZManVector<int> vecind(0);
    vecind.Resize(numvec);
    vecind.Fill(-1,0,numvec);
    for(auto n=0; n<numvec; n++) {
      vecind[n] = matp->VariableIndex(fVecNames[n]);
      if (vecind[n] == -1) {
        std::cout << "Post processing vector name " << fVecNames[n]
                  << " not found\n";
      }
    }
    for(auto n=0; n<numvec; n++){
      if(vecind[n] != -1){
        TPZSimpleTimer timer("vec_each");
        (fOutFile) << "VECTORS " << fVecNames[n] << " float" << std::endl;
        const int64_t nnod = fNodeMap.NElements();
        for(auto i=0;i<nnod;i++) {
          TPZGraphNode *node = &fNodeMap[i];
          if(node) node->DrawSolution(vecind[n], EVTKStyle);
        }
        (fOutFile) << std::endl;
      }
    }
  }
  if(numtens){
    TPZSimpleTimer timer("tens");
    TPZVec<int> tensind(0);
    tensind.Resize(numtens);
    tensind.Fill(-1,0,numtens);
    for(auto n=0; n<numtens; n++) {
      tensind[n] = matp->VariableIndex(fTensorNames[n]);
      if(tensind[n] == -1){
        std::cout << "Post processing name " << fTensorNames[n] << " not found\n";
        DebugStop();
      }
    }
    for(auto n=0; n<numtens; n++){
      TPZSimpleTimer timer("tens_each");
      (fOutFile) << "TENSORS " << fTensorNames[n] << " float" << std::endl;
      const int64_t nnod = fNodeMap.NElements();
      for(auto i=0;i<nnod;i++) {
        TPZGraphNode *node = &fNodeMap[i];
        if(node) node->DrawSolution(tensind[n], EVTKStyle);
      }
      (fOutFile) << std::endl;
    }
  }
	fOutFile.close();
	fNumSteps++;
}


void TPZVTKGraphMesh::SequenceNodes(){
	TPZGraphMesh::SequenceNodes();
}

void TPZVTKGraphMesh::DrawNodes(){

	const int64_t nnod = fNodeMap.NElements();//total number of nodes, valid or not
  int64_t nn = 0L;//?
  
	for(auto i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n){
      nn += n->NPoints();
    }
	}
	(fOutFile) << "POINTS ";
	(fOutFile) << nn <<  " float" << endl;

  int64_t count{0};
	for(auto i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n){
      n->DrawCo(EVTKStyle);
    }
	}
}

void TPZVTKGraphMesh::DrawConnectivity(MElementType type) {
	
	int64_t nel = fElementList.NElements();
	if(!nel) return;
	TPZGraphEl *el = (TPZGraphEl *) fElementList[0];
	(fOutFile) << std::endl;
	(fOutFile) << "CELLS ";
	int64_t nelem = 0;
	int64_t nint = 0;
	int64_t i;
	for(i=0;i<nel;i++) {
		el = (TPZGraphEl *) fElementList[i];
		if(el)
		{
			nelem += el->NElements();
			nint += (el->NNodes() + 1) * el->NElements();
		}
 	}
	(fOutFile) << nelem << " " << nint << std::endl;
	for(i=0;i<nel;i++) {
		el = (TPZGraphEl *) fElementList[i];
		if(el) el->Connectivity(EVTKStyle);
 	}
	(fOutFile) << std::endl;
	(fOutFile) << "CELL_TYPES " << nelem << std::endl;
	for(i=0;i<nel;i++) {
		el = (TPZGraphEl *) fElementList[i];
		if(el)
		{
			int64_t j, nElPerEl = el->NElements();
			for(j = 0; j < nElPerEl; j++)(fOutFile) << el->ExportType(EVTKStyle) << std::endl;
		}
 	}
	(fOutFile) << std::endl;
}

void TPZVTKGraphMesh::DrawSolution(TPZBlock &/*Sol*/) {
    cout << "TPZMVGraphMesh::DrawSolution not Implemented\n";
}

void TPZVTKGraphMesh::DrawSolution(char * /*var = 0*/) {
    cout << "TPZMVGraphMesh::DrawSolution not Implemented\n";
}
