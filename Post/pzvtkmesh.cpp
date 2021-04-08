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
	
    std::set<int> matid = MaterialIds();
    std::set<int> matids = MaterialIds(); /// partial solution
    if(matids.size() == 0) {
        cout << "TPZMVGraphMesh no material found\n";
        return;
    }
    set<int>::iterator it = matids.begin();
    TPZMaterial * matp = fCompMesh->FindMaterial(*it);
	if(fOutFile.is_open())
	{
		fOutFile.close();
	}
	int n;
	{
		std::stringstream sout;
		sout << fFileName.substr(0,fFileName.size()-4) << ".scal_vec." << step << ".vtk";
		fOutFile.open(sout.str().c_str());
	}
	(fOutFile) << "# vtk DataFile Version 3.0" << endl;
	(fOutFile) << "Arquivo gerado por PZ" << endl;
	(fOutFile) << "ASCII\n" << endl;
	(fOutFile) << "DATASET UNSTRUCTURED_GRID" << endl;
	DrawNodes();
	DrawConnectivity(ECube);
	int numscal = fScalarNames.NElements();
	int numvec = fVecNames.NElements();
    int numtens = fTensorNames.NElements();
	if(numscal || numvec)(fOutFile) << "POINT_DATA " << NPoints() << endl;
	if(numscal)
	{
		TPZManVector<int> scalind(0);
		scalind.Fill(-1,0,numscal);
		scalind.Resize(numscal);
		for(n=0; n<numscal; n++) {
			scalind[n] = matp->VariableIndex( fScalarNames[n]);
            if (scalind[n] == -1) {
                std::cout << fScalarNames[n] << " not recognized as post processing name\n";
            }
		}
		for(n=0; n<numscal; n++)
		{
            if (scalind[n] != -1)
            {
                (fOutFile) << "SCALARS " << fScalarNames[n] << " float" << endl << "LOOKUP_TABLE default\n";
                int64_t nnod = fNodeMap.NElements(),i;
                for(i=0;i<nnod;i++) {
                    TPZGraphNode *node = &fNodeMap[i];
                    if(node) node->DrawSolution(scalind[n], EVTKStyle);
                }
                (fOutFile) << std::endl;
            }
		}
	}
	if(numvec)
	{
		TPZManVector<int> vecind(0);
		vecind.Resize(numvec);
		vecind.Fill(-1,0,numvec);
		for(n=0; n<numvec; n++) {
			vecind[n] = matp->VariableIndex(fVecNames[n]);
            if(vecind[n] == -1)
            {
                std::cout << "Post processing vector name " << fVecNames[n] << " not found\n";
                
            }
		}
		for(n=0; n<numvec; n++)
		{
            if(vecind[n] != -1)
            {
                (fOutFile) << "VECTORS " << fVecNames[n] << " float" << std::endl;
                int64_t nnod = fNodeMap.NElements(), i;
                for(i=0;i<nnod;i++) {
                    TPZGraphNode *node = &fNodeMap[i];
                    if(node) node->DrawSolution(vecind[n], EVTKStyle);
                }
                (fOutFile) << std::endl;
            }
		}
	}
    if(numtens)
    {
        TPZVec<int> tensind(0);
        tensind.Resize(numtens);
        tensind.Fill(-1,0,numtens);
        for(n=0; n<numtens; n++) {
            tensind[n] = matp->VariableIndex(fTensorNames[n]);
            if(tensind[n] == -1)
            {
                std::cout << "Post processing name " << fTensorNames[n] << " not found\n";
                DebugStop();
            }
        }
        for(n=0; n<numtens; n++)
        {
            (fOutFile) << "TENSORS " << fTensorNames[n] << " float" << std::endl;
            int64_t nnod = fNodeMap.NElements(), i;
            for(i=0;i<nnod;i++) {
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
	
	int64_t nn = 0L;
	int64_t nnod = fNodeMap.NElements();
	int64_t i;
	for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) nn += n->NPoints();
	}
	(fOutFile) << "POINTS ";
	(fOutFile) << nn <<  " float" << endl;
	
	for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) n->DrawCo(EVTKStyle);
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
