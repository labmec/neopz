/**
 * @file
 * @brief Contains the implementation of the TPZVTKGraphMesh methods. 
 */
/*
 *  pzvtkmesh.cpp
 *  NeoPZ
 *
 *  Created by Philippe Devloo on 04/12/08.
 *  Copyright 2008 UNICAMP. All rights reserved.
 *
 */

#include <sstream>
#include "pzvtkmesh.h"
#include "pzcmesh.h"
#include "pzmaterial.h"
#include "pzgraphnode.h"
#include "pzgraphel.h"

using namespace std;

TPZVTKGraphMesh::TPZVTKGraphMesh(TPZCompMesh *cmesh, int dimension, TPZAutoPointer<TPZMaterial> mat,
								 const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames) : TPZGraphMesh(cmesh, dimension, mat) {
	fNumCases = 0;
	fNumSteps = 0;
	fStyle = EVTKStyle;
	fVecNames = vecnames;
	fScalarNames = scalnames;
}

TPZVTKGraphMesh::TPZVTKGraphMesh(TPZCompMesh *cmesh, int dimension, TPZVTKGraphMesh *graph,TPZAutoPointer<TPZMaterial> mat) :
TPZGraphMesh(cmesh, dimension,mat) {
	if(!mat) fMaterial = graph->fMaterial;
	fNumCases = graph->fNumCases;
	fNumSteps = graph->fNumSteps;
	fStyle = EVTKStyle;
}

void TPZVTKGraphMesh::DrawMesh(int numcases) {
	
	fNumCases = numcases;
	fNumSteps = 0;
}

void TPZVTKGraphMesh::DrawSolution(int step, REAL time){
	
	TPZAutoPointer<TPZMaterial> matp = Material();
	if(!matp) {
		cout << "TPZMVGraphMesh no material found\n";
		return;
	}
	if(fOutFile.is_open())
	{
		fOutFile.close();
	}
	int n;
	{
		std::stringstream sout;
		sout << fFileName.substr(0,fFileName.size()-4) << ".scal_vec" << step << ".vtk";
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
	if(numscal || numvec)(fOutFile) << "POINT_DATA " << NPoints() << endl;
	if(numscal)
	{
		TPZVec<int> scalind(0);
		scalind.Fill(-1,0,numscal);
		scalind.Resize(numscal);
		for(n=0; n<numscal; n++) {
			scalind[n] = matp->VariableIndex( fScalarNames[n]);
		}
		for(n=0; n<numscal; n++)
		{
			(fOutFile) << "SCALARS " << fScalarNames[n] << " float" << endl << "LOOKUP_TABLE default\n";
			int nnod = fNodeMap.NElements(),i;
			for(i=0;i<nnod;i++) {
				TPZGraphNode *node = &fNodeMap[i];
				if(node) node->DrawSolution(scalind[n], EVTKStyle);
			}
			(fOutFile) << std::endl;
		}
	}
	if(numvec)
	{
		TPZVec<int> vecind(0);
		vecind.Resize(numvec);
		vecind.Fill(-1,0,numvec);
		for(n=0; n<numvec; n++) {
			vecind[n] = matp->VariableIndex(fVecNames[n]);
		}
		for(n=0; n<numvec; n++)
		{
			(fOutFile) << "VECTORS " << fVecNames[n] << " float" << std::endl;
			int nnod = fNodeMap.NElements(), i;
			for(i=0;i<nnod;i++) {
				TPZGraphNode *node = &fNodeMap[i];
				if(node) node->DrawSolution(vecind[n], EVTKStyle);
			}
			(fOutFile) << std::endl;
		}
	}
	fNumSteps++;
}


void TPZVTKGraphMesh::SequenceNodes(){
	TPZGraphMesh::SequenceNodes();
}

void TPZVTKGraphMesh::DrawNodes(){
	
	long nn = 0L;
	int nnod = fNodeMap.NElements();
	int i;
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
	
	int nel = fElementList.NElements();
	if(!nel) return;
	TPZGraphEl *el = (TPZGraphEl *) fElementList[0];
	(fOutFile) << std::endl;
	(fOutFile) << "CELLS ";
	int nelem = 0;
	int nint = 0;
	int i;
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
			int j, nElPerEl = el->NElements();
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
