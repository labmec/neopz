/**
 * @file
 * @brief Contains the implementation of the TPZReadTetGen methods. 
 */

#include "pzreadtetgen.h"
#include "pzgmesh.h"
#include "TPZTimer.h"

#include <iostream>
#include <fstream>

TPZReadTetGen::TPZReadTetGen(){
	
}//method

TPZReadTetGen::~TPZReadTetGen(){
	
}//method

TPZGeoMesh * TPZReadTetGen::Process(std::string NodeFileName, std::string FaceFileName, std::string TetraFileName){
	TPZGeoMesh * gmesh = new TPZGeoMesh();
	int nnodes, nfaces, nvols;
	bool check;
	
	TPZMultiTimer time(4);
	time.processName(0) = "Reading nodes";
	time.processName(1) = "Reading faces";
	time.processName(2) = "Reading tetrahedras";
	time.processName(3) = "BuildConnectivity";
	
	time.start(0);
	check = this->ProcessNodes(NodeFileName, *gmesh, nnodes);
	time.stop(0);
	if (check == false){
		delete gmesh;
		return NULL;
	}
	
	time.start(1);
	check = this->ProcessFaces(FaceFileName, *gmesh, nfaces);
	time.stop(1);
	if (check == false){
		delete gmesh;
		return NULL;
	}
	
	time.start(2);
	check = this->ProcessTetra(TetraFileName, *gmesh, nvols);
	time.stop(2);
	if (check == false){
		delete gmesh;
		return NULL;
	}
    
	time.start(3);
	gmesh->BuildConnectivity();
	time.stop(3);
	
	std::cout << time << std::endl;
	
	return gmesh;
	
}//method

bool TPZReadTetGen::ProcessNodes(std::string NodeFileName,  TPZGeoMesh &gmesh, int & numbernodes){
	std::ifstream NodeFile(NodeFileName.c_str());
	this->fNodeIndices.clear();
	
	int ID, index, nnodes, dim, nattributes, BoundMarkers, bcmark;
	REAL X, Y, Z, attr;
	TPZManVector<REAL,3> coord(3);
	
	NodeFile >> nnodes >> dim >> nattributes >> BoundMarkers;  
	numbernodes = nnodes;
	
	if (dim != 3){
		std::cout << __PRETTY_FUNCTION__ << " - dim must be equal to 3" << std::endl;
		return false;  
	}//if
	
	for(int i = 0; i < nnodes; i++){
		NodeFile >> ID >> X >> Y >> Z;
		if (nattributes) NodeFile >> attr;
		if (BoundMarkers) NodeFile >> bcmark;
		
		index = gmesh.NodeVec().AllocateNewElement();
		coord[0] = X;
		coord[1] = Y;
		coord[2] = Z;
		gmesh.NodeVec()[index].Initialize(ID, coord, gmesh);
		
		this->fNodeIndices[ ID ] = index;        
	}//for i
	
	return true;
	
}//method

bool TPZReadTetGen::ProcessFaces(std::string FaceFileName,  TPZGeoMesh &gmesh, int & numberfaces){
	
	std::ifstream FaceFile(FaceFileName.c_str());
	
	int nfaces, BoundMarker, ID, mat;
	int NO1, NO2, NO3, index;
	TPZManVector<int,3> nodind(3);
	//TPZGeoEl * gel;
	
	FaceFile >> nfaces >> BoundMarker;
	
	numberfaces = nfaces;
	
	for(int i = 0; i < nfaces; i++){
		FaceFile >> ID >> NO1 >> NO2 >> NO3;
		if (BoundMarker) FaceFile >> mat;
		else mat = -1;
		nodind[0] = this->fNodeIndices[ NO1 ];
		nodind[1] = this->fNodeIndices[ NO2 ];
		nodind[2] = this->fNodeIndices[ NO3 ];
		/*gel = */gmesh.CreateGeoElement(ETriangle, nodind, mat, index);
	}//for i
	
	return true;
	
}//method

bool TPZReadTetGen::ProcessTetra(std::string TetraFileName, TPZGeoMesh &gmesh, int & numbervols){
	std::ifstream TetraFile(TetraFileName.c_str());
	int nvols, n, mat, nattr, ID, NO1, NO2, NO3, NO4, index;
	TPZManVector<int,4> nodind(4);
	//TPZGeoEl * gel;
	TetraFile >> nvols >> n >> nattr;
	
	numbervols = nvols;
	
	if (n != 4) std::cout << __PRETTY_FUNCTION__ << " - tetrahedra must have only four nodes" << std::endl;
	
	for(int i = 0; i < nvols; i++){
		TetraFile >> ID >> NO1 >> NO2 >> NO3 >> NO4;
		if (nattr) TetraFile >> mat;
		else mat = 9;
		nodind[0] = this->fNodeIndices[ NO1 ];
		nodind[1] = this->fNodeIndices[ NO2 ];
		nodind[2] = this->fNodeIndices[ NO3 ];
		nodind[3] = this->fNodeIndices[ NO4 ];
		/*gel = */gmesh.CreateGeoElement(ETetraedro, nodind, mat, index);    
		
	}//for i
	return true;
	
}//method
