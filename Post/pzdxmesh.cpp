/**
 * @file
 * @brief Contains the implementation of the TPZDXGraphMesh methods. 
 */
//$Id: pzdxmesh.cpp,v 1.16 2010-01-12 12:13:18 caju Exp $

#include "pzdxmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzgraphel.h"
#include "pztrigraph.h"
#include "pzgraphnode.h"
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"

using namespace std;

TPZDXGraphMesh::TPZDXGraphMesh(TPZCompMesh *cmesh, int dimension, TPZAutoPointer<TPZMaterial> mat, const TPZVec<std::string> &scalarnames, const TPZVec<std::string> &vecnames) :
TPZGraphMesh(cmesh,dimension,mat) {
	SetNames(scalarnames,vecnames);
	fNextDataField = 1;
	fStyle = EDXStyle;
	//	int index = 0;
	TPZCompEl *ce = FindFirstInterpolatedElement(cmesh,dimension);
	fElementType = "noname";
	if(ce) {
		int type = ce->Type();
		if(type == EOned)               fElementType = "lines";
		if(type == ETriangle)           fElementType = "triangles";
		if(type == EQuadrilateral)      fElementType = "quads";
		if(type == ECube)               fElementType = "cubes";
		if(type == EPrisma)             fElementType = "cubes";
		TPZGeoEl *gel = ce->Reference();
		if( type == EDiscontinuous || (type == EAgglomerate && gel) ){
			int nnodes = gel->NNodes();
			if(nnodes==4 && dimension==2) fElementType = "quads";
			if(nnodes==3 && dimension==2) fElementType = "triangles";
			if(dimension==3)              fElementType = "cubes";
		}
	}
	fNumCases = 0;
	fNumConnectObjects[0] = 1;
	fNumConnectObjects[1] = 1;
	fNumConnectObjects[2] = 1;
	fNormalObject = 0;
	
	for(int i = 0; i < 8; i++) fElConnectivityObject[i] = -1;
	
}

TPZDXGraphMesh::TPZDXGraphMesh(TPZCompMesh *cmesh,int dim,TPZDXGraphMesh *graph,TPZAutoPointer<TPZMaterial> mat) :
TPZGraphMesh(cmesh,dim,mat) {
	if(!mat) fMaterial = graph->fMaterial;
	fNextDataField = graph->fNextDataField;
	fStyle = EDXStyle;
	fElementType = "noname";
	fElementType[0] = '\0';
	fElementType = graph->fElementType;
	fNumCases = graph->fNumCases;
	int i;
	for(i=0;i<3;i++) {
		fNumConnectObjects[i] = 1;
		fFirstFieldValues[i] = graph->fFirstFieldValues[i];
	}
	for(;i<8;i++) fNumConnectObjects[i] = 0;
	fTimes = graph->fTimes;
	fFileName = graph->fFileName;
	/**independing the transfered output file from old graphmesh */
	graph->fOutFile.close();
	/** Load output file from old graphmesh to new graphmesh */
	fOutFile.open(fFileName.c_str());
	
	for(int i = 0; i < 8; i++) fElConnectivityObject[i] = -1;
	
}

TPZDXGraphMesh::~TPZDXGraphMesh() {
	Close();
}

int TPZDXGraphMesh::NNodes() {
	if(fElementType == "noname")    return 0;
	if(fElementType == "lines")     return 2;
	if(fElementType == "triangles") return 3;
	if(fElementType == "quads")     return 4;
	if(fElementType == "cubes")     return 8;
	return -1;
}

std::string TPZDXGraphMesh::ElementName() {
	return fElementType;
}

void TPZDXGraphMesh::DrawMesh(int numcases) {
	
	int dim = Material()->Dimension();
	int dim1 = dim-1;
	MElementType eltypes[] = {ENoType, EOned, EQuadrilateral, ETriangle,ECube};
	int firsttype[4] = {1,2,4,0}, lasttype[4] = {2,4,5,0};
	int numnod[] = {0,2,4,3,8};
	std::string elname[] = {"noname","lines","quads","triangles","cubes"};
	int object = fNextDataField;
	fNodePosObject[dim1] = object;
	long nn = NPoints();
	//      int imax = (1<<fResolution);
	
	//field nodes / positions
	(fOutFile) <<  "object " << object << " class array type float rank 1 shape 3 items "
	<< nn << " data follows"<< endl;
	DrawNodes();
	(fOutFile) << "attribute \"dep\" string \"positions\"" << endl;
	(fOutFile) << "#" << endl;
	object++;
	
	//   if(dim==1) {
	// 	//normal vectors
	// 	fNormalObject = object;
	// 	(fOutFile) <<  "object " << object << " class array type float rank 1 shape 3 items "
	// 	      << nn << " data follows"<< endl;
	// 	DrawNormals(nn);
	// 	(fOutFile) << "attribute \"dep\" string \"positions\"" << endl;
	// 	(fOutFile) << "#" << endl;
	//   object++;
	//   } else {
	fNormalObject = -1;
	//  }
	
	//  field connectivity
	int it;
	for(it=firsttype[dim1]; it<lasttype[dim1]; it++) {
		MElementType type = eltypes[it];
		int nel = NElements(type);
		if(nel) {
			//      fNumConnectObjects[dim1]++;
			fElConnectivityObject[type] = object;
			(fOutFile) << "object " << object << " class array type int rank 1 shape "
			<< numnod[it];
			(fOutFile) << " items "
			<< NElements(type) << " data follows"<<endl;
			object++;
			DrawConnectivity(type);
			
			(fOutFile) << "attribute \"element type\" string ";
			(fOutFile) << "\"" << elname[it] << "\"" << endl;
			(fOutFile) << "attribute \"ref\" string \"positions\"" << endl;
			(fOutFile) << "#" << endl;
		} else {
			fElConnectivityObject[type] = 0;
		}
	}
	fNumCases = numcases;
	//  long zero = 0;
	fNextDataField = object;
	//fFirstFieldValues[0] = new TPZVec<int>(numcases,zero);
	//fFirstFieldValues[1] = new TPZVec<int>(numcases,zero);
	//fFirstFieldValues[2] = new TPZVec<int>(numcases,zero);
	//fTimes = new TPZVec<REAL>(numcases,0.);
}

void TPZDXGraphMesh::DrawSolution(int step, REAL time){//0,
	//								  TPZVec<char *> &scalarnames, TPZVec<char *> &vectornames) {
	
	TPZAutoPointer<TPZMaterial> matp = Material();
	int i,nel;
	if(!matp) return;
	int dim = matp->Dimension();
	int dim1 = dim-1;
	
	int numscal = fScalarNames.NElements();
	int numvec = fVecNames.NElements();
	TPZVec<int> scalind(0);
	TPZVec<int> vecind(0);
	scalind.Resize(numscal);
	vecind.Resize(numvec);
	scalind.Fill(-1);
	vecind.Fill(-1);
	int n;
	for(n=0; n<numscal; n++) {
		scalind[n] = matp->VariableIndex( fScalarNames[n]);
	}
	for(n=0; n<numvec; n++) {
		vecind[n] = matp->VariableIndex( fVecNames[n]);
	}
	fFirstFieldValues[dim1].Push(fNextDataField+1);
	fTimes.Push(time);
	cout << "TPZDXGraphMesh.DrawSolution step = " << step << " time = " << time << endl;
	long numpoints = NPoints();
	TPZVec<int> scal(1);
	for(n=0; n<numscal; n++) {
		(fOutFile) << "object " << fNextDataField << " class array type float rank 0 items "
		<< numpoints << " data follows " << endl;
		scal[0] = scalind[n];
		nel = fNodeMap.NElements();
		for(i=0;i<nel;i++) {
			TPZGraphNode *np = &fNodeMap[i];
			if(np) np->DrawSolution(scal, fStyle);
		}
		(fOutFile) << "attribute \"dep\" string \"positions\"" << endl;
		(fOutFile) << "#" << endl;
		fNextDataField ++;
		
		
		//+++++ Cubes    
		if(dim == 3 && (NNodes() == 8 || NNodes() == 6)) {
			(fOutFile) << "object " << (fNextDataField) << " class field" << endl;
			(fOutFile) << "component \"data\" value " << (fNextDataField-1) << endl;
			(fOutFile) << "component \"positions\" value " << fNodePosObject[dim1] << endl;
			(fOutFile) << "component \"connections\" value " << fElConnectivityObject[ECube] << endl;
			(fOutFile) << "attribute \"name\" string \"" << (std::string) fScalarNames[n]
			<< step << (0) << "\"" << endl;
			(fOutFile) << "#" << endl;
			fNextDataField ++;
		}
		
		if(dim == 2 && (NNodes() == 4 || NNodes() == 3)) {
			(fOutFile) << "object " << (fNextDataField) << " class field" << endl;
			(fOutFile) << "component \"data\" value " << (fNextDataField-1) << endl;
			(fOutFile) << "component \"positions\" value " << fNodePosObject[dim1] << endl;
			(fOutFile) << "component \"connections\" value " << fElConnectivityObject[EQuadrilateral] << endl;
			(fOutFile) << "attribute \"name\" string \"" << (std::string) fScalarNames[n]; (fOutFile).flush();
			(fOutFile) << step; (fOutFile).flush();
			(fOutFile) << (0); (fOutFile).flush();
			(fOutFile) << "\""; (fOutFile).flush();
			(fOutFile) << endl; (fOutFile).flush();
			(fOutFile) << "#" << endl; (fOutFile).flush();
			fNextDataField ++;
		}
		if(dim == 1) {
			(fOutFile) << "object " << (fNextDataField) << " class field" << endl;
			(fOutFile) << "component \"data\" value " << (fNextDataField-1) << endl;
			(fOutFile) << "component \"positions\" value " << fNodePosObject[dim1] << endl;
			(fOutFile) << "component \"connections\" value " << fElConnectivityObject[EOned] << endl;
			if(fNormalObject > 0) (fOutFile) << "component \"normals\" value " << fNormalObject << endl;
			(fOutFile) << "attribute \"name\" string \"" << (std::string) fScalarNames[n]
			<< step << (0) << "\"" << endl;
			(fOutFile) << "#" << endl;
			fNextDataField ++;
		}
	}
	TPZVec<int> vec(1);
	for(n=0; n<numvec; n++) {
		(fOutFile) << "object " << fNextDataField << " class array type float rank 1 shape " <<
		matp->NSolutionVariables(vecind[n]) << " items "
		<< numpoints << " data follows " << endl;
		vec[0] = vecind[n];
		nel = fNodeMap.NElements();
		for(i=0;i<nel;i++) {
			TPZGraphNode *np = &fNodeMap[i];
			if(np) np->DrawSolution(vec, fStyle);
		}
		(fOutFile) << "attribute \"dep\" string \"positions\"" << endl;
		(fOutFile) << "#" << endl;
		fNextDataField ++;
		//+++++ Cubes        
		if(dim == 3 && fElConnectivityObject[ECube]) {
			(fOutFile) << "object " << (fNextDataField) << " class field" << endl;
			(fOutFile) << "component \"data\" value " << (fNextDataField-1) << endl;
			(fOutFile) << "component \"positions\" value " << fNodePosObject[dim1] << endl;
			(fOutFile) << "component \"connections\" value " << fElConnectivityObject[ECube] << endl;
			(fOutFile) << "attribute \"name\" string \"" << (std::string) fVecNames[n]
			<< step << (0) << "\"" << endl;
			(fOutFile) << "#" << endl;
			fNextDataField ++;
		}
		//|\|\|\|\|\|\|\|\|\|\|
		if(dim == 2 && fElConnectivityObject[EQuadrilateral]) {
			(fOutFile) << "object " << (fNextDataField) << " class field" << endl;
			(fOutFile) << "component \"data\" value " << (fNextDataField-1) << endl;
			(fOutFile) << "component \"positions\" value " << fNodePosObject[dim1] << endl;
			(fOutFile) << "component \"connections\" value " << fElConnectivityObject[EQuadrilateral] << endl;
			(fOutFile) << "attribute \"name\" string \"" << (std::string) fVecNames[n]
			<< step << (0) << "\"" << endl;
			(fOutFile) << "#" << endl;
			fNextDataField ++;
		}
		if(dim == 2 && fElConnectivityObject[ETriangle]) {
			(fOutFile) << "object " << (fNextDataField) << " class field" << endl;
			(fOutFile) << "component \"data\" value " << (fNextDataField-1) << endl;
			(fOutFile) << "component \"positions\" value " << fNodePosObject[dim1] << endl;
			(fOutFile) << "component \"connections\" value " << fElConnectivityObject[ETriangle] << endl;
			(fOutFile) << "attribute \"name\" string \"" << (std::string) fVecNames[n]
			<< step << (1) << "\"" << endl;
			(fOutFile) << "#" << endl;
			fNextDataField ++;
		}
		if(dim == 1) {
			(fOutFile) << "object " << (fNextDataField) << " class field" << endl;
			(fOutFile) << "component \"data\" value " << (fNextDataField-1) << endl;
			(fOutFile) << "component \"positions\" value " << fNodePosObject[dim1] << endl;
			(fOutFile) << "component \"connections\" value " << fElConnectivityObject[EOned] << endl;
			(fOutFile) << "attribute \"name\" string \"" << (std::string) fVecNames[n]
			<< step << (0) << "\"" << endl;
			(fOutFile) << "#" << endl;
			fNextDataField ++;
		}
	}
}


void TPZDXGraphMesh::Close(){
	
	int n;
	int dim = fDimension;
	int dim1 = dim-1;
	fNumCases = fFirstFieldValues[dim1].NElements();
	int numscal = fScalarNames.NElements();
	int numvec = fVecNames.NElements();
	int ist;
	//  TPZStack<int> *FV = &fFirstFieldValues[dim1];
	TPZStack<int> &FVR = fFirstFieldValues[dim1];
	cout << "DxMesh finalizing\n";
	//cout << " fTimes = " << (*fTimes)[0] << ' ' << (*fTimes)[1] << endl;
	for(n=0; n<numscal; n++) {
		(fOutFile) << "object \"" <<  fScalarNames[n] << "\" class series" << endl;
		for(ist =0; ist<fNumCases; ist++) {
			int icon;
			for(icon=0; icon< fNumConnectObjects[dim1]; icon++) {
				long FVRval = FVR[ist]+icon;
				(fOutFile) << "member " << ist << " position " << (fTimes)[ist]
				<< " value " << FVRval << endl;
			}
			fFirstFieldValues[dim1][ist] += 1+fNumConnectObjects[dim1];
		}
		(fOutFile) << "#" << endl;
	}
	for(n=0; n<numvec; n++) {
		(fOutFile) << "object \"" << fVecNames[n] << "\" class series" << endl;
		for(ist =0; ist<fNumCases; ist++) {
			int icon;
			for(icon=0; icon< fNumConnectObjects[dim1]; icon++) {
				long FVRval = FVR[ist]+icon;
				(fOutFile) << "member " << ist << " position " << (fTimes)[ist]
				<< " value " << FVRval << endl;
			}
			fFirstFieldValues[dim1][ist] += 1+fNumConnectObjects[dim1];
		}
		(fOutFile) << "#" << endl;
	}
	(fOutFile) << "object \"ALL\" class group\n";
	for(n=0; n<numscal; n++) {
		(fOutFile) << "member \"" << fScalarNames[n] << "\" value \"" << fScalarNames[n] << '\"' << endl;
	}
	for(n=0; n<numvec; n++) {
		(fOutFile) << "member \"" << fVecNames[n] << "\" value \"" << fVecNames[n] << '\"' << endl;
	}
	(fOutFile) << "end\n";
	fNumCases = 0;
	fNumConnectObjects[0] = 1;
	fNumConnectObjects[1] = 1;
	fNumConnectObjects[2] = 1;
	fNormalObject = 0;
	fTimes.Resize(0);
	fFirstFieldValues[0].Resize(0);
	fFirstFieldValues[1].Resize(0);
	fFirstFieldValues[2].Resize(0);
	fNextDataField = 1;
	
}


void TPZDXGraphMesh::DrawSolution(char * var)
{
	
    //int nmat = fCompMesh->MaterialVec().NElements();
    TPZAutoPointer<TPZMaterial> matp = Material();
    int i,varind;
    varind = matp->VariableIndex(var);
    TPZVec<int> vec(1);
    (fOutFile) << "object " << fNextDataField << " class array type float rank 1 shape " <<
	matp->NSolutionVariables(varind) << " items " << NPoints() << " data follows " << endl;
    vec[0] = varind;
    int nel = fCompMesh->ConnectVec().NElements();
    for(i=0;i<nel;i++) {
		TPZGraphNode n = fNodeMap[i];
		n.DrawSolution(vec,fStyle);
    }
    (fOutFile) << "attribute \"dep\" string \"positions\"" << endl;
    (fOutFile) << "#" << endl;
    (fOutFile) << "object " << (fNextDataField+1) << " class field" << endl;
    (fOutFile) << "component \"data\" value " << fNextDataField << endl;
	//    (fOutFile) << "component \"positions\" value " << fNodeCoField << endl;
	//    (fOutFile) << "component \"connections\" value " << fConnectField << endl;
    (fOutFile) << "attribute \"name\" string \"" << var << fNextDataField <<
	"\"" << endl;
    (fOutFile) << "#" << endl;
    fNextDataField += 2;
	
}

void TPZDXGraphMesh::DrawSolution(TPZBlock<REAL> &/*bl*/)
{
	/*
	 Pix i = fCompMesh->MaterialVec().first();
	 TPZMaterial *matp = (TPZMaterial *) fCompMesh->MaterialVec().contents(i);
	 (fOutFile) << "object " << fNextDataField << " class array type float rank 1 shape " <<
	 matp->NumVariables() << " items "
	 << NPoints() << " data follows " << endl;
	 i = fConnectVec.first();
	 while(i) {
	 TGraphNode *n = (TGraphNode *) fConnectVec.contents(i);
	 n->DrawSolution(bl,fStyle);
	 fConnectVec.next(i);
	 }
	 (fOutFile) << "attribute \"dep\" string \"positions\"" << endl;
	 (fOutFile) << "#" << endl;
	 (fOutFile) << "object " << (fNextDataField+1) << " class field" << endl;
	 (fOutFile) << "component \"data\" value " << fNextDataField << endl;
	 (fOutFile) << "component \"positions\" value " << fNodeCoField << endl;
	 (fOutFile) << "component \"connections\" value " << fConnectField << endl;
	 (fOutFile) << "attribute \"name\" string \"" << fNextDataField <<
	 "\"" << endl;
	 (fOutFile) << "#" << endl;
	 fNextDataField += 2;
	 */
}

void TPZDXGraphMesh::DrawNormals(int numnorm) {
	int i;
	for(i=0; i<numnorm; i++) {
		(fOutFile) << "0. 1. 0." << endl;
	}
}

void TPZDXGraphMesh::SetFileName(const std::string &filename)
{
	TPZGraphMesh::SetFileName(filename);
	fOutFile.open(filename.c_str());
}

