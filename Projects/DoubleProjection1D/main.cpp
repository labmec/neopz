#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"

#include "tpzgeolinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"

#include <iostream>
#include <string>

#include <math.h>

using namespace std;

const int matId = 1;
const int bc1 = -1;
const int bc2 = -2;

const int dirichlet = 0;
const int neumann = 1;

//
 /*programa para resolver o problema de advecao-dufusao 1D*/
/*    - a.u''(x) + b.u'(x) = f(x), em  0< x <1*/
/* com : u(0) = g*, du/dn = 0 em x =1*/
//

TPZGeoMesh * MalhaGeom(int h, REAL xL, REAL  xR);
TPZGeoMesh * MalhaComp(TPZGeoMesh * gmesh, int p, REAL a, REAL b);

int main(int argc, char *argv[])
{	
	std::string logs("log4cxx.doubleprojection1d");
	InitializePZLOG(logs);
	// gRefDBase.InitializeRefPatterns();
	
	cout << "\n\n ===== bla bla bla ==="<<endl;
	MalhaGeom(2, 0., 20.);
	 
	return EXIT_SUCCESS;
}

TPZGeoMesh * MalhaGeom(int h, REAL xL, REAL  xR){
	
	int Qnodes = pow(2., h) + 1;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolPoint(1);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	//indice dos nos
	int nel = pow(2., h);
	REAL dx = fabs(xR - xL)/nel;
	REAL valx;
	int id = 0;
	for(int xi = 0; xi < Qnodes; xi++)
	{
		valx = xL+xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id =0;
	TopolPoint[0] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc1, *gmesh);
	id++;
	for (int eli=0; eli<nel; eli++) {
		TopolLine[0] = eli;
		TopolLine[1] = eli+1;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
		id++;
	}
	
	TopolPoint[0] = Qnodes-1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc2,*gmesh);
	
	gmesh->BuildConnectivity();
	
	ofstream arg("gmesh.txt");
	gmesh->Print(arg);
	
	return gmesh;
}
