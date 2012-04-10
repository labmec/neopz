/*
 * @file
 * @brief Contains the implementation of the functions grouping the processes to: \n
 * - Generate one-dimensional geometric mesh (gmesh)
 * - Generate a computational mesh (cmesh) associated with gmesh,
 * - Assembling and solving linear system for 1D linear diff. eq. on \f$[0.0, 1.0]\f$
 */

#include "PZ_Process.h"

using namespace std;

TPZGeoMesh *GeomMesh1D(int h,TPZVec<int> &matId,TPZVec<int> &bc,TPZVec<REAL> &xL,TPZVec<REAL> &xR)
{
	if(!matId.NElements() || bc.NElements() < 2)
		return NULL;
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	int Qnodes = 2;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolPoint(1);
	
	// Creating geometric nodes
	int id = 0;
	int ndiv = 1;
	REAL dx = fabs(xR[0] - xL[0])/ndiv;
	REAL pointco;
	for(int xi = 0; xi < Qnodes; xi++) {
		pointco = xL[0]+xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,pointco);//coord X
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	// Creating geometric elements (point and linear)
	id=0;
	TopolPoint[0] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc[0], *gmesh);
	id++;
	for (int eli=0; eli<ndiv; eli++) {
		TopolLine[0] = eli;
		TopolLine[1] = eli+1;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId[0],*gmesh);
		id++;
	}
	TopolPoint[0] = Qnodes-1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc[1],*gmesh);
	
	// Building the connectivity between elements as (element, side)
	gmesh->BuildConnectivity();
	
	// If necessary it refines level by level initial mesh
	UniformRefinement(h,gmesh);
	
	return gmesh;
}

void UniformRefinement(int h,TPZGeoMesh *gmesh) {
	//uniform refinement
	TPZVec<TPZGeoEl *> filhos;
	for (int ref = 0; ref < h; ref++) {
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++ ) {
			TPZGeoEl *gel = gmesh->ElementVec()[i];
			if(gel->Type() != EPoint) gel->Divide(filhos);  // You can to divide point element but it do nothing (only error message)
		}
	}
}

TPZCompMesh *CompMesh1D(TPZGeoMesh *gmesh,int p, TPZMaterial *material,TPZVec<int> &bc,TPZVec<int> &bcType) {
	if(!material || bc.NElements()<2 || bcType.NElements() != bc.NElements()) return NULL;
	int dim = 1;
	
	
	TPZAutoPointer<TPZMaterial> mat(material);
	
	// related to interpolation space
	TPZCompEl::SetgOrder(p);
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	
	// Related to boundary conditions
	//	REAL uN=1-cosh(1.)/sinh(1.);
	TPZFMatrix<REAL> val1(1,1,0.), val2(1,1,0.);
	if(!bcType[0])  // dirichlet
		val2.PutVal(0,0,0.0);
	TPZAutoPointer<TPZMaterial> BCond1 = material->CreateBC(mat, bc[0],bcType[0], val1, val2);
	cmesh->InsertMaterialObject(BCond1);
	
	if(!bcType[1])  // dirichlet
		val2.PutVal(0,0,0.0);
	TPZAutoPointer<TPZMaterial> BCond2 = material->CreateBC(mat, bc[1],bcType[1], val1, val2);
	cmesh->InsertMaterialObject(BCond2);
	
	//Adjusting data
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements(); 
	cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;
}

void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	// Symmetric case	
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
	an.Solution().Zero();
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt);
	an.SetSolver(step);
	an.Run();
	
	//	Nonsymmetric case
	//	TPZBandStructMatrix full(fCmesh);
	//	an.SetStructuralMatrix(full);
	//	TPZStepSolver step;
	//	step.SetDirect(ELU);
	//	an.SetSolver(step);
	//	an.Run();
}


// CASE 2D

#include "pzgengrid.h"

TPZGeoMesh *GeomMesh2D(int h,TPZVec<int> &matId,TPZVec<int> &bc,TPZVec<REAL> &x0,TPZVec<REAL> &x1,int elementType)
{
	if(!matId.NElements() || bc.NElements() < 2 || elementType < 0)
		return NULL;
	if(elementType)
		elementType = 1;
	
    TPZGeoMesh *gmesh = new TPZGeoMesh;
	
	TPZManVector<int> nx(2,2);   // subdivisions in X and in Y. 
	TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.5 and hy = 0.5  
	gen.SetElementType(elementType);       // type = 0 means rectangular elements  --> type = 1 means triangular elements
	gen.Read(gmesh);             // generating grid in gmesh

	// setting bc condition -2 (u = 0) from (0.0, 0.0) until (1.0, 0.0)
	x1[1] = 0.0;
	gen.SetBC(gmesh,x0,x1,-2);
	// setting bc condition -1 from (0.0, 1.0) until (1.0, 1.0)
	x0[0] = x0[1] = 1.;
	gen.SetBC(gmesh,x1,x0,-1);
	// setting bc condition -1 from (1.0, 1.0) until (0.0, 1.0)
	x1[0] = 0.;	x1[1] = 1.;
	gen.SetBC(gmesh,x0,x1,-1);
	// setting bc condition -1 from (0.0, 1.0) until (0.0, 0.0)
	x0[0] = x0[1] = 0.;
	gen.SetBC(gmesh, x1, x0, -1);

	// Building the connectivity between elements as (element, side)
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	// If necessary it refines level by level initial mesh
	UniformRefinement(h,gmesh);
	
	return gmesh;
}

TPZCompMesh *CompMesh(TPZGeoMesh *gmesh,int p, TPZMaterial *material,TPZVec<int> &bc,TPZVec<int> &bcType) {	
	if(!material || bc.NElements()<2 || bcType.NElements() != bc.NElements()) return NULL;
	int dim = 2;
	
	TPZAutoPointer<TPZMaterial> mat(material);
	
	// related to interpolation space
	TPZCompEl::SetgOrder(p);
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	
	// Related to boundary conditions
	//	REAL uN=1-cosh(1.)/sinh(1.);
	TPZFMatrix<REAL> val1(1,1,0.), val2(1,1,0.);
	if(!bcType[0])  // dirichlet
		val2.PutVal(0,0,0.);
	TPZAutoPointer<TPZMaterial> BCond1 = material->CreateBC(mat, bc[0],bcType[0], val1, val2);
	cmesh->InsertMaterialObject(BCond1);
	
	if(!bcType[1])  // dirichlet
		val2.PutVal(0,0,0.25);
	TPZAutoPointer<TPZMaterial> BCond2 = material->CreateBC(mat, bc[1],bcType[1], val1, val2);
	cmesh->InsertMaterialObject(BCond2);
	
	//Adjusting data
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements(); 
	cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;
}


// OUTPUT FORMAT

// Output as Mathematica format
void OutputMathematica(std::ofstream &outMath,int var,int pointsByElement,TPZCompMesh *cmesh) {
	int i, j, k, nnodes;
	int nelem = cmesh->ElementVec().NElements();
	int dim = cmesh->Dimension();   // Dimension of the model
	REAL w;
	if(var-1 < 0) var = 1;
	// Map to store the points and values 
	map<REAL,TPZVec<REAL> > Graph;
	TPZVec<REAL> tograph(4,0.);
	map<TPZVec<REAL>,REAL> Graphics;
	
	for(i=0;i<nelem;i++) {
		TPZCompEl *cel = cmesh->ElementVec()[i];
		TPZGeoEl *gel = cel->Reference();
		TPZInterpolationSpace * sp = dynamic_cast <TPZInterpolationSpace*>(cel);
		int nstates = cel->Material()->NStateVariables();
		// If var is higher than nstates of the element, go to next element
		if(var > nstates)
			continue;
		TPZVec<REAL> qsi(3,0.), sol(nstates,0.), outfem(3,0.);
		nnodes = gel->NNodes();
		if(pointsByElement < nnodes) pointsByElement = nnodes;
		for(j=0;j<gel->NNodes();j++) {
			// Get corners points to compute solution on
			gel->CenterPoint(j,qsi);
			sp->Solution(qsi,0,sol);
			cel->Reference()->X(qsi,outfem);
			// Jointed point coordinates and solution value on			
			for(k=0;k<3;k++) tograph[k] = outfem[k];
			tograph[k] = sol[var-1];
			Graph.insert(pair<REAL,TPZVec<REAL> >(outfem[0],tograph));
			Graphics.insert(pair<TPZVec<REAL>,REAL>(outfem,sol[var-1]));
			// If cel is point gets one point value
			if(cel->Type() == EPoint) {
				break;
			}
		}
		// If cel is point gets one point value
		if(cel->Type() == EPoint) continue;
		// Print another points using integration points
		TPZIntPoints *rule = NULL;
		int order = 1, npoints = 0;
		while(pointsByElement-(npoints+nnodes) > 0) {
			if(rule) delete rule;   // Cleaning unnecessary allocation
			int nsides = gel->NSides();
			// Get the integration rule to compute internal points to print, not to print
			rule = gel->CreateSideIntegrationRule(nsides-1,order);
			if(!rule) break;
			npoints = rule->NPoints();
			order += 2;
		}
		for(j=0;j<npoints;j++) {
			// Get integration points to get internal points
			rule->Point(j,qsi,w);
			sp->Solution(qsi,0,sol);
			cel->Reference()->X(qsi,outfem);
			// Jointed point coordinates and solution value on
			for(k=0;k<3;k++) tograph[k] = outfem[k];
			tograph[k] = sol[var-1];
			Graph.insert(pair<REAL,TPZVec<REAL> >(outfem[0],tograph));
			Graphics.insert(pair<TPZVec<REAL>,REAL>(outfem,sol[var-1]));
		}
	}
	
	// Printing the points and values into the Mathematica file
	outMath << "Saida = { ";
	if(dim<2) {
		map<REAL,TPZVec<REAL> >::iterator it;
		for(it=Graph.begin();it!=Graph.end();it++) {
			if(it!=Graph.begin()) outMath << ",";
			outMath << "{";
			for(j=0;j<dim;j++)
				outMath << (*it).second[j] << ",";
			outMath << (*it).second[3] << "}";
		}
		outMath << "}" << std::endl;
		// Choose Mathematica command depending on model dimension
		outMath << "ListPlot[Saida,Joined->True]"<< endl;
	}
	else {
		map<TPZVec<REAL>,REAL>::iterator it;
		for(it=Graphics.begin();it!=Graphics.end();it++) {
			if(it!=Graphics.begin()) outMath << ",";
			outMath << "{";
			for(j=0;j<dim;j++)
				outMath << (*it).first[j] << ",";
			outMath << (*it).second << "}";
		}
		outMath << "}" << std::endl;
		outMath << "ListPlot3D[Saida]"<< endl;
	}
}

// Output as VTK (Visualization Tool Kit) format
void OutputVTK(std::string &outVTK, TPZCompMesh *cmesh,TPZAnalysis &an) {
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
	vecnames[0] = "Derivate";
	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,outVTK);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}
