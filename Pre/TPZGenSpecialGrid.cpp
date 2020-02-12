/**
 * @file
 * @brief Contains the implementation of the methods to generate polygonal grids approximating three dimensional special surfaces. 
 */

#include "pzmanvector.h"
#include "pzstack.h"
#include "pzgeoel.h"

#include "TPZGenSpecialGrid.h"

#include <math.h>

/** 
 * Function to generate a polygonal mesh approximating a sphere from a octahedron mesh (polygonal) based on number of uniform refinements
 */
TPZGeoMesh *TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(TPZVec<REAL> &center, REAL radius, int nUniformRefs) {
	// Initial mesh data
	// In this case is a octahedron
	// Octahedron has 6 nodes. We considerer that the nodes is on axes and radius is the distance between any node and origin.
	const int64_t nnode = 6;
	// Octahedron has 8 (triangular) faces
	const int64_t nelem = 8;
	
	// Initial nodes and initial triangular faces of the octahedron
    REAL initialcoord[nnode][3] = {{center[0]-radius,center[1],center[2]},{center[0],center[1]+radius,center[2]},{center[0],center[1],center[2]-radius},{center[0],center[1],center[2]+radius},{center[0],center[1]-radius,center[2]},{center[0]+radius,center[1],center[2]}};
    int indices[nelem][nnode] = {{3,4,5},{3,5,1},{3,1,0},{3,0,4},{4,0,2},{4,2,5},{2,0,1},{5,2,1}};
	
	// Geometric element vector
    TPZGeoEl *elvec[nelem], *gel;
	
	// Creating geometric initial mesh
    TPZGeoMesh *gmesh = new TPZGeoMesh();
	
    // Initializing nodes of the polygonal initial mesh
    int64_t node;
    for(node=0; node<nnode; node++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL> coord(3);
        coord[0] = initialcoord[node][0];
        coord[1] = initialcoord[node][1];
        coord[2] = initialcoord[node][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(node,coord,*gmesh);
    }
    // Creating triangular elements
    int64_t el, index;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(3);
        for(node=0; node<3; node++) nodind[node]=indices[el][node];
        elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
    }
	// Building connectivity for initial mesh - OCTAHEDRON
	gmesh->BuildConnectivity();
	
	//Loop making uniform refinement and changing coordinates of the nodes (projecting into the sphere) until tolerance is reached
	TPZManVector<REAL> baryparam(3,0.), barycenter(3,0.);
	for(int ii=0;ii<nUniformRefs;ii++) {
		// Make a uniform refinement
		UniformRefinement(1,gmesh,2);
		
		// Projecting all nodes into the sphere
		TPZVec<REAL> coordinates(3,0.);
		REAL norm;
		int i;
		for(node=0;node<gmesh->NNodes();node++) {
			TPZGeoNode *gnode = &gmesh->NodeVec()[node];
			gnode->GetCoordinates(coordinates);
			for(i=0;i<3;i++)
				coordinates[i] -= center[i];
			norm = sqrt(coordinates[0]*coordinates[0] + coordinates[1]*coordinates[1] + coordinates[2]*coordinates[2]);
			for(i=0;i<3;i++) coordinates[i] /= norm;
			for(i=0;i<3;i++) coordinates[i] += center[i];
			gnode->SetCoord(coordinates);
		}
	}
	
	// Before return, all the "father" elements are deleted, but the "son" elements are not
	for(el=0;el<gmesh->NElements();el++) {
		gel = gmesh->ElementVec()[el];
		if(!gel) continue;
		if(gel->HasSubElement()) {
			gel->ResetSubElements();
			gmesh->DeleteElement(gel);
		}
		else
			gel->SetFatherIndex(-1);
	}
	// The connectivity are updated
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	return gmesh;
}

/** 
 * Function to generate a polygonal mesh approximating a sphere from a octahedron mesh (polygonal) 
 */
TPZGeoMesh *TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(TPZVec<REAL> &center,REAL radius, REAL tol) {
	// Initial mesh data
	// In this case is a octahedron
	// Octahedron has 6 nodes. We considerer that the nodes is on axes and radius is the distance between any node and origin.
	const int64_t nnode = 6;
	// Octahedron has 8 (triangular) faces
	const int64_t nelem = 8;
	
	// Initial nodes and initial triangular faces of the octahedron
    REAL initialcoord[nnode][3] = {{center[0]-radius,center[1],center[2]},{center[0],center[1]+radius,center[2]},{center[0],center[1],center[2]-radius},{center[0],center[1],center[2]+radius},{center[0],center[1]-radius,center[2]},{center[0]+radius,center[1],center[2]}};
    int indices[nelem][nnode] = {{3,4,5},{3,5,1},{3,1,0},{3,0,4},{4,0,2},{4,2,5},{2,0,1},{5,2,1}};
	
	// Geometric element vector
    TPZGeoEl *elvec[nelem], *gel;
	
	// Creating geometric initial mesh
    TPZGeoMesh *gmesh = new TPZGeoMesh();
	
    // Initializing nodes of the polygonal initial mesh
    int64_t node;
    for(node=0; node<nnode; node++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL> coord(3);
        coord[0] = initialcoord[node][0];
        coord[1] = initialcoord[node][1];
        coord[2] = initialcoord[node][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(node,coord,*gmesh);
    }
    // Creating triangular elements
    int64_t el, index;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(3);
        for(node=0; node<3; node++) nodind[node]=indices[el][node];
        elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
    }
	// Building connectivity for initial mesh - OCTAHEDRON
	gmesh->BuildConnectivity();
	
	//Loop making uniform refinement and changing coordinates of the nodes (projecting into the sphere) until tolerance is reached
	TPZManVector<REAL> baryparam(3,0.), barycenter(3,0.);
	REAL dist = 1.;
	while(tol < dist) {
		// Make a uniform refinement
		UniformRefinement(1,gmesh,2);
		
		// Projecting all nodes into the sphere
		TPZVec<REAL> coordinates(3,0.);
		REAL norm;
		int i;
		for(node=0;node<gmesh->NNodes();node++) {
			TPZGeoNode *gnode = &gmesh->NodeVec()[node];
			gnode->GetCoordinates(coordinates);
			for(i=0;i<3;i++)
				coordinates[i] -= center[i];
			norm = sqrt(coordinates[0]*coordinates[0] + coordinates[1]*coordinates[1] + coordinates[2]*coordinates[2]);
			for(i=0;i<3;i++) coordinates[i] /= norm;
			for(i=0;i<3;i++) coordinates[i] += center[i];
			gnode->SetCoord(coordinates);
		}
		
		// Find element no null to calculate distance between barycenter and projection into sphere
		for(el=0;el<gmesh->NElements();el++) {
			gel = gmesh->ElementVec()[el];
			if(!gel || gel->HasSubElement()) continue;
			gel->CenterPoint(6,baryparam);
			gel->X(baryparam,barycenter);
			break;
		}
		dist = radius - sqrt((barycenter[0]-center[0])*(barycenter[0]-center[0]) + (barycenter[1]-center[1])*(barycenter[1]-center[1]) + (barycenter[2]-center[2])*(barycenter[2]-center[2]));
		if(dist < 0.0) {
			delete gmesh;
			gmesh = 0;
			return gmesh;
		} 
	}
	
	// Before return, all the "father" elements are deleted, but the "son" elements are not
	for(el=0;el<gmesh->NElements();el++) {
		gel = gmesh->ElementVec()[el];
		if(!gel) continue;
		if(gel->HasSubElement()) {
			gel->ResetSubElements();
			gmesh->DeleteElement(gel);
		}
		else
			gel->SetFatherIndex(-1);
	}
	// The connectivity are updated
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();

	return gmesh;
}

/** 
 * Make uniform refinement of the geometric mesh. This method finishs making the new connectivities for the elements (ResetConnectivities() and BuildConnectivity()) 
 */
void TPZGenSpecialGrid::UniformRefinement(int nUniformRefs,TPZGeoMesh *gmesh, const int dimelfordivision, bool allmaterial, const int matidtodivided) {
	TPZManVector<TPZGeoEl*> filhos;
	for(int Division=0; Division<nUniformRefs; Division++)
	{
		int64_t nels = gmesh->NElements();
		for(int64_t elem = 0; elem < nels; elem++)
		{    
			TPZGeoEl * gel = gmesh->ElementVec()[elem];
			if(!gel || gel->HasSubElement())
				continue;
			if(dimelfordivision > 0 && gel->Dimension() != dimelfordivision) continue;
			if(!allmaterial){
				if(gel->MaterialId() == matidtodivided){
					gel->Divide(filhos);
				}
			}
			else{
				gel->Divide(filhos);
			}
		}
	}
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
