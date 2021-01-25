//
//  TPZGeoMeshBuilder.cpp
//  pz
//
//  Created by Omar Durán on 2/12/19.
//

#include "TPZGeoMeshBuilder.h"

void TPZGeoMeshBuilder::InsertNodes(TPZGeoMesh * gmesh, std::vector<std::size_t> & node_identifiers, std::vector<double> & coord) {
    
    int64_t n_nodes = node_identifiers.size();
    gmesh->NodeVec().Resize(n_nodes);
    gmesh->SetMaxNodeId(n_nodes-1);
    
    int64_t node_id;
    REAL nodecoordX,nodecoordY,nodecoordZ;
    
    /// Inserting nodes
    TPZGeoNode node_obj;
    for (int64_t inode = 0; inode < n_nodes; inode++) {
    
        node_id = node_identifiers[inode]-1; //  because pz is zero based.
        int pos = inode*3; // because the model is dwell in R3
        nodecoordX = coord[pos];
        nodecoordY = coord[pos+1];
        nodecoordZ = coord[pos+2];
        
        
        node_obj.SetNodeId(node_id);
        node_obj.SetCoord(0,nodecoordX);
        node_obj.SetCoord(1,nodecoordY);
        node_obj.SetCoord(2,nodecoordZ);
        gmesh->NodeVec()[node_id] = node_obj;
    }
}

TPZGeoEl* TPZGeoMeshBuilder::InsertElement(TPZGeoMesh * gmesh, int & physical_identifier, int & el_type, int & el_identifier, std::vector<int> & node_identifiers){
    
    TPZManVector <int64_t,15> Topology;
    int n_nodes = node_identifiers.size();
    Topology.Resize(n_nodes, 0);
    for (int k_node = 0; k_node<n_nodes; k_node++) {
        Topology[k_node] = node_identifiers[k_node]-1;
    }
    TPZGeoEl* gel = nullptr;
    el_identifier--;
    switch (el_type) {
        case 1:
        {   // Line
            gel = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 2:
        {
            // Triangle
            gel = new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (el_identifier, Topology, physical_identifier, *gmesh);
            
        }
            break;
        case 3:
        {
            // Quadrilateral
            gel = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (el_identifier, Topology, physical_identifier, *gmesh);
            
        }
            break;
        case 4:
        {
            // Tetrahedron
            gel = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (el_identifier, Topology, physical_identifier, *gmesh);
            
        }
            break;
        case 5:
        {
            // Hexahedra
            gel = new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 6:
        {
            // Prism
            gel = new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 7:
        {
            // Pyramid
            gel = new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 8:
        {
            // Quadratic Line
            gel = new TPZGeoElRefPattern< pzgeom::TPZQuadraticLine> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 9:
        {
            // Triangle
            gel = new TPZGeoElRefPattern< pzgeom::TPZQuadraticTrig> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 10:
        {
            // Quadrilateral
            gel = new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 11:
        {
            // Tetrahedron
            gel = new TPZGeoElRefPattern< pzgeom::TPZQuadraticTetra> (el_identifier, Topology, physical_identifier, *gmesh);
            
        }
            break;
        case 12:
        {
            // Hexahedra
            gel = new TPZGeoElRefPattern< pzgeom::TPZQuadraticCube> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 13:
        {
            // Prism
            gel = new TPZGeoElRefPattern< pzgeom::TPZQuadraticPrism> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        case 15:{
            // Point
            gel = new TPZGeoElement< pzgeom::TPZGeoPoint, pzrefine::TPZRefPoint> (el_identifier, Topology, physical_identifier, *gmesh);
        }
            break;
        default:
        {
            std::cout << "Element not impelemented." << std::endl;
            DebugStop();
        }
            break;
    }
    return gel;
}

int TPZGeoMeshBuilder::GetNumberofNodes(int & el_type){
    
    int n_nodes;
    switch (el_type) {
        case 1:
        {   // Line
            n_nodes = 2;
        }
            break;
        case 2:
        {
            // Triangle
            n_nodes = 3;
        }
            break;
        case 3:
        {
            // Quadrilateral
            n_nodes = 4;
        }
            break;
        case 4:
        {
            // Tetrahedron
            n_nodes = 4;
        }
            break;
        case 5:
        {
            // Hexahedra
            n_nodes = 8;
        }
            break;
        case 6:
        {
            // Prism
            n_nodes = 6;
        }
            break;
        case 7:
        {
            // Pyramid
            n_nodes = 5;
        }
            break;
        case 8:
        {
            // Quadratic Line
            n_nodes = 3;
        }
            break;
        case 9:
        {
            // Quadratic Triangle
            n_nodes = 6;
        }
            break;
        case 10:
        {
            // Quadratic Quadrilateral
            n_nodes = 8;
        }
            break;
        case 11:
        {
            // Quadratic Tetrahedron
            n_nodes = 10;
            
        }
            break;
        case 12:
        {
            // Quadratic Hexahedra
            n_nodes = 20;
        }
            break;
        case 13:
        {
            // Quadratic Prism
            n_nodes = 15;
        }
            break;
        case 15:{
            // Point
            n_nodes = 1;
        }
            break;
        default:
        {
            std::cout << "Element not impelemented." << std::endl;
            n_nodes = 0;
            DebugStop();
        }
            break;
    }
    
    return n_nodes;
}

void TPZGeoMeshBuilder::PrintGeometry(TPZGeoMesh * gmesh, std::string  & name){
    
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name   << name.c_str() << ".txt";
    vtk_name    << name.c_str() << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    gmesh->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}
