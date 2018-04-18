//
//  TRMGmshReader.cpp
//  PZ
//
//  Created by Omar on 1/15/17.
//
//

#include "TRMGmshReader.h"

#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "tpzcube.h"


#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticcube.h"
#include "tpzquadratictetra.h"
#include "tpzgeoblend.h"

#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include <tpzarc3d.h>

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

TRMGmshReader::TRMGmshReader() {
    fVolNumber = 0;
    fBCNumber = 0;
    fProblemDimension = 0;
    fDimensionlessL = 1.0;
}//method

TRMGmshReader::~TRMGmshReader() {
    
}//method


TPZGeoMesh * TRMGmshReader::GeometricGmshMesh(std::string file_name)
{
    
    std::string string_temp;
    
    //  Mesh Creation
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    {
        
        // reading a general mesh information by filter
        std::ifstream read (file_name.c_str());
        
        while(read)
        {
            char buf[1024];
            read.getline(buf, 1024);
            std::string str(buf);
            
            if(str == "$MeshFormat" || str == "$MeshFormat\r")
            {
                read.getline(buf, 1024);
                std::string str(buf);
                std::cout << "Reading mesh format = " << str << std::endl;
                
            }
            
            if(str == "$PhysicalNames" || str == "$PhysicalNames\r" )
            {

                int64_t n_entities;
                read >> n_entities;
                int max_dimension = 0;
                
                int dimension, id;
                std::string name;
                std::pair<int, std::string> chunk;
                
                for (int64_t inode = 0; inode < n_entities; inode++) {

                    read.getline(buf, 1024);
                    read >> dimension;
                    read >> id;
                    read >> name;
                    fMaterialDataVec.fMatID.Push(id);
                    chunk.first = id;
                    chunk.second = name;
                    fMaterialDataVec.fMaterial.Push(chunk);
                    
                    if (max_dimension < dimension) {
                        max_dimension = dimension;
                    }
                }
                gmesh->SetDimension(max_dimension);
                
                char buf_end[1024];
                read.getline(buf_end, 1024);
                read.getline(buf_end, 1024);
                std::string str_end(buf_end);
                if(str_end == "$EndPhysicalNames" || str_end == "$EndPhysicalNames\r")
                {
                    std::cout << "Read mesh physical entities = " << n_entities << std::endl;
                }
                continue;
            }
            
            if(str == "$Nodes" || str == "$Nodes\r")
            {
                
                int64_t n_nodes;
                read >> n_nodes;
                
                int64_t node_id;
                double nodecoordX , nodecoordY , nodecoordZ ;
                gmesh -> NodeVec().Resize(n_nodes);
                gmesh->SetMaxNodeId(n_nodes-1);
                
                // needed for node insertion
                const int64_t Tnodes = n_nodes;
                TPZVec <TPZGeoNode> Node(Tnodes);
                
                for (int64_t inode = 0; inode < n_nodes; inode++) {

                    read.getline(buf, 1024);
                    read >> node_id;
                    read >> nodecoordX;
                    read >> nodecoordY;
                    read >> nodecoordZ;
                
                    Node[node_id-1].SetNodeId(node_id-1);
                    Node[node_id-1].SetCoord(0,nodecoordX/fDimensionlessL);
                    Node[node_id-1].SetCoord(1,nodecoordY/fDimensionlessL);
                    Node[node_id-1].SetCoord(2,nodecoordZ/fDimensionlessL);
                    gmesh->NodeVec()[node_id-1] = Node[node_id-1];
                    
                }

                
                char buf_end[1024];
                read.getline(buf_end, 1024);
                read.getline(buf_end, 1024);
                std::string str_end(buf_end);
                if(str_end == "$EndNodes" || str_end == "$EndNodes\r")
                {
                    std::cout << "Read mesh nodes = " <<  gmesh->NNodes() << std::endl;
                }
                continue;
            }
            
            if(str == "$Elements" || str == "$Elements\r")
            {
                
                int64_t n_elements;
                read >> n_elements;
                gmesh->SetMaxElementId(n_elements-1);
                
                for (int64_t iel = 0; iel < n_elements; iel++) {
                    this->InsertElement(gmesh, read);
                }

                char buf_end[1024];
                read.getline(buf_end, 1024);
                read.getline(buf_end, 1024);
                std::string str_end(buf_end);
                if(str_end == "$EndElements" || str_end == "$EndElements\r")
                {
                    std::cout << "Read mesh elements = " << gmesh->NElements() << std::endl;
                }
                continue;
            }
            
        }
        
    }

    std::cout << "Read General Mesh Data -> done!" << std::endl;
    gmesh->BuildConnectivity();
    std::cout << "Geometric Mesh Connectivity -> done!" << std::endl;
    return gmesh;
    
}// End Method

void TRMGmshReader::SetfDimensionlessL(REAL dimensionlessL)
{
    fDimensionlessL = dimensionlessL;
}

/** @brief Insert elements following msh file format */
bool TRMGmshReader::InsertElement(TPZGeoMesh * gmesh, std::ifstream & line){
    
    // first implementation based on linear elements: http://gmsh.info/doc/texinfo/gmsh.html#File-formats
    TPZManVector <int64_t,1> TopolPoint(1);
    TPZManVector <int64_t,2> TopolLine(2);
    TPZManVector <int64_t,3> TopolTriangle(3);
    TPZManVector <int64_t,4> TopolQuad(4);
    TPZManVector <int64_t,4> TopolTet(4);
    TPZManVector <int64_t,5> TopolPyr(5);
    TPZManVector <int64_t,8> TopolHex(8);
    
    TPZManVector <int64_t,2> TopolLineQ(3);
    TPZManVector <int64_t,6> TopolTriangleQ(6);
    TPZManVector <int64_t,8> TopolQuadQ(8);
    TPZManVector <int64_t,10> TopolTetQ(10);
    TPZManVector <int64_t,5> TopolPyrQ(14);
    TPZManVector <int64_t,8> TopolHexQ(20);
    

    int64_t element_id, type_id, div_id, physical_id, elementary_id;
    
    char buf[1024];
    line.getline(buf, 1024);
    line >> element_id;
    line >> type_id;
    line >> div_id;
    line >> physical_id;
    line >> elementary_id;
    
    switch (type_id) {
        case 1:
        {
            // Triangle
            line >> TopolLine[0]; //node 1
            line >> TopolLine[1]; //node 2
            element_id--;
            TopolLine[0]--;
            TopolLine[1]--;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, TopolLine, physical_id, *gmesh);
        }
            break;
        case 2:
        {
            // Triangle
            line >> TopolTriangle[0]; //node 1
            line >> TopolTriangle[1]; //node 2
            line >> TopolTriangle[2]; //node 3
            element_id--;
            TopolTriangle[0]--;
            TopolTriangle[1]--;
            TopolTriangle[2]--;
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (element_id, TopolTriangle, physical_id, *gmesh);
           
        }
            break;
        case 3:
        {
            // Quadrilateral
            line >> TopolQuad[0]; //node 1
            line >> TopolQuad[1]; //node 2
            line >> TopolQuad[2]; //node 3
            line >> TopolQuad[3]; //node 4
            element_id--;
            TopolQuad[0]--;
            TopolQuad[1]--;
            TopolQuad[2]--;
            TopolQuad[3]--;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (element_id, TopolQuad, physical_id, *gmesh);
            
        }
            break;
        case 4:
        {
            // Tetrahedron
            line >> TopolTet[0]; //node 1
            line >> TopolTet[1]; //node 2
            line >> TopolTet[2]; //node 3
            line >> TopolTet[3]; //node 4
            element_id--;
            TopolTet[0]--;
            TopolTet[1]--;
            TopolTet[2]--;
            TopolTet[3]--;
            new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (element_id, TopolTet, physical_id, *gmesh);
            
        }
            break;
        case 5:
        {
            // Hexahedra
            line >> TopolHex[0]; //node 1
            line >> TopolHex[1]; //node 2
            line >> TopolHex[2]; //node 3
            line >> TopolHex[3]; //node 4
            line >> TopolHex[4]; //node 5
            line >> TopolHex[5]; //node 6
            line >> TopolHex[6]; //node 7
            line >> TopolHex[7]; //node 8
            element_id--;
            TopolHex[0]--;
            TopolHex[1]--;
            TopolHex[2]--;
            TopolHex[3]--;
            TopolHex[4]--;
            TopolHex[5]--;
            TopolHex[6]--;
            TopolHex[7]--;
            new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (element_id, TopolHex, physical_id, *gmesh);
        }
            break;
        case 8:
        {
            // Triangle
            line >> TopolLineQ[0]; //node 1
            line >> TopolLineQ[1]; //node 2
            line >> TopolLineQ[2]; //node 2
            element_id--;
            TopolLineQ[0]--;
            TopolLineQ[1]--;
            TopolLineQ[2]--;
            new TPZGeoElRefPattern< pzgeom::TPZQuadraticLine> (element_id, TopolLineQ, physical_id, *gmesh);
        }
            break;
        case 9:
        {
            // Triangle
            line >> TopolTriangleQ[0]; //node 1
            line >> TopolTriangleQ[1]; //node 2
            line >> TopolTriangleQ[2]; //node 3
            line >> TopolTriangleQ[3]; //node 4
            line >> TopolTriangleQ[4]; //node 5
            line >> TopolTriangleQ[5]; //node 6
            element_id--;
            TopolTriangleQ[0]--;
            TopolTriangleQ[1]--;
            TopolTriangleQ[2]--;
            TopolTriangleQ[3]--;
            TopolTriangleQ[4]--;
            TopolTriangleQ[5]--;
            new TPZGeoElRefPattern< pzgeom::TPZQuadraticTrig> (element_id, TopolTriangleQ, physical_id, *gmesh);
        }
            break;
        case 10:
        {
            // Quadrilateral
            line >> TopolQuadQ[0]; //node 1
            line >> TopolQuadQ[1]; //node 2
            line >> TopolQuadQ[2]; //node 3
            line >> TopolQuadQ[3]; //node 4
            line >> TopolQuadQ[4]; //node 5
            line >> TopolQuadQ[5]; //node 6
            line >> TopolQuadQ[6]; //node 7
            line >> TopolQuadQ[7]; //node 8
            element_id--;
            TopolQuadQ[0]--;
            TopolQuadQ[1]--;
            TopolQuadQ[2]--;
            TopolQuadQ[3]--;
            TopolQuadQ[4]--;
            TopolQuadQ[5]--;
            TopolQuadQ[6]--;
            TopolQuadQ[7]--;
            new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (element_id, TopolQuadQ, physical_id, *gmesh);
        }
            break;
        case 111:
        {
            // Tetrahedron
            line >> TopolTetQ[0]; //node 1
            line >> TopolTetQ[1]; //node 2
            line >> TopolTetQ[2]; //node 3
            line >> TopolTetQ[3]; //node 4
            
            line >> TopolTetQ[4]; //node 5
            line >> TopolTetQ[5]; //node 6
            line >> TopolTetQ[6]; //node 7
            line >> TopolTetQ[7]; //node 8
            line >> TopolTetQ[8]; //node 9
            line >> TopolTetQ[9]; //node 10
            
            element_id--;
            TopolTetQ[0]--;
            TopolTetQ[1]--;
            TopolTetQ[2]--;
            TopolTetQ[3]--;
            
            TopolTetQ[4]--;
            TopolTetQ[5]--;
            TopolTetQ[6]--;
            TopolTetQ[7]--;
            TopolTetQ[8]--;
            TopolTetQ[9]--;
            DebugStop();            
            new TPZGeoElRefPattern< pzgeom::TPZQuadraticTetra> (element_id, TopolTetQ, physical_id, *gmesh);
            
        }
            break;
        case 122:
        {
            // Hexahedra
            line >> TopolHexQ[0]; //node 1
            line >> TopolHexQ[1]; //node 2
            line >> TopolHexQ[2]; //node 3
            line >> TopolHexQ[3]; //node 4
            line >> TopolHexQ[4]; //node 5
            line >> TopolHexQ[5]; //node 6
            line >> TopolHexQ[6]; //node 7
            line >> TopolHexQ[7]; //node 8
            
            line >> TopolHexQ[8]; //node 9
            line >> TopolHexQ[9]; //node 10
            line >> TopolHexQ[10]; //node 11
            line >> TopolHexQ[11]; //node 12
            line >> TopolHexQ[12]; //node 13
            line >> TopolHexQ[13]; //node 14
            line >> TopolHexQ[14]; //node 15
            line >> TopolHexQ[15]; //node 16
            line >> TopolHexQ[16]; //node 17
            line >> TopolHexQ[17]; //node 18
            line >> TopolHexQ[18]; //node 19
            line >> TopolHexQ[19]; //node 20
            
            element_id--;
            TopolHexQ[0]--;
            TopolHexQ[1]--;
            TopolHexQ[2]--;
            TopolHexQ[3]--;
            TopolHexQ[4]--;
            TopolHexQ[5]--;
            TopolHexQ[6]--;
            TopolHexQ[7]--;
            
            TopolHexQ[8]--;
            TopolHexQ[9]--;
            TopolHexQ[10]--;
            TopolHexQ[11]--;
            TopolHexQ[12]--;
            TopolHexQ[13]--;
            TopolHexQ[14]--;
            TopolHexQ[15]--;
            TopolHexQ[16]--;
            TopolHexQ[17]--;
            TopolHexQ[18]--;
            TopolHexQ[19]--;
            DebugStop();
            new TPZGeoElRefPattern< pzgeom::TPZQuadraticCube> (element_id, TopolHexQ, physical_id, *gmesh);
        }
            break;
        default:
        {
            std::cout << "Element not impelemented." << std::endl;
            DebugStop();
        }
            break;
    }
   
    return true;
}