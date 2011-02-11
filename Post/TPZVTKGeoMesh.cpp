/*
 *  TPZVTKGeoMesh.cpp
 *  Crack
 *
 *  Created by Cesar Lucci on 16/08/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZVTKGeoMesh.h"

#include <sstream>

/**
 * Generate an output of all geomesh to VTK
 */
void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, bool matColor)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	std::stringstream node, connectivity, type, material;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{				
		if(gmesh->ElementVec()[el]->Type() == EOned && !gmesh->ElementVec()[el]->IsLinearMapping())//Exclude Arc3D and Ellipse3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->HasSubElement())
		{
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}			
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType = TPZVTKGeoMesh::GetVTK_ElType(gmesh->ElementVec()[el]);
		type << elType << std::endl;
		
		if(matColor == true)
		{
			material << gmesh->ElementVec()[el]->MaterialId() << std::endl;
		}
		else
		{
			material << elType << std::endl;
		}
		
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str() << std::endl;
	
	file << "CELL_DATA" << " " << nVALIDelements << std::endl;
	file << "FIELD FieldData 1" << std::endl;
	if(matColor == true)
	{
		file << "material 1 " << nVALIDelements << " int" << std::endl;
	}
	else
	{
		file << "ElementType 1 " << nVALIDelements << " int" << std::endl;
	}
	file << material.str();
	
	file.close();
}

/**
 * Generate an output of all geomesh to VTK, associating to each one the given data
 */
void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, TPZVec<int> &elData)
{
	if(gmesh->NElements() != elData.NElements())
	{
		std::cout << "Wrong vector size of elements data!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
	}
	file.clear();
	int nelements = gmesh->NElements();
	
	std::stringstream node, connectivity, type, material;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{				
		if(gmesh->ElementVec()[el]->Type() == EOned && !gmesh->ElementVec()[el]->IsLinearMapping())//Exclude Arc3D and Ellipse3D
		{
			continue;
		}
		//			if(gmesh->ElementVec()[el]->HasSubElement())
		if (elData[el] == -999) {
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}			
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType = TPZVTKGeoMesh::GetVTK_ElType(gmesh->ElementVec()[el]);
		type << elType << std::endl;
		
		material << elData[el] << std::endl;
		
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str() << std::endl;
	
	file << "CELL_DATA" << " " << nVALIDelements << std::endl;
	file << "FIELD FieldData 1" << std::endl;
	
	file << "Substructure 1 " << nVALIDelements << " int" << std::endl;
	
	file << material.str();
	
	file.close();
}

/**
 * Based on a given geomesh, just the elements that have an neighbour with a given material id will be exported to an VTK file
 */
void TPZVTKGeoMesh::PrintGMeshVTKneighbour_material(TPZGeoMesh * gmesh, std::ofstream &file, int neighMaterial, bool matColor)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	std::stringstream node, connectivity, type, material;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{				
		if(gmesh->ElementVec()[el]->Type() == EOned && !gmesh->ElementVec()[el]->IsLinearMapping())//Exclude Arc3D and Ellipse3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->HasSubElement())
		{
			continue;
		}
		
		bool matFound = false;
		for(int s = 0; s < gmesh->ElementVec()[el]->NSides(); s++)
		{
			TPZGeoElSide thisSide(gmesh->ElementVec()[el], s);
			TPZGeoElSide neighSide = thisSide.Neighbour();
			
			while(thisSide != neighSide)
			{
				if(neighSide.Element()->MaterialId() == neighMaterial)
				{
					matFound = true;
					break;
				}
				neighSide = neighSide.Neighbour();
			}
			if(matFound)
			{
				break;
			}
		}
		if(!matFound)
		{
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}			
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType = TPZVTKGeoMesh::GetVTK_ElType(gmesh->ElementVec()[el]);
		type << elType << std::endl;
		
		if(matColor == true)
		{
			material << gmesh->ElementVec()[el]->MaterialId() << std::endl;
		}
		else
		{
			material << elType << std::endl;
		}
		
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();
	
	file << "CELL_DATA" << " " << nVALIDelements << std::endl;
	file << "FIELD FieldData 1" << std::endl;
	if(matColor == true)
	{
		file << "material 1 " << nVALIDelements << " int" << std::endl;
	}
	else
	{
		file << "ElementType 1 " << nVALIDelements << " int" << std::endl;
	}
	file << material.str();
	
	file.close();
}

/**
 * Based on a given geomesh, just the elements that have the given material id will be exported to an VTK file
 */
void TPZVTKGeoMesh::PrintGMeshVTKmy_material(TPZGeoMesh * gmesh, std::ofstream &file, std::set<int> myMaterial, bool matColor)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	std::stringstream node, connectivity, type, material;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{				
		if(gmesh->ElementVec()[el]->Type() == EOned && !gmesh->ElementVec()[el]->IsLinearMapping())//Exclude Arc3D and Ellipse3D
		{
			continue;
		}
		int mat = gmesh->ElementVec()[el]->MaterialId();
		bool found = !(myMaterial.find(mat) == myMaterial.end() );
		if(gmesh->ElementVec()[el]->HasSubElement() || !found)
		{
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}			
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType = TPZVTKGeoMesh::GetVTK_ElType(gmesh->ElementVec()[el]);
		type << elType << std::endl;
		
		if(matColor == true)
		{
			material << gmesh->ElementVec()[el]->MaterialId() << std::endl;
		}
		else
		{
			material << elType << std::endl;
		}
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();
	
	file << "CELL_DATA" << " " << nVALIDelements << std::endl;
	file << "FIELD FieldData 1" << std::endl;
	if(matColor == true)
	{
		file << "material 1 " << nVALIDelements << " int" << std::endl;
	}
	else
	{
		file << "ElementType 1 " << nVALIDelements << " int" << std::endl;
	}
	file << material.str();
	
	file.close();
}

int TPZVTKGeoMesh::GetVTK_ElType(TPZGeoEl * gel)
{
	MElementType pzElType = gel->Type();
	
	int elType = -1;
	switch (pzElType)
	{
		case(EPoint):
		{
			elType = 1;
			break;
		}
		case(EOned):
		{
			if(gel->NNodes() == 2)
			{
				elType = 3;	
			}
			break;
		}
		case (ETriangle):
		{
			if(gel->NNodes() == 3)
			{
				elType = 5;	
			}
			
			break;				
		}
		case (EQuadrilateral):
		{
			if(gel->NNodes() == 4)
			{
				elType = 9;					
			}
			
			break;				
		}
		case (ETetraedro):
		{
			if(gel->NNodes() == 8)
			{
				elType = 10;	
			}
			break;				
		}
		case (EPiramide):
		{
			if(gel->NNodes() == 5)
			{
				elType = 14;	
			}
			break;				
		}
		case (EPrisma):
		{
			if(gel->NNodes() == 6)
			{
				elType = 13;	
			}
			break;				
		}
		case (ECube):
		{
			if(gel->NNodes() == 8)
			{
				elType = 12;	
			}
			break;				
		}
		default:
		{
			std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
			DebugStop();
			break;	
		}
	}
	
	return elType;
}