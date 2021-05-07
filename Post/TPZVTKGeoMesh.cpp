/**
 * @file
 * @brief Contains the implementation of the TPZVTKGeoMesh methods. 
 */

#include "TPZVTKGeoMesh.h"
#include "pzgeoel.h"
#include "pzintel.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeopoint.h"
#include <sstream>
#include "pzcmesh.h"

/**
 * Generate an output of all geomesh to VTK
 */
void TPZVTKGeoMesh::PrintCMeshVTK(TPZCompMesh * cmesh, std::ofstream &file, bool matColor) {
    cmesh->LoadReferences();

    file.clear();
    int64_t nelements = cmesh->NElements();

    std::stringstream node, connectivity, type, material, index, referenceIndex;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, size = 0, nVALIDelements = 0;

    for (int64_t el = 0; el < nelements; el++) {
        TPZCompEl * cel = cmesh->ElementVec()[el];
        if (cel == NULL || cel->Reference() == NULL) {
            continue;
        }
        TPZGeoEl * gel = cel->Reference();
        if (!gel) continue;
        if (!gel->Reference()) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = cmesh->Reference()->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;

            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;

        if (matColor == true) {
            material << gel->MaterialId() << std::endl;
        } else {
            material << elType << std::endl;
        }
        index << cel->Index() << std::endl;
        referenceIndex << gel->Index() << std::endl;

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

    file << "CELL_DATA " << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData 3 " << std::endl;
    if (matColor == true) {
        file << "material 1 " << nVALIDelements << " int " << std::endl;
    } else {
        file << "ElementType 1 " << nVALIDelements << " int " << std::endl;
    }
    file << material.str();

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    file << "geoElIndex 1 " << nVALIDelements << " int" << std::endl;
    file << referenceIndex.str();

    file.close();
}

/*
* Generate an output of all geomesh to VTK
*/
void TPZVTKGeoMesh::PrintCMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, bool matColor) {
    
    file.clear();
    int64_t nelements = gmesh->NElements();
    
    std::stringstream node, connectivity, type, material, index, referenceIndex;
    
    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;
    
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";
    
    int64_t actualNode = -1, size = 0, nVALIDelements = 0;
    
    for (int64_t el = 0; el < nelements; el++) {
        TPZGeoEl * gel = gmesh->ElementVec()[el];
        if (gel == NULL || gel->Reference() == NULL) {
            continue;
        }
        TPZCompEl *cel = gel->Reference();
        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);
        
        size += (1 + elNnodes);
        connectivity << elNnodes;
        
        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = gmesh->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;
            
            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;
        
        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;
        
        if (matColor == true) {
            material << gel->MaterialId() << std::endl;
        } else {
            material << elType << std::endl;
        }
        index << cel->Index() << std::endl;
        referenceIndex << gel->Index() << std::endl;
        
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
    
    file << "CELL_DATA " << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData 3 " << std::endl;
    if (matColor == true) {
        file << "material 1 " << nVALIDelements << " int " << std::endl;
    } else {
        file << "ElementType 1 " << nVALIDelements << " int " << std::endl;
    }
    file << material.str();
    
    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();
    
    file << "geoElIndex 1 " << nVALIDelements << " int" << std::endl;
    file << referenceIndex.str();
    
    file.close();
}

/**
 * Generate an output of all geomesh to VTK
 */
void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, bool matColor, bool dimension) {
    file.clear();
    int64_t nelements = gmesh->NElements();
    TPZGeoEl *gel;

    std::stringstream node, connectivity, type, material, index, eldimension;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, size = 0, nVALIDelements = 0;

    for (int64_t el = 0; el < nelements; el++) {
        gel = gmesh->ElementVec()[el];
        if (!gel)//|| (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (gel->HasSubElement()) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = gmesh->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;

            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;

        if (matColor == true) {
            material << gel->MaterialId() << std::endl;
        }
        if (dimension == true) {
            eldimension << gel->Dimension() << std::endl;
        }
        index << gel->Index() << std::endl;
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
    file << "FIELD FieldData " << 1+int(matColor)+int(dimension) << std::endl;
    if (matColor == true) {
        file << "material 1 " << nVALIDelements << " int" << std::endl;
        file << material.str();
    }
    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();
    if (dimension == true) {
        file << "Dimension 1 " << nVALIDelements << " int" << std::endl;
        file << eldimension.str();
    }

    file.close();
}


// Generate an output of all geomesh to VTK, associating to each one the given data

void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, TPZVec<int> &elData) {
    if (gmesh->NElements() != elData.NElements()) {
        std::cout << "Wrong vector size of elements data!" << std::endl;
        std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
    }
    file.clear();
    int64_t nelements = gmesh->NElements();

    std::stringstream node, connectivity, type, material, eldat, eldim, index;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, size = 0, nVALIDelements = 0;
    TPZGeoEl *gel;

    for (int64_t el = 0; el < nelements; el++) {
        gel = gmesh->ElementVec()[el];
        if (!gel || (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (elData[el] == -999) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = gmesh->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;

            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;

        eldat  << elData[el] << std::endl;
        material << gel->MaterialId() << std::endl;
        eldim << gel->Dimension() << std::endl;

        index << gel->Index() << std::endl;

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
    file << "FIELD FieldData 4" << std::endl;

    file << "ElData 1 " << nVALIDelements << " int" << std::endl;
    file << eldat.str();

    file << "Material 1 " << nVALIDelements << " int" << std::endl;
    file << material.str();

    file << "Dimension 1 " << nVALIDelements << " int" << std::endl;
    file << eldim.str();
    
    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    file.close();
}


// Generate an output of all geomesh to VTK, associating to each one the given data

void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, TPZVec<REAL> &elData) {
    if (gmesh->NElements() != elData.NElements()) {
        std::cout << "Wrong vector size of elements data!" << std::endl;
        std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
    }
    file.clear();
    int64_t nelements = gmesh->NElements();

    std::stringstream node, connectivity, type, material, eldat, eldim, index;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, size = 0, nVALIDelements = 0;
    TPZGeoEl *gel;

    for (int64_t el = 0; el < nelements; el++) {
        gel = gmesh->ElementVec()[el];
        if (!gel || (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (elData[el] == -999.) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = gmesh->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;

            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;

        eldat  << elData[el] << std::endl;
        material << gel->MaterialId() << std::endl;
        eldim << gel->Dimension() << std::endl;


        index << gel->Index() << std::endl;

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
    file << "FIELD FieldData 4" << std::endl;

    file << "ElData 1 " << nVALIDelements << " float" << std::endl;
    file << eldat.str();

    file << "Material 1 " << nVALIDelements << " int" << std::endl;
    file << material.str();

    file << "Dimension 1 " << nVALIDelements << " int" << std::endl;
    file << eldim.str();

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    file.close();
}

/** @brief Generate an output of all geometric elements that have a computational counterpart to VTK */
void TPZVTKGeoMesh::PrintCMeshVTK(TPZCompMesh * cmesh, std::ofstream &file, TPZVec<REAL> &elData, std::string dataName) {
    cmesh->LoadReferences();

    file.clear();
    int64_t nelements = cmesh->NElements();

    std::stringstream node, connectivity, type, material, index, referenceIndex;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, size = 0, nVALIDelements = 0;

    for (int64_t el = 0; el < nelements; el++) {
        TPZCompEl * cel = cmesh->ElementVec()[el];
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl * gel = cel->Reference();
        if (!gel) continue;
        if (!gel || (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (elData[el] == -999.) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = cmesh->Reference()->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;

            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;

        material << elData[el] << std::endl;

        index << cel->Index() << std::endl;
        referenceIndex << gel->Index() << std::endl;

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
    file << "FIELD FieldData 3" << std::endl;

    file << dataName << " 1 " << nVALIDelements << " float" << std::endl;

    file << material.str();

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    file << "geoElIndex 1 " << nVALIDelements << " int" << std::endl;
    file << referenceIndex.str();

    file.close();
}

// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given

void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, char *filename, TPZChunkVector<int> &elData) {
    std::ofstream file(filename);
#ifdef PZDEBUG
    if (!file.is_open())
        DebugStop();
#endif

    if (gmesh->NElements() != elData.NElements()) {
        std::cout << "Wrong vector size of elements data!" << std::endl;
        std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
    }
    file.clear();
    int64_t nelements = gmesh->NElements();

    std::stringstream connectivity, type, material, index;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t t, c, el;
    int64_t actualNode = -1, size = 0, nVALIDelements = 0;
    int64_t counternodes = gmesh->NNodes();
    TPZGeoEl *gel;

    for (el = 0; el < nelements; el++) {
        gel = gmesh->ElementVec()[el];
        if (!gel || (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (elData[el] == -999) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (t = 0; t < elNnodes; t++) {
            actualNode = gel->NodeIndex(t);
            if (actualNode < 0)
                DebugStop();

            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;

        material << elData[el] << std::endl;

        index << gel->Index() << std::endl;

        nVALIDelements++;
    }

    // Printing all nodes of the mesh
    file << counternodes << " float" << std::endl;
    for (t = 0; t < counternodes; t++) {
        TPZGeoNode *node = &(gmesh->NodeVec()[t]);
        for (c = 0; c < 3; c++) {
            double coord = node->Coord(c);
            file << coord << " ";
        }
        file << std::endl;
    }

    file << std::endl << "CELLS " << nVALIDelements << " ";

    file << size << std::endl;
    file << connectivity.str() << std::endl;

    file << "CELL_TYPES " << nVALIDelements << std::endl;
    file << type.str() << std::endl;

    file << "CELL_DATA" << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData 2" << std::endl;

    file << "Substructure 1 " << nVALIDelements << " int" << std::endl;

    file << material.str();

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    file.close();
}

// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given

void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, char *filename, TPZVec<REAL> &elData) {
    std::ofstream file(filename);
#ifdef PZDEBUG
    if (!file.is_open())
        DebugStop();
#endif

    if (gmesh->NElements() != elData.NElements()) {
        std::cout << "Wrong vector size of elements data!" << std::endl;
        std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
    }
    file.clear();
    int64_t nelements = gmesh->NElements();

    std::stringstream connectivity, type, datael, index;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t t, c, el;
    int64_t actualNode = -1, size = 0, nVALIDelements = 0;
    int64_t counternodes = gmesh->NNodes();
    TPZGeoEl *gel;

    for (el = 0; el < nelements; el++) {
        gel = gmesh->ElementVec()[el];
        if (!gel || (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (elData[el] == -999) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (t = 0; t < elNnodes; t++) {
            actualNode = gel->NodeIndex(t);
            if (actualNode < 0)
                DebugStop();

            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;

        datael << elData[el] << std::endl;
        index << gel->Index() << std::endl;

        nVALIDelements++;
    }

    // Printing all nodes of the mesh
    file << counternodes << " float" << std::endl;
    for (t = 0; t < counternodes; t++) {
        TPZGeoNode *node = &(gmesh->NodeVec()[t]);
        for (c = 0; c < 3; c++) {
            double coord = node->Coord(c);
            file << coord << " ";
        }
        file << std::endl;
    }

    file << std::endl << "CELLS " << nVALIDelements << " ";

    file << size << std::endl;
    file << connectivity.str() << std::endl;

    file << "CELL_TYPES " << nVALIDelements << std::endl;
    file << type.str() << std::endl;

    file << "CELL_DATA" << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData 2" << std::endl;

    file << "Substructure 1 " << nVALIDelements << " float" << std::endl;
    file << datael.str();

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    file.close();
}
// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given

void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, char *filename, TPZVec<TPZVec<REAL> > &elData) {
    std::ofstream file(filename);
#ifdef PZDEBUG
    if (!file.is_open())
        DebugStop();
#endif

    if (gmesh->NElements() != elData.NElements()) {
        std::cout << "Wrong vector size of elements data!" << std::endl;
        std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
    }
    file.clear();
    int64_t nelements = gmesh->NElements();
    int64_t ndatas = elData[0].NElements();
    if (!ndatas) {
        file.close();
        return;
    }

    std::stringstream connectivity, type, index;
    std::stringstream *datael = new std::stringstream[ndatas];

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t dat, t, c, el;
    int64_t actualNode = -1, size = 0, nVALIDelements = 0;
    int64_t counternodes = gmesh->NNodes();
    TPZGeoEl *gel;

    for (el = 0; el < nelements; el++) {
        gel = gmesh->ElementVec()[el];
        if (!gel || (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (elData[el][0] < 0.) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (t = 0; t < elNnodes; t++) {
            actualNode = gel->NodeIndex(t);
            if (actualNode < 0)
                DebugStop();

            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;
        index << gel->Index() << std::endl;

        for (dat = 0; dat < ndatas; dat++)
            *(datael + dat) << elData[el][dat] << std::endl;

        nVALIDelements++;
    }

    // Printing all nodes of the mesh
    file << counternodes << " float" << std::endl;
    for (t = 0; t < counternodes; t++) {
        TPZGeoNode *node = &(gmesh->NodeVec()[t]);
        for (c = 0; c < 3; c++) {
            double coord = node->Coord(c);
            file << coord << " ";
        }
        file << std::endl;
    }

    file << std::endl << "CELLS " << nVALIDelements << " ";

    file << size << std::endl;
    file << connectivity.str() << std::endl;

    file << "CELL_TYPES " << nVALIDelements << std::endl;
    file << type.str() << std::endl;

    file << "CELL_DATA" << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData " << ndatas + 1 << std::endl;

    for (dat = 0; dat < ndatas; dat++) {
        file << "Substructure" << dat + 1 << " 1 " << nVALIDelements << " float" << std::endl;
        file << datael[dat].str();
        file << std::endl;
    }

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    if (datael)
        delete[] datael;
    file.close();
}

// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given

void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, const char *filename, int var) {
    std::ofstream file(filename);
#ifdef PZDEBUG
    if (!file.is_open())
        DebugStop();
#endif

    file.clear();
    int64_t nelements = gmesh->NElements();

    std::stringstream connectivity, type, material, index;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t t, c, el;
    int64_t actualNode = -1, size = 0, nVALIDelements = 0;
    int64_t counternodes = gmesh->NNodes();
    TPZGeoEl *gel;
    TPZCompEl *cel;
    TPZVec<REAL> qsi(3, 0.);
    TPZVec<STATE> sol;

    for (el = 0; el < nelements; el++) {
        gel = gmesh->ElementVec()[el];
        // Print only to elements not refines (computational elements actives?)
        //		if(!gel || gel->HasSubElement())
        //			continue;
        if (!gel) continue;
        cel = gel->Reference();
        //		if(!cel) continue;
        if (gel->HasSubElement()) {
            continue;
        }
        if (gel->Type() == EOned && !gel->IsLinearMapping())//Exclude Arc3D and Ellipse3D
            continue;

        int elNnodes = gel->NCornerNodes();
        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (t = 0; t < elNnodes; t++) {
            actualNode = gel->NodeIndex(t);
            if (actualNode < 0)
                DebugStop();

            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
        type << elType << std::endl;

        // calculando o valor da solucao para o elemento
        //		cel->Solution(qsi,var,sol);
        //		if(sol.NElements()) material << sol[0] << std::endl;
        //		else material << 0.0 << std::endl;
        material << gel->MaterialId() << std::endl;

        index << gel->Index() << std::endl;

        nVALIDelements++;
    }

    // Printing all nodes of the mesh
    file << counternodes << " float" << std::endl;
    for (t = 0; t < counternodes; t++) {
        TPZGeoNode *node = &(gmesh->NodeVec()[t]);
        for (c = 0; c < 3; c++) {
            double coord = node->Coord(c);
            file << coord << " ";
        }
        file << std::endl;
    }

    file << std::endl << "CELLS " << nVALIDelements << " ";

    file << size << std::endl;
    file << connectivity.str() << std::endl;

    file << "CELL_TYPES " << nVALIDelements << std::endl;
    file << type.str() << std::endl;

    file << "CELL_DATA" << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData 2" << std::endl;

    file << "Substructure 1 " << nVALIDelements << " double" << std::endl;

    file << material.str();

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    file.close();
}

/**
 * Based on a given geomesh, just the elements that have an neighbour with a given material id will be exported to an VTK file
 */
void TPZVTKGeoMesh::PrintGMeshVTKneighbour_material(TPZGeoMesh * gmesh, std::ofstream &file, int neighMaterial, bool matColor) {
    file.clear();
    int64_t nelements = gmesh->NElements();

    std::stringstream node, connectivity, type, material, index;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, size = 0, nVALIDelements = 0;

    for (int64_t el = 0; el < nelements; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        if (gel->Type() == EOned && !gel->IsLinearMapping())//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (gel->HasSubElement()) {
            continue;
        }

        bool matFound = false;
        for (int s = 0; s < gel->NSides(); s++) {
            TPZGeoElSide thisSide(gmesh->ElementVec()[el], s);
            TPZGeoElSide neighSide = thisSide.Neighbour();

            while (thisSide != neighSide) {
                if (neighSide.Element()->MaterialId() == neighMaterial) {
                    matFound = true;
                    break;
                }
                neighSide = neighSide.Neighbour();
            }
            if (matFound) {
                break;
            }
        }
        if (!matFound) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = gmesh->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;

            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gmesh->ElementVec()[el]);
        type << elType << std::endl;

        if (matColor == true) {
            material << gel->MaterialId() << std::endl;
        } else {
            material << elType << std::endl;
        }
        index << gel->Index() << std::endl;
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
    file << "FIELD FieldData 2" << std::endl;
    if (matColor == true) {
        file << "material 1 " << nVALIDelements << " int" << std::endl;
    } else {
        file << "ElementType 1 " << nVALIDelements << " int" << std::endl;
    }
    file << material.str();

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    file.close();
}

void TPZVTKGeoMesh::PrintGMeshVTKneighbourhood(TPZGeoMesh * gmesh, int64_t elIndex, std::ofstream &file) {
    int elMat = 999;
    int surrMat = 555;
    std::set<int> myMaterial;
    myMaterial.insert(elMat);
    myMaterial.insert(surrMat);

    TPZGeoMesh * gmeshCP = new TPZGeoMesh(*gmesh);

    TPZGeoEl * gel = gmeshCP->ElementVec()[elIndex];
    SetMaterialVTK(gel, elMat);

    int nsides = gel->NSides();
    for (int s = 0; s < nsides; s++) {
        TPZGeoElSide thisSide(gel, s);
        TPZGeoElSide neighSide = thisSide.Neighbour();

        while (thisSide != neighSide) {
            TPZGeoEl * neighEl = neighSide.Element();
            if (thisSide.IsAncestor(neighSide) == false) {
                SetMaterialVTK(neighEl, surrMat);
            }
            neighSide = neighSide.Neighbour();
        }

    }
    PrintGMeshVTKmy_material(gmeshCP, file, myMaterial, true);
}

void TPZVTKGeoMesh::PrintGMeshVTK(TPZGeoMesh * gmesh, std::set<int64_t> & elIndex, std::ofstream &file) {
    file.clear();
    int64_t nelements = gmesh->NElements();

    std::stringstream node, connectivity, type, material;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, size = 0, nVALIDelements = 0;

    for (int64_t el = 0; el < nelements; el++) {
        if (elIndex.find(el) == elIndex.end()) {
            continue;
        }
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        if (gel->Type() == EOned && !gel->IsLinearMapping())//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        if (gel->HasSubElement()) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = gmesh->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;

            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gmesh->ElementVec()[el]);
        type << elType << std::endl;

        material << el << std::endl;

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

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;

    file << material.str();

    file.close();
}

void TPZVTKGeoMesh::SetMaterialVTK(TPZGeoEl * gel, int mat) {
    gel->SetMaterialId(mat);

    int64_t nnodes = gel->NNodes();
    for (int64_t nd = 0; nd < nnodes; nd++) {
        TPZVec<REAL> NodeCoord(3);
        TPZVec<int64_t> Topol(1);

        int64_t elIndex = 0;

        Topol[0] = gel->NodeIndex(nd);

        gel->Mesh()->NodeVec()[gel->NodeIndex(nd)].GetCoordinates(NodeCoord);
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elIndex, Topol, mat, *(gel->Mesh()));
    }

    if (gel->HasSubElement()) {
        int nSons = gel->NSubElements();
        for (int s = 0; s < nSons; s++) {
            TPZGeoEl * son = gel->SubElement(s);
            SetMaterialVTK(son, mat);
        }
    }
}

/**
 * Based on a given geomesh, just the elements that have the given material id will be exported to an VTK file
 */
void TPZVTKGeoMesh::PrintGMeshVTKmy_material(TPZGeoMesh * gmesh, std::ofstream &file, std::set<int> myMaterial, bool matColor, bool dimension) {
    file.clear();
    int64_t nelements = gmesh->NElements();

    std::stringstream node, connectivity, type, material, index, eldimension;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, size = 0, nVALIDelements = 0;

    for (int64_t el = 0; el < nelements; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (!gel->IsLinearMapping() && gel->Dimension() == 1)//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
        int mat = gel->MaterialId();
        bool found = !(myMaterial.find(mat) == myMaterial.end());
        if (gel->HasSubElement() || !found) {
            continue;
        }

        MElementType elt = gel->Type();
        int elNnodes = MElementType_NNodes(elt);

        size += (1 + elNnodes);
        connectivity << elNnodes;

        for (int t = 0; t < elNnodes; t++) {
            for (int c = 0; c < 3; c++) {
                double coord = gmesh->NodeVec()[gel->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;

            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;

        int elType = TPZVTKGeoMesh::GetVTK_ElType(gmesh->ElementVec()[el]);
        type << elType << std::endl;

        if (matColor == true) {
            material << gel->MaterialId() << std::endl;
        } else {
            material << elType << std::endl;
        }
        if (dimension == true) {
            eldimension << gel->Dimension() << std::endl;
        }
        index << gel->Index() << std::endl;
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
    file << "FIELD FieldData " << 2+int(dimension) << std::endl;
    if (matColor == true) {
        file << "material 1 " << nVALIDelements << " int" << std::endl;
    } else {
        file << "ElementType 1 " << nVALIDelements << " int" << std::endl;
    }
    file << material.str();

    file << "elIndex 1 " << nVALIDelements << " int" << std::endl;
    file << index.str();

    if (dimension == true) {
        file << "Dimension 1 " << nVALIDelements << " int" << std::endl;
        file << eldimension.str();
    }
    
    file.close();
}

int TPZVTKGeoMesh::GetVTK_ElType(TPZGeoEl * gel) {
    MElementType pzElType = gel->Type();

    int elType = -1;
    switch (pzElType) {
        case(EPoint):
        {
            elType = 1;
            break;
        }
        case(EOned):
        {
            elType = 3;
            break;
        }
        case (ETriangle):
        {
            elType = 5;
            break;
        }
        case (EQuadrilateral):
        {
            elType = 9;
            break;
        }
        case (ETetraedro):
        {
            elType = 10;
            break;
        }
        case (EPiramide):
        {
            elType = 14;
            break;
        }
        case (EPrisma):
        {
            elType = 13;
            break;
        }
        case (ECube):
        {
            elType = 12;
            break;
        }
        default:
        {
            std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
            break;
        }
    }
    if (elType == -1) {
        std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
        std::cout << "MIGHT BE CURVED ELEMENT (quadratic or quarter point)" << std::endl;
        DebugStop();
    }

    return elType;
}

/** Print a pointmesh whose values are the polynomial orders */
void TPZVTKGeoMesh::PrintPOrderPoints(TPZCompMesh &cmesh, std::set<int> dimensions, std::ofstream &file) {
    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";
    std::stringstream points;
    int64_t numpoints = 0;
    std::stringstream celldata;
    std::stringstream fielddata;
    int pointtype = 1;
    int64_t nel = cmesh.NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = cmesh.ElementVec()[iel];
        if (!cel) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        int nsides = gel->NSides();
        for (int side = 0; side < nsides; side++) {
            if (intel->NSideConnects(side) == 0) {
                continue;
            }
            int sidedim = gel->SideDimension(side);
            if (dimensions.find(sidedim) == dimensions.end()) {
                continue;
            }

            TPZManVector<REAL, 4> coord(3, 0.), xi(gel->Dimension(), 0.);
            gel->CenterPoint(side, xi);
            gel->X(xi, coord);
            points << coord << std::endl;
            celldata << "1 " << numpoints << std::endl;
            numpoints++;
            TPZConnect &c = intel->MidSideConnect(side);
            int corder = c.Order();
            fielddata << corder << std::endl;

        }
    }
    file << numpoints << " float" << std::endl << points.str();

    file << "CELLS " << numpoints << " ";

    file << 2 * numpoints << std::endl;
    file << celldata.str() << std::endl;

    file << "CELL_TYPES " << numpoints << std::endl;
    for (int64_t i = 0; i < numpoints; i++) {
        file << pointtype << std::endl;
    }

    file << "CELL_DATA" << " " << numpoints << std::endl;
    file << "FIELD FieldData 1" << std::endl;
    file << "porder 1 " << numpoints << " int" << std::endl;
    file << fielddata.str();
    file << std::endl;

    file.close();

}
