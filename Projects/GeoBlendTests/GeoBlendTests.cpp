
#include <pzgeoelrefless.h>
#include <tpzarc3d.h>
#include <TPZVTKGeoMesh.h>
#include <pzanalysis.h>
#include "tpzgeoblend.h"
#include "pzgeoquad.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

void CreateGeoMesh(TPZGeoMesh * &gmesh, bool printGMesh, int nDiv);

/**
 * Routine for testing the TPZGeoBlend refactoring
 * @return
 */
int main()
{

    bool printGMesh = true;
    int nDiv = 4;


    TPZGeoMesh * gmesh = nullptr;

    pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>::fUseNewX = false;
    CreateGeoMesh(gmesh, printGMesh, nDiv);
    delete gmesh;

//    pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>::fUseNewX = true;
//    CreateGeoMesh(gmesh, printGMesh, nDiv);
//    delete gmesh;


}

void CreateGeoMesh(TPZGeoMesh * &gmesh, bool printGMesh, int nDiv) {
    gmesh = new TPZGeoMesh();

    int64_t nnodes = 5;
    gmesh->NodeVec().Resize(nnodes);

    TPZVec<REAL> coord(3, 0.);
    for (int64_t i = 0; i < 4; i++) {
        coord[0] = ((i + 1) / 2) % 2;
        coord[1] = i / 2;
        gmesh->NodeVec()[i].SetCoord(coord);
        gmesh->NodeVec()[i].SetNodeId(i);
    }
    coord[0] = 0.5 + M_SQRT2 / 2;
    coord[1] = 0.5;

    gmesh->NodeVec()[4].SetCoord(coord);
    gmesh->NodeVec()[4].SetNodeId(4);

    int matIdQuad = 1, matIdArc = 2;
    TPZVec<int64_t> nodesIdVec(4);
    int64_t elId = 0;

    {
        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1;
        nodesIdVec[2] = 2;
        nodesIdVec[3] = 3;

        TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > *quadEl =
                new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(nodesIdVec, matIdQuad,
                                                                                 *gmesh);

    }
    nodesIdVec.Resize(3);
    {
        nodesIdVec[0] = 1;
        nodesIdVec[1] = 2;
        nodesIdVec[2] = 4;
        TPZGeoElRefPattern<pzgeom::TPZArc3D> *arc =
                new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    }
    gmesh->BuildConnectivity();


    {
        TPZVec<TPZGeoEl *> sons;
        for (int iDiv = 0; iDiv < nDiv; iDiv++) {
            int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *geo = gmesh->Element(iel);
                TPZGeoEl *current = geo;
                int nFather = 0;
                while (current->Father()) {
                    current = current->Father();
                    nFather++;
                }
                if (nFather <= iDiv && !geo->HasSubElement()) {
                    geo->Divide(sons);
                    nel = gmesh->NElements();
                }
            }
        }
    }

    if (printGMesh) {
        std::string meshFileName = "gmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        gmesh->Print(outTXT);
        outTXT.close();
        outVTK.close();
    }
}
