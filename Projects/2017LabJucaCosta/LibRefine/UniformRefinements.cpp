//
//  SimilarUniformRefinements.cpp
//  PZ
//
//  Created by labmec on 08/10/17.
//
//

#include <stdio.h>
#include "CreateAndRefineMeshes.h"


/* 2. Functions to uniform refinement of the geometric meshes.
 Projects:
 Poisson3D_Shock
 */
// The function refines all geometric elements (without subelements) for any material or for one matidtodivided material - Only for one material
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh) {
    if(nDiv < 1 || !gmesh || gmesh->Dimension()<0) {
        std::cout << "It is nothing done! (UniformRefinement)." << std::endl;
        return;
    }
    
    TPZManVector<TPZGeoEl*> filhos;
    // Loop over (number of refinements) nDiv - It refines all geometric elements without subelements
    for(int D=0; D<nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            if(!gel || gel->HasSubElement() || !gel->Dimension())
                continue;
            else
                gel->Divide(filhos);
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

// The function refines all geometric elements (without subelements) for any material or for material ids in MatIdsVec. If dim < 0 will be refining all elements of the mesh (all dimension not zero)
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, TPZVec<int> *MatIdsVec) {
    if(!nDiv || !gmesh || !dim || dim>gmesh->Dimension()) {
        std::cout << "It is nothing done! (UniformRefinement)." << std::endl;
        return;
    }
    // If dim is negative, the user want to refine geometrical elements not points
    if(dim<0) {
        UniformRefinement(nDiv,gmesh);
        return;
    }
    
    // If MatIdsVector is null, the refinement will be made for all materials
    bool allmaterial = false;
    if(!MatIdsVec || !MatIdsVec->NElements()) {
        allmaterial = true;
    }
    
    int matidtodivided;
    TPZManVector<TPZGeoEl*> filhos;
    // Loop over (number of refinements) nDiv - It refines all geometric elements without subelements
    for(int D=0; D<nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            int geldim = gel->Dimension();
            if(!gel || gel->HasSubElement() || !geldim || geldim != dim)
                continue;
            if(!allmaterial) {
                int matidcurrent = gel->MaterialId();
                for(int count=0;count<MatIdsVec->NElements();count++) {
                    matidtodivided = (*MatIdsVec)[count];
                    if(matidcurrent == matidtodivided) {
                        gel->Divide(filhos);
                        break;
                    }
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

void RegularizeMesh(TPZGeoMesh *gmesh, int dimension)
{
    //Control flag
    bool changed = true;
    // If exists something wrong
    if(!gmesh || gmesh->Dimension() < 0)
        DebugStop();

    while (changed)
    {
        changed = false;
        int nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->ElementVec()[el];
            if (gel->HasSubElement()) {
                continue;
            }
            int dim = gel->Dimension();
            if (dim != dimension)
            {
                continue;
            }
            int nsides = gel->NSides();
            int nrefined = 0;
            int nsidedim = 0;
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                nsidedim++;
            }
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                TPZGeoElSide thisside(gel,is);
                TPZGeoElSide neighbour = thisside.Neighbour();
                if (neighbour != thisside) {
                    TPZStack<TPZGeoElSide> subelements;
                    neighbour.GetSubElements2(subelements);
                    int nsub = subelements.size();
                    if (nsub > 0) {
                        nrefined++;
                    }
                    for (int isub=0; isub<nsub; isub++) {
                        TPZGeoElSide sub = subelements[isub];
                        if (sub.Dimension() != dim-1) {
                            continue;
                        }
                        if (sub.HasSubElement()) {
                            TPZManVector<TPZGeoEl *> newsub;
                            gel->Divide(newsub);
                            changed = true;
                            break;
                        }
                    }
                }
                if (gel->HasSubElement()) {
                    break;
                }
            }
            if (nrefined >= nsidedim-1) {
                TPZManVector<TPZGeoEl *> newsub;
                gel->Divide(newsub);
                changed = true;
            }
        }
    }
	gmesh->CleanUp();
	gmesh->BuildConnectivity();
}
