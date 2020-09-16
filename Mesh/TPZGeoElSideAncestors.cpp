//
// Created by Gustavo A. Batistela on 15/09/2020.
//

#include "TPZGeoElSideAncestors.h"

void TPZGeoElSideAncestors::AddAncestor(TPZGeoElSide gelside) {
    
    if(!gelside) DebugStop();
    
    TPZGeoElSide lower = gelside.StrictFather();
    if(lower)
    {
        fAncestors.Push({gelside,lower});
    }
    TPZGeoElSide neighbour(gelside.Neighbour());
    while(neighbour != gelside)
    {
        lower = neighbour.StrictFather();
        if(lower)
        {
            fAncestors.Push({neighbour,lower});
            AddAncestor(lower);
            break;
        }
        neighbour = neighbour.Neighbour();
    }
}

void TPZGeoElSideAncestors::BuildAncestors() {

    fAncestors.resize(0);
    if (!fCurrent) {
        return;
    }
    AddAncestor(fCurrent);
}

TPZTransform<REAL> TPZGeoElSideAncestors::BuildTransform(TPZGeoElSide larger)
{
    if(!fCurrent) DebugStop();
    // check if larger is an ancestor of the current element
    int64_t num_ancestors = fAncestors.size();
    int ilevel = -1;
    for (int il = 0; il<num_ancestors; il++) {
        if(larger.IsNeighbour(fAncestors[il].second))
        {
            ilevel = il;
            break;
        }
    }
    if(ilevel == -1) DebugStop();
    // use the TPZGeoElSide method
    TPZTransform<REAL> result(fCurrent.Dimension());
    fCurrent.SideTransform3(larger, result);
    return result;
}

/// return true is a (strict) larger element with matid exists
bool TPZGeoElSideAncestors::HasLarger(int matid)
{
    int64_t num_ancestors = fAncestors.size();
    for (int il = 0; il<num_ancestors; il++) {
        TPZGeoElSide larger = fAncestors[il].second;
        if(larger.HasNeighbour(matid))
        {
            return true;
        }
    }
    return false;
}
