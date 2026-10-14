#include "TPZMHMApproxCreator.h"
#include "pzgmesh.h"


TPZMHMApproxCreator::TPZMHMApproxCreator() {
}

TPZMHMApproxCreator::TPZMHMApproxCreator(TPZVec<int64_t>& elementPartition) {
    fGeoToMHMDomain = elementPartition;
}

TPZMHMApproxCreator::~TPZMHMApproxCreator() {
}

void TPZMHMApproxCreator::SetMHMDomain(TPZGeoEl *gel, int domain) {
    if (!gel) {
        DebugStop();
    }
    const int64_t index = gel->Index();
    if (index >= fGeoToMHMDomain.size()) {
        fGeoToMHMDomain.Resize(index + 1, -1);
    }
    fGeoToMHMDomain[index] = domain;
}

/// Compute a map from the subdomain index to the set of geometric element indices that belong to the subdomain.
void TPZMHMApproxCreator::ComputeDomainMap() {
    fDomainMap.clear();
    for (int64_t i = 0; i < fGeoToMHMDomain.size(); i++) {
        int64_t domain = fGeoToMHMDomain[i];
        if (domain < 0) {
            continue;
        }
        fDomainMap[domain].insert(i);
    }
}

void TPZMHMApproxCreator::SetSubelementDomain(TPZGeoMesh *gmesh) {
    if (!gmesh) {
        DebugStop();
    }
    const int64_t nel = gmesh->NElements();
    if (fGeoToMHMDomain.size() < nel) {
        fGeoToMHMDomain.Resize(nel, -1);
    }
    for (int64_t el = 0; el < nel; ++el) {
        const int64_t domain = fGeoToMHMDomain[el];
        if (domain == -1) {
            continue;
        }
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || !gel->HasSubElement()) {
            continue;
        }
        for (int isub = 0; isub < gel->NSubElements(); ++isub) {
            TPZGeoEl *subel = gel->SubElement(isub);
            if (subel) {
                SetMHMDomain(subel, domain);
            }
        }
    }
}
