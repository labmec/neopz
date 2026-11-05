#ifndef PZ_TPZMHMCREATOR_H
#define PZ_TPZMHMCREATOR_H

#include "pzvec.h"
class TPZGeoMesh;
#include <map>
#include <set>

class TPZGeoEl;

class TPZMHMApproxCreator
{
protected:
  /// @brief MHM domain associated with each geometric element index.
  TPZVec<int64_t> fGeoToMHMDomain;

  /// @brief Vector of integers specifying the indices of the geometric elements that will be used as Lagrange multipliers.
  TPZVec<int64_t> fLagrangeElements;

  /// @brief Data structure identifying the geometric elements associated with each domain. The key is the domain index, and the value is a set of geometric element indices that belong to the domain.
  std::map<int64_t, std::set<int64_t>> fDomainMap;

public:
  TPZMHMApproxCreator();
  TPZMHMApproxCreator(TPZVec<int64_t> &elementPartition);

  ~TPZMHMApproxCreator();

  /// Set the MHM domain of a geometric element.
  void SetMHMDomain(TPZGeoEl *gel, int domain);

  /// Compute a map from the subdomain index to the set of geometric element indices that belong to the subdomain.
  void ComputeDomainMap();

  /// Set the MHM domain index of all subelements of every geometric element in a mesh.
  void SetSubelementDomain(TPZGeoMesh *gmesh);

  /// @brief Get the MHM domain associated with a geometric element index.
  const std::map<int64_t, std::set<int64_t>> &GetMHMDomainMap() const
  {
    return fDomainMap;
  }
};

#endif
