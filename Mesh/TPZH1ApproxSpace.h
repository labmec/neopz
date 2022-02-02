#ifndef _TPZH1APPROXSPACE_H_
#define _TPZH1APPROXSPACE_H_

#include "TPZApproxSpace.h"

/// Enum stating which family of H1 spaces will be created
enum class H1Family {
  Standard,//< Standard
  EnrichedPrism//<Prism is enriched such that it is De Rham compatible with HDiv els
};

class TPZH1ApproxSpace : public TPZApproxSpace{
private:
  H1Family fFamily{H1Family::Standard};
public:
  //! Default constructor
  TPZH1ApproxSpace() = default;
  //! Copy constructor
  TPZH1ApproxSpace(const TPZH1ApproxSpace &copy) = default;
  //! Copy assignment operator
  TPZH1ApproxSpace &operator=(const TPZH1ApproxSpace &copy) = default;
  //! Move constructor
  TPZH1ApproxSpace(TPZH1ApproxSpace &&copy) = default;
  //! Move assignment operator
  TPZH1ApproxSpace &operator=(TPZH1ApproxSpace &&copy) = default;
  //! Destructor
  ~TPZH1ApproxSpace() = default;
  //! Gets which H1 family will be built
  H1Family GetFamily() const { return fFamily;}
  //! Sets which H1 family will be built
  void SetFamily(H1Family fam) { fFamily = fam;}

  //the following methods override methods from TPZApproxSpace
  
  //! Class identifier
  int ClassId() const override;
  //! Read method
  void Read(TPZStream &buf, void *context) override;
  //! Write method
  void Write(TPZStream &buf, int withclassid) const override;
  //! Creates a computational element using the function pointer for the topology
  TPZCompEl *CreateCompEl(TPZGeoEl *gel, TPZCompMesh &mesh) const override;
  //@}
};

#endif /* _TPZH1APPROXSPACE_H_ */
