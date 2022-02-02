#ifndef _TPZAPPROXSPACE_H_
#define _TPZAPPROXSPACE_H_

#include "TPZSavable.h"
#include "pzvec.h"

#include <set>

class TPZCompEl;
class TPZCompMesh;
class TPZGeoEl;

/*
 * @brief Administer the creation of approximation spaces
 * @author Philippe Devloo
 * @since 2009
 * @ingroup interpolation
 */
class TPZApproxSpace : public TPZSavable {
private:
  //! Allows to create hybrid meshes
  bool fCreateHybridMesh{false};
    
  //! Indicates whether each element should have an aditional lagrange multiplier
  bool fCreateLagrangeMultiplier{false};
    
  //! Indicates that the elements need to be created with memory (at int point)
  bool fCreateWithMemory{false};
    
public:
  //! Default constructor
  TPZApproxSpace() = default;
  //! Copy constructor
  TPZApproxSpace(const TPZApproxSpace &copy) = default;
  //! Copy assignment operator
  TPZApproxSpace &operator=(const TPZApproxSpace &copy) = default;
  //! Move constructor
  TPZApproxSpace(TPZApproxSpace &&copy) = default;
  //! Move assignment operator
  TPZApproxSpace &operator=(TPZApproxSpace &&copy) = default;
  //! Destructor
  ~TPZApproxSpace() = default;
  /**
   * Functions to be overriden in child classes
   */
  //@{
  //! Class identifier
  int ClassId() const override;
  //! Read method
  void Read(TPZStream &buf, void *context) override;
  //! Write method
  void Write(TPZStream &buf, int withclassid) const override;
  //! Creates a computational element using the function pointer for the topology
  virtual TPZCompEl *CreateCompEl(TPZGeoEl *gel, TPZCompMesh &mesh) const = 0;
  //@}

  //! Sets whether to create lagrange multiplier
  void SetCreateLagrange(bool flag)
  {
    fCreateLagrangeMultiplier = flag;
  }
  //! Sets whether to create with memory  
  void CreateWithMemory(bool flag)
  {
    fCreateWithMemory = flag;
  }
  /** @brief Determine if the mesh will be created with disconnected elements
   * After the mesh is created, interface elements need to be created "by hand"
   */
  void CreateDisconnectedElements(bool create)
  {
    fCreateHybridMesh = create;
  }
  //! Gets whether the approximation space needs memory at the integration points
  bool NeedsMemory()
  {
    return fCreateWithMemory;
  }
	//! Creates the computational elements and the degrees of freedom for specified materials
	void BuildMesh(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs) const;
	
	//! Creates the computational elements and the degrees of freedom 
	void BuildMesh(TPZCompMesh &cmesh) const;
    
  //! Creates the computational elements and the degrees of freedom for specified elements
  void BuildMesh(TPZCompMesh &cmesh, const TPZVec<int64_t> &gelindexes) const;
    
  /** @brief Creates the interface elements */
	/** Only element of material id in the set<int> will be created */
	static void CreateInterfaces(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs);
	
	//! Creates the interface elements
	static void CreateInterfaces(TPZCompMesh &cmesh);
    
  //! Encapsulate the elements in condensed computational elements
  static void CondenseLocalEquations(TPZCompMesh &cmesh);
    
  //! Undo the encapsulate elements
  static void UndoCondenseLocalEquations(TPZCompMesh &cmesh);

  //! Create interface elements between the computational elements
  static void CreateInterfaceElements(TPZCompMesh *mesh, bool onlydiscontinuous = true, bool multiphysics = false);  
  /// this method will substitute all interface elements with materialid within the set by three elements : one H1 element and two interface elements
  static void Hybridize(TPZCompMesh &cmesh,const std::set<int> &matids, bool isconnectedElem = false);
};

#endif /* _TPZAPPROXSPACE_H_ */
