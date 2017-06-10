#ifndef TPZEulerBernoulliBeamDataH
#define TPZEulerBernoulliBeamDataH

#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <map>

/** Data structure of material and section data required
 * by TPZEulerBernoulliBeam computations
 * @author Tiago
 * @since Feb 2017
 */
class TPZEulerBernoulliBeamData {

  private:

  //1 means stiffness, 0 means mass
  int fCurrenttlyComputingStiffness;

  //Gravity vector: by default Z is adopted as the vertical upward direction, i.e. gravity = -Z, but user can change it
  TPZFMatrix<STATE> fGravity;

  public:

  TPZEulerBernoulliBeamData();

  ~TPZEulerBernoulliBeamData();

    //Returns if stiffness matrix is to be computed
  bool ComputeStiffness() const{ return (this->fCurrenttlyComputingStiffness); }

  //Returns if mass matrix is to be computed
  bool ComputeMass() const{ return !(this->fCurrenttlyComputingStiffness); }

  //Sets stiffness matrix is to be computed
  void SetComputeStiffness() { this->fCurrenttlyComputingStiffness = 1; }

  //Set mass matrix is to be computed
  void SetComputeMass() { this->fCurrenttlyComputingStiffness = 0; }

  ///Material properties: Young's module, transversal module and density
  struct MaterialProperties{
    REAL fE, fG, fRho;
    void Print(std::ostream & out) const{
      out << fE << " fE\n";
      out << fG << " fG\n";
      out << fRho << " fRho\n";
    }
  };

  /** Section properties:
   * fA: cross section area
   * fIy and fIz: moments of inertia
   * fJ = fIp for circular sections for which the axis of twist and the centroidal axis coincide.
   */
  struct SectionProperties{
    REAL fA, fIy, fIz, fJt, fIp;
    void Print(std::ostream & out) const{
      out << fA << " fA\n";
      out << fIy << " fIy\n";
      out << fIz << " fIz\n";
      out << fJt << " fJt\n";
      out << fIp << " fIp\n";
    }
  };

  //Material properties: MaterialId -> data
  std::map<int, MaterialProperties> fMaterialProp;

  ///Section properties: SectionId -> data
  std::map<int, SectionProperties> fSectionProp;

  void Print(std::ostream &out) const;

  void GetData(int matId, int sectionId, MaterialProperties& mat, SectionProperties& sec) const;

  ///Gravity vector. Z is adopted as the vertical upward direction, i.e. gravity = -Z
  void g(TPZFMatrix<STATE> & g) const;

  void SetGravity(TPZFMatrix<STATE> &g);

};
#endif

