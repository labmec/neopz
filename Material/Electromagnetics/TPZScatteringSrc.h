#ifndef _TPZSCATTERINGSRC_H_
#define _TPZSCATTERINGSRC_H_

#include <Electromagnetics/TPZScattering.h>
#include <TPZMatWithMem.h>


/*
  desired memory:
  sol = solution at the integration point
*/

//! Data to be stored at each integration point
struct TPZScatteredSol3D{
  TPZManVector<CSTATE,2> sol = {0,0};//< solution value (tg component)
  TPZManVector<REAL,3> x = {0,0,0};//<for debugging purposes
  //! Unique identifier for serialization purposes
  [[nodiscard]] int ClassId() const;
  //! Write to stream(serialization method)
  void Write(TPZStream &buf, int withclassid) const;
  //! Read from stream(serialization method)
  void Read(TPZStream &buf, void *context);
  //! Print contents(debugging method)
  void Print(std::ostream &out) const;
};


inline std::ostream& operator<<( std::ostream& out, const TPZScatteredSol3D& t ){
  t.Print(out);
  return out;
}

//! Implements a source for scattering problems of 2D waveguides
class TPZScatteringSrc : public TPZScattering,
                                 public TPZMatWithMem<TPZScatteredSol3D>{
public:
  //! All constructors from base class shall be available
  using TPZScattering::TPZScattering;
  //! Contribution to the matrix and rhs at the integration point
  void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                  TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override
  {
    Contribute(data,weight,ef);
  }
  //! Contribution to the rhs at the integration point
  void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                  TPZFMatrix<CSTATE> &ef) override;
  
  void FillDataRequirements(TPZMaterialData &data) const override;
  //! Gets propagation constant associated with source
  [[nodiscard]] CSTATE Beta() const{return fBeta;}
  //! Set propagation constant associated with source
  void SetBeta(const CSTATE val){
    fBeta = val;
  }
  //! Returns the integrable dimension of the material
  int Dimension() const override {return 2;}
  //! Creates a copy of this instance
  TPZScatteringSrc * NewMaterial() const override;
  //! Returns name of the class
  std::string Name() const override { return "TPZScatteringSrc"; }
  //! Unique identifier for serialization purposes
  [[nodiscard]] int ClassId() const override;
  //! Write to stream(serialization method)
  void Write(TPZStream &buf, int withclassid) const override;
  //! Read from stream(serialization method)
  void Read(TPZStream &buf, void *context) override;
private:
  CSTATE fBeta{0};
};


#endif /* _TPZSCATTERINGSRC_H_ */
