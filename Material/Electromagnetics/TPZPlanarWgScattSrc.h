#ifndef _TPZPLANARWGSCATTERINGSRC_H_
#define _TPZPLANARWGSCATTERINGSRC_H_

#include <Electromagnetics/TPZPlanarWgScatt.h>
#include <TPZMatWithMem.h>


/*
  desired memory:
  sol = solution at the integration point
  dsol = solution's derivatives
*/

//! Data to be stored at each integration point
struct TPZScatteredSol2D{
  CSTATE sol = 0;//< solution value
  TPZManVector<CSTATE,3> dsol = {0,0,0};//<derivative value
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


inline std::ostream& operator<<( std::ostream& out, const TPZScatteredSol2D& t ){
  t.Print(out);
  return out;
}

//! Implements a source for scattering problems of planar waveguides
class TPZPlanarWgScattSrc : public TPZPlanarWgScatt,
                                 public TPZMatWithMem<TPZScatteredSol2D>{
public:
  //! All constructors from base class shall be available
  using TPZPlanarWgScatt::TPZPlanarWgScatt;
  //! Contribution to the integration point
  void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                  TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;
  void FillDataRequirements(TPZMaterialData &data) const override;
  //! Gets propagation constant associated with source
  [[nodiscard]] CSTATE Beta() const{return fBeta;}
  //! Set propagation constant associated with source
  void SetBeta(const CSTATE val){fBeta = val;}
  //! Creates a copy of this instance
  TPZPlanarWgScattSrc * NewMaterial() const override;
  //! Returns name of the class
  std::string Name() const override { return "TPZPlanarWgScattSrc"; }
  //! Unique identifier for serialization purposes
  [[nodiscard]] int ClassId() const override;
  //! Write to stream(serialization method)
  void Write(TPZStream &buf, int withclassid) const override;
  //! Read from stream(serialization method)
  void Read(TPZStream &buf, void *context) override;
private:
  CSTATE fBeta{0};
};


#endif /* _TPZPLANARWGSCATTERINGSRC_H_ */
