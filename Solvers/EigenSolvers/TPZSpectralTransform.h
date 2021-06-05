#ifndef TPZSPECTRALTRANSFORM_H
#define TPZSPECTRALTRANSFORM_H

#include "TPZSavable.h"
#include "pzreal.h"

template<class T>
class TPZKrylovEigenSolver;
template<class T>
class TPZMatrix;
template<class T>
class TPZVec;
template<class T>
class TPZAutoPointer;
/**
   @brief Deifnes the interface for the spectral transformations used in the 
   TPZKrylovEigenSolver class.

   These are based on transforming the spectrum of the original problem as 
   to accelerate convergence. They are specially important when computing 
   internal eigenvalues.
 */
template<class TVar>
class TPZSpectralTransform : public virtual TPZSavable{
  friend class TPZKrylovEigenSolver<TVar>;
public:
  /** @name BasicUsage */
  /** @{*/
  //! Default constructor
  TPZSpectralTransform() = default;
  //! Creates a clone of the spectral transformation
  virtual TPZSpectralTransform<TVar> * Clone() const = 0;
  /** @}*/
protected:
  /** @name Protected */
  /** @{*/
  //! Calculates the appropriate transformation for a generalised EVP
  virtual TPZAutoPointer<TPZMatrix<TVar>> CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A,
                                     TPZAutoPointer<TPZMatrix<TVar>>B) const = 0;
  //! Calculates the appropriate transformation for an EVP
  virtual TPZAutoPointer<TPZMatrix<TVar>> CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A) const = 0;
  //! Calculates the original eigenvalues from the mapped ones
  virtual void TransformEigenvalues(TPZVec<CTVar> &w) const = 0;
  /** @}*/
public:
  /** @name ReadWrite */
  /** @{*/
  int ClassId() const override;
  inline void Write(TPZStream &buf, int withclassid) const override {}
  inline void Read(TPZStream &buf, void *context) override {}
  /** @} */
};

//! Defines the shift of origin spectral transformation
template<class TVar>
class TPZSTShiftOrigin : public TPZSpectralTransform<TVar>{
  friend class TPZKrylovEigenSolver<TVar>;
public:
  /** @name BasicUsage */
  /** @{*/
  //!Default constructor
  TPZSTShiftOrigin() = default;
  //!Creates a shift of origin spectral transformation with a given shift
  inline TPZSTShiftOrigin(const TVar s);
  //! Creates a clone of the spectral transformation
  inline TPZSTShiftOrigin<TVar> *Clone() const override;
  //! Sets shift
  inline void SetShift(TVar s);
  //! Gets shift
  inline TVar Shift() const;
  /** @} */
protected:
  /**@name Protected*/
  /**@{*/
  //! Calculates the appropriate transformation for a generalised EVP
  TPZAutoPointer<TPZMatrix<TVar>> CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A, TPZAutoPointer<TPZMatrix<TVar>>B) const override;
  //! Calculates the appropriate transformation for an EVP
  TPZAutoPointer<TPZMatrix<TVar>> CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A) const override;
  //! Calculates the original eigenvalues from the mapped ones
  void TransformEigenvalues(TPZVec<CTVar> &w) const override;
  /** @} */
public:
  /** @name ReadWrite */
  /** @{*/
  int ClassId() const override;
  void Write(TPZStream &buf, int withclassid) const override;
  void Read(TPZStream &buf, void *context) override;
  /** @} */
protected:
  //! Shift to be applied in the system
  TVar fShift{0};
};

template<class TVar>
TPZSTShiftOrigin<TVar>::TPZSTShiftOrigin(const TVar s)
{
  fShift = s;
}

template<class TVar>
TPZSTShiftOrigin<TVar> * TPZSTShiftOrigin<TVar>::Clone() const
{
  return new TPZSTShiftOrigin<TVar>(*this);
}

template<class TVar>
void TPZSTShiftOrigin<TVar>::SetShift(TVar s)
{
  fShift = s;
}
  
template<class TVar>
TVar TPZSTShiftOrigin<TVar>::Shift() const
{
  return fShift;
}

//! Defines the shift and invert spectral transformation
template<class TVar>
class TPZSTShiftAndInvert : public TPZSTShiftOrigin<TVar>{
  friend class TPZKrylovEigenSolver<TVar>;
  /**@name BasicUsage*/
  /**@{*/
public:
  //! Inherits constructors from base class
  using TPZSTShiftOrigin<TVar>::TPZSTShiftOrigin;
  //! Creates a clone of the spectral transformation
  inline TPZSTShiftAndInvert<TVar> * Clone() const override;
  /**@}*/
protected:
  /**@name BasicUsage*/
  /**@{*/
  //! Calculates the appropriate transformation for a generalised EVP
  TPZAutoPointer<TPZMatrix<TVar>> CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A, TPZAutoPointer<TPZMatrix<TVar>>B) const override;
  //! Calculates the appropriate transformation for an EVP
  TPZAutoPointer<TPZMatrix<TVar>> CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A) const override;
  //! Calculates the original eigenvalues from the mapped ones
  void TransformEigenvalues(TPZVec<CTVar> &w) const override;
  /**@}*/
  /** @name ReadWrite */
  /** @{*/
  int ClassId() const override;
  inline void Write(TPZStream &buf, int withclassid) const override {}
  inline void Read(TPZStream &buf, void *context) override {}
  /** @} */
};

template<class TVar>
TPZSTShiftAndInvert<TVar> * TPZSTShiftAndInvert<TVar>::Clone() const
{
  return new TPZSTShiftAndInvert<TVar>(*this);
}

#define INSTANTIATE_TEMPLATES(TCLASS)                 \
  extern template class TCLASS<float>;                \
  extern template class TCLASS<double>;               \
  extern template class TCLASS<std::complex<float>>;  \
  extern template class TCLASS<std::complex<double>>; \


INSTANTIATE_TEMPLATES(TPZSTShiftOrigin)
INSTANTIATE_TEMPLATES(TPZSTShiftAndInvert)
#undef INSTANTIATE_TEMPLATES
#endif