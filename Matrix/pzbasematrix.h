/**
 * @file
 * @brief Contains TPZBaseMatrix class, root matrix class.
 * Since this class is not a template class, its will only
 * define methods that take no other matrices.
 */

#ifndef TPZBASEMATRIX_H
#define TPZBASEMATRIX_H

#include "TPZSavable.h"
#include "pzreal.h"

/** @brief Defines output format */
enum MatrixOutputFormat {EFormatted, EInputFormat, EMathematicaInput, EMatlabNonZeros, EMatrixMarket, ECSV, EFixedColumn};

class TPZStream;
/** \addtogroup matrix
 * @{
 */
/**
 * @enum DecomposeType
 * @brief Defines decomposition type for any matrix classes
 * @param ENoDecompose Not decomposed
 * @param ELU Decomposed using LU method
 * @param ECholesky Decomposed using Cholesky method
 * @param ELDLt Decomposed using LDLt method
 */
enum DecomposeType { ENoDecompose, ELU, ELUPivot, ECholesky, ELDLt };

/** @brief Root matrix class (abstract). \ref matrix "Matrix"
 * Abstract class TPZBaseMatrix which is agnostic with respect to
 * the arithmetic type being used. */
class TPZBaseMatrix : public TPZSavable{
public:
  /** @brief Default constructor */
  TPZBaseMatrix() : TPZRegisterClassId(&TPZBaseMatrix::ClassId) {
    fDecomposed = 0;
    fDefPositive = 0;
    fRow = 0;
    fCol = 0;
  }
  /** @brief Copy constructor*/
  TPZBaseMatrix(const TPZBaseMatrix &cp) = default;
  /** @brief Move constructor*/
  TPZBaseMatrix(TPZBaseMatrix &&cp) = default;
  /** @brief Simple destructor */
  virtual ~TPZBaseMatrix() = default;

  /** @brief Copy assignment operator*/
  TPZBaseMatrix &operator=(const TPZBaseMatrix &) = default;
  /** @brief Move assignment operator*/
  TPZBaseMatrix &operator=(TPZBaseMatrix &&) = default;

  /**
   * @brief Returns the approximate size of the memory footprint (amount
   * of memory required to store this object).
   */
  virtual int64_t MemoryFootprint() const = 0;

  /** @brief Fill matrix storage with randomic values */
  virtual void AutoFill(int64_t nrow, int64_t ncol, int symmetric) = 0;
    
  /** @brief Checks if current matrix value is symmetric */
  virtual int VerifySymmetry(REAL tol) const = 0;

  /** @brief Returns number of rows */
  inline int64_t Rows() const { return fRow; }
  /** @brief Returns number of cols */
  inline int64_t Cols() const { return fCol; }

  /** @brief Returns the dimension of the matrix if the matrix is square.*/
  /** If the matrix is not square, returns an error */
  inline virtual int64_t Dim() const {
    if (IsSquare())
      return Rows();
    PZError << "Cannot call TPZBaseMatrix::Dim() for a non-square ";
    PZError << "matrix! Aborting...";
    DebugStop();
    return 0; // compiler wont complain
  }

  /**
   * @brief Redimensions a matriz keeping the previous values
   * @param newRows Specifies the new number of rows in matrix
   * @param newCols Specifies the new number of Columns in matrix
   */
  virtual int Resize(const int64_t newRows, const int64_t newCols) = 0;
  /**
   * @brief Redimensions the matrix reinitializing it with zero
   * @param newRows Specifies the new number of rows in matrix.
   * @param newCols Specifies the new number of Columns in matrix.
   */
  virtual int Redim(const int64_t newRows, const int64_t newCols) = 0;

  /** @brief Zeroes the matrix */
  virtual int Zero() = 0;


  /** @brief Checks if the current matrix is symmetric */
  virtual int IsSymmetric() const{ return 0;}
  /** @brief Checks if current matrix is square */
  inline int IsSquare() const {
    return fRow == fCol; }

  /** @brief Simetrizes copies upper plan to the lower plan, making its data
   * simetric */
  virtual void Simetrize() = 0;

  /** @brief Checks if current matrix is definite positive */
  virtual int IsDefPositive() const {
    return fDefPositive; }
  /** @brief Checks if current matrix is already decomposed */
  int IsDecomposed() const {
    return fDecomposed; }



  /** @brief Sets current matrix to decomposed state */
  virtual void SetIsDecomposed(int val) {
    fDecomposed = (char)val; }

  /** @brief decompose the system of equations acording to the decomposition
   * scheme */
  virtual int Decompose(const DecomposeType dt, std::list<int64_t> &singular) {
    switch (dt) {
    case ELU:
      return Decompose_LU(singular);
      break;
    case ELDLt:
      return Decompose_LDLt(singular);
      break;
    case ECholesky:
      return Decompose_Cholesky(singular);
      break;
    default:
      DebugStop();
      break;
    }
    return -1;
  }

  /** @brief decompose the system of equations acording to the decomposition
   * scheme */
  virtual int Decompose(const DecomposeType dt) {
    switch (dt) {
    case ELU:
      return Decompose_LU();
      break;
    case ELDLt:
      return Decompose_LDLt();
      break;
    case ECholesky:
      return Decompose_Cholesky();
      break;
    default:
      DebugStop();
      break;
    }
    return -1;
  }
  /**
   * @name Factorization
   * @brief Those member functions perform the matrix factorization
   * @{
   */

  /** @brief Decomposes the current matrix using LU decomposition. */
  virtual int Decompose_LU(std::list<int64_t> &singular) = 0;
  /** @brief Decomposes the current matrix using LU decomposition. */
  virtual int Decompose_LU() = 0;

  
  /**
   * @brief Decomposes the current matrix using Cholesky method.
   *
   * The curent matrix has to be hermitian.
   */
  virtual int Decompose_Cholesky(std::list<int64_t> &singular) = 0;
  /** @brief Decomposes the current matrix using Cholesky*/
  virtual int Decompose_Cholesky() = 0;
  /**
   * @brief Decomposes the current matrix using LDLt.
   *
   * The current matrix has to be symmetric.
   * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal
   * matrix.
   */
  virtual int Decompose_LDLt(std::list<int64_t> &singular) = 0;
  /** @brief Decomposes the current matrix using LDLt. */
  virtual int Decompose_LDLt() = 0;
  
  /** @} */

  /** @brief It prints the matrix data in a MatrixFormat Rows X Cols */
  virtual void Print(const char *name, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const = 0;
  int ClassId() const override;

  /**
   * @brief Unpacks the object structure from a stream of bytes
   * @param buf The buffer containing the object in a packed form
   * @param context
   */
  void Read(TPZStream &buf, void *context) override;

  /**
   * @brief Packs the object structure in a stream of bytes
   * @param buf Buffer which will receive the bytes
   * @param withclassid
   */
  void Write(TPZStream &buf, int withclassid) const override;

protected:
  /** @brief Number of rows in matrix */
  int64_t fRow;
  /** @brief Number of cols in matrix */
  int64_t fCol;
  /** @brief Decomposition type used to decompose the current matrix */
  char fDecomposed;
  /** @brief Definite Posistiveness of current matrix */
  char fDefPositive;
};

/** @} */

#endif
