/**
 * @file
 * @brief Contains TPZBlockDiagonal class which defines block diagonal matrices.
 */

#ifndef _TBLOCKDIAG_
#define _TBLOCKDIAG_


#include <iostream>
#include "pzmatrix.h"
#include "pzvec.h"

/**
 * @brief Defines block diagonal matrices. \ref matrix "Matrix"
 * @ingroup matrix
 * @author Philippe Devloo.
 * @since 01/2001
 */
template <class TVar>
class TPZBlockDiagonal : public TPZMatrix<TVar>
{
	
public:
	/** @brief Simple constructor */
	TPZBlockDiagonal ();
	/**
     * @brief Constructor with initialization parameters
     * @param blocksizes Size of blocks on Block Diagonal matrix
     * @param glob Global matrix which will be blocked
	 */
	TPZBlockDiagonal (const TPZVec<int> &blocksizes, const TPZFMatrix<TVar> &glob );
	/**
     * @brief Constructor with initialization parameters
     * @param blocksizes Size of blocks on Block Diagonal matrix
	 */
	TPZBlockDiagonal (const TPZVec<int> &blocksizes);
	/** @brief Copy constructor */
	TPZBlockDiagonal (const TPZBlockDiagonal & );

  TPZBlockDiagonal* NewMatrix() const override{ return new TPZBlockDiagonal{};}
	CLONEDEF(TPZBlockDiagonal)

	/** @brief Creates a copy from another TPZBlockDiagonal*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override
  {                                                           
    auto *from = dynamic_cast<const TPZBlockDiagonal<TVar> *>(mat);                
    if (from) {                                               
      *this = *from;                                          
    }                                                         
    else                                                      
      {                                                       
        PZError<<__PRETTY_FUNCTION__;                         
        PZError<<"\nERROR: Called with incompatible type\n."; 
        PZError<<"Aborting...\n";                             
        DebugStop();                                          
      }                                                       
  }
	/** @brief Simple destructor */
	~TPZBlockDiagonal();
	
	int    Put(const int64_t row,const int64_t col,const TVar& value ) override;
	const TVar Get(const int64_t row,const int64_t col ) const override;
	
	TVar &operator()(const int64_t row, const int64_t col);
	virtual TVar &s(const int64_t row, const int64_t col) override;

	/** @brief This method don't make verification if the element exist. It is fast than Put */
	int    PutVal(const int64_t row,const int64_t col,const TVar& value ) override;
	/** @brief This method don't make verification if the element exist. It is fast than Get */
	const  TVar GetVal(const int64_t row,const int64_t col ) const override;
	
	/** @brief Computes z = alpha * opt(this)*x + beta * y */
	/** @note z and x cannot overlap in memory */
	void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
				 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
	
	int64_t Dim() const  override { return this->Rows(); }
	
	/** @brief Zeroes all the elements of the matrix. */
	int Zero() override;
	
	/**
	 * @brief Return the choosen block size
	 * @param blockid - block index
	 */
	int GetSizeofBlock(int64_t blockid) {return fBlockSize[blockid];}
	
	void Transpose(TPZMatrix<TVar> *const T) const override;
	virtual int Decompose_LU() override;
	virtual int Decompose_LU(std::list<int64_t> &singular) override;
	
	/** @brief Makes the backward and forward substitutions whether the matrix was LU decomposed */
	virtual int Substitution( TPZFMatrix<TVar> * B ) const override;
	
	/** @brief Updates the values of the matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat) override;
	
	/** @brief This method checks the working of the class */
	static int main();
    
    /** Fill the matrix with random values (non singular matrix) */
    void AutoFill(int64_t dim, int64_t dimj, int symmetric) override;
	
private:
	
	/** @brief Clean data matrix. Zeroes number of columns and rows. */
	int Clear() override;

public:
	/**
     * @brief Initializes current matrix based on blocksize
     * @param blocksize Used to initialize current matrix
	 */
	void Initialize(const TPZVec<int> &blocksize);
	/**
     * @brief Adds a block to current matrix
     * @param i Adds in ith position
     * @param block Block to be added
	 */
	void AddBlock(int64_t i, TPZFMatrix<TVar> &block);
	/**
     * @brief Sets a block in the current matrix
     * @param i Adds in ith position
     * @param block Block to be added
	 */
	void SetBlock(int64_t i, TPZFMatrix<TVar> &block);
	
	/**
     * @brief Gets a block from current matrix
     * @param i Returns teh ith block
     * @param block Contains returned block
	 */
	void GetBlock(int64_t i, TPZFMatrix<TVar> &block);
	
	/**
     @brief Builds a block from matrix
     @param matrix Matrix to build from
	 */
	void BuildFromMatrix(TPZMatrix<TVar> &matrix);
	/**
	 * @brief Prints current matrix data
     * @param message Message to be printed
     * @param out Output device
	 * @param format Output format to print
	 */
	virtual void Print(const char *message, std::ostream &out = std::cout, const MatrixOutputFormat format =EFormatted) const override;
	
	int64_t NumberofBlocks() {return fBlockSize.NElements();}
    public:
int ClassId() const override;

protected:
	inline TVar *&Elem() override
  {
    return fStorage.begin();
  }
  inline const TVar *Elem() const override
  {
    return fStorage.begin();
  }
  inline int64_t Size() const override
  {
    return fStorage.size();
  }
	/** @brief Stores matrix data */
	TPZVec<TVar> fStorage;
	/** @brief Stores blocks data */
	TPZVec<int64_t> fBlockPos;
	/** @brief Stores block sizes data */
	TPZVec<int> fBlockSize;
};

template<class TVar>
inline TVar &TPZBlockDiagonal<TVar>::s(const int64_t row, const int64_t col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return this->operator()(row,col);
}

#endif

