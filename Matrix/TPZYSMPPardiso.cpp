/**
 * @file
 * @brief Contains the implementation of the TPZFYsmpMatrix methods.
 */

#ifdef USING_MKL
#include "TPZYSMPPardiso.h"
#include "pzfmatrix.h"

#ifdef USING_MKL
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl_spblas.h>
#endif
// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZFYsmpMatrixPardiso<TVar>::TPZFYsmpMatrixPardiso() : TPZRegisterClassId(&TPZFYsmpMatrixPardiso::ClassId),
TPZFYsmpMatrix<TVar>() {

}

template<class TVar>
void TPZFYsmpMatrixPardiso<TVar>::CopyFrom(const TPZMatrix<TVar> *  mat)
{                                                           
  auto *from = dynamic_cast<const TPZFYsmpMatrixPardiso<TVar> *>(mat);
  if (from) {
    *this = *from;
  }
  else
  {
    auto *from = dynamic_cast<const TPZFYsmpMatrix<TVar> *>(mat);
    if (from && from->IsDecomposed() == ENoDecompose) {
      *this = *from;
    }
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: Called with incompatible type\n.";
    PZError<<"Aborting...\n";
    DebugStop();
  }
}

template<class TVar>
int TPZFYsmpMatrixPardiso<TVar>::ClassId() const{
    return Hash("TPZFYsmpMatrixPardiso") ^ TPZFYsmpMatrix<TVar>::ClassId() << 1;
}


template<class TVar>
void
TPZFYsmpMatrixPardiso<TVar>::MultAdd(const TPZFMatrix<TVar> &x,
                                     const TPZFMatrix<TVar> &y,
                                     TPZFMatrix<TVar> &z,
                                     const TVar alpha,
                                     const TVar beta,const int opt) const {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
this->MultAddChecks(x,y,z,alpha,beta,opt);

		//suported MKL types
		if constexpr ((
										(std::is_same_v<TVar,float>) ||
										(std::is_same_v<TVar,double>) ||
										(std::is_same_v<TVar,std::complex<float>>) ||
										(std::is_same_v<TVar,std::complex<double>>)
									 )){
			const int64_t m_rows = this->Rows();
			const int64_t m_cols = this->Cols();
			const int64_t x_rows = x.Rows();

			if(beta==(TVar)0.0){
        const auto zr = opt ? m_cols : m_rows;
        const auto zc = x.Cols();
        z.Redim(zr,zc);
      }else{
        z = y;
        const auto  r = (opt) ? this->Rows() : this->Cols();
        if(r == 0){
          z*=beta;
          return;
        }
      }
      
			const int64_t z_cols = z.Cols();
			const int64_t z_rows = z.Rows();
			
			sparse_status_t status; 
			sparse_operation_t op =
				opt ?
        (opt==1 ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_CONJUGATE_TRANSPOSE)
        : SPARSE_OPERATION_NON_TRANSPOSE;
			sparse_index_base_t idx = SPARSE_INDEX_BASE_ZERO;
			sparse_matrix_t A;
			matrix_descr descr;
      switch(this->GetSymmetry()){
      case SymProp::NonSym:
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        break;
      case SymProp::Sym:
        descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        break;
      case SymProp::Herm:
        descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
        break;
      }
			descr.mode = SPARSE_FILL_MODE_FULL;
			descr.diag = SPARSE_DIAG_NON_UNIT;

			auto CheckStatus = [] (auto status){
				switch(status){
				case SPARSE_STATUS_SUCCESS:
					return;
					break;
				case SPARSE_STATUS_NOT_INITIALIZED:
					std::cout<<"The routine encountered an empty handle or matrix array. "<<std::endl;
					DebugStop();
					break;
				case SPARSE_STATUS_ALLOC_FAILED:
					std::cout<<"Internal memory allocation failed. "<<std::endl;
					DebugStop();
					break;
				case SPARSE_STATUS_INVALID_VALUE:
					std::cout<<"The input parameters contain an invalid value. "<<std::endl;
					DebugStop();
					break;
				case SPARSE_STATUS_EXECUTION_FAILED:
					std::cout<<"Execution failed. "<<std::endl;
					DebugStop();
					break;
				case SPARSE_STATUS_INTERNAL_ERROR:
					std::cout<<"An error in algorithm implementation occurred. "<<std::endl;
					DebugStop();
					break;
				case SPARSE_STATUS_NOT_SUPPORTED:
					std::cout<<"The requested operation is not supported. "<<std::endl;
					DebugStop();
					break;
				}
			};

      static_assert(sizeof(int64_t)==sizeof(long long), "incompatible sizes for MKL");
      long long * ia_b = (long long*) this->fIA.begin();
      long long * ja_b = (long long*) this->fJA.begin();
      
			//create A mat
			if constexpr (std::is_same_v<TVar,double>){
				status = mkl_sparse_d_create_csr(&A,idx, m_rows, m_cols,
                                         ia_b, ia_b+1,
																				 ja_b,this->fA.begin());
				CheckStatus(status);
			}
			else if constexpr (std::is_same_v<TVar,float>){

				status = mkl_sparse_s_create_csr(&A,idx, m_rows, m_cols, ia_b, ia_b+1,
																				 ja_b,this->fA.begin());
				CheckStatus(status);
			}
			else if constexpr (std::is_same_v<TVar,std::complex<double>>){
				status = mkl_sparse_z_create_csr(&A,idx, m_rows, m_cols, ia_b, ia_b+1,
																				 ja_b,this->fA.begin());
				CheckStatus(status);
			}
			else if constexpr (std::is_same_v<TVar,std::complex<float>>){
				status = mkl_sparse_c_create_csr(&A,idx, m_rows, m_cols, ia_b, ia_b+1,
																				 ja_b,this->fA.begin());
				CheckStatus(status);
			}

			for(int c = 0; c < z_cols; c++){
				const TVar *x_ptr = x.Elem() + x_rows*c;
				TVar *z_ptr = z.Elem() + z_rows*c;
				if constexpr (std::is_same_v<TVar,double>){
					status = mkl_sparse_d_mv(op,alpha,A,descr,x_ptr,beta,z_ptr);
					CheckStatus(status);
				}
				else if constexpr (std::is_same_v<TVar,float>){
					status = mkl_sparse_s_mv(op,alpha,A,descr,x_ptr,beta,z_ptr);
					CheckStatus(status);
				}
				else if constexpr (std::is_same_v<TVar,std::complex<double>>){
					status = mkl_sparse_z_mv(op,alpha,A,descr,x_ptr,beta,z_ptr);
					CheckStatus(status);
				}
				else if constexpr (std::is_same_v<TVar,std::complex<float>>){
					status = mkl_sparse_c_mv(op, alpha,A,descr,x_ptr,beta,z_ptr);
					CheckStatus(status);
				}
			}

      status = mkl_sparse_destroy(A);
      CheckStatus(status);
			return;
		}else{
      //unsupported type
      DebugStop();
    }

}

template<class TVar>
void TPZFYsmpMatrixPardiso<TVar>::SetIsDecomposed(DecomposeType val){
  TPZBaseMatrix::SetIsDecomposed(val);
  if(val){fPardisoControl.fDecomposed = true;}
}

template<class TVar>
int TPZFYsmpMatrixPardiso<TVar>::Decompose(const DecomposeType dt)
{
  if(this->fDecomposed && this->fDecomposed != dt){
    this->Error(__PRETTY_FUNCTION__,"matrix is already decomposed with other scheme");
  }

  if(!fPardisoControl.HasCustomSettings()){
    const auto sysType = this->GetSymmetry();
    typename TPZPardisoSolver<TVar>::MProperty prop =
      this->IsDefPositive() ?
      TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
      TPZPardisoSolver<TVar>::MProperty::EIndefinite;
    fPardisoControl.SetMatrixType(sysType,prop);
  }
  fPardisoControl.Decompose(this);
  this->SetIsDecomposed(dt);
  return 0;
}

template<class TVar>
int TPZFYsmpMatrixPardiso<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt)
{
  if(this->fDecomposed && this->fDecomposed != dt){
    this->Error(__PRETTY_FUNCTION__,"matrix is already decomposed with other scheme");
  }
  if(this->fDecomposed == ENoDecompose) this->Decompose(dt);
  const TPZFYsmpMatrixPardiso<TVar>* this_ct = const_cast<const TPZFYsmpMatrixPardiso<TVar>*>(this);
  this_ct->SolveDirect(F,dt);
  return 0;
}

template<class TVar>
int TPZFYsmpMatrixPardiso<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const
{
  if(this->fDecomposed && this->fDecomposed != dt){
    this->Error(__PRETTY_FUNCTION__,"matrix is already decomposed with other scheme");
  }
  if(!this->fDecomposed){
    this->Error(__PRETTY_FUNCTION__,"matrix should've been decomposed already");
  }
  TPZFMatrix<TVar> x(F);
  fPardisoControl.Solve(this,x,F);
  return 0;
}

template class TPZFYsmpMatrixPardiso<double>;
template class TPZFYsmpMatrixPardiso<float>;
template class TPZFYsmpMatrixPardiso<std::complex<float>>;
template class TPZFYsmpMatrixPardiso<std::complex<double>>;


template class TPZFYsmpMatrixPardiso<long double>;
template class TPZFYsmpMatrixPardiso<std::complex<long double>>;
#endif
