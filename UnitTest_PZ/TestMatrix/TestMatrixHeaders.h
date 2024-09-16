#ifndef TESTMATRIXINC
#define TESTMATRIXINC
#include "pzfmatrix.h"
#include "pzsfulmat.h"
#include "pzbndmat.h"
#include "pzsbndmat.h"
#include "pzskylnsymmat.h"
#include "pzskylmat.h"
#include "TPZYSMPMatrix.h"
#include "TPZSYSMPMatrix.h"
#ifdef PZ_USING_MKL
#include "TPZYSMPPardiso.h"
#include "TPZSYSMPPardiso.h"
#endif
#include "pzblockdiag.h"
#include "tpzsparseblockdiagonal.h"
#include "fad.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>

using namespace std::complex_literals;

template<>
struct Catch::StringMaker<long double> {
  static std::string convert(long double ref);
  static int precision;
};

int Catch::StringMaker<long double>::precision = 10;

std::string Catch::StringMaker<long double>::convert(long double value) {
  std::ostringstream out;
  out.precision(precision);
  out << std::fixed << value;
  return out.str();
}

template<class MAT>
struct SymmetricStorage : std::false_type {};
template<class TVar>
struct SymmetricStorage<TPZSFMatrix<TVar>> : std::true_type {};
template<class TVar>
struct SymmetricStorage<TPZSBMatrix<TVar>> : std::true_type {};
template<class TVar>
struct SymmetricStorage<TPZSkylMatrix<TVar>> : std::true_type {};
template<class TVar>
struct SymmetricStorage<TPZSYsmpMatrix<TVar>> : std::true_type {};
#ifdef PZ_USING_MKL
template<class TVar>
struct SymmetricStorage<TPZSYsmpMatrixPardiso<TVar>> : std::true_type {};
#endif

template<class MAT>
inline constexpr bool IsSymmetricStorage = SymmetricStorage<MAT>::value;

#endif
