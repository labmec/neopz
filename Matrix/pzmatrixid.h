/**
 * @file
 * @brief Contains TPZMatrix<TVar>class, root matrix class.
 */

#ifndef _TMATRIXIDHH_
#define _TMATRIXIDHH_

/**
 * @addtogroup matrix
 * @{
 */

/** @brief Id of TPZMATRED with matrix very sparse */
const int TPZMATRED_VERYSPARSE_ID = 500;
/** @brief If of TPZMATRED with full matrix */
const int TPZMATRED_FMATRIX_ID = 501;
/** @brief Id of full matrix */
const int TPZFMATRIX_DOUBLE_ID = 100;
const int TPZFMATRIX_FLOAT_ID = 101;
const int TPZFMATRIX_COMPLEX_ID = 102;
const int TPZFMATRIX_LONG_DOUBLE_ID = 103;

/** @brief Id of the DohrMatrix condense */
const int TPZDOHRMATRIXSUBSTRUCTCONDENSEFLOAT = 505;
const int TPZDOHRMATRIXSUBSTRUCTCONDENSEDOUBLE = 506;
const int TPZDOHRMATRIXSUBSTRUCTCONDENSELONGDOUBLE = 507;
/** @brief Id of the DohrMatrix */
const int TPZDOHRMATRIXSUBSTRUCTFLOAT = 510;
const int TPZDOHRMATRIXSUBSTRUCTDOUBLE = 511;
const int TPZDOHRMATRIXSUBSTRUCTLONGDOUBLE = 512;

/** @brief Id of the DohrMatrix */
const int TPZDOHRPRECOND_FLOAT_ID = 515;
const int TPZDOHRPRECOND_DOUBLE_ID = 516;
const int TPZDOHRPRECOND_LONGDOUBLE_ID = 517;

const int TPZDOHRPRECONDCONDENSE_FLOAT_ID = 518;
const int TPZDOHRPRECONDCONDENSE_DOUBLE_ID = 519;
const int TPZDOHRPRECONDCONDENSE_LONGDOUBLE_ID = 520;

const int TPZSTEPSOLVERFLOAT_ID = 525;
const int TPZSTEPSOLVERDOUBLE_ID = 526;
const int TPZSTEPSOLVERLONGDOUBLE_ID = 527;

const int TPZSTEPSOLVERCOMPLEXFLOAT_ID = 530;

const int TPZSTEPSOLVERCOMPLEXDOUBLE_ID = 535;
const int TPZSTEPSOLVERCOMPLEXLONGDOUBLE_ID = 536;

/** @brief Id of the Skyline Matrix */
const int TSKYLMATRIX_DOUBLE_ID = 540;
const int TSKYLMATRIX_FLOAT_ID = 541;

/** @} */

#endif

