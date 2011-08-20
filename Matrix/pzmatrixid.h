/**
 * @file
 * @brief Contains TPZFMATRIXID and TPZBLOCKID ids.
 */

#ifndef PZMATRIXIDH
#define PZMATRIXIDH

/** \addtogroup matrix 
 * @{ */
/** @brief Id of full matrix */
const int TPZFMATRIXID = 100;
/** @brief Id of the diagonal block matrix */
const int TPZBLOCKID = 101;

/** @brief Function to register matrix classes */
void RegisterMatrixClasses();

/** @} */

#endif
