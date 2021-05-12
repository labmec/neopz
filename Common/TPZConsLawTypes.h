/**
 * @file TPZConsLawTypes.h
 * @brief Contains types used in meshes and materials related to conservation laws.
 */
#ifndef TPZCONSLAWTYPES_H
#define TPZCONSLAWTYPES_H
/**
 * @enum TPZTimeDiscr
 * @brief Indicates the type of time discretization
 * @var None_TD
 * @brief No time discretization.
 * @var Explicit_TD
 * @brief Explicit time discretization. Can to be Euler method, Runge Kutta method, etc.
 * @var ApproxImplicit_TD
 * @brief Semi implicit time discretization
 * @var Implicit_TD
 * @brief Implicit time discretization.
 * @var Unknown_TD
 * @brief Unknown time discretization.
 */
enum TPZTimeDiscr
{
	None_TD = -1,
	Explicit_TD = 0,
	ApproxImplicit_TD = 1,
	Implicit_TD = 2,
	Unknown_TD = 3
};

/**
 * @enum TPZContributeTime
 * @brief Indicates which term is put in the right hand side and tangent matrix
 */
enum TPZContributeTime
{
	None_CT = -1,
	Last_CT = 0,
	Advanced_CT = 1
};

/**
 * @enum TPZResidualType
 * @brief Which terms are being contributed
 */
enum TPZResidualType
{
	None_RT = -1,
	Residual_RT = 0,
	Flux_RT = 1
};
#endif