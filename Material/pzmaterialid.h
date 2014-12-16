/**
 * @file
 * @brief Contains IDs for materials.
 */

#ifndef PZMATERIALIDH
#define PZMATERIALIDH

/** \addtogroup material
 * @{
 */

/** @brief Id of material class */
const int TPZMATERIALID = 300;
/** @brief Id of Discontinuous Galerkin material */
const int TPZDISCONTINUOUSGALERKIN = 301;
/** @brief Id of Bi-dimensional linear material */
const int TPZMAT2DLINID = 302;
/** @brief Id of Bi-dimensional conservation law material */
const int TPZCONSERVATIONLAW2ID = 303;
/** @brief Id of Euler conservation law material */
const int TPZEULERCONSLAW2ID = 304;
/** @brief Id of artificial diffusivity */
const int TPZARTDIFFID = 305;
/** @brief Id of Boundary condition */
const int TPZBNDCONDID = 306;
/** @brief Id of Elasticity material */
const int TPZELASTICITYMATERIALID = 307;
/** @brief Id of Three-dimensional elasticity material */
const int TPZELASTICITY3DMATERIALID = 308;
/** @brief Id of Three-dimensional material for tests */
const int TPZMATTEST3DID = 309;
/** @brief Id of one-dimensional linear material */
const int TPZNLMAT1D = 310;
/** @brief Id of three-dimensional poisson material */
const int TPZMATPOISSON3D = 311;
/** @brief Id of Elasticity axi material */
const int TPZELASTICITYAXIMATERIALID = 312;
/** @brief Id of Material Data */
const int TPZMATERIALDATAID = 313;
/** @brief Id of material void flux */
const int TPZMATERIALVOIDFLUX = 314;
/** @brief Id of post processing material */
const int TPZPOSTPROCMAT_ID = 315;
const int TPBRFICTICIOUSHORIZONTALWELLFLOWID = 316;
const int TPBR3DDARCYFLOWID = 317;
const int TPBRSKINFUNCTIONID = 318;
const int TPBRCASEWELLSKINDATA = 319;

const int TPBRCONSTANTSKINID = 320;

const int TPBRADJUSTEDPERMEABILITYFUNCTIONID = 321;

const int TPBRPERMFUNCTIONID = 322;

const int TPBRRESERVOIRPERMEABILITYFUNCTIONID = 323;

const int TPZVISCOELASTICITYMATERIALID = 324;

const int TPZUncoupledMultiPhysicsID = 325;

/** @brief Id of Elasticity material for hybrid formulation*/
const int TPZELASTICITYHYBRIDMATERIALID = 326;

/** @brief Id of mhm formulation of Darcy problems  material */
const int TPZMATDARCYMHM = 327;

/** @brief Id of Three-dimensional elasticity material */
const int TPZELASTICITY3DNLINEARMATERIALID = 328;

/** @brief Id of a simple laplace equation */
const int TPZMatLaplacianID = 329;

/** @brief Id of a Lagrange Multiplier Material object */
const int TPZLagrangeMultiplierID = 330;

/** @brief Id of a forcing function representing a poro elastic drag force */
const int TPZBiotForceID = 331;

/** @brief Id of a simple laplace equation with lagrange multiplier */
const int TPZMatLaplacianLagrangeID = 332;

/** @brief Id of a material to the double projection method */
const int TPZMDPMaterialID = 333;

/** @brief Id of the elasto plastic material included in Sest2D */
const int TPZMatElastoPlasticSest2DID = 334;

/** @brief Id of the elastic material included in Sest2D */
const int TPZElasticityMaterialSest2DID = 335;

/** @brief Id of a forcing function representing a poro elastic drag force */
const int TPBrAcidFuncID = 336;



/** @} */

#endif //PZMATERIALIDH
