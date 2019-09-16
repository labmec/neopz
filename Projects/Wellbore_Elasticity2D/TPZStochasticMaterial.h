//
/**
 * @author NATHALIA BATALHA
 * @since 09/16/2019.
 * @brief Class to generate stochastic field
 * @brief Variables: Young's modulus and Poisson's ratio
 */

#ifndef TPZSTOCHASTICMATERIAL_H
#define TPZSTOCHASTICMATERIAL_H

#include <iostream>
#include "pzfmatrix.h"
#include "pzvec.h"
#include "GeometricMesh.hpp"
#include "ComputationalMesh.hpp"

#include <math.h>
#include <complex>
#include <string>
#include <random>
#include <chrono>


class TPZStochasticMaterial  {
    
protected:
    
    /** @brief Number of square elements in the mesh */
    int fnSquareElements;
    
    /** @brief Number of square elements in the mesh for inclined wellbores */
    int fmatsize;
    
    /** @brief Defines if the stochastic field is inclined or not
     * @note \f$fstochasticInclined = 1\f$ => Inclined field
     * @note \f$fstochasticInclined != 1\f$ => Vertical field
     */
    int fstochasticInclined;
    
    /** @brief If inclined field, wellbore direction/azimuth */
    REAL fdirection;
    
    /** @brief If inclined field, wellbore inclination */
    REAL finclination;
    
    /** @brief wellbore radius */
    REAL frw;
    
    /** @brief wellbore external radius */
    REAL frext;
    
    /** @brief If inclined, 3D field total height */
    REAL fH;
    
    /** @brief If inclined, elements heigh in the 3D field (implemented for square elements so far) */
    REAL fh;
    
    /** @brief Correlation lenght, given by a scale of fluctuation */
    REAL fscale;
    
    /** @brief Geometric mesh */
    TPZGeoMesh* fgmesh;
    
    /** @brief Correlation Matrix */
    TPZFMatrix<REAL> fK;
    
    /** @brief Random Distribution */
    TPZFMatrix<REAL> fRand_U;
    
    /** @brief Random correlated* distribution */
    TPZFMatrix<REAL> fU;
    
    /** @brief Decomposed matrix from Mathematica */
    TPZFMatrix<REAL> fM;
    
    /** @brief Random correlated* distribution of E */
    TPZFMatrix<REAL> fU_E;
    
    /** @brief Random correlated* distribution of nu */
    TPZFMatrix<REAL> fU_nu;
    
    /** @brief Random distribution type of E
     * @note \f$fE_dist = 1\f$ => Normal distribution
     * @note \f$fE_dist = 2\f$ => Lognormal Distribution
     */
    int fE_dist;
    
    /** @brief Random distribution type of nu
     * @note \f$fnu_dist = 1\f$ => Normal distribution
     * @note \f$fnu_dist = 2\f$ => Lognormal Distribution
     */
    int fnu_dist;
    
    /** @brief Correlation function or kernel
     * @note \f$fE_funct = 1\f$ => Exponential kernel
     * @note \f$fE_funct = 2\f$ => Spherical kernel
     */
    int fE_funct;
    
    /** @brief Correlation function or kernel
     * @note \f$fnu_funct = 1\f$ => Exponential kernel
     * @note \f$fnu_funct = 2\f$ => Spherical kernel
     */
    int fnu_funct;
    
    /** @brief Correlation Matrix of E */
    TPZFMatrix<REAL> fKE;
    
    /** @brief Correlation Matrix of nu */
    TPZFMatrix<REAL> fKnu;
    
public:
    
     TPZStochasticMaterial();
    
    /**
     * @brief Creates an object with:
     * @param geometricMesh geometric mesh
     * @param numSquareElems number of square elements
     * @param stochasticInclined inclined field
     * @param direction wellbore direction
     * @param inclination wellbore inclination
     * @param rw wellbore radius
     * @param rext external radius
     * @param M read correlation matrix decomposed
     * @param scale correlation length
     * @param funcE E function
     * @param funcnu nu function
     * @param distribE E distribution
     * @param distribnu nu distribution
     */
    TPZStochasticMaterial(TPZGeoMesh* geometricMesh, int numSquareElems, int stochasticInclined, REAL direction,
                       REAL inclination, REAL rw, REAL rext, const TPZFMatrix<REAL> &M, REAL scale, int funcE,
                          int funcnu, int distribE, int distribnu);
    
    /** @brief Copy constructor */
    TPZStochasticMaterial(const TPZStochasticMaterial &cp);
    
    /** @brief Destructor  */
    virtual ~TPZStochasticMaterial();
    
    /** @brief Set distribution type and function of E */
    virtual void SetYoungField(int distribution, int function);
    
    /** @brief Set distribution type and function of nu */
    virtual void SetPoissonField(int distribution, int function);
    
    /** @brief Set distribution type: normal or lognormal */
    virtual void SetReadMatrix(const TPZFMatrix<REAL> &M);
    
    /** @brief Set inclined parameters if wellbore is inclined */
    virtual void SetInclinedField(int stochasticInclined,REAL direction, REAL inclination);
    
    /** @brief Set geometric parameters  */
    virtual void SetFieldGeometry();
    
    /** @brief Set inclined geometry */
    virtual void InclinedFieldGeometry();
    
    /** @brief Calculates Correlation either vertical or inclined */
    TPZFMatrix<REAL>  EvaluateCorrelation(int function);
    
    /** @brief Calculates Correlation either vertical or inclined */
    TPZFMatrix<REAL> GetCorrelatedVector(int distribution);
    
    /** @brief Get vector of distribution type */
    TPZFMatrix<REAL> GetDistribution(int matrixSize, int distribution);

    /** @brief Calculates Correlation either vertical or inclined */
    virtual void GetStochasticField( TPZFMatrix<REAL> f_E, TPZFMatrix<REAL> f_nu);
    
    /** @brief Calculates correlation matrix for vertical well */
    TPZFMatrix<REAL> calcCorrelationMatrix(int function);
    
    /** @brief Calculates correlation matrix for inclined well */
    TPZFMatrix<REAL> calcCorrelationMatrixInclined(int function);
    
    /** @brief Print correlation matrix to be decomposed at Mathematica - This should be removed after SVD decomposition!!! */
    virtual void PrintCorrelation();
    
    // Calc Correlation Matrix
    TPZVec<REAL> calcStochasticField();

};

#endif /* defined(TPZSTOCHASTICMATERIAL_H) */
