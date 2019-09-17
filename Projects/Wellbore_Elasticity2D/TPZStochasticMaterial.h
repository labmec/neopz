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
#include "pzfunction.h"
#include "pzerror.h"
#include "tpzverysparsematrix.h"
#include "pzsfulmat.h"

#include <math.h>
#include <complex>
#include <string>
#include <random>
#include <chrono>

class TPZStochasticMaterial {
    
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
    TPZFMatrix<STATE> fK;
    
    /** @brief Random Distribution */
    TPZFMatrix<STATE> fRand_U;
    
    /** @brief Random correlated* distribution */
    TPZFMatrix<STATE> fU;
    
    /** @brief Decomposed matrix from Mathematica */
    TPZFMatrix<STATE> fM;
    
    /** @brief Random distribution of E */
    TPZFMatrix<STATE> fU_E;
    
    /** @brief Random distribution of nu */
    TPZFMatrix<STATE> fU_nu;
    
    /** @brief Random correlated* distribution of E */
    TPZFMatrix<STATE> f_E;
    
    /** @brief Random correlated* distribution of nu */
    TPZFMatrix<STATE> f_nu;
    
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
    TPZFMatrix<STATE> fKE;
    
    /** @brief Correlation Matrix of nu */
    TPZFMatrix<STATE> fKnu;
    
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
                       REAL inclination, REAL rw, REAL rext, const TPZFMatrix<STATE> &M, REAL scale, int funcE,
                          int funcnu, int distribE, int distribnu);
    
    /** @brief Copy constructor */
    TPZStochasticMaterial(const TPZStochasticMaterial &cp);
    
    /** @brief Simple destructor  */
    virtual ~TPZStochasticMaterial();
    
    /** @brief Calculate Stochastic Field for both material properties: Young's modulsu and nu
     * It returns f_E and f_nu, with size: 2D elements number
     * These Matrices(Elnumber,1) should be acessed in TPZMaterial for each elements id */
    TPZFMatrix<STATE> calcStochasticField();
    
    /** @brief Set distribution type and function of E */
    virtual void SetYoungField(int distribution, int function);
    
    /** @brief Set distribution type and function of nu */
    virtual void SetPoissonField(int distribution, int function);
    
    /** @brief Set Matrix to be read - should be deleted after SVD decomposition implementation */
    virtual void SetReadMatrix(const TPZFMatrix<STATE> &M);
    
    /** @brief Set inclined parameters if wellbore is inclined */
    virtual void SetInclinedField(REAL rw, REAL rext, int stochasticInclined,REAL direction, REAL inclination);
    
    /** @brief Set inclined geometry */
    virtual void InclinedFieldGeometry();
    
    /** @brief Calculates Correlation for either vertical or inclined well */
    TPZFMatrix<STATE>  EvaluateCorrelation(int function);
    
    /** @brief Gets random distribution type */
    TPZFMatrix<STATE> GetRandomDistribution(int distribution);
    
    /** @brief Get distribution for either vertical or inclined well */
    TPZFMatrix<STATE> GetDistribution(int matrixSize, int distribution);

    /** @brief Gets stochastic field for Young's modulus and Poisson's ratio */
    virtual void GetStochasticField( TPZFMatrix<STATE> f_E, TPZFMatrix<STATE> f_nu);
    
    /** @brief Calculates correlation matrix for vertical well */
    TPZFMatrix<STATE> calcCorrelationMatrix(int function);
    
    /** @brief Calculates correlation matrix for inclined well */
    TPZFMatrix<STATE> calcCorrelationMatrixInclined(int function);
    
    /** @brief Print correlation matrix to be decomposed at Mathematica - This should be removed after SVD decomposition!!! */
    virtual void PrintCorrelation();

};

#endif /* defined(TPZSTOCHASTICMATERIAL_H) */
