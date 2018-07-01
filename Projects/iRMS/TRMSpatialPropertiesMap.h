//
//  TRMSpatialPropertiesMap.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMSpatialPropertiesMap__
#define __PZ__TRMSpatialPropertiesMap__

#include <stdio.h>
#include "pzmanvector.h"
#include "pzfmatrix.h"

#include "pzgmesh.h"
#include "tpzhierarquicalgrid.h"
#include "pzcmesh.h"
#include "pzl2projection.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "TRMSpatialMap.h"
#include "pzanalysis.h"
#include "pzinterpolationspace.h"

class TRMSpatialPropertiesMap{
    
    
private:
    
    /** @brief L2 computational mesh for compute spatial properties of SPE10 */
    TPZCompMesh * fSPE10Cmesh;

    /** @brief Gmsh grid file */
    std::string fGridName;
    
    /** @brief Set SPE10 fields file */
    std::pair< std::string , std::string > fPermPorFields;
    
    /** @brief number of blocks i, j and k  */
    TPZStack<int> fNBlocks;
    
    /** @brief size of blocks dx, dy and dz  */
    TPZStack<REAL> fBlocks_sizes;
    
    /** @brief spatial properties model */
    int fMap_model; // map_model = {0 -> constan map, 1 -> linear map, 2 -> kriged map}

    // Constant case
    
    void Kappa_c(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    void lambda_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda, TPZManVector<STATE,10> &state_vars);
    
    void S_e_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &S_e, TPZManVector<STATE,10> &state_vars);
    
    void mu_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars);
    
    void alpha_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &alpha, TPZManVector<STATE,10> &state_vars);
    
    // Omar function
    
    void Kappa_f(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi_f(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    // Omar Interpolation
    
    void Kappa_int(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi_int(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    double phi_cdf(double x);
    
public:
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap();
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap(TRMSpatialPropertiesMap &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap &operator=(TRMSpatialPropertiesMap &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief default destructor */
    ~TRMSpatialPropertiesMap();
    
    /**
     * @defgroup Set and get methods
     * @{
     */
    
    /** @brief Set spatial properties model */
    void SetMapModel(int model){
        fMap_model = model;
    }
    
    /** @brief Get spatial properties model */
    int MapModel()
    {
        return fMap_model;
    }


    /** @brief Geological Stress $\sigma_{0}$ */    
    void S_0(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &s_0);
    
    void Kappa(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    void lambda(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda, TPZManVector<STATE,10> &state_vars);
    
    void S_e(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &S_e, TPZManVector<STATE,10> &state_vars);
    
    void mu(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars);
    
    void alpha(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &alpha, TPZManVector<STATE,10> &state_vars);
    
    /////////////////////////////////////// SPE10 utilities ///////////////////////////////////////////////
    
    
    /** @brief Set SPE10 fields file */
    void SetSpatialFields(TPZStack<int> NBlocks, TPZStack<REAL> Blocks_sizes, std::pair< std::string , std::string > PermPorFields){
        fNBlocks = NBlocks;
        fBlocks_sizes = Blocks_sizes;
        fPermPorFields = PermPorFields;
    }
    
    /** @brief Set SPE10 fields file */
    std::pair< std::string , std::string > & SpatialFields(){
        return fPermPorFields;
    }
    
    /** @brief Get SPE10 fields file */
    TPZStack<int> & NBlocks(){
        return fNBlocks;
    }
    
    /** @brief Get SPE10 fields file */
    TPZStack<REAL> & Blocks_sizes(){
        return fBlocks_sizes;
    }
    
    /** @brief Load spatial properties from SPE10 cartesian intact fields Kx, ky, kz and phi */
    void LoadSPE10Map(bool PrintMapQ);
    
    /** @brief Load spatial properties from SPE10 cartesian intact fields Kx, ky, kz and phi */
    bool ComputePropertieSPE10Map(long & index, TPZVec<STATE> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, REAL & phi);

    /** @brief Insert spatial properties from SPE10 on pz mesh of order zero */
    bool Insert_Inside_Map(int n_data);
    
    /** @brief Get dof for spatial properties from SPE10 on pz mesh with connect solution (kx,ky,kz,phi) */
    void ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes);
    
    /** @brief Create a reservoir-box geometry */
    TPZGeoMesh * CreateGeometricBoxMesh(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy, TPZManVector<REAL,2> dz);
    
    static void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    static void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    static void ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    void ExpandGeomesh(TPZGeoMesh *gmesh, REAL sx, REAL  sy, REAL  sz);
    
    void TraslateGeomesh(TPZGeoMesh *gmesh, TPZVec<REAL> t_vec);
    
    void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);
    
};

#endif /* defined(__PZ__TRMSpatialPropertiesMap__) */
