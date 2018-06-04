//
//  TPZSimulationData.h
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#ifndef TPZSimulationData_h
#define TPZSimulationData_h

#include <stdio.h>
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include "tinystr.h"
#include "tinyxml.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"

class TPZSimulationData {
    
protected:
    
    /** @brief Flag that states is the poroperm coupling is mixed */
    bool m_is_mixed_formulation_Q;
    
    /** @brief Spatial refinemenet level */
    int m_h_level;
    
    /** @brief Polynomial order for elasticity component */
    int m_elasticity_order;
    
    /** @brief Polynomial order for diffusion component */
    int m_diffusion_order;
    
    /** @brief Physical dimension of the domain */
    int m_dimesion;

    /** @brief Gravity field */
    TPZManVector<REAL,3> m_g;

    /** @brief Prestress state */
    TPZFMatrix<REAL> m_sigma_0;
    
    /** @brief Number of time steps */
    int m_n_steps;
    
    /** @brief Time step size */
    REAL m_dt;
    
    /** @brief Store time values to be reported */
    TPZStack< REAL , 500 > m_reporting_times;
    
    /** @brief Number of iteration */
    int m_n_iteraions;
    
    /** @brief Residue overal tolerance */
    REAL m_epsilon_res;
    
    /** @brief Correction overal tolerance */
    REAL m_epsilon_cor;
    
    /** @brief Name for the Gmsh geometry file being used */
    std::string m_geometry_file;
    
    /** @brief Neopz geometry description */
    TPZGeoMesh * m_geometry;
    
    /** @brief Name for the vtk files being postprocessed */
    std::string m_vtk_file;
    
    /** @brief Number of vtk resolution during postprocessing */
    int m_vtk_resolution;
    
    /** @brief Vector that storage scalar names for postprocessing */
    TPZManVector<std::string,50> m_scalnames;
    
    /** @brief Vector that storage vector names for postprocessing */
    TPZManVector<std::string,50> m_vecnames;
    
    /** @brief Integer that define the number of regions presented in the geometry */
    int m_n_regions;
    
    /** @brief Material and boundaries identifiers sorted per region */
    TPZManVector<std::pair<int, TPZManVector<int,8>>,8> m_mat_ids;
    
    /** @brief Material properties sorted per region */
    TPZManVector<TPZManVector<REAL,8>,8> m_mat_props;

    // Controled by the kernel
    
    /** @brief Initial state directive */
    bool m_is_initial_state_Q;
    
    /** @brief Current time directive */
    bool m_is_current_state_Q;
    
    /** @brief Current time value */
    REAL m_time;
    
    /** @brief Map that storage all the boundary conditions supported  */
    std::map< std::string,std::pair<int,std::vector<std::string> > >  m_condition_type_to_index_value_names;
    
    /** @brief Map that storage the boundary condition identifier with the numerical values provided  */
    std::map<int, std::vector<REAL> > m_bc_id_to_values;
    
    /** @brief Map that storage the provided bc identifiers with the type of boundary condition  */
    std::map<int, std::string> m_bc_id_to_type;
    

    
public:
    
    
    /** @brief default constructor */
    TPZSimulationData();
    
    /** @brief default constructor */
    TPZSimulationData(const TPZSimulationData & other)
    {
        m_is_mixed_formulation_Q            = other.m_is_mixed_formulation_Q;
        m_h_level                           = other.m_h_level;
        m_elasticity_order                  = other.m_elasticity_order;
        m_diffusion_order                   = other.m_diffusion_order;
        m_dimesion                          = other.m_dimesion;
        m_g                                 = other.m_g;
        m_sigma_0                           = other.m_sigma_0;
        m_n_steps                           = other.m_n_steps;
        m_dt                                = other.m_dt;
        m_reporting_times                   = other.m_reporting_times;
        m_n_iteraions                       = other.m_n_iteraions;
        m_epsilon_res                       = other.m_epsilon_res;
        m_epsilon_cor                       = other.m_epsilon_cor;
        m_geometry_file                     = other.m_geometry_file;
        m_geometry                          = other.m_geometry;
        m_vtk_file                          = other.m_vtk_file;
        m_vtk_resolution                    = other.m_vtk_resolution;
        m_scalnames                         = other.m_scalnames;
        m_vecnames                          = other.m_vecnames;
        m_n_regions                         = other.m_n_regions;
        m_mat_ids                           = other.m_mat_ids;
        m_mat_props                         = other.m_mat_props;
        m_is_initial_state_Q                = other.m_is_initial_state_Q;
        m_is_current_state_Q                = other.m_is_current_state_Q;
        m_time                              = other.m_time;
    }
    
    /** @brief default constructor */
    TPZSimulationData &operator=(const TPZSimulationData &other)
    {
        if (this != & other) // prevent self-assignment
        {
            m_is_mixed_formulation_Q            = other.m_is_mixed_formulation_Q;
            m_h_level                           = other.m_h_level;
            m_elasticity_order                  = other.m_elasticity_order;
            m_diffusion_order                   = other.m_diffusion_order;
            m_dimesion                          = other.m_dimesion;
            m_g                                 = other.m_g;
            m_sigma_0                           = other.m_sigma_0;
            m_n_steps                           = other.m_n_steps;
            m_dt                                = other.m_dt;
            m_reporting_times                   = other.m_reporting_times;
            m_n_iteraions                       = other.m_n_iteraions;
            m_epsilon_res                       = other.m_epsilon_res;
            m_epsilon_cor                       = other.m_epsilon_cor;
            m_geometry_file                     = other.m_geometry_file;
            m_geometry                          = other.m_geometry;
            m_vtk_file                          = other.m_vtk_file;
            m_vtk_resolution                    = other.m_vtk_resolution;
            m_scalnames                         = other.m_scalnames;
            m_vecnames                          = other.m_vecnames;
            m_n_regions                         = other.m_n_regions;
            m_mat_ids                           = other.m_mat_ids;
            m_mat_props                         = other.m_mat_props;
            m_is_initial_state_Q                = other.m_is_initial_state_Q;
            m_is_current_state_Q                = other.m_is_current_state_Q;
            m_time                              = other.m_time;
        }
        return *this;
    }
    
    /** @brief destructor */
    ~TPZSimulationData();
    
    /** @brief Read the xml input file */
    void ReadSimulationFile(char *simulation_file);
    
    /** @brief Set the flag that states is the poroperm coupling is mixed */
    void SetMixedFormulationQ(bool flag) { m_is_mixed_formulation_Q = flag; }
    
    /** @brief Set initial state */
    void SetInitialStateQ(bool state) { m_is_initial_state_Q = state; }
    
    /** @brief Get initial state */
    bool IsInitialStateQ() {return m_is_initial_state_Q;}
    
    /** @brief Set current time state */
    void SetCurrentStateQ(bool state) { m_is_current_state_Q = state; }
    
    /** @brief Get current time state */
    bool IsCurrentStateQ() {return m_is_current_state_Q;}

    /** @brief Setup for reporting times and time step size */
    void SetTimeControls(int n_times, REAL dt);
    
    /** @brief Setup for Newton method controls */
    void SetNumericControls(int n_iterations, REAL epsilon_res, REAL epsilon_cor);
    
    /** @brief Get time values being reported */
    TPZStack< REAL , 500 > ReportingTimes(){
        return m_reporting_times;
    }
    
    /** @brief Time step size */
    REAL dt() { return m_dt; }
    
    /** @brief Set the current time value */
    void SetTime(REAL time) { m_time = time; }
    
    /** @brief Get the current time value */
    REAL t() { return m_time; }
    
    /** @brief Get the  number of time steps */
    int n_steps() { return m_n_steps; }
    
    /** @brief Get the number of iterations steps */
    int n_iterations() { return m_n_iteraions; }
    
    /** @brief Get the residue overal tolerance */
    REAL epsilon_res() { return m_epsilon_res; }
    
    /** @brief Get the correction overal tolerance */
    REAL epsilon_cor() { return m_epsilon_cor; }
    
    
    
    
    /** @brief Get Number of vtk resolution during postprocessing */
    int n_div() { return m_vtk_resolution; }
    
    
    
    
    /** @brief Get the gravity field */
    TPZVec<REAL> & Gravity()
    {
        return m_g;
    }
    
    /** @brief Get prestress state */
    TPZFMatrix<REAL> & Sigma_0()
    {
        return m_sigma_0;
    }
    
    /** @brief Get the neopz geometry description */
    TPZGeoMesh * Geometry()
    {
        return m_geometry;
    }

    /** @brief Get the number of regions presented in the geometry */
    int NumberOfRegions() { return m_n_regions; }
    
    /** @brief Get the material and boundaries identifiers sorted per region */
    TPZManVector<std::pair<int, TPZManVector<int,8>>,8> & MaterialIds() { return m_mat_ids; }
    
    /** @brief Get the material properties sorted per region */
    TPZManVector<TPZManVector<REAL,8>,8> & MaterialProps() { return m_mat_props; }
    
    /** @brief Get the physical dimension of the domain */
    int Dimension() { return m_dimesion; }
    
    /** @brief Get the spatial refinemenet level */
    int HLevel() { return m_h_level; }
    
    /** @brief Get the polynomial order for elasticity component */
    int ElasticityOrder() { return m_elasticity_order; }
    
    /** @brief Get the polynomial order for diffusion component */
    int DiffusionOrder() { return m_diffusion_order; }
    
    /** @brief Print the all members */
    void Print();
    
    /** @brief Print the geometry member */
    void PrintGeometry();
    
    /** @brief Get the map that storage all the boundary conditions supported  */
    std::map< std::string,std::pair<int,std::vector<std::string> > > & ConditionTypeToBCIndex() { return m_condition_type_to_index_value_names; }
    
    /** @brief Get the map that storage the type of boundary condition with the numerical values provided  */
    std::map< int , std::vector<REAL> > & BCIdToBCValues() { return m_bc_id_to_values; }
    
    /** @brief Get the map that storage the provided bc identifiers with the type of boundary condition  */
    std::map<int, std::string> & BCIdToConditionType() { return m_bc_id_to_type; }
    
private:
    
    /** @brief Read the Gmsh file and set the geometry member */
    void ReadGeometry();
    
    /** @brief Fillup the map that storage all the boundary conditions supported */
    void LoadBoundaryConditions();
    
    
    
};


#endif /* TPZSimulationData_h */
