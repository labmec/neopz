//
//  TPZSimulationData.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPZSimulationData.h"

/** @brief costructor */
TPZSimulationData::TPZSimulationData()
{
    
    m_is_mixed_formulation_Q = false;
    m_h_level = 0;
    m_elasticity_order = 0;
    m_diffusion_order = 0;
    m_dimesion = 0;
    m_g.Resize(0);
    m_sigma_0.Resize(0, 0);
    m_n_steps = 0;
    m_dt = 0.0;
    m_reporting_times.resize(0);
    m_n_iteraions = 0;
    m_epsilon_res = 1.0;
    m_epsilon_cor = 1.0;
    m_n_threads = 0;
    m_geometry_file = "";
    m_geometry = NULL;
    m_vtk_file = "";
    m_vtk_resolution = 0;
    m_scalnames.Resize(0);
    m_vecnames.Resize(0);
    m_n_regions = 0;
    m_mat_ids.Resize(0);
    m_mat_props.Resize(0);
    m_is_initial_state_Q = false;
    m_is_current_state_Q = false;
    m_time = 0.0;
    
}

/** @brief destructor */
TPZSimulationData::~TPZSimulationData()
{
    
}

/** @brief simulation file reader */
void TPZSimulationData::ReadSimulationFile(char *simulation_file)
{
    
    TiXmlDocument document(simulation_file);
    bool file_ok_Q = false;
    
    file_ok_Q = document.LoadFile();
    if (file_ok_Q)
    {
        std::cout << "This Xml is ok! -> " << simulation_file << std::endl;
    }
    else
    {
        std::cout << "Failed to load file \n" << simulation_file << std::endl;
        std::cout <<    "Check the given path or your xml structure. \n" << std::endl;
        DebugStop();
    }
    
    
    // TiXmlElement dummy object
    TiXmlElement * container;
    const char * char_container;
    TiXmlHandle doc_handler( & document );
    
    // Begin:: Geometry description
    container = doc_handler.FirstChild("CaseData").FirstChild("Mesh").FirstChild("MeshFile").ToElement();
    const char * geometry_file = container->Attribute("mesh_file");
    m_geometry_file = geometry_file;
    
    this->ReadGeometry();
    int dimension = m_geometry->Dimension();
    m_dimesion = dimension;
    // End:: Geometry description
    
    
    // Begin:: Time controls
    container = doc_handler.FirstChild("CaseData").FirstChild("TimeControls").FirstChild("StepSize").ToElement();
    char_container = container->Attribute("dt");
    REAL dt = std::atof(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("TimeControls").FirstChild("StepNumber").ToElement();
    char_container = container->Attribute("n_time_steps");
    int n_stpes = std::atoi(char_container);
    
    SetTimeControls(n_stpes,dt);
    // End:: Time controls
    
    
    // Begin:: Newton method controls
    container = doc_handler.FirstChild("CaseData").FirstChild("NewtonControls").FirstChild("Iterations").ToElement();
    char_container = container->Attribute("n_iterations");
    int n_iterations = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("NewtonControls").FirstChild("Residue").ToElement();
    char_container = container->Attribute("res_tolerance");
    REAL epsilon_res = std::atof(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("NewtonControls").FirstChild("Correction").ToElement();
    char_container = container->Attribute("cor_tolerance");
    REAL epsilon_cor = std::atof(char_container);
    
    SetNumericControls(n_iterations,epsilon_res,epsilon_cor);
    // End:: Newton method controls
    
    // Begin:: Parallel controls
    container = doc_handler.FirstChild("CaseData").FirstChild("ParallelControls").FirstChild("Numthreads").ToElement();
    char_container = container->Attribute("n_threads");
    int n_threads = std::atoi(char_container);
    m_n_threads = n_threads;
    // End:: Parallel controls
    
    // Begin:: Finite elements
    container = doc_handler.FirstChild("CaseData").FirstChild("FEM").FirstChild("MixedFormulationQ").ToElement();
    char_container = container->Attribute("useQ");
    bool is_mixed_formulation_Q = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("FEM").FirstChild("HRefine").ToElement();
    char_container = container->Attribute("h_level");
    int h_level = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("FEM").FirstChild("Elasticity").ToElement();
    char_container = container->Attribute("p_order");
    int elasticity_order = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("FEM").FirstChild("Diffusion").ToElement();
    char_container = container->Attribute("p_order");
    int diffusion_order = std::atoi(char_container);
    
    m_is_mixed_formulation_Q    = is_mixed_formulation_Q;
    m_h_level                   = h_level;
    m_elasticity_order          = elasticity_order;
    m_diffusion_order           = diffusion_order;
    // End:: Finite elements
    
    
    // Begin:: Outputs
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("OutputFolder").ToElement();
    char_container = container->Attribute("name");
    system(char_container);
    std::cout << "Creating directory using : "<< char_container << "\n" << std::endl;
    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("LogFolder").ToElement();
    char_container = container->Attribute("name");
    system(char_container);
    std::cout << "Creating directory using : "<< char_container << "\n" << std::endl;
    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("PostProcessing").ToElement();
    const char * vtk_file = container->Attribute("vtk_file");
    m_vtk_file = vtk_file;
    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("PostProcessing").ToElement();
    char_container = container->Attribute("n_divisions");
    int vtk_resolution = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("PostProcessing").ToElement();
    char_container = container->Attribute("n_scalars");
    int n_scalars = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("PostProcessing").ToElement();
    char_container = container->Attribute("n_vectors");
    int n_vectors = std::atoi(char_container);
    
    TPZManVector<std::string,50> scalnames(n_scalars), vecnames(n_vectors);

    int iscalar = 0;
    container = doc_handler.FirstChild( "CaseData" ).FirstChild( "OutputControls" ).FirstChild( "Scalars" ).FirstChild( "Var" ).ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        char_container = container->Attribute("name");
        scalnames[iscalar] = char_container;
        iscalar++;
    }

    int ivectorial = 0;
    container = doc_handler.FirstChild( "CaseData" ).FirstChild( "OutputControls" ).FirstChild( "Vectors" ).FirstChild( "Var" ).ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        char_container = container->Attribute("name");
        vecnames[ivectorial] = char_container;
        ivectorial++;
    }
    
    m_vtk_file = vtk_file;
    m_vtk_resolution = vtk_resolution;
    m_scalnames = scalnames;
    m_vecnames = vecnames;
    // End:: Outputs
    
    
    // Begin:: Physics
    container = doc_handler.FirstChild("CaseData").FirstChild("Physics").FirstChild("GravityConstant").ToElement();
    char_container = container->Attribute("gravity");
    REAL g_c = std::atof(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("Physics").FirstChild("GravityDirection").ToElement();
    char_container = container->Attribute("x_direction");
    REAL x_dir = std::atof(char_container);
    char_container = container->Attribute("y_direction");
    REAL y_dir = std::atof(char_container);
    m_g.Resize(dimension,0.0);
    m_sigma_0.Resize(3, 3);
    m_sigma_0.Zero();
    m_g[0] = g_c * x_dir;
    m_g[1] = g_c * y_dir;
    
    if (dimension == 3)
    {
        char_container = container->Attribute("z_direction");
        REAL z_dir = std::atof(char_container);
        m_g[2] = g_c * z_dir;
    }
    
    container = doc_handler.FirstChild("CaseData").FirstChild("Physics").FirstChild("Prestress").ToElement();
    char_container = container->Attribute("sxx");
    REAL sxx = std::atof(char_container);
    char_container = container->Attribute("syy");
    REAL syy = std::atof(char_container);
    char_container = container->Attribute("sxy");
    REAL sxy = std::atof(char_container);
    
    m_sigma_0(0,0) = sxx;
    m_sigma_0(1,1) = syy;
    m_sigma_0(1,0) = sxy;
    m_sigma_0(0,1) = sxy;
    
    if (dimension == 3)
    {
        
        char_container = container->Attribute("szz");
        REAL szz = std::atof(char_container);
        char_container = container->Attribute("sxz");
        REAL sxz = std::atof(char_container);
        char_container = container->Attribute("syz");
        REAL syz = std::atof(char_container);
        
        m_sigma_0(2,2) = szz;
        m_sigma_0(0,2) = sxz;
        m_sigma_0(2,0) = sxz;
        m_sigma_0(1,2) = syz;
        m_sigma_0(2,1) = syz;
    
    }
    // End:: Physics
    
    
    // Begin:: Regions and materials parameters
    container = doc_handler.FirstChild("CaseData").FirstChild("ReservoirRegions").FirstChild("RegionNumber").ToElement();
    char_container = container->Attribute("n_regions");
    int n_regions = std::atoi(char_container);
    
    m_n_regions = n_regions;
    m_mat_ids.Resize(n_regions);
    m_mat_props.Resize(n_regions);
    
    // Block that will increase with the number of parameters that the PDE require
    // This follow the Parameters structure in the input xml file
    TPZStack<std::string> par_names;
    par_names.Push("eyoung");
    par_names.Push("nu");
    par_names.Push("phi");
    par_names.Push("kappa");
    par_names.Push("alpha");
    par_names.Push("m");
    par_names.Push("rho");
    par_names.Push("mu");

    
    
    int iregion = 0;
    TiXmlElement * sub_container;
    container = doc_handler.FirstChild( "CaseData" ).FirstChild( "RegionsDefinition" ).FirstChild("RegionData").ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        
        char_container = container->Attribute("mat_id");
        int mat_id = std::atoi(char_container);
        char_container = container->Attribute("n_boundaries");
        int n_boundaries = std::atoi(char_container);
        char_container = container->Attribute("n_parameters");
        int n_parameters = std::atoi(char_container);
        
        m_mat_ids[iregion].first = mat_id;
        m_mat_ids[iregion].second.Resize(n_boundaries);
        
        sub_container = container->FirstChild("Boundaries")->FirstChild("Boundary")->ToElement();
        int iboundary = 0;
        for( ; sub_container; sub_container=sub_container->NextSiblingElement())
        {
            char_container = sub_container->Attribute("bc_id");
            int bc_id = std::atoi(char_container);
            m_mat_ids[iregion].second [iboundary] = bc_id;
            iboundary++;
        }
        
        
        m_mat_props[iregion].Resize(n_parameters);
        sub_container = container->FirstChild("Parameters")->ToElement();
        
#ifdef PZDEBUG
        if (n_parameters != par_names.size())
        {
            std::cout << "The list of parameters does not conicides with the impelemented during reading the xml." << std::endl;
            DebugStop();
        }
#endif
    
        for (int ipar = 0; ipar < n_parameters; ipar++)
        {
            char_container = sub_container->Attribute(par_names[ipar].c_str());
            REAL par = std::atof(char_container);
            m_mat_props[iregion][ipar] = par;
        }
        
        iregion++;
    }
    // End:: Regions and materials parameters
    
    
    // Begin::  Block that define the material parameters
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
    int eyoung = 0, nu = 1, phi = 2, kappa = 3, alpha = 4, m = 5, rho = 6, mu = 7;
    m_young = m_mat_props[iregion][eyoung];
    m_nu = m_mat_props[iregion][nu];
    m_porosity = m_mat_props[iregion][phi];
    m_k = m_mat_props[iregion][kappa];
    m_alpha = m_mat_props[iregion][alpha];
    m_Se = 1.0/m_mat_props[iregion][m];
    m_eta = m_mat_props[iregion][mu];
    m_rho_f = m_mat_props[iregion][rho];
    }
    // End::  Block that define the material parameters

    
    
    // Begin:: Regions and materials parameters
    this->LoadBoundaryConditions();
    
    std::pair<int, std::string > bc_id_to_type_chunk;
    std::pair<int , std::vector<REAL> > bc_id_to_values_chunk;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator chunk;
    container = doc_handler.FirstChild("CaseData").FirstChild("BoundaryConditions").FirstChild("BoundaryCondition").ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        
        char_container = container->Attribute("bc_id");
        int bc_id = std::atoi(char_container);
        
        char_container = container->Attribute("type");
        std::string condition(char_container);
        chunk = m_condition_type_to_index_value_names.find(condition);
        
        // Association bc type with numerical values
        bc_id_to_values_chunk.first = bc_id;
        bc_id_to_values_chunk.second.resize(0);
        int n_data = chunk->second.second.size();
        for (int i = 0; i < n_data; i++)
        {
            char_container = container->Attribute(chunk->second.second[i].c_str());
#ifdef PZDEBUG
            if (!char_container)
            {
                std::cout << " the boundary " << condition << "  needs the value " << chunk->second.second[i] << std::endl;
                std::cout << " Please review your boundary condition definitions. " << std::endl;
                DebugStop();
            }
#endif
            REAL bc_value = std::atof(char_container);
            bc_id_to_values_chunk.second.push_back(bc_value);
        }
        m_bc_id_to_values.insert(bc_id_to_values_chunk);

        // Association bc identifier with bc type
        bc_id_to_type_chunk.first = bc_id;
        bc_id_to_type_chunk.second = condition;
        m_bc_id_to_type.insert(bc_id_to_type_chunk);
        
        
    }
    // End:: Regions and materials parameters
    
}

/** @brief Setup reporting times and time step size */
void TPZSimulationData::SetTimeControls(int n_times, REAL dt)
{
    
    m_n_steps    = n_times;
    m_dt         = dt;
    m_reporting_times.Resize(n_times, 0.0);
    for (int it = 0; it < n_times; it++)
    {
        m_reporting_times[it] = it*dt;
    }
    
}

/** @brief Setup reporting times and time step size */
void TPZSimulationData::SetNumericControls(int n_iterations, REAL epsilon_res, REAL epsilon_cor)
{
    
    m_n_iteraions  =   n_iterations;
    m_epsilon_res    =   epsilon_res;
    m_epsilon_cor    =   epsilon_cor;
    
}

/** @brief Print the all members */
void TPZSimulationData::Print()
{
    
    std::cout << " TPZSimulationData class members : " << std::endl;
    std::cout << std::endl;
    std::cout << " m_is_mixed_formulation_Q = " << m_is_mixed_formulation_Q << std::endl;
    std::cout << " m_h_level = " << m_h_level << std::endl;
    std::cout << " m_elasticity_order = " << m_elasticity_order << std::endl;
    std::cout << " m_diffusion_order = " << m_diffusion_order << std::endl;
    std::cout << " m_dimesion = " << m_dimesion << std::endl;
    std::cout << " m_g = " << m_g << std::endl;
    std::cout << " m_sigma_0 = " << m_sigma_0 << std::endl;
    std::cout << " m_n_steps = " << m_n_steps << std::endl;
    std::cout << " m_dt = " << m_dt << std::endl;
    std::cout << " m_reporting_times = " << m_reporting_times << std::endl;
    std::cout << " m_n_iteraions = " << m_n_iteraions << std::endl;
    std::cout << " m_epsilon_res = " << m_epsilon_res << std::endl;
    std::cout << " m_epsilon_cor = " << m_epsilon_cor << std::endl;
    std::cout << " m_n_threads = " << m_n_threads << std::endl;
    std::cout << " m_geometry_file = " << m_geometry_file << std::endl;
    std::cout << " m_geometry = " << m_geometry << std::endl;
    std::cout << " m_vtk_file = " << m_vtk_file << std::endl;
    std::cout << " m_vtk_resolution = " << m_vtk_resolution << std::endl;
    std::cout << " m_scalnames = " << m_scalnames << std::endl;
    std::cout << " m_vecnames = " << m_vecnames << std::endl;
    std::cout << " m_n_regions = " << m_n_regions << std::endl;
    
    std::cout << " m_mat_ids = " << std::endl;
    int n_data = m_mat_ids.size();
    for (int i = 0; i < n_data; i++)
    {
        std::cout << " region material id = " << m_mat_ids[i].first << std::endl;
        int n_bc = m_mat_ids[i].second.size();
        std::cout << " bc material ids = " << m_mat_ids[i].first << std::endl;
        for (int j = 0; j <n_bc; j++)
        {
            std::cout << " " << m_mat_ids[i].second [j];
        }
        std::cout << std::endl;
    }
    
    std::cout << " m_mat_props = " << std::endl;
    n_data = m_mat_props.size();
    for (int i = 0; i < n_data; i++) {
        std::cout << " region number = " << i << std::endl;
        int n_bc = m_mat_props[i].size();
        for (int j = 0; j <n_bc; j++) {
            std::cout << " " << m_mat_props[i][j];
        }
        std::cout << std::endl;
    }
    
    std::cout << " m_is_initial_state_Q = " <<m_is_initial_state_Q << std::endl;
    std::cout << " m_is_current_state_Q = " << m_is_current_state_Q << std::endl;
    std::cout << " m_time = " << m_time << std::endl;
    std::cout << " ---------------------- " << std::endl;
    std::cout << std::endl;
    
}


/** @brief read the geometry */
void TPZSimulationData::ReadGeometry()
{
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    m_geometry = Geometry.GeometricGmshMesh(m_geometry_file);
    
#ifdef PZDEBUG
    if (!m_geometry)
    {
        std::cout << "The mesh was not generated." << std::endl;
        DebugStop();
    }
#endif
    
}

/** @brief print the geometry */
void TPZSimulationData::PrintGeometry()
{
    
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name  << "geometry" << ".txt";
    vtk_name   << "geometry"  << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    m_geometry->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(m_geometry, vtkfile, true);

#ifdef PZDEBUG
    TPZCheckGeom checker(m_geometry);
    checker.CheckUniqueId();
    if(checker.PerformCheck())
    {
        DebugStop();
    }
#endif
    
}


/** @brief applying the boundary conditions */
void TPZSimulationData::LoadBoundaryConditions()
{
 
    std::pair<std::string,std::pair<int,std::vector<std::string> > > chunk;
    if (m_dimesion!=3)
    {
        
        // 2D conditions
        
        // Dirichlet for elasticity and Dirichlet for diffusion
        chunk.first = "Du_Dp"; // name
        chunk.second.first = 0; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in x_direction and Dirichlet for diffusion
        chunk.first = "Dux_Dp"; // name
        chunk.second.first = 1; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in y_direction and Dirichlet for diffusion
        chunk.first = "Duy_Dp"; // name
        chunk.second.first = 2; // index
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Neumann for elasticity and Dirichlet for diffusion
        chunk.first = "Nt_Dp"; // name
        chunk.second.first = 3; // index
        chunk.second.second.push_back("tx");
        chunk.second.second.push_back("ty");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Neumann for elasticity and Dirichlet for diffusion (Wellbore Boundary)
        chunk.first = "Ntn_Dp"; // name
        chunk.second.first = 4; // index
        chunk.second.second.push_back("tn");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity and Neumann for diffusion
        chunk.first = "Du_Nq"; // name
        chunk.second.first = 5; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in x_direction and Neumann for diffusion
        chunk.first = "Dux_Nq"; // name
        chunk.second.first = 6; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in y_direction and Neumann for diffusion
        chunk.first = "Duy_Nq"; // name
        chunk.second.first = 7; // index
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Neumann for elasticity and Neumann for diffusion
        chunk.first = "Nt_Nq"; // name
        chunk.second.first = 8; // index
        chunk.second.second.push_back("tx");
        chunk.second.second.push_back("ty");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Neumann for elasticity and Neumann for diffusion (Wellbore Boundary)
        chunk.first = "Ntn_Nq"; // name
        chunk.second.first = 9; // index
        chunk.second.second.push_back("tn");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);       

        // Dirichlet for elasticity in y_direction and Dirichlet for diffusion (time dependent)
        chunk.first = "Duy_time_Dp"; // name
        chunk.second.first = 10; // index
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Neumann for elasticity and Dirichlet for diffusion (time dependent)
        chunk.first = "Ntn_time_Dp"; // name
        chunk.second.first = 11; // index
        chunk.second.second.push_back("tx");
        chunk.second.second.push_back("ty");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        
    }
    else
    {
        // 3D conditions
        
        // Dirichlet for elasticity and Dirichlet for diffusion
        chunk.first = "Du_Dp"; // name
        chunk.second.first = 0; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("uz");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in x & y direction and Dirichlet for diffusion
        chunk.first = "Duxy_Dp"; // name
        chunk.second.first = 1; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in x & z direction and Dirichlet for diffusion
        chunk.first = "Duxz_Dp"; // name
        chunk.second.first = 2; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("uz");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in y & z direction and Dirichlet for diffusion
        chunk.first = "Duyz_Dp"; // name
        chunk.second.first = 3; // index
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("uz");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in x_direction and Dirichlet for diffusion
        chunk.first = "Dux_Dp"; // name
        chunk.second.first = 4; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in y_direction and Dirichlet for diffusion
        chunk.first = "Duy_Dp"; // name
        chunk.second.first = 5; // index
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in z_direction and Dirichlet for diffusion
        chunk.first = "Duz_Dp"; // name
        chunk.second.first = 6; // index
        chunk.second.second.push_back("uz");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        
        // Neumann for elasticity and Dirichlet for diffusion
        chunk.first = "Nt_Dp"; // name
        chunk.second.first = 7; // index
        chunk.second.second.push_back("tx");
        chunk.second.second.push_back("ty");
        chunk.second.second.push_back("tz");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Neumann for elasticity and Dirichlet for diffusion (Wellbore Boundary)
        chunk.first = "Ntn_Dp"; // name
        chunk.second.first = 8; // index
        chunk.second.second.push_back("tn");
        chunk.second.second.push_back("p");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity and Neumann for diffusion
        chunk.first = "Du_Nq"; // name
        chunk.second.first = 9; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("uz");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in x & y direction and Dirichlet for diffusion
        chunk.first = "Duxy_Nq"; // name
        chunk.second.first = 10; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in x & z direction and Dirichlet for diffusion
        chunk.first = "Duxz_Nq"; // name
        chunk.second.first = 11; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("uz");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in y & z direction and Dirichlet for diffusion
        chunk.first = "Duyz_Nq"; // name
        chunk.second.first = 12; // index
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("uz");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in x_direction and Neumann for diffusion
        chunk.first = "Dux_Nq"; // name
        chunk.second.first = 13; // index
        chunk.second.second.push_back("ux");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in y_direction and Neumann for diffusion
        chunk.first = "Duy_Nq"; // name
        chunk.second.first = 14; // index
        chunk.second.second.push_back("uy");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Dirichlet for elasticity in z_direction and Neumann for diffusion
        chunk.first = "Duz_Nq"; // name
        chunk.second.first = 15; // index
        chunk.second.second.push_back("uz");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
        
        // Neumann for elasticity and Neumann for diffusion
        chunk.first = "Nt_Nq"; // name
        chunk.second.first = 16; // index
        chunk.second.second.push_back("tx");
        chunk.second.second.push_back("ty");
        chunk.second.second.push_back("tz");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);

        // Neumann for elasticity and Neumann for diffusion (Wellbore Boundary)
        chunk.first = "Ntn_Nq"; // name
        chunk.second.first = 17; // index
        chunk.second.second.push_back("tn");
        chunk.second.second.push_back("qn");
        m_condition_type_to_index_value_names.insert(chunk);
        chunk.second.second.resize(0);
    }
    

    return;
}
