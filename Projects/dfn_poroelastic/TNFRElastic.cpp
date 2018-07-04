//
//  TNFRElastic.cpp
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#include "TNFRElastic.h"
#include "TNRFElasticMemory.h"
#include "pzbndcond.h"
#include "TPZMatWithMem.h"

// Constructor
TNFRElastic::TNFRElastic(){
    
    m_geometry  = nullptr;
    m_cmesh     = nullptr;
    m_Operator  = nullptr;
    m_boundary_data.clear();
}

// Descructor
TNFRElastic::~TNFRElastic(){
    
}


// setup the internal members
void TNFRElastic::BuildOperator(TPZGeoMesh *  geometry, int p_order){
    
    m_geometry = geometry;
    m_cmesh = new TPZCompMesh(m_geometry);
    m_cmesh->SetDefaultOrder(p_order);
    m_cmesh->SetDimModel(m_geometry->Dimension());
    m_cmesh->SetAllCreateFunctionsContinuous();
    
    // Inserting a single material and related boundary conditions
    int rock_id = 1;
    TNFRElasticMaterial * rock_material = new TNFRElasticMaterial(rock_id);
    m_cmesh->InsertMaterialObject(rock_material);

    unsigned int n_bc = m_boundary_data.size();
    TPZFMatrix<REAL> val1(1,1,0.0),val2(1,1,0.0);
    for (unsigned int i = 0; i < n_bc; i++) {
        TNFRBoundaryDescription bc = m_boundary_data[i];
        unsigned int n_values = bc.GetBCValues().size();
        val2.Resize(n_values, 1);
        for (unsigned int j = 0; j < n_values; j++) {
            val2(i,0) = bc.GetBCValues()[j];
        }
        
        TPZMatWithMem<TNRFElasticMemory,TPZBndCond> * boundary_c = new TPZMatWithMem<TNRFElasticMemory,TPZBndCond>;
        boundary_c->SetNumLoadCases(1);
        boundary_c->SetMaterial(rock_material);
        boundary_c->SetId(bc.GetBCId());
        boundary_c->SetType(bc.GetBCType());
        boundary_c->SetValues(val1, val2);
        m_cmesh->InsertMaterialObject(boundary_c);
    }
    
    m_cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshDeformation.txt");
    m_cmesh->Print(out);
#endif
    
    
}
