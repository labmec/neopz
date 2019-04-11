//
//  TElastoPlasticData.h
//  IntegrationPointExperiments
//
//  Created by Omar Durán on 3/26/19.
//

#ifndef TElastoPlasticData_h
#define TElastoPlasticData_h

#include <stdio.h>
#include "TBCData.h"
#include "TPZElasticResponse.h"

class TElastoPlasticData {
    
protected:
    
    /// Stands for boundary data
    std::vector<TBCData> m_gamma_data;
    
    /// material identifier
    int m_id = -1;
    
    TPZElasticResponse m_LER;
    
    /// Friction angle
    REAL m_MC_phi = -1;
    
    /// Cohesion
    REAL m_MC_c = -1;

public:
    
    /// Default constructor
    TElastoPlasticData();
    
    /// Copy constructor
    TElastoPlasticData(const TElastoPlasticData &  other);
    
    /// Assignmet constructor
    TElastoPlasticData & operator=(const TElastoPlasticData &  other);
    
    /// Default destructor
    ~TElastoPlasticData();

    void SetBoundaryData (std::vector<TBCData> bcdata);

    std::vector<TBCData> BoundaryData ();

    void SetId (int id);

    int Id ();

    void SetMaterialParameters (TPZElasticResponse ER, REAL phi, REAL c);

    TPZElasticResponse ElasticResponse ();

    REAL FrictionAngle();

    REAL Cohesion();


};

#endif /* TElastoPlasticData_h */
