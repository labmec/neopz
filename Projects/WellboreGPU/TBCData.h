//
//  TBCData.h
//  IntegrationPointExperiments
//
//  Created by Omar Durán on 3/26/19.
//

#ifndef TBCData_h
#define TBCData_h

#include <stdio.h>
#include <vector>
#include "pzreal.h"

class TBCData {

protected:
    
    /// material identifier
    int m_id = -1;
    
    /// bc type - 0 -> Dirichlet and 1 -> Neumann
    int m_type = -1;
    
    /// The boundary data
    std::vector<REAL> m_value;

    REAL m_initial_value = -1;

public:
    
    /// Default constructor
    TBCData();
    
    /// Copy constructor
    TBCData(const TBCData &  other);
    
    /// Assignmet constructor
    TBCData & operator=(const TBCData &  other);
    
    /// Default destructor
    ~TBCData();

    void SetId (int id);

    int Id ();

    void SetType (int type);

    int Type ();

    void SetValue (std::vector<REAL> m_value);

    std::vector<REAL> Value ();

    void SetInitialValue (REAL initv);

    REAL InitialValue();
};

#endif /* TBCData_h */
