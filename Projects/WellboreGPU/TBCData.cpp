//
//  TBCData.cpp
//  IntegrationPointExperiments
//
//  Created by Omar Durán on 3/26/19.
//

#include "TBCData.h"

TBCData::TBCData(){
    m_value.resize(0);
}

TBCData::TBCData(const TBCData &  other){
    
    m_id            = other.m_id;
    m_type          = other.m_type;
    m_value         = other.m_value;
    m_initial_value = other.m_initial_value;
}

TBCData & TBCData::operator=(const TBCData &  other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_id            = other.m_id;
    m_type          = other.m_type;
    m_value         = other.m_value;
    m_initial_value = other.m_initial_value;
    
    return *this;
}

TBCData::~TBCData(){
    
}

void TBCData::SetId (int id) {
    m_id = id;
}

int TBCData::Id () {
    return m_id;
}

void TBCData::SetType (int type) {
    m_type = type;
}

int TBCData::Type () {
    return m_type;
}

void TBCData::SetValue (std::vector<REAL> value) {
    m_value = value;
}

std::vector<REAL> TBCData::Value () {
    return m_value;
}

void TBCData::SetInitialValue (REAL initv) {
    m_initial_value = initv;
}

REAL TBCData::InitialValue() {
    return m_initial_value;
}
