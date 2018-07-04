//
//  TNRFElasticMemory.cpp
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#include "TNRFElasticMemory.h"


TNRFElasticMemory::TNRFElasticMemory(){
    
}

TNRFElasticMemory::TNRFElasticMemory(const TNRFElasticMemory & other){
    
}

const TNRFElasticMemory & TNRFElasticMemory::operator=(const TNRFElasticMemory & other){
    
    if (this != & other) // prevent self-assignment
    {

    }
    return *this;
}

TNRFElasticMemory::~TNRFElasticMemory(){
    
}

void TNRFElasticMemory::Print(std::ostream &out) const {
    
}

int TNRFElasticMemory::ClassId() const{;
    return Hash("TNRFElasticMemory");
}
