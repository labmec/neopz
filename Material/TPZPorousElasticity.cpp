//
//  TPZPorousElasticity.cpp
//  pz
//
//  Created by Omar Durán on 4/13/18.
//

#include "TPZPorousElasticity.h"


int TPZPorousElasticity::ClassId() const {
    return Hash("TPZPorousElasticity");
}

TPZPorousElasticity::TPZPorousElasticity(){
    DebugStop();
}

TPZPorousElasticity::TPZPorousElasticity(int matid){
    DebugStop();
}

TPZPorousElasticity & TPZPorousElasticity::operator=(const TPZPorousElasticity &copy){
    DebugStop();
}


TPZPorousElasticity::~TPZPorousElasticity(){
    DebugStop();
}
