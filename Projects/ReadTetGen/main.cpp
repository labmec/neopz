//$Id: main.cpp,v 1.1 2006-03-04 15:39:02 tiago Exp $

/**
 * Validation test of TPZReadtetGen class wich implements the interface between tetgen and PZ.
 * March 03, 2006
 */

#include "pzgmesh.h"
#include "pzreadtetgen.h"
int main(){

  TPZReadTetGen read;
  TPZGeoMesh * gmesh = read.Process("malha10.1.node", "malha10.1.face", "malha10.1.ele");
//  gmesh->Print(std::cout);


}
