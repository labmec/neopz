/***************************************************************************
                          pzidentifyrefpattern.cpp  -  description
                             -------------------
    begin                : Mon Mar 8 2004
    copyright            : (C) 2004 by cesar
    email                : cesar@becks.fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "pzidentifyrefpattern.h"
#include "pzgeoel.h"
#include <pzeltype.h>
#include <TPZRefPattern.h>

#include <pzrefpoint.h>
#include <TPZRefLinear.h>
#include <pzreftriangle.h>
#include <pzrefquad.h>
#include <pzrefpyram.h>
#include <pzrefprism.h>
#include <pzreftetrahedra.h>
#include <TPZRefCube.h>

#include <iterator>

TPZIdentifyRefPattern::TPZIdentifyRefPattern(string &path){
  fPath = path;  
}
TPZIdentifyRefPattern::~TPZIdentifyRefPattern(){
}
/** Returns the refinement pattern that generates the given refinement */
TPZRefPattern * TPZIdentifyRefPattern::GetRefPattern (TPZGeoEl *father, TPZVec<TPZGeoEl *> subelem){
  int eltype =  father->Type();
  int nelem = subelem.NElements();
  TPZRefPattern *rp;
  switch (eltype) {
    case (EPoint) : {
      if (nelem != 1){
        PZError << "TPZIdentifyRefPattern::GetRefPattern ERROR : point partition detected!" << endl;
        exit (-1);
      }
      return 0;// return refpattern for point!!
    }
    case (EOned) :{
      if (nelem != 2){
        PZError << "TPZIdentifyRefPattern::GetRefPattern ERROR : wrong linear partition detected!" << endl;
        exit (-1);
      }
      rp = new TPZRefPattern ("/home/pos/cesar/RefPattern/Unif_Linear.rpt");
    }
    default : {
      //return a uniform refinement pattern
      if (nelem == UniformSubElem(eltype)) {
        rp = GetUniform(eltype);
      }else{
        //identify the side refinement pattern
        int side = IdentifySide(father,subelem);
        rp = GetSideRefPattern(eltype,side);
      }
    }
  }
  TPZRefPattern *mesh_rp = father->Mesh()->GetRefPattern (eltype,rp->GetName());
  //If the refinement pattern is already defined delete the created refinement pattern
  if (mesh_rp) {
    delete rp;
    return mesh_rp;
  }
  //Insert a new refinement pattern into mesh
  father->Mesh()->InsertRefPattern(rp);
  return rp;
}

/** Identify the side of the refinement pattern */
int TPZIdentifyRefPattern::IdentifySide(TPZGeoEl *father, TPZVec<TPZGeoEl *> subelem){
  set< TSide > Father;
  set< TSide > Sons;
  set< TSide > Result;
  int iside,isub;

  for (iside=0;iside<father->NSides();iside++){
    if (father->SideDimension(iside) != 1) continue;
    TPZGeoElSide gelside (father,iside);
    TSide sidefat(gelside);
    Father.insert(sidefat);
  }
  for (isub=0;isub<subelem.NElements();isub++){
    for (iside=0;iside<subelem[isub]->NSides();iside++){
      if (subelem[isub]->SideDimension(iside) != 1) continue;
      TPZGeoElSide gelside (subelem[isub],iside);
      TSide sideson(gelside);
      Sons.insert(sideson);
    }
  }

  set_difference(Father.begin(),Father.end(),Sons.begin(),Sons.end(),inserter(Result,Result.begin()));
  if (Result.size() != 1) return -1;

  TSide side_result = *Result.begin();
  int side  = side_result.fSide;
  return side;
}

using namespace pzrefine;

int TPZIdentifyRefPattern::UniformSubElem( int eltype){
  switch (eltype) {
    case (EPoint)         : return TPZRefPoint::NSubEl;
    case (EOned)          : return TPZRefLinear::NSubEl;
    case (ETriangle)      : return TPZRefTriangle::NSubEl;
    case (EQuadrilateral) : return TPZRefQuad::NSubEl;
    case (ETetraedro)     : return TPZRefTetrahedra::NSubEl;
    case (EPiramide)      : return TPZRefPyramid::NSubEl;
    case (EPrisma)        : return TPZRefPrism::NSubEl;
    case (ECube)          : return TPZRefCube::NSubEl;
  }
  PZError << "TPZIdentifyRefPattern::UniformSubElem ERROR unknown eltype : " << eltype << endl;
  return (-1);
}

TPZRefPattern * TPZIdentifyRefPattern::GetUniform(int eltype){
  string fullfilename = fPath;
  fullfilename += "/";
  switch (eltype) {
    case (EPoint) : {
      fullfilename += "Point_";
      break;
    }
    case (EOned) : {
      fullfilename += "Linear_";
      break;
    }
    case (ETriangle) : {
      fullfilename += "Triang_";
      break;
    }
    case (EQuadrilateral) : {
      fullfilename += "Quad_";
      break;
    }
    case (ETetraedro) : {
      fullfilename += "Tetra_";
      break;
    }
    case (EPiramide)  : {
      fullfilename += "Piram_";
      break;
    }
    case (EPrisma) : {
      fullfilename += "Prism_";
      break;
    }
    case (ECube) : {
      fullfilename += "Hexa_";
      break;
    }
    default:{
      PZError << "TPZIdentifyRefPattern::GetUniform ERROR unknown eltype : " << eltype << endl;
      return (0);
    }
  }
  fullfilename += "Unif.rpt";
  TPZRefPattern *rp = new TPZRefPattern(fullfilename);  
  return rp;
}

TPZRefPattern * TPZIdentifyRefPattern::GetSideRefPattern(int eltype,int side){
  string fullfilename = fPath;
  fullfilename += "/";
  switch (eltype) {
    case (EPoint) : {
      fullfilename += "Point_";
      break;
    }
    case (EOned) : {
      fullfilename += "Linear_";
      break;
    }
    case (ETriangle) : {
      fullfilename += "Triang_";
      break;
    }
    case (EQuadrilateral) : {
      fullfilename += "Quad_";
      break;
    }
    case (ETetraedro) : {
      fullfilename += "Tetra_";
      break;
    }
    case (EPiramide)  : {
      fullfilename += "Piram_";
      break;
    }
    case (EPrisma) : {
      fullfilename += "Prism_";
      break;
    }
    case (ECube) : {
      fullfilename += "Hexa_";
      break;
    }
    default:{
      PZError << "TPZIdentifyRefPattern::GetUniform ERROR unknown eltype : " << eltype << endl;
      return (0);
    }
  }
  fullfilename += "Unif.rpt";
  TPZRefPattern *rp = new TPZRefPattern(fullfilename);
  return rp;
}


