#ifndef _TPZVTKGENERATOR_H_
#define _TPZVTKGENERATOR_H_

#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzeltype.h"
#include <iostream>

class TPZCompMesh;

/*********************************************************************/
/* File:   TPZVTKGenerator.h                                         */
/* Author: Francisco Orlandini                                       */
/* Date:   5. July 2022                                              */
/* Adapted from: NGSolve's vtkoutput.hpp                             */
/*********************************************************************/

namespace TPZVTK{//useful constants for compile time
  //hexahedron: nnodes
  static constexpr int MAX_PTS{8};
  //pyramid: nsub els
  static constexpr int MAX_SUBEL{10};
  static int CellType(const MElementType el){
    switch (el){
    case EPoint:
      return 1;
    case EOned:
      return 3;
    case ETriangle:
      return 5;
    case EQuadrilateral:
      return 9;
    case ETetraedro:
      return 10;
    case EPiramide:
      return 14;
    case EPrisma:
      return 13;
    case ECube:
      return 12;
    default:
      std::cout << "TPZVTKGenerator Element Type "
                << MElementType_Name(el) << " not supported!"
                << std::endl;
      DebugStop();
      unreachable();
    }
  }
};

class TPZVTKField : public TPZVec<STATE>
{
public:
  enum class Type{scal,vec,tens};
  TPZVTKField() { ; };
  TPZVTKField(Type t, int id, std::string name);
  void SetType(Type t) {fType = t;}
  Type GetType() const {return fType;}
  std::string TypeName() const{
    switch (fType){
    case Type::scal: return "SCALARS";
    case Type::vec: return "VECTORS";
    case Type::tens: return "TENSORS";
    }
  }
  
  int Dimension() const {
    switch (fType){
    case Type::scal: return 1;
    case Type::vec: return 3;
    case Type::tens: return 9;
    }
  }
  void SetName(std::string name) { fName = name; }
  std::string Name() const{ return fName; }
  void SetId(int id) {fId = id;}
  int Id() {return fId;}
private:
  Type fType = Type::scal;
  int fId = -1;
  std::string fName = "none";
};

class TPZVTKGenerator{
protected:
  TPZAutoPointer<TPZCompMesh> fCMesh = nullptr;

  std::string fFilename = "";
  int fSubdivision{0};
  TPZVec<TPZAutoPointer<TPZVTKField>> fFields;
  TPZVec<std::array<REAL,3>> fPoints;
  int fNPtsPerPt = fPoints.size() * 5;
  TPZVec<std::array<int,TPZVTK::MAX_PTS+2>> fCells;//max 

  int fOutputCount = 0;
  TPZVec<REAL> fTimes = {0};
  TPZAutoPointer<std::ofstream> fFileout{nullptr};
  std::string fLastOutputName = "";

  //! Resets fFields, fPoints, fCells
  void ResetArrays();

  void FillReferenceLine(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems);
  void FillReferenceTrig(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems);
  void FillReferenceQuad(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems);
  void FillReferenceTet(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems);
  void FillReferenceHex(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems);
  void FillReferencePrism(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems);
  
  //! Print all points in VTK Legacy format
  void PrintPointsLegacy();
  //! Print all cells in VTK Legacy format
  void PrintCellsLegacy();
  //! Print all cell types in VTK Legacy format
  void PrintCellTypesLegacy(int meshdim);
  //! Print all fields in VTK Legacy format
  void PrintFieldDataLegacy();
  //! Check if element should be processed
  bool IsValidEl(TPZCompEl *el);
public:
  //! Creates instance for generating .vtk results for a given mesh
  TPZVTKGenerator(TPZAutoPointer<TPZCompMesh> cmesh,
                  const TPZVec<std::string> fields,
                  std::string filename,
                  int vtkres);
  //!Generates .vtk file for a given current solution of the mesh
  void Do(REAL time = -1);
};

#endif /* _TPZVTKGENERATOR_H_ */