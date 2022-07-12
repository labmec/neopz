#ifndef _TPZVTKGENERATOR_H_
#define _TPZVTKGENERATOR_H_

#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzeltype.h"
#include <iostream>
#include <map>

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
  TPZCompMesh* fCMesh = nullptr;

  std::string fFilename = "";
  int fSubdivision{0};
  TPZVec<TPZAutoPointer<TPZVTKField>> fFields;

  TPZVec<TPZManVector<REAL,3>> fPoints;
  TPZVec<std::array<int,TPZVTK::MAX_PTS+2>> fCells;//max 

  int fOutputCount = 0;
  TPZVec<REAL> fTimes = {0};
  TPZAutoPointer<std::ofstream> fFileout{nullptr};
  std::string fLastOutputName = "";

  std::map<MElementType,TPZVec<TPZManVector<REAL,3>>> fRefVertices;
  std::map<MElementType,TPZVec<std::array<int,TPZVTK::MAX_PTS+2>>> fRefEls;
  //! (el,first pos) list of elements and the initial position of their solution in the field vec
  TPZVec<std::pair<TPZCompEl*,int>> fElementVec;

  template<class TOPOL>
  void FillReferenceEl(TPZVec<TPZManVector<REAL,3>> &ref_coords,
                       TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems);

  void FillRefEls(const int meshdim);
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
  //! Compute all post-processing points and create list of valid elements
  void ComputePoints();
public:
  /**
     @brief Creates instance for generating .vtk results for a given mesh
     @param[in] cmesh Computational mesh
     @param[in] fields names of fields to be post-processed
     @param[in] filename filename without extension
     @param[in] vtkres resolution of vtk post-processing (number of el subdivision)
  */
  TPZVTKGenerator(TPZCompMesh* cmesh,
                  const TPZVec<std::string> fields,
                  std::string filename,
                  int vtkres);
  /**
     @brief Creates instance for generating .vtk results for a given mesh
     @param[in] cmesh Computational mesh
     @param[in] fields names of fields to be post-processed
     @param[in] filename filename without extension
     @param[in] vtkres resolution of vtk post-processing (number of el subdivision)
  */
  TPZVTKGenerator(TPZAutoPointer<TPZCompMesh> cmesh,
                  const TPZVec<std::string> fields,
                  std::string filename,
                  int vtkres);
  //!Generates .vtk file for a given current solution of the mesh
  void Do(REAL time = -1);

  /** @brief Resets post-processing data structure. Call this function if
   the mesh has had changes between Do() calls (refinement, etc)*/
  void ResetArrays();
};

#endif /* _TPZVTKGENERATOR_H_ */