#include "TPZVTKGenerator.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "TPZMaterial.h"
#include "TPZSimpleTimer.h"
#include "pzinterpolationspace.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMaterialDataT.h"
#include "pzmultiphysicselement.h"

/*********************************************************************/
/* File:   TPZVTKGenerator.cpp                                       */
/* Author: Francisco Orlandini                                       */
/* Date:   5. July 2022                                              */
/* Adapted from: NGSolve's vtkoutput.hpp                             */
/*********************************************************************/


template<class TVar>
class TPZPostProcEl{
public:
  TPZPostProcEl(TPZCompEl *cel);
  void InitData();
  void ComputeRequiredData(TPZVec<REAL> &qsi);
  void Solution(const TPZVec<REAL> &qsi, const int id, TPZVec<TVar> &sol);
private:
  bool fIsMultiphysics{false};
  TPZCompEl *fCel{nullptr};
  TPZMaterialDataT<TVar> fMatdata;
  TPZManVector<TPZMaterialDataT<TVar>,10> fDatavec; 
};

template<class TVar>
void ComputeFieldAtEl(TPZCompEl *cel,
                      const TPZVec<std::array<REAL,3>> &ref_vertices,
                      TPZVec<TPZAutoPointer<TPZVTKField>>& fields){
  TPZManVector<TVar,9> sol;

  const auto celdim = cel->Dimension();

  TPZPostProcEl<TVar> graphel(cel);
  graphel.InitData();
  TPZManVector<REAL,3> qsi(celdim,0);
  for (const auto &ip : ref_vertices){
    //copy to tpzvec with appropriate size
    for(int ix = 0; ix < celdim; ix++){qsi[ix] = ip[ix];}

    graphel.ComputeRequiredData(qsi);
    for (int i = 0; i < fields.size(); i++){
      auto &field = *(fields[i]);
      auto fdim = field.Dimension();
      sol.Resize(fdim);
      graphel.Solution(qsi, field.Id(), sol);
      const auto sz = sol.size();
      if constexpr (std::is_same_v<TVar,CSTATE>){
        for (int d = 0; d < sz; ++d){
          AppendToVec(field, std::real(sol[d]));
        }
        for (int d = sz; d < fdim; ++d){
          AppendToVec(field, 0.0);
        }
      }else{
        for (int d = 0; d < fdim; ++d){
          AppendToVec(field, sol[d]);
        }
        for (int d = sz; d < fdim; ++d){
          AppendToVec(field, 0.0);
        }
      }
    }
  }
}


TPZVTKField::TPZVTKField(TPZVTKField::Type type, int id, std::string aname) :
  fType( type), fId(id), fName(aname) { ; }


TPZVTKGenerator::TPZVTKGenerator(TPZAutoPointer<TPZCompMesh> cmesh,
                                    const TPZVec<std::string> fields,
                                    std::string filename,
                                    int vtkres)
  : fCMesh(cmesh), fFilename(filename), fSubdivision(vtkres)
{

  //let us init the field ids
  const int nvars = fields.size();

  fFields.resize(nvars);
  auto * matp = (fCMesh->MaterialVec().begin())->second;

  for(int i = 0; i < nvars; i++){
    const auto &name = fields[i];
    const auto index = matp->VariableIndex(name);
    const auto dim = matp->NSolutionVariables(index);
    const TPZVTKField::Type type = [dim](){
      switch(dim){
      case 1:
        return TPZVTKField::Type::scal;
      case 2:
      case 3:
        return TPZVTKField::Type::vec;
      case 9:
        return TPZVTKField::Type::tens;
      default:
        PZError<<"TPZVTKGenerator:\n"
               <<"field dim "<<dim
               <<" not supported. It should be either "
               <<"\t1 (scalar)\n"
               <<"\t3 (vector)\n"
               <<"\t9 (tensor)\n";
        DebugStop();
        return TPZVTKField::Type::scal;
      }
    }();
    fFields[i] = new TPZVTKField(type,index,name);
  }
}

/// Empty all fields, points and cells
void TPZVTKGenerator::ResetArrays()
{
  fPoints.Resize(0);
  fCells.Resize(0);
  for(auto field : fFields){
    field->Resize(0);
  }
}


void TPZVTKGenerator::FillReferenceLine(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems)
{
  const int r = 1 << fSubdivision;

  const REAL h = 2.0 / r;
  constexpr REAL init = -1;
  for (int i = 0; i <= r; ++i){
    AppendToVec(ref_coords,std::array<REAL,3>{init + i * h, 0});
    AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(EOned), 2, i, i+1});
  }
}


void TPZVTKGenerator::FillReferenceTrig(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems)
{
  if (fSubdivision == 0){
    AppendToVec(ref_coords,std::array<REAL,3>{0.0, 0.0, 0.0});
    AppendToVec(ref_coords,std::array<REAL,3>{1.0, 0.0, 0.0});
    AppendToVec(ref_coords,std::array<REAL,3>{0.0, 1.0, 0.0});
    AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETriangle),3, 0, 1, 2});
  }
  else{
    const int r = 1 << fSubdivision;
    const int s = r + 1;

    const REAL h = 1.0 / r;

    int pidx = 0;
    for (int i = 0; i <= r; ++i)
      for (int j = 0; i + j <= r; ++j){
        AppendToVec(ref_coords,std::array<REAL,3>{j * h, i * h});
      }

    pidx = 0;

    for (int i = 0; i <= r; ++i)
      for (int j = 0; i + j <= r; ++j, pidx++)
      {
        // int pidx_curr = pidx;
        if (i + j == r)
          continue;
        int pidx_incr_i = pidx + 1;
        int pidx_incr_j = pidx + s - i;

        AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETriangle), 3, pidx, pidx_incr_i, pidx_incr_j});
        int pidx_incr_ij = pidx_incr_j + 1;

        if (i + j + 1 < r){
          AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETriangle), 3, pidx_incr_i, pidx_incr_ij, pidx_incr_j});
        }
      }
  }
}
/// Fill principil lattices (points and connections on subdivided reference simplex) in 2D
void TPZVTKGenerator::FillReferenceQuad(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems)
{
  if (fSubdivision == 0)
    {
      AppendToVec(ref_coords,std::array<REAL,3>{-1.0, -1.0, 0.0});
      AppendToVec(ref_coords,std::array<REAL,3>{ 1.0, -1.0, 0.0});
      AppendToVec(ref_coords,std::array<REAL,3>{ 1.0,  1.0, 0.0});
      AppendToVec(ref_coords,std::array<REAL,3>{-1.0,  1.0, 0.0});
      AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS+2>{TPZVTK::CellType(EQuadrilateral), 4, 0, 1, 2, 3});
    }
  else
    {
      const int r = 1 << fSubdivision;
      // const int s = r + 1;

      const REAL h = 2.0 / r;
      constexpr REAL init{-1};
      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; j <= r; ++j)
          {
            AppendToVec(ref_coords,std::array<REAL,3>{init + j * h, init + i * h});
          }

      for (int i = 0; i < r; ++i)
        {
          int incr_i = r + 1;
          pidx = i * incr_i;
          for (int j = 0; j < r; ++j, pidx++)
            {
              AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(EQuadrilateral),4, pidx, pidx + 1, pidx + incr_i + 1, pidx + incr_i});
            }
        }
    }
}

/// Fill principil lattices (points and connections on subdivided reference simplex) in 3D
void TPZVTKGenerator::FillReferenceTet(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems)
{
  if (fSubdivision == 0)
    {
      AppendToVec(ref_coords,std::array<REAL,3>{0.0, 0.0, 0.0});
      AppendToVec(ref_coords,std::array<REAL,3>{1.0, 0.0, 0.0});
      AppendToVec(ref_coords,std::array<REAL,3>{0.0, 1.0, 0.0});
      AppendToVec(ref_coords,std::array<REAL,3>{0.0, 0.0, 1.0});
      AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETetraedro),4, 0, 1, 2, 3});
    }
  else
    {
      const int r = 1 << fSubdivision;
      const int s = r + 1;

      const REAL h = 1.0 / r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i + j <= r; ++j)
          for (int k = 0; i + j + k <= r; ++k)
            {
              AppendToVec(ref_coords,std::array<REAL,3>{i * h, j * h, k * h});
            }

      for (int i = 0; i <= r; ++i)
        for (int j = 0; i + j <= r; ++j)
          for (int k = 0; i + j + k <= r; ++k, pidx++)
            {
              if (i + j + k == r)
                continue;
              // int pidx_curr = pidx;
              int pidx_incr_k = pidx + 1;
              int pidx_incr_j = pidx + s - i - j;
              int pidx_incr_i = pidx + (s - i) * (s + 1 - i) / 2 - j;

              int pidx_incr_kj = pidx_incr_j + 1;

              int pidx_incr_ij = pidx + (s - i) * (s + 1 - i) / 2 - j + s - (i + 1) - j;
              int pidx_incr_ki = pidx + (s - i) * (s + 1 - i) / 2 - j + 1;
              int pidx_incr_kij = pidx + (s - i) * (s + 1 - i) / 2 - j + s - (i + 1) - j + 1;

              AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETetraedro), 4, pidx, pidx_incr_k, pidx_incr_j, pidx_incr_i});
              if (i + j + k + 1 == r)
                continue;

              AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETetraedro), 4, pidx_incr_k, pidx_incr_kj, pidx_incr_j, pidx_incr_i});
              AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETetraedro), 4, pidx_incr_k, pidx_incr_kj, pidx_incr_ki, pidx_incr_i});

              AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETetraedro), 4, pidx_incr_j, pidx_incr_i, pidx_incr_kj, pidx_incr_ij});
              AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETetraedro), 4, pidx_incr_i, pidx_incr_kj, pidx_incr_ij, pidx_incr_ki});

              if (i + j + k + 2 != r)
                AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ETetraedro), 4, pidx_incr_kj, pidx_incr_ij, pidx_incr_ki, pidx_incr_kij});
            }
    }
}

/// Fill principil lattices (points and connections on subdivided reference hexahedron) in 3D
void TPZVTKGenerator::FillReferenceHex(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems)
{
  if (fSubdivision == 0)
    {
      AppendToVec(ref_coords,std::array<REAL,3>{-1.0, -1.0, -1.0});
      AppendToVec(ref_coords,std::array<REAL,3>{ 1.0, -1.0, -1.0});
      AppendToVec(ref_coords,std::array<REAL,3>{ 1.0,  1.0, -1.0});
      AppendToVec(ref_coords,std::array<REAL,3>{-1.0,  1.0, -1.0});
      AppendToVec(ref_coords,std::array<REAL,3>{-1.0, -1.0,  1.0});
      AppendToVec(ref_coords,std::array<REAL,3>{ 1.0, -1.0,  1.0});
      AppendToVec(ref_coords,std::array<REAL,3>{ 1.0,  1.0,  1.0});
      AppendToVec(ref_coords,std::array<REAL,3>{-1.0,  1.0,  1.0});
      AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(ECube),8,0,1,2,3,4,5,6,7});
    }
  else
    {
      const int r = 1 << fSubdivision;
      // const int s = r + 1;

      const REAL h = 2.0 / r;
      constexpr REAL init{-1};
      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; j <= r; ++j)
          for (int k = 0; k <= r; ++k)
            {
              AppendToVec(ref_coords,std::array<REAL,3>{init + k * h,
                                                        init + j * h,
                                                        init + i * h});
            }

      for (int i = 0; i < r; ++i)
        {
          int incr_i = (r + 1) * (r + 1);
          for (int j = 0; j < r; ++j)
            {
              int incr_j = r + 1;
              pidx = i * incr_i + j * incr_j;
              for (int k = 0; k < r; ++k, pidx++)
                {
                  AppendToVec(ref_elems,
                              std::array<int,TPZVTK::MAX_PTS + 2>
                              {TPZVTK::CellType(ECube), 8, pidx, pidx + 1, pidx + incr_j + 1, pidx + incr_j,
                               pidx + incr_i, pidx + incr_i + 1, pidx + incr_i
                               + incr_j + 1, pidx + incr_j + incr_i});
                }
            }
        }
    }
}


void TPZVTKGenerator::FillReferencePrism(TPZVec<std::array<REAL,3>> &ref_coords, TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems)
{
  if (fSubdivision == 0){
    AppendToVec(ref_coords,std::array<REAL,3>{0.0, 0.0, -1.0});
    AppendToVec(ref_coords,std::array<REAL,3>{1.0, 0.0, -1.0});
    AppendToVec(ref_coords,std::array<REAL,3>{0.0, 1.0, -1.0});
    AppendToVec(ref_coords,std::array<REAL,3>{0.0, 0.0,  1.0});
    AppendToVec(ref_coords,std::array<REAL,3>{1.0, 0.0,  1.0});
    AppendToVec(ref_coords,std::array<REAL,3>{0.0, 1.0,  1.0});
    std::array<int,TPZVTK::MAX_PTS + 2> elem;
    elem[0] = TPZVTK::CellType(EPrisma);
    elem[1] = 6;
    for (int i = 0; i < MElementType_NNodes(EPrisma); i++){
      elem[i + 2] = i;
    }
    AppendToVec(ref_elems,elem);
  }
  else{
    const int r = 1 << fSubdivision;
    const int s = r + 1;

    const REAL h = 1.0 / r;
    const REAL hz = 2.0 /r;
    constexpr REAL initz = -1;
    int pidx = 0;
    for (int k = 0; k <= r; k++)
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i + j <= r; ++j)
          {
            AppendToVec(ref_coords,std::array<REAL,3>{j * h,
                                                      i * h,
                                                      initz + k * hz});
          }

    pidx = 0;
    for (int k = 0; k < r; k++)
      {
        int incr_k = (r + 2) * (r + 1) / 2;
        pidx = k * incr_k;
        for (int i = 0; i <= r; ++i)
          for (int j = 0; i + j <= r; ++j, pidx++)
            {
              // int pidx_curr = pidx;
              if (i + j == r)
                continue;
              int pidx_incr_i = pidx + 1;
              int pidx_incr_j = pidx + s - i;

              AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(EPrisma), 6,
                                                                        pidx, pidx_incr_i, pidx_incr_j, pidx + incr_k, pidx_incr_i + incr_k, pidx_incr_j + incr_k, 0, 0});

              int pidx_incr_ij = pidx_incr_j + 1;

              if (i + j + 1 < r)
                AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS + 2>{TPZVTK::CellType(EPrisma), 6,
                                                                          pidx_incr_i, pidx_incr_ij, pidx_incr_j, pidx_incr_i + incr_k, pidx_incr_ij + incr_k, pidx_incr_j + incr_k, 0, 0});
            }
      }
  }
}


void TPZVTKGenerator::PrintPointsLegacy()
{
  TPZSimpleTimer timer("PrintPts");
  
  (*fFileout) << "POINTS " << fPoints.size() << " float" << std::endl;
  for(const auto &p : fPoints){
    for(const auto &x : p){
      *fFileout << x<<'\t';
    }
    *fFileout << std::endl;
  }
}

/// output of cells in form vertices
void TPZVTKGenerator::PrintCellsLegacy()
{
  TPZSimpleTimer timer("PrintCells");
  
  // count number of data for cells, one + number of vertices
  int ndata = 0;
  for (auto &c : fCells){
    ndata++;
    ndata += c[1];
  }
  *fFileout << "CELLS " << fCells.size() << " " << ndata << std::endl;
  for (const auto &c : fCells){
    const int nv = c[1];
    *fFileout << nv << '\t';
    for (int i = 0; i < nv; i++)
      *fFileout << c[i + 2] << '\t';
    *fFileout << std::endl;
  }
}

/// output of cell types
void TPZVTKGenerator::PrintCellTypesLegacy(int meshdim)
{
  TPZSimpleTimer timer("PrintCellTypes");
  
  *fFileout << "CELL_TYPES " << fCells.size() << std::endl;

  for (const auto &c : fCells){
    *fFileout << c[0] << '\n';
  }
  *fFileout << "CELL_DATA " << fCells.size() << std::endl;
  *fFileout << "POINT_DATA " << fPoints.size() << std::endl;
}

/// output of field data (coefficient values)
void TPZVTKGenerator::PrintFieldDataLegacy()
{
  TPZSimpleTimer timer("PrintField");
  
  for (auto field : fFields){
    *fFileout << field->TypeName() <<' ' << field->Name()
              << " float " << std::endl;
    if(field->GetType() == TPZVTKField::Type::scal){
      *fFileout << "LOOKUP_TABLE default" << std::endl;
      for (const auto &v : *field)
      {*fFileout << v << ' ';}
      *fFileout << '\n';
    }else{
      //both tensor and vector should be written as
      //(1,2,3)
      //(1,2,3)
      //etc
      const auto nfs = field->size();
      for (int i = 0; i < nfs; i++){
        *fFileout << (*field)[i] << ' ';
        if((i+1)%3 == 0){
          *fFileout << '\n';
        }
      }
    }
    
    
  }
}

bool TPZVTKGenerator::IsValidEl(TPZCompEl *cel)
{
  if(!cel || ! cel->Reference()){return false;}
  
  const auto gel = cel->Reference();
  const auto geldim = gel->Dimension();
  const auto meshdim = gel->Mesh()->Dimension();
  if(geldim != meshdim){return false;}
  if(gel->HasSubElement()){return false;}
  if(gel->Type() == EPiramide){
    std::cout<<__PRETTY_FUNCTION__
             <<"\n pyramid element not supported yet! Aborting..."
             <<std::endl;
    DebugStop();
  }
  switch (gel->Type()){
    case EPoint:
    case EOned:
    case ETriangle:
    case EQuadrilateral:
    case ETetraedro:
    case EPiramide:
    case EPrisma:
    case ECube:
      return true;
    case EPolygonal: /*8*/	
    case EInterface: /*9*/	
    case EInterfacePoint: /*10*/	
    case EInterfaceLinear: /*11*/	
    case EInterfaceSurface: /*12*/	
    case ESubstructure: /*13*/	
    case EGlobLoc: /*14*/	
    case EDiscontinuous: /*15*/	
    case EAgglomerate: /*16*/	
    case ENoType: /*17*/	
      return false;
    }
}

void TPZVTKGenerator::Do(REAL time)
{
  TPZSimpleTimer timer("Do");

  const bool isCplxMesh = fCMesh->GetSolType() == ESolType::EComplex;

  std::ostringstream filenamefinal;
  std::stringstream appended;
  std::vector<int> datalength;
  int offs = 0;

  filenamefinal << fFilename;
  filenamefinal << "_step" << std::setw(5) << std::setfill('0')
                  << fOutputCount;

  fLastOutputName = filenamefinal.str();

  filenamefinal << ".vtk";
  
  if (fOutputCount > 0) {
    // cout << IM(4) << " ( " << fOutputCount << " )";
    const auto currt = fTimes.size();
    fTimes.resize(currt + 1);
    if (time == -1) {
      AppendToVec(fTimes, fOutputCount);
    } else {
      AppendToVec(fTimes, time);
    }
  } else {
    if (time != -1) {
      fTimes[0] = time;
    }
  }
  fOutputCount++;

  fFileout = new std::ofstream(filenamefinal.str());

  ResetArrays();

  TPZManVector<std::array<REAL, 3>,200> ref_vertices_line(0),
    ref_vertices_trig(0), ref_vertices_quad(0),
    ref_vertices_tet(0),  ref_vertices_hex(0),
    ref_vertices_prism(0);

  TPZVec<std::array<int, TPZVTK::MAX_PTS + 2>> ref_lines,
    ref_trigs, ref_quads, ref_tets, ref_hexes, ref_prisms;

  TPZManVector<std::array<REAL, 3>,200> ref_vertices;
  TPZVec<std::array<int, TPZVTK::MAX_PTS + 2>> ref_elems;


  const auto meshdim = fCMesh->Dimension();
  {
    TPZSimpleTimer timer("FillRefEls");
    if(meshdim == 3){
      FillReferenceTet(ref_vertices_tet, ref_tets);
      FillReferenceHex(ref_vertices_hex, ref_hexes);
      FillReferencePrism(ref_vertices_prism, ref_prisms);
    }else if(meshdim == 2){
      FillReferenceTrig(ref_vertices_trig, ref_trigs);
      FillReferenceQuad(ref_vertices_quad, ref_quads);
    }else{
      FillReferenceLine(ref_vertices_line, ref_lines);
    }
  }

  // header:
  *fFileout << "# vtk DataFile Version 3.0" << std::endl;
  *fFileout << "vtk output" << std::endl;
  *fFileout << "ASCII" << std::endl;
  *fFileout << "DATASET UNSTRUCTURED_GRID" << std::endl;

  
  
  for (auto cel : fCMesh->ElementVec()) {
    if (! IsValidEl(cel)){continue;}
    const auto eltype = cel->Reference()->Type();

    switch (eltype) {
    case EOned:
      ref_vertices = ref_vertices_line;
      ref_elems = ref_lines;
      break;
    case ETriangle:
      ref_vertices = ref_vertices_trig;
      ref_elems = ref_trigs;
      break;
    case EQuadrilateral:
      ref_vertices = ref_vertices_quad;
      ref_elems = ref_quads;
      break;
    case ETetraedro:
      ref_vertices = ref_vertices_tet;
      ref_elems = ref_tets;
      break;
    case ECube:
      ref_vertices = ref_vertices_hex;
      ref_elems = ref_hexes;
      break;
    case EPrisma:
      ref_vertices = ref_vertices_prism;
      ref_elems = ref_prisms;
      break;
    default:
      PZError << "VTK output for element type" << MElementType_Name(eltype)
              << "not supported";
      DebugStop();
    }

    const int offset = fPoints.size();
    const int eldim = cel->Dimension();
    TPZManVector<REAL, 3> qsi(eldim, 0), pt(3, 0.);
    std::array<REAL, 3> ptx;
    for (const auto &ip : ref_vertices) {
      for (auto x = 0; x < eldim; x++) {
        qsi[x] = ip[x];
      }
      cel->Reference()->X(qsi, pt);
      for (auto x = 0; x < 3; x++) {
        ptx[x] = pt[x];
      }
      AppendToVec(fPoints, ptx);
      // does it make a difference if printBound?
    }

    if (isCplxMesh) {
      ComputeFieldAtEl<CSTATE>(cel, ref_vertices, fFields);
    } else {
      ComputeFieldAtEl<STATE>(cel, ref_vertices, fFields);
    }

    for (const auto &elem : ref_elems) {
      std::array<int, TPZVTK::MAX_PTS + 2> new_elem = elem;
      for (int i = 2; i <= new_elem[0]; ++i)
        new_elem[i] += offset;
      AppendToVec(fCells, new_elem);
    }
  }

  PrintPointsLegacy();
  PrintCellsLegacy();
  PrintCellTypesLegacy(meshdim);
  PrintFieldDataLegacy();

  
  std::cout << " Done." << std::endl;
}


template<class TVar>
TPZPostProcEl<TVar>::TPZPostProcEl(TPZCompEl *cel) : fCel(cel){
  auto test = dynamic_cast<TPZMultiphysicsElement*>(cel);
  if(test){fIsMultiphysics = true;}
}

template<class TVar>
void TPZPostProcEl<TVar>::InitData(){
  if(fIsMultiphysics){
    auto mfcel = dynamic_cast<TPZMultiphysicsElement*>(fCel);
    const int64_t nref = mfcel->NMeshes();
    fDatavec.resize(nref);
    mfcel->InitMaterialData(fDatavec);
  }else{
    auto intel = dynamic_cast<TPZInterpolationSpace*>(fCel);
    intel->InitMaterialData(fMatdata);
  }
}

template<class TVar>
void TPZPostProcEl<TVar>::ComputeRequiredData(TPZVec<REAL> &qsi){
  if(fIsMultiphysics){
    auto mfcel = dynamic_cast<TPZMultiphysicsElement*>(fCel);
    const int64_t nref = mfcel->NMeshes();
    for(int ir = 0; ir < nref; ir++){
      fDatavec[ir].fNeedsSol= true;
    }
    TPZManVector<TPZTransform<> > trvec;
    mfcel->AffineTransform(trvec);
    mfcel->ComputeRequiredData(qsi, trvec, fDatavec);
  }else{
    auto intel = dynamic_cast<TPZInterpolationSpace*>(fCel);
    fMatdata.fNeedsSol = true;
		intel->ComputeRequiredData(fMatdata, qsi);
  }
}


template<class TVar>
void TPZPostProcEl<TVar>::Solution(const TPZVec<REAL> &qsi, const int id, TPZVec<TVar> &sol)
{
  TPZMaterial * material = fCel->Material();
  if(fIsMultiphysics){
    auto mfcel = dynamic_cast<TPZMultiphysicsElement*>(fCel);
    auto *matCombined =
       dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(material);
    matCombined->Solution(fDatavec,id,sol);
  }else{
    auto intel = dynamic_cast<TPZInterpolationSpace*>(fCel);
    auto *matSingle =
       dynamic_cast<TPZMatSingleSpaceT<TVar>*>(material);
    matSingle->Solution(fMatdata,id,sol);
  }
}
