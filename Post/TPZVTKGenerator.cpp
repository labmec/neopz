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

#include "TPZRefPatternDataBase.h"


#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"
#include "tpzpyramid.h"

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
                      const TPZVec<TPZManVector<REAL,3>> &ref_vertices,
                      TPZVec<TPZAutoPointer<TPZVTKField>>& fields){
  TPZManVector<TVar,9> sol;

  const auto celdim = cel->Dimension();

  TPZPostProcEl<TVar> graphel(cel);
  graphel.InitData();
  for (auto &ip : ref_vertices){
    ip.Resize(celdim);
    graphel.ComputeRequiredData(ip);
    for (int i = 0; i < fields.size(); i++){
      auto &field = *(fields[i]);
      auto fdim = field.Dimension();
      sol.Resize(fdim);
      graphel.Solution(ip, field.Id(), sol);
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

template <class TOPOL>
void TPZVTKGenerator::FillReferenceEl(TPZVec<TPZManVector<REAL,3>> &ref_coords,
                                      TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems){
  if(fSubdivision == 0)
  {
    static constexpr auto nnodes = TOPOL::NCornerNodes;
    ref_coords.Resize(nnodes,{0,0,0});
    //create array with all entries = nnodes
    AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS+2>{});
    ref_elems[0][0] = TPZVTK::CellType(TOPOL::Type());
    ref_elems[0][1] = nnodes;
    for(int i = 0; i < nnodes; i++){
      TOPOL::ParametricDomainNodeCoord(i,ref_coords[i]);
      ref_elems[0][i+2] = i;
    }
  }else
  {
    auto refp = gRefDBase.GetUniformRefPattern(TOPOL::Type());
    auto refpmesh = refp->RefPatternMesh();


    const int nrefs = fSubdivision - 1;

    TPZManVector<TPZGeoEl*,TPZVTK::MAX_SUBEL> sons(TPZVTK::MAX_SUBEL);
    for(int i = 0; i < nrefs; i++){
      for(auto gel : refpmesh.ElementVec()){
        if(gel->HasSubElement()==false){
          gel->Divide(sons);
        }
      }
    }
    
    
    //fill nodesg
    const auto nnodes = refpmesh.NNodes();
    ref_coords.Resize(nnodes,{0,0,0});
    for(auto in = 0; in < nnodes; in++){
      auto &node = refpmesh.NodeVec()[in];
      node.GetCoordinates(ref_coords[in]);
    }


    //fill cells
    int ic = 0;//i-th cell
    for(TPZGeoEl* gel : refpmesh.ElementVec()){
      if(gel->HasSubElement()==false){
        const int gnnodes = gel->NNodes();
        AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS+2>{});
        ref_elems[ic][0] = TPZVTK::CellType(gel->Type());
        ref_elems[ic][1] = gnnodes;
        for(int in = 0; in < gnnodes; in++){
          ref_elems[ic][in+2] = gel->NodeIndex(in);
        }
        ic++;
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

  fFileout = new std::ofstream(filenamefinal.str());

  ResetArrays();

  TPZManVector<TPZManVector<REAL, 3>,200> ref_vertices_line(0),
    ref_vertices_trig(0), ref_vertices_quad(0),ref_vertices_tet(0),
    ref_vertices_prism(0), ref_vertices_hex(0), ref_vertices_pyr(0) ;

  TPZVec<std::array<int, TPZVTK::MAX_PTS + 2>> ref_lines, ref_trigs,
    ref_quads, ref_tets, ref_prisms, ref_hexes, ref_pyrs;

  TPZManVector<TPZManVector<REAL, 3>,200> ref_vertices;
  TPZVec<std::array<int, TPZVTK::MAX_PTS + 2>> ref_elems;


  const auto meshdim = fCMesh->Dimension();
  {
    TPZSimpleTimer timer("FillRefEls");
    if(fSubdivision && !fOutputCount){//just need to do it once
      TPZManVector<TPZManVector<MElementType,3>,4> eltypes{
        {},
        {EOned},
        {ETriangle, EQuadrilateral},
        {ETetraedro, EPrisma, ECube, EPiramide}
      };
      for( auto type : eltypes[meshdim]){
        auto refp = gRefDBase.GetUniformRefPattern(type);
        if(!refp){gRefDBase.InitializeUniformRefPattern(type);}
      }
    }
    if (meshdim == 3){
      FillReferenceEl<pztopology::TPZPyramid>(ref_vertices_pyr, ref_pyrs);
      FillReferenceEl<pztopology::TPZTetrahedron>(ref_vertices_tet, ref_tets);
      FillReferenceEl<pztopology::TPZCube>(ref_vertices_hex, ref_hexes);
      FillReferenceEl<pztopology::TPZPrism>(ref_vertices_prism, ref_prisms);
    }
    else if (meshdim == 2){
      FillReferenceEl<pztopology::TPZTriangle>(ref_vertices_trig, ref_trigs);
      FillReferenceEl<pztopology::TPZQuadrilateral>(ref_vertices_quad, ref_quads);
    }else{
      FillReferenceEl<pztopology::TPZLine>(ref_vertices_line, ref_lines);
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
    TPZManVector<REAL, 3> pt(3, 0.);
    for (auto &ip : ref_vertices) {
      cel->Reference()->X(ip, pt);
      AppendToVec(fPoints, pt);
    }

    if (isCplxMesh) {
      ComputeFieldAtEl<CSTATE>(cel, ref_vertices, fFields);
    } else {
      ComputeFieldAtEl<STATE>(cel, ref_vertices, fFields);
    }

    for (const auto &elem : ref_elems) {
      std::array<int, TPZVTK::MAX_PTS + 2> new_elem = elem;
      const int npts = new_elem[1];
      for (int i = 0; i <= npts; ++i)
        new_elem[i+2] += offset;
      AppendToVec(fCells, new_elem);
    }
  }

  PrintPointsLegacy();
  PrintCellsLegacy();
  PrintCellTypesLegacy(meshdim);
  PrintFieldDataLegacy();

  fOutputCount++;
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


#define FILLREF(T)\
  template \
  void TPZVTKGenerator::FillReferenceEl<T>(TPZVec<TPZManVector<REAL,3>> &, \
                                           TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &);

FILLREF(pztopology::TPZTriangle)
FILLREF(pztopology::TPZQuadrilateral)
FILLREF(pztopology::TPZTetrahedron)
FILLREF(pztopology::TPZCube)
FILLREF(pztopology::TPZPrism)
#undef FILLREF
