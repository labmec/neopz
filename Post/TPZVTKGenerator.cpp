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
                      TPZVec<TPZAutoPointer<TPZVTKField>>& fields,
                      const TPZVec<int> &init_pos){
  TPZManVector<TVar,9> sol;

  const auto celdim = cel->Dimension();

  TPZPostProcEl<TVar> graphel(cel);
  graphel.InitData();

  int iv = 0;
  for (auto &ip : ref_vertices){
    ip.Resize(celdim);
    graphel.ComputeRequiredData(ip);
    const int nfields = fields.size();
    for (int i = 0; i < nfields; i++){
      auto &field = *(fields[i]);
      auto fdim = field.Dimension();
      auto pos = init_pos[i] + fdim*iv;
      sol.Resize(fdim);
      graphel.Solution(ip, field.Id(), sol);
      const auto sz = sol.size();
      if constexpr (std::is_same_v<TVar,CSTATE>){
        for (int d = 0; d < sz; ++d){
          field[pos++] =  std::real(sol[d]);
        }
        for (int d = sz; d < fdim; ++d){
          field[pos++] =  0.0;
        }
      }else{
        for (int d = 0; d < fdim; ++d){
          field[pos++] = sol[d];
        }
        for (int d = sz; d < fdim; ++d){
          field[pos++] = sol[d];
        }
      }
    }
    iv++;
  }
}


TPZVTKField::TPZVTKField(TPZVTKField::Type type, int id, std::string aname) :
  fType( type), fId(id), fName(aname) { ; }

TPZVTKGenerator::TPZVTKGenerator(TPZAutoPointer<TPZCompMesh> cmesh,
                                 const TPZVec<std::string> fields,
                                 std::string filename,
                                 int vtkres)
  : TPZVTKGenerator(cmesh.operator->(),fields,filename,vtkres)//delegates to other ctor
{
}



TPZVTKGenerator::TPZVTKGenerator(TPZCompMesh* cmesh,
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
  const int meshdim = fCMesh->Dimension();
  FillRefEls(meshdim);
  ComputePoints();
}

void TPZVTKGenerator::FillRefEls(const int meshdim)
{
  TPZSimpleTimer timer("FillRefEls");
  
  TPZManVector<TPZManVector<MElementType,3>,4> eltypes{
    {},
    {EOned},
    {ETriangle, EQuadrilateral},
    {ETetraedro, EPrisma, ECube, EPiramide}
  };

  for( auto type : eltypes[meshdim]){
    auto refp = gRefDBase.GetUniformRefPattern(type);
    if(!refp){gRefDBase.InitializeUniformRefPattern(type);}
    TPZVec<TPZManVector<REAL,3>> ref_coords;
    TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> ref_elems;
    switch(type){
    case EOned:
      FillReferenceEl<pztopology::TPZLine>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case ETriangle:
      FillReferenceEl<pztopology::TPZTriangle>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case EQuadrilateral:
      FillReferenceEl<pztopology::TPZQuadrilateral>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case ETetraedro:
      FillReferenceEl<pztopology::TPZTetrahedron>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case EPrisma:
      FillReferenceEl<pztopology::TPZPrism>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case ECube:
      FillReferenceEl<pztopology::TPZCube>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case EPiramide:
      FillReferenceEl<pztopology::TPZPyramid>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    default:
      unreachable();
    }
  }
}

void TPZVTKGenerator::ComputePoints()
{
  fPoints.resize(0);
  fElementVec.resize(0);
  TPZManVector<TPZCompEl*,10> compelvec;
  for (auto orig_cel : fCMesh->ElementVec()) {
    if(!orig_cel){continue;}
    /**
       orig_cel might be a TPZSubCompMesh, TPZElementGroup, etc.
       So we need to check exactly how many computational elements are "inside" it.
       note: the same logic does not apply for multiphysics elements, 
       where post-processing each element individually wouldn't make sense.
    */
    TPZStack<TPZCompEl*> cellist;
    orig_cel->GetCompElList(cellist);
    for(auto cel : cellist){
      if (! IsValidEl(cel)){continue;}
      const auto type = cel->Reference()->Type();
      //add to valid elements for post-processing

      const int offset = fPoints.size();
      AppendToVec(fElementVec,std::make_pair(cel,offset));
      const int eldim = cel->Dimension();
      TPZManVector<REAL, 3> pt(3, 0.);
      for (auto &ip : fRefVertices[type]) {
        cel->Reference()->X(ip, pt);
        AppendToVec(fPoints, pt);
      }
      for (const auto &elem : fRefEls[type]) {
        std::array<int, TPZVTK::MAX_PTS + 2> new_elem = elem;
        const int npts = new_elem[1];
        for (int i = 0; i <= npts; ++i)
          new_elem[i+2] += offset;
        AppendToVec(fCells, new_elem);
      }
    }
    //at this point we know how many points there are
    const int npts = fPoints.size();
    for(auto &f : fFields){
      f->Resize(0);
      const int fdim = f->Dimension();
      f->Resize(npts*fdim);
      //we first resize to zero to ensure that we allocate the exact size
    }
  }
}

/// Empty all fields, points and cells
void TPZVTKGenerator::ResetArrays()
{
  fElementVec.Resize(0);
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

  if(fPoints.size() == 0){
    //perhaps the mesh has changed
    ComputePoints();
  }

  const auto meshdim = fCMesh->Dimension();

  // header:
  *fFileout << "# vtk DataFile Version 3.0" << std::endl;
  *fFileout << "vtk output" << std::endl;
  *fFileout << "ASCII" << std::endl;
  *fFileout << "DATASET UNSTRUCTURED_GRID" << std::endl;

  const int nfields = fFields.size();
  TPZVec<int> posvec(nfields);
  
  for (auto [cel,pos] : fElementVec) {
    const auto eltype = cel->Reference()->Type();

    for(auto f = 0; f < nfields; f++){
      posvec[f] = pos * fFields[f]->Dimension();
    }
    if (isCplxMesh) {
      ComputeFieldAtEl<CSTATE>(cel, fRefVertices[eltype], fFields,posvec);
    } else {
      ComputeFieldAtEl<STATE>(cel, fRefVertices[eltype], fFields,posvec);
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
