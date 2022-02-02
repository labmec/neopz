#include "TPZH1ApproxSpace.h"
#include "TPZCompElH1.h"

#include "pzgeoel.h"

/*
  The following functions are defined at
  pzcreateapproxspace.cpp (without the H1 in the name)
  and should be moved to TPZCompElH1
*/
TPZCompEl* CreateH1PointEl(TPZGeoEl*,TPZCompMesh&,const H1Family);
TPZCompEl* CreateH1LineEl(TPZGeoEl*,TPZCompMesh&,const H1Family);
TPZCompEl* CreateH1TriangleEl(TPZGeoEl*,TPZCompMesh&,const H1Family);
TPZCompEl* CreateH1QuadEl(TPZGeoEl*,TPZCompMesh&,const H1Family);
TPZCompEl* CreateH1TetraEl(TPZGeoEl*,TPZCompMesh&,const H1Family);
TPZCompEl* CreateH1PyramEl(TPZGeoEl*,TPZCompMesh&,const H1Family);
TPZCompEl* CreateH1PrismEl(TPZGeoEl*,TPZCompMesh&,const H1Family);
TPZCompEl* CreateH1HexaEl(TPZGeoEl*,TPZCompMesh&,const H1Family);



TPZCompEl *TPZH1ApproxSpace::CreateCompEl(TPZGeoEl *gel, TPZCompMesh &mesh) const{

#ifdef PZDEBUG
  if(!gel){
    PZError<<__PRETTY_FUNCTION__
           <<"\nERROR: geometric element does not exist.\nAborting..."
           <<std::endl;
    DebugStop();
  }
#endif // PZDEBUG
  
  const H1Family h1fam = this->fFamily;

  switch(gel->Type()){
  case EPoint:
    return CreateH1PointEl(gel,mesh,h1fam);
    break;
	case EOned:
    return CreateH1LineEl(gel,mesh,h1fam);
    break;
	case ETriangle:
    return CreateH1TriangleEl(gel,mesh,h1fam);
    break;
	case EQuadrilateral:
    return CreateH1QuadEl(gel,mesh,h1fam);
    break;
	case ETetraedro:
    return CreateH1TetraEl(gel,mesh,h1fam);
    break;
	case EPiramide:
    return CreateH1PyramEl(gel,mesh,h1fam);
    break;
	case EPrisma:
    return CreateH1PrismEl(gel,mesh,h1fam);
    break;
	case ECube:
    return CreateH1HexaEl(gel,mesh,h1fam);
    break;
	case EPolygonal:
	case EInterface:
	case EInterfacePoint:
	case EInterfaceLinear:
	case EInterfaceSurface:
	case ESubstructure:
	case EGlobLoc:
	case EDiscontinuous:
	case EAgglomerate:
	case ENoType:
    PZError<<__PRETTY_FUNCTION__
           <<"\nInvalid element type.Aborting.."
           <<std::endl;
    DebugStop();
  }
}


int TPZH1ApproxSpace::ClassId() const{
  return Hash("TPZH1ApproxSpace") ^ TPZApproxSpace::ClassId() << 1;
}

//for read and write method
constexpr static int64_t tag_standard{1};
constexpr static int64_t tag_enrichedprism{2};

void TPZH1ApproxSpace::Read(TPZStream &buf, void *context){
  TPZApproxSpace::Read(buf,context);
  int64_t val{0};
  buf.Read(&val);
  if(val == tag_standard){
    this->fFamily = H1Family::Standard;
  }else if(val == tag_enrichedprism){
    this->fFamily = H1Family::EnrichedPrism;
  }
}
void TPZH1ApproxSpace::Write(TPZStream &buf, int withclassid) const{
  TPZApproxSpace::Write(buf,withclassid);
  
  switch(this->fFamily){
  case H1Family::Standard:
    buf.Write(tag_standard);
    break;
  case H1Family::EnrichedPrism:
    buf.Write(tag_enrichedprism);
    break;
  }
}