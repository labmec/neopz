#include "tmbheidi.h"
#include "tmbxubinshih.h"
#include "tmbadapinterface.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPattern.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"

TMBAdaptInterface::TMBAdaptInterface(	TPZCompMesh *cmesh, int nstate, TPZVec<int> &dimstate,
  										TPZVec<int> &statetoanalyse,TPZVec<REAL> &sol)/* :
            					fMaxError(0),fMinError(0),fErInd(0),fErAnType(0) */{
  fMesh = cmesh;
  fSolution = &sol;
  if (dimstate.NElements() != nstate) {
  	cout << "TMBAdaptInterface::TMBAdaptInterface ERROR - incorrect number of state variables\n";
   	exit (-1);
  }
  fDimState = dimstate;
  fStateAnalyse = statetoanalyse;
  fMaxLevel = 100;
  fMark.Resize(0);
}

TMBAdaptInterface:: ~TMBAdaptInterface (){
}

//void TMBAdaptInterface::SetStateVariables(TPZVec<int> &state){
//  fStateAnalyse = state;
//}

void TMBAdaptInterface::SetMaxMinError (TPZVec<REAL> &maxervec, TPZVec<REAL> &minervec){
  fMaxError = maxervec;
  fMinError = minervec;
}

TPZCompMesh *TMBAdaptInterface::GetAdaptedMesh(TPZVec<int> &erind, TPZVec<int> &erantype,
                                               bool level_check, int reftype,
                                               int side_ref_dim, int side_ref_state){
  if (erind.NElements() != fStateAnalyse.NElements() || erantype.NElements()!= fStateAnalyse.NElements()) {
    cout << "TMBAdaptInterface::GetAdaptedMesh ERROR - incorrect number of state variables\n";
    exit (-1);
  }
  fErInd = erind;
  fErAnType = erantype;
  fMark.Resize(fMesh->ElementVec().NElements());
  fMark.Fill(0);
  
  TPZVec<int> side_ref;
  MarkedElListH (fMark,side_ref,side_ref_dim,side_ref_state);
//  cout << fMark << endl;
  int i;
  int nel = fMesh->ElementVec().NElements();
  switch (reftype){
    case 1:{
      //The side based has not aglomerate function (how refinement pattern was used?)
      for (i=0;i<nel;i++){
        TPZCompEl * el = fMesh->ElementVec()[i];
        if (!el || !fMark[i]) continue;
        if (fMark[i] > 0) SideRefine(el,side_ref[i]);
      }
      break;
    }
    default:{
      //Uniform refinement pattern
      for (i=0;i<nel;i++){
        TPZCompEl * el = fMesh->ElementVec()[i];
        if (!el || !fMark[i]) continue;
        if (fMark[i] > 0) Refine(el,level_check);
//        if (fMark[i] < 0) Aglomerate(el);
      }
      break;
    }
  }
  fMesh->AdjustBoundaryElements();
  return fMesh;
}

void TMBAdaptInterface::MarkedElListH(TPZVec<int> &elindex, TPZVec<int> &side, int sidedim, int sidestate){
  int i;
  TPZStack <int> state_heidi;
  TPZStack <int> state_XuShih;
  TPZStack <int> ertype_heidi;
  TPZStack <int> ertype_XuShih;
  TPZStack <REAL> max_heidi;
  TPZStack <REAL> min_heidi;
  TPZStack <REAL> max_XuShih;
  TPZStack <REAL> min_XuShih;

  elindex.Fill(0);
  
  for (i=0;i<fStateAnalyse.NElements();i++){
    switch (fErInd[i]){
      case 0:{
        state_heidi.Push(fStateAnalyse[i]);
        ertype_heidi.Push(fErAnType[i]);
        max_heidi.Push(fMaxError[i]);
        min_heidi.Push(fMinError[i]);
        break;
      }
      case 1:{
        state_XuShih.Push(fStateAnalyse[i]);
        ertype_XuShih.Push(fErAnType[i]);
        max_XuShih.Push(fMaxError[i]);
        min_XuShih.Push(fMinError[i]);
        break;
      }
      default:{
        state_heidi.Push(fStateAnalyse[i]);
        ertype_heidi.Push(fErAnType[i]);
        max_heidi.Push(fMaxError[i]);
        min_heidi.Push(fMinError[i]);
        break;
      }
    }
  }
  if (state_heidi.NElements() > 0){
    TMBHeidi heidi(fDimState.NElements(), fDimState,state_heidi,*fSolution);
    heidi.SetError(max_heidi,min_heidi,ertype_heidi);
    heidi.MarkedElListH(elindex,side,sidedim,sidestate);
  }
  if (state_XuShih.NElements() > 0){
    int level = MaxLevel(fMesh);
    if (level < 2 || !(level%2)){
      TMBXubinShih xushi(fMesh,fDimState.NElements(), fDimState,state_XuShih,*fSolution);
      xushi.SetError(max_XuShih,min_XuShih,ertype_XuShih);

      xushi.MarkedElListH(elindex,side,sidedim,sidestate);
    }
    else{
      TMBXubinShih xushi(fMesh,fDimState.NElements(), fDimState,state_XuShih,*fSolution,1);
      xushi.SetError(max_XuShih,min_XuShih,ertype_XuShih);
      xushi.MarkedElListH(elindex,side,sidedim,sidestate);
    }
  }
}

void TMBAdaptInterface::Refine(TPZCompEl *el, bool level_check){
  if (!el) return;
  TPZVec<int> subelindex;
  TPZGeoEl *gel = el->Reference();
  int level = gel->Level();
  //If the element is already in the maximum level...
  if (level >= fMaxLevel) return;
  //The element will be refined...
  if (level_check){
    level ++;
    int i;
    int nsides = gel->NSides();
    for (i=0;i<nsides;i++){
      TPZGeoElSide neighside = gel->Neighbour(i);
      TPZGeoEl *neigh = neighside.Element();
      while (neigh && neigh != gel){
        int neighlevel = neigh->Level();
        if (!neigh->Reference()) {
          neighside = neigh->Neighbour(neighside.Side());
          neigh = neighside.Element();
          continue;
        }
        int neighindex = neigh->Reference()->Index();
        if ((level-neighlevel) >= 2){
          Refine(neigh->Reference(),level_check);
          fMark[neighindex] = 0;
        }
        neighside = neigh->Neighbour(neighside.Side());
        neigh = neighside.Element();
      }
    }
  }
  el->Divide(el->Index(),subelindex,0);
  return;
}

void TMBAdaptInterface::Aglomerate(TPZCompEl *el){

	//verificar possibilidade de aglomerar
	TPZGeoEl *father = el->Reference()->Father();
	if (!father) return;
	int level = father->Level();
	int i;
	TPZStack<TPZGeoElSide> sidesubelem;
	father->GetSubElements2(father->NSides()-1, sidesubelem, father->Dimension());
	int nsubel = sidesubelem.NElements();
	int sum = 0;
	int index;
	for (i=0;i<nsubel;i++){
		TPZGeoEl *subel = sidesubelem[i].Element();
		index = subel->Reference()->Index();
		if (fMark[index] == -1) sum++;	
	}

	if ( ((REAL)sum / (REAL)nsubel ) < 0.5) {
		for (i=0;i<nsubel;i++){
			TPZGeoEl *subel = sidesubelem[i].Element();
			index = subel->Reference()->Index();
			if (fMark[index] == -1 ) fMark[index] = 0;
		}
		return;
	}
	
	for (i=0;i<father->NSides();i++){
		TPZGeoElSide neighside = father->Neighbour(i);
		TPZGeoEl *neigh = neighside.Element();
		if (!neigh) continue;
		int neighlevel = neigh->Level();
		int neighindex = neigh->Reference()->Index();
		if ((level-neighlevel) > 2 && fMark[neighindex] != -1) return;
	}
	
	TPZStack< TPZGeoElSide > subel;
	int side = father->NSides()-1;
	int dimension = father->Dimension();
	
	father->GetSubElements2 (side,subel,dimension);
	TPZVec<int> subelindex (subel.NElements(),-1);
	for (i=0;i<subel.NElements();i++){
		subelindex[i] = subel[i].Element()->Reference()->Index();
	}
	fMesh->Coarsen(subelindex,index);
}

int TMBAdaptInterface::MaxLevel(TPZCompMesh *mesh) {
  int nel = mesh->NElements();
  int el;
  int level = 0;
  for(el=0; el<nel; el++) {
    TPZCompEl *cel = mesh->ElementVec()[el];
    if(!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if(!gel) continue;
    int gellev = gel->Level();
    level = (level <gellev) ? gellev : level;
  }
  return level;
}


/*void TMBAdaptInterface::SetPZMeshSolution(){
  int nvolcont = fSolution.GetNumberVolumeContainers();
  int nvols, i, j, globalindex;
  const int nsize = 10;
  const int vsize = 20;
  int nstate = fStateVec.NElements();
  TPZManVector<int, nstate> statepos(nstate);
  TPZManVector<REAL,nstate> values(nstate);
  TPZStack<int> volindex;
  FillVolIndex(volindex);
  //  OutVars[0] = RHO; OutVars[1] = RHOU; OutVars[2] = RHOV; OutVars[3] = RHOW; OutVars[4] = E;

  TMBLocalData *localdata;

  for(i=0;i<NumVolConts;i++){
    nvols = fSolution.GetVolumeContainer(i).GetNumberVolumes();
    localdata = fState.GetLocalData(i);
    localdata->SetExternalIndices(fState, statepos);

    for(j=0;j<nvols;j++){
      LocalState->GetDataValue(values,j,statepos);
      globaindex =  fSolution.GetVolumeContainer(i).LocalToGlobal(j);
      TPZCompEl * el = fMesh->ElementVec()[volindex[globalindex]];

      TPZConnect connect = el->Connect(0);
      int sn = connect.SequenceNumber();
      int bl_pos = fMesh->Block().Position(sn);
      int bl_siz = fMesh->Block().Size(sn);

      if (bl_siz != nstate){
	cout << "TMBAdaptInterface::Error - solution block size inconsistancy\n";
	exit(-1);
      }

      for (k=0;k<bl_siz;k++){
	fMesh->Solution()(bl_pos+k,0) = values[k];
      }
    }
  }
  
}*/

/*void TMBAdaptInterface::FillVolIndex(TPZStack<int> &volindex){
  int i;
  int nel = fMesh->ElementVec().NElements();
  for (i=0;i<nel;i++){
    TPZCompEl *cel = fMesh->ElementVec()[i];
    if (!cel) continue;
    int eltype = cel->Type();
    if (eltype == EDiscontinuous){
      if (cel->Dimension() == 3) volindex.Push(i);
    }
  }
  
}*/
/** Implements side based refinement for the specified element. */
void TMBAdaptInterface::SideRefine(TPZCompEl *cel, int side){
  int i,j;
  TPZCompElSide celside(cel,side);
  TPZStack<TPZCompElSide> elsidevec;
  elsidevec.Push(celside);
  TPZVec<int> subelindex;
  celside.EqualLevelElementList(elsidevec,0,0);
  TPZStack<TPZGeoElSide> geosidevec;
  for (j=0;j<elsidevec.NElements();j++){
    geosidevec.Push(elsidevec[j].Reference());
  }
  int nel = elsidevec.NElements();
  TPZManVector<TPZConnect *,3> con(3,0),conc(3,0);
  for(i=0;i<nel;i++) {
    TPZCompEl *celneig = elsidevec[i].Element();
    TPZGeoEl *gel = celneig->Reference();
    int refside = elsidevec[i].Side();
    bool chg = ChangeRefPattern(gel,refside);
    if (!chg) return;
    celneig->Divide(celneig->Index(),subelindex,0);
    TPZStack<TPZGeoElSide> subels;
    gel->GetSubElements2(refside,subels);
    if(con[0] == 0) {
      con[0] = &subels[0].Element()->Reference()->Connect(subels[0].Side());
      con[1] = &subels[1].Element()->Reference()->Connect(subels[1].Side());
      con[2] = &subels[2].Element()->Reference()->Connect(subels[2].Side());
    } else {
      conc[0] = &subels[0].Element()->Reference()->Connect(subels[0].Side());
      conc[1] = &subels[1].Element()->Reference()->Connect(subels[1].Side());
      conc[2] = &subels[2].Element()->Reference()->Connect(subels[2].Side());
      if(con[0] != conc[0]) {
        cout << "fodeu1\n";
        cout << geosidevec << endl;
        for (j=0;j<geosidevec.NElements();j++){
          int nsub = geosidevec[j].Element()->NSubElements();
          for (int isub=0;isub<nsub;isub++){
            geosidevec[j].Element()->SubElement(isub)->Print(cout);
          }
        }
      }
      if(con[1] != conc[1] && con[1] != conc[2]) {
        cout << "fodeu2\n";
        cout << geosidevec << endl;
      }
      if(con[2] != conc[1] && con[2] != conc[2]) {
        cout << "fodeu3\n";
        cout << geosidevec << endl;
      }
    }
    int index = celneig->Index();
    if (index < fMark.NElements()) fMark[index] = 0;
  }
}

bool TMBAdaptInterface::ChangeRefPattern(TPZGeoEl * gel, int side){
 
  TPZRefPattern *refpat = fMesh->Reference()->GetRefPattern(gel,side);
  int eltype = gel->Type();
  switch (eltype){
    case (EPoint) :{
      TPZGeoElRefPattern<TPZShapePoint,TPZGeoPoint> *gelref =
                               dynamic_cast <TPZGeoElRefPattern<TPZShapePoint,TPZGeoPoint>* >(gel);
      if (!gelref) return false;        
      gelref->SetRefPattern(refpat);
      break;
    }
    case (EOned) :{
      TPZGeoElRefPattern<TPZShapeLinear,TPZGeoLinear> *gelref =
                               dynamic_cast <TPZGeoElRefPattern<TPZShapeLinear,TPZGeoLinear>* >(gel);
      if (!gelref) return false;
      gelref->SetRefPattern(refpat);
      break;
    }
    case (ETriangle) :{
      TPZGeoElRefPattern<TPZShapeTriang,TPZGeoTriangle> *gelref =
                               dynamic_cast <TPZGeoElRefPattern<TPZShapeTriang,TPZGeoTriangle>* >(gel);
      if (!gelref) return false;
      gelref->SetRefPattern(refpat);
      break;
    }
    case (EQuadrilateral) :{
      TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad> *gelref =
                               dynamic_cast <TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad>* >(gel);
      if (!gelref) return false;
      gelref->SetRefPattern(refpat);
      break;
    }
    case (ETetraedro) :{
      TPZGeoElRefPattern<TPZShapeTetra,TPZGeoTetrahedra> *gelref =
                               dynamic_cast <TPZGeoElRefPattern<TPZShapeTetra,TPZGeoTetrahedra>* >(gel);
      if (!gelref) return false;
      gelref->SetRefPattern(refpat);
      break;
    }
    case (EPiramide) :{
      TPZGeoElRefPattern<TPZShapePiram,TPZGeoPyramid> *gelref =
                               dynamic_cast <TPZGeoElRefPattern<TPZShapePiram,TPZGeoPyramid>* >(gel);
      if (!gelref) return false;
      gelref->SetRefPattern(refpat);
      break;
    }
    case (EPrisma) :{
      TPZGeoElRefPattern<TPZShapePrism,TPZGeoPrism> *gelref =
                               dynamic_cast <TPZGeoElRefPattern<TPZShapePrism,TPZGeoPrism>* >(gel);
      if (!gelref) return false;
      gelref->SetRefPattern(refpat);
      break;
    }
    case (ECube) :{
      TPZGeoElRefPattern<TPZShapeCube,TPZGeoCube> *gelref =
                               dynamic_cast <TPZGeoElRefPattern<TPZShapeCube,TPZGeoCube>* >(gel);
      if (!gelref) return false;
      gelref->SetRefPattern(refpat);
      break;
    }
    default: {
      cout << "TMBAdaptInterface::ChangeRefPattern ERROR : Unidentied element! "  << endl;
      return false;
    }
  }
  return true;
}
