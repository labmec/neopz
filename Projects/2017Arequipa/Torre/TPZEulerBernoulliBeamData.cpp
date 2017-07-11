#include "TPZEulerBernoulliBeamData.h"

TPZEulerBernoulliBeamData::TPZEulerBernoulliBeamData(){
  fGravity.Redim(3,1);
  fGravity(2,0) = -9.81;
  this->fCurrenttlyComputingStiffness = 1;
}

TPZEulerBernoulliBeamData::~TPZEulerBernoulliBeamData(){
  //nothing here
}

void TPZEulerBernoulliBeamData::Print(std::ostream &out) const{
  out << "fMaterialProp:\n";
  for(std::map<int, MaterialProperties>::const_iterator w = fMaterialProp.begin(); w != fMaterialProp.end(); w++){
    out << "MaterialId: " << w->first << "\n";
    w->second.Print(out);
  }
  out << "fSectionProp:\n";
  for(std::map<int, SectionProperties>::const_iterator w = fSectionProp.begin(); w != fSectionProp.end(); w++){
    out << "SectionId: " << w->first << "\n";
    w->second.Print(out);
  }
  out << "fCurrenttlyComputingStiffness: " << fCurrenttlyComputingStiffness << "\n";
}//Print

void TPZEulerBernoulliBeamData::GetData(int matId, int sectionId, MaterialProperties& mat, SectionProperties& sec) const{

  {//escope only
    std::map<int, TPZEulerBernoulliBeamData::MaterialProperties>::const_iterator wm = this->fMaterialProp.find( matId );
    if(wm == this->fMaterialProp.end()){
      DebugStop();
    }
    mat = wm->second;
  }

  {//escope only
    std::map<int, TPZEulerBernoulliBeamData::SectionProperties>::const_iterator wm = this->fSectionProp.find( sectionId );
    if(wm == this->fSectionProp.end()) DebugStop();
    sec = wm->second;
  }

}//GetData

void TPZEulerBernoulliBeamData::g(TPZFMatrix<STATE> & g) const{
  g = this->fGravity;
}

void TPZEulerBernoulliBeamData::SetGravity(TPZFMatrix<STATE> &g){
  if(g.Rows() != 3 || g.Cols() != 1) DebugStop();
  this->fGravity = g;
}


