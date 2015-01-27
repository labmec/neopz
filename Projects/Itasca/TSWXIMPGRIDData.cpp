//---------------------------------------------------------------------------

#include "TSWXIMPGRIDData.h"

void TSWXGridData::Print(std::ostream &myfile) {
  for(unsigned int i = 0; i < fGridPoints.size(); i++){
    myfile << "G " << fGridPoints[i].fId << " " << fGridPoints[i].fX[0]
                                         << " " << fGridPoints[i].fX[1]
                                         << " " << fGridPoints[i].fX[2] << "\n";
  }

  for(unsigned int i = 0; i < fB8Elements.size(); i++){
    myfile << "Z B8 " << fB8Elements[i].fId;
    for(int k = 0; k < 8; k++) myfile << " " << fB8Elements[i].fIndices[k];
    myfile << "\n";
  }

  for(unsigned int i = 0; i < fZoneGroups.size(); i++){
    myfile << "ZGROUP " << fZoneGroups[i].fGroupName;
    for(unsigned int k = 0; k < fZoneGroups[i].fEls.size(); k++){
      if(k%15 == 0) myfile << "\n";
      myfile << " " << fZoneGroups[i].fEls[k];
    }
    myfile << "\n";
  }
}///void



