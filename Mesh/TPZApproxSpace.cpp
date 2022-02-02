#include "TPZApproxSpace.h"
#include "pzcmesh.h"
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzcreateapproximationspace");
#endif

/**
   WE STILL NEED TO IMPLEMENT ALL THE REMAINING METHODS
**/

/** @brief Creates the computational elements, and the degree of freedom nodes */
void TPZApproxSpace::BuildMesh(TPZCompMesh &cmesh) const {
  std::set<int> materialids;
  std::map<int,TPZMaterial *>::iterator it;
  for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
    materialids.insert(it->first);
  }
  BuildMesh(cmesh,materialids);
}

void TPZApproxSpace::BuildMesh(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs) const {
  TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh.Reference()->ElementVec();
  int64_t i, nelem = elvec.NElements();
  int64_t neltocreate = 0;
  int64_t index;
  for(i=0; i<nelem; i++) {
    TPZGeoEl *gel = elvec[i];
    if(!gel) continue;
    if (gel->Reference()) {
      continue;
    }
    if(!gel->HasSubElement()) {
      neltocreate++;
    }
  }
  std::set<int> matnotfound;
  int64_t nbl = cmesh.Block().NBlocks();
  if(neltocreate > nbl)
    {
      cmesh.Block().SetNBlocks(neltocreate);
    }
  cmesh.Block().SetNBlocks(nbl);
    
  for(i=0; i<nelem; i++) {
    TPZGeoEl *gel = elvec[i];
    if(!gel || gel->Reference()) continue;
    if(!gel->HasSubElement()) {
      int matid = gel->MaterialId();
      TPZMaterial * mat = cmesh.FindMaterial(matid);
      if(!mat)
        {
          matnotfound.insert(matid);
          continue;
        }
      int printing = 0;
      if (printing) {
        gel->Print(std::cout);
      }
            
      //checking material in MaterialIDs
      std::set<int>::const_iterator found = MaterialIDs.find(matid);
      if (found == MaterialIDs.end())
        {
          continue;
        }
            
      if(!gel->Reference() && gel->NumInterfaces() == 0)
        {
          index = CreateCompEl(gel,cmesh)->Index();
          if (fCreateHybridMesh) {
            cmesh.ElementVec()[index]->Reference()->ResetReference();
          }
#ifdef PZ_LOG
          if (logger.isDebugEnabled())
            {
                    
              std::stringstream sout;
              if (index < 0) {
                if(gel->Dimension() == 0){
                  sout << "Zero dimensional element, is your approximation space Hdiv?. " << std::endl;
                  gel->Print(sout);
                  sout << "No computational element was created. " << std::endl;
                }
              }else{
                cmesh.ElementVec()[index]->Print(sout);
              }
              LOGPZ_DEBUG(logger, sout.str())
                }
#endif
        }
    }
  }
  cmesh.InitializeBlock();
}
