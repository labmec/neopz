#include "TPZMHMH1ApproxCreator.h"
#include "TPZInterfaceEl.h"
#include "TPZNullMaterial.h"
#include "TPZCompElDisc.h"
#include "pzgeoelbc.h"
#include "Projection/TPZL2Projection.h"
#include "pzsubcmesh.h"
#include "TPZMultiphysicsCompMesh.h"

TPZMHMH1ApproxCreator::TPZMHMH1ApproxCreator(/* args */) : TPZMHMApproxCreator(), TPZH1ApproxCreator()
{
}

TPZMHMH1ApproxCreator::TPZMHMH1ApproxCreator(TPZVec<int64_t>& elementPartition, TPZGeoMesh * gmesh) : TPZMHMApproxCreator(elementPartition), TPZH1ApproxCreator(gmesh)
{
  SetHybridType(HybridizationType::EStandard);
}

TPZMHMH1ApproxCreator::~TPZMHMH1ApproxCreator()
{
}

TPZCompMesh *TPZMHMH1ApproxCreator::CreateL2Space()
{
  fGeoMesh->ResetReference();
  TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
  InsertL2MaterialObjects(cmesh);
  cmesh->SetDefaultOrder(fDefaultPOrder+fExtraInternalPOrder);
  cmesh->SetDimModel(fGeoMesh->Dimension());
  cmesh->SetAllCreateFunctionsContinuous();
  const std::map<int64_t, std::set<int64_t>> &domainMap = fDomainMap;
  int64_t nel = fGeoMesh->NElements();
  // create the computational elements for each domain
  for(auto &it : domainMap)
  {
    int64_t domain = it.first;
    auto &elements = it.second;
    for (auto el : elements)
    {
      TPZGeoEl *gel = fGeoMesh->Element(el);
      if(gel->HasSubElement())
      {
        continue;
      }
      int matid = gel->MaterialId();
      if(cmesh->FindMaterial(matid) == nullptr)
      {
        continue;
      }
      cmesh->ApproxSpace().CreateCompEl(gel,*cmesh);
    }
    for (auto el : elements)
    {
      TPZGeoEl *gel = fGeoMesh->Element(el);
      gel->ResetReference();
    }
  }
  cmesh->InitializeBlock();
  return cmesh;
}

void TPZMHMH1ApproxCreator::GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh)
{
  std::map<int64_t, TPZSubCompMesh *> submeshes;
  for(auto &it : fDomainMap)
  {
    int64_t domain = it.first;
    TPZSubCompMesh *submesh = new TPZSubCompMesh(*mcmesh);
    submeshes[domain] = submesh;
  }
  mcmesh->LoadReferences();
  {
    std::ofstream out("gmesh.txt");
    this->fGeoMesh->Print(out);
  }
  for(auto &it : fDomainMap)
  {
    int64_t domain = it.first;
    auto &elements = it.second;
    TPZSubCompMesh *submesh = submeshes[domain];
    for(auto el : elements)
    {
      TPZGeoEl *gel = fGeoMesh->Element(el);
      TPZCompEl *cel = gel->Reference();
      if(!cel) continue;
      submesh->TransferElement(mcmesh, cel->Index());
    }
    submesh->InitializeBlock();
  }
  mcmesh->ComputeNodElCon();
  int64_t ncon = mcmesh->NConnects();
  for(int64_t icon = 0; icon < ncon; icon++)
  {
    TPZConnect &c = mcmesh->ConnectVec()[icon];
    int lagrangeLevel = c.LagrangeMultiplier();
    if(lagrangeLevel >= 3) c.IncrementElConnected();
  }
  for(auto &it : submeshes)
  {
    TPZSubCompMesh *submesh = it.second;
    submesh->MakeAllInternal();
    submesh->CleanUpUnconnectedNodes();
    submesh->SetAnalysisSparse(0);

  }
  mcmesh->ComputeNodElCon();
  mcmesh->CleanUpUnconnectedNodes();
#ifdef PZDEBUG
  {
    std::ofstream out("mcmesh.txt");
    mcmesh->Print(out);
  }
#endif
}

TPZCompMesh *TPZMHMH1ApproxCreator::CreateBoundaryHDivSpace()
{
  return TPZH1ApproxCreator::CreateBoundaryHDivSpace();
}

TPZCompMesh *TPZMHMH1ApproxCreator::CreateConstantSpace(const int &lagMult)
{
  TPZCompMesh *constMesh = new TPZCompMesh(fGeoMesh);
  constMesh->SetAllCreateFunctionsDiscontinuous();
  constMesh->SetDefaultOrder(0);
  int nstate = 1;
  int meshDim = fGeoMesh->Dimension();
  if(fProbType == ProblemType::EElastic) {
      nstate = meshDim*(meshDim+1)/2;
  }

  std::set<int> volumetricMatIds = GetVolumetricMatIds();
  for(auto it : volumetricMatIds){
      TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim,nstate);
      constMesh->InsertMaterialObject(mat);
  }

  const std::map<int64_t,std::set<int64_t>> &domainMap = fDomainMap;
  int64_t nel = fGeoMesh->NElements();
  for(auto &it : domainMap)
  {
    int64_t domain = it.first;
    auto &elements = it.second;
    REAL DomainArea = 0.;
    TPZManVector<REAL,3> DomainCenter(3,0.);
    for (auto el : elements)
    {
      TPZGeoEl *gel = fGeoMesh->Element(el);
      int matid = gel->MaterialId();
      if(volumetricMatIds.find(matid) == volumetricMatIds.end())
      {
        continue;
      }
      TPZGeoElSide gelside(gel);
      REAL area = gelside.Area();
      TPZManVector<REAL,3> xcenter(3,0.);
      gelside.CenterX(xcenter);
      DomainArea += area;
      for(int i = 0; i < 3; i++)
      {
        DomainCenter[i] += area*xcenter[i];
      }
    }
    for(int i = 0; i < 3; i++)
    {
      DomainCenter[i] /= DomainArea;
    }
    int64_t domainconnect = -1;
    for (auto el : elements)
    {
      TPZGeoEl *gel = fGeoMesh->Element(el);
      int matid = gel->MaterialId();
      if(volumetricMatIds.find(matid) == volumetricMatIds.end())
      {
        continue;
      }
      auto *cel = constMesh->ApproxSpace().CreateCompEl(gel,*constMesh);
      TPZCompElDisc *discel = dynamic_cast<TPZCompElDisc *>(cel);
      if(!discel)
      {
        DebugStop();
      }
      if(domainconnect == -1)
      {
        domainconnect = discel->ConnectIndex(0);
        TPZConnect &c = constMesh->ConnectVec()[domainconnect];
        c.SetLagrangeMultiplier(lagMult);
      }
      else
      {
        discel->SetConnectIndex(0,domainconnect);
      }
      discel->SetCenterPoint(DomainCenter);
      discel->SetFalseUseQsiEta();
      gel->ResetReference();
    }
  }
  constMesh->InitializeBlock();
  constMesh->ComputeNodElCon();
  constMesh->CleanUpUnconnectedNodes();
  return constMesh;
}

void TPZMHMH1ApproxCreator::AddHybridizationGeoElements()
{
  CreateSkeletonElements();
}

void TPZMHMH1ApproxCreator::CreateSkeletonElements()
{
  if (fHybridizationData.fInterfaces.size())
  {
    DebugStop();
  }

  if (fHybridizationData.fLagrangeMatId < 0)
    DebugStop();

  std::set<int64_t> lagrangeElements;
  auto *cmesh = CriaMalhaTemporaria();
  int64_t nel = fGeoMesh->NElements();
  int dimension = fGeoMesh->Dimension();
  int ninterf;

  for (int64_t iel = 0; iel < nel; iel++)
  {
    TPZGeoEl *gel = fGeoMesh->ElementVec()[iel];
    if (!gel)
      continue;

    ninterf = gel->NumInterfaces();
    if (ninterf > 1)
      DebugStop();
    if (ninterf == 1)
    {
      TPZCompEl *cel = gel->Reference();
      TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
      if (!intface)
      {
        DebugStop();
      }
      lagrangeElements.insert(iel);
      gel->SetMaterialId(fHybridizationData.fLagrangeMatId);
      gel->DecrementNumInterfaces();
    }
  }

  fGeoMesh->ResetReference();
  delete cmesh;

  CreateWrapElements();
  ComputeDomainMap();
}

void TPZMHMH1ApproxCreator::CreateWrapElements()
{
  int64_t nel = fGeoMesh->NElements();
  int meshdim = fGeoMesh->Dimension();
  std::set<int> bcMatIds = GetBCMatIds();
  SetSubelementDomain(fGeoMesh);
  for (int64_t el = 0; el < nel; el++)
  {
    TPZGeoEl *gel = fGeoMesh->Element(el);
    if (!gel || gel->Dimension() != meshdim || gel->HasSubElement())
    {
      continue;
    }
    int64_t domain = -1;
    if (el < fGeoToMHMDomain.size())
    {
      domain = fGeoToMHMDomain[el];
    }
    if (domain == -1)
    {
      DebugStop();
    }
    int firstside = gel->FirstSide(meshdim - 1);
    int lastside = gel->FirstSide(meshdim);
    for (int side = firstside; side < lastside; side++)
    {
      TPZGeoElSide gelside(gel, side);
      TPZGeoElSide boundary = gelside.HasNeighbour(bcMatIds);
      if(boundary && fHybridizationData.fHybridizeBCLevel == 0) {
        SetMHMDomain(boundary.Element(), domain);
        continue;
      }
      if (boundary)
      {
        TPZGeoElBC gbc(gelside, fHybridizationData.fWrapMatId);
        TPZGeoEl *wrap = gbc.CreatedElement();
        SetMHMDomain(wrap, domain);
        int interfaceMatId = fHybridizationData.fLeftInterfaceMatId;
        TPZGeoElBC gbc2(gbc, interfaceMatId);
        SetMHMDomain(gbc2.CreatedElement(), domain);
        int64_t interfaceIndex = gbc2.CreatedElement()->Index();
        fHybridizationData.fInterfaces[interfaceIndex] = TPZHybrid(interfaceIndex, wrap->Index(), boundary.Element()->Index());
        continue;
      }
      // if there is no boundary, there must be a neighbour of the same dimension.
      TPZGeoElSide neighbour = gelside.Neighbour();
      while (neighbour != gelside && neighbour.Element()->Dimension() != meshdim)
      {
        neighbour = neighbour.Neighbour();
      }
      if (neighbour == gelside)
      {
        DebugStop();
      }
      int64_t neighbourIndex = neighbour.Element()->Index();
      int64_t neighbourDomain = -1;
      if (neighbourIndex < fGeoToMHMDomain.size())
      {
        neighbourDomain = fGeoToMHMDomain[neighbourIndex];
      }
      // if the neighbour is in a different domain, we need to create a wrap element and an interface element
      if (neighbourDomain != domain)
      {
        TPZGeoElBC gbc(gelside, fHybridizationData.fWrapMatId);
        TPZGeoEl *wrap = gbc.CreatedElement();
        SetMHMDomain(wrap, domain);

        TPZGeoElSide skelside = gelside.HasNeighbour(fHybridizationData.fLagrangeMatId);
        if (!skelside)
        {
          skelside = gelside.HasLowerLevelNeighbour(fHybridizationData.fLagrangeMatId);
        }
        if (!skelside)
        {
          DebugStop();
        }
        if(neighbourDomain < domain) {
          int interfaceMatId = fHybridizationData.fLeftInterfaceMatId;
          TPZGeoElBC gbc2(gelside, interfaceMatId);
          SetMHMDomain(gbc2.CreatedElement(), domain);
          int64_t interfaceIndex = gbc2.CreatedElement()->Index();
          fHybridizationData.fInterfaces[interfaceIndex] = TPZHybrid(interfaceIndex, wrap->Index(), skelside.Element()->Index());
        } else {
          int interfaceMatId = fHybridizationData.fRightInterfaceMatId;
          TPZGeoElBC gbc2(gelside, interfaceMatId);
          SetMHMDomain(gbc2.CreatedElement(), domain);
          int64_t interfaceIndex = gbc2.CreatedElement()->Index();
          fHybridizationData.fInterfaces[interfaceIndex] = TPZHybrid(interfaceIndex, skelside.Element()->Index(), wrap->Index());
        }
      }
    }
  }
}

TPZCompMesh *TPZMHMH1ApproxCreator::CriaMalhaTemporaria()
{
  fGeoMesh->ResetReference();
  int dim = fGeoMesh->Dimension();
  std::set<int> matids, bcids;
  TPZGeoMesh &gmesh = *fGeoMesh;
  int64_t nel = gmesh.NElements();
  auto belongsToMHM = [this](int64_t el) {
    return el < fGeoToMHMDomain.size() && fGeoToMHMDomain[el] != -1;
  };

  for (int64_t el = 0; el < nel; el++)
  {
    TPZGeoEl *gel = gmesh.ElementVec()[el];
    if (!gel)
    {
      continue;
    }
    if (gel->Dimension() != dim || !belongsToMHM(el))
    {
      continue;
    }
    int materialid = gel->MaterialId();
    matids.insert(materialid);
    int ns = gel->NSides();
    for (int is = 0; is < ns; is++)
    {
      if (gel->SideDimension(is) != dim - 1)
      {
        continue;
      }
      TPZGeoElSide gelside(gel, is);
      TPZGeoElSide neighbour = gelside.Neighbour();
      while (neighbour != gelside)
      {
        if (neighbour.Element()->Dimension() == dim - 1)
        {
          bcids.insert(neighbour.Element()->MaterialId());
        }
        if (neighbour.Element()->Dimension() == dim)
        {
          int neighmatid = neighbour.Element()->MaterialId();
          if (materialid != neighmatid && materialid < neighmatid)
          {
            fGeoMesh->AddInterfaceMaterial(materialid, neighmatid, materialid);
            fGeoMesh->AddInterfaceMaterial(neighmatid, materialid, materialid);
          }
        }
        neighbour = neighbour.Neighbour();
      }
    }
  }

  TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
  cmesh->SetDimModel(dim);
  TPZManVector<STATE, 1> sol(1, 0.);
  int nstate = 1;
  std::set<int>::iterator it = matids.begin();
  TPZMaterialT<STATE> *meshmat = nullptr;
  while (it != matids.end())
  {
    TPZL2Projection<STATE> *material = new TPZL2Projection<STATE>(*it, dim, nstate, sol);
    cmesh->InsertMaterialObject(material);
    if (!meshmat)
    {
      meshmat = material;
    }
    it++;
  }

  if (!meshmat)
  {
    DebugStop();
  }
  {
    TPZL2Projection<STATE> *material = new TPZL2Projection<STATE>(fHybridizationData.fWrapMatId, dim, nstate, sol);
    cmesh->InsertMaterialObject(material);
  }
  if(fHybridizationData.fHybridizeBCLevel != 0)
  {
    it = bcids.begin();
    while (it != bcids.end())
    {
      TPZFMatrix<STATE> val1(2, 2, 0.);
      TPZVec<STATE> val2(2, 0.);
      TPZBndCond *BCondD1 = meshmat->CreateBC(meshmat, *it, 0, val1, val2);
      cmesh->InsertMaterialObject(BCondD1);
      it++;
    }
  }
  cmesh->SetAllCreateFunctionsDiscontinuous();
  TPZCreateApproximationSpace &create = cmesh->ApproxSpace();

  for (int64_t el = 0; el < nel; el++)
  {
    TPZGeoEl *gel = gmesh.ElementVec()[el];
    if (!gel)
    {
      continue;
    }
    if (gel->Dimension() != dim || !belongsToMHM(el))
    {
      continue;
    }
    create.CreateCompEl(gel, *cmesh);
    if(fHybridizationData.fHybridizeBCLevel == 0)
    {
      continue;
    }
    int ns = gel->NSides();
    for (int is = 0; is < ns; is++)
    {
      if (gel->SideDimension(is) != dim - 1)
      {
        continue;
      }
      TPZGeoElSide gelside(gel, is);
      TPZGeoElSide neighbour = gelside.Neighbour();
      while (neighbour != gelside)
      {
        if (neighbour.Element()->Dimension() == dim - 1)
        {
          int matid = neighbour.Element()->MaterialId();
          if (bcids.find(matid) != bcids.end())
          {
            create.CreateCompEl(neighbour.Element(), *cmesh);
          }
        }
        neighbour = neighbour.Neighbour();
      }
    }
  }

  TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
  cmesh->ExpandSolution();
  return cmesh;
}
