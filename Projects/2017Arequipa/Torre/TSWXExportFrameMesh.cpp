#include "TSWXExportFrameMesh.h"
#include "TPZEulerBernoulliBeam.h"

void TSWXExportFrameMesh::GenerateGraphMesh(TPZCompMesh &cmesh, double time){

  //nodal solution displacements x,y,z
  TSWXGraphSol displSol;
  displSol.fSolType = ENodeSolution;
  displSol.fTitle = "Displacements";

  //normal forces
  TSWXGraphSol normalForces;
  normalForces.fSolType = ECellSolution;
  normalForces.fTitle = "NormalForce";

  //section ids
  TSWXGraphSol sectionIds;
  sectionIds.fSolType = ECellSolution;
  sectionIds.fTitle = "SectionId";

  for(int iel = 0; iel < cmesh.NElements(); iel++){
    TPZEulerBernoulliBeam * beam = dynamic_cast<TPZEulerBernoulliBeam*>(cmesh.ElementVec()[iel]);
    if(!beam) continue;
    TPZGeoEl * gel = beam->Reference();
    if(gel->Dimension() != 1) continue;

    const int nnodesEl = 3;
    TSWXGraphNode nodeI(3), nodeJ(3);
    for(int i = 0; i < 3; i++) nodeI[i] = gel->NodePtr(0)->Coord(i);
    for(int i = 0; i < 3; i++) nodeJ[i] = gel->NodePtr(1)->Coord(i);
    TSWXGraphNode nodeMiddle(3);
    for(int i = 0; i < 3; i++) nodeMiddle[i] = (gel->NodePtr(0)->Coord(i)+gel->NodePtr(1)->Coord(i))/2.;
    const int Node0 = this->fGraphMesh.fNodes.size();
    this->fGraphMesh.fNodes.push_back(nodeI);
    this->fGraphMesh.fNodes.push_back(nodeMiddle);
    this->fGraphMesh.fNodes.push_back(nodeJ);

    TSWXGraphEl locEl;
    locEl.fElType = EvtkQuadraticEdge;
    locEl.fMatId = beam->SectionId();
    locEl.fIncid.resize(nnodesEl);
    for(int i = 0; i < 3; i++){
      locEl.fIncid[i] = Node0 + i;
    }
//    const int elementId = this->fGraphMesh.fElem.size();
    this->fGraphMesh.fElem.push_back(locEl);

    //nodal solution displacements x,y,z
    TPZFNMatrix< 12 > NodalDispl(12,1,0.);
    beam->GetSolutionVector(NodalDispl);
    TSWXGraphSingleSol u(3);
    //node I
    for(int i = 0; i < 3; i++){
      u[i] = NodalDispl(i,0);
    }
    displSol.fData.push_back(u);

    //node J
    for(int i = 0; i < 3; i++){
      u[i] = (NodalDispl(i,0)+NodalDispl(6+i,0))/2.; //we could improve visualization by actually computing nodeMiddle solution values by Hermitian interpolation
    }
    displSol.fData.push_back(u);

    //middle node
    for(int i = 0; i < 3; i++){
      u[i] = NodalDispl(6+i,0);
    }
    displSol.fData.push_back(u);

    //normal force
    TSWXGraphSingleSol N(1);
    N[0] = beam->GetStaticNormalForce();
    normalForces.fData.push_back( N );

    //section ids
    TSWXGraphSingleSol secId(1);
    secId[0] = beam->SectionId();
    sectionIds.fData.push_back(secId);
  }///for iel

  //solution data
  TPZStack< TSWXGraphSol > sol;
  sol.push_back(displSol);
  sol.push_back(normalForces);
  sol.push_back(sectionIds);
  this->fGraphMesh.fSol.AddSol(time,sol);

}//void




