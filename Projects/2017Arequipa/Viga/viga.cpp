#include "viga.h"

const int materialId = 1;
const int apoioMovelId = -2;
const int apoioFixoId = -22;
const int cargaMatId = +7;
const int prismaApoios = +17;

TPZGeoMesh * GeoMesh();
TPZCompMesh * CompMesh(TPZGeoMesh * gmesh, int pOrder);
void IncidHexa(int nx, int ny, int iel, TPZVec<int64_t> &incid);
int CreateMidNode(int nodeA, int nodeB, TPZGeoMesh * gmesh, REAL Yoffset);

int main() {

  TPZMaterial::gBigNumber = 1e12;
  const int pOrder = 1;
  const int nthreads = 0;//0 eh serial

  TPZGeoMesh * gmesh = GeoMesh();
  TPZCompMesh * cmesh = CompMesh(gmesh,pOrder);

  std::ofstream myfile("SENBAnalysis.txt"); //single edge notched bend
  TPZPullOutTestAnalysis an(cmesh);

  TPZSkylineStructMatrix StrMatrix(an.Mesh());
  StrMatrix.SetNumThreads(nthreads);
  an.SetStructuralMatrix(StrMatrix);

  TPZStepSolver<STATE> stepP;
  stepP.SetDirect( ECholesky  );
  an.SetSolver(stepP);

  an.Solution().Zero();
  an.LoadSolution();

  int substeps = 100;

  bool linesearch = true;
  an.TimeSteps(myfile, 1e0,100,linesearch,substeps,"SENBTest");
  myfile.flush();

  delete cmesh;
  delete gmesh;

  return 0;
}//void

TPZGeoMesh * GeoMesh(){

  //malha nova
  const int nx = 17;
  //em mm
  REAL xcoord[nx] = {0.,2.,4.,6.,12., 20.,28.,60.,90.,120.,150.,180.,210.,245.,255.,300.,350.};
  const int ny = 23;
  REAL ycoord[ny] = {0,25,31.25,37.5,43.75,50,56.25,62.5,68.75,75,81.25,87.5,93.75,100,106.25,112.5,118.75,125,131.25,137.5,143.75,149.,150};

  TPZVec<REAL> xc, yc;
  xc.Resize(nx);
  for(int i = 0; i < nx; i++) xc[i] = xcoord[i];
  yc.Resize(ny);
  for(int i = 0; i < ny; i++) yc[i] = ycoord[i];

  TPZGeoMesh * gmesh = new TPZGeoMesh();

  //criando nos
  const int nnodesPlano = nx*ny;
  REAL zc[2] = {0,150};
  for(int iz = 0; iz < 2; iz++){
    for(int iy = 0; iy < ny; iy++){
      for(int ix = 0; ix < nx; ix++){
        const int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL,3> nodeCoord(3);
        nodeCoord[0] = xc[ix];
        nodeCoord[1] = yc[iy];
        nodeCoord[2] = zc[iz];
        const int nodeID = nodind;
        gmesh->NodeVec()[nodind].Initialize(nodeID,nodeCoord,*gmesh);
      }
    }
  }

  //criando hexas
  const int nelHexa = (nx-1)*(ny-1);
  for(int iel = 0; iel < nelHexa; iel++){
  	TPZManVector<int64_t,8> incid(8);
    int64_t index;
    IncidHexa(nx, ny, iel, incid);
    gmesh->CreateGeoElement(ECube,incid,materialId,index);
  }

  //apoio vertical - penultima e antepenultima fila de nos a direita
  {
    //criando nos extras para o prisma do apoio
    TPZManVector<int64_t,6> nodindPrism(6);
    nodindPrism[0] = nx - 2;
    nodindPrism[1] = (nx - 2) - 1;
    nodindPrism[2] = CreateMidNode(  nodindPrism[0],  nodindPrism[1], gmesh, -20. );
    nodindPrism[3] = nx - 2 + nnodesPlano;
    nodindPrism[4] = (nx - 2 + nnodesPlano) - 1;
    nodindPrism[5] = CreateMidNode(  nodindPrism[3],  nodindPrism[4], gmesh, -20. );
    int64_t index;
    gmesh->CreateGeoElement(EPrisma,nodindPrism,prismaApoios,index);
    TPZManVector<int64_t,2> lineInd(2);
    lineInd[0] = nodindPrism[2];
    lineInd[1] = nodindPrism[5];
    gmesh->CreateGeoElement(EOned,lineInd,apoioMovelId,index);
  }

  //carga aplicada
  {
    TPZManVector<int64_t,4> nodind(4);
    nodind[0] = nx*(ny-1);
    nodind[1] = nx*(ny-1) + nnodesPlano;
    nodind[2] = (nx*(ny-1) + nnodesPlano) + 1;
    nodind[3] = (nx*(ny-1)) + 1;
    int64_t index;
    gmesh->CreateGeoElement(EQuadrilateral,nodind,cargaMatId,index);
  }

	gmesh->BuildConnectivity();

  //criando a outra metade da viga

  const int nnodes = gmesh->NNodes();
  TPZVec<int> nodesRef(nnodes);
  TPZVec<REAL> coord(3);
  for(int in = 0; in < nnodes; in++){
    gmesh->NodeVec()[in].GetCoordinates(coord);
    bool ranhura = false;
    if( fabs(coord[0]) < 1e-3 && fabs(coord[1]) < 1e-3 ) ranhura = true;
    if( fabs(coord[0]) < 1e-3 && !ranhura){
      //nao sera duplicado, porque esta no plano x = 0, exceto o (0,0,z) que eh a ranhura aberta na viga
      nodesRef[in] = in;
    }
    else{
      coord[0] *= -1.;//espelhando no plano x = 0
      const int nodind = gmesh->NodeVec().AllocateNewElement();
      const int nodeID = nodind;
      gmesh->NodeVec()[nodind].Initialize(nodeID,coord,*gmesh);
      nodesRef[in] = nodind;
    }
  }//nos copiados

  const int nel = gmesh->NElements();
  TPZManVector<int64_t> incid;
  for(int iel = 0; iel < nel; iel++){
    TPZGeoEl * gel = gmesh->ElementVec()[iel];
    if(!gel) continue;
    const int nnodesOfEl = gel->NNodes();
    incid.Resize(nnodesOfEl);
    for(int ic = 0; ic < nnodesOfEl; ic++){
      incid[ic] = nodesRef[ gel->NodeIndex(ic) ];
    }
    int64_t index;
    int matid = gel->MaterialId();
    if(matid == apoioMovelId) matid = apoioFixoId;
    gmesh->CreateGeoElement(gel->Type(),incid,matid,index);
  }//elementos copiados

  gmesh->BuildConnectivity();

  return gmesh;
}//method

TPZCompMesh * CompMesh(TPZGeoMesh * gmesh, int pOrder){

  if(!gmesh) DebugStop();

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  TPZVec<REAL> nullforce(3,0.);
  const REAL poisson = 0.2;
  REAL Ey = 33488.;

  REAL fc, ft;
  fc = 40.;
  ft = 4.0654;

  //material elasticidade
  TPZElasticity3DGD * elastMat = new TPZElasticity3DGD(materialId, Ey, poisson, nullforce,0.,0.,0.);

  //material para apoios apenas
  cmesh->InsertMaterialObject(new TPZElasticity3DGD(prismaApoios, Ey, poisson, nullforce,0.,0.,0.));

  TPZMaterial* material = elastMat;
  cmesh->InsertMaterialObject(material);

  {//apoio vertical movel
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    val1.Zero();
    val1(1,1) = TPZMaterial::gBigNumber;
    val1(2,2) = TPZMaterial::gBigNumber;//para tirar mov corpo rigido
    val2.Zero();
    TPZMaterial* bcmat = material->CreateBC(material,apoioMovelId,2,val1,val2);
	bcmat->SetLinearContext(false);
    cmesh->InsertMaterialObject(bcmat);
  }
  {//apoio vertical fixo
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    val1.Zero();
    val1(0,0) = TPZMaterial::gBigNumber;
    val1(1,1) = TPZMaterial::gBigNumber;
    val1(2,2) = TPZMaterial::gBigNumber;//para tirar mov corpo rigido
    val2.Zero();
    TPZMaterial* bcmat = material->CreateBC(material,apoioFixoId,2,val1,val2);
	bcmat->SetLinearContext(false);
    cmesh->InsertMaterialObject(bcmat);
  }

  REAL deslocamentoVertical = -30.;// em mm

  {//deslocamento
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    val1.Zero();
    val1(1,1) = TPZMaterial::gBigNumber;
    val2.Zero();
    val2(1,0) = deslocamentoVertical;
    TPZMaterial* bcmatDeslocamento = material->CreateBC(material,cargaMatId,2,val1,val2);
	bcmatDeslocamento->SetLinearContext(false);
    cmesh->InsertMaterialObject(bcmatDeslocamento);
  }

	cmesh->SetDefaultOrder( pOrder );
  cmesh->SetDimModel(3);

  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();

	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();

  return cmesh;

}///method

void IncidHexa(int nx, int ny, int iel, TPZVec<int64_t> &incid){
  const int linha = iel/(nx-1);//linhas debaixo pra cima
  const int coluna = iel % (nx-1);//colunas da esquerda pra direita
  const int noinicial = linha*nx+coluna;

  int nos[4] = {noinicial,noinicial+1,noinicial+1+nx,noinicial+nx};
  incid.Resize(8);
  for(int i = 0; i < 4; i++){
    incid[i] = nos[i];
  }
  const int nnodesPlano = nx*ny;
  for(int i = 0; i < 4; i++){
    incid[i+4] = nos[i]+nnodesPlano;
  }
}

int CreateMidNode(int nodeA, int nodeB, TPZGeoMesh * gmesh, REAL Yoffset){
  TPZManVector<REAL,3> coordA(3), coordB(3);
  gmesh->NodeVec()[ nodeA ].GetCoordinates(coordA);
  gmesh->NodeVec()[ nodeB ].GetCoordinates(coordB);
  TPZManVector<REAL,3> nodeCoord(3);
  for(int i = 0; i < 3; i++) nodeCoord[i] = (coordA[i]+coordB[i])/2.;
  nodeCoord[1] += Yoffset;

  const int nodind = gmesh->NodeVec().AllocateNewElement();
  const int nodeID = nodind;
  gmesh->NodeVec()[nodind].Initialize(nodeID,nodeCoord,*gmesh);
  return nodind;
}




