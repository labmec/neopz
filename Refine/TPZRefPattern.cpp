/* Generated by Together */

#include "TPZRefPattern.h"
//#include "pzgmesh.h"
#include "pztrnsform.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pzelg1d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
#include "pzelgt3d.h"
#include "pzelgpi3d.h"
#include "pzelgpr3d.h"
#include "pzelgc3d.h"
#include "pzquad.h"
#include "pzvec.h"
//#include "pzvecc.h"

#include <fstream>

TPZRefPattern::TPZRefPattern() : fSubElRefPattern(0),fBCRefPattern(0) {
  fFileRefPatt = '\0';
  fRefineType = -1;
  fMesh = 0;
  fBCRefPattern.Resize(0);
  fName = "";
}

TPZRefPattern::TPZRefPattern(char * file ) : fSubElRefPattern(0),fBCRefPattern(0) {
  fMesh = new TPZGeoMesh;
  fMesh->SetName("* * * Refinement Pattern * * *");
  fFileRefPatt = file;/**arquivo contendo o padr�o de refinamento*/
  ReadPattern();/**l� o arquivo contendo o refinamento e cria a malha fMesh*/
  fMesh->BuildConnectivity();/**conectividades entre sub-elementos*/
//  MeshPrint();//TETSTE
  ComputeTransforms();/**calcula as transforma��es entre filhos e pai*/
  ComputePartition();/**efetua a parti��o do elemento pai de acordo com os lados dos sub-elementos*/
  fBCRefPattern.Resize(0);
}

void TPZRefPattern::ReadPattern(){
  ifstream in(fFileRefPatt);
  int nnodes,nelems,inode;
  in >> nnodes >> nelems;
  in >> fRefineType;
  in >> fName;
  TPZVec<REAL> coord(3);
  fMesh->NodeVec().Resize(nnodes);
  //criac�o dos objetos nodais
  for(inode=0;inode<nnodes;inode++){
    in >> coord[0];
    in >> coord[1];
    in >> coord[2];
//    cout << coord << endl;
    fMesh->NodeVec()[inode].Initialize(inode,coord,*fMesh);
  }
  TPZGeoEl *father;
  //criac�o dos elementos geom�tricos que definem a partic�o
  int ntype,nummat,ncorners,incid,el;
  for(el=0;el<nelems;el++) {/**os sub-elementos podem n�o ter de uma mesma geometria*/
    in >> ntype >> nummat;
    ncorners = ntype;
    if(ntype == 7) ncorners = 4;/**tetraedro*/
    TPZVec<int> nodes(ncorners);
    for(incid=0;incid<ncorners;incid++){
      in >> nodes[incid];
    }
    TPZGeoEl *subel = CreateGeoEl(ntype,nummat,nodes,fMesh);
    if(el==0){
      father = subel;
      //father->ResizeSubEl(nelems-1);//n�mero de filhos
    }
    if(el>0){
      father->SetSubElement(el-1,subel);
      subel->SetFather(father);
    }
  }
 // MeshPrint();
}

/** Preenche-se a estrutura elemento/lado com aqueles objetos que tem
 *  transforma��o associada.Calcula-se a transforma��o entre o lado
 *  do sub-elemento e o lado do pai respectivo.
 */
void TPZRefPattern::ComputePartition(){
  int sizeinit = fTransforms.fInitSonSides.NElements()-1;/**igual ao n�mero de filhos*/
  int init,iside;
  TPZGeoEl *fat = Element(0);/**elemento pai da divis�o*/
  DefinitionOfSizePartition();/**calcula o tamanho da partic�o dos lados do pai e inicializa fInitSide*/
  int sidestot = fFatherSides.fPartitionSubSide.NElements();
  /**inicializando o vetor de transformac�es*/
  for(int p=0;p<sidestot;p++) fFatherSides.fPartitionSubSide[p] = TPZGeoElSide();/**zera as entradas*/
  /**calcular a parti��o*/
  for(init=0;init<sizeinit;init++){/**percorre os filhos*/
    int initsideson = fTransforms.fInitSonSides[init];/**come�o dos lados do filho*/
//    int lastsideson = fTransforms.fInitSonSides[init+1];/**comeco dos lados do pr�ximo filho*/
//    int nsonsd = lastsideson-initsideson;/**comprimento da partic�o*/
    TPZGeoEl *son = Element(init+1);/**filho*/
    int nsides = son->NSides();
    for(iside=0;iside<nsides;iside++){/**percorre-se os lados do filho*/
      int sdf = fTransforms.fFatherSide[initsideson+iside];/**lado do pai associado ao lado do sub-elemento*/
      if(sdf < fat->NNodes()) continue;/**cantos do pai n�o entram na parti��o do lado*/
      int pos = 0;
      int sdf2 = sdf - fat->NNodes();/**fInitSide n�o contempla os cantos, s� a partir de aresta para cima*/
      int initss = fFatherSides.fInitSide[sdf2];/**inicio da partic�o do lado sdf do pai*/
      TPZGeoElSide sonside = fFatherSides.fPartitionSubSide[initss+pos];/**1o sub/lado da parti��o do lado do pai*/
      TPZGeoElSide empty(sonside);/**copia de sonside*/
      int sizept = fFatherSides.fInitSide[sdf2+1]-initss;/**comprimento da faixa do lado*/
      while(empty.Element() && pos < sizept){/**percorre-se a faixa do lado do pai*/
        pos++;/**procura-se o primeiro lugar desocupado na faixa do lado do pai*/
        empty = fFatherSides.fPartitionSubSide[initss+pos];
      }
      if(pos > sizept-1){/**a faixa toda esta ocupada com sub/lado n�o nulo*/
        PZError << "TPZRefPattern::ComputePartition erro de dimensionamento\n";
        //exit(-1);
      }
      fFatherSides.fPartitionSubSide[initss+pos] = TPZGeoElSide(son,iside);
    }
  }
  /**conferindo a consist�ncia da partic�o*/
  for(int ss=0;ss<sidestot;ss++){
    if(!fFatherSides.fPartitionSubSide[ss].Element()){
      PZError << "TPZRefPattern::ComputePartition particao inconsistente";
    }
  }
  /**saida do arquivo de dados da particao*/
  ofstream out("partition.out");
  fFatherSides.Print(*fMesh,out);
  /**extraindo sub/side quando o side � repetido dentro da partic�o do lado*/
  int init2;
  sizeinit = fFatherSides.fInitSide.NElements()-1;
  for(iside=0;iside<sizeinit;iside++){/**percorrendo fInitSide*/
    init = fFatherSides.fInitSide[iside];/**posic�o inicial*/
    init2 = fFatherSides.fInitSide[iside+1];/**posic�o final*/
    for(int sd=init;sd<init2;sd++){/**percorrendo a partic�o do lado do pai*/
      TPZGeoElSide gs = fFatherSides.fPartitionSubSide[sd];/**1o sub/side do lado*/
      if(!gs.Element()) continue;
      for(int sd2=sd+1;sd2<init2;sd2++){
        TPZGeoElSide gs2 = fFatherSides.fPartitionSubSide[sd2];/**sub/side sub-seguinte do lado*/
        if(gs2.Element() && gs.Element()->NeighbourExists(gs.Side(),gs2)){
           fFatherSides.fPartitionSubSide[sd2] = TPZGeoElSide();/**apagando o sub/side com side repetido*/
        }
      }
    }
  }
  /**tirando os buracos da partic�o*/
  TPZVec<int> newinit(fFatherSides.fInitSide.NElements());/**capacidade m�xima*/
  TPZVec<TPZGeoElSide> newpartition(fFatherSides.fPartitionSubSide.NElements());/**capacidade m�xima*/
  newinit[0] = 0;
  int count = 0;
  for(iside=0;iside<sizeinit;iside++){/**percorrendo fInitSide*/
    init = fFatherSides.fInitSide[iside];/**posic�o inicial*/
    init2 = fFatherSides.fInitSide[iside+1];/**posic�o final*/
    newinit[iside+1] = newinit[iside];
    for(int sd=init;sd<init2;sd++){/**percorrendo a partic�o do lado do pai*/
      TPZGeoElSide gs = fFatherSides.fPartitionSubSide[sd];/**1o sub/side do lado*/
      if(!gs.Element()) continue;
      newinit[iside+1]++;
      newpartition[count++] = gs;
    }
  }
  fFatherSides.fInitSide = newinit;
  fFatherSides.fPartitionSubSide = newpartition;
  NSideSubElements();/**preenche fNSubSideFather com o n�mero de elementos associados a cada lado do pai*/
  out << "\n\n               *** Particao enxuta ***\n\n";
  fFatherSides.Print(*fMesh,out);
  out.flush();
  out.close();
  cout << "\nTPZRefPattern::ComputePartition arquivo contendo a particao dos lados do pai: partition.out\n";
  ofstream out1("fathersides.out");
  Print1(*fMesh,out1);
  out1.flush();
  out1.close();
  cout << "\nTPZRefPattern::ComputePartition arquivo contendo o lado do pai associado ao lado do filho: fathersides.out\n";
}

void TPZRefPattern::TPZPartitionFatherSides::Print(TPZGeoMesh &gmesh,ostream &out){
  out << "TPZRefPattern::TPZSideTransform::Print partic�o dos lados do pai pelos lados dos sub-elementos\n\n";
  int iss,iside;//,count=0;
  TPZGeoEl *father = gmesh.ElementVec()[0];/**elemento pai*/
  int nsfat = father->NSides();
  int nnod = father->NNodes();
  int ntot = nsfat-nnod;
  for(iside=0;iside<ntot;iside++){
    out << "\nLado do elemento pai = " << (iside+nnod) << endl;
    int initsideson = fInitSide[iside];
    int lastsideson = fInitSide[iside+1];
    for(iss=initsideson;iss<lastsideson;iss++){/**percorre-se os sub/lado da parti��o*/
      TPZGeoElSide subside = fPartitionSubSide[iss];
      if(!subside.Element()){
        cout <<  "ERRO : Elemento da partic�o nulo\n ";
        continue;
      }
      out << "Sub id = " << subside.Element()->Id() << "  Lado do sub " << subside.Side()  << endl;
    }
    out << endl;  
  }  
}

void TPZRefPattern::Print1(TPZGeoMesh &gmesh,ostream &out){
  out << "TPZRefPattern::TPZSideTransform::Print lado do pai associados aos lados dos sub-elementos\n\n";
  int iside,isub;
//  int nnod = Element(0)->NNodes();
  int nsubs = Element(0)->NSubElements();
  out << "Refinement Pattern id " << fRefineType << " named " << fName << endl;
  for(isub=0;isub<nsubs;isub++){
    TPZGeoEl *sub = Element(isub+1);
    int nsides = sub->NSides();
    for(iside=0;iside<nsides;iside++){
      out << "sub/lado = " << isub << "/" << iside << "  Lado do pai = " << FatherSide(iside,isub)  << endl;
    }
  }
  this->fFatherSides.Print(gmesh,out);
  this->fTransforms.Print(gmesh,out);
}

/** Determina-se, para todos e cada um dos lados dos sub-elementos, o lado do pai no
 *  qual est� contido. Caso o lado do elemento � igual ao lado do pai asigna-se
 *  o valor -1. Neste caso n� existe depend�ncia. 
 */
void TPZRefPattern::ComputeTransforms(){
/**calcula transformacoes entre lado de filho e lado de pai*/
//  int ind = 0;
//  int sizemesh = NSubElements()+1;
  TPZGeoEl *fath = Element(0);/**Elemento pai deve ser o primeiro elemento da lista*/
  if(!fath){
    PZError << "TPZRefPattern::ComputePartition Father not exists?!\n";
    exit(-1);
  }
  int isub,nsubs = NSubElements();/**total filhos*/
  /**preenchendo a estrutura TPZFatherSides*/
  int initside=0,cont=0,side,fatside;  
  TPZVec<REAL> masscent(3,0.),xpoint(3,0.);
  xpoint[0] = 0.;
  xpoint[1] = 0.;
  xpoint[2] = 0.; 
  fTransforms.fInitSonSides.Resize(nsubs+1);/** +1 para incluir a posi��o final de fSideFather*/
  int size = SizeOfSubsSides(nsubs);
  fTransforms.fFatherSide.Resize(size);/**vale -1 quando o lado n�o � contido propriamente*/
  fTransforms.fSideTransform.Resize(size);
  for(isub=0;isub<nsubs;isub++){
    fTransforms.fInitSonSides[isub] = initside;
    TPZGeoEl *son = Element(isub+1);/**o pai n�o sabe quem s�o os filhos*/
    int nsides = son->NSides();/**um dos elementos implementados no PZ*/
    //cout << "\n------------------------> ELemento de id = " << son->Id() << endl << endl;
    for(side=0;side<nsides;side++){/**cantos + arestas + faces, o interior n�o � compartilhado*/
      son->CenterPoint(side,masscent);/**percorre todos os lados do elemento filho*/
      son->X(masscent,xpoint);
      //cout << "\nxpoint = " << xpoint[0] << " , " << xpoint[1] << " , " << xpoint[2] << endl << endl;
      fatside = fath->WhichSide(xpoint);/**lado do pai contendo o lado do filho*/
      fTransforms.fSideTransform[cont] = son->ComputeParamTrans(fath,fatside,side);
      //cout << "sub/side -> fath side: " << son->Id() << "/" << side << " -> " << fatside << endl; 
      //fTransforms.fSideTransform[cont].Mult().Print(" T ");
      //fTransforms.fSideTransform[cont].Sum().Print(" b ");
      fTransforms.fFatherSide[cont++] = fatside;
      initside++;
    }
  }
  fTransforms.fInitSonSides[isub] = initside;/**posi��o final em fSideFather*/
  ofstream out("transformacoes.out");
  fTransforms.Print(*fMesh,out);
  out.flush();
  out.close();
  cout << "\nTPZRefPattern::ComputeTransforms lados do pai asssociados aos subs, fathersides.out\n";
}

void TPZRefPattern::TPZSideTransform::Print(TPZGeoMesh &gmesh,ostream &out){
  out << "TPZRefPattern::TPZSideTransform::Print transformacoes parametricas\n\n";
  int isub,iside,count=0;
  int nsubs = gmesh.ElementVec().NElements()-1;/**a malha cont�m um elemento pai e filhos*/
  for(isub=0;isub<nsubs;isub++){
    int elid = gmesh.ElementVec()[isub+1]->Id();
    out << "\n>>>>>>>>>>>> Sub-element id = " << elid << " <<<<<<<<<<<< " << endl << endl;/**a informa��o da malha � est�tica*/
    int nsides = gmesh.ElementVec()[isub+1]->NSides();
    for(iside=0;iside<nsides;iside++){/**percorre-se os lados de cada filho*/
      out << "> Sub-element/side = " << elid << "/" << iside << "  Side of father = " << fFatherSide[count] << endl;
      fSideTransform[count].Mult().Print("Transformacao  T :",out);
      fSideTransform[count++].Sum().Print("Vetor b : ",out);
    }
    out << endl;  
  }
}

int TPZRefPattern::FatherSide(int side, int sub ){
  int nsides = fMesh->ElementVec()[sub+1]->NSides();
  int nsubs = NSubElements();
  if(side<0 || side>nsides ||  sub <0 ||  sub >nsubs){
    PZError << "TPZRefPattern::FatherSide arguments wrong argument\n";
    PZError << "side = " << side << " sub = " << sub;
    return -1;
  }
  int pos = fTransforms.fInitSonSides[sub];
  return ( fTransforms.fFatherSide[pos+side]  );
}

TPZTransform TPZRefPattern::Transform(int side, int sub){
  int nsides = Element(sub+1)->NSides();
  int nsubs = NSubElements();
  if(side<0 || side>nsides ||  sub <0 ||  sub >nsubs){
    PZError << "TPZRefPattern::Transform wrong arguments\n";
    PZError << "side = " << side << " sub = " << sub;
    return TPZTransform(0);
  }
  int pos = fTransforms.fInitSonSides[sub];
  return ( fTransforms.fSideTransform[pos+side]  );
}

void TPZRefPattern::SideNodes(int side, TPZVec<int> &vecnodes){
  /**side � um lado � do pai, � preciso achar todos os n�s internos ao lado*/
  TPZGeoEl *father = Element(0);
  int nsdfat = father->NSides();
  int nnod = father->NNodes();
  if(side<0 || side>nsdfat){
    PZError << "TPZRefPattern::SideNodes wrong side, side = " << side << endl;
    vecnodes.Resize(0);
    return;
  }
  if (side < nnod){
    vecnodes.Resize(1);
    vecnodes[0] = Element(0)->NodeIndex(side);
    return ;
  }
  //Aparentemente a estrutura nao contempla os nos...
  side -= nnod;

#ifdef HUGE_DEBUG
  cout << "************" << endl;
  cout << fFatherSides.fInitSide << endl;
  cout << "************" << endl;
  cout << fFatherSides.fNSubSideFather << endl;
  cout << "************" << endl;
  int nss = fFatherSides.fPartitionSubSide.NElements()-1;
  for (int i =0; i<nss ;i++){
    cout << "\nsubelement " << i << endl;
    cout << "side " << fFatherSides.fPartitionSubSide[i].Side() << endl;
    if (fFatherSides.fPartitionSubSide[i].Element())
      fFatherSides.fPartitionSubSide[i].Element()->Print(cout);
    else cout <<  "No associated element " << endl;
  }
#endif  
    
  int pos = fFatherSides.fInitSide[side];
  int pos2 = fFatherSides.fInitSide[side+1];
  vecnodes.Resize(pos2-pos);/**o n�mero de n�s no lado � menor que isto*/
  int count = 0;
  for(int par=pos;par<pos2;par++){/**intervalo do lado side*/
    TPZGeoElSide subs = fFatherSides.fPartitionSubSide[par];//[pos];
    if(!subs.Element()){
      PZError << "TPZRefPattern::NSideSubElements puncture in the partition\n";
      vecnodes.Resize(0);
      return;
    }
//    subs.Element()->Print(cout);
    int sd = subs.Side();
    if(sd < nnod){
      TPZGeoEl *el = subs.Element();
      int node = el->NodeIndex(sd);
      //vecnodes[par-pos] = sd;/**cada n� aparece uma �nica ves na partic�o do lado*/
      vecnodes[par-pos] = node;
      count++;
    }
  }
  vecnodes.Resize(count);
//  cout << "VecNodes " << vecnodes << endl;
}

int TPZRefPattern::NNodes(){
  return ( fMesh->NodeVec().NElements() );
}

int TPZRefPattern::NSideNodes(int side){
  int nsdfat = Element(0)->NSides();
  if(side<0 || side>nsdfat){
    PZError << "TPZRefPattern::NSideNodes wrong side, side = " << side << endl;
    return -1;
  }
  TPZVec<int> vec;
  SideNodes(side,vec);  
  return ( vec.NElements() );
}

int TPZRefPattern::NSubElements(){
  return ( fMesh->ElementVec().NElements() - 1 );
}

void TPZRefPattern::SideSubElement(int sidein, int position, int & sub, int & sideout){
  if(sidein<0 || sidein>Element(0)->NSides()){
    PZError << "TPZRefPattern::SideSubElement null side, side = " << sidein << endl;  
  }
  //Para contemplar os lados de dimensao 0
/*  if (side < Element(0).NCornerNodes()){
    TPZGeoElSide thisside (Element(0),sidein);
    TPZGeoElSide neighbour = thisside.Neighbour();
    while (neighbour != thisside && neighbour.Exists()){
      sub.Push (neighbour.Element().Id() -1);
      sideout.Push (neighbour.Side());
      neighbour = neighbour.Neighbour();
    }
    return;
  }*/
  int insd = fFatherSides.fInitSide[sidein];
  int insd2 = fFatherSides.fInitSide[sidein+1];
/*
  for (int i = 0; i<fFatherSides.fPartitionSubSide.NElements(); i++){
    if (!fFatherSides.fPartitionSubSide[i].Exists()) cout << "error for element side " << i << endl;
    else cout << "Fathersides.PartitionSubSide " << fFatherSides.fPartitionSubSide[i].Element()->Id() -1
         << " side " << fFatherSides.fPartitionSubSide[i].Side() << endl;
  }
*/  
  if(position < 0 || position > (insd2-insd)){/**subs2-subsd � o n�mero de sub's do lado*/
     PZError << "TPZRefPattern::SideSubElement wrong position, position = " << position << endl;
     sub = sideout = -1;
     return;
  }
  /**a ordem � determinada pela partic�o dos lados do pai, � uma ordem fixa*/
  sub = fFatherSides.fPartitionSubSide[insd].Element()->Id()-1;/**id contemplado como filho*/
  sideout = fFatherSides.fPartitionSubSide[insd].Side();
}

void TPZRefPattern::NSideSubElements(){
  /**procura-se o n�mero de sub-elementos da partic�o, isto n�o � igual ao n�mero de
     elementos da partic�o*/
  /**side � um lado do pai, � preciso achar o n�mero de sub-elementos ligados a este lado*/
  TPZGeoEl *father = Element(0);
  int nsdfat = father->NSides();
  int nnod = father->NNodes();
  fFatherSides.fNSubSideFather.Resize(nsdfat); 
  for(int side=0;side<nsdfat;side++){/**percorre-se os lados do pai*/
    if(side<nnod){/**os cantos n�o fazem parte da partic�o*/
      TPZGeoElSide fat(Element(0),side),neigh;/**elemento pai e seu vizinho*/
      if(!fat.Element()){
        PZError <<  "TPZRefPattern::NSideSubElements null father\n"; 
        exit(-1);/**radicalizou*/
      }
      neigh = fat.Neighbour();
      //isso  um crime... pq conta aqui en no conta para os outros lados ??? int count = 1;/**o pai est� contado*/
      int count = 0;
      while(neigh.Element() && neigh.Element()!=fat.Element()){
        neigh = neigh.Neighbour(); 
        count++;
      }
      if(!neigh.Element()){
        PZError <<  "TPZRefPattern::NSideSubElements null neighbour\n"; 
        exit(-1);/**terrorismo, extremismo*/
      }   
      fFatherSides.fNSubSideFather[side] = count;
      continue;
    }
    /**arestas, faces e interior: side > nnod*/
    /**tomo todos os elemento distintos da partic�o do lado*/
    int sidepos = side-nnod;
    int pos = fFatherSides.fInitSide[sidepos];
    int pos2 = fFatherSides.fInitSide[sidepos+1];
    int p,size = pos2-pos;
    //    TPZVec<int> subs(size,0);/**tamanho m�ximo*/
    //    for(p=pos;p<pos2;p++) subs[p-pos] = fFatherSides.fPartitionSubSide[p].Element()->Id();/**todos positivos > 0*/
    //    for(int i=0;i<size;i++) for(int j=i+1;j<size;j++) if(subs[j]==subs[i]) subs[j] = 0;/**anulando os repetidos*/
    //    int count = 0;
    //    for(int k=0;k<size;k++) if(subs[k]) count++;
    //    fFatherSides.fNSubSideFather[side] = count;
    fFatherSides.fNSubSideFather[side] = size;
  }
}

int TPZRefPattern::NSideSubElements(int side){
  if(side<0 || side>Element(0)->NSides()){
    PZError << "TPZRefPattern::NSideSubElement null side, side = " << side << endl;
    return -1;/**danou-se*/
  }
//  cout << "NSideSubElements " <<  fFatherSides.fNSubSideFather << endl;
  return (fFatherSides.fNSubSideFather[side]);
}


void TPZRefPattern::MeshPrint(){
  ofstream out("meshrefpatt.out");
  cout << "\nTPZRefPattern::Print imprime a malha do padrao de refinamento, arquivo de saida: meshrefpatt.out\n";
  fMesh->Print(out);
  out.flush();
  out.close();
}

TPZGeoEl *TPZRefPattern::CreateGeoEl(int ntype, int mat,TPZVec<int> &nodes,TPZGeoMesh *gmesh){
  switch(ntype) {//tipo de elemento
    case 2://unidimensional ; elg1d =
      return new TPZGeoEl1d(nodes,mat,*gmesh);
      //return;
    case 3://tri�ngulo ; elgt2d =
      return new TPZGeoElT2d(nodes,mat,*gmesh);
      //return;
    case 4://quadril�tero ; elgq2d
      return new TPZGeoElQ2d(nodes,mat,*gmesh);
      //return;
    case 7://tetraedro ; elgt3d =
      return new TPZGeoElT3d(nodes,mat,*gmesh);
      //return;
    case 5://pir�mide ; elgpi3d =
      return new TPZGeoElPi3d(nodes,mat,*gmesh);
      //return;
    case 6://prisma ; elgpi3d =
      return new TPZGeoElPr3d(nodes,mat,*gmesh);
      //return;
    case 8://cubo ; elgc3d =
      return new TPZGeoElC3d(nodes,mat,*gmesh);
      //return;
    default:
      for(int i=0;i<1;i++)
        PZError << "\nTPZRefPattern::CreateGeoEl -> Elemento nao conhecido\n";
        PZError << "\nGood Bye::Program Aborted\n";
        exit(1);
  }
  PZError <<  "\nTPZRefPattern::CreateGeoEl ntype error, ntype = " << ntype << endl;
  PZError << "\nAborted program\n";
  exit(-1);/**acabou a festa*/
}

int TPZRefPattern::SizeOfSubsSides(int ison){/**ison � o n�mero do sub-elemento*/
  int nsubs = NSubElements();
  if(ison==nsubs) ison--;
  if(ison < 0 || ison > nsubs-1){
    PZError <<  "TPZRefPattern::SizeOfSubsSides filho nao existe, filho = " << ison << endl;
  }
  int count = 0,isub;
  for(isub=0;isub<nsubs;isub++){    
    count += Element(isub+1)->NSides();
    if(isub == ison) return count;
  }
  return 0;
}

int TPZRefPattern::IsFatherNeighbour(TPZGeoElSide fathside,TPZGeoEl *son){
  /**se elementos filho e pai tem um lado em comum eles compartilham a vizinhan�a*/
  TPZGeoElSide neighbour = fathside.Neighbour();
  while(neighbour.Element() && neighbour.Element()->Id() != fathside.Element()->Id()){
    if(neighbour.Element()->Id() == son->Id()) return 1;/**elementos pai e filho s�o vizinhos*/
    neighbour = neighbour.Neighbour();/**percorre-se a vizinhan�a do lado do pai*/
  }
  return 0;//**se n�o � vizinho ent�o nenhum sub-elemento compartilha esse lado*/
}

void TPZRefPattern::DefinitionOfSizePartition(){
  TPZGeoEl *fat = Element(0);/**elemento pai da divis�o*/
  int nsidefat = fat->NSides();
  int nnodes = fat->NNodes();
  int nsides = nsidefat-nnodes;
  fFatherSides.fInitSide.Resize(nsides+1);/**n�mero de lados do pai + final*/  
  int maxsize = 0,sidefat,sf;
//  int nsubs = NSubElements();
  int size = fTransforms.fFatherSide.NElements();  
  fFatherSides.fInitSide[0] = 0;
  for(sf=nnodes;sf<nsidefat;sf++){/**os cantos do pai n�o s�o particionados*/
    int sd = 0;
    int nsidestot = 0;
    while(sd < size){
      sidefat = fTransforms.fFatherSide[sd];/**percorre fFatherSide: lado do pai contendo lado de filho*/
      if(sidefat == sf) nsidestot++;/**conta todos os lados repetidos igual a sf*/      
      sd++;
    }
    maxsize += nsidestot;/**n�mero total de lados que particionam os lados do pai (exceto cantos do pai)*/
    fFatherSides.fInitSide[sf-nnodes + 1] = maxsize;/**comeco do pr�ximo lado do pai*/    
  }
  fFatherSides.fPartitionSubSide.Resize(maxsize);
  /**PARA TESTES*/
  ofstream out("dimofpartition.out");
  out << "Valores de fFatherSides.fInitSide[sd]\n\n";
  for(int sd=0;sd<nsides;sd++){
    out << "sd : fInitSide[sd] = " << sd << " : " << fFatherSides.fInitSide[sd] << endl;
  }
  out << "Tamanho de fPartitionSubSide =  " << maxsize; 
  cout << "\nTPZRefPattern::DefinitionOfSizePartition arquivo da particao: dimofpartition.out\n";
}

TPZGeoEl *TPZRefPattern::Element(int iel){
  int nel = NSubElements()+1;/**filhos mais o pai*/
  if(iel<0 || iel>nel){
    PZError <<  "TPZRefPattern::Element elemento nao existe, elemento de id = " << iel << endl;
  }
  return ( fMesh->ElementVec()[iel]  );
}

int TPZRefPattern::SidePartition(TPZVec<TPZGeoElSide> &gelvec, int side){
  int nsides = Element(0)->NSides();/**n�mero de lados do pai*/
  int nnodes = Element(0)->NNodes();
  if(side<nnodes || side>nsides){
    PZError <<  "TPZRefaPattern::SidePartition side error: side = " << side << endl;
    gelvec.Resize(0);
    return 0;
  }
  side -= Element(0)->NNodes();/**n�o inclui cantos*/
  int firstpos = fFatherSides.fInitSide[side];/**posic�o inicial da partic�o do lado*/
  int lastpos = fFatherSides.fInitSide[side+1];/**posic�o final da partic�o do lado*/
  int size = lastpos-firstpos;
  gelvec.Resize(size);/**tamanho: n�mero de elementos da partic�o do lado*/
  int pos;
  for(pos = firstpos;pos<lastpos;pos++){
    int pos0 = pos - firstpos;
    gelvec[pos0] = fFatherSides.fPartitionSubSide[pos];
  }
  return size;/**n�mero de elementos da partic�o do lado*/;
}

/**verifica as transformac�es do lado do sub-elemento para o lado do pai*/
void TPZRefPattern::TransformationTest(){
  int isub,sd,ip;
  //x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> x1(3),pf(3),x2(3),xpf(3,0.);
  REAL weight;
  TPZGeoEl *father = Element(0);/**pai*/
  int dimfatside,fatside,nsides;
  int dimfat = father->Dimension();
  TPZGeoEl *subel;
  int nsubs = NSubElements();
  for(isub=0;isub<nsubs;isub++){    
    subel  = Element(isub+1);
    nsides = subel->NSides();
    for(sd=0;sd<nsides;sd++){
      if( sd<subel->NNodes() && IsFatherNeighbour(TPZGeoElSide(father,sd),subel) ) continue;
      int dims = subel->SideDimension(sd);
      int dimsub = subel->Dimension();
      /**regra de integrac�o para o espaco param�trico do lado do sub-elemento*/
      TPZIntPoints *rule = subel->CreateSideIntegrationRule(sd,5);
      TPZVec<int> order(dims,5);
      TPZVec<REAL> point(dims,0.),point2(dimsub);
      rule->SetOrder(order);
      for(ip=0;ip<rule->NPoints();ip++){
        /**ponto no espaco param�trico do lado do filho*/
        rule->Point(ip,point,weight);        
        TPZTransform sidet(dimsub);/**transformac�o unit�ria*/
        if(dims < dimsub){          
          sidet = subel->SideToSideTransform(sd,nsides-1);/**transf. no elemento metre*/
       //   TPZTransform told(sidet);
      //    sidet =  subel->SideToElemShapeT(sd);
      //    if(IsNotEqual(sidet,told)) 
       //     PZError <<  "TPZRefPattern::TransformationTest as transformacoes nao sao iguais -> 1\n\n";
        }
        sidet.Apply(point,point2);//**transformac�o para o interior do mestre do sub-elemento*
        subel->X(point2,x1);/**ponto no lado do filho deformado*/
        /**transformac�o: espaco param�trico do filho/lado  -> espaco param�trico do pai/lado*/
        fatside = father->WhichSide(x1);/**no elemento mestre do pai, father � o deformado*/;
        dimfatside = father->SideDimension(fatside);
        TPZTransform trans = Transform(sd,isub);/**transforma��o calculada pelo TPZRefPattern*/
        /**------------teste das transforma��es-----------*/
        TPZTransform sdtosd(dims);
        TPZGeoElSide(subel,sd).SideTransform3(TPZGeoElSide(father,fatside),sdtosd);
        if(IsNotEqual(trans,sdtosd)){
          PZError <<  "TPZRefPattern::TransformationTest as transformacoes nao sao iguais -> 2\n\n";
        }
        /**-----------fim teste das transforma��es--------*/
        TPZVec<REAL> pf(dimfatside);
        trans.Apply(point,pf);/**ponto pf no espaco param�trico do lado do pai*/
        TPZTransform sidetf(dimfatside);/**unit�ria do lado do pai*/
        if(dimfatside < dimfat){
          /**para o interior do sub-elemento*/  
          sidetf = father->SideToSideTransform(fatside,father->NSides()-1);
          //sidetf =  father->SideToElemShapeT(fatside);
        }
        sidetf.Apply(pf,xpf);/**do espaco param�trico do lado do pai para o interior do pai*/
        father->X(xpf,x2);//ponto no lado do pai deformado	
        if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
          PZError << "\nTransformacao errada\n";
          PZError << "son    = " << (subel->Id()) << endl;
          PZError << "father = " << (father->Id()) << endl;
          PZError << "side   = " << sd << endl << endl;
          int ok;
          cin >> ok;
        } else {
/*           cout << "Transformacao OK!\n"; */
/*           cout << "Filho/lado : " << subel->Id() << "/" << sd << endl; */
/*           cout << "Pai : " << father->Id() << endl << endl; */
        }//fim if sqrt..
      }//fim rule
    }//fim for sd
  }//fim for isub
}

int TPZRefPattern::IsNotEqual(TPZTransform &Told, TPZTransform &Tnew){
  int nrows = Told.Mult().Rows();
  int ncols = Told.Mult().Cols();
  if(Tnew.Mult().Rows()!=nrows || Tnew.Mult().Cols()!=ncols) return 1;
  if(Tnew.Sum().Rows()!=Told.Sum().Rows() || Tnew.Sum().Cols()!=Told.Sum().Cols()) return 1;
  int c,r;
   for(r=0;r<nrows;r++){
     for(c=0;c<ncols;c++){         
      if( fabs(Tnew.Mult()(r,c)-Told.Mult()(r,c)) > 1.e-12 ) return 1;
    }
    if( fabs(Tnew.Sum()(r,0)-Told.Sum()(r,0)) > 1.e-12 ) return 1; 
  }
  return 0;
}


void TPZRefPattern::CreateNewNodes (TPZGeoEl * gel, TPZVec<int> &newnodeindexes){
  int side;
  int nnodes = gel->NCornerNodes();
  int totalnodes = fMesh->NNodes();
  newnodeindexes.Resize(totalnodes);
  
  if (gel->HasSubElement()){
    cout << "CreateNewNodes called for an element which is already divided. \n";
    cout.flush();
    return;
// Este trecho nao teria funcionado mesmo. Nada garante que esta ordem seria
// respeitada pelo padrao de refinamento
//    for (side=gel->NCornerNodes();side<gel->NSides();side++){
//      gel->MidSideNodeIndex(side,index);
//      newnodeindexes[side-sum] = index;
//    }
//    return;
  }

  int nsides = gel->NSides();
  //Nao estou mais iterando em i... melhorou??
  for (side = nnodes;side<nsides;side++){
    CreateMidSideNodes(gel,side,newnodeindexes);
  }
}

void TPZRefPattern::CreateMidSideNodes (TPZGeoEl * gel, int side, TPZVec<int> &newnodeindexes){
  
  int i,j,k,index;
  TPZGeoMesh *gmesh = gel->Mesh();
  //SideNodes retorna um vetor com os indices dos nos internos da malha refpatern
  //com ele eu sei quantos nos internos ou, na linguagem antiga , quantos MidSideNodes tem
  TPZManVector<int> sidenodes;
  SideNodes(side,sidenodes);
  TPZGeoElSide gelside(gel,side);
  TPZGeoElSide neighbour(gelside.Neighbour());
  TPZManVector<int> sideindices(0);
  while(neighbour.Element() && neighbour != gelside) {
    //if(!neighbour.HasSubElement()) {
    if(neighbour.HasSubElement()) {      
      neighbour.Element()->MidSideNodeIndices(neighbour.Side(),sideindices);
      break;
    }
    neighbour = neighbour.Neighbour();
  }
  for (j=0;j<sidenodes.NElements();j++){
    index = sidenodes[j];
    //coordenadas do novo no na malha ref pattern
    TPZVec<REAL> refnodecoord(3,0.);
    TPZManVector<REAL,3> neighbourcoord(3,0.);
    for (k=0;k<3;k++) refnodecoord[k] = Mesh()->NodeVec()[index].Coord(k);
    //passando para as coordenadas do elemento da malha real...
    TPZManVector<REAL,3> newnodecoord(Element(0)->Dimension(),0.);
    //coordenada no espaco do elemento mestre do elemento de
    //referencia da malha refpattern
    Element(0)->ComputeXInverse(refnodecoord,newnodecoord);
    //coordenada espacial do no na malha real
    gel->X(newnodecoord,refnodecoord);
    newnodeindexes[index] = -1;
    //verificar se um vizinho ja criou o no
    for(i=0; i< sideindices.NElements(); i++) {
      for(k=0; k<3; k++) neighbourcoord[k] = gmesh->NodeVec()[sideindices[i]].Coord(k);
      REAL dif = 0.;
      for (k=0;k<3;k++) {
        dif += (refnodecoord[k] - neighbourcoord[k]) * (refnodecoord[k] - neighbourcoord[k]);
      }
      if (dif < 1e-12) {
        newnodeindexes[index] = sideindices[i];
        break;
      }
    }
    if (newnodeindexes[index] == -1) {
      //Caso o no nao exista nos vizinhos sera necessario cria-lo...
      int newindex = gmesh->NodeVec().AllocateNewElement();
      gmesh->NodeVec()[newindex].Initialize(refnodecoord,*gmesh);
      newnodeindexes[index] = newindex;
    }
  }
}
/** Returns the refinement pattern identifier */
string TPZRefPattern::GetName(){
  return fName;
}
