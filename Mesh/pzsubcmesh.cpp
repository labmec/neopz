//$Id: pzsubcmesh.cpp,v 1.22 2007-08-14 12:36:05 phil Exp $

// subcmesh.cpp: implementation of the TPZSubCompMesh class.
//
//////////////////////////////////////////////////////////////////////

#include "pzsubcmesh.h"
#include "pzgmesh.h"
#include "pzcompel.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "pzmathyperelastic.h"
#include "pznonlinanalysis.h"
#include "pzskylmat.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontStructMatrix.h"
#include "pzsmfrontalanal.h"
#include "pzbndcond.h"

#include <stdio.h>

#include <sstream>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.subcmesh"));
#endif

const int numel=1;

static REAL angle = 0.2;

static void Forcing(TPZVec<REAL> &x, TPZVec<REAL> &disp){
	disp[0] = -(x[1]-0.5)*sin(angle)+(x[0]-0.5)*cos(angle)-(x[0]-0.5);
	disp[1] = (x[1]-0.5)*cos(angle)+(x[0]-0.5)*sin(angle)-(x[1]-0.5);
	disp[2] = 0.;
}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
int TPZSubCompMesh::main() {
//	int index;

	//Create the Geometric Mesh
	TPZGeoMesh geo;

	//Define the output file name
	std::ofstream output("output.dat");\

	//Set the basic coordinate nodes
	double coordstore[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};

	TPZVec<REAL> coord(3,0.);
	int i,j;

	//Set the node coordinates
	for(i=0; i<4*(numel+1); i++) {
		for (j=0; j<2; j++) {
			coord[j] = coordstore[i%4][j];
			coord[2] = i/4;
		}
		//   	int nodeindex = geo.NodeVec().AllocateNewElement();
	geo.NodeVec()[i].Initialize(i,coord,geo);
	}

	// create the elements
	TPZGeoEl *gel[numel];
	TPZVec<int> indices(8);

	// Set the connectivities
	for(i=0; i<numel; i++) {
		// initialize node indexes
		for(j=0; j<8; j++) indices[j] = 4*i+j;
    int index;
		gel[i] = geo.CreateGeoElement(ECube,indices,1, index);
	}
	//	TPZGeoElBC t3(gel[0],20,-1,geo);
	//	TPZGeoElBC t4(gel[numel-1],25,-2,geo);
	geo.BuildConnectivity();

	//Create the computacional mesh
	TPZCompMesh mesh(&geo);

	// Insert the materials
	TPZAutoPointer<TPZMaterial> meumat = new TPZMatHyperElastic(1,1.e5,0.25);
	mesh.InsertMaterialObject(meumat);

	int numeq;
	TPZVec<int> skyline;

	// Insert the boundary conditions
	TPZFMatrix val1(3,3,0.),val2(3,1,0.);
	TPZAutoPointer<TPZMaterial> bnd = meumat->CreateBC (meumat,-1,0,val1,val2);
	mesh.InsertMaterialObject(bnd);
        bnd = TPZAutoPointer<TPZMaterial>(meumat->CreateBC (meumat,-2,0,val1,val2));
	bnd->SetForcingFunction(Forcing);
	mesh.InsertMaterialObject(bnd);

	mesh.AutoBuild();
	mesh.InitializeBlock();


	numeq = mesh.NEquations();

	// Teste 1 colocar os elementos, inclusive intermedi�rios, como sub elementos
	TPZSubCompMesh *sub[numel];

	int index = -1;
	for (i=0;i<numel;i++){
		sub[i] = new TPZSubCompMesh(mesh,index);
	}

	//Teste 2 - Passar todos os sub elementos para os subelementos

	for (i=0;i<numel;i++){
		sub[i]->TransferElement(&mesh,i);
	}

	for (i=0;i<numel;i++){
		sub[i]->MakeAllInternal();
//		sub[i]->Prints(output);
	}

//	mesh.ComputeNodElCon();
//	mesh.Print(output);
	TPZNonLinearAnalysis an(&mesh,output);


//	mesh.Print(output);
	output.flush();
//	TPZFMatrix *rhs = new TPZFMatrix(skyline);
	TPZSkylineStructMatrix strskyl(&mesh);
	an.SetStructuralMatrix(strskyl);
	an.Solution().Zero();
	TPZStepSolver sol;
//	sol.ShareMatrix(an.Solver());
	sol.SetDirect(ELDLt);
	an.SetSolver(sol);
	//	an.Solver().SetDirect(ELDLt);
	an.IterativeProcess(output,0.00001,5);

	//mesh.Print(output);
	sub[0]->LoadSolution();
	sub[0]->SetName("sub[0]");
	sub[0]->Print(output);
	output.flush();


//	int a=1;
	/*int m1, m2;
	while (a != 0){
		std::cout << "mesh 1 e 2\n";
		cin >> m1 >> m2;

		for (int i=0; i<11; i++){
			if (&mesh == sub[m1]->CommonMesh(sub[m2])){
				std::cout << 10 << "\n";
				break;
			}
			if (sub[i] == sub[m1]->CommonMesh(sub[m2])){
				std::cout << i << "\n";
				break;
			}
			else {
				std::cout << sub[i] << "\t" << sub[m1]->CommonMesh(sub[m2]) << "\n";
			}
		}
		std::cout << "Digite 0 para sair \n";
		cin >> a;
	}  */
//	sub[4]->TransferElementFrom(&mesh,0);//gel[0]->Index());
//	sub[5]->TransferElementFrom(sub[4],1);

//	sub[5]->TransferElementFrom(&mesh,0);

//	sub[5]->TransferElementTo(sub[4],0);//gel[0]->Index());
//	sub[0]->TransferElementFrom(&mesh,1);//gel[1]->Index());
//	sub[1]->TransferElementFrom(&mesh,1);//gel[1]->Index());
//	sub[2]->TransferElementFrom(&mesh,1);//gel[1]->Index());
//	sub[2]->TransferElementTo(sub[8],0);
/*	sub[8]->TransferElementFrom(&mesh,0);
	sub[8]->TransferElementTo(&mesh,0);
	sub[4]->TransferElementFrom(&mesh,0);
	sub[4]->TransferElementTo(&mesh,1);
*/

	//sub[2]->TransferElementTo(sub[4],0);//gel[0]->Index());

	/*TPZVec<int> indices(3,0);
	int nodeindex = geo.NodeVec().AllocateNewElement();
	geo.NodeVec().Initialize(nodeindex,indices,sub[9]);
	sub[9]->NodeIndex(nodeindex,sub[2]);

	TPZFMatrix dep(2,2,2.);
	int is1 = sub[9]->AllocateNewConnect();
	int is2 = sub[9]->AllocateNewConnect();
	int is3 = sub[9]->AllocateNewConnect();
	(sub[9]->fConnectVec[is1]).AddDependency(is1,is2,dep,0,0,2,2);
	sub[9]->AllocateNewConnect();
	int lastcreated = sub[9]->AllocateNewConnect();
	sub[9]->MakeExternal(lastcreated);
	sub[9]->MakeExternal(is1);
	sub[9]->NodeIndex(lastcreated,sub[2]);
	sub[2]->NodeIndex(sub[2]->AllocateNewConnect(),sub[9]);	*/

/*	sub[0]->Prints(output);
	output.flush();
	sub[1]->Prints(output);
	output.flush();
	sub[2]->Prints(output);
	output.flush();
	sub[3]->Prints(output);
	output.flush();
	sub[4]->Prints(output);
	output.flush();
	sub[5]->Prints(output);
	output.flush();
	sub[6]->Prints(output);
	output.flush();
	sub[7]->Prints(output);
	output.flush();
	sub[8]->Prints(output);
	output.flush();
	sub[9]->Prints(output);
	output.flush();
	//sub[10]->Prints(out);
	//out.flush();  */


	return 0;
}



TPZSubCompMesh::TPZSubCompMesh(TPZCompMesh &mesh, int &index) : TPZCompMesh(mesh.Reference()), TPZCompEl(mesh,0,index)  {

	fAnalysis = NULL;

}

TPZSubCompMesh::TPZSubCompMesh() : TPZCompMesh(), TPZCompEl()  {

	fAnalysis = NULL;
}

TPZSubCompMesh::~TPZSubCompMesh(){
}


TPZCompMesh * TPZSubCompMesh::FatherMesh(){
	return Mesh();
}


TPZCompMesh * TPZSubCompMesh::CommonMesh(TPZCompMesh *mesh){

	TPZStack<TPZCompMesh *> s1, s2;
	int pos1=0, pos2, comind;
	TPZCompMesh *father = FatherMesh();
	s1.Push(this);
	while (father){
		s1.Push((father));
		pos1++;
		father = s1[pos1]->FatherMesh();
	}
	pos2 = 0;
	s2.Push(mesh);
	father = mesh->FatherMesh();
	while (father){
		s2.Push(father);
		pos2++;
		father = s2[pos2]->FatherMesh();
	}
	if (s1[pos1] != s2[pos2]) return 0;
	comind=0; //The first mesh is common for all submeshes
	for (; pos1>=0 && pos2>=0; pos1--, pos2--) {
		if((s1[pos1])!=(s2[pos2])) {
			comind=pos1+1;
			return (s1[comind]);
		}
	}
	return (pos1 >=0 ) ? (s1[pos1+1]) : s2[pos2+1];
}

int TPZSubCompMesh::NConnects() const{
	return fConnectIndex.NElements();
}

int TPZSubCompMesh::ConnectIndex(int i) const{
	return fConnectIndex[i];
}

int TPZSubCompMesh::Dimension() const {
	return -1;
}


//void TPZSubCompMesh::SetMaterial(TPZAutoPointer<TPZMaterial> mat){
//}

int TPZSubCompMesh::NodeIndex(int nolocal, TPZCompMesh *super)
{
  if(super == this) return nolocal;
  TPZCompMesh *root = CommonMesh(super);
  if(!root || fExternalLocIndex[nolocal] == -1) return -1;
  int result = fConnectIndex[fExternalLocIndex[nolocal]];

  if(root == FatherMesh())
  {
    return result;
  }
  else
  {
    TPZCompMesh *father =  FatherMesh();
    TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *> (father);
    if(sub)
    {
      return sub->NodeIndex(result,super);
    }
    else
    {
      return -1;
    }
  }
/*	int rootindex = PutinSuperMesh(nolocal,root);
	return neighbour->GetFromSuperMesh(rootindex,root);*/
}

int TPZSubCompMesh::AllocateNewConnect(int blocksize, int order){
/*	int connectindex = fConnectIndex.AllocateNewElement();
	int blocknum = fBlock.NBlocks();
	fBlock.SetNBlocks(blocknum+1);
	fConnectVec[connectindex].SetSequenceNumber(blocknum);
	return connectindex;
	*/

	int connectindex = TPZCompMesh::AllocateNewConnect();
	int seqnum = fConnectVec[connectindex].SequenceNumber();
	fBlock.Set(seqnum,blocksize);
        fConnectVec[connectindex].SetOrder(order);
	int i,oldsize = fExternalLocIndex.NElements();

	if(oldsize <= connectindex) {
		fExternalLocIndex.Resize(connectindex+1);
		for(i=oldsize; i<=connectindex;i++) fExternalLocIndex[i] = -1;
	} else {
		fExternalLocIndex[connectindex] = -1;
	}
/*
	int nodeindex;

	if (fConnectIndex.NElements()<connectindex){
		fConnectIndex.Resize(fConnectIndex.NElements()+1);
		fExternalLocIndex.Resize(fExternalLocIndex.NElements()+1);
		fExternalLocIndex[fExternalLocIndex.NElements()-1]=-1;
		nodeindex = NodeIndex (connectindex,
		MakeExternal(connectindex);
	}
*/
return connectindex;
}


void TPZSubCompMesh::MakeExternal(int local){
	if(fExternalLocIndex[local] == -1) {
	//Allocate the dependent nodes of the selected local node in father mesh
		int extconnect;
		int lastext = fConnectIndex.NElements();
		fConnectIndex.Resize(lastext+1);
	//Allocate the selected local node in father mesh
		int blocksize = fConnectVec[local].NDof(*(TPZCompMesh *)this);
		extconnect = FatherMesh()->AllocateNewConnect(blocksize);

//		int lastext = fConnectIndex.NElements();
//		fConnectIndex.Resize(lastext+1);

		fConnectIndex[lastext] = extconnect;
		fExternalLocIndex[local] = lastext;
		TPZConnect::TPZDepend *listdepend = fConnectVec[local].FirstDepend();
		while(listdepend) {
			int depindex = listdepend->fDepConnectIndex;
			MakeExternal(listdepend->fDepConnectIndex);
			int depextind = fConnectIndex[fExternalLocIndex[depindex]];
			int r = listdepend->fDepMatrix.Rows();
			int c = listdepend->fDepMatrix.Cols();
			FatherMesh()->ConnectVec()[extconnect].AddDependency(extconnect,depextind,listdepend->fDepMatrix,0,0,r,c);
			fConnectVec[local].RemoveDepend(local,depindex);
			listdepend = fConnectVec[local].FirstDepend();
		}
	} else {
		if(fConnectVec[local].FirstDepend() ) {
			std::cout << "TPZSubCompMesh iconsistent data structure !";
		}
	}
}

int TPZSubCompMesh::PutinSuperMesh(int local, TPZCompMesh *super){
	if(super == this) return local;
	if(fExternalLocIndex[local] == -1) MakeExternal(local);
	return FatherMesh()->PutinSuperMesh(fConnectIndex[fExternalLocIndex[local]],super);
}

int TPZSubCompMesh::GetFromSuperMesh(int superind, TPZCompMesh *super){
	if(super == this) return superind;
	if(super != FatherMesh()) superind = FatherMesh()->GetFromSuperMesh(superind,super);
	int i,nc = fConnectIndex.NElements();
	for(i=0; i<nc; i++) if(fConnectIndex[i] == superind) break;
	if(i== nc) {
		int blocksize=super->ConnectVec()[superind].NDof(*(TPZCompMesh *)super);
                int order = super->ConnectVec()[superind].Order();
		int gl = AllocateNewConnect(blocksize,order);
		fConnectIndex.Resize(fConnectIndex.NElements()+1);
		fConnectIndex[fConnectIndex.NElements()-1] = superind;
		fExternalLocIndex[gl] = fConnectIndex.NElements()-1;
		return gl;
	} else {
		int j;
		for(j=0; j<fExternalLocIndex.NElements(); j++) {
			if(fExternalLocIndex[j] == i) return j;
		}
		int blocksize=super->ConnectVec()[superind].NDof(*(TPZCompMesh *)super);
                int order = super->ConnectVec()[superind].Order();
		j = AllocateNewConnect(blocksize,order);
		fExternalLocIndex[j] = i;
		return j;
	}
}

void TPZSubCompMesh::Print(std::ostream &out){

	out << "Sub Mesh" << (void *) this;
 TPZCompEl::Print(out);
	TPZCompMesh::Print(out);
	out.flush();
	int i;
	for (i=0; i<fConnectVec.NElements(); i++){
		out << "Node[" << i <<"]\t" << fExternalLocIndex[i];
		if (fExternalLocIndex[i] != -1) out << "Index in father mesh:\t" << fConnectIndex[fExternalLocIndex[i]];
		out << std::endl;
	}
}

  /**
 * Transfer the dependency list of a connect. This will
 * make the dependency disappear for the corresponding father mesh
 * It is necessary that the number of elements connected to the connect be equal one
   */
void TPZSubCompMesh::TransferDependencies(int local)
{
  if (fExternalLocIndex[local] == -1) return;
  TPZCompMesh *father = FatherMesh();
  int superind = fConnectIndex[fExternalLocIndex[local]];
  if(father->ConnectVec()[superind].NElConnected() != 1)
  {
    std::cout << __PRETTY_FUNCTION__ << " number of elements connected to connect " << superind <<
        " = " << father->ConnectVec()[superind].NElConnected() << std::endl;
  }
  if(father && RootMesh(local) != father) {
    std::cout << "ERROR";
  }
  TPZConnect::TPZDepend *listdepend = father->ConnectVec()[superind].FirstDepend();
  while(listdepend) {
    int depfatherindex = listdepend->fDepConnectIndex;
    int depindexlocal = GetFromSuperMesh(depfatherindex,father);
    int r = listdepend->fDepMatrix.Rows();
    int c = listdepend->fDepMatrix.Cols();
    ConnectVec()[local].AddDependency(local,depindexlocal,listdepend->fDepMatrix,0,0,r,c);
    father->ConnectVec()[depfatherindex].RemoveDepend(superind,depfatherindex);
    listdepend = father->ConnectVec()[superind].FirstDepend();
  }
}


void TPZSubCompMesh::MakeInternal(int local){
  TransferDependencies(local);
  int i;
  int localindex = fExternalLocIndex[local];
  for (i=fExternalLocIndex[local]; i<fConnectIndex.NElements()-1; i++){
          fConnectIndex[i]= fConnectIndex[i+1];
  }
  for(i=0; i<fConnectVec.NElements(); i++) {
          if(fExternalLocIndex[i] != -1 && fExternalLocIndex[i] > localindex) fExternalLocIndex[i]--;
  }
  fConnectIndex.Resize(fConnectIndex.NElements()-1);
  fExternalLocIndex[local]= -1;
}


TPZCompMesh * TPZSubCompMesh::RootMesh(int local){
	if (fExternalLocIndex[local] == -1) return this;
	else return (FatherMesh()->RootMesh(fConnectIndex[fExternalLocIndex[local]]));
	return NULL;
}

/**
 * Este m�todo deve estar errado. Primeiro tem que por os connects que tem dependencias
 * caso contrario n�s com dependencias serao duplicados
 *
 * talvez primeiro copiar a estrutura dos n�s dependentes e DEPOIS tir� los da malha pai
 */
void TPZSubCompMesh::MakeAllInternal(){
  TPZStack<int> stack;
  int i,j;
  TPZCompMesh *father = FatherMesh();
  father->ComputeNodElCon();
  //TPZCompMesh::Print();
  //father->Print();
  for (i=0;i<fConnectVec.NElements();i++){
    if (fExternalLocIndex[i]==-1) continue;
    // put the candidate nodes in the stack
    if (father->ConnectVec()[fConnectIndex[fExternalLocIndex[i]]].NElConnected() == 1) stack.Push(i);
  }
  // put the independent connects first
  while(stack.NElements()) {
    int locind = stack.Pop();
    TPZConnect &coni = father->ConnectVec()[fConnectIndex[fExternalLocIndex[locind]]];
    int can = 0;
    // special procedure when the node has dependencies
    if (coni.FirstDepend()){
      TPZConnect::TPZDepend *listdepend = coni.FirstDepend();
      for(j=0;j<stack.NElements(); j++){
        int jlocind = stack[j];
        if (jlocind == locind) continue;
        // if the node upon which locind is dependent is already on the stack, no further analysis required
        if (listdepend->HasDepend(fConnectIndex[fExternalLocIndex[jlocind]])) break;
      }
      // no element on the stack is listed as dependent from the current node
      if (j == stack.NElements())
      {
        can=1;
      }
      // we found an element in the dependency list. Let s check it first
      else

      {
        // put the node upon which the current node depends in the current position and the dependent node at the end
        int jlocind = stack[j];
        stack[j] = locind;
        stack.Push(jlocind);
      }
    }
    // the node has no dependencies
    else {
            can=1;
    }
    // if the node is not internal to the fathermesh, don't put it on the stack
    if(can && RootMesh(locind) != FatherMesh()) can = 0;
    if (can) {
            MakeInternal(locind);
    }
  }
  //TPZCompMesh::Print();
  //father->Print();
  //std::cout.flush();
}

void TPZSubCompMesh::PotentialInternal(std::list<int> &connectindices){
  int i;
  TPZCompMesh *father = FatherMesh();
  father->ComputeNodElCon();
  //TPZCompMesh::Print();
  //father->Print();
  for (i=0;i<fConnectVec.NElements();i++){
    if (fExternalLocIndex[i]==-1)
    {
      connectindices.push_back(i);
    }
    else
    {
      int extcon = this->fConnectIndex[fExternalLocIndex[i]];
      if(father->ConnectVec()[extcon].NElConnected() == 1) connectindices.push_back(i);
    }
  }
}


void TPZSubCompMesh::SetConnectIndex(int inode, int index){
	fConnectIndex[inode] = index;
}

int TPZSubCompMesh::TransferElementFrom(TPZCompMesh *mesh, int elindex){
	if(mesh == this) return elindex;
		if (! IsAllowedElement(mesh,elindex)) {
		std::cout <<"TPZSubCompMesh::TransferElementFrom ERROR: trying to transfer an element not allowed" << std::endl;
		return -1;
	}
	if (mesh != FatherMesh()){
		elindex = FatherMesh()->TransferElementFrom(mesh,elindex);
	}
	if (CommonMesh(mesh) != mesh){
		std::cout <<"TPZSubCompMesh::TransferElementFrom ERROR: mesh is not supermesh" << std::endl;
		return -1;
	}
	TPZCompMesh *father = FatherMesh();
	TPZCompEl *cel = father->ElementVec()[elindex];
	if (!cel) {
		std::cout <<"TPZSubCompMesh::TransferElementFrom ERROR: element not existing" << std::endl;
		return -1;
	}
	int i,ncon = cel->NConnects();
	for (i=0; i<ncon; i++){
		int superind = cel->ConnectIndex(i);
		int subindex = GetFromSuperMesh(superind,father);
		cel->SetConnectIndex(i,subindex);
	}
	cel->SetMesh(this);
        if(cel->Reference())
        {
          TPZAutoPointer<TPZMaterial> mat = cel->Material();
          if(!mat)
          {
            father->CopyMaterials(*this);
          }
        }

//	int blocksize=mesh->ConnectVec()[elindex].NDof((TPZCompMesh *)mesh);
	int newelind = fElementVec.AllocateNewElement();
	fElementVec[newelind] = cel;
        cel->SetIndex(newelind);
	father->ElementVec()[elindex] = 0;
	father->ElementVec().SetFree(elindex);
	return newelind;
}

int TPZSubCompMesh::TransferElementTo(TPZCompMesh *mesh, int elindex){
	TPZCompMesh *common = CommonMesh(mesh);
	if ( common!= mesh){
		std::cout <<"TPZSubCompMesh::TransferElementTo ERROR: mesh is not supermesh" << std::endl;
		return -1;
	}
	if(mesh == this) return elindex;

	if (mesh != FatherMesh()){
		FatherMesh()->TransferElementTo(mesh,elindex);
	}


	TPZCompMesh *father = FatherMesh();
	if(elindex >= ElementVec().NElements()){
		std::cout <<"TPZSubCompMesh::TransferElementTo ERROR: not possible transfer non existing element" << std::endl;
		return -1;
	}
	TPZCompEl *cel = ElementVec()[elindex];
	if (!cel) {
		std::cout <<"TPZSubCompMesh::TransferElementTo ERROR: not possible transfer null element" << std::endl;
		return -1;
	}
	int i,ncon = cel->NConnects();
	for (i=0; i<ncon; i++){
		int subindex = cel->ConnectIndex(i);
		MakeExternal(subindex);
		int superind = fConnectIndex[fExternalLocIndex[subindex]];
		cel->SetConnectIndex(i,superind);
	}
//	int blocksize=father->ConnectVec()[elind].NDof(father);
	int newelind = father->ElementVec().AllocateNewElement();
	father->ElementVec()[newelind] = cel;
	cel->SetMesh(father);
        cel->SetIndex(newelind);
	ElementVec()[elindex] = 0;
	ElementVec().SetFree(elindex);
	return newelind;
}

int TPZSubCompMesh::TransferElement(TPZCompMesh *mesh, int elindex){
	TPZCompMesh *comm = CommonMesh(mesh);
	int newelind = mesh->TransferElementTo(comm,elindex);
	int ell=TransferElementFrom(comm,newelind);
	InitializeBlock();
	return ell;
}

int TPZSubCompMesh::IsAllowedElement(TPZCompMesh *mesh, int elindex){
	if (CommonMesh(mesh) == mesh){
		TPZCompMesh *father = this;
		while(father->FatherMesh() != mesh) {
			father = father->FatherMesh();
		}
		TPZSubCompMesh *sub = (TPZSubCompMesh *) father;
		int index = sub->Index();
		return (elindex != index);
	}
	return 1;
}

void TPZSubCompMesh::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){
	if(!fAnalysis) this->SetAnalysis();
	int i=0;
	CleanUpUnconnectedNodes();
	PermuteExternalConnects();

	TPZBlock &block = Mesh()->Block();
	//	TPZFMatrix &MeshSol = Mesh()->Solution();
	// clean ek and ef

	int nmeshnodes = fConnectVec.NElements();
	int numeq=0;
	//??

	for (i=0; i< nmeshnodes; i++){
		if(fExternalLocIndex[i] == -1) {
				TPZConnect &df = fConnectVec[i];
				if(df.HasDependency() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
				int seqnum = df.SequenceNumber();
				numeq += Block().Size(seqnum);
		}
	}
	numeq = (TPZCompMesh::NEquations()) - numeq;
//??

	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,1);

	int nelemnodes = NConnects();

	ek.fBlock.SetNBlocks(nelemnodes);
	ef.fBlock.SetNBlocks(nelemnodes);
	for (i = 0; i < nelemnodes ; i++)	{
		int nodeindex = ConnectIndex(i);
  		ek.fBlock.Set(i,block.Size(nodeindex));
  		ef.fBlock.Set(i,block.Size(nodeindex));
	  }
	  ek.fConnect.Resize(nelemnodes);
	  ef.fConnect.Resize(nelemnodes);

	  for(i=0; i<nelemnodes; ++i){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	  }
	if (! fAnalysis){
		TPZStructMatrix::Assemble(ek.fMat,ef.fMat,*this);
	}
	else{
		fAnalysis->Run(std::cout);
		fAnalysis->CondensedSolution(ek.fMat,ef.fMat);
//		ek.fMat->Print("ek reduzido");
//		std::cout.flush();
	}
	//ek.fMat->Print();
}

void TPZSubCompMesh::SetAnalysis(){
	if(fAnalysis) delete fAnalysis;
	fAnalysis = new TPZSubMeshFrontalAnalysis(this);
	//	int numint = NumInternalEquations();
	TPZFrontStructMatrix<TPZFrontSym> fstr(this);
	fAnalysis->SetStructuralMatrix(fstr);
	TPZStepSolver solver;
	fAnalysis->SetSolver(solver);
	PermuteExternalConnects();
//	ofstream out("subcmesh.dat");
//	Prints(out);
}

  /**
 * Permute the potentially internal connects to the first on the list
 * Respect the previous order of the connects
   */
void TPZSubCompMesh::PermuteInternalFirst()
{
  // map from sequence number of the pontentially internal nodes to the node indices
  // first the independent nodes, then the dependent nodes
  std::map<int,int> independent;
  std::list<int> internal;
  this->PotentialInternal(internal);
/*#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << "Internal connects ic/seqnum";
    std::list<int>::iterator it;
    for(it=internal.begin(); it!= internal.end(); it++)
    {
      sout << *it << "/" << ConnectVec()[*it].SequenceNumber() << " ";
    }
    LOGPZ_DEBUG(logger,sout.str())
  }
#endif*/
  TPZCompMesh *father = this->FatherMesh();
  std::list<int>::iterator it;
  for(it=internal.begin(); it!= internal.end(); it++)
  {
    int locind = *it;
    int superind = fConnectIndex[this->fExternalLocIndex[locind]];
    if(father->ConnectVec()[superind].FirstDepend())
    {
    }
    else
    {
      independent[ConnectVec()[locind].SequenceNumber()] = locind;
    }
  }
  TPZManVector<int> permute(fConnectVec.NElements(),-1);

  int count = 0;
  std::map<int,int>::iterator mapit;
  for(mapit=independent.begin(); mapit!=independent.end(); mapit++)
  {
    permute[mapit->first] = count++;
  }
  std::map<int,int> seqmap;
  int ind;
  for(ind=0; ind < fConnectVec.NElements(); ind++)
  {
    int seqnum = fConnectVec[ind].SequenceNumber();
    if(seqnum == -1) continue;
    seqmap[seqnum]=ind;
  }
  for(mapit=seqmap.begin(); mapit!=seqmap.end(); mapit++)
  {
    if(permute[mapit->first] == -1) permute[mapit->first] = count++;
  }
/*#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << "Permutation vector " << permute;
    LOGPZ_DEBUG(logger,sout.str())
  }
#endif*/
  Permute(permute);
}

void TPZSubCompMesh::PermuteExternalConnects(){
	//compute number of internal nodes -> numinternal
//	TPZCompMesh::Print();

	int i=0, numinternal=0;
	int nconnects = fConnectVec.NElements();
//std::cout << "fExternalLocIndex\n";
//for(i=0; i<nconnects; i++) std::cout << fExternalLocIndex[i] << ' ';
//std::cout << std::endl;
	for(i=0;i<nconnects; i++){
		if (fExternalLocIndex[i]==-1){
			// which are not free and which are not constrained
			TPZConnect &no = fConnectVec[i];

			if(no.NElConnected() == 0) continue;
			//se n�o tiver elemento conectado tambe'm
			numinternal+= 1;
		}
	}
	// initialize a counter for internal nodes
	i=0;
	int seqnum=0;
	int countint=0;
	TPZManVector<int> permute(nconnects);
	for (i=0;i<nconnects;i++) permute[i] = i;

	// loop over all nodes
	for (i=0;i<fConnectVec.NElements();i++){
		// take seqnum = the sequencenumber of the node
		//TPZConnect &no = fConnectIndex[fExternalLocIndex[i]];
		TPZConnect &no = fConnectVec[i];
		seqnum = no.SequenceNumber();
		// if the node is free or constrained
		if (no.HasDependency() || no.NElConnected() == 0) {
			//->set permute[sequnum] to itself
			continue;
		}
		// if the node is internal
		if (fExternalLocIndex[i] == -1){
			//-> set permute[sequnum] to counter
			permute[seqnum] = countint;
			//-> increment counter
			countint += 1;
		}
		// if the node is external
		else{
			// ->set permute[seqnum] = fExternalConnectIndex+numinternal
			permute [seqnum] = fExternalLocIndex[i]+numinternal;
		// end loop
		}
	}
	//for (i=0;i<NConnects();i++){
	//	std::cout << "Permute [" <<i <<"]\n";
	//}
	Permute(permute);
}

void TPZSubCompMesh::LoadSolution(){
  //	int count = 0;
	int i=0;
	int seqnumext;
	int seqnumint;
	//	int numinteq = NumInternalEquations();
	int size;
	TPZFMatrix &sol = Mesh()->Solution();

	for (i=0;i<fConnectVec.NElements(); i++) {
		if (fExternalLocIndex[i] != -1) {
			TPZConnect &noext = Mesh()->ConnectVec()[fConnectIndex[fExternalLocIndex[i]]];
			TPZConnect &noint = fConnectVec[i];
			seqnumext = noext.SequenceNumber();
			size = (Mesh()->Block()).Size(seqnumext);
			seqnumint = noint.SequenceNumber();
			int posext = Mesh()->Block().Position(seqnumext);
			int posint = fBlock.Position(seqnumint);
			int l;
			for(l=0;l<size;l++) {
				fSolution(posint+l,0) = sol(posext+l,0);
			}
		}
	}
	if(fAnalysis) fAnalysis->LoadSolution(fSolution);
	TPZCompMesh::LoadSolution(fSolution);
}


void TPZSubCompMesh::Skyline(TPZVec<int> &skyline) {
	TPZCompMesh::Skyline(skyline);
	skyline.Resize(NumInternalEquations());
}

int TPZSubCompMesh::NumInternalEquations() {
	int nmeshnodes = fConnectVec.NElements();
	int numeq=0;
	//??

	int i;
	for (i=0; i< nmeshnodes; i++){
		if(fExternalLocIndex[i] == -1) {
				TPZConnect &df = fConnectVec[i];
				if(df.HasDependency() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
				int seqnum = df.SequenceNumber();
				numeq += Block().Size(seqnum);
		}
	}
	return numeq;

}


REAL TPZSubCompMesh::CompareElement(int var, char *matname){
	return CompareMesh(var,matname);
}


void TPZSubCompMesh::LoadElementReference()
{
	TPZCompMesh::LoadReferences();
}


void TPZSubCompMesh::GetExternalConnectIndex (TPZVec<int> &extconn){
	int i;
	int count = 0;
	for(i=0; i<fExternalLocIndex.NElements(); i++){
		if (fExternalLocIndex[i] != -1) count++;
	}
	extconn.Resize(count);

	count=0;
	for(i=0; i<fExternalLocIndex.NElements(); i++){
		if (fExternalLocIndex[i] != -1){
			extconn[count] = i;
			count++;
		}
	}
	return;
}


TPZAnalysis * TPZSubCompMesh::GetAnalysis()
{
	return fAnalysis;
}

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
int TPZSubCompMesh::ClassId() const
{
  return TPZSUBCOMPMESHID;
}
template class
    TPZRestoreClass< TPZSubCompMesh, TPZSUBCOMPMESHID>;

  /**
  Save the element data to a stream
  */
void TPZSubCompMesh::Write(TPZStream &buf, int withclassid)
{
  TPZCompEl::Write(buf,withclassid);
  TPZCompMesh::Write(buf,0);
  WriteObjects(buf,fConnectIndex);
  WriteObjects(buf,fExternalLocIndex);
}

  /**
  Read the element data from a stream
  */
void TPZSubCompMesh::Read(TPZStream &buf, void *context)
{
  TPZCompEl::Read(buf,context);
  TPZCompMesh::Read(buf,Mesh()->Reference());
  ReadObjects(buf,fConnectIndex);
  ReadObjects(buf,fExternalLocIndex);
}

void TPZSubCompMesh::ComputeSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes){
  PZError << __PRETTY_FUNCTION__ << " - ERROR! This method is not implemented\n";
}

void TPZSubCompMesh::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                                const TPZFMatrix &axes,  TPZVec<REAL> &sol, TPZFMatrix &dsol){
  PZError << __PRETTY_FUNCTION__ << " - ERROR! This method is not implemented\n";
}

void TPZSubCompMesh::ComputeSolution(TPZVec<REAL> &qsi,
                                TPZVec<REAL> &normal,
                                TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                                TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes){
  PZError << __PRETTY_FUNCTION__ << " - ERROR! This method is not implemented\n";
}


