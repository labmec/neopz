/**
 * @file
 * @brief Contains the implementation of the TPZSubCompMesh methods.
 */

#include "pzsubcmesh.h"
#include "pzgmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "pzmathyperelastic.h"
#include "pznonlinanalysis.h"
#include "pzskylmat.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzsmfrontalanal.h"
#include "pzsmanal.h"
#include "pzbndcond.h"
#include "pzvisualmatrix.h"

#include <stdio.h>

#include <sstream>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.subcmesh"));
static LoggerPtr logger2(Logger::getLogger("pz.mesh.tpzcompmesh"));
#endif

/// Number of elements to test 
const int numel=1;
/// Angle in radians to test
static REAL angle = 0.2;

/// Defining function force (external to material) \f$ F(x,y) = (0.5-y)seno(angle) + (x-0.5)[coseno(angle) - 1] \f$ \f$ (F = disp) \f$
static void Forcing(const TPZVec<REAL> &x, TPZVec<REAL> &disp){
	disp[0] = -(x[1]-0.5)*sin(angle)+(x[0]-0.5)*cos(angle)-(x[0]-0.5);
	disp[1] = (x[1]-0.5)*cos(angle)+(x[0]-0.5)*sin(angle)-(x[1]-0.5);
	disp[2] = 0.;
}


// Construction/Destruction

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
	
	//int numeq;
	TPZVec<int> skyline;
	
	// Insert the boundary conditions
	TPZFMatrix val1(3,3,0.),val2(3,1,0.);
	TPZAutoPointer<TPZMaterial> bnd = meumat->CreateBC (meumat,-1,0,val1,val2);
	mesh.InsertMaterialObject(bnd);
	bnd = TPZAutoPointer<TPZMaterial>(meumat->CreateBC (meumat,-2,0,val1,val2));
    
	bnd->SetForcingFunction(new TPZDummyFunction(Forcing));
	mesh.InsertMaterialObject(bnd);
	
	mesh.AutoBuild();
	mesh.InitializeBlock();
	
	
	//numeq = mesh.NEquations();
	
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
	
	return 0;
}


TPZSubCompMesh::TPZSubCompMesh(TPZCompMesh &mesh, int &index) : TPZCompMesh(mesh.Reference()), TPZCompEl(mesh,0,index),
fSingularConnect(-1) {
	
	fAnalysis = NULL;
	
}

TPZSubCompMesh::TPZSubCompMesh() : TPZCompMesh(), TPZCompEl(), fSingularConnect(-1)  {
	
	fAnalysis = NULL;
}

TPZSubCompMesh::~TPZSubCompMesh(){
}


TPZCompMesh * TPZSubCompMesh::FatherMesh() const{
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

/**
 * Compute the number of elements connected to each connect object
 */
void TPZSubCompMesh::ComputeNodElCon()
{
	TPZCompMesh::ComputeNodElCon();
	std::map<int,int>::iterator it;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) {
		fConnectVec[it->second].IncrementElConnected();
	}
    if(fSingularConnect >= 0)
    {
        fConnectVec[fSingularConnect].IncrementElConnected();
    }
	/*
	 int ic;
	 for(ic = 0; ic< fConnectVec.NElements(); ic++)
	 {
	 if(fExternalLocIndex[ic] != -1)
	 {
	 fConnectVec[ic].IncrementElConnected();
	 }
	 }
	 */
}

/**
 * Compute the number of elements connected to each connect object
 */
void TPZSubCompMesh::ComputeNodElCon(TPZVec<int> &nelconnected) const
{
	TPZCompMesh::ComputeNodElCon(nelconnected);
	std::map<int,int>::const_iterator it;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) {
		nelconnected[it->second]++;
	}
	/*	int ic;
	 for(ic = 0; ic< fConnectVec.NElements(); ic++)
	 {
	 if(fExternalLocIndex[ic] != -1)
	 {
	 nelconnected[ic]++;
	 }
	 }
	 */
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

int TPZSubCompMesh::AllocateNewConnect(int nshape, int nstate, int order){
	
	int connectindex = TPZCompMesh::AllocateNewConnect(nshape, nstate, order);
	int seqnum = fConnectVec[connectindex].SequenceNumber();
    int blocksize = nshape*nstate;
	fBlock.Set(seqnum,blocksize);
	fConnectVec[connectindex].SetOrder(order);
	int i,oldsize = fExternalLocIndex.NElements();
	
	if(oldsize <= connectindex) {
		fExternalLocIndex.Resize(connectindex+1);
		for(i=oldsize; i<=connectindex;i++) fExternalLocIndex[i] = -1;
	} else {
		fExternalLocIndex[connectindex] = -1;
	}
	return connectindex;
}


void TPZSubCompMesh::MakeExternal(int local){
	if(fExternalLocIndex[local] == -1) {
		//Allocate the dependent nodes of the selected local node in father mesh
		int extconnect;
		int lastext = fConnectIndex.NElements();
		fConnectIndex.Resize(lastext+1);
		//Allocate the selected local node in father mesh
        TPZConnect &c = fConnectVec[local];
        int nshape = c.NShape();
        int nstate = c.NState();
        int order  = c.Order();
		extconnect = FatherMesh()->AllocateNewConnect(nshape,nstate,order);
		
		fConnectIndex[lastext] = extconnect;
		fExternalLocIndex[local] = lastext;
		fFatherToLocal[extconnect] = local;
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
	if(super != FatherMesh()) 
	{
		superind = FatherMesh()->GetFromSuperMesh(superind,super);
		super = FatherMesh();
	}
	std::map<int,int>::iterator it = fFatherToLocal.find(superind);
	//	int i,nc = fConnectIndex.NElements();
	//	for(i=0; i<nc; i++) if(fConnectIndex[i] == superind) break;
	//	if(i== nc) {
	if(it == fFatherToLocal.end())
	{
        TPZConnect &c = super->ConnectVec()[superind];
		int gl = AllocateNewConnect(c.NShape(),c.NState(),c.Order());
		fConnectIndex.Resize(fConnectIndex.NElements()+1);
		fConnectIndex[fConnectIndex.NElements()-1] = superind;
		fExternalLocIndex[gl] = fConnectIndex.NElements()-1;
		fFatherToLocal[superind] = gl;
		return gl;
	} else {
		int j;
		j = it->second;
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Connect in fathermesh " << superind << "  existing connect found : corresponds to connect " << j << " in subcompmesh";
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif

		return j;
	}
}

void TPZSubCompMesh::Print(std::ostream &out) const {
	
	out << "Sub Mesh" << (void *) this;
	TPZCompEl::Print(out);
	TPZCompMesh::Print(out);
	out.flush();
	int i;
	for (i=0; i<fConnectVec.NElements(); i++){
		out << "Node[" << i <<"]\t" << fExternalLocIndex[i];
		if (fExternalLocIndex[i] != -1) out << " Index in father mesh:\t" << fConnectIndex[fExternalLocIndex[i]];
		out << std::endl;
	}
	std::map<int,int>::const_iterator it;
	out << "Global to Local map ";
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) {
		out << it->first << '|' << it->second << ' ';
	}
	out << std::endl;
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
#ifdef DEBUG 
	if(father->ConnectVec()[superind].NElConnected() != 1)
	{
		std::cout << __PRETTY_FUNCTION__ << " number of elements connected to connect " << superind <<
        " = " << father->ConnectVec()[superind].NElConnected() << std::endl;
	}
	if(father && RootMesh(local) != father) {
		std::cout << "ERROR";
	}
#endif
	TPZConnect::TPZDepend *listdepend = father->ConnectVec()[superind].FirstDepend();
	while(listdepend) {
		int depfatherindex = listdepend->fDepConnectIndex;
		int depindexlocal = GetFromSuperMesh(depfatherindex,father);
		int r = listdepend->fDepMatrix.Rows();
		int c = listdepend->fDepMatrix.Cols();
		ConnectVec()[local].AddDependency(local,depindexlocal,listdepend->fDepMatrix,0,0,r,c);
		father->ConnectVec()[superind].RemoveDepend(superind,depfatherindex);
		listdepend = father->ConnectVec()[superind].FirstDepend();
	}
}

void TPZSubCompMesh::MakeInternal(int local){
	TransferDependencies(local);
	int i;
	int localindex = fExternalLocIndex[local];
	int fatherindex = fConnectIndex[localindex];
	for (i=fExternalLocIndex[local]; i<fConnectIndex.NElements()-1; i++){
		fConnectIndex[i]= fConnectIndex[i+1];
	}
	for(i=0; i<fConnectVec.NElements(); i++) {
		if(fExternalLocIndex[i] != -1 && fExternalLocIndex[i] > localindex) fExternalLocIndex[i]--;
	}
	fConnectIndex.Resize(fConnectIndex.NElements()-1);
	fFatherToLocal.erase(fatherindex);
	fExternalLocIndex[local]= -1;
}

void TPZSubCompMesh::MakeInternalFast(int local){
	TransferDependencies(local);
	int localindex = fExternalLocIndex[local];
	int fatherindex = fConnectIndex[localindex];
	fConnectIndex[localindex] = -1;
	fFatherToLocal.erase(fatherindex);
	fExternalLocIndex[local]= -1;
}

TPZCompMesh * TPZSubCompMesh::RootMesh(int local){
	if (fExternalLocIndex[local] == -1) return this;
	else return (FatherMesh()->RootMesh(fConnectIndex[fExternalLocIndex[local]]));
	//return NULL;
}

/**
 * Este m�todo deve estar errado. Primeiro tem que por os connects que tem dependencias
 * caso contrario n�s com dependencias serao duplicados
 *
 * talvez primeiro copiar a estrutura dos n�s dependentes e DEPOIS tir� los da malha pai
 */
void TPZSubCompMesh::MakeAllInternal(){
	//	TPZStack<int> stack;
	//TPZVec<int> nelcon;
	TPZCompMesh *father = FatherMesh();
	TPZAdmChunkVector<TPZConnect> &connectvec = father->ConnectVec();
#ifdef DEBUG
	//father->ComputeNodElCon();
#endif
	//father->ComputeNodElCon(nelcon);
	//#ifdef DEBUG 
	//	int in;
	//	int nn = nelcon.NElements();
	//	for (in=0; in<nn; in++) {
	//		if(father->ConnectVec()[in].NElConnected() != nelcon[in])
	//		{
	//			std::cout << "NelConnected " << in << " " << father->ConnectVec()[in].NElConnected() << " != " << nelcon[in] << std::endl;
	//		}
	//	}
	//#endif
	//TPZCompMesh::Print();
	//father->Print();
	std::set<int> cantransfer;
	std::set<int> delaytransfer;
	std::map<int,int>::iterator it;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) {
		// put the candidate nodes in the stack
		if (father->ConnectVec()[it->first].NElConnected() == 1) 
		{
			cantransfer.insert(it->second);
			//			stack.Push(it->second);
			//#ifdef DEBUG 
			//			if(father->ConnectVec()[it->first].NElConnected() != 1)
			//			{
			//				int in = it->first;
			//				std::cout << "NelConnected " << in << " " << father->ConnectVec()[in].NElConnected() << " != " << nelcon[in] << std::endl;
			//			}
			//#endif
		}
	}
	// look for dependent nodes
	while (cantransfer.size() || delaytransfer.size()) 
	{
		std::set<int>::iterator itset;
		for (itset = cantransfer.begin(); itset != cantransfer.end(); itset++) {
			TPZConnect &con = connectvec[*itset];
			TPZConnect::TPZDepend *listdepend = con.FirstDepend();
			while (listdepend) {
				if (cantransfer.find(listdepend->fDepConnectIndex) != cantransfer.end()) {
					delaytransfer.insert(listdepend->fDepConnectIndex);
				}
				listdepend = listdepend->fNext;
			}
		}
		for (itset=delaytransfer.begin(); itset != delaytransfer.end(); itset++) {
			cantransfer.erase(*itset);
		}
		
		for (itset=cantransfer.begin(); itset!=cantransfer.end(); itset++) 
		{
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Making the connect index " << *itset << " internal";
				LOGPZ_DEBUG(logger,sout.str())				
			}
#endif
			MakeInternalFast(*itset);
		}
		cantransfer = delaytransfer;
		delaytransfer.clear();
	}
	/*
	 for (i=0;i<fConnectVec.NElements();i++){
	 if (fExternalLocIndex[i]==-1) continue;
	 // put the candidate nodes in the stack
	 if (father->ConnectVec()[fConnectIndex[fExternalLocIndex[i]]].NElConnected() == 1) stack.Push(i);
	 }
	 */
	// put the independent connects first
	
	/*
	 while(stack.NElements()) {
	 int locind = stack.Pop();
	 int can = 0;
	 for(j=0; j<stack.NElements();j++)
	 {
	 int jlocind = stack[j];
	 if (jlocind == locind) continue;
	 TPZConnect &conj = father->ConnectVec()[fConnectIndex[fExternalLocIndex[stack[j]]]];
	 // special procedure when the node has dependencies
	 if (conj.FirstDepend())
	 {
	 TPZConnect::TPZDepend *listdepend = conj.FirstDepend();
	 // if the node upon which locind is dependent is already on the stack, no further analysis required
	 if (listdepend->HasDepend(fConnectIndex[fExternalLocIndex[locind]])) 
	 {
	 #ifdef LOG4CXX
	 {
	 std::stringstream sout;
	 sout << "Connect " << locind << " cannot be made internal because of " << jlocind;
	 LOGPZ_DEBUG(logger,sout.str())
	 }
	 #endif
	 break;
	 }
	 }
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
	 // if the node is not internal to the fathermesh, don't put it on the stack
	 if(can && RootMesh(locind) != FatherMesh()) can = 0;
	 if (can) 
	 {
	 #ifdef LOG4CXX
	 {
	 std::stringstream sout;
	 sout << "Making the connect index " << locind << " internal";
	 LOGPZ_DEBUG(logger,sout.str())				
	 }
	 #endif
	 MakeInternal(locind);
	 }
	 }
	 */
	
	fConnectIndex.Resize(fFatherToLocal.size());
	int count = 0;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) 
	{
		fConnectIndex[count] = it->first;
		fExternalLocIndex[it->second] = count;
		count++;
	}	
#ifdef DEBUG 
	if (count != (int)fFatherToLocal.size()) {
		DebugStop();
	}
#endif
	TPZCompMesh::ExpandSolution();
	//TPZCompMesh::Print();
	//father->Print();
	//std::cout.flush();
}

void TPZSubCompMesh::PotentialInternal(std::list<int> &connectindices) const {
	int i;
	TPZCompMesh *father = FatherMesh();
	TPZVec<int> nelconnected;
	father->ComputeNodElCon(nelconnected);
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
			if(father->ConnectVec()[extcon].NElConnected() == 1) 
			{
				connectindices.push_back(i);
			}
		}
	}
}


void TPZSubCompMesh::SetConnectIndex(int inode, int index){
	fConnectIndex[inode] = index;
}

int TPZSubCompMesh::TransferElementFrom(TPZCompMesh *mesh, int elindex){
	if(mesh == this) return elindex;
#ifdef DEBUG
	if (! IsAllowedElement(mesh,elindex)) {
		std::cout <<"TPZSubCompMesh::TransferElementFrom ERROR: trying to transfer an element not allowed" << std::endl;
		DebugStop();
		return -1;
	}
#endif
	if (mesh != FatherMesh()){
		elindex = FatherMesh()->TransferElementFrom(mesh,elindex);
		mesh = FatherMesh();
	}
#ifdef DEBUG 
	if (CommonMesh(mesh) != mesh){
		std::cout <<"TPZSubCompMesh::TransferElementFrom ERROR: mesh is not supermesh" << std::endl;
		DebugStop();
		return -1;
	}
#endif
	TPZCompMesh *father = FatherMesh();
	TPZCompEl *cel = father->ElementVec()[elindex];
	if (!cel) {
		std::cout <<"TPZSubCompMesh::TransferElementFrom ERROR: element not existing" << std::endl;
		DebugStop();
		return -1;
	}
	int i,ncon = cel->NConnects();
	for (i=0; i<ncon; i++){
		int superind = cel->ConnectIndex(i);
		int subindex = GetFromSuperMesh(superind,father);
		TPZInterfaceElement *interf = dynamic_cast<TPZInterfaceElement *> (cel);
		if(!interf)
		{
			cel->SetConnectIndex(i,subindex);
		}
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
#ifdef DEBUG 
	TPZCompMesh *common = CommonMesh(mesh);
	if ( common!= mesh){
		std::cout <<"TPZSubCompMesh::TransferElementTo ERROR: mesh is not supermesh" << std::endl;
		return -1;
	}
#endif
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
	if(fAnalysis)
	{
	}
	else
	{
		DebugStop();//this->SetAnalysis();
	}
	std::set<int> matids = fAnalysis->StructMatrix()->MaterialIds();
	if(!NeedsComputing(matids))
	{
		ek.Reset();
		ef.Reset();
		return;
	}	
	int i=0;
	CleanUpUnconnectedNodes();
	PermuteExternalConnects();
	
	
	
	TPZBlock &block = Mesh()->Block();
	//	TPZFMatrix &MeshSol = Mesh()->Solution();
	// clean ek and ef
	
	//	int nmeshnodes = fConnectVec.NElements();
	int numeq=0, numeq2=0;
	//??
	int ic;
	for (ic=0; ic<fConnectIndex.NElements(); ic++) {
		int conindex = fConnectIndex[ic];
		TPZConnect &cn = Mesh()->ConnectVec()[conindex];
		if (cn.SequenceNumber()<0 || cn.HasDependency()) {
			DebugStop();
		}
		int seqnum = cn.SequenceNumber();
		int blsize = Mesh()->Block().Size(seqnum);
		numeq2 += blsize;
	}
	int numeq3 = Mesh()->NEquations();
    {
        int ftlsize = fFatherToLocal.size();
        int ncon = fConnectIndex.NElements();
		
        if(ftlsize != ncon)
        {
            std::cout << "Number of connects of the submesh out of sink\n";
            DebugStop();
        }
    }
	std::map<int,int>::iterator it;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) 
	{
		i = it->second;
		//	for (i=0; i< nmeshnodes; i++){
		//		if(fExternalLocIndex[i] == -1) {
		TPZConnect &df = fConnectVec[i];
		if (fExternalLocIndex[i] == -1) {
			LOGPZ_ERROR(logger,"fExternalLocIndex and fFatherToLocal are out of sink")
			DebugStop();
		}
		if(df.HasDependency() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
		int seqnum = df.SequenceNumber();
		numeq += Block().Size(seqnum);
		//		}
	}
	int nconstrconnects = 0;
	int globeq2 = 0;
	for (ic=0; ic<fConnectVec.NElements(); ic++) {
		TPZConnect &cn = fConnectVec[ic];
		if (cn.HasDependency()) {
			nconstrconnects++;
		}
		else if(cn.SequenceNumber() >= 0) {
			globeq2 += Block().Size(cn.SequenceNumber());
		}
		
	}
	
	if (numeq != numeq2 || numeq > numeq3) {
		DebugStop();
	}
	int numextconnects = fConnectIndex.NElements();
	int nconnects = fConnectVec.NElements();
	int numintconnects = nconnects-numextconnects-nconstrconnects;
	
	// check whether the connects are properly enumerated
#ifdef DEBUG 
	{
		int globeq = TPZCompMesh::NEquations();
		if (globeq2 != globeq) {
			DebugStop();
		}
		int numinteq = globeq - numeq;
		int in;
		// verify whether the block structure is resequenced...
		for (in=0; in<nconnects-1; in++) {
			int blsize = Block().Size(in);
			int pos1 = Block().Position(in);
			int pos2 = Block().Position(in+1);
			if (pos2-pos1 != blsize) {
				DebugStop();
			}
		}
		int numinteq2 = Block().Position(numintconnects);
		if (numinteq != numinteq2) {
			DebugStop();
		}
		//int globeq = TPZCompMesh::NEquations();
		
		for(in=0; in<nconnects; in++)
		{
			TPZConnect &df = ConnectVec()[in];
			if( ! df.NElConnected() || df.SequenceNumber() == -1) continue;
			int seqnum = df.SequenceNumber();
			int eq = Block().Position(seqnum);
			int blsize = Block().Size(seqnum);
			if((eq < numinteq || seqnum < numintconnects) && fExternalLocIndex[in] != -1 )
			{
				std::stringstream sout;
				sout << "Connect " << in << " has equation " << eq << " but is external";
				LOGPZ_ERROR(logger,sout.str())
				DebugStop();
			}
			if ((eq >= numinteq || seqnum >= numintconnects) && blsize && !df.HasDependency() && fExternalLocIndex[in] == -1) {
				std::stringstream sout;
				sout << "Connect " << in << " has equation " << eq << " but is internal and has no dependencies ";
				df.Print(*this,sout);
				LOGPZ_ERROR(logger,sout.str())
				DebugStop();
			}
			if((eq < globeq || seqnum < nconnects-nconstrconnects) && df.HasDependency())
			{
				std::stringstream sout;
				sout << "Connect " << in << " with dependency was not put at the end of the stack equation " << eq << " global equations " << globeq;
				LOGPZ_ERROR(logger,sout.str())
				DebugStop();
			}
		}
	}
#endif
	//??
	
	TPZAutoPointer<TPZMaterial> mat = MaterialVec().begin()->second;
	int nstate = mat->NStateVariables();
	ek.fNumStateVars = nstate;
	ef.fNumStateVars = nstate;
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,1);
	
	int nelemnodes = NConnects();
	
	ek.fBlock.SetNBlocks(nelemnodes);
	ef.fBlock.SetNBlocks(nelemnodes);
	for (i = 0; i < nelemnodes ; i++)	{
		//int nodeindex = ConnectIndex(i);
		int seqnum = Connect(i).SequenceNumber();
  		ek.fBlock.Set(i,block.Size(seqnum));
  		ef.fBlock.Set(i,block.Size(seqnum));
	}
	ek.fConnect.Resize(nelemnodes);
	ef.fConnect.Resize(nelemnodes);
	
	for(i=0; i<nelemnodes; ++i){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
	if (! fAnalysis){
		TPZFStructMatrix local(this);
		TPZAutoPointer<TPZMatrix> stiff = local.CreateAssemble(ef.fMat,NULL);
		ek.fMat = *(stiff.operator->());
		//		TPZStructMatrix::Assemble(ek.fMat,ef.fMat,*this,-1,-1);
	}
	else{
		//if(!fAnalysis->Solver().Matrix())
		{
			fAnalysis->Run(std::cout);
			if(fAnalysis->AmIKilled()){
				return;
			}
		}
		
		TPZSubMeshFrontalAnalysis *sman = dynamic_cast<TPZSubMeshFrontalAnalysis *> (fAnalysis.operator->());
		if(sman)
		{
			TPZAbstractFrontMatrix *frontmat = dynamic_cast<TPZAbstractFrontMatrix *> (fAnalysis->Solver().Matrix().operator->());
			if(frontmat)
			{
				sman->SetFront(frontmat->GetFront());
			}
		}
		
		//Trying to get a derived Analysis which is a SubMeshAnalysis.
		//It could be better done with an abstract class SubMeshAnalysis which defines CondensedSolution method
		TPZSubMeshAnalysis * castedAnal = dynamic_cast<TPZSubMeshAnalysis *>(fAnalysis.operator->());
		if(castedAnal){
			castedAnal->CondensedSolution(ek.fMat,ef.fMat);
		}
		
		TPZSubMeshFrontalAnalysis * castedAnalFrontal = dynamic_cast<TPZSubMeshFrontalAnalysis *>(fAnalysis.operator->());
		if(castedAnalFrontal){
			castedAnalFrontal->CondensedSolution(ek.fMat,ef.fMat);
		}
		
		if(!castedAnal && !castedAnalFrontal){
			DebugStop();
		}
		
	}
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Substructure stiffness matrix\n";
		ek.Print(sout);
		sout << "Substructure right hand side\n";
		ef.Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//ek.fMat->Print();
}

void TPZSubCompMesh::SetAnalysisSkyline(int numThreads, int preconditioned, TPZAutoPointer<TPZGuiInterface> guiInterface){
	fAnalysis = new TPZSubMeshAnalysis(this);
	fAnalysis->SetGuiInterface(guiInterface);
	
	TPZAutoPointer<TPZStructMatrix> str = NULL;
	
	if(numThreads > 0){
		str = new TPZParSkylineStructMatrix(this,numThreads);
		str->AssembleOnlyInternalEquations();
	}
	else{
		str = new TPZSkylineStructMatrix(this);
		str->AssembleOnlyInternalEquations();
	}
	PermuteExternalConnects();
    SaddlePermute();
	
	
	str->SetNumThreads(numThreads);
    TPZAutoPointer<TPZMatrix> mat = str->Create();
    TPZAutoPointer<TPZMatrix> mat2 = mat->Clone();
	
	fAnalysis->SetStructuralMatrix(str);
	TPZStepSolver *step = new TPZStepSolver(mat);
    TPZStepSolver *gmrs = new TPZStepSolver(mat2);
    step->SetReferenceMatrix(mat2);
	step->SetDirect(ELDLt);
    gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
	TPZAutoPointer<TPZMatrixSolver> autostep = step;
    TPZAutoPointer<TPZMatrixSolver> autogmres = gmrs;
    if(preconditioned)
    {
        fAnalysis->SetSolver(autogmres);
    }
    else
    {
        fAnalysis->SetSolver(autostep);
    }
	
#ifdef DEBUG 
	{
		TPZFMatrix fillin;
		int resolution = 100;
		ComputeFillIn(resolution,fillin);		
#ifdef USING_BOOST
		std::string out("matrix_boost.vtk");
#else
		std::string out("matrix_native.vtk");
#endif
		VisualMatrix(fillin,out);
	}
#endif
	
}

void TPZSubCompMesh::SetAnalysisFrontal(int numThreads, TPZAutoPointer<TPZGuiInterface> guiInterface){
	
	fAnalysis = new TPZSubMeshFrontalAnalysis(this);
	fAnalysis->SetGuiInterface(guiInterface);
	
#ifdef DEBUG
	{
		TPZFMatrix fillin;
		int resolution = 100;
		ComputeFillIn(resolution,fillin);		
#ifdef USING_BOOST
		std::string out("matrix_boost.vtk");
#else
		std::string out("matrix_native.vtk");
#endif
		VisualMatrix(fillin,out);
	}
#endif
	
	TPZAutoPointer<TPZStructMatrix> fstr = NULL;
	if(numThreads){
		fstr = new TPZParFrontStructMatrix<TPZFrontNonSym>(this);
		static_cast<TPZParFrontStructMatrix<TPZFrontNonSym> *>(fstr.operator->())
		->SetNumberOfThreads(numThreads+2); //o frontal tem dois threads auxiliares
	}
	else{
		fstr = new TPZFrontStructMatrix<TPZFrontNonSym>(this);
	}
	
	fstr->SetNumThreads(numThreads);
	fAnalysis->SetStructuralMatrix(fstr);
	
	TPZStepSolver solver;
    solver.SetDirect(ELU);
	fAnalysis->SetSolver(solver);
	
	LOGPZ_DEBUG(logger2,__PRETTY_FUNCTION__)
	PermuteExternalConnects();
}

/**
 * Compute the permutation vector which puts the internal connects to the first on the list
 * Respect the previous order of the connects
 */
void TPZSubCompMesh::ComputePermutationInternalFirst(TPZVec<int> &permute) const
{
	// map from sequence number of the pontentially internal nodes to the node indices
	// first the independent nodes, then the dependent nodes
	std::map<int,int> independent;
	std::list<int> internal;
	this->PotentialInternal(internal);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Index = " << Index() << " Internal connects ic/seqnum";
		std::list<int>::iterator it;
		for(it=internal.begin(); it!= internal.end(); it++)
		{
			sout << *it << "/" << ConnectVec()[*it].SequenceNumber() << " ";
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	TPZCompMesh *father = this->FatherMesh();
	std::list<int>::iterator it;
	for(it=internal.begin(); it!= internal.end(); it++)
	{
		int locind = *it;
		int externallocindex = this->fExternalLocIndex[locind];
		if(externallocindex > 0)
		{
			int superind = fConnectIndex[externallocindex];
			if(father->ConnectVec()[superind].FirstDepend())
			{
			}
			else
			{
				independent[ConnectVec()[locind].SequenceNumber()] = locind;
			}
		}
		else if (!ConnectVec()[locind].FirstDepend())
		{
			independent[ConnectVec()[locind].SequenceNumber()] = locind;
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Mesh Address " << (void *) this << " Index = " << Index() << " \nIndependent connect sequence numbers and indices ";
		std::map<int,int>::iterator mapit;
		for(mapit=independent.begin(); mapit!= independent.end(); mapit++)
		{
			sout << "[" << mapit->first << " , " << mapit->second << "] ";
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	permute.Resize(0);
	permute.Resize(fConnectVec.NElements(),-1);
	
	int count = 0;
	std::map<int,int>::iterator mapit;
	for(mapit=independent.begin(); mapit!=independent.end(); mapit++)
	{
		permute[mapit->first] = count++;
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Index = " << Index() << " Permutation vector 1 " << permute;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
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
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Index = " << Index() << " Permutation vector 2 " << permute;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

/**
 * Permute the potentially internal connects to the first on the list
 * Respect the previous order of the connects
 */
void TPZSubCompMesh::PermuteInternalFirst(TPZVec<int> &permute)
{
	this->ComputePermutationInternalFirst(permute);
	LOGPZ_DEBUG(logger,"Permuting")
	Permute(permute);
}

void TPZSubCompMesh::PermuteExternalConnects(){
	//compute number of internal nodes -> numinternal
	//	TPZCompMesh::Print();
	
	int i=0, numinternal=0, numconstraints = 0, numexternal=0;
	//int countinternal=0
	int countconstraint=0;
	int nconnects = fConnectVec.NElements();
	std::set<int> internalseqnum;
	//std::cout << "fExternalLocIndex\n";
	//for(i=0; i<nconnects; i++) std::cout << fExternalLocIndex[i] << ' ';
	//std::cout << std::endl;
	for(i=0;i<nconnects; i++){
		if (fExternalLocIndex[i]==-1){
			// which are not free and which are not constrained
			TPZConnect &no = fConnectVec[i];
			
			if(no.NElConnected() == 0) continue;
			//se n�o tiver elemento conectado tambe'm
			if(no.HasDependency())
			{
				numconstraints++;
			}
			else
			{
				numinternal+= 1;
				internalseqnum.insert(no.SequenceNumber());
			}
		}
		else
		{
			numexternal++;
		}
	}
	countconstraint = numinternal+numexternal;
	// initialize a counter for internal nodes
	i=0;
	TPZManVector<int> permute(nconnects);
	for (i=0;i<nconnects;i++) permute[i] = i;
	std::set<int>::iterator it;
	int seqnum = 0;
	for(it=internalseqnum.begin(); it!=internalseqnum.end(); it++)
	{
		permute[*it] = seqnum++;
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << " numinternal " << numinternal << " numconstraints " << numconstraints << " numexternal " << numexternal << std::endl;
		//	sout << " permute so far " << permute;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	// loop over all nodes
	for (i=0;i<fConnectVec.NElements();i++){
		// take seqnum = the sequencenumber of the node
		//TPZConnect &no = fConnectIndex[fExternalLocIndex[i]];
		TPZConnect &no = fConnectVec[i];
		seqnum = no.SequenceNumber();
		// if the node is free or constrained
		//		if (no.HasDependency() || no.NElConnected() == 0) {
		//			//->set permute[sequnum] to itself
		//			continue;
		//		}
		// if the node is internal
		if (fExternalLocIndex[i] == -1){
			//-> set permute[sequnum] to counter
			// ->set permute[seqnum] = fExternalConnectIndex+numinternal
			if(no.HasDependency())
			{
				permute[seqnum] = countconstraint;
				countconstraint++;
			}
			else
			{
			}
		}
		// if the node is external
		else
		{
			permute [seqnum] = fExternalLocIndex[i]+numinternal;
			// end loop
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Index = " << " Permutations " << permute;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	Permute(permute);
}

void TPZSubCompMesh::LoadSolution() {
// #warning ME TIRE DAQUI
	
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


void TPZSubCompMesh::SkylineInternal(TPZVec<int> &skyline) {
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
			if(df.HasDependency() || df.IsCondensed() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
            int dfsize = df.NShape()*df.NState();
#ifdef DEBUG
			int seqnum = df.SequenceNumber();
			int blsize = Block().Size(seqnum);
            if (blsize != dfsize) {
                DebugStop();
            }
#endif
            numeq += dfsize;
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

/*
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
 */


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
	WriteObjects(buf, fFatherToLocal);
    buf.Write(&fSingularConnect,1);
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
	ReadObjects(buf, fFatherToLocal);
    buf.Read(&fSingularConnect,1);
}

void TPZSubCompMesh::ComputeSolution(TPZVec<REAL> &qsi,
                                     TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix &axes){
	PZError << __PRETTY_FUNCTION__ << " - ERROR! This method is not implemented\n";
}

void TPZSubCompMesh::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
									 const TPZFMatrix &axes,  TPZSolVec &sol, TPZGradSolVec &dsol){
	PZError << __PRETTY_FUNCTION__ << " - ERROR! This method is not implemented\n";
}

void TPZSubCompMesh::ComputeSolution(TPZVec<REAL> &qsi,
									 TPZVec<REAL> &normal,
									 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix &leftaxes,
									 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix &rightaxes){
	PZError << __PRETTY_FUNCTION__ << " - ERROR! This method is not implemented\n";
}


/**
 * Creates corresponding graphical element(s) if the dimension matches
 * graphical elements are used to generate output files
 * @param graphmesh graphical mesh where the element will be created
 * @param dimension target dimension of the graphical element
 */
void TPZSubCompMesh::CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension)
{
	int nel = fElementVec.NElements();
	int iel;
	for(iel=0; iel<nel; iel++)
	{
		TPZCompEl *cel = fElementVec[iel];
		if(!cel) continue;
		cel->CreateGraphicalElement(graphmesh, dimension);
	}
}

/**
 * Verifies if any element needs to be computed corresponding to the material ids
 */
bool TPZSubCompMesh::NeedsComputing(const std::set<int> &matids)
{
	std::set<int> meshmatids;
	if(! matids.size())
	{
		return true;
	}
	int numtrue=0,numfalse=0;
	// loop over the elements
	int iel, nelem = ElementVec().NElements();
	for(iel=0; iel<nelem; iel++)
	{
		TPZCompEl *cel = ElementVec()[iel];
		if(!cel) continue;
		TPZAutoPointer<TPZMaterial> mat = cel->Material();
		if(!mat)
		{
			TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (cel);
			if(submesh)
			{
				bool result = submesh->NeedsComputing(matids);
				if(result) numtrue++;
				else numfalse++;
			}
		}
		else 
		{
			int matid = mat->Id();
			meshmatids.insert(matid);
			if(matids.find(matid) != matids.end())
			{
				numtrue++;
			}
			else {
				numfalse++;
			}
			
		}
	}
	{
		std::stringstream sout;
		sout << "Material ids contained in the mesh ";
		std::set<int>::iterator it;
		for(it = meshmatids.begin(); it != meshmatids.end(); it++)
		{
			sout << *it << " ";
		}
		sout << std::endl << "Material ids which should be computed ";
		std::set<int>::const_iterator it2;
		for(it2= matids.begin(); it2 != matids.end(); it2++)
		{
			sout << *it2 << " ";
		}
		sout << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
	if(numtrue && numfalse)
	{
		std::stringstream sout;
		sout << "A substructure should have either all elements computable or not numtrue " << numtrue << " numfalse " << numfalse;
		LOGPZ_WARN(logger,sout.str())
	}
	if(numtrue)
	{
		return true;
	}
	else {
		return false;
	}
	return false;
}

bool TPZSubCompMesh::VerifyDatastructureConsistency()
{
	// all elements of fConnectIndex should be found in fFatherToSub map
	if (fConnectIndex.NElements() != (int)fFatherToLocal.size()) {
		DebugStop();
	}
	int numberexternal = fConnectIndex.NElements();
	int i;
	for (i=0; i<numberexternal; i++) {
		if (fFatherToLocal.find(fConnectIndex[i]) == fFatherToLocal.end()) {
			DebugStop();
		}
	}
	// the number of external connects in the fExternalLocIndex should be size also
	int nel = fExternalLocIndex.NElements();
	int numext = 0;
	for (i=0; i<nel; i++) {
		if (fExternalLocIndex[i] != -1) {
			numext++;
		}
	}
	if (numext != numberexternal) {
		DebugStop();
	}
	std::map<int,int>::iterator it;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) {
		if (fExternalLocIndex[it->second] == -1) {
			DebugStop();
		}
	}
	return true;
}

int TPZSubCompMesh::NumberRigidBodyModes()
{
	if (fSingularConnect == -1) {
		return 0;
	}
	int seqnum = fConnectVec[fSingularConnect].SequenceNumber();
	return fBlock.Size(seqnum);
	
}


/// Set the number of rigid body modes associated with the internal degrees of freedom
void TPZSubCompMesh::SetNumberRigidBodyModes(int nrigid)
{
	if (fSingularConnect == -1) {
        int nshape = nrigid;
        int nstate = 1;
        int order = 1;
		fSingularConnect = AllocateNewConnect(nshape,nstate,order);
		fConnectVec[fSingularConnect].IncrementElConnected();
		int extind = FatherMesh()->AllocateNewConnect(nshape,nstate,order);
		FatherMesh()->ConnectVec()[extind].IncrementElConnected();
		int next = fConnectIndex.NElements();
		fConnectIndex.Resize(next+1);
		fConnectIndex[next] = extind;
		fExternalLocIndex[fSingularConnect] = next;
        fFatherToLocal[extind] = fSingularConnect;
        ExpandSolution();
	}
	else if(fSingularConnect != -1 && nrigid >0 ) {
		int seqnum = fConnectVec[fSingularConnect].SequenceNumber();
		fBlock.Set(seqnum,nrigid);
        ExpandSolution();
		int extind = fExternalLocIndex[fSingularConnect];
		TPZCompMesh *fathermesh = FatherMesh();
		if (fathermesh && extind < 0) {
			DebugStop();
		}
		while (fathermesh && extind > 0) {
			seqnum = fathermesh->ConnectVec()[extind].SequenceNumber();
			fathermesh->Block().Set(seqnum, nrigid);
            fathermesh->ExpandSolution();
			TPZSubCompMesh *subfather = dynamic_cast<TPZSubCompMesh *> (fathermesh);
			if (subfather) {
				extind = subfather->fExternalLocIndex[extind];
			}
			fathermesh = fathermesh->FatherMesh();
		}
	}
	else {
		// not implemented yet
		DebugStop();
	}	
}
