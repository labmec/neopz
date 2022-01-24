/**
 * @file
 * @brief Contains the implementation of the TPZSubCompMesh methods.
 */

#include "pzsubcmesh.h"
#include "pzgmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzcmesh.h"
#include "TPZElementMatrixT.h"
#include "pznonlinanalysis.h"
#include "pzskylmat.h"
#include "TPZMatrixSolver.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzsmfrontalanal.h"
#include "pzsmanal.h"
#include "TPZMaterial.h"
#include "TPZMatLoadCases.h"
#include "TPZBndCond.h"
#include "pzvisualmatrix.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


#include <stdio.h>

#include <sstream>
#include <iterator>
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.subcmesh");
static TPZLogger logger2("pz.mesh.tpzcompmesh");
#else
static int logger;
#endif

TPZSubCompMesh::TPZSubCompMesh(TPZCompMesh &mesh) : TPZRegisterClassId(&TPZSubCompMesh::ClassId), TPZCompMesh(mesh.Reference()), TPZCompEl(mesh,0),
fSingularConnect(-1) {
    SetDimModel(mesh.Dimension());
	fAnalysis = NULL;
	
}

TPZSubCompMesh::TPZSubCompMesh() : TPZRegisterClassId(&TPZSubCompMesh::ClassId),TPZCompMesh(), TPZCompEl(), fSingularConnect(-1)  {
	
	fAnalysis = NULL;
}

TPZSubCompMesh::~TPZSubCompMesh(){

	// THIS ROUTINE NEEDS TO INCLUDE THE DELETION OF THE LIST POINTERS
	TPZGeoMesh *ref = TPZCompMesh::Reference();
	if (ref){
		ref->ResetReference();
		this->LoadReferences();
	}
#ifdef PZDEBUG
    ComputeNodElCon();
#endif
	int64_t i, nelem = this->NElements();
    
	//deleting subcompmesh
	for(i=0; i<nelem; i++){
		TPZCompEl *el = fElementVec[i];
		TPZSubCompMesh * subc = dynamic_cast<TPZSubCompMesh*>(el);
		if(subc){
			delete subc;
		}
	}
	
    
	//unwrapping condensed compel
	for(i=0; i<nelem; i++){
		TPZCompEl *el = fElementVec[i];
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(el);
		if(cond){
			cond->Unwrap();
		}
	}
	
	//unwrapping element groups
	for(i=0; i<nelem; i++){
		TPZCompEl *el = fElementVec[i];
		TPZElementGroup * group = dynamic_cast<TPZElementGroup*>(el);
		if(group){
            group->Unwrap();
		}
	}
	
	//deleting interface elements
	for(i=0; i<nelem; i++){
		TPZCompEl *el = fElementVec[i];
		TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(el);
		if(face){
            delete el;
		}
	}
	
	//deleting other elements
	for(i=0; i<nelem; i++) {
		TPZCompEl *el = fElementVec[i];
		if(!el) continue;
		delete el;
	}
	
	fElementVec.Resize(0);
	fElementVec.CompactDataStructure(fElementVec.NOW);
	fConnectVec.Resize(0);
	fConnectVec.CompactDataStructure(fConnectVec.NOW);
	
	MaterialVec().clear();
}


TPZCompMesh * TPZSubCompMesh::FatherMesh() const{
	return Mesh();
}


TPZCompMesh * TPZSubCompMesh::CommonMesh(TPZCompMesh *mesh){
	
	TPZStack<TPZCompMesh *> s1, s2;
	int64_t pos1=0, pos2, comind;
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
	std::map<int64_t,int64_t>::iterator it;
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
	std::map<int64_t,int64_t>::const_iterator it;
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

int64_t TPZSubCompMesh::ConnectIndex(int i) const{
	return fConnectIndex[i];
}

int TPZSubCompMesh::Dimension() const {
	return -1;
}


//void TPZSubCompMesh::SetMaterial(TPZMaterial * mat){
//}

int64_t TPZSubCompMesh::NodeIndex(int64_t nolocal, TPZCompMesh *super)
{
	if(super == this) return nolocal;
	TPZCompMesh *root = CommonMesh(super);
	if(!root || fExternalLocIndex[nolocal] == -1) return -1;
	int64_t result = fConnectIndex[fExternalLocIndex[nolocal]];
	
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

int64_t TPZSubCompMesh::AllocateNewConnect(int nshape, int nstate, int order){
	
	int64_t connectindex = TPZCompMesh::AllocateNewConnect(nshape, nstate, order);
	int64_t seqnum = fConnectVec[connectindex].SequenceNumber();
    int blocksize = nshape*nstate;
	fBlock.Set(seqnum,blocksize);
	fConnectVec[connectindex].SetOrder(order,connectindex);
	int64_t i,oldsize = fExternalLocIndex.NElements();
	
	if(oldsize <= connectindex) {
		fExternalLocIndex.Resize(connectindex+1);
		for(i=oldsize; i<=connectindex;i++) fExternalLocIndex[i] = -1;
	} else {
		fExternalLocIndex[connectindex] = -1;
	}
	return connectindex;
}

int64_t TPZSubCompMesh::AllocateNewConnect(const TPZConnect &connect){
	
	int64_t connectindex = TPZCompMesh::AllocateNewConnect(connect);
	int64_t seqnum = fConnectVec[connectindex].SequenceNumber();
    int nshape = connect.NShape();
    int nstate = connect.NState();
    int blocksize = nshape*nstate;
	fBlock.Set(seqnum,blocksize);
	int64_t i,oldsize = fExternalLocIndex.NElements();
	
	if(oldsize <= connectindex) {
		fExternalLocIndex.Resize(connectindex+1);
		for(i=oldsize; i<=connectindex;i++) fExternalLocIndex[i] = -1;
	} else {
		fExternalLocIndex[connectindex] = -1;
	}
	return connectindex;
}


void TPZSubCompMesh::MakeExternal(int64_t local){
	if(fExternalLocIndex[local] == -1) {
		//Allocate the dependent nodes of the selected local node in father mesh
		int64_t extconnect;
		int64_t lastext = fConnectIndex.NElements();
		fConnectIndex.Resize(lastext+1);
		//Allocate the selected local node in father mesh
        TPZConnect &c = fConnectVec[local];
		extconnect = FatherMesh()->AllocateNewConnect(c);
		
		fConnectIndex[lastext] = extconnect;
		fExternalLocIndex[local] = lastext;
		fFatherToLocal[extconnect] = local;
		TPZConnect::TPZDepend *listdepend = fConnectVec[local].FirstDepend();
		while(listdepend) {
			int64_t depindex = listdepend->fDepConnectIndex;
			MakeExternal(listdepend->fDepConnectIndex);
			int64_t depextind = fConnectIndex[fExternalLocIndex[depindex]];
			int64_t r = listdepend->fDepMatrix.Rows();
			int64_t c = listdepend->fDepMatrix.Cols();
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

int64_t TPZSubCompMesh::PutinSuperMesh(int64_t local, TPZCompMesh *super){
	if(super == this) return local;
	if(fExternalLocIndex[local] == -1) MakeExternal(local);
	return FatherMesh()->PutinSuperMesh(fConnectIndex[fExternalLocIndex[local]],super);
}

int64_t TPZSubCompMesh::GetFromSuperMesh(int64_t superind, TPZCompMesh *super){
	if(super == this) return superind;
	if(super != FatherMesh()) 
	{
		superind = FatherMesh()->GetFromSuperMesh(superind,super);
		super = FatherMesh();
	}
	std::map<int64_t,int64_t>::iterator it = fFatherToLocal.find(superind);
	//	int i,nc = fConnectIndex.NElements();
	//	for(i=0; i<nc; i++) if(fConnectIndex[i] == superind) break;
	//	if(i== nc) {
	if(it == fFatherToLocal.end())
	{
        TPZConnect &c = super->ConnectVec()[superind];
		int64_t gl = AllocateNewConnect(c);
		fConnectIndex.Resize(fConnectIndex.NElements()+1);
		fConnectIndex[fConnectIndex.NElements()-1] = superind;
		fExternalLocIndex[gl] = fConnectIndex.NElements()-1;
		fFatherToLocal[superind] = gl;
		return gl;
	} else {
		int64_t j;
		j = it->second;
		
#ifdef PZ_LOG2
        if(logger.isDebugEnabled())
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
    
    {
        int64_t nc = NConnects();
        for(int ic=0; ic<nc; ic++)
        {
            TPZCompEl::Connect(ic).Print(*Mesh(),out);
        }
    }
    
	TPZCompMesh::Print(out);
	out.flush();
	int64_t i;
	for (i=0; i<fConnectVec.NElements(); i++){
		out << "Node[" << i <<"]\t" << fExternalLocIndex[i];
		if (fExternalLocIndex[i] != -1) out << " Index in father mesh:\t" << fConnectIndex[fExternalLocIndex[i]];
		out << std::endl;
	}
	std::map<int64_t,int64_t>::const_iterator it;
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
void TPZSubCompMesh::TransferDependencies(int64_t local)
{
	if (fExternalLocIndex[local] == -1) return;
	TPZCompMesh *father = FatherMesh();
	int64_t superind = fConnectIndex[fExternalLocIndex[local]];
#ifdef PZDEBUG 
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
		int64_t depfatherindex = listdepend->fDepConnectIndex;
		int64_t depindexlocal = GetFromSuperMesh(depfatherindex,father);
		int64_t r = listdepend->fDepMatrix.Rows();
		int64_t c = listdepend->fDepMatrix.Cols();
		ConnectVec()[local].AddDependency(local,depindexlocal,listdepend->fDepMatrix,0,0,r,c);
		//father->ConnectVec()[superind].RemoveDepend(superind,depfatherindex);
        listdepend = listdepend->fNext;
	}
}

void TPZSubCompMesh::MakeInternal(int64_t local){
	TransferDependencies(local);
	int64_t i;
	int64_t localindex = fExternalLocIndex[local];
	int64_t fatherindex = fConnectIndex[localindex];
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

void TPZSubCompMesh::MakeInternalFast(int64_t local){
	TransferDependencies(local);
	int64_t localindex = fExternalLocIndex[local];
	int64_t fatherindex = fConnectIndex[localindex];
//    Mesh()->ConnectVec()[fatherindex].RemoveDepend();
	fConnectIndex[localindex] = -1;
	fFatherToLocal.erase(fatherindex);
	fExternalLocIndex[local]= -1;
}

static void GatherDependency(TPZCompMesh &cmesh, TPZConnect &start, std::set<int64_t> &dependency)
{
    if (!start.HasDependency()) {
        return;
    }
    TPZConnect::TPZDepend *depend = start.FirstDepend();
    while (depend) {
        dependency.insert(depend->fDepConnectIndex);
        TPZConnect &c = cmesh.ConnectVec()[depend->fDepConnectIndex];
        GatherDependency(cmesh, c, dependency);
        depend = depend->fNext;
    }
}

TPZCompMesh * TPZSubCompMesh::RootMesh(int64_t local){
	if (fExternalLocIndex[local] == -1) return this;
	else return (FatherMesh()->RootMesh(fConnectIndex[fExternalLocIndex[local]]));
	//return NULL;
}

/**
 * Este metodo deve estar errado. Primeiro tem que por os connects que tem dependencias
 * caso contrario nos com dependencias serao duplicados
 *
 * talvez primeiro copiar a estrutura dos nos dependentes e DEPOIS tira los da malha pai
 */
void TPZSubCompMesh::MakeAllInternal(){
	//	TPZStack<int> stack;
	//TPZVec<int> nelcon;
	TPZCompMesh *father = FatherMesh();
	TPZAdmChunkVector<TPZConnect> &connectvec = father->ConnectVec();
#ifdef PZDEBUG
	//father->ComputeNodElCon();
#endif
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Connect indexes " << fConnectIndex;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif

	//father->ComputeNodElCon(nelcon);
	//#ifdef PZDEBUG
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
    
    // cantransfer contains indices in the local mesh
	std::set<int64_t> cantransfer;
	std::set<int64_t> delaytransfer;
	std::map<int64_t,int64_t>::iterator it;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) {
		// put the candidate nodes in the stack
        TPZConnect &fatherconnect = father->ConnectVec()[it->first];
#ifdef PZ_LOG
        if (logger.isDebugEnabled() && fatherconnect.HasDependency())
        {
            std::stringstream sout;
            sout << "Father connect indexes " << it->first;
            sout << "Submesh connect index " << it->second;
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif

		if (fatherconnect.NElConnected() == 1 || fatherconnect.HasDependency())
		{
			cantransfer.insert(it->second);
		}
	}
	// look for dependent nodes
	while (cantransfer.size() || delaytransfer.size()) 
	{
		std::set<int64_t>::iterator itset;
		for (itset = cantransfer.begin(); itset != cantransfer.end(); itset++) {
            // this doesnt make sense : connectvec is a vector of the father mesh
            int64_t csubmeshindex = *itset;
            int64_t elindex = this->fExternalLocIndex[csubmeshindex];
            int64_t fatherindex = fConnectIndex[elindex];
#ifdef PZDEBUG
            if (fFatherToLocal[fatherindex] != csubmeshindex) {
                DebugStop();
            }
#endif
			TPZConnect &con = connectvec[fatherindex];
            
            std::set<int64_t> dependset;
            GatherDependency(*father, con, dependset);
            for (std::set<int64_t>::iterator it = dependset.begin(); it != dependset.end(); it++) {
                int64_t connectindex = *it;
                if (fFatherToLocal.find(connectindex) != fFatherToLocal.end()) {
                    int64_t submeshconnectindex = fFatherToLocal[connectindex];
                    if (cantransfer.find(submeshconnectindex) != cantransfer.end())
                    {
                        delaytransfer.insert(submeshconnectindex);
                        cantransfer.erase(submeshconnectindex);
                    }
                }
                else
                {
                    std::cout << "The dependency of a connect should be transferred with " <<
                    " the element \n";
//                    DebugStop();
                }
            }
		}
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << " connect indexes that will be transferred ";
            std::copy(cantransfer.begin(), cantransfer.end(), std::ostream_iterator<int64_t>(sout, " "));
            sout << std::endl;
            sout << " connect indexes that are delayed ";
            std::copy(delaytransfer.begin(), delaytransfer.end(), std::ostream_iterator<int64_t>(sout, " "));
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
		for (itset=cantransfer.begin(); itset!=cantransfer.end(); itset++)
		{
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
				std::stringstream sout;
                int64_t localindex = fExternalLocIndex[*itset];
                int64_t fatherindex = fConnectIndex[localindex];
                father->ConnectVec()[fatherindex].Print(*father,sout);
                sout << "Making the connect index " << *itset << " internal " << " index in the father mesh " << fatherindex << std::endl;
                sout << "Connect indexes " << fConnectIndex;
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
	 #ifdef PZ_LOG
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
	 #ifdef PZ_LOG
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
	int64_t count = 0;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) 
	{
		fConnectIndex[count] = it->first;
		fExternalLocIndex[it->second] = count;
		count++;
	}	
#ifdef PZDEBUG 
	if (count != (int64_t)fFatherToLocal.size()) {
		DebugStop();
	}
#endif
	TPZCompMesh::ExpandSolution();
	//TPZCompMesh::Print();
	//father->Print();
	//std::cout.flush();
}

void TPZSubCompMesh::PotentialInternal(std::list<int64_t> &connectindices) const {
	int64_t i;
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
			int64_t extcon = this->fConnectIndex[fExternalLocIndex[i]];
			if(father->ConnectVec()[extcon].NElConnected() == 1) 
			{
				connectindices.push_back(i);
			}
		}
	}
}


void TPZSubCompMesh::SetConnectIndex(int inode, int64_t index){
	fConnectIndex[inode] = index;
}

int64_t TPZSubCompMesh::TransferElementFrom(TPZCompMesh *mesh, int64_t elindex){
	if(mesh == this) return elindex;
#ifdef PZDEBUG
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
#ifdef PZDEBUG 
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
    TPZInterfaceElement *interf = dynamic_cast<TPZInterfaceElement *> (cel);
    TPZMultiphysicsInterfaceElement *multinterf = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
    TPZCompEl *left = 0;
    TPZCompEl *right = 0;
    if (interf) {
        left = interf->LeftElement();
        right = interf->RightElement();
    }
    if (multinterf) {
        left = multinterf->LeftElement();
        right = multinterf->RightElement();
    }
    
//    if(!interf && !multinterf)
    {
        TPZStack<int64_t> allc;
        std::map<int64_t,int64_t> fathertolocal;
        cel->BuildConnectList(allc);
        int64_t ncon = allc.size();
        for (int64_t i=0; i<ncon; i++){
            int64_t superind = allc[i];
            int64_t subindex = GetFromSuperMesh(superind,father);
            fathertolocal[superind] = subindex;
        }
        int64_t ncon2 = cel->NConnects();
        for (int64_t ic = 0; ic<ncon2; ic++) {
            int64_t cindex = cel->ConnectIndex(ic);
            if(fathertolocal.find(cindex) == fathertolocal.end()) DebugStop();
            int64_t subindex = fathertolocal[cindex];
            cel->SetConnectIndex(ic,subindex);
        }
    }

//    else
//    {
//        int nleftcon = left->NConnects();
//        {
//            TPZCompMesh *comm = CommonMesh(left->Mesh());
//            int ncon = nleftcon;
//            for (int ic=0; ic<ncon ; ic++) {
//                int64_t superind = left->ConnectIndex(ic);
//                int64_t commind = left->Mesh()->PutinSuperMesh(superind, comm);
//                int64_t subindex = GetFromSuperMesh(commind, comm);
//                if (multinterf) {
//                    cel->SetConnectIndex(ic, subindex);
//                }
//            }
//        }
//        {
//            TPZCompMesh *comm = CommonMesh(right->Mesh());
//            int ncon = right->NConnects();
//            for (int ic=0; ic<ncon ; ic++) {
//                int64_t superind = right->ConnectIndex(ic);
//                int64_t commind = right->Mesh()->PutinSuperMesh(superind, comm);
//                int64_t subindex = GetFromSuperMesh(commind, comm);
//                if (multinterf) {
//                    cel->SetConnectIndex(ic+nleftcon, subindex);
//                }
//            }
//        }
//    }
    
    if(cel->Reference())
    {
        TPZMaterial * matfather;
        matfather = cel->Material();
        if(!matfather)
        {
            // I don't know what to do...
            DebugStop();
        }
        int matid = matfather->Id();
        TPZMaterial * matthis = FindMaterial(matid);
        
        // perform a "shallow copy" of the material
        if (!matthis) {
            MaterialVec()[matfather->Id()] = matfather;
        }
        
        // for boundary conditions we need to copy the referred material too
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(matfather);
        if (bnd) {
            TPZMaterial *ref = bnd->Material();
            int refid = ref->Id();
            TPZMaterial *matthis = FindMaterial(refid);
            if(!matthis)
            {
                MaterialVec()[refid] = ref;
            }
        }
    }
	cel->SetMesh(this);
    /*
     if(cel->Reference())
     {
     TPZMaterial * mat = cel->Material();
     if(!mat)
     {
     father->CopyMaterials(*this);
     }
     }
     */
	//	int blocksize=mesh->ConnectVec()[elindex].NDof((TPZCompMesh *)mesh);
	int64_t newelind = fElementVec.AllocateNewElement();
	fElementVec[newelind] = cel;
	cel->SetIndex(newelind);
	father->ElementVec()[elindex] = 0;
	father->ElementVec().SetFree(elindex);
	return newelind;
}

int64_t TPZSubCompMesh::TransferElementTo(TPZCompMesh *mesh, int64_t elindex){
#ifdef PZDEBUG 
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
		int64_t subindex = cel->ConnectIndex(i);
		MakeExternal(subindex);
		int64_t superind = fConnectIndex[fExternalLocIndex[subindex]];
		cel->SetConnectIndex(i,superind);
	}
	//	int blocksize=father->ConnectVec()[elind].NDof(father);
	int64_t newelind = father->ElementVec().AllocateNewElement();
	father->ElementVec()[newelind] = cel;
	cel->SetMesh(father);
	cel->SetIndex(newelind);
	ElementVec()[elindex] = 0;
	ElementVec().SetFree(elindex);
	return newelind;
}

int64_t TPZSubCompMesh::TransferElement(TPZCompMesh *mesh, int64_t elindex){
	TPZCompMesh *comm = CommonMesh(mesh);
	int64_t newelind = mesh->TransferElementTo(comm,elindex);
	int64_t ell=TransferElementFrom(comm,newelind);
	return ell;
}

int TPZSubCompMesh::IsAllowedElement(TPZCompMesh *mesh, int64_t elindex){
	if (CommonMesh(mesh) == mesh){
		TPZCompMesh *father = this;
		while(father->FatherMesh() != mesh) {
			father = father->FatherMesh();
		}
		TPZSubCompMesh *sub = (TPZSubCompMesh *) father;
		int64_t index = sub->Index();
		return (elindex != index);
	}
	return 1;
}

/// Assemble the stiffness matrix in locally kept datastructure
void TPZSubCompMesh::Assemble()
{
    if(fAnalysis)
    {
    }
    else
    {
        std::cout << "The SubCompMesh needs a configured analysis\n";
        DebugStop();//this->SetAnalysis();
    }
    std::set<int> matids = fAnalysis->StructMatrix()->MaterialIds();
    if(!NeedsComputing(matids))
    {
        return;
    }
    int i=0;
    CleanUpUnconnectedNodes();
    PermuteExternalConnects();
    fAnalysis->Assemble();

    //Trying to get a derived Analysis which is a SubMeshAnalysis.
    //It could be better done with an abstract class SubMeshAnalysis which defines CondensedSolution method
    TPZSubMeshAnalysis * castedAnal = dynamic_cast<TPZSubMeshAnalysis *>(fAnalysis.operator->());
    if(!castedAnal)
    {
        DebugStop();
    }

    TPZAutoPointer<TPZMatrix<STATE> > ReducableStiff = castedAnal->Matrix();
    if (!ReducableStiff) {
        DebugStop();
    }
    TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (ReducableStiff.operator->());
    if(!matred) DebugStop();
    
    matred->SetF(fAnalysis->Rhs());

}

/// Initialize the datastructure of ef
void TPZSubCompMesh::InitializeEF(TPZElementMatrix &ef)
{
    TPZBlock &block = Mesh()->Block();
    TPZMaterial * mat = MaterialVec().begin()->second;
    int nstate = mat->NStateVariables();
    const int numloadcases = [mat](){
        if (auto *tmp = dynamic_cast<TPZMatLoadCasesBase*>(mat); tmp){
            return tmp->NumLoadCases();
        }else{
            return 1;
        }
    }();
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
    
    int nelemnodes = NConnects();
    ef.Block().SetNBlocks(nelemnodes);
    for (int i = 0; i < nelemnodes ; i++)    {
        //int nodeindex = ConnectIndex(i);
        int64_t seqnum = Connect(i).SequenceNumber();
        ef.Block().Set(i,block.Size(seqnum));
    }
    ef.fConnect.Resize(nelemnodes);
    
    for(int i=0; i<nelemnodes; ++i){
        (ef.fConnect)[i] = ConnectIndex(i);
    }
}


template<class TVar>
void TPZSubCompMesh::CalcStiffInternal(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef){
	if(!fAnalysis)
	{
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe SubCompMesh needs a configured analysis\n";
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
	//	TPZFMatrix<REAL> &MeshSol = Mesh()->Solution();
	// clean ek and ef
	
	//	int nmeshnodes = fConnectVec.NElements();
	int64_t numeq=0, numeq2=0;
	//??
	int64_t ic;
	for (ic=0; ic<fConnectIndex.NElements(); ic++) {
		int64_t conindex = fConnectIndex[ic];
		TPZConnect &cn = Mesh()->ConnectVec()[conindex];
		if (cn.SequenceNumber()<0 || cn.HasDependency()) {
			DebugStop();
		}
		int64_t seqnum = cn.SequenceNumber();
		int blsize = Mesh()->Block().Size(seqnum);
		numeq2 += blsize;
	}
	int64_t numeq3 = Mesh()->NEquations();
    {
        int ftlsize = fFatherToLocal.size();
        int ncon = fConnectIndex.NElements();
		
        if(ftlsize != ncon)
        {
            std::cout << "Number of connects of the submesh out of sink\n";
            DebugStop();
        }
    }
	std::map<int64_t,int64_t>::iterator it;
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
		int64_t seqnum = df.SequenceNumber();
		numeq += Block().Size(seqnum);
		//		}
	}
	int64_t nconstrconnects = 0;
	int globeq2 = 0;
	for (ic=0; ic<fConnectVec.NElements(); ic++) {
		TPZConnect &cn = fConnectVec[ic];
		if (cn.HasDependency() || cn.IsCondensed()) {
			nconstrconnects++;
		}
		else if(cn.SequenceNumber() >= 0) {
			globeq2 += Block().Size(cn.SequenceNumber());
		}
		
	}
	
	if (numeq != numeq2 || numeq > numeq3) {
		DebugStop();
	}
	
	// check whether the connects are properly enumerated
#ifdef PZDEBUG 
	int64_t numextconnects = fConnectIndex.NElements();
	int64_t nconnects = fConnectVec.NElements();
	int64_t numintconnects = nconnects-numextconnects-nconstrconnects;
	{
		int64_t globeq = TPZCompMesh::NEquations();
		if (globeq2 != globeq) {
			DebugStop();
		}
		int64_t numinteq = globeq - numeq;
		int64_t in;
		// verify whether the block structure is resequenced...
		for (in=0; in<nconnects-1; in++) {
			int blsize = Block().Size(in);
			int64_t pos1 = Block().Position(in);
			int64_t pos2 = Block().Position(in+1);
			if (pos2-pos1 != blsize) {
				DebugStop();
			}
		}
        
        int64_t numinteq2 = 0;
        if(numintconnects != 0) numinteq2 = Block().Position(numintconnects-1)+Block().Size(numintconnects-1);
		if (numinteq != numinteq2) {
			DebugStop();
		}
		//int globeq = TPZCompMesh::NEquations();
		
		for(in=0; in<nconnects; in++)
		{
			TPZConnect &df = ConnectVec()[in];
			if( ! df.NElConnected() || df.SequenceNumber() == -1) continue;
			int64_t seqnum = df.SequenceNumber();
			int64_t eq = Block().Position(seqnum);
			int blsize = Block().Size(seqnum);
			if((eq < numinteq || seqnum < numintconnects) && fExternalLocIndex[in] != -1 )
			{
				std::stringstream sout;
				sout << "Connect " << in << " has equation " << eq << " but is external";
				LOGPZ_ERROR(logger,sout.str())
				DebugStop();
			}
			if ((eq >= numinteq || seqnum >= numintconnects) && blsize && !df.HasDependency() && !df.IsCondensed() && fExternalLocIndex[in] == -1) {
				std::stringstream sout;
				sout << "Connect " << in << " has equation " << eq << " but is internal and has no dependencies ";
				df.Print(*this,sout);
				LOGPZ_ERROR(logger,sout.str())
				DebugStop();
			}
			if((eq < globeq || seqnum < nconnects-nconstrconnects) && (df.HasDependency() || df.IsCondensed()) )
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
	
	TPZMaterial * mat = MaterialVec().begin()->second;
	int nstate = mat->NStateVariables();
    const int numloadcases = [mat](){
        if (auto *tmp = dynamic_cast<TPZMatLoadCasesBase*>(mat); tmp){
            return tmp->NumLoadCases();
        }else{
            return 1;
        }
    }();
    ek.fMesh = Mesh();
    ef.fMesh = Mesh();

	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
    ek.fType = TPZElementMatrix::EK;
    ef.fType = TPZElementMatrix::EF;
	
	int nelemnodes = NConnects();
	
	ek.fBlock.SetNBlocks(nelemnodes);
	ef.fBlock.SetNBlocks(nelemnodes);
	for (i = 0; i < nelemnodes ; i++)	{
		//int nodeindex = ConnectIndex(i);
		int64_t seqnum = Connect(i).SequenceNumber();
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
		TPZFStructMatrix<TVar> local(this);
		TPZAutoPointer<TPZMatrix<TVar> > stiff =
            dynamic_cast<TPZMatrix<TVar>*>(local.CreateAssemble(ef.fMat,nullptr));
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
			TPZAbstractFrontMatrix<TVar> *frontmat = dynamic_cast<TPZAbstractFrontMatrix<TVar> *> (fAnalysis->MatrixSolver<TVar>().Matrix().operator->());
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
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
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

/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
template<class TVar>
void TPZSubCompMesh::CalcResidualInternal(TPZElementMatrixT<TVar> &ef)
{
    TPZFMatrix<TVar> rhs;
    fAnalysis->AssembleResidual();
    TPZSubMeshAnalysis * castedAnal = dynamic_cast<TPZSubMeshAnalysis *>(fAnalysis.operator->());

    if (!castedAnal) {
        PZError<<__PRETTY_FUNCTION__;
        DebugStop();
    }
    InitializeEF(ef);
    castedAnal->ReducedRightHandSide(ef.fMat);
//    TPZCompMesh::CalcResidual(ef);
//    ef.PermuteGather(fIndexes);
//    fCondensed.SetF(ef.fMat);
//    //const TPZFMatrix<REAL> &f1 = fCondensed.F1Red();
//    TPZFNMatrix<100,TVar> f1(fCondensed.Dim1(),ef.fMat.Cols());
//    fCondensed.F1Red(f1);
//    int64_t dim1 = f1.Rows();
//    int64_t dim = ef.fMat.Rows();
//    int64_t dim0 = dim-dim1;
//    for (int64_t i= dim0; i<dim; i++) {
//        ef.fMat(i,0) = f1.GetVal(i-dim0,0);
//    }
}

/** @brief Sets the analysis type. */
void TPZSubCompMesh::SetAnalysisSparse(int numThreads)
{
    fAnalysis = new TPZSubMeshAnalysis(this);
    TPZAutoPointer<TPZStructMatrix> str = NULL;
    str = new TPZSSpStructMatrix<STATE>(this);
    if(numThreads > 0){
        str->SetNumThreads(numThreads);
    }
    
    SaddlePermute();
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    PermuteExternalConnects();
    str->SetNumThreads(numThreads);
    int64_t numinternal = NumInternalEquations();
    str->EquationFilter().SetMinMaxEq(0, numinternal);
    TPZAutoPointer<TPZMatrix<STATE> > mat =
        dynamic_cast<TPZMatrix<STATE>*>(str->Create());
    if(!mat){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"ERROR\n: Incompatible types. Aborting...\n";
        DebugStop();
    }
    str->EquationFilter().Reset();
    fAnalysis->SetStructuralMatrix(str);
    TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
    step->SetDirect(ELDLt);
    TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
    fAnalysis->SetSolver(autostep);

}

/** @brief Sets the analysis type. */
void TPZSubCompMesh::SetAnalysisNonSymSparse(int numThreads)
{
    fAnalysis = new TPZSubMeshAnalysis(this);
    TPZAutoPointer<TPZStructMatrix> str = NULL;

    if(numThreads > 0){
        str = new TPZSpStructMatrix<STATE>(this);
        str->SetNumThreads(numThreads);
    }
    else{
        str = new TPZSpStructMatrix<STATE>(this);
    }

    SaddlePermute();
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    PermuteExternalConnects();
    str->SetNumThreads(numThreads);
    int64_t numinternal = NumInternalEquations();
    str->EquationFilter().SetMinMaxEq(0, numinternal);
    TPZAutoPointer<TPZMatrix<STATE> > mat =
        dynamic_cast<TPZMatrix<STATE>*>(str->Create());
    if(!mat){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"ERROR\n: Incompatible types. Aborting...\n";
        DebugStop();
    }
    str->EquationFilter().Reset();
    fAnalysis->SetStructuralMatrix(str);
    TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
    step->SetDirect(ELDLt);
    TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
    fAnalysis->SetSolver(autostep);

}


/** @brief Sets the analysis type. */
void TPZSubCompMesh::SetAnalysisFStruct(int numThreads)
{
    fAnalysis = new TPZSubMeshAnalysis(this);
    TPZAutoPointer<TPZStructMatrix> str = NULL;

    str = new TPZFStructMatrix<STATE>(this);
    if(numThreads > 0){
        str->SetNumThreads(numThreads);
    }

    SaddlePermute();
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    PermuteExternalConnects();
    str->SetNumThreads(numThreads);
    int64_t numinternal = NumInternalEquations();
    str->EquationFilter().SetMinMaxEq(0, numinternal);
    TPZAutoPointer<TPZMatrix<STATE> > mat =
        dynamic_cast<TPZMatrix<STATE>*>(str->Create());
    str->EquationFilter().Reset();
    fAnalysis->SetStructuralMatrix(str);
    TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
    step->SetDirect(ELU);
    TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
    fAnalysis->SetSolver(autostep);

}


void TPZSubCompMesh::SetAnalysisSkyline(int numThreads, int preconditioned, TPZAutoPointer<TPZGuiInterface> guiInterface){
	fAnalysis = new TPZSubMeshAnalysis(this);
	fAnalysis->SetGuiInterface(guiInterface);
	TPZAutoPointer<TPZStructMatrix> str = NULL;
	str = new TPZSkylineStructMatrix<STATE>(this);
	if(numThreads > 0){
        str->SetNumThreads(numThreads);
	}
    
    SaddlePermute();
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	PermuteExternalConnects();
    
	
	
	str->SetNumThreads(numThreads);
    int64_t numinternal = NumInternalEquations();
    str->EquationFilter().SetMinMaxEq(0, numinternal);
    TPZMatrix<STATE> *skylmat = dynamic_cast<TPZMatrix<STATE> *>(str->Create());
    TPZAutoPointer<TPZMatrix<STATE> > mat = skylmat;
    str->EquationFilter().Reset();
    TPZAutoPointer<TPZMatrix<STATE> > mat2 = mat->Clone();
	
	fAnalysis->SetStructuralMatrix(str);
	TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
    TPZStepSolver<STATE> *gmrs = new TPZStepSolver<STATE>(mat2);
    step->SetReferenceMatrix(mat2);
	step->SetDirect(ELDLt);
    gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
	TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
    TPZAutoPointer<TPZMatrixSolver<STATE> > autogmres = gmrs;
    if(preconditioned)
    {
        fAnalysis->SetSolver(autogmres);
    }
    else
    {
        fAnalysis->SetSolver(autostep);
    }
	
    
#ifdef PZDEBUG2
	{
		TPZFMatrix<REAL> fillin;
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

void TPZSubCompMesh::SetAnalysisSkyline(int numThreads, int preconditioned, TPZAutoPointer<TPZRenumbering> renumber){
    fAnalysis = new TPZSubMeshAnalysis;
    fAnalysis->SetRenumber(renumber);
    fAnalysis->SetCompMesh(this, true);
    TPZAutoPointer<TPZStructMatrix> str = NULL;
    
    str = new TPZSkylineStructMatrix<STATE>(this);
    if(numThreads > 0){
        str->SetNumThreads(numThreads);
    }
    
    SaddlePermute();
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    PermuteExternalConnects();
    
    
    
    str->SetNumThreads(numThreads);
    int64_t numinternal = NumInternalEquations();
    str->EquationFilter().SetMinMaxEq(0, numinternal);
    TPZAutoPointer<TPZMatrix<STATE> > mat =
        dynamic_cast<TPZMatrix<STATE>*>(str->Create());
    str->EquationFilter().Reset();
    TPZAutoPointer<TPZMatrix<STATE> > mat2 = mat->Clone();
    
    fAnalysis->SetStructuralMatrix(str);
    TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
    TPZStepSolver<STATE> *gmrs = new TPZStepSolver<STATE>(mat2);
    step->SetReferenceMatrix(mat2);
    step->SetDirect(ELDLt);
    gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
    TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
    TPZAutoPointer<TPZMatrixSolver<STATE> > autogmres = gmrs;
    if(preconditioned)
    {
        fAnalysis->SetSolver(autogmres);
    }
    else
    {
        fAnalysis->SetSolver(autostep);
    }
    
    
#ifdef PZDEBUG
    {
        TPZFMatrix<REAL> fillin;
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
	
#ifdef PZDEBUG2
	{
		TPZFMatrix<REAL> fillin;
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
		fstr = new TPZParFrontStructMatrix<TPZFrontNonSym<STATE> >(this);
		static_cast<TPZParFrontStructMatrix<TPZFrontNonSym<STATE> > *>(fstr.operator->())
		->SetNumThreads(numThreads+2); //o frontal tem dois threads auxiliares
	}
	else{
		fstr = new TPZFrontStructMatrix<TPZFrontNonSym<STATE> >(this);
	}
	
	fstr->SetNumThreads(numThreads);
	fAnalysis->SetStructuralMatrix(fstr);
	
	TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
	fAnalysis->SetSolver(solver);
	
#ifdef PZ_LOG
	LOGPZ_DEBUG(logger2, __PRETTY_FUNCTION__)
#endif
	PermuteExternalConnects();
}

/**
 * Compute the permutation vector which puts the internal connects to the first on the list
 * Respect the previous order of the connects
 */
void TPZSubCompMesh::ComputePermutationInternalFirst(TPZVec<int64_t> &permute) const
{
	// map from sequence number of the pontentially internal nodes to the node indices
	// first the independent nodes, then the dependent nodes
	std::map<int64_t,int64_t> independent;
	std::list<int64_t> internal;
	this->PotentialInternal(internal);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Index = " << Index() << " Internal connects ic/seqnum";
		std::list<int64_t>::iterator it;
		for(it=internal.begin(); it!= internal.end(); it++)
		{
			sout << *it << "/" << ConnectVec()[*it].SequenceNumber() << " ";
		}
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
	}
#endif
	TPZCompMesh *father = this->FatherMesh();
	std::list<int64_t>::iterator it;
	for(it=internal.begin(); it!= internal.end(); it++)
	{
		int64_t locind = *it;
		int64_t externallocindex = this->fExternalLocIndex[locind];
		if(externallocindex > 0)
		{
			int64_t superind = fConnectIndex[externallocindex];
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
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Mesh Address " << (void *) this << " Index = " << Index() << " \nIndependent connect sequence numbers and indices ";
		std::map<int64_t,int64_t>::iterator mapit;
		for(mapit=independent.begin(); mapit!= independent.end(); mapit++)
		{
			sout << "[" << mapit->first << " , " << mapit->second << "] ";
		}
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
	}
#endif
	permute.Resize(0);
	permute.Resize(fConnectVec.NElements(),-1);
	
	int64_t count = 0;
	std::map<int64_t,int64_t>::iterator mapit;
	for(mapit=independent.begin(); mapit!=independent.end(); mapit++)
	{
		permute[mapit->first] = count++;
	}
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Index = " << Index() << " Permutation vector 1 " << permute;
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
	}
#endif
	std::map<int64_t,int64_t> seqmap;
	int64_t ind;
	for(ind=0; ind < fConnectVec.NElements(); ind++)
	{
		int64_t seqnum = fConnectVec[ind].SequenceNumber();
		if(seqnum == -1) continue;
		seqmap[seqnum]=ind;
	}
	for(mapit=seqmap.begin(); mapit!=seqmap.end(); mapit++)
	{
		if(permute[mapit->first] == -1) permute[mapit->first] = count++;
	}
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Index = " << Index() << " Permutation vector 2 " << permute;
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
	}
#endif
}

/**
 * Permute the potentially internal connects to the first on the list
 * Respect the previous order of the connects
 */
void TPZSubCompMesh::PermuteInternalFirst(TPZVec<int64_t> &permute)
{
	this->ComputePermutationInternalFirst(permute);
#ifdef PZ_LOG
	if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger, "Permuting")
#endif
	Permute(permute);
}

void TPZSubCompMesh::PermuteExternalConnects(){
	//compute number of internal nodes -> numinternal
	//	TPZCompMesh::Print();
    
    ComputeNodElCon();
	
	int64_t i=0, numinternal=0, numconstraints = 0, numexternal=0;
	//int countinternal=0
	int64_t countconstraint=0;
	int64_t nconnects = fConnectVec.NElements();
	std::set<int64_t> internalseqnum;
	//std::cout << "fExternalLocIndex\n";
	//for(i=0; i<nconnects; i++) std::cout << fExternalLocIndex[i] << ' ';
	//std::cout << std::endl;
	for(i=0;i<nconnects; i++){
		if (fExternalLocIndex[i]==-1){
			// which are not free and which are not constrained
			TPZConnect &no = fConnectVec[i];
			
			if(no.NElConnected() == 0) continue;
			//se n�o tiver elemento conectado tambe'm
			if(no.HasDependency() || no.IsCondensed())
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
	TPZManVector<int64_t> permute(nconnects);
	for (i=0;i<nconnects;i++) permute[i] = i;
	std::set<int64_t>::iterator it;
	int64_t seqnum = 0;
	for(it=internalseqnum.begin(); it!=internalseqnum.end(); it++)
	{
		permute[*it] = seqnum++;
	}
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
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
			if(no.HasDependency() || no.IsCondensed())
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
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Index = " << " Permutations " << permute << std::endl;
        std::set<int64_t> permval;
        permval.insert(&permute[0], (&permute[permute.size()-1]+1));
        sout << " Number of distinct values in permute " << permval.size();
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	Permute(permute);
}

template<class TVar>
void TPZSubCompMesh::LoadSolutionInternal(TPZFMatrix<TVar> &mysol) {
	
	int64_t i=0;
	int64_t seqnumext;
	int64_t seqnumint;
	//	int numinteq = NumInternalEquations();
	int size;
	TPZFMatrix<TVar> &extsol = *Mesh()->Block().Matrix<TVar>();
	
	for (i=0;i<fConnectVec.NElements(); i++) {
		if (fExternalLocIndex[i] != -1) {
			TPZConnect &noext = Mesh()->ConnectVec()[fConnectIndex[fExternalLocIndex[i]]];
			TPZConnect &noint = fConnectVec[i];
			seqnumext = noext.SequenceNumber();
			size = (Mesh()->Block()).Size(seqnumext);
			seqnumint = noint.SequenceNumber();
			int64_t posext = Mesh()->Block().Position(seqnumext);
			int64_t posint = fBlock.Position(seqnumint);
			int l;
			for(l=0;l<size;l++) {
				mysol(posint+l,0) = extsol(posext+l,0);
			}
		}
	}
//    fSolution.Print(std::cout);
    
	if(fAnalysis) fAnalysis->LoadSolution(mysol);
	else TPZCompMesh::LoadSolution(fSolution);
}

void TPZSubCompMesh::LoadSolution(){
    LoadSolutionInternal<STATE>(fSolution);
}

/**
 * @brief Compute the integral of a variable defined by the string if the material id is included in matids
 */
TPZVec<STATE> TPZSubCompMesh::IntegrateSolution(const std::string &varname, const std::set<int> &matids)
{
    return TPZCompMesh::Integrate(varname,matids);
}


void TPZSubCompMesh::TransferMultiphysicsElementSolution()
{
    int64_t nel = this->NElements();
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        fSolution.Print("SubMeshSol",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = this->Element(iel);
        if (!cel) {
            continue;
        }
        cel->TransferMultiphysicsElementSolution();
    }
}


void TPZSubCompMesh::SkylineInternal(TPZVec<int64_t> &skyline) {
	TPZCompMesh::Skyline(skyline);
	skyline.Resize(NumInternalEquations());
}

int64_t TPZSubCompMesh::NumInternalEquations() {
	int64_t nmeshnodes = fConnectVec.NElements();
	int64_t numeq=0;
	//??
	
	int64_t i;
	for (i=0; i< nmeshnodes; i++){
		if(fExternalLocIndex[i] == -1) {
			TPZConnect &df = fConnectVec[i];
			if(df.HasDependency() || df.IsCondensed() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
            int dfsize = df.NShape()*df.NState();
#ifdef PZDEBUG
			int64_t seqnum = df.SequenceNumber();
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
int TPZSubCompMesh::ClassId() const{
    return Hash("TPZSubCompMesh") ^ TPZCompMesh::ClassId() << 1 ^ TPZCompEl::ClassId() << 2;
}

#ifndef BORLAND
template class TPZRestoreClass< TPZSubCompMesh>;
#endif

/**
 Save the element data to a stream
 */
void TPZSubCompMesh::Write(TPZStream &buf, int withclassid) const
{
    //std::map<int, TPZMaterial *> matmap = MaterialVec();
    //MaterialVec().clear();
	TPZCompEl::Write(buf,withclassid);
	TPZCompMesh::Write(buf,0);
    //MaterialVec() = matmap;//AQUIFRAN
    const std::map<int, TPZMaterial *> &matmap = fMaterialVec;
    TPZManVector<int> matindex(matmap.size(),-1);
    int count=0;
    for (std::map<int,TPZMaterial *>::const_iterator it = matmap.begin(); it != matmap.end(); it++) {
        matindex[count++] = it->first;
    }
    buf.Write( matindex);
	buf.Write(fConnectIndex);
	buf.Write(fExternalLocIndex);
	buf.Write(fFatherToLocal);
    buf.Write(&fSingularConnect,1);
}

/**
 Read the element data from a stream
 */
void TPZSubCompMesh::Read(TPZStream &buf, void *context)
{
	TPZCompEl::Read(buf,context);
	TPZCompMesh::Read(buf,Mesh()->Reference());
    TPZCompMesh *mesh = (TPZCompMesh *) context;
    TPZManVector<int> matindex;
    buf.Read( matindex);
    int sz = matindex.size();
    for (int im=0; im<sz; im++) {
        MaterialVec()[matindex[im]] = mesh->MaterialVec()[matindex[im]];
    }
	buf.Read(fConnectIndex);
	buf.Read(fExternalLocIndex);
	buf.Read( fFatherToLocal);
    buf.Read(&fSingularConnect,1);
}

/**
 * Creates corresponding graphical element(s) if the dimension matches
 * graphical elements are used to generate output files
 * @param graphmesh graphical mesh where the element will be created
 * @param dimension target dimension of the graphical element
 */
void TPZSubCompMesh::CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension)
{
	int64_t nel = fElementVec.NElements();
	int64_t iel;
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
	int64_t iel, nelem = ElementVec().NElements();
	for(iel=0; iel<nelem; iel++)
	{
		TPZCompEl *cel = ElementVec()[iel];
		if(!cel) continue;
		TPZMaterial * mat = cel->Material();
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
#ifdef PZ_LOG
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
		LOGPZ_DEBUG(logger, sout.str())
	}
#endif
	if(numtrue && numfalse)
	{
		std::stringstream sout;
		sout << "A substructure should have either all elements computable or not numtrue " << numtrue << " numfalse " << numfalse;
		LOGPZ_ERROR(logger,sout.str())
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
	if (fConnectIndex.NElements() != fFatherToLocal.size()) {
		DebugStop();
	}
	int64_t numberexternal = fConnectIndex.NElements();
	int64_t i;
	for (i=0; i<numberexternal; i++) {
		if (fFatherToLocal.find(fConnectIndex[i]) == fFatherToLocal.end()) {
			DebugStop();
		}
	}
	// the number of external connects in the fExternalLocIndex should be size also
	int64_t nel = fExternalLocIndex.NElements();
	int64_t numext = 0;
	for (i=0; i<nel; i++) {
		if (fExternalLocIndex[i] != -1) {
			numext++;
		}
	}
	if (numext != numberexternal) {
		DebugStop();
	}
	std::map<int64_t,int64_t>::iterator it;
	for (it=fFatherToLocal.begin(); it!=fFatherToLocal.end(); it++) {
		if (fExternalLocIndex[it->second] == -1) {
			DebugStop();
		}
	}
	return true;
}

/** Verify if the material associated with the element is contained in the set */
bool TPZSubCompMesh::HasMaterial(const std::set<int> &materialids) const
{
    int nel = NElements();
    for (int el=0; el<nel ; el++) {
        TPZCompEl *cel = ElementVec()[el];        
        if (!cel) {
            continue;
        }
        bool has_material_Q = cel->HasMaterial(materialids);
        if (has_material_Q) {
            return true;
        }
    }
    return false;
}

int TPZSubCompMesh::NumberRigidBodyModes()
{
	if (fSingularConnect == -1) {
		return 0;
	}
	int64_t seqnum = fConnectVec[fSingularConnect].SequenceNumber();
	return fBlock.Size(seqnum);
	
}


/// Set the number of rigid body modes associated with the internal degrees of freedom
void TPZSubCompMesh::SetNumberRigidBodyModes(int nrigid, unsigned char lagrange)
{
	if (fSingularConnect == -1) {
        int nshape = nrigid;
        int nstate = 1;
        int order = 1;
		fSingularConnect = AllocateNewConnect(nshape,nstate,order);
		fConnectVec[fSingularConnect].IncrementElConnected();
        fConnectVec[fSingularConnect].SetLagrangeMultiplier(lagrange);
		int64_t extind = FatherMesh()->AllocateNewConnect(nshape,nstate,order);
		FatherMesh()->ConnectVec()[extind].IncrementElConnected();
		FatherMesh()->ConnectVec()[extind].SetLagrangeMultiplier(lagrange);
		int64_t next = fConnectIndex.NElements();
		fConnectIndex.Resize(next+1);
		fConnectIndex[next] = extind;
		fExternalLocIndex[fSingularConnect] = next;
        fFatherToLocal[extind] = fSingularConnect;
        ExpandSolution();
	}
	else if(fSingularConnect != -1 && nrigid >0 ) {
		int64_t seqnum = fConnectVec[fSingularConnect].SequenceNumber();
		fConnectVec[fSingularConnect].SetLagrangeMultiplier(lagrange);
		fBlock.Set(seqnum,nrigid);
        ExpandSolution();
		int64_t extind = fExternalLocIndex[fSingularConnect];
		TPZCompMesh *fathermesh = FatherMesh();
		if (fathermesh && extind < 0) {
			DebugStop();
		}
		while (fathermesh && extind > 0) {
			seqnum = fathermesh->ConnectVec()[extind].SequenceNumber();
			fathermesh->ConnectVec()[extind].SetLagrangeMultiplier(lagrange);
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

/** @brief adds the connect indexes associated with base shape functions to the set */
void TPZSubCompMesh::BuildCornerConnectList(std::set<int64_t> &connectindexes) const
{
    int nel = NElements();
    for (int el=0; el<nel ; el++) {
        TPZCompEl *cel = ElementVec()[el];
        if (!cel) {
            continue;
        }
        std::set<int64_t> locconind;
        cel->BuildCornerConnectList(locconind);
        std::set<int64_t>::iterator it;
        for (it=locconind.begin(); it != locconind.end(); it++) {
            int64_t index = *it;
            int64_t extlocindex = fExternalLocIndex[index];
            if (extlocindex != -1) {
                int64_t cornerind = fConnectIndex[extlocindex];
                connectindexes.insert(cornerind);
            }
        }
    }
}

/// return the index in the subcompmesh of a connect with index within the father
int64_t TPZSubCompMesh::InternalIndex(int64_t IndexinFather)
{
    if (fFatherToLocal.find(IndexinFather) == fFatherToLocal.end()) {
        DebugStop();
    }
    return fFatherToLocal[IndexinFather];
}

void TPZSubCompMesh::EvaluateError(TPZVec<REAL> &errors, bool store_errors){
    errors.Fill(0.);
    fAnalysis->PostProcessError(errors,store_errors);
    int NErrors = errors.size();
    if(store_errors)
    {
        int64_t index = Index();
        TPZFMatrix<STATE> &elvals = Mesh()->ElementSolution();
        if (elvals.Cols() < NErrors) {
            std::cout << "The element solution of the mesh should be resized before EvaluateError\n";
            DebugStop();
        }
        for (int ier=0; ier <NErrors; ier++) {
            elvals(index,ier) = errors[ier];
        }
    }

}


/**
 * Compute the residual norm of the internal equation
 * This method gives accurate results after CalcStiff or CalcResidual has been called
 */
REAL TPZSubCompMesh::InternalResidualNorm()
{
    REAL norm = 0.;
    if(!fAnalysis) DebugStop();
    TPZFMatrix<STATE> &rhs = fAnalysis->Rhs();
    // identify the internal equations
    int64_t ncon = ConnectVec().NElements();
    for (int64_t ic = 0; ic<ncon; ic++) {
        if(fExternalLocIndex[ic] != -1) continue;
        TPZConnect &c = ConnectVec()[ic];
        int nvar = c.NState()*c.NShape();
        if(c.HasDependency() || c.IsCondensed() || nvar == 0) continue;
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = Block().Position(seqnum);
        for(int iv=0; iv<nvar; iv++)
        {
            REAL rhsabs = abs(rhs(pos+iv,0));
            norm += rhsabs*rhsabs;
        }
    }
    return sqrt(norm);
    // sum the square of the residuals of the internal equations
    
}
