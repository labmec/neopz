//METHODS DEFINITIONS FOR CLASS COMPUTATIONAL MESH
// _*_ c++ _*_
#include "pzerror.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzconnect.h"
#include "pzbndcond.h"
#include "pzmaterial.h"
#include "pzsolve.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzelmat.h"
#include "TPZCompElDisc.h"
//#include "TPZEulerConsLaw.h"
#include "pztrnsform.h"
#include "pztransfer.h"
#include "pzavlmap.h"
#include "pzvec.h"
#include "pzadmchunk.h"

#include "pzmetis.h"

TPZCompMesh::TPZCompMesh (TPZGeoMesh* gr) : fElementVec(0),
  fMaterialVec(0), fConnectVec(0),fSolution(0,1),
  fBCConnectVec(0) {
  
  //Initializing class members
  
  fReference = gr;
  fChecked = 0;
  fName[0] = '\0';
  fName[126] = '\0';
  if(gr) {
    SetName( gr->Name() );
    gr->SetReference(this);
  }
  fBlock.SetMatrix(&fSolution);
  fSolutionBlock.SetMatrix(&fSolution);
}


TPZCompMesh::~TPZCompMesh() {

  // THIS NEEDS TO INCLUDE THE DELETION ROUTINES OF ALL ITEMS
  CleanUp();
}

void TPZCompMesh::CleanUp() {
  
  // THIS ROUTINE NEEDS TO INCLUDE THE DELETION OF THE LIST POINTERS
  TPZGeoMesh *ref = Reference();
  ref->ResetReference();
  LoadReferences();
  int i, nelem = NElements();
  for(i=0; i<nelem; i++) {
    TPZCompEl *el = fElementVec[i];
    if(el) delete el;
    fElementVec[i] = 0;
  }
  fElementVec.Resize(0);
  fElementVec.CompactDataStructure(1);
  fConnectVec.Resize(0);
  fConnectVec.CompactDataStructure(1);
  nelem = NMaterials();
  for(i=0; i<nelem; i++) {
    if(fMaterialVec[i]) delete fMaterialVec[i];
  }
  fMaterialVec.Resize(0);
  fMaterialVec.CompactDataStructure(1);
  fBCConnectVec.Resize(0);
  fBCConnectVec.CompactDataStructure(1);
  
  fBlock.SetNBlocks(0);
  fSolutionBlock.SetNBlocks(0);
  fSolution.Redim(0,0);
}

void TPZCompMesh::SetName (char *nm) {
  
  if (nm != NULL) {
    strncpy(fName,nm,126);
  } else {
    fName[0] = '\0';
  }
}

void TPZCompMesh::Print (ostream & out) {

	ComputeNodElCon();
  out << "\n\t\tCOMPUTABLE GRID INFORMATIONS:\n\n";
  out << "TITLE-> " << fName << "\n\n";

  out << "number of connects            = " << NConnects() << endl;
  out << "number of elements            = " << NElements() << endl;
  out << "number of materials           = " << NMaterials() << endl;
  out << "number of nodal bound cond    = " << NBCConnects() << endl;

  out << "\n\t Connect Information:\n\n";
  int i, nelem = NConnects();
  for(i=0; i<nelem; i++) {
    if(fConnectVec[i].SequenceNumber() == -1) {
    	if(fConnectVec[i].HasDependency()) {
      	cout << "TPZCompMesh::Print inconsistency of connect\n";
         cout << "Index " << i << ' ';
         fConnectVec[i].Print(*this,cout);
      }
    	continue;
    }
    out << "Index " << i << ' ';
    fConnectVec[i].Print(*this,out);
  }
  out << "\n\t Computable Element Information:\n\n";
  nelem = NElements();
  for(i=0; i<nelem; i++) {
    if(!fElementVec[i]) continue;
    TPZCompEl *el = fElementVec[i];
    out << "Index " << i << ' ';
    el->Print(out);
    if(!el->Reference()) continue;
    out << "\tReference Id = " << el->Reference()->Id() << endl;
  }
  out << "\n\tMaterial Information:\n\n";
  nelem = NMaterials();
  for(i=0; i<nelem; i++) {
    if(fMaterialVec[i] == 0) continue;
    TPZMaterial *mt = fMaterialVec[i];
    mt->Print(out);
  }
  out << "\nNodal boundary conditions\n";
  nelem = NBCConnects();
  for(i=0; i<nelem; i++) {
    if(!fBCConnectVec[i].fConnect) continue;
    fBCConnectVec[i].Print(*this,out);
  }
}

/**Insert a material object in the datastructure*/
int TPZCompMesh::InsertMaterialObject(TPZMaterial *mat) {
	if(!mat) return -1;
   TPZMaterial *othermat;
	int index;
   othermat = FindMaterial(mat->Id());
   if(!othermat) {
   	index = fMaterialVec.AllocateNewElement();
      fMaterialVec[index] = mat;
   } else {
		int nelem= NMaterials();
  		for(index=0; index<nelem; index++) {
	    	TPZMaterial *localmat = fMaterialVec[index];
   	 	if(!localmat) continue;
    		if(localmat == othermat) break;
      }
      if(index == nelem) {
      	PZError << "TPZCompMesh::InsertMaterialObject I don't understand\n";
      	index = fMaterialVec.AllocateNewElement();
      }
      if(othermat != mat) delete othermat;
      fMaterialVec[index] = mat;
   }
   return index;
}

TPZMaterial* TPZCompMesh::FindMaterial(int matid){	// find the material object with id matid
  int i, nelem= NMaterials();
  for(i=0; i<nelem; i++) {
    TPZMaterial *mat = fMaterialVec[i];
    if(!mat) continue;
    if(mat->Id() == matid) return mat;
  }
  return 0;
}

void TPZCompMesh::AutoBuild() {
  TPZAdmChunkVector<TPZGeoEl *> &elvec = Reference()->ElementVec();
  int i, nelem = elvec.NElements();
  int index;
  for(i=0; i<nelem; i++) {
    TPZGeoEl *gel = elvec[i];
    if(!gel) continue;
    if(!gel->HasSubElement()) {
      gel->CreateCompEl(*this,index);
    }
  }
  TPZCompEl *cel;
  TPZAdmChunkVector<TPZGeoElBC> &elbcvec = Reference()->BCElementVec();
  nelem = elbcvec.NElements();
  // InitializeBlock();
  for(i=0; i<nelem; i++) {
    if(!elbcvec[i].fBCElement) { 
      cel = elbcvec[i].fElement->CreateBCCompEl(elbcvec[i].fSide,elbcvec[i].fId,*this);
      if(cel) elbcvec[i].fBCElement = cel->Reference();
    }
  }
  InitializeBlock();
}

void TPZCompMesh::InitializeBlock() {
  ExpandSolution();
  CleanUpUnconnectedNodes();
}

void TPZCompMesh::ExpandSolution() {
  fBlock.Resequence();
  int ibl,nblocks = fBlock.NBlocks();
  
  TPZFMatrix OldSolution(fSolution);
  
  int cols = fSolution.Cols();
  fSolution.Redim(fBlock.Dim(),cols);
  int minblocks = nblocks < fSolutionBlock.NBlocks() ? nblocks : fSolutionBlock.NBlocks();

  int ic;
  for(ic=0; ic<cols; ic++) {
     for(ibl = 0;ibl<minblocks;ibl++) {
       int oldsize = fSolutionBlock.Size(ibl);
       int oldposition = fSolutionBlock.Position(ibl);
       int newsize = fBlock.Size(ibl);
       int newposition = fBlock.Position(ibl);
       int minsize = (oldsize < newsize) ? oldsize : newsize;
       int ieq;
       for(ieq=0; ieq<minsize; ieq++) {
         fSolution(newposition+ieq,ic) = OldSolution(oldposition+ieq,ic);
       }
     }
  }
  fSolutionBlock = fBlock;
}

void TPZCompMesh::LoadSolution(TPZFMatrix &mat){
  
  int nrow = mat.Rows();
  int ncol = mat.Cols();
  int i,j;
  for(j=0;j<ncol;j++) for(i=0;i<nrow;i++) fSolution(i,j) = (mat)(i,j);
  int nelem = NElements();
  TPZCompEl *cel;
  for(i=0; i<nelem; i++) {
    cel = fElementVec[i];
    if(!cel || !cel->IsInterpolated()) continue;
    cel->LoadSolution();
  }
}

void TPZCompMesh::LoadReferences() {
  
//	Reference()->ResetReference();
	Reference()->SetReference(this);
	int i, nelem = NElements();
	for(i=0; i<nelem; i++) {
		TPZCompEl *el = fElementVec[i];
		if(!el || !el->IsInterpolated()) continue;
		el->LoadElementReference();
/*  TPZGeoEl *gel = el->Reference();
    if(!gel) continue;
    gel->SetReference(el);
	*/
	}
}

void TPZCompMesh::CleanUpUnconnectedNodes() {
  ComputeNodElCon();
  int i, nelem = NConnects();
  int need = 0;
  for (i=0;i<nelem;i++) {
    TPZConnect &no = fConnectVec[i];
    if (no.SequenceNumber() == -1) continue;
    if (no.HasDependency() && no.NElConnected() == 0) {
      PZError << "TPZCompMesh::CleanUpUnconnectedNodes node has dependency\n";
      continue;
    }
    if (no.HasDependency()) break;
  }
  for (;i<nelem;i++) {
    TPZConnect &no = fConnectVec[i];
    if (!no.HasDependency() && no.SequenceNumber() != -1) need = 1;
  }
  int nblocks = fBlock.NBlocks();
  TPZManVector<int> permute(nblocks);
  for (i=0;i<nblocks;i++) permute[i] = i;

  if (need) {
    for(i=nelem-1; i>=0; i--) {
      TPZConnect &no = fConnectVec[i];
      if(no.SequenceNumber() == -1) continue;
      if(no.HasDependency()) {
	int seq = no.SequenceNumber();
	for (int j=0;j<nblocks;j++) if (permute[j] > permute[seq]) permute[j]--;
	permute[seq] = nblocks-1;
      }
    }

  }
//  cout << "permute to put the dependent connects to the back." << endl;
//  for (i=0;i<nblocks;i++) cout << permute[i] << " ";
//  cout << endl;
  int numfree = 0;
  for(i=0; i<nelem; i++) {
    TPZConnect &no = fConnectVec[i];
    if(no.SequenceNumber() == -1) continue;
    if(no.NElConnected() == 0) {
		if(no.HasDependency()) {
	      PZError << "TPZCompMesh::CleanUpUnconnectedNodes node has dependency\n";
			continue;
		}
      int seq = no.SequenceNumber();
      need = 1;
      for (int j=0;j<nblocks;j++) if (permute[j] > permute[seq]) permute[j]--;
      permute[seq] = nblocks-1;
      fBlock.Set(seq,0);
      no.SetSequenceNumber(-1);
      fConnectVec.SetFree(i);
      numfree++;
    }
  }
//  cout << "permute to put the free connects to the back\n";
//  for (i=0;i<nblocks;i++) cout << permute[i] << ' ';
//  cout << "need = " << need << endl;

  if (need) {
    Permute(permute);
    fBlock.SetNBlocks(nblocks-numfree);
  }
}

void TPZCompMesh::ComputeNodElCon() {
  
  int i, nelem = NConnects();
  for(i=0; i<nelem; i++) {
    TPZConnect &no = fConnectVec[i];
    if(no.SequenceNumber() == -1) continue;
    no.ResetElConnected();
  }
  
  
  TPZStack<int> nodelist;
  int numnod;
  // modified Philippe 22/7/97
  // in order to account for constrained nodes
  nelem = NElements();
  for(i=0; i<nelem; i++) {
    TPZCompEl *el = fElementVec[i];
    if(!el || !el->IsInterpolated()) continue;
//    if(!el) continue;
    nodelist.Resize(0);
    el->BuildConnectList(nodelist);
    numnod = nodelist.NElements();
    for (int in=0; in<numnod; ++in) {
      int dfnindex = nodelist[in];
      TPZConnect *dfn = &fConnectVec[dfnindex];
      dfn->IncrementElConnected();
    }
  }
}

//ofstream check("check.dat");
void TPZCompMesh::Assemble(TPZMatrix &stiffness,TPZFMatrix &rhs) {

  int iel;
  int numel = 0, nelem = NElements();
  TPZElementMatrix ek,ef;
  TPZManVector<int> destinationindex(0);
  TPZManVector<int> sourceindex(0);
  REAL stor1[1000],stor2[1000],stor3[100],stor4[100];
  ek.fMat = new TPZFMatrix(0,0,stor1,1000);
  ek.fConstrMat = new TPZFMatrix(0,0,stor2,1000);
  ef.fMat = new TPZFMatrix(0,0,stor3,100);
  ef.fConstrMat = new TPZFMatrix(0,0,stor4,100);

  for(iel=0; iel < nelem; iel++) {
    TPZCompEl *el = fElementVec[iel];
    if(!el) continue;
    //		int dim = el->NumNodes();
    el->CalcStiff(ek,ef);
  	 //ek.fMat->Print(out);
    //ef.fMat->Print();
    if(!(numel%20)) cout << endl << numel;
    cout << '*';
	cout.flush();
    numel++;

    if(!el->HasDependency()) {
      //ek.fMat->Print("stiff has no constraint",test);
      //ef.fMat->Print("rhs has no constraint",test);
      //test.flush();
      destinationindex.Resize(ek.fMat->Rows());
      int destindex = 0;
      int numnod = ek.NConnects();
      for(int in=0; in<numnod; in++) {
         int npindex = ek.ConnectIndex(in);
         TPZConnect &np = fConnectVec[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = fBlock.Position(blocknumber);
         int ndf = fBlock.Size(blocknumber);
         for(int idf=0; idf<ndf; idf++) {
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      stiffness.AddKel(*ek.fMat,destinationindex);
      rhs.AddFel(*ef.fMat,destinationindex);                 //  ??????????? Erro
    } else {
      // the element has dependent nodes
      el->ApplyConstraints(ek,ef);
      //ek.fMat->Print("stif no constraint",test);
      //ek.fConstrMat->Print("stif constrained",test);
      //ef.fMat->Print("rhs no constraint",test);
      //ef.fConstrMat->Print("rhs constrained",test);
      //test.flush();
      //test << "sum of columns\n";
      int destindex = 0;
      int fullmatindex = 0;
      destinationindex.Resize(ek.fConstrMat->Rows());
      sourceindex.Resize(ek.fConstrMat->Rows());
      int numnod = ek.fConstrConnect.NElements();
      for(int in=0; in<numnod; in++) {
         int npindex = ek.fConstrConnect[in];
         TPZConnect &np = fConnectVec[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = fBlock.Position(blocknumber);
      	int ndf = fBlock.Size(blocknumber);
         if(np.HasDependency()) {
           fullmatindex += ndf;
           continue;
         }
         for(int idf=0; idf<ndf; idf++) {
           sourceindex[destindex] = fullmatindex++;
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      sourceindex.Resize(destindex);
      destinationindex.Resize(destindex);
      stiffness.AddKel(*ek.fConstrMat,sourceindex,destinationindex);
      rhs.AddFel(*ef.fConstrMat,sourceindex,destinationindex);
/*
if(ek.fConstrMat->Decompose_LU() != -1) {
    el->ApplyConstraints(ek,ef);
    ek.Print(*this,check);
    check.flush();
}
*/
    }
  }//fim for iel
  // Add contribution of nodal boundary conditions
  // Philippe 8/4/97
  //test.flush();
  static REAL BigNumber = 1.e9;
  nelem = fBCConnectVec.NElements();
  for(iel=0; iel<nelem; iel++) {
    TPZConnectBC &nodbnd = fBCConnectVec[iel];
    TPZConnect *mp = nodbnd.fConnect;
    if(!mp) continue;
    if(mp->HasDependency()) {
      PZError << "TPZCompMesh::Assemble a boundary node has dependency!\n";
      continue;
    }
    int seqnum = mp->SequenceNumber();
    if(seqnum < 0) continue;//Cedric 11/04/99
    int blocksize = fBlock.Size(seqnum);
    TPZBndCond *bc = nodbnd.fBC;
    int firsteq = fBlock.Position(seqnum);
    int numvar = (bc->Val1()).Rows();
    if(numvar != blocksize) {
      PZError << "TPZCompMesh::Assemble the block size and size of the "
	"boundary condition are incompatible\n";
      continue;
    }
    int type = bc->Type();
    switch(type) {
       case 0: {	// condicao dirichlet
         REAL prevval;
         int ieq;
         for(ieq = firsteq; ieq<firsteq+numvar; ieq++) {
            prevval = stiffness.GetVal(ieq,ieq);
            prevval += BigNumber;
            stiffness.PutVal(ieq,ieq,prevval);
            rhs(ieq,0) += BigNumber*(bc->Val2())(ieq-firsteq,0);
         }
         break;
       }
       case 1: { 	// condicao Neumann
         int ieq;
         for(ieq = firsteq; ieq<firsteq+numvar; ieq++) {
            rhs(ieq,0) += (bc->Val2())(ieq-firsteq,0);
         }
         break;
       }
       case 2: {	// condicao mixta
            REAL prevval;
            int ieq;
            for(ieq = 0; ieq<numvar; ieq++) {
               int ieqnumber = firsteq+ieq;
               rhs(ieqnumber,0) += (bc->Val2())(ieq,0);
               int jeq;
               for(jeq = 0; jeq < numvar; jeq++) {
                 int jeqnumber = firsteq+jeq;
                 prevval = stiffness.GetVal(ieqnumber,jeqnumber);
                 prevval += bc->Val1()(ieq,jeq);
                 stiffness.PutVal(ieqnumber,jeqnumber,prevval);
               }
            }
            break;
       }
    }
  }

  cout << endl;
}



//ofstream check("check.dat");
void TPZCompMesh::Assemble(TPZFMatrix &rhs) {

  int iel;
  int numel = 0, nelem = NElements();
  TPZElementMatrix ef;
  TPZManVector<int> destinationindex(0);
  TPZManVector<int> sourceindex(0);
  REAL stor3[100],stor4[100];
  ef.fMat = new TPZFMatrix(0,0,stor3,100);
  ef.fConstrMat = new TPZFMatrix(0,0,stor4,100);

  for(iel=0; iel < nelem; iel++) {
    TPZCompEl *el = fElementVec[iel];
    if(!el) continue;
    //		int dim = el->NumNodes();
    el->CalcResidual(ef);
  	 //ek.fMat->Print(out);
    //ef.fMat->Print();
//    if(!(numel%20)) cout << endl << numel;
//    cout << '*';
//	cout.flush();
    numel++;

    if(!el->HasDependency()) {
      //ek.fMat->Print("stiff has no constraint",test);
      //ef.fMat->Print("rhs has no constraint",test);
      //test.flush();
      destinationindex.Resize(ef.fMat->Rows());
      int destindex = 0;
      int numnod = ef.NConnects();
      for(int in=0; in<numnod; in++) {
         int npindex = ef.ConnectIndex(in);
         TPZConnect &np = fConnectVec[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = fBlock.Position(blocknumber);
         int ndf = fBlock.Size(blocknumber);
         for(int idf=0; idf<ndf; idf++) {
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      rhs.AddFel(*ef.fMat,destinationindex);                 //  ??????????? Erro
    } else {
      // the element has dependent nodes
      el->ApplyConstraints(ef);
      //ek.fMat->Print("stif no constraint",test);
      //ek.fConstrMat->Print("stif constrained",test);
      //ef.fMat->Print("rhs no constraint",test);
      //ef.fConstrMat->Print("rhs constrained",test);
      //test.flush();
      //test << "sum of columns\n";
      int destindex = 0;
      int fullmatindex = 0;
      destinationindex.Resize(ef.fConstrMat->Rows());
      sourceindex.Resize(ef.fConstrMat->Rows());
      int numnod = ef.fConstrConnect.NElements();
      for(int in=0; in<numnod; in++) {
         int npindex = ef.fConstrConnect[in];
         TPZConnect &np = fConnectVec[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = fBlock.Position(blocknumber);
      	int ndf = fBlock.Size(blocknumber);
         if(np.HasDependency()) {
           fullmatindex += ndf;
           continue;
         }
         for(int idf=0; idf<ndf; idf++) {
           sourceindex[destindex] = fullmatindex++;
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      sourceindex.Resize(destindex);
      destinationindex.Resize(destindex);
      rhs.AddFel(*ef.fConstrMat,sourceindex,destinationindex);
/*
if(ek.fConstrMat->Decompose_LU() != -1) {
    el->ApplyConstraints(ek,ef);
    ek.Print(*this,check);
    check.flush();
}
*/
    }
  }//fim for iel
  // Add contribution of nodal boundary conditions
  // Philippe 8/4/97
  //test.flush();
  static REAL BigNumber = 1.e9;
  nelem = fBCConnectVec.NElements();
  for(iel=0; iel<nelem; iel++) {
    TPZConnectBC &nodbnd = fBCConnectVec[iel];
    TPZConnect *mp = nodbnd.fConnect;
    if(!mp) continue;
    if(mp->HasDependency()) {
      PZError << "TPZCompMesh::Assemble a boundary node has dependency!\n";
      continue;
    }
    int seqnum = mp->SequenceNumber();
    if(seqnum < 0) continue;//Cedric 11/04/99
    int blocksize = fBlock.Size(seqnum);
    TPZBndCond *bc = nodbnd.fBC;
    int firsteq = fBlock.Position(seqnum);
    int numvar = (bc->Val1()).Rows();
    if(numvar != blocksize) {
      PZError << "TPZCompMesh::Assemble the block size and size of the "
	"boundary condition are incompatible\n";
      continue;
    }
    int type = bc->Type();
    switch(type) {
       case 0: {	// condicao dirichlet
         int ieq;
         for(ieq = firsteq; ieq<firsteq+numvar; ieq++) {
            rhs(ieq,0) += BigNumber*(bc->Val2())(ieq-firsteq,0);
         }
         break;
       }
       case 1: { 	// condicao Neumann
         int ieq;
         for(ieq = firsteq; ieq<firsteq+numvar; ieq++) {
            rhs(ieq,0) += (bc->Val2())(ieq-firsteq,0);
         }
         break;
       }
       case 2: {	// condicao mixta
            int ieq;
            for(ieq = 0; ieq<numvar; ieq++) {
               int ieqnumber = firsteq+ieq;
               rhs(ieqnumber,0) += (bc->Val2())(ieq,0);
            }
            break;
       }
    }
  }

//  cout << endl;
}

int TPZCompMesh::NEquations() {

  int neq = 0;
  int i, ncon = NConnects();
  for(i=0; i<ncon; i++) {
    TPZConnect &df = fConnectVec[i];
    if(df.HasDependency() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
    int seqnum = df.SequenceNumber();
    neq += fBlock.Size(seqnum);
  }
  return neq;
}

/** Este metodo nao e apropriado para aproximantes descontinuas */
int TPZCompMesh::BandWidth() {
  
  int bw = 0;
  TPZStack<int> connectlist;
  // modified Philippe 24/7/97
  // in order to take dependent nodes into account
  
  int i, nelem = NElements();
  for(i=0; i<nelem; i++) {
    connectlist.Resize(0);
    TPZCompEl *el = fElementVec[i];
    if(!el || !(el->IsInterpolated())) continue;
    el->BuildConnectList(connectlist);
    int nnod = connectlist.NElements();
    if(!nnod) continue;
    int ifirstnode = 0;
    TPZConnect *np = &fConnectVec[connectlist[ifirstnode++]];
    while(ifirstnode < nnod && (np->HasDependency() || !fBlock.Size(np->SequenceNumber()))) {
      np = &fConnectVec[connectlist[ifirstnode++]];
    }
    int ibl = np->SequenceNumber();
    int loweq = fBlock.Position(ibl);
    int higheq = loweq+fBlock.Size(ibl)-1;
    for(int n=ifirstnode;n<nnod;n++) {
      np = &fConnectVec[connectlist[n]];
      if(np->HasDependency()) continue;
      int ibl = np->SequenceNumber();
      if(!fBlock.Size(ibl)) continue;
      int leq = fBlock.Position(ibl);
      int heq = leq+fBlock.Size(ibl)-1;
      loweq = (loweq > leq) ? leq : loweq;
      higheq = (higheq < heq) ? heq : higheq;
    }
    int elbw = higheq - loweq;
    bw = (bw < elbw) ? elbw : bw;
  }
  return bw;
}

void TPZCompMesh::Skyline(TPZVec<int> &skyline) {

   TPZStack<int> connectlist;
   // modified Philippe 24/7/97
   // in order to take dependent nodes into account

   int neq = NEquations();
   skyline.Resize(neq);
   skyline.NElements();
//cout << "Element Equations";
//int eleq=0;
   int i, n, l, nelem = NElements();
   for(i=0; i<neq; i++) skyline[i] = i;
   for(i=0; i<nelem; i++) {
      TPZCompEl *el = fElementVec[i];
      if(!el || !el->IsInterpolated()) continue;
//      if(!el) continue;
      connectlist.Resize(0);
      el->BuildConnectList(connectlist);
      int nnod = connectlist.NElements();
      if(!nnod) continue;
      int ifirstnode = 0;
      TPZConnect *np = &fConnectVec[connectlist[0]];
      while(ifirstnode < nnod && np->HasDependency()) {
         ifirstnode++;
         np = &fConnectVec[connectlist[ifirstnode]];
      }
      int ibl = np->SequenceNumber();
      int loweq = fBlock.Position(ibl);
      int higheq = loweq+fBlock.Size(ibl)-1;
      for(n=ifirstnode;n<nnod;n++) {
         np = &fConnectVec[connectlist[n]];
         if(np->HasDependency()) continue;
         int ibl = np->SequenceNumber();
         int leq = fBlock.Position(ibl);
         int heq = leq+fBlock.Size(ibl)-1;
//for(int _eq=leq; _eq<= heq; _eq++) {
//   if((eleq%20==0)) cout << endl;
//   cout << _eq << ' ';
//   eleq++;
//}
         loweq = (loweq > leq) ? leq : loweq;
         higheq = (higheq < heq) ? heq : higheq;
      }
//cout << endl;
      for(n=ifirstnode;n<nnod;n++) {
         np = &fConnectVec[connectlist[n]];
         if(np->HasDependency()) continue;
         int ibl = np->SequenceNumber();
         int leq = fBlock.Position(ibl);
         int heq = leq+fBlock.Size(ibl);
         for(l=leq;l<heq;l++) {
            skyline[l] = skyline[l] < loweq ? skyline[l] : loweq;
         }
      }
   }
}

void TPZCompMesh::BuildTransferMatrix(TPZCompMesh &coarsemesh, TPZTransfer &transfer) {

  TPZBlock &localblock = Block();
  TPZBlock &coarseblock = coarsemesh.Block();
  // adapt the block size of the blocks, dividing by the number of variables
  //  of the material
  int i, nmat = NMaterials();
  TPZMaterial *mat = 0;
  for(i=0; i<nmat; i++) {
    mat = fMaterialVec[i];
    if(mat) break;
  }
  if(!mat) {
    PZError << "TPZCompMesh::BuildTransferMatrix, no material object found\n";
    return;
  }
  int nvar = mat->NStateVariables();
  int dim = mat->Dimension();
  
  transfer.SetBlocks(localblock,coarseblock,nvar,NIndependentConnects(),coarsemesh.NIndependentConnects());
  Reference()->ResetReference();
  coarsemesh.LoadReferences();
  int nelem = NElements();
  for(i=0; i<nelem; i++) {
    TPZInterpolatedElement *locel = 0;
    if(fElementVec[i] && fElementVec[i]->IsInterpolated())
      locel = (TPZInterpolatedElement *) fElementVec[i];
    if(!locel || !locel->IsInterpolated()) continue;
    if(locel->Dimension() != dim) continue;
    TPZGeoEl *locgel = locel->Reference();
    TPZGeoEl *coarsegel = locgel;
    if(!locgel) {
      cout << "TPZCompMesh::BuildTransferMatrix is not implemented for super elements\n";
      continue;
    }
    TPZInterpolatedElement *coarsel = 0;
    while(coarsegel && !coarsegel->Reference()) {
      coarsegel = coarsegel->Father();
    }
    if(!coarsegel) {
      cout << "TPZCompMesh::BuildTransferMatrix corresponding coarse element not found\n";
      locel->Print(cout);
      continue;
    }
    if(coarsegel->Reference()->IsInterpolated())
      coarsel = (TPZInterpolatedElement *) coarsegel->Reference();
    
    if(!coarsel) continue;
    
    if(coarsel->Mesh() != &coarsemesh) {
      cout << "TPZCompMesh::BuildTransferMatrix is not implemented for transfers"
	" between superelements\n";
      continue;
    }
    TPZTransform t(coarsel->Dimension());
    locgel->BuildTransform2(locel->NConnects()-1,coarsegel,t);
    locel->BuildTransferMatrix(*coarsel,t,transfer);
  }
}

void TPZCompMesh::CreateConnectBC() {
  TPZGeoMesh *geo = Reference();
  TPZAdmChunkVector<TPZGeoNodeBC> *geobndcondvec = &geo->BCNodeVec();
  int ibc, nnodebc = geobndcondvec->NElements();
  for(ibc=0; ibc<nnodebc; ibc++) {
    TPZGeoNodeBC &gbc = (*geobndcondvec)[ibc];
    int bcnumber = gbc.fBCId;
    TPZMaterial *bc = FindMaterial(bcnumber);
    if(!bc) continue;
    // find a geometric element which has a
    //      corresponding computational element
    TPZGeoElSide gel(gbc.fGeoEl,gbc.fGeoElSide);
    TPZStack<TPZGeoElSide> neighbourset;
    gel.AllNeighbours(neighbourset);
    int in=0,nneigh;
    nneigh = neighbourset.NElements();
    while(in<nneigh && !neighbourset[in].Reference().Exists()) in++;
//     neighbour = gel.Neighbour();
//     while(neighbour.Exists() && neighbour != gel && !neighbour.Reference().Exists()) {
//       neighbour = neighbour.Neighbour();
//     }
//     if(!neighbour.Exists() || ! neighbour.Reference().Exists()) continue;
    if(in == nneigh) continue;
    TPZCompElSide cel = neighbourset[in].Reference();
    TPZConnect *df = &cel.Element()->Connect(neighbourset[in].Side());
    if(!df) continue;
    TPZConnectBC dfbc(df,(TPZBndCond *) bc);
    int key = BCConnectVec().AllocateNewElement();
    BCConnectVec()[key] = dfbc;
  }
}

/*
void TPZCompMesh::ComputeConnecttoElGraph(TPZVec<int> &firstel, TPZVec<int> &connectelgraph){
  int connectstackstore[50];
  TPZStack<int> connectstack(connectstackstore,50);
  int i, ncon = NConnects();
  firstel.Resize(ncon+1);
  firstel[0] = 0;
  for(i=0; i<ncon; i++) {
    TPZConnect &c = fConnectVec[i];
    int seqnum = c.SequenceNumber();
    if(seqnum == -1) {
      firstel[i+1] = firstel[i];
    } else {
      firstel[i+1] = firstel[i] + c.NElConnected();
    }
  }
  connectelgraph.Resize(firstel[ncon]);
  connectelgraph.Fill(-1);
  int nelem = NElements();
  for(i=0; i<nelem; i++) {
    TPZCompEl *el = fElementVec[i];
    if(!el) continue;
    connectstack.Resize(0);
    el->BuildConnectList(connectstack);
    int in;
    ncon = connectstack.NElements();
    for(in=0; in<ncon; in++) {
      int ic = connectstack[in];
      int first = firstel[ic];
      int last = firstel[ic+1];
      while(connectelgraph[first] != -1 && first < last) first++;
      if(first == last) {
	PZError << "TPZCompMesh::ComputeConnecttoElGraph wrong data structure\n";
	continue;
      }
      connectelgraph[first] = i;
    }
  }
}
*/
int TPZCompMesh::NIndependentConnects() {
  int i, ncon = NConnects();
  int NIndependentConnects = 0;
  for(i=0; i<ncon; i++) {
      TPZConnect &c = fConnectVec[i];
      if(c.HasDependency() || c.SequenceNumber() == -1) continue;
      NIndependentConnects++;
  }
  return NIndependentConnects;
}

void TPZCompMesh::ComputeElGraph(TPZStack<int> &elgraph, TPZVec<int> &elgraphindex){
  int i, ncon;
  TPZStack<int> connectstack;
  int nelem = NElements();
  elgraphindex.Resize(nelem+1);
  elgraphindex[0] = 0;
  elgraph.Resize(0);
  int curel=0;
  for(i=0; i<nelem; i++) {
    TPZCompEl *el = fElementVec[i];
    if(!el || !el->IsInterpolated()){
      elgraphindex[curel+1]=elgraph.NElements();
      curel++;
      continue;
    }
    connectstack.Resize(0);
    el->BuildConnectList(connectstack);
    int in;
    ncon = connectstack.NElements();
    for(in=0; in<ncon; in++) {
      int ic = connectstack[in];
      TPZConnect &c = fConnectVec[ic];
      if(c.HasDependency()) continue;
      elgraph.Push(c.SequenceNumber());
    }
    elgraphindex[curel+1]=elgraph.NElements();
    curel++;
  }
  elgraphindex.Resize(curel+1);
}

void TPZCompMesh::Divide(int index,TPZVec<int> &subindex,int interpolate) {

  TPZCompEl * el = fElementVec[index];
  if (!el) {
    PZError << "TPZCompMesh::Divide element not found index = " << index << endl;
    subindex.Resize(0);
    return;
  }
  el->Divide(index,subindex,interpolate);
}

void TPZCompMesh::CoarsenDisc(TPZVec<int> &elements, int &index) {
  int i,nelem = elements.NElements();
  if(!nelem) {
    index = -1;
    return;
  }

  TPZGeoEl *father = 0;
  TPZCompEl *cel = fElementVec[elements[0]];

  if(!cel->IsInterpolated() ) {
    index = -1;
    return;
  }

  TPZCompElDisc *el;
  el = (TPZCompElDisc *) fElementVec[elements[0]];
  if(el) father = el->Reference()->Father();
  if(!father) {
    index = -1;
    return;
  }

  for(i=1;i<nelem;i++) {
    if(!fElementVec[elements[i]]) {
      index = -1;
      return;
    }
    TPZGeoEl *father2 = fElementVec[elements[i]]->Reference()->Father();
    if(!father2 || father != father2) {
      index = -1;
      return;
    }
  }

  if(nelem != father->NSubElements())
    cout << "TPZCompEl::Coarsen : incomplete list of elements sons\n";

  for(i=0; i<nelem; i++) {
    el = (TPZCompElDisc *) fElementVec[elements[i]];
    el->RemoveInterfaces();
    el->Reference()->ResetReference();
  }

  for(i=0; i<nelem; i++) {
    el = (TPZCompElDisc *) fElementVec[elements[i]];
    delete el;
  }
  father->CreateCompEl(*this,index);
}

void TPZCompMesh::Coarsen(TPZVec<int> &elements, int &index) {
  int i,nelem = elements.NElements();
  if(!nelem) {
    index = -1;
    return;
  }

  TPZGeoEl *father = 0;
  TPZCompEl *cel = fElementVec[elements[0]];
  if(!cel->IsInterpolated()) {
    index = -1;
    return;
  }
  TPZInterpolatedElement *el;
  el = (TPZInterpolatedElement *)fElementVec[elements[0]];
  if(el) father = el->Reference()->Father();
  if(!father) {
    index = -1;
    return;
  }

  for(i=1;i<nelem;i++) {
    if(!fElementVec[elements[i]]) {
      index = -1;
      return;
    }
    TPZGeoEl *father2 = fElementVec[elements[i]]->Reference()->Father();
    if(!father2 || father != father2) {
      index = -1;
      return;
    }
  }

  if(nelem != father->NSubElements())
    cout << "TPZCompEl::Coarsen : incomplete list of elements sons\n";

  for(i=0; i<nelem; i++) {
    el = (TPZInterpolatedElement *)fElementVec[elements[i]];
    if(el->CanToHaveInterface()) el->DeleteInterfaces();
//
    el->RemoveSideRestraintsII(TPZInterpolatedElement::EDelete);
    el->Reference()->ResetReference();
  }

  for(i=0; i<nelem; i++) {
    el = (TPZInterpolatedElement *)fElementVec[elements[i]];
    delete el;
  }
  father->CreateCompEl(*this,index);
}

/**ExpandSolution must be called before calling this*/
void TPZCompMesh::Permute(TPZVec<int> &permute) {

  ExpandSolution();
//   if (permute.NElements() != fBlock.NBlocks()) {
//     PZError << "TPZCompMesh::Permute : permute vector size not equal to fBlock size\n";
//   }
  int i,j;
  int permutenel = permute.NElements();
  for (i = 0; i < permutenel; i++) fBlock.Set(permute[i],fSolutionBlock.Size(i));
  fBlock.Resequence();
  if (fSolution.Rows() != 0) {
    TPZFMatrix	newsol(fSolution);
    for (i=0;i<fBlock.NBlocks();i++) {
      int oldpos = fSolutionBlock.Position(i);
      int newpos;
      if(i < permutenel) {
	      newpos = fBlock.Position(permute[i]);
      } else {
	      newpos = fBlock.Position(i);
      }
      for (j=0;j<fSolutionBlock.Size(i);j++) fSolution(newpos+j,0) = newsol(oldpos+j,0);
    }    //a sol. inicial esta em newsol
  }

  fSolutionBlock = fBlock;
  int ncon = NConnects();
  for(i=0; i<ncon; i++) {
    TPZConnect &df = fConnectVec[i];
    int seqnum = df.SequenceNumber();
    if(seqnum == -1) continue;
    if(seqnum < permutenel) df.SetSequenceNumber(permute[seqnum]);
  }
}

void TPZCompMesh::ConnectSolution(ostream & out) {

  out << "\n\t\tCONNECT INDEX SOLUTION:\n\n";
  int i;
  int ncon = NConnects();
  for(i=0; i<ncon; i++) {
        out << i << ") ";
        TPZConnect &df = ConnectVec()[i];
        int seqnum = df.SequenceNumber();
		  if(df.NElConnected()==0) {
          out << "free node" << endl;
        } else if (seqnum < 0 || Block().Size(seqnum)==0) {
           out << "non solution connect" << endl;
        } else {
        		int pos = Block().Position(seqnum);
           for(int j=0;j<Block().Size(seqnum);j++)
	           	out << Solution()(pos+j,0) << "  ";
			  	out << endl;
        }
  }
}

void TPZCompMesh::EvaluateError(
       void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
       REAL &true_errorSum,REAL &L2_errorSum,REAL &estimateSum) {

	L2_errorSum = 0.;
   true_errorSum = 0.;
   estimateSum = 0.;
   REAL true_error = 0.;
   REAL estimate = 0.;
   REAL L2_error = 0.;
   TPZBlock *flux = 0;
   TPZCompEl *cel;
   //soma de erros sobre os elementos
   for(int el=0;el< fElementVec.NElements();el++) {
      cel = fElementVec[el];
      if(!cel || !cel->IsInterpolated() || cel->Material()->Id() < 0) continue;
		cel->EvaluateError(fp,true_error,L2_error,flux,estimate);
      if(fabs(true_error) > 100. || fabs(L2_error) > 100. || fabs(estimate) > 100.) {
      	int key;
         cout << "\n Error overflow \n";
      	cin >> key;
      }
      true_errorSum += true_error*true_error;
      L2_errorSum += L2_error*L2_error;
      estimateSum += estimate*estimate;
   }
   true_errorSum = sqrt(true_errorSum);
   L2_errorSum = sqrt(L2_errorSum);
   estimateSum = sqrt(estimateSum);
}

void TPZCompMesh::AdjustBoundaryElements() {
  int changed = 1;
  while(changed) {
    changed = 0;
    int nel = fElementVec.NElements();
    int el;
    TPZVec<int> subelindex;
    for(el=0; el<nel; el++) {
      TPZStack<TPZCompElSide> elvec;
      TPZCompEl *elp = fElementVec[el];
      if(!elp || !elp->IsInterpolated()) continue;
      TPZMaterial *mat = elp->Material();
      if(mat && mat->Id() >= 0) continue;
      if(!elp->IsInterpolated()) continue;
      int nsides = elp->Reference()->NSides();
      int is;
      int nc = elp->Reference()->NCornerNodes();
      for(is=nc; is<nsides; is++) {
	TPZCompElSide elpside(elp,is);
	if(elp->Reference()->SideDimension(is) == 0) continue;
	elvec.Resize(0);
	elpside.HigherLevelElementList(elvec,0,0);
	if(elvec.NElements()) {
	  TPZGeoElSide fatherside = elvec[0].Reference();//el BC
	  int face = elp->Reference()->NSides() - 1;
	  while(fatherside.Exists()) {
	    if(elp->Reference()->NeighbourExists(face,fatherside)) break;
	    fatherside = fatherside.Father2();
	  }
	  if(fatherside.Exists()) {
	    Divide(el,subelindex,0);
	    changed = 1;
	    break;
	  }
	}
	if(is == nsides-1) {
	  TPZInterpolatedElement *elpint = dynamic_cast <TPZInterpolatedElement *> (elp);
	  if(!elpint) continue;
	  int porder = elpint->PreferredSideOrder(is);
	  int maxorder = 0;
	  elvec.Resize(0);
	  elpside.EqualLevelElementList(elvec,0,0);
	  int eq;
	  for(eq=0; eq<elvec.NElements(); eq++) {
	    TPZInterpolatedElement *eqel = dynamic_cast<TPZInterpolatedElement *> (elvec[eq].Element());
	    int eqside = elvec[eq].Side();
	    if(!eqel) continue;
	    if(maxorder < eqel->PreferredSideOrder(eqside)) maxorder =  eqel->PreferredSideOrder(eqside);
	  }
	  if(porder < maxorder) {
	    elpint->PRefine(maxorder);
	    changed = 1;
	  }
	}
      }
    }
  }
}

int TPZCompMesh::PutinSuperMesh (int local, TPZCompMesh *super){
	if (super != this) return -1;
	else return local;
}

int TPZCompMesh::GetFromSuperMesh (int superind, TPZCompMesh *super){
	if (super != this) return -1;
	else return superind;
}



/*void TPZCompMesh::AdjustBoundaryElements() {
	int changed = 1;
	while(changed) {
		changed = 0;
		int nel = fElementVec.NElements();
		int el;
		TPZVec<int> subelindex;
		for(el=0; el<nel; el++) {
			TPZStack<TPZCompElSide> elvec(10);
			TPZCompEl *elp = fElementVec[el];
			if(!elp) continue;
			TPZMaterial *mat = elp->Material();
			if(mat && mat->Id() >= 0) continue;
			if(!elp->IsInterpolated()) continue;
         //TPZInterpolatedElement *cel = (TPZInterpolatedElement *) elp;
//         cel->PRefine(99);
// coucou!!!
			int nsides = elp->Reference()->NSides();
			int is;
         int nc = elp->Reference()->NCornerNodes();//Cedric 21/05/99
			//for(is=0; is<nsides; is++) {
			for(is=nc; is<nsides; is++) {//Cedric 02/04/00
				TPZCompElSide elpside(elp,is);
				elvec.Resize(0);
				elpside.HigherLevelElementList(elvec,0,0);
				if(elvec.NElements()) {
					Divide(el,subelindex,0);
					changed = 1;
					break;
				}
			}
		}
	}
}*/
/*
void TPZCompMesh::hp_Adaptive_Mesh_Design(void (*fExact)(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv),
               		REAL nadmerror,int niter,TPZVec<REAL> &singularpoint,TPZAnalysis &an,ostream &out) {
   int index;
	int Q = 2;//valor dado > 1 na procura de hn para elementos proximos ao ponto singular
   int Nc = 4;//number of layers
   int iter = 0;//iteração atual
   REAL tol = .05;//tolerancia na aproximação do nova ordem do elemento
   int maxiter = 20;//maximo numero de iterações do algoritmo do ponto fixo
   REAL w = 0.3;//fator de relaxação
   cout << "\n\nStep  0\n";
	an.Run(out);//solução malha inicial
   Print(out);
   Reference()->Print(out);
   an.Print("FEM SOLUTION ",out);

	while(++iter < niter || MaximLocalError(fExact) > AdmLocalError(fExact,nadmerror) ) {

   	TPZCompEl *singel = SingularElement(singularpoint,ElementVec());//elemento que contem o ponto singular
      TPZStack<int> indexlist(0);
      Step3(singel,singularpoint,Q,Nc,indexlist);//distribui a ordem dos els. proximos ao ponto singular
      Step4(fExact,w,nadmerror,tol,maxiter,indexlist);//processa os restantes elementos
      InitializeBlock();
      TPZAnalysis step_an(this,out);
      int numeq = NEquations();
      TPZFMatrix *stiff = new TPZFMatrix(numeq,numeq);
      step_an.SetMatrix(stiff);
      step_an.Solver().SetDirect(ELU);
      cout << "\n\nStep " << (iter+1) << endl;
      step_an.Run(out);
      Print(out);
      Reference()->Print(out);
      an.Print("FEM SOLUTION ",out);
   }
}

TPZCompEl *TPZCompMesh::SingularElement(TPZVec<REAL> &point,TPZAdmChunkVector<TPZCompEl *> &listel) {
//point deve ser interior a algum elemento da malha
	int nel = listel.NElements();
   TPZInterpolatedElement *intel,*intelkeep=0;
   for(int index=0;index<nel;index++) {
   	intel = (TPZInterpolatedElement *) listel[index];
      TPZVec<REAL> coord(3,0.),result(3);
      REAL dist = 0.,mindist = 100.0;
      if(intel) {
         for(int con=0;con<intel->NCornerConnects();con++) {
            intel->Reference()->X(coord,result);
            for(int i=0;i<3;i++)
               dist += pow( result[i] - intel->Reference()->NodePtr(con)->Coord(i) , 2.0);
            if(dist < mindist) {
               mindist = dist;//dist = sqrt(dist);
               intelkeep = intel;
            }//if
         }//for con
      }
   }//for index
   return intelkeep;
}

REAL TPZCompMesh::h_Parameter(TPZCompEl *cel) {

 	REAL h = 0.,cicjdist;
   TPZGeoEl *gel = cel->Reference();
   int nconn = gel->NCornerNodes();
   for(int conni=0;conni<nconn;conni++) {
   	for(int connj=conni;connj<nconn;connj++) {
         cicjdist = 0.;
         for(int coordi=0;coordi<3;coordi++) {
         	REAL coor1 = gel->NodePtr(conni)->Coord(coordi);
            REAL coor2 = gel->NodePtr(connj)->Coord(coordi);
            cicjdist += pow( coor1 - coor2, 2.0);
         }
         cicjdist = sqrt(cicjdist);
         if(h < cicjdist) h = cicjdist;
      }
   }
   return h;
}

void TPZCompMesh::Step3(TPZCompEl *cel,TPZVec<REAL> &point,int Q,int Nc,TPZStack<int> &indexlist) {

      int sum = (1+Q);
      int nc = 2;
      while(nc < Nc+1) {
          sum += pow(Q,nc);
          nc++;
      }
      REAL h = h_Parameter(cel);
      REAL hn = h/sum;
      if(hn > 1.3*h) hn = 2.0*h*hn / (h + hn);
      REAL hsub = 100.0;
      TPZCompEl *locel = cel;
      //obter um subelemento que contem o ponto singular e tem tamanho <= hn
      while(hsub > hn) {
      	TPZVec<int> indexsubs;
         int index = locel->Index();
         locel->Divide(index,indexsubs,1);
         int nsub = indexsubs.NElements();
         TPZAdmChunkVector<TPZCompEl *> listsub(nsub);
         for(int k=0;k<nsub;k++) {
         	index = listsub.AllocateNewElement();
         	listsub[index] = fElementVec[indexsubs[k]];
         }
         locel = SingularElement(point,listsub);
         hsub = h_Parameter(locel);
      }
      TPZInterpolatedElement *intel = (TPZInterpolatedElement *) locel;
      intel->SetSideOrder(2,Nc+1);
      indexlist.Push(intel->Index());
      //os elemento viz devem ter ordens menores a cel quanto mais longe de point
      TPZInterpolatedElement *neighkeep,*neigh;
      //feito só para o caso 1d , extender para o caso geral
      int dim = intel->Dimension();
      if(dim != 1) {
      	cout << "TPZCompMesh::Step3 not dimension implemented , dimension = " << cel->Dimension() << endl;
         exit(1);
      }
      for(int side=0;side<2;side++) {
         int ly = Nc;
         TPZGeoElSide neighside = intel->Reference()->Neighbour(side);
         TPZGeoElSide neighsidekeep = neighside;
         while(ly > 0 && neighsidekeep.Exists() && neighsidekeep.Element()->Id() > 0) {
            neigh = (TPZInterpolatedElement *) neighsidekeep.Element()->Reference();
            neigh->SetSideOrder(2,ly);
            int otherside = (neighsidekeep.Side()+1)%2;
            neighsidekeep.SetSide(otherside);
            indexlist.Push(neighsidekeep.Reference().Element()->Index());
            neighside = neighsidekeep.Neighbour();
            while(!neighside.Exists()) {
               TPZStack<TPZCompElSide> elvec(0);
               TPZCompElSide neighsidecomp = neighsidekeep.Reference();
               neighsidecomp.HigherLevelElementList(elvec,1,1);
               if(elvec.NElements()) {
               	neighside = elvec[0].Reference();
                  break;
               }
               neighside = neighsidecomp.LowerLevelElementList(1).Reference();
            }
            if(!neighside.Exists()) break;
            neighsidekeep = neighside;
            ly--;
         }
      }
}

void TPZCompMesh::Step4(
				void (*fExact)(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv),
         REAL w,REAL nadmerror,REAL tol,int maxiter,TPZStack<int> &indexlist) {

   REAL H1_error,L2_error,estimate;
   TPZVec<REAL> flux;
   int nlist = indexlist.NElements();

	for(int iel=0;iel<NElements();iel++) {
      int count = 0;
      //descartando os elemento próximos ao ponto singular processados com Step3(..)
      while(count < nlist) if(iel == indexlist[count++]) break;
      if(count == nlist) continue;
      fElementVec[iel]->EvaluateError(fExact,H1_error,L2_error,flux,estimate);
      REAL csi = H1_error / AdmLocalError(fExact,nadmerror);
      //calculo da ordem pn do elemento atual
      TPZInterpolatedElement *elem = (TPZInterpolatedElement *) fElementVec[iel];
      int p = elem->SideOrder(2);
      REAL pnj = p;
      REAL p_ln_csi = p + log(csi);
      REAL hdivp = h_Parameter(elem) / p;
      int iter=0;
      //algoritmo do ponto fixo
      REAL pnjmais1 = 0;
      while(iter < maxiter) {
	      REAL fipn_j = p_ln_csi - p*log(pnj * hdivp);
         pnjmais1 = (1.-w)*pnj + w*fipn_j;
         if(fabs(pnjmais1 - pnj) < tol) break;
         pnj = pnjmais1;
         iter++;
      }
      //escolha do inteiro mais proximo a pnjmais1
      int k=0;
      while(pnjmais1 - k > 1) k++;
      REAL pn;
      if(fabs(pnjmais1 - k) > .5) pn = k+1;
      else pn = k;
      elem->SetSideOrder(2,pn);
   }
}

REAL TPZCompMesh::MaximLocalError(void (*fExact)(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv)) {

   REAL H1_error,L2_error,estimate;
   TPZVec<REAL> flux(0);
   int nel = NElements();
   REAL maxerror;
	for(int iel=0;iel<nel;iel++) {
      fElementVec[iel]->EvaluateError(fExact,H1_error,L2_error,flux,estimate);
      if(maxerror < H1_error) maxerror = H1_error;
   }
	return maxerror;
}

void NullFunction(TPZVec<REAL> &point,TPZVec<REAL>&val,TPZFMatrix &deriv);
REAL TPZCompMesh::AdmLocalError(void (*fExact)(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv),REAL nadmerror) {

	REAL H1_e,H1_uhp,L2_error,estimate;
   int m = NElements();
	EvaluateError(fExact,H1_e,L2_error,estimate);
   EvaluateError(NullFunction,H1_uhp,L2_error,estimate);
   REAL admerror = nadmerror*sqrt(H1_uhp*H1_uhp + H1_e*H1_e) / sqrt(m);
	return admerror;
}

void NullFunction(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv) {

    val = 0*point[0];
    deriv(0,0) = 0.;
}
*/

REAL TPZCompMesh::CompareMesh(int var, char *matname){

	REAL error = 0.;
	int i=0;
	for (i=0;i<fElementVec.NElements();i++){
		TPZCompEl *el = fElementVec[i];
		if(el) error+= el->CompareElement(var,matname);
	}
	return (error);
}

void TPZCompMesh::SetElementSolution(int i, TPZVec<REAL> &sol) {
  if(sol.NElements() != NElements()) {
    cout << "TPZCompMesh::SetElementSolution size of the vector doesn't match\n";
  }
  if(fElementSolution.Cols() <= i) fElementSolution.Resize(NElements(),i+1);
  int el,nel= NElements();
  for(el=0; el<nel; el++) {
    fElementSolution(el,i) = sol[el];
  }
}

void TPZCompMesh::GetRefPatches(TPZStack<TPZGeoEl *> &grpatch){
	int i,j;
	bool test = 0;
	int nel = NElements();
	//	cout << "GetRefPatches\n" << nel << endl;
	for (i=0; i<nel; i++){
		if (fElementVec[i]){
			TPZGeoEl *gel = fElementVec[i]->Reference();
			if (gel)
			 //gel = fElementVec[i]->Reference();
			  //	cout << "Iniciando procura do elemento de referência do elemento " << fElementVec[i]->Index() << endl;
			gel = fElementVec[i]->GetRefElPatch();
		//	gel->Print();
			if (gel){
				test = false;
				int npat =grpatch.NElements();
				for (j=0; j<npat; j++){
				//	grpatch[j]->Print();
					if (grpatch[j] == gel) {
						test = true;
						break;
					}
				}
				if (!test){
					grpatch.Push(gel);
				}
			}
		}
	}
	return;
}


void  TPZCompMesh::GetNodeToElGraph(TPZVec<int> &nodtoelgraph, TPZVec<int> &nodtoelgraphindex, TPZStack<int> &elgraph, TPZVec<int> &elgraphindex){

  ComputeElGraph(elgraph,elgraphindex);
  
//  int i,j;
//  for (i=0; i<elgraphindex.NElements()-1; i++){
//    cout << "Block: " << i << "\t";
//    for (j = elgraphindex[i]; j<elgraphindex[i+1]; j++){
//      cout << elgraph[j] << "\t";
//    }
//   cout << endl;
//  }

  TPZMetis renum(elgraphindex.NElements() -1 ,NIndependentConnects());
  renum.SetElementGraph(elgraph,elgraphindex);
  renum.NodeToElGraph(elgraph,elgraphindex,nodtoelgraph, nodtoelgraphindex);	
/*   TPZRenumbering *re = &renum; */
/*   re->Print(nodtoelgraph, nodtoelgraphindex , "Grapho de nós para elementos "); */


}


void TPZCompMesh::GetElementPatch(TPZVec<int> nodtoelgraph, TPZVec<int> nodtoelgraphindex, TPZStack<int> &elgraph, TPZVec<int> &elgraphindex,int elind ,TPZStack<int> &patch){

  int aux =0;
  TPZAVLMap<int,int> elconmap(aux);
  int i,j;
  for (i= elgraphindex[elind]; i<elgraphindex[elind+1];i++){
    int node = elgraph[i];
    for (j = nodtoelgraphindex[node];  j<nodtoelgraphindex[node+1]; j++){
      elconmap[nodtoelgraph[j]] = elind; 
      
    }
  }
  patch.Resize(0);

  TPZPix iter = elconmap.First();
  while(iter){
    patch.Push(elconmap.Key(iter));
    elconmap.Next(iter);
  }	
}


TPZCompMesh* TPZCompMesh::Clone(){

  TPZGeoMesh *gmesh = Reference();
  gmesh->ResetReference();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  CopyMaterials(cmesh);

  int iel, nel = fElementVec.NElements();
  for (iel=0;iel<nel;iel++){
    if (!fElementVec[iel]) continue;
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (fElementVec[iel]);

    TPZGeoEl *gel = fElementVec[iel]->Reference();
    if (!gel) continue;
    int indexcomp;
    gel->CreateCompEl(*cmesh,indexcomp);
    
    TPZInterpolatedElement *clintel = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[indexcomp]);
    int side, nsides = intel->NConnects();
    for (side =0;side<nsides;side++){
      int porder = intel->SideOrder(side);
      clintel->PRefine(side,porder);
	  int clorder = clintel->SideOrder(side);
	  if(porder != clorder) {
		  cout << "PZCompMesh::Clone order was not honoured side " << 
			  side << " porder " << porder << " clorder " << clorder <<
			  endl;
	  }
      int elndof = intel->Connect(side).NDof(*this);
      int clelndof = clintel->Connect(side).NDof(*cmesh);
	  if(elndof != clelndof) {
		  cout << "PZCompMesh::Clone ndof different side " << 
			  side << " porder " << porder << " clorder " << clorder <<
			  " elndof = "<< elndof << " clelndof " << clelndof << endl;
	  }
    }
  }

  cmesh->InitializeBlock();

  gmesh->ResetReference();
  LoadReferences();

 
  int clnel = cmesh->ElementVec().NElements();
//  if (clnel != nel){
//    cout << "Clone mesh inconsistancy: clone elements failed!\n";
//    return 0;
//  }

  for (iel=0;iel<clnel;iel++){
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    TPZGeoEl *gel = cel->Reference();
    TPZCompEl *el = gel->Reference();
    int ic,nc = el->NConnects();
    for (ic=0;ic<nc;ic++){
      int elseqnum = el->Connect(ic).SequenceNumber();
      int elndof = el->Connect(ic).NDof(*this);
      int clelseqnum = cel->Connect(ic).SequenceNumber();
      int clelndof = cel->Connect(ic).NDof(*cmesh);
      if (elndof != clelndof){
		cout << "Number of degree of freedom incompatible between clone and original mesh!\n";
		cout << "Clone connect id: " << ic << "  Clone dof: " << clelndof << "  Original dof: " << clelndof << endl;
		continue;
      }
      int idof;
      for (idof=0; idof<elndof; idof++){
	cmesh->fSolutionBlock(clelseqnum,0,idof,0) = Block()(elseqnum,0,idof,0); 
      }
    }
  }
  return cmesh;
}


void TPZCompMesh::CopyMaterials(TPZCompMesh *mesh){
  int nmat = MaterialVec().NElements();
  int m;
  for(m=0; m<nmat; m++) {
    TPZMaterial *mat = MaterialVec()[m];
    if(!mat) continue;
    mat->Clone(mesh->MaterialVec());
  }

}

REAL TPZCompMesh::DeltaX(){

  int nel = ElementVec().NElements(),i,j;
  if(nel == 0) cout << "\nTPZCompMesh::DeltaX nenhum elemento computacional foi criado\n";
  REAL maxdist = 0.0,dist=0.0;
  TPZVec<REAL> point0(3),point1(3);
  TPZGeoNode *node0,*node1;
  for(i=0;i<nel;i++){
    TPZCompEl *com = ElementVec()[i];
    if(!com) continue;
    if(com->Type()==16 || com->Material()->Id() < 0) continue;
    node0 = com->Reference()->NodePtr(0);
    node1 = com->Reference()->NodePtr(1);
    for(j=0;j<3;j++){
      point0[j] = node0->Coord(j);
      point1[j] = node1->Coord(j);
    }
    dist = TPZGeoEl::Distance(point0,point1);
    if(dist > maxdist) maxdist = dist;
  }
  return maxdist;
}

REAL TPZCompMesh::MaximumRadiusOfMesh(){

  int nel = ElementVec().NElements(),i;
  if(nel == 0) cout << "\nTPZCompMesh::MaximumRadiusOfMesh malha vazia\n";
  REAL maxdist = 0.0,dist=0.0;
  TPZVec<REAL> point0(3),point1(3);
  for(i=0;i<nel;i++){
    TPZCompEl *com = ElementVec()[i];
    if(!com) continue;
    if(com->Type()==16 || com->Material()->Id() < 0) continue;
    dist = com->MaximumRadiusOfEl();
    if(dist > maxdist) maxdist = dist;
  }
  return maxdist;
}

REAL TPZCompMesh::LesserEdgeOfMesh(){

  int nel = ElementVec().NElements(),i;
  if(nel == 0) cout << "\nTPZCompMesh::MaximumRadiusOfMesh malha vazia\n";
  REAL mindist =10000.0,dist=0.0;
  for(i=0;i<nel;i++){
    TPZCompEl *com = ElementVec()[i];
    if(!com) continue;
    if(com->Type()==16 || com->Material()->Id() < 0) continue;
    dist = com->LesserEdgeOfEl();
    if(dist < mindist) mindist = dist;
  }
  return mindist;
}

