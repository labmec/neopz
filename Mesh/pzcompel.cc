//METHODS DEFINITION FOR CLASS ELBAS

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzcmesh.h"
#include "pzbndcond.h"
#include "pzelmat.h"
#include "pzconnect.h"
#include "pzblockdiag.h"

//#include "pzshapelinear.h"
//#include "pzshapequad.h"
//#include "pzshapetriang.h"
//#include "pzshapetetra.h"
//#include "pzshapepiram.h"
//#include "pzshapeprism.h"
//#include "pzshapecube.h"

#include "pzerror.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzsolve.h"

#include "pzmanvector.h"
#include "pzstack.h"

#include "pztransfer.h"
#include "pztrnsform.h"
#include "pzquad.h"

#include <math.h>
#include <stdlib.h>

void TPZCompEl::CalcBlockDiagonal(TPZStack<int> &connectlist, TPZBlockDiagonal & blockdiag) {

  TPZElementMatrix ek,ef;
  int b;
  CalcStiff(ek,ef);
  if(HasDependency()) {
    ApplyConstraints(ek,ef);
    int numblock = ek.fConstrConnect.NElements();
    TPZVec<int> blocksize(numblock);
    for(b=0; b<numblock; b++) blocksize[b] = ek.fConstrBlock->Size(b);
    blockdiag.Initialize(blocksize);
    connectlist = ek.fConstrConnect;
    for(b=0; b<numblock; b++) {
        int blsize = blocksize[b];
        int conind = ek.fConstrConnect[b];
        if(Mesh()->ConnectVec()[conind].HasDependency()) continue;
        TPZFMatrix ekbl(blsize,blsize);
        int r,c;
    	TPZBlock *mbl = ek.fConstrBlock;
        for(r=0; r<blsize; r++) {
            for(c=0; c<blsize; c++) {
                ekbl(r,c) = (*mbl)(b,b,r,c);
            }
        }
        blockdiag.AddBlock(b,ekbl);
    }
  } else {
    int numblock = ek.fConnect.NElements();
    TPZVec<int> blocksize(numblock);
    for(b=0; b<numblock; b++) blocksize[b] = ek.fBlock->Size(b);
    blockdiag.Initialize(blocksize);
    connectlist = ek.fConnect;
    for(b=0; b<numblock; b++) {
        int blsize = blocksize[b];
        TPZFMatrix ekbl(blsize,blsize);
        int r,c;
	    TPZBlock *mbl = ek.fBlock;
        for(r=0; r<blsize; r++) {
            for(c=0; c<blsize; c++) {
                ekbl(r,c) = (*mbl)(b,b,r,c);
            }
        }
        blockdiag.AddBlock(b,ekbl);
    }
  }
}

int TPZCompEl::gOrder = 2;

/*
void (*TPZCompEl::fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi) = TPZCompEl::Chebyshev;
void (*TPZShapeLinear::fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi) = TPZCompEl::fOrthogonal;
void (*TPZShapeQuad::fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi) = TPZCompEl::fOrthogonal;
void (*TPZShapeTriang::fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi) = TPZCompEl::fOrthogonal;
void (*TPZShapeTetra::fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi) = TPZCompEl::fOrthogonal;
void (*TPZShapePiram::fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi) = TPZCompEl::fOrthogonal;
void (*TPZShapePrism::fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi) = TPZCompEl::fOrthogonal;
void (*TPZShapeCube::fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi) = TPZCompEl::fOrthogonal;
*/

TPZCompEl::TPZCompEl() {
  fMesh = 0;
  fIndex = -1;
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, int &index) {
  fMesh = &mesh;
  index = mesh.ElementVec().AllocateNewElement();
  mesh.ElementVec()[index] = this;
  fIndex = index;
}

TPZCompEl::~TPZCompEl() {
  int index = Index();
	fMesh->ElementVec()[index] = 0;
	fMesh->ElementVec().SetFree(index);
}

MElementType TPZCompEl::Type() {
  cout << "TPZCompEl::Type unknown\n";
  return ENoType;
}

void TPZCompEl::LoadSolution() {
  // an element without mesh is a free element
  if(!Mesh() || !HasDependency()) return;
  TPZStack<int> connectlist;
  int totalconnects;
  BuildConnectList(connectlist);
  totalconnects = connectlist.NElements();
  TPZManVector<int> dependenceorder(totalconnects);
  BuildDependencyOrder(connectlist,dependenceorder);
  TPZMaterial *mat = Material();
  if(!mat) return;
  int numstate = mat->NStateVariables();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  int maxdep = 0;
  int in;
  int iv,jv,idf;
  REAL coef;
  for(in=0;in<totalconnects;in++)
    maxdep = (maxdep < dependenceorder[in]) ? dependenceorder[in] : maxdep;
  int current_order = maxdep-1;
  while(current_order >= 0) {
    for(in=0; in<totalconnects; in++) {
      if(dependenceorder[in] != current_order) continue;
      TPZConnect *dfn = &Mesh()->ConnectVec()[connectlist[in]];
      if(!dfn->HasDependency()) continue;
      int bl = dfn->SequenceNumber();
      int nvar = block.Size(bl);
      //         int numshape = nvar/numstate;
      TPZConnect::TPZDepend *dep = dfn->FirstDepend();
	  int blpos = block.Position(bl);
      for(iv=0; iv<nvar; iv++) MeshSol(blpos+iv, 0) = 0.;
      while(dep) {
	int depconindex = dep->fDepConnectIndex;
	TPZConnect &depcon = Mesh()->ConnectVec()[depconindex];
	int depseq = depcon.SequenceNumber();
	int numdepvar = block.Size(depseq);
	int depseqpos = block.Position(depseq);
	for(iv=0; iv<nvar; iv+=numstate) {
	  for(jv=0; jv<numdepvar; jv+=numstate) {
	    coef = dep->fDepMatrix(iv/numstate,jv/numstate);
	    for(idf=0; idf<numstate; idf++) MeshSol(blpos+iv+idf,0) += coef*MeshSol(depseqpos+jv+idf,0);
	  }
	}
	dep = dep->fNext;
      }
    }
    current_order--;
  }
}

void TPZCompEl::SetMesh(TPZCompMesh *mesh) {
  // The mesh can be null, indicating that the element is now unused
  //	if(!mesh)
  //		PZError << "TPZCompEl.SetMesh called with zero pointer\n";

  fMesh = mesh;
}

TPZCompMesh *TPZCompEl::Mesh() {
  //	if(!fMesh) PZError << "TPZCompEl.Mesh returned NULL pointer\n";
  return fMesh;
}

TPZConnect &TPZCompEl::Connect(int i) {
  if(fMesh) {
    int connectindex = ConnectIndex(i);
    if(connectindex >= 0) {
      return fMesh->ConnectVec()[connectindex];
    } else {
      PZError << "TPZCompEl::Connect called for noninitialized connect\n";
    }
  } else {
    PZError << "Connect called for an element without mesh\n";
    exit(-1);
  }
  static TPZConnect dummy;
  return dummy;
}


void TPZCompEl::Print(ostream & out) {
  out << "output for a computable element\n";
  out << "Number of connects = " << NConnects() << " Node indexes : ";
  int nod;
  for(nod=0; nod< NConnects(); nod++)
    out << ConnectIndex(nod) <<  ' ' ;
  out << endl;

}

ostream& operator<<(ostream &s,TPZCompEl & el){
  el.Print(s);
  return s;
}

void TPZCompEl::PrintSolution(TPZVec<REAL> &point,char *varname,ostream &s) {
  TPZGeoEl *georef = Reference();
  TPZMaterial *mat = Material();
  if(!georef || !mat) {
    cout << "TPZCompEl::PrintSolution should not be called for an element"
      " which doesnt have a geometric reference or material\n";
    Print();
    return;
  }
  int varindex = mat->VariableIndex(varname);
  int numvar = mat->NSolutionVariables(varindex);
  if(varindex == -1) return;
  TPZManVector<REAL> sol(numvar);
  sol.Fill(0.);
  Solution(point,varindex,sol);
  for(int i=0; i<sol.NElements(); i++) {
    s << sol[i] << '\t';
  }
}

void TPZCompEl::PrintCoordinate(TPZVec<REAL> &point,int CoordinateIndex,ostream &s) {
  TPZGeoEl *georef = Reference();
  if(!georef) {
    cout << "TPZCompEl::PrintCoordinate should not be called on an"
      " element which doesnt have a geometric reference\n";
    return;
  }
  TPZManVector<REAL> X(3);
  X.Fill(0.);
  georef->X(point,X);
  s << X[CoordinateIndex] << '\t';
}

void TPZCompEl::PrintTitle(char *varname,ostream &s) {
  TPZMaterial *mat = Material();
  TPZGeoEl *georef = Reference();
  if(!georef || !mat) {
    cout << "TPZCompEl::PrintTitle should not be called for an element"
      " which doesnt have a material\n";
    return;
  }
  int varindex = mat->VariableIndex(varname);
  if(varindex == -1) return;
  int numvar = mat->NSolutionVariables(varindex);
  if(numvar == 0) return;
  if(numvar == 1) {
    s << varname << '\t';
    return;
  }
  for(int i=0; i<numvar; i++) s << varname << '_' << i << '\t';
}

void TPZCompEl::Chebyshev(REAL x,int num,TPZFMatrix &phi,TPZFMatrix &dphi){
  // Quadratic or higher shape functions
  if(num <= 0) return;
  phi.Put(0,0,1.0);
  dphi.Put(0,0, 0.0);
  if(num == 1) return;
  phi.Put(1,0, x);
  dphi.Put(0,1, 1.0);
  int ord;
  for(ord = 2;ord<num;ord++) {
    phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
    dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
  }
}

inline void TPZCompEl::Divide(int index, TPZVec<int> &subindex, int interpolate) {
  subindex.Resize(0);
  cout<< "TPZCompEl::Divide called\n";
}

/**The TElementMatrix objects will be called ekmat and efmat respectively
vector containing all nodes of the constrained element will be kept into fConstrNod*/
void TPZCompEl::ApplyConstraints(TPZElementMatrix &ekmat,TPZElementMatrix &efmat) {
  int totalnodes= ekmat.fConnect.NElements();
  ekmat.fConstrConnect.Resize(totalnodes);
  if(totalnodes) ekmat.fConstrConnect.Fill(0,0);
  int in;
  for(in=0; in<totalnodes; in++) ekmat.fConstrConnect[in] = ekmat.fConnect[in];
  // total number of nodes of the constrained element
  BuildConnectList(ekmat.fConstrConnect);
  totalnodes = ekmat.fConstrConnect.NElements();

  // compute the list of nodes and their proper order of processing
  TPZVec<int> DependenceOrder(0);
  // ekmat.fConstrNod, totalnodes and DependenceOrder
  // are initialized using codes documented above
  BuildDependencyOrder(ekmat.fConstrConnect,DependenceOrder);

  // compute the number of statevariables
  // the number of state variables is the number of unknowns associated with
  // each shapefunction
  // numstate is best initialized during computation of the stiffness matrix
  TPZMaterial *mat = Material();
  int numstate = mat->NStateVariables();

  // initialize the block structure
  if(!ekmat.fConstrBlock) ekmat.fConstrBlock = new TPZBlock(0,totalnodes);
  else ekmat.fConstrBlock->SetNBlocks(totalnodes);

  // toteq contains the total number of equations of the constrained matrix
  int toteq = 0;
  for(in=0; in<totalnodes; in++) {
    int dfnindex = ekmat.fConstrConnect[in];
    TPZConnect &dfn = Mesh()->ConnectVec()[dfnindex];
    int ndf = dfn.NDof(*Mesh());
    ekmat.fConstrBlock->Set(in,ndf);
    toteq += ndf;
  }

  // initialize the constrained matrices
  if(!ekmat.fConstrMat) ekmat.fConstrMat = new TPZFMatrix;
  if(!efmat.fConstrMat) efmat.fConstrMat = new TPZFMatrix;

  ekmat.fConstrBlock->Resequence();
  ekmat.fConstrBlock->SetMatrix(ekmat.fConstrMat);

  int nrhs = efmat.fMat->Cols();
  ekmat.fConstrMat->Redim(toteq,toteq);
  efmat.fConstrMat->Redim(toteq,nrhs);

  // copy the original matrix to the constrained matrix
  int numnod = ekmat.fConnect.NElements();
  for(in=0; in<numnod; in++) {
    int irnode =0;
    int idfn = ekmat.fConnect[in];
    // find the index of the node in the destination (constrained) matrix
    while(irnode < totalnodes && ekmat.fConstrConnect[irnode] != idfn) irnode++;

    // first and last rows in the original matrix
    int ifirst = ekmat.fBlock->Position(in);
    int ilast = ifirst+ekmat.fBlock->Size(in);

    // first and last rows in the desination (reception) matrix
    int irfirst = ekmat.fConstrBlock->Position(irnode);
    //	   int irlast = irfirst+ekmat.fConstrBlock->Size(irnode);

    int i,ir,ieq;
    for(i=ifirst,ir=irfirst;i<ilast;i++,ir++) {
      for(ieq=0; ieq<nrhs; ieq++) {
	(*efmat.fConstrMat)(ir,ieq) = (*efmat.fMat)(i,ieq);
      }
    }
    int jn;
    for(jn=0; jn<numnod; jn++) {
      int jrnode = 0;
      int jdfn = ekmat.fConnect[jn];
      // find the index of the node in the destination (constrained) matrix
      while(jrnode < totalnodes && ekmat.fConstrConnect[jrnode] != jdfn) jrnode++;
      if(jrnode == totalnodes) {
	cout << "TPZCompEl::ApplyConstraints node not found in node list\n";
      }
      // first and last columns in the original matrix
      int jfirst = ekmat.fBlock->Position(jn);
      int jlast = jfirst+ekmat.fBlock->Size(jn);
      // first and last columns in the desination (reception) matrix
      int jrfirst = ekmat.fConstrBlock->Position(jrnode);
      //			int jrlast = irfirst+ekmat.fConstrBlock->Size(jrnode);
      int j,jr;
      for(i=ifirst,ir=irfirst;i<ilast; i++,ir++) {
	for(j=jfirst,jr=jrfirst;j<jlast; j++,jr++) {
	  (*ekmat.fConstrMat)(ir,jr) = (*ekmat.fMat)(i,j);
	}
      }
    }
  }

  int numnodes_processed = 0;
  int current_order = 0;
  while(numnodes_processed < totalnodes) {
    int in;
    for(in=0; in<totalnodes; in++) {
      int dfnindex = ekmat.fConstrConnect[in];
      TPZConnect *dfn = &Mesh()->ConnectVec()[dfnindex];
      if(DependenceOrder[in] != current_order) continue;

      // only nodes which have dependency order equal to the
      // current order are processed
      numnodes_processed++;

      int inpos = ekmat.fConstrBlock->Position(in);
      int insize = ekmat.fConstrBlock->Size(in);
      // inpos : position of the dependent equation
      // insize : number of equations processed

      // loop over the nodes from which dfn depends
      TPZConnect::TPZDepend *dep = dfn->FirstDepend();
      while(dep) {
	int depnodeindex = dep->fDepConnectIndex;
	// look for the index where depnode is found
	int depindex=0;
	while(depindex < totalnodes && ekmat.fConstrConnect[depindex] != depnodeindex) depindex++;
	if(depindex == totalnodes) {
	  cout << "TPZCompEl::ApplyConstraints node not found in node list\n";
	}

	int deppos = ekmat.fConstrBlock->Position(depindex);
	int depsize = ekmat.fConstrBlock->Size(depindex);
				// deppos : position of the receiving equation
				// depsize : number of receiving equations

	// process the rows of the constrained matrix
	int send;
	int receive;
	int ieq;
	REAL coef;
	int idf;

	for(send=inpos; send<inpos+insize; send += numstate) {
	  for(receive=deppos; receive<deppos+depsize; receive += numstate) {
	    coef = dep->fDepMatrix((send-inpos)/numstate,(receive-deppos)/numstate);
	    for(ieq=0; ieq<toteq; ieq++) for(idf=0; idf<numstate; idf++)  {
	      (*ekmat.fConstrMat)(receive+idf,ieq) += coef*(*ekmat.fConstrMat)(send+idf,ieq);
	    }
	    for(ieq=0; ieq<nrhs; ieq++) for(idf=0; idf<numstate; idf++) {
	      (*efmat.fConstrMat)(receive+idf,ieq) += coef*(*efmat.fConstrMat)(send+idf,ieq);
	    }
	  }
	}

	for(send=inpos; send<inpos+insize; send += numstate) {
	  for(receive=deppos; receive<deppos+depsize; receive += numstate) {
	    coef = dep->fDepMatrix((send-inpos)/numstate,(receive-deppos)/numstate);
	    for(ieq=0; ieq<toteq; ieq++) for(idf=0; idf<numstate; idf++) {
	      (*ekmat.fConstrMat)(ieq,receive+idf) += coef*(*ekmat.fConstrMat)(ieq,send+idf);
	    }
	  }
	}

	dep = dep->fNext;
      } // end of while
    } // end of loop over all nodes
    current_order++;
  } // end of while loop
}

/**The TElementMatrix objects will be called ekmat and efmat respectively
vector containing all nodes of the constrained element will be kept into fConstrNod*/
void TPZCompEl::ApplyConstraints(TPZElementMatrix &efmat) {
  int totalnodes= efmat.fConnect.NElements();
  efmat.fConstrConnect.Resize(totalnodes);
  if(totalnodes) efmat.fConstrConnect.Fill(0,0);
  int in;
  for(in=0; in<totalnodes; in++) efmat.fConstrConnect[in] = efmat.fConnect[in];
  // total number of nodes of the constrained element
  BuildConnectList(efmat.fConstrConnect);
  totalnodes = efmat.fConstrConnect.NElements();

  // compute the list of nodes and their proper order of processing
  TPZVec<int> DependenceOrder(0);
  // ekmat.fConstrNod, totalnodes and DependenceOrder
  // are initialized using codes documented above
  BuildDependencyOrder(efmat.fConstrConnect,DependenceOrder);

  // compute the number of statevariables
  // the number of state variables is the number of unknowns associated with
  // each shapefunction
  // numstate is best initialized during computation of the stiffness matrix
  TPZMaterial *mat = Material();
  int numstate = mat->NStateVariables();

  // initialize the block structure
  if(!efmat.fConstrBlock) efmat.fConstrBlock = new TPZBlock(0,totalnodes);
  else efmat.fConstrBlock->SetNBlocks(totalnodes);

  // toteq contains the total number of equations of the constrained matrix
  int toteq = 0;
  for(in=0; in<totalnodes; in++) {
    int dfnindex = efmat.fConstrConnect[in];
    TPZConnect &dfn = Mesh()->ConnectVec()[dfnindex];
    int ndf = dfn.NDof(*Mesh());
    efmat.fConstrBlock->Set(in,ndf);
    toteq += ndf;
  }

  // initialize the constrained matrices
  if(!efmat.fConstrMat) efmat.fConstrMat = new TPZFMatrix;

  efmat.fConstrBlock->Resequence();
  efmat.fConstrBlock->SetMatrix(efmat.fConstrMat);

  int nrhs = efmat.fMat->Cols();
  efmat.fConstrMat->Redim(toteq,nrhs);

  // copy the original matrix to the constrained matrix
  int numnod = efmat.fConnect.NElements();
  for(in=0; in<numnod; in++) {
    int irnode =0;
    int idfn = efmat.fConnect[in];
    // find the index of the node in the destination (constrained) matrix
    while(irnode < totalnodes && efmat.fConstrConnect[irnode] != idfn) irnode++;

    // first and last rows in the original matrix
    int ifirst = efmat.fBlock->Position(in);
    int ilast = ifirst+efmat.fBlock->Size(in);

    // first and last rows in the desination (reception) matrix
    int irfirst = efmat.fConstrBlock->Position(irnode);
    //	   int irlast = irfirst+ekmat.fConstrBlock->Size(irnode);

    int i,ir,ieq;
    for(i=ifirst,ir=irfirst;i<ilast;i++,ir++) {
      for(ieq=0; ieq<nrhs; ieq++) {
	(*efmat.fConstrMat)(ir,ieq) = (*efmat.fMat)(i,ieq);
      }
    }
    int jn;
    for(jn=0; jn<numnod; jn++) {
      int jrnode = 0;
      int jdfn = efmat.fConnect[jn];
      // find the index of the node in the destination (constrained) matrix
      while(jrnode < totalnodes && efmat.fConstrConnect[jrnode] != jdfn) jrnode++;
      if(jrnode == totalnodes) {
	cout << "TPZCompEl::ApplyConstraints node not found in node list\n";
      }
    }
  }

  int numnodes_processed = 0;
  int current_order = 0;
  while(numnodes_processed < totalnodes) {
    int in;
    for(in=0; in<totalnodes; in++) {
      int dfnindex = efmat.fConstrConnect[in];
      TPZConnect *dfn = &Mesh()->ConnectVec()[dfnindex];
      if(DependenceOrder[in] != current_order) continue;

      // only nodes which have dependency order equal to the
      // current order are processed
      numnodes_processed++;

      int inpos = efmat.fConstrBlock->Position(in);
      int insize = efmat.fConstrBlock->Size(in);
      // inpos : position of the dependent equation
      // insize : number of equations processed

      // loop over the nodes from which dfn depends
      TPZConnect::TPZDepend *dep = dfn->FirstDepend();
      while(dep) {
	int depnodeindex = dep->fDepConnectIndex;
	// look for the index where depnode is found
	int depindex=0;
	while(depindex < totalnodes && efmat.fConstrConnect[depindex] != depnodeindex) depindex++;
	if(depindex == totalnodes) {
	  cout << "TPZCompEl::ApplyConstraints node not found in node list\n";
	}

	int deppos = efmat.fConstrBlock->Position(depindex);
	int depsize = efmat.fConstrBlock->Size(depindex);
				// deppos : position of the receiving equation
				// depsize : number of receiving equations

	// process the rows of the constrained matrix
	int send;
	int receive;
	int ieq;
	REAL coef;
	int idf;

	for(send=inpos; send<inpos+insize; send += numstate) {
	  for(receive=deppos; receive<deppos+depsize; receive += numstate) {
	    coef = dep->fDepMatrix((send-inpos)/numstate,(receive-deppos)/numstate);
	    for(ieq=0; ieq<nrhs; ieq++) for(idf=0; idf<numstate; idf++) {
	      (*efmat.fConstrMat)(receive+idf,ieq) += coef*(*efmat.fConstrMat)(send+idf,ieq);
	    }
	  }
	}
	dep = dep->fNext;
      } // end of while
    } // end of loop over all nodes
    current_order++;
  } // end of while loop
}

void TPZCompEl::EvaluateError(void (* /*fp*/)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
			      REAL &/*true_error*/,REAL &/*0L2_error*/,TPZBlock * /*flux*/,REAL &/*estimate*/) {
  cout << "TPZCompEl::EvaluateError is called." << endl;
}

void TPZCompEl::Solution(TPZVec<REAL> &/*qsi*/,int var,TPZManVector<REAL> &sol){
  if(var >= 100) {
    int ind = Index();
    if(fMesh->ElementSolution().Cols() > var-100) {
      sol[0] = fMesh->ElementSolution()(ind,var-100);
    } else {
      sol[0] = 0.;
    }
  } else {
    sol.Resize(0);
  }
}


void TPZCompEl::BuildConnectList(TPZStack<int> &connectlist) {
	TPZConnect *dfn;
	int dfnindex;
	int nconnects = NConnects();
	int in;
	for(in=0; in<nconnects; in++) {
		dfnindex = ConnectIndex(in);
		dfn = &Connect(in);
		dfn->AddToList(dfnindex,*Mesh(),connectlist);
	}
}

int TPZCompEl::HasDependency() {
  int nconnects = NConnects();
  int in;
  for(in=0; in<nconnects; in++) if(Connect(in).HasDependency()) return 1;
  return 0;
}

void TPZCompEl::BuildDependencyOrder(TPZVec<int> &connectlist,
				     TPZVec<int> &DependenceOrder) {
  // nodelist (input) : vector which contains pointers to all nodes which
  // are in the dependency chain of the nodes of the element
  int totalnodes = connectlist.NElements();
  DependenceOrder.Resize(totalnodes);
  DependenceOrder.Fill(0,0);
  // initialize the vector which contains the
  // dependency order to zero
  int CurrentOrder = 0;
  // order which is currently processed
  int numnodes_processed = totalnodes;
  // number of nodes processed during the current cycle

  while(numnodes_processed) {

    numnodes_processed = 0;
    int i;
    for(i=0; i<totalnodes; i++) {
      int dfnindex = connectlist[i];
      TPZConnect *dfn = &Mesh()->ConnectVec()[dfnindex];
      if(dfn->HasDependency() && DependenceOrder[i] == CurrentOrder) {
	dfn->SetDependenceOrder(dfnindex,*Mesh(),CurrentOrder,connectlist,DependenceOrder);
	// this method will fill in the DependenceOrder vector by recursively
	// calling SetDependenceOrder for the nodes upon which dfn depends
	numnodes_processed++;
      }
    }
    CurrentOrder++;
  }
}

int TPZCompEl::Index() {

/*   int i=0; */
/*   int numel = fMesh->NElements(); */
/*   TPZAdmChunkVector<TPZCompEl *> &vec = fMesh->ElementVec(); */
/*   while ( i < numel ) { */
/*     if (vec[i] == this) break; */
/*     i++; */
/*   } */
/*   return i; */
  return fIndex;
}

int TPZCompEl::NEquations(){
	int i=0;
	int numeq=0;
	for (;i<NConnects(); i++){
		TPZConnect &df = Connect(i);
		if(df.HasDependency() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
		int seqnum = df.SequenceNumber();
		numeq += Mesh()->Block().Size(seqnum);
	}
	return numeq;
}

REAL TPZCompEl::CompareElement(int var, char *matname){
	cout << "TPZCompEl::CompareElement called!\n";
	return 0.;
}


// TPZCompElSide ///////////////////////////////////////////////////////////////
////////////// TPZCompElSide ///////////////////////////////////////////////////
////////////////////////// TPZCompElSide ///////////////////////////////////////
//TPZCompElSide::~TPZCompElSide() {}

TPZCompElSide::TPZCompElSide() {
  fEl = 0;
  fSide = -1;
  //  georeftest = 0;
}

TPZCompElSide::TPZCompElSide(const TPZCompElSide &celside) {
  fEl = celside.fEl;
  fSide   = celside.fSide;
  //  georeftest = celside.georeftest;
}

TPZCompElSide::TPZCompElSide(TPZCompEl *cel,int side) {
  fEl = cel;
  fSide   = side;
  //  if(cel) georeftest = cel->Reference();
  //  else georeftest = 0;
}

TPZGeoElSide TPZCompElSide::Reference() const {
  if(!fEl) return TPZGeoElSide();
  TPZGeoElSide sideel(fEl->Reference(),fSide);
  return sideel;
}

void TPZCompElSide::SetSide(int side) {
  fSide = side;
  Reference().SetSide(side);
}

void TPZCompElSide::HigherLevelElementList(TPZStack<TPZCompElSide> &elvec,
					   int onlyinterpolated, int removeduplicates) {
  TPZGeoElSide georef = Reference();
  if(!georef.Element()) return;
  georef.HigherLevelCompElementList2(elvec,onlyinterpolated,removeduplicates);
}

void TPZCompElSide::EqualLevelElementList(TPZStack<TPZCompElSide> &elsidevec,
					  int onlyinterpolated, int removeduplicates) {

  TPZGeoElSide georef = Reference();
  if(!georef.Exists()) return;
  georef.EqualLevelCompElementList(elsidevec,onlyinterpolated,removeduplicates);
}

void TPZCompElSide::HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates) {
  TPZGeoElSide georef = Reference();
  if(!georef.Exists()) return;
  georef.HigherDimensionElementList(elsidevec,onlyinterpolated);
  if(removeduplicates) RemoveDuplicates(elsidevec);

}


void TPZCompElSide::RemoveDuplicates(TPZStack<TPZCompElSide> &elvec) {

  RemoveConnectDuplicates(elvec);
  int i, nelem = elvec.NElements();
  for(i=0; i<nelem; i++) {
    TPZGeoElSide geli = elvec[i].Reference();
    if(!geli.Exists()) continue;
    int j;
    for(j=i+1; j<nelem; j++) {
      TPZGeoElSide gelj = elvec[j].Reference();
      if(!gelj.Exists()) continue;
      if(geli.NeighbourExists(gelj)) {
	int k;
	cout << "TPZCompElSide::RemoveDuplicates : case not identified by RemoveConnectDuplicates\n";
	for(k=j;k<nelem-1;k++) elvec[k] = elvec[k+1];
	elvec.Resize(nelem-1);
	nelem--;
	j--;
      }
    }
  }
}

void TPZCompElSide::ConnectedElementList(TPZStack<TPZCompElSide> &ellist,
					 int onlyinterpolated, int removeduplicates) {
  TPZCompElSide father = LowerLevelElementList(onlyinterpolated);
  if(father.Exists()) ellist.Push(father);
  EqualLevelElementList(ellist,onlyinterpolated,removeduplicates);
  HigherLevelElementList(ellist,onlyinterpolated,removeduplicates);
}

//retorna o elem. de menor nivel ao qual estou restrito
TPZCompElSide TPZCompElSide::LowerLevelElementList(int onlyinterpolated) {
  TPZGeoElSide georef = Reference();
  if(!georef.Exists()) return TPZCompElSide();

  return georef.LowerLevelCompElementList2(onlyinterpolated);
}

void TPZCompElSide::ExpandConnected(TPZStack<TPZCompElSide> &expandvec,int onlyinterpolated) {
  //insidevec contem todos os elem. de maior nivel que dependem de mim,
  //obtidos com HigherLevelElementList(..)
  TPZGeoElSide georef = Reference();
  if(!georef.Exists()) return;
  TPZStack<TPZCompElSide> highdimsidevec;
  TPZStack<int> smallsides;
  TPZCompElSide compside,compsidenext;
  int exnel = expandvec.NElements();
  for (int i=0;i<exnel;i++) {

    TPZGeoElSide ref = expandvec[i].Reference();
    smallsides.Resize(0);
    ref.Element()->LowerDimensionSides(ref.Side(),smallsides);
    for (int k=0;k<smallsides.NElements();k++) {
      compside = TPZCompElSide(Element(),smallsides[k]);
      compsidenext = compside.LowerLevelElementList(onlyinterpolated);
      //      TPZGeoElSide smallgeoside = compside.Reference();
      while(compsidenext.Exists()) {
	//         TPZGeoElSide geoside = compsidenext.Reference();//linha para teste
         if (compsidenext.Reference().NeighbourExists(this->Reference())) {
            TPZCompElSide lowid = LowerIdElementList(compside,onlyinterpolated);
            int l=0;
            while (l<expandvec.NElements() && lowid.Reference() != expandvec[l].Reference()) l++;
            if(l == expandvec.NElements()) {
                expandvec.Push(lowid);
            }
            break;
         }
         compside = compsidenext;
         compsidenext = compside.LowerLevelElementList(onlyinterpolated);
      }
    }
  }
}

TPZCompElSide TPZCompElSide::LowerIdElementList(TPZCompElSide &expandvec,int onlyinterpolated) {

  if(!expandvec.Element()) return TPZCompElSide();
  TPZGeoElSide gelside = expandvec.Reference();
  TPZStack<TPZGeoElSide> neighbourset;
  gelside.AllNeighbours(neighbourset);
  neighbourset.Push(gelside);
  int lowidindex = neighbourset.NElements()-1;
  //  TPZGeoElSide neighbour = gelside.Neighbour(),
  TPZGeoElSide lowidneigh(gelside);
//   if(!neighbour.Element()) return expandvec;
  int lowid = gelside.Id();
  int in, nneigh = neighbourset.NElements()-1;
  while(in < nneigh) {
    TPZCompEl *ref = neighbourset[in].Reference().Element();
    if(neighbourset[in].Id() < lowid && ref && (!onlyinterpolated || ref->IsInterpolated())) {
      lowidneigh = neighbourset[in];
      lowid = lowidneigh.Id();
      lowidindex = in;
    }
    in++;
  }
    
//   while(neighbour!=gelside) {     //pode ser um viz. inativo
//     TPZCompEl *ref = neighbour.Reference().Element();
//     if(neighbour.Id() < lowid && ref && (!onlyinterpolated || ref->IsInterpolated())) {
//       lowidneigh = neighbour;
//       lowid = lowidneigh.Id();
//     }
//     neighbour = neighbour.Neighbour();
//   }//for
  return lowidneigh.Reference();
}

void TPZCompElSide::RemoveConnectDuplicates(TPZStack<TPZCompElSide> &expandvec){
   int i,k;
	int nelems = expandvec.NElements();
   TPZStack<TPZCompElSide> locexpand;
   for(i=0;i<nelems;i++) locexpand.Push(expandvec[i]);
   //TPZCompElSide store[100];
   //TPZStack<TPZCompElSide> locexpand(store,nelems);
   //for(i=0;i<nelems;i++) locexpand[i] = expandvec[i];
   expandvec.Resize(0);
   for(k=0;k<nelems;k++){
      TPZCompEl *kel = locexpand[k].Element();
      if(!kel) continue;
      int kside = locexpand[k].Side();
      i=k+1;
      while(i<nelems){
         TPZCompEl *iel = locexpand[i].Element();
         int iside = locexpand[i].Side();
      	if(iel && kel->ConnectIndex(kside) == iel->ConnectIndex(iside))
            locexpand[i] = TPZCompElSide();
         i++;
      }
   }
   for(i=0;i<nelems;i++)
   	if(locexpand[i].Element()) expandvec.Push(locexpand[i]);
}


void TPZCompEl::LoadElementReference()
{

	TPZGeoEl *ref = Reference();
	if(ref) ref->SetReference(this);
}

void TPZCompEl::CalcResidual(TPZElementMatrix &ef){
	TPZElementMatrix ek;
	CalcStiff(ek,ef);
  cout << "TPZCompEl::CalcResidual(*) is called." << endl;
}

TPZGeoEl * TPZCompEl::GetRefElPatch(){
  // cout << "Obtendo elemento geométrico de referência " << Index() << endl;
//  Print();
  TPZGeoEl *ref = Reference();
  if (!ref) return (0);
//  ref->Print();

//  cout << "1";
  TPZStack <TPZGeoEl *> father;
  father.Push(ref);
  //		cout << ",2";
  //  int count=0;
  while(ref->Father()){
	  father.Push(ref->Father());
	  ref = ref->Father();
//	  ref->Print();
	}
  int j;
  while(father.NElements()) {
    TPZGeoEl *aux = father.Pop();
    for (j=0; j<aux->NSides(); j++){
    	int sidedimension = aux->SideDimension(j);
    	if(!sidedimension){
    	 continue;
    	}
    	
      TPZGeoElSide side(aux,j);
      TPZStack <TPZCompElSide> stack;
      side.EqualLevelCompElementList(stack,1,1);
      if(stack.NElements()){
	//      	cout << " \n \n \n ==================================\n ================================\nElemento PatchReference\n";
	//      	aux->Print();
      	return aux;
      }
    }
  }
  // 	cout << " \n \n \n ==================================\n ================================\nElemento PatchReference falho\n";
	// 	Reference()->Print();
  return (Reference());
}

/*
void TPZCompel::GetPatch(TPZStack <TPZCompel *> patch){
*/
  
REAL TPZCompEl::MaximumRadiusOfEl(){

  //O elemento deve ser a envoltura convexa dos seus vértices

  if(!this) cout << "TPZCompMesh::MaximumRadiusOfEl() null element";

  int i,j,k;
  REAL maxdist = 0.0,dist=0.0;
  TPZVec<REAL> point0(3),point1(3);
  TPZGeoNode *node0,*node1;
  int nvertices = Reference()->NNodes();
  for(i=0;i<nvertices;i++){
    for(j=i+1;j<nvertices;j++){
      node0 = Reference()->NodePtr(i);
      node1 = Reference()->NodePtr(j);
      for(k=0;k<3;k++){
	point0[k] = node0->Coord(k);
	point1[k] = node1->Coord(k);
      }
      dist = TPZGeoEl::Distance(point0,point1);
      if(dist > maxdist) maxdist = dist;
    }
  }
  return maxdist;
}

REAL TPZCompEl::LesserEdgeOfEl(){

  if(!this) cout << "TPZCompMesh::LesserEdgeOfEl null element";

  int i,j,k;
  REAL mindist = 1000.0,dist=0.0;
  TPZVec<REAL> point0(3),point1(3);
  TPZGeoNode *node0,*node1;
  int nvertices = Reference()->NNodes();
  for(i=0;i<nvertices;i++){
    for(j=i+1;j<nvertices;j++){
      node0 = Reference()->NodePtr(i);
      node1 = Reference()->NodePtr(j);
      for(k=0;k<3;k++){
	point0[k] = node0->Coord(k);
	point1[k] = node1->Coord(k);
      }
      dist = TPZGeoEl::Distance(point0,point1);
      if(dist < mindist) mindist = dist;
    }
  }
  return mindist;
}
