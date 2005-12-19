//$Id: pzcompel.cpp,v 1.21 2005-12-19 11:56:57 tiago Exp $

//METHODS DEFINITION FOR CLASS ELBAS

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzcmesh.h"
#include "pzbndcond.h"
#include "pzelmat.h"
#include "pzconnect.h"
#include "pzblockdiag.h"

#include "pzshapelinear.h"
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

#include "TPZInterfaceEl.h"
#include "TPZCompElDisc.h"
#include "pzintel.h"

#include <math.h>
#include <stdlib.h>

#include <sstream>
using namespace std;

// LOGPZ
//#include <LOGPZ/logger.h>
// #include <LOGPZ/basicconfigurator.h>
// #include <LOGPZ/propertyconfigurator.h>
// #include <LOGPZ/helpers/exception.h>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompel"));
static LoggerPtr loggerSide(Logger::getLogger("pz.mesh.tpzcompelside"));
#endif

void TPZCompEl::CalcBlockDiagonal(TPZStack<int> &connectlist, TPZBlockDiagonal & blockdiag) {
  LOGPZ_INFO(logger, "Entering CalcBlockDiagonal.");
  TPZElementMatrix ek,ef;
  int b;
  CalcStiff(ek,ef);
  if(HasDependency()) {
    ApplyConstraints(ek,ef);
    int numblock = ek.fConstrConnect.NElements();
    TPZVec<int> blocksize(numblock);
    for(b=0; b<numblock; b++) blocksize[b] = ek.fConstrBlock.Size(b);
    blockdiag.Initialize(blocksize);
    connectlist = ek.fConstrConnect;
    for(b=0; b<numblock; b++) {
      int blsize = blocksize[b];
      int conind = ek.fConstrConnect[b];
      if(Mesh()->ConnectVec()[conind].HasDependency()) continue;
      TPZFMatrix ekbl(blsize,blsize);
      int r,c;
      TPZBlock &mbl = ek.fConstrBlock;
      for(r=0; r<blsize; r++) {
        for(c=0; c<blsize; c++) {
          ekbl(r,c) = (mbl)(b,b,r,c);
        }
      }
      blockdiag.AddBlock(b,ekbl);
    }
  } else {
    int numblock = ek.fConnect.NElements();
    TPZVec<int> blocksize(numblock);
    for(b=0; b<numblock; b++) blocksize[b] = ek.fBlock.Size(b);
    blockdiag.Initialize(blocksize);
    connectlist = ek.fConnect;
    for(b=0; b<numblock; b++) {
      int blsize = blocksize[b];
      TPZFMatrix ekbl(blsize,blsize);
      int r,c;
      TPZBlock &mbl = ek.fBlock;
      for(r=0; r<blsize; r++) {
        for(c=0; c<blsize; c++) {
          ekbl(r,c) = (mbl)(b,b,r,c);
        }
      }
      blockdiag.AddBlock(b,ekbl);
    }
  }
  LOGPZ_INFO(logger, "Exiting CalcBlockDiagonal.");
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
  LOGPZ_INFO(logger, "Entering empty constructor.");
  fMesh = 0;
  fIndex = -1;
  LOGPZ_INFO(logger, "Exiting empty constructor.");
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) {
  LOGPZ_INFO(logger, "Entering constructor (TPZCompMesh &mesh, TPZGeoEl *ref, int &index).");
  fMesh = &mesh;
  index = mesh.ElementVec().AllocateNewElement();
  mesh.ElementVec()[index] = this;
  fIndex = index;
  fReference = ref;
  fReferenceIndex = (ref == 0) ? -1 : ref->Index();
  LOGPZ_INFO(logger, "Exiting empty constructor (TPZCompMesh &mesh, TPZGeoEl *ref, int &index).");
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy) {
  LOGPZ_INFO(logger, "Entering constructor (TPZCompMesh &mesh, const TPZCompEl &copy).");
  fMesh = &mesh;
  int index = copy.fIndex;
  if(index >= 0) mesh.ElementVec()[index] = this;
  fIndex = index;
  fReference = copy.Reference();
  fReferenceIndex = copy.fReferenceIndex;
  LOGPZ_INFO(logger, "Exiting empty constructor (TPZCompMesh &mesh, const TPZCompEl &copy).");
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, int &index) {
  LOGPZ_INFO(logger, "Entering constructor (TPZCompMesh &mesh, const TPZCompEl &copy, int &index).");
  fMesh = &mesh;
  index = mesh.ElementVec().AllocateNewElement();
  if(index >= 0) mesh.ElementVec()[index] = this;
  fIndex = index;
  fReference = copy.fReference;
  fReferenceIndex = copy.fReferenceIndex;
  LOGPZ_INFO(logger, "Exiting empty constructor (TPZCompMesh &mesh, const TPZCompEl &copy, int &index).");
}

TPZCompEl::~TPZCompEl() {
  LOGPZ_INFO(logger, "Entering destructor.");
  int index = Index();
  fMesh->ElementVec()[index] = 0;
  fMesh->ElementVec().SetFree(index);
  LOGPZ_INFO(logger, "Exiting destructor.");
}

MElementType TPZCompEl::Type() {
  LOGPZ_INFO(logger, "Entering Type.");
  LOGPZ_WARN(logger, "Type unknown");
  LOGPZ_INFO(logger, "Exiting Type.");
  return ENoType;
}

void TPZCompEl::LoadSolution() {
  LOGPZ_INFO(logger, "Entering LoadSolution.");
  // an element without mesh is a free element
  if(!Mesh() || !HasDependency()) return;
  TPZStack<int> connectlist;
  int totalconnects;
  BuildConnectList(connectlist);
  totalconnects = connectlist.NElements();
  TPZManVector<int> dependenceorder(totalconnects);
  BuildDependencyOrder(connectlist,dependenceorder);
  TPZMaterial *mat = Material();
  if(!mat) {
    LOGPZ_WARN(logger, "Exiting LoadSolution because a null material was reached.");
    return;
  }
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
  LOGPZ_INFO(logger, "Exiting normally from LoadSolution.");
}

void TPZCompEl::SetMesh(TPZCompMesh *mesh) {
  LOGPZ_INFO(logger, "Entering SetMesh.");
  // The mesh can be null, indicating that the element is now unused
  if(!mesh)
    LOGPZ_ERROR(logger, "TPZCompEl.SetMesh called with zero pointer.");
  fMesh = mesh;
  LOGPZ_INFO(logger, "Exiting SetMesh.");
}

TPZCompMesh *TPZCompEl::Mesh() {
  LOGPZ_INFO(logger, "Entering Mesh.");
  if(!fMesh) 
    LOGPZ_ERROR(logger, "TPZCompEl.Mesh called for a uninitialized element.");
  LOGPZ_INFO(logger, "Exiting Mesh.");
  return fMesh;
}

TPZConnect &TPZCompEl::Connect(int i) {
  LOGPZ_INFO(logger, "Entering Connect.");
#ifndef NODEBUG
  if(fMesh) {
    int connectindex = ConnectIndex(i);
    if(connectindex >= 0) {
      return fMesh->ConnectVec()[connectindex];
    } else {
      LOGPZ_ERROR(logger, "TPZCompEl::Connect called for noninitialized connect\n");
    }
  } else {
    LOGPZ_FATAL(logger, "Connect called for an element without mesh\n");
    exit(-1);
  }
  static TPZConnect dummy;
  LOGPZ_INFO(logger, "Exiting Connect.")
  return dummy;
#else
  LOGPZ_INFO(logger, "Exiting Connect.")
  return fMesh->ConnectVec()[ConnectIndex(i)];
#endif
}

void TPZCompEl::Print(std::ostream & out) {
  LOGPZ_INFO(logger, "Entering Print.");
  out << "output for a computable element\n";
  out << "Number of connects = " << NConnects() << " Node indexes : ";
  int nod;
  for(nod=0; nod< NConnects(); nod++)
    out << ConnectIndex(nod) <<  ' ' ;
  out << endl;
  LOGPZ_INFO(logger, "Exiting Print.");
}

std::ostream& operator<<(std::ostream &s,TPZCompEl & el){
  LOGPZ_INFO(logger, "Entering operator <<.");
  el.Print(s);
  LOGPZ_INFO(logger, "Exiting operator <<.");
  return s;
}

void TPZCompEl::PrintSolution(TPZVec<REAL> &point,char *varname,std::ostream &s) {
  LOGPZ_INFO(logger, "Entering PrintSolution.");
  TPZGeoEl *georef = Reference();
  TPZMaterial *mat = Material();
  if(!georef || !mat) {
    LOGPZ_WARN(logger, "Exiting PrintSolution should not be called for an element which doesnt have a geometric reference or material");
    Print();
    return;
  }
  int varindex = mat->VariableIndex(varname);
  int numvar = mat->NSolutionVariables(varindex);
  if(varindex == -1) {
    LOGPZ_WARN(logger, "Exiting PrintSolution should not be called for an element which has unknown variable index");
    return;
  }
  TPZManVector<REAL> sol(numvar);
  sol.Fill(0.);
  Solution(point,varindex,sol);
  for(int i=0; i<sol.NElements(); i++) {
    s << sol[i] << '\t';
  }
  LOGPZ_INFO(logger, "Exiting PrintSolution normally.");
}

void TPZCompEl::PrintCoordinate(TPZVec<REAL> &point,int CoordinateIndex,std::ostream &s) {
  LOGPZ_INFO(logger, "Entering PrintCoordinate.");
  TPZGeoEl *georef = Reference();
  if(!georef) {
    LOGPZ_WARN(logger,"Exiting PrintCoordinate should not be called on an element which doesnt have a geometric reference");
    return;
  }
  TPZManVector<REAL> X(3);
  X.Fill(0.);
  georef->X(point,X);
  s << X[CoordinateIndex] << '\t';
  LOGPZ_INFO(logger, "Exiting PrintCoordinate normally.");
}

void TPZCompEl::PrintTitle(char *varname,std::ostream &s) {
  LOGPZ_INFO(logger, "Entering PrintTitle.");
  TPZMaterial *mat = Material();
  TPZGeoEl *georef = Reference();
  if(!georef || !mat) {
    LOGPZ_WARN(logger,"Exiting PrintTitle should not be called for an element which doesnt have a material");
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
  LOGPZ_INFO(logger, "Exiting PrintTitle.");
}

void TPZCompEl::Chebyshev(REAL x,int num,TPZFMatrix &phi,TPZFMatrix &dphi){
  LOGPZ_INFO(logger, "Entering Chebyshev");
  // Quadratic or higher shape functions
  if(num <= 0) {
    LOGPZ_INFO(logger, "Exiting Chebyshev - num <= 0");
    return;
  }
  phi.Put(0,0,1.0);
  dphi.Put(0,0, 0.0);
  if(num == 1) {
    LOGPZ_INFO(logger, "Exiting Chebyshev num == 1");
    return;
  }
  phi.Put(1,0, x);
  dphi.Put(0,1, 1.0);
  int ord;
  for(ord = 2;ord<num;ord++) {
    phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
    dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
  }
  LOGPZ_INFO(logger, "Exiting Chebyshev");
}

inline void TPZCompEl::Divide(int index, TPZVec<int> &subindex, int interpolate) {
  LOGPZ_INFO(logger, "Entering Divide");
  subindex.Resize(0);
  LOGPZ_WARN(logger,"TPZCompEl::Divide called");
  LOGPZ_INFO(logger, "Exiting Divide");
}

/**The TElementMatrix objects will be called ekmat and efmat respectively
vector containing all nodes of the constrained element will be kept into fConstrNod*/
void TPZCompEl::ApplyConstraints(TPZElementMatrix &ekmat,TPZElementMatrix &efmat) {
  LOGPZ_INFO(logger, "Entering ApplyConstraints - ek ef");
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
  ekmat.fConstrBlock.SetNBlocks(totalnodes);

  // toteq contains the total number of equations of the constrained matrix
  int toteq = 0;
  for(in=0; in<totalnodes; in++) {
    int dfnindex = ekmat.fConstrConnect[in];
    TPZConnect &dfn = Mesh()->ConnectVec()[dfnindex];
    int ndf = dfn.NDof(*Mesh());
    ekmat.fConstrBlock.Set(in,ndf);
    toteq += ndf;
  }

  ekmat.fConstrBlock.Resequence();
  ekmat.fConstrBlock.SetMatrix(&ekmat.fConstrMat);

  int nrhs = efmat.fMat.Cols();
  ekmat.fConstrMat.Redim(toteq,toteq);
  efmat.fConstrMat.Redim(toteq,nrhs);

  // copy the original matrix to the constrained matrix
  int numnod = ekmat.fConnect.NElements();
  for(in=0; in<numnod; in++) {
    int irnode =0;
    int idfn = ekmat.fConnect[in];
    // find the index of the node in the destination (constrained) matrix
    while(irnode < totalnodes && ekmat.fConstrConnect[irnode] != idfn) irnode++;

    // first and last rows in the original matrix
    int ifirst = ekmat.fBlock.Position(in);
    int ilast = ifirst+ekmat.fBlock.Size(in);

    // first and last rows in the desination (reception) matrix
    int irfirst = ekmat.fConstrBlock.Position(irnode);
    //	   int irlast = irfirst+ekmat.fConstrBlock->Size(irnode);

    int i,ir,ieq;
    for(i=ifirst,ir=irfirst;i<ilast;i++,ir++) {
      for(ieq=0; ieq<nrhs; ieq++) {
        (efmat.fConstrMat)(ir,ieq) = (efmat.fMat)(i,ieq);
      }
    }
    int jn;
    for(jn=0; jn<numnod; jn++) {
      int jrnode = 0;
      int jdfn = ekmat.fConnect[jn];
      // find the index of the node in the destination (constrained) matrix
      while(jrnode < totalnodes && ekmat.fConstrConnect[jrnode] != jdfn) jrnode++;
      if(jrnode == totalnodes) {
        LOGPZ_WARN(logger, "node not found in node list");
      }
      // first and last columns in the original matrix
      int jfirst = ekmat.fBlock.Position(jn);
      int jlast = jfirst+ekmat.fBlock.Size(jn);
      // first and last columns in the desination (reception) matrix
      int jrfirst = ekmat.fConstrBlock.Position(jrnode);
      //int jrlast = irfirst+ekmat.fConstrBlock->Size(jrnode);
      int j,jr;
      for(i=ifirst,ir=irfirst;i<ilast; i++,ir++) {
        for(j=jfirst,jr=jrfirst;j<jlast; j++,jr++) {
          (ekmat.fConstrMat)(ir,jr) = (ekmat.fMat)(i,j);
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

      int inpos = ekmat.fConstrBlock.Position(in);
      int insize = ekmat.fConstrBlock.Size(in);
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
          LOGPZ_WARN(logger,"node not found in node list");
        }

        int deppos = ekmat.fConstrBlock.Position(depindex);
        int depsize = ekmat.fConstrBlock.Size(depindex);
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
              (ekmat.fConstrMat)(receive+idf,ieq) += coef*(ekmat.fConstrMat)(send+idf,ieq);
            }
            for(ieq=0; ieq<nrhs; ieq++) for(idf=0; idf<numstate; idf++) {
              (efmat.fConstrMat)(receive+idf,ieq) += coef*(efmat.fConstrMat)(send+idf,ieq);
            }
          }
        }

        for(send=inpos; send<inpos+insize; send += numstate) {
          for(receive=deppos; receive<deppos+depsize; receive += numstate) {
            coef = dep->fDepMatrix((send-inpos)/numstate,(receive-deppos)/numstate);
            for(ieq=0; ieq<toteq; ieq++) for(idf=0; idf<numstate; idf++) {
              (ekmat.fConstrMat)(ieq,receive+idf) += coef*(ekmat.fConstrMat)(ieq,send+idf);
            }
          }
        }

        dep = dep->fNext;
      } // end of while
    } // end of loop over all nodes
    current_order++;
  } // end of while loop
  LOGPZ_INFO(logger, "Exiting ApplyConstraints - ek ef");
}

/**The TElementMatrix objects will be called ekmat and efmat respectively
vector containing all nodes of the constrained element will be kept into fConstrNod*/
void TPZCompEl::ApplyConstraints(TPZElementMatrix &efmat) {
  LOGPZ_INFO(logger, "Entering ApplyConstraints - ef");
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
  efmat.fConstrBlock.SetNBlocks(totalnodes);

  // toteq contains the total number of equations of the constrained matrix
  int toteq = 0;
  for(in=0; in<totalnodes; in++) {
    int dfnindex = efmat.fConstrConnect[in];
    TPZConnect &dfn = Mesh()->ConnectVec()[dfnindex];
    int ndf = dfn.NDof(*Mesh());
    efmat.fConstrBlock.Set(in,ndf);
    toteq += ndf;
  }

  // initialize the constrained matrices

  efmat.fConstrBlock.Resequence();
  efmat.fConstrBlock.SetMatrix(&efmat.fConstrMat);

  int nrhs = efmat.fMat.Cols();
  efmat.fConstrMat.Redim(toteq,nrhs);

  // copy the original matrix to the constrained matrix
  int numnod = efmat.fConnect.NElements();
  for(in=0; in<numnod; in++) {
    int irnode =0;
    int idfn = efmat.fConnect[in];
    // find the index of the node in the destination (constrained) matrix
    while(irnode < totalnodes && efmat.fConstrConnect[irnode] != idfn) irnode++;

    // first and last rows in the original matrix
    int ifirst = efmat.fBlock.Position(in);
    int ilast = ifirst+efmat.fBlock.Size(in);

    // first and last rows in the desination (reception) matrix
    int irfirst = efmat.fConstrBlock.Position(irnode);
    //	   int irlast = irfirst+ekmat.fConstrBlock->Size(irnode);

    int i,ir,ieq;
    for(i=ifirst,ir=irfirst;i<ilast;i++,ir++) {
      for(ieq=0; ieq<nrhs; ieq++) {
        (efmat.fConstrMat)(ir,ieq) = (efmat.fMat)(i,ieq);
      }
    }
    int jn;
    for(jn=0; jn<numnod; jn++) {
      int jrnode = 0;
      int jdfn = efmat.fConnect[jn];
      // find the index of the node in the destination (constrained) matrix
      while(jrnode < totalnodes && efmat.fConstrConnect[jrnode] != jdfn) jrnode++;
      if(jrnode == totalnodes) {
        LOGPZ_WARN(logger, "node not found in node list");
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

      int inpos = efmat.fConstrBlock.Position(in);
      int insize = efmat.fConstrBlock.Size(in);
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
          LOGPZ_WARN(logger, "node not found in node list");
        }
        int deppos = efmat.fConstrBlock.Position(depindex);
        int depsize = efmat.fConstrBlock.Size(depindex);
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
              (efmat.fConstrMat)(receive+idf,ieq) += coef*(efmat.fConstrMat)(send+idf,ieq);
            }
          }
        }
        dep = dep->fNext;
      } // end of while
    } // end of loop over all nodes
    current_order++;
  } // end of while loop
  LOGPZ_INFO(logger, "Exiting ApplyConstraints - ef");
}

void TPZCompEl::EvaluateError(void (* /*fp*/)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
                              TPZVec<REAL> &/*errors*/,TPZBlock * /*flux*/) {
  LOGPZ_INFO(logger, "Entering EvaluateError");
  LOGPZ_WARN(logger, "EvaluateError is called.");
  LOGPZ_INFO(logger, "Exiting EvaluateError");
}

void TPZCompEl::Solution(TPZVec<REAL> &/*qsi*/,int var,TPZManVector<REAL> &sol){
  LOGPZ_INFO(logger, "Entering Solution");
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
  LOGPZ_INFO(logger, "Exiting Solution");
}

void TPZCompEl::BuildConnectList(TPZStack<int> &connectlist) {
  LOGPZ_INFO(logger, "Entering BuildConnectList");
  TPZConnect *dfn;
  int dfnindex;
  int nconnects = NConnects();
  int in;
  for(in=0; in<nconnects; in++) {
    dfnindex = ConnectIndex(in);
    dfn = &Connect(in);
    dfn->AddToList(dfnindex,*Mesh(),connectlist);
  }
  LOGPZ_INFO(logger, "Exiting BuildConnectList");  
}

int TPZCompEl::HasDependency() {
  LOGPZ_INFO(logger, "Entering HasDependency");
  int nconnects = NConnects();
  int in;
  for(in=0; in<nconnects; in++) if(Connect(in).HasDependency()){
    LOGPZ_INFO(logger, "Exiting True HasDependency");
    return 1;
  }
  LOGPZ_INFO(logger, "Exiting False HasDependency");
  return 0;
}

void TPZCompEl::BuildDependencyOrder(TPZVec<int> &connectlist,
                                     TPZVec<int> &DependenceOrder) {
  LOGPZ_INFO(logger, "Entering BuildDependencyOrder");
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
  LOGPZ_INFO(logger, "Exiting BuildDependencyOrder");
}

int TPZCompEl::Index() {
  LOGPZ_INFO(logger, "Entering Index");
/*   int i=0; */
/*   int numel = fMesh->NElements(); */
/*   TPZAdmChunkVector<TPZCompEl *> &vec = fMesh->ElementVec(); */
/*   while ( i < numel ) { */
/*     if (vec[i] == this) break; */
/*     i++; */
/*   } */
/*   return i; */
  LOGPZ_INFO(logger, "Exiting Index");
  return fIndex;
}

int TPZCompEl::NEquations(){
  LOGPZ_INFO(logger, "Entering NEquations");
  int i=0;
  int numeq=0;
  for (;i<NConnects(); i++){
    TPZConnect &df = Connect(i);
    if(df.HasDependency() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
    int seqnum = df.SequenceNumber();
    numeq += Mesh()->Block().Size(seqnum);
  }
  LOGPZ_INFO(logger, "Exiting NEquations");
  return numeq;
}

REAL TPZCompEl::CompareElement(int var, char *matname){
	cout << "TPZCompEl::CompareElement called!\n";
	LOGPZ_INFO(logger, "Entering CompareElement");
	LOGPZ_WARN(logger, "CompareElement called!");
	LOGPZ_INFO(logger, "Exiting CompareElement");
	return 0.;
}


// TPZCompElSide ///////////////////////////////////////////////////////////////
////////////// TPZCompElSide ///////////////////////////////////////////////////
////////////////////////// TPZCompElSide ///////////////////////////////////////
//TPZCompElSide::~TPZCompElSide() {}

void TPZCompEl::LoadElementReference()
{
  LOGPZ_INFO(logger, "Entering LoadElementReference");
  TPZGeoEl *ref = Reference();
  if(ref) ref->SetReference(this);
  LOGPZ_INFO(logger, "Exiting LoadElementReference");
}

void TPZCompEl::CalcResidual(TPZElementMatrix &ef){
  LOGPZ_INFO(logger, "Entering CalcResidual");
  TPZElementMatrix ek;
  CalcStiff(ek,ef);
//  cout << "TPZCompEl::CalcResidual(*) is called." << endl;
  LOGPZ_INFO(logger, "Exiting CalcResidual");
}

TPZGeoEl * TPZCompEl::GetRefElPatch(){
  LOGPZ_INFO(logger, "Entering GetRefElPatch");
  std::stringstream sout;
  sout << "Obtendo elemento geom�rico de refer�cia " << Index() << endl;
  Print(sout);
  TPZGeoEl *ref = Reference();
  if (!ref) {
    LOGPZ_WARN(logger, "reached a null reference");
    return (0);
  }
  ref->Print(sout);
  TPZStack <TPZGeoEl *> father;
  father.Push(ref);
  while(ref->Father()){
    father.Push(ref->Father());
    ref = ref->Father();
    ref->Print(sout);
  }
  int j;
  LOGPZ_DEBUG(logger, sout.str());
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
        //cout << " \n \n \n ==================================\n ================================\nElemento PatchReference\n";
        //aux->Print();
        LOGPZ_INFO(logger, "Exing GetRefElPatch");
        return aux;
      }
    }
  }
  //cout << " \n \n \n ==================================\n ================================\nElemento PatchReference falho\n";
  //Reference()->Print();
  LOGPZ_WARN(logger, "Exing GetRefElPatch - Elemento PatchReference falho");
  return (Reference());
}

REAL TPZCompEl::MaximumRadiusOfEl(){
  LOGPZ_INFO(logger, "Entering MaximumRadiusOfEl");
  //O elemento deve ser a envoltura convexa dos seus v�tices
  if(!this) {
    LOGPZ_ERROR(logger,"TPZCompMesh::MaximumRadiusOfEl() null element");
  }

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
  LOGPZ_INFO(logger, "Exiting MaximumRadiusOfEl");
  return maxdist;
}

REAL TPZCompEl::LesserEdgeOfEl(){
  LOGPZ_INFO(logger, "Entering LesserEdgeOfEl");
  
  if(!this) LOGPZ_INFO(logger,"TPZCompMesh::LesserEdgeOfEl null element");

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
  LOGPZ_INFO(logger, "Exiting LesserEdgeOfEl");
  return mindist;
}


  /**
  Save the element data to a stream
  */
void TPZCompEl::Write(TPZStream &buf, int withclassid) 
{
  LOGPZ_INFO(logger, "Entering Write");
  TPZSaveable::Write(buf,withclassid);
  buf.Write(&fIndex,1);
  buf.Write(&fReferenceIndex,1);
  LOGPZ_INFO(logger, "Exiting Write");
}  
  
  /**
  Read the element data from a stream
  */
void TPZCompEl::Read(TPZStream &buf, void *context)
{
  LOGPZ_INFO(logger, "Entering Read");
  TPZSaveable::Read(buf,context);
  fMesh = (TPZCompMesh *) context;
  buf.Read(&fIndex,1);
  buf.Read(&fReferenceIndex,1);
  LOGPZ_INFO(logger, "Exiting Read");
}

void TPZCompEl::SetOrthogonalFunction(void (*orthogonal)(REAL x,int num,
                                      TPZFMatrix & phi,TPZFMatrix & dphi)) {
  LOGPZ_INFO(logger, "Entering SetOrthogonalFunction");
  pzshape::TPZShapeLinear::fOrthogonal = orthogonal;
  LOGPZ_INFO(logger, "Exiting SetOrthogonalFunction");
}

TPZInterfaceElement * TPZCompEl::CreateInterface(int side)
{
  LOGPZ_INFO(logger, "Entering CreateInterface");
  TPZInterfaceElement * newcreatedinterface = NULL;

  TPZGeoEl *ref = Reference();
  if(!ref) {
    LOGPZ_WARN(logger, "Exiting CreateInterface Null reference reached - NULL interface used");
    return newcreatedinterface;
  }

  TPZCompElSide thisside(this,side);
  TPZStack<TPZCompElSide> list;
  list.Resize(0);
  thisside.EqualLevelElementList(list,0,1);//retorna distinto ao atual ou nulo
  int size = list.NElements();
  //espera-se ter os elementos computacionais esquerdo e direito 
  //j�criados antes de criar o elemento interface
  if(size){
    //Interface has the same material of the neighbour with lesser dimension.
    //It makes the interface has the same material of boundary conditions (TPZCompElDisc with interface dimension)
    LOGPZ_DEBUG(logger, "There are size");
    int matid;
    int thisdim = this->Dimension();
    int neighbourdim = list[0].Element()->Dimension();
    if (thisdim == neighbourdim) matid = this->Material()->Id();
    else 
    {
      if (thisdim < neighbourdim) matid = this->Material()->Id();
      else matid = list[0].Element()->Material()->Id();
    }
    
    TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid);
    //isto acertou as vizinhanas da interface geom�rica com o atual
    int index;

    TPZCompEl *list0 = list[0].Element();
    int list0side = list[0].Side();
    TPZCompElDisc * thisdisc  = dynamic_cast<TPZCompElDisc*>(this);
    TPZCompElDisc * neighdisc = dynamic_cast<TPZCompElDisc*>(list0);
    int thisside = -1;//discontinuous elements
    if (!thisdisc) thisside = side;
    int neighside = - 1;//discontinuous elements
    if (!neighdisc) neighside = list0side;



    if(Dimension() > list0->Dimension()){
       //o de volume �o direito caso um deles seja BC
       //a normal aponta para fora do contorno
       newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,this,list0, thisside, neighside);
    } else {
       //caso contr�io ou caso ambos sejam de volume 
       newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,list0,this, neighside, thisside);
    }
    LOGPZ_WARN(logger, "Exiting CreateInterface - NULL interface used");
    return newcreatedinterface;
  }

//If there is no equal level element, we try the lower elements.
//Higher elements will not be considered by this method. In that case the interface must be created by the neighbour.
  TPZCompElSide lower = thisside.LowerLevelElementList(0);
  if(lower.Exists()){
    LOGPZ_DEBUG(logger, "lower exists");
    //Interface has the same material of the neighbour with lesser dimension.
    //It makes the interface has the same material of boundary conditions (TPZCompElDisc with interface dimension)
    int matid;
    int thisdim = this->Dimension();
    int neighbourdim = lower.Element()->Dimension();
    if (thisdim == neighbourdim) matid = this->Material()->Id();
    else 
    {
      if (thisdim < neighbourdim) matid = this->Material()->Id();
      else matid = lower.Element()->Material()->Id();
     }

    //existem esquerdo e direito: this e lower
    TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid);
    int index;



    TPZCompEl *lowcel = lower.Element();
    int lowside = lower.Side();
    TPZCompElDisc * thisdisc  = dynamic_cast<TPZCompElDisc*>(this);
    TPZCompElDisc * neighdisc = dynamic_cast<TPZCompElDisc*>(lowcel);
    int thisside = -1;//discontinuous elements
    if (!thisdisc) thisside = side;
    int neighside = - 1;//discontinuous elements
    if (!neighdisc) neighside = lowside;

    if(Dimension() > lowcel->Dimension()){
       //para que o elemento esquerdo seja de volume
       newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,this,lowcel,thisside,neighside);
    } else {
       newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,lowcel,this,neighside,thisside);
    }
    LOGPZ_INFO(logger, "Exiting CreateInterface");
    return newcreatedinterface;
  }
  LOGPZ_INFO(logger, "Exiting CreateInterface");
  return newcreatedinterface;
}


void TPZCompEl::RemoveInterfaces(){

  int nsides = Reference()->NSides();
  if (!this->Material()){
    std::stringstream mess;
    mess << __PRETTY_FUNCTION__ << " - this->Material() == NULL, I can't RemoveInterfaces()";
    PZError << mess.str() << std::endl;
    LOGPZ_ERROR(logger, mess.str());
    return;  
  }
  int InterfaceDimension = this->Material()->Dimension() - 1;
  int is;
  TPZStack<TPZCompElSide> list,equal;
  for(is=0;is<nsides;is++){
    TPZCompElSide thisside(this,is);
    if(thisside.Reference().Dimension() != InterfaceDimension) continue;
    // procurar na lista de elementos iguais
    list.Resize(0);// o lado atual �uma face
    //thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
    RemoveInterface(is);// chame remove interface do elemento atual (para o side atual)
    thisside.HigherLevelElementList(list,0,0);// procurar na lista de elementos menores (todos)
    int size = list.NElements(),i;            // 'isto pode incluir elementos interfaces'
    //tirando os elementos de interface da lista
    for(i=0;i<size;i++){
      if(list[i].Element()->Type() == EInterface) list[i] = TPZCompElSide();//tirando interface
    }
    for(i=0;i<size;i++){// percorre os elementos menores
      if(!list[i].Element()) continue;
      TPZGeoElSide geolist = list[i].Reference();//TESTE
      if(geolist.Dimension() != InterfaceDimension) continue;
      equal.Resize(0);// para cada elemento menor e' preciso verificar a dimensao,
      list[i].EqualLevelElementList(equal,0,0);//montar a lista de elementos iguais (todos)
      equal.Push(list[i]);//n� �incorporado no m�odo anterior
      int neq = equal.NElements(),k=-1;
      while(++k < neq) if(equal[k].Element()->Type() != EInterface) break;//procurando elemento descont�uo cujo
      if(!neq || k == neq){                               //lado faz parte da parti� do lado side do this
	PZError << __PRETTY_FUNCTION__ << " - inconsistency of data";
	exit(-1);//elemento descont�uo n� achado: ERRO
      }// chame removeinterface do elemento menor

      equal[k].Element()->RemoveInterface(equal[k].Side());
    }
  }

}

void TPZCompEl::RemoveInterface(int side) {
  
  TPZStack<TPZCompElSide> list;
  list.Resize(0);
  TPZCompElSide thisside(this,side);
  thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
  int size = list.NElements(),i=-1;
  while(++i < size) if(list[i].Element()->Type() == EInterface) break;// procura aquele que e derivado de TPZInterfaceEl
  if(!size || i == size){
    return;// nada a ser feito
  }
  // aqui existe a interface
  TPZCompEl *cel = list[i].Element();
  TPZGeoEl *gel = cel->Reference();
  gel->RemoveConnectivities();// deleta o elemento das vizinhancas
  TPZGeoMesh *gmesh = Mesh()->Reference();
  int index = gmesh->ElementIndex(gel);// identifica o index do elemento
  gmesh->ElementVec()[index] = NULL;
  delete cel;
  delete gel;// deleta o elemento
  gmesh->ElementVec().SetFree(index);// Chame SetFree do vetor de elementos da malha geometrica para o index

}



// TPZCompElSide ///////////////////////////////////////////////////////////////
////////////// TPZCompElSide ///////////////////////////////////////////////////
////////////////////////// TPZCompElSide ///////////////////////////////////////
//TPZCompElSide::~TPZCompElSide() {}

TPZCompElSide::TPZCompElSide() {
  LOGPZ_INFO(loggerSide, "Entering empty constructor");
  fEl = 0;
  fSide = -1;
  //  georeftest = 0;
  LOGPZ_INFO(loggerSide, "Exiting empty constructor");
}

TPZCompElSide::TPZCompElSide(const TPZCompElSide &celside) {
  LOGPZ_INFO(loggerSide, "Entering copy constructor");
  fEl = celside.fEl;
  fSide   = celside.fSide;
  //  georeftest = celside.georeftest;
  LOGPZ_INFO(loggerSide, "Exiting copy constructor");
}

TPZCompElSide::TPZCompElSide(TPZCompEl *cel,int side) {
  LOGPZ_INFO(loggerSide, "Entering cel,side constructor");
  fEl = cel;
  fSide   = side;
  //  if(cel) georeftest = cel->Reference();
  //  else georeftest = 0;
  LOGPZ_INFO(loggerSide, "Exiting cel,side constructor");
}

TPZGeoElSide TPZCompElSide::Reference() const {
  LOGPZ_INFO(loggerSide, "Entering Reference");
  if(!fEl) {
    LOGPZ_WARN(loggerSide, "Exiting Reference - non initialized side element reached");
    return TPZGeoElSide();
  }
  TPZGeoElSide sideel(fEl->Reference(),fSide);
  LOGPZ_INFO(loggerSide, "Exiting Reference");
  return sideel;
}

void TPZCompElSide::SetSide(int side) {
  LOGPZ_INFO(loggerSide, "Entering SetSide");
  fSide = side;
  Reference().SetSide(side);
  LOGPZ_INFO(loggerSide, "Exiting SetSide");
}

void TPZCompElSide::HigherLevelElementList(TPZStack<TPZCompElSide> &elvec,
                                           int onlyinterpolated, int removeduplicates) {
  LOGPZ_INFO(loggerSide, "Entering HigherLevelElementList");
  TPZGeoElSide georef = Reference();
  if(!georef.Element()) {
    LOGPZ_WARN(loggerSide, "Exiting HigherLevelElementList - null reference reached");
    return;
  }
  LOGPZ_INFO(loggerSide, "Exiting HigherLevelElementList");
  georef.HigherLevelCompElementList2(elvec,onlyinterpolated,removeduplicates);
}

void TPZCompElSide::EqualLevelElementList(TPZStack<TPZCompElSide> &elsidevec,
                                          int onlyinterpolated, int removeduplicates) {
  LOGPZ_INFO(loggerSide, "Entering EqualLevelElementList");
  TPZGeoElSide georef = Reference();
  if(!georef.Exists()) {
    LOGPZ_WARN(loggerSide, "Exiting EqualLevelElementList - null reference reached");
    return;
  }
  LOGPZ_INFO(loggerSide, "Exiting EqualLevelElementList");
  georef.EqualLevelCompElementList(elsidevec,onlyinterpolated,removeduplicates);
}

void TPZCompElSide::HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, 
                                               int onlyinterpolated, int removeduplicates) {
  LOGPZ_INFO(loggerSide, "Entering HigherDimensionElementList");
  TPZGeoElSide georef = Reference();
  if(!georef.Exists()) {
    LOGPZ_INFO(loggerSide, "Entering HigherDimensionElementList - null reference reached");
    return;
  }
  georef.HigherDimensionElementList(elsidevec,onlyinterpolated);
  if(removeduplicates) RemoveDuplicates(elsidevec);
  LOGPZ_INFO(loggerSide, "Exiting HigherDimensionElementList");
}


void TPZCompElSide::RemoveDuplicates(TPZStack<TPZCompElSide> &elvec) {
  LOGPZ_INFO(loggerSide, "Entering RemoveDuplicates");
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
        LOGPZ_ERROR(loggerSide, "case not identified by RemoveConnectDuplicates");
        for(k=j;k<nelem-1;k++) elvec[k] = elvec[k+1];
        elvec.Resize(nelem-1);
        nelem--;
        j--;
      }
    }
  }
  LOGPZ_INFO(loggerSide, "Exiting RemoveDuplicates");
}

void TPZCompElSide::ConnectedElementList(TPZStack<TPZCompElSide> &ellist,
                                         int onlyinterpolated, int removeduplicates) {
  LOGPZ_INFO(loggerSide, "Entering ConnectedElementList");
  TPZCompElSide father = LowerLevelElementList(onlyinterpolated);
  if(father.Exists()) ellist.Push(father);
  EqualLevelElementList(ellist,onlyinterpolated,removeduplicates);
  HigherLevelElementList(ellist,onlyinterpolated,removeduplicates);
  LOGPZ_INFO(loggerSide, "Exiting ConnectedElementList");
}

//retorna o elem. de menor nivel ao qual estou restrito
TPZCompElSide TPZCompElSide::LowerLevelElementList(int onlyinterpolated) {
  LOGPZ_INFO(loggerSide, "Entering LowerLevelElementList");
  TPZGeoElSide georef = Reference();
  if(!georef.Exists()) {
    LOGPZ_WARN(loggerSide, "Exiting LowerLevelElementList - null reference reached");
    return TPZCompElSide();
  }
  LOGPZ_INFO(loggerSide, "Exiting LowerLevelElementList");
  return georef.LowerLevelCompElementList2(onlyinterpolated);
}

void TPZCompElSide::ExpandConnected(TPZStack<TPZCompElSide> &expandvec,int onlyinterpolated) {
  LOGPZ_INFO(loggerSide, "Entering ExpandConnected");
  //insidevec contem todos os elem. de maior nivel que dependem de mim,
  //obtidos com HigherLevelElementList(..)
  TPZGeoElSide georef = Reference();
  if(!georef.Exists()) {
    LOGPZ_WARN(loggerSide, "Exiting ExpandConnected - null reference reached");
    return;
  }
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
  LOGPZ_INFO(loggerSide, "Exiting ExpandConnected");
}

TPZCompElSide TPZCompElSide::LowerIdElementList(TPZCompElSide &expandvec,int onlyinterpolated) {
  LOGPZ_INFO(loggerSide, "Entering LowerIdElementList");
  if(!expandvec.Element()) {
    LOGPZ_WARN(loggerSide, "Exiting LowerIdElementList - empty list");
    return TPZCompElSide();
  }
  TPZGeoElSide gelside = expandvec.Reference();
  TPZStack<TPZGeoElSide> neighbourset;
  gelside.AllNeighbours(neighbourset);
  neighbourset.Push(gelside);
  int lowidindex = neighbourset.NElements()-1;
  //  TPZGeoElSide neighbour = gelside.Neighbour(),
  TPZGeoElSide lowidneigh(gelside);
  //if(!neighbour.Element()) return expandvec;
  int lowid = gelside.Id();
  int in, nneigh = neighbourset.NElements()-1;
  while(in < nneigh) {
    TPZCompEl *ref = neighbourset[in].Reference().Element();    
    if(neighbourset[in].Id() < lowid && ref && (!onlyinterpolated || dynamic_cast<TPZInterpolatedElement*>(ref)    )) {
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
  LOGPZ_INFO(loggerSide, "Exiting LowerIdElementList");
  return lowidneigh.Reference();
}

void TPZCompElSide::RemoveConnectDuplicates(TPZStack<TPZCompElSide> &expandvec){
  LOGPZ_INFO(loggerSide, "Entering RemoveConnectDuplicates");
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
  LOGPZ_INFO(loggerSide, "Exiting RemoveConnectDuplicates");
}

