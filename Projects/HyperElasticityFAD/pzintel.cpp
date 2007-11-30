// -*- c++ -*-
#include "pzintel.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pztrnsform.h"
#include "pztransfer.h"
#include "pzbndcond.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzquad.h"
#include "pzelmat.h"
#include "pzmat1dlin.h"
//#include "pztempmat.h"
#include "time.h"
#include "pzmanvector.h"
#include "pzblockdiag.h"
#include "pzcheckrestraint.h"
#include "pzdebug.h"

TPZInterpolatedElement::TPZInterpolatedElement(TPZCompMesh &mesh, TPZGeoEl *reference, int &index) :
  TPZCompEl(mesh,index) {
  fReference = reference;
  int materialid = reference->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
}

TPZInterpolatedElement::~TPZInterpolatedElement() {
}

int TPZInterpolatedElement::MaterialId() {
  if(!fMaterial) {
    PZError << "TPZIntEl::MaterialIndex. This element has no material.\n";
    return -1;
  }

  return fMaterial->Id();
}

int TPZInterpolatedElement::NShapeF() {
  int nn = NConnects();
  int in,res=0;
  for(in=0;in<nn;in++) res += NConnectShapeF(in);
  return res;
}

int TPZInterpolatedElement::NSideShapeF(int side) {
  int ns = NSideConnects(side);
  int in,res=0;
  for(in=0;in<ns;in++) res+= NConnectShapeF(SideConnectLocId(in,side));
  return res;
}


int TPZInterpolatedElement::MidSideConnectLocId(int side) {
  int il = 1 + NConnects() - (fReference->NSides());
  int nodloc = SideConnectLocId(NSideConnects(side)-il,side);
  return nodloc;
}


int TPZInterpolatedElement::SideConnectIndex(int connect, int side) {
  return ConnectIndex(SideConnectLocId(connect,side));
}

TPZConnect *TPZInterpolatedElement::SideConnect(int connect,int side) {
  if(side<0 || connect<0 || side>fReference->NSides()) {
    PZError << "TPZIntEl::SideConnect has bad first or second parameter.\n";
    return 0;
  }
  return &(fMesh->ConnectVec()[SideConnectIndex(connect,side)]);
}



void TPZInterpolatedElement::IdentifySideOrder(int side){
  TPZCompElSide thisside(this,side);
  TPZCompElSide large = thisside.LowerLevelElementList(1);
  int sideorder = SideOrder(side);
  int neworder;
  TPZStack<TPZCompElSide> elvec;
  thisside.EqualLevelElementList(elvec,1,0);
  elvec.Push(thisside);
  int cap,il,computorder;
  TPZInterpolatedElement *equal;
  int equalside;
  if(large.Exists()) {
    // There is a larger element
    // identify the order of the larger element and set the interpolation order
    // of the current element and all its neighbours of equal level to this order
    TPZInterpolatedElement *largel = (TPZInterpolatedElement *) large.Element();
    neworder = largel->SideOrder(large.Side());
    // We assume the datastructure of the elements is consistent in the sense
    // that if the current element has the same side order than the large
    //  element, then all its neighbours will also have the same order
    if(neworder != sideorder) {
      RemoveSideRestraintWithRespectTo(side,large);
      cap = elvec.NElements();
      for(il = 0; il<cap; il++) {
	equal = (TPZInterpolatedElement *) elvec[il].Element();
	equalside = elvec[il].Side();
	if(equal->ConnectIndex(equalside) != -1) {
	  equal->SetSideOrder(equalside,neworder);
	}
      }
      if(largel->ConnectIndex(large.Side()) != -1) {
	RestrainSide(side,largel,large.Side());
      }
    }
    computorder = neworder;
  } else {
    // There is no larger element connected to the side
    // identify the new side order by comparing the orders of the equal level elements
    neworder = ComputeSideOrder(elvec);
    // Verify is the side order of all elements is equal to neworder
    cap = elvec.NElements();
    il = 0;
    computorder = neworder;//Cedric
    while(il<cap) {//SideOrder(int side)
      equal = (TPZInterpolatedElement *) elvec[il].Element();
      equalside = elvec[il].Side();
      int equalorder = equal->SideOrder(equalside);
      if(equalorder != neworder) {
	computorder = neworder+1;//Cedric
      }
      il++;
    }
  }
  if(neworder != sideorder || neworder != computorder) {//Cedric :  || neworder != computorder
    // The order of the current element changed
    // Therefore the constraints of all smaller elements connected to the current
    // element will need to be adapted
    // The following loop will happen twice, but doesn't affect anything

    elvec.Resize(0);
    thisside.EqualLevelElementList(elvec,1,0);
    elvec.Push(thisside);
    cap = elvec.NElements();
    for(il=0; il<cap; il++) {
      equal = (TPZInterpolatedElement *) elvec[il].Element();
      equalside = elvec[il].Side();
      if(equal->ConnectIndex(equalside) != -1) {
	equal->SetSideOrder(equalside,neworder);
      }
    }
    TPZStack<TPZCompElSide> highdim;
    if(thisside.Reference().Dimension() == 1) {
      for(il=0; il<cap; il++) {
	elvec[il].HigherDimensionElementList(highdim,1,1);
      }
    }
    elvec.Resize(0);
    // Adapt the restraints of all smaller elements connected to the current side
    thisside.HigherLevelElementList(elvec,1,1);
//    thisside.ExpandConnected(elvec,1);
    cap = elvec.NElements();
    TPZInterpolatedElement *small;
    int smallside;
    int dimension;
    for(dimension = 0; dimension < 4; dimension++) {
      for(il=0; il<cap; il++) {
	if(elvec[il].Reference().Dimension() != dimension) continue;
	small = (TPZInterpolatedElement *) elvec[il].Element();
	smallside = elvec[il].Side();
	// Identify Side Order is used because it will call itself recursively
	// for its smaller elements too.
	small->IdentifySideOrder(smallside);
      }
    }
	// look for a large element with dimension 2
	while(large.Exists() && large.Reference().Dimension() <2) {
		large = large.LowerLevelElementList(1);
	}
    cap = highdim.NElements();
    for(il=0; il<cap; il++) {
      TPZInterpolatedElement *el;
      int highside;
      el = dynamic_cast<TPZInterpolatedElement *> (highdim[il].Element());
      highside = highdim[il].Side();
      int order, comporder;
      order = el->SideOrder(highside);
      TPZStack<TPZCompElSide> equallist;
      highdim[il].EqualLevelElementList(equallist,1,0);
      equallist.Push(highdim[il]);
      TPZCompElSide highlarge = highdim[il].LowerLevelElementList(1);
      // when the  element highdim is restricted by the same element as the original side, do nothing
      // the restriction will be taken care of in the future
	  //EXPERIMENTAL
	  if(highlarge.Exists()) continue;
      if(highlarge.Element() == large.Element() && highlarge.Side() == large.Side() && large.Element()) continue;
      if(highlarge.Exists()) {
	TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (highlarge.Element());
	comporder = cel->SideOrder(highlarge.Side());
      } else {
	comporder = el->ComputeSideOrder(equallist);
      }
      if(order != comporder) {
	el->IdentifySideOrder(highside);
      } else {
	el->RecomputeRestraints(highside);
      }
    }
  }
}

void TPZInterpolatedElement::RecomputeRestraints(int side) {
  TPZStack<TPZCompElSide> elvec;
  elvec.Resize(0);
  TPZCompElSide thisside(this,side);
  // Adapt the restraints of all smaller elements connected to the current side
  thisside.HigherLevelElementList(elvec,1,1);
//  thisside.ExpandConnected(elvec,1);
  int cap,il;
  cap = elvec.NElements();
  TPZInterpolatedElement *small;
  int smallside;
  int dimension;
  for(dimension =0; dimension <4; dimension++) {
    for(il=0; il<cap; il++) {
      if(elvec[il].Reference().Dimension() != dimension) continue;
      small = (TPZInterpolatedElement *) elvec[il].Element();
      smallside = elvec[il].Side();
      // Identify Side Order is used because it will call itself recursively
      // for its smaller elements too.
      small->RemoveSideRestraintWithRespectTo(smallside,thisside);
      small->RestrainSide(smallside,this,side);
    }
  }
}

void TPZInterpolatedElement::UpdateNeighbourSideOrder(int side, TPZVec<TPZCompElSide> &elvec) {
  int orside = SideOrder(side);
  int neighbourside;
  TPZInterpolatedElement *neighbour;
  int il;
  int cap = elvec.NElements();
  for(il=0; il<cap; il++) {
    neighbour = dynamic_cast<TPZInterpolatedElement *> ( elvec[il].Element() );
    if(!neighbour) continue;
    neighbourside = elvec[il].Side();
    int neighord = neighbour->SideOrder(neighbourside);
    if(orside != neighord) neighbour->IdentifySideOrder(neighbourside);
  }
}

void TPZInterpolatedElement::BuildTransferMatrix(TPZInterpolatedElement &coarsel, TPZTransform &t, TPZTransfer &transfer){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t
  int locnod = NConnects();
  int cornod = coarsel.NConnects();
  int locmatsize = NShapeF();
  int cormatsize = coarsel.NShapeF();

  // compare interpolation orders
  // the minimum interpolation order of this needs to be larger than the maximum interpolation order of coarse

  int myminorder = SideOrder(locnod-1);
  int ncorners = NCornerConnects();
  int ic;
  for(ic=ncorners; ic<locnod; ic++) {
    if(SideOrder(ic) < myminorder) myminorder = SideOrder(ic);
  }
  ncorners = coarsel.NCornerConnects();
  int coarsemaxorder = 0;
  for(ic=ncorners; ic<cornod; ic++) {
    if(coarsel.SideOrder(ic) > coarsemaxorder) coarsemaxorder = coarsel.SideOrder(ic);
  }
  if(coarsemaxorder > myminorder) {
    //    cout << "TPZInterpolatedElement compute the transfer matrix coarse " << coarsemaxorder << " me " << myminorder << endl;
    return;
  }
  TPZStack<int> connectlistcoarse,dependencyordercoarse, corblocksize;
  connectlistcoarse.Resize(0);
  dependencyordercoarse.Resize(0);
  corblocksize.Resize(0);
  for(ic=0; ic<cornod; ic++) connectlistcoarse.Push(coarsel.ConnectIndex(ic));
  coarsel.BuildConnectList(connectlistcoarse);
  coarsel.BuildDependencyOrder(connectlistcoarse,dependencyordercoarse);
  cornod = connectlistcoarse.NElements();
  int nvar = coarsel.Material()->NStateVariables();
  TPZBlock corblock(0,cornod);
  int in;
  for(in = 0; in < NConnects(); in++) {
    int blsize = coarsel.NConnectShapeF(in);
    corblock.Set(in,blsize);
    corblocksize.Push(blsize);
  }
  int c;
  for(;in<cornod; in++) {
    c = connectlistcoarse[in];
    int blsize = coarsel.Mesh()->ConnectVec()[c].NDof(*(coarsel.Mesh()))/nvar;
    corblock.Set(in,blsize);
    corblocksize.Push(blsize);
    cormatsize += blsize;
  }
  corblock.Resequence();

  REAL loclocmatstore[500] = {0.};
  TPZFMatrix loclocmat(locmatsize,locmatsize,loclocmatstore,500);
  TPZFMatrix loccormat(locmatsize,cormatsize);
  loclocmat.Zero();
  loccormat.Zero();

  TPZIntPoints &intrule = GetIntegrationRule();
  int dimension = Dimension();

  TPZManVector<int> prevorder(dimension),
    order(dimension);
  intrule.GetOrder(prevorder);

  TPZManVector<int> interpolation(dimension);
  GetInterpolationOrder(interpolation);

  // compute the interpolation order of the shapefunctions squared
  int dim;
  int maxorder = interpolation[0];
  for(dim=0; dim<interpolation.NElements(); dim++) {
    maxorder = interpolation[dim] < maxorder ? maxorder : interpolation[dim];
  }
  for(dim=0; dim<dimension; dim++) {
    order[dim] = maxorder*2+2;
  }
  intrule.SetOrder(order);

  TPZBlock locblock(0,locnod);

  for(in = 0; in < locnod; in++) {
    locblock.Set(in,NConnectShapeF(in));
  }
  locblock.Resequence();

  REAL locphistore[50]={0.},locdphistore[150]={0.};
  TPZFMatrix locphi(locmatsize,1,locphistore,50);
  TPZFMatrix locdphi(dimension,locmatsize,locdphistore,150);
  locphi.Zero();
  locdphi.Zero();
  // derivative of the shape function
  // in the master domain

  TPZFMatrix corphi(cormatsize,1);
  TPZFMatrix cordphi(dimension,cormatsize);
  // derivative of the shape function
  // in the master domain

  REAL jacobianstore[9],
    axesstore[9];
  TPZManVector<REAL> int_point(dimension),
    coarse_int_point(dimension);
  TPZFMatrix jacobian(dimension,dimension,jacobianstore,9),jacinv(dimension,dimension);
  TPZFMatrix axes(3,3,axesstore,9);
  TPZManVector<REAL> x(3);
  int_point.Fill(0.,0);
  REAL jac_det = 1.;
  fReference->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
  REAL multiplier = 1./jac_det;

  int numintpoints = intrule.NPoints();
  REAL weight;
  int lin,ljn,cjn;

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule.Point(int_ind,int_point,weight);
    fReference->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
    fReference->X(int_point, x);
    Shape(int_point,locphi,locdphi);
    weight *= jac_det;
    t.Apply(int_point,coarse_int_point);
    corphi.Zero();
    cordphi.Zero();
    coarsel.Shape(coarse_int_point,corphi,cordphi);

    coarsel.ExpandShapeFunctions(connectlistcoarse,dependencyordercoarse,corblocksize,corphi,cordphi);

    for(lin=0; lin<locmatsize; lin++) {
      for(ljn=0; ljn<locmatsize; ljn++) {
	loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0)*multiplier;
      }
      for(cjn=0; cjn<cormatsize; cjn++) {
	loccormat(lin,cjn) += weight*locphi(lin,0)*corphi(cjn,0)*multiplier;
      }
    }
    jacobian.Zero();
  }
  loclocmat.SolveDirect(loccormat,ELDLt);


  for(in=0; in<locnod; in++) {
    //    int cind = connectlistcoarse[in];
    if(Connect(in).HasDependency()) continue;
    int locblocknumber = Connect(in).SequenceNumber();
    int locblocksize = locblock.Size(in);
    int locblockpos = locblock.Position(in);
    TPZStack<int> locblockvec;
    TPZStack<int> globblockvec;
    int numnonzero = 0,jn;
    //      if(transfer.HasRowDefinition(locblocknumber)) continue;

    for(jn = 0; jn<cornod; jn++) {
      int corblocksize = corblock.Size(jn);
      int corblockpos = corblock.Position(jn);
      int cind = connectlistcoarse[jn];
      TPZConnect &con = coarsel.Mesh()->ConnectVec()[cind];
      if(con.HasDependency()) continue;
      int corblocknumber = con.SequenceNumber();
      if(locblocksize == 0 || corblocksize == 0) continue;
      TPZFMatrix small(locblocksize,corblocksize,0.);
      loccormat.GetSub(locblockpos,corblockpos,
		       locblocksize,corblocksize,small);
      REAL tol = Norm(small);
      if(tol >= 1.e-10) {
	locblockvec.Push(jn);
	globblockvec.Push(corblocknumber);
	numnonzero++;
      }
    }
    if(transfer.HasRowDefinition(locblocknumber)) continue;
    transfer.AddBlockNumbers(locblocknumber,globblockvec);
    int jnn;
    for(jnn = 0; jnn<numnonzero; jnn++) {
      jn = locblockvec[jnn];
      int corblocksize = corblock.Size(jn);
      int corblockpos = corblock.Position(jn);
      if(corblocksize == 0 || locblocksize == 0) continue;
      TPZFMatrix small(locblocksize,corblocksize,0.);
      loccormat.GetSub(locblockpos,corblockpos,locblocksize,corblocksize,small);
      transfer.SetBlockMatrix(locblocknumber,globblockvec[jnn],small);
    }
  }
  intrule.SetOrder(prevorder);

}

int TPZInterpolatedElement::CreateMidSideConnect(int side) {
  TPZCompMesh *cmesh = Mesh();
  TPZMaterial *mat = Material();
  int nvar = 1;
  if(mat) nvar = mat->NStateVariables();
  int newnodeindex;
  int il;
  int nodloc = MidSideConnectLocId(side);


  TPZStack<TPZCompElSide> elvec;
  TPZCompElSide thisside(this,side);
  if(side < NCornerConnects()) {
    thisside.EqualLevelElementList(elvec,0,0);
    int nel = elvec.NElements();                // (1)
    if(nel && elvec[nel-1].Reference().Dimension() == thisside.Reference().Dimension()){
      newnodeindex =  elvec[nel-1].ConnectIndex();
      SetConnectIndex(nodloc,newnodeindex);
    } else {
      newnodeindex = cmesh->AllocateNewConnect();
      TPZConnect &newnod = cmesh->ConnectVec()[newnodeindex];
      if(newnod.HasDependency()) {
	cout << "TPZInterpolatedElement::CreateMidSideConnect, new node has dependency\n";
	newnod.Print(*cmesh);
      }
      int seqnum = newnod.SequenceNumber();
      cmesh->Block().Set(seqnum,nvar*NConnectShapeF(nodloc));
      SetConnectIndex(nodloc,newnodeindex);   //Is true only to one-dimensional case
      // We created a new node, check whether the node needs to be constrained
      TPZCompElSide father = thisside.LowerLevelElementList(1);
      if(father.Exists()) {
      	int side_neig = father.Side();
      	TPZInterpolatedElement *cel = (TPZInterpolatedElement *) father.Element();
      	RestrainSide(side,cel,side_neig);
      }
    }
    return newnodeindex;
  }

  // Connect looks for a connecting element of equal or lower level
  TPZInterpolatedElement *cel = 0;
  int side_neig = 0;
  thisside.EqualLevelElementList(elvec,1,1);
  int nelem = elvec.NElements();
  // find an element in the list which is interpolated
  if(nelem) {
    cel = (TPZInterpolatedElement *) elvec[0].Element();
    side_neig = elvec[0].Side();
  }
  int newnodecreated = 0;
  if(cel) {
    newnodeindex = cel->ConnectIndex(side_neig);
    SetConnectIndex(nodloc,newnodeindex);
  } else {
    newnodeindex = cmesh->AllocateNewConnect();
    TPZConnect &newnod = cmesh->ConnectVec()[newnodeindex];
    int seqnum = newnod.SequenceNumber();
    cmesh->Block().Set(seqnum,nvar*NConnectShapeF(nodloc));
    SetConnectIndex(nodloc,newnodeindex);
    newnodecreated = 1;
  }
  TPZCompElSide father = thisside.LowerLevelElementList(1);
  if(father.Exists() && newnodecreated) {
    // There is a larger element connected along the side that is being created
    // We will adopt the interpolation order of the father and compute the
    // restraints
    // IT IS ASSUMED THAT ALL NEIGHBOURING ELEMENTS AND SMALLER ELEMENTS CONNECTED
    // TO THIS SIDE HAVE CONSISTENT INTERPOLATION ORDERS
    side_neig = father.Side();
    TPZInterpolatedElement *elfather = (TPZInterpolatedElement *) father.Element();
    int sideorder = elfather->SideOrder(side_neig);
    SetSideOrder(side,sideorder);
    RestrainSide(side,elfather,side_neig);
    elvec.Resize(0);
    thisside.HigherLevelElementList(elvec,1,1);
//    thisside.ExpandConnected(elvec,1);
//    thisside.RemoveDuplicates(elvec);
    int cap = elvec.NElements();
    int dim;
    for(dim=0; dim<3; dim++) {
      for(il=0; il<cap; il++) {
	if(elvec[il].Reference().Dimension() != dim) continue;
	cel = (TPZInterpolatedElement *) elvec[il].Element();
	Reference()->ResetReference();
	cel->RemoveSideRestraintWithRespectTo(elvec[il].Side(),father);
	Reference()->SetReference(this);

	cel->RestrainSide(elvec[il].Side(),this,side);
      }
    }
  } else if(!father.Exists() && !newnodecreated){
    // The element does not have a larger element connected along the side
    // The insertion of the new element may have an effect on the interpolation
    // order of all equal level elements which are connected
    elvec.Resize(0);
    thisside.EqualLevelElementList(elvec,1,0);
    int oldorder = -1;
    if(elvec.NElements()) {
      TPZInterpolatedElement *cel = (TPZInterpolatedElement *) elvec[0].Element();
      oldorder = cel->SideOrder(elvec[0].Side());
    }
    elvec.Push(thisside);
    int sideorder = ComputeSideOrder(elvec);
    if(sideorder != oldorder) {
      TPZInterpolatedElement *cel = (TPZInterpolatedElement *) elvec[0].Element();
      cel->IdentifySideOrder(elvec[0].Side());
    } else {
      SetSideOrder(side,sideorder);
    }
  } else if(!father.Exists() && newnodecreated) {
    elvec.Resize(0);
    thisside.HigherLevelElementList(elvec,1,1);
//    thisside.ExpandConnected(elvec,1);
//    thisside.RemoveDuplicates(elvec);
    int cap = elvec.NElements();
    int sideorder = PreferredSideOrder(side);
    SetSideOrder(side,sideorder);
    // We check on all smaller elements connected to the current element
    // whether their interpolation order is consistent with the interpolation
    // order of the current element/side
    int dim;
    for(dim=0; dim<3; dim++) {
      for(il=0; il<cap; il++) {
	if(elvec[il].Reference().Dimension() != dim) continue;
	TPZInterpolatedElement *cel = (TPZInterpolatedElement *) elvec[il].Element();
	// the restraint will only be applied if the first
	// element in Connect is equal to the current element
	int celsideorder = cel->SideOrder(elvec[il].Side());
	if(celsideorder == sideorder) {
	  // if the interpolation order remains unchanged, just restrain the side
	  cel->RestrainSide(elvec[il].Side(),this,side);
	} else {
	  // if the interpolation order changed, recompute the side order of all
	  // smaller elements and recompute the constraints
	  cel->IdentifySideOrder(elvec[il].Side());
	}
      }
    }
  }
  return newnodeindex;
}

void TPZInterpolatedElement::RestrainSide(int side, TPZInterpolatedElement *large, int neighbourside) {
  if(!IsConnectContinuous(side)) return;
  TPZCompElSide thisside(this,side);
  TPZCompElSide locallarge = thisside.LowerLevelElementList(1);
  TPZCompElSide largecompside(large,neighbourside);
  TPZGeoElSide locallargeref = locallarge.Reference();
  TPZGeoElSide largecompsideref = largecompside.Reference();
  if(!locallarge.Exists() || !locallargeref.NeighbourExists(largecompsideref)) {
    PZError << "TPZInterpolatedElement::RestrainSide called for a wrong large element\n";
    return;
  }
  TPZInterpolatedElement *cel = 0;
  if(locallarge.Exists()) cel = (TPZInterpolatedElement *) locallarge.Element();
  if(!cel) {
    cout << "TPZInterpolatedElement::RestrainSide, I don't understand cel = 0\n";
    return;
/*  } else if(cel->Reference()->Level() < large->Reference()->Level()) {
    cout << "TPZInterpolatedElement::RestrainSide, I don't understand cel is larger\n";
    cout << "my level " << Reference()->Level() << " cel level " << cel->Reference()->Level()
	 << " large level " << large->Reference()->Level() << endl;
    cout << "my id " << Reference()->Id() << " cel id " << cel->Reference()->Id()
	 << " large id " << large->Reference()->Id() << endl;
    return;

  } else if(cel->Reference()->Level() > large->Reference()->Level()) {
    cout << "TPZInterpolatedElement::RestrainSide no restriction between " << Reference()->Id()
	 << " and " << large->Reference()->Id() << " because of " << cel->Reference()->Id() << endl;
    return;
*/
  }
  TPZConnect &myconnect = Connect(side);
  if(myconnect.HasDependency() && locallargeref.Dimension() > 0) {
    cout << "TPZInterpolatedElement unnecessary call to restrainside\n";
  }
  if (cel->ConnectIndex(locallarge.Side()) == -1){
    cout <<  "TPZInterpolatedElement::RestrainSide : Side of large element not initialized\n ";
    return;
  }
  if(locallargeref.Dimension() == 0) return;
  TPZGeoElSide thisgeoside(Reference(),side);
  TPZGeoElSide largeside = largecompside.Reference();
  TPZTransform t(thisgeoside.Dimension());
  thisgeoside.SideTransform3(largeside,t);
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(side,SideOrder(side)*2);
  if(!intrule) return;
  int numint = intrule->NPoints();
  int numshape = NSideShapeF(side);
  int numshapel = large->NSideShapeF(neighbourside);
  TPZFMatrix phis(numshape,1),dphis(2,numshape),phil(numshapel,1),dphil(2,numshapel);
  TPZFMatrix M(numshape,numshape,0.),MSL(numshape,numshapel,0.);
  TPZVec<REAL> par(3),pointl(3),point(3);
  int in,jn;
  REAL weight;
  for(int it=0; it<numint; it++) {
    intrule->Point(it,par,weight);
    SideShapeFunction(side,par,phis,dphis);
    t.Apply(par,pointl);
    large->SideShapeFunction(neighbourside,pointl,phil,dphil);
    for(in=0; in<numshape; in++) {
      for(jn=0; jn<numshape; jn++) {
	M(in,jn) += phis(in,0)*phis(jn,0)*weight;
      }
      for(jn=0; jn<numshapel; jn++) {
	MSL(in,jn) += phis(in,0)*phil(jn,0)*weight;
      }
    }
  }
  TPZStepSolver MSolve(&M);
  MSolve.SetDirect(ELU);
  MSolve.Solve(MSL,MSL);
  MSolve.ResetMatrix();
  int numsidenodes_small = NSideConnects(side);
  // Philippe 12/3/99
  //  int numsidenodes_large = NSideConnects(neighbourside);
  int numsidenodes_large = large->NSideConnects(neighbourside);
  TPZBlock MBlocksmall(0,numsidenodes_small), MBlocklarge(0,numsidenodes_large);
  for(in = 0; in<numsidenodes_small; in++) {
    int locid = SideConnectLocId(in,side);
    MBlocksmall.Set(in,NConnectShapeF(locid));
  }
  for(in = 0; in<numsidenodes_large; in++) {
    int locid = large->SideConnectLocId(in,neighbourside);
    MBlocklarge.Set(in,large->NConnectShapeF(locid));
  }

  MBlocksmall.Resequence();
  MBlocklarge.Resequence();
  TPZFMatrix blocknorm(numsidenodes_small,numsidenodes_large,0.);
  for(in = 0; in< numsidenodes_small; in++) {
    int ibl = MBlocksmall.Size(in);
    if(!ibl) continue;
    for(jn = 0; jn<numsidenodes_large; jn++) {
      int jbl = MBlocklarge.Size(jn);
      if(!jbl) continue;
      int i,j;
      int ipos = MBlocksmall.Position(in);
      int jpos = MBlocklarge.Position(jn);
      for(i=0; i<ibl; i++) for(j=0; j<jbl; j++) blocknorm(in,jn) += MSL(ipos+i,jpos+j)*MSL(ipos+i,jpos+j);
      blocknorm(in,jn) /= (ibl*jbl);
      blocknorm(in,jn) = sqrt(blocknorm(in,jn));
    }
  }
  CheckConstraintConsistency(side);
  TPZConnect &inod = Connect(side);
  int inodindex = ConnectIndex(side);
  int ndepend = 0;
  in = numsidenodes_small-1;
  for(jn = 0; jn<numsidenodes_large; jn++) {
    if(blocknorm(in,jn) < 1.e-8) continue;
    int jnodindex = large->SideConnectIndex(jn,neighbourside);
    inod.AddDependency(inodindex,jnodindex,MSL,MBlocksmall.Position(in),MBlocklarge.Position(jn),
		       MBlocksmall.Size(in),MBlocklarge.Size(jn));
    ndepend++;
  }

  if (! ndepend){
    //  cout << "Caso esquisito!!! Chame o Boss que vc receber� um pr�mio\n";
    for(jn = 0; jn<numsidenodes_large; jn++) {
      int jnodindex = large->SideConnectIndex(jn,neighbourside);
      inod.AddDependency(inodindex,jnodindex,MSL,MBlocksmall.Position(in),MBlocklarge.Position(jn),
			 MBlocksmall.Size(in),MBlocklarge.Size(jn));
      ndepend++;
    }
  }
  delete intrule;
  // a matriz frestraint deveria ser igual a MSL
  TPZCheckRestraint test(thisside,largecompside);
  //test.Print(cout);

  int imsl, jmsl;
  int rmsl = MSL.Rows();
  int cmsl = MSL.Cols();

  int rtest = test.RestraintMatrix().Rows();
  int ctest = test.RestraintMatrix().Cols();

  if (rtest!=rmsl || ctest!=cmsl){
    cout << "TPZInterpolatedElement::Error::Restraint matrix side incompatibility: MSL (rows,cols): ( " << rmsl
	 << " , " << cmsl << " )" << " RestraintMatrix (rows,cols): (" << rtest << " , "  << ctest << " )\n";
    int a;
    cin >> a;
    return;
  }

  TPZFMatrix mslc (MSL);
  mslc -= test.RestraintMatrix();

  REAL normmsl = 0.;
  for (imsl=0; imsl<rmsl; imsl++){
    for (jmsl=0; jmsl<cmsl; jmsl++){
      normmsl += sqrt(mslc(imsl,jmsl)*mslc(imsl,jmsl));
    }
  }
  if (normmsl > 1.E-6){
    cout << "TPZInterpolatedElement::Error::MSL matrix has non zero norm " << normmsl << "\n";
    mslc.Print("Difference Matrix ",cout);
    for (imsl=0; imsl<rmsl; imsl++){
      for (jmsl=0; jmsl<cmsl; jmsl++){
	if (fabs(MSL(imsl,jmsl) - test.RestraintMatrix()(imsl,jmsl)) > 1.E-6){
	  cout << "msl[ " << imsl << " , " << jmsl << " ] = " << MSL(imsl,jmsl) << "\t " << test.RestraintMatrix()(imsl,jmsl) << endl;
	}
      }
    }    int a;
    gDebug = 1;
    cin >> a;
  }


  // verificar a norma de MSL
  if(test.CheckRestraint()) {
    cout << "TPZInterpolatedElement::Error::Bad restraints detected\n";// recado de erro.
    test.Print(cout);
    int a;
    gDebug = 1;
    cin >> a;
    test.Diagnose();
    TPZCheckRestraint test2(thisside,largecompside);
  }

}

void TPZInterpolatedElement::CheckConstraintConsistency() {
  int nc = fReference->NSides();
  //  int a;
  for(int c=0; c<nc; c++) CheckConstraintConsistency(c);
}


int TPZInterpolatedElement::CheckElementConsistency(){
  int dimel = Dimension();
  int iside;
  int a;
  for (iside = 0; iside<(fReference->NSides()-1); iside++){
    TPZCompElSide celside(this,iside);
    int dimsmall = celside.Reference().Dimension();

    TPZIntPoints *sirule = Reference()->CreateSideIntegrationRule (iside,SideOrder(iside)*2);

    if (dimsmall >= dimel){
      cout << "TPZInterpolatedElement::CheckConstraintConsistency : dismall >= dimel: " << dimsmall << " >= " <<  dimel << endl;
      cin >> a;
    }

    int idim;

    int nshapes = NSideShapeF(iside);
    TPZFMatrix phis(nshapes,1);
    TPZFMatrix dphis(dimsmall,nshapes);

    for (idim = (dimsmall + 1); idim <= dimel; idim++){
      TPZStack <TPZGeoElSide> geoelsidevec;
      celside.Reference().Element()->AllHigherDimensionSides(iside,idim,geoelsidevec);

      int nelhigh = geoelsidevec.NElements();
      int inh;
      for (inh = 0; inh < nelhigh; inh++){

	int sidel =  geoelsidevec[inh].Side();
	int nshapel =  NSideShapeF(sidel);

	TPZFMatrix phil(nshapel, 1);
	TPZFMatrix dphil(idim,nshapel);

	int npts = sirule->NPoints();
	int ipt;

	TPZTransform transform (Reference()->SideToSideTransform(iside,sidel));

	for (ipt = 0; ipt<npts; ipt++){
	  TPZVec <REAL> pts(dimsmall);
	  TPZVec <REAL> ptl(idim);
	  REAL w;
	  sirule->Point(ipt,pts,w);
	  SideShapeFunction(iside,pts,phis,dphis);
	  transform.Apply(pts,ptl);
	  SideShapeFunction(iside,ptl,phil,dphil);
	  int check = CompareShapeF(iside,sidel,phis,dphis,phil,dphil,transform);
	  if (!check) return check;
	}
      }
    }
  }
  return 1;
}


int TPZInterpolatedElement::CompareShapeF(int sides, int sidel, TPZFMatrix &phis, TPZFMatrix &dphis, TPZFMatrix &phil, TPZFMatrix &dphil, TPZTransform &transform){
  int ncons = NSideConnects(sides);
  int nconl = NSideConnects(sidel);
  TPZVec<int> posl(nconl+1), poss(ncons+1);
  int icon;
  if(nconl) posl[0] =0;
  if(ncons) poss[0] = 0;
  for(icon=0; icon<nconl; icon++) {
    posl[icon+1] = posl[icon] + NConnectShapeF(SideConnectLocId(icon,sidel));
  }
  for(icon=0; icon<ncons; icon++) {
    poss[icon+1] = poss[icon] + NConnectShapeF(SideConnectLocId(icon,sides));
  }

  for (icon=0; icon<nconl; icon++){
    int conlindl = SideConnectLocId(icon,sidel);
    int conscounter = -1;
    //    int i;
    TPZCompElSide cellsmallside(this,sides);

    for(conscounter=0; conscounter<ncons; conscounter++) if(SideConnectLocId(conscounter,sides) == conlindl) break;
    if(conscounter != ncons) {
      int firsts = poss[conscounter];
      int firstl = posl[icon];

      int dimsmall = cellsmallside.Reference().Dimension();
      int idim;
      for (idim = 0; idim < dimsmall; idim++){
	if (phis(firsts+idim,0) != phil(firstl+idim,0)){
	  cout << "TPZInterpolatedElement::CompareShapeF : Inconsistent transformation side detected\n";
	  cout << "small side connect shape [ " << firsts+idim << " ] = " << phis(firsts+idim,0)
	       << "large side connect shape [ " << firstl+idim << " ] = " << phil(firstl+idim,0) << endl;
	  return 0;
	}
      }

      //como aplicar o transform a matrizes??
      TPZFMatrix left = transform.Mult()*dphil;
      left -= dphis;

      int m,n;
      for (m=0;m<left.Rows();m++)
	for (n=0;n<left.Cols();n++)
	  if (fabs(left(m,n)) > 1E-6){
	    cout << "TPZInterpolatedElement::CompareShapeF : Erro TrT * dphil not equal dphis\n"
		 << "dphis =  [ " << m << " , " << n << " ] = " << dphis(m,n) << " TrT*dphil = " << left(m,n);
	  }


    }
    else {
      int idim;
      int dimsmall = cellsmallside.Reference().Dimension();
      int firsts = poss[conscounter];

      for (idim = 0; idim < dimsmall; idim++){
	if (phis(firsts+idim,0) != 0){
	  cout << "TPZInterpolatedElement::CompareShapeF : Inconsistent transformation side detected\n";
	  cout << "small side connect shape [ " << firsts+idim << " ] = " << phis(firsts+idim,0) << "  Must be 0\n";
	  return 0;
	}
      }

      TPZFMatrix left = transform.Mult()*dphil;

      int m,n;
      for (m=0;m<left.Rows();m++)
	for (n=0;n<left.Cols();n++)
	  if (fabs(left(m,n)) > 1E-6){
	    cout << "TPZInterpolatedElement::CompareShapeF : Erro TrT * dphil not equal dphis\n"
		 << "dphis =  [ " << m << " , " << n << " ] = " << dphis(m,n) << " TrT*dphil = " << left(m,n);
	  }
    }
  }
  return 1;
}

void TPZInterpolatedElement::CheckConstraintConsistency(int side) {
  TPZCompElSide thisside(this,side);
  TPZCompElSide large = thisside.LowerLevelElementList(1);
  if(large.Exists()) {
    int largeside = large.Side();
    TPZInterpolatedElement *largel = (TPZInterpolatedElement *) large.Element();
    int nlargesideconnects = largel->NSideConnects(largeside);
    TPZConnect &thiscon = Connect(side);
    int jn;
    TPZConnect::TPZDepend *dep = thiscon.FirstDepend();
    while(dep) {
      for(jn = 0; jn < nlargesideconnects; jn++) {
	int largeindex = largel->SideConnectIndex(jn,largeside);
	if(dep->fDepConnectIndex == largeindex || largeindex == ConnectIndex(side)) break;
      }
      if(jn == nlargesideconnects) {
	cout << "TPZInterpolatedElement::Dependency inconsistency\n";
	thiscon.Print(*Mesh(),cout);
	largel->Print(cout);
      }
      dep = dep->fNext;
    }
  } else {
    TPZConnect &thiscon = Connect(side);
    if(thiscon.HasDependency()) {
      large = thisside.LowerLevelElementList(1);
      cout << "Dependency inconsistency\n";
      thiscon.Print(*Mesh(),cout);
    }
  }

}

void TPZInterpolatedElement::RemoveSideRestraintWithRespectTo(int side,
							      const TPZCompElSide & neighbour) {
  TPZCompElSide thisside(this,side);
  //  TPZGeoElSide thisgeoside = thisside.Reference();
  //  TPZGeoElSide largegeoside = neighbour.Reference();
  TPZCompElSide largeset = thisside.LowerLevelElementList(1);
  if(!largeset.Exists()) {
    PZError << "TPZInterpolatedElement::RemoveSideRestraintWithRespectTo inconsistent 1\n";
    return;
  }
//  int smallfatherlevel = largeset.Element()->Reference()->Level();
  // look for the first larger element connected to cel
//  if(smallfatherlevel > neighbour.Reference().Level() ) {
//    return;
//  } else if(smallfatherlevel < neighbour.Reference().Level() ) {
//    PZError << "TPZInterpolatedElement::RemoveSideRestraintWithRespectTo inconsistent 2\n";
//  }
  TPZInterpolatedElement *smallfather = (TPZInterpolatedElement *) largeset.Element();
  int smallfatherside = largeset.Side();
  if(!neighbour.Reference().NeighbourExists(largeset.Reference())) {
    PZError << "TPZInterpolatedElement::RemoveSideRestraintWithRespectTo inconsistent 3\n";
  }
  int js;
  int nsfather = smallfather->NSideConnects(smallfatherside);
  int dfiindex = ConnectIndex(side);
  TPZConnect *dfi = &Connect(side);
  for(js=0; js<nsfather; js++) {
    int dfjindex = smallfather->SideConnectIndex(js,smallfatherside);
    dfi->RemoveDepend(dfiindex,dfjindex);
  }
  if(dfi->HasDependency()) {
    int dim1,dim2;
    dim1 =  Reference()->SideDimension(side);
    dim2 =  neighbour.Element()->Reference()->SideDimension(neighbour.Side());
    if(dim1 || dim2) {
      PZError << "TPZInterpolatedElement::RemoveSideRestraintWithRespectTo fishy " << dfiindex << endl;
    }
  }
}

void TPZInterpolatedElement::RemoveSideRestraintsII(MInsertMode mode) {
  TPZCompElSide large;//elemento grande
  TPZStack<TPZCompElSide> elemset;//elementos pequenos
  int numsides = Reference()->NSides();
  int side,nelem,iel;

  for(side = 0; side < numsides; side++) {

    TPZCompElSide thisside(this,side);
    if(mode == EInsert && side >= NCornerConnects()) {//modo insercao
      PZError << "RemoveSideRestraintsII with mode insert should not be called\n";
      return;
      elemset.Resize(0);
      thisside.EqualLevelElementList(elemset,1,0);
      nelem = elemset.NElements();
      large = thisside.LowerLevelElementList(1);
      if(!nelem && large.Exists()) {//existe grande e nao existe igual
	elemset.Resize(0);
	thisside.HigherLevelElementList(elemset,1,1);
//	large.ExpandConnected(elemset,1);
	nelem = elemset.NElements();// se existem pequenos
	for(iel=0; iel<nelem; iel++) {//remover restricoes dos pequenos para o grande
	  TPZInterpolatedElement *cels = (TPZInterpolatedElement *) elemset[iel].Element();
	  //cels->RemoveSideRestraintWithRespectTo(elemset[iel].Side(),large,&thisside);
	  cels->RemoveSideRestraintWithRespectTo(elemset[iel].Side(),large);
	}
	//no construtor � feito a identificacao da ordem : IdentifySideOrder(i)
      }
      continue;
    }//fim mode Einsert

    if(mode == EDelete) {//modo remo��o
      elemset.Resize(0);
      thisside.EqualLevelElementList(elemset,1,0);//iguais
      if(side < NCornerConnects()) thisside.HigherLevelElementList(elemset,1,0);
      nelem = elemset.NElements();
      TPZStack<TPZCompElSide> smallset;
      thisside.HigherLevelElementList(smallset,1,1);//menores
      int nsmall = smallset.NElements();
      large = thisside.LowerLevelElementList(1);
      if(nelem && !large.Exists()) {//iguais e nao grande
	//se existem iguais s� deve-se recalcular a ordem de todos os iguais que restaram
	/* 	int oldorder = SideOrder(side); */
	/* 	int order = ComputeSideOrder(elemset); */
	/* 	if(order != oldorder) { */
	Reference()->ResetReference();
	TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (elemset[0].Element());
	if(cel) {
	  cel->IdentifySideOrder(elemset[0].Side());
	}
	Reference()->SetReference(this);
	/* 	  for(iel=0; iel<nelem; iel++) { */
	/* 	    TPZInterpolatedElement *cel = (TPZInterpolatedElement *) elemset[iel].Element(); */
	/* 	    cel->SetSideOrder(elemset[iel].Side(),order); */
	/* 	  } */
	/* 	  SetSideOrder(side,order); */
	/* 	  if(nsmall) {//existem pequenos : 12c */
	/* 	    thisside.ExpandConnected(smallset,1); */
	/* 	    nsmall = smallset.NElements(); */
	/* 	    int dim; */
	/* 	    for(dim=0; dim<4; dim++) { */
	/* 	      for(iel=0; iel<nsmall; iel++) { */
	/* 		if(smallset[iel].Reference().Dimension() == dim) { */
	/* 		  TPZInterpolatedElement *cel = (TPZInterpolatedElement *) smallset[iel].Element(); */
	/* 		  cel->IdentifySideOrder(smallset[iel].Side()); */
	/* 		} */
	/* 	      } */
	/* 	    } */
	/*  	  } */
	//	}//12c e 1bc
      }
      else if(!nelem) {//nao existem iguais
	if(large.Exists()) {//existe grande : a23 , ab3
	  RemoveSideRestraintWithRespectTo(side,large);
	}//fim large
	if(nsmall) {//existem pequenos : a23
//	  thisside.ExpandConnected(smallset,1);
	  nsmall = smallset.NElements();
	  for(iel=0; iel<nsmall; iel++) {
	    TPZInterpolatedElement *cel = (TPZInterpolatedElement *) smallset[iel].Element();
	    cel->RemoveSideRestraintWithRespectTo(smallset[iel].Side(),thisside);
	  }
	}//ab3
	if(large.Exists() && nsmall) {
	  Reference()->ResetReference();
	  int dim;
	  for(dim=0; dim<4; dim++) {
	    for(iel=0; iel<nsmall; iel++) {
	      if(smallset[iel].Reference().Dimension() == dim) {
		TPZInterpolatedElement *cel = (TPZInterpolatedElement *) smallset[iel].Element();
		TPZInterpolatedElement *largel = (TPZInterpolatedElement *) large.Element();
		cel->RestrainSide(smallset[iel].Side(),largel,large.Side());
	      }
	    }
	  }
	  Reference()->SetReference(this);
	}
      }//fim nelem e large
    }//fim mode EDelete
  }//fim for
}//fim todos

int TPZInterpolatedElement::ComputeSideOrder(TPZVec<TPZCompElSide> &smallset) {
  int nelem = smallset.NElements();
  if(!nelem) {
    cout << "TPZInterpolatedElement::ComputeSideOrder called for empty list, -1 returned\n";
    return -1;
  }
  TPZInterpolatedElement *cel = (TPZInterpolatedElement *) smallset[0].Element();
  int minorder = cel->PreferredSideOrder(smallset[0].Side());
  int iel;
  for(iel=1; iel<nelem; iel++) {
    cel = (TPZInterpolatedElement *) smallset[iel].Element();
    int celorder = cel->PreferredSideOrder(smallset[iel].Side());
    minorder = minorder < celorder ? minorder : celorder;
  }
  return minorder;
}

void TPZInterpolatedElement::InterpolateSolution(TPZInterpolatedElement &coarsel){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t
  TPZTransform t(Dimension());

  //Cedric 16/03/99
  //  Reference()->BuildTransform(NConnects(),coarsel.Reference(),t);
  t = Reference()->BuildTransform2(fReference->NSides()-1,coarsel.Reference(),t);

  int locnod = NConnects();
  int cornod = coarsel.NConnects();
  int locmatsize = NShapeF();
  int cormatsize = coarsel.NShapeF();
  int nvar = fMaterial->NStateVariables();
  int dimension = Dimension();

  TPZFMatrix loclocmat(locmatsize,locmatsize,0.);
  TPZFMatrix projectmat(locmatsize,nvar,0.);

  TPZVec<int> prevorder(dimension),order(dimension);
  TPZIntPoints &intrule = GetIntegrationRule();
  intrule.GetOrder(prevorder);

  TPZVec<int> interpolation(dimension);
  GetInterpolationOrder(interpolation);
  int printing = 0;
  if (printing){
    cout << "locnod = " << locnod << endl;
    cout << "coarnod = " << cornod << endl;
    cout << "local n shape f " << locmatsize << endl;
    cout << "cormatsize " << cormatsize << endl;
    cout << "nvar" << nvar << endl;
    cout << "dimension " << dimension << endl;
    int ccc;
    cout << "interpolation={";
    for (ccc = 0; ccc < interpolation.NElements(); ccc++)
      cout << "  " << interpolation[ccc];
    cout << "  }" << endl;
  }

  // compute the interpolation order of the shapefunctions squared
  int nel = interpolation.NElements();
  int dim,maxorder = interpolation[0];
  for(dim=1;dim<nel;dim++) maxorder = (interpolation[dim] > maxorder) ? interpolation[dim] : maxorder;

  for(dim=0; dim<dimension; dim++) {
    order[dim] = maxorder*2;
  }
  intrule.SetOrder(order);

  TPZFMatrix locphi(locmatsize,1);
  TPZFMatrix locdphi(dimension,locmatsize);	// derivative of the shape function
  // in the master domain

  TPZFMatrix corphi(cormatsize,1);
  TPZFMatrix cordphi(dimension,cormatsize);	// derivative of the shape function
  // in the master domain

  TPZVec<REAL> int_point(dimension),coarse_int_point(dimension);
  TPZFMatrix jacobian(dimension,dimension),jacinv(dimension,dimension);
  TPZFMatrix axes(3,3,0.);
  REAL zero = 0.;
  TPZVec<REAL> x(3,zero);
  TPZVec<REAL> u(nvar);

  int numintpoints = intrule.NPoints();
  REAL weight;
  int lin,ljn,cjn;
  TPZConnect *df;
  TPZBlock &coarseblock = coarsel.Mesh()->Block();

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule.Point(int_ind,int_point,weight);
    REAL jac_det = 1.;
    Reference()->Jacobian( int_point, jacobian , axes,jac_det,jacinv);
    Reference()->X(int_point, x);
    Shape(int_point,locphi,locdphi);
    weight *= jac_det;
    t.Apply(int_point,coarse_int_point);
    coarsel.Shape(coarse_int_point,corphi,cordphi);
    u.Fill(0.);
    int iv = 0;
    for(lin=0; lin<cornod; lin++) {
      df = &coarsel.Connect(lin);
      int dfseq = df->SequenceNumber();

      int dfvar = coarseblock.Size(dfseq);
      int nconshapef = coarsel.NConnectShapeF(lin);
      dfvar = dfvar < nconshapef*nvar ? dfvar : nconshapef;
      for(ljn=0; ljn<dfvar; ljn++) {
	u[iv%nvar] += corphi(iv/nvar,0)*coarseblock(dfseq,0,ljn,0);
	iv++;
      }
      if(dfvar < nconshapef*nvar) iv += nconshapef*nvar- dfvar;
    }
    for(lin=0; lin<locmatsize; lin++) {
      for(ljn=0; ljn<locmatsize; ljn++) {
	loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0);
      }
      for(cjn=0; cjn<nvar; cjn++) {
	projectmat(lin,cjn) += weight*locphi(lin,0)*u[cjn];
      }
    }
    jacobian.Zero();
  }

  if (printing){
    int matiter;
    for (matiter = 0; matiter < loclocmat.Cols();matiter++)
      cout << loclocmat(25,matiter) << "\t";
    cout << endl;
    projectmat.Print("projecmat");
  }

  loclocmat.SolveDirect(projectmat,ELU);
  // identify the non-zero blocks for each row
  TPZBlock &fineblock = Mesh()->Block();
  int iv=0,in;
  for(in=0; in<locnod; in++) {
    df = &Connect(in);
    int dfseq = df->SequenceNumber();
    int dfvar = fineblock.Size(dfseq);
    for(ljn=0; ljn<dfvar; ljn++) {
      fineblock(dfseq,0,ljn,0) = projectmat(iv/nvar,iv%nvar);
      iv++;
    }
  }
  intrule.SetOrder(prevorder);

}

/**calculate the element stiffness matrix*/
void TPZInterpolatedElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef) {

  int i;

  if(fMaterial == NULL){
    PZError << "TPZInterpolatedElement::CalcStiff : no material for this element\n";
    Print(PZError);
    return;
  }
  int numdof = fMaterial->NStateVariables();
  int ncon = NConnects();
  int dim = Dimension();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  // clean ek and ef
  if(!ek.fMat) ek.fMat = new TPZFMatrix();
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ek.fBlock) ek.fBlock = new TPZBlock(ek.fMat);
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  int numeq = nshape*numdof;
  ek.fMat->Redim(numeq,numeq);
  ef.fMat->Redim(numeq,1);
  ek.fBlock->SetNBlocks(ncon);
  ef.fBlock->SetNBlocks(ncon);
  TPZVec<REAL> sol(numdof,0.);
  for (i = 0; i < ncon ; i++)	{
    ek.fBlock->Set(i,NConnectShapeF(i)*numdof);
    ef.fBlock->Set(i,NConnectShapeF(i)*numdof);
  }

  if( !ek.fMat || !ef.fMat || !ek.fBlock || !ef.fBlock){
    cout << "TPZInterpolatedElement.calc_stiff : not enough storage for local stifness"
      " matrix \n";
    Print(cout);
    if(ek.fMat)   delete ek.fMat;
    if(ek.fBlock) delete ek.fBlock;
    if(ef.fMat)   delete ef.fMat;
    if(ef.fBlock) delete ef.fBlock;
    ek.fMat=  NULL;
    ek.fBlock = NULL;
    ef.fMat = NULL;
    ef.fBlock = NULL;
    return;
  }
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);

  for(i=0; i<ncon; ++i){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }
  //suficiente para ordem 5 do cubo
  REAL phistore[220],dphistore[660],dphixstore[660];
  TPZFMatrix phi(nshape,1,phistore,220);
  TPZFMatrix dphi(dim,nshape,dphistore,660),dphix(dim,nshape,dphixstore,660);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL weight = 0.;

  REAL dsolstore[90];
  TPZFMatrix dsol(dim,numdof,dsolstore,90);


  TPZIntPoints &intrule = GetIntegrationRule();
  TPZVec<int> order(3,1);
  intrule.SetOrder(order);
  for(int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind){

    intrule.Point(int_ind,intpoint,weight);

    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);

    fReference->X(intpoint, x);

    weight *= fabs(detjac);

    Shape(intpoint,phi,dphi);

    int ieq;
    switch(dim) {
    case 0:
      break;
    case 1:
      dphix = dphi*(1./detjac);
      break;
    case 2:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
      }
      break;
    case 3:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
	dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
      }
      break;
    default:
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }
    int iv=0,d;

    sol.Fill(0.);
    dsol.Zero();
    for(int in=0; in<ncon; in++) {
      TPZConnect *df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      for(int jn=0; jn<dfvar; jn++) {
	sol[iv%numdof] += phi(iv/numdof,0)*MeshSol(pos+jn,0);
	for(d=0; d<dim; d++)
	  dsol(d,iv%numdof) += dphix(d,iv/numdof)*MeshSol(pos+jn,0);
	iv++;
      }
    }

    cout << "\nCalcStiff sol\n" << sol;
    cout << "\nCalcStiff dsol\n" << dsol;
    cout << "\nCalcStiff phi\n" << phi;
    cout << "\nCalcStiff dphix\n" << dphix;
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phi,dphix,*ek.fMat,*ef.fMat);
    ek.fMat->Print("Correct stiffness matrix");
    ef.fMat->Print("Correct right hand side");
  }
}

void TPZInterpolatedElement::ProjectFlux(TPZElementMatrix &ek, TPZElementMatrix &ef) {

  if(fMaterial == NULL){
    cout << "TPZCompEl1d.project_flux : no material for this element\n";
    Print(cout);
    return;
  }

  int numdof = fMaterial->NStateVariables();//misael
  int num_flux = fMaterial->NFluxes();
  int dim = Dimension();
  int nshape = NShapeF();
  int ncon = NConnects();
  TPZBlock &block = Mesh()->Block();
  TPZIntPoints &intrule = GetIntegrationRule();

  if(!ek.fMat) ek.fMat = new TPZFMatrix();
  if(!ek.fBlock) ek.fBlock = new TPZBlock(ek.fMat);
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  int numeq = nshape;
  ek.fMat->Resize(numeq,numeq);
  ek.fBlock->SetNBlocks(ncon);
  ef.fMat->Resize(numeq,num_flux);
  ef.fBlock->SetNBlocks(ncon);

  if( !ek.fMat || !ef.fMat || !ek.fBlock || !ef.fBlock){
    PZError << "TPZInterpolatedElement.projectflux : not enough storage for local stifne"
      " matrix \n";
    Print(cout);
    if(ek.fMat)   delete ek.fMat;
    if(ek.fBlock) delete ek.fBlock;
    if(ef.fMat)   delete ef.fMat;
    if(ef.fBlock) delete ef.fBlock;
    ek.fMat=  NULL;
    ek.fBlock = NULL;
    ef.fMat = NULL;
    ef.fBlock = NULL;
    return;
  }

  for(int i=0; i<ncon; ++i){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }

  TPZFMatrix phi(nshape,1);
  TPZFMatrix dphi(dim,nshape);
  TPZFMatrix dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3);
  TPZVec<REAL> sol(numdof);
  TPZFMatrix dsol(dim,numdof);
  TPZVec<REAL> flux(num_flux,1);
  TPZVec<REAL> intpoint(dim);
  double weight = 0.;

  for(int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind){

    intrule.Point(int_ind,intpoint,weight);

    fReference->Jacobian( intpoint , jacobian, axes , detjac , jacinv);

    weight *= fabs(detjac);

    fReference->X( intpoint , x);

    Shape(intpoint,phi,dphi);
    int ieq;

    switch(dim) {
    case 0:
      break;
    case 1:
      dphix = dphi*(1./detjac);
      break;
    case 2://Cedric
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
      }
      break;
    case 3:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
	dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
      }
      break;
    default:
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }

    int iv=0,in,jn,d;
    TPZConnect *df;
    for(in=0; in<ncon; in++) {
      df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      for(jn=0; jn<dfvar; jn++) {
	sol[iv%numdof] += phi(iv/numdof,0)*block(dfseq,0,jn,0);
	for(d=0; d<dim; d++) dsol(d,iv%numdof) += dphix(d,iv/numdof)*block(dfseq,0,jn,0);
	iv++;
      }
    }
    fMaterial->Flux(x,sol,dsol,axes,flux);
    for(in=0; in<nshape; in++){
      for(int ifl=0; ifl<num_flux; ifl++){
	(*ef.fMat)(in,ifl) += flux[ifl]*phi(in,0)*weight;
      }
      for(int jn = 0; jn<nshape; jn++){
	(*ek.fMat)(in,jn) += phi(in,0)*phi(jn,0)*weight;
      }
    }
  }
}


/**Implement the refinement of an interpolated element*/
void TPZInterpolatedElement::Divide(int index,TPZVec<int> &sub,int interpolatesolution) {
  if (fMesh->ElementVec()[index] != this) {
    PZError << "TPZInterpolatedElement::Divide index error";
    sub.Resize(0);
    return;
  }
  int nsubelements = fReference->NSubElements();
  sub.Resize(nsubelements);

  TPZManVector<TPZGeoEl *> pv(nsubelements);
  fReference->Divide(pv);//o elemento geometrico correspondente ao atual elemento computacional � dividido
  if(!pv.NElements()) {
    sub.Resize(0);
    return;
  }

  int i;
  RemoveSideRestraintsII(EDelete);//Cedric 25/03/99

  TPZGeoEl *ref;
  TPZInterpolatedElement *cel;
  Reference()->ResetReference();
  int ncon = fReference->NSides();
  TPZInterpolatedElement::gOrder = PreferredSideOrder(ncon-1);
  for (i=0;i<nsubelements;i++)	{
    ref = pv[i];//ponteiro para subelemento i
    ref->CreateCompEl(*fMesh,sub[i]);
    cel = (TPZInterpolatedElement *) fMesh->ElementVec()[sub[i]];
    cel->CheckConstraintConsistency();
    // e' assumido que CreateCompEl inseri o elemento comp no vetor de elementos da malha
  }
  if(interpolatesolution) {
    Mesh()->ExpandSolution();
    for(i=0; i<nsubelements; i++) {
      cel = (TPZInterpolatedElement *) fMesh->ElementVec()[sub[i]];
      cel->CheckConstraintConsistency();
      cel->InterpolateSolution(*this);
      // e' assumido que CreateCompEl inseri o elemento comp no vetor de elementos da malha
    }
  }
  // we assume that the calling program will delete the element from
  // the data structure
  delete this;// deve ser relegado para o Refine
}

REAL TPZInterpolatedElement::CompareElement(int var, char *matname) {
  if (strcmp(matname,Material()->Name())) return(0.);
  REAL error=0.;
  int dim = Dimension();
  int numdof = fMaterial->NSolutionVariables(var);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZManVector<REAL> sol(numdof,0.);
  TPZManVector<REAL> othersol(numdof,0.);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL weight = 0.;

  // At this point we assume grids are identical
  TPZCompEl *otherelement = Reference()->Reference();

  TPZIntPoints &intrule = GetIntegrationRule();
  //GetIntegrationRule().NPoints()
  for(int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind){
    //GetIntegrationRule().Point(int_ind,intpoint,weight);
    intrule.Point(int_ind,intpoint,weight);

    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);

    //fReference->X(intpoint, x);
    weight *= fabs(detjac);
    sol.Fill(0.);
    Solution(intpoint,var,sol);
    otherelement->Solution(intpoint,var,othersol);
    int i=0;

    //Compare the solutions on each point coordinate.
    for (i=0; i<sol.NElements(); i++){
      error += (sol[i]-othersol[i])*(sol[i]-othersol[i]);
    }

  }
  return error;
}

void TPZInterpolatedElement::Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol) {

  if(var >= 100) {
    TPZCompEl::Solution(qsi,var,sol);
    return;
  }
  int nshape = NShapeF();
  int dim = Dimension();
  int ncon = NConnects();
  if(var == 99) {
    sol[0] = SideOrder(ncon-1);
    return;
  }
  TPZBlock &block = fMesh->Block();
  TPZFMatrix &Sol = fMesh->Solution();

  if(fMaterial == NULL){
    PZError << "TPZIntEl::Solution : no Material for this element\n";
    Print(PZError);
    return;
  }

  int numdof = fMaterial->NStateVariables();
  REAL phistore[220],dphistore[660],dphixstore[660];
  TPZFMatrix phi(nshape,1,phistore,220);
  TPZFMatrix dphi(dim,nshape,dphistore,660),dphix(dim,nshape,dphixstore,660);
  TPZManVector<REAL> u(numdof);
  TPZFMatrix du(dim,numdof,0.);
  TPZFMatrix axes(3,3,0.);
  REAL jacstore[10],jacinvstore[10];
  TPZFMatrix jacobian(dim,dim,jacstore,10);
  TPZFMatrix jacinv(dim,dim,jacinvstore,10);
  //TPZVec<REAL> x(3);
  TPZManVector<REAL> x(3);
  REAL detjac;
  int ieq;
  fReference->Jacobian(qsi,jacobian,axes,detjac,jacinv);
  Shape(qsi,phi,dphi);
  fReference->X(qsi,x);
  switch(dim) {
  case 0:
    //dphix.Redim(1,1);
    //dphix(0,0) = dphi(0,0);
    break;
  case 1:
    dphix = dphi*(1./detjac);
    break;
  case 2:
    for(ieq = 0; ieq < nshape; ieq++) {
      dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
      dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
    }
    break;
  case 3:
    for(ieq = 0; ieq < nshape; ieq++) {
      dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
      dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
      dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
    }
    break;

  default:
    PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
    PZError.flush();
  }

  int iv=0,in,jn,d;
  TPZConnect *df;
  u.Fill(0.);
  du.Zero();
  for(in=0; in<ncon; in++) {
    df = &Connect(in);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    for(jn=0; jn<dfvar; jn++) {
      u[iv%numdof] += phi(iv/numdof,0)*Sol(pos+jn,0);
      for(d=0; d<dim; d++){
	du(d,iv%numdof) += dphix(d,iv/numdof)*Sol(pos+jn,0);
      }
      iv++;
    }
  }
  fMaterial->Solution(u,du,axes,var,sol);
}

void TPZInterpolatedElement::Print(ostream &out) {
  out << "Number of connects = " << NConnects() << " Node indexes : ";
  int nod;
  for(nod=0; nod< NConnects(); nod++)
    out << ConnectIndex(nod) <<  '/' << NConnectShapeF(nod) << ' ' ;
  out << endl;
  out << "Side orders = ";
  for (nod=0; nod< NConnects(); nod++) out << SideOrder(nod) << ' ';
  out << endl;
  if (fMaterial) {
    out << "material id " << fMaterial->Id() << endl;
  } else if(Reference()) {
    Reference()->Print(out);
  }
  int id;
  TPZIntPoints &intrule = GetIntegrationRule();
  int dim = Dimension();

  TPZManVector<int> prevorder(dim);
  intrule.GetOrder(prevorder);

  int norders = prevorder.NElements();
  out << "Integration orders : \t";
  for (id=0;id<norders;id++){
    out << prevorder[id] << "\t" ;
  }
  out << endl;
}

void TPZInterpolatedElement::PRefine(int side, int order) {
  SetPreferredOrder(side,order);
  IdentifySideOrder(side);
//   if (side == NConnects()-1){
//     int trueorder = SideOrder(side);
//     SetIntegrationRule(2*trueorder+2);
//   }
}

void TPZInterpolatedElement::PRefine(int order) {
  int side = NCornerConnects();
  int nconnects = NConnects();
  for(;side < nconnects; side++) {
    PRefine(side,order);
  }
  //  int maxorder = SideOrder(nconnects-1);
  //  SetIntegrationRule(2*maxorder);
}

void TPZInterpolatedElement::EvaluateError(
					   void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
					   REAL &true_error,REAL &L2_error,TPZBlock * /*flux */,REAL &estimate) {

  true_error=0.;
  L2_error=0.;
  estimate=0.;
  if(fMaterial == NULL){
    PZError << "TPZInterpolatedElement::EvaluateError : no material for this element\n";
    Print(PZError);
    return;
  }
  // Adjust the order of the integration rule
  TPZIntPoints &intrule = GetIntegrationRule();
  int dim = Dimension();

  TPZManVector<int> prevorder(dim),
    order(dim);
  intrule.GetOrder(prevorder);

  TPZManVector<int> interpolation(0);
  GetInterpolationOrder(interpolation);

  // compute the interpolation order of the shapefunctions squared
  // This is wrong!
  int maxorder = interpolation[0];
  int d;
  for(d=0; d<interpolation.NElements(); d++) {
    maxorder = maxorder < interpolation[d] ? interpolation[d] : maxorder;
  }
  for(d=0; d<dim; d++) {
    order[d] = 2*maxorder+2;
  }
  intrule.SetOrder(order);


  int ndof = fMaterial->NStateVariables();
  int nflux = fMaterial->NFluxes();
  int nshape = NShapeF();
  //suficiente para ordem 5 do cubo
  REAL phistore[220],dphistore[660],dphixstore[660];
  TPZFMatrix phi(nshape,1,phistore,220);
  TPZFMatrix dphi(dim,nshape,dphistore,660),dphix(dim,nshape,dphixstore,660);
  REAL jacobianstore[9],axesstore[9];
  TPZFMatrix jacobian(dim,dim,jacobianstore,9);
  TPZFMatrix axes(3,3,axesstore,9);
  TPZManVector<REAL> x(3);//TPZVec<REAL> x(3,0.);
  REAL duexactstore[90];
  TPZManVector<REAL> u_exact(ndof);
  TPZFMatrix du_exact(dim,ndof,duexactstore,90);
  TPZManVector<REAL> intpoint(3),values(3);
  REAL detjac,weight;
  TPZManVector<REAL> u(ndof);
  REAL dudxstore[90];//,jacinvstore[9];
  TPZFMatrix dudx(dim,ndof,dudxstore,90);
  TPZManVector<REAL> flux_el(nflux);
  TPZMaterial *matp = (TPZMaterial *) fMaterial;
  int ncon = NConnects();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix jacinv(dim,dim);
  int ieq;

  for(int nint=0; nint<GetIntegrationRule().NPoints(); nint++) {

    GetIntegrationRule().Point(nint,intpoint,weight);
    fReference->Jacobian( intpoint , jacobian, axes, detjac , jacinv);
    Shape(intpoint,phi,dphi);
    fReference->X( intpoint , x);
    weight *= fabs(detjac);
    switch(dim) {
    case 0:
      //dphix.Redim(1,1);
      //dphix(0,0) = dphi(0,0);
      break;
    case 1:
      dphix = dphi*(1./detjac);
      break;
    case 2:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
      }
      break;
    case 3:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
	dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
      }
      break;
    default:
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }
    int iv=0,in,jn,d;
    TPZConnect *df;
    u.Fill(0.);
    dudx.Zero();
    for(in=0; in<ncon; in++) {
      df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      for(jn=0; jn<dfvar; jn++) {
	u[iv%ndof] += phi(iv/ndof,0)*block(dfseq,0,jn,0);
	for(d=0; d<dim; d++)
	  dudx(d,iv%ndof) += dphix(d,iv/ndof)*block(dfseq,0,jn,0);
	iv++;
      }
    }//solucao calculculada no sistema local : elementos 2d
    //contribu��es dos erros
    if(fp) {
      fp(x,u_exact,du_exact);
      matp->Errors(x,u,dudx,axes,flux_el,u_exact,du_exact,values);
      true_error += values[0]*weight;
      L2_error += values[1]*weight;
      estimate += values[2]*weight;
    }
  }//fim for : integration rule
   //Norma sobre o elemento
  true_error = sqrt(true_error);
  L2_error = sqrt(L2_error);
  estimate = sqrt(estimate);
  intrule.SetOrder(prevorder);
}
void TPZInterpolatedElement::CalcResidual(TPZElementMatrix &ef) {

  int i;

  if(fMaterial == NULL){
    PZError << "TPZCompEl1d::CalcResidual : no material for this element\n";
    Print(PZError);
    return;
  }
  int numdof = fMaterial->NStateVariables();
  int ncon = NConnects();
  int dim = Dimension();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  // clean ek and ef
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  int numeq = nshape*numdof;
  ef.fMat->Redim(numeq,1);
  ef.fBlock->SetNBlocks(ncon);
  TPZVec<REAL> sol(numdof,0.);
  for (i = 0; i < ncon ; i++)	{
    ef.fBlock->Set(i,NConnectShapeF(i)*numdof);
  }

  if( !ef.fMat || !ef.fBlock){
    cout << "TPZCompEl1d.calc_stiff : not enough storage for local stifness"
      " matrix \n";
    Print(cout);
    if(ef.fMat)   delete ef.fMat;
    if(ef.fBlock) delete ef.fBlock;
    ef.fMat = NULL;
    ef.fBlock = NULL;
    return;
  }
  ef.fConnect.Resize(ncon);

  for(i=0; i<ncon; ++i){
    (ef.fConnect)[i] = ConnectIndex(i);
  }
  //suficiente para ordem 5 do cubo
  REAL phistore[220],dphistore[660],dphixstore[660];
  TPZFMatrix phi(nshape,1,phistore,220);
  TPZFMatrix dphi(dim,nshape,dphistore,660),dphix(dim,nshape,dphixstore,660);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL weight = 0.;

  REAL dsolstore[90];
  TPZFMatrix dsol(dim,numdof,dsolstore,90);


  TPZIntPoints &intrule = GetIntegrationRule();
  for(int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind){
    intrule.Point(int_ind,intpoint,weight);

    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);

    fReference->X(intpoint, x);

    weight *= fabs(detjac);

    Shape(intpoint,phi,dphi);

    int ieq;
    switch(dim) {
    case 0:
      break;
    case 1:
      dphix = dphi*(1./detjac);
      break;
    case 2:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
      }
      break;
    case 3:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
	dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
      }
      break;
    default:
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }
    int iv=0,d;

    sol.Fill(0.);
    dsol.Zero();
    for(int in=0; in<ncon; in++) {
      TPZConnect *df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      for(int jn=0; jn<dfvar; jn++) {
	sol[iv%numdof] += phi(iv/numdof,0)*MeshSol(pos+jn,0);
	for(d=0; d<dim; d++)
	  dsol(d,iv%numdof) += dphix(d,iv/numdof)*MeshSol(pos+jn,0);
	iv++;
      }
    }

    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phi,dphix,*ef.fMat);
  }
}
void TPZInterpolatedElement::CalcBlockDiagonal(TPZStack<int> &connectlist, TPZBlockDiagonal & blockdiag) {

  int i;

  if(fMaterial == NULL){
    PZError << "TPZCompEl1d::CalcDiagonal : no material for this element\n";
    Print(PZError);
    return;
  }
  int ncon = NConnects();
  TPZCompMesh &mesh = *Mesh();
  int numdof = fMaterial->NStateVariables();
  TPZVec<int> dependencyorder;
  connectlist.Resize(0);
  for(i=0; i<ncon; i++) connectlist.Push(ConnectIndex(i));
  BuildConnectList(connectlist);
  BuildDependencyOrder(connectlist,dependencyorder);
  int dim = Dimension();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();

  int numblocks = connectlist.NElements();
  TPZVec<int> blocksizes(numblocks,0);
  int b,c, nexpandedshape = 0,blsize;
  for(b=0; b<numblocks; b++) {
    c = connectlist[b];
    blsize = mesh.ConnectVec()[c].NDof(mesh)/numdof;
    nexpandedshape += blsize;
    blocksizes[b] = blsize*numdof;
  }
  blockdiag.Initialize(blocksizes);
  for(b=0; b<numblocks; b++) {
    blocksizes[b] /= numdof;
  }
  //  int numeq = nshape*numdof;
  TPZVec<REAL> sol(numdof,0.);
  //suficiente para ordem 5 do cubo
  REAL phistore[220],dphistore[660],dphixstore[660];
  TPZFMatrix phi(nexpandedshape,1,phistore,220);
  TPZFMatrix dphi(dim,nexpandedshape,dphistore,660),dphix(dim,nexpandedshape,dphixstore,660);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL weight = 0.;

  REAL dsolstore[90];
  TPZFMatrix dsol(dim,numdof,dsolstore,90);


  TPZIntPoints &intrule = GetIntegrationRule();
  for(int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind){

    intrule.Point(int_ind,intpoint,weight);

    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);

    fReference->X(intpoint, x);

    weight *= fabs(detjac);

    phi.Zero();
    dphi.Zero();
    dphix.Zero();
    Shape(intpoint,phi,dphi);

    int ieq;
    switch(dim) {
    case 0:
      //dphix.Redim(1,1);
      //dphix(0,0) = dphi(0,0);
      break;
    case 1:
      dphix = dphi*(1./detjac);
      break;
    case 2:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
      }
      break;
    case 3:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
	dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
      }
      break;
    default:
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }
    int iv=0,d;

    sol.Fill(0.);
    dsol.Zero();
    for(int in=0; in<ncon; in++) {
      TPZConnect *df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      for(int jn=0; jn<dfvar; jn++) {
	sol[iv%numdof] += phi(iv/numdof,0)*MeshSol(pos+jn,0);
	for(d=0; d<dim; d++)
	  dsol(d,iv%numdof) += dphix(d,iv/numdof)*MeshSol(pos+jn,0);
	iv++;
      }
    }
    // Expand the values of the shape functions and their derivatives
    ExpandShapeFunctions(connectlist,dependencyorder,blocksizes,phi,dphix);
    // Make the contribution in small blocks
    int eq = 0;
    for(b=0; b<numblocks; b++) {
      int cind = connectlist[b];
      int i,j;
      TPZConnect &con = mesh.ConnectVec()[cind];
      blsize = blocksizes[b];
      if(!blsize || con.HasDependency()) {
	eq+= blsize;
	continue;
      }
      TPZFMatrix phil(blsize,1);
      TPZFMatrix dphil(dim,blsize);
      for(i=0; i<blsize; i++) {
	phil(i,0) = phi(eq+i,0);
	for(j=0; j<dim; j++) {
	  dphil(j,i) = dphix(j,eq+i);
	}
      }
      eq += blsize;
      TPZFMatrix ekl(blsize*numdof,blsize*numdof,0.), efl(blsize*numdof,1,0.);
      fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phil,dphil,ekl,efl);
      blockdiag.AddBlock(b,ekl);
    }
  }
}

void TPZInterpolatedElement::ExpandShapeFunctions(TPZVec<int> &connectlist, TPZVec<int> &dependencyorder, TPZVec<int> &blocksizes, TPZFMatrix &phi, TPZFMatrix &dphix) {


  int numblocks =  connectlist.NElements();
  TPZCompMesh &mesh = *Mesh();
  int nhandled=0;
  int current_order = 0;
  int current_block =0;
  while(nhandled < numblocks) {
    if(dependencyorder[current_block] == current_order) {
      nhandled++;
      int cind = connectlist[current_block];
      TPZConnect &con = mesh.ConnectVec()[cind];
      con.ExpandShape(cind,connectlist,blocksizes,phi,dphix);
    }
    current_block++;
    if(current_block == numblocks) {
      current_block = 0;
      current_order++;
    }
  }
}

REAL TPZInterpolatedElement::MeanSolution(int var) {
  int dim = Dimension(), nvars;
  if(!fMaterial) return 0.;

  nvars = fMaterial->NSolutionVariables(var);
  if(nvars!=1) PZError << "TPZInterpolatedElement::MeanSolution is not implemented to nvars != 1." << endl;
  TPZManVector<REAL> sol(nvars,0.);

  int i;
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> intpoint(dim,0.);
  double weight = 0., meanvalue = 0.;
  REAL area = 0.;
  for(i=0;i<GetIntegrationRule().NPoints();i++){
    GetIntegrationRule().Point(i,intpoint,weight);
    fReference->Jacobian(intpoint,jacobian,axes,detjac,jacinv);

    /** Compute the solution value at point integration*/
    Solution(intpoint,var,sol);
	area += weight*fabs(detjac);
    meanvalue += (weight*fabs(detjac)*sol[0]);   //  meanvalue += (weight*fabs(detjac)*sol[j]);
  }
  return (meanvalue/area);
}
/**Compute the contribution to stiffness matrix and load vector on the element*/
void TPZInterpolatedElement::CalcIntegral(TPZElementMatrix &ef) {
  int i;
  if(fMaterial == NULL) return;

  int numdof = fMaterial->NStateVariables();
  int ncon = NConnects();
  int dim = Dimension();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  int numeq = nshape*numdof;
  ef.fMat->Redim(numeq,1);
  ef.fBlock->SetNBlocks(ncon);
  TPZVec<REAL> sol(numdof,0.);
  for(i=0;i<ncon;i++)
    ef.fBlock->Set(i,NConnectShapeF(i)*numdof);

  int in,jn;
  TPZConnect *df;

  if(!ef.fMat || !ef.fBlock){
    if(ef.fMat)   delete ef.fMat;
    if(ef.fBlock) delete ef.fBlock;
    ef.fBlock = NULL;
    ef.fMat = NULL;
    return;
  }
  ef.fConnect.Resize(ncon);

  for(i=0; i<ncon; ++i)
    (ef.fConnect)[i] = ConnectIndex(i);

  TPZFMatrix phi(nshape,1);
  TPZFMatrix dphi(dim,nshape);
  TPZFMatrix dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  double weight = 0.;

  for(int int_ind = 0; int_ind < GetIntegrationRule().NPoints(); ++int_ind){
    GetIntegrationRule().Point(int_ind,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    fReference->X(intpoint, x);
    weight *= fabs(detjac);
    Shape(intpoint,phi,dphi);

    int l, iv=0;
    REAL coef;
    for(in=0;in<numdof;in++) sol[in] = 0.;
    for(in=0; in<ncon; in++) {
      df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      for(jn=0;jn<dfvar;jn++) {
        coef = block(dfseq,0,jn,0);
        sol[iv%numdof] += phi(iv/numdof,0)*coef;
        iv++;
      }
    }
    for(in=0;in<nshape;in++)
      for(l=0;l<numdof;l++)
        (*ef.fMat)(in*numdof+l,0) += weight*phi(in,0)*sol[l];
  }
}

int TPZInterpolatedElement::AdjustPreferredSideOrder(int side, int order) {
  TPZGeoEl *gel = Reference();
  if(!gel) return order;
  int dim = gel->SideDimension(side);
  if(dim != 2) return order;
  TPZStack<int> elsides;
  gel->LowerDimensionSides(side,elsides);
  int maxorder = order;
  int nsides = elsides.NElements();
  int is;
  for(is=0; is<nsides; is++) {
	  TPZGeoElSide gelside(gel,elsides[is]);
    if(gelside.Dimension() != 1) continue;
    int s = elsides[is];
    int sorder = SideOrder(s);
    maxorder = maxorder < sorder ? sorder : maxorder;
  }
  return maxorder;

}


#ifdef _AUTODIFF

/**calculate the element Energy*/
void TPZInterpolatedElement::CalcEnergy(TPZElementMatrix &ek, TPZElementMatrix &ef) {
	cout << "\nCalcEnergy Called\n";
  int i;

  if(fMaterial == NULL){
    PZError << "TPZInterpolatedElement::CalcStiff : no material for this element\n";
    Print(PZError);
    return;
  }
  int numdof = fMaterial->NStateVariables();
  int ncon = NConnects();
  int dim = Dimension();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  TPZFMatrix &MeshSol = Mesh()->Solution();
  // clean ek and ef
  if(!ek.fMat) ek.fMat = new TPZFMatrix();
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ek.fBlock) ek.fBlock = new TPZBlock(ek.fMat);
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  int numeq = nshape*numdof;
  ek.fMat->Redim(numeq,numeq);
  ef.fMat->Redim(numeq,1);
  ek.fBlock->SetNBlocks(ncon);
  ef.fBlock->SetNBlocks(ncon);

  for (i = 0; i < ncon ; i++)	{
    ek.fBlock->Set(i,NConnectShapeF(i)*numdof);
    ef.fBlock->Set(i,NConnectShapeF(i)*numdof);
  }

  if( !ek.fMat || !ef.fMat || !ek.fBlock || !ef.fBlock){
    cout << "TPZInterpolatedElement.calc_stiff : not enough storage for local stifness"
      " matrix \n";
    Print(cout);
    if(ek.fMat)   delete ek.fMat;
    if(ek.fBlock) delete ek.fBlock;
    if(ef.fMat)   delete ef.fMat;
    if(ef.fBlock) delete ef.fBlock;
    ek.fMat=  NULL;
    ek.fBlock = NULL;
    ef.fMat = NULL;
    ef.fBlock = NULL;
    return;
  }
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);

  for(i=0; i<ncon; ++i){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }
  //suficiente para ordem 5 do cubo
  REAL phistore[220],dphistore[660],dphixstore[660];
  TPZFMatrix phi(nshape,1,phistore,220);
  TPZFMatrix dphi(dim,nshape,dphistore,660),dphix(dim,nshape,dphixstore,660);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL weight = 0.;

  TPZVec<FADFADREAL> sol(numdof);
  TPZVec<FADFADREAL> dsol(numdof * dim);// x, y and z data aligned

  FADREAL defaultFAD(numeq, 0., 0.);
  if(defaultFAD.dx(0)==1.)PZError << "\nError: FAD doesn't have default constructor for parameters: (number of derivatives, default value, default derivative value) !";
  FADFADREAL defaultFADFAD(numeq, defaultFAD, defaultFAD);

  FADFADREAL U(defaultFADFAD); // Zeroed Energy Value -> ready for contribution


  TPZIntPoints &intrule = GetIntegrationRule();
  TPZVec<int> order(3,1);
  intrule.SetOrder(order);
  for(int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind){

    intrule.Point(int_ind,intpoint,weight);

    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);

    fReference->X(intpoint, x);

    weight *= fabs(detjac);

    Shape(intpoint,phi,dphi);

    int ieq;
    switch(dim) {
    case 0:
      break;
    case 1:
      dphix = dphi*(1./detjac);
      break;
    case 2:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
      }
      break;
    case 3:
      for(ieq = 0; ieq < nshape; ieq++) {
	dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
	dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
	dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
      }
      break;
    default:
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }
    int iv=0,d;

    sol.Fill(defaultFADFAD);
    dsol.Fill(defaultFADFAD);

    for(int in=0; in<ncon; in++) {
      TPZConnect *df = &Connect(in);
      int dfseq = df->SequenceNumber();
      int dfvar = block.Size(dfseq);
      int pos = block.Position(dfseq);
      for(int jn=0; jn<dfvar; jn++) {
        /*FADFADREAL upos(numeq, iv, FADREAL(numeq, iv, MeshSol(pos+jn,0)));

	sol[iv%numdof] += upos * FADREAL(phi(iv/numdof,0));*/
	// Using direct access to the fad derivatives to enhance performance

	sol[iv%numdof].val().val() += MeshSol(pos+jn,0) * phi(iv/numdof,0);
	sol[iv%numdof].val().fastAccessDx(iv) += phi(iv/numdof,0);
	sol[iv%numdof].fastAccessDx(iv).val() += phi(iv/numdof,0);
	for(d=0; d<dim; d++)
	{
	//dsol[d+(iv%numdof)*dim] += upos * FADREAL (dphix(d, iv/numdof));
	// Using direct access to the fad derivatives to enhance performance

	   dsol[d+(iv%numdof)*dim].val().val() += MeshSol(pos+jn,0) * dphix(d, iv/numdof);
	   dsol[d+(iv%numdof)*dim].val().fastAccessDx(iv) += dphix(d, iv/numdof);
	   dsol[d+(iv%numdof)*dim].fastAccessDx(iv).val() += dphix(d, iv/numdof);
	}

	iv++;
      }
    }
/*
    cout << "\nCalcEnergy sol\n" << sol;
    cout << "\nCalcEnergy dsol\n" << dsol;
    cout << "\nCalcEnergy phi\n" << phi;
    cout << "\nCalcEnergy dphix\n" << dphix;
*/

    fMaterial->ContributeEnergy(x,sol,dsol,U,weight);
    FADToMatrix(U, *ek.fMat, *ef.fMat);
    ef.fMat->Print("Right hand side");
  }

  FADToMatrix(U, *ek.fMat, *ef.fMat);
}

void TPZInterpolatedElement::FADToMatrix(FADFADREAL &U, TPZFMatrix & ek, TPZFMatrix & ef)
{

  int efsz = ef.Rows();
  int ekrows = ek.Rows();
  int ekcols = ek.Cols();

  int Ucols = U.size();
  int Urows = U.val().size();

  if(efsz != Urows)PZError << "Energy Fad type and ef vectors are of different sizes\n";
  if(ekrows != Urows || ekcols != Ucols)PZError << "Energy Fad type and ek matrix are of different sizes\n";

  FADREAL * pBufferFAD;
  int i,j;
  for(j = 0; j < Urows; j++)
  {
     pBufferFAD = & U.fastAccessDx(j);
     ef(j,0) = - pBufferFAD->val();
     // U.val().fastAccessDx(i); must be the same as U.fastAccessDx(i).val();
     for(i = 0; i < Ucols; i++)
     {
        ek(i,j) = pBufferFAD->fastAccessDx(i);
     }
  }
}

#endif
