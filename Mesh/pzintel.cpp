/**
 * @file
 * @brief Contains the implementation of the TPZInterpolatedElement methods.
 */

#include "pzintel.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pztrnsform.h"
#include "pztransfer.h"
#include "TPZMatrixSolver.h"
#include "pzstepsolver.h"
#include "pzquad.h"
#include "TPZElementMatrixT.h"
#include "time.h"
#include "pzmanvector.h"
#include "pzblockdiag.h"
#include "pzcheckrestraint.h"

#include "pzcheckmesh.h"
#include "TPZCompElDisc.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzinterpolatedelement");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#else
static int logger;
#endif


#include "fadType.h"

#include <sstream>
using namespace std;

TPZInterpolatedElement::TPZInterpolatedElement(TPZCompMesh &mesh, TPZGeoEl *reference) :
TPZInterpolationSpace(mesh, reference) {
}

TPZInterpolatedElement::TPZInterpolatedElement(TPZCompMesh &mesh, const TPZInterpolatedElement &copy) :
TPZInterpolationSpace(mesh, copy) {
}

TPZInterpolatedElement::TPZInterpolatedElement(TPZCompMesh &mesh,
        const TPZInterpolatedElement &copy,
        std::map<int64_t, int64_t> & gl2lcElMap) :
TPZInterpolationSpace(mesh, copy, gl2lcElMap) {
}

TPZInterpolatedElement::TPZInterpolatedElement() :
TPZInterpolationSpace() {
}

TPZInterpolatedElement::~TPZInterpolatedElement() {
}

int TPZInterpolatedElement::NShapeF() const {
    int nn = NConnects();
    int in, res = 0;
    for (in = 0; in < nn; in++) {
        TPZConnect &c = Connect(in);
        res += c.NShape();
    }
    return res;
}

int TPZInterpolatedElement::NSideShapeF(int side) const {
    int ns = NSideConnects(side);
    int in, res = 0;
    for (in = 0; in < ns; in++) {
        TPZConnect &c = Connect(SideConnectLocId(in, side));
        res += NConnectShapeF(SideConnectLocId(in, side), c.Order());
    }
    return res;
}

int TPZInterpolatedElement::MidSideConnectLocId(int side) const {
    int il = 1;
    if (NSideConnects(side) == 0) return -1;
    int nodloc = SideConnectLocId(NSideConnects(side) - il, side);
    return nodloc;
}

/**
 * returns a reference to the connect in the middle of the side
 * @param is side which is being queried
 */
TPZConnect &TPZInterpolatedElement::MidSideConnect(int is) const {
    if (NSideConnects(is) == 0) {
        DebugStop();
    }
    return Connect(MidSideConnectLocId(is));
}

int64_t TPZInterpolatedElement::SideConnectIndex(int connect, int side) const {
    return ConnectIndex(SideConnectLocId(connect, side));
}

TPZConnect &TPZInterpolatedElement::SideConnect(int connect, int side) {
    if (side < 0 || connect < 0 || side > Reference()->NSides()) {
        LOGPZ_ERROR(logger, "Exiting SideConnect - has bad first or second parameter.");
        static TPZConnect con;
        return con;
    }
    return (fMesh->ConnectVec()[SideConnectIndex(connect, side)]);
}

void TPZInterpolatedElement::ForceSideOrder(int side, int order) {
    if ((side < Reference()->NCornerNodes() && order != 0) || side >= Reference()->NSides()) {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " setting an order for a corner side " << side << " order " << order;
        LOGPZ_ERROR(logger, sout.str())
        DebugStop();
    }
    TPZCompElSide thisside(this, side);
    TPZCompElSide large = thisside.LowerLevelElementList(1);
    if (large.Exists()) {
#ifdef PZ_LOG
        LOGPZ_INFO(logger, "Exiting ForceSideOrder - large exists.");
#endif
        return;
    }
    TPZConnect &c = MidSideConnect(side);
    int sideorder = c.Order();
    int neworder = order;
    int orderchanged = (sideorder != neworder);
    if (!orderchanged) return;
    TPZStack<TPZCompElSide> elvec;
    thisside.EqualLevelElementList(elvec, 1, 0);
    elvec.Push(thisside);
    int64_t cap, il;
    TPZInterpolatedElement *equal;
    int equalside;
    if (orderchanged == 1) {
        //		elvec.Push(thisside);
        cap = elvec.NElements();
        for (il = 0; il < cap; il++) {
            equal = dynamic_cast<TPZInterpolatedElement *> (elvec[il].Element());
            if (!equal) continue;
            equalside = elvec[il].Side();
            int cindex = equal->MidSideConnectLocId(equalside);
            if (equal->ConnectIndex(cindex) != -1) {
                equal->SetSideOrder(equalside, neworder);
            }
        }
        // Put the accumulated higher dimension sides in highdim (neighbours will not be included because of the third parameter)
        // The higher dimension sides are only analysed for one dimensional sides
        TPZStack<TPZCompElSide> highdim;
        if (thisside.Reference().Dimension() == 1) {
            for (il = 0; il < cap; il++) {
                elvec[il].HigherDimensionElementList(highdim, 1, 1);
            }
        }

        // reuse elvec to put the smaller elements (element/sides restrained by the current element side

        elvec.Resize(0);
        // Adapt the restraints of all smaller elements connected to the current side
        thisside.HigherLevelElementList(elvec, 1, 1);

        // analyse their side order because they are restrained by me
        cap = elvec.NElements();
        TPZInterpolatedElement *smalll;
        int smallside;
        int dimension;
        for (dimension = 0; dimension < 4; dimension++) {
            for (il = 0; il < cap; il++) {
                if (elvec[il].Reference().Dimension() != dimension) continue;
                smalll = dynamic_cast<TPZInterpolatedElement *> (elvec[il].Element());
                if (!smalll) continue;
                smallside = elvec[il].Side();
                // Identify Side Order is used because it will call itself recursively
                // for its smaller elements too.
                smalll->IdentifySideOrder(smallside);
            }
        }

        // Loop over all higher dimension sides
        cap = highdim.NElements();
        for (il = 0; il < cap; il++) {
            // verify if the higher dimension element/side is restrained.
            // if it is restrained then its order will not have changed
            TPZCompElSide highlarge = highdim[il].LowerLevelElementList(1);
            if (highlarge.Exists()) continue;

            // verify if the order of the higher dimension side has changed
            TPZInterpolatedElement *el;
            int highside;
            el = dynamic_cast<TPZInterpolatedElement *> (highdim[il].Element());
            if (!el) continue;
            highside = highdim[il].Side();
            int order, comporder;
            TPZStack<TPZCompElSide> equallist;
            highdim[il].EqualLevelElementList(equallist, 1, 0);
            equallist.Push(highdim[il]);
            // when the  element highdim is restricted by the same element as the original side, do nothing
            // the restriction will be taken care of in the future
            // verify if the order of the higher dimension side changed due to the change in order of the side being studied
            order = el->MidSideConnect(highside).Order();
            comporder = el->ComputeSideOrder(equallist);
            if (order != comporder) {
                // the order has changed
                el->IdentifySideOrder(highside);
            } else {
                // the order hasnt changed
                el->RecomputeRestraints(highside);
            }
        }
    }
}

void TPZInterpolatedElement::IdentifySideOrder(int side) {
    if (MidSideConnectLocId(side) == -1) {
        // there is no connect to adjust the order for
        return;
    }
    TPZCompElSide thisside(this, side);
    TPZCompElSide large = thisside.LowerLevelElementList(1);
    //	int sideorder = MidSideConnect(side).Order();
    int dimension = Reference()->SideDimension(side);

    int neworder;
    int orderchanged = 0;
    TPZStack<TPZCompElSide> elvecall, elvecequal;
    elvecall.Push(thisside);
    thisside.EqualLevelElementList(elvecall, 1, 0);
    int index = MidSideConnectLocId(side);
    int64_t sideconnectindex = ConnectIndex(index);
    int i;
    for (i = 0; i < elvecall.NElements(); i++) {
        int64_t elvecconnectindex = elvecall[i].ConnectIndex();
        // skip the element/sides which have a different connect than my own
        if (sideconnectindex != -1 && sideconnectindex != elvecconnectindex) {
            continue;
        }
        if (elvecall[i].ConnectIndex() != -1) elvecequal.Push(elvecall[i]);
    }

    TPZInterpolatedElement *equal;
    int equalside;
    if (large.Exists()) {
        // There is a larger element
        // identify the order of the larger element and set the interpolation order
        // of the current element and all its neighbours of equal level to this order
        TPZInterpolatedElement *largel = dynamic_cast<TPZInterpolatedElement *> (large.Element());
        neworder = largel->EffectiveSideOrder(large.Side());
        // improve the check on the connect dependencies
        // We assume the datastructure of the elements is consistent in the sense
        // that if the current element has the same side order than the large
        //  element, then all its neighbours will also have the same order
        if (largel->MidSideConnectLocId(large.Side()) != -1 && VerifyConstraintConsistency(side, large) == false) {
            orderchanged = 1;
            RemoveSideRestraintWithRespectTo(side, large);
            SetSideOrder(side, neworder);
            RestrainSide(side, largel, large.Side());
            int64_t cap = elvecequal.NElements();
            for (int64_t il = 0; il < cap; il++) {
                equal = dynamic_cast<TPZInterpolatedElement *> (elvecequal[il].Element());
                if (!equal) continue;
                equalside = elvecequal[il].Side();
                //if(equal->ConnectIndex(equalside) != -1)
                if (equal->MidSideConnectLocId(equalside) != -1) {
                    equal->SetSideOrder(equalside, neworder);
                }
            }
#ifdef PZDEBUG
            if (VerifyConstraintConsistency(side, large) == false) {
                std::cout << "I don't understand\n";
                VerifyConstraintConsistency(side, large);
            }
#endif
        }
    } else {
        // There is no larger element connected to the side
        // identify the new side order by comparing the orders of the equal level elements
        if (MidSideConnect(side).HasDependency()) {
            TPZConnect &c = MidSideConnect(side);
            c.Print(*Mesh());
            LOGPZ_WARN(logger, "TPZInterpolatedElement SetSideOrder fodeu");
            DebugStop();
            large = thisside.LowerLevelElementList(1);
        }
        /// comparte the fPreferred Orders of the neighbouring elements
        neworder = ComputeSideOrder(elvecequal);
        // Verify is the side order of all elements is equal to neworder
        int64_t connectindex = ConnectIndex(MidSideConnectLocId(side));
        int connectorder = MidSideConnect(side).Order();
        if (neworder != connectorder) {
            orderchanged = 1;
            SetSideOrder(side, neworder);
        }
#ifdef PZDEBUG
        int64_t cap = elvecequal.NElements();
        int64_t il = 0;
        while (il < cap) {//SideOrder(int side)
            equal = dynamic_cast<TPZInterpolatedElement *> (elvecequal[il].Element());
            equalside = elvecequal[il].Side();
            TPZConnect connect = equal->Connect(equal->MidSideConnectLocId(equalside));

            int64_t equalindex = equal->ConnectIndex(equal->MidSideConnectLocId(equalside));
            if (equalindex != connectindex) {

                if (connect.LagrangeMultiplier() == 1) {
                    il++;
                    continue;
                }

                DebugStop();
            }
            il++;
        }
#endif

    }
    /// if the dimension of the side = 0 there will be no smaller connects depending on it
    if (dimension == 0) {
        return;
    }

    // reuse elvec to put the smaller elements (element/sides restrained by the current element side
    int sideorder = EffectiveSideOrder(side);
    TPZStack<TPZCompElSide> elvechighlevel;
    // Adapt the restraints of all smaller elements connected to the current side
    thisside.HigherLevelElementList(elvechighlevel, 1, 1);
    for (int dim = 0; dim < 4; dim++) {
        for (int64_t il = 0; il < elvechighlevel.size(); il++) {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*> (elvechighlevel[il].Element());
            TPZGeoEl *ref = intel->Reference();
            if (ref->SideDimension(elvechighlevel[il].Side()) != dim) {
                continue;
            }
            if (intel->VerifyConstraintConsistency(elvechighlevel[il].Side(), thisside) == false) {
                orderchanged = 1;
                intel->RemoveSideRestraintWithRespectTo(elvechighlevel[il].Side(), thisside);
                intel->SetSideOrder(elvechighlevel[il].Side(), sideorder);
                intel->RestrainSide(elvechighlevel[il].Side(), this, side);
                intel->IdentifySideOrder(elvechighlevel[il].Side());
#ifdef PZDEBUG
                if (intel->VerifyConstraintConsistency(elvechighlevel[il].Side(), thisside) == false) {
                    std::cout << "Verify\n";
                    intel->VerifyConstraintConsistency(elvechighlevel[il].Side(),thisside);
                }
#endif
            }
        }
    }
    /// Now we will look at the smaller connects
    if (orderchanged == 1) {//Cedric :  || neworder != computorder
        // The order of the current element changed
        // Therefore the constraints of all smaller elements connected to the current
        // element will need to be adapted
        // The following loop will happen twice, but doesn't affect anything

        //    elvec.Resize(0);
        //    thisside.EqualLevelElementList(elvec,1,0);
        // Put the accumulated higher dimension sides in highdim (neighbours will not be included because of the third parameter)
        // The higher dimension sides are only analysed for one dimensional sides
        TPZStack<TPZCompElSide> highdim;
        {
            int64_t cap = elvecequal.size();
            if (thisside.Reference().Dimension() == 1) {
                for (int64_t il = 0; il < cap; il++) {
                    elvecequal[il].HigherDimensionElementList(highdim, 1, 1);
                }
            }
        }

        for (int64_t il = 0; il < highdim.size(); il++) {

            // verify if the higher dimension element/side is restrained.
            // if it is restrained then its order will not have changed
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (highdim[il].Element());
            if (!intel || intel->ConnectIndex(intel->MidSideConnectLocId(highdim[il].Side())) == -1) {
                continue;
            }
            TPZCompElSide highlarge = highdim[il].LowerLevelElementList(1);
            if (highlarge.Exists()) {
                TPZInterpolatedElement *highlargel = dynamic_cast<TPZInterpolatedElement *> (highlarge.Element());
                if (!highlargel) {
                    DebugStop();
                }
                int sideorder = highlargel->EffectiveSideOrder(highlarge.Side());
                if (intel->VerifyConstraintConsistency(highdim[il].Side(), highlarge) == false) {
                    intel->RemoveSideRestraintWithRespectTo(highdim[il].Side(), highlarge);
                    intel->SetSideOrder(highdim[il].Side(), sideorder);
                    intel->RestrainSide(highdim[il].Side(), highlargel, highlarge.Side());
                    intel->IdentifySideOrder(highdim[il].Side());
                    if (intel->VerifyConstraintConsistency(highdim[il].Side(), highlarge) == false) {
                        std::cout << "Verify\n";
                        intel->VerifyConstraintConsistency(highdim[il].Side(), highlarge);
                    }
                }
            } else {
                intel->IdentifySideOrder(highdim[il].Side());
            }
        }
    }
}

void TPZInterpolatedElement::RecomputeRestraints(int side) {
    TPZStack<TPZCompElSide> elvec;
    elvec.Resize(0);
    TPZCompElSide thisside(this, side);
    //Adapt the restraints of all smaller elements connected to the current side
    thisside.HigherLevelElementList(elvec, 1, 1);
    //thisside.ExpandConnected(elvec,1);
    int64_t cap, il;
    cap = elvec.NElements();
    TPZInterpolatedElement *smalll;
    int order = EffectiveSideOrder(side);
    int smallside;
    int dimension;
    for (dimension = 0; dimension < 4; dimension++) {
        for (il = 0; il < cap; il++) {
            if (elvec[il].Reference().Dimension() != dimension) continue;
            smalll = dynamic_cast<TPZInterpolatedElement *> (elvec[il].Element());
            if (!smalll) continue;
            smallside = elvec[il].Side();
            // Identify Side Order is used because it will call itself recursively
            // for its smaller elements too.
            smalll->RemoveSideRestraintWithRespectTo(smallside, thisside);
            smalll->SetSideOrder(smallside, order);
            smalll->RestrainSide(smallside, this, side);
        }
    }
}

void TPZInterpolatedElement::UpdateNeighbourSideOrder(int side, TPZVec<TPZCompElSide> &elvec) {
    int orside = MidSideConnect(side).Order();
    int neighbourside;
    TPZInterpolatedElement *neighbour;
    int64_t il;
    int64_t cap = elvec.NElements();
    for (il = 0; il < cap; il++) {
        neighbour = dynamic_cast<TPZInterpolatedElement *> (elvec[il].Element());
        if (!neighbour) continue;
        neighbourside = elvec[il].Side();
        int neighord = neighbour->MidSideConnect(neighbourside).Order();
        if (orside != neighord) neighbour->IdentifySideOrder(neighbourside);
    }
}

void TPZInterpolatedElement::BuildTransferMatrix(TPZInterpolatedElement &coarsel, TPZTransform<> &t, TPZTransfer<STATE> &transfer) {
    // accumulates the transfer coefficients between the current element and the
    // coarse element into the transfer matrix, using the transformation t
    TPZGeoEl *ref = Reference();
    int locnod = NConnects();
    int cornod = coarsel.NConnects();
    int locmatsize = NShapeF();
    int cormatsize = coarsel.NShapeF();

    // compare interpolation orders
    // the minimum interpolation order of this needs to be larger than the maximum interpolation order of coarse

    int myminorder = MidSideConnect(locnod - 1).Order();
    int ncorners = NCornerConnects();
    int ic;
    for (ic = ncorners; ic < locnod; ic++) {
        if (MidSideConnect(ic).Order() < myminorder) myminorder = MidSideConnect(ic).Order();
    }
    ncorners = coarsel.NCornerConnects();
    int coarsemaxorder = 0;
    for (ic = ncorners; ic < cornod; ic++) {
        if (coarsel.MidSideConnect(ic).Order() > coarsemaxorder) coarsemaxorder = coarsel.MidSideConnect(ic).Order();
    }
    if (coarsemaxorder > myminorder) {
        stringstream sout;
        sout << "Exiting BuildTransferMatrix - compute the transfer matrix coarse "
                << coarsemaxorder << " me " << myminorder << endl;
        LOGPZ_ERROR(logger, sout.str());
        return;
    }
    TPZStack<int64_t> connectlistcoarse;
    TPZStack<int> dependencyordercoarse, corblocksize;
    connectlistcoarse.Resize(0);
    dependencyordercoarse.Resize(0);
    corblocksize.Resize(0);
    for (ic = 0; ic < cornod; ic++) connectlistcoarse.Push(coarsel.ConnectIndex(ic));
    coarsel.BuildConnectList(connectlistcoarse);
    TPZConnect::BuildDependencyOrder(connectlistcoarse, dependencyordercoarse, *coarsel.Mesh());

    // cornod = number of connects associated with the coarse element
    cornod = connectlistcoarse.NElements();
    int nvar = coarsel.Material()->NStateVariables();

    // number of blocks is cornod
    TPZBlock corblock(0, cornod);
    int in;
    // corblock and corblocksize have the same function
    for (in = 0; in < NConnects(); in++) {
        TPZConnect &c = Connect(in);
        int blsize = c.NShape();
        corblock.Set(in, blsize);
        corblocksize.Push(blsize);
    }

    // cormatsize is , so far, the number of shape functions of the coarse element
    int c;
    for (; in < cornod; in++) {
        c = connectlistcoarse[in];
        int nshape = coarsel.Mesh()->ConnectVec()[c].NShape();
        int blsize = coarsel.Mesh()->ConnectVec()[c].NDof(*(coarsel.Mesh())) / nvar;
        if (nshape != blsize) {
            DebugStop();
        }
        corblock.Set(in, blsize);
        corblocksize.Push(blsize);
        cormatsize += blsize;
    }
    corblock.Resequence();

    //  REAL loclocmatstore[500] = {0.};
    // loclocmat is the inner product of the shape functions of the local element
    // loccormat is the inner product of the shape functions with the shape functions
    //    of the coarse element, both dependent and independent
    TPZFNMatrix<500, STATE> loclocmat(locmatsize, locmatsize);
    TPZFNMatrix<500, STATE> loccormat(locmatsize, cormatsize);
    loclocmat.Zero();
    loccormat.Zero();

    TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
    int dimension = Dimension();

    TPZManVector<int> prevorder(dimension), order(dimension);
    intrule->GetOrder(prevorder);

    TPZManVector<int> interpolation(dimension);
    GetInterpolationOrder(interpolation);

    // compute the interpolation order of the shapefunctions squared
    int dim;
    int maxorder = interpolation[0];
    for (dim = 0; dim < interpolation.NElements(); dim++) {
        maxorder = interpolation[dim] < maxorder ? maxorder : interpolation[dim];
    }
    for (dim = 0; dim < dimension; dim++) {
        order[dim] = maxorder * 2 + 2;
    }
    intrule->SetOrder(order);

    TPZBlock locblock(0, locnod);

    for (in = 0; in < locnod; in++) {
        TPZConnect &c = Connect(in);
        locblock.Set(in, c.NShape());
    }
    locblock.Resequence();

    REAL locphistore[50] = {0.}, locdphistore[150] = {0.};
    TPZFMatrix<REAL> locphi(locmatsize, 1, locphistore, 50);
    TPZFMatrix<REAL> locdphi(dimension, locmatsize, locdphistore, 150);
    locphi.Zero();
    locdphi.Zero();
    // derivative of the shape function
    // in the master domain

    TPZFMatrix<REAL> corphi(cormatsize, 1);
    TPZFMatrix<REAL> cordphi(dimension, cormatsize);
    // derivative of the shape function
    // in the master domain

    REAL jacobianstore[9],
            axesstore[9];
    TPZManVector<REAL> int_point(dimension),
            coarse_int_point(dimension);
    TPZFMatrix<REAL> jacobian(dimension, dimension, jacobianstore, 9), jacinv(dimension, dimension);
    TPZFMatrix<REAL> axes(3, 3, axesstore, 9);
    TPZManVector<REAL> x(3);
    int_point.Fill(0., 0);
    REAL jac_det = 1.;
    ref->Jacobian(int_point, jacobian, axes, jac_det, jacinv);
    REAL multiplier = 1. / jac_det;

    int numintpoints = intrule->NPoints();
    REAL weight;
    int lin, ljn, cjn;

    for (int int_ind = 0; int_ind < numintpoints; ++int_ind) {

        intrule->Point(int_ind, int_point, weight);
        ref->Jacobian(int_point, jacobian, axes, jac_det, jacinv);
        ref->X(int_point, x);
        Shape(int_point, locphi, locdphi);
        weight *= jac_det;
        t.Apply(int_point, coarse_int_point);
        corphi.Zero();
        cordphi.Zero();
        coarsel.Shape(coarse_int_point, corphi, cordphi);

        coarsel.ExpandShapeFunctions(connectlistcoarse, dependencyordercoarse, corblocksize, corphi, cordphi);

        for (lin = 0; lin < locmatsize; lin++) {
            for (ljn = 0; ljn < locmatsize; ljn++) {
                loclocmat(lin, ljn) += weight * locphi(lin, 0) * locphi(ljn, 0) * multiplier;
            }
            for (cjn = 0; cjn < cormatsize; cjn++) {
                loccormat(lin, cjn) += weight * locphi(lin, 0) * corphi(cjn, 0) * multiplier;
            }
        }
        jacobian.Zero();
    }
    loclocmat.SolveDirect(loccormat, ELDLt);


    for (in = 0; in < locnod; in++) {
        //    int cind = connectlistcoarse[in];
        if (Connect(in).HasDependency()) continue;
        int64_t locblocknumber = Connect(in).SequenceNumber();
        int locblocksize = locblock.Size(in);
        int64_t locblockpos = locblock.Position(in);
        TPZStack<int> locblockvec;
        TPZStack<int> globblockvec;
        int numnonzero = 0, jn;
        //      if(transfer.HasRowDefinition(locblocknumber)) continue;

        for (jn = 0; jn < cornod; jn++) {
            int corblocksize = corblock.Size(jn);
            int64_t corblockpos = corblock.Position(jn);
            int cind = connectlistcoarse[jn];
            TPZConnect &con = coarsel.Mesh()->ConnectVec()[cind];
            if (con.HasDependency()) continue;
            int64_t corblocknumber = con.SequenceNumber();
            if (locblocksize == 0 || corblocksize == 0) continue;
            TPZFMatrix<STATE> smalll(locblocksize, corblocksize, 0.);
            loccormat.GetSub(locblockpos, corblockpos,
                    locblocksize, corblocksize, smalll);
            REAL tol = Norm(smalll);
            if (tol >= 1.e-10) {
                locblockvec.Push(jn);
                globblockvec.Push(corblocknumber);
                numnonzero++;
            }
        }
        if (transfer.HasRowDefinition(locblocknumber)) continue;
        transfer.AddBlockNumbers(locblocknumber, globblockvec);
        int jnn;
        for (jnn = 0; jnn < numnonzero; jnn++) {
            jn = locblockvec[jnn];
            int corblocksize = corblock.Size(jn);
            int64_t corblockpos = corblock.Position(jn);
            if (corblocksize == 0 || locblocksize == 0) continue;
            TPZFMatrix<STATE> smalll(locblocksize, corblocksize, 0.);
            loccormat.GetSub(locblockpos, corblockpos, locblocksize, corblocksize, smalll);
            transfer.SetBlockMatrix(locblocknumber, globblockvec[jnn], smalll);
        }
    }
    intrule->SetOrder(prevorder);
}

int64_t TPZInterpolatedElement::CreateMidSideConnect(int side) {
    TPZCompMesh *cmesh = Mesh();
    TPZMaterial * mat = Material();
#ifdef PZDEBUG
    if (!mat) {
        std::cout << __PRETTY_FUNCTION__ << " no material associated with matid " << Reference()->MaterialId() << std::endl;
    }
#endif
    int nvar = 1;
    if (mat) nvar = mat->NStateVariables();
    int64_t newnodeindex;
    int64_t il;
    int64_t nodloc = MidSideConnectLocId(side);


    TPZStack<TPZCompElSide> elvec;
    TPZCompElSide thisside(this, side);
    if (side < NCornerConnects()) {
        thisside.EqualLevelElementList(elvec, 1, 0);
        int64_t nel = elvec.NElements(); // (1)
        if (nel && elvec[nel - 1].Reference().Dimension() == thisside.Reference().Dimension()) {
            newnodeindex = elvec[nel - 1].ConnectIndex();
            SetConnectIndex(nodloc, newnodeindex);
        } else {
            // corner nodes have order 1
            int order = 1; //cmesh->GetDefaultOrder();
            int nshape = NConnectShapeF(nodloc, order);
            // create the connects with order 1
            newnodeindex = cmesh->AllocateNewConnect(nshape, nvar, order);
            TPZConnect &newnod = cmesh->ConnectVec()[newnodeindex];
            if (newnod.HasDependency()) {
                LOGPZ_ERROR(logger, "CreateMidSideConnect, new node has dependency");
                newnod.Print(*cmesh);
                DebugStop();
            }
            int64_t seqnum = newnod.SequenceNumber();
            SetConnectIndex(nodloc, newnodeindex); //Is true only to one-dimensional case
            cmesh->Block().Set(seqnum, nvar * nshape);
            // We created a new node, check whether the node needs to be constrained
            TPZCompElSide father = thisside.LowerLevelElementList(1);
            if (father.Exists()) {
                int side_neig = father.Side();
                TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (father.Element());
                if (cel) {
                    // corner nodes to not change order
                    //newnod.SetOrder(cel->SideOrder(side_neig));
                    RestrainSide(side, cel, side_neig);
                }
            }
        }
        return newnodeindex;
    }

    // Connect looks for a connecting element of equal or lower level
    TPZInterpolatedElement *cel = 0;
    int side_neig = 0;
    thisside.EqualLevelElementList(elvec, 1, 1);
    int64_t nelem = elvec.NElements();
    // find an element in the list which is interpolated
    if (nelem) {
        cel = dynamic_cast<TPZInterpolatedElement *> (elvec[0].Element());
        side_neig = elvec[0].Side();
    }
    int64_t newnodecreated = 0;
    if (cel) {
        newnodeindex = cel->ConnectIndex(cel->MidSideConnectLocId(side_neig));
        SetConnectIndex(nodloc, newnodeindex);
    } else {
        int nshape = 0;
        int order = cmesh->GetDefaultOrder();
        newnodeindex = cmesh->AllocateNewConnect(nshape, nvar, order);
        TPZConnect &newnod = cmesh->ConnectVec()[newnodeindex];
        int64_t seqnum = newnod.SequenceNumber();
        SetConnectIndex(nodloc, newnodeindex);
        nshape = NConnectShapeF(nodloc, order);
        newnod.SetNShape(nshape);
        cmesh->Block().Set(seqnum, nvar * nshape);
        newnodecreated = 1;
    }
    TPZCompElSide father = thisside.LowerLevelElementList(1);
    if (father.Exists() && newnodecreated) {
        // There is a larger element connected along the side that is being created
        // We will adopt the interpolation order of the father and compute the
        // restraints
        // IT IS ASSUMED THAT ALL NEIGHBOURING ELEMENTS AND SMALLER ELEMENTS CONNECTED
        // TO THIS SIDE HAVE CONSISTENT INTERPOLATION ORDERS
        side_neig = father.Side();
        TPZInterpolatedElement *elfather = dynamic_cast<TPZInterpolatedElement *> (father.Element());
        int sideorder = elfather->EffectiveSideOrder(side_neig);
        SetSideOrder(side, sideorder);
        RestrainSide(side, elfather, side_neig);
        elvec.Resize(0);
        thisside.HigherLevelElementList(elvec, 1, 1);
        //    thisside.ExpandConnected(elvec,1);
        //    thisside.RemoveDuplicates(elvec);
        int64_t cap = elvec.NElements();
        int dim;
        for (dim = 0; dim < 3; dim++) {
            for (il = 0; il < cap; il++) {
                if (elvec[il].Reference().Dimension() != dim) continue;
                cel = dynamic_cast<TPZInterpolatedElement *> (elvec[il].Element());
                Reference()->ResetReference();
                cel->RemoveSideRestraintWithRespectTo(elvec[il].Side(), father);
                Reference()->SetReference(this);

                cel->RestrainSide(elvec[il].Side(), this, side);
            }
        }
    } else if (father.Exists() && !newnodecreated) {
        side_neig = father.Side();
        TPZInterpolatedElement *elfather = dynamic_cast<TPZInterpolatedElement *> (father.Element());
        int sideorder = elfather->EffectiveSideOrder(side_neig);
        SetSideOrder(side, sideorder);
    } else if (!father.Exists() && !newnodecreated) {
        // The element does not have a larger element connected along the side
        // The insertion of the new element may have an effect on the interpolation
        // order of all equal level elements which are connected
        if (Connect(MidSideConnectLocId(side)).HasDependency()) {
            stringstream sout;
            sout << "TPZInterpolatedElement fodeu side " << side << endl;
            LOGPZ_ERROR(logger, sout.str());
            father = thisside.LowerLevelElementList(1);
            DebugStop();
        }
        elvec.Resize(0);
        thisside.EqualLevelElementList(elvec, 1, 0);
        int oldorder = -1;
        if (elvec.NElements()) {
            TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (elvec[0].Element());
            oldorder = cel->MidSideConnect(elvec[0].Side()).Order();
        }
        elvec.Push(thisside);
        int sideorder = ComputeSideOrder(elvec);
        if (sideorder != oldorder) {
            TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (elvec[0].Element());
            cel->IdentifySideOrder(elvec[0].Side());
        } else {
            SetSideOrder(side, sideorder);
        }
    } else if (!father.Exists() && newnodecreated) {
        elvec.Resize(0);
        thisside.HigherLevelElementList(elvec, 1, 1);
        //    thisside.ExpandConnected(elvec,1);
        //    thisside.RemoveDuplicates(elvec);
        int64_t cap = elvec.NElements();
        int sideorder = PreferredSideOrder(side);
        SetSideOrder(side, sideorder);
        sideorder = EffectiveSideOrder(side);
        // We check on all smaller elements connected to the current element
        // whether their interpolation order is consistent with the interpolation
        // order of the current element/side
        int dim;
        for (dim = 0; dim < 3; dim++) {
            for (il = 0; il < cap; il++) {
                if (elvec[il].Reference().Dimension() != dim) continue;
                TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (elvec[il].Element());
                if (!cel) continue;
                // the restraint will only be applied if the first
                // element in Connect is equal to the current element
                int locindex = cel->MidSideConnectLocId(elvec[il].Side());
                if (locindex < 0) {
                    continue;
                }
                int celsideorder = cel->Connect(locindex).Order();
#ifdef PZDEBUG
                {
                    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (elvec[il].Element());
                    TPZGeoEl *gel = Reference();
                    TPZConnect &c = intel->MidSideConnect(elvec[il].Side());
                    if (c.HasDependency()) {
                        int64_t cindex = intel->ConnectIndex(intel->MidSideConnectLocId(elvec[il].Side()));
                        gel->ResetReference();
                        TPZCompElSide locallarge;
                        TPZGeoElSide localgelside(intel->Reference(), elvec[il].Side());
                        locallarge = localgelside.LowerLevelCompElementList2(1);
                        std::cout << "I dont understand\n";
                    }
                }
#endif
                if (celsideorder != sideorder) {
                    cel->SetSideOrder(elvec[il].Side(), sideorder);
                }
                // if the interpolation order remains unchanged, just restrain the side
                cel->RestrainSide(elvec[il].Side(), this, side);
                if (celsideorder != sideorder) {
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
    TPZCompElSide thisside(this, side);
    TPZGeoElSide thisgeoside = thisside.Reference();
    TPZCompElSide largecompside(large, neighbourside);
    TPZGeoElSide largegeoside = largecompside.Reference();
    TPZInterpolatedElement *cel = nullptr;
    if (largecompside.Exists()) cel = dynamic_cast<TPZInterpolatedElement *> (largecompside.Element());
    if (!cel) {
        LOGPZ_ERROR(logger, "Exiting RestrainSide - null computational element.");
        return;
    }
    int locind = MidSideConnectLocId(side);
    if (locind < 0) {
        DebugStop();
    }
    TPZConnect &myconnect = Connect(locind);
    if (myconnect.NShape() == 0) {
        /// no shape functions to restrain
        return;
    }
    if (myconnect.HasDependency() && largegeoside.Dimension() > 0) {
        LOGPZ_WARN(logger, "RestrainSide - unnecessary call to restrainside");
        DebugStop();
    }
    if (cel->ConnectIndex(cel->MidSideConnectLocId(largecompside.Side())) == -1) {
        LOGPZ_ERROR(logger, "Exiting RestrainSide - Side of large element not initialized");
        DebugStop();
        return;
    }
    if (largegeoside.Dimension() < thisgeoside.Dimension()) {
        LOGPZ_ERROR(logger, "Exiting RestrainSide - Trying to restrict to a side of smaller dimension");
        DebugStop();
        return;
    }
    TPZTransform<> t(thisgeoside.Dimension());
    thisgeoside.SideTransform3(largegeoside, t);
    int nsideconnects = NSideConnects(side);
    int maxord = 1;
    for (int sidecon = 0; sidecon < nsideconnects; sidecon++) {
        TPZConnect &c = SideConnect(sidecon, side);
        int sideord = c.Order();
        maxord = maxord < sideord ? sideord : maxord;
    }
    int largeorder = large->EffectiveSideOrder(neighbourside);
    int sideorder = MidSideConnect(side).Order();
    if (sideorder < largeorder && thisgeoside.Dimension() && largegeoside.Dimension()) {
        DebugStop();
    }
    TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(side, maxord * 2);
    if (!intrule) {
        LOGPZ_ERROR(logger, "Exiting RestrainSide - cannot create side integration rule");
        return;
    }
    int numint = intrule->NPoints();
    int numshape = NSideShapeF(side);
    int sidedimension = thisgeoside.Dimension();
    int largesidedimension = largegeoside.Dimension();
    int numshapel = large->NSideShapeF(neighbourside);
    TPZFNMatrix<100, REAL> phis(numshape, 1), dphis(2, numshape), phil(numshapel, 1), dphil(2, numshapel);
    TPZFNMatrix<1000, REAL> MSL(numshape, numshapel, 0.);
    TPZFNMatrix<1000, REAL> *M = new TPZFNMatrix<1000, REAL>(numshape, numshape, 0.);
    TPZManVector<REAL, 3> par(sidedimension), pointl(largesidedimension); //,point(3);
    int64_t in, jn;
    REAL weight;
    for (int it = 0; it < numint; it++) {
        intrule->Point(it, par, weight);
        SideShapeFunction(side, par, phis, dphis);
        t.Apply(par, pointl);
        large->SideShapeFunction(neighbourside, pointl, phil, dphil);
        for (in = 0; in < numshape; in++) {
            for (jn = 0; jn < numshape; jn++) {
                (*M)(in, jn) += phis(in, 0) * phis(jn, 0) * weight;
            }
            for (jn = 0; jn < numshapel; jn++) {
                MSL(in, jn) += phis(in, 0) * phil(jn, 0) * weight;
            }
        }
    }
#ifdef PZ_LOG_keep
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        M->Print("MSS = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZStepSolver<REAL> MSolve(M);
    MSolve.SetDirect(ELU);

    MSolve.Solve(MSL, MSL);

    //MSolve.ResetMatrix();
    int numsidenodes_small = NSideConnects(side);
    // Philippe 12/3/99
    //  int numsidenodes_large = NSideConnects(neighbourside);
    int numsidenodes_large = large->NSideConnects(neighbourside);
    TPZBlock MBlocksmall(0, numsidenodes_small), MBlocklarge(0, numsidenodes_large);
    for (in = 0; in < numsidenodes_small; in++) {
        int locid = SideConnectLocId(in, side);
        TPZConnect &c = Connect(locid);
#ifdef PZDEBUG
        if (NConnectShapeF(locid, c.Order()) != c.NShape()) {
            DebugStop();
        }
#endif
        unsigned int nshape = c.NShape();
        MBlocksmall.Set(in, nshape);
    }
    for (in = 0; in < numsidenodes_large; in++) {
        int locid = large->SideConnectLocId(in, neighbourside);
        TPZConnect &c = large->Connect(locid);
        unsigned int nshape = c.NShape();
#ifdef PZDEBUG
        if (large->NConnectShapeF(locid, c.Order()) != nshape) {
            DebugStop();
        }
#endif
        MBlocklarge.Set(in, nshape);
    }

    MBlocksmall.Resequence();
    MBlocklarge.Resequence();
    TPZFNMatrix<20, REAL> blocknorm(numsidenodes_small, numsidenodes_large, 0.);
    for (in = 0; in < numsidenodes_small; in++) {
        int ibl = MBlocksmall.Size(in);
        if (!ibl) continue;
        for (jn = 0; jn < numsidenodes_large; jn++) {
            int jbl = MBlocklarge.Size(jn);
            if (!jbl) continue;
            int i, j;
            int64_t ipos = MBlocksmall.Position(in);
            int64_t jpos = MBlocklarge.Position(jn);
            for (i = 0; i < ibl; i++) for (j = 0; j < jbl; j++) blocknorm(in, jn) += fabs(MSL(ipos + i, jpos + j)) * fabs(MSL(ipos + i, jpos + j));
            blocknorm(in, jn) /= (ibl * jbl);
            blocknorm(in, jn) = sqrt(blocknorm(in, jn));
        }
    }
#ifdef PZDEBUG
    CheckConstraintConsistency(side);
#endif
    TPZConnect &inod = Connect(MidSideConnectLocId(side));
    int64_t inodindex = ConnectIndex(MidSideConnectLocId(side));
    int64_t ndepend = 0;
    in = numsidenodes_small - 1;
#ifdef PZDEBUG
    if (MBlocksmall.Size(in) == 0) {
        DebugStop();
    }
#endif
    for (jn = 0; jn < numsidenodes_large; jn++) {
        if (MBlocksmall.Size(in) == 0 || MBlocklarge.Size(jn) == 0) {
            continue;
        }
        int64_t jnodindex = large->SideConnectIndex(jn, neighbourside);
        TPZConnect::TPZDepend *depend = inod.AddDependency(inodindex, jnodindex, MSL, MBlocksmall.Position(in), MBlocklarge.Position(jn),
                MBlocksmall.Size(in), MBlocklarge.Size(jn));
        if (blocknorm(in, jn) < 1.e-8) {
            depend->fDepMatrix.Zero();
        }
        ndepend++;
    }

    if (!ndepend) {
        //cout << "Caso esquisito!!! Chame o Boss que vc receberem premio\n";
        for (jn = 0; jn < numsidenodes_large; jn++) {
            int64_t jnodindex = large->SideConnectIndex(jn, neighbourside);
            if (MBlocklarge.Size(jn)) {
                inod.AddDependency(inodindex, jnodindex, MSL, MBlocksmall.Position(in), MBlocklarge.Position(jn),
                        MBlocksmall.Size(in), MBlocklarge.Size(jn));
            }
            ndepend++;
        }
    }
    delete intrule;

#ifdef HUGE_DEBUG
    // a matriz frestraint deveria ser igual a MSL
    TPZCheckRestraint test(thisside, largecompside);
    //test.Print(cout);
    int64_t imsl, jmsl;
    int64_t rmsl = MSL.Rows();
    int64_t cmsl = MSL.Cols();
    int64_t rtest = test.RestraintMatrix().Rows();
    int64_t ctest = test.RestraintMatrix().Cols();

    if (rtest != rmsl || ctest != cmsl) {
        stringstream sout;
        sout << "Exiting - Restraint matrix side incompatibility: MSL (rows,cols): ( " << rmsl
                << " , " << cmsl << " )" << " RestraintMatrix (rows,cols): (" << rtest << " , " << ctest << " )\n"
                << "press any key to continue";
        LOGPZ_ERROR(logger, sout.str());
        int a;
        cin >> a;
        return;
    }

    TPZFMatrix<REAL> mslc(MSL);
    mslc -= test.RestraintMatrix();

    REAL normmsl = 0.;
    for (imsl = 0; imsl < rmsl; imsl++) {
        for (jmsl = 0; jmsl < cmsl; jmsl++) {
            normmsl += sqrt(mslc(imsl, jmsl) * mslc(imsl, jmsl));
        }
    }
    if (normmsl > 1.E-6) {
        stringstream sout;
        sout << "TPZInterpolatedElement::Error::MSL matrix has non zero norm " << normmsl << "\n";
        mslc.Print("Difference Matrix ", sout);
        for (imsl = 0; imsl < rmsl; imsl++) {
            for (jmsl = 0; jmsl < cmsl; jmsl++) {
                if (fabs(MSL(imsl, jmsl) - test.RestraintMatrix()(imsl, jmsl)) > 1.E-6) {
                    sout << "msl[ " << imsl << " , " << jmsl << " ] = " << MSL(imsl, jmsl) << "\t "
                            << test.RestraintMatrix()(imsl, jmsl) << endl;
                }
            }
        }
        LOGPZ_ERROR(logger, sout.str());
        //     int a;
        //     gDebug = 1;
        //     cin >> a;
    }

    // verificar a norma de MSL
    if (test.CheckRestraint()) {
        stringstream sout
        sout << "TPZInterpolatedElement::Error::Bad restraints detected\n"; // recado de erro.
        test.Print(sout);
        //     int a;
        //     gDebug = 1;
        //     cin >> a;
        test.Diagnose();
        LOGPZ_ERROR(logger, sout.str());
        TPZCheckRestraint test2(thisside, largecompside);
    }
#endif
}

void TPZInterpolatedElement::CheckConstraintConsistency() {
    int nc = Reference()->NSides();
    //  int a;
    for (int c = 0; c < nc; c++) CheckConstraintConsistency(c);
}

int TPZInterpolatedElement::CheckElementConsistency() {
    int dimel = Dimension();
    int iside;
    int a;
    int nstate = 1;
    nstate = Material()->NStateVariables();

    for (iside = 0; iside < NConnects(); iside++) {
        int nshape = Connect(iside).NShape();
        if (Connect(iside).CheckDependency(nshape, Mesh(), nstate) == -1) {
            LOGPZ_WARN(logger, "CheckElementConsistency detected inconsistency 1");
            return 0;
        }
        if (Connect(iside).NDof(*Mesh()) != nshape * nstate) {
            LOGPZ_WARN(logger, "CheckElementConsistency detected inconsistency 2");
            return 0;
        }
    }
    for (iside = 0; iside < (Reference()->NSides() - 1); iside++) {
        TPZCompElSide celside(this, iside);
        int dimsmall = celside.Reference().Dimension();

        TPZIntPoints *sirule = Reference()->CreateSideIntegrationRule(iside, EffectiveSideOrder(iside)*2);

        if (dimsmall >= dimel) {
            stringstream sout;
            sout << "TPZInterpolatedElement::CheckConstraintConsistency : dismall >= dimel: "
                    << dimsmall << " >= " << dimel << endl
                    << "press any key to continue";
            LOGPZ_ERROR(logger, sout.str());
            DebugStop();
            delete sirule;
            return 0;
        }
        int idim;
        int nshapes = NSideShapeF(iside);
        TPZFMatrix<REAL> phis(nshapes, 1);
        TPZFMatrix<REAL> dphis(dimsmall, nshapes);
        for (idim = (dimsmall + 1); idim <= dimel; idim++) {
            TPZStack <TPZGeoElSide> geoelsidevec;
            celside.Reference().Element()->AllHigherDimensionSides(iside, idim, geoelsidevec);

            int nelhigh = geoelsidevec.NElements();
            int inh;
            for (inh = 0; inh < nelhigh; inh++) {
                int sidel = geoelsidevec[inh].Side();
                int nshapel = NSideShapeF(sidel);
                TPZFMatrix<REAL> phil(nshapel, 1);
                TPZFMatrix<REAL> dphil(idim, nshapel);
                int npts = sirule->NPoints();
                int ipt;
                TPZTransform<> transform(Reference()->SideToSideTransform(iside, sidel));
                for (ipt = 0; ipt < npts; ipt++) {
                    TPZVec <REAL> pts(dimsmall);
                    TPZVec <REAL> ptl(idim);
                    REAL w;
                    sirule->Point(ipt, pts, w);
                    SideShapeFunction(iside, pts, phis, dphis);
                    transform.Apply(pts, ptl);
                    SideShapeFunction(sidel, ptl, phil, dphil);
                    int check = CompareShapeF(iside, sidel, phis, dphis, phil, dphil, transform);
                    if (!check) {
#ifdef PZ_LOG
                        LOGPZ_INFO(logger, "Exiting CheckElementConsistency - don't compare shapefunctions.");
#endif
                        delete sirule;
                        return check;
                    }
                }
            }
        }
        delete sirule;
    }
    return 1;
}

int TPZInterpolatedElement::CompareShapeF(int sides, int sidel, TPZFMatrix<REAL> &phis, TPZFMatrix<REAL> &dphis, TPZFMatrix<REAL> &phil, TPZFMatrix<REAL> &dphil, TPZTransform<> &transform) {
    int ncons = NSideConnects(sides);
    int nconl = NSideConnects(sidel);
    TPZVec<int> posl(nconl + 1), poss(ncons + 1);
    int icon, icons, iconl;
    if (nconl) posl[0] = 0;
    if (ncons) poss[0] = 0;
    for (icon = 0; icon < nconl; icon++) {
        posl[icon + 1] = posl[icon] + Connect(SideConnectLocId(icon, sidel)).NShape();
    }
    for (icon = 0; icon < ncons; icon++) {
        poss[icon + 1] = poss[icon] + Connect(SideConnectLocId(icon, sides)).NShape();
    }
    REAL diff = 0.;
    for (icons = 0; icons < ncons; icons++) {
        int consind = SideConnectLocId(icons, sides);
        for (iconl = 0; iconl < nconl; iconl++) {
            int conlind = SideConnectLocId(iconl, sidel);
            if (consind == conlind) {
                int nshape = poss[icons + 1] - poss[icons];
                int dims = Reference()->SideDimension(sides);
                int diml = Reference()->SideDimension(sidel);
                int ishape, idim, jdim;
                for (ishape = 0; ishape < nshape; ishape++) {
                    int shapel = posl[iconl] + ishape;
                    int shapes = poss[icons] + ishape;
                    diff += (phis(shapes, 0) - phil(shapel, 0))*(phis(shapes, 0) - phil(shapel, 0));
                    REAL derivcomp[3];
                    for (idim = 0; idim < dims; idim++) {
                        derivcomp[idim] = 0.;
                        for (jdim = 0; jdim < diml; jdim++) {
                            derivcomp[idim] += dphil(jdim, shapel) * transform.Mult()(jdim, idim);
                        }
                        diff += (dphis(idim, shapes) - derivcomp[idim])*(dphis(idim, shapes) - derivcomp[idim]);
                    }
                }
            }
        }
    }

    if (diff >= 1.e-6) {
        stringstream sout;
        sout << "TPZInterpolatedElement::CompareShapeF sides " << sides << " sidel " << sidel << " do not compare diff " << diff << endl;
        LOGPZ_ERROR(logger, sout.str());
    }
    return 1;
}

void TPZInterpolatedElement::CheckConstraintConsistency(int side) {
    TPZCompElSide thisside(this, side);
    TPZCompElSide large = thisside.LowerLevelElementList(1);
    if (large.Exists()) {
        int largeside = large.Side();
        TPZInterpolatedElement *largel = dynamic_cast<TPZInterpolatedElement *> (large.Element());
        int nlargesideconnects = largel->NSideConnects(largeside);
        int n_small_side_connect = NSideConnects(side)-1;
        TPZConnect &thiscon = SideConnect(n_small_side_connect, side);
        int jn;
        TPZConnect::TPZDepend *dep = thiscon.FirstDepend();
        while (dep) {
            for (jn = 0; jn < nlargesideconnects; jn++) {
                int largeindex = largel->SideConnectIndex(jn, largeside);
                if (dep->fDepConnectIndex == largeindex || largeindex == ConnectIndex(side)) break;
            }
            if (jn == nlargesideconnects) {
                stringstream sout;
                sout << "Dependency inconsistency";
                thiscon.Print(*Mesh(), sout);
                largel->Print(sout);
                LOGPZ_ERROR(logger, sout.str());
            }
            dep = dep->fNext;
        }
    } else {
        int n_small_side_connect = NSideConnects(side)-1;
        TPZConnect &thiscon = SideConnect(n_small_side_connect, side);
        if (thiscon.HasDependency()) {
            stringstream sout;
            large = thisside.LowerLevelElementList(1);
            sout << "Dependency inconsistency\n";
            thiscon.Print(*Mesh(), sout);
            LOGPZ_ERROR(logger, sout.str());
        }
    }
}

void TPZInterpolatedElement::RemoveSideRestraintWithRespectTo(int side,
        const TPZCompElSide & neighbour) {
    TPZCompElSide thisside(this, side);
    TPZCompElSide largeset = thisside.LowerLevelElementList(1);
    if (!largeset.Exists()) {
        LOGPZ_WARN(logger, "Exiting RemoveSideRestraintWithRespectTo inconsistent 1");
        return;
    }
    TPZInterpolatedElement *smallfather = dynamic_cast<TPZInterpolatedElement *> (largeset.Element());
    int smallfatherside = largeset.Side();
    if (!neighbour.Reference().NeighbourExists(largeset.Reference())) {
        LOGPZ_WARN(logger, "Exiting RemoveSideRestraintWithRespectTo inconsistent 3");
    }
    int js;
    int nsfather = smallfather->NSideConnects(smallfatherside);
    int dfiindex = ConnectIndex(MidSideConnectLocId(side));
    TPZConnect *dfi = &Connect(MidSideConnectLocId(side));
    for (js = 0; js < nsfather; js++) {
        int dfjindex = smallfather->SideConnectIndex(js, smallfatherside);
        dfi->RemoveDepend(dfiindex, dfjindex);
    }
    if (dfi->HasDependency()) {
        int dim1, dim2;
        dim1 = Reference()->SideDimension(side);
        dim2 = neighbour.Element()->Reference()->SideDimension(neighbour.Side());
        if (dim1 || dim2) {
            stringstream sout;
            sout << "TPZInterpolatedElement::RemoveSideRestraintWithRespectTo fishy " << dfiindex;
            LOGPZ_ERROR(logger, sout.str());
        }
    }
}

void TPZInterpolatedElement::RemoveSideRestraintsII(MInsertMode mode) {
    if (mode == EInsert) {//modo insercao
        LOGPZ_WARN(logger, "Exiting RemoveSideRestraintsII with mode insert should not be called");
        return;
    }
    TPZCompElSide large; //elemento grande
    TPZStack<TPZCompElSide> elemset; //elementos pequenos
    int numsides = Reference()->NSides();
    int side, nelem, iel;
    if (mode == EDelete) {
        for (side = 0; side < numsides; side++) {
            if (NSideConnects(side) == 0) {
                continue;
            }
            TPZCompElSide thisside(this, side);
            elemset.Resize(0);
            thisside.EqualLevelElementList(elemset, 1, 0); //iguais
            if (side < NCornerConnects()) thisside.HigherLevelElementList(elemset, 1, 0);
            nelem = elemset.NElements();
            TPZStack<TPZCompElSide> smallset;
            thisside.HigherLevelElementList(smallset, 1, 1); //menores
            int nsmall = smallset.NElements();
            large = thisside.LowerLevelElementList(1);
            if (!nelem) {//nao existem iguais
                if (large.Exists()) {//existe grande : a23 , ab3
                    RemoveSideRestraintWithRespectTo(side, large);
                }//fim large
                if (nsmall) {//existem pequenos : a23
                    //	    thisside.ExpandConnected(smallset,1);
                    nsmall = smallset.NElements();
                    for (iel = 0; iel < nsmall; iel++) {
                        TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (smallset[iel].Element());
                        int smallside = smallset[iel].Side();
                        if (cel && cel->NSideConnects(smallside) > 0) cel->RemoveSideRestraintWithRespectTo(smallset[iel].Side(), thisside);
                    }
                }//ab3
            }
        }
    }
    TPZCompEl *refloaded = Reference()->Reference();
    Reference()->ResetReference();
    for (side = 0; side < numsides; side++) {
        if (NSideConnects(side) == 0) {
            continue;
        }
        if (mode == EDelete) {//modo remocao
            TPZCompElSide thisside(this, side);
            elemset.Resize(0);
            thisside.EqualLevelElementList(elemset, 1, 0); //iguais
            if (side < NCornerConnects()) thisside.HigherLevelElementList(elemset, 1, 0);
            nelem = elemset.NElements();
            TPZStack<TPZCompElSide> smallset;
            thisside.HigherLevelElementList(smallset, 1, 1); //menores
            int nsmall = smallset.NElements();
            large = thisside.LowerLevelElementList(1);
            if (nelem && !large.Exists()) {//iguais e nao grande
                //se existem iguais s�deve-se recalcular a ordem de todos os iguais que restaram
                /* 	int oldorder = SideOrder(side); */
                /* 	int order = ComputeSideOrder(elemset); */
                /* 	if(order != oldorder) { */
                //        Reference()->ResetReference();
                TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (elemset[0].Element());
                if (cel) {
                    cel->IdentifySideOrder(elemset[0].Side());
                }
            } else if (!nelem) {//nao existem iguais
                if (large.Exists() && nsmall) {
                    //        Reference()->ResetReference();
                    int dim;
                    for (dim = 0; dim < 4; dim++) {
                        for (iel = 0; iel < nsmall; iel++) {
                            if (smallset[iel].Reference().Dimension() == dim) {
                                TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (smallset[iel].Element());
                                TPZInterpolatedElement *largel = dynamic_cast<TPZInterpolatedElement *> (large.Element());
                                cel->RestrainSide(smallset[iel].Side(), largel, large.Side());
                            }
                        }
                    }
                }
            }//fim nelem e large
        }//fim mode EDelete
    }//fim for
    Reference()->SetReference(refloaded);
}//fim todos

int TPZInterpolatedElement::ComputeSideOrder(TPZVec<TPZCompElSide> &smallset) {
    int nelem = smallset.NElements();
    if (!nelem) {
        LOGPZ_ERROR(logger, "Exiting ComputeSideOrder called for empty list, -1 returned");
        return -1;
    }
    TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (smallset[0].Element());
    int minorder = cel->PreferredSideOrder(smallset[0].Side());
#ifdef PZ_LOG2
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Order of first side " << minorder;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int iel;
    for (iel = 1; iel < nelem; iel++) {
        cel = dynamic_cast<TPZInterpolatedElement *> (smallset[iel].Element());
        int celorder = cel->PreferredSideOrder(smallset[iel].Side());
        minorder = minorder < celorder ? minorder : celorder;
    }
    return minorder;
}

/**Implement the refinement of an interpolated element*/
void TPZInterpolatedElement::Divide(int64_t index,TPZVec<int64_t> &sub,int interpolatesolution) {
	
	Mesh()->SetAllCreateFunctions(*this);
	
	//necessary to allow continuous and discontinuous elements in same simulation
	this->RemoveInterfaces();
	
	if (fMesh->ElementVec()[index] != this) {
		LOGPZ_ERROR(logger,"Exiting Divide: index error");
		sub.Resize(0);
		return;
	}
	TPZGeoEl *ref = Reference();
	
	TPZManVector<TPZGeoEl *,10> pv;
	ref->Divide(pv);//o elemento geometrico correspondente ao atual elemento computacional �dividido
	if(!pv.NElements()) {
		sub.Resize(0);
		if(Type() != EPoint)
			LOGPZ_ERROR(logger,"Exiting Divide: no subelements acquired");
		return;
	}
	// The refinement pattern may be adjusted during the division process
	int nsubelements = ref->NSubElements();
	sub.Resize(nsubelements);
	
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << (void*) Mesh() << " Divide " << Index() << " " << Reference()->Index();
        LOGPZ_DEBUG(loggerdiv, sout.str())
    }
#endif

    int i;

#ifdef HUGE_DEBUG
    {
        stringstream sout;
        TPZCheckMesh chk(Mesh(), &sout);
        if (chk.CheckConnectOrderConsistency() != -1) {
            LOGPZ_WARN(logger, "CheckConnectOrderConsistency failed");
            LOGPZ_WARN(logger, sout.str());
        }
    }
#endif

    RemoveSideRestraintsII(EDelete); //Cedric 25/03/99

    fMesh->ElementVec()[index] = NULL; //Cesar 2003-11-19
    Reference()->ResetReference();

    TPZGeoEl *cref;
    TPZInterpolatedElement *cel;

#ifdef HUGE_DEBUG
    if (NConnects() == 7) {
        stringstream sout;
        TPZCheckMesh chk(Mesh(), &sout);
        if (chk.CheckConnectOrderConsistency() != -1) {
            LOGPZ_WARN(logger, "CheckConnectOrderConsistency failed after RemoveSideRestraintsII");
            LOGPZ_WARN(logger, sout.str());
        }
    }
#endif

    int ncon = ref->NSides();
    fMesh->SetDefaultOrder(PreferredSideOrder(ncon - 1));
    for (i = 0; i < nsubelements; i++) {
        cref = pv[i]; //ponteiro para subelemento i
        TPZCompEl* cel = fMesh->CreateCompEl(cref);
        sub[i] = cel->Index();
        // e' assumido que CreateCompEl inseri o elemento comp no vetor de elementos da malha
    }
    if (interpolatesolution) {
        Mesh()->ExpandSolution();
        for (i = 0; i < nsubelements; i++) {
            cel = dynamic_cast<TPZInterpolatedElement *> (fMesh->ElementVec()[sub[i]]);
            if (!cel) {
                LOGPZ_WARN(logger, "Divide interpolate cast error");
                continue;
            }
            cel->CheckConstraintConsistency();
            cel->InterpolateSolution(*this);
            // e' assumido que CreateCompEl inseri o elemento comp no vetor de elementos da malha
        }
    }
    // we assume that the calling program will delete the element from
    // the data structure
    delete this; // deve ser relegado para o Refine
}

REAL TPZInterpolatedElement::CompareElement(int var, char *matname) {
    TPZMaterial * material = Material();
    if (!material) {
        PZError << "\nError at " << __PRETTY_FUNCTION__ << " - no material " << std::endl;
        LOGPZ_ERROR(logger, "InterpolateSolution no material");
        return 0.;
    }
    if (matname != material->Name()) return (0.);
    STATE error = 0.;
    int dim = Dimension();
    int numdof = material->NSolutionVariables(var);

    TPZFNMatrix<9> axes(3, 3, 0.);
    TPZFNMatrix<9> jacobian(dim, dim), jacinv(dim, dim);

    TPZManVector<STATE> sol(numdof, 0.), othersol(numdof, 0.);
    TPZVec<REAL> intpoint(dim, 0.);
    REAL weight = 0.;
    REAL detjac;

    // At this point we assume grids are identical
    TPZCompEl *otherelement = Reference()->Reference();

    const TPZIntPoints &intrule = GetIntegrationRule();
    TPZGeoEl *ref = Reference();

    for (int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind) {
        intrule.Point(int_ind, intpoint, weight);
        this->Solution(intpoint, var, sol);
        otherelement->Solution(intpoint, var, othersol);

        ref->Jacobian(intpoint, jacobian, axes, detjac, jacinv);
        weight *= fabs(detjac);

        for (int i = 0; i < sol.NElements(); i++) {
            error += (sol[i] - othersol[i])*(sol[i] - othersol[i])*(STATE) weight;
        }
    }
    return fabs(error);
}

void TPZInterpolatedElement::Print(std::ostream &out) const {
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZInterpolationSpace::Print(out);
    TPZGeoEl *gel = Reference();
    out << "Index = " << fIndex;
    out << " - Center coordinate: ";
    for (int i = 0; i < gel->NCornerNodes(); i++) {
        TPZVec< REAL > center(3, 0.);
        for (int j = 0; j < 3; j++) center[j] = gel->NodePtr(i)->Coord(j);
        out << "[" << i << "]" << center << " ";
    }
    out << std::endl;
    int nconects = this->NConnects();
    out << "Number of connects = " << this->NConnects() << " Connect indexes/NShape : ";
    int nod;
    for (nod = 0; nod < nconects/*NConnects()*/; nod++) {
        TPZConnect &c = Connect(nod);
#ifdef PZDEBUG
        int ncon_compute = NConnectShapeF(nod, c.Order());
        if (c.NShape() != ncon_compute) {
            DebugStop();
        }
#endif
        out << ConnectIndex(nod) << '/' << c.NShape() << ' ';
    }
    out << endl;
    out << "Side orders = ";
    int NSides = Reference()->NSides();
    for (int is = 0; is < NSides; is++) {
        int nsconnect = NSideConnects(is);
        out << " side " << is << " orders ";
        for (int ic = 0; ic < nsconnect; ic++) {
            int sloc = SideConnectLocId(ic, is);
            int order = Connect(sloc).Order();
            out << order << " ";
        }
        out << std::endl;
    }
    //	for (nod=0; nod< NConnects(); nod++) out << SideOrder(nod) << ' ';
    //	out << endl;
    TPZMaterial * material = Material();
    if (!material) {
        out << " no material " << std::endl;
    }

    if (material) {
        out << "material id " << material->Id() << endl;
    }
    int id;
    const TPZIntPoints &intrule = GetIntegrationRule();
    int dim = Dimension();

    TPZManVector<int> prevorder(dim);
    intrule.GetOrder(prevorder);

    int norders = prevorder.NElements();
    out << "Integration orders : \t";
    for (id = 0; id < norders; id++) {
        out << prevorder[id] << "\t";
    }

    intrule.Print(out);

    out << endl;
}

void TPZInterpolatedElement::ShortPrint(std::ostream &out) const {
    out << __PRETTY_FUNCTION__ << std::endl;
    out << std::endl;
    out << "Index = " << fIndex;
    out << " - Center coordinate: ";
    for (int i = 0; i < Reference()->NCornerNodes(); i++) {
        TPZVec< REAL > center(3, 0.);
        for (int j = 0; j < 3; j++) center[j] = Reference()->NodePtr(i)->Coord(j);
        out << "[" << i << "]" << center << " ";
    }
    int nconects = this->NConnects();
    out << "Number of connects = " << this->NConnects() << " Connect indexes/NShape : ";
    int nod;
    for (nod = 0; nod < nconects/*NConnects()*/; nod++) {
        TPZConnect &c = Connect(nod);
#ifdef PZDEBUG
        if (c.NShape() != NConnectShapeF(nod, c.Order())) {
            DebugStop();
        }
#endif
        out << ConnectIndex(nod) << '/' << c.NShape() << ' ';
    }
    out << endl;
}

void TPZInterpolatedElement::PRefine(int order) {
#ifdef PZDEBUG
    auto elementGMesh = this->fMesh->Reference();
    auto referenceGMesh = this->Reference()->Mesh();
    if(elementGMesh != referenceGMesh){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Computational mesh of this element seems not to be associated\n";
        PZError<<"with the same mesh as the one this element's reference belongs to.\n";
        PZError<<"You may want to call TPZGeoMesh::ResetReference()\n";
        PZError<<" and TPZCompMesh::LoadReferences() before PRefine\n";
    }

#endif
    SetPreferredOrder(order);

#ifdef PZ_LOG
    if (loggerdiv.isDebugEnabled()) {
        std::stringstream sout;
        sout << (void*) Mesh() << " PRefine elindex " << Index() << " gel index " << Reference()->Index() << " " << order;
        LOGPZ_DEBUG(loggerdiv, sout.str())
    }
#endif
    TPZGeoEl *gel = Reference();
    int ns = gel->NSides();
    for (int is = 0; is < ns; is++) {
        IdentifySideOrder(is);
    }
#ifdef PZ_LOG
    if (loggerdiv.isDebugEnabled()) {
        std::stringstream sout;
        sout << " PRefine connect orders ";
        int nc = NConnects();
        for(int ic=0; ic<nc; ic++) sout << (int)Connect(ic).Order() << " ";
        LOGPZ_DEBUG(loggerdiv, sout.str())
    }
#endif
}

REAL TPZInterpolatedElement::MeanSolution(int var) {
    int dim = Dimension(), nvars;
    TPZMaterial * material = Material();
    if (!material) {
        cout << __PRETTY_FUNCTION__ << " no material " << std::endl;
        LOGPZ_ERROR(logger, "Meansolution no material");
        return 0.;
    }

    nvars = material->NSolutionVariables(var);
    if (nvars != 1) {
        LOGPZ_ERROR(logger, "Exiting MeanSolution: is not implemented to nvars != 1.");
        return 0.;
    }
    TPZManVector<STATE> sol(nvars, 0.);

    int i;
    TPZFMatrix<REAL> axes(3, 3, 0.);
    TPZFMatrix<REAL> jacobian(dim, dim);
    TPZFMatrix<REAL> jacinv(dim, dim);
    REAL detjac;
    TPZVec<REAL> intpoint(dim, 0.);
    REAL weight = 0.;
    REAL meanvalue = 0.;
    REAL area = 0.;
    TPZGeoEl *ref = Reference();

    for (i = 0; i < GetIntegrationRule().NPoints(); i++) {
        GetIntegrationRule().Point(i, intpoint, weight);
        ref->Jacobian(intpoint, jacobian, axes, detjac, jacinv);

        /** Compute the solution value at point integration*/
        Solution(intpoint, var, sol);
        area += weight * fabs(detjac);
#ifdef STATE_COMPLEX
        meanvalue += (weight * fabs(detjac)) * sol[0].real(); //  meanvalue += (weight*fabs(detjac)*sol[j]);
#else
        meanvalue += (weight * fabs(detjac)) * sol[0];
#endif
    }
    return (meanvalue / area);
}

/**Compute the contribution to stiffness matrix and load vector on the element*/
template<class TVar>
void TPZInterpolatedElement::CalcIntegral(TPZElementMatrixT<TVar> &ef) {
    int i;
    TPZMaterial * material = Material();
    if (!material) {
        cout << __PRETTY_FUNCTION__ << " no material " << std::endl;
        LOGPZ_ERROR(logger, "CalcIntegral no material");
        ef.Reset();
        return;
    }

    int numdof = material->NStateVariables();
    int ncon = NConnects();
    int dim = Dimension();
    int nshape = NShapeF();
    TPZBlock &block = Mesh()->Block();
    TPZFMatrix<TVar> &solution = Mesh()->Solution();
    int numloadcases = solution.Cols();

    int numeq = nshape*numdof;
    ef.fMat.Redim(numeq, numloadcases);
    ef.fBlock.SetNBlocks(ncon);
    TPZVec<TVar> sol(numdof, 0.);
    for (i = 0; i < ncon; i++) {
        TPZConnect &c = Connect(i);
        unsigned int nshape = c.NShape();
#ifdef PZDEBUG
        if (NConnectShapeF(i, c.Order()) != nshape || c.NState() != numdof) {
            DebugStop();
        }
#endif
        ef.fBlock.Set(i, nshape * numdof);
    }
    int in, jn;
    TPZConnect *df;

    ef.fConnect.Resize(ncon);

    for (i = 0; i < ncon; ++i)
        (ef.fConnect)[i] = ConnectIndex(i);

    TPZFMatrix<REAL> phi(nshape, 1);
    TPZFMatrix<REAL> dphi(dim, nshape);
    TPZFMatrix<REAL> dphix(dim, nshape);
    TPZFMatrix<REAL> axes(3, 3, 0.);
    TPZFMatrix<REAL> jacobian(dim, dim);
    TPZFMatrix<REAL> jacinv(dim, dim);
    REAL detjac;
    TPZVec<REAL> x(3, 0.);
    TPZVec<REAL> intpoint(dim, 0.);
    REAL weight = 0.;
    TPZGeoEl *ref = Reference();

    for (int int_ind = 0; int_ind < GetIntegrationRule().NPoints(); ++int_ind) {
        GetIntegrationRule().Point(int_ind, intpoint, weight);
        ref->Jacobian(intpoint, jacobian, axes, detjac, jacinv);
        ref->X(intpoint, x);
        weight *= fabs(detjac);
        Shape(intpoint, phi, dphi);

        int64_t l, iv = 0;
        TVar coef;
        for (int ist = 0; ist < numloadcases; ist++) {
            for (in = 0; in < numdof; in++) sol[in] = 0.;
            for (in = 0; in < ncon; in++) {
                df = &Connect(in);
                int64_t dfseq = df->SequenceNumber();
                int dfvar = block.Size(dfseq);
                for (jn = 0; jn < dfvar; jn++) {
                    int64_t pos = block.Position(dfseq);
                    coef = solution(pos + jn, ist);
                    sol[iv % numdof] += (TVar) phi(iv / numdof, 0) * coef;
                    iv++;
                }
            }
            for (in = 0; in < nshape; in++)
                for (l = 0; l < numdof; l++)
                    (ef.fMat)(in * numdof + l, 0) += (TVar) weight * (TVar) phi(in, 0) * sol[l];
        }
    }
}

int TPZInterpolatedElement::AdjustPreferredSideOrder(int side, int order) {
    TPZGeoEl *gel = Reference();
    if (!gel) {
        LOGPZ_ERROR(logger, "Exiting AdjustPreferredSideOrder: null reference element.");
        return order;
    }
    int dim = gel->SideDimension(side);
    if (dim != 2 || ConnectIndex(MidSideConnectLocId(side)) == -1) {
        //    std::stringstream sout;
        //    sout << "Exiting AdjustPreferredSideOrder: dimension != 2 " << dim;
        //    LOGPZ_ERROR(logger,sout.str().c_str());
        return order;
    }
    TPZGeoElSide gelside(gel, side);
    TPZStack<TPZCompElSide> celsidestack;
    gelside.HigherLevelCompElementList2(celsidestack, 1, 1);
    // if there are no smaller elements do nothing
    if (celsidestack.size() == 0) {
        return order;
    }
    TPZStack<int> elsides;
    gel->LowerDimensionSides(side, elsides);
    int maxorder = order;
    int nsides = elsides.NElements();
    int is;
    for (is = 0; is < nsides; is++) {
        TPZGeoElSide gelside(gel, elsides[is]);
        if (gelside.Dimension() != 1) continue;
        int s = elsides[is];
        int sorder = MidSideConnect(s).Order();
        maxorder = maxorder < sorder ? sorder : maxorder;
    }
    return maxorder;
}

int TPZInterpolatedElement::ClassId() const {
    return Hash("TPZInterpolatedElement") ^ TPZInterpolationSpace::ClassId() << 1;
}

/** Save the element data to a stream */
void TPZInterpolatedElement::Write(TPZStream &buf, int withclassid) const {
    TPZInterpolationSpace::Write(buf, withclassid);
}

/** Read the element data from a stream */
void TPZInterpolatedElement::Read(TPZStream &buf, void *context) {
    TPZInterpolationSpace::Read(buf, context);
}


/** @brief adds the connect indexes associated with base shape functions to the set */
void TPZInterpolatedElement::BuildCornerConnectList(std::set<int64_t> &connectindexes) const {
    int ncorner = NCornerConnects();
    for (int ic = 0; ic < ncorner; ic++) {
        connectindexes.insert(ConnectIndex(ic));
    }
}

/// return true if the connects associated with the side have dependency with large and if the dependency dimensions match

bool TPZInterpolatedElement::VerifyConstraintConsistency(int side, TPZCompElSide large) const {
    int midsideconnectlocindex = MidSideConnectLocId(side);
    // there is no connect associated with the side
    if (midsideconnectlocindex < 0) {
        return true;
    }
    TPZGeoEl *gelref = Reference();
    TPZGeoElSide gelside(gelref, side);
    TPZConnect &c = Connect(midsideconnectlocindex);
    TPZInterpolatedElement *largel = dynamic_cast<TPZInterpolatedElement *> (large.Element());
    if (!largel) {
        DebugStop();
    }
    int nsideconnectslarge = largel->NSideConnects(large.Side());
    /// the number of sides with nonzero number of shapefunctions
    int nonzerosides = 0;
    TPZConnect &largec = largel->MidSideConnect(large.Side());
    for (int ic = 0; ic < nsideconnectslarge; ic++) {
        TPZConnect &cl = largel->SideConnect(ic, large.Side());
        if (cl.NShape() != 0) {
            nonzerosides++;
        }
    }
    int largeorder = largel->EffectiveSideOrder(large.Side());
    if (gelside.Dimension() != 0 && c.Order() != largeorder) {
        return false;
    }
    // if we optimized the dependency list to exclude zero norm dependencies, then we will pay the price of recomputing every time
    int cnumdepend = c.NumDepend();
    if (cnumdepend != nonzerosides && c.NShape() != 0) {
        return false;
    }
    int nshapeconnect = c.NShape();
    if (cnumdepend == 0 && nshapeconnect != 0) {
        PZError<<"Error at: "<<__PRETTY_FUNCTION__<<std::endl;
        PZError<<"cnumdepend: "<<cnumdepend<<std::endl;
        PZError<<"nshapeconnect: "<<nshapeconnect<<std::endl;
        DebugStop();
    }
    std::map<int64_t, int> largenshapeconnect;
    for (int ic = 0; ic < nsideconnectslarge; ic++) {
        TPZConnect &clarge = largel->SideConnect(ic, large.Side());
        int64_t clargeindex = largel->SideConnectIndex(ic, large.Side());
        largenshapeconnect[clargeindex] = clarge.NShape();
    }
    TPZConnect::TPZDepend *depend = c.FirstDepend();
    while (depend) {
        int nrow = depend->fDepMatrix.Rows();
        int ncol = depend->fDepMatrix.Cols();
        /// the number of rows of the dependency matrix is equal to the number of shape functions of the connect associated with the small side
        if (nrow != c.NShape()) {
            return false;
        }
#ifdef PZDEBUG
        if (largenshapeconnect.find(depend->fDepConnectIndex) == largenshapeconnect.end()) {
            DebugStop();
        }
#endif
        /// the number of columns of the dependency matrix is equal to the number of shape functions of the connect of the large element
        if (ncol != largenshapeconnect[depend->fDepConnectIndex]) {
            return false;
        }
        depend = depend->fNext;
    }
    return true;
}

void TPZInterpolatedElement::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsContinuous();
}


template
void TPZInterpolatedElement::
CalcIntegral<STATE>(TPZElementMatrixT<STATE>&ef);
