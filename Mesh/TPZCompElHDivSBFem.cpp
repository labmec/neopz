/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivSBFem methods.
 */


#include "TPZCompElHDivSBFem.h"
#include "TPZShapeHDivCollapsed.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZMaterialDataT.h"
#include "pzelchdiv.h"
#include "pzvec_extras.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivSBFem"));
#endif

// Initialize with the geometry of the SBFemVolume
template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel, TPZGeoElSide &gelside) :
TPZCompElHDivCollapsed<TSHAPE>(mesh, gel), fGeoElVolSide(gelside)
{
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel) :
TPZCompElHDivCollapsed<TSHAPE>(mesh, gel), fGeoElVolSide(0)
{
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem(TPZCompMesh &mesh, const TPZCompElHDivSBFem<TSHAPE> &copy) :
TPZCompElHDivCollapsed<TSHAPE>(mesh,copy)
{
    fGeoElVolSide = copy.fGeoElVolSide;
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem() : TPZCompElHDivCollapsed<TSHAPE>()
{
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::~TPZCompElHDivSBFem()
{
}

// Set the GeoElSide for the volumetric SBFEM element
template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::SetGelVolumeSide(TPZGeoElSide &gelside)
{
    fGeoElVolSide = gelside;
}

// Get the GeoElSide for the volumetric SBFEM element
template<class TSHAPE>
TPZGeoElSide & TPZCompElHDivSBFem<TSHAPE>::GetGelVolumeSide()
{
    return fGeoElVolSide;
}

template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi)
{
    // Compute data for the 1d functions first
    bool needssol = data.fNeedsSol;
    data.fNeedsSol = false;
    TPZCompElHDivCollapsed<TSHAPE>::ComputeRequiredData(data, qsi);
    data.fNeedsSol = needssol;

    TPZShapeData &shapedata = data;

    int porderboundbottom = this->ConnectOrder(TSHAPE::NFacets+1); // to be implemented
    int porderboundtop = this->ConnectOrder(TSHAPE::NFacets+2);

    int nshape = TPZShapeHDivCollapsed<TSHAPE>::NShapeF(data);
    int nshapeboundleft = TPZCompElHDivCollapsed<TSHAPE>::NConnectShapeF(TSHAPE::NFacets+1, porderboundbottom);
    int nshapeboundright = TPZCompElHDivCollapsed<TSHAPE>::NConnectShapeF(TSHAPE::NFacets+2, porderboundtop);
    int nshape1d = nshape - nshapeboundleft - nshapeboundright;

    ComputeSBFemVolumeHdivData(data, nshape1d);

    // Adjusting divergence and phi values
    for (int i = 0; i < nshape1d; i++)
    {
        data.divphi(i) = data.fDPhi(0,i);
        data.phi(i) = data.fPhi(i);
    }
    for (int i = 0; i < nshapeboundleft; i++)
    {
        data.divphi(i+nshape1d) = 0;
    }
    for (int i = 0; i < nshapeboundright; i++)
    {
        data.divphi(i+nshape1d+nshapeboundleft) = data.fPhi(i+nshape1d);
        data.phi(i+nshape1d+nshapeboundleft) = 0;
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        stringstream sout;
        sout << "qsi = " << data.xParametric << endl;
        data.phi.Print("phi = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger,sout.str());
    }
#endif
}

// This function will compute the directions for the HDiv collapsed based on the information of neighbourhood
template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::ComputeSBFemVolumeHdivData(TPZMaterialDataT<STATE> &data, int64_t nshape1d)
{
    // Computing the deformed directions for the 2d functions using the information of the neighbourhood
    // Inspired in TPZSBFemVolume::ComputeKMatrices

    // Before entering in this method, the material data was computed
    // (ComputeRequiredData calls this method and before it called ComputeRequiredData from TPZCompElHDivSBFem)

    // The Reference element will be the skeleton
    TPZGeoEl *Ref1D = this->Reference();
    int dim1 = TSHAPE::Dimension;
    auto detjac1d = data.detjac;
    auto axes1d = data.axes;

    // The Volume element will be the one passed as an argument by the TPZSBFemVolumeHDiv
    // The SBFemVolumeHdiv will compute the Contribute for the SBFemHdiv element
    TPZGeoEl * gelvolume = fGeoElVolSide.Element();
    int dim2 = gelvolume->Dimension();

    // Find the higher side of the skeleton element
    TPZGeoElSide SkeletonSide(Ref1D, Ref1D->NSides() - 1);

    // Create a transformation between the sides: Skeleton and Volume
    TPZTransform<REAL> tr(dim2, dim1);
    tr = SkeletonSide.NeighbourSideTransform(fGeoElVolSide);
    TPZTransform<REAL> t2 = gelvolume->SideToSideTransform(fGeoElVolSide.Side(), gelvolume->NSides() - 1);
    tr = t2.Multiply(tr);

    // Applying the transformation to obtain the 2d parametric coordinate
    TPZManVector<REAL, 3> xi(dim1), xiquad(dim2), xivol(dim2);
    xi = data.xParametric;
    tr.Apply(xi, xiquad);
    xivol = xiquad;
    xivol[dim2 - 1] = -0.5;

    // Computing data for the volumetric element
    gelvolume->X(xivol, data.x);
    gelvolume->Jacobian(xiquad, data.jacobian, data.axes, data.detjac, data.jacinv);

    auto axes = data.axes;
    if (dim2 == 3)
    {
        PZError << __PRETTY_FUNCTION__ << "is not ready for 3D examples yet. \n";
        DebugStop();
        AdjustAxes3D(axes, data.axes, data.jacobian, data.jacinv, data.detjac);
    }

    // Adjusting derivatives
    ExtendShapeFunctions(data, detjac1d);

    auto ndir = data.fDeformedDirections.Cols()-1;
    
    TPZFNMatrix<9,REAL> grad(dim2,dim2,0);
    gelvolume->GradX(xivol,grad);
    REAL detjac = sqrt(2*data.detjac);
    
    TPZFNMatrix<9,REAL> jaccollapsed;
    data.axes.Resize(dim2,dim2);
    // data.fDeformedDirections.Resize(3,data.fMasterDirections.Cols());
    data.fDeformedDirections.Zero();
    TPZFMatrix<REAL> internaldir(data.fDeformedDirections);
    TPZIntelGen<TSHAPE>::Reference()->HDivDirections(data.xParametric, internaldir);

    auto signal = 0, signal0 = 0;
    TPZManVector<REAL,2> signalvec(2,0);
    // for (int i = 0; i < dim2; i++)
    // {
    //     signal += -data.axes(1,i);
    // }
    TPZFMatrix<REAL> directions(2,2,0.);
    data.axes.Multiply(data.jacobian,directions);
    
    for (auto j = 0; j < nshape1d; j++)
    {
        for (auto i = 0; i < dim2; i++)
        {
            data.fDeformedDirections(i,j) = -internaldir(i,0);
        }
    }
    for (int j = 0; j < data.fHDivNumConnectShape[3]; j++)
    {
        for (auto i = 0; i < dim2; i++)
        {
            data.fDeformedDirections(i,nshape1d+j) = this->fbottom_side_orient*grad(i,1)*2./fabs(detjac1d);
        }
    }
    for (int j = 0; j < data.fHDivNumConnectShape[4]; j++)
    {
        for (auto i = 0; i < dim2; i++)
        {
            data.fDeformedDirections(i,nshape1d+data.fHDivNumConnectShape[3]+j) = this->ftop_side_orient*grad(i,1)*2./fabs(detjac1d);
        }
    }
}

template <class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::AdjustAxes3D(const TPZFMatrix<REAL> &axes2D, TPZFMatrix<REAL> &axes3D, TPZFMatrix<REAL> &jac3D, TPZFMatrix<REAL> &jacinv3D, REAL detjac)
{
    TPZManVector<REAL, 3> ax1(3), ax2(3), ax3(3);
    for (int i = 0; i < 3; i++) {
        ax1[i] = axes2D.g(0, i);
        ax2[i] = axes2D.g(1, i);
        Cross(ax1, ax2, ax3);
    }
    for (int i = 0; i < 3; i++) {
        axes3D(0, i) = ax1[i];
        axes3D(1, i) = ax2[i];
        axes3D(2, i) = ax3[i];
        if (detjac < 0.) {
            axes3D(2, i) *= -1.;
        }
    }
    TPZFNMatrix<9, REAL> jacnew(3, 3), axest(3, 3), jacinv(3, 3);
    axes3D.Transpose(&axest);
    axes3D.Multiply(jac3D, jacnew);
    jacinv3D.Multiply(axest, jacinv);
    jac3D = jacnew;
    jacinv3D = jacinv;
#ifdef PZDEBUG
    // check whether the axes are orthogonal and whether the jacobian is still the inverse of jacinv
    {
        TPZFNMatrix<9, REAL> ident1(3, 3, 0.), ident2(3, 3, 0.), identity(3, 3);
        identity.Identity();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    ident1(i, j) += axes3D(i, k) * axes3D(j, k);
                    ident2(i, j) += jac3D(i, k) * jacinv3D(k, j);
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (fabs(ident1(i, j) - identity(i, j)) > 1.e-6) {
                    DebugStop();
                }
                if (fabs(ident2(i, j) - identity(i, j)) > 1.e-6) {
                    DebugStop();
                }
            }
        }
    }
#endif
}

template <class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::ExtendShapeFunctions(TPZMaterialDataT<STATE> &data, REAL detjac1d)
{
    int dim = fGeoElVolSide.Element()->Dimension();

    TPZShapeData& shapedata = data;
    int porderboundbottom = this->ConnectOrder(TSHAPE::NFacets+1);
    int porderboundtop = this->ConnectOrder(TSHAPE::NFacets+2);

    auto nshape = TPZShapeHDivCollapsed<TSHAPE>::NShapeF(shapedata);
    int nshapeboundleft = TPZCompElHDivCollapsed<TSHAPE>::NConnectShapeF(TSHAPE::NFacets+1, porderboundbottom);
    int nshapeboundright = TPZCompElHDivCollapsed<TSHAPE>::NConnectShapeF(TSHAPE::NFacets+2, porderboundtop);
    int nshape1d = nshape - nshapeboundleft - nshapeboundright;

    shapedata.fDPhi.Resize(dim, nshape);

    for (int ish = 0; ish < nshapeboundleft; ish++)
    {
        for (int d = 0; d < dim - 1; d++)
        {
            shapedata.fDPhi(d, ish + nshape1d) = 0.;
            shapedata.fDPhi(d, ish + nshape1d+nshapeboundleft) = shapedata.fPhi(ish);
        }
        shapedata.fDPhi(dim-1, ish+nshape1d) = -shapedata.fPhi(ish) / 2.;
        shapedata.fDPhi(dim-1, ish+nshape1d+nshapeboundleft) = 0.;
    }
   
    TPZFNMatrix<50,REAL> philoc(data.fPhi.Rows(),data.fPhi.Cols(),0.),
        dphiloc(data.fDPhi.Rows(),data.fDPhi.Cols(),0.);
    TPZManVector<int,TSHAPE::NSides> ord(TSHAPE::NFacets+3, 0);

    const int nconnects = shapedata.fHDivConnectOrders.size();
    ord = shapedata.fHDivConnectOrders;

    int nc = this->Reference()->NCornerNodes();
    TPZManVector<int64_t,8> id(nc);
    for (int ic=0; ic<nc; ic++)
    {
        id[ic] = this->Reference()->Node(ic).Id();
    }

    TSHAPE::Shape(data.xParametric, id, ord, philoc, dphiloc);

    TPZManVector<int,9> permutegather(TSHAPE::NSides);
    int transformid = TSHAPE::GetTransformId(id);
    TSHAPE::GetSideHDivPermutation(transformid, permutegather);
    
    TPZManVector<int64_t,27> FirstIndex(TSHAPE::NSides+1);
    fCelFlux->FirstShapeIndex(FirstIndex);
    int order = shapedata.fHDivConnectOrders[nconnects-2];
    
    for (int side=0; side < TSHAPE::NSides; side++)
    {
        int ifirst = FirstIndex[side];
        int kfirst = FirstIndex[permutegather[side]];
        int nshapeloc = TSHAPE::NConnectShapeF(side,order);
        for (int i=0; i<nshapeloc; i++)
        {
            data.fPhi(nshape1d + ifirst+i,0) = philoc(kfirst+i,0);
            for (int d=0; d< TSHAPE::Dimension; d++)
            {
                data.fDPhi(d,nshape1d + ifirst+i) = dphiloc(d,kfirst+i);
            }
        }
    }

    data.divsol.Resize(1);
    data.divsol[0].Resize(1,0.);

    TPZInterpolationSpace::Convert2Axes(data.fDPhi, data.jacinv, data.dphix);
}

template <class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialDataT<STATE> &data)
{
    DebugStop();
}

#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;

template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeQuad>>;

template class TPZCompElHDivSBFem<TPZShapeTriang>;
template class TPZCompElHDivSBFem<TPZShapeLinear>;
template class TPZCompElHDivSBFem<TPZShapeQuad>;
