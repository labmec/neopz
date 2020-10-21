/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivCollapsed methods.
 */


#include "TPZCompElHDivCollapsed.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzelchdiv.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivCollapsed"));
#endif

template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,gel,index), fBottom(mesh,gel,index), fTop(mesh,gel,index)
{
    index = this->fIndex;
    mesh.ElementVec().SetFree(fTop.Index());
    mesh.ElementVec().SetFree(fBottom.Index());
    int64_t bottom_c_index = fBottom.ConnectIndex(0);
    int64_t top_c_index = fTop.ConnectIndex(0);
    fBottom.SetIndex(-1);
    fTop.SetIndex(-1);
    

#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	 {
         std::stringstream sout;
         sout << "Finalizando criacao do elemento ";
         this->Print(sout);
         LOGPZ_DEBUG(logger,sout.str())
	 }
#endif
	 
}

template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh, const TPZCompElHDivCollapsed<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy),fBottom(copy.fBottom), fTop(copy.fTop)
{
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh,
												 const TPZCompElHDivCollapsed<TSHAPE> &copy,
												 std::map<int64_t,int64_t> & gl2lcConMap,
												 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap),fBottom(mesh,copy.fBottom,gl2lcConMap,gl2lcElMap),
fTop(mesh,copy.fBottom,gl2lcConMap,gl2lcElMap)
{
	
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<=TSHAPE::NFacets;i++)
	{
		int64_t lcIdx = -1;
		int64_t glIdx = copy.fConnectIndexes[i];
		if(glIdx == -1)
		{
			// nothing to clone
			this->fConnectIndexes[i] = -1;
			continue;
		}
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end())
		{
			lcIdx = gl2lcConMap[glIdx];
		}
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			this-> fConnectIndexes[i] = -1;
			return;
		}
		this-> fConnectIndexes[i] = lcIdx;
	}
    // write the code when this constructor is called
    DebugStop();
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed() :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(),fBottom(),fTop()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NFacets+1;i++) {
		this-> fConnectIndexes[i] = -1;
	}
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::~TPZCompElHDivCollapsed(){
    TPZGeoEl *gel = this->Reference();
    if(!gel) return;
    if (gel && gel->Reference() != this) {
        return;
    }
    int side = TSHAPE::NSides-1;
    TPZGeoElSide gelside(this->Reference(),side);
    TPZStack<TPZCompElSide> celstack;
    TPZCompElSide largecel = gelside.LowerLevelCompElementList2(0);
    if (largecel) {
        int cindex = this->SideConnectLocId(0, side);
        TPZConnect &c = this->Connect(cindex);
        c.RemoveDepend();
    }
    gelside.HigherLevelCompElementList3(celstack, 0, 1);
    int64_t ncel = celstack.size();
    for (int64_t el=0; el<ncel; el++) {
        TPZCompElSide celsidesmall = celstack[el];
        TPZGeoElSide gelsidesmall = celsidesmall.Reference();
        if (gelsidesmall.Dimension() != gel->Dimension()) {
            continue;
        }
        TPZCompEl *cel = celsidesmall.Element();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            DebugStop();
        }
        int cindex = intel->SideConnectLocId(0, celsidesmall.Side());
        TPZConnect &c = intel->Connect(cindex);
        c.RemoveDepend();
    }
    if (gel){
        gel->ResetReference();
    }

}

// NAO TESTADO
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetSideOrient(int side, int sideorient)
{
    if (side < TSHAPE::NFacets || side >= TSHAPE::NSides) {
        DebugStop();
    }
    if(side < TSHAPE::NSides-1) TPZCompElHDiv<TSHAPE>::SetSideOrient(side, sideorient);
    else fBottom.SetSideOrient(side, sideorient);
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::GetSideOrient(int side)
{
    if (side >= TSHAPE::NFacets && side < TSHAPE::NSides - 1) {
        return TPZCompElHDiv<TSHAPE>::GetSideOrient(side);
    }
    else if(side == TSHAPE::NSides-1) return fBottom.GetSideOrient(side);
    DebugStop();
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::NConnects() const {
	
	return TPZCompElHDiv<TSHAPE>::NConnects()+2;
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
	if(i < TSHAPE::NFacets)
	{
        SetConnectIndex(i, connectindex);
	}
    else if(i == TSHAPE::NFacets)
    {
        fBottom.SetConnectIndex(0, connectindex);
    }
    else if(i == TSHAPE::NFacets+1)
    {
        fTop.SetConnectIndex(0, connectindex);
    }
    else
    {
        DebugStop();
    }
	
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::NConnectShapeF(int connect, int connectorder) const
{
    if(connect <= TSHAPE::NFacets)
    {
        return TPZCompElHDiv<TSHAPE>::NConnectShapeF(connect, connectorder);
    }
    else if(connect <= TSHAPE::NFacets+2)
    {
        return fBottom.NConnectShapeF(0, connectorder);
    }
    else DebugStop();
    return -1;
}


/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZCompElHDiv<TSHAPE>::InitMaterialData(data);
    if(data.fUserData) DebugStop();
    auto datapair = new std::pair<TPZMaterialData,TPZMaterialData>;
    data.fUserData = datapair;
    TPZMaterialData &datatop = datapair->second, &databottom = datapair->first;
    fTop.InitMaterialData(datatop);
    fBottom.InitMaterialData(databottom);
    // expand the shape vector and normal vector
    int nvecshape = data.fVecShapeIndex.size();
    int nscalar = data.phi.Rows();
    int nscalartop = datatop.phi.Rows();
    int nscalarbottom = databottom.phi.Rows();
    int nvec = data.fDeformedDirections.Cols();
    const int dim = TSHAPE::Dimension;
    data.fMasterDirections.Resize(dim+1, nvec+2);
    data.fMasterDirections(dim,nvec) = 1.;
    data.fMasterDirections(dim,nvec+1) = -1.;
    data.fDeformedDirections.Resize(3, nvec+2);
    data.fVecShapeIndex.Resize(nvecshape+nscalartop+nscalarbottom, {0,0});
    for(int i=0; i<nscalartop; i++) data.fVecShapeIndex[nvecshape+i] = std::pair<int,int64_t>(nvec,nscalar+i);
    for(int i=0; i<nscalarbottom; i++) data.fVecShapeIndex[nvecshape+nscalartop+i] = std::pair<int,int64_t>(nvec+1,nscalar+nscalartop+i);
    data.phi.Resize(nscalar+nscalartop+nscalarbottom, 1);
    data.dphi.Resize(dim+1,nscalar+nscalartop+nscalarbottom);
    data.divphi.Resize(nscalar+nscalartop+nscalarbottom,1);
#ifdef _AUTODIFF
    if(data.fNeedsDeformedDirectionsFad) DebugStop();
#endif
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
        std::stringstream sout;
        sout << "After InitMaterialData\n";
        data.Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
    
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
    if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension || point.size() != TSHAPE::Dimension ){
		DebugStop() ;
	}
    fBottom.SideShapeFunction(side, point, phi, dphi);
    
    return;
}

/** Compute the shape function at the integration point */
//template<class TSHAPE>
//void TPZCompElHDivCollapsed<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
//{
//    fBottom.Shape(pt, phi, dphi);
//}

/** Read the element data from a stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZIntelGen<TSHAPE>::Read(buf,context);
    fBottom.Read(buf,context);
    fTop.Read(buf, context);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	TPZIntelGen<TSHAPE>::Write(buf,withclassid);
    fBottom.Write(buf, false);
    fTop.Write(buf, false);
}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZCompElHDiv<TSHAPE>::Print(out);
    fBottom.Print(out);
    fTop.Print(out);
    
}

#include "pzvec_extras.h"

static void ExpandAxes(TPZFMatrix<REAL> &axinput, TPZMatrix<REAL> &axout)
{
    int dim = axinput.Rows();
    switch (dim) {
        case 1:
        {
            REAL norms[3];
            TPZManVector<REAL,3> v1(3),try1[3];
            for(int i=0; i<3; i++) v1[i] = axinput(0,i);
            for(int i=0;i<3;i++)
            {
                TPZManVector<REAL,3> uni(3,0.);
                uni[i] = 1.;
                try1[i].resize(3);
                Cross(v1,uni,try1[i]);
                norms[i] = Norm<REAL>(try1[i]);
            }
            int maxi = 0;
            if(norms[1]>norms[0]) maxi = 1;
            if(norms[2]>norms[maxi]) maxi = 2;
            for(int i=0; i<3; i++)
            {
                axout(0,i) = axinput(0,i);
                axout(1,i) = try1[maxi][i]/norms[maxi];
            }
        }
            break;
        case 2:
        {
            TPZManVector<REAL,3> v1(3),v2(3),v3(3);
            for(int i=0; i<3; i++)
            {
                v1[i] = axinput(0,i);
                v2[i] = axinput(1,i);
            }
            Cross<REAL>(v1,v2,v3);
            for(int i=0; i<3; i++)
            {
                axout(0,i) = v1[i];
                axout(1,i) = v2[i];
                axout(2,i) = v3[i];
            }

        }
            break;
        default:
            DebugStop();
            break;
    }
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
    
//    TPZManVector<int,TSHAPE::NSides*TSHAPE::Dimension> normalsidesDG(TSHAPE::Dimension*TSHAPE::NSides);

    TPZCompElHDiv<TSHAPE>::ComputeRequiredData(data, qsi);
    TPZFNMatrix<9,REAL> axeslocal(TSHAPE::Dimension+1,3);
    ExpandAxes(data.axes, axeslocal);
    data.axes = axeslocal;
    int dim = TSHAPE::Dimension+1;
    // compute the deformed directions for the two additional vectors
    {
        int nvec = TSHAPE::NSides*TSHAPE::Dimension;
        TPZFNMatrix<3,REAL> masterdir(TSHAPE::Dimension+1,2);
        for(int i=0; i<3; i++){
            for(int k=0; k<2; k++){
                data.fDeformedDirections(i,nvec+k) = 0.;
                for(int l=0; l<dim; l++){
                    data.fDeformedDirections(i,nvec+k) += data.axes(l,i)*data.fMasterDirections(l,nvec+k);
                }
            }
        }
    }
    data.ComputeFunctionDivergence();
    std::pair<TPZMaterialData,TPZMaterialData> *datapair = (std::pair<TPZMaterialData,TPZMaterialData> *) data.fUserData;
    TPZMaterialData &datatop = datapair->second, &databottom = datapair->first;

    int nsides = this->Reference()->NSides();
    // compute the divergence of the top and bottom elements
    // the value is the value of the shape function times the sign of the vector in master direction
    {
        fTop.ComputeRequiredData(datatop, qsi);
        fBottom.ComputeRequiredData(databottom, qsi);
        int64_t numvec = data.divphi.Rows();
        int64_t nvec_top = datatop.phi.Rows();
        int64_t nvec_bottom = databottom.phi.Rows();
        int64_t nvec_hdiv = numvec-nvec_top-nvec_bottom;
        //
        for (int64_t i= nvec_hdiv; i<numvec-nvec_top; i++) {
            data.divphi(i) = databottom.phi(i-nvec_hdiv);
        }
        for (int64_t i= nvec_hdiv+nvec_bottom; i<numvec; i++) {
            data.divphi(i) = databottom.phi(i-nvec_hdiv-nvec_bottom);
        }
    }
    
    if (data.fNeedsSol) {
        TPZCompElHDiv<TSHAPE>::ComputeSolution(qsi, data);
    }


#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        data.fDeformedDirections.Print("Normal Vectors " , sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    

}//void

template<class TSHAPE>
int64_t TPZCompElHDivCollapsed<TSHAPE>::ConnectIndex(int con) const
{
    if(con <= TSHAPE::NFacets) return TPZCompElHDiv<TSHAPE>::ConnectIndex(con);
    if(con > TSHAPE::NFacets + 2) DebugStop();
    if(con == TSHAPE::NFacets+1) return fBottom.ConnectIndex(0);
    if(con == TSHAPE::NFacets+2) return fTop.ConnectIndex(0);
    DebugStop();
    return -1;
}

/**
 * @brief Destroy internally allocated data structures
 */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::CleanupMaterialData(TPZMaterialData &data)
{
    std::pair<TPZMaterialData,TPZMaterialData> *userdata = (std::pair<TPZMaterialData,TPZMaterialData> *) data.fUserData;
    delete userdata;
    data.fUserData = 0;
}




template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsHDiv();
}

#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;

template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeQuad>>;


template class TPZCompElHDivCollapsed<TPZShapeTriang>;
template class TPZCompElHDivCollapsed<TPZShapeLinear>;
template class TPZCompElHDivCollapsed<TPZShapeQuad>;

