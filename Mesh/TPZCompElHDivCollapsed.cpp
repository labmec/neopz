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
#include "TPZMaterialDataT.h"
#include "pzelchdiv.h"
#include "TPZShapeHDivCollapsed.h"


#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDivCollapsed");
#else
static int logger;
#endif

template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh, TPZGeoEl *gel) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,gel)
{
//    fbottom_c_index = -1; set at constructor in .h
//    ftop_c_index = -1; set at constructor in .h
//    fbottom_side_orient = -1; in .h
//    ftop_side_orient = 1; in .h
    this->Reference()->SetReference(this);
    

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
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
TPZCompElHDiv<TSHAPE>(mesh,copy),fbottom_c_index(copy.fbottom_c_index), ftop_c_index(copy.ftop_c_index),
fbottom_side_orient(copy.fbottom_side_orient), ftop_side_orient(copy.ftop_side_orient)
{
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh,
												 const TPZCompElHDivCollapsed<TSHAPE> &copy,
												 std::map<int64_t,int64_t> & gl2lcConMap,
												 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap)
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
    // Initialize top and bottom connects properly if needed be!
    DebugStop();
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed() :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>()
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
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides + 1) {
        DebugStop();
    }

    if(side < TSHAPE::NSides-1) TPZCompElHDiv<TSHAPE>::SetSideOrient(side, sideorient);
    else if(side == TSHAPE::NSides-1 ) fbottom_side_orient = sideorient;
    else ftop_side_orient = sideorient;
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::GetSideOrient(int side)
{
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides + 1) {
        DebugStop();
    }

    if (side < TSHAPE::NSides - 1) {
        return TPZCompElHDiv<TSHAPE>::GetSideOrient(side);
    }
    else if(side == TSHAPE::NSides-1) return fbottom_side_orient;
    else return ftop_side_orient;
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::NConnects() const {
    int nconnects = TPZCompElHDiv<TSHAPE>::NConnects();
    if (fbottom_c_index != -1) {
        nconnects++;
    }
    if (ftop_c_index != -1) {
        nconnects++;
    }

	return nconnects;
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
	if(i <= TSHAPE::NFacets)
	{
        TPZCompElHDiv<TSHAPE>::SetConnectIndex(i, connectindex);
	}
    else if(i == TSHAPE::NFacets+1)
    {
        fbottom_c_index = connectindex;
    }
    else if(i == TSHAPE::NFacets+2)
    {
        ftop_c_index = connectindex;
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
    else if(connect == TSHAPE::NFacets+1)
    {
        TPZManVector<int,22> order(TSHAPE::NSides-TSHAPE::NCornerNodes,connectorder);
        return TSHAPE::NShapeF(order);
    }
    else if(connect == TSHAPE::NFacets+2)
    {
        TPZManVector<int,22> order(TSHAPE::NSides-TSHAPE::NCornerNodes,connectorder);
        return TSHAPE::NShapeF(order);
    }
    else DebugStop();
    return -1;
}


/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::InitMaterialData(TPZMaterialData &data){
    auto *tmp =
        dynamic_cast<TPZMaterialDataT<STATE>*>(&data);
    if(tmp){
        InitMaterialDataT(*tmp);
    }else{
        auto *tmp =
        dynamic_cast<TPZMaterialDataT<CSTATE>*>(&data);
        if(tmp){
            InitMaterialDataT(*tmp);
        }
    }
}
template<class TSHAPE>
template<class TVar>
void TPZCompElHDivCollapsed<TSHAPE>::InitMaterialDataT(TPZMaterialDataT<TVar> &data)
{
    TPZShapeData& shapedata = data;
    
    // Order of top and bot connect
    TPZConnect& connbot = this->Mesh()->ConnectVec()[fbottom_c_index];
    TPZConnect& conntop = this->Mesh()->ConnectVec()[ftop_c_index];
    const int toporder = conntop.Order();
    const int bottomorder = connbot.Order();
    
    // Number of shape and dim
    const int64_t nvecshapecollpased = this->NShapeF();
    const int dim = TSHAPE::Dimension+1;
    
    // Order of connects
    int ncon = NConnects();
    TPZManVector<int> connectorders(ncon,0);
    for(int i=0; i<TSHAPE::NFacets+1; i++) connectorders[i] = this->Connect(i).Order();
    connectorders[ncon-2] = bottomorder;
    connectorders[ncon-1] = toporder;
    
    // Side orient
    TPZManVector<int> sideorient(TSHAPE::NFacets+2,0);
    for(int i=0; i<TSHAPE::NFacets; i++) sideorient[i] = this->SideOrient(i);
    sideorient[TSHAPE::NFacets] = GetSideOrient(TSHAPE::NSides-1);
    sideorient[TSHAPE::NFacets+1] = GetSideOrient(TSHAPE::NSides);
    
    // Node ids of geoel
    TPZGeoEl *gel = this->Reference();
    TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes);
    for(int i=0; i<TSHAPE::NCornerNodes; i++) ids[i] = gel->NodeIndex(i);
    
    // Init shape data
    TPZShapeHDivCollapsed<TSHAPE>::Initialize(ids, connectorders, sideorient, shapedata);

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
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
	
    DebugStop(); // is this function ever used? If so, check if it is ok without fBottom
    if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension || point.size() != TSHAPE::Dimension ){
		DebugStop() ;
	}
//    fBottom.SideShapeFunction(side, point, phi, dphi);
    
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
    buf.Read(&fbottom_c_index,1);
    buf.Read(&ftop_c_index,1);
    buf.Read(&fbottom_side_orient,1);
    buf.Read(&ftop_side_orient,1);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	TPZIntelGen<TSHAPE>::Write(buf,withclassid);
    buf.Write(&fbottom_c_index,1);
    buf.Write(&ftop_c_index,1);
    buf.Write(&fbottom_side_orient,1);
    buf.Write(&ftop_side_orient,1);
}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZCompElHDiv<TSHAPE>::Print(out);
    out << "\nfbottom_c_index = " << fbottom_c_index << std::endl;
    out << "ftop_c_index = " << ftop_c_index << std::endl;
    out << "fbottom_side_orient = " << fbottom_side_orient << std::endl;
    out << "ftop_side_orient = " << ftop_side_orient << std::endl;
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

/*
template<class TSHAPE>
template<class TVar>
void TPZCompElHDivCollapsed<TSHAPE>::ComputeRequiredDataT(TPZMaterialDataT<TVar> &data,
                                                TPZVec<REAL> &qsi){
    
//    TPZManVector<int,TSHAPE::NSides*TSHAPE::Dimension> normalsidesDG(TSHAPE::Dimension*TSHAPE::NSides);

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = false;
    TPZCompElHDiv<TSHAPE>::ComputeRequiredData(data, qsi);
    data.fNeedsSol = needsol;
    TPZFNMatrix<9,REAL> axeslocal(TSHAPE::Dimension+1,3);
    ExpandAxes(data.axes, axeslocal);
    data.axes = axeslocal;
    const int dim = TSHAPE::Dimension+1;
    const int nvecshapestd = data.fDeformedDirections.Cols();
    TPZManVector<REAL> topdir(3,0.), botdir(3,0.); // top and bot directions in the deformed element
    TPZManVector<REAL,3> vecup={0,0,0}, vecdown={0,0,0};
    vecup[dim-1] = 1.;
    vecdown[dim-1] = -1.;
    // compute the deformed directions for the two additional vectors
    {
        for(int i=0; i<3; i++){
                topdir[i] = data.axes(dim-1,i);
                botdir[i] = -data.axes(dim-1,i);
        }
    }
    
    std::pair<TPZMaterialDataT<TVar>,TPZMaterialDataT<TVar>> *datapair = (std::pair<TPZMaterialDataT<TVar>,TPZMaterialDataT<TVar>> *) data.fUserData;
    TPZMaterialDataT<TVar> &datatop = datapair->second, &databottom = datapair->first;
    
    // compute the divergence of the top and bottom elements
    // the value is the value of the shape function times the sign of the vector in master direction
    {
        fTop.ComputeRequiredData(datatop, qsi);
        fBottom.ComputeRequiredData(databottom, qsi);
        int64_t numvec = data.divphi.Rows();
        int64_t numphi = data.phi.Rows();
        int64_t nvec_top = datatop.phi.Rows();
        int64_t nvec_bottom = databottom.phi.Rows();
        int64_t nvec_hdiv = numvec;
      
        // fDeformedDirections (for now) represents the H1 shape functions
        // times the element vectors. So, it is already the hdiv shape function itself.
        // Its size is, therefore, the size for a standard 2d hdiv element, plus
        // the shape functions related to the top and bottom connect that communicate
        // with the adjacent 3D elements
        const int64_t nvecshapecollpased = nvecshapestd+nvec_top+nvec_bottom;
        data.fDeformedDirections.Resize(3,nvecshapecollpased);
        data.fVecShapeIndex.Resize(nvecshapecollpased);
        data.divphi.Resize(nvecshapecollpased,1);
        data.phi.Resize(nvecshapecollpased,1);
        for(int i=numvec; i<nvecshapecollpased; i++)
        {
            data.fVecShapeIndex[i] = std::pair<int,int64_t>(i,i);
            data.phi(i) = 1.;
        }

        // First we append the bottom shapes and then the top shapes
        for (int64_t i= nvec_hdiv; i<nvecshapecollpased-nvec_top; i++) {
            for (int d = 0; d < 3; d++) {
                data.fDeformedDirections(d,i) = databottom.phi(i-nvec_hdiv)*botdir[d];
            }
        }
        for (int64_t i= nvec_hdiv+nvec_bottom; i<nvecshapecollpased; i++) {
            for (int d = 0; d < 3; d++) {
                data.fDeformedDirections(d,i) = datatop.phi(i-nvec_hdiv-nvec_bottom)*topdir[d];
            }
        }
        // Same for divphi
        for (int64_t i= nvec_hdiv; i<nvecshapecollpased-nvec_top; i++) {
            data.divphi(i,0) = -databottom.phi(i-nvec_hdiv);
        }
        for (int64_t i= nvec_hdiv+nvec_bottom; i<nvecshapecollpased; i++) {
            data.divphi(i,0) = datatop.phi(i-nvec_hdiv-nvec_bottom);
        }
    }
    
    if (data.fNeedsSol) {
        TPZCompElHDiv<TSHAPE>::ReallyComputeSolution(data);
    }


#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        data.fDeformedDirections.Print("Normal Vectors " , sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    

}//void
*/

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data) {


    // Some variable definitions
    TPZConnect& connbot = this->Mesh()->ConnectVec()[fbottom_c_index];
    TPZConnect& conntop = this->Mesh()->ConnectVec()[ftop_c_index];
    const int toporder = conntop.Order();
    const int bottomorder = connbot.Order();
    const int64_t nvec_top = conntop.NShape();
    const int64_t nvec_bottom = connbot.NShape();
    const int64_t nvecshapecollpased = this->NShapeF();
    const int64_t numvec = nvecshapecollpased - nvec_top - nvec_bottom;
    const int dim = TSHAPE::Dimension+1;
    
    // Expanding element axes
    TPZFNMatrix<9,REAL> axeslocal(TSHAPE::Dimension+1,3);
    ExpandAxes(data.axes, axeslocal);
    
    // adding a column do gradx pointing to the normal direction (direction that communicates with the 3d neighbor element)
    TPZFMatrix<REAL> gradx(3,TSHAPE::Dimension,0.);
    this->Reference()->GradX(qsi, gradx);
    TPZFNMatrix<9> gradxlocal(3,dim);
    for (int i=0; i<3; i++) {
        for (int d=0; d<dim-1; d++) {
            gradxlocal(i,d) = gradx(i,d);
        }
        gradxlocal(i,dim-1) = axeslocal(dim-1,i);
    }
    
    // Computing shape functions using TPZShapeHDivCollapsed structure
    TPZShapeData& shapedata = data;
    TPZShapeHDivCollapsed<TSHAPE>::Shape(qsi,shapedata,data.phi,data.divphi);
    data.divphi *= 1./data.detjac;
    gradxlocal.Multiply(data.phi,data.fDeformedDirections);
    data.fDeformedDirections *= 1./data.detjac;
    
    // Filling phi with 1 and fVecShapeIndex with 1,1 to make actual materials work
    data.fVecShapeIndex.Resize(nvecshapecollpased);
    data.phi.Resize(nvecshapecollpased,1);
    for(int i=0; i<nvecshapecollpased; i++)
    {
        data.fVecShapeIndex[i] = std::pair<int,int64_t>(i,i);
        data.phi(i) = 1.;
    }
}



template<class TSHAPE>
int64_t TPZCompElHDivCollapsed<TSHAPE>::ConnectIndex(int con) const
{
    if(con <= TSHAPE::NFacets) return TPZCompElHDiv<TSHAPE>::ConnectIndex(con);
    if(con > TSHAPE::NFacets + 2) DebugStop();
    if(con == TSHAPE::NFacets+1) return fbottom_c_index;
    if(con == TSHAPE::NFacets+2) return ftop_c_index;
    DebugStop();
    return -1;
}

/**
 * @brief Destroy internally allocated data structures
 */
//auto datapair = new std::pair<TPZMaterialDataT<TVar>,TPZMaterialDataT<TVar>>;
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::CleanupMaterialData(TPZMaterialData &data)
{
    TPZMaterialDataT<STATE> *dataS = dynamic_cast<TPZMaterialDataT<STATE> *>(&data);
    if(dataS)
    {
        std::pair<TPZMaterialDataT<STATE>,TPZMaterialDataT<STATE>> *userdataS = (std::pair<TPZMaterialDataT<STATE>,TPZMaterialDataT<STATE>> *) data.fUserData;
        if(userdataS) delete userdataS;
    }
    TPZMaterialDataT<CSTATE> *dataC = dynamic_cast<TPZMaterialDataT<CSTATE> *>(&data);
    if(dataC)
    {
        std::pair<TPZMaterialDataT<CSTATE>,TPZMaterialDataT<CSTATE>> *userdataC = (std::pair<TPZMaterialDataT<CSTATE>,TPZMaterialDataT<CSTATE>> *) data.fUserData;
        if(userdataC) delete userdataC;
    }
    data.fUserData = nullptr;
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::ConnectOrder(int connect) const {
	if (connect < 0 || connect >= NConnects()){
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "Connect index out of range connect " << connect <<
			" nconnects " << NConnects();
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return -1;
	}

    TPZConnect &c = this->Connect(connect);
    return c.Order();
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

