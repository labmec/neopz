/**
 * @file
 * @brief Contains the implementation of the TPZGeoBlend methods. 
 */
#include "tpzgeoblend.h"

#include "pzgeoelside.h"
#include "tpzgeoelmapped.h"
#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"
#include "pzextractval.h"
#include "pzlog.h"
#include <sstream>

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.geoblend");
#endif

/** @brief Copy constructor */
template<class TGeo>
pzgeom::TPZGeoBlend<TGeo>::TPZGeoBlend(const TPZGeoBlend<TGeo> &cp, TPZGeoMesh &destmesh) : TPZRegisterClassId(&TPZGeoBlend<TGeo>::ClassId), TGeo(cp,destmesh)
{
    int64_t gelindex = cp.fGeoEl->Index();
    for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
        fNeighbours[is] =  cp.fNeighbours[is];
        fTrans[is] =  cp.fTrans[is];
    }
    fGeoEl = destmesh.Element(gelindex);
    if(!fGeoEl) DebugStop();
}


template<class TGeo>
bool pzgeom::TPZGeoBlend<TGeo>::IsLinearMapping(int side) const
{
    TPZStack<int> LowAllSides;
    TGeo::LowerDimensionSides(side,LowAllSides);
    if(side < 0 || side > TGeo::NSides-1)
    {
        DebugStop();
        return 0;
    }
    bool straight = true;
    
    if(side >= TGeo::NNodes && fNeighbours[side-TGeo::NNodes].ElementIndex() != -1)
    {
        straight = false;
        return straight;
    }
    for(int lowside = 0; lowside < LowAllSides.NElements(); lowside++)
    {
        if(LowAllSides[lowside] >= TGeo::NNodes && fNeighbours[LowAllSides[lowside]-TGeo::NNodes].ElementIndex() != -1)
        {
            straight = false;
            return straight;
        }
    }
    return straight;
}

template <class TGeo>
bool pzgeom::TPZGeoBlend<TGeo>::ResetBlendConnectivity(const int64_t &side, const int64_t &index){
    if(fNeighbours[side-TGeo::NNodes].ElementIndex() == index){
        fNeighbours[side-TGeo::NNodes].SetElementIndex(-1);
        return true;
    }
    return false;
}

template <class TGeo>
void pzgeom::TPZGeoBlend<TGeo>::SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform<> &trans)
{
	if(!(fNeighbours[side-TGeo::NNodes].ElementIndex() != -1))
	{
		fNeighbours[side-TGeo::NNodes] = neigh;
		fTrans[side - TGeo::NNodes] = trans;
	}
	else
	{
#ifdef PZ_LOG
        if(logger.isWarnEnabled())
        {
        
            std::stringstream mess;
            mess << "Trying to SetNeighbourInfo for an already set element\n";
            mess << "* this * = " << __PRETTY_FUNCTION__ << "\n";
            this->Print(mess);
            mess << "* neigh * = \n";
            neigh.Element()->Print(mess);
            LOGPZ_WARN(logger,mess.str());
        }
#endif
	}
}

template <class TGeo>
template<class T>
void pzgeom::TPZGeoBlend<TGeo>::GradX(TPZFMatrix<REAL> &coord, TPZVec<T> &xiInterior, TPZFMatrix<T> &gradx) const {

    TPZGeoEl &gel = *fGeoEl;
    TPZGeoMesh *gmesh = gel.Mesh();
#ifdef PZ_LOG
    const auto VAL_WIDTH = 10;
    std::ostringstream soutLogDebug;
    if(logger.isDebugEnabled())
    {   soutLogDebug<<"======================_______REF_1"<<std::endl;
        soutLogDebug << "TPZGeoBlend<" <<MElementType_Name(TGeo::Type())<<">::GradX"<<std::endl;
        soutLogDebug << "element id " <<gel.Id()<<std::endl;
        soutLogDebug << "xi: ";
        for(int i = 0; i < xiInterior.size(); i++) soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<xiInterior[i]<<"\t";
        soutLogDebug<<std::endl;
    }
#endif
    const REAL zero = 1e-14;
    gradx.Redim(3,TGeo::Dimension);

    if (fNeighbours[TGeo::NSides - 1 -TGeo::NNodes].ElementIndex() != -1){
        TPZManVector<T, 3> neighXi;
        TPZFNMatrix<9, T> dNeighXiDXi;
        if (!MapToNeighSide(TGeo::NSides-1, TGeo::Dimension, xiInterior, neighXi, dNeighXiDXi)) {
#ifdef PZ_LOG2
            if(logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "MapToNeighSide is singular for par " << xi << " and side " << TGeo::NSides-1 << ". Aborting...";
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            DebugStop();
        }
        TPZFNMatrix<9, T> gradNeigh;
        Neighbour(TGeo::NSides-1, gmesh).GradX(neighXi, gradNeigh);
        gradNeigh.Multiply(dNeighXiDXi,gradx);//gradNonLinSide = gradNeigh.dNeighXiDXi
        return;
    }

    /**
     * The non-linear mapping is calculated from deviations of the linear mapping.
     * The linear mapping of an element (or of any of its side) can be calculated with the barycentric coordinates
     * of the nodes contained in it
     */
    TPZFNMatrix<9, T> phi(TGeo::NNodes, 1), dPhiDxi(TGeo::Dimension, TGeo::NNodes);
    TGeo::TShape(xiInterior, phi, dPhiDxi);//gets the barycentric coordinates


    TPZFNMatrix<45,T> gradXLin(3, TGeo::Dimension,(T)0);
    for (int iNode = 0; iNode < TGeo::NNodes; iNode++) {//calculates the linear mapping
        for(int x = 0; x < 3; x++) {
            for (int xi = 0; xi < TGeo::Dimension; xi++) {
                gradXLin(x, xi) += coord(x, iNode) * dPhiDxi(xi, iNode);
            }
        }
    }
    for(int i = 0; i < gradx.Rows(); i++){
        for(int j = 0; j < gradx.Cols(); j++){
            gradx(i,j) = gradXLin(i,j);
        }
    }

#ifdef PZ_LOG
    soutLogDebug<<"FIRST TERM"<<std::endl;
    soutLogDebug << "gradient of linear mapping:\n";
    if (logger.isDebugEnabled()) {
        for (int i = 0; i < gradXLin.Rows(); i++) {
            for (int j = 0; j < gradXLin.Cols(); j++) {
                soutLogDebug << std::setw(VAL_WIDTH) << std::right << gradXLin(i, j);
                if (j != gradXLin.Cols() - 1) soutLogDebug << "\t";
            }
            soutLogDebug << std::endl;
        }
    }

#endif
    /**
     * Now, the deviation for any non-linearity of the sides' mappings must be taken into account.
     */

    TPZVec<TPZFMatrix<T> > gradNonLinSideVec(TGeo::NSides - TGeo::NNodes,
                                             TPZFNMatrix<27,T>(3, TGeo::Dimension, 0));
    TPZVec<TPZFMatrix<T> > gradLinSideVec(TGeo::NSides - TGeo::NNodes,
                                          TPZFNMatrix<27,T>(3, TGeo::Dimension, 0));
    TPZManVector<T, 20> blendFactor(TGeo::NSides - TGeo::NNodes, (T)0);
    TPZFNMatrix<27,T> dCorrFactorDxi(TGeo::NSides - TGeo::NNodes, TGeo::Dimension, (T) 0);
    TPZFNMatrix<27, T> linearSideMappings(TGeo::NSides - TGeo::NNodes, 3, 0.);
    TPZFNMatrix<27, T> nonLinearSideMappings(TGeo::NSides - TGeo::NNodes, 3, 0.);
    for (int sideIndex = 0; sideIndex < TGeo::NSides - TGeo::NNodes - 1; sideIndex++) {
        TPZManVector<T, 3> xiProjectedOverSide(TGeo::Dimension, 0);
        TPZFMatrix<T> &gradNonLinSide = gradNonLinSideVec[sideIndex];
        TPZFMatrix<T> &gradLinSide = gradLinSideVec[sideIndex];
        gradLinSide.Redim(3,TGeo::Dimension);
        gradNonLinSide.Redim(3,TGeo::Dimension);

        int side = TGeo::NNodes + sideIndex;
        TPZGeoElSide gelside(fNeighbours[sideIndex], gmesh);
        #ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            soutLogDebug << "================================" << std::endl;
            soutLogDebug << "side: " << side << " is linear: ";
        }
        #endif
        if (IsLinearMapping(side) || !gelside.Exists()) {
            blendFactor[sideIndex] = 0;
            #ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                if (IsLinearMapping(side)) soutLogDebug << "true" << std::endl;
                else soutLogDebug << "false (no gelside) " << std::endl;
            }
            #endif
            continue;
        }
        #ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            soutLogDebug << "false (gelside exists) " << std::endl;
        }
        #endif
        /**
     * Calculates the linear mapping of the side sideIndex, and the projected point on sideIndex
     */
        TPZManVector<T, 3> sideXi;
        TPZFNMatrix<9, T> transfXiToSideXi;
        bool regularMap = TGeo::CheckProjectionForSingularity(side, xiInterior);
        if (!regularMap) {
            #ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                soutLogDebug << "mapping is not regular. skipping side... ";

            }
            #endif
            continue;
        }

        this->MapToSide(side, xiInterior, sideXi, transfXiToSideXi);
        MElementType sideType = TGeo::Type(side);
        const int nSideNodes = MElementType_NNodes(sideType);
        const int sideDim = TGeo::SideDimension(side);
        TPZFNMatrix<9, T> sidePhiVec(nSideNodes, 1),
                dSidePhiDSideXiVec(sideDim, nSideNodes);
        GetSideShapeFunction<TGeo>(side, sideXi, sidePhiVec, dSidePhiDSideXiVec);

        TPZFNMatrix<9, T> dXprojDxiSide(TGeo::Dimension,sideDim,(T)0),
                dXiProjectedOverSideDxi(TGeo::Dimension,TGeo::Dimension,(T)0);
        TPZManVector<REAL, 3> nodeCoord;
        //calculation of transformation for the derivatives of sidephi
        
        TPZFNMatrix<9, T> gradLinSideXiSide(3,sideDim,(T)0);
        for (int iNode = 0; iNode < nSideNodes; iNode++) {
            const int currentNode = TGeo::SideNodeLocId(side, iNode);

            for (int x = 0; x < 3; x++) {
                linearSideMappings(sideIndex, x) +=
                        coord(x, currentNode) * sidePhiVec(iNode, 0);
            }

            TGeo::ParametricDomainNodeCoord(currentNode, nodeCoord);
            for (int x = 0; x < TGeo::Dimension; x++) {
                xiProjectedOverSide[x] += nodeCoord[x] * sidePhiVec(iNode, 0);
            }

            for(int x = 0; x < 3; x++){
                for(int xi = 0; xi < sideDim; xi++) {
                    gradLinSideXiSide(x,xi) += coord(x, currentNode) * dSidePhiDSideXiVec(xi,iNode);
                }
            }

            for (int i = 0; i < TGeo::Dimension; i++) {
                for (int j = 0; j < sideDim; j++) {
                    dXprojDxiSide(i,j) += nodeCoord[i] * dSidePhiDSideXiVec(j,iNode);
                }
            }

        }

        gradLinSideXiSide.Multiply(transfXiToSideXi,gradLinSide);//gradLinSideTemp = gradLinSideXiSide.transfXiToSideXi
        dXprojDxiSide.Multiply(transfXiToSideXi,dXiProjectedOverSideDxi);//dXprojDxiTemp = dXprojDxiSide.transfXiToSideXi

        #ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            soutLogDebug << "xi projection over side: ";
            for (int x = 0; x < TGeo::Dimension; x++) soutLogDebug << xiProjectedOverSide[x] << "\t";
            soutLogDebug << "\nGrad of projected point to side:\n";
            for(int i = 0; i < dXiProjectedOverSideDxi.Rows(); i++){
                for(int j = 0; j < dXiProjectedOverSideDxi.Cols(); j++){
                    soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<dXiProjectedOverSideDxi(i,j);
                    if(j != dXiProjectedOverSideDxi.Cols() - 1) soutLogDebug<<"\t";
                }
                soutLogDebug<<std::endl;
            }
            soutLogDebug << std::endl << "xi projection in side domain:";
            for (int x = 0; x < sideXi.size(); x++) soutLogDebug << sideXi[x] << "\t";
            soutLogDebug << "\nTransformation from x side to xi interior:\n";
            for(int i = 0; i < transfXiToSideXi.Rows(); i++){
                for(int j = 0; j < transfXiToSideXi.Cols(); j++){
                    soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<transfXiToSideXi(i,j);
                    if(j != transfXiToSideXi.Cols() - 1) soutLogDebug<<"\t";
                }
                soutLogDebug<<std::endl;
            }
            soutLogDebug << std::endl << "\nLinear mapping of projected point: ";
            for (int x = 0; x < 3; x++) soutLogDebug << linearSideMappings(sideIndex, x) << "\t";
            soutLogDebug << "\nGrad of linear mapping of point projected to side:\n";
            for(int i = 0; i < gradLinSide.Rows(); i++){
                for(int j = 0; j < gradLinSide.Cols(); j++){
                    soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<gradLinSide(i,j);
                    if(j != gradLinSide.Cols() - 1) soutLogDebug<<"\t";
                }
                soutLogDebug<<std::endl;
            }
        }
        #endif
        /**
         * Calculates the non-linear mapping of the side sideIndex
         */
        {
            TPZManVector<T,3> dCorrFactor(TGeo::Dimension,(T)0);
            TGeo::BlendFactorForSide(side, xiInterior, blendFactor[sideIndex], dCorrFactor);
            for(int iXi = 0; iXi < TGeo::Dimension; iXi++){
                dCorrFactorDxi(sideIndex,iXi) = dCorrFactor[iXi];
            }
        }

        int sidedim = gelside.Dimension();
        TPZManVector<T, 3> neighXi;
        TPZFNMatrix<9, T> dNeighXiDXi;
        if (!MapToNeighSide(side, sidedim, xiInterior, neighXi, dNeighXiDXi)) {
#ifdef PZ_LOG2
            if(logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "MapToNeighSide is singular for par " << par << " and side " << byside << " skipping the side ";
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            continue;
        }
        TPZManVector<T, 3> Xside(3, 0.);
        Neighbour(side, gmesh).X(neighXi, Xside);
        TPZFNMatrix<9, T> gradNeigh;
        Neighbour(side, gmesh).GradX(neighXi,gradNeigh);
        for (int x = 0; x < 3; x++) {
            nonLinearSideMappings(sideIndex, x) = Xside[x];
        }
        gradNeigh.Multiply(dNeighXiDXi,gradNonLinSide);//gradNonLinSide = gradNeigh.dNeighXiDXi

//        else{
        TPZStack<int> containedNodesInSide;
        TGeo::LowerDimensionSides(side, containedNodesInSide, 0);

        TPZStack<int> allContainedSides;
        TGeo::LowerDimensionSides(side, allContainedSides);
        for (int subSideIndex = containedNodesInSide.NElements();
             subSideIndex < allContainedSides.NElements(); subSideIndex++) {
//            TPZFNMatrix<9, T> gradNonLinSubSide(3,TGeo::Dimension,(T)0);
            const int subSide = allContainedSides[subSideIndex];
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                soutLogDebug << "\tSubside " << subSideIndex << " (global: " << subSide << ") is linear: ";
                if (IsLinearMapping(subSide)) soutLogDebug << "true" << std::endl;
                else soutLogDebug << "false" << std::endl;
            }
#endif
            if (IsLinearMapping(subSide)) continue;

            T blendFactorSide = -1;
            TPZManVector<T,3> dCorrFactorSideDxiProj(TGeo::Dimension,(T)0);
            TPZFNMatrix<3, T> dCorrFactorSideDxi(TGeo::Dimension,1,(T)0);
            TGeo::BlendFactorForSide(subSide, xiProjectedOverSide, blendFactorSide, dCorrFactorSideDxiProj);

#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                soutLogDebug << "\tSubside influence :" << blendFactorSide << std::endl;
            }
#endif
            TPZFNMatrix<3, T> dCorrFactorSideDxiProjMat(1,TGeo::Dimension,(T)0);
            for(int xi = 0; xi < TGeo::Dimension; xi++){
                dCorrFactorSideDxiProjMat(0,xi) = dCorrFactorSideDxiProj[xi];
            }

//                 dCorrFactorSideDxi = dXiProjectedOverSideDxi. dCorrFactorSideDxiProjMat
//                dXiProjectedOverSideDxi.Multiply(dCorrFactorSideDxiProjMat,dCorrFactorSideDxi);//mostrarphil

//                dCorrFactorSideDxi = dCorrFactorSideDxiProjMat . dXiProjectedOverSideDxi;
            dCorrFactorSideDxiProjMat.Multiply(dXiProjectedOverSideDxi,dCorrFactorSideDxi);
            #ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                soutLogDebug << "\n\tGrad of correction factor of point projected to side:\n";
                for(int i = 0; i < dCorrFactorSideDxi.Rows(); i++){
                    soutLogDebug<<"\t";
                    for(int j = 0; j < dCorrFactorSideDxi.Cols(); j++){
                        soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<dCorrFactorSideDxi(i,j);
                        if(j != dCorrFactorSideDxi.Cols() - 1) soutLogDebug<<"\t";
                    }
                    soutLogDebug<<std::endl;
                }
            }
            #endif
            TPZFMatrix<T> &gradNonLinSubSide = gradNonLinSideVec[subSide - TGeo::NNodes];
            TPZFMatrix<T> &gradLinSubSide = gradLinSideVec[subSide - TGeo::NNodes];
            for (int x = 0; x < 3; x++) {
                nonLinearSideMappings(sideIndex, x) -=
                        blendFactorSide *
                        (nonLinearSideMappings(subSide - TGeo::NNodes, x) -
                         linearSideMappings(subSide - TGeo::NNodes, x));
                for (int j = 0; j < TGeo::Dimension; j++) {
                    gradNonLinSide(x,j) -= blendFactorSide *  (gradNonLinSubSide(x,j)-gradLinSubSide(x,j));
                    gradNonLinSide(x,j) -= dCorrFactorSideDxi(0,j) *
                                           (nonLinearSideMappings(subSide - TGeo::NNodes, x) -
                                            linearSideMappings(subSide - TGeo::NNodes, x));
                }
            }
        }
//        }
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            soutLogDebug << "Non-linear mapping: ";
            for (int x = 0; x < 3; x++) soutLogDebug << nonLinearSideMappings(sideIndex, x) << "\t";
            soutLogDebug << "\nGrad of non-linear mapping:\n";
            for(int i = 0; i < gradNonLinSide.Rows(); i++){
                for(int j = 0; j < gradNonLinSide.Cols(); j++){
                    soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<gradNonLinSide(i,j);
                    if(j != gradNonLinSide.Cols() - 1) soutLogDebug<<"\t";
                }
                soutLogDebug<<std::endl;
            }
            soutLogDebug << std::endl;
            soutLogDebug << "adding to result mapping of side: " << side << std::endl;
            soutLogDebug << "\t\tcorrection factor: " << blendFactor[sideIndex] << std::endl;
            soutLogDebug << "\t\tdCorrFactorDxi: ";
            for(int i = 0; i < TGeo::Dimension; i++) soutLogDebug<<dCorrFactorDxi(sideIndex,i)<<"\t";
            soutLogDebug << "\n";

            soutLogDebug<<"SECOND TERM (SIDE "<<sideIndex<<") :"<<std::endl;
            soutLogDebug << "\t\tblendFactor *  (gradNonLinSide-gradLinSide): " << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < TGeo::Dimension; j++) {
                    const T val = blendFactor[sideIndex] *  (gradNonLinSide(i,j)-gradLinSide(i,j));
                    soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<val;
                    if(j != TGeo::Dimension - 1) soutLogDebug<<"\t";
                }
                soutLogDebug << "\n";
            }
            soutLogDebug<<"THIRD TERM (SIDE "<<sideIndex<<") :"<<std::endl;
            soutLogDebug << "\t\tdCorrFactorDxi *\n"
                            "                              (nonLinearSideMappings -\n"
                            "                               linearSideMappings): " << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < TGeo::Dimension; j++) {
                    const T val = dCorrFactorDxi(sideIndex,j) *
                                  (nonLinearSideMappings(sideIndex, i) -
                                   linearSideMappings(sideIndex,i));
                    soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<val;
                    if(j != TGeo::Dimension - 1) soutLogDebug<<"\t";
                }
                soutLogDebug << "\n";
            }

            soutLogDebug<<std::endl;
        }
#endif

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < TGeo::Dimension; j++) {
                gradx(i,j) += blendFactor[sideIndex] *  (gradNonLinSide(i,j)-gradLinSide(i,j));
                gradx(i,j) += dCorrFactorDxi(sideIndex,j) *
                              (nonLinearSideMappings(sideIndex, i) -
                               linearSideMappings(sideIndex,i));
            }
        }
    }
#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
        soutLogDebug << "================================" <<std::endl;
        soutLogDebug << "=============result=============" <<std::endl;
        soutLogDebug <<"================================"<<std::endl;
        for(int i = 0; i < gradx.Rows(); i++){
            for(int j = 0; j < gradx.Cols(); j++){
                soutLogDebug<<std::setw(VAL_WIDTH) << std::right<<gradx(i,j);
                if(j != gradx.Cols() - 1) soutLogDebug<<"\t";
            }
            soutLogDebug<<std::endl;
        }
        soutLogDebug << "================================" <<std::endl;
        soutLogDebug << "================================" <<std::endl;
        soutLogDebug << "================================" <<std::endl;
        LOGPZ_DEBUG(logger,soutLogDebug.str())
        soutLogDebug.str("");
    }
#endif
#ifdef PZDEBUG
    for(int i = 0; i < gradx.Rows();i++){
        for(int j = 0; j < gradx.Cols(); j++) {
            if (std::isnan(TPZExtractVal::val(gradx(i,j)))) {
                DebugStop();
            }
        }
    }
#endif

}

template<class TGeo>
template<class T>
void pzgeom::TPZGeoBlend<TGeo>::X(TPZFMatrix<REAL> &coord, TPZVec<T> &xi, TPZVec<T> &result) const {
//    TPZGeoEl &gel = *fGeoEl;
    if(!fGeoEl) DebugStop();
    TPZGeoMesh *gmesh = fGeoEl->Mesh();
    TPZManVector<T,3> notUsedHereVec(3,(T)0);
    TPZFNMatrix<9,T> notUsedHereMat(TGeo::NSides, TGeo::NSides,(T)0);//since some methods dont resize the matrix, it is
    // bigger than needed.

    #ifdef PZ_LOG
    std::ostringstream soutLogDebug;
    if(logger.isDebugEnabled())
    {
        soutLogDebug << "TPZGeoBlend<" <<MElementType_Name(TGeo::Type())<<">::X2_______REF_1"<<std::endl;
        soutLogDebug << "element id " <<fGeoEl->Id()<<std::endl;
        soutLogDebug << "xi: ";
        for(int i = 0; i < xi.size(); i++) soutLogDebug<<xi[i]<<"\n";
        soutLogDebug<<std::endl;
    }
    #endif
    const REAL zero = 1e-14;
    result.Resize(3);
    result.Fill(0);

    if (fNeighbours[TGeo::NSides - 1 -TGeo::NNodes].ElementIndex() != -1){
        TPZManVector<T, 3> neighXi;
        if (!MapToNeighSide(TGeo::NSides-1, TGeo::Dimension, xi, neighXi, notUsedHereMat)) {
#ifdef PZ_LOG2
            if(logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "MapToNeighSide is singular for par " << xi << " and side " << TGeo::NSides-1 << ". Aborting...";
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            DebugStop();
        }
        Neighbour(TGeo::NSides-1, gmesh).X(neighXi, result);
        return;
    }

    /**
     * The non-linear mapping is calculated from deviations of the linear mapping.
     * The linear mapping of an element (or of any of its side) can be calculated with the barycentric coordinates
     * of the nodes contained in it
     */
    TPZFNMatrix<9, T> phi(TGeo::NNodes, 1);
    TGeo::TShape(xi, phi, notUsedHereMat);//gets the barycentric coordinates


    for (int iNode = 0; iNode < TGeo::NNodes; iNode++) {//calculates the linear mapping
        for (int x = 0; x < 3; x++) {
            result[x] += coord(x, iNode) * phi(iNode, 0);
        }
    }
    #ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        soutLogDebug << "Linear mapping:\n"<<std::endl;
        for(int i = 0; i < result.size(); i++) soutLogDebug<<result[i]<<"\n";
        soutLogDebug<<std::endl;
    }
    #endif

    /**
     * Now, the deviation for any non-linearity of the sides' mappings must be taken into account.
     */
    TPZManVector<bool, 20> isRegularMapping(TGeo::NSides - TGeo::NNodes, 0.);
    for (int sideIndex = 0; sideIndex < TGeo::NSides - TGeo::NNodes - 1; sideIndex++) {
        int side = TGeo::NNodes + sideIndex;
        isRegularMapping[sideIndex] = TGeo::CheckProjectionForSingularity(side,xi);
    }

    TPZManVector<T, 20> blendFactor(TGeo::NSides - TGeo::NNodes, 0.);
    TPZFNMatrix<27, T> projectedPointOverSide(TGeo::NSides - TGeo::NNodes, TGeo::Dimension, 0.);
    TPZFNMatrix<27, T> linearSideMappings(TGeo::NSides - TGeo::NNodes, 3, 0.);
    TPZFNMatrix<27, T> nonLinearSideMappings(TGeo::NSides - TGeo::NNodes, 3, 0.);
    for (int sideIndex = 0; sideIndex < TGeo::NSides - TGeo::NNodes - 1; sideIndex++) {
        int side = TGeo::NNodes + sideIndex;
        if(!isRegularMapping[sideIndex]) {
            #ifdef PZ_LOG
            if(logger.isDebugEnabled()){
                soutLogDebug <<"mapping is not regular. skipping side... ";

            }
            #endif
            continue;
        }
        TPZGeoElSide gelside(fNeighbours[sideIndex], gmesh);
        #ifdef PZ_LOG
        if(logger.isDebugEnabled())
        {
            soutLogDebug << "================================" <<std::endl;
            soutLogDebug << "side: "<<side<<" is linear: ";
        }
        #endif
        if (IsLinearMapping(side) || !gelside.Exists()) {
            blendFactor[sideIndex] = 0;
            #ifdef PZ_LOG
            if(logger.isDebugEnabled()){
                if( IsLinearMapping(side) )  soutLogDebug <<"true"<<std::endl;
                else    soutLogDebug <<"false (no gelside) "<<std::endl;
            }
            #endif
            continue;
        }
        #ifdef PZ_LOG
        if(logger.isDebugEnabled()){
            soutLogDebug <<"false (gelside exists) "<<std::endl;
        }
        #endif
        /**
     * Calculates the linear mapping of the side sideIndex, and the projected point on sideIndex
     */
        TPZManVector<T, 3>  sideXi;
        this->MapToSide(side, xi, sideXi, notUsedHereMat);

        MElementType sideType = TGeo::Type(side);
        const int nSideNodes = MElementType_NNodes(sideType);
        TPZFNMatrix<9, T> sidePhi(nSideNodes, 1);
        notUsedHereMat.Resize(TGeo::Dimension,nSideNodes);
        GetSideShapeFunction<TGeo>(side, sideXi, sidePhi, notUsedHereMat);

        TPZManVector<REAL,3> nodeCoord;
        for (int iNode = 0; iNode < nSideNodes; iNode++) {
            const int currentNode = TGeo::SideNodeLocId(side, iNode);

            for (int x = 0; x < 3; x++) {
                linearSideMappings(sideIndex, x) +=
                        coord(x, currentNode) * sidePhi(iNode, 0);
            }
            TGeo::ParametricDomainNodeCoord(currentNode,nodeCoord);
            for (int x = 0; x < TGeo::Dimension; x++) {
                projectedPointOverSide(sideIndex,x) += nodeCoord[x] * sidePhi(iNode,0);
            }
        }
        /**
         * Calculates the non-linear mapping of the side sideIndex
         */
        #ifdef PZ_LOG
        if(logger.isDebugEnabled()){
            soutLogDebug <<"xi projection over side:\n";
            for (int x = 0; x < TGeo::Dimension; x++) soutLogDebug<<projectedPointOverSide(sideIndex,x)<<"\n";
            soutLogDebug<<std::endl<<"xi projection in side domain:\n";
            for (int x = 0; x < sideXi.size(); x++) soutLogDebug<<sideXi[x]<<"\n";
            soutLogDebug<<std::endl<<"linear mapping of projected point:\n";
            for (int x = 0; x < 3; x++) soutLogDebug<<linearSideMappings(sideIndex, x)<<"\n";
        }
        #endif

        TGeo::BlendFactorForSide(side, xi, blendFactor[sideIndex], notUsedHereVec);
        int sidedim = gelside.Dimension();
        TPZManVector<T, 3> neighXi;
        if (!MapToNeighSide(side, sidedim, xi, neighXi, notUsedHereMat)) {
#ifdef PZ_LOG2
            if(logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "MapToNeighSide is singular for par " << par << " and side " << byside << " skipping the side ";
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            continue;
        }
        TPZManVector<T, 3> Xside(3, 0.);
        TPZGeoElSide neighbour = Neighbour(side, gmesh);
        neighbour.X(neighXi, Xside);
        for (int x = 0; x < 3; x++) {
            nonLinearSideMappings(sideIndex, x) = Xside[x];
        }

//        else{
        TPZStack<int> containedNodesInSide;
        TGeo::LowerDimensionSides(side, containedNodesInSide, 0);

        TPZStack<int> allContainedSides;
        TGeo::LowerDimensionSides(side, allContainedSides);
        for (int subSideIndex = containedNodesInSide.NElements();
             subSideIndex < allContainedSides.NElements(); subSideIndex++) {
            const int subSide = allContainedSides[subSideIndex];
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                soutLogDebug << "\tSubside " << subSideIndex << " (global: " << subSide << ") is linear: ";
                if (IsLinearMapping(subSide)) soutLogDebug << "true" << std::endl;
                else soutLogDebug << "false" << std::endl;
            }
#endif
            if(!isRegularMapping[subSide - TGeo::NNodes]) {
                #ifdef PZ_LOG
                if(logger.isDebugEnabled()){
                    soutLogDebug <<"mapping of subside is not regular. skipping side... ";

                }
                #endif
                continue;
            }
            if (IsLinearMapping(subSide)) continue;
            TPZManVector<T, 3> projectedPoint(TGeo::Dimension, -1);
            for (int x = 0; x < TGeo::Dimension; x++) {
                projectedPoint[x] = projectedPointOverSide(sideIndex, x);
            }

            T blendFactorSide = -1;
            TGeo::BlendFactorForSide(subSide, projectedPoint, blendFactorSide,notUsedHereVec);

#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                soutLogDebug << "\tSubside influence :" << blendFactorSide << std::endl;
            }
#endif
            for (int x = 0; x < 3; x++) {
                nonLinearSideMappings(sideIndex, x) -=
                        blendFactorSide *
                        (nonLinearSideMappings(subSide - TGeo::NNodes, x) -
                         linearSideMappings(subSide - TGeo::NNodes, x));
            }

        }
//        }
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            soutLogDebug << "Non-linear mapping\n";
            for (int x = 0; x < 3; x++) soutLogDebug << nonLinearSideMappings(sideIndex, x) << "\n";
            soutLogDebug << std::endl;

            soutLogDebug << "adding to result mapping of side: " << side << std::endl;
            soutLogDebug << "\t\tblend factor: " << blendFactor[sideIndex] << std::endl;
//            LOGPZ_DEBUG(logger,soutLogDebug.str())
//            soutLogDebug.str("");
        }
#endif
        for (int x = 0; x < 3; x++) {
            result[x] += blendFactor[sideIndex] *
                         (nonLinearSideMappings(sideIndex, x) - linearSideMappings(sideIndex, x));
        }
    }
#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
        soutLogDebug << "================================" <<std::endl;
        soutLogDebug << "=============result=============" <<std::endl;
        soutLogDebug <<"================================"<<std::endl;
        for (int x = 0; x < 3; x++) soutLogDebug<<result[x]<<"\n";
        soutLogDebug << "\n================================" <<std::endl;
        soutLogDebug << "================================" <<std::endl;
        soutLogDebug << "================================" <<std::endl;
        LOGPZ_DEBUG(logger,soutLogDebug.str())
        soutLogDebug.str("");
    }
#endif
}

template <class TGeo>
void pzgeom::TPZGeoBlend<TGeo>::Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
{

    TPZGeoEl &gel = *fGeoEl;
    TPZManVector<REAL,3> NeighPar, SidePar, Xside(3,0.), XNode(3,0.);
    int majorSide = TGeo::NSides - 1;
	
    TPZManVector<REAL> SidesCounter(TGeo::NSides,0);
    TPZStack<int> LowNodeSides, LowAllSides;
	
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "input parameter par " << par;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
    TPZFNMatrix<24> blend(TGeo::NNodes,1), Dblend(TGeo::Dimension,TGeo::NNodes);
    TGeo::Shape(par,blend,Dblend);
	
    TPZFNMatrix<9> J1, J2, Ax, JacTemp(3,TGeo::Dimension, 0.), Jneighbourhood;
    REAL Det;
    TPZGeoMesh *gmesh = gel.Mesh();
    for(int byside = majorSide; byside >= TGeo::NNodes; byside--)
    {
        TPZGeoElSide neighbyside = Neighbour(byside, gmesh);
        if(neighbyside.Exists())
        {
            TGeo::LowerDimensionSides(byside,LowNodeSides,0);
            TGeo::LowerDimensionSides(byside,LowAllSides);
            int dim = Neighbour(byside, gmesh).Dimension();
            TPZFNMatrix<9> Inv(dim,dim);
            int sidedim = neighbyside.Dimension();
            if(!MapToNeighSide(byside,sidedim,par,NeighPar, Jneighbourhood))
			{
				continue;
			}
            Neighbour(byside,gmesh).X(NeighPar,Xside);
            Neighbour(byside,gmesh).Jacobian(NeighPar,J1,Ax,Det,Inv);
            Ax.Transpose(); 
            Ax.Multiply(J1,J2);
			
#ifdef PZ_LOG
			if(logger.isDebugEnabled())
			{
				std::stringstream sout;
				sout << "byside " << byside << std::endl;
        		sout << "side of the neighbour " << Neighbour(byside,gmesh).Side() << std::endl;
        		Neighbour(byside,gmesh).Element()->Print(sout);
				sout << "neighbour parameter(NeighPar) " << NeighPar << std::endl;
				sout << "Jacobian of the map(Jneighborhood) " << Jneighbourhood << std::endl;
				sout << "Xside " << Xside << std::endl;
				sout << "jacobian neighbour(J1) " << J1 << std::endl;
        		Ax.Print("Ax of the neighbour (Ax) ",sout);
				sout << "jacobian of the neighbour multiplied by the axes(J2) " << J2 << std::endl;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			
            J2.Multiply(Jneighbourhood,J1);
			
#ifdef PZ_LOG
			if(logger.isDebugEnabled())
			{
				std::stringstream sout;
				sout << "acumulated jacobian(J1) " << J1 << std::endl;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			
            REAL blendTemp = 0.; 
            TPZManVector<REAL,3> DblendTemp(TGeo::Dimension,0.);
            for(int a = 0; a < LowNodeSides.NElements(); a++)
            {
                TPZManVector<REAL> parChanged(par.NElements());
                TGeo::Shape(parChanged,blend,Dblend);
				
                blendTemp += blend(LowNodeSides[a],0);
                for(int b = 0; b < TGeo::Dimension; b++)
                {
                    DblendTemp[b] += Dblend(b,LowNodeSides[a]);
                }
            }
			
            for(int a = 0; a < 3; a++)
            {
                for(int b = 0; b < TGeo::Dimension; b++)
                {
                    JacTemp(a,b) += (1 - SidesCounter[byside]) * (J1(a,b)*blendTemp + Xside[a]*DblendTemp[b]);
                }
            }
            for(int a = 0; a < LowAllSides.NElements(); a++) 
            {
                SidesCounter[LowAllSides[a]] += (1 - SidesCounter[byside]);
            }
        }
    }
	
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		JacTemp.Print("Jabobian before contributing the nodes",sout);
		sout << "SidesCounter " << SidesCounter << std::endl;
		sout << "DBlend " << Dblend << std::endl;
		sout << "NodeCoord " << coord << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
    for(int a = 0; a < TGeo::NNodes; a++)
    {
        for(int b = 0; b < 3; b++) 
        { 
            for(int c = 0; c < TGeo::Dimension; c++)
            {
            	JacTemp(b,c) += (1 - SidesCounter[a]) * coord(b,a)*Dblend(c,a);
            }
        }
    }
	
    if(TGeo::Dimension == 1)
    {
        TPZFNMatrix<9> axest;
        JacTemp.GramSchmidt(axest,jacobian);
        axest.Transpose(&axes);
        jacinv.Resize(jacobian.Rows(), jacobian.Cols());
        detjac = jacobian(0,0);
        if(IsZero(detjac)){
            detjac = ZeroTolerance();
        }
        jacinv(0,0) = 1./detjac;
    } else if(TGeo::Dimension == 2)
    {
        TPZFNMatrix<9> axest;
        JacTemp.GramSchmidt(axest,jacobian);
        axest.Transpose(&axes);
        
        jacinv.Resize(jacobian.Rows(), jacobian.Cols());
        
        detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
        if(IsZero(detjac)){
            detjac = ZeroTolerance();
        }
        jacinv(0,0) =  jacobian(1,1) / detjac;
        jacinv(1,1) =  jacobian(0,0) / detjac;
        jacinv(0,1) = -jacobian(0,1) / detjac;
        jacinv(1,0) = -jacobian(1,0) / detjac;
    }
    else
    {
        jacobian = JacTemp;
        axes.Resize(3,3); axes.Zero();
        
        jacinv.Resize(jacobian.Rows(), jacobian.Cols());
        
        axes(0,0) = 1.; axes(1,1) = 1.; axes(2,2) = 1.;
        detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0);//- a02 a11 a20
        detjac += jacobian(0,1)*jacobian(1,2)*jacobian(2,0);//+ a01 a12 a20
        detjac += jacobian(0,2)*jacobian(1,0)*jacobian(2,1);//+ a02 a10 a21
        detjac -= jacobian(0,0)*jacobian(1,2)*jacobian(2,1);//- a00 a12 a21
        detjac -= jacobian(0,1)*jacobian(1,0)*jacobian(2,2);//- a01 a10 a22
        detjac += jacobian(0,0)*jacobian(1,1)*jacobian(2,2);//+ a00 a11 a22
        
        if(IsZero(detjac))
		{
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Singular Jacobian " << detjac;
                LOGPZ_ERROR(logger, sout.str())
            }
#endif
			detjac = ZeroTolerance();
		}
        
        jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2)) / detjac;//-a12 a21 + a11 a22
        jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2)) / detjac;// a02 a21 - a01 a22
        jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2)) / detjac;//-a02 a11 + a01 a12
        jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2)) / detjac;// a12 a20 - a10 a22
        jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2)) / detjac;//-a02 a20 + a00 a22
        jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2)) / detjac;// a02 a10 - a00 a12
        jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1)) / detjac;//-a11 a20 + a10 a21
        jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1)) / detjac;// a01 a20 - a00 a21
        jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1)) / detjac;//-a01 a10 + a00 a11
    }
}

template <class TGeo>
void pzgeom::TPZGeoBlend<TGeo>::Print(std::ostream &out) const
{
	TGeo::Print(out);
	out << "Neighbours/transformations used for mapping the sides :\n";
	int is;
	for(is=TGeo::NNodes; is<TGeo::NSides; is++)
	{
        out << "Side: " << is << " El/side: " << fNeighbours[is-TGeo::NNodes].ElementIndex() << ":" <<
        fNeighbours[is-TGeo::NNodes].Side() << '\n';
	}
}

template <class TGeo>
void pzgeom::TPZGeoBlend<TGeo>::Initialize(TPZGeoEl *refel)
{
    fGeoEl = refel;
    for(int byside = TGeo::NNodes; byside < (TGeo::NSides); byside++)
    {
        if (refel->SideIsUndefined(byside)) {
            continue;
        }
        TPZGeoElSide ElemSide(refel,byside);
        TPZGeoElSide NextSide(ElemSide.Neighbour());
        if(!NextSide.Element()) continue;
        while(NextSide.Element() != ElemSide.Element())
        {
            if(NextSide.Exists() && !NextSide.Element()->IsLinearMapping() && !NextSide.Element()->IsGeoBlendEl())
            {
                TPZGeoElSide NeighSide = NextSide;
                TPZTransform<> NeighTransf(NeighSide.Dimension(),NeighSide.Dimension());
                ElemSide.SideTransform3(NeighSide,NeighTransf);
                SetNeighbourInfo(byside,NeighSide,NeighTransf);
                break;
            }
            NextSide = NextSide.Neighbour();
        }
    }
}

/*
 //inseri este método pois o compilador não encontrou
 template <class TGeo>
 void TPZGeoBlend<TGeo>::Initialize(TPZVec<int> &nodeindexes)
 {
 TGeo::Initialize(nodeindexes);
 }
 */

#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeoprism.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzcmesh.h"


using namespace pzgeom;


/// create an example element based on the topology
/* @param gmesh mesh in which the element should be inserted
 @param matid material id of the element
 @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
 @param size (in) size of space where the element should be created
 */

#include "TPZWavyLine.h"

template <class TGeo>
void pzgeom::TPZGeoBlend<TGeo>::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    
    TGeo::InsertExampleElement(gmesh, -1, lowercorner, size);
    int64_t elid = gmesh.ElementVec().NElements()-1;
    TPZManVector<int64_t,3> nodeindexes(8);
    TPZGeoEl *gel = gmesh.Element(elid);
    int NNodes = TGeo::NCornerNodes;
    for (int i=0; i<NNodes; i++) {
        nodeindexes[i] = gel->NodeIndex(i);
    }
    TPZGeoElRefPattern<TPZWavyLine> *gelwave = new TPZGeoElRefPattern<TPZWavyLine>(nodeindexes,matid,gmesh);
    TPZManVector<REAL,3> wavedir(3,0.02);
    wavedir[0] = 0.;
    gelwave->Geom().SetData(wavedir, 2);
    delete gel;
    int64_t index;
    gmesh.CreateGeoBlendElement(TGeo::Type(), nodeindexes, matid, index);
}




template class pzgeom::TPZGeoBlend<TPZGeoCube>;
template class pzgeom::TPZGeoBlend<TPZGeoTriangle>;
template class pzgeom::TPZGeoBlend<TPZGeoPrism>;
template class pzgeom::TPZGeoBlend<TPZGeoPyramid>;
template class pzgeom::TPZGeoBlend<TPZGeoTetrahedra>;
template class pzgeom::TPZGeoBlend<TPZGeoQuad>;
template class pzgeom::TPZGeoBlend<TPZGeoLinear>;
template class pzgeom::TPZGeoBlend<TPZGeoPoint>;


///CreateGeoElement -> TPZGeoBlend
#define IMPLEMENTBLEND(TGEO,CREATEFUNCTION) \
\
template void pzgeom::TPZGeoBlend<TGEO>::X(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZVec<REAL> &result) const; \
template void pzgeom::TPZGeoBlend<TGEO>::X(TPZFMatrix<REAL> &coord, TPZVec<Fad<REAL>> &par, TPZVec<Fad<REAL>> &result) const; \
template void pzgeom::TPZGeoBlend<TGEO>::GradX(TPZFMatrix<REAL> &coord, TPZVec<REAL> &xiInterior, TPZFMatrix<REAL> &gradx) const;\
template void pzgeom::TPZGeoBlend<TGEO>::GradX(TPZFMatrix<REAL> &coord, TPZVec<Fad<REAL>> &xiInterior, TPZFMatrix<Fad<REAL>> &gradx) const;\
template class \
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoBlend<TGEO> >>; \
\
template class TPZGeoElRefLess<TPZGeoBlend<TGEO> >;\
template class TPZGeoElRefPattern<TPZGeoBlend<TGEO> >;

IMPLEMENTBLEND(pzgeom::TPZGeoPoint,CreatePointEl)
IMPLEMENTBLEND(pzgeom::TPZGeoLinear,CreateLinearEl)
IMPLEMENTBLEND(pzgeom::TPZGeoQuad,CreateQuadEl)
IMPLEMENTBLEND(pzgeom::TPZGeoTriangle,CreateTriangleEl)
IMPLEMENTBLEND(pzgeom::TPZGeoCube,CreateCubeEl)
IMPLEMENTBLEND(pzgeom::TPZGeoPrism,CreatePrismEl)
IMPLEMENTBLEND(pzgeom::TPZGeoPyramid,CreatePyramEl)
IMPLEMENTBLEND(pzgeom::TPZGeoTetrahedra,CreateTetraEl)

#undef IMPLEMENTBLEND
