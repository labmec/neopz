//
//  TPBRWellBBox.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/27/15.
//
//

#include "TPBRWellBBox.h"
#include "tpzgeoblend.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoCube.h"
#include "TRMFlowConstants.h"
#include "TPZRefPatternTools.h"
#include "tpzarc3d.h"

#include "pzcheckgeom.h"
#include "tpzgeoblend.h"



/// add geometric nodes
void TPBRWellBBox::GenerateNodeLayers()
{
    REAL zouter[] = {fZBottom,fZBottom,fZTop,fZTop};
    REAL xouter[] = {-fLx/2.,fLx/2.,fLx/2.,-fLx/2.};
    
    REAL zwell[] = {-fWellDiameter/(2.*M_SQRT2),-fWellDiameter/(2.*M_SQRT2),fWellDiameter/(2.*M_SQRT2),fWellDiameter/(2.*M_SQRT2)};
    REAL xwell[] = {-fWellDiameter/(2.*M_SQRT2),fWellDiameter/(2.*M_SQRT2),fWellDiameter/(2.*M_SQRT2),-fWellDiameter/(2.*M_SQRT2)};
    
    REAL zarc[] = {-fWellDiameter/2.,0.,fWellDiameter/2.,0.};
    REAL xarc[] = {0.,fWellDiameter/2.,0.,-fWellDiameter/2.};
    
    int nlayers = fWellPosY.size();
    
    fGMesh->NodeVec().Resize(12*nlayers);
//    TPZGeoMesh *gmesh = fGMesh.operator->();
    int nodeindex = 0;
    
    fNodeLayers.Resize(nlayers);
    for (int ilayer = 0; ilayer<nlayers; ilayer++) {
        REAL youter = fOuterPosY[ilayer];
        REAL ywell = fWellPosY[ilayer];
        for (int inode = 0; inode<4; inode++) {
            TPZManVector<REAL,3> co(3);
            co[0] = xouter[inode];
            co[1] = youter;
            co[2] = zouter[inode];
            fGMesh->NodeVec()[nodeindex].Initialize(co, fGMesh);
            fNodeLayers[ilayer].fOuterNodeIndices[inode] = nodeindex;
            nodeindex++;
            co[0] = xwell[inode];
            co[1] = ywell;
            co[2] = zwell[inode];
            fGMesh->NodeVec()[nodeindex].Initialize(co, fGMesh);
            fNodeLayers[ilayer].fWellNodeIndices[inode] = nodeindex;
            nodeindex++;
            co[0] = xarc[inode];
            co[1] = ywell;
            co[2] = zarc[inode];
            fGMesh->NodeVec()[nodeindex].Initialize(co, fGMesh);
            fNodeLayers[ilayer].fArcNodes[inode] = nodeindex;
            nodeindex++;
        }
    }
}

/// generate 4 reservoir elements and one well element between both faces
void TPBRWellBBox::GenerateElementLayer(TFace &NodesLeft, TFace &NodesRight, bool createArc3DAtLeft, bool createArc3DAtRight)
{
//    int64_t elementid = 0;
    for (int volel=0; volel<4; volel++) {
        TPZGeoMesh *gmesh = fGMesh.operator->();
        
        TPZManVector<int64_t,8> nodeindices(8);
        nodeindices[0] = NodesLeft.fOuterNodeIndices[volel];
        nodeindices[1] = NodesLeft.fOuterNodeIndices[(volel+1)%4];
        nodeindices[2] = NodesLeft.fWellNodeIndices[(volel+1)%4];
        nodeindices[3] = NodesLeft.fWellNodeIndices[volel];
        nodeindices[4] = NodesRight.fOuterNodeIndices[volel];
        nodeindices[5] = NodesRight.fOuterNodeIndices[(volel+1)%4];
        nodeindices[6] = NodesRight.fWellNodeIndices[(volel+1)%4];
        nodeindices[7] = NodesRight.fWellNodeIndices[volel];
        new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(nodeindices,_ReservMatId,*gmesh);
        
    }
    {
        TPZGeoMesh *gmesh = fGMesh.operator->();
        
        TPZManVector<int64_t,8> nodeindices(8);
        nodeindices[0] = NodesLeft.fWellNodeIndices[0];
        nodeindices[1] = NodesLeft.fWellNodeIndices[1];
        nodeindices[2] = NodesLeft.fWellNodeIndices[2];
        nodeindices[3] = NodesLeft.fWellNodeIndices[3];
        nodeindices[4] = NodesRight.fWellNodeIndices[0];
        nodeindices[5] = NodesRight.fWellNodeIndices[1];
        nodeindices[6] = NodesRight.fWellNodeIndices[2];
        nodeindices[7] = NodesRight.fWellNodeIndices[3];
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(nodeindices,_WellMatId3D,*gmesh);
        
        TPZManVector<int64_t,3> arc3dIndices(3);
        if (createArc3DAtLeft)
        {
            arc3dIndices[0] = NodesLeft.fWellNodeIndices[0];
            arc3dIndices[1] = NodesLeft.fWellNodeIndices[1];
            arc3dIndices[2] = NodesLeft.fArcNodes[0];
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(arc3dIndices, _BlendArcMatId, *gmesh);
            
            arc3dIndices[0] = NodesLeft.fWellNodeIndices[1];
            arc3dIndices[1] = NodesLeft.fWellNodeIndices[2];
            arc3dIndices[2] = NodesLeft.fArcNodes[1];
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(arc3dIndices, _BlendArcMatId, *gmesh);
            
            arc3dIndices[0] = NodesLeft.fWellNodeIndices[2];
            arc3dIndices[1] = NodesLeft.fWellNodeIndices[3];
            arc3dIndices[2] = NodesLeft.fArcNodes[2];
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(arc3dIndices, _BlendArcMatId, *gmesh);
            
            arc3dIndices[0] = NodesLeft.fWellNodeIndices[3];
            arc3dIndices[1] = NodesLeft.fWellNodeIndices[0];
            arc3dIndices[2] = NodesLeft.fArcNodes[3];
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(arc3dIndices, _BlendArcMatId, *gmesh);
        }
        
        if(createArc3DAtRight)
        {
            arc3dIndices[0] = NodesRight.fWellNodeIndices[0];
            arc3dIndices[1] = NodesRight.fWellNodeIndices[1];
            arc3dIndices[2] = NodesRight.fArcNodes[0];
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(arc3dIndices, _BlendArcMatId, *gmesh);
            
            arc3dIndices[0] = NodesRight.fWellNodeIndices[1];
            arc3dIndices[1] = NodesRight.fWellNodeIndices[2];
            arc3dIndices[2] = NodesRight.fArcNodes[1];
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(arc3dIndices, _BlendArcMatId, *gmesh);
            
            arc3dIndices[0] = NodesRight.fWellNodeIndices[2];
            arc3dIndices[1] = NodesRight.fWellNodeIndices[3];
            arc3dIndices[2] = NodesRight.fArcNodes[2];
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(arc3dIndices, _BlendArcMatId, *gmesh);
            
            arc3dIndices[0] = NodesRight.fWellNodeIndices[3];
            arc3dIndices[1] = NodesRight.fWellNodeIndices[0];
            arc3dIndices[2] = NodesRight.fArcNodes[3];
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(arc3dIndices, _BlendArcMatId, *gmesh);
        }
    }
}

/// add to left and right elements
void TPBRWellBBox::AddElementClosure()
{
    {
        TFace &NodesLeft = fNodeLayers[0];
        TPZManVector<int64_t,8> nodeindices(8);
        nodeindices[0] = NodesLeft.fOuterNodeIndices[0];
        nodeindices[1] = NodesLeft.fOuterNodeIndices[1];
        nodeindices[2] = NodesLeft.fOuterNodeIndices[2];
        nodeindices[3] = NodesLeft.fOuterNodeIndices[3];
        nodeindices[4] = NodesLeft.fWellNodeIndices[0];
        nodeindices[5] = NodesLeft.fWellNodeIndices[1];
        nodeindices[6] = NodesLeft.fWellNodeIndices[2];
        nodeindices[7] = NodesLeft.fWellNodeIndices[3];
        TPZGeoMesh *gmesh = fGMesh.operator->();
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(nodeindices,_ReservMatId,*gmesh);
    }
    {
        TFace &NodesLeft = fNodeLayers[fNodeLayers.size()-1];
        TPZManVector<int64_t,8> nodeindices(8);
        nodeindices[0] = NodesLeft.fWellNodeIndices[0];
        nodeindices[1] = NodesLeft.fWellNodeIndices[1];
        nodeindices[2] = NodesLeft.fWellNodeIndices[2];
        nodeindices[3] = NodesLeft.fWellNodeIndices[3];
        nodeindices[4] = NodesLeft.fOuterNodeIndices[0];
        nodeindices[5] = NodesLeft.fOuterNodeIndices[1];
        nodeindices[6] = NodesLeft.fOuterNodeIndices[2];
        nodeindices[7] = NodesLeft.fOuterNodeIndices[3];
        TPZGeoMesh *gmesh = fGMesh.operator->();
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(nodeindices,_ReservMatId,*gmesh);
    }
}

/// Set the well partitions and the rib indices
void TPBRWellBBox::SetWellDivision(TPZVec<REAL> &position, int divisionstowardswell)
{
    fNumDiv = divisionstowardswell;
    fWellPosY = position;
    
    int sz = fWellPosY.size();
    fOuterPosY.Resize(sz, 0.);
    
    fOuterPosY[0] = -fLy/2.;
    fOuterPosY[sz-1] = +fLy/2.;
    
    for(int p = 1; p < sz-1; p++)
    {
        fOuterPosY[p] = fWellPosY[p];
    }
}

/// Generate the mesh based on the internal data structure
TPZAutoPointer<TPZGeoMesh> TPBRWellBBox::GenerateMesh()
{
    fGMesh = new TPZGeoMesh;
    GenerateNodeLayers();
    int nlayer = this->fNodeLayers.size()-1;
    
    int firstlayer = 0; // 0;
    for (int ilayer=firstlayer; ilayer < nlayer; ilayer++)
    {
        bool createArc3DAtRight = (ilayer == (nlayer-1));
        createArc3DAtRight = false;
        bool createArc3DAtLeft = true;
        if (ilayer == 0) {
            createArc3DAtLeft = false;
        }
        GenerateElementLayer(fNodeLayers[ilayer], fNodeLayers[ilayer+1],createArc3DAtLeft, createArc3DAtRight);
    }
    AddElementClosure();
    fGMesh->BuildConnectivity();
    
    this->ApplyDivisions();
    if (fPlaneZcoordinate != 0.)
    {
        this->CutMeshatPlane();
    }
    
    return fGMesh;
}

/// apply the successive refinements
void TPBRWellBBox::ApplyDivisions()
{
    std::set<int> matids;
    matids.insert(_WellMatId3D);
    for (int idiv = 0; idiv < fNumDiv; idiv++) {
        int nel = fGMesh->NElements();
        for (int iel = 0; iel<nel; iel++) {
            TPZGeoEl *gel = fGMesh->ElementVec()[iel];
            if (!gel || gel->HasSubElement()) {
                continue;
            }
            TPZRefPatternTools::RefineDirectional(gel, matids);
        }
    }
}

/**
 * This methos is used by RefineDirectional method!!!
 * Returns the refpattern that matches the sides refinement intensity and midnodes coordinates with respect to sidestorefine vector
 * @param gel - input data: geometric element for which the perfect match refpattern will be returned
 * @param sidestorefine - input data: vector filled with sides refinement intensity
 */
//static TPZAutoPointer<TPZRefPattern> PerfectMatchRefPattern(TPZGeoEl *gel, TPZVec<int> &sidestorefine);
/// determine at which height the nodes are which need to be moved
REAL TPBRWellBBox::HeightOfNearestPlane()
{
    if (fabs(fPlaneZcoordinate) < fWellDiameter) {
        DebugStop();
    }
    REAL mindist = fZTop-fZBottom;
    REAL bestheight = -1;
    int numpoints = 0;
    TPZStack<int> pointindices;
    /// look for elements which are divided and for which the plane crosses their edges
    int nnodes = fGMesh->NNodes();
    for (int in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        fGMesh->NodeVec()[in].GetCoordinates(xco);
        if (fabs(xco[2]) < fWellDiameter) {
            continue;
        }
        REAL dist = fabs(xco[2]-fPlaneZcoordinate);
        if (fabs(dist - mindist) < 1.e-6) {
            numpoints++;
            pointindices.Push(in);
        }
        if (dist < mindist) {
            mindist = dist;
            bestheight = xco[2];
            numpoints = 1;
            pointindices.Resize(0);
            pointindices.Push(in);
        }
    }
    
    std::cout << " Numpoints at mindist " << numpoints << " " << pointindices << std::endl;
    return bestheight;
}

/// move nodes at a certain height to a given height along their ribs
void TPBRWellBBox::MovePoints(REAL fromheight, REAL toheight)
{
    int nel = fGMesh->NElements();
    for (int el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->ElementVec()[el];
        if (!gel || !gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is = gel->NCornerNodes(); is < nsides; is++) {
            if (gel->SideDimension(is) != 1) continue;
            int64_t index = -1;
            gel->MidSideNodeIndex(is, index);
            if (index >= 0) {
                TPZManVector<REAL,3> xco(3);
                fGMesh->NodeVec()[index].GetCoordinates(xco);
                if (fabs(xco[2]-fromheight) < 1.e-6) {
                    MoveMidSideNode(gel, is);
                }
            }
            
        }
        
    }
}

/// divide elements whose edges cross the plane Z coordinate. return the set of newly created node indices
void TPBRWellBBox::DirectionalDivide(std::list<TPZGeoElSide> &dividedsides)
{
    int nelem = fGMesh->NElements();
    for (int el=0; el<nelem; el++) {
        TPZGeoEl *gel = fGMesh->ElementVec()[el];
        if (!gel ||gel->HasSubElement()) {
            continue;
        }
        std::list<TPZGeoElSide> gelsides;
        for (int is=0; is<gel->NSides(); is++) {
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            int nodindex1 = gel->SideNodeIndex(is, 0);
            int nodindex2 = gel->SideNodeIndex(is, 1);
            TPZManVector<REAL,3> co1(3),co2(3);
            fGMesh->NodeVec()[nodindex1].GetCoordinates(co1);
            fGMesh->NodeVec()[nodindex2].GetCoordinates(co2);
            if ((co1[2]-fPlaneZcoordinate)*(co2[2]-fPlaneZcoordinate) < -1.e-2) {
                gelsides.push_back(TPZGeoElSide(gel,is));
            }
        }
        if (gelsides.size()) {
            // identify the refinement pattern
            TPZManVector<int,30> sides(gel->NSides(),0);
            std::list<TPZGeoElSide>::iterator it = gelsides.begin();
            int i=0;
            while (it != gelsides.end()) {
                sides[it->Side()] = 1;
                dividedsides.push_back(*it);
                i++;
                it++;
            }
            TPZAutoPointer<TPZRefPattern> refpat = TPZRefPatternTools::PerfectMatchRefPattern(gel, sides);
            gel->SetRefPattern(refpat);
            TPZManVector<TPZGeoEl*,10> subgel;
            gel->Divide(subgel);
        }
    }
}

/// move the newly created nodes to the plane Z coordinate
void TPBRWellBBox::MoveNodes(std::set<int> &newnodes)
{
    int nel = fGMesh->NElements();
    for (int el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->ElementVec()[el];
        if (!gel) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is = gel->NCornerNodes(); is < nsides; is++) {
            if (gel->SideDimension(is) != 1) continue;
            int64_t index = -1;
            gel->MidSideNodeIndex(is, index);
            if (index >= 0) {
                if (newnodes.find(index) != newnodes.end()) {
                    MoveMidSideNode(gel, is);
                }
            }
            
        }
        
    }
}

/// move the midsidenode of the element to the fPlaneZcoordinate
void TPBRWellBBox::MoveMidSideNode(TPZGeoEl *gel, int side)
{
    int node1 = gel->SideNodeIndex(side, 0);
    int node2 = gel->SideNodeIndex(side, 1);
    int64_t midnode;
    gel->MidSideNodeIndex(side, midnode);
    
    if(midnode >= 0)
    {
        TPZManVector<REAL,3> co1(3),co2(3),coM(3);
        fGMesh->NodeVec()[node1].GetCoordinates(co1);
        fGMesh->NodeVec()[node2].GetCoordinates(co2);
        fGMesh->NodeVec()[midnode].GetCoordinates(coM);
        if (fabs(coM[2]-fPlaneZcoordinate) < 1.e-6)
        {
            return;
        }
        REAL factor = (fPlaneZcoordinate-co1[2])/(co2[2]-co1[2]);
        for (int i=0; i<3; i++)
        {
            coM[i] = co1[i] + factor*(co2[i]-co1[i]);
        }
        fGMesh->NodeVec()[midnode].SetCoord(coM);
    }
}

/// cut the elements at the height of fPlaneZcoordinate
void TPBRWellBBox::CutMeshatPlane()
{
    if (fPlaneZcoordinate == 0.) {
        return;
    }
    REAL closestplane;
    closestplane = HeightOfNearestPlane();
    MovePoints(closestplane, fPlaneZcoordinate);
    std::list<TPZGeoElSide> sidesdivided;
    DirectionalDivide(sidesdivided);
    std::list<TPZGeoElSide>::iterator it = sidesdivided.begin();
    while (it != sidesdivided.end()) {
        TPZGeoElSide gelside = *it;
        MoveMidSideNode(gelside.Element(), gelside.Side());
        it++;
    }
}


