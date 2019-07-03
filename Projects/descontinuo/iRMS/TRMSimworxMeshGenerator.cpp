//
//  TRMSimworxMeshGenerator.cpp
//  PZ
//
//  Created by Nathan Shauer on 6/15/15.
//
//

#include "TRMSimworxMeshGenerator.h"
#include "TRMRawData.h"
#include "TRMFlowConstants.h"
#include "TPBRWellBBox.h"

#include "TPZVTKGeoMesh.h"
#include "TPZRefPatternTools.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "tpzellipse3d.h"
#include "pzgeoelbc.h"

const int arcEllip = -1;
const int mioloMat  = 1;
const int reservMat = 2;
const int reservTriangleMat = -2;
const int reservMatTemp = -3;

StructMioloData::StructMioloData()
{
    m_Lx = 0.;
    m_Ly = 0.;
    m_LzBottom = 0.;
    m_LzTop = 0.;
    m_thereIsCutPlane = false;
    m_CutPlaneZ = 0.;
    m_wellPositionsY.Resize(0);
    m_ribsIndexes.Resize(0);
    m_reservInterfaceNodeIndices.clear();
}

TRMSimworxMeshGenerator::TRMSimworxMeshGenerator()
{
    m_auxGMesh = NULL;
    fRibs.clear();
}

TRMSimworxMeshGenerator::~TRMSimworxMeshGenerator()
{
    m_auxGMesh = NULL;
    delete m_auxGMesh;
}

void TRMSimworxMeshGenerator::FillLinerHolesData(TPZFMatrix<> &linerHolesData)
{
    linerHolesData.Redim(5, 3);
    linerHolesData(0,0) = -250.;
    linerHolesData(1,0) = -100.;
    linerHolesData(2,0) = -70.;
    linerHolesData(3,0) = 0.;
    linerHolesData(4,0) = 250.;
    
    linerHolesData(0,1) = 0.;
    linerHolesData(1,1) = 164.042;
    linerHolesData(2,1) = 0.;
    linerHolesData(3,1) = 65.6168;
    linerHolesData(4,1) = 131.234;
    
    linerHolesData(0,2) = 0.;
    linerHolesData(1,2) = 0.00635;
    linerHolesData(2,2) = 0.;
    linerHolesData(3,2) = 0.00635;
    linerHolesData(4,2) = 0.003175;
    
}

void TRMSimworxMeshGenerator::FillCaseHolesData(TPZFMatrix<> &caseHolesData)
{
    caseHolesData.Redim(4, 3);
    caseHolesData(0,0) = -250.;
    caseHolesData(1,0) = -150.;
    caseHolesData(2,0) = -50.;
    caseHolesData(3,0) = 250.;
    
    caseHolesData(0,1) = 0.;
    caseHolesData(1,1) = 65.6168;
    caseHolesData(2,1) = 0.;
    caseHolesData(3,1) = 65.6168;
    
    caseHolesData(0,2) = 0.00635;
    caseHolesData(1,2) = 0.00635;
    caseHolesData(2,2) = 0.00635;
    caseHolesData(3,2) = 0.00635;
    
}

TPZGeoMesh * TRMSimworxMeshGenerator::CreateSimworxGeoMesh(TRMRawData &rawdata, bool withwellbc)
{
    const REAL reservoirWidth = rawdata.fReservoirWidth;
    const REAL reservoirLength = rawdata.fReservoirLength;
    const REAL MyreservoirHeight = rawdata.fReservoirHeight;
    const REAL MywellLength = rawdata.fLw;
    
    const REAL reservoirSemiAxeX = reservoirWidth/2.;
    const REAL reservoirSemiAxeY = reservoirLength/2.;
    const REAL mioloWidth = reservoirWidth/5.;
    const REAL reservoirheight = MyreservoirHeight;
    
    const int ndiv = 12;
    TPZManVector<REAL,ndiv> espacamentoReservY(ndiv);
    
    REAL prop = 0.70;
    REAL wellength = MywellLength;
    espacamentoReservY[0] = -prop*wellength;
    espacamentoReservY[ndiv-1] = prop *wellength;
    for (int i=1; i<ndiv-1; i++) {
        // 8 -> 1.397
        // 12 -> 1.2
        espacamentoReservY[i] = 1.25*(-wellength/2.+i*(wellength/(ndiv-1)));
    }
    espacamentoReservY[1] = -wellength/2.+0.3;
    espacamentoReservY[ndiv-2] = wellength/2.-0.3;
    
    TPZManVector<REAL,6> espacamentoZ(5);
    espacamentoZ[0] = -reservoirheight/2.;
    espacamentoZ[1] = -20.;
    espacamentoZ[2] = -5.;
    espacamentoZ[3] = +20.;
    espacamentoZ[4] = reservoirheight/2.;
    bool thereIsCutPlane = true;
    
    TPZGeoMesh * reservoirGMesh = ReallyGenerateGeoMesh(reservoirSemiAxeX,reservoirSemiAxeY,mioloWidth,espacamentoReservY,espacamentoZ,thereIsCutPlane,rawdata);
    if (withwellbc) {
        CreateWellBoundaries(reservoirGMesh);
    }
    
    return reservoirGMesh;
}

TPZGeoMesh * TRMSimworxMeshGenerator::ReallyGenerateGeoMesh(const REAL semiAxeX,
                                                 const REAL semiAxeY,
                                                 const REAL mioloLx,
                                                 const TPZVec<REAL> & espacamentoMioloY,
                                                 const TPZVec<REAL> & espacamentoZ,
                                                 const bool thereIsCutPlane,
                                                 TRMRawData &rawdata)
{
    REAL Lw = rawdata.fLw;
    
    //verificacoes
    {
        if(espacamentoMioloY.NElements() < 2)
        {
            DebugStop();
        }
        if(espacamentoZ.NElements() < 2)
        {
            DebugStop();
        }
        REAL miolo1st = espacamentoMioloY[0];
        REAL mioloLast = espacamentoMioloY[espacamentoMioloY.NElements()-1];
        if(fabs(miolo1st+mioloLast) > 1.e-3 || miolo1st > mioloLast)//Velificando se (miolo1st != -mioloLast || miolo1st > mioloLast)
        {
            DebugStop();
        }
        
        if(espacamentoMioloY.NElements() > 2)
        {
            REAL miolo2nd = espacamentoMioloY[1];
            REAL mioloPenultimate = espacamentoMioloY[espacamentoMioloY.NElements()-2];
            
            if(fabs(miolo2nd) > Lw/2.)//Somente o miolo1st e mioloLast devem estar alem da extensao do poco.
            {
                DebugStop();
            }
            if(fabs(mioloPenultimate) > Lw/2.)
            {
                DebugStop();
            }
        }
    }
    
    StructMioloData mioloData;
    
    //Reservoir GeoMesh
    TPZGeoMesh * reservoirGMesh = CreateEllipticalReservoirGeoMesh(semiAxeX,
                                                                   semiAxeY,
                                                                   mioloLx,
                                                                   espacamentoMioloY,
                                                                   espacamentoZ,
                                                                   Lw,
                                                                   thereIsCutPlane,
                                                                   mioloData,
                                                                   rawdata);
    
    
    //Miolo GeoMesh
    REAL welldiam = rawdata.fWellDiam;
    REAL Lx = mioloData.m_Lx;
    REAL Ly = mioloData.m_Ly;
    REAL ZBottom = -mioloData.m_LzBottom;
    REAL ZTop = mioloData.m_LzTop;
    
   
    TPBRWellBBox box(welldiam,Lx,Ly,ZBottom,ZTop);
    int numdiv = 8;
    
    box.SetWellDivision(mioloData.m_wellPositionsY, numdiv);
    if(thereIsCutPlane)
    {
        box.SetCutPlane(mioloData.m_CutPlaneZ);
    }
    
    
    TPZAutoPointer<TPZGeoMesh> boxGMesh = box.GenerateMesh();
    
    
    //Matching nodes between reserv and miolo GeoMeshes
    std::map<int,int> miolo_reserv_nodeIndices;
    Pair_Miolo_Reserv_Nodes(reservoirGMesh, mioloData.m_reservInterfaceNodeIndices,
                            boxGMesh, miolo_reserv_nodeIndices);
    
    
    TPZGeoMesh * completeGMesh = MergeGeoMeshes(reservoirGMesh,boxGMesh,miolo_reserv_nodeIndices);
    
    HangFacesTreatment(completeGMesh);
    return completeGMesh;
}

TPZGeoMesh * TRMSimworxMeshGenerator::CreateEllipticalReservoirGeoMesh(const REAL semiAxeX,
                                              const REAL semiAxeY,
                                              const REAL mioloLx,
                                              const TPZVec<REAL> & espacamentoMioloY,
                                              const TPZVec<REAL> & espacamentoZ,
                                              const REAL Lw,
                                              const bool thereIsCutPlane,
                                              StructMioloData & mioloData,
                                              TRMRawData &rawdata)
{
    //Construindo malha eliptica auxiliar (para mapear os pontos da malha cartesiana para a elipse)
    REAL innerRectangleLx, innerRectangleLy;
    BuildAuxiliary2DGeoMesh(semiAxeX,semiAxeY,
                            innerRectangleLx, innerRectangleLy);
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    int nnodesX = 4;
    int nnodesY = 2 + espacamentoMioloY.NElements();
    int nnodes = nnodesX*nnodesY;
    
    gmesh->NodeVec().Resize(nnodes);
    
    int actNode = 0;
    for(int ny = 0; ny < nnodesY; ny++)
    {
        for(int nx = 0; nx < nnodesX; nx++)
        {
            TPZManVector<REAL,3> coordIn(3,0.);
            TPZManVector<REAL,3> coordOut(3,0.);
            
            ///Abscissa
            if(nx == 0)
            {
                coordIn[0] = -innerRectangleLx/2.;
            }
            else if(nx == 1)
            {
                coordIn[0] = -mioloLx/2.;
            }
            else if(nx == 2)
            {
                coordIn[0] = +mioloLx/2.;
            }
            else if(nx == 3)
            {
                coordIn[0] = +innerRectangleLx/2.;
            }
            else
            {
                DebugStop();
            }
            
            ///Ordenada
            if(ny == 0)
            {
                coordIn[1] = -innerRectangleLy/2.;
            }
            else if(ny == nnodesY-1)
            {
                coordIn[1] = +innerRectangleLy/2.;
            }
            else
            {
                coordIn[1] = espacamentoMioloY[ny-1];
            }
            
            if(IsNodeOn2DCentralCore(nnodes,actNode) == false)
            {
                MapCoordIn2DEllipse(coordIn,innerRectangleLx,innerRectangleLy,coordOut);
            }
            else
            {
                coordOut = coordIn;
            }
            
            coordOut[2] = espacamentoZ[0];
            gmesh->NodeVec()[actNode].SetCoord(coordOut);
            actNode++;
        }//nx
    }//ny
    
    for(int ny = 0; ny < nnodesY-1; ny++)
    {
        if(ny == 0)
        { //Criacao de elementos abaixo
            TPZManVector<int64_t,3> Topol3(3);
            TPZManVector<int64_t,4> Topol4(4);
            
            Topol3[0] = 4; Topol3[1] = 0; Topol3[2] = 5;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(Topol3, reservTriangleMat, *(gmesh));
            
            Topol3[0] = 1; Topol3[1] = 0; Topol3[2] = 5;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(Topol3, reservTriangleMat, *(gmesh));
            
            Topol4[0] = 2; Topol4[1] = 6; Topol4[2] = 5; Topol4[3] = 1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(Topol4, reservMatTemp, *(gmesh));
            
            Topol3[0] = 2; Topol3[1] = 3; Topol3[2] = 6;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(Topol3, reservTriangleMat, *(gmesh));
            Topol3[0] = 7; Topol3[1] = 3; Topol3[2] = 6;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(Topol3, reservTriangleMat, *(gmesh));
        }
        else if(ny == nnodesY-2)
        { //Criacao de elementos acima
            TPZManVector<int64_t,3> Topol3(3);
            TPZManVector<int64_t,4> Topol4(4);
            
            Topol3[0] = nnodes-8; Topol3[1] = nnodes-4; Topol3[2] = nnodes-7;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(Topol3, reservTriangleMat, *(gmesh));
            
            Topol3[0] = nnodes-3; Topol3[1] = nnodes-4; Topol3[2] = nnodes-7;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(Topol3, reservTriangleMat, *(gmesh));
            
            Topol4[0] = nnodes-3; Topol4[1] = nnodes-7; Topol4[2] = nnodes-6; Topol4[3] = nnodes-2;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(Topol4, reservMatTemp, *(gmesh));
            
            Topol3[0] = nnodes-2; Topol3[1] = nnodes-1; Topol3[2] = nnodes-6;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(Topol3, reservTriangleMat, *(gmesh));
            
            Topol3[0] = nnodes-5; Topol3[1] = nnodes-1; Topol3[2] = nnodes-6;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(Topol3, reservTriangleMat, *(gmesh));
        }
        else
        {//Criacao de elementos intermediarios (miolo e laterais)
            TPZManVector<int64_t,4> Topol4(4);
            
            Topol4[0] = ny*nnodesX+0; Topol4[1] = ny*nnodesX+1; Topol4[2] = (ny+1)*nnodesX+1; Topol4[3] = (ny+1)*nnodesX+0;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(Topol4, reservMat, *(gmesh));
            
            Topol4[0] = ny*nnodesX+1; Topol4[1] = ny*nnodesX+2; Topol4[2] = (ny+1)*nnodesX+2; Topol4[3] = (ny+1)*nnodesX+1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(Topol4, mioloMat, *(gmesh));
            
            Topol4[0] = ny*nnodesX+2; Topol4[1] = ny*nnodesX+3; Topol4[2] = (ny+1)*nnodesX+3; Topol4[3] = (ny+1)*nnodesX+2;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(Topol4, reservMat, *(gmesh));
        }
    }
    
    gmesh->BuildConnectivity();
    AdjustPrismCoordinates(gmesh,semiAxeX,semiAxeY);
    
    InsertEllipseArcs(gmesh, semiAxeX, semiAxeY, espacamentoZ[0]);
    
    REAL mioloLast = espacamentoMioloY[espacamentoMioloY.NElements()-1];
    int ndiv = std::max(1.,log((semiAxeY-mioloLast)/mioloLx)/log(2.));
    
    gmesh = ElementsCrusher2D(gmesh,ndiv);
    gmesh = Extrude2DMesh(gmesh,espacamentoZ,thereIsCutPlane);
    
    FillStructMiolo(mioloLx, espacamentoMioloY, espacamentoZ, Lw, thereIsCutPlane, mioloData, rawdata);
    CreateReservoirFaces(*(gmesh),mioloData);
    
    return gmesh;
    
}

void TRMSimworxMeshGenerator::BuildAuxiliary2DGeoMesh(const REAL semiX, const REAL semiY,
                             REAL & innerRectangleLx, REAL & innerRectangleLy)
{
    innerRectangleLx = 2.*semiX/sqrt(2.);
    innerRectangleLy = 2.*semiY/sqrt(2.);
    
    
    
    //////////////////////////////////////////////////
    m_auxGMesh = new TPZGeoMesh;
    m_auxGMesh->NodeVec().Resize(4);
    TPZManVector<REAL,3> NodesCoords(3,0.);
    
    m_auxGMesh->NodeVec()[0].SetNodeId(0);
    NodesCoords[0] = -innerRectangleLx/2.;
    NodesCoords[1] = -innerRectangleLy/2.;
    m_auxGMesh->NodeVec()[0].SetCoord(NodesCoords);
    
    m_auxGMesh->NodeVec()[1].SetNodeId(1);
    NodesCoords[0] = +innerRectangleLx/2.;
    NodesCoords[1] = -innerRectangleLy/2.;
    m_auxGMesh->NodeVec()[1].SetCoord(NodesCoords);
    
    m_auxGMesh->NodeVec()[2].SetNodeId(2);
    NodesCoords[0] = +innerRectangleLx/2.;
    NodesCoords[1] = +innerRectangleLy/2.;
    m_auxGMesh->NodeVec()[2].SetCoord(NodesCoords);
    
    m_auxGMesh->NodeVec()[3].SetNodeId(3);
    NodesCoords[0] = -innerRectangleLx/2.;
    NodesCoords[1] = +innerRectangleLy/2.;
    m_auxGMesh->NodeVec()[3].SetCoord(NodesCoords);
    
    TPZManVector<int64_t,4> Topol(4);
    for(int i = 0; i < 4; i++) Topol[i] = i;
    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(Topol, 0, *(m_auxGMesh));
    
    TPZManVector<REAL,3> ellipOrigin(3,0.);
    TPZManVector<REAL,3> semiAxeX(3,0.);
    TPZManVector<REAL,3> semiAxeY(3,0.);
    semiAxeX[0] = semiX; semiAxeY[1] = semiY;
    
    Topol.Resize(2);
    
    Topol[0] = 0;
    Topol[1] = 1;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge4 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,1, *(m_auxGMesh));
    ellipEdge4->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    Topol[0] = 1;
    Topol[1] = 2;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge5 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,1, *(m_auxGMesh));
    ellipEdge5->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    Topol[0] = 2;
    Topol[1] = 3;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge6 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,1, *(m_auxGMesh));
    ellipEdge6->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    Topol[0] = 3;
    Topol[1] = 0;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge7 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,1, *(m_auxGMesh));
    ellipEdge7->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    m_auxGMesh->BuildConnectivity();
}

void TRMSimworxMeshGenerator::Pair_Miolo_Reserv_Nodes(TPZGeoMesh * reservGMesh, const std::set<int> & reservNodeIndices,
                             TPZAutoPointer<TPZGeoMesh> mioloGMesh,
                             std::map<int,int> & miolo_reserv_nodeIndices)
{
    miolo_reserv_nodeIndices.clear();
    
    std::set<int>::const_iterator it_miolo, it_reserv;
    
    for(int nodeMiolo = 0; nodeMiolo < mioloGMesh->NNodes(); nodeMiolo++)
    {
        TPZManVector<REAL,3> mioloCoord(3,0.);
        mioloGMesh->NodeVec()[nodeMiolo].GetCoordinates(mioloCoord);
        
        bool nodeMatched = false;
        for(it_reserv = reservNodeIndices.begin(); it_reserv != reservNodeIndices.end(); it_reserv++)
        {
            int nodeReserv = *it_reserv;
            TPZManVector<REAL,3> reservCoord(3,0.);
            reservGMesh->NodeVec()[nodeReserv].GetCoordinates(reservCoord);
            
            if(DistFromPoints(mioloCoord,reservCoord) < 1.e-3)
            { //nohs pareados
                nodeMatched = true;
                miolo_reserv_nodeIndices[nodeMiolo] = nodeReserv;
                break;
            }
        }
        if(nodeMatched == false)
        { //criacao de nohs na malha reservGMesh, seguido do pareamento
            int lastNode = reservGMesh->NodeVec().NElements();
            reservGMesh->NodeVec().Resize(lastNode+1);
            reservGMesh->NodeVec()[lastNode].SetNodeId(lastNode);
            reservGMesh->NodeVec()[lastNode].SetCoord(mioloCoord);
            miolo_reserv_nodeIndices[nodeMiolo] = lastNode;
        }
    }
}

REAL TRMSimworxMeshGenerator::DistFromPoints(TPZVec<REAL> & ptA, TPZVec<REAL> & ptB)
{
    if(ptA.NElements() != ptB.NElements())
    {
        DebugStop();
    }
    
    REAL dist = 0.;
    for(int i = 0; i < ptA.NElements(); i++)
    {
        dist += (ptA[i] - ptB[i]) * (ptA[i] - ptB[i]);
    }
    dist = sqrt(dist);
    
    return dist;
}


TPZGeoMesh * TRMSimworxMeshGenerator::MergeGeoMeshes(TPZGeoMesh * reservGMesh,
                            TPZAutoPointer<TPZGeoMesh> mioloGMesh,
                            const std::map<int,int> & miolo_reserv_nodeIndices)
{
    
    for(int el = 0; el < mioloGMesh->NElements(); el++)
    {
        TPZGeoEl * gel = mioloGMesh->ElementVec()[el];
        
        if(!gel || gel->HasSubElement()) continue;
        
        TPZManVector<int64_t> Topol(gel->NNodes());
        
        //pareamento dos nohs
        for(int n = 0; n < gel->NNodes(); n++)
        {
            std::map<int,int>::const_iterator it = miolo_reserv_nodeIndices.find(gel->NodeIndex(n));
            if(it != miolo_reserv_nodeIndices.end())
            {
                Topol[n] = it->second;
            }
            else
            {
                DebugStop();
            }
        }
        
        
        if(gel->Dimension() == 1)//Arc3D
        {
            new TPZGeoElRefPattern<pzgeom::TPZArc3D>(Topol, _BlendArcMatId, *reservGMesh);
            continue;
        }
        else
        {
            bool createGeoBlend = false;
            for(int s = gel->NNodes(); s < gel->NSides(); s++)
            {
                //Verificando se eh vizinho de arc3D
                TPZGeoElSide geoSide(gel,s);
                TPZGeoElSide neighSide(geoSide.Neighbour());
                while(neighSide != geoSide)
                {
                    if(neighSide.Element()->Dimension() == 1)//vizinho de arc3D.
                    {
                        createGeoBlend = true;
                        break;
                    }
                    neighSide = neighSide.Neighbour();
                }
            }
            
            int64_t index;
            if(createGeoBlend)
            {
                reservGMesh->CreateGeoBlendElement(gel->Type(),Topol,gel->MaterialId(),index);
            }
            else
            {
                reservGMesh->CreateGeoElement(gel->Type(),Topol,gel->MaterialId(),index);
            }
        }
    }
    
    reservGMesh->BuildConnectivity();
    
    //Criando faces externas que envolvem a malha.
    int nelements = reservGMesh->NElements();
    for(int el = 0; el < nelements; el++)
    {
        TPZGeoEl * gel = reservGMesh->ElementVec()[el];
        if(!gel) continue;
        if(gel->HasSubElement()) continue;
        
        for(int f = gel->NNodes(); f < gel->NSides()-1; f++)
        {
            TPZGeoElSide side(gel,f);
            if(side.Dimension() == 2)
            {
                if(side.Neighbour() == side)
                {
                    gel->CreateBCGeoEl(f, _mioloBoundary);
                }
            }
        }
    }
    
    int nnodes = reservGMesh->NNodes();
    for (int in=0; in<nnodes; in++) {
        reservGMesh->NodeVec()[in].SetNodeId(in);
    }
    reservGMesh->SetMaxNodeId(nnodes);
    
    return reservGMesh;
}

void TRMSimworxMeshGenerator::HangFacesTreatment(TPZGeoMesh * mergedGeoMesh)
{
    int nelem = mergedGeoMesh->NElements();
    
    std::list<TPZGeoElSide> reservFaces;
    std::list<TPZGeoElSide> mioloFaces;
    
    for(int el = 0; el < nelem; el++)
    {
        TPZGeoEl * gel = mergedGeoMesh->ElementVec()[el];
        if(!gel)
        {
            continue;
        }
        if(gel->MaterialId() == _mioloReservFaces)
        {
            if(gel->Dimension() != 2)
            {
                DebugStop();
            }
            TPZGeoElSide face(gel,gel->NSides()-1);
            reservFaces.push_back(face);
        }
        else if(gel->MaterialId() == _mioloBoundary)
        {
            if(gel->Dimension() != 2)
            {
                DebugStop();
            }
            TPZGeoElSide face(gel,gel->NSides()-1);
            mioloFaces.push_back(face);
        }
    }
    
    std::list<TPZGeoElSide>::iterator itReserv;
    for(itReserv = reservFaces.begin(); itReserv != reservFaces.end(); itReserv++)
    {
        TPZGeoElSide reservFace = *itReserv;
        
        TPZStack<TPZGeoElSide> allneigh;
        reservFace.AllNeighbours(allneigh);
        if (allneigh.NElements() == 2) {
            //Se face tem 2 vizinhos (miolo3D, reservatorio3D), nao precisa fazer nada.
            continue;
        }
        
        //Caso contrario, eh para ter apenas 1 vizinho (reservatorio3D). Neste caso eh hang face.
        if(allneigh.NElements() != 1)
        {
            DebugStop();
        }
        
        TPZManVector<TPZGeoEl*> sons(0);
        std::map<REAL, TPZGeoEl *> sonorder;
        std::list<TPZGeoElSide>::iterator itMiolo;
        for(itMiolo = mioloFaces.begin(); itMiolo != mioloFaces.end(); itMiolo++)
        {
            TPZGeoElSide mioloFace = *itMiolo;
            
            TPZManVector<REAL,2> xCenterMiolo(3,0.);
            mioloFace.CenterX(xCenterMiolo);
            
            if(PointIsInsideElementDomain(reservFace,xCenterMiolo))
            {
                sonorder[xCenterMiolo[1]] = mioloFace.Element();
            }
        }
        
        
        if (sonorder.size() == 0) {
            TPZManVector<REAL,3> xCenterFace(3,0.);
            reservFace.CenterX(xCenterFace);
            std::cout << "Center face " << xCenterFace << std::endl;
            DebugStop();
        }
        sons.Resize(sonorder.size(), 0);
        int index = 0;
        for (std::map<REAL,TPZGeoEl *>::iterator it = sonorder.begin(); it != sonorder.end(); it++,index++) {
            sons[index] = it->second;
        }
        RefineReservFaceToMatchMiolo(reservFace.Element(),sons);
    }
    mergedGeoMesh->BuildConnectivity();
}



bool TRMSimworxMeshGenerator::PointIsInsideElementDomain(TPZGeoElSide gelside, TPZVec<REAL> & coord)
{
    TPZManVector<REAL,3> XCenter(3);
    gelside.CenterX(XCenter);
    if (fabs(coord[0]-XCenter[0]) > 1.e-3 || fabs(coord[2]-XCenter[2]) > 1.e-3 ) {
        return false;
    }
    REAL Y0,Y2;
    REAL ymin, ymax;
    Y0 = gelside.Element()->NodePtr(0)->Coord(1);
    Y2 = gelside.Element()->NodePtr(2)->Coord(1);
    ymin = std::min(Y0, Y2);
    ymax = std::max(Y0, Y2);
    if (coord[1] < ymin || coord[1] > ymax) {
        return false;
    }
    return true;
}


void TRMSimworxMeshGenerator::RefineReservFaceToMatchMiolo(TPZGeoEl * reservGel, const TPZVec<TPZGeoEl*> & sons)
{
    TPZVec<TPZGeoEl *> elementVec(sons.NElements()+1);
    elementVec[0] = reservGel;
    for(int s = 1; s < elementVec.NElements(); s++)
    {
        elementVec[s] = sons[s-1];
    }
    TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::GetRefPatternBasedOnRealMeshElements(elementVec);
    if(!refp)
    {
        DebugStop();
    }
    reservGel->SetRefPattern(refp);
    
    for(int s = 1; s < elementVec.NElements(); s++)
    {
        elementVec[0]->SetSubElement(s-1,elementVec[s]);
        elementVec[s]->SetFather(elementVec[0]);
        elementVec[s]->SetMaterialId(elementVec[0]->MaterialId());
    }
}


void TRMSimworxMeshGenerator::MapCoordIn2DEllipse(const TPZVec<REAL> &coordIn,
                         const REAL innerRectangleLx,
                         const REAL innerRectangleLy,
                         TPZVec<REAL> &coordOut)
{
    if(!m_auxGMesh)
    {
        DebugStop();
    }
    if(m_auxGMesh->NElements() == 0)
    {
        DebugStop();
    }
    if(!m_auxGMesh->ElementVec()[0])
    {
        DebugStop();
    }
    
    TPZManVector<REAL,2> qsi(2);
    qsi[0] = coordIn[0]/(innerRectangleLx/2.);
    qsi[1] = coordIn[1]/(innerRectangleLy/2.);
    m_auxGMesh->ElementVec()[0]->X(qsi,coordOut);
}


bool TRMSimworxMeshGenerator::IsNodeOn2DCentralCore(int nnodes, int nodeIndex)
{
    if(nodeIndex < 5 || nodeIndex > (nnodes-6))
    {
        return false;
    }
    else
    {
        if(nodeIndex%4 == 1 || nodeIndex%4 == 2)
        {
            return true;
        }
    }
    
    return false;
}




void TRMSimworxMeshGenerator::InsertEllipseArcs(TPZGeoMesh * gmesh, REAL semiX, REAL semiY, REAL zCoord)
{
    int nnodes = gmesh->NNodes();
    
    TPZManVector<REAL,3> ellipOrigin(3,0.); ellipOrigin[2] = zCoord;
    TPZManVector<REAL,3> semiAxeX(3,0.), semiAxeY(3,0.);
    semiAxeX[0] = semiX; semiAxeY[1] = semiY;
    TPZManVector<int64_t,2> Topol(2);
    
    ///////////////////////////bottom
    Topol[0] = 0;
    Topol[1] = 1;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge0 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,arcEllip, *(gmesh));
    ellipEdge0->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    Topol[0] = 1;
    Topol[1] = 2;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge1 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,arcEllip, *(gmesh));
    ellipEdge1->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    Topol[0] = 2;
    Topol[1] = 3;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge2 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,arcEllip, *(gmesh));
    ellipEdge2->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    ///////////////////////////top
    Topol[0] = nnodes-1;
    Topol[1] = nnodes-2;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge3 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,arcEllip, *(gmesh));
    ellipEdge3->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    Topol[0] = nnodes-2;
    Topol[1] = nnodes-3;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge4 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,arcEllip, *(gmesh));
    ellipEdge4->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    Topol[0] = nnodes-3;
    Topol[1] = nnodes-4;
    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge5 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,arcEllip, *(gmesh));
    ellipEdge5->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
    
    ///////////////////////////laterals
    int node0right = 3;
    int node1right = 7;
    int node0left = nnodes-4;
    int node1left = nnodes-8;
    
    while(node1left >= 0)
    {
        Topol[0] = node0right;
        Topol[1] = node1right;
        TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdgeRight = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,arcEllip, *(gmesh));
        ellipEdgeRight->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
        
        Topol[0] = node0left;
        Topol[1] = node1left;
        TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdgeLeft = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (Topol,arcEllip, *(gmesh));
        ellipEdgeLeft->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
        
        node0right = node1right;
        node1right += 4;
        
        node0left = node1left;
        node1left -= 4;
    }
    
    gmesh->BuildConnectivity();
}



TPZGeoMesh * TRMSimworxMeshGenerator::ElementsCrusher2D(TPZGeoMesh * gmesh, int ndiv)
{
    TPZAutoPointer<TPZRefPattern> refTri = gRefDBase.FindRefPattern("Tri0001110Local");
    if(!refTri)
    {
        std::stringstream refTriSTR;
        refTriSTR << "6 4\n"
        << "-50 Tri0001110Local\n"
        << "0. 0. 0.\n"
        << "1. 0. 0.\n"
        << "0. 1. 0.\n"
        << "0.5 0. 0.\n"
        << "0.5 0.5 0.\n"
        << "0. 0.5 0.\n"
        << "2 3 0 1 2\n"
        << "3 4 0 3 4 5\n"
        << "2 3 3 1 4\n"
        << "2 3 5 4 2";
        refTri = new TPZRefPattern(refTriSTR);
        gRefDBase.InsertRefPattern(refTri);
    }
    
    TPZAutoPointer<TPZRefPattern> refQuad2edges = gRefDBase.FindRefPattern("Qua000010100Local");
    if(!refQuad2edges)
    {
        std::stringstream refQuad2edgesSTR;
        refQuad2edgesSTR << "6 3\n"
        << "-51 Qua000010100Local\n"
        << "-1. -1. 0.\n"
        << "1. -1. 0.\n"
        << "1. 1. 0.\n"
        << "-1. 1. 0.\n"
        << "0. -1. 0.\n"
        << "0. 1. 0.\n"
        << "3 4 0 1 2 3\n"
        << "3 4 0 4 5 3\n"
        << "3 4 4 1 2 5";
        refQuad2edges = new TPZRefPattern(refQuad2edgesSTR);
        gRefDBase.InsertRefPattern(refQuad2edges);
    }
    
    TPZAutoPointer<TPZRefPattern> refQuadVertical = gRefDBase.FindRefPattern("Qua000001010Local");
    if(!refQuadVertical)
    {
        std::stringstream refQuadVerticalSTR;
        refQuadVerticalSTR << "6 3\n"
        << "-52 Qua000001010Local\n"
        << "-1. -1. 0.\n"
        << "1. -1. 0.\n"
        << "1. 1. 0.\n"
        << "-1. 1. 0.\n"
        << "-1. 0. 0.\n"
        << "1. 0. 0.\n"
        << "3 4 0 1 2 3\n"
        << "3 4 0 1 5 4\n"
        << "3 4 4 5 2 3";
        refQuadVertical = new TPZRefPattern(refQuadVerticalSTR);
        gRefDBase.InsertRefPattern(refQuadVertical);
    }
    
    TPZAutoPointer<TPZRefPattern> refQuadFlag = gRefDBase.FindRefPattern("Qua000001001Local");
    if(!refQuadFlag)
    {
        std::stringstream refQuadFlagSTR;
        refQuadFlagSTR << "6 4\n"
        << "-52 Qua000001001Local\n"
        << "-1. -1. 0.\n"
        << "1. -1. 0.\n"
        << "1. 1. 0.\n"
        << "-1. 1. 0.\n"
        << "-1. 0. 0.\n"
        << "0. 0. 0.\n"
        << "3 4 0 1 2 3\n"
        << "3 4 0 1 5 4\n"
        << "3 4 4 5 2 3\n"
        << "2 3 5 1 2";
        refQuadFlag = new TPZRefPattern(refQuadFlagSTR);
        gRefDBase.InsertRefPattern(refQuadFlag);
    }
    
    //Refinando elementos, exceto o prolongamento em Y do miolo para farfield.
    for(int d = 0; d < ndiv; d++)
    {
        int nelem = gmesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[el];
            if(gel->HasSubElement() || gel->MaterialId() == arcEllip || gel->MaterialId() == mioloMat)
            {
                continue;
            }
            if(gel->Type() == ETriangle)
            {
                gel->SetRefPattern(refTri);
                TPZManVector<TPZGeoEl*> sons(0);
                gel->Divide(sons);
            }
            else if(gel->Type() == EQuadrilateral)
            {
                if(gel->MaterialId() == reservTriangleMat)
                {
                    gel->SetRefPattern(gRefDBase.GetUniformRefPattern(EQuadrilateral));
                }
                else
                {
                    gel->SetRefPattern(refQuad2edges);
                }
                TPZManVector<TPZGeoEl*> sons(0);
                gel->Divide(sons);
            }
        }
    }
    
    //Refinando prolongamento em Y do miolo para farfield.
    int nelem = gmesh->NElements();
    for(int el = 0; el < nelem; el++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[el];
        if(gel->HasSubElement() || gel->MaterialId() != reservMatTemp)
        {
            continue;
        }
        if(gel->Type() != EQuadrilateral)
        {
            DebugStop();
        }
        
        bool mioloNeigh = false;
        for(int s = 4; s <= 7; s++)
        {
            if(gel->Neighbour(s).Element()->MaterialId() == mioloMat)
            {
                mioloNeigh = true;
                break;
            }
        }
        if(mioloNeigh)
        {
            gel->SetRefPattern(refQuadFlag);
            TPZManVector<TPZGeoEl*> sons(0);
            gel->Divide(sons);
        }
        else
        {
            gel->SetRefPattern(refQuadVertical);
            TPZManVector<TPZGeoEl*> sons(0);
            gel->Divide(sons);
        }
    }
    
    RefineEllipseArcs(gmesh);
    gmesh = ExtractLeaf2DMesh(gmesh);
    
    return gmesh;
}

void TRMSimworxMeshGenerator::RefineEllipseArcs(TPZGeoMesh * gmesh)
{
    bool couldStop = false;
    
    while(couldStop == false)
    {
        couldStop = true;
        
        int nelem = gmesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[el];
            if(gel->HasSubElement() || gel->Dimension() != 1)
            {
                continue;
            }
            
            TPZGeoElSide edgeSide(gel,2);
            TPZGeoElSide neighSide(edgeSide.Neighbour());
            if(edgeSide == neighSide) DebugStop();//tem que ter vizinho
            
            //Procurando vizinho que foi refinado pelo lado igual a mim (gel)
            TPZAutoPointer< TPZRefPattern > sideRefPat = NULL;
            while( !sideRefPat && neighSide != edgeSide ){
                if(neighSide.Element()->GetRefPattern()){
                    sideRefPat = neighSide.Element()->GetRefPattern()->SideRefPattern(neighSide.Side());
                }
                if(!sideRefPat)
                {
                    neighSide = neighSide.Neighbour();
                }
            }
            if(!sideRefPat) continue;
            if(neighSide.Element()->HasSubElement() == false) continue;
            
            gel->SetRefPattern( sideRefPat );
            TPZManVector<TPZGeoEl*> sons(0);
            gel->Divide(sons);
            
            couldStop = false;
        }
    }
}

TPZGeoMesh * TRMSimworxMeshGenerator::Extrude2DMesh(TPZGeoMesh * gmesh2D, const TPZVec<REAL> & espacamentoZ, const bool thereIsCutPlane)
{
    int nnodesPerPlane = gmesh2D->NNodes();
    int nPlanes = espacamentoZ.NElements();
    
    TPZGeoMesh * gmesh3D = new TPZGeoMesh;
    gmesh3D->NodeVec().Resize(nnodesPerPlane * nPlanes);
    
    for(int plane = 0; plane < nPlanes; plane++)
    {
        REAL z = espacamentoZ[plane];
        
        for(int n = 0; n < nnodesPerPlane; n++)
        {
            TPZManVector<REAL,3> coord(3);
            gmesh2D->NodeVec()[n].GetCoordinates(coord);
            coord[2] = z;
            gmesh3D->NodeVec()[plane*nnodesPerPlane + n].SetCoord(coord);
        }
    }
    
    for(int el = 0; el < gmesh2D->NElements(); el++)
    {
        TPZGeoEl * gel = gmesh2D->ElementVec()[el];
        for(int plane = 0; plane < nPlanes-1; plane++)
        {
            bool mioloPlane = (espacamentoZ[plane] < 0. && espacamentoZ[plane+1] > 0.);
            if(thereIsCutPlane && mioloPlane == false && plane < (nPlanes-2))
            {
                mioloPlane = (espacamentoZ[plane+1] < 0. && espacamentoZ[plane+2] > 0.);
            }
            if(mioloPlane && gel->MaterialId() == mioloMat)
            {
                //Elementos no miolo serao gerados por um gerador de malha exclusivo
                //e incluidos nesta malha posteriormente.
                continue;
            }
            if(gel->NCornerNodes() > 2)
            {
                //Criacao de hexaedro ou elipse.
                TPZManVector<int64_t,8> Topol3D(2*gel->NNodes());
                for(int n = 0; n < gel->NNodes(); n++)
                {
                    Topol3D[n] = plane*nnodesPerPlane + gel->NodeIndex(n);
                    Topol3D[n+gel->NNodes()] = (plane+1)*nnodesPerPlane + gel->NodeIndex(n);
                }
                if(gel->IsGeoBlendEl())
                {
                    //Elementos que encostam no farfield sao GeoBlend.
                    if(gel->NCornerNodes() == 3)
                    {
                        new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoPrism> >(Topol3D, _ReservMatId, *(gmesh3D));
                    }
                    else if(gel->NCornerNodes() == 4)
                    {
                        new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(Topol3D, _ReservMatId, *(gmesh3D));
                    }
                }
                else
                {
                    //Elementos no interior do dominio sao Lineares.
                    if(gel->NCornerNodes() == 3)
                    {
                        new TPZGeoElRefPattern<pzgeom::TPZGeoPrism>(Topol3D, _ReservMatId, *(gmesh3D));
                    }
                    else if(gel->NCornerNodes() == 4)
                    {
                        new TPZGeoElRefPattern<pzgeom::TPZGeoCube>(Topol3D, _ReservMatId, *(gmesh3D));
                    }
                }
            }
            else
            {
                //Criacao de arco de elipse.
                TPZManVector<int64_t,2> Topol3D(2);
                Topol3D[0] = plane*nnodesPerPlane + gel->NodeIndex(0);
                Topol3D[1] = plane*nnodesPerPlane + gel->NodeIndex(1);
                TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * newEllip = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D>(Topol3D, _ellipseArcMat, *(gmesh3D));
                
                TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellip = NULL;
                if(gel->Father())
                {
                    ellip = dynamic_cast< TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * > (gel->EldestAncestor());
                }
                else
                {
                    ellip = dynamic_cast< TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * > (gel);
                }
                if(!ellip)
                {
                    DebugStop();
                }
                TPZManVector<REAL,3> origin = ellip->Geom().Origin();
                origin[2] = espacamentoZ[plane];
                TPZManVector<REAL,3> semiX = ellip->Geom().SemiAxeX();
                TPZManVector<REAL,3> semiY = ellip->Geom().SemiAxeY();
                newEllip->Geom().SetAxes(origin,semiX,semiY);
                
                if(plane == nPlanes-2)
                {
                    //A camada mais acima de elipses eh contemplada aqui.
                    Topol3D[0] = (nPlanes-1)*nnodesPerPlane + gel->NodeIndex(0);
                    Topol3D[1] = (nPlanes-1)*nnodesPerPlane + gel->NodeIndex(1);
                    newEllip = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D>(Topol3D, _ellipseArcMat, *(gmesh3D));
                    
                    origin[2] = espacamentoZ[plane+1];
                    newEllip->Geom().SetAxes(origin,semiX,semiY);
                }
            }
        }
    }
    
    gmesh3D->BuildConnectivity();
    
    return gmesh3D;
}

TPZGeoMesh * TRMSimworxMeshGenerator::ExtractLeaf2DMesh(TPZGeoMesh * gmesh)
{
    TPZGeoMesh * leafGeoMesh = new TPZGeoMesh;
    
    leafGeoMesh->NodeVec() = gmesh->NodeVec();
    
    leafGeoMesh->ElementVec().Resize(0);
    
    for(int el = 0; el < gmesh->NElements(); el++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[el];
        
        if(!gel || gel->HasSubElement()) continue;
        
        TPZManVector<int64_t> nodeindexes(gel->NNodes());
        for(int n = 0; n < gel->NNodes(); n++)
        {
            nodeindexes[n] = gel->NodeIndex(n);
        }
        CheckNodesSequenceForExtrusion(gmesh,nodeindexes);
        
        if(gel->Dimension() == 1)
        {
            //Elipse
            TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * newEllip = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D>(nodeindexes, gel->MaterialId(), *(leafGeoMesh));
            
            TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellip = NULL;
            if(gel->Father())
            {
                ellip = dynamic_cast< TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * > (gel->EldestAncestor());
            }
            else
            {
                ellip = dynamic_cast< TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * > (gel);
            }
            if(!ellip)
            {
                DebugStop();
            }
            TPZManVector<REAL,3> origin = ellip->Geom().Origin();
            TPZManVector<REAL,3> semiX = ellip->Geom().SemiAxeX();
            TPZManVector<REAL,3> semiY = ellip->Geom().SemiAxeY();
            newEllip->Geom().SetAxes(origin,semiX,semiY);
        }
        else
        {
            //Triangulo ou quadrilatero
            bool createGeoBlend = false;
            
            for(int s = gel->NNodes(); s < gel->NSides(); s++)
            {
                //Verificando se eh vizinho de elipse
                TPZGeoElSide geoSide(gel,s);
                TPZGeoElSide neighSide(geoSide.Neighbour());
                while(neighSide != geoSide)
                {
                    if(neighSide.Element()->Dimension() == 1)
                    {
                        createGeoBlend = true;
                        break;
                    }
                    neighSide = neighSide.Neighbour();
                }
            }
            
            if(gel->Type() == ETriangle)
            {
                if(createGeoBlend)
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(nodeindexes, gel->MaterialId(), *(leafGeoMesh));
                }
                else
                {
                    new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeindexes, gel->MaterialId(), *(leafGeoMesh));
                }
            }
            else if(gel->Type() == EQuadrilateral)
            {
                if(createGeoBlend)
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(nodeindexes, gel->MaterialId(), *(leafGeoMesh));
                }
                else
                {
                    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeindexes, gel->MaterialId(), *(leafGeoMesh));
                }
            }
        }
    }
    
    leafGeoMesh->BuildConnectivity();
    
    return leafGeoMesh;
}

void TRMSimworxMeshGenerator::CheckNodesSequenceForExtrusion(TPZGeoMesh * gmesh, TPZVec<int64_t> &nodeIndexes)
{
    if(nodeIndexes.NElements() != 3 && nodeIndexes.NElements() != 4)
    {
        return;
    }
    
    TPZManVector<REAL,3> coord0(3), coord1(3), coord2(3);
    gmesh->NodeVec()[nodeIndexes[0]].GetCoordinates(coord0);
    gmesh->NodeVec()[nodeIndexes[1]].GetCoordinates(coord1);
    gmesh->NodeVec()[nodeIndexes[2]].GetCoordinates(coord2);
    
    REAL zComponentVectorialProduct = coord0[0]*coord1[1]   //componente z do produto vetorial (coord1-coord0)x(coord2-coord0)
    + coord1[0]*coord2[1]
    + coord0[1]*coord2[0]
    - coord0[1]*coord1[0]
    - coord1[1]*coord2[0]
    - coord0[0]*coord2[1];
    
    if(zComponentVectorialProduct < 0.)
    {
        //Corrigir sequencia topologica.
        int64_t indexTemp = nodeIndexes[1];
        nodeIndexes[1] = nodeIndexes[nodeIndexes.NElements()-1];  //Deixe assim pois funciona para triangulo e quadrilatero simultaneamente!
        nodeIndexes[nodeIndexes.NElements()-1] = indexTemp;
    }
}

void TRMSimworxMeshGenerator::FillStructMiolo(const REAL mioloLx,
                     const TPZVec<REAL> & espacamentoMioloY,
                     const TPZVec<REAL> & espacamentoZ,
                     const REAL Lw,
                     const bool thereIsCutPlane,
                     StructMioloData & mioloData,
                     TRMRawData &rawdata)
{
    mioloData.m_Lx = mioloLx;
    mioloData.m_Ly = espacamentoMioloY[espacamentoMioloY.NElements()-1] - espacamentoMioloY[0];
    mioloData.m_thereIsCutPlane = thereIsCutPlane;
    mioloData.m_CutPlaneZ = 0.;//serah definido mais ao final deste metodo.
    
    std::set<int> wellPosInMilimeters;//wellPosInMilimeters intervalo [0;+Lw] em milimetros
    wellPosInMilimeters.insert(-Lw/2.*1000.);
    wellPosInMilimeters.insert(+Lw/2.*1000.);
    
    for(int p = 0; p < espacamentoMioloY.NElements(); p++)
    {
        REAL pos = espacamentoMioloY[p];
        if(fabs(pos) < Lw/2.-1.e-2)
        {
            wellPosInMilimeters.insert(espacamentoMioloY[p]*1000.);
        }
    }
    bool LzTopBottomFound = false;
    for(int p = 0; p < espacamentoZ.NElements()-1; p++)
    {
        REAL pos = espacamentoZ[p];
        if(pos < 0. && espacamentoZ[p+1] > 0.)
        {
            mioloData.m_LzBottom = fabs(pos);
            mioloData.m_LzTop = espacamentoZ[p+1];
            
            if(thereIsCutPlane)
            {
                if(p == 0)
                {
                    //Quer plano de corte mas nao criou camada para aquifero logo abaixo do poco?
                    DebugStop();
                }
                mioloData.m_LzBottom = fabs(espacamentoZ[p-1]);
                mioloData.m_CutPlaneZ = espacamentoZ[p];
            }
            LzTopBottomFound = true;
            
            break;
        }
    }
    if(LzTopBottomFound == false)
    {
        DebugStop();
    }
    
    bool hasLiner = rawdata.fHasLiner;
    if(hasLiner)
    {
        //Segundo documentacao TPBRDataKernel.h --> class TLiner.
        /** Dados de furacao do liner. Matriz com tres colunas
         * (0) Posicao do final do trecho (em m), dado em WellPath (-Lw/2 a + Lw/2)
         * (1) Densidade de furos no liner: numero de furos / comprimento do tubo para cada trecho do liner
         * (2) Diametro aparente dos furos em m
         */
        TPZFMatrix<REAL> linerHolesData;// = DataKernel->fCompletion.fLiner.fFuracao;
        FillLinerHolesData(linerHolesData);
        int nrows = linerHolesData.Rows();
        for(int r = 0; r < nrows-1; r++)
        {
            double pos_meterBegin = linerHolesData(r,0); //pos_m intervalo [-Lw/2;+Lw/2]
            double pos_meterEnd = linerHolesData(r+1,0); //pos_m intervalo [-Lw/2;+Lw/2]
            double furos_metro = linerHolesData(r,1);
            if(furos_metro < 1.e-2)
            {
                int pos_milimeterBegin = pos_meterBegin*1000.;
                wellPosInMilimeters.insert(pos_milimeterBegin);
                
                int pos_milimeterEnd = pos_meterEnd*1000.;
                wellPosInMilimeters.insert(pos_milimeterEnd);
            }
        }
    }
    
    bool hasCasing = rawdata.fHasCasing;
    if(hasCasing)
    {
        //Segundo documentacao TPBRDataKernel.h --> class TCaseData.
        /** Dados de furacao do liner. Matriz com tres colunas
         * (0) Posicao do final do trecho (em m), dado em WellPath (-Lw/2 a + Lw/2)
         * (1) Densidade de furos no liner: numero de furos / comprimento do tubo para cada trecho do liner
         * (2) Diametro aparente dos furos em m
         */
        TPZFMatrix<REAL> caseHolesData;// = DataKernel->fCompletion.fCaseData.fCanhoneado;
        FillCaseHolesData(caseHolesData);
        int nrows = caseHolesData.Rows();
        for(int r = 0; r < nrows-1; r++)
        {
            double pos_meterBegin = caseHolesData(r,0); //pos_m intervalo [-Lw/2;+Lw/2]
            double pos_meterEnd = caseHolesData(r+1,0); //pos_m intervalo [-Lw/2;+Lw/2]
            double furos_metro = caseHolesData(r,1);
            if(furos_metro < 1.e-2)
            {
                int pos_milimeterBegin = pos_meterBegin*1000.;
                wellPosInMilimeters.insert(pos_milimeterBegin);
                
                int pos_milimeterEnd = pos_meterEnd*1000.;
                wellPosInMilimeters.insert(pos_milimeterEnd);
            }
        }
    }
    
    mioloData.m_wellPositionsY.Resize(wellPosInMilimeters.size());
    mioloData.m_ribsIndexes.Resize(0);
    
    std::set<int>::iterator it;
    int p = 0;
    for(it = wellPosInMilimeters.begin(); it != wellPosInMilimeters.end(); it++, p++)
    {
        REAL pos_meter = (*it)/1000.;
        for(int i = 0; i < espacamentoMioloY.NElements(); i++)
        {
            REAL pos = espacamentoMioloY[i];
            if(fabs(pos - pos_meter) < 1.e-1)
            {
                int oldSize = mioloData.m_ribsIndexes.NElements();
                mioloData.m_ribsIndexes.Resize(oldSize+1);
                mioloData.m_ribsIndexes[oldSize] = p;
            }
        }
        mioloData.m_wellPositionsY[p] = pos_meter;
    }
}

void TRMSimworxMeshGenerator::CreateReservoirFaces(TPZGeoMesh & gmesh, StructMioloData & mioloData)
{
    const int nel = gmesh.NElements();
    mioloData.m_reservInterfaceNodeIndices.clear();
    
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl * gel = gmesh.ElementVec()[iel];
        if (!gel || gel->HasSubElement() || gel->Dimension() != 3)
        {
            continue;
        }
        
        const int nsides = gel->NSides();
        for (int is = gel->NNodes(); is < nsides; is++)
        {
            TPZGeoElSide myself(gel, is);
            if (myself.Dimension() != 2) continue;
            
            TPZGeoElSide neighbour = myself.Neighbour();
            if (neighbour == myself)
            {
                TPZGeoEl * face = NULL;
                
                if(gel->IsGeoBlendEl())
                {
                    face = CreateBCGeoBlendEl(gel, is, _LateralReservBC);
                }
                else
                {
                    face = gel->CreateBCGeoEl(is, _LateralReservBC);
                }
                
                ///verificar se eh face do miolo
                TPZManVector<REAL,2> centerQsi(2);
                TPZManVector<REAL,3> centerX(3);
                
                face->CenterPoint(face->NSides()-1, centerQsi);
                face->X(centerQsi,centerX);
                
                REAL tol = 1.e-2;
                if( ( fabs(centerX[0]) < (mioloData.m_Lx/2. + tol) ) &&
                   ( fabs(centerX[1]) < (mioloData.m_Ly/2. + tol) ) &&
                   ( centerX[2] < (mioloData.m_LzTop + tol) ) &&
                   ( centerX[2] > (-mioloData.m_LzBottom - tol) ) )
                {
                    face->SetMaterialId(_mioloReservFaces);
                    
                    for(int n = 0; n < face->NCornerNodes(); n++)
                    {
                        mioloData.m_reservInterfaceNodeIndices.insert(face->NodeIndex(n));
                    }
                }
                else
                {
                    ///verificar se a face nao eh de camada confinante
                    TPZFNMatrix<9> jac(2,2);
                    TPZFNMatrix<9> axes(2,3);
                    TPZFNMatrix<9> jacinv(2,2);
                    REAL detjac;
                    
                    face->Jacobian(centerQsi, jac, axes, detjac, jacinv);
                    
                    if (fabs(axes(0,2)) + fabs(axes(1,2)) < 1.e-10)
                    {
                        ///entao a normal eh vertical
                        if(face->Mesh()->NodeVec()[face->NodeIndex(0)].Coord(2) > 0.)
                        {
                            face->SetMaterialId(_ConfinementReservBCtop);
                        }
                        else
                        {
                            face->SetMaterialId(_ConfinementReservBCbottom);
                        }
                    }
                }
            }
        }///for is
    }///for iel
}


TPZGeoEl * TRMSimworxMeshGenerator::CreateBCGeoBlendEl(TPZGeoEl *orig, int side, int bc)
{
    int ns = orig->NSideNodes(side);
    TPZManVector<int64_t> nodeindices(ns);
    int in;
    for(in=0; in<ns; in++)
    {
        nodeindices[in] = orig->SideNodeIndex(side,in);
    }
    int64_t index;
    
    TPZGeoMesh * mesh = orig->Mesh();
    MElementType type = orig->Type(side);
    
    TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
    TPZGeoElSide me(orig,side);
    TPZGeoElSide newelside(newel,newel->NSides()-1);
    
    newelside.InsertConnectivity(me);
    newel->Initialize();
    
    return newel;
}

/// adiciona elementos unidimensionais no poco tri dimensional
void TRMSimworxMeshGenerator::AddRibElements(TPZGeoMesh *gmesh, int WellMatId1D, int WellMatFake1D)
{
#ifdef LOGANDO
    std::ofstream check("../testribs.txt");
#endif
    int nel = gmesh->NElements();
    TPZStack<TPZGeoEl *> createdwell1d;
    // loop over the 3d well elements
    for(int el=0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        if(!gel) continue;
        if(gel->MaterialId() != _WellMatId3D)
        {
            continue;
        }
        
        // data structure to keep track of the 4 ribs and associated faces
        TWellRibs newrib(el);
        // data structure to keep track of faces that have been considered
        std::set<int> facesides;
        int nsides = gel->NSides();
        int ribcount = 0;
        for(int is=0; is<nsides; is++)
        {
            int dim = gel->SideDimension(is);
            if(dim != 1)
            {
                continue;
            }
            TPZManVector<REAL,3> x1(3),x2(3),diff(3);
            int nod1 = gel->SideNodeIndex(is,0);
            int nod2 = gel->SideNodeIndex(is,1);
            gmesh->NodeVec()[nod1].GetCoordinates(x1);
            gmesh->NodeVec()[nod2].GetCoordinates(x2);
            // look for the ribs along the y axis
            // create elements along the ribs
            // along the main rib the material id = WellMatId1d
            // along the three other ribs it is WellMatFake1D
            // Computational elements of this type will have their connects restrained
            for(int i=0; i<3; i++) diff[i] = x2[i]-x1[i];
            if(std::abs(diff[0]) < 1.e-3 && std::abs(diff[2]) < 1.e-3)
            {
#ifdef LOGANDO
                check << "Element " << gel->Index() << " co1 " << x1 << " co2 " << x2 << std::endl;
#endif
                if(ribcount == 0)
                {
                    
                    TPZGeoElBC gbc(gel,is,WellMatId1D);
#ifdef LOGANDO
                    gbc.CreatedElement()->Print(check);
#endif
                    createdwell1d.Push(gbc.CreatedElement());
                    newrib.fRibElements[ribcount] = gbc.CreatedElement()->Index();
                    TPZStack<TPZGeoElSide> faces;
                    gel->AllHigherDimensionSides(is,2,faces);
                    if(faces.NElements() != 2) DebugStop();
                    // find a face that hasnt been found yet
                    int correctface = -1;
                    for(int i=0; i<2; i++)
                    {
                        if(facesides.find(faces[i].Side()) == facesides.end())
                        {
                            correctface = faces[i].Side();
                            break;
                        }
                    }
                    if(correctface == -1) DebugStop();
                    facesides.insert(correctface);
                    TPZGeoElSide gelside(gel,correctface);
                    TPZGeoElSide neighbour = gelside.Neighbour();
                    if(!neighbour.Element()) DebugStop();
                    newrib.fReservoirSides[ribcount] = neighbour;
                    ribcount++;
                }
                else
                {
                    TPZGeoElBC gbc(gel,is,WellMatFake1D);
#ifdef LOGANDO
                    gbc.CreatedElement()->Print(check);
#endif
                    newrib.fRibElements[ribcount] = gbc.CreatedElement()->Index();
                    TPZStack<TPZGeoElSide> faces;
                    gel->AllHigherDimensionSides(is,2,faces);
                    if(faces.NElements() != 2) DebugStop();
                    // find a face that hasnt been found yet
                    int correctface = -1;
                    for(int i=0; i<2; i++)
                    {
                        if(facesides.find(faces[i].Side()) == facesides.end())
                        {
                            correctface = faces[i].Side();
                            break;
                        }
                    }
                    if(correctface == -1) DebugStop();
                    facesides.insert(correctface);
                    TPZGeoElSide gelside(gel,correctface);
                    TPZGeoElSide neighbour = gelside.Neighbour();
                    if(!neighbour.Element()) DebugStop();
                    newrib.fReservoirSides[ribcount] = neighbour;
                    ribcount++;
                }
            }
        }
        fRibs.Push(newrib);
    }
    int nwell1d = createdwell1d.NElements();
    int firstwell = fRibs[0].fRibElements[0];
    TPZGeoEl *firstgel = gmesh->ElementVec()[firstwell];
    TPZGeoElSide heelside(firstgel, 0);
    int lastwell = fRibs[nwell1d-1].fRibElements[0];
    TPZGeoEl *lastgel = gmesh->ElementVec()[lastwell];
    TPZGeoElSide toeside(lastgel, 1);
    TPZGeoElBC heelbc(heelside,_WellHeelMatId);
    TPZGeoElBC toebc(toeside,_WellToeMatId);
    
#ifdef LOGANDO
    int nr = fRibs.NElements();
    check << "Output for " << nr << " detected well elements\n";
    for(int ir=0; ir<nr; ir++)
    {
        fRibs[ir].Print(check);
    }
#endif
}

void TRMSimworxMeshGenerator::CreateWellBoundaries(TPZAutoPointer<TPZGeoMesh> reservoirGMesh){
    
    const int64_t nel =reservoirGMesh->NElements();
    
    for (int64_t iel = 0 ; iel < nel ; iel++) {
        
        // filters
        
        TPZGeoEl * gel = reservoirGMesh->Element(iel);
        if (!gel) {
            continue;
        }
        
        if (gel->MaterialId() != _WellMatId3D ) {
            continue;
        }
        
#ifdef PZDEBUG
        if (gel->Type() != ECube) {
            DebugStop();
        }
#endif
        
        const int nc = gel->NCornerNodes();
        const int nedges = 12;
        const int nfaces = 6;
        
        
        for (int ic = nc + nedges; ic < nc + nedges + nfaces; ic++) {
            bool hasinterface = false;
            bool hasendpoint = false;
            bool isinternal = false;
            bool isinterfaceside = false;
            TPZGeoElSide gelside(gel,ic);
            TPZGeoElSide neigh = gelside.Neighbour();
            while (gelside != neigh) {
                int neighmatid = neigh.Element()->MaterialId();
                if (neigh.Element()->Dimension() == 2 && (neighmatid == _WellToeMatId || neighmatid == _WellHeelMatId)) {
                    hasendpoint = true;
                }
                if (neigh.Element()->Dimension() == 3 && neigh.Element()->MaterialId() == _ReservMatId) {
                    isinterfaceside = true;
                }
                if (neigh.Element()->Dimension() == 3 && neigh.Element()->MaterialId() == _WellMatId3D) {
                    isinternal = true;
                }
                if (neighmatid == _Well3DReservoirFaces) {
                    hasinterface = true;
                }
                neigh = neigh.Neighbour();
            }
            if (hasinterface == true)
            {
                std::cout << __PRETTY_FUNCTION__ << " called twice?\n";
            }
            if (hasinterface == true && (ic < 21 || ic > 25)) {
                std::cout << __LINE__ << " I don't understand\n";
            }
            if (hasinterface == false && (ic >20 && ic < 25)) {
                TPZGeoElBC(gelside, _Well3DReservoirFaces);
            }
            if (hasendpoint == true && ic != 20 && ic != 25) {
                std::cout << __LINE__ << " I don't understand\n";
            }
            if (ic == 20 && ! hasendpoint && !isinternal) {
                TPZGeoElBC(gelside,_WellHeelMatId);
            }
            if (ic == 25 && !hasendpoint && !isinternal) {
                TPZGeoElBC(gelside,_WellToeMatId);
            }
        }
    }    
}

void TRMSimworxMeshGenerator::AdjustPrismCoordinates(TPZGeoMesh * gmesh, REAL semiX, REAL semiY)
{
    const int nnodes = gmesh->NNodes();
    
    const REAL xnn2 = gmesh->NodeVec()[nnodes-2].Coord(0);
    const REAL ynn2 = gmesh->NodeVec()[nnodes-2].Coord(1);
    
    const REAL xnn5 = gmesh->NodeVec()[nnodes-5].Coord(0);
    const REAL ynn5 = gmesh->NodeVec()[nnodes-5].Coord(1);
    
    const REAL tempX = (semiY)*(semiY)*(xnn2 + xnn5)*(xnn2 + xnn5) + (semiX)*(semiX)*(ynn2 + ynn5)*(ynn2 + ynn5);
    const REAL newX = (semiX*semiY*(xnn2 + xnn5))/sqrt(std::max(1.e-10,tempX));
    
    const REAL tempY = 1. - newX*newX/(semiX*semiX);
    const REAL newY = semiY*sqrt(std::max(0.,tempY));
    
    gmesh->NodeVec()[0].SetCoord(0,-newX);
    gmesh->NodeVec()[0].SetCoord(1,-newY);
    
    gmesh->NodeVec()[3].SetCoord(0,+newX);
    gmesh->NodeVec()[3].SetCoord(1,-newY);
    
    gmesh->NodeVec()[nnodes-1].SetCoord(0,+newX);
    gmesh->NodeVec()[nnodes-1].SetCoord(1,+newY);
    
    gmesh->NodeVec()[nnodes-4].SetCoord(0,-newX);
    gmesh->NodeVec()[nnodes-4].SetCoord(1,+newY);
}