//
//  TRMSimworxMeshGenerator.h
//  PZ
//
//  Created by Nathan Shauer on 6/15/15.
//
//

#ifndef __PZ__TRMSimworxMeshGenerator__
#define __PZ__TRMSimworxMeshGenerator__

#include <stdio.h>
#include "pzgmesh.h"
#include "pzfmatrix.h"

class TRMRawData;
class TPZGeoElSide; //a

/**
 * @brief Struct to store data from the miolo of the mesh
 * @author Geo-Adapta
 */
struct StructMioloData
{
public:
    
    StructMioloData();
    ~StructMioloData(){}
    
    /// largura do miolo (dimensao em planta na direcao X
    REAL m_Lx;
    /// comprimento do miolo (dimensao em planta na direcao Y
    REAL m_Ly;
    /// distancia em Z da base do miolo ao poco
    REAL m_LzBottom;
    /// distancia em Z do topo do miolo ao poco
    REAL m_LzTop;
    /// Flag que diz se tem plano de corte horizontal no miolo (motivacao: aquifero)
    bool m_thereIsCutPlane;
    /// Cota Z do plano de corte (quando houver)
    REAL m_CutPlaneZ;
    
    /// espacamento do poco
    TPZManVector<REAL,20> m_wellPositionsY;
    /// indices do vetor m_wellPositionsY que compreendem ribs
    TPZManVector<int,20> m_ribsIndexes;
    
    std::set<int> m_reservInterfaceNodeIndices;
};

/**
 * @brief Class that manages the generation of the simworx mesh
 * @author Geo-Adapta (Philippe Devloo more precisely)
 */
class TRMSimworxMeshGenerator
{
private:

    /// Auxilary geomesh
    TPZGeoMesh *m_auxGMesh = NULL;
    
public:
    
    TRMSimworxMeshGenerator();
    ~TRMSimworxMeshGenerator();
    
    TPZGeoEl * CreateBCGeoBlendEl(TPZGeoEl *orig, int side, int bc);
    
    void CreateReservoirFaces(TPZGeoMesh & gmesh, StructMioloData & mioloData);
    
    void CheckNodesSequenceForExtrusion(TPZGeoMesh * gmesh, TPZVec<long> &nodeIndexes);
    
    TPZGeoMesh * ExtractLeaf2DMesh(TPZGeoMesh * gmesh);
    
    TPZGeoMesh * Extrude2DMesh(TPZGeoMesh * gmesh2D, const TPZVec<REAL> & espacamentoZ, const bool thereIsCutPlane);
    
    TPZGeoMesh * ElementsCrusher2D(TPZGeoMesh * gmesh, int ndiv);
    
    void RefineEllipseArcs(TPZGeoMesh * gmesh);
    
    void MapCoordIn2DEllipse(const TPZVec<REAL> &coordIn,
                             const REAL innerRectangleLx,
                             const REAL innerRectangleLy,
                             TPZVec<REAL> &coordOut);
    
    void FillStructMiolo(const REAL mioloLx,
                         const TPZVec<REAL> & espacamentoMioloY,
                         const TPZVec<REAL> & espacamentoZ,
                         const REAL Lw,
                         const bool thereIsCutPlane,
                         StructMioloData & mioloData,
                         TRMRawData &rawdata);
    
    void InsertEllipseArcs(TPZGeoMesh * gmesh, REAL semiX, REAL semiY, REAL zCoord);
    
    bool IsNodeOn2DCentralCore(int nnodes, int nodeIndex);
    
    bool PointIsInsideElementDomain(TPZGeoElSide gelside, TPZVec<REAL> & coord);
    
    void HangFacesTreatment(TPZGeoMesh * mergedGeoMesh);
    
    void RefineReservFaceToMatchMiolo(TPZGeoEl * reservGel, const TPZVec<TPZGeoEl*> & sons);
    
    void BuildAuxiliary2DGeoMesh(const REAL semiX, const REAL semiY, REAL & innerRectangleLx, REAL & innerRectangleLy);
    
    TPZGeoMesh * MergeGeoMeshes(TPZGeoMesh * reservGMesh,
                                TPZAutoPointer<TPZGeoMesh> mioloGMesh,
                                const std::map<int,int> & miolo_reserv_nodeIndices);
    
    REAL DistFromPoints(TPZVec<REAL> & ptA, TPZVec<REAL> & ptB);
    
    void FillLinerHolesData(TPZFMatrix<> &linerHolesData);
    
    void FillCaseHolesData(TPZFMatrix<> &caseHolesData);
    
    TPZGeoMesh * CreateEllipticalReservoirGeoMesh(const REAL semiAxeX,
                                                  const REAL semiAxeY,
                                                  const REAL mioloLx,
                                                  const TPZVec<REAL> & espacamentoMioloY,
                                                  const TPZVec<REAL> & espacamentoZ,
                                                  const REAL Lw,
                                                  const bool thereIsCutPlane,
                                                  StructMioloData & mioloData,
                                                  TRMRawData &rawdata);
    
    TPZAutoPointer<TPZGeoMesh> ReallyGenerateGeoMesh(const REAL semiAxeX,
                                                     const REAL semiAxeY,
                                                     const REAL mioloLx,
                                                     const TPZVec<REAL> & espacamentoMioloY,
                                                     const TPZVec<REAL> & espacamentoZ,
                                                     const bool thereIsCutPlane,
                                                     TRMRawData &rawdata);
    TPZAutoPointer<TPZGeoMesh> CreateSimworxGeoMesh(TRMRawData &rawdata);
    
    void Pair_Miolo_Reserv_Nodes(TPZGeoMesh * reservGMesh, const std::set<int> & reservNodeIndices,
                                 TPZAutoPointer<TPZGeoMesh> mioloGMesh,
                                 std::map<int,int> & miolo_reserv_nodeIndices);
    
};


#endif /* defined(__PZ__TRMSimworxMeshGenerator__) */
