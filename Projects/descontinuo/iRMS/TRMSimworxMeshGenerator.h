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
#include "pzgeoelside.h"
#include "TPZInterfaceEl.h"

class TRMRawData;
class TPZGeoElSide;

/**
 * @brief Struct to store data from the miolo of the mesh
 * @author Geo-Adapta
 */
struct StructMioloData
{
public:
    
    StructMioloData();
    StructMioloData(const StructMioloData &copy)
    {
        DebugStop();
    }
    StructMioloData &operator=(const StructMioloData &copy)
    {
        DebugStop();
        return *this;
    }
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
 * @brief Stores and groups the elements along the sides of the well elements
 * @author Geo-Adapta
 */
struct TWellRibs
{
    /// Elemento do poco 3D
    int fWellEl3D;
    /// Elementos nas arrestas
    TPZManVector<int,4> fRibElements;
    
    /// Elemento unidimensional de liner
    int fLinerRibElement;
    
    /// Lados dos elementos de reservatorio
    TPZManVector<TPZGeoElSideIndex,4> fReservoirSides;
    
    /// Elementos de interface associados
    TPZManVector<TPZInterfaceElement *,4> fInterfaces;
    
    TWellRibs(int elindex = -1) : fWellEl3D(elindex), fRibElements(4,-1), fLinerRibElement(-1), fReservoirSides(4), fInterfaces(4,0)
    {
    }
    TWellRibs(const TWellRibs &cp) : fWellEl3D(cp.fWellEl3D), fRibElements(cp.fRibElements), fLinerRibElement(cp.fLinerRibElement),
    fReservoirSides(cp.fReservoirSides), fInterfaces(cp.fInterfaces)
    {
    }
    TWellRibs &operator=(const TWellRibs &cp)
    {
        fWellEl3D = cp.fWellEl3D;
        fRibElements = cp.fRibElements;
        fLinerRibElement = cp.fLinerRibElement;
        fReservoirSides = cp.fReservoirSides;
        fInterfaces = cp.fInterfaces;
        return *this;
    }
    
    void Print(std::ostream &out) const
    {
        out << "3D Element index " << fWellEl3D << std::endl;
        out << "Ribs Element Indices\n";
        for(int i=0; i<4;i++) out << fRibElements[i] << " ";
        out << "Liner Rib Element " << fLinerRibElement << std::endl;
        out << "\nReservoir Element/sides\n";
        //for(int i=0; i<4; i++) fReservoirSides[i].Print(out);
        out << "Interface element indexes ";
        for(int i=0; i<4; i++) if(fInterfaces[i]) out << fInterfaces[i]->Index() << " ";
        out << std::endl;
    }
};


/**
 * @brief Class that manages the generation of the simworx mesh
 * @author Geo-Adapta (Philippe Devloo more precisely)
 */
class TRMSimworxMeshGenerator
{
private:
    
    /// Auxilary geomesh
    TPZGeoMesh *m_auxGMesh;
    
    /// Ribs of the mesh
    TPZStack<TWellRibs> fRibs;
    
public:
    
    TRMSimworxMeshGenerator();
    TRMSimworxMeshGenerator(const TRMSimworxMeshGenerator &copy)
    {
        DebugStop();
    }
    TRMSimworxMeshGenerator &operator=(const TRMSimworxMeshGenerator &copy)
    {
        DebugStop();
        return *this;
    }
    ~TRMSimworxMeshGenerator();
    
    TPZGeoEl * CreateBCGeoBlendEl(TPZGeoEl *orig, int side, int bc);
    
    void CreateReservoirFaces(TPZGeoMesh & gmesh, StructMioloData & mioloData);
    
    void CheckNodesSequenceForExtrusion(TPZGeoMesh * gmesh, TPZVec<int64_t> &nodeIndexes);
    
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
    
    TPZGeoMesh * ReallyGenerateGeoMesh(const REAL semiAxeX,
                                                     const REAL semiAxeY,
                                                     const REAL mioloLx,
                                                     const TPZVec<REAL> & espacamentoMioloY,
                                                     const TPZVec<REAL> & espacamentoZ,
                                                     const bool thereIsCutPlane,
                                                     TRMRawData &rawdata);
    TPZGeoMesh * CreateSimworxGeoMesh(TRMRawData &rawdata, bool withwellbc);
    
    /**
     * Os cornerNodes do retangulo inscrito na elipse sao movidos para um ponto que
     * gera malha com melhor aspecto dos prismas para situacoes diversas de geometrias.
     *
     */
    void AdjustPrismCoordinates(TPZGeoMesh * gmesh, REAL semiX, REAL semiY);
    
    void Pair_Miolo_Reserv_Nodes(TPZGeoMesh * reservGMesh, const std::set<int> & reservNodeIndices,
                                 TPZAutoPointer<TPZGeoMesh> mioloGMesh,
                                 std::map<int,int> & miolo_reserv_nodeIndices);
    
    void AddRibElements(TPZGeoMesh *gmesh, int WellMatId1D, int WellMatFake1D);
    
    void CreateWellBoundaries(TPZAutoPointer<TPZGeoMesh> reservoirGMesh);
    
};


#endif /* defined(__PZ__TRMSimworxMeshGenerator__) */
