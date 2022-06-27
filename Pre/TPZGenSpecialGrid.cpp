/**
 * @file
 * @brief Contains the implementation of the methods to generate polygonal grids approximating three dimensional special surfaces. 
 */

#include "TPZGenSpecialGrid.h"

#include "pzmanvector.h"
#include "pzstack.h"
#include "pzgeoel.h"
#include "tpzgeoblend.h"
#include "tpzgeoelrefpattern.h"
#include "tpzarc3d.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"
#include "TPZVTKGeoMesh.h"

#include <math.h>
#include <fstream>
//#define BLEND_VERBOSE //outputs x and grad comparisons
#define BLEND_OUTPUT_TXT//prints all elements in .txt format
#define BLEND_OUTPUT_VTK//prints all elements in .vtk format


/** 
 * Function to generate a polygonal mesh approximating a sphere from a octahedron mesh (polygonal) based on number of uniform refinements
 */
TPZGeoMesh *TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(TPZVec<REAL> &center, REAL radius, int nUniformRefs) {
    // Initial mesh data
    // In this case is a octahedron
    // Octahedron has 6 nodes. We considerer that the nodes is on axes and radius is the distance between any node and origin.
    const int64_t nnode = 6;
    // Octahedron has 8 (triangular) faces
    const int64_t nelem = 8;
    
    // Initial nodes and initial triangular faces of the octahedron
    REAL initialcoord[nnode][3] = {{center[0]-radius,center[1],center[2]},{center[0],center[1]+radius,center[2]},{center[0],center[1],center[2]-radius},{center[0],center[1],center[2]+radius},{center[0],center[1]-radius,center[2]},{center[0]+radius,center[1],center[2]}};
    int indices[nelem][nnode] = {{3,4,5},{3,5,1},{3,1,0},{3,0,4},{4,0,2},{4,2,5},{2,0,1},{5,2,1}};
    
    // Geometric element vector
    TPZGeoEl *elvec[nelem], *gel;
    
    // Creating geometric initial mesh
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    // Initializing nodes of the polygonal initial mesh
    int64_t node;
    for(node=0; node<nnode; node++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL> coord(3);
        coord[0] = initialcoord[node][0];
        coord[1] = initialcoord[node][1];
        coord[2] = initialcoord[node][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(node,coord,*gmesh);
    }
    // Creating triangular elements
    int64_t el, index;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(3);
        for(node=0; node<3; node++) nodind[node]=indices[el][node];
        elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
    }
    // Building connectivity for initial mesh - OCTAHEDRON
    gmesh->BuildConnectivity();
    
    //Loop making uniform refinement and changing coordinates of the nodes (projecting into the sphere) until tolerance is reached
    TPZManVector<REAL> baryparam(3,0.), barycenter(3,0.);
    for(int ii=0;ii<nUniformRefs;ii++) {
        // Make a uniform refinement
        UniformRefinement(1,gmesh,2);
        
        // Projecting all nodes into the sphere
        TPZVec<REAL> coordinates(3,0.);
        REAL norm;
        int i;
        for(node=0;node<gmesh->NNodes();node++) {
            TPZGeoNode *gnode = &gmesh->NodeVec()[node];
            gnode->GetCoordinates(coordinates);
            for(i=0;i<3;i++)
                coordinates[i] -= center[i];
            norm = sqrt(coordinates[0]*coordinates[0] + coordinates[1]*coordinates[1] + coordinates[2]*coordinates[2]);
            for(i=0;i<3;i++) coordinates[i] /= norm;
            for(i=0;i<3;i++) coordinates[i] += center[i];
            gnode->SetCoord(coordinates);
        }
    }
    
    // Before return, all the "father" elements are deleted, but the "son" elements are not
    for(el=0;el<gmesh->NElements();el++) {
        gel = gmesh->ElementVec()[el];
        if(!gel) continue;
        if(gel->HasSubElement()) {
            gel->ResetSubElements();
            gmesh->DeleteElement(gel);
        }
        else
            gel->SetFatherIndex(-1);
    }
    // The connectivity are updated
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    
    return gmesh;
}

/** 
 * Function to generate a polygonal mesh approximating a sphere from a octahedron mesh (polygonal) 
 */
TPZGeoMesh *TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(TPZVec<REAL> &center,REAL radius, REAL tol) {
    // Initial mesh data
    // In this case is a octahedron
    // Octahedron has 6 nodes. We considerer that the nodes is on axes and radius is the distance between any node and origin.
    const int64_t nnode = 6;
    // Octahedron has 8 (triangular) faces
    const int64_t nelem = 8;
    
    // Initial nodes and initial triangular faces of the octahedron
    REAL initialcoord[nnode][3] = {{center[0]-radius,center[1],center[2]},{center[0],center[1]+radius,center[2]},{center[0],center[1],center[2]-radius},{center[0],center[1],center[2]+radius},{center[0],center[1]-radius,center[2]},{center[0]+radius,center[1],center[2]}};
    int indices[nelem][nnode] = {{3,4,5},{3,5,1},{3,1,0},{3,0,4},{4,0,2},{4,2,5},{2,0,1},{5,2,1}};
    
    // Geometric element vector
    TPZGeoEl *elvec[nelem], *gel;
    
    // Creating geometric initial mesh
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    // Initializing nodes of the polygonal initial mesh
    int64_t node;
    for(node=0; node<nnode; node++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL> coord(3);
        coord[0] = initialcoord[node][0];
        coord[1] = initialcoord[node][1];
        coord[2] = initialcoord[node][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(node,coord,*gmesh);
    }
    // Creating triangular elements
    int64_t el, index;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(3);
        for(node=0; node<3; node++) nodind[node]=indices[el][node];
        elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
    }
    // Building connectivity for initial mesh - OCTAHEDRON
    gmesh->BuildConnectivity();
    
    //Loop making uniform refinement and changing coordinates of the nodes (projecting into the sphere) until tolerance is reached
    TPZManVector<REAL> baryparam(3,0.), barycenter(3,0.);
    REAL dist = 1.;
    while(tol < dist) {
        // Make a uniform refinement
        UniformRefinement(1,gmesh,2);
        
        // Projecting all nodes into the sphere
        TPZVec<REAL> coordinates(3,0.);
        REAL norm;
        int i;
        for(node=0;node<gmesh->NNodes();node++) {
            TPZGeoNode *gnode = &gmesh->NodeVec()[node];
            gnode->GetCoordinates(coordinates);
            for(i=0;i<3;i++)
                coordinates[i] -= center[i];
            norm = sqrt(coordinates[0]*coordinates[0] + coordinates[1]*coordinates[1] + coordinates[2]*coordinates[2]);
            for(i=0;i<3;i++) coordinates[i] /= norm;
            for(i=0;i<3;i++) coordinates[i] += center[i];
            gnode->SetCoord(coordinates);
        }
        
        // Find element no null to calculate distance between barycenter and projection into sphere
        for(el=0;el<gmesh->NElements();el++) {
            gel = gmesh->ElementVec()[el];
            if(!gel || gel->HasSubElement()) continue;
            gel->CenterPoint(6,baryparam);
            gel->X(baryparam,barycenter);
            break;
        }
        dist = radius - sqrt((barycenter[0]-center[0])*(barycenter[0]-center[0]) + (barycenter[1]-center[1])*(barycenter[1]-center[1]) + (barycenter[2]-center[2])*(barycenter[2]-center[2]));
        if(dist < 0.0) {
            delete gmesh;
            gmesh = 0;
            return gmesh;
        } 
    }
    
    // Before return, all the "father" elements are deleted, but the "son" elements are not
    for(el=0;el<gmesh->NElements();el++) {
        gel = gmesh->ElementVec()[el];
        if(!gel) continue;
        if(gel->HasSubElement()) {
            gel->ResetSubElements();
            gmesh->DeleteElement(gel);
        }
        else
            gel->SetFatherIndex(-1);
    }
    // The connectivity are updated
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    
    return gmesh;
}

/** 
 * Make uniform refinement of the geometric mesh. This method finishs making the new connectivities for the elements (ResetConnectivities() and BuildConnectivity()) 
 */
void TPZGenSpecialGrid::UniformRefinement(int nUniformRefs,TPZGeoMesh *gmesh, const int dimelfordivision, bool allmaterial, const int matidtodivided) {
    TPZManVector<TPZGeoEl*> filhos;
    for(int Division=0; Division<nUniformRefs; Division++)
    {
        int64_t nels = gmesh->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {    
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            if(!gel || gel->HasSubElement())
                continue;
            if(dimelfordivision > 0 && gel->Dimension() != dimelfordivision) continue;
            if(!allmaterial){
                if(gel->MaterialId() == matidtodivided){
                    gel->Divide(filhos);
                }
            }
            else{
                gel->Divide(filhos);
            }
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

/// create da 2D mesh corresponding to a circle divided in quads and triangles
TPZGeoMesh *TPZGenSpecialGrid::CreateGeoMesh2D_Circle(int nDiv)
{
#ifdef BLEND_VERBOSE
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    TPZVec<REAL> coord(3, 0.);
    for (int64_t i = 0; i < 4; i++) {
        coord[0] = i / 2 == i % 2 ? -1 : 1;
        coord[1] = -1 + 2 * (i / 2);
        //            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
        const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
    }//quad nodes
    
    TPZManVector<int64_t, 4> nodesIdArcsVec(4);
    for (int64_t i = 0; i < 4; i++) {
        const int xOff = i % 2;
        const int yOff = 1 - xOff;
        const int xSign = 1 - 2 * (i / 2);
        const int ySign = -1 * xSign;
        coord[0] = xSign * xOff * M_SQRT2;
        coord[1] = ySign * yOff * M_SQRT2;
        //            std::cout<<"xOff:\t"<<xOff<<"\txSign:"<<xSign<<std::endl;
        //            std::cout<<"yOff:\t"<<yOff<<"\tySign:"<<ySign<<std::endl;
        //            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
        const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
        nodesIdArcsVec[i] = newindex;
    }//quad nodes
    
    int matIdQuad = 1, matIdArc = 2;
    TPZManVector<int64_t, 4> nodesIdVec(1);
    int64_t elId = 0;
    
    nodesIdVec.Resize(4);
    TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > *quadEl = nullptr;
    {
        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1;
        nodesIdVec[2] = 2;
        nodesIdVec[3] = 3;
        quadEl = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(nodesIdVec, matIdQuad,
                                                                                  *gmesh);
    }
    
    
    TPZAutoPointer<TPZRefPattern> refp;
    char buf[] =
    "6 4 "
    "-98 QuadIntoAll "
    "-1 -1   0 "//0
    " 1 -1   0 "//1
    " 1  1   0 "//2
    "-1  1   0 "//3
    " 0 -1   0 "//4
    " 0  1   0 "//5
    "3 4  0 1 2 3 " //master quad
    "3 4  0 4 5 3 " //left half of the cube
    "2 3  4 1 2 "//triang1
    "2 3  4 2 5 "; //triang2
    
    std::istringstream str(buf);
    refp = new TPZRefPattern(str);
    refp->GenerateSideRefPatterns();
    gRefDBase.InsertRefPattern(refp);
    const int firstEdge = pztopology::TPZQuadrilateral::NumSides(0);
    const int nEdges = pztopology::TPZQuadrilateral::NumSides(1);
    quadEl->SetRefPattern(refp);
    TPZGeoElRefPattern<pzgeom::TPZArc3D> *arc = nullptr;
    for (int edge = firstEdge; edge < firstEdge + nEdges; edge++) {
        const int nNodes = pztopology::TPZQuadrilateral::NSideNodes(edge);
        nodesIdVec.Resize(nNodes + 1);
        for (int node = 0; node < nNodes; node++) {
            nodesIdVec[node] = pztopology::TPZQuadrilateral::SideNodeLocId(edge, node);
            //                std::cout<<nodesIdVec[node]<<"\t";
        }
        nodesIdVec[nNodes] = nodesIdArcsVec[edge - firstEdge];
        //            std::cout<<nodesIdVec[nNodes]<<std::endl;
        arc = new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
        auto arcRefp = refp->SideRefPattern(edge);
        arc->SetRefPattern(arcRefp);
    }
    gmesh->BuildConnectivity();
    
    {
        TPZVec<TPZGeoEl *> sons;
        const int nel = gmesh->NElements();
        quadEl->Divide(sons);
        for (int iel = 1; iel < nel; iel++) {
            TPZGeoEl *geo = gmesh->Element(iel);
            auto geoElSide = geo->Neighbour(geo->NSides() - 1);
            if (geoElSide.NSubElements() > 1) {
                geo->Divide(sons);
            }
        }
    }
    
    //Create GeoBlend elements
    {
        TPZGeoMesh *newgmesh = new TPZGeoMesh();
        for (int64_t i = 0; i < gmesh->NNodes(); i++) {
            coord[0] = gmesh->NodeVec()[i].Coord(0);
            coord[1] = gmesh->NodeVec()[i].Coord(1);
            coord[2] = gmesh->NodeVec()[i].Coord(2);
            const int64_t newindex = newgmesh->NodeVec().AllocateNewElement();
            newgmesh->NodeVec()[newindex].Initialize(coord, *newgmesh);
        }
        //TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >
        const int nel = gmesh->NElements();
        newgmesh->ElementVec().Resize(nel);
        for (int iel = 0; iel < nel; iel++) {
            newgmesh->ElementVec()[iel] = nullptr;
            TPZGeoEl *geo = gmesh->Element(iel);
            nodesIdVec.resize(geo->NNodes());
            for (int i = 0; i < nodesIdVec.size(); i++) {
                nodesIdVec[i] = geo->NodeIndex(i);
            }
            const int matId = geo->MaterialId();
            const int elId = geo->Index();
            //                geo->Ind
            switch (geo->Type()) {
                case ETriangle:
                    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(elId, nodesIdVec, matId,
                                                                                        *newgmesh);
                    break;
                case EQuadrilateral:
                    if (geo->HasSubElement()) continue;//the original cube should not be copied
                    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>(elId, nodesIdVec, matId,
                                                                                    *newgmesh);
                    break;
                default:
                    geo->Clone(*newgmesh);
                    break; //the element should not be deleted
            }
            //TPZGeoElRefPattern(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
        }
        {
            TPZVec<TPZGeoEl *> sons;
            std::vector<std::string> loading = {"-","/","|","\\"};
            for (int iDiv = 0; iDiv < nDiv; iDiv++) {
                std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
                const int nel = gmesh->NElements();
                for (int iel = 0; iel < nel; iel++) {
                    std::cout<<"\b"<<loading[iel%4]<<std::flush;
                    TPZGeoEl *geo = gmesh->ElementVec()[iel];
                    if (geo && !geo->HasSubElement()) {
                        geo->Divide(sons);
                    }
                }
                std::cout<<"\b";
            }
        }

        delete gmesh;
        gmesh = newgmesh;
    }
    
    gmesh->SetDimension(2);
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName =
        "blendmesh2D_partial.txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef BLEND_OUTPUT_VTK
    {
        const std::string meshFileName =
        "blendmesh2D_partial.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
    return gmesh;
}

/// create a 2D circle segment mesh made of a single quadrilateral
TPZGeoMesh *TPZGenSpecialGrid::CreateArcMesh(REAL inner, REAL outer, REAL alpha0, REAL alpha1)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    const int matIdQuad = 1;
    const int matIdArc = 2;
    TPZManVector<REAL,3> coords(3);
    int64_t newindex = -1;
    const int64_t nNodes = 6;
    for(int i = 0; i < nNodes; i++){
        const REAL r = ((i+1)/2) % 2 ? outer : inner;
        const REAL aux1 = i/2 ? 1.0 : 0.0;
        const REAL aux2 = (i/2)/2 ? 1.0 : 0.0;
        const REAL theta = alpha0+(aux1 - 0.5 * aux2)*(alpha1-alpha0);
        //            std::cout <<"aux1 = "<<aux1<<"aux2 = "<<aux2<<std::endl;
        coords[0] = r * cos(theta);
        coords[1] = r * sin(theta);
        coords[2] = 0.0;
        newindex = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newindex].Initialize(coords, *gmesh);
    }
    
    TPZManVector<int64_t,3> nodesIdVec(4);
    //nodesIdVec={0,1,2,3};//TODO:use this assignment when merged with master branch
    nodesIdVec[0]=0;
    nodesIdVec[1]=1;
    nodesIdVec[2]=2;
    nodesIdVec[3]=3;
    auto *quad = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>(nodesIdVec,matIdQuad, *gmesh);
    TPZGeoElRefPattern<pzgeom::TPZArc3D> *arc = nullptr;
    nodesIdVec.Resize(3);
    //nodesIdVec = {0,3,4};//TODO:use this assignment when merged with master branch
    nodesIdVec[0]=0;
    nodesIdVec[1]=3;
    nodesIdVec[2]=4;
    arc = new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    //nodesIdVec = {1,2,5};//TODO:use this assignment when merged with master branch
    nodesIdVec[0]=1;
    nodesIdVec[1]=2;
    nodesIdVec[2]=5;
    arc = new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    gmesh->BuildConnectivity();
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName = "gmesh_circle.txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef BLEND_OUTPUT_VTK
    {
        const std::string meshFileName = "gmesh_circle.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
    return gmesh;

}

/// create mesh of a sphere composed of all 3D topologies and mapped to a sphere
TPZGeoMesh *TPZGenSpecialGrid::CreateGeoMesh3D_DividedSphere(int nDiv)
{
#ifdef BLEND_VERBOSE
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    const int64_t nNodesHexa = 8;
    const int64_t nNodes= nNodesHexa;//+nNodesArc;
    
    const REAL radius = std::sqrt(3);
    const REAL cubeSide= 2;
    TPZVec<REAL> coord(3, 0.);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///                                                                                                         ////
    ////the face containing the nodes 0,1,2 and 3 will be joined with a half-sphere                             ////
    ///                                                                                                         ////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int64_t i = 0; i < nNodesHexa; i++) {
        int iAux = i%4;
        bool xFactor = (bool) !( (iAux % 2) == (iAux / 2) );
        bool yFactor = (bool) ( i / 4 );
        bool zFactor = (bool) ( ( iAux /2 ) % 2 );
        coord[0] = 0.5 * cubeSide * (2 * (int)xFactor - 1) ;
        coord[1] = 0.5 * cubeSide  * (1 - 2 * (int)yFactor);
        coord[2] = 0.5 * cubeSide * (2 * (int)zFactor - 1) ;
        const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
        //            std::cout<<"x : "<<xFactor<<"\ty : "<<yFactor<<"\tz : "<<zFactor<<std::endl;
        //            std::cout<<"x : "<<coord[0]<<"\ty : "<<coord[1]<<"\tz : "<<coord[2]<<std::endl;
    }
    
    
    int matIdVol = 1, matIdSphere = 2;
    TPZVec<int64_t> nodesIdVec(1);
    int64_t elId = 0;
    
    nodesIdVec.Resize(8);
    
    nodesIdVec[0] = 0;
    nodesIdVec[1] = 1;
    nodesIdVec[2] = 2;
    nodesIdVec[3] = 3;
    nodesIdVec[4] = 4;
    nodesIdVec[5] = 5;
    nodesIdVec[6] = 6;
    nodesIdVec[7] = 7;
    TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > *hexaEl =
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(nodesIdVec, matIdVol, *gmesh);
    TPZAutoPointer<TPZRefPattern> refp;
    char buf[] =
    "12 5 "
    "-99 CubeIntoAll "
    "-1 -1  -1 "
    " 1 -1  -1 "
    " 1  1  -1 "
    "-1  1  -1 "
    "-1 -1   0 "
    " 1 -1   0 "
    " 1  1   0 "
    "-1  1   0 "
    "-1 -1   1 "
    " 1 -1   1 "
    " 1  1   1 "
    "-1  1   1 "
    "7 8  0 1 2 3 8 9 10 11 " //master cube
    "7 8  4 5 6 7 8 9 10 11 " //upper half of the cube
    "6 6  0 5 4 3 6 7 "//prism
    "5 5  0 3 6 5 2 " //pyramid
    "4 4  1 2 5 0 "; //tetra
    
    std::istringstream str(buf);
    refp = new TPZRefPattern(str);
    refp->GenerateSideRefPatterns();
    gRefDBase.InsertRefPattern(refp);
    const int firstFace = pztopology::TPZCube::NumSides(0) + pztopology::TPZCube::NumSides(1);
    const int nFaces = pztopology::TPZCube::NumSides(2);
    hexaEl->SetRefPattern(refp);
    TPZGeoElRefPattern<pzgeom::TPZQuadSphere<>> *sphere = nullptr;
    for(int face = firstFace; face < firstFace + nFaces; face++){
        const int nNodes = pztopology::TPZCube::NSideNodes(face);
        nodesIdVec.Resize(nNodes);
        for (int node = 0; node < nNodes; node++){
            nodesIdVec[node] = pztopology::TPZCube::SideNodeLocId(face,node);
        }
        sphere = new TPZGeoElRefPattern<pzgeom::TPZQuadSphere<>>(nodesIdVec, matIdSphere, *gmesh);
        coord[0] = coord[1] = coord[2] = 0;
        sphere->Geom().SetData(radius, coord);
        auto sphereRefp = refp->SideRefPattern(face);
        sphere->SetRefPattern(sphereRefp);
    }
    gmesh->BuildConnectivity();
    
    {
        TPZVec<TPZGeoEl *> sons;
        const int nel = gmesh->NElements();
        // divide element zero (the cube)
        hexaEl->Divide(sons);
        for (int iel = 1; iel < nel; iel++) {
            TPZGeoEl *geo = gmesh->Element(iel);
            auto geoElSide = geo->Neighbour(geo->NSides()-1);
            // geoElSide is the cube
            if(geoElSide.NSubElements() > 1)
            {
                // divide the sphere elements in accordance with the divided cube
                geo->Divide(sons);
            }
        }
    }
    // the subelements of the cube and spheres are geoelmapped elements
    // force the 3d subelements to blend elements, copy the surface elements
    //Create GeoBlend elements
    {
        TPZGeoMesh *newgmesh = new TPZGeoMesh();
        for (int64_t i = 0; i < gmesh->NNodes(); i++) {
            coord[0] = gmesh->NodeVec()[i].Coord(0);
            coord[1] = gmesh->NodeVec()[i].Coord(1);
            coord[2] = gmesh->NodeVec()[i].Coord(2);
            const int64_t newindex = newgmesh->NodeVec().AllocateNewElement();
            newgmesh->NodeVec()[newindex].Initialize(coord,*newgmesh);
        }
        //TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >
        const int nel = gmesh->NElements();
        newgmesh->ElementVec().Resize(nel);
        for (int iel = 0; iel < nel; iel++) {
            newgmesh->ElementVec()[iel] = nullptr;
            TPZGeoEl *geo = gmesh->Element(iel);
            nodesIdVec.resize(geo->NNodes());
            for(int i = 0; i < nodesIdVec.size(); i++){
                nodesIdVec[i] = geo->NodeIndex(i);
            }
            const int matId = geo->MaterialId();
            const int elId = geo->Index();
            //                geo->Ind
            switch(geo->Type()){
                case ETetraedro:
                    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTetrahedra>>(elId,nodesIdVec,matId,*newgmesh);
                    break;
                case EPiramide:
                    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoPyramid>>(elId,nodesIdVec,matId,*newgmesh);
                    break;
                case EPrisma:
                    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoPrism>>(elId,nodesIdVec,matId,*newgmesh);
                    break;
                case ECube:
                    if(geo->HasSubElement()) continue;//the original cube should not be copied
                    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube>>(elId,nodesIdVec,matId,*newgmesh);
                    break;
                default:
                    auto geoclone = geo->Clone(*newgmesh);
                    break; //the element should not be deleted
            }
            //TPZGeoElRefPattern(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
        }
        delete gmesh;
        gmesh = newgmesh;
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName =
        "blendmesh3D_partial.txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef BLEND_OUTPUT_VTK
    {
        const std::string meshFileName =
        "blendmesh3D_partial.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
    {
        TPZVec<TPZGeoEl *> sons;
        std::vector<std::string> loading = {"-","/","|","\\"};
        for (int iDiv = 0; iDiv < nDiv; iDiv++) {
#ifdef BLEND_OUTPUT_TXT
            {
                const std::string meshFileName =
                "blendmesh3D_ref" + std::to_string(iDiv) + ".txt";
                std::ofstream outTXT(meshFileName.c_str());
                gmesh->Print(outTXT);
                outTXT.close();
            }
#endif
            std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
            const int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                std::cout<<"\b"<<loading[iel%4]<<std::flush;
                TPZGeoEl *geo = gmesh->ElementVec()[iel];
                if (geo && !geo->HasSubElement()) {
                    geo->Divide(sons);
                }
            }
            std::cout<<"\b";
        }
    }
#ifdef BLEND_OUTPUT_TXT
    {
        const std::string meshFileName =
        "blendmesh3D.txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef BLEND_OUTPUT_VTK
    {
        const std::string meshFileName =
        "blendmesh3D.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
    
    return gmesh;
}

