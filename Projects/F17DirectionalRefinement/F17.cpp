/***********
 * This project used a F17 aircraft as a model to demonstrate
 * the capabilities of automatic directional refinement in NeoPZ.
 * These capabilities might be useful, for instance, in the case of
 * boundary layer analysis, hence the chosen example.
 * The program first builds a mesh modelling the **OUTSIDE** of the aircraft,
 * and the aircraft is part of the boundary of the created mesh.
 *
 * Then, through sucessive iterations, it models the boundary layer of the jet.
 */

#include "pzgmesh.h"//Mesh/pzgmesh.h
#include "TPZRefPatternDataBase.h"//Refine/TPZRefPatternDataBase.h
#include "TPZRefPatternTools.h"//Refine/TPZRefPatternTools.h
#include "TPZVTKGeoMesh.h"//Mesh/TPZVTKGeoMesh.h


#include <chrono>
/**
 * This function will create a TPZGeoMesh from a set of files containing the nodes, triangles and
 * tetrahedra forming the original F17 Mesh
 * @return
 */
TPZGeoMesh * CreateF17Mesh();

/**
 * This function will identify which bc elements are actually part of a reentrant portion of the mesh(F17 shell)
 * and which ones belong to the bounding box. It will assign, respectively, the material identifiers -1 and -4.
 * This method will only work if the bounding box is in the form $[a,b] \times [c,d] \times [e,f$, i.e., if
 * all it's faces are aligned with one of the directions of our coordinate system
 * @param geomesh
 */
void  FilterBoundingBox(TPZGeoMesh *geomesh);

int main()
{
    static const std::chrono::time_point<std::chrono::system_clock> wall_time_start = std::chrono::system_clock::now();

    gRefDBase.InitializeRefPatterns();

    TPZGeoMesh *geomesh = CreateF17Mesh();
    std::cout << "Filtering bounding box...\n";
    FilterBoundingBox(geomesh);

    const int nref = 5;
    std::cout << "Number of refinement steps: "<<nref<<std::endl;
    //cin >> nref;
    const std::string meshFileName("F17model_");
    std::ofstream outt(meshFileName+std::to_string(0)+".vtk");
    std::set<int> matids;
    matids.insert(-1);//material id associated with the F17 shell

    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outt, 1);


    int destMatId = 1;
    for(auto i=0; i<nref; i++)
    {
        int nelements = geomesh->NElements();
        std::cout << "\tEntering refinement step: " << i << " number of elements = " << nelements << '\n';
        std::cout << "\tRefining..."<<std::endl;
        //directional refinement will be performed in the direction of elements with material ids
        //contained in the set matids. elements that are refined will have the material id set to destmatid

        const auto reftimer_b = std::chrono::system_clock::now();
        for (auto el=0; el<nelements; el++)
        {
            TPZGeoEl *element = geomesh->ElementVec()[el];
            if(!element || element->MaterialId() != destMatId) continue;
            TPZRefPatternTools::RefineDirectional(element, matids, destMatId);
        }
        const auto reftimer_e = std::chrono::system_clock::now();
        std::cout<<"\tRefined!"<<'\n';
        std::cout <<"\tTotal time: ";
        std::cout<<std::chrono::duration_cast<std::chrono::milliseconds>(reftimer_e-reftimer_b).count();
        std::cout<<std::endl;
        //we want to differentiate which elements are refined in which round
        destMatId += 1;
        std::ofstream out(meshFileName+std::to_string(i+1)+".vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, 1);
    }
    std::cout << "Finalizing...\n";
    return 0 ;
}

TPZGeoMesh * CreateF17Mesh()
{
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    std::string path = PZSOURCEDIR;
    path += "/Projects/F17DirectionalRefinement/Files/";

    std::cout << "===============================================================\n"
              << "Reading F17 mesh\n";

    //First step: reading the nodes coordinates
    std::string nodesFileName = path;
    nodesFileName += "yf17.xyz";
    std::cout << "\t\tInput file for nodes = " << nodesFileName.c_str()
              << "\n\t\t\tprocessing nodes...\n";
    std::ifstream nodesFile (nodesFileName.c_str());
    if(!nodesFile.is_open()){
        std::cout<<"Could not find nodes at file "<<nodesFileName.c_str()<<'\n';
        std::cout<<"Aborting..."<<std::endl;
        DebugStop();
    }

    TPZManVector< REAL,3 > nodes(3);
    int i;
    while(nodesFile) {
        nodesFile >> nodes[0];
        nodesFile >> nodes[1];
        nodesFile >> nodes[2];
        if(!nodesFile) continue;
        int nodind = gmesh->NodeVec().AllocateNewElement();

        gmesh->NodeVec()[nodind] = TPZGeoNode(i,nodes,*gmesh);
    }

    // Second step: Generation of the triangles.
    std::string faceFileName = path;
    faceFileName += "yf17.tri";
    std::cout << "\t\tInput file for faces = " << faceFileName.c_str()
              << "\n\t\t\tprocessing faces...\n";
    std::ifstream facesFile (faceFileName.c_str());
    if(!facesFile.is_open()){
        std::cout<<"Could not find triangles at file "<<faceFileName.c_str()<<'\n';
        std::cout<<"Aborting..."<<std::endl;
        DebugStop();
    }
    TPZManVector<int64_t,3> indices(3);

    while (facesFile){
        facesFile >> indices[0] ;
        facesFile >> indices[1] ;
        facesFile >> indices[2] ;
        int64_t index;
        indices[0]--;
        indices[1]--;
        indices[2]--;
        if(!facesFile) continue;
        gmesh->CreateGeoElement(ETriangle,indices,-1,index,1);
    }

    // Third step: Generation of the tetrahedra.
    std::string tetraFileName = path;
    tetraFileName += "yf17.tet";
    std::cout << "\t\tInput file for volumes = " << tetraFileName.c_str()
              << "\n\t\t\tprocessing volumes...\n";
    std::ifstream tetraFile (tetraFileName.c_str());
    TPZManVector<int64_t,4> indices2(4);
    if(!tetraFile.is_open()){
        std::cout<<"Could not find volumes at file "<<tetraFileName.c_str()<<'\n';
        std::cout<<"Aborting..."<<std::endl;
        DebugStop();
    }
    //  int p;
    while (tetraFile){
        tetraFile >> indices2[0] ;
        tetraFile >> indices2[1] ;
        tetraFile >> indices2[2] ;
        tetraFile >> indices2[3] ;
        indices2[0]--;
        indices2[1]--;
        indices2[2]--;
        indices2[3]--;
        if(!tetraFile) continue;
        int64_t index;
        gmesh->CreateGeoElement(ETetraedro,indices2,1,index,1);
    }

    // Fourth step: Building the mesh.
    std::cout << "\n\tGenerating the connectivities and neighborhood information...\n";
    gmesh->BuildConnectivity();
    std::cout << "F17 mesh read\n"
              << "===============================================================\n";
    return gmesh;
}

void FilterBoundingBox(TPZGeoMesh *geomesh)
{
    constexpr REAL geoTol = 1e-4;

    const auto nnode = geomesh->NodeVec().NElements();
    TPZFNMatrix<6> xminmax(3,2,0);

    //identifying the maximum/minimum coordinates of the bounding box
    for(auto no=0; no<nnode; no++)
    {
        if(no == 0)
        {
            for(auto idf=0; idf<3; idf++)   xminmax(idf,0) = geomesh->NodeVec()[no].Coord(idf);
            for(auto idf=0; idf<3; idf++)   xminmax(idf,1) = geomesh->NodeVec()[no].Coord(idf);
        }
        for(auto idf=0; idf<3; idf++)
        {
            if(geomesh->NodeVec()[no].Coord(idf) < xminmax(idf,0)) xminmax(idf,0) = geomesh->NodeVec()[no].Coord(idf);
            if(geomesh->NodeVec()[no].Coord(idf) > xminmax(idf,1)) xminmax(idf,1) = geomesh->NodeVec()[no].Coord(idf);
        }
    }

    const auto nelem = geomesh->ElementVec().NElements();
    int num4 = 0;
    int num1 = 0;
    for(auto el=0; el<nelem; el++)
    {
        TPZGeoEl *gel = geomesh->ElementVec()[el];
        //just boundary elements for now
        if(!geomesh->ElementVec()[el] || gel->MaterialId() != -1) continue;

        //geting the min X and max X for one element
        TPZFNMatrix<9> xminmaxloc(3,2,0.);
        for(auto no=0; no<gel->NNodes(); no++)
        {
            TPZGeoNode *gno = gel->NodePtr(no);
            for(auto idf=0; idf<3; idf++)
            {
                if(no == 0)
                {
                    xminmaxloc(idf,0) = gno->Coord(idf);
                    xminmaxloc(idf,1) = gno->Coord(idf);
                }
                if(gno->Coord(idf) < xminmaxloc(idf,0)) xminmaxloc(idf,0) = gno->Coord(idf);
                if(gno->Coord(idf) > xminmaxloc(idf,1)) xminmaxloc(idf,1) = gno->Coord(idf);
            }
        }

        //calculates the difference between minX_loc (maxX_loc) and minX (maxX). this way it can detect if an
        //element has a node touching the bounding box
        xminmaxloc -= xminmax;
        REAL mindif = 1.;
        for(auto no=0; no<2; no++)
            for(auto idf=0; idf<3; idf++)
                mindif = (mindif < fabs(xminmaxloc(idf,no)))? mindif : fabs(xminmaxloc(idf,no));


        // Elements that does not contain a node in the bounding box are certainly part of the F17 Mesh.
        // So, they are skipped. Elements close to the boundary must be analysed
        if(mindif < geoTol)
        {
            //now we will find in which plane the element is contained to determine if it belongs to the
            //bounding box.
            TPZManVector<bool,3> dirMinDif(3,false);
            for(auto no=0; no<2; no++)
                for(auto idf=0; idf<3; idf++)
                    if(fabs(xminmaxloc(idf,no)) < geoTol) dirMinDif[idf] = true;
            bool isInFace{false};
            for(auto idf=0; idf<3; idf++)
            {
                if(!dirMinDif[idf]) continue;
                bool isInFaceIdf{true};
                REAL avgCoord = 0;
                for(auto no=0; no<gel->NNodes(); no++)
                {
                    if(!dirMinDif[idf]) continue;
                    TPZGeoNode *gno = gel->NodePtr(no);
                    const auto coord = gno->Coord(idf);
                    avgCoord = (avgCoord*no + coord )/(no+1);
                    isInFaceIdf = fabs(avgCoord-coord) < geoTol;
                    if(!isInFaceIdf) break;
                }
                isInFace = isInFace || isInFaceIdf;
            }
            if(isInFace)
            {
                gel->SetMaterialId(-4);
                num4++;
            }
            else
            {
                num1++;
            }

            //old way: using the element axes for the same purpose
//            TPZFNMatrix<9> axes(3,3),jac(2,2),jacinv(2,2);
//            REAL detjac;
//            TPZManVector<REAL,3> coor(2,0.3333);
//            gel->Jacobian(coor,jac,axes,detjac,jacinv);
//
//            TPZManVector<REAL,3> vectorialProd(3,0.);
//            vectorialProd[0] = -axes(0,2)*axes(1,1) + axes(0,1)*axes(1,2);
//            vectorialProd[1] = axes(0,2)*axes(1,0) - axes(0,0)*axes(1,2);
//            vectorialProd[2] = -axes(0,1)*axes(1,0) + axes(0,0)*axes(1,1);
//
//
//            TPZManVector<bool,3> dirMinDif(3,false);
//            for(auto no=0; no<2; no++)
//                for(auto idf=0; idf<3; idf++)
//                    if(fabs(xminmaxloc(idf,no)) < geoTol) dirMinDif[idf] = true;
//            //if they are really close to the boundary and their axes are aligned with a cartesian plane,
//            //they are not part of the F17 mesh
//            if(
//                    ((fabs (fabs(vectorialProd[0])-1.) < geoTol) && dirMinDif[0]) ||
//                    ((fabs (fabs(vectorialProd[1])-1.) < geoTol) && dirMinDif[1]) ||
//                    ((fabs (fabs(vectorialProd[2])-1.) < geoTol) && dirMinDif[2]) )
//            {
//                gel->SetMaterialId(-4);
//                num4++;
//            }
//            else
//            {
//                num1++;
//            }
        }
    }
    std::cout << "Number of elements with -4 condition " << num4 << " with -1 condition " << num1 << std::endl;
}