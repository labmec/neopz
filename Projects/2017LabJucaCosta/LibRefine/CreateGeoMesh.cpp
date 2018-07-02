//
//  SimilarFunctionsIntoProjects.cpp
//  PZ
//
//  Created by labmec on 08/10/17.
//
//

#include "CreateAndRefineMeshes.h"

/* 1. Functions contructing geometrical meshes.
 Projects:
    Poisson3D_Shock
 */

/*************************************************/
/**** Creating Geometric Mesh in cube or square **/
/*************************************************/
TPZGeoMesh *CreateGeomMesh(int typeel,int materialId,int id_bc0,int id_bc1,int id_bc2) {
    TPZManVector<REAL> point(3,0.), pointlast(3,0.);
    TPZGeoMesh* gmesh;
    switch (typeel) {
        case EOned:
        {
            pointlast[0] = 1.;
            gmesh = new TPZGeoMesh;
            int Qnodes = 2;
            
            gmesh->SetMaxNodeId(Qnodes-1);
            gmesh->NodeVec().Resize(Qnodes);
            TPZVec<TPZGeoNode> Node(Qnodes);
            
            TPZVec <int64_t> TopolLine(2);
            TPZVec <int64_t> TopolPoint(1);
            
            int64_t id = 0;
            for (int j=0; j<2;j++) {
                Node[id].SetNodeId(id);
                if(!j) Node[id].SetCoord(point);//coord x
                else Node[id].SetCoord(pointlast);
                gmesh->NodeVec()[id] = Node[id];
                id++;
            }
            
            //indice dos elementos
            id = 0;
            
            TopolLine[0] = 0;
            TopolLine[1] = 1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,materialId,*gmesh);
            id++;
            
            TopolPoint[0] = 0;
            new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,id_bc0,*gmesh);
            id++;
            TopolPoint[0] = 1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,id_bc0,*gmesh);
            
            gmesh->SetDimension(1);
            gmesh->ResetConnectivities();
            gmesh->BuildConnectivity();
        }
            break;
        case EQuadrilateral:
        {
            const int nnode = 4;
            REAL co[nnode][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};
            
            TPZGeoEl *elvec[1];
            gmesh = new TPZGeoMesh();
            gmesh->SetDimension(2);
            int nod;
            for(nod=0; nod<nnode; nod++) {
                int nodind = gmesh->NodeVec().AllocateNewElement();
                TPZVec<REAL> coord(3);
                coord[0] = co[nod][0];
                coord[1] = co[nod][1];
                coord[2] = co[nod][2];
                gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
            }
            
            int el = 0;
            TPZVec<int64_t> nodind(nnode);
            for(nod=0; nod<nnode; nod++) nodind[nod]=nod;
            int64_t index;
            elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,materialId,index);
            
            gmesh->BuildConnectivity();
            
            // bc -1 -> Dirichlet
            TPZGeoElBC gbc1(elvec[0],4,id_bc0);
            TPZGeoElBC gbc2(elvec[0],5,id_bc0);
            TPZGeoElBC gbc3(elvec[0],6,id_bc0);
            TPZGeoElBC gbc4(elvec[0],7,id_bc0);
            gmesh->ResetConnectivities();
            gmesh->BuildConnectivity();
        }
            break;
        case ETriangle:
        {
            const int nelem = 4;
            const int nnode = 5;
            
            REAL co[nnode][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.},{0.5,0.5}};
            int indices[nelem][nnode] = {{0,1,4},{1,2,4},{2,3,4},{0,3,4}};
            
            TPZGeoEl *elvec[nelem];
            gmesh = new TPZGeoMesh();
            gmesh->SetDimension(2);
            
            int nod;
            for(nod=0; nod<nnode; nod++) {
                int nodind = gmesh->NodeVec().AllocateNewElement();
                TPZVec<REAL> coord(3,0.0);
                coord[0] = co[nod][0];
                coord[1] = co[nod][1];
                gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
            }
            
            int el;
            for(el=0; el<nelem; el++) {
                TPZVec<int64_t> nodind(3);
                for(nod=0; nod<3; nod++) nodind[nod]=indices[el][nod];
                int64_t index;
                elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
            }
            
            gmesh->BuildConnectivity();
            
            // bc -1 -> Dirichlet
            TPZGeoElBC gbc1(elvec[0],3,id_bc0);
            TPZGeoElBC gbc2(elvec[1],3,id_bc0);
            TPZGeoElBC gbc3(elvec[2],3,id_bc0);
            TPZGeoElBC gbc4(elvec[3],3,id_bc0);
            gmesh->ResetConnectivities();
            gmesh->BuildConnectivity();
        }
            break;
        case ETetraedro:
            gmesh = ConstructingTetrahedraInCube(1.,materialId,id_bc0,id_bc1);
            break;
        case EPrisma:
            gmesh = ConstructingPrismsInCube(1.,materialId,id_bc0,id_bc1);
            break;
        case EPiramide:
            gmesh = ConstructingPyramidsInCube(1.,materialId,id_bc0,id_bc1);
            break;
        case ECube:
            gmesh = ConstructingPositiveCube(1.,typeel,materialId,id_bc0,id_bc1);
            break;
        default:
            gmesh = 0;
            break;
    }
    
    return gmesh;
}

#include "TPZRefPatternDataBase.h"
TPZGeoMesh *ConstructingPositiveCube(REAL InitialL,int typeel,int materialId,int id_bc0,int id_bc1,int id_bc2) {
    // CREATING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // Dependig on dimension of the typeel
    const int nelem = 1;
    const int nnode = 8;
    REAL co[nnode][3] = {
        {0.,0.,0.},
        {InitialL,0.,0.},
        {InitialL,InitialL,0.},
        {0.,InitialL,0.},
        {0.,0.,InitialL},
        {InitialL,0.,InitialL},
        {InitialL,InitialL,InitialL},
        {0.,InitialL,InitialL}
    };
    TPZVec<TPZVec<int64_t> > indices(nelem);
    indices[0].Resize(nnode);
    int nod;
    for(nod=0;nod<nnode;nod++)
        indices[0][nod] = nod;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    gmesh->SetDimension(3);
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(ECube,indices[el],1,index);
    }
    gmesh->BuildConnectivity();
    
    // Introduzing boundary condition for cube - ALL DIRICHLET
    // Boundary condition on one dimensional sides
    //	TPZGeoElBC gbc00(gmesh->ElementVec()[0],8,id_bc0);
    //	TPZGeoElBC gbc01(gmesh->ElementVec()[0],9,id_bc0);
    //	TPZGeoElBC gbc02(gmesh->ElementVec()[0],10,id_bc0);
    //	TPZGeoElBC gbc03(gmesh->ElementVec()[0],11,id_bc0);
    //	TPZGeoElBC gbc04(gmesh->ElementVec()[0],12,id_bc0);
    //	TPZGeoElBC gbc05(gmesh->ElementVec()[0],13,id_bc0);
    //	TPZGeoElBC gbc06(gmesh->ElementVec()[0],14,id_bc0);
    //	TPZGeoElBC gbc07(gmesh->ElementVec()[0],15,id_bc0);
    //	TPZGeoElBC gbc08(gmesh->ElementVec()[0],16,id_bc0);
    //	TPZGeoElBC gbc09(gmesh->ElementVec()[0],17,id_bc0);
    //	TPZGeoElBC gbc10(gmesh->ElementVec()[0],18,id_bc0);
    //	TPZGeoElBC gbc11(gmesh->ElementVec()[0],19,id_bc0);
    // face 0 (20) bottom XY - face 1 (21) lateral left XZ - face 4 (24) lateral back YZ : Dirichlet
    TPZGeoElBC gbc20(gmesh->ElementVec()[0],20,id_bc0);
    TPZGeoElBC gbc21(gmesh->ElementVec()[0],21,id_bc0);
    TPZGeoElBC gbc22(gmesh->ElementVec()[0],24,id_bc0);
    // face 2 (22) Neumann - Partial derivative (du/dx) - lateral front
    TPZGeoElBC gbc23(gmesh->ElementVec()[0],22,id_bc0);
    // face 3 (23) Neumann - Partial derivative (du/dy) - lateral right
    TPZGeoElBC gbc24(gmesh->ElementVec()[0],23,id_bc0);
    // face 5 (25) Neumann - Partial derivative (du/dz) - top
    TPZGeoElBC gbc25(gmesh->ElementVec()[0],25,id_bc0);
    
    TPZVec<TPZGeoEl *> sub;
    std::string filename = "D:\\";
    switch (typeel) {
        case ENoType:
        {
            char buf[1024];
            std::istringstream str(buf);
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
            TPZAutoPointer<TPZRefPattern> refpatFound = gRefDBase.FindRefPattern(refpat);
            if(!refpatFound)
            {
                gRefDBase.InsertRefPattern(refpat);
            }
            else
            {
                refpatFound->SetName(refpat->Name());
            }
            refpat->InsertPermuted();
        }
            break;
            
        case EPrisma:   // hexahedron -> four prisms
        {
            // Dividing hexahedron in four prisms (anymore)
            filename += "/3D_Hexa_directional_2faces.rpt";
            
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
            TPZAutoPointer<TPZRefPattern> refpatFound = gRefDBase.FindRefPattern(refpat);
            if(!refpatFound)
            {
                gRefDBase.InsertRefPattern(refpat);
            }
            else
            {
                refpatFound->SetName(refpat->Name());
            }
            refpat->InsertPermuted();
            TPZGeoEl *gel = gmesh->ElementVec()[0];
            TPZGeoElRefPattern <TPZGeoCube> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
            gelrp->SetRefPattern(refpat);
            gel->Divide(sub);
        }
            break;
        case EPiramide:
        {
            // Dividing hexahedron in four pyramids (anymore)
            filename += "/3D_Hexa_Rib_Side_08.rpt";
            
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
            if(!gRefDBase.FindRefPattern(refpat))
            {
                gRefDBase.InsertRefPattern(refpat);
            }
            TPZGeoEl *gel = gmesh->ElementVec()[0];
            TPZGeoElRefPattern <TPZGeoCube> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
            gelrp->SetRefPattern(refpat);
            gel->Divide(sub);
        }
            break;
        case ETetraedro:
        {
            // Dividing hexahedron in two tetrahedras, two prisms and one pyramid
            filename += "/3D_Hexa_Rib_Side_16_17_18.rpt";
            
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
            if(!gRefDBase.FindRefPattern(refpat))
            {
                gRefDBase.InsertRefPattern(refpat);
            }
            TPZGeoEl *gel = gmesh->ElementVec()[0];
            TPZGeoElRefPattern <TPZGeoCube> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
            gelrp->SetRefPattern(refpat);
            gel->Divide(sub);
        }
            break;
        case ECube:
            break;
        default:
        {
            // hexahedron -> three prisms (anymore)
            // Dividing hexahedron in prisms
            filename += "/3D_Hexa_Rib_Side_16_18.rpt";
            
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
            if(!gRefDBase.FindRefPattern(refpat))
            {
                gRefDBase.InsertRefPattern(refpat);
            }
            TPZGeoEl *gel = gmesh->ElementVec()[0];
            TPZGeoElRefPattern <TPZGeoCube> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
            gelrp->SetRefPattern(refpat);
            gel->Divide(sub);
        }
            break;
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    
    return gmesh;
}

TPZGeoMesh *ConstructingTetrahedraInCube(REAL InitialL,int materialId,int id_bc0,int id_bc1,int id_bc2) {
    // CONSIDERING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // And dividing into five tetrahedras
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(3);
    const int nelem = 5;
    const int nnode = 8;
    REAL co[nnode][3] = {
        {0.,0.,0.},
        {InitialL,0.,0.},
        {InitialL,InitialL,0.},
        {0.,InitialL,0.},
        {0.,0.,InitialL},
        {InitialL,0.,InitialL},
        {InitialL,InitialL,InitialL},
        {0.,InitialL,InitialL},
    };
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    TPZVec<TPZVec<int64_t> > indices(nelem);
    int nnodebyelement = 4;
    int el;
    for(el=0;el<nelem;el++)
        indices[el].Resize(nnodebyelement);
    // nodes to first element
    indices[0][0] = 0;
    indices[0][1] = 1;
    indices[0][2] = 3;
    indices[0][3] = 4;
    // nodes to second element
    indices[1][0] = 1;
    indices[1][1] = 2;
    indices[1][2] = 3;
    indices[1][3] = 6;
    // nodes to third element
    indices[2][0] = 4;
    indices[2][1] = 5;
    indices[2][2] = 6;
    indices[2][3] = 1;
    // nodes to fourth element
    indices[3][0] = 6;
    indices[3][1] = 7;
    indices[3][2] = 4;
    indices[3][3] = 3;
    // nodes to fifth element
    indices[4][0] = 1;
    indices[4][1] = 4;
    indices[4][2] = 6;
    indices[4][3] = 3;
    
    TPZGeoEl *elvec[nelem];
    for(el=0; el<nelem; el++) {
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(ETetraedro,indices[el],materialId,index);
    }
    gmesh->BuildConnectivity();
    
    // Introduzing boundary condition for cube - ALL DIRICHLET
    // boundary condition over one dimensional sides
    TPZGeoElBC gbc01(gmesh->ElementVec()[0],4,id_bc0);
    TPZGeoElBC gbc02(gmesh->ElementVec()[0],5,id_bc0);
    TPZGeoElBC gbc03(gmesh->ElementVec()[0],6,id_bc0);
    TPZGeoElBC gbc04(gmesh->ElementVec()[0],7,id_bc0);
    TPZGeoElBC gbc05(gmesh->ElementVec()[0],8,id_bc0);
    TPZGeoElBC gbc06(gmesh->ElementVec()[0],9,id_bc0);
    TPZGeoElBC gbc07(gmesh->ElementVec()[1],5,id_bc0);
    TPZGeoElBC gbc08(gmesh->ElementVec()[1],6,id_bc0);
    TPZGeoElBC gbc09(gmesh->ElementVec()[1],7,id_bc0);
    TPZGeoElBC gbc10(gmesh->ElementVec()[1],8,id_bc0);
    TPZGeoElBC gbc11(gmesh->ElementVec()[1],9,id_bc0);
    TPZGeoElBC gbc12(gmesh->ElementVec()[2],4,id_bc0);
    TPZGeoElBC gbc13(gmesh->ElementVec()[2],5,id_bc0);
    TPZGeoElBC gbc14(gmesh->ElementVec()[2],6,id_bc0);
    TPZGeoElBC gbc15(gmesh->ElementVec()[2],8,id_bc0);
    TPZGeoElBC gbc16(gmesh->ElementVec()[3],4,id_bc0);
    TPZGeoElBC gbc17(gmesh->ElementVec()[3],5,id_bc0);
    TPZGeoElBC gbc18(gmesh->ElementVec()[3],8,id_bc0);
    
    // face 0 (20) bottom XY
    TPZGeoElBC gbc20(gmesh->ElementVec()[0],10,id_bc0);
    TPZGeoElBC gbc21(gmesh->ElementVec()[0],11,id_bc0);
    TPZGeoElBC gbc22(gmesh->ElementVec()[0],13,id_bc0);
    TPZGeoElBC gbc23(gmesh->ElementVec()[1],10,id_bc0);
    TPZGeoElBC gbc24(gmesh->ElementVec()[1],11,id_bc0);
    TPZGeoElBC gbc25(gmesh->ElementVec()[1],12,id_bc0);
    TPZGeoElBC gbc26(gmesh->ElementVec()[2],10,id_bc0);
    TPZGeoElBC gbc27(gmesh->ElementVec()[2],11,id_bc0);
    TPZGeoElBC gbc28(gmesh->ElementVec()[2],12,id_bc0);
    TPZGeoElBC gbc29(gmesh->ElementVec()[3],10,id_bc0);
    TPZGeoElBC gbc30(gmesh->ElementVec()[3],11,id_bc0);
    TPZGeoElBC gbc31(gmesh->ElementVec()[3],12,id_bc0);
    return gmesh;
}

TPZGeoMesh *ConstructingPyramidsInCube(REAL InitialL,int materialId,int id_bc0,int id_bc1,int id_bc2) {
    // CONSIDERING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // And dividing into six pyramids
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    const int nelem = 6;
    const int nnode = 9;
    REAL co[nnode][3] = {
        {0.,0.,0.},
        {InitialL,0.,0.},
        {InitialL,InitialL,0.},
        {0.,InitialL,0.},
        {0.,0.,InitialL},
        {InitialL,0.,InitialL},
        {InitialL,InitialL,InitialL},
        {0.,InitialL,InitialL},
        {0.5*InitialL,0.5*InitialL,0.5*InitialL}
    };
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    TPZVec<TPZVec<int64_t> > indices(nelem);
    int nnodebyelement = 5;
    int el;
    for(el=0;el<nelem;el++)
        indices[el].Resize(nnodebyelement);
    // nodes to first element
    indices[0][0] = 0;
    indices[0][1] = 1;
    indices[0][2] = 2;
    indices[0][3] = 3;
    indices[0][4] = 8;
    // nodes to second element
    indices[1][0] = 0;
    indices[1][1] = 1;
    indices[1][2] = 5;
    indices[1][3] = 4;
    indices[1][4] = 8;
    // nodes to third element
    indices[2][0] = 1;
    indices[2][1] = 2;
    indices[2][2] = 6;
    indices[2][3] = 5;
    indices[2][4] = 8;
    // nodes to fourth element
    indices[3][0] = 2;
    indices[3][1] = 3;
    indices[3][2] = 7;
    indices[3][3] = 6;
    indices[3][4] = 8;
    // nodes to fifth element
    indices[4][0] = 0;
    indices[4][1] = 3;
    indices[4][2] = 7;
    indices[4][3] = 4;
    indices[4][4] = 8;
    // nodes to sixth element
    indices[5][0] = 4;
    indices[5][1] = 5;
    indices[5][2] = 6;
    indices[5][3] = 7;
    indices[5][4] = 8;
    
    TPZGeoEl *elvec[nelem];
    for(el=0; el<nelem; el++) {
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(EPiramide,indices[el],materialId,index);
    }
    gmesh->BuildConnectivity();
    
    // Introduzing boundary condition for cube - ALL DIRICHLET
    // boundary condition on one dimensional sides
    TPZGeoElBC gbc01(gmesh->ElementVec()[0],5,id_bc0);
    TPZGeoElBC gbc02(gmesh->ElementVec()[0],6,id_bc0);
    TPZGeoElBC gbc03(gmesh->ElementVec()[0],7,id_bc0);
    TPZGeoElBC gbc04(gmesh->ElementVec()[0],8,id_bc0);
    TPZGeoElBC gbc05(gmesh->ElementVec()[1],6,id_bc0);
    TPZGeoElBC gbc06(gmesh->ElementVec()[1],7,id_bc0);
    TPZGeoElBC gbc07(gmesh->ElementVec()[1],8,id_bc0);
    TPZGeoElBC gbc08(gmesh->ElementVec()[2],6,id_bc0);
    TPZGeoElBC gbc09(gmesh->ElementVec()[2],7,id_bc0);
    TPZGeoElBC gbc10(gmesh->ElementVec()[3],6,id_bc0);
    TPZGeoElBC gbc11(gmesh->ElementVec()[3],7,id_bc0);
    TPZGeoElBC gbc12(gmesh->ElementVec()[4],7,id_bc0);
    
    // face 0 (20) bottom XY
    TPZGeoElBC gbc21(gmesh->ElementVec()[0],13,id_bc0);
    TPZGeoElBC gbc22(gmesh->ElementVec()[1],13,id_bc0);
    TPZGeoElBC gbc23(gmesh->ElementVec()[2],13,id_bc0);
    TPZGeoElBC gbc24(gmesh->ElementVec()[3],13,id_bc0);
    TPZGeoElBC gbc25(gmesh->ElementVec()[4],13,id_bc0);
    TPZGeoElBC gbc26(gmesh->ElementVec()[5],13,id_bc0);
    
    return gmesh;
}

TPZGeoMesh *ConstructingPrismsInCube(REAL InitialL,int materialId,int id_bc0,int id_bc1,int id_bc2) {
    // CREATING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // And dividing into four prisms
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    gmesh->SetDimension(3);
    const int nelem = 4;
    const int nnode = 12;
    REAL co[nnode][3] = {
        {0.,0.,0.},
        {InitialL,0.,0.},
        {InitialL,InitialL,0.},
        {0.,InitialL,0.},
        {0.,0.,InitialL},
        {InitialL,0.,InitialL},
        {InitialL,InitialL,InitialL},
        {0.,InitialL,InitialL},
        {0.,0.,0.5*InitialL},
        {InitialL,0.,0.5*InitialL},
        {InitialL,InitialL,0.5*InitialL},
        {0.,InitialL,0.5*InitialL},
    };
    for(int nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    TPZVec<TPZVec<int64_t> > indices(nelem);
    int nnodebyelement = 6;
    int el;
    for(el=0;el<nelem;el++)
        indices[el].Resize(nnodebyelement);
    // nodes to first element
    indices[0][0] = 0;
    indices[0][1] = 1;
    indices[0][2] = 2;
    indices[0][3] = 8;
    indices[0][4] = 9;
    indices[0][5] = 10;
    // nodes to second element
    indices[1][0] = 0;
    indices[1][1] = 2;
    indices[1][2] = 3;
    indices[1][3] = 8;
    indices[1][4] = 10;
    indices[1][5] = 11;
    // nodes to third element
    indices[2][0] = 8;
    indices[2][1] = 9;
    indices[2][2] = 10;
    indices[2][3] = 4;
    indices[2][4] = 5;
    indices[2][5] = 6;
    // nodes to fourth element
    indices[3][0] = 8;
    indices[3][1] = 10;
    indices[3][2] = 11;
    indices[3][3] = 4;
    indices[3][4] = 6;
    indices[3][5] = 7;
    
    TPZGeoEl *elvec[nelem];
    for(el=0; el<nelem; el++) {
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(EPrisma,indices[el],materialId,index);
    }
    gmesh->BuildConnectivity();
    
    // Introduzing boundary condition for cube - ALL DIRICHLET
    // boundary condition on one dimensional sides
    TPZGeoElBC gbc00(gmesh->ElementVec()[0],6,id_bc0);
    TPZGeoElBC gbc01(gmesh->ElementVec()[0],7,id_bc0);
    TPZGeoElBC gbc02(gmesh->ElementVec()[0],8,id_bc0);
    TPZGeoElBC gbc03(gmesh->ElementVec()[0],9,id_bc0);
    TPZGeoElBC gbc04(gmesh->ElementVec()[0],10,id_bc0);
    TPZGeoElBC gbc05(gmesh->ElementVec()[0],11,id_bc0);
    TPZGeoElBC gbc06(gmesh->ElementVec()[0],12,id_bc0);
    TPZGeoElBC gbc07(gmesh->ElementVec()[0],13,id_bc0);
    TPZGeoElBC gbc08(gmesh->ElementVec()[1],7,id_bc0);
    TPZGeoElBC gbc09(gmesh->ElementVec()[1],8,id_bc0);
    TPZGeoElBC gbc10(gmesh->ElementVec()[1],11,id_bc0);
    TPZGeoElBC gbc11(gmesh->ElementVec()[1],13,id_bc0);
    TPZGeoElBC gbc12(gmesh->ElementVec()[1],14,id_bc0);
    TPZGeoElBC gbc13(gmesh->ElementVec()[2],9,id_bc0);
    TPZGeoElBC gbc14(gmesh->ElementVec()[2],10,id_bc0);
    TPZGeoElBC gbc15(gmesh->ElementVec()[2],11,id_bc0);
    TPZGeoElBC gbc16(gmesh->ElementVec()[2],12,id_bc0);
    TPZGeoElBC gbc17(gmesh->ElementVec()[2],13,id_bc0);
    TPZGeoElBC gbc18(gmesh->ElementVec()[2],14,id_bc0);
    TPZGeoElBC gbc19(gmesh->ElementVec()[3],11,id_bc0);
    TPZGeoElBC gbc20(gmesh->ElementVec()[3],13,id_bc0);
    TPZGeoElBC gbc21(gmesh->ElementVec()[3],14,id_bc0);
    
    // face 0 (20) bottom XY
    TPZGeoElBC gbc30(gmesh->ElementVec()[0],15,id_bc0);
    TPZGeoElBC gbc31(gmesh->ElementVec()[0],16,id_bc0);
    TPZGeoElBC gbc32(gmesh->ElementVec()[0],17,id_bc0);
    TPZGeoElBC gbc33(gmesh->ElementVec()[1],15,id_bc0);
    TPZGeoElBC gbc34(gmesh->ElementVec()[1],17,id_bc0);
    TPZGeoElBC gbc35(gmesh->ElementVec()[1],18,id_bc0);
    TPZGeoElBC gbc36(gmesh->ElementVec()[2],16,id_bc0);
    TPZGeoElBC gbc37(gmesh->ElementVec()[2],17,id_bc0);
    TPZGeoElBC gbc38(gmesh->ElementVec()[2],19,id_bc0);
    TPZGeoElBC gbc39(gmesh->ElementVec()[3],17,id_bc0);
    TPZGeoElBC gbc40(gmesh->ElementVec()[3],18,id_bc0);
    TPZGeoElBC gbc41(gmesh->ElementVec()[3],19,id_bc0);
    
    return gmesh;
}
TPZGeoMesh *ConstructingSeveral3DElementsInCube(REAL InitialL,MElementType typeel,int id_bc0,int id_bc1,int id_bc2) {
    // CREATING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // Dependig on dimension of the typeel
    const int nelem = 1;
    const int nnode = 8;
    REAL co[nnode][3] = {
        {0.,0.,0.},
        {InitialL,0.,0.},
        {InitialL,InitialL,0.},
        {0.,InitialL,0.},
        {0.,0.,InitialL},
        {InitialL,0.,InitialL},
        {InitialL,InitialL,InitialL},
        {0.,InitialL,InitialL}
    };
    TPZVec<TPZVec<int64_t> > indices(nelem);
    indices[0].Resize(nnode);
    int nod;
    for(nod=0;nod<nnode;nod++)
        indices[0][nod] = nod;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(ECube,indices[el],1,index);
    }
    gmesh->BuildConnectivity();
    
    // Introduzing boundary condition for cube - ALL DIRICHLET
    // face 0 (20) bottom XY - face 1 (21) lateral left XZ - face 4 (24) lateral back YZ : Dirichlet
    TPZGeoElBC gbc10(gmesh->ElementVec()[0],20,id_bc0);
    TPZGeoElBC gbc11(gmesh->ElementVec()[0],21,id_bc0);
    TPZGeoElBC gbc12(gmesh->ElementVec()[0],24,id_bc0);
    // face 2 (22) Neumann - Partial derivative (du/dx) - lateral front
    TPZGeoElBC gbc13(gmesh->ElementVec()[0],22,id_bc0);
    // face 3 (23) Neumann - Partial derivative (du/dy) - lateral right
    TPZGeoElBC gbc14(gmesh->ElementVec()[0],23,id_bc0);
    // face 5 (25) Neumann - Partial derivative (du/dz) - top
    TPZGeoElBC gbc15(gmesh->ElementVec()[0],25,id_bc0);
    
    return gmesh;
}

//*******Shell to deforming************
TPZGeoMesh *CreateGeomMesh(std::string &archivo) {
    
    // Generacion de una malla utilizando o GID previamente 
    TPZReadGIDGrid grid;
    TPZGeoMesh *meshgrid = grid.GeometricGIDMesh(archivo);
    if(!meshgrid->NElements())
        return 0;
    
    return meshgrid;
}

void PrintNRefinementsByType(int nref, int64_t nels,int64_t newnels,TPZVec<int64_t> &counter,std::ostream &out) {
    out << "\n HP Refinement done, on  " << nels << " elements, given " << newnels << " elements. "<< std::endl;
    out << " NRef = " << nref << std::endl;
    for(int j=0;j<counter.NElements();j++)
        if(counter[j]) {
            out << " Refinement type " << j << " : " << counter[j] << std::endl;
        }
    out << " Processed elements " << (nels-counter[0]);
}

int MaxLevelReached(TPZCompMesh *cmesh) {
    int levelreached = 0, level;
    TPZCompEl *el = 0;
    for(int i=0;i<cmesh->NElements();i++) {
        el = cmesh->ElementVec()[i];
        if(!el || el->Dimension() != cmesh->Dimension()) continue;
        level = el->Reference()->Level();
        levelreached = (level > levelreached) ? level : levelreached;
    }
    return levelreached;
}

