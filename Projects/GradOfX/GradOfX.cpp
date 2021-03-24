#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZRefPatternTools.h"


#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"

#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticpyramid.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#include "pzgmesh.h"

#include "fad.h"

#include "TPZVTKGeoMesh.h"

//Linear shape functions
TPZGeoMesh * OneDimensional();
TPZGeoMesh * TwoDimensionalT();
TPZGeoMesh * TwoDimensionalQ();
TPZGeoMesh * ThreeDimensionalT();
TPZGeoMesh * ThreeDimensionalH();
TPZGeoMesh * ThreeDimensionalPr();
TPZGeoMesh * ThreeDimensionalPy();

//Quadratic shape functions
TPZGeoMesh * OneDimensionalQuadratic();
TPZGeoMesh * TwoDimensionalTQuadratic();
TPZGeoMesh * TwoDimensionalQQuadratic();
TPZGeoMesh * ThreeDimensionalTQuadratic();
TPZGeoMesh * ThreeDimensionalHQuadratic();
TPZGeoMesh * ThreeDimensionalPrQuadratic();
TPZGeoMesh * ThreeDimensionalPyQuadratic();




void ComputeGradofX(TPZGeoMesh * mesh, std::string file_name, TPZFMatrix<REAL> &triplets);


int main()
{
    std::string file_name;
    TPZFMatrix<REAL> triplets;
    
    
    {
        TPZGeoMesh * mesh = OneDimensional();
        file_name = "x_and_grad_1D.txt";
        triplets.Resize(3,3);
        triplets.Zero();
        triplets(0,0) = -1.0;
        triplets(1,0) = +0.0;
        triplets(2,0) = +1.0;
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = OneDimensionalQuadratic();
        file_name = "x_and_grad_1D_Quadratic.txt";
        triplets.Resize(3,3);
        triplets.Zero();
        triplets(0,0) = -1.0;
        triplets(1,0) = +0.0;
        triplets(2,0) = +1.0;
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = TwoDimensionalT();
        file_name = "x_and_grad_2D_T.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = -0.0;
        triplets(0,1) = +0.0;
        triplets(1,0) = +1.0;
        triplets(1,1) = +0.0;
        triplets(2,0) = -0.0;
        triplets(2,1) = +1.0;
        triplets(3,0) = +0.25;
        triplets(3,1) = +0.25;
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = TwoDimensionalTQuadratic();
        file_name = "x_and_grad_2D_T_Quadratic.txt";
        triplets.Resize(6,3);
        triplets.Zero();
        triplets(0,0) = -0.0;
        triplets(0,1) = +0.0;
        triplets(1,0) = +1.0;
        triplets(1,1) = +0.0;
        triplets(2,0) = -0.0;
        triplets(2,1) = +1.0;
        triplets(3,0) = 0.5;
        triplets(3,1) = 0.0;
        triplets(4,0) = 0.5;
        triplets(4,1) = 0.5;
        triplets(5,0) = 0.0;
        triplets(5,1) = 0.5;
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = TwoDimensionalQ();
        file_name = "x_and_grad_2D_Q.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = +0.25;
        triplets(0,1) = +0.25;
        triplets(1,0) = -0.25;
        triplets(1,1) = +0.25;
        triplets(2,0) = -0.25;
        triplets(2,1) = -0.25;
        triplets(3,0) = +0.25;
        triplets(3,1) = -0.25;
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = TwoDimensionalQQuadratic();
        file_name = "x_and_grad_2D_Q_Quadratic.txt";
        triplets.Resize(8,3);
        triplets.Zero();
        triplets(0,0) = -1.0;
        triplets(0,1) = -1.0;
        triplets(1,0) = +1.0;
        triplets(1,1) = -1.0;
        triplets(2,0) = +1.0;
        triplets(2,1) = +1.0;
        triplets(3,0) = -1.0;
        triplets(3,1) = +1.0;
        triplets(4,0) = 0.0;
        triplets(4,1) = -1.0;
        triplets(5,0) = 1.0;
        triplets(5,1) = 0.0;
        triplets(6,0) = 0.0;
        triplets(6,1) = 1.0;
        triplets(7,0) = -1.0;
        triplets(7,1) = 0.0;
        
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = ThreeDimensionalT();
        file_name = "x_and_grad_3D_T.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = +0.25;
        triplets(0,1) = +0.25;
        triplets(0,2) = +0.00;
        
        triplets(1,0) = +0.15;
        triplets(1,1) = +0.15;
        triplets(1,2) = +0.25;
        
        triplets(2,0) = +0.05;
        triplets(2,1) = +0.05;
        triplets(2,2) = +0.25;
        
        triplets(3,0) = +0.00;
        triplets(3,1) = +0.00;
        triplets(3,2) = +1.00;
        
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = ThreeDimensionalTQuadratic();
        file_name = "x_and_grad_3D_T_Quadratic.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = +0.25;
        triplets(0,1) = +0.25;
        triplets(0,2) = +0.00;
        
        triplets(1,0) = +0.15;
        triplets(1,1) = +0.15;
        triplets(1,2) = +0.25;
        
        triplets(2,0) = +0.05;
        triplets(2,1) = +0.05;
        triplets(2,2) = +0.25;
        
        triplets(3,0) = +0.00;
        triplets(3,1) = +0.00;
        triplets(3,2) = +1.00;
        
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = ThreeDimensionalH();
        file_name = "x_and_grad_3D_H.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = +0.25;
        triplets(0,1) = +0.25;
        triplets(0,2) = +0.00;
        
        triplets(1,0) = -0.15;
        triplets(1,1) = -0.15;
        triplets(1,2) = +0.25;
        
        triplets(2,0) = +0.05;
        triplets(2,1) = +0.05;
        triplets(2,2) = -0.25;
        
        triplets(3,0) = +1.00;
        triplets(3,1) = +1.00;
        triplets(3,2) = +1.00;
        
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = ThreeDimensionalHQuadratic();
        file_name = "x_and_grad_3D_H_Quadratic.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = +0.25;
        triplets(0,1) = +0.25;
        triplets(0,2) = +0.00;
        
        triplets(1,0) = -0.15;
        triplets(1,1) = -0.15;
        triplets(1,2) = +0.25;
        
        triplets(2,0) = +0.05;
        triplets(2,1) = +0.05;
        triplets(2,2) = -0.25;
        
        triplets(3,0) = +1.00;
        triplets(3,1) = +1.00;
        triplets(3,2) = +1.00;
        
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = ThreeDimensionalPr();
        file_name = "x_and_grad_3D_Pr.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = +0.25;
        triplets(0,1) = +0.25;
        triplets(0,2) = +0.00;
        
        triplets(1,0) = +0.15;
        triplets(1,1) = +0.15;
        triplets(1,2) = +0.25;
        
        triplets(2,0) = +0.05;
        triplets(2,1) = +0.05;
        triplets(2,2) = +0.25;
        
        triplets(3,0) = +0.00;
        triplets(3,1) = +0.00;
        triplets(3,2) = +1.00;
        
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = ThreeDimensionalPrQuadratic();
        file_name = "x_and_grad_3D_Pr_Quadratic.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = +0.25;
        triplets(0,1) = +0.25;
        triplets(0,2) = +0.00;
        
        triplets(1,0) = +0.15;
        triplets(1,1) = +0.15;
        triplets(1,2) = +0.25;
        
        triplets(2,0) = +0.05;
        triplets(2,1) = +0.05;
        triplets(2,2) = +0.25;
        
        triplets(3,0) = +0.00;
        triplets(3,1) = +0.00;
        triplets(3,2) = +1.00;
        
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = ThreeDimensionalPy();
        file_name = "x_and_grad_3D_Py.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = -0.25;
        triplets(0,1) = -0.25;
        triplets(0,2) = +0.00;
        
        triplets(1,0) = +0.15;
        triplets(1,1) = +0.15;
        triplets(1,2) = +0.25;
        
        triplets(2,0) = -0.05;
        triplets(2,1) = -0.05;
        triplets(2,2) = +0.25;
        
        triplets(3,0) = +0.00;
        triplets(3,1) = +0.00;
        triplets(3,2) = +0.99999;
        
        ComputeGradofX(mesh,file_name,triplets);
    }
    
    {
        TPZGeoMesh * mesh = ThreeDimensionalPyQuadratic();
        file_name = "x_and_grad_3D_Py_Quadratic.txt";
        triplets.Resize(4,3);
        triplets.Zero();
        triplets(0,0) = -0.25;
        triplets(0,1) = -0.25;
        triplets(0,2) = +0.00;
        
        triplets(1,0) = +0.15;
        triplets(1,1) = +0.15;
        triplets(1,2) = +0.25;
        
        triplets(2,0) = -0.05;
        triplets(2,1) = -0.05;
        triplets(2,2) = +0.25;
        
        triplets(3,0) = +0.00;
        triplets(3,1) = +0.00;
        triplets(3,2) = +0.99999;
        
        ComputeGradofX(mesh,file_name,triplets);
    }

    return 0;
}

void ComputeGradofX(TPZGeoMesh * mesh, std::string file_name, TPZFMatrix<REAL> &triplets){
    
    std::ofstream file(file_name.c_str());
    
    TPZManVector<REAL,3> triplet_xi_eta_zeta(3,0);
    TPZManVector<REAL,3> x(3,0);
    
    TPZVec<Fad<REAL> > Fad_triplet_xi_eta_zeta(3);
    
    
    TPZFMatrix<REAL> gradx(3,3,0.0);
    TPZFMatrix<Fad<REAL>> gradxFad(3,3,0.0);
    
    TPZFMatrix<REAL> coordinates;
    TPZFMatrix<REAL> jac;
    TPZFMatrix<REAL> axes;
    REAL detjac;
    TPZFMatrix<REAL> jacinv;
    
    int nplot_points = triplets.Rows();

    
    int nel = mesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel  = mesh->Element(iel);
        if (!gel) {
            std::cout << "not geometric element found." << std::endl;
            DebugStop();
        }
        
        gel->NodesCoordinates(coordinates);
        
        file << " ------------------------------------------------------ " << std::endl;
        file << " Geometric Element identifier = " << gel->Id() << std::endl;
        int iel_dim = gel->Dimension();
        
        
        for (int ip = 0; ip <  nplot_points; ip++) {
            triplet_xi_eta_zeta[0] = triplets(ip,0);
            triplet_xi_eta_zeta[1] = triplets(ip,1);
            triplet_xi_eta_zeta[2] = triplets(ip,2);
            
            for(int i = 0; i < iel_dim; i++){
                REAL val = triplets(ip,i);
                Fad<REAL> a(iel_dim,i,val);
                Fad_triplet_xi_eta_zeta[i] = a;
             
            }

            file << " triplet_xi_eta_zeta coordinate = " << triplet_xi_eta_zeta[0] << " ,  " << triplet_xi_eta_zeta[1] << " ,  " << triplet_xi_eta_zeta[2] << std::endl;
            file << " Fad_triplet_xi_eta_zeta coordinate = " << Fad_triplet_xi_eta_zeta[0] << " ,  " << Fad_triplet_xi_eta_zeta[1] << " ,  " << Fad_triplet_xi_eta_zeta[2] << std::endl;
            
            gel->X(triplet_xi_eta_zeta, x);
            file << " Mapped x coordinate = " << x[0] << " ,  " << x[1] << " ,  " << x[2] << std::endl;
            
            gel->GradX(triplet_xi_eta_zeta, gradx);
            file << " Grad of x = " << gradx << std::endl;
            
            int r = gradx.Rows();
            int c = gradx.Cols();
            gel->GradX(Fad_triplet_xi_eta_zeta, gradxFad);
            for(int i = 0; i < r; i++ ){
                for(int j = 0; j < c; j++ ){
                    gradx(i,j)=gradxFad(i,j).val();
                }
            }
            file << " Grad of x using Fad (GradXFad) = " << gradx << std::endl;
            
            TPZFMatrix<REAL> X_dx(r,c,0.0);
            TPZVec<Fad<REAL> > x(3);
            gel->X(Fad_triplet_xi_eta_zeta, x);
            for(int i = 0; i < r; i++ ){
                for(int j = 0; j < c; j++ ){
                    X_dx(i,j)=x[i].dx(j);
                }
            }
            file << " Grad of x using Fad (X.dx) = " << X_dx << std::endl;
        
            gel->Jacobian(triplet_xi_eta_zeta, jac, axes, detjac, jacinv);
            file << " axes = " << axes << std::endl;
            file << " jacobian = " << jac << std::endl;
            file << " detjac = " << detjac << std::endl;
            file << " *********************************** " << std::endl;
//            gel->JacobianXYZ(triplet_xi_eta_zeta, jac, axes, detjac, jacinv);
//            file << " axes = " << axes << std::endl;
//            file << " jacobian = " << jac << std::endl;
//            file << " detjac = " << detjac << std::endl;
//            file  << std::endl;
            
        }
        
    }
    
    file.flush();
    
    return;
    
}


TPZGeoMesh * OneDimensional(){
    
    int mat_id = 1;
    int64_t nodes  = 2;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolLine(2);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);     //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 2.0;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);     //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh1D.txt";
    std::string out_name_vtk = "GeometricMesh1D.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * OneDimensionalQuadratic(){
    
    int mat_id = 1;
    int64_t nodes  = 3;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolLine(2);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);     //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.25;
    z = 0.0;
    Node[id].SetCoord(0 , x);     //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);     //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elementid = 0;
    TPZVec < int64_t > nodeindex(3,0);
    
    nodeindex[0] = 0;
    nodeindex[1] = 2;
    nodeindex[2] = 1;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid , nodeindex, mat_id, *geo_mesh);
    elementid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh1D_Quadratic.txt";
    std::string out_name_vtk = "GeometricMesh1D_Quadratic.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}


TPZGeoMesh * TwoDimensionalT(){
    
    int mat_id = 1;
    int64_t nodes  = 3;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolTriangle(3);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolTriangle[0] = 0;
    TopolTriangle[1] = 1;
    TopolTriangle[2] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elid,TopolTriangle,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh2DT.txt";
    std::string out_name_vtk = "GeometricMesh2DT.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * TwoDimensionalTQuadratic(){
    
    int mat_id = 1;
    int64_t nodes  = 6;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolTriangle(6);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = -0.207107;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 1.0;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = -0.207107;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolTriangle[0] = 0;
    TopolTriangle[1] = 1;
    TopolTriangle[2] = 2;
    TopolTriangle[3] = 3;
    TopolTriangle[4] = 4;
    TopolTriangle[5] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZQuadraticTrig > (elid,TopolTriangle,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh2DT_Quadratic.txt";
    std::string out_name_vtk = "GeometricMesh2DT_Quadratic.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}


TPZGeoMesh * TwoDimensionalQ(){
    
    int mat_id = 1;
    int64_t nodes  = 4;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolQuadrilateral(4);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 0.0;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolQuadrilateral[0] = 0;
    TopolQuadrilateral[1] = 1;
    TopolQuadrilateral[2] = 2;
    TopolQuadrilateral[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elid,TopolQuadrilateral,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh2DQ.txt";
    std::string out_name_vtk = "GeometricMesh2DQ.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * TwoDimensionalQQuadratic(){
    
    int mat_id = 1;
    int64_t nodes  = 8;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolQuadrilateral(8);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 0.0;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = -0.207107;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.20711;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 1.20711;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = -0.207107;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolQuadrilateral[0] = 0;
    TopolQuadrilateral[1] = 1;
    TopolQuadrilateral[2] = 2;
    TopolQuadrilateral[3] = 3;
    TopolQuadrilateral[4] = 4;
    TopolQuadrilateral[5] = 5;
    TopolQuadrilateral[6] = 6;
    TopolQuadrilateral[7] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad > (elid,TopolQuadrilateral,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh2DQ_Quadratic.txt";
    std::string out_name_vtk = "GeometricMesh2DQ_Quadratic.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * ThreeDimensionalT(){
    
    int mat_id = 1;
    int64_t nodes  = 4;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolTetrahedron(4);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolTetrahedron[0] = 0;
    TopolTetrahedron[1] = 1;
    TopolTetrahedron[2] = 2;
    TopolTetrahedron[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra > (elid,TopolTetrahedron,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh3DT.txt";
    std::string out_name_vtk = "GeometricMesh3DT.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * ThreeDimensionalTQuadratic(){
    
    int mat_id = 1;
    int64_t nodes  = 10;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolTetrahedron(10);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.25;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.25;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
   
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.0;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.25;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    
    //  Geometric Elements
    int elid = 0;
    
    TopolTetrahedron[0] = 0;
    TopolTetrahedron[1] = 1;
    TopolTetrahedron[2] = 2;
    TopolTetrahedron[3] = 3;
    TopolTetrahedron[4] = 4;
    TopolTetrahedron[5] = 5;
    TopolTetrahedron[6] = 6;
    TopolTetrahedron[7] = 7;
    TopolTetrahedron[8] = 8;
    TopolTetrahedron[9] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZQuadraticTetra > (elid,TopolTetrahedron,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh3DT_Quadratic.txt";
    std::string out_name_vtk = "GeometricMesh3DT_Quadratic.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}


TPZGeoMesh * ThreeDimensionalH(){
    
    int mat_id = 1;
    int64_t nodes  = 8;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolHexahedron(8);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.5;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolHexahedron[0] = 0;
    TopolHexahedron[1] = 1;
    TopolHexahedron[2] = 2;
    TopolHexahedron[3] = 3;
    TopolHexahedron[4] = 4;
    TopolHexahedron[5] = 5;
    TopolHexahedron[6] = 6;
    TopolHexahedron[7] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (elid,TopolHexahedron,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh3DH.txt";
    std::string out_name_vtk = "GeometricMesh3DH.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * ThreeDimensionalHQuadratic(){
    
    int mat_id = 1;
    int64_t nodes  = 20;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolHexahedron(20);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);     //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.5;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.25;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.25;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.5;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.25;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.5;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.25;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;


    
    //  Geometric Elements
    int elid = 0;
    
    TopolHexahedron[0] = 0;
    TopolHexahedron[1] = 1;
    TopolHexahedron[2] = 2;
    TopolHexahedron[3] = 3;
    TopolHexahedron[4] = 4;
    TopolHexahedron[5] = 5;
    TopolHexahedron[6] = 6;
    TopolHexahedron[7] = 7;
    TopolHexahedron[8] = 8;
    TopolHexahedron[9] = 9;
    TopolHexahedron[10] = 10;
    TopolHexahedron[11] = 11;
    TopolHexahedron[12] = 12;
    TopolHexahedron[13] = 13;
    TopolHexahedron[14] = 14;
    TopolHexahedron[15] = 15;
    TopolHexahedron[16] = 16;
    TopolHexahedron[17] = 17;
    TopolHexahedron[18] = 18;
    TopolHexahedron[19] = 19;

    new TPZGeoElRefPattern< pzgeom::TPZQuadraticCube > (elid,TopolHexahedron,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh3DH_Quadratic.txt";
    std::string out_name_vtk = "GeometricMesh3DH_Quadratic.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * ThreeDimensionalPr(){
    
    int mat_id = 1;
    int64_t nodes  = 6;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolPrism(6);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolPrism[0] = 0;
    TopolPrism[1] = 1;
    TopolPrism[2] = 2;
    TopolPrism[3] = 3;
    TopolPrism[4] = 4;
    TopolPrism[5] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPrism > (elid,TopolPrism,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh3DPr.txt";
    std::string out_name_vtk = "GeometricMesh3DPr.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * ThreeDimensionalPrQuadratic(){
    
    int mat_id = 1;
    int64_t nodes  = 15;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolPrism(15);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.25;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.25;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.25;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.0;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.25;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.25;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolPrism[0] = 0;
    TopolPrism[1] = 1;
    TopolPrism[2] = 2;
    TopolPrism[3] = 3;
    TopolPrism[4] = 4;
    TopolPrism[5] = 5;
    TopolPrism[6] = 6;
    TopolPrism[7] = 7;
    TopolPrism[8] = 8;
    TopolPrism[9] = 9;
    TopolPrism[10] = 10;
    TopolPrism[11] = 11;
    TopolPrism[12] = 12;
    TopolPrism[13] = 13;
    TopolPrism[14] = 14;

    new TPZGeoElRefPattern< pzgeom::TPZQuadraticPrism > (elid,TopolPrism,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh3DPr_Quadratic.txt";
    std::string out_name_vtk = "GeometricMesh3DPr_Quadratic.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}


TPZGeoMesh * ThreeDimensionalPy(){
    
    int mat_id = 1;
    int64_t nodes  = 5;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolPyramid(5);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.25;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolPyramid[0] = 0;
    TopolPyramid[1] = 1;
    TopolPyramid[2] = 2;
    TopolPyramid[3] = 3;
    TopolPyramid[4] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid> (elid,TopolPyramid,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh3DPy.txt";
    std::string out_name_vtk = "GeometricMesh3DPy.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}

TPZGeoMesh * ThreeDimensionalPyQuadratic(){
    
    int mat_id = 1;
    int64_t nodes  = 13;
    
    TPZGeoMesh * geo_mesh= new TPZGeoMesh;
    
    geo_mesh->SetMaxNodeId(nodes-1);
    geo_mesh->NodeVec().Resize(nodes);
    TPZVec<TPZGeoNode> Node(nodes);
    
    TPZVec <int64_t> TopolPyramid(13);
    REAL x, y, z;
    
    // Nodes
    int64_t id = 0;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.25;
    y = 0.25;
    z = 1.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 0.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 1.0;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.5;
    y = 1.0;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.0;
    y = 0.5;
    z = 0.0;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.125;
    y = 0.125;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;

    Node[id].SetNodeId(id);
    x = 0.625;
    y = 0.125;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.625;
    y = 0.625;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    x = 0.125;
    y = 0.625;
    z = 0.5;
    Node[id].SetCoord(0 , x);         //coord x
    Node[id].SetCoord(1 , y);     //coord y
    Node[id].SetCoord(2 , z);     //coord z
    geo_mesh->NodeVec()[id] = Node[id];
    id++;
    
    //  Geometric Elements
    int elid = 0;
    
    TopolPyramid[0] = 0;
    TopolPyramid[1] = 1;
    TopolPyramid[2] = 2;
    TopolPyramid[3] = 3;
    TopolPyramid[4] = 4;
    TopolPyramid[5] = 5;
    TopolPyramid[6] = 6;
    TopolPyramid[7] = 7;
    TopolPyramid[8] = 8;
    TopolPyramid[9] = 9;
    TopolPyramid[10] = 10;
    TopolPyramid[11] = 11;
    TopolPyramid[12] = 12;

    new TPZGeoElRefPattern< pzgeom::TPZQuadraticPyramid> (elid,TopolPyramid,mat_id,*geo_mesh);
    elid++;
    
    geo_mesh->BuildConnectivity();
    
    std::string out_name_text = "GeometicMesh3DPy_Quadratic.txt";
    std::string out_name_vtk = "GeometricMesh3DPy_Quadratic.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    geo_mesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(geo_mesh,Dummyfile, true);
    
    return geo_mesh;
    
}
