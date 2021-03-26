#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgnode.h"
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzmatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "time.h"
#include "pzconvectionproblem.h"
#include "pzmultiphase.h"
#include "pzl2projection.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"


#include <time.h>
#include <stdio.h>

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"
#include <cmath>



void Parametricfunction(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void Parametricfunction2(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void Parametricfunction3(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void BasicForm();
void GIDExtrusion();

bool ftriang = false;
REAL angle = 0.0*M_PI/4.0;

int main()
{   
//    BasicForm();
    GIDExtrusion();

 
    return 0;
}

void BasicForm(){
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
    GeoMesh1->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh1->NodeVec()[0]=Node;
    
    TPZVec<int64_t> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
    GeoMesh1->BuildConnectivity();
    GeoMesh1->SetDimension(0);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew1.txt");
        GeoMesh1->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh1,Dummyfile, true);
    }
    
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh1);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(Parametricfunction, 5);
    CreateGridFrom.SetParametricFunction(ParFunc);
    REAL t=0.0;
    REAL dt = 0.04;
    int n = 25.0;
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh2 = CreateGridFrom.ComputeExtrusion(t, dt, n);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew2.txt");
        GeoMesh2->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew2.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh2,Dummyfile, true);
    }
    
    
    
    TPZHierarquicalGrid CreateGridFrom2(GeoMesh2);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(Parametricfunction2, 5);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dt, n);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew3.txt");
        GeoMesh3->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew3.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh3,Dummyfile, true);
    }
    
    
    
    TPZHierarquicalGrid CreateGridFrom3(GeoMesh3);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc3 = new TPZDummyFunction<REAL>(Parametricfunction3, 5);
    CreateGridFrom3.SetParametricFunction(ParFunc3);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh4 = CreateGridFrom3.ComputeExtrusion(t, dt, n);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew4.txt");
        GeoMesh4->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew4.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh4,Dummyfile, true);
    }
    
}
void GIDExtrusion(){
    
   	std::string dirname = PZSOURCEDIR;
    //  Reading mesh
    std::string GridFileName;
    GridFileName = dirname + "/Projects/Gridgeneration/";
    //    GridFileName += "OilWaterSystemUnit.dump";
    GridFileName += "BaseGeometryDakeThin.dump";

    TPZReadGIDGrid GeometryInfo;
    GeometryInfo.SetfDimensionlessL(100.0);
    TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    gmesh->SetDimension(2);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicGIDMesh1.txt");
        gmesh->Print(argument);
        std::ofstream Dummyfile("GeometricGIDMesh1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    REAL t=0.0;
    REAL dt = 0.1;
    int n = 10.0;
    
    TPZHierarquicalGrid CreateGridFrom3(gmesh);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc3 = new TPZDummyFunction<REAL>(Parametricfunction3, 5);
    CreateGridFrom3.SetParametricFunction(ParFunc3);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh2 = CreateGridFrom3.ComputeExtrusion(t, dt, n);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicGIDMesh2.txt");
        GeoMesh2->Print(argument);
        std::ofstream Dummyfile("GeometricGIDMesh2.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh2,Dummyfile, true);
    }
    
}

void Parametricfunction(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = cos(par[0]);
    X[1] = sin(par[0]);
    X[2] = 0.0;
}

void Parametricfunction2(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = par[0];
    X[1] = par[0];
    X[2] = 0.0;
}

void Parametricfunction3(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;//cos(par[0]);
    X[2] = par[0];//par[0];
}
