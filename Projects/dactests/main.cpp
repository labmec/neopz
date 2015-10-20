
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"

#include "tpzchangeel.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPBSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "TPZReadGIDGrid.h"
#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

#include "pzhdivfull.h"
#include "pzelchdiv.h"

#include "pzgeopyramid.h"

#include "pznumeric.h"

#include "TPZExtendGridDimension.h"
#include "pzelchdivbound2.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"

#include "TPZLagrangeMultiplier.h"
#include "pzmatmixedpoisson3d.h"

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


#include "pyramidalmesh.h"

#include "tools.h"
#include "LaplaceInCylinder.h"
#include "LaplaceInCircle.h"
#include "LaplaceInSphere.h"
#include "LaplaceInQuadrilateral.h"
#include "LaplaceInCube.h"

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;
using namespace pzshape;

//int matId = 1;
//
//int dirichlet = 0;
//int neumann = 1;
//
//int bc0 = -1;
//int bc1 = -2;
//int bc2 = -3;
//int bc3 = -4;
//int bc4 = -5;
//int bc5 = -6;
//int matskeleton = -7;

//// just for print data
///** @brief Map used norms */
std::map<REAL,REAL> fDebugMapL2, fDebugMapHdiv;
///** @brief Map used Degrees of Freedom */
std::map<int,int> fDebugDoF;

//int tetraedra_2[6][4]=
//{
//    {1,2,5,4},
//    {4,7,3,2},
//    {0,1,2,4},
//    {0,2,3,4},
//    {4,5,6,2},
//    {4,6,7,2}
//};
//
//int piramide_2[6][5]=
//{
//    {0,1,2,3,8},
//    {0,1,5,4,8},
//    {1,2,6,5,8},
//    {3,2,6,7,8},
//    {0,3,7,4,8},
//    {4,5,6,7,8}
//};


//bool MyDoubleComparer(REAL a, REAL b);
//
//void GenerateNodes(TPZGeoMesh *gmesh, long nelem);
//void GenerateNodesPyramid(TPZGeoMesh *gmesh, long nelem);
//
//TPZGeoMesh *GMesh(int dimensao, bool ftriang, int ndiv);
//
//TPZGeoMesh *GMeshCirculoGeob(int dimensao, int ndiv);
//TPZGeoMesh *GMeshCirculoQuad(int dimensao, int ndiv);
//
//TPZGeoMesh *CreateOneCubo(int nref=0);
//TPZGeoMesh *CreateOneCuboWithTetraedrons(long nelem=1, int MaterialId=1);
//TPZGeoMesh *GMeshCubeWithPyramids(long nelem=1, int MaterialId=1);
//
//
//TPZGeoMesh * BasicForm(int n, REAL t, REAL dt);
//void Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X);
//void Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X);
//void Parametricfunction3(const TPZVec<STATE> &par, TPZVec<STATE> &X);

int dim = 3;
REAL aa = 0.0;
REAL bb = 0.0;
REAL cc = 0.0;
REAL Epsilon = 0.4;


//REAL const Pi = M_PI;//4.*atan(1.);



//bool ftriang = false;//true;//
bool IsCubedomain = true;
bool IsPrism = false;
bool IsTetra = false;
bool IsPiram = false;

//bool isspheredomain = true, iscircledomain = false, iscylinderdomain = false, isquaddomain = false;
//bool iscircledomain = true, isspheredomain = false, iscylinderdomain = false, isquaddomain = false;
bool iscylinderdomain = true, iscircledomain = false, isspheredomain = false, isquaddomain = false;
//bool isquaddomain = true, iscircledomain = false, isspheredomain = false, iscylinderdomain = false;


#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

#include "pztransfer.h"

int main(int argc, char *argv[])
{
//#ifdef LOG4CXX
//    InitializePZLOG();
//#endif
    
//  gRefDBase.InitializeAllUniformRefPatterns();
//	gRefDBase.InitializeRefPatterns();

    int p = 1;
    int ndiv = 0;
    HDivPiola = 0;
    ofstream saidaerros("../ErroNormas.txt",ios::app);
    
    for(p=1;p<4;p++)
    {
        saidaerros << "\nPARA p = " << p << endl;
        saidaerros << "ndiv " << setw(6) << "DoFT" << setw(20) << "DofCond" << setw(20) << "ErroL2Primal" << setw(20) << "ErroL2Dual" << setw(20) << "ErroL2Div" << setw(20) << "ErroHDivDual"  << endl;
        
        for (ndiv=1; ndiv<5/*7-p*/; ndiv++)
        {
            
            if (dim==2)
            {
                //TPZGeoMesh *gmesh2d = GMesh(2, ftriang, ndiv);
                if (iscircledomain) {
                    LaplaceInCircle  * circ = new LaplaceInCircle();
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    circ->Run( k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else if (iscylinderdomain)
                {
                    LaplaceInCylinder * cilind = new LaplaceInCylinder( );
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    cilind->Run( k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else  if(isspheredomain)
                {
                    LaplaceInSphere * sphere = new LaplaceInSphere( );
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    sphere->Run( k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else if(isquaddomain)
                {
                    LaplaceInQuadrilateral * quad = new LaplaceInQuadrilateral();
                    quad->setTriangTrue();
                    //quad->setH1True();
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    quad->Run(k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);
                }
                else
                {
                    //Ops! Nao esta fazendo nada.
                    DebugStop();
                }
            }
            else // dim == 3
            {
                if (IsCubedomain)
                {
                    LaplaceInCube * cubo = new LaplaceInCube();
                    //cubo->setTetraTrue();
                    //cubo->setPrismaTrue();
//                    cubo->setH1True();
                    bool HdivMaisMais = false;
                    int k = HdivMaisMais ? p+1 : p;
                    cubo->Run(k, ndiv, fDebugMapL2, fDebugMapHdiv, saidaerros, HdivMaisMais);

                }
                else if(IsTetra)
                {
                    DebugStop();
//                    REAL dndiv = ndiv;
//                    int nref = (int) pow(2., dndiv);
//                    
//                    gmesh = CreateOneCuboWithTetraedrons(nref, matId);
                }
                else if (IsPrism)
                {
                    // Tem que escolher a malha certa.
                    DebugStop();
                }
                else if(IsPiram)
                {
                    DebugStop();
//                    REAL dndiv = ndiv;
//                    int nref = (int) pow(2., dndiv);
//                    gmesh = GMeshCubeWithPyramids( nref);
//                    PyramidalMesh(gmesh,ndiv);
                }
                else
                {
                    // Nenhuma malha foi escolhida
                    DebugStop();
                }
            }
            

//            std::cout<< " grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
//            std::string filename("InputDataMeta");
//            std::string L2("L2.txt");
//            std::string Hdiv("Hdiv.txt");
//            std::string HdivData,L2Data;
//            HdivData = filename+str+Hdiv;
//            L2Data = filename+str+L2;
//            
//            PrintDebugMapForMathematica(HdivData,L2Data);
            
            std::cout<< " FIM - grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
            
        }
        saidaerros << "\n ----------------------------------------------------------------------------- " << endl;
        fDebugMapHdiv.clear();
        fDebugMapL2.clear();
    }
    
    std::cout<< " fim " << std::endl;
    
	return EXIT_SUCCESS;
}







//
//TPZGeoMesh * BasicForm(int n, REAL t, REAL dt)
//{
//    
//    // Creating a 0D element to be extruded
//    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
//    GeoMesh1->NodeVec().Resize(1);
//    TPZGeoNode Node;
//    TPZVec<REAL> coors(3,0.0);
//    Node.SetCoord(coors);
//    Node.SetNodeId(0);
//    GeoMesh1->NodeVec()[0]=Node;
//    
//    TPZVec<long> Topology(1,0);
//    int elid=0;
//    int matid=1;
//    
//    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
//    GeoMesh1->BuildConnectivity();
//    GeoMesh1->SetDimension(0);
//    {
//        //  Print Geometrical Base Mesh
//        std::ofstream argument("GeometicMeshNew1.txt");
//        GeoMesh1->Print(argument);
//        std::ofstream Dummyfile("GeometricMeshNew1.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh1,Dummyfile, true);
//    }
//    
//    
//    TPZHierarquicalGrid CreateGridFrom(GeoMesh1);
//    TPZAutoPointer<TPZFunction<STATE> > ParFunc = new TPZDummyFunction<STATE>(Parametricfunction);
//    CreateGridFrom.SetParametricFunction(ParFunc);
//    
//    // Computing Mesh extruded along the parametric curve Parametricfunction
//    TPZGeoMesh * GeoMesh2 = CreateGridFrom.ComputeExtrusion(t, dt, n);
//    
//    {
//        //  Print Geometrical Base Mesh
//        std::ofstream argument("GeometicMeshNew2.txt");
//        GeoMesh2->Print(argument);
//        std::ofstream Dummyfile("GeometricMeshNew2.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh2,Dummyfile, true);
//    }
//    
//    
//    
//    TPZHierarquicalGrid CreateGridFrom2(GeoMesh2);
//    TPZAutoPointer<TPZFunction<STATE> > ParFunc2 = new TPZDummyFunction<STATE>(Parametricfunction2);
//    CreateGridFrom2.SetParametricFunction(ParFunc2);
//    
//    // Computing Mesh extruded along the parametric curve Parametricfunction2
//    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dt, n);
//    {
//        //  Print Geometrical Base Mesh
//        std::ofstream argument("GeometicMeshNew3.txt");
//        GeoMesh3->Print(argument);
//        std::ofstream Dummyfile("GeometricMeshNew3.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh3,Dummyfile, true);
//    }
//    
//    
//    
//    TPZHierarquicalGrid CreateGridFrom3(GeoMesh3);
//    TPZAutoPointer<TPZFunction<STATE> > ParFunc3 = new TPZDummyFunction<STATE>(Parametricfunction3);
//    CreateGridFrom3.SetParametricFunction(ParFunc3);
//    
//    // Computing Mesh extruded along the parametric curve Parametricfunction2
//    TPZGeoMesh * GeoMesh4 = CreateGridFrom3.ComputeExtrusion(t, dt, n);
//    {
//        //  Print Geometrical Base Mesh
//        std::ofstream argument("GeometicMeshNew4.txt");
//        GeoMesh4->Print(argument);
//        std::ofstream Dummyfile("GeometricMeshNew4.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh4,Dummyfile, true);
//    }
//    return GeoMesh4;
//}
//
//TPZGeoMesh *GMesh(int d, bool ftriang, int ndiv)
//{
//    
//    if(dim!=2)
//    {
//        std::cout << "dimensao errada" << std::endl;
//        dim = 2;
//        DebugStop();
//    }
//    
//    int Qnodes =  4;
//    
//    TPZGeoMesh * gmesh = new TPZGeoMesh;
//    gmesh->SetMaxNodeId(Qnodes-1);
//    gmesh->NodeVec().Resize(Qnodes);
//    TPZVec<TPZGeoNode> Node(Qnodes);
//    
//    gmesh->SetDimension(dim);
//    
//    TPZVec <long> TopolQuad(4);
//    TPZVec <long> TopolTriang(3);
//    TPZVec <long> TopolLine(2);
//    TPZVec <long> TopolPoint(1);
//    
//    //indice dos nos
//    long id = 0;
//    //    REAL valx;
//    //    for(int xi = 0; xi < Qnodes/2; xi++)
//    //    {
//    //        valx = xi*Lx;
//    //        Node[id].SetNodeId(id);
//    //        Node[id].SetCoord(0 ,valx );//coord X
//    //        Node[id].SetCoord(1 ,0. );//coord Y
//    //        gmesh->NodeVec()[id] = Node[id];
//    //        id++;
//    //    }
//    //
//    //    for(int xi = 0; xi < Qnodes/2; xi++)
//    //    {
//    //        valx = Lx - xi*Lx;
//    //        Node[id].SetNodeId(id);
//    //        Node[id].SetCoord(0 ,valx );//coord X
//    //        Node[id].SetCoord(1 ,Ly);//coord Y
//    //        gmesh->NodeVec()[id] = Node[id];
//    //        id++;
//    //    }
//    //
//    TPZManVector<REAL,3> coord(2,0.);
//    int in = 0;
//    //c0
//    coord[0] = -1.0;//0.0; //
//    coord[1] = -1.0;//0.0; //
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c1
//    coord[0] =   1.0;//1.0;//
//    coord[1] =  -1.0;//0.0;//
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c2
//    coord[0] =   1.0;//0.0;//
//    coord[1] =   1.0;//1.0;//
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c3
//    coord[0] =  -1.0;//1.0;//
//    coord[1] =   1.0;//1.0;//
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    //indice dos elementos
//    id = 0;
//    
//    if(ftriang==true) // triangulo
//    {
//        TopolTriang[0] = 0;
//        TopolTriang[1] = 1;
//        TopolTriang[2] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
//        id++;
//        
//        TopolTriang[0] = 0;
//        TopolTriang[1] = 2;
//        TopolTriang[2] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
//        id++;
//        
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//        id++;
//        
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//        id++;
//        
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//        id++;
//        
//        TopolLine[0] = 3;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
//        id++;
//    }
//    else{
//        
//        TopolQuad[0] = 0;
//        TopolQuad[1] = 1;
//        TopolQuad[2] = 2;
//        TopolQuad[3] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//        id++;
//        
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//        id++;
//        
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//        id++;
//        
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//        id++;
//        
//        TopolLine[0] = 3;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
//    }
//    
//    gmesh->BuildConnectivity();
//    
//    /// gmesh para aqui
//    
//    TPZVec<TPZGeoEl *> sons;
//    for (int iref = 0; iref < ndiv; iref++) {
//        int nel = gmesh->NElements();
//        for (int iel = 0; iel < nel; iel++) {
//            TPZGeoEl *gel = gmesh->ElementVec()[iel];
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            gel->Divide(sons);
//        }
//    }
//    
//    
//    //#ifdef LOG4CXX
//    //	if(logdata->isDebugEnabled())
//    //	{
//    //        std::stringstream sout;
//    //        sout<<"\n\n Malha Geometrica Inicial\n ";
//    //        gmesh->Print(sout);
//    //        LOGPZ_DEBUG(logdata,sout.str())
//    //	}
//    //#endif
//    
//    return gmesh;
//}
//
//TPZGeoMesh *CreateOneCuboWithTetraedrons(long nelem, int MaterialId)
//{
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    GenerateNodes(gmesh,nelem);
//    
//    for (long i=0; i<nelem; i++) {
//        for (long j=0; j<nelem; j++) {
//            for (long k=0; k<nelem; k++) {
//                TPZManVector<long,8> nodes(8,0);
//                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
//                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
//                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
//                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
//                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
//                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
//                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
//                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
//#ifdef LOG4CXX
//                if(logdata->isDebugEnabled())
//                {
//                    std::stringstream sout;
//                    sout << "Tetrahedral nodes " << nodes;
//                    LOGPZ_DEBUG(logdata, sout.str())
//                }
//#endif
//                for (int el=0; el<6; el++)
//                {
//                    TPZManVector<long,4> elnodes(4);
//                    long index;
//                    for (int il=0; il<4; il++) {
//                        elnodes[il] = nodes[tetraedra_2[el][il]];
//                    }
//                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index);
//                }
//            }
//        }
//    }
//    gmesh->BuildConnectivity();
//    
//    // Boundary Conditions
//    const int numelements = gmesh->NElements();
//    //    const int bczMinus = -3, bczplus = -2, bcids = -1;
//    //    const int bczMinus = -1, bczplus = -1, bcids = -1;
//    
//    for(int el=0; el<numelements; el++)
//    {
//        TPZManVector <TPZGeoNode,4> Nodefinder(4);
//        TPZManVector <REAL,3> nodecoord(3);
//        TPZGeoEl *tetra = gmesh->ElementVec()[el];
//        TPZVec<long> ncoordVec(0); long sizeOfVec = 0;
//        
//        // na face z = 0
//        for (int i = 0; i < 4; i++)
//        {
//            long pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[2],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc0);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face y = 0
//        for (int i = 0; i < 4; i++)
//        {
//            long pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[1],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc1);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face x = 1
//        for (int i = 0; i < 4; i++)
//        {
//            long pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[0],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc2);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face y = 1
//        for (int i = 0; i < 4; i++)
//        {
//            long pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[1],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc3);
//        }
//        
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face x = 0
//        for (int i = 0; i < 4; i++)
//        {
//            long pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[0],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc4);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face z = 1
//        for (int i = 0; i < 4; i++)
//        {
//            long pos = tetra->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[2],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 3)
//        {
//            int lado = tetra->WhichSide(ncoordVec);
//            TPZGeoElSide tetraSide(tetra, lado);
//            TPZGeoElBC(tetraSide,bc5);
//        }
//        
//        
//        
//    }
//    
//    return gmesh;
//}
//
//bool MyDoubleComparer(REAL a, REAL b)
//{
//    if (IsZero(a-b)){
//        return true;
//    }
//    else{
//        return false;
//    }
//}
//
//void GenerateNodes(TPZGeoMesh *gmesh, long nelem)
//{
//    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
//    for (long i=0; i<=nelem; i++) {
//        for (long j=0; j<=nelem; j++) {
//            for (long k=0; k<=nelem; k++) {
//                TPZManVector<REAL,3> x(3);
//                x[0] = k*1./nelem;
//                x[1] = j*1./nelem;
//                x[2] = i*1./nelem;
//                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
//            }
//        }
//    }
//}
//
//void GenerateNodesPyramid(TPZGeoMesh *gmesh, long nelem)
//{
//    long sizenodevec = (nelem+1)*(nelem+1)*(nelem+1)+(nelem*nelem*nelem);
//    gmesh->NodeVec().Resize(sizenodevec);
//    int posicao = 0;
//    for (long i=0; i<=nelem; i++) {
//        for (long j=0; j<=nelem; j++) {
//            for (long k=0; k<=nelem; k++) {
//                TPZManVector<REAL,3> x(3);
//                x[0] = k*1./nelem;
//                x[1] = j*1./nelem;
//                x[2] = i*1./nelem;
//                posicao = i*(nelem+1)*(nelem+1)+j*(nelem+1)+k;
//                gmesh->NodeVec()[posicao].Initialize(x, *gmesh);
//                 
//            }
//        }
//    }
//    for (long i=0; i<nelem; i++) {
//        for (long j=0; j<nelem; j++) {
//            for (long k=0; k<nelem; k++) {
//                TPZManVector<REAL,3> xc(3);
//                xc[0] = k*1./nelem + (0.5)*1./nelem;
//                xc[1] = j*1./nelem + (0.5)*1./nelem;
//                xc[2] = i*1./nelem + (0.5)*1./nelem;
//                posicao++;
//                gmesh->NodeVec()[posicao].Initialize(xc, *gmesh);
//            }
//        }
//    }
//}
//
//
//TPZGeoMesh *CreateOneCubo(int nref)
//{
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    int nnodes = 8;
//    int idf0=-1;
//    int idf1=-2;
//    int idf2=-3;
//    int idf3=-4;
//    int idf4=-5;
//    int idf5=-6;
//    
//    gmesh->SetDimension(3);
//    gmesh->NodeVec().Resize(nnodes);
//    
//    TPZManVector<REAL,3> coord(3,0.);
//    int in = 0;
//    
//    //cubo [0,1]ˆ3
//    //    //c0
//    //    coord[0] = 0.0;
//    //    coord[1] = 0.0;
//    //    coord[2] = 0.0;
//    //    gmesh->NodeVec()[in].SetCoord(coord);
//    //    gmesh->NodeVec()[in].SetNodeId(in++);
//    //    //c1
//    //    coord[0] =  1.0;
//    //    coord[1] = 0.0;
//    //    coord[2] = 0.0;
//    //    gmesh->NodeVec()[in].SetCoord(coord);
//    //    gmesh->NodeVec()[in].SetNodeId(in++);
//    //    //c2
//    //    coord[0] =  1.0;
//    //    coord[1] =  1.0;
//    //    coord[2] = 0.0;
//    //    gmesh->NodeVec()[in].SetCoord(coord);
//    //    gmesh->NodeVec()[in].SetNodeId(in++);
//    //    //c3
//    //    coord[0] = 0.0;
//    //    coord[1] =  1.0;
//    //    coord[2] = 0.0;
//    //    gmesh->NodeVec()[in].SetCoord(coord);
//    //    gmesh->NodeVec()[in].SetNodeId(in++);
//    //
//    //    //c4
//    //    coord[0] = 0.0;
//    //    coord[1] = 0.0;
//    //    coord[2] =  1.0;
//    //    gmesh->NodeVec()[in].SetCoord(coord);
//    //    gmesh->NodeVec()[in].SetNodeId(in++);
//    //    //c5
//    //    coord[0] =  1.0;
//    //    coord[1] = 0.0;
//    //    coord[2] =  1.0;
//    //    gmesh->NodeVec()[in].SetCoord(coord);
//    //    gmesh->NodeVec()[in].SetNodeId(in++);
//    //    //c6
//    //    coord[0] =  1.0;
//    //    coord[1] =  1.0;
//    //    coord[2] =  1.0;
//    //    gmesh->NodeVec()[in].SetCoord(coord);
//    //    gmesh->NodeVec()[in].SetNodeId(in++);
//    //    //c7
//    //    coord[0] = 0.0;
//    //    coord[1] =  1.0;
//    //    coord[2] =  1.0;
//    //    gmesh->NodeVec()[in].SetCoord(coord);
//    //    gmesh->NodeVec()[in].SetNodeId(in++);
//    
//    // cubo [-1,1]^3
//    //c0
//    coord[0] = -1.0;
//    coord[1] = -1.0;
//    coord[2] = -1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c1
//    coord[0] =  1.0;
//    coord[1] = -1.0;
//    coord[2] = -1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c2
//    coord[0] =  1.0;
//    coord[1] =  1.0;
//    coord[2] = -1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c3
//    coord[0] = -1.0;
//    coord[1] =  1.0;
//    coord[2] = -1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    
//    //c4
//    coord[0] = -1.0;
//    coord[1] = -1.0;
//    coord[2] =  1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c5
//    coord[0] =  1.0;
//    coord[1] = -1.0;
//    coord[2] =  1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c6
//    coord[0] =  1.0;
//    coord[1] =  1.0;
//    coord[2] =  1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    //c7
//    coord[0] = -1.0;
//    coord[1] =  1.0;
//    coord[2] =  1.0;
//    gmesh->NodeVec()[in].SetCoord(coord);
//    gmesh->NodeVec()[in].SetNodeId(in);
//    in++;
//    
//    
//    int index = 0;
//    
//    TPZVec<long> TopologyQuad(4);
//    
//    // bottom
//    TopologyQuad[0]=0;
//    TopologyQuad[1]=1;
//    TopologyQuad[2]=2;
//    TopologyQuad[3]=3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf0,*gmesh);
//    index++;
//    
//    // Front
//    TopologyQuad[0]=0;
//    TopologyQuad[1]=1;
//    TopologyQuad[2]=5;
//    TopologyQuad[3]=4;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf1,*gmesh);
//    index++;
//    
//    // Rigth
//    TopologyQuad[0]=1;
//    TopologyQuad[1]=2;
//    TopologyQuad[2]=6;
//    TopologyQuad[3]=5;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf2,*gmesh);
//    index++;
//    // Back
//    TopologyQuad[0]=3;
//    TopologyQuad[1]=2;
//    TopologyQuad[2]=6;
//    TopologyQuad[3]=7;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf3,*gmesh);
//    index++;
//    
//    // Left
//    TopologyQuad[0]=0;
//    TopologyQuad[1]=3;
//    TopologyQuad[2]=7;
//    TopologyQuad[3]=4;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf4,*gmesh);
//    index++;
//    
//    // Top
//    TopologyQuad[0]=4;
//    TopologyQuad[1]=5;
//    TopologyQuad[2]=6;
//    TopologyQuad[3]=7;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf5,*gmesh);
//    index++;
//    
//    TPZManVector<long,8> TopolCubo(8,0);
//    TopolCubo[0] = 0;
//    TopolCubo[1] = 1;
//    TopolCubo[2] = 2;
//    TopolCubo[3] = 3;
//    TopolCubo[4] = 4;
//    TopolCubo[5] = 5;
//    TopolCubo[6] = 6;
//    TopolCubo[7] = 7;
//    
//    
//    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index, TopolCubo, matId, *gmesh);
//    index++;
//    
//    gmesh->BuildConnectivity();
//    
//    /// gmesh para aqui
//    
//    TPZVec<TPZGeoEl *> sons;
//    for (int iref = 0; iref < nref; iref++) {
//        int nel = gmesh->NElements();
//        for (int iel = 0; iel < nel; iel++) {
//            TPZGeoEl *gel = gmesh->ElementVec()[iel];
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            gel->Divide(sons);
//        }
//    }
//    
//    
//    std::ofstream out("SingleCubeWithBcs.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
//    
//    return gmesh;
//}
//
//TPZGeoMesh *GMeshCubeWithPyramids(long nelem, int MaterialId)
//{
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    GenerateNodesPyramid(gmesh,nelem);
//    
//    int posicao = 0;
//    for (long i=0; i<=nelem; i++) {
//        for (long j=0; j<=nelem; j++) {
//            for (long k=0; k<=nelem; k++) {
//                posicao = i*(nelem+1)*(nelem+1)+j*(nelem+1)+k;
//            }
//        }
//    }
//
//    
//    for (long i=0; i<nelem; i++) {
//        for (long j=0; j<nelem; j++) {
//            for (long k=0; k<nelem; k++) {
//                TPZManVector<long,9> nodes(9,0);
//                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
//                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
//                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
//                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
//                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
//                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
//                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
//                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
//                nodes[8] = posicao + (k)*(nelem)*(nelem)+(j)*(nelem)+i+1;
//#ifdef LOG4CXX
//                if(logdata->isDebugEnabled())
//                {
//                    std::stringstream sout;
//                    sout << "Pyramidal nodes " << nodes;
//                    LOGPZ_DEBUG(logdata, sout.str())
//                }
//#endif
//                for (int el=0; el<6; el++)
//                {
//                    TPZManVector<long,5> elnodes(5);
//                    long index;
//                    for (int il=0; il<5; il++) {
//                        elnodes[il] = nodes[piramide_2[el][il]];
//                    }
//                    gmesh->CreateGeoElement(EPiramide, elnodes, MaterialId, index);
//                }
//            }
//        }
//    }
//    gmesh->BuildConnectivity();
//    
//    // Boundary Conditions
//    const int numelements = gmesh->NElements();
//    //    const int bczMinus = -3, bczplus = -2, bcids = -1;
//    //    const int bczMinus = -1, bczplus = -1, bcids = -1;
//    
//    for(int el=0; el<numelements; el++)
//    {
//        TPZManVector <TPZGeoNode,5> Nodefinder(5);
//        TPZManVector <REAL,3> nodecoord(3);
//        TPZGeoEl *piramide = gmesh->ElementVec()[el];
//        TPZVec<long> ncoordVec(0); long sizeOfVec = 0;
//        
//        // na face z = 0
//        for (int i = 0; i < 5; i++)
//        {
//            long pos = piramide->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[2],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 4)
//        {
//            int lado = piramide->WhichSide(ncoordVec);
//            TPZGeoElSide piramideSide(piramide, lado);
//            TPZGeoElBC(piramideSide,bc0);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face y = 0
//        for (int i = 0; i < 5; i++)
//        {
//            long pos = piramide->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[1],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 4)
//        {
//            int lado = piramide->WhichSide(ncoordVec);
//            TPZGeoElSide piramideSide(piramide, lado);
//            TPZGeoElBC(piramideSide,bc1);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face x = 1
//        for (int i = 0; i < 5; i++)
//        {
//            long pos = piramide->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[0],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 4)
//        {
//            int lado = piramide->WhichSide(ncoordVec);
//            TPZGeoElSide piramideSide(piramide, lado);
//            TPZGeoElBC(piramideSide,bc2);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face y = 1
//        for (int i = 0; i < 5; i++)
//        {
//            long pos = piramide->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[1],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 4)
//        {
//            int lado = piramide->WhichSide(ncoordVec);
//            TPZGeoElSide piramideSide(piramide, lado);
//            TPZGeoElBC(piramideSide,bc3);
//        }
//        
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face x = 0
//        for (int i = 0; i < 5; i++)
//        {
//            long pos = piramide->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[0],0.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 4)
//        {
//            int lado = piramide->WhichSide(ncoordVec);
//            TPZGeoElSide piramideSide(piramide, lado);
//            TPZGeoElBC(piramideSide,bc4);
//        }
//        
//        ncoordVec.clear();
//        sizeOfVec = 0;
//        // na face z = 1
//        for (int i = 0; i < 5; i++)
//        {
//            long pos = piramide->NodeIndex(i);
//            Nodefinder[i] = gmesh->NodeVec()[pos];
//            Nodefinder[i].GetCoordinates(nodecoord);
//            if (MyDoubleComparer(nodecoord[2],1.))
//            {
//                sizeOfVec++;
//                ncoordVec.Resize(sizeOfVec);
//                ncoordVec[sizeOfVec-1] = pos;
//            }
//        }
//        if(ncoordVec.NElements() == 4)
//        {
//            int lado = piramide->WhichSide(ncoordVec);
//            TPZGeoElSide piramideSide(piramide, lado);
//            TPZGeoElBC(piramideSide,bc5);
//        }
//        
//        
//        
//    }
//    
//    
//    std::ofstream out("CubeWithBcs.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
//    
//    return gmesh;
//}
//
//
//void Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X)
//{
//    X[0] = par[0];
//    X[1] = 0.0;
//    X[2] = 0.0;
//}
//
//void Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X)
//{
//    X[0] = 0.0;
//    X[1] = par[0];
//    X[2] = 0.0;
//}
//
//void Parametricfunction3(const TPZVec<STATE> &par, TPZVec<STATE> &X)
//{
//    X[0] = 0.0;
//    X[1] = 0.0;
//    X[2] = par[0];
//}


