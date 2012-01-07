#include <iostream>
#include "pzfmatrix.h"
#include "pzelasmat.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzlog.h"
#include "tpzgeoblend.h"
#include "tpzellipse3d.h"
#include "tpzarc3d.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeopoint.h"
#include "pzgeopyramid.h"
#include "pzgeoquad.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "pzgmesh.h"
#include "tpzmathtools.h"

#include "TPZPlaneFracture.h"
#include "TPZPoligonalChain.h"

#include "pzgengrid.h"

using namespace std;

const int matPoint = -3;

//** just for visualize given dots in vtk */
void InsertDots4VTK(TPZGeoMesh * gmesh, TPZVec<REAL> &fractureDots);

void FillFractureDotsExampleEllipse(TPZVec<REAL> &fractureDots);
void FillFractureDotsExampleCrazy(TPZVec<REAL> &fractureDots);

//----------------------------------------------------------------------------------------------------------------------------------

#include "TPZTimer.h"

int main(int argc, char * const argv[])
{	
    TPZTimer readRef("ReadingRefPatterns");
    readRef.start();
    
    //#define writeAgain
    #ifdef writeAgain
        gRefDBase.InitializeRefPatterns();
        
        std::ofstream outRefP("RefPatternsUsed.txt");
        gRefDBase.WriteRefPatternDBase(outRefP);
    #else
        std::ifstream inRefP("RefPatternsUsed.txt");
        gRefDBase.ReadRefPatternDBase("RefPatternsUsed.txt");
    #endif
    readRef.stop();
    std::cout << "DeltaT leitura refpatterns = " << readRef.seconds() << " s" << std::endl;
    
    double lw = 200.;
    double bulletDepthIni = 20.;
    double bulletDepthFin = 180.;
    TPZVec< std::map<double,double> > pos_stress(5);
    pos_stress[0][12.]  = 0.;
    pos_stress[1][24.]  = 0.;
    pos_stress[2][36.]  = 0.;
    pos_stress[3][64.]  = 0.;
    pos_stress[4][176.] = 0.;
    
    TPZPlaneFracture plfrac(lw, bulletDepthIni, bulletDepthFin, pos_stress);
    
    TPZVec<REAL> fractureDots;
    
    
    
    
//    FillFractureDotsExampleCrazy(fractureDots);
//    
//    TPZTimer clockIni1("PartyBegins1");
//    clockIni1.start();
//
//    TPZGeoMesh * fractureMesh = plfrac.GetFractureMesh(fractureDots);
//    clockIni1.stop();
//    std::cout << "DeltaT get fracture mesh = " << clockIni1.seconds() << " s" << std::endl;
//
//    InsertDots4VTK(fractureMesh, fractureDots);
//
//    std::ofstream out1("FractureZprofile1.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(fractureMesh, out1, true);
    
    
    
    
    FillFractureDotsExampleEllipse(fractureDots);
    
    TPZTimer clockIni2("PartyBegins2");
    clockIni2.start();
    
    TPZGeoMesh * fractureMesh2 = plfrac.GetFractureMesh(fractureDots);
    clockIni2.stop();
    std::cout << "DeltaT get fracture mesh = " << clockIni2.seconds() << " s" << std::endl;
    
    InsertDots4VTK(fractureMesh2, fractureDots);
    
    std::ofstream out2("FractureZprofile2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fractureMesh2, out2, true);
    
    
    #ifndef writeAgain
    std::ofstream outRefP("RefPatternsUsed.txt");
    gRefDBase.WriteRefPatternDBase(outRefP);
    #endif
    
    return 0;
}

//** just for visualize given dots in vtk */
void InsertDots4VTK(TPZGeoMesh * gmesh, TPZVec<REAL> &fractureDots)
{
    int nDots = fractureDots.size() / 2;
    int nnodesOriginal = gmesh->NNodes();
	int Qnodes = nnodesOriginal + nDots;
    	
	//initializing gmesh->NodeVec()
	gmesh->NodeVec().Resize(Qnodes);
	TPZGeoNode Node;
    TPZVec<REAL> NodeCoord(3);
    TPZVec<int> Topol(1);
    
    int elId = gmesh->NElements();
	for(int n = nnodesOriginal; n < Qnodes; n++)
	{
        Topol[0] = n;

        int actDot = n - nnodesOriginal;
        
        NodeCoord[0] = fractureDots[2*actDot];
        NodeCoord[1] = 0.;
        NodeCoord[2] = fractureDots[2*actDot + 1];
        
		Node.SetNodeId(n);
		Node.SetCoord(NodeCoord);
		gmesh->NodeVec()[n] = Node; 
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol, matPoint,*gmesh);
        elId++;
	}
}


void FillFractureDotsExampleEllipse(TPZVec<REAL> &fractureDots)
{
    int nnodes = 80;
    
    fractureDots.Resize(2*nnodes, 0.);
    int node;
    
    node = 0;
    
    fractureDots[2*node] = 0.5; fractureDots[2*node+1] = 95.;
    
    node = 1;
    
    fractureDots[2*node] = 5.; fractureDots[2*node+1] = 94.9852;
    
    node = 2;
    
    fractureDots[2*node] = 10.; fractureDots[2*node+1] = 94.9408;
    
    node = 3;
    
    fractureDots[2*node] = 15.; fractureDots[2*node+1] = 94.8667;
    
    node = 4;
    
    fractureDots[2*node] = 20.; fractureDots[2*node+1] = 94.7627;
    
    node = 5;
    
    fractureDots[2*node] = 25.; fractureDots[2*node+1] = 94.6286;
    
    node = 6;
    
    fractureDots[2*node] = 30.; fractureDots[2*node+1] = 94.4643;
    
    node = 7;
    
    fractureDots[2*node] = 35.; fractureDots[2*node+1] = 94.2692;
    
    node = 8;
    
    fractureDots[2*node] = 40.; fractureDots[2*node+1] = 94.0431;
    
    node = 9;
    
    fractureDots[2*node] = 45.; fractureDots[2*node+1] = 93.7854;
    
    node = 10;
    
    fractureDots[2*node] = 50.; fractureDots[2*node+1] = 93.4956;
    
    node = 11;
    
    fractureDots[2*node] = 55.; fractureDots[2*node+1] = 93.173;
    
    node = 12;
    
    fractureDots[2*node] = 60.; fractureDots[2*node+1] = 92.8169;
    
    node = 13;
    
    fractureDots[2*node] = 65.; fractureDots[2*node+1] = 92.4264;
    
    node = 14;
    
    fractureDots[2*node] = 70.; fractureDots[2*node+1] = 92.0006;
    
    node = 15;
    
    fractureDots[2*node] = 75.; fractureDots[2*node+1] = 91.5385;
    
    node = 16;
    
    fractureDots[2*node] = 80.; fractureDots[2*node+1] = 91.0387;
    
    node = 17;
    
    fractureDots[2*node] = 85.; fractureDots[2*node+1] = 90.4998;
    
    node = 18;
    
    fractureDots[2*node] = 90.; fractureDots[2*node+1] = 89.9204;
    
    node = 19;
    
    fractureDots[2*node] = 95.; fractureDots[2*node+1] = 89.2986;
    
    node = 20;
    
    fractureDots[2*node] = 100.; fractureDots[2*node+1] = 88.6323;
    
    node = 21;
    
    fractureDots[2*node] = 105.; fractureDots[2*node+1] = 87.9193;
    
    node = 22;
    
    fractureDots[2*node] = 110.; fractureDots[2*node+1] = 87.1567;
    
    node = 23;
    
    fractureDots[2*node] = 115.; fractureDots[2*node+1] = 86.3416;
    
    node = 24;
    
    fractureDots[2*node] = 120.; fractureDots[2*node+1] = 85.4702;
    
    node = 25;
    
    fractureDots[2*node] = 125.; fractureDots[2*node+1] = 84.5384;
    
    node = 26;
    
    fractureDots[2*node] = 130.; fractureDots[2*node+1] = 83.541;
    
    node = 27;
    
    fractureDots[2*node] = 135.; fractureDots[2*node+1] = 82.4721;
    
    node = 28;
    
    fractureDots[2*node] = 140.; fractureDots[2*node+1] = 81.3243;
    
    node = 29;
    
    fractureDots[2*node] = 145.; fractureDots[2*node+1] = 80.0886;
    
    node = 30;
    
    fractureDots[2*node] = 150.; fractureDots[2*node+1] = 78.7537;
    
    node = 31;
    
    fractureDots[2*node] = 155.; fractureDots[2*node+1] = 77.305;
    
    node = 32;
    
    fractureDots[2*node] = 160.; fractureDots[2*node+1] = 75.7233;
    
    node = 33;
    
    fractureDots[2*node] = 165.; fractureDots[2*node+1] = 73.9822;
    
    node = 34;
    
    fractureDots[2*node] = 170.; fractureDots[2*node+1] = 72.0442;
    
    node = 35;
    
    fractureDots[2*node] = 175.; fractureDots[2*node+1] = 69.8515;
    
    node = 36;
    
    fractureDots[2*node] = 180.; fractureDots[2*node+1] = 67.3077;
    
    node = 37;
    
    fractureDots[2*node] = 185.; fractureDots[2*node+1] = 64.2256;
    
    node = 38;
    
    fractureDots[2*node] = 190.; fractureDots[2*node+1] = 60.125;
    
    node = 39;
    
    fractureDots[2*node] = 195.; fractureDots[2*node+1] = 50.;
    
    node = 40;
    
    fractureDots[2*node] = 195.; fractureDots[2*node+1] = 42.;
    
    node = 41;
    
    fractureDots[2*node] = 190.; fractureDots[2*node+1] = 39.875;
    
    node = 42;
    
    fractureDots[2*node] = 185.; fractureDots[2*node+1] = 35.7744;
    
    node = 43;
    
    fractureDots[2*node] = 180.; fractureDots[2*node+1] = 32.6923;
    
    node = 44;
    
    fractureDots[2*node] = 175.; fractureDots[2*node+1] = 30.1485;
    
    node = 45;
    
    fractureDots[2*node] = 170.; fractureDots[2*node+1] = 27.9558;
    
    node = 46;
    
    fractureDots[2*node] = 165.; fractureDots[2*node+1] = 26.0178;
    
    node = 47;
    
    fractureDots[2*node] = 160.; fractureDots[2*node+1] = 24.2767;
    
    node = 48;
    
    fractureDots[2*node] = 155.; fractureDots[2*node+1] = 22.695;
    
    node = 49;
    
    fractureDots[2*node] = 150.; fractureDots[2*node+1] = 21.2463;
    
    node = 50;
    
    fractureDots[2*node] = 145.; fractureDots[2*node+1] = 19.9114;
    
    node = 51;
    
    fractureDots[2*node] = 140.; fractureDots[2*node+1] = 18.6757;
    
    node = 52;
    
    fractureDots[2*node] = 135.; fractureDots[2*node+1] = 17.5279;
    
    node = 53;
    
    fractureDots[2*node] = 130.; fractureDots[2*node+1] = 16.459;
    
    node = 54;
    
    fractureDots[2*node] = 125.; fractureDots[2*node+1] = 15.4616;
    
    node = 55;
    
    fractureDots[2*node] = 120.; fractureDots[2*node+1] = 14.5298;
    
    node = 56;
    
    fractureDots[2*node] = 115.; fractureDots[2*node+1] = 13.6584;
    
    node = 57;
    
    fractureDots[2*node] = 110.; fractureDots[2*node+1] = 12.8433;
    
    node = 58;
    
    fractureDots[2*node] = 105.; fractureDots[2*node+1] = 12.0807;
    
    node = 59;
    
    fractureDots[2*node] = 100.; fractureDots[2*node+1] = 11.3677;
    
    node = 60;
    
    fractureDots[2*node] = 95.; fractureDots[2*node+1] = 10.7014;
    
    node = 61;
    
    fractureDots[2*node] = 90.; fractureDots[2*node+1] = 10.0796;
    
    node = 62;
    
    fractureDots[2*node] = 85.; fractureDots[2*node+1] = 9.50016;
    
    node = 63;
    
    fractureDots[2*node] = 80.; fractureDots[2*node+1] = 8.96134;
    
    node = 64;
    
    fractureDots[2*node] = 75.; fractureDots[2*node+1] = 8.46154;
    
    node = 65;
    
    fractureDots[2*node] = 70.; fractureDots[2*node+1] = 7.99937;
    
    node = 66;
    
    fractureDots[2*node] = 65.; fractureDots[2*node+1] = 7.57359;
    
    node = 67;
    
    fractureDots[2*node] = 60.; fractureDots[2*node+1] = 7.18313;
    
    node = 68;
    
    fractureDots[2*node] = 55.; fractureDots[2*node+1] = 6.82703;
    
    node = 69;
    
    fractureDots[2*node] = 50.; fractureDots[2*node+1] = 6.50444;
    
    node = 70;
    
    fractureDots[2*node] = 45.; fractureDots[2*node+1] = 6.21462;
    
    node = 71;
    
    fractureDots[2*node] = 40.; fractureDots[2*node+1] = 5.95692;
    
    node = 72;
    
    fractureDots[2*node] = 35.; fractureDots[2*node+1] = 5.73079;
    
    node = 73;
    
    fractureDots[2*node] = 30.; fractureDots[2*node+1] = 5.53573;
    
    node = 74;
    
    fractureDots[2*node] = 25.; fractureDots[2*node+1] = 5.37135;
    
    node = 75;
    
    fractureDots[2*node] = 20.; fractureDots[2*node+1] = 5.23731;
    
    node = 76;
    
    fractureDots[2*node] = 15.; fractureDots[2*node+1] = 5.13333;
    
    node = 77;
    
    fractureDots[2*node] = 10.; fractureDots[2*node+1] = 5.05921;
    
    node = 78;
    
    fractureDots[2*node] = 5.; fractureDots[2*node+1] = 5.0148;
    
    node = 79;
    
    fractureDots[2*node] = 0.5; fractureDots[2*node+1] = 5.;
}


void FillFractureDotsExampleCrazy(TPZVec<REAL> &fractureDots)
{
    int nnodes = 62;
    
    fractureDots.Resize(2*nnodes, 0.);
    int node;
    
    ///
    node = 0;
    double desloc = 0.;
    
    fractureDots[2*node] = 5.; fractureDots[2*node+1] = 80.;
    
    node = 1;
    
    fractureDots[2*node] = 9.02368; fractureDots[2*node+1] = 82.1856;
    
    node = 2;
    
    fractureDots[2*node] = 25.6151; fractureDots[2*node+1] = 82.9122;
    
    node = 3;
    
    fractureDots[2*node] = 48.4505; fractureDots[2*node+1] = 82.5833;
    
    node = 4;
    
    fractureDots[2*node] = 74.1186; fractureDots[2*node+1] = 81.5259;
    
    node = 5;
    
    fractureDots[2*node] = 100. + desloc; fractureDots[2*node+1] = 80.;//<----------------------
    
    node = 6;
    
    fractureDots[2*node] = 124.157; fractureDots[2*node+1] = 78.2079;
    
    node = 7;
    
    fractureDots[2*node] = 145.23; fractureDots[2*node+1] = 76.3027;
    
    node = 8;
    
    fractureDots[2*node] = 162.349; fractureDots[2*node+1] = 74.3954;
    
    node = 9;
    
    fractureDots[2*node] = 175.043; fractureDots[2*node+1] = 72.5625;
    
    node = 10;
    
    fractureDots[2*node] = 183.172; fractureDots[2*node+1] = 70.8516;
    
    node = 11;
    
    fractureDots[2*node] = 186.85; fractureDots[2*node+1] = 69.2873;
    
    node = 12;
    
    fractureDots[2*node] = 186.391; fractureDots[2*node+1] = 67.8763;
    
    node = 13;
    
    fractureDots[2*node] = 182.253; fractureDots[2*node+1] = 66.6118;
    
    node = 14;
    
    fractureDots[2*node] = 174.99; fractureDots[2*node+1] = 65.4771;
    
    node = 15;
    
    fractureDots[2*node] = 165.212; fractureDots[2*node+1] = 64.4495;
    
    node = 16;
    
    fractureDots[2*node] = 153.552; fractureDots[2*node+1] = 63.5025;
    
    node = 17;
    
    fractureDots[2*node] = 140.633; fractureDots[2*node+1] = 62.609;
    
    node = 18;
    
    fractureDots[2*node] = 127.05; fractureDots[2*node+1] = 61.7426;
    
    node = 19;
    
    fractureDots[2*node] = 113.346; fractureDots[2*node+1] = 60.8796;
    
    node = 20;
    
    fractureDots[2*node] = 100. + desloc; fractureDots[2*node+1] = 60.;//<----------------------
    
    node = 21;
    
    fractureDots[2*node] = 87.418; fractureDots[2*node+1] = 59.0884;
    
    node = 22;
    
    fractureDots[2*node] = 75.9253; fractureDots[2*node+1] = 58.1348;
    
    node = 23;
    
    fractureDots[2*node] = 65.7644; fractureDots[2*node+1] = 57.1342;
    
    node = 24;
    
    fractureDots[2*node] = 57.0956; fractureDots[2*node+1] = 56.0873;
    
    node = 25;
    
    fractureDots[2*node] = 50.; fractureDots[2*node+1] = 55.;
    
    node = 26;
    
    fractureDots[2*node] = 44.4857; fractureDots[2*node+1] = 53.8827;
    
    node = 27;
    
    fractureDots[2*node] = 40.495; fractureDots[2*node+1] = 52.75;
    
    node = 28;
    
    fractureDots[2*node] = 37.9143; fractureDots[2*node+1] = 51.6198;
    
    node = 29;
    
    fractureDots[2*node] = 36.5847; fractureDots[2*node+1] = 50.5124;
    
    node = 30;
    
    fractureDots[2*node] = 36.3137; fractureDots[2*node+1] = 49.4496;
    
    node = 31;
    
    fractureDots[2*node] = 36.888; fractureDots[2*node+1] = 48.4535;
    
    node = 32;
    
    fractureDots[2*node] = 38.0857; fractureDots[2*node+1] = 47.5455;
    
    node = 33;
    
    fractureDots[2*node] = 39.6892; fractureDots[2*node+1] = 46.7455;
    
    node = 34;
    
    fractureDots[2*node] = 41.4972; fractureDots[2*node+1] = 46.0702;
    
    node = 35;
    
    fractureDots[2*node] = 43.3363; fractureDots[2*node+1] = 45.5328;
    
    node = 36;
    
    fractureDots[2*node] = 45.0707; fractureDots[2*node+1] = 45.1417;
    
    node = 37;
    
    fractureDots[2*node] = 46.6112; fractureDots[2*node+1] = 44.8996;
    
    node = 38;
    
    fractureDots[2*node] = 47.9218; fractureDots[2*node+1] = 44.8032;
    
    node = 39;
    
    fractureDots[2*node] = 49.0243; fractureDots[2*node+1] = 44.8424;
    
    node = 40;
    
    fractureDots[2*node] = 50.; fractureDots[2*node+1] = 45.;
    
    node = 41;
    
    fractureDots[2*node] = 50.9889; fractureDots[2*node+1] = 45.2519;
    
    node = 42;
    
    fractureDots[2*node] = 52.1852; fractureDots[2*node+1] = 45.5668;
    
    node = 43;
    
    fractureDots[2*node] = 53.8292; fractureDots[2*node+1] = 45.9072;
    
    node = 44;
    
    fractureDots[2*node] = 56.1951; fractureDots[2*node+1] = 46.2297;
    
    node = 45;
    
    fractureDots[2*node] = 59.5747; fractureDots[2*node+1] = 46.4863;
    
    node = 46;
    
    fractureDots[2*node] = 64.2563; fractureDots[2*node+1] = 46.626;
    
    node = 47;
    
    fractureDots[2*node] = 70.4977; fractureDots[2*node+1] = 46.5964;
    
    node = 48;
    
    fractureDots[2*node] = 78.4951; fractureDots[2*node+1] = 46.3461;
    
    node = 49;
    
    fractureDots[2*node] = 88.3449; fractureDots[2*node+1] = 45.8275;
    
    node = 50;
    
    fractureDots[2*node] = 100.; fractureDots[2*node+1] = 45.;
    
    node = 51;
    
    fractureDots[2*node] = 113.22; fractureDots[2*node+1] = 43.8338;
    
    node = 52;
    
    fractureDots[2*node] = 127.511; fractureDots[2*node+1] = 42.3139;
    
    node = 53;
    
    fractureDots[2*node] = 142.068; fractureDots[2*node+1] = 40.4454;
    
    node = 54;
    
    fractureDots[2*node] = 155.692; fractureDots[2*node+1] = 38.2585;
    
    node = 55;
    
    fractureDots[2*node] = 166.721; fractureDots[2*node+1] = 35.8147;
    
    node = 56;
    
    fractureDots[2*node] = 172.932; fractureDots[2*node+1] = 33.2139;
    
    node = 57;
    
    fractureDots[2*node] = 171.452; fractureDots[2*node+1] = 30.601;
    
    node = 58;
    
    fractureDots[2*node] = 158.644; fractureDots[2*node+1] = 28.175;
    
    node = 59;
    
    fractureDots[2*node] = 129.998; fractureDots[2*node+1] = 26.1969;
    
    node = 60;
    
    fractureDots[2*node] = 80.; fractureDots[2*node+1] = 25.;
    
    node = 61;
    
    fractureDots[2*node] = 22.; fractureDots[2*node+1] = 25.;
}
