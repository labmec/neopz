#include <iostream>

#include "pzgmesh.h"
#include "pzgeopoint.h"
#include "TPZRefPatternDataBase.h"
#include "tpzgeoelrefpattern.h"
#include "TPZPoligonalChain.h"
#include "TPZPlaneFracture.h"

//to delete
#include "tpzchangeel.h"
#include "pzmaterial.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzstrmatrix.h"
#include "pzskylmat.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeotetrahedra.h"


using namespace std;

const int matPoint = -3;

//** just for visualize given dots in vtk */

void FillFractureDotsExampleEllipse(TPZVec<REAL> &fractureDots);
void FillFractureDotsExampleCrazy(TPZVec<REAL> &fractureDots);

//----------------------------------------------------------------------------------------------------------------------------------

#include "TPZTimer.h"

/** Exemplo da utilizacao da quadratura adaptativa (integral adaptativa) */
/*

#include "adapt.h"

class funcao
{
public:
    
    funcao()
    {
        
    }
    ~funcao()
    {
        
    }
    double operator()(double x)
    {
        double val = this->func(x);
        return val;
    }
    
    static REAL func(REAL x)
    {
        REAL val = sin(log10(x) + x*4.*atan(1.))*x*x;
        return val;
    }
};

int main(int argc, char * const argv[])
{
    Adapt intRule(1.E-100);
    
    funcao funcaoOBJ;
    
    const REAL a = 2.;
    const REAL b = 10.;
    double integr = intRule.integrate(funcaoOBJ,a,b);
    
    std::cout.precision(16);
    std::cout << std::endl << integr << std::endl;
    
    return 0;
}
*/

int main(int argc, char * const argv[])
{	
    std::cout << "\e";
    TPZTimer readRef("ReadingRefPatterns");
    readRef.start();
    
    //#define writeAgain
    #ifdef writeAgain
        gRefDBase.InitializeRefPatterns();
    #else
        std::ifstream inRefP("RefPatternsUsed.txt");
        gRefDBase.ReadRefPatternDBase("RefPatternsUsed.txt");
    #endif
    
    readRef.stop();
    std::cout << "DeltaT leitura refpatterns = " << readRef.seconds() << " s" << std::endl;
    
    double lw = 202.;
    double bulletDepthIni = 20.;
    double bulletDepthFin = 80.;    
    TPZVec< std::map<double,double> > pos_stress(5);
    //stretch #0
    pos_stress[0][0.]  = 2.;
    pos_stress[0][11.]  = 3.;
    //stretch #1
    pos_stress[1][11.]  = 4.;
    pos_stress[1][25.]  = 7.;
    //stretch #2
    pos_stress[2][25.]  = 3.;
    pos_stress[2][37.]  = 6.;
    //stretch #3
    pos_stress[3][37.]  = 6.;
    pos_stress[3][63.]  = 8.;
    //stretch #4
    pos_stress[4][63.] = 8.;
    pos_stress[4][210.] = 10.;
    TPZPlaneFracture plfrac(lw, bulletDepthIni, bulletDepthFin, pos_stress);
    TPZVec<REAL> fractureDots(0);
    FillFractureDotsExampleEllipse(fractureDots);
    
    TPZTimer clockIni2("PartyBegins2");
    clockIni2.start();    
    
    std::string vtkFile = "fracturePconstant.vtk";
    plfrac.RunThisFractureGeometry(fractureDots, vtkFile);
    
    clockIni2.stop();
    std::cout << "DeltaT get fracture cmesh = " << clockIni2.seconds() << " s" << std::endl;
  
    std::ofstream outRefP("RefPatternsUsed.txt");
    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}

int mainTestes(int argc, char * const argv[])
{
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ECube);
    
    int Qnodes = 8;
	
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(Qnodes);
    
    int nid = 0;
	TPZGeoNode Node;
    
    {    
        Node.SetNodeId(nid);
        Node.SetCoord(0,0.);
        Node.SetCoord(1,0.);
        Node.SetCoord(2,0.);
        gmesh->NodeVec()[nid] = Node;
        nid++;
        
        Node.SetNodeId(nid);
        Node.SetCoord(0,1.);
        Node.SetCoord(1,0.);
        Node.SetCoord(2,0.);
        gmesh->NodeVec()[nid] = Node;
        nid++;
        
        Node.SetNodeId(nid);
        Node.SetCoord(0,1.);
        Node.SetCoord(1,1.);
        Node.SetCoord(2,0.);
        gmesh->NodeVec()[nid] = Node;
        nid++;
        
        Node.SetNodeId(nid);
        Node.SetCoord(0,0.);
        Node.SetCoord(1,1.);
        Node.SetCoord(2,0.);
        gmesh->NodeVec()[nid] = Node;
        nid++;
        
        Node.SetNodeId(nid);
        Node.SetCoord(0,0.);
        Node.SetCoord(1,0.);
        Node.SetCoord(2,1.);
        gmesh->NodeVec()[nid] = Node;
        nid++;
        
        Node.SetNodeId(nid);
        Node.SetCoord(0,1.);
        Node.SetCoord(1,0.);
        Node.SetCoord(2,1.);
        gmesh->NodeVec()[nid] = Node;
        nid++;
        
        Node.SetNodeId(nid);
        Node.SetCoord(0,1.);
        Node.SetCoord(1,1.);
        Node.SetCoord(2,1.);
        gmesh->NodeVec()[nid] = Node;
        nid++;
        
        Node.SetNodeId(nid);
        Node.SetCoord(0,0.);
        Node.SetCoord(1,1.);
        Node.SetCoord(2,1.);
        gmesh->NodeVec()[nid] = Node;
        nid++;
	}
    
	int matId, elId = 0;
	TPZVec <int> Topol(Qnodes);
    
    matId = 1;
    for(int n = 0; n < Qnodes; n++) Topol[n] = n;
    TPZGeoEl * gel = new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (elId,Topol, matId,*gmesh);
    elId++;
    
    gmesh->BuildConnectivity();
    
    gel->CreateBCGeoEl(20, -1);
    gel->CreateBCGeoEl(25, -2);
    
    int nels = gmesh->NElements();
    for(int el = 0; el < nels; el++)
    {
        TPZVec<TPZGeoEl * > sons;
        gmesh->ElementVec()[el]->Divide(sons);
    }
    
    gmesh->Print();
    
    int foundN;
    for(int n = 0; n < gmesh->NNodes(); n++)
    {
        if(gmesh->NodeVec()[n].Id() != 18) continue;
        
        foundN = n;
    
        gmesh->NodeVec()[n].SetCoord(0, 0.7);
        gmesh->NodeVec()[n].SetCoord(1, 0.3);
        //gmesh->NodeVec()[n].SetCoord(2, 0.4);
        
    }
    
    nels = gmesh->NElements();
    for(int el = 0; el < nels; el++)
    {
        int target = -1;
        for(int eln = 0; eln < gmesh->ElementVec()[el]->NNodes(); eln++)
        {
            if(foundN == gmesh->ElementVec()[el]->NodeIndex(eln))
            {
                target = eln;
                break;
            }
        }
        if(target > -1)
        {
            TPZChangeEl::ChangeToQuarterPoint(gmesh, el, target);
        }
    }

    
    gmesh->Print();
//    
//    int elIndex = 0;
//    int targetSide = 0;
//    TPZChangeEl::ChangeToQuarterPoint(gmesh, elIndex, targetSide);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    ////////////////////////////////////
    
    int ddim = 3;
    cmesh->SetDimModel(ddim);
    cmesh->SetAllCreateFunctionsContinuous();
    
    REAL young = 1000.;
    REAL poisson = 0.;
    TPZVec<REAL> force(3,0.);
    
    TPZAutoPointer<TPZMaterial> materialQpoint = new TPZElasticity3D(1, young, poisson, force);
    cmesh->InsertMaterialObject(materialQpoint);
    
    ////BCs    
    TPZFMatrix<REAL> k(3,3,0.), f(3,1,0.);
    int dirichlet = 0;
    
    //    //teste 2: mov. corpo rigido (rotacao)
    {
        f(0,0) = 0.;
        f(1,0) = 0.;
        f(2,0) = 0.;
        TPZAutoPointer<TPZMaterial> materialMixedPoint1 = new TPZElasticity3D(-301, young, poisson, force);
        TPZBndCond * pontoDeApoio1 = new TPZBndCond(materialMixedPoint1,-1, dirichlet, k, f);
        cmesh->InsertMaterialObject(pontoDeApoio1);
        
        f(0,0) = 0.;
        f(1,0) = 0.;
        f(2,0) = 1.;
        TPZAutoPointer<TPZMaterial> materialMixedPoint2 = new TPZElasticity3D(-302, young, poisson, force);
        TPZBndCond * pontoDeApoio2 = new TPZBndCond(materialMixedPoint2,-2, dirichlet, k, f);
        cmesh->InsertMaterialObject(pontoDeApoio2);
    }
    cmesh->SetDefaultOrder(2);
    cmesh->AutoBuild();
    cmesh->Print();
    ///////////////////////////////////////
    
    ////Analysis
	TPZAnalysis an(cmesh);
    
	TPZSkylineStructMatrix skylin(cmesh); //caso simetrico
	TPZStepSolver<REAL> step;
	step.SetDirect(ECholesky);
    
    an.SetStructuralMatrix(skylin);
	an.SetSolver(step);
	an.Run();
    
    ////Post Processing
    TPZManVector<std::string,10> scalnames(6), vecnames(1);
    scalnames[0] = "DisplacementX";
    scalnames[1] = "DisplacementY";
    scalnames[2] = "DisplacementZ";
    scalnames[3] = "StressX";
    scalnames[4] = "StressY";
    scalnames[5] = "StressZ";
    
    vecnames[0] = "PrincipalStress";
    
    int div = 1;
    const int dim = 3;
    std::string vtkFile = "validandoRotQPts.vtk";
    an.DefineGraphMesh(dim,scalnames,vecnames,vtkFile);
    an.PostProcess(div,dim);
    
    
    return 0;
}

void FillFractureDotsExampleEllipse(TPZVec<REAL> &fractureDots)
{
    int nnodes = 80;
    
    fractureDots.Resize(2*nnodes, 0.);
    int node;
    double shiftZ = -100.;
    
    node = 0;
    
    fractureDots[2*node] = 0.5; fractureDots[2*node+1] = shiftZ + 95.;
    
    node = 1;
    
    fractureDots[2*node] = 5.; fractureDots[2*node+1] = shiftZ + 94.9852;
    
    node = 2;
    
    fractureDots[2*node] = 10.; fractureDots[2*node+1] = shiftZ + 94.9408;
    
    node = 3;
    
    fractureDots[2*node] = 15.; fractureDots[2*node+1] = shiftZ + 94.8667;
    
    node = 4;
    
    fractureDots[2*node] = 20.; fractureDots[2*node+1] = shiftZ + 94.7627;
    
    node = 5;
    
    fractureDots[2*node] = 25.; fractureDots[2*node+1] = shiftZ + 94.6286;
    
    node = 6;
    
    fractureDots[2*node] = 30.; fractureDots[2*node+1] = shiftZ + 94.4643;
    
    node = 7;
    
    fractureDots[2*node] = 35.; fractureDots[2*node+1] = shiftZ + 94.2692;
    
    node = 8;
    
    fractureDots[2*node] = 40.; fractureDots[2*node+1] = shiftZ + 94.0431;
    
    node = 9;
    
    fractureDots[2*node] = 45.; fractureDots[2*node+1] = shiftZ + 93.7854;
    
    node = 10;
    
    fractureDots[2*node] = 50.; fractureDots[2*node+1] = shiftZ + 93.4956;
    
    node = 11;
    
    fractureDots[2*node] = 55.; fractureDots[2*node+1] = shiftZ + 93.173;
    
    node = 12;
    
    fractureDots[2*node] = 60.; fractureDots[2*node+1] = shiftZ + 92.8169;
    
    node = 13;
    
    fractureDots[2*node] = 65.; fractureDots[2*node+1] = shiftZ + 92.4264;
    
    node = 14;
    
    fractureDots[2*node] = 70.; fractureDots[2*node+1] = shiftZ + 92.0006;
    
    node = 15;
    
    fractureDots[2*node] = 75.; fractureDots[2*node+1] = shiftZ + 91.5385;
    
    node = 16;
    
    fractureDots[2*node] = 80.; fractureDots[2*node+1] = shiftZ + 91.0387;
    
    node = 17;
    
    fractureDots[2*node] = 85.; fractureDots[2*node+1] = shiftZ + 90.4998;
    
    node = 18;
    
    fractureDots[2*node] = 90.; fractureDots[2*node+1] = shiftZ + 89.9204;
    
    node = 19;
    
    fractureDots[2*node] = 95.; fractureDots[2*node+1] = shiftZ + 89.2986;
    
    node = 20;
    
    fractureDots[2*node] = 100.; fractureDots[2*node+1] = shiftZ + 88.6323;
    
    node = 21;
    
    fractureDots[2*node] = 105.; fractureDots[2*node+1] = shiftZ + 87.9193;
    
    node = 22;
    
    fractureDots[2*node] = 110.; fractureDots[2*node+1] = shiftZ + 87.1567;
    
    node = 23;
    
    fractureDots[2*node] = 115.; fractureDots[2*node+1] = shiftZ + 86.3416;
    
    node = 24;
    
    fractureDots[2*node] = 120.; fractureDots[2*node+1] = shiftZ + 85.4702;
    
    node = 25;
    
    fractureDots[2*node] = 125.; fractureDots[2*node+1] = shiftZ + 84.5384;
    
    node = 26;
    
    fractureDots[2*node] = 130.; fractureDots[2*node+1] = shiftZ + 83.541;
    
    node = 27;
    
    fractureDots[2*node] = 135.; fractureDots[2*node+1] = shiftZ + 82.4721;
    
    node = 28;
    
    fractureDots[2*node] = 140.; fractureDots[2*node+1] = shiftZ + 81.3243;
    
    node = 29;
    
    fractureDots[2*node] = 145.; fractureDots[2*node+1] = shiftZ + 80.0886;
    
    node = 30;
    
    fractureDots[2*node] = 150.; fractureDots[2*node+1] = shiftZ + 78.7537;
    
    node = 31;
    
    fractureDots[2*node] = 155.; fractureDots[2*node+1] = shiftZ + 77.305;
    
    node = 32;
    
    fractureDots[2*node] = 160.; fractureDots[2*node+1] = shiftZ + 75.7233;
    
    node = 33;
    
    fractureDots[2*node] = 165.; fractureDots[2*node+1] = shiftZ + 73.9822;
    
    node = 34;
    
    fractureDots[2*node] = 170.; fractureDots[2*node+1] = shiftZ + 72.0442;
    
    node = 35;
    
    fractureDots[2*node] = 175.; fractureDots[2*node+1] = shiftZ + 69.8515;
    
    node = 36;
    
    fractureDots[2*node] = 180.; fractureDots[2*node+1] = shiftZ + 67.3077;
    
    node = 37;
    
    fractureDots[2*node] = 185.; fractureDots[2*node+1] = shiftZ + 64.2256;
    
    node = 38;
    
    fractureDots[2*node] = 190.; fractureDots[2*node+1] = shiftZ + 60.125;
    
    node = 39;
    
    fractureDots[2*node] = 195.; fractureDots[2*node+1] = shiftZ + 50.;
    
    node = 40;
    
    fractureDots[2*node] = 195.; fractureDots[2*node+1] = shiftZ + 42.;
    
    node = 41;
    
    fractureDots[2*node] = 190.; fractureDots[2*node+1] = shiftZ + 39.875;
    
    node = 42;
    
    fractureDots[2*node] = 185.; fractureDots[2*node+1] = shiftZ + 35.7744;
    
    node = 43;
    
    fractureDots[2*node] = 180.; fractureDots[2*node+1] = shiftZ + 32.6923;
    
    node = 44;
    
    fractureDots[2*node] = 175.; fractureDots[2*node+1] = shiftZ + 30.1485;
    
    node = 45;
    
    fractureDots[2*node] = 170.; fractureDots[2*node+1] = shiftZ + 27.9558;
    
    node = 46;
    
    fractureDots[2*node] = 165.; fractureDots[2*node+1] = shiftZ + 26.0178;
    
    node = 47;
    
    fractureDots[2*node] = 160.; fractureDots[2*node+1] = shiftZ + 24.2767;
    
    node = 48;
    
    fractureDots[2*node] = 155.; fractureDots[2*node+1] = shiftZ + 22.695;
    
    node = 49;
    
    fractureDots[2*node] = 150.; fractureDots[2*node+1] = shiftZ + 21.2463;
    
    node = 50;
    
    fractureDots[2*node] = 145.; fractureDots[2*node+1] = shiftZ + 19.9114;
    
    node = 51;
    
    fractureDots[2*node] = 140.; fractureDots[2*node+1] = shiftZ + 18.6757;
    
    node = 52;
    
    fractureDots[2*node] = 135.; fractureDots[2*node+1] = shiftZ + 17.5279;
    
    node = 53;
    
    fractureDots[2*node] = 130.; fractureDots[2*node+1] = shiftZ + 16.459;
    
    node = 54;
    
    fractureDots[2*node] = 125.; fractureDots[2*node+1] = shiftZ + 15.4616;
    
    node = 55;
    
    fractureDots[2*node] = 120.; fractureDots[2*node+1] = shiftZ + 14.5298;
    
    node = 56;
    
    fractureDots[2*node] = 115.; fractureDots[2*node+1] = shiftZ + 13.6584;
    
    node = 57;
    
    fractureDots[2*node] = 110.; fractureDots[2*node+1] = shiftZ + 12.8433;
    
    node = 58;
    
    fractureDots[2*node] = 105.; fractureDots[2*node+1] = shiftZ + 12.0807;
    
    node = 59;
    
    fractureDots[2*node] = 100.; fractureDots[2*node+1] = shiftZ + 11.3677;
    
    node = 60;
    
    fractureDots[2*node] = 95.; fractureDots[2*node+1] = shiftZ + 10.7014;
    
    node = 61;
    
    fractureDots[2*node] = 90.; fractureDots[2*node+1] = shiftZ + 10.0796;
    
    node = 62;
    
    fractureDots[2*node] = 85.; fractureDots[2*node+1] = shiftZ + 9.50016;
    
    node = 63;
    
    fractureDots[2*node] = 80.; fractureDots[2*node+1] = shiftZ + 8.96134;
    
    node = 64;
    
    fractureDots[2*node] = 75.; fractureDots[2*node+1] = shiftZ + 8.46154;
    
    node = 65;
    
    fractureDots[2*node] = 70.; fractureDots[2*node+1] = shiftZ + 7.99937;
    
    node = 66;
    
    fractureDots[2*node] = 65.; fractureDots[2*node+1] = shiftZ + 7.57359;
    
    node = 67;
    
    fractureDots[2*node] = 56.; fractureDots[2*node+1] = shiftZ + 7.18313;
    
    node = 68;
    
    fractureDots[2*node] = 55.; fractureDots[2*node+1] = shiftZ - 7.;//+ 6.82703;
    
    node = 69;
    
    fractureDots[2*node] = 50.; fractureDots[2*node+1] = shiftZ + 6.50444;
    
    node = 70;
    
    fractureDots[2*node] = 45.; fractureDots[2*node+1] = shiftZ + 6.21462;
    
    node = 71;
    
    fractureDots[2*node] = 40.; fractureDots[2*node+1] = shiftZ + 5.95692;
    
    node = 72;
    
    fractureDots[2*node] = 31.; fractureDots[2*node+1] = shiftZ + 5.73079;
    
    node = 73;
    
    fractureDots[2*node] = 29.; fractureDots[2*node+1] = shiftZ - 6.;//5.53573;
    
    node = 74;
    
    fractureDots[2*node] = 25.; fractureDots[2*node+1] = shiftZ + 5.37135;
    
    node = 75;
    
    fractureDots[2*node] = 20.; fractureDots[2*node+1] = shiftZ + 5.23731;
    
    node = 76;
    
    fractureDots[2*node] = 15.; fractureDots[2*node+1] = shiftZ + 5.13333;
    
    node = 77;
    
    fractureDots[2*node] = 10.; fractureDots[2*node+1] = shiftZ + 5.05921;
    
    node = 78;
    
    fractureDots[2*node] = 5.; fractureDots[2*node+1] = shiftZ + 5.0148;
    
    node = 79;
    
    fractureDots[2*node] = 0.5; fractureDots[2*node+1] = shiftZ + 5.;
}


void FillFractureDotsExampleCrazy(TPZVec<REAL> &fractureDots)
{
    int nnodes = 62;
    
    fractureDots.Resize(2*nnodes, 0.);
    int node;
    
    double shiftZ = -100.;
    
    ///
    node = 0;
    double desloc = 0.;
    
    fractureDots[2*node] = 5.; fractureDots[2*node+1] = shiftZ + 80.;
    
    node = 1;
    
    fractureDots[2*node] = 9.02368; fractureDots[2*node+1] = shiftZ + 82.1856;
    
    node = 2;
    
    fractureDots[2*node] = 25.6151; fractureDots[2*node+1] = shiftZ + 82.9122;
    
    node = 3;
    
    fractureDots[2*node] = 48.4505; fractureDots[2*node+1] = shiftZ + 82.5833;
    
    node = 4;
    
    fractureDots[2*node] = 74.1186; fractureDots[2*node+1] = shiftZ + 81.5259;
    
    node = 5;
    
    fractureDots[2*node] = 100. + desloc; fractureDots[2*node+1] = shiftZ + 80.;//<----------------------
    
    node = 6;
    
    fractureDots[2*node] = 124.157; fractureDots[2*node+1] = shiftZ + 78.2079;
    
    node = 7;
    
    fractureDots[2*node] = 145.23; fractureDots[2*node+1] = shiftZ + 76.3027;
    
    node = 8;
    
    fractureDots[2*node] = 162.349; fractureDots[2*node+1] = shiftZ + 74.3954;
    
    node = 9;
    
    fractureDots[2*node] = 175.043; fractureDots[2*node+1] = shiftZ + 72.5625;
    
    node = 10;
    
    fractureDots[2*node] = 183.172; fractureDots[2*node+1] = shiftZ + 70.8516;
    
    node = 11;
    
    fractureDots[2*node] = 186.85; fractureDots[2*node+1] = shiftZ + 69.2873;
    
    node = 12;
    
    fractureDots[2*node] = 186.391; fractureDots[2*node+1] = shiftZ + 67.8763;
    
    node = 13;
    
    fractureDots[2*node] = 182.253; fractureDots[2*node+1] = shiftZ + 66.6118;
    
    node = 14;
    
    fractureDots[2*node] = 174.99; fractureDots[2*node+1] = shiftZ + 65.4771;
    
    node = 15;
    
    fractureDots[2*node] = 165.212; fractureDots[2*node+1] = shiftZ + 64.4495;
    
    node = 16;
    
    fractureDots[2*node] = 153.552; fractureDots[2*node+1] = shiftZ + 63.5025;
    
    node = 17;
    
    fractureDots[2*node] = 140.633; fractureDots[2*node+1] = shiftZ + 62.609;
    
    node = 18;
    
    fractureDots[2*node] = 127.05; fractureDots[2*node+1] = shiftZ + 61.7426;
    
    node = 19;
    
    fractureDots[2*node] = 113.346; fractureDots[2*node+1] = shiftZ + 60.8796;
    
    node = 20;
    
    fractureDots[2*node] = 100. + desloc; fractureDots[2*node+1] = shiftZ + 60.;//<----------------------
    
    node = 21;
    
    fractureDots[2*node] = 87.418; fractureDots[2*node+1] = shiftZ + 59.0884;
    
    node = 22;
    
    fractureDots[2*node] = 75.9253; fractureDots[2*node+1] = shiftZ + 58.1348;
    
    node = 23;
    
    fractureDots[2*node] = 65.7644; fractureDots[2*node+1] = shiftZ + 57.1342;
    
    node = 24;
    
    fractureDots[2*node] = 57.0956; fractureDots[2*node+1] = shiftZ + 56.0873;
    
    node = 25;
    
    fractureDots[2*node] = 50.; fractureDots[2*node+1] = shiftZ + 55.;
    
    node = 26;
    
    fractureDots[2*node] = 44.4857; fractureDots[2*node+1] = shiftZ + 53.8827;
    
    node = 27;
    
    fractureDots[2*node] = 40.495; fractureDots[2*node+1] = shiftZ + 52.75;
    
    node = 28;
    
    fractureDots[2*node] = 37.9143; fractureDots[2*node+1] = shiftZ + 51.6198;
    
    node = 29;
    
    fractureDots[2*node] = 36.5847; fractureDots[2*node+1] = shiftZ + 50.5124;
    
    node = 30;
    
    fractureDots[2*node] = 36.3137; fractureDots[2*node+1] = shiftZ + 49.4496;
    
    node = 31;
    
    fractureDots[2*node] = 36.888; fractureDots[2*node+1] = shiftZ + 48.4535;
    
    node = 32;
    
    fractureDots[2*node] = 38.0857; fractureDots[2*node+1] = shiftZ + 47.5455;
    
    node = 33;
    
    fractureDots[2*node] = 39.6892; fractureDots[2*node+1] = shiftZ + 46.7455;
    
    node = 34;
    
    fractureDots[2*node] = 41.4972; fractureDots[2*node+1] = shiftZ + 46.0702;
    
    node = 35;
    
    fractureDots[2*node] = 43.3363; fractureDots[2*node+1] = shiftZ + 45.5328;
    
    node = 36;
    
    fractureDots[2*node] = 45.0707; fractureDots[2*node+1] = shiftZ + 45.1417;
    
    node = 37;
    
    fractureDots[2*node] = 46.6112; fractureDots[2*node+1] = shiftZ + 44.8996;
    
    node = 38;
    
    fractureDots[2*node] = 47.9218; fractureDots[2*node+1] = shiftZ + 44.8032;
    
    node = 39;
    
    fractureDots[2*node] = 49.0243; fractureDots[2*node+1] = shiftZ + 44.8424;
    
    node = 40;
    
    fractureDots[2*node] = 50.; fractureDots[2*node+1] = shiftZ + 45.;
    
    node = 41;
    
    fractureDots[2*node] = 50.9889; fractureDots[2*node+1] = shiftZ + 45.2519;
    
    node = 42;
    
    fractureDots[2*node] = 52.1852; fractureDots[2*node+1] = shiftZ + 45.5668;
    
    node = 43;
    
    fractureDots[2*node] = 53.8292; fractureDots[2*node+1] = shiftZ + 45.9072;
    
    node = 44;
    
    fractureDots[2*node] = 56.1951; fractureDots[2*node+1] = shiftZ + 46.2297;
    
    node = 45;
    
    fractureDots[2*node] = 59.5747; fractureDots[2*node+1] = shiftZ + 46.4863;
    
    node = 46;
    
    fractureDots[2*node] = 64.2563; fractureDots[2*node+1] = shiftZ + 46.626;
    
    node = 47;
    
    fractureDots[2*node] = 70.4977; fractureDots[2*node+1] = shiftZ + 46.5964;
    
    node = 48;
    
    fractureDots[2*node] = 78.4951; fractureDots[2*node+1] = shiftZ + 46.3461;
    
    node = 49;
    
    fractureDots[2*node] = 88.3449; fractureDots[2*node+1] = shiftZ + 45.8275;
    
    node = 50;
    
    fractureDots[2*node] = 100.; fractureDots[2*node+1] = shiftZ + 45.;
    
    node = 51;
    
    fractureDots[2*node] = 113.22; fractureDots[2*node+1] = shiftZ + 43.8338;
    
    node = 52;
    
    fractureDots[2*node] = 127.511; fractureDots[2*node+1] = shiftZ + 42.3139;
    
    node = 53;
    
    fractureDots[2*node] = 142.068; fractureDots[2*node+1] = shiftZ + 40.4454;
    
    node = 54;
    
    fractureDots[2*node] = 155.692; fractureDots[2*node+1] = shiftZ + 38.2585;
    
    node = 55;
    
    fractureDots[2*node] = 166.721; fractureDots[2*node+1] = shiftZ + 35.8147;
    
    node = 56;
    
    fractureDots[2*node] = 172.932; fractureDots[2*node+1] = shiftZ + 33.2139;
    
    node = 57;
    
    fractureDots[2*node] = 171.452; fractureDots[2*node+1] = shiftZ + 30.601;
    
    node = 58;
    
    fractureDots[2*node] = 158.644; fractureDots[2*node+1] = shiftZ + 28.175;
    
    node = 59;
    
    fractureDots[2*node] = 129.998; fractureDots[2*node+1] = shiftZ + 26.1969;
    
    node = 60;
    
    fractureDots[2*node] = 80.; fractureDots[2*node+1] = shiftZ + 25.;
    
    node = 61;
    
    fractureDots[2*node] = 22.; fractureDots[2*node+1] = shiftZ + 25.;
}
