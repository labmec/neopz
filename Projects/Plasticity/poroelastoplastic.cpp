#include "poroelastoplastic.h"
#include "pzelctemp.h" // TPZIntelGen<TSHAPE>
#include "pzshapecube.h" // TPZShapeCube
#include "TPZLadeKim.h"
#include "pzmat2dlin.h"
#include "pzporoanalysis.h"
#include "pzbfilestream.h"
#include <sstream>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>
#include "pzelasmat.h"
#include "TPZVTKGeoMesh.h"
#include "BrazilianTestGeoMesh.h"
#include "TPZProjectEllipse.h"

#include "WellBoreAnalysis.h"

#define PV

void VisualizeSandlerDimaggio(std::stringstream &FileName, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD);

void SolverSetUp2(TPZAnalysis &an, TPZCompMesh *fCmesh);
void SetUPPostProcessVariables2(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames );
void ManageIterativeProcess2(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost);
void ManageIterativeProcessPesoProprio2(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost);
void RotationMatrix(TPZFMatrix<REAL> &R, double thetaRad, int axis);
void RotateMatrix(TPZFMatrix<REAL> &Mat, double thetaRad,int rotateaboutaxes);
void RotateMesh(TPZGeoMesh &geomesh, REAL angleDegree,int rotateabout);
void calcSDBar();
TPZGeoMesh * BarMesh(int h);
void Cmesh(TPZCompMesh *CMesh, TPZMaterial * mat,REAL theta,int axes);
void wellcmesh();
void wellboreanalysis();

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

using namespace pzshape; // needed for TPZShapeCube and related classes

#include <math.h>



//#define MACOS
#ifdef MACOS

#include <iostream>
#include <math.h>
#include <signal.h>
#include <fenv.h>
#include <xmmintrin.h>

#define ENABLE_FPO_EXCEPTIONS _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);


#define DECLARE_FPO_HANDLER_FUNC void InvalidFPOHandler(int signo) {\
switch(signo) {\
case SIGFPE: std::cout << "ERROR : Invalid Arithmetic operation." << std::endl; break;\
}\
exit(signo);\
}

#define ATTACH_FPO_SIGNAL struct sigaction act = {};\
act.sa_handler = InvalidFPOHandler;\
sigaction(SIGFPE, &act, NULL);


DECLARE_FPO_HANDLER_FUNC;
#endif

int startfrom = 1;
#include "TPZTimer.h"


int main2 ()
{
    
#ifdef MACOS
    
    ENABLE_FPO_EXCEPTIONS;
    ATTACH_FPO_SIGNAL;
    
#endif
    
    TPZTimer time1,time2;
    time1.start();
    InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();
    
    time2.start();
    TPZWellBoreAnalysis well;
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;

    REAL computedquarter = 7.05761678496926;
    REAL sqj2_refine;
    const int nsubsteps = 5;
    std::cout << std::setprecision(15);
    if (startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        TPZManVector<STATE,3> confinementEffective(3,0.), confinementTotal(3,0.);
        REAL SH,Sh,SV;
        Sh=-62.1;
        SH=-45.9;
        SV=-48.2;
        
        confinementEffective[0] = Sh;
        confinementEffective[1] = SH;
        confinementEffective[2] = SV;
        REAL effectiveWellPressure = 19.5; // 19.5 ou 23.4 ou 28.9
        STATE biotcoef = 0.659;
        well.SetBiotCoefficient(biotcoef);
        
        STATE WellPressure = effectiveWellPressure/(1.-biotcoef);
        for (int i=0; i<3; i++) {
            confinementTotal[i] = confinementEffective[i]-biotcoef*WellPressure;
        }
        
        well.SetConfinementTotalStresses(confinementTotal, WellPressure);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        

      
        bool modelMC =true;
				
        if (modelMC)
        {
//            REAL cohesion = 7.;
//            REAL Phi = 0.25;
            REAL cohesion = A - C;
            REAL Phi = B*C;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            sqj2_refine=0.001;
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            sqj2_refine=0.0001;
            
        }
			
        int divisions = 20;
        REAL delx = 0.2*innerradius*M_PI_2/divisions;
        TPZManVector<int,2> numdiv(2,divisions);
        numdiv[1] = 40;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->CreateMesh();
        int porder = 2;
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        //well.LinearConfiguration(1);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Wellbore0.bin");
        well.Write(save);
    }
    
    
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        
    }
    
    
    time2.stop();
    std::cout << "\n tempo 0 = "<< time2.seconds()<< std::endl;
    time2.reset();
    time2.start();
    
    if (startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore0.bin");
        well.Read(read);
    }
    return 0;
    if (startfrom <= 1)
    {
        
        int nsteps = 10;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        TPZBFileStream save;
        save.OpenWrite("Wellbore1.bin");
        well.Write(save);
        
    }
    
    
    
    
    
    
    extern std::map<int,long> gF1Stat;
    extern std::map<int,long> gF2Stat;
    
    std::map<int,long>::iterator it;
    for (it=gF1Stat.begin(); it!=gF1Stat.end(); it++) {
        std::cout << "Numero de iteracoes F1 " << it->first << " numero de ocorrencias " << it->second << std::endl;
    }
    for (it=gF2Stat.begin(); it!=gF2Stat.end(); it++) {
        std::cout << "Numero de iteracoes F2 " << it->first << " numero de ocorrencias " << it->second << std::endl;
    }
    
    gF1Stat.clear();
    gF2Stat.clear();
    time2.stop();
    std::cout << "\n tempo 1 = "<< time2.seconds()<< std::endl;
    time2.reset();
    time2.start();
    
    if (startfrom == 2)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore1.bin");
        well.Read(read);
    }
    if (startfrom <= 2) {
        
        //well.PRefineElementAbove(0.0001, 3);
        well.DivideElementsAbove(0.0001);
        well.PRefineElementAbove(0.0001, 3);
        well.ExecuteSimulation(nsubsteps);
        REAL analyticarea = M_PI*(outerradius*outerradius-innerradius*innerradius)/4.;
        REAL originalarea = well.GetCurrentConfig()->ComputeTotalArea();
        REAL openingangle = well.GetCurrentConfig()->OpeningAngle(0.00000001);
        REAL plastifiedarea = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(sqj2_refine);
        std::cout << "Analytical area " << analyticarea << " computed area " << originalarea << std::endl;
        std::cout << "Analytical - computed area " << analyticarea-originalarea << std::endl;
        std::cout << "Plastified area " << plastifiedarea << std::endl;
        std::cout << "Opening angle " << openingangle << std::endl;
        std::cout << "Saving Wellbore2.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore2.bin");
        well.Write(save);
    }
    time2.stop();
    std::cout << "\n tempo 2 = "<< time2.seconds()<< std::endl;
    time2.reset();
    time2.start();
    if (startfrom == 3) {
        TPZBFileStream read;
        read.OpenRead("Wellbore2.bin");
        well.Read(read);
    }
    
    for (it=gF1Stat.begin(); it!=gF1Stat.end(); it++) {
        std::cout << "Numero de iteracoes F1 " << it->first << " numero de ocorrencias " << it->second << std::endl;
    }
    for (it=gF2Stat.begin(); it!=gF2Stat.end(); it++) {
        std::cout << "Numero de iteracoes F2 " << it->first << " numero de ocorrencias " << it->second << std::endl;
    }
    
    gF1Stat.clear();
    gF2Stat.clear();
    
    if (startfrom <= 3)
    {
        // valor de a e b para sqJ2 = 0.00025 E USANDO Pef = 23.4
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.010;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.79;
        //vvalor de a e b para sqJ2 = 0.0005 p 19.5
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.03409;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.829545;
        //vvalor de a e b para sqJ2 = 0.0007 p 19.5
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.014;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.90;
        
        
        std::multimap<REAL, REAL> polygonalChainbase, polygonalChain;
        well.GetJ2Isoline(sqj2_refine, polygonalChainbase);
        REAL maxy = well.GetCurrentConfig()->MaxYfromLastBreakout();
        for (std::multimap<REAL, REAL>::iterator it = polygonalChainbase.begin(); it != polygonalChainbase.end(); it++) {
            if (it->second < maxy) {
                polygonalChain.insert(std::make_pair(it->first,it->second));
            }
        }
        if (polygonalChain.size()!=0)
        {
            TPZProjectEllipse ellips(polygonalChain);
            TPZManVector<REAL,2> center(2),ratios(2),verify(2);
            ellips.StandardFormatForSimpleEllipse(center, ratios);
            verify[0] = ratios[0]/well.GetCurrentConfig()->fInnerRadius;
            verify[1] = ratios[1]/well.GetCurrentConfig()->fInnerRadius;
            REAL a = ratios[0];
            REAL b = ratios[1];
            well.AddEllipticBreakout(a, b);
            
        }
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        //well.PostProcess(1);
        well.ExecuteSimulation(nsubsteps);
        well.DivideElementsAbove(0.0001);
        well.PRefineElementAbove(0.0001, 3);
        well.ExecuteSimulation(nsubsteps);
        
        REAL analyticarea = M_PI*(outerradius*outerradius-innerradius*innerradius)/4.;
        REAL originalarea = well.GetCurrentConfig()->ComputeTotalArea();
        REAL openingangle = well.GetCurrentConfig()->OpeningAngle(0.0001);
        REAL plastifiedarea = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(sqj2_refine);
        std::cout << "Analytical area " << analyticarea << " computed area " << originalarea << std::endl;
        std::cout << "computed quarter - domain area " << computedquarter-originalarea << std::endl;
        std::cout << "Plastified area " << plastifiedarea << std::endl;
        std::cout << "Opening angle " << openingangle << std::endl;
        
        
        std::cout << "Saving Wellbore3.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore3.bin");
        well.Write(save);
        
        if(startfrom == 3)
        {
            TPZFileStream write;
            write.OpenWrite("afterread.txt");
            well.GetCurrentConfig()->Write(write);
        }
        else {
            TPZFileStream write;
            write.OpenWrite("afterrun.txt");
            well.GetCurrentConfig()->Write(write);
        }
        
    }
    time2.stop();
    std::cout << "\n tempo 3 = "<< time2.seconds()<< std::endl;
    time2.reset();
    time2.start();
    if (startfrom == 4)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore3.bin");
        well.Read(read);
    }
    for (it=gF1Stat.begin(); it!=gF1Stat.end(); it++) {
        std::cout << "Numero de iteracoes F1 " << it->first << " numero de ocorrencias " << it->second << std::endl;
    }
    for (it=gF2Stat.begin(); it!=gF2Stat.end(); it++) {
        std::cout << "Numero de iteracoes F2 " << it->first << " numero de ocorrencias " << it->second << std::endl;
    }
    
    gF1Stat.clear();
    gF2Stat.clear();
    
    if (startfrom <= 4)
    {
        // valor de a e b para sqJ2 = 0.00025 E USANDO Pef = 23.4
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.040;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.46;
        //vvalor de a e b para sqJ2 = 0.0005 p 19.5
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.10227;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.596591;
        //vvalor de a e b para sqJ2 = 0.0007 Pef = 19.5
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.044;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.66;
        std::multimap<REAL, REAL> polygonalChainbase, polygonalChain;
        well.GetJ2Isoline(sqj2_refine, polygonalChainbase);
        REAL innerradius = well.GetCurrentConfig()->fInnerRadius;
        REAL maxy = well.GetCurrentConfig()->MaxYfromLastBreakout();
        int numgood = 0;
        int numbad = 0;
        for (std::multimap<REAL, REAL>::iterator it = polygonalChainbase.begin(); it != polygonalChainbase.end(); it++) {
            TPZManVector<REAL,2> co(2);
            co[0] = 0.;
            co[1] = it->second;
            well.GetCurrentConfig()->ProjectNode(co);
            std::cout << "distance " << (it->first-co[0])/innerradius << endl;
            if (it->second < maxy && it->first-co[0] > 0.01*innerradius) {
                polygonalChain.insert(std::make_pair(it->first,it->second));
                numgood++;
            }
            else{
                numbad++;
            }
        }
        TPZFMatrix<REAL> polybase(polygonalChainbase.size(),2);
        cout << "base chain\n";
        int i=0;
        for (std::multimap<REAL, REAL>::iterator it = polygonalChainbase.begin(); it != polygonalChainbase.end(); it++,i++) {
            std::cout << it->first << " " << it->second << endl;
            polybase(i,0) = it->first;
            polybase(i,1) = it->second;
        }
        polybase.Print("Polybase = " , cout , EMathematicaInput);
        TPZFMatrix<REAL> poly(polygonalChain.size(),2);
        cout << "filtered chain\n";
        
        i=0;
        for (std::multimap<REAL, REAL>::iterator it = polygonalChain.begin(); it != polygonalChain.end(); it++,i++) {
            std::cout << it->first << " " << it->second << endl;
            poly(i,0) = it->first;
            poly(i,1) = it->second;
        }
        poly.Print("Poly = " , cout , EMathematicaInput);
        TPZProjectEllipse ellips(polygonalChain);
        TPZManVector<REAL,2> center(2),ratios(2),verify(2);
        ellips.StandardFormatForSimpleEllipse(center, ratios);
        cout << "chain fitting\n";
        for (std::multimap<REAL, REAL>::iterator it = polygonalChain.begin(); it != polygonalChain.end(); it++) {
            REAL xellips = ratios[0]*sqrt(1-it->second*it->second/ratios[1]/ratios[1]);
            std::cout << it->first << " " << xellips << " error " <<  (it->first-xellips)/innerradius << endl;
        }
        well.AddEllipticBreakout(ratios[0], ratios[1]);
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        well.PostProcess(1);
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteSimulation(nsubsteps);
        well.DivideElementsAbove(0.0001);
        well.PRefineElementAbove(0.0001, 3);
        well.ExecuteSimulation(nsubsteps);
        
        REAL analyticarea = M_PI*(outerradius*outerradius-innerradius*innerradius)/4.;
        REAL originalarea = well.GetCurrentConfig()->ComputeTotalArea();
        REAL openingangle = well.GetCurrentConfig()->OpeningAngle(0.000001);
        REAL plastifiedarea = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(sqj2_refine);
        std::cout << "Analytical area " << analyticarea << " computed area " << originalarea << std::endl;
        std::cout << "computed quarter - domain area " << computedquarter-originalarea << std::endl;
        std::cout << "Plastified area " << plastifiedarea << std::endl;
        std::cout << "Opening angle " << openingangle << std::endl;
        
        
        std::cout << "Saving Wellbore4.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore4.bin");
        well.Write(save);
    }
    time2.stop();
    std::cout << "\n tempo 4 = "<< time2.seconds()<< std::endl;
    time2.reset();
    time2.start();
    if (startfrom == 5) {
        TPZBFileStream read;
        read.OpenRead("Wellbore4.bin");
        well.Read(read);
    }
    if (startfrom <= 5)
    {
        // valor de a e b para sqJ2 = 0.00025 E USANDO Pef = 23.4
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.082;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.30;
        // valor de a e b para sqJ2 = 0.0005 P = 19.5
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.17045;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.465909;
        // valor de a e b para sqJ2 = 0.0007 P = 19.5
        //        REAL a = well.GetCurrentConfig()->fInnerRadius*1.099;
        //        REAL b = well.GetCurrentConfig()->fInnerRadius*0.419;
        std::multimap<REAL, REAL> polygonalChainbase, polygonalChain;
        well.GetJ2Isoline(sqj2_refine, polygonalChainbase);
        REAL maxy = well.GetCurrentConfig()->MaxYfromLastBreakout();
        
        for (std::multimap<REAL, REAL>::iterator it = polygonalChainbase.begin(); it != polygonalChainbase.end(); it++) {
            TPZManVector<REAL,2> co(2);
            co[0] = 0.;
            co[1] = it->second;
            well.GetCurrentConfig()->ProjectNode(co);
            std::cout << "distance " << (it->first-co[0])/innerradius << endl;
            
            if (it->second < maxy && it->first-co[0] > 0.01*innerradius) {
                polygonalChain.insert(std::make_pair(it->first,it->second));
            }
        }
        TPZProjectEllipse ellips(polygonalChain);
        TPZManVector<REAL,2> center(2),ratios(2),verify(2);
        ellips.StandardFormatForSimpleEllipse(center, ratios);
        well.AddEllipticBreakout(ratios[0], ratios[1]);
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        well.PostProcess(1);
        well.ExecuteSimulation(nsubsteps);
        well.DivideElementsAbove(0.0001);
        well.PRefineElementAbove(0.0001, 3);
        well.ExecuteSimulation(nsubsteps);
        
        REAL analyticarea = M_PI*(outerradius*outerradius-innerradius*innerradius)/4.;
        REAL originalarea = well.GetCurrentConfig()->ComputeTotalArea();
        REAL plastifiedarea = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(sqj2_refine);
        REAL openingangle = well.GetCurrentConfig()->OpeningAngle(0.0001);
        std::cout << "Analytical area " << analyticarea << " computed area " << originalarea << std::endl;
        std::cout << "computed quarter - domain area " << computedquarter-originalarea << std::endl;
        std::cout << "Plastified area " << plastifiedarea << std::endl;
        std::cout << "Opening angle " << openingangle << std::endl;
        
        
        std::cout << "Saving Wellbore5.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore5.bin");
        well.Write(save);
    }
    time2.stop();
    std::cout << "\n tempo 5 = "<< time2.seconds()<< std::endl;
    time2.reset();
    time2.start();
    if (startfrom == 6)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore5.bin");
        well.Read(read);
        
    }
    time1.stop();
    cout << "tempo de simulação em segundos = "<<time1.seconds() <<endl;
    return 0;
    if (startfrom <= 6)
    {
        //        well.VerifyGlobalEquilibrium();
        if(0)
        {
            
            TPZStack<std::string> postprocess;
            postprocess.Push("I1J2Stress");
            TPZFMatrix<STATE> valuetable;
            //TPZManVector<REAL,3> x(3,0.);
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(well.GetCurrentConfig()->fCMesh.ElementVec()[0]);
            if (!intel) {
                DebugStop();
            }
            TPZIntPoints &intpoints = intel->GetIntegrationRule();
            TPZManVector<REAL> ksi(2,0.),xco(3,0.);
            REAL weight;
            intpoints.Point(0, ksi, weight);
            TPZMaterialData data;
            intel->InitMaterialData(data);
            intel->ComputeRequiredData(data, ksi);
            TPZManVector<long> memindices(intpoints.NPoints());
            intel->GetMemoryIndices(memindices);
            data.intGlobPtIndex = memindices[0];
            TPZMaterial *mat = intel->Material();
            int varindex = mat->VariableIndex("I1J2Stress");
            int nvar = mat->NSolutionVariables(varindex);
            TPZManVector<STATE> post(nvar);
            mat->Solution(data, varindex, post);
            std::cout << "Post processed " << post << std::endl;
            intel->Reference()->X(ksi, xco);
            //x[0] = 1.1;
            well.PostProcessedValues(xco , postprocess, valuetable);
            valuetable.Print("Post processed I1=J2",std::cout);
            xco.Fill(0.);
            well.PostProcessedValues(xco , postprocess, valuetable);
            valuetable.Print("Post processed I1=J2",std::cout);
        }
        
        //well.ChangeMaterialId(-2, -6);
        well.DeleteElementsAbove(0.0001);
        well.ChangeMaterialId(-6, -2);
        well.ExecuteSimulation(nsubsteps);
        std::cout << "Saving Wellbore7.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore7.bin");
        well.Write(save);
        
    }
    time2.stop();
    std::cout << "\n tempo 6 = "<< time2.seconds()<< std::endl;
    time2.reset();
    time2.start();
    if (startfrom == 7)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore7.bin");
        well.Read(read);
    }
    time1.stop();
    cout << "tempo de simulação em segundos = "<<time1.seconds() <<endl;
    return 0;
    	
    return EXIT_SUCCESS;
    
    
}

#include "TPZGenSpecialGrid.h"

void BuildPlasticSurface(TPZCompMesh *cmesh, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD);

void VisualizeSandlerDimaggio(std::stringstream &fileName, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD)
{
    TPZVec<REAL> coords(0);
    TPZGeoMesh *gmesh = TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(coords, 0.001,1);
    TPZCompMesh *cgrid = new TPZCompMesh(gmesh);
    TPZManVector<STATE> force(3,0.);
    
    for (int imat=1; imat<4; imat++) {
        TPZMat2dLin *mat = new TPZMat2dLin(imat);
        TPZFNMatrix<9,STATE> xk(3,3,0.),xc(3,3,0.),xf(3,1,0.);
        mat->SetMaterial(xk,xc,xf);
        //    TPZMaterial * mat = new TPZElasticity3D(1,1.e5,0.2,force);
        cgrid->InsertMaterialObject(mat);
    }
    cgrid->AutoBuild();
    TPZFMatrix<REAL> elsol(cgrid->NElements(),1,0.);
    
    cgrid->ElementSolution() = elsol;
    
    TPZAnalysis an(cgrid);
    std::stringstream vtkfilename;
    vtkfilename << fileName.str();
    vtkfilename << ".vtk";
    std::ofstream meshout(vtkfilename.str().c_str());
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh,meshout);
    BuildPlasticSurface(cgrid,pSD);
    TPZStack<std::string> scalnames,vecnames;
    vecnames.Push("state");
    scalnames.Push("Error");
    
    an.DefineGraphMesh(2, scalnames, vecnames, "plot.vtk");
    an.PostProcess(0);
    TPZPlasticState<REAL> state = pSD->GetState();
    for (REAL alfa = 1.e-5; alfa<1.e-4; alfa+=1.e-5) {
        state.fAlpha = alfa;
        pSD->SetState(state);
        BuildPlasticSurface(cgrid, pSD);
        an.PostProcess(0);
    }
    delete cgrid;
    delete gmesh;
    
    
}

int ComputeMultiplier(TPZVec<REAL> &stress, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD, TPZVec<REAL> &stressresult);

void BuildPlasticSurface(TPZCompMesh *cmesh, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD)
{
    int ncon = cmesh->NConnects();
    TPZVec<int> computed(ncon,0);
    for (int el=0; el<cmesh->NElements(); el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        TPZManVector<REAL,3> centerksi(2,0.),xcenter(3,0.);
        gel->CenterPoint(gel->NSides()-1, centerksi);
        gel->X(centerksi, xcenter);
        TPZManVector<REAL> stress(3,0.);
        int matid = ComputeMultiplier(xcenter, pSD,stress);
        cmesh->ElementSolution()(cel->Index(),0) = matid;
        //        gel->SetMaterialId(matid);
        for (int icon=0; icon<gel->NCornerNodes(); icon++)
        {
            TPZConnect &c = cel->Connect(icon);
            TPZGeoNode &gnod = *gel->NodePtr(icon);
            TPZManVector<REAL> co(3,0.),stress(3,0.);
            gnod.GetCoordinates(co);
            
            ComputeMultiplier(co, pSD,stress);
            int seqnum = c.SequenceNumber();
            for (int idf=0; idf<3; idf++) {
                cmesh->Block()(seqnum,0,idf,0) = -co[idf]+stress[idf];
            }
        }
    }
}

int ComputeMultiplier(TPZVec<REAL> &stress, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD,TPZVec<REAL> &stressresult)
{
    REAL mult = 1.;
    REAL incr = 1.;
    TPZTensor<REAL> stresstensor,epsilon,stresscenter, epscenter;
    stresstensor.XX() = stress[0];
    stresstensor.YY() = stress[1];
    stresstensor.ZZ() = stress[2];
    stresscenter.XX() = 0.1;
    stresscenter.YY() = 0.1;
    stresscenter.ZZ() = 0.1;
    SANDLERDIMAGGIOSTEP1 *pSDP = dynamic_cast<SANDLERDIMAGGIOSTEP1 *>(pSD);
    pSDP->fER.ComputeDeformation(stresstensor,epsilon);
    pSDP->fER.ComputeDeformation(stresscenter, epscenter);
    TPZTensor<REAL> epsstart(epsilon);
    TPZManVector<REAL,3> phi(2,0.);
    pSD->Phi(epscenter, phi);
    pSD->Phi(epsilon, phi);
    while((phi[0]) > 0. || (phi[1]) > 0.)
    {
        mult *= 0.5;
        epsilon = epsstart;
        epsilon.Scale(mult);
        epsilon.Add(epscenter, 1.);
        pSD->Phi(epsilon, phi);
    }
    while((phi[0]) < 0. && (phi[1]) < 0.)
    {
        mult *= 2.;
        epsilon = epsstart;
        epsilon.Scale(mult);
        epsilon.Add(epscenter, 1.);
        pSD->Phi(epsilon, phi);
    }
    REAL tol = mult * 1.e-4;
    mult *= 0.5;
    incr = mult/2.;
    while (incr > tol) {
        mult += incr;
        epsilon = epsstart;
        epsilon.Scale(mult);
        epsilon.Add(epscenter, 1.);
        pSD->Phi(epsilon, phi);
        if ((phi[0]) > 0. || (phi[1]) > 0) {
            mult -= incr;
            incr /= 2.;
        }
    }
    epsilon = epsstart;
    epsilon.Scale(mult);
    epsilon.Add(epscenter, 1.);
    pSD->Phi(epsilon, phi);
    int result = 3;
    if (fabs(phi[0]) < fabs(phi[1])) {
        result = 1;
    }
    else {
        result = 2;
    }
    pSD->fER.Compute(epsilon, stresstensor);
    stressresult[0] = stresstensor.XX();
    stressresult[1] = stresstensor.YY();
    stressresult[2] = stresstensor.ZZ();
    return result;
}
