/**
 * \file
 * @brief Contains implementations of the TPZPlaneFracture methods.
 * @author Cesar Lucci
 * @since 09/08/2010
 */

#include "TPZPlaneFracture.h"

#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZVTKGeoMesh.h"
#include "tpzgeoelrefpattern.h"
#include "tpzchangeel.h"
#include "pzanalysis.h"
#include "pzelast3d.h"
#include "pzelasmat.h"
#include "pzbndcond.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzreducedspace.h"
#include "TPZElast3Dnlinear.h"


#include "adapt.h"
#include "TPZJIntegral.h"
#include "TPZReynoldsFlow.h"

using namespace pztopology;


/** PUBLIC METHODS */
TPZPlaneFracture::TPZPlaneFracture()
{
    std::cout << "Default constructor would not be used in this class!\n";
    DebugStop();
}
//------------------------------------------------------------------------------------------------------------


TPZPlaneFracture::TPZPlaneFracture(REAL lw, REAL bulletDepthTVDIni, REAL bulletDepthTVDFin,
                                   TPZVec< std::map<REAL,REAL> > & posTVD_stress, REAL xLength, REAL yLength, REAL Lmax)
{
    fposTVD_stress = posTVD_stress;
    
    fInitialElIndex = 0;
    
    fPreservedMesh = new TPZGeoMesh;
    
    fLmax = Lmax;
    
    std::set<REAL> espacamentoVerticalTVD;
    std::list<REAL> espacamentoVerticalDEPTH;
    
    std::map<REAL,REAL>::iterator itM;
    
    espacamentoVerticalTVD.insert(((REAL)0.));
    espacamentoVerticalTVD.insert(lw);
    espacamentoVerticalTVD.insert(bulletDepthTVDIni);
    espacamentoVerticalTVD.insert(bulletDepthTVDFin);
    
    int nstretches = posTVD_stress.NElements();
    for(int s = 0; s < nstretches; s++)
    {
        for(itM = posTVD_stress[s].begin(); itM != posTVD_stress[s].end(); itM++)
        {
            REAL pos = itM->first;
            espacamentoVerticalTVD.insert(pos);
        }
    }
    
    REAL pos0 = 0.;
    std::set<REAL>::iterator itS = espacamentoVerticalTVD.begin(); itS++;
    for(; itS != espacamentoVerticalTVD.end(); itS++)
    {
        REAL pos1 = *itS;
        REAL deltaZ = fabs(pos1 - pos0);
        
        int nrows = 1;
        REAL deltaZused = deltaZ/nrows;
        while(deltaZused > Lmax)
        {
            nrows++;
            deltaZused = deltaZ/nrows;
        }
        
        for(int r = 1; r <= nrows; r++)
        {
            REAL z = pos0 + r*deltaZused;
            espacamentoVerticalTVD.insert(z);
        }
        pos0 = pos1;
    }
    
    for(itS = espacamentoVerticalTVD.begin(); itS != espacamentoVerticalTVD.end(); itS++)
    {
        REAL posDepth = *itS;
        espacamentoVerticalDEPTH.push_back(-posDepth);//Converting positions (TVD) in depth.
    }
    
    GeneratePreservedMesh(espacamentoVerticalDEPTH, xLength, yLength);
}
//------------------------------------------------------------------------------------------------------------


TPZPlaneFracture::~TPZPlaneFracture()
{
    delete fPreservedMesh;
    fposTVD_stress.Resize(0);
}
//------------------------------------------------------------------------------------------------------------


//#define elastLinear
void TPZPlaneFracture::RunThisFractureGeometry(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                               REAL pressureInsideCrack,
                                               std::string vtkFile,
                                               bool printVTKfile)
{
    int porder = 1;
    
    //--------------------------------------------------------------------------------------------
    REAL pressureInsideCrackRef = 1.;//Pressao da malha computacional referenciada
    
#ifdef elastLinear
    TPZCompMesh * fractureCMeshRef = this->GetFractureCompMesh(poligonalChain, porder, pressureInsideCrackRef);
#else
    TPZCompMesh * fractureCMeshRef = this->GetFractureCompMeshNLinear(poligonalChain, porder, pressureInsideCrackRef);
#endif
    
    int neq = fractureCMeshRef->NEquations();
    std::cout << "Numero de equacoes = " << neq << std::endl;
    
    ////Analysis 1
    TPZAnalysis anRef(fractureCMeshRef);
    
    TPZSkylineStructMatrix skylinRef(fractureCMeshRef);
    TPZStepSolver<STATE> stepRef;
    stepRef.SetDirect(ECholesky);
    
    anRef.SetStructuralMatrix(skylinRef);
    anRef.SetSolver(stepRef);
    
#ifdef elastLinear
    anRef.Run();
#else
    {
        REAL chute = 10.;
        for(int r = 0; r < anRef.Solution().Rows(); r++)
            for(int c = 0; c < anRef.Solution().Cols(); c++)
                anRef.Solution()(r,c) = chute;
        anRef.Run();
        for(int r = 0; r < anRef.Solution().Rows(); r++)
            for(int c = 0; c < anRef.Solution().Cols(); c++)
                anRef.Solution()(r,c) += chute;
    }
#endif
    
    if(printVTKfile)
    {
        std::ofstream out0("0out.txt");
        out0 << "sol0={";
        int nelem = fractureCMeshRef->NElements();
        for(int el = 0; el < nelem; el++)
        {
            TPZCompEl * cel = fractureCMeshRef->ElementVec()[el];
            if(cel->Reference()->MaterialId() == 10)
            {
                TPZVec<REAL> centerQsi(3,0.);
				TPZVec<STATE> solTemp(3,0.);
                cel->Reference()->CenterPoint(cel->Reference()->NSides()-1, centerQsi);
                
                cel->Solution(centerQsi, 0, solTemp);
                for(int s = 0; s < solTemp.NElements(); s++)
                {
                    out0 << solTemp[s] << ",";
                }
            }
        }
        out0 << "};\n";
    }
    
//    if(printVTKfile)
//    {
//        TPZManVector<std::string,10> scalnames(3), vecnames(1);
//        
//        //        scalnames[0] = "EDisplacementX";
//        //        scalnames[1] = "EDisplacementY";
//        scalnames[0] = "StressX";
//        scalnames[1] = "StressY";
//        scalnames[2] = "StressZ";
//        
//        vecnames[0] = "Displacement";
//        
//        const int dim = 3;
//        int div = 0;
//        anRef.DefineGraphMesh(dim,scalnames,vecnames,vtkFile);
//        anRef.PostProcess(div);
//    }

    //--------------------------------------------------------------------------------------------
#ifdef elastLinear
    TPZCompMeshReferred * fractureCMesh = this->GetFractureCompMeshReferred(poligonalChain,porder, pressureInsideCrack, fractureCMeshRef);
#else
    TPZCompMeshReferred * fractureCMesh = this->GetFractureCompMeshReferredNLinear(poligonalChain,porder, pressureInsideCrack, fractureCMeshRef);
#endif
    
    ////Analysis 2
    TPZAnalysis an(fractureCMesh);
    
    TPZSkylineStructMatrix skylin(fractureCMesh); //caso simetrico
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    
    an.SetStructuralMatrix(skylin);
    an.SetSolver(step);
    
#ifdef elastLinear
    an.Run();
#else
    {
        REAL chute = 10.;
        for(int r = 0; r < an.Solution().Rows(); r++)
            for(int c = 0; c < an.Solution().Cols(); c++)
                an.Solution()(r,c) = chute;
        an.Run();
        for(int r = 0; r < an.Solution().Rows(); r++)
            for(int c = 0; c < an.Solution().Cols(); c++)
                an.Solution()(r,c) += chute;
    }
#endif
    
    if(printVTKfile)
    {
        std::ofstream out1("1out.txt");
        out1 << "sol1={";
        int nelem = fractureCMesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            TPZCompEl * cel = fractureCMesh->ElementVec()[el];
            if(cel->Reference()->MaterialId() == 10)
            {
                TPZVec<REAL> centerQsi(3,0.);
				TPZVec<STATE> solTemp(3,0.);
                cel->Reference()->CenterPoint(cel->Reference()->NSides()-1, centerQsi);
                
                cel->Solution(centerQsi, 0, solTemp);
                for(int s = 0; s < solTemp.NElements(); s++)
                {
                    out1 << solTemp[s] << ",";
                }
            }
        }
        out1 << "};\n";
    }
    
//    if(printVTKfile)
//    {
//        TPZManVector<std::string,10> scalnames(3), vecnames(1);
//    
//    //        scalnames[0] = "EDisplacementX";
//    //        scalnames[1] = "EDisplacementY";
//        scalnames[0] = "StressX";
//        scalnames[1] = "StressY";
//        scalnames[2] = "StressZ";
//        
//        vecnames[0] = "Displacement";
//
//        std::string vtkFile1 = "fracturePconstant1.vtk";
//        const int dim = 3;
//        int div = 0;
//        an.DefineGraphMesh(dim,scalnames,vecnames,vtkFile1);
//        an.PostProcess(div);
//    }
    
    
    /////// Example of J-Integral
    
    //    TPZVec<REAL> originXYZ(3,0.), direction(3,0.);
    //
    //    int POSmiddle1D = int(REAL(fcrackBoundaryElementsIds.NElements())/2. + 0.5);
    //    int middle1DId = fcrackBoundaryElementsIds[POSmiddle1D];
    //
    //    TPZVec<REAL> originQSI(1,0.);
    //    TPZGeoEl * gel1D = fractureCMesh->Reference()->ElementVec()[middle1DId];
    //    gel1D->X(originQSI, originXYZ);
    //
    //    direction[0] = 0.;
    //    direction[1] = 0.;
    //    direction[2] = 1.;
    //
    //    originXYZ[2] = -5.;
    //
    //    REAL radius = 0.6;
    //    Path3D * Path3DMiddle = new Path3D(fractureCMesh, originXYZ, direction, radius, pressureInsideCrack);
    //
    //    JIntegral3D jInt;
    //    jInt.PushBackPath3D(Path3DMiddle);
    //    TPZVec<REAL> Jvector(3);
    //
    //    Jvector = jInt.IntegratePath3D(0);
}

/** PRIVATE METHODS */
//------------------------------------------------------------------------------------------------------------


TPZCompMesh * TPZPlaneFracture::GetFractureCompMesh(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                    int porder, REAL pressureInsideCrack)
{
    ////GeoMesh
    TPZGeoMesh * gmesh = this->GetFractureGeoMesh(poligonalChain);
    
    ////CompMesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDimModel(3);
    
    STATE young = 0.29E5;
    STATE poisson = 0.25;
    TPZVec<STATE> force(3,0.);
    STATE prestressXX = 0.;
    STATE prestressYY = 0.;
    STATE prestressZZ = 0.;
    
    TPZMaterial * materialLin = new TPZElasticity3D(__3DrockMat_linear, young, poisson, force, prestressXX, prestressYY, prestressZZ);
    cmesh->InsertMaterialObject(materialLin);
    
    TPZMaterial * materialQpoint = new TPZElasticity3D(__3DrockMat_quarterPoint, young, poisson, force, prestressXX, prestressYY, prestressZZ);
    cmesh->InsertMaterialObject(materialQpoint);
    
    ////BCs
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, dirichDir = 3;
    //    int mixed = 2;
    
    {
        f.Zero();
        f(0,0) = 1.;
        TPZMaterial * materialMixedLeft = new TPZElasticity3D(-300, young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,__2DleftMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedLeft);
        
        f.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialMixedOutFracture = new TPZElasticity3D(-301, young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,__2DfractureMat_outside, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);
        
        f.Zero();
        f(2,0) = 1.;
        TPZMaterial * materialMixedTop = new TPZElasticity3D(-302, young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,__2DtopMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        //
        TPZMaterial * materialMixedBottom = new TPZElasticity3D(-303, young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,__2DbottomMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedBottom);
        
        ///////////farField
        k.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialNewmannFarField = new TPZElasticity3D(-304, young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,__2DfarfieldMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f.Zero();
        f(1,0) = pressureInsideCrack;
        TPZMaterial * materialNewmannInsideFract = new TPZElasticity3D(-305, young, poisson, force);
        TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,__2DfractureMat_inside, newmann, k, f);
        cmesh->InsertMaterialObject(newmannInsideFract);
    }
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    cmesh->SetDefaultOrder(porder);
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

TPZCompMesh * TPZPlaneFracture::GetFractureCompMeshNLinear(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                           int porder, REAL pressureInsideCrack)
{
    ////GeoMesh
    TPZGeoMesh * gmesh = this->GetFractureGeoMesh(poligonalChain);
    
    ////CompMesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDimModel(3);
    
    STATE young = 0.29E5;
    STATE poisson = 0.25;
    TPZVec<STATE> force(3,0.);
    STATE prestressXX = 0.;
    STATE prestressYY = 0.;
    STATE prestressZZ = 0.;
    
    TPZMaterial * materialLin = new TPZElast3Dnlinear(__3DrockMat_linear, young, poisson, force, prestressXX, prestressYY, prestressZZ);
    cmesh->InsertMaterialObject(materialLin);
    
    TPZMaterial * materialQpoint = new TPZElast3Dnlinear(__3DrockMat_quarterPoint, young, poisson, force, prestressXX, prestressYY, prestressZZ);
    cmesh->InsertMaterialObject(materialQpoint);
    
    ////BCs
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, dirichDir = 3;
    //    int mixed = 2;
    
    {
        f.Zero();
        f(0,0) = 1.;
        TPZMaterial * materialMixedLeft = new TPZElast3Dnlinear(-300, young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,__2DleftMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedLeft);
        
        f.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialMixedOutFracture = new TPZElast3Dnlinear(-301, young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,__2DfractureMat_outside, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);
        
        f.Zero();
        f(2,0) = 1.;
        TPZMaterial * materialMixedTop = new TPZElast3Dnlinear(-302, young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,__2DtopMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        //
        TPZMaterial * materialMixedBottom = new TPZElast3Dnlinear(-303, young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,__2DbottomMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedBottom);
        
        ///////////farField
        k.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialNewmannFarField = new TPZElast3Dnlinear(-304, young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,__2DfarfieldMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f.Zero();
        f(1,0) = pressureInsideCrack;
        TPZMaterial * materialNewmannInsideFract = new TPZElast3Dnlinear(-305, young, poisson, force);
        TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,__2DfractureMat_inside, newmann, k, f);
        cmesh->InsertMaterialObject(newmannInsideFract);
    }
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    cmesh->SetDefaultOrder(porder);
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}


TPZCompMeshReferred * TPZPlaneFracture::GetFractureCompMeshReferred(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                                    int porder, REAL pressureInsideCrack,
                                                                    TPZCompMesh * cmeshRef)
{
    ////GeoMesh
    TPZGeoMesh * gmesh = cmeshRef->Reference();
    gmesh->ResetReference();
    
    ////CompMeshReferred
    TPZCompMeshReferred * cmesh = new TPZCompMeshReferred(gmesh);
    
    cmesh->SetDimModel(3);
    
    STATE young = 0.29E5;
    STATE poisson = 0.25;
    TPZVec<STATE> force(3,0.);
    STATE prestressXX = 0.;
    STATE prestressYY = 0.;
    STATE prestressZZ = 0.;
    
    TPZMaterial * materialLin = new TPZElasticity3D(__3DrockMat_linear, young, poisson, force, prestressXX, prestressYY, prestressZZ);
    cmesh->InsertMaterialObject(materialLin);
    
    TPZMaterial * materialQpoint = new TPZElasticity3D(__3DrockMat_quarterPoint, young, poisson, force, prestressXX, prestressYY, prestressZZ);
    cmesh->InsertMaterialObject(materialQpoint);
    
    ////BCs
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, dirichDir = 3;
    //    int mixed = 2;
    
    {
        f.Zero();
        f(0,0) = 1.;
        TPZMaterial * materialMixedLeft = new TPZElasticity3D(-300, young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,__2DleftMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedLeft);
        
        f.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialMixedOutFracture = new TPZElasticity3D(-301, young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,__2DfractureMat_outside, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);
        
        f.Zero();
        f(2,0) = 1.;
        TPZMaterial * materialMixedTop = new TPZElasticity3D(-302, young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,__2DtopMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        //
        TPZMaterial * materialMixedBottom = new TPZElasticity3D(-303, young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,__2DbottomMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedBottom);
        
        ///////////farField
        k.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialNewmannFarField = new TPZElasticity3D(-304, young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,__2DfarfieldMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f.Zero();
        f(1,0) = pressureInsideCrack;
        TPZMaterial * materialNewmannInsideFract = new TPZElasticity3D(-305, young, poisson, force);
        TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,__2DfractureMat_inside, newmann, k, f);
        cmesh->InsertMaterialObject(newmannInsideFract);
    }
    
    int numsol = cmesh->Solution().Cols();
    cmesh->AllocateNewConnect(numsol, 1, 1);
    
    TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmesh);
    
    cmesh->SetDefaultOrder(porder);
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    cmesh->LoadReferred(cmeshRef);
    
    return cmesh;
}


TPZCompMeshReferred * TPZPlaneFracture::GetFractureCompMeshReferredNLinear(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                                           int porder, REAL pressureInsideCrack,
                                                                           TPZCompMesh * cmeshRef)
{
    ////GeoMesh
    TPZGeoMesh * gmesh = cmeshRef->Reference();
    gmesh->ResetReference();
    
    ////CompMeshReferred
    TPZCompMeshReferred * cmesh = new TPZCompMeshReferred(gmesh);
    
    cmesh->SetDimModel(3);
    
    STATE young = 0.29E5;
    STATE poisson = 0.25;
    TPZVec<STATE> force(3,0.);
    STATE prestressXX = 0.;
    STATE prestressYY = 0.;
    STATE prestressZZ = 0.;
    
    TPZMaterial * materialLin = new TPZElast3Dnlinear(__3DrockMat_linear, young, poisson, force, prestressXX, prestressYY, prestressZZ);
    cmesh->InsertMaterialObject(materialLin);
    
    TPZMaterial * materialQpoint = new TPZElast3Dnlinear(__3DrockMat_quarterPoint, young, poisson, force, prestressXX, prestressYY, prestressZZ);
    cmesh->InsertMaterialObject(materialQpoint);
    
    ////BCs
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, dirichDir = 3;
    //    int mixed = 2;
    
    {
        f.Zero();
        f(0,0) = 1.;
        TPZMaterial * materialMixedLeft = new TPZElast3Dnlinear(-300, young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,__2DleftMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedLeft);
        
        f.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialMixedOutFracture = new TPZElast3Dnlinear(-301, young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,__2DfractureMat_outside, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);
        
        f.Zero();
        f(2,0) = 1.;
        TPZMaterial * materialMixedTop = new TPZElast3Dnlinear(-302, young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,__2DtopMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        //
        TPZMaterial * materialMixedBottom = new TPZElast3Dnlinear(-303, young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,__2DbottomMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedBottom);
        
        ///////////farField
        k.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialNewmannFarField = new TPZElast3Dnlinear(-304, young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,__2DfarfieldMat, dirichDir, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f.Zero();
        f(1,0) = pressureInsideCrack;
        TPZMaterial * materialNewmannInsideFract = new TPZElast3Dnlinear(-305, young, poisson, force);
        TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,__2DfractureMat_inside, newmann, k, f);
        cmesh->InsertMaterialObject(newmannInsideFract);
    }
    
    int numsol = cmesh->Solution().Cols();
    cmesh->AllocateNewConnect(numsol, 1, 1);
    
    TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmesh);
    
    cmesh->SetDefaultOrder(porder);
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    cmesh->LoadReferred(cmeshRef);
    
    return cmesh;
}


TPZGeoMesh * TPZPlaneFracture::GetFractureGeoMesh(const TPZVec<std::pair<REAL,REAL> > &poligonalChain)
{
    TPZGeoMesh * refinedMesh = new TPZGeoMesh(*fPreservedMesh);
    
	long nelem = refinedMesh->NElements();
    
	std::map< long, std::set<REAL> > elId_TrimCoords;
	std::list< std::pair<long,REAL> > elIdSequence;
	
	DetectEdgesCrossed(poligonalChain, refinedMesh, elId_TrimCoords, elIdSequence);
	
	//Refining auxiliar 1D elements
	TPZVec<TPZGeoEl*> sons;
	std::map< long, std::set<REAL> >::iterator it;
	for(it = elId_TrimCoords.begin(); it != elId_TrimCoords.end(); it++)
	{
		int el1DId = it->first;
        TPZGeoEl * el1D = refinedMesh->ElementVec()[el1DId];
        
		TPZAutoPointer<TPZRefPattern> linRefp = Generate1DRefPatt(it->second);
        
		el1D->SetRefPattern(linRefp);
		el1D->Divide(sons);
	}
    
	//Refining 2D and 3D elements with the intention to match the geometry of the crack boundary
	for(long el = 0; el < nelem; el++)
	{
		TPZGeoEl * gel = refinedMesh->ElementVec()[el];//2D element in 2D mesh
        
		TPZAutoPointer<TPZRefPattern> elRefp = TPZRefPatternTools::PerfectMatchRefPattern(gel);
		if(elRefp)
		{
            gel = refinedMesh->ElementVec()[el];//2D element in 3D mesh
            if(gel)
            {
                gel->SetRefPattern(elRefp);
                gel->Divide(sons);
                
                int innerSide = gel->NSides() - 1;
                gel = gel->Neighbour(innerSide).Element();//3D element in 3D mesh
                
                elRefp = TPZRefPatternTools::PerfectMatchRefPattern(gel);
                if(elRefp)
                {
                    gel->SetRefPattern(elRefp);
                    gel->Divide(sons);
                }
                int firstFace = 20;
                int lastFace = 25;
                for(int f = firstFace; f < lastFace; f++)//2D BC element in 3D mesh
                {
                    TPZGeoElSide hexaFace(gel,f);
                    TPZGeoElSide quadriFace = hexaFace.Neighbour();
                    if( quadriFace != hexaFace && IsBoundaryMaterial(quadriFace.Element()) )
                    {
                        elRefp = TPZRefPatternTools::PerfectMatchRefPattern(quadriFace.Element());
                        if(elRefp)
                        {
                            quadriFace.Element()->SetRefPattern(elRefp);
                            quadriFace.Element()->Divide(sons);
                        }
                    }
                }
            }
            else
            {
                DebugStop();
            }
		}
	}
	
	GenerateCrackBoundary(refinedMesh, elIdSequence);
    SeparateElementsInMaterialSets(refinedMesh);
    TurnIntoQuarterPoint(refinedMesh);
    
    //    std::ofstream cc("cc.vtk");
    //    TPZVTKGeoMesh::PrintGMeshVTK(refinedMesh, cc, true);
    
	return refinedMesh;
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::GeneratePreservedMesh(std::list<REAL> & espacamento, REAL xLength, REAL yLength)
{
    TPZVec< TPZVec<REAL> > NodeCoord(0);
    long nrows, ncols;
    
    std::list<REAL>::iterator it = espacamento.end(); it--;
    
    int nDirRef = int(log(fabs(yLength)/fLmax)/log(2.));
    int nLayers = nDirRef + 2;
    
    REAL Y = 0.;
    GenerateNodesAtPlaneY(espacamento, xLength, NodeCoord, nrows, ncols, Y);
    
    Y = yLength/pow(2.,nDirRef);
    for(int lay = 1; lay < nLayers; lay++)
    {
        GenerateNodesAtPlaneY(espacamento, xLength, NodeCoord, nrows, ncols, Y);
        Y *= 2.;
    }
    long nNodesByLayer = nrows*ncols;
    long Qnodes = nNodesByLayer * nLayers;
	
	//initializing gmesh->NodeVec()
	fPreservedMesh->NodeVec().Resize(Qnodes);
    
    long pos = 0;
	TPZGeoNode Node;
    for(long l = 0; l < nLayers; l++)
    {
        for(long n = 0; n < nNodesByLayer; n++)
        {
            Node.SetNodeId(pos);
            Node.SetCoord(NodeCoord[pos]);
            fPreservedMesh->NodeVec()[pos] = Node;
            pos++;
        }
    }
	
	//inserting quadrilaterals
	long elId = 0;
	TPZVec <long> Topol(4);
	for(long r = 0; r < (nrows-1); r++)
	{
		for(long c = 0; c < (ncols-1); c++)
		{
			Topol[0] = ncols*r+c; Topol[1] = ncols*r+c+1; Topol[2] = ncols*(r+1)+c+1; Topol[3] = ncols*(r+1)+c;
			new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DfractureMat_outside,*fPreservedMesh);
			elId++;
		}
	}
    
    //inserting hexaedrons
    for(long l = 0; l < (nLayers-1); l++)
    {
        for(long r = 0; r < (nrows-1); r++)
        {
            for(long c = 0; c < (ncols-1); c++)
            {
               	Topol.Resize(8);
                Topol[0] = ncols*r+c + l*nNodesByLayer;
                Topol[1] = ncols*r+c+1 + l*nNodesByLayer;
                Topol[2] = ncols*(r+1)+c+1 + l*nNodesByLayer;
                Topol[3] = ncols*(r+1)+c + l*nNodesByLayer;
                //
                Topol[4] = ncols*r+c + (l+1)*nNodesByLayer;
                Topol[5] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                Topol[6] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                Topol[7] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (elId,Topol,__3DrockMat_linear,*fPreservedMesh);
                elId++;
                
                Topol.Resize(4);
                if(l == (nLayers-2))//farfield cc
                {
                    Topol[0] = ncols*r+c + (l+1)*nNodesByLayer;
                    Topol[1] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                    Topol[2] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DfarfieldMat,*fPreservedMesh);
                    elId++;
                }
                if(c == 0)//left cc
                {
                    Topol[0] = ncols*r+c + l*nNodesByLayer;
                    Topol[1] = ncols*r+c + (l+1)*nNodesByLayer;
                    Topol[2] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*(r+1)+c + l*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DleftMat,*fPreservedMesh);
                    elId++;
                }
                else if(c == (ncols - 2))//right cc
                {
                    Topol[0] = ncols*r+c+1 + l*nNodesByLayer;
                    Topol[1] = ncols*(r+1)+c+1 + l*nNodesByLayer;
                    Topol[2] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DrightMat,*fPreservedMesh);
                    elId++;
                }
                if(r == 0)//top cc
                {
                    Topol[0] = ncols*r+c + l*nNodesByLayer;
                    Topol[1] = ncols*r+c+1 + l*nNodesByLayer;
                    Topol[2] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*r+c + (l+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DtopMat,*fPreservedMesh);
                    elId++;
                }
                else if(r == (nrows - 2))//bottom cc
                {
                    Topol[0] = ncols*(r+1)+c + l*nNodesByLayer;
                    Topol[1] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                    Topol[2] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*(r+1)+c+1 + l*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DbottomMat,*fPreservedMesh);
                    elId++;
                }
            }
        }
    }
    
	fPreservedMesh->BuildConnectivity();
    fPreservedMesh->SetMaxElementId(fPreservedMesh->NElements()-1);
    fPreservedMesh->SetMaxNodeId(fPreservedMesh->NNodes()-1);
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::GenerateNodesAtPlaneY(std::list<REAL> & espacamento, REAL xLength,
                                             TPZVec< TPZVec<REAL> > & NodeCoord, long & nrows, long & ncols,
                                             REAL Y)
{
    nrows = espacamento.size();
    
    std::list<REAL>::iterator it = espacamento.end(); it--;
    
    int nStretches = 1;
    REAL deltaX = xLength/nStretches;
    while(deltaX > fLmax)
    {
        nStretches++;
        deltaX = xLength/nStretches;
    }
    ncols = nStretches+1;
    
	const long nNodesByLayer = nrows*ncols;
    long oldNNodes = NodeCoord.NElements();
    long newNNodes = oldNNodes + nNodesByLayer;
    
	NodeCoord.Resize(newNNodes);
	
	long nodeId = oldNNodes;
	
	//setting nodes coords
	for(it = espacamento.begin(); it != espacamento.end(); it++)
	{
		for(long c = 0; c < ncols; c++)
		{
            NodeCoord[nodeId].Resize(3,0.);
			NodeCoord[nodeId][0] = c*deltaX;
            NodeCoord[nodeId][1] = Y;
			NodeCoord[nodeId][2] = *it;
			nodeId++;
		}
	}
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::DetectEdgesCrossed(const TPZVec<std::pair<REAL,REAL> > &poligonalChain, TPZGeoMesh * planeMesh,
                                          std::map< long, std::set<REAL> > &elId_TrimCoords,
                                          std::list< std::pair<long,REAL> > &elIdSequence)
{
	long npoints = (poligonalChain.NElements());
	long nelem = planeMesh->NElements();
	
    TPZVec<REAL> coord(3,0.);
    int firstPt = 0;
    coord[0] = poligonalChain[firstPt].first;
    coord[2] = poligonalChain[firstPt].second;
    
    int axe0 = 0;//axe X
    int axe1 = 2;//axe Z
    int axeNormal = 1;//axe Y
    TPZVec<REAL> qsi(2,0.);
    TPZGeoEl * firstGel = planeMesh->FindElement(coord, qsi, fInitialElIndex, 2);
    
	TPZGeoEl * gel = firstGel;
	TPZGeoEl * nextGel = NULL;
	
#ifdef DEBUG
	if(!gel)
	{
		std::cout << "first point of crack tip boundary does NOT belong to any 2D element" << std::endl;
		std::cout << "(or was given a mesh with zero quantity of 2D elements)!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
#endif
	
	REAL alphaMin;
	bool reachNextPoint;
	long nElsCrossed, thispoint, nextpoint;
	std::map< int, std::set<REAL> >::iterator it;
	std::set<REAL> trim;
	TPZVec< TPZVec<REAL> > intersectionPoint;
	TPZManVector<REAL,3> x(3,0.), dx(3,0.);
    TPZManVector<REAL,3> xNext(3,0.), qsi2D(2,0.);
	TPZVec <long> Topol(2), edgeVec;
	
	long p;
	for(p = 0; p < (npoints-1); p++)
	{
		nElsCrossed = 0;
		alphaMin = 0.;
		thispoint = p;
		nextpoint = p+1;
        
        x[0] = poligonalChain[thispoint].first;
        x[2] = poligonalChain[thispoint].second;
        
        dx[0] = poligonalChain[nextpoint].first - poligonalChain[thispoint].first;
        dx[2] = poligonalChain[nextpoint].second - poligonalChain[thispoint].second;
        
        xNext[0] = poligonalChain[nextpoint].first;
        xNext[2] = poligonalChain[nextpoint].second;
        
		REAL norm = 0.;
		for(int c = 0; c < 3; c++)
		{
			norm += dx[c]*dx[c];
		}
		norm = sqrt(norm);
		for(int c = 0; c < 3; c++)
		{
			dx[c] = dx[c]/norm;
		}
		
		REAL Tol;
		ZeroTolerance(Tol);
		reachNextPoint = gel->ComputeXInverse(xNext, qsi2D,Tol);
        int axe0 = 0;//axe X
        int axe1 = 2;//axe Z
        int axeNormal = 1;//axe Y
		while(reachNextPoint == false && nElsCrossed < nelem)
		{
			nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, true, axe0, axe1, axeNormal, false);
			
			alphaMin = __smallNum;//qdo vai de vizinho a vizinho ateh chegar no proximo ponto, nao deve-se incluir os (alphaX=0)
			if(nextGel->ComputeXInverse(xNext, qsi2D,Tol))
			{
				reachNextPoint = true;
			}
			gel = nextGel;
			nElsCrossed++;
		}
		if(nElsCrossed == nelem)
		{
			//Deve ter alternado entre vizinhos!!!
			DebugStop();
		}
	}
    
    TPZGeoEl * lastGel = gel;
	
	dx[0] = -1.;//direcao oposta ao eixo x do sistema de coordenadas
    dx[1] =  0.;//direcao oposta ao eixo x do sistema de coordenadas
    dx[2] =  0.;//direcao oposta ao eixo x do sistema de coordenadas
	
	p = 0;//Fechando o inicio da fratura
	thispoint = p;
	nextpoint = p+1;
    
    x[0] = poligonalChain[thispoint].first;
    x[1] = 0.;
    x[2] = poligonalChain[thispoint].second;
    
	nextGel = firstGel;
	gel = NULL;
	alphaMin = 0.;
	while(gel != nextGel)
	{
		gel = nextGel;
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, false, axe0, axe1, axeNormal, true);
        if(!nextGel)
        {
            DebugStop();
        }
		alphaMin = __smallNum;
	}
	
	p = npoints;//Fechando o final da fratura
	thispoint = p-1;
    
    x[0] = poligonalChain[thispoint].first;
    x[1] = 0.;
    x[2] = poligonalChain[thispoint].second;
    
	nextGel = lastGel;
	gel = NULL;
	alphaMin = 0.;
	while(gel != nextGel)
	{
		gel = nextGel;
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, true, axe0, axe1, axeNormal, true);
        if(!nextGel)
        {
            DebugStop();
        }
		alphaMin = __smallNum;
	}
}
//------------------------------------------------------------------------------------------------------------


TPZGeoEl * TPZPlaneFracture::CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, REAL alphaMin,
												  std::map< long, std::set<REAL> > &elId_TrimCoords,
                                                  std::list< std::pair<long,REAL> > &elIdSequence,
                                                  bool pushback, int planeAxe0, int planeAxe1, int planeNormal,
                                                  bool closingFracture)
{
	bool thereIsAn1DElemAlready;
	int edge;
	std::map< long, std::set<REAL> >::iterator it;
	std::set<REAL> trim;
	TPZVec< TPZVec<REAL> > ExactIntersectionPoint, ModulatedIntersectionPoint;
	TPZVec<REAL> qsi1Dvec(1), xCrackBoundary(3);
	TPZVec<long> Topol(2);
	TPZVec<int>  edgeVec;
	REAL Tol;
	ZeroTolerance(Tol);
	
	TPZGeoMesh * planeMesh = gel->Mesh();
	
	bool haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint, alphaMin, planeAxe0, planeAxe1, planeNormal);
	
	if(haveIntersection == false)
	{
		haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint, 0., planeAxe0, planeAxe1, planeNormal);
        
#ifdef DEBUG
		if(haveIntersection == false)
		{
            TPZManVector<REAL,3> qsi2D(2,0.);
            if(gel->ComputeXInverse(x, qsi2D,Tol))
            {
                std::cout << "The point inside element does NOT intersect its edges! EdgeIntersection method face an exeption!" << std::endl;
                std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
                
                std::cout << "Element " << gel->Id() << std::endl;
                for(long n = 0; n < gel->NNodes(); n++)
                {
                    std::cout << "Node " << n << std::endl;
                    gel->NodePtr(n)->Print(std::cout);
                    std::cout << std::endl;
                }
                std::cout << "x: " << x[0] << " , " << x[1] << " , " << x[2] << std::endl;
                std::cout << "dx: " << dx[0] << " , " << dx[1] << " , " << dx[2] << std::endl;
                std::cout << "alphaMin = " << alphaMin << std::endl << std::endl;
                
                DebugStop();
            }
            else
            {
                std::cout << "Trying to find edge intersection from an point that doesnt belong to given element!" << std::endl;
                std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
                
                DebugStop();
            }
		}
#endif
	}
	TPZVec<REAL> xLin(3);
	for(long edg = 0; edg < edgeVec.NElements(); edg++)
	{
		edge = edgeVec[edg];
		TPZVec<REAL> n0(3), n1(3);
		planeMesh->NodeVec()[gel->SideNodeIndex(edge, 0)].GetCoordinates(n0);
		planeMesh->NodeVec()[gel->SideNodeIndex(edge, 1)].GetCoordinates(n1);
		for(int c = 0; c < 3; c++)
		{
            REAL coordM = ModulatedIntersectionPoint[edg][c];
			xLin[c] = coordM;
		}
		REAL qsi1D = LinearComputeXInverse(xLin, n0, n1);
		
		TPZGeoElSide gelEdge(gel, edge);
		TPZGeoElSide neighEdge = gelEdge.Neighbour();
		thereIsAn1DElemAlready = false;
		while(neighEdge != gelEdge)
		{
			if(neighEdge.Element()->Dimension() == 1)//jah existe um elemento 1D inserido nesta aresta!
			{
				thereIsAn1DElemAlready = true;
				long neighEdgeId = neighEdge.Element()->Id();
				it = elId_TrimCoords.find(neighEdgeId);
				TPZTransform transBetweenNeigh = neighEdge.NeighbourSideTransform(gelEdge);
				qsi1D *= transBetweenNeigh.Mult()(0,0);
				it->second.insert(qsi1D);
				if(pushback)
				{
					elIdSequence.push_back(std::make_pair(neighEdgeId, qsi1D));
				}
				else // push_FRONT
				{
					elIdSequence.push_front(std::make_pair(neighEdgeId, qsi1D));
				}
                
				break;
			}
			neighEdge = neighEdge.Neighbour();
		}
		if(thereIsAn1DElemAlready == false)//nao existe um elemento 1D nesta aresta!
		{
			trim.clear();
			trim.insert(qsi1D);
			Topol[0] = gel->SideNodeIndex(edge, 0);
			Topol[1] = gel->SideNodeIndex(edge, 1);
			TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * linGeo =
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, __aux1DEl_Mat, *planeMesh);
            
			planeMesh->BuildConnectivity();
			
			long linGeoId = linGeo->Id();
			elId_TrimCoords[linGeoId] = trim;
			
			if(pushback) // push_BACK
			{
				elIdSequence.push_back(std::make_pair(linGeoId, qsi1D));
			}
			else // push_FRONT
			{
				elIdSequence.push_front(std::make_pair(linGeoId, qsi1D));
			}
		}
	}
	x = ExactIntersectionPoint[ExactIntersectionPoint.NElements() - 1];
	TPZGeoElSide gelEdge(gel, edge);
	TPZGeoElSide neighEdge = gelEdge.Neighbour();
	bool neighFound = false;
	while(neighEdge != gelEdge)
	{
		if(neighEdge.Element()->Dimension() == 2 &&
           neighEdge.Element()->MaterialId() == __2DfractureMat_outside &&
           neighEdge.Element()->Father() == NULL)
		{
            neighFound = true;
			break;
		}
		neighEdge = neighEdge.Neighbour();
	}
    if(neighFound)
    {
        return neighEdge.Element();
    }
    else if(closingFracture == false)
    {
        DebugStop();//fratura saiu para fora do dominio
    }
    else if(gelEdge.Element()->Dimension() == 2 &&
            gelEdge.Element()->MaterialId() == __2DfractureMat_outside &&
            gelEdge.Element()->Father() == NULL)
    {
        //fechando fratura no inicio ou no final, portanto o
        //ultimo elemento nao encontra vizinho do tipo 2d com
        //materialId == __2DfractureMat_outside, retornando a si mesmo
        return gelEdge.Element();
    }
	
	return NULL;
}
//------------------------------------------------------------------------------------------------------------


bool TPZPlaneFracture::EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
                                        TPZVec< TPZVec<REAL> > &ExactIntersect, REAL alphaMin,
                                        int planeAxe0, int planeAxe1, int planeNormal)
{
    long nearNode;
	bool IsNearNode = TPZChangeEl::NearestNode(gel, x, nearNode, __smallNum);
	
	edge.Resize(0);
	ExactIntersect.Resize(0);
	
	int ncnodes = gel->NCornerNodes();
	
	TPZVec< TPZVec<REAL> > node(ncnodes), dnode(ncnodes);
	TPZVec<REAL> nodeCoord(3);
	for(int n = 0; n < ncnodes; n++)
	{
		node[n].Resize(3,1);
		dnode[n].Resize(3,1);
		gel->NodePtr(n)->GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			(node[n])[c] = nodeCoord[c];
		}
	}
	std::list<REAL> edgeNorm;
	/** <REAL: alpha , <REAL: alphaNodemod, REAL: alphaNodesmooth, REAL: norm, int: edge that intersect> > */
    std::map<REAL, TPZVec<REAL> > alpha;
    
	REAL alphaX, alphaNodemod, alphaNodesmooth, norm;
	for(int n = 0; n < ncnodes; n++)
	{
		norm = 0.;
		for(int c = 0; c < 3; c++)//computing the direction from node N to node N+1
		{
			if(n != (ncnodes-1) )
			{
				(dnode[n])[c] = node[n+1][c] - node[n][c];
			}
			else
			{
				dnode[n][c] = node[0][c] - node[n][c];
			}
			norm += dnode[n][c]*dnode[n][c];
		}
		norm = sqrt(norm);
		for(int c = 0; c < 3; c++)//normalizing computed direction
		{
			dnode[n][c] /= norm;
		}
		
		alphaX = ComputeAlphaX(x, dx, node[n], dnode[n], planeAxe0, planeAxe1, planeNormal);
		alphaNodemod = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, true, false, planeAxe0, planeAxe1, planeNormal);
		alphaNodesmooth = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, IsNearNode, true, planeAxe0, planeAxe1, planeNormal);
		if(alphaX >= alphaMin && alphaNodesmooth >= 0. && alphaNodesmooth <= norm)
		{
			int thisEdge = n+ncnodes;
			TPZVec<REAL> someData(4);
			someData[0] = alphaNodemod; someData[1] = alphaNodesmooth; someData[2] = norm; someData[3] = thisEdge;
			alpha[alphaX] = someData;
		}
	}
	
	// a aresta que serah interseccionada serah a que tiver menor alphaX não negativo, i.e.: o primeiro par do mapa alpha!
	std::map<REAL, TPZVec<REAL> >::iterator it = alpha.begin();
	if(it != alpha.end())
	{
		alphaX = it->first;
		alphaNodemod = (it->second)[0];
		alphaNodesmooth = (it->second)[1];
		norm = (it->second)[2];
		edge.Resize(1); edge[0] = int((it->second)[3]);
		
		ExactIntersect.Resize(1); ExactIntersect[0].Resize(3,1);
		for(int c = 0; c < 3; c++)
		{
			ExactIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodesmooth*dnode[edge[0]-ncnodes][c];
		}
		
		if(alphaX <= alphaMin && alpha.size() > 1)
		{
			it++;
			alphaX = it->first;
			alphaNodemod = (it->second)[0];
			alphaNodesmooth = (it->second)[1];
			norm = (it->second)[2];
			edge.Resize(2); edge[1] = int((it->second)[3]);
			
			ExactIntersect.Resize(2); ExactIntersect[1].Resize(3,1);
			for(int c = 0; c < 3; c++)
			{
				ExactIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodesmooth*dnode[edge[1]-ncnodes][c];
			}
		}
		return true;
	}
	
	//se o ponto p estah sobre um noh e a direcao dp aponta para fora do elemento...
	if(IsNearNode)
	{
		std::cout << "Estah no noh " << nearNode << " e nao foi encontrada interseccao!" << std::endl;
		std::cout << "Tratar este caso!" << std::endl << std::endl;
	}
	
	return false;
}
//------------------------------------------------------------------------------------------------------------


bool TPZPlaneFracture::EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
										TPZVec< TPZVec<REAL> > &ExactIntersect,
                                        TPZVec< TPZVec<REAL> > &ModulatedIntersect, REAL alphaMin,
                                        int planeAxe0, int planeAxe1, int planeNormal)
{
	long nearNode;
	bool IsNearNode = TPZChangeEl::NearestNode(gel, x, nearNode, __smallNum);
	
	edge.Resize(0);
	ExactIntersect.Resize(0);
	ModulatedIntersect.Resize(0);
	
	int ncnodes = gel->NCornerNodes();
	
	TPZVec< TPZVec<REAL> > node(ncnodes), dnode(ncnodes);
	TPZVec<REAL> nodeCoord(3);
	for(int n = 0; n < ncnodes; n++)
	{
		node[n].Resize(3,1);
		dnode[n].Resize(3,1);
		gel->NodePtr(n)->GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			(node[n])[c] = nodeCoord[c];
		}
	}
	std::list<REAL> edgeNorm;
	/** <REAL: alpha , <REAL: alphaNodemod, REAL: alphaNodesmooth, REAL: norm, int: edge that intersect> > */
    std::map<REAL, TPZVec<REAL> > alpha;
    
	REAL alphaX, alphaNodemod, alphaNodesmooth, norm;
	for(int n = 0; n < ncnodes; n++)
	{
		norm = 0.;
		for(int c = 0; c < 3; c++)//computing the direction from node N to node N+1
		{
			if(n != (ncnodes-1) )
			{
				(dnode[n])[c] = node[n+1][c] - node[n][c];
			}
			else
			{
				dnode[n][c] = node[0][c] - node[n][c];
			}
			norm += dnode[n][c]*dnode[n][c];
		}
		norm = sqrt(norm);
		for(int c = 0; c < 3; c++)//normalizing computed direction
		{
			dnode[n][c] /= norm;
		}
		
		alphaX = ComputeAlphaX(x, dx, node[n], dnode[n], planeAxe0, planeAxe1, planeNormal);
		alphaNodemod = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, true, false, planeAxe0, planeAxe1, planeNormal);
		alphaNodesmooth = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, IsNearNode, true, planeAxe0, planeAxe1, planeNormal);
		if(alphaX >= alphaMin && alphaNodesmooth >= 0. && alphaNodesmooth <= norm)
		{
			int thisEdge = n+ncnodes;
			TPZVec<REAL> someData(4);
			someData[0] = alphaNodemod; someData[1] = alphaNodesmooth; someData[2] = norm; someData[3] = thisEdge;
			alpha[alphaX] = someData;
		}
	}
	
	// a aresta que serah interseccionada serah a que tiver menor alphaX não negativo, i.e.: o primeiro par do mapa alpha!
	std::map<REAL, TPZVec<REAL> >::iterator it = alpha.begin();
	if(it != alpha.end())
	{
		alphaX = it->first;
		alphaNodemod = (it->second)[0];
		alphaNodesmooth = (it->second)[1];
		norm = (it->second)[2];
		edge.Resize(1); edge[0] = int((it->second)[3]);
		
		ModulatedIntersect.Resize(1); ModulatedIntersect[0].Resize(3,1);
		for(int c = 0; c < 3; c++)
		{
			ModulatedIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodemod*dnode[edge[0]-ncnodes][c];
		}
		ExactIntersect.Resize(1); ExactIntersect[0].Resize(3,1);
		for(int c = 0; c < 3; c++)
		{
			ExactIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodesmooth*dnode[edge[0]-ncnodes][c];
		}
		
		if(alphaX <= alphaMin && alpha.size() > 1)
		{
			it++;
			alphaX = it->first;
			alphaNodemod = (it->second)[0];
			alphaNodesmooth = (it->second)[1];
			norm = (it->second)[2];
			edge.Resize(2); edge[1] = int((it->second)[3]);
			
			ModulatedIntersect.Resize(2); ModulatedIntersect[1].Resize(3,1);
			for(int c = 0; c < 3; c++)
			{
				ModulatedIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodemod*dnode[edge[1]-ncnodes][c];
			}
			ExactIntersect.Resize(2); ExactIntersect[1].Resize(3,1);
			for(int c = 0; c < 3; c++)
			{
				ExactIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodesmooth*dnode[edge[1]-ncnodes][c];
			}
		}
		return true;
	}
	
	//se o ponto p estah sobre um noh e a direcao dp aponta para fora do elemento...
	if(IsNearNode)
	{
		std::cout << "Estah no noh " << nearNode << " e nao foi encontrada interseccao!" << std::endl;
		std::cout << "Tratar este caso!" << std::endl << std::endl;
	}
	
	return false;
}
//------------------------------------------------------------------------------------------------------------


REAL TPZPlaneFracture::ComputeAlphaNode(TPZVec<REAL> &x, TPZVec<REAL> &dx,
                                        TPZVec<REAL> &node, TPZVec<REAL> &dnode,
                                        REAL norm, bool modulate, bool smooth,
                                        int planeAxe0, int planeAxe1, int planeNormal)
{
	REAL fractionNumQ =	dx[planeAxe1]*node[planeAxe0] - dx[planeAxe0]*node[planeAxe1] - dx[planeAxe1]*x[planeAxe0] + dx[planeAxe0]*x[planeAxe1];
	
	REAL fractionDenomQ =	dnode[planeAxe1]*dx[planeAxe0] - dnode[planeAxe0]*dx[planeAxe1];
	
	REAL alphaNode = -1.;
	if(fabs(fractionDenomQ) > __smallNum)
	{
		alphaNode = fractionNumQ/fractionDenomQ;
		if(fabs(alphaNode) < __smallNum) alphaNode = 0.;
		if(fabs(alphaNode - norm) < __smallNum) alphaNode = norm;
	}
	
	if(modulate)
	{
		int stretchesQTD = __EdgeStretchesQTD;
		if(smooth)
		{
			stretchesQTD *= __TrimQTDmultiplier;
		}
		int nsegm = int(alphaNode/(norm/stretchesQTD) + 0.5);
		if(nsegm == 0)//Tirando a interseccao do noh inicial da aresta
		{
			nsegm = 1;
		}
		else if(nsegm == stretchesQTD)//Tirando a interseccao do noh final da aresta
		{
			nsegm = stretchesQTD - 1;
		}
		alphaNode = nsegm*(norm/stretchesQTD);//modulando o ponto para multiplo de (norm/__EdgeStretchesQTD)
	}
	
	return alphaNode;
}
//------------------------------------------------------------------------------------------------------------


REAL TPZPlaneFracture::ComputeAlphaX(TPZVec<REAL> &x, TPZVec<REAL> &dx,
                                     TPZVec<REAL> &node, TPZVec<REAL> &dnode,
                                     int planeAxe0, int planeAxe1, int planeNormal)
{
	//computing alpha (dx vector multiplier to intersect element edges)
	REAL fractionNumP   =	dnode[planeAxe1]*node[planeAxe0] - dnode[planeAxe0]*node[planeAxe1] -
	dnode[planeAxe1]*x[planeAxe0] + dnode[planeAxe0]*x[planeAxe1];
	
	REAL fractionDenomP =	dnode[planeAxe1]*dx[planeAxe0] - dnode[planeAxe0]*dx[planeAxe1];
	
	REAL alphaX = -1.;
	if(fabs(fractionDenomP) > __smallNum)
	{
		alphaX = fractionNumP/fractionDenomP;
		if(fabs(alphaX) < __smallNum)
		{
			alphaX = 0.;
		}
	}
	
	return alphaX;
}
//------------------------------------------------------------------------------------------------------------


REAL TPZPlaneFracture::LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1)
{
	REAL dL = 0., L = 0.;
	for(int c = 0; c < 3; c++)
	{
		L += (n1[c]-n0[c])*(n1[c]-n0[c]);
		dL += (x[c]-n0[c])*(x[c]-n0[c]);
	}
	L = sqrt(L);
	dL = sqrt(dL);
	
#ifdef DEBUG
	if(fabs(L) < __smallNum || fabs(dL) < __smallNum || fabs(L - dL) < __smallNum)
	{
		std::cout << "n0 and n1 are coincident nodes!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
#endif
	
	REAL qsi = -1. + dL/L*2.;
	
	return qsi;
}
//------------------------------------------------------------------------------------------------------------


TPZAutoPointer<TPZRefPattern> TPZPlaneFracture::Generate1DRefPatt(std::set<REAL> &TrimCoord)
{
	if(TrimCoord.size() == 0)
	{
		return NULL;
	}
	
	long Qnodes = TrimCoord.size() + 2;
	long Qelements = TrimCoord.size() + 2;
	long QsubElements = Qelements - 1;
	TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
	
	//1. setting initial and final nodes coordinates of 1D element mesh
	NodeCoord[0][0] = -1.;
	NodeCoord[1][0] =  1.;
	//.
	
	//2. setting intermediate nodes coordinates of 1D element mesh
	int c = 2;
	std::set<REAL>::iterator it;
	for(it = TrimCoord.begin(); it != TrimCoord.end(); it++)
	{
		REAL coord = *it;
		
		(NodeCoord[c])[0] = coord;
		c++;
	}
	//.
	
	//3. initializing internal mesh of refPattern
	TPZGeoMesh internalMesh;
	internalMesh.SetMaxNodeId(Qnodes-1);
	internalMesh.SetMaxElementId(Qelements-1);
	internalMesh.NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(NodeCoord[n]);
		internalMesh.NodeVec()[n] = Node[n];
	}
	//.
	
	//4. inserting 1D elements on internal mesh of refPattern
	long elId = 0;
	TPZVec <long> Topol(2);
	
	//4.1 inserting father element
	Topol[0] = 0; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__aux1DEl_Mat,internalMesh);
	elId++;
	//.
	
	//4.2 inserting subelements
	//first subelement
	Topol[0] = 0; Topol[1] = 2;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son1 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__aux1DEl_Mat,internalMesh);
	son1->SetFather(father);
	son1->SetFather(father->Index());
	elId++;
	//
	
	//last subelement
	Topol[0] = TrimCoord.size() + 1; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son2 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__aux1DEl_Mat,internalMesh);
	son2->SetFather(father);
	son2->SetFather(father->Index());
	elId++;
	//
	
	for(long el = 2; el < QsubElements; el++)
	{
		Topol[0] = el; Topol[1] = el+1;
		TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__aux1DEl_Mat,internalMesh);
		son->SetFather(father);
		son->SetFather(father->Index());
		elId++;
	}
	//.
	//.
	
	internalMesh.BuildConnectivity();
	
	TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(internalMesh);
	TPZAutoPointer<TPZRefPattern> Found = gRefDBase.FindRefPattern(refPattern);
	if(!Found)
	{
		gRefDBase.InsertRefPattern(refPattern);
		refPattern->InsertPermuted();
		
		return refPattern;
	}
	else
	{
		return Found;
	}
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::UpdatePoligonalChain(TPZGeoMesh * gmesh,
                                            std::list< std::pair<long,REAL> > &elIdSequence,
                                            TPZVec<std::pair<REAL,REAL> > &poligonalChainUpdated)
{
	int nptos = elIdSequence.size();
	poligonalChainUpdated.Resize(nptos);
    
	TPZVec<REAL> qsi1Dvec(1), ptoCoord(3);
	int el1Did, p = 0;
	REAL qsi1D;
	
	std::list< std::pair<long,REAL> >::iterator it;
	for(it = elIdSequence.begin(); it != elIdSequence.end(); it++)
	{
		el1Did = it->first;
		TPZGeoEl * el1D = gmesh->ElementVec()[el1Did];
		
#ifdef DEBUG
		int elDim = el1D->Dimension();
		if(elDim != 1)
		{
			std::cout << "The elIdSequence supposedly would contains ids of elements exclusively 1D!" << std::endl;
			std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
			DebugStop();
		}
#endif
		
		qsi1D = it->second;
		qsi1Dvec[0] = qsi1D;
		
		el1D->X(qsi1Dvec, ptoCoord);
		
		poligonalChainUpdated[p] = std::make_pair(ptoCoord[0],ptoCoord[2]);
		
		p++;
	}
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::GenerateCrackBoundary(TPZGeoMesh * refinedMesh,
                                             std::list< std::pair<long,REAL> > &elIdSequence)
{
    fcrackBoundaryElementsIds.Resize(0);
	TPZVec<REAL> qsi0vec(1), qsi1vec(1), node0coord(3), node1coord(3);
	TPZVec<long> Topol(2);
	long el0id, el1id, n0, n1;
	REAL qsi0, qsi1;
	std::list< std::pair<long,REAL> >::iterator crackit0, crackit1, crackitEnd;
	crackitEnd = elIdSequence.end(); crackitEnd--;
	for(crackit0 = elIdSequence.begin(); crackit0 != crackitEnd; crackit0++)
	{
		crackit1 = crackit0; crackit1++;
		
		el0id = crackit0->first;
		TPZGeoEl * el0 = refinedMesh->ElementVec()[el0id];
		qsi0 = crackit0->second;
		qsi0vec[0] = qsi0;
		el0->X(qsi0vec, node0coord);
		n0 = TPZChangeEl::NearestNode(refinedMesh, node0coord, __smallNum);
		Topol[0] = n0;
		
		el1id = crackit1->first;
		TPZGeoEl * el1 = refinedMesh->ElementVec()[el1id];
		qsi1 = crackit1->second;
		qsi1vec[0] = qsi1;
		el1->X(qsi1vec, node1coord);
		n1 = TPZChangeEl::NearestNode(refinedMesh, node1coord, __smallNum);
		Topol[1] = n1;
		
		TPZGeoEl * crack1D = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, __1DcrackTipMat, *refinedMesh);
        int oldSize = fcrackBoundaryElementsIds.NElements();
        fcrackBoundaryElementsIds.Resize(oldSize+1);
        fcrackBoundaryElementsIds[oldSize] = crack1D->Id();
	}
    
    refinedMesh->BuildConnectivity();
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::SeparateElementsInMaterialSets(TPZGeoMesh * refinedMesh)
{
    int n1Dels = fcrackBoundaryElementsIds.NElements();
    std::map<int,TPZFracture2DEl> fracturedElems;
    
    //Capturando subelementos que encostam no contorno da fratura
    for(int el = 0; el < n1Dels; el++)
    {
        int id = fcrackBoundaryElementsIds[el];
        TPZGeoEl * gel = refinedMesh->ElementVec()[id];//1D element of crach boundary
        
#ifdef DEBUG
        if(gel->Dimension() != 1)
        {
            DebugStop();
        }
#endif
        
        for(int sd = 0; sd < gel->NSides(); sd++)
        {
            TPZGeoElSide side1D(gel,sd);
            TPZGeoElSide sideNeigh = side1D.Neighbour();
            while(sideNeigh != side1D)
            {
                int sideNeighbyside = sideNeigh.Side();
                TPZGeoEl * neigh = sideNeigh.Element();
                
                if(neigh->HasSubElement() || neigh->Dimension() == 1)
                {
                    sideNeigh = sideNeigh.Neighbour();
                    continue;
                }
                
                //Primeiro estagio da identificacao dos elementos 2D que estao no interior da fratura.
                //brief: primeiramente sao identificados os elementos no interior da fratura que encostam no crack tip.
                //Obs.: A condicao (IsBoundaryMaterial(neigh) == false) eh para excluir os quadrilateros de condicao de contorno.
                if(sd == 2 && neigh->Dimension() == 2 && IsBoundaryMaterial(neigh) == false)
                {
                    TPZVec<REAL> neighCenterQSI(neigh->Dimension()), neighCenterX(3);
                    neigh->CenterPoint(neigh->NSides()-1, neighCenterQSI);
                    neigh->X(neighCenterQSI, neighCenterX);
                    
                    TPZVec<REAL> n0(3), n1(3);
                    refinedMesh->NodeVec()[gel->SideNodeIndex(sd, 0)].GetCoordinates(n0);
                    refinedMesh->NodeVec()[gel->SideNodeIndex(sd, 1)].GetCoordinates(n1);
                    
                    //Como o contorno da fratura foi construido no sentido antihorario no plano x,z (normal Y > 0),
                    //interessam os elementos aa direita do elemento 1D. Portanto eh feito produto vetorial entre os vetores
                    //frac=(n1-n0) e cg_neigh=(cg-n0). O vizinho aa direita apresentarah componente em Y positiva.
                    REAL crossYcomp = n0[2]*n1[0] - n0[0]*n1[2] -
                    n0[2]*neighCenterX[0] + n1[2]*neighCenterX[0] +
                    n0[0]*neighCenterX[2] - n1[0]*neighCenterX[2];
                    
                    if(crossYcomp > 0.)
                    {
                        neigh->SetMaterialId(__2DfractureMat_inside);
                        
                        std::map<int,TPZFracture2DEl>::iterator edgIt_temp = fracturedElems.find(neigh->Id());
                        if(edgIt_temp != fracturedElems.end())
                        {
                            TPZFracture2DEl neighFractEl = edgIt_temp->second;
                            neighFractEl.RemoveThisEdge(sideNeighbyside);
                            edgIt_temp->second = neighFractEl;
                        }
                        else
                        {
                            TPZFracture2DEl fractEl(neigh);
                            fractEl.RemoveThisEdge(sideNeighbyside);
                            fracturedElems[fractEl.Id()] = fractEl;
                        }
                    }
                }
                
                sideNeigh = sideNeigh.Neighbour();
            }
        }
    }
    
    //capturanto demais elementos no interior da fratura
    std::set<int> finishedFracturedElems;
    while(fracturedElems.size() > 0)
    {
        std::map<int,TPZFracture2DEl>::iterator edgIt = fracturedElems.begin();
        
        TPZFracture2DEl actEl = edgIt->second;
        
        std::set<int>::iterator sideIt;
        for(sideIt = actEl.fEdge.begin(); sideIt != actEl.fEdge.end(); sideIt++)
        {
            int side = *sideIt;
            TPZGeoElSide actElEdge(actEl.fElem2D,side);
            TPZGeoElSide neighElSide = actElEdge.Neighbour();
            
            bool wellDone = false;
            while(actElEdge != neighElSide && wellDone == false)
            {
                int sideNeighbyside = neighElSide.Side();
                TPZGeoEl * neighEl = neighElSide.Element();
                if(neighEl->Dimension() == 2 && neighEl->MaterialId() == __2DfractureMat_outside && !neighEl->HasSubElement())
                {
                    int neighElId = neighEl->Id();
                    if(finishedFracturedElems.find(neighElId) != finishedFracturedElems.end())
                    {
                        wellDone = true;
                        continue;
                    }
                    std::map<int,TPZFracture2DEl>::iterator edgIt_temp = fracturedElems.find(neighElId);
                    if(edgIt_temp != fracturedElems.end())
                    {
                        TPZFracture2DEl neighFractEl = edgIt_temp->second;
                        neighFractEl.RemoveThisEdge(sideNeighbyside);
                        
                        if(neighFractEl.IsOver())
                        {
                            finishedFracturedElems.insert(neighFractEl.Id());
                            fracturedElems.erase(edgIt_temp);
                        }
                        else
                        {
                            edgIt_temp->second = neighFractEl;
                        }
                    }
                    else
                    {
                        neighEl->SetMaterialId(__2DfractureMat_inside);
                        
                        TPZFracture2DEl fractEl(neighEl);
                        fractEl.RemoveThisEdge(sideNeighbyside);
                        fracturedElems[fractEl.Id()] = fractEl;
                    }
                    
                    wellDone = true;
                }
                neighElSide = neighElSide.Neighbour();
            }
        }
        
        finishedFracturedElems.insert(actEl.Id());
        fracturedElems.erase(edgIt);
    }
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::TurnIntoQuarterPoint(TPZGeoMesh * refinedMesh)
{
    RefinementProceedings(refinedMesh);
    
    for(int i = 0; i < fcrackBoundaryElementsIds.NElements(); i++)
    {
        TPZGeoEl * gel1D = refinedMesh->ElementVec()[fcrackBoundaryElementsIds[i]];
        
#ifdef DEBUG
        if(gel1D->Dimension() != 1)
        {
            DebugStop();
        }
#endif
        
        for(int s = 0; s < 2; s++)
        {
            TPZGeoElSide edge(gel1D,s);
            TPZGeoElSide neigh(edge.Neighbour());
            
            while(edge != neigh)
            {
                if(neigh.Element()->HasSubElement() == false && neigh.Element()->IsLinearMapping())
                {
                    if(neigh.Element()->Dimension() == 3)
                    {
                        neigh.Element()->SetMaterialId(__3DrockMat_quarterPoint);
                    }
                    int neighSide = neigh.Side();
                    TPZGeoEl * neighEl = TPZChangeEl::ChangeToQuarterPoint(refinedMesh, neigh.Element()->Id(), neighSide);
                    neigh = neighEl->Neighbour(neighSide);
                }
                else
                {
                    neigh = neigh.Neighbour();
                }
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::RefinementProceedings(TPZGeoMesh * refinedMesh)
{return;
    REAL desiredSize = fLmax;
    int ndiv = log((fLmax/2.)/desiredSize)/log(2.);
    if(ndiv < 1)
    {
        ndiv = 1;
    }
    
    std::set<int> fracturePlaneMat;
    fracturePlaneMat.insert(__2DfractureMat_inside);
    
    for(int div = 0; div < ndiv; div++)
    {
        int nelem = refinedMesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            TPZGeoEl * gel = refinedMesh->ElementVec()[el];
            if(gel->HasSubElement() == false)
            {
                TPZRefPatternTools::RefineDirectional(gel, fracturePlaneMat);
            }
        }
    }
    
    std::set<int> crackTipMat;
    crackTipMat.insert(__1DcrackTipMat);
    for(int div = 0; div < ndiv; div++)
    {
        int nelem = refinedMesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            TPZGeoEl * gel = refinedMesh->ElementVec()[el];
            if(gel->HasSubElement() == false)
            {
                TPZRefPatternTools::RefineDirectional(gel, crackTipMat);
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------


bool TPZPlaneFracture::IsBoundaryMaterial(TPZGeoEl * gel)
{
    int materialId = gel->MaterialId();
    
    if(materialId == __2DfarfieldMat)
    {
        return true;
    }
    else if(materialId == __2DleftMat)
    {
        return true;
    }
    else if(materialId == __2DrightMat)
    {
        return true;
    }
    else if(materialId == __2DtopMat)
    {
        return true;
    }
    else if(materialId == __2DbottomMat)
    {
        return true;
    }
    
    return false;
}
//------------------------------------------------------------------------------------------------------------


//** just for visualize given dots in vtk */
void TPZPlaneFracture::InsertDots4VTK(TPZGeoMesh * gmesh, const TPZVec<REAL> &fractureDots)
{
    int nDots = fractureDots.size() / 2;
    long nnodesOriginal = gmesh->NNodes();
	long Qnodes = nnodesOriginal + nDots;
    
	//initializing gmesh->NodeVec()
	gmesh->NodeVec().Resize(Qnodes);
	TPZGeoNode Node;
    TPZVec<REAL> NodeCoord(3);
    TPZVec<long> Topol(1);
    
    long elId = gmesh->NElements();
	for(long n = nnodesOriginal; n < Qnodes; n++)
	{
        Topol[0] = n;
        
        long actDot = n - nnodesOriginal;
        
        NodeCoord[0] = fractureDots[2*actDot];
        NodeCoord[1] = 0.;
        NodeCoord[2] = fractureDots[2*actDot + 1];
        
		Node.SetNodeId(n);
		Node.SetCoord(NodeCoord);
		gmesh->NodeVec()[n] = Node;
        
        int matPoint = -100;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol, matPoint,*gmesh);
        elId++;
	}
}
//------------------------------------------------------------------------------------------------------------







