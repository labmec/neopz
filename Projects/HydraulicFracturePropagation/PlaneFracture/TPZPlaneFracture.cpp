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


TPZPlaneFracture::TPZPlaneFracture(TPZVec<TPZLayerProperties> & layerVec, REAL bulletDepthTVDIni, REAL bulletDepthTVDFin,
                                   REAL xLength, REAL yLength, REAL Lmax, int nstripes)

{
    fLayerVec = layerVec;
    
    fInitialElIndex = 0;
    
    fPreservedMesh = new TPZGeoMesh;
    
    fLmax = Lmax;
    fLfrac = 0.;
    fnstripes = nstripes;
    
    std::set<REAL> espacamentoVerticalTVD;
    
    espacamentoVerticalTVD.insert(bulletDepthTVDIni);
    espacamentoVerticalTVD.insert(bulletDepthTVDFin);
    
    //Inserindo TVDs impostos (TVDs das camadas fornecidas)
    //>>>>>>>> Obs.: Eh considerado que o TVDfin da camada (s) corresponde ao TVDini da camada (s+1)
    int nstretches = layerVec.NElements();
    for(int s = 0; s < nstretches; s++)
    {
        REAL pos = layerVec[s].fTVDini;
        espacamentoVerticalTVD.insert(pos);
        if(s == (nstretches-1))
        {
            pos = layerVec[s].fTVDfin;
            espacamentoVerticalTVD.insert(pos);
        }
    }
    
    //Inserindo TVDs intermediarios para geracao da malha geometrica (baseados no Lmax)
    std::set<REAL>::iterator itS = espacamentoVerticalTVD.begin();
    REAL pos0 = *itS;
    itS++;
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
    
    
    //Conversao de TVD para profundidade (MD)
    std::list<REAL> espacamentoVerticalDEPTH;
    for(itS = espacamentoVerticalTVD.begin(); itS != espacamentoVerticalTVD.end(); itS++)
    {
        REAL posDepth = *itS;
        espacamentoVerticalDEPTH.push_back(-posDepth);
    }
    
    GeneratePreservedMesh(espacamentoVerticalDEPTH, bulletDepthTVDIni, bulletDepthTVDFin, xLength, yLength);
}
//------------------------------------------------------------------------------------------------------------


TPZPlaneFracture::~TPZPlaneFracture()
{
    delete fPreservedMesh;
}
//------------------------------------------------------------------------------------------------------------


#define elastLinear
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
    
    if(printVTKfile)
    {
        TPZManVector<std::string,10> scalnames(3), vecnames(1);
        
        //        scalnames[0] = "EDisplacementX";
        //        scalnames[1] = "EDisplacementY";
        scalnames[0] = "StressX";
        scalnames[1] = "StressY";
        scalnames[2] = "StressZ";
        
        vecnames[0] = "Displacement";
        
        const int dim = 3;
        int div = 0;
        anRef.DefineGraphMesh(dim,scalnames,vecnames,vtkFile);
        anRef.PostProcess(div);
    }

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
    
    if(printVTKfile)
    {
        TPZManVector<std::string,10> scalnames(3), vecnames(1);
    
    //        scalnames[0] = "EDisplacementX";
    //        scalnames[1] = "EDisplacementY";
        scalnames[0] = "StressX";
        scalnames[1] = "StressY";
        scalnames[2] = "StressZ";
        
        vecnames[0] = "Displacement";

        std::string vtkFile1 = "fracturePconstant1.vtk";
        const int dim = 3;
        int div = 0;
        an.DefineGraphMesh(dim,scalnames,vecnames,vtkFile1);
        an.PostProcess(div);
    }
    
    
    /////// Example of J-Integral
    
    //    TPZVec<REAL> originXYZ(3,0.), direction(3,0.);
    //
    //    int POSmiddle1D = int(REAL(fcrackBoundaryElementsIndexes.NElements())/2. + 0.5);
    //    int middle1DId = fcrackBoundaryElementsIndexes[POSmiddle1D];
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
    
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, dirichDir = 3;
    //int mixed = 2;
    
    for(int lay = 0; lay < fLayerVec.NElements(); lay++)
    {
        STATE young = fLayerVec[lay].fYoung;
        STATE poisson = fLayerVec[lay].fPoisson;
        TPZVec<STATE> force(3,0.);
        STATE prestressXX = fLayerVec[lay].fSigmaMax;
        STATE prestressYY = fLayerVec[lay].fSigmaMin;
        STATE prestressZZ = 0.;
//        STATE prestressXX = 0.;
//        STATE prestressYY = 0.;
//        STATE prestressZZ = 0.;
        
        ////Rock
        TPZMaterial * materialLin = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force, prestressXX, prestressYY, prestressZZ);
        cmesh->InsertMaterialObject(materialLin);
        
        ////BCs
        f.Zero();
        f(0,0) = 1.;
        TPZMaterial * materialMixedLeft = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,globMaterialIdGen.LeftMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedLeft);
        
        TPZMaterial * materialMixedRight = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedRight = new TPZBndCond(materialMixedRight,globMaterialIdGen.RightMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedRight);
        
        f.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialMixedOutFracture = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,globMaterialIdGen.OutSideFractMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);
        
        f.Zero();
        f(2,0) = 1.;
        TPZMaterial * materialMixedTop = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,globMaterialIdGen.TopMatId(), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        
        TPZMaterial * materialMixedBottom = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,globMaterialIdGen.BottomMatId(), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedBottom);
        
        k.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialNewmannFarField = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,globMaterialIdGen.FarfieldMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f.Zero();
        f(1,0) = pressureInsideCrack;
        for(int stripe = 0; stripe < fnstripes; stripe++)
        {
            TPZMaterial * materialNewmannInsideFract = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
            TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,globMaterialIdGen.InsideFractMatId(lay, stripe), newmann, k, f);
            cmesh->InsertMaterialObject(newmannInsideFract);
        }
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
    
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, dirichDir = 3;
    //int mixed = 2;
    
    for(int lay = 0; lay < fLayerVec.NElements(); lay++)
    {
        STATE young = fLayerVec[lay].fYoung;
        STATE poisson = fLayerVec[lay].fPoisson;
        TPZVec<STATE> force(3,0.);
        STATE prestressXX = fLayerVec[lay].fSigmaMax;
        STATE prestressYY = fLayerVec[lay].fSigmaMin;
        STATE prestressZZ = 0.;
//        STATE prestressXX = 0.;
//        STATE prestressYY = 0.;
//        STATE prestressZZ = 0.;

        ////Rock
        TPZMaterial * materialLin = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force, prestressXX, prestressYY, prestressZZ);
        cmesh->InsertMaterialObject(materialLin);
    
        ////BCs
        f.Zero();
        f(0,0) = 1.;
        TPZMaterial * materialMixedLeft = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,globMaterialIdGen.LeftMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedLeft);
        
        TPZMaterial * materialMixedRight = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedRight = new TPZBndCond(materialMixedRight,globMaterialIdGen.RightMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedRight);
        
        f.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialMixedOutFracture = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,globMaterialIdGen.OutSideFractMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);
        
        f.Zero();
        f(2,0) = 1.;
        TPZMaterial * materialMixedTop = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,globMaterialIdGen.TopMatId(), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        //
        TPZMaterial * materialMixedBottom = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,globMaterialIdGen.BottomMatId(), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedBottom);
        
        ///////////farField
        k.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialNewmannFarField = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,globMaterialIdGen.FarfieldMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f.Zero();
        f(1,0) = pressureInsideCrack;
        for(int stripe = 0; stripe < fnstripes; stripe++)
        {
            TPZMaterial * materialNewmannInsideFract = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
            TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,globMaterialIdGen.InsideFractMatId(lay, stripe), newmann, k, f);
            cmesh->InsertMaterialObject(newmannInsideFract);
        }
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
    
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, dirichDir = 3;
    //    int mixed = 2;

    for(int lay = 0; lay < fLayerVec.NElements(); lay++)
    {
        STATE young = fLayerVec[lay].fYoung;
        STATE poisson = fLayerVec[lay].fPoisson;
        TPZVec<STATE> force(3,0.);
        STATE prestressXX = fLayerVec[lay].fSigmaMax;
        STATE prestressYY = fLayerVec[lay].fSigmaMin;
        STATE prestressZZ = 0.;

        ////Rock
        TPZMaterial * materialLin = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force, prestressXX, prestressYY, prestressZZ);
        cmesh->InsertMaterialObject(materialLin);

        ////BCs
        f.Zero();
        f(0,0) = 1.;
        TPZMaterial * materialMixedLeft = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,globMaterialIdGen.LeftMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedLeft);

        TPZMaterial * materialMixedRight = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedRight = new TPZBndCond(materialMixedRight,globMaterialIdGen.RightMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedRight);
        
        f.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialMixedOutFracture = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,globMaterialIdGen.OutSideFractMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);
        
        f.Zero();
        f(2,0) = 1.;
        TPZMaterial * materialMixedTop = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,globMaterialIdGen.TopMatId(), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        //
        TPZMaterial * materialMixedBottom = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,globMaterialIdGen.BottomMatId(), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedBottom);
        
        ///////////farField
        k.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialNewmannFarField = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,globMaterialIdGen.FarfieldMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f.Zero();
        f(1,0) = pressureInsideCrack;
        for(int stripe = 0; stripe < fnstripes; stripe++)
        {
            TPZMaterial * materialNewmannInsideFract = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force);
            TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,globMaterialIdGen.InsideFractMatId(lay, stripe), newmann, k, f);
            cmesh->InsertMaterialObject(newmannInsideFract);
        }
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
    
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, dirichDir = 3;
    //    int mixed = 2;
    
    for(int lay = 0; lay < fLayerVec.NElements(); lay++)
    {
        STATE young = fLayerVec[lay].fYoung;
        STATE poisson = fLayerVec[lay].fPoisson;
        TPZVec<STATE> force(3,0.);
        STATE prestressXX = fLayerVec[lay].fSigmaMax;
        STATE prestressYY = fLayerVec[lay].fSigmaMin;
        STATE prestressZZ = 0.;
        
        TPZMaterial * materialLin = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force, prestressXX, prestressYY, prestressZZ);
        cmesh->InsertMaterialObject(materialLin);
    
        ////BCs
        f.Zero();
        f(0,0) = 1.;
        TPZMaterial * materialMixedLeft = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,globMaterialIdGen.LeftMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedLeft);

        TPZMaterial * materialMixedRight = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedRight = new TPZBndCond(materialMixedRight,globMaterialIdGen.RightMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedRight);
        
        f.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialMixedOutFracture = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,globMaterialIdGen.OutSideFractMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);
        
        f.Zero();
        f(2,0) = 1.;
        TPZMaterial * materialMixedTop = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,globMaterialIdGen.TopMatId(), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        //
        TPZMaterial * materialMixedBottom = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,globMaterialIdGen.BottomMatId(), dirichDir, k, f);
        cmesh->InsertMaterialObject(mixedBottom);
        
        ///////////farField
        k.Zero();
        f(1,0) = 1.;
        TPZMaterial * materialNewmannFarField = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,globMaterialIdGen.FarfieldMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f.Zero();
        f(1,0) = pressureInsideCrack;
        for(int stripe = 0; stripe < fnstripes; stripe++)
        {
            TPZMaterial * materialNewmannInsideFract = new TPZElast3Dnlinear(globMaterialIdGen.RockMatId(lay), young, poisson, force);
            TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,globMaterialIdGen.InsideFractMatId(lay, stripe), newmann, k, f);
            cmesh->InsertMaterialObject(newmannInsideFract);
        }
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
    
	std::map< long, std::set<REAL> > elIndex_TrimCoords;
	std::list< std::pair<long,REAL> > elIndexSequence;
	
	DetectEdgesCrossed(poligonalChain, refinedMesh, elIndex_TrimCoords, elIndexSequence);
	
	//Refining auxiliar 1D elements
	TPZVec<TPZGeoEl*> sons;
	std::map< long, std::set<REAL> >::iterator it;
	for(it = elIndex_TrimCoords.begin(); it != elIndex_TrimCoords.end(); it++)
	{
		int el1DIndex = it->first;
        TPZGeoEl * el1D = refinedMesh->ElementVec()[el1DIndex];
        
		TPZAutoPointer<TPZRefPattern> linRefp = Generate1DRefPatt(it->second);
        
		el1D->SetRefPattern(linRefp);
		el1D->Divide(sons);
	}
    
	//Refining 2D and 3D elements with the intention to match the geometry of the crack boundary
	for(long el = 0; el < nelem; el++)
	{
		TPZGeoEl * gel = refinedMesh->ElementVec()[el];//2D element in 2D mesh
        if(globMaterialIdGen.IsOutsideFractMat(gel->MaterialId()) == false)
        {
            continue;
        }
        
#ifdef DEBUG
        if(gel->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
		TPZAutoPointer<TPZRefPattern> elRefp = TPZRefPatternTools::PerfectMatchRefPattern(gel);
		if(elRefp)
		{
            gel->SetRefPattern(elRefp);
            gel->Divide(sons);
            
            int innerSide = gel->NSides() - 1;
            gel = gel->Neighbour(innerSide).Element();//3D element in 3D mesh
            
#ifdef DEBUG
            if(gel->Dimension() != 3)
            {
                DebugStop();
            }
#endif

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
	}
	
	GenerateCrackBoundary(refinedMesh, elIndexSequence);
    SeparateElementsInMaterialSets(refinedMesh);
    //TurnIntoQuarterPoint(refinedMesh);
    
//    std::ofstream cc("cc.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(refinedMesh, cc, true);
    
	return refinedMesh;
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::GeneratePreservedMesh(std::list<REAL> & espacamento,
                                             REAL bulletDepthTVDIni, REAL bulletDepthTVDFin,
                                             REAL xLength, REAL yLength)
{
    REAL bulletDepthDEPTHIni = -bulletDepthTVDIni;
    REAL bulletDepthDEPTHFin = -bulletDepthTVDFin;
    
    TPZVec< TPZVec<REAL> > NodeCoord(0);
    long nrows, ncols;
    
    //std::list<REAL>::iterator it = espacamento.end(); it--;
    
    int nDirRef = int(log(fabs(yLength)/fLmax)/log(2.));
    int nLayersY = nDirRef + 2;
    
    REAL Y = 0.;
    GenerateNodesAtPlaneY(espacamento, xLength, NodeCoord, nrows, ncols, Y);
    
    Y = yLength/pow(2.,nDirRef);
    for(int lay = 1; lay < nLayersY; lay++)
    {
        GenerateNodesAtPlaneY(espacamento, xLength, NodeCoord, nrows, ncols, Y);
        Y *= 2.;
    }
    long nNodesByLayer = nrows*ncols;
    long Qnodes = nNodesByLayer * nLayersY;
	
	//initializing gmesh->NodeVec()
	fPreservedMesh->NodeVec().Resize(Qnodes);
    
    long pos = 0;
	TPZGeoNode Node;
    for(long l = 0; l < nLayersY; l++)
    {
        for(long n = 0; n < nNodesByLayer; n++)
        {
            Node.SetNodeId(pos);
            Node.SetCoord(NodeCoord[pos]);
            fPreservedMesh->NodeVec()[pos] = Node;
            pos++;
        }
    }
	
	//inserting quadrilaterals and 1D bullet region
	TPZVec<long> Topol2(2), Topol4(4), Topol8(8);
	for(long r = 0; r < (nrows-1); r++)
	{
		for(long c = 0; c < (ncols-1); c++)
		{
			Topol4[0] = ncols*r+c; Topol4[1] = ncols*r+c+1; Topol4[2] = ncols*(r+1)+c+1; Topol4[3] = ncols*(r+1)+c;
            
            REAL z0 = fPreservedMesh->NodeVec()[Topol4[0]].Coord(2);
            REAL z3 = fPreservedMesh->NodeVec()[Topol4[3]].Coord(2);
            REAL zMed = (z0 + z3)/2.;
            int whatLayer = GetLayer(zMed);
            
			new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.OutSideFractMatId(whatLayer),*fPreservedMesh);
            
            //detecting bullet region
            if(c == 0)
            {
                REAL tol = 1.E-3;
                if(zMed < (bulletDepthDEPTHIni + tol) && zMed > (bulletDepthDEPTHFin - tol))
                {
                    Topol2[0] = ncols*r;
                    Topol2[1] = ncols*(r+1);
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol2,globMaterialIdGen.BulletMatId(whatLayer),*fPreservedMesh);
                }
            }
		}
	}
    
    //inserting hexaedrons
    for(long ly = 0; ly < (nLayersY-1); ly++)
    {
        for(long r = 0; r < (nrows-1); r++)
        {
            for(long c = 0; c < (ncols-1); c++)
            {
                Topol8[0] = ncols*r+c + ly*nNodesByLayer;
                Topol8[1] = ncols*r+c+1 + ly*nNodesByLayer;
                Topol8[2] = ncols*(r+1)+c+1 + ly*nNodesByLayer;
                Topol8[3] = ncols*(r+1)+c + ly*nNodesByLayer;
                //
                Topol8[4] = ncols*r+c + (ly+1)*nNodesByLayer;
                Topol8[5] = ncols*r+c+1 + (ly+1)*nNodesByLayer;
                Topol8[6] = ncols*(r+1)+c+1 + (ly+1)*nNodesByLayer;
                Topol8[7] = ncols*(r+1)+c + (ly+1)*nNodesByLayer;
             
                REAL z0 = fPreservedMesh->NodeVec()[Topol8[0]].Coord(2);
                REAL z3 = fPreservedMesh->NodeVec()[Topol8[3]].Coord(2);
                REAL zMed = (z0+z3)/2.;
                int whatLayer = GetLayer(zMed);
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (Topol8,globMaterialIdGen.RockMatId(whatLayer),*fPreservedMesh);
                
                if(ly == (nLayersY-2))//farfield cc
                {
                    Topol4[0] = ncols*r+c + (ly+1)*nNodesByLayer;
                    Topol4[1] = ncols*r+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[2] = ncols*(r+1)+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*(r+1)+c + (ly+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.FarfieldMatId(whatLayer),*fPreservedMesh);
                }
                if(c == 0)//left cc
                {
                    Topol4[0] = ncols*r+c + ly*nNodesByLayer;
                    Topol4[1] = ncols*r+c + (ly+1)*nNodesByLayer;
                    Topol4[2] = ncols*(r+1)+c + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*(r+1)+c + ly*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.LeftMatId(whatLayer),*fPreservedMesh);
                }
                else if(c == (ncols - 2))//right cc
                {
                    Topol4[0] = ncols*r+c+1 + ly*nNodesByLayer;
                    Topol4[1] = ncols*(r+1)+c+1 + ly*nNodesByLayer;
                    Topol4[2] = ncols*(r+1)+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*r+c+1 + (ly+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.RightMatId(whatLayer),*fPreservedMesh);
                }
                if(r == 0)//top cc
                {
                    Topol4[0] = ncols*r+c + ly*nNodesByLayer;
                    Topol4[1] = ncols*r+c+1 + ly*nNodesByLayer;
                    Topol4[2] = ncols*r+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*r+c + (ly+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.TopMatId(),*fPreservedMesh);
                }
                else if(r == (nrows - 2))//bottom cc
                {
                    Topol4[0] = ncols*(r+1)+c + ly*nNodesByLayer;
                    Topol4[1] = ncols*(r+1)+c + (ly+1)*nNodesByLayer;
                    Topol4[2] = ncols*(r+1)+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*(r+1)+c+1 + ly*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.BottomMatId(),*fPreservedMesh);
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
    std::list<REAL>::iterator it;
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
                                          std::map< long, std::set<REAL> > &elIndex_TrimCoords,
                                          std::list< std::pair<long,REAL> > &elIndexSequence)
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
        
        fLfrac = std::max(fLfrac,xNext[0]);
        
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
			nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elIndex_TrimCoords, elIndexSequence, true, axe0, axe1, axeNormal, false);
			
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
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elIndex_TrimCoords, elIndexSequence, false, axe0, axe1, axeNormal, true);
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
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elIndex_TrimCoords, elIndexSequence, true, axe0, axe1, axeNormal, true);
        if(!nextGel)
        {
            DebugStop();
        }
		alphaMin = __smallNum;
	}
}
//------------------------------------------------------------------------------------------------------------


TPZGeoEl * TPZPlaneFracture::CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, REAL alphaMin,
												  std::map< long, std::set<REAL> > &elIndex_TrimCoords,
                                                  std::list< std::pair<long,REAL> > &elIndexSequence,
                                                  bool pushback, int planeAxe0, int planeAxe1, int planeNormal,
                                                  bool closingFracture)
{
	bool thereIsAn1DElemAlready;
	int edge;
	std::map< long, std::set<REAL> >::iterator it;
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
                
                std::cout << "Element " << gel->Index() << std::endl;
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
				TPZTransform transBetweenNeigh = neighEdge.NeighbourSideTransform(gelEdge);
				qsi1D *= transBetweenNeigh.Mult()(0,0);
                long neighEdgeIndex = neighEdge.Element()->Index();
				it = elIndex_TrimCoords.find(neighEdgeIndex);
                if(it == elIndex_TrimCoords.end())//Deve ser elemento 1D de materialId == globMaterialIdGen.__1DbulletMat
                {
                    std::set<REAL> trim;
                    trim.insert(qsi1D);
                    elIndex_TrimCoords[neighEdgeIndex] = trim;
                }
                else
                {
                    it->second.insert(qsi1D);
                }
				if(pushback)
				{
					elIndexSequence.push_back(std::make_pair(neighEdgeIndex, qsi1D));
				}
				else // push_FRONT
				{
					elIndexSequence.push_front(std::make_pair(neighEdgeIndex, qsi1D));
				}
                
				break;
			}
			neighEdge = neighEdge.Neighbour();
		}
		if(thereIsAn1DElemAlready == false)//nao existe um elemento 1D nesta aresta!
		{
			std::set<REAL> trim;
			trim.insert(qsi1D);
			Topol[0] = gel->SideNodeIndex(edge, 0);
			Topol[1] = gel->SideNodeIndex(edge, 1);
			TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * linGeo =
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, globMaterialIdGen.Aux1DMatId(), *planeMesh);
            
			planeMesh->BuildConnectivity();
			
			long linGeoIndex = linGeo->Index();
			elIndex_TrimCoords[linGeoIndex] = trim;
			
			if(pushback) // push_BACK
			{
				elIndexSequence.push_back(std::make_pair(linGeoIndex, qsi1D));
			}
			else // push_FRONT
			{
				elIndexSequence.push_front(std::make_pair(linGeoIndex, qsi1D));
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
           globMaterialIdGen.IsOutsideFractMat(neighEdge.Element()->MaterialId()) &&
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
            globMaterialIdGen.IsOutsideFractMat(gelEdge.Element()->MaterialId()) &&
            gelEdge.Element()->Father() == NULL)
    {
        //fechando fratura no inicio ou no final, portanto o
        //ultimo elemento nao encontra vizinho do tipo 2d com
        //materialId == globMaterialIdGen.__2DfractureMat_outside, retornando a si mesmo
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
	TPZVec<long> Topol(2);
	
	//4.1 inserting father element
	Topol[0] = 0; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol,globMaterialIdGen.Aux1DMatId(),internalMesh);
	//.
	
	//4.2 inserting subelements
	//first subelement
	Topol[0] = 0; Topol[1] = 2;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son1 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol,globMaterialIdGen.Aux1DMatId(),internalMesh);
	son1->SetFather(father);
	son1->SetFather(father->Index());
	//
	
	//last subelement
	Topol[0] = TrimCoord.size() + 1; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son2 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol,globMaterialIdGen.Aux1DMatId(),internalMesh);
	son2->SetFather(father);
	son2->SetFather(father->Index());
	//
	
	for(long el = 2; el < QsubElements; el++)
	{
		Topol[0] = el; Topol[1] = el+1;
		TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol,globMaterialIdGen.Aux1DMatId(),internalMesh);
		son->SetFather(father);
		son->SetFather(father->Index());
	}
	//.
	//.
	
	internalMesh.BuildConnectivity();
    internalMesh.SetMaxNodeId(internalMesh.NNodes()-1);
	internalMesh.SetMaxElementId(internalMesh.NElements()-1);
	
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
                                            std::list< std::pair<long,REAL> > &elIndexSequence,
                                            TPZVec<std::pair<REAL,REAL> > &poligonalChainUpdated)
{
	int nptos = elIndexSequence.size();
	poligonalChainUpdated.Resize(nptos);
    
	TPZVec<REAL> qsi1Dvec(1), ptoCoord(3);
	int el1Did, p = 0;
	REAL qsi1D;
	
	std::list< std::pair<long,REAL> >::iterator it;
	for(it = elIndexSequence.begin(); it != elIndexSequence.end(); it++)
	{
		el1Did = it->first;
		TPZGeoEl * el1D = gmesh->ElementVec()[el1Did];
		
#ifdef DEBUG
		int elDim = el1D->Dimension();
		if(elDim != 1)
		{
			std::cout << "The elIndexSequence supposedly would contains ids of elements exclusively 1D!" << std::endl;
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
                                             std::list< std::pair<long,REAL> > &elIndexSequence)
{
    fcrackBoundaryElementsIndexes.Resize(0);
	TPZVec<REAL> qsi0vec(1), qsi1vec(1), node0coord(3), node1coord(3);
	TPZVec<long> Topol(2);
	long el0id, el1id, n0, n1;
	REAL qsi0, qsi1;
	std::list< std::pair<long,REAL> >::iterator crackit0, crackit1, crackitEnd;
	crackitEnd = elIndexSequence.end(); crackitEnd--;
	for(crackit0 = elIndexSequence.begin(); crackit0 != crackitEnd; crackit0++)
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
		
		TPZGeoEl * crack1D = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, globMaterialIdGen.CrackTipMatId(), *refinedMesh);
        int oldSize = fcrackBoundaryElementsIndexes.NElements();
        fcrackBoundaryElementsIndexes.Resize(oldSize+1);
        fcrackBoundaryElementsIndexes[oldSize] = crack1D->Index();
	}
    
    refinedMesh->BuildConnectivity();
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::SeparateElementsInMaterialSets(TPZGeoMesh * refinedMesh)
{
    int n1Dels = fcrackBoundaryElementsIndexes.NElements();
    std::map<int,TPZFracture2DEl> fracturedElems;
    
    //Capturando subelementos que encostam no contorno da fratura
    for(int el = 0; el < n1Dels; el++)
    {
        int cracktipIndex = fcrackBoundaryElementsIndexes[el];
        TPZGeoEl * gel = refinedMesh->ElementVec()[cracktipIndex];//1D element of crach boundary
        
#ifdef DEBUG
        if(gel->Dimension() != 1)
        {
            DebugStop();
        }
#endif
        int sd = 2;
        TPZGeoElSide side1D(gel,sd);
        TPZGeoElSide sideNeigh = side1D.Neighbour();
        while(sideNeigh != side1D)
        {
            int sideNeighbyside = sideNeigh.Side();
            TPZGeoEl * neigh = sideNeigh.Element();
            
            if(neigh->HasSubElement() || neigh->Dimension() != 2)
            {
                sideNeigh = sideNeigh.Neighbour();
                continue;
            }
            
            //Primeiro estagio da identificacao dos elementos 2D que estao no interior da fratura.
            //brief: primeiramente sao identificados os elementos no interior da fratura que encostam no crack tip.
            //Obs.: A condicao (IsBoundaryMaterial(neigh) == false) eh para excluir os quadrilateros de condicao de contorno.
            if(IsBoundaryMaterial(neigh) == false)
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
                    if(neigh->LowestFather() != neigh)
                    {
                        neigh->LowestFather()->CenterPoint(neigh->LowestFather()->NSides()-1, neighCenterQSI);
                        neigh->LowestFather()->X(neighCenterQSI, neighCenterX);
                    }
                    REAL Xc = neighCenterX[0];
                    REAL Zc = neighCenterX[2];
                    int stripe = std::min( fnstripes-1 , (int)(Xc/(fLfrac/fnstripes)) );
                    int layer = GetLayer(Zc);
                    
                    neigh->SetMaterialId(globMaterialIdGen.InsideFractMatId(layer, stripe));
                    
                    std::map<int,TPZFracture2DEl>::iterator edgIt_temp = fracturedElems.find(neigh->Index());
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
                        fracturedElems[fractEl.Index()] = fractEl;
                    }
                }
            }
            
            sideNeigh = sideNeigh.Neighbour();
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
                if(neighEl->Dimension() == 2 && globMaterialIdGen.IsOutsideFractMat(neighEl->MaterialId()) && !neighEl->HasSubElement())
                {
                    int neighElIndex = neighEl->Index();
                    if(finishedFracturedElems.find(neighElIndex) != finishedFracturedElems.end())
                    {
                        wellDone = true;
                        continue;
                    }
                    std::map<int,TPZFracture2DEl>::iterator edgIt_temp = fracturedElems.find(neighElIndex);
                    if(edgIt_temp != fracturedElems.end())
                    {
                        TPZFracture2DEl neighFractEl = edgIt_temp->second;
                        neighFractEl.RemoveThisEdge(sideNeighbyside);
                        
                        if(neighFractEl.IsOver())
                        {
                            finishedFracturedElems.insert(neighFractEl.Index());
                            fracturedElems.erase(edgIt_temp);
                        }
                        else
                        {
                            edgIt_temp->second = neighFractEl;
                        }
                    }
                    else
                    {
                        TPZVec<REAL> neighCenterQSI(neighEl->Dimension()), neighCenterX(3);
                        neighEl->LowestFather()->CenterPoint(neighEl->LowestFather()->NSides()-1, neighCenterQSI);
                        neighEl->LowestFather()->X(neighCenterQSI, neighCenterX);
                        
                        REAL Xc = neighCenterX[0];
                        REAL Zc = neighCenterX[2];
                        int stripe = std::min( fnstripes-1 , (int)(Xc/(fLfrac/fnstripes)) );
                        int layer = GetLayer(Zc);
                        
                        neighEl->SetMaterialId(globMaterialIdGen.InsideFractMatId(layer, stripe));
                        
                        TPZFracture2DEl fractEl(neighEl);
                        fractEl.RemoveThisEdge(sideNeighbyside);
                        fracturedElems[fractEl.Index()] = fractEl;
                    }
                    
                    wellDone = true;
                }
                neighElSide = neighElSide.Neighbour();
            }
        }
        
        finishedFracturedElems.insert(actEl.Index());
        fracturedElems.erase(edgIt);
    }
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::TurnIntoQuarterPoint(TPZGeoMesh * refinedMesh)
{
//    RefinementProceedings(refinedMesh);
//    
//    for(int i = 0; i < fcrackBoundaryElementsIndexes.NElements(); i++)
//    {
//        TPZGeoEl * gel1D = refinedMesh->ElementVec()[fcrackBoundaryElementsIndexes[i]];
//        
//#ifdef DEBUG
//        if(gel1D->Dimension() != 1)
//        {
//            DebugStop();
//        }
//#endif
//        
//        for(int s = 0; s < 2; s++)
//        {
//            TPZGeoElSide edge(gel1D,s);
//            TPZGeoElSide neigh(edge.Neighbour());
//            
//            while(edge != neigh)
//            {
//                if(neigh.Element()->HasSubElement() == false && neigh.Element()->IsLinearMapping())
//                {
//                    if(neigh.Element()->Dimension() == 3)
//                    {
//                        neigh.Element()->SetMaterialId(globMaterialIdGen.__3DrockMat_quarterPoint);
//                    }
//                    int neighSide = neigh.Side();
//                    TPZGeoEl * neighEl = TPZChangeEl::ChangeToQuarterPoint(refinedMesh, neigh.Element()->Index(), neighSide);
//                    neigh = neighEl->Neighbour(neighSide);
//                }
//                else
//                {
//                    neigh = neigh.Neighbour();
//                }
//            }
//        }
//    }
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::RefinementProceedings(TPZGeoMesh * refinedMesh)
{
//    REAL desiredSize = fLmax;
//    int ndiv = log((fLmax/2.)/desiredSize)/log(2.);
//    if(ndiv < 1)
//    {
//        ndiv = 1;
//    }
//    
//    std::set<int> fracturePlaneMat;
//    fracturePlaneMat.insert(globMaterialIdGen.__2DfractureMat_inside);
//    
//    for(int div = 0; div < ndiv; div++)
//    {
//        int nelem = refinedMesh->NElements();
//        for(int el = 0; el < nelem; el++)
//        {
//            TPZGeoEl * gel = refinedMesh->ElementVec()[el];
//            if(gel->HasSubElement() == false)
//            {
//                TPZRefPatternTools::RefineDirectional(gel, fracturePlaneMat);
//            }
//        }
//    }
//    
//    std::set<int> crackTipMat;
//    crackTipMat.insert(globMaterialIdGen.__1DcrackTipMat);
//    for(int div = 0; div < ndiv; div++)
//    {
//        int nelem = refinedMesh->NElements();
//        for(int el = 0; el < nelem; el++)
//        {
//            TPZGeoEl * gel = refinedMesh->ElementVec()[el];
//            if(gel->HasSubElement() == false)
//            {
//                TPZRefPatternTools::RefineDirectional(gel, crackTipMat);
//            }
//        }
//    }
}
//------------------------------------------------------------------------------------------------------------


bool TPZPlaneFracture::IsBoundaryMaterial(TPZGeoEl * gel)
{
    int materialId = gel->MaterialId();
 
    if(fabs(materialId) >= 3010 && fabs(materialId) <= 7010)
    {//is left, right, farfield, top or bottom
        return true;
    }
    else
    {
        return false;
    }
}
//------------------------------------------------------------------------------------------------------------


int TPZPlaneFracture::GetLayer(REAL zMed)
{
    for(int lay = 0; lay < fLayerVec.NElements(); lay++)
    {
        if(fabs(zMed) > fLayerVec[lay].fTVDini && fabs(zMed) < fLayerVec[lay].fTVDfin)
        {
            return lay;
        }
    }
    
    DebugStop();//nao achou o layer
    return -1;
}

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
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (Topol, matPoint,*gmesh);
	}
}
//------------------------------------------------------------------------------------------------------------







