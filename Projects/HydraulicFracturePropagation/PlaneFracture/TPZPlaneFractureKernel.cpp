//
//  TPZPlaneFractureKernel.cpp
//  PZ
//
//  Created by Cesar Lucci on 18/11/13.
//
//

#include "TPZPlaneFractureKernel.h"

#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

TPZPlaneFractureKernel::TPZPlaneFractureKernel()
{
    fPlaneFractureMesh = NULL;
}

TPZPlaneFractureKernel::TPZPlaneFractureKernel(TPZVec<TPZLayerProperties> & layerVec, REAL bulletDepthTVDIni, REAL bulletDepthTVDFin,
                                               REAL xLength, REAL yLength, REAL Lmax, int nstripes)
{
    fPlaneFractureMesh = new TPZPlaneFractureMesh(layerVec, bulletDepthTVDIni, bulletDepthTVDFin, xLength, yLength, Lmax, nstripes);
}

TPZPlaneFractureKernel::~TPZPlaneFractureKernel()
{
    delete fPlaneFractureMesh;
}

#define elastLinear
void TPZPlaneFractureKernel::RunThisFractureGeometry(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                     REAL pressureInsideCrack,
                                                     std::string vtkFile,
                                                     bool printVTKfile)
{
    int porder = 1;
    
    //--------------------------------------------------------------------------------------------
    REAL pressureInsideCrackRef = 1.;//Pressao da malha computacional referenciada
    
#ifdef elastLinear
    TPZCompMesh * fractureCMeshRef = fPlaneFractureMesh->GetFractureCompMesh(poligonalChain, porder, pressureInsideCrackRef);
#else
    TPZCompMesh * fractureCMeshRef = fPlaneFractureMesh->GetFractureCompMeshNLinear(poligonalChain, porder, pressureInsideCrackRef);
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
    TPZCompMeshReferred * fractureCMesh = fPlaneFractureMesh->GetFractureCompMeshReferred(porder, pressureInsideCrack, fractureCMeshRef);
#else
    TPZCompMeshReferred * fractureCMesh = fPlaneFractureMesh->GetFractureCompMeshReferredNLinear(porder, pressureInsideCrack, fractureCMeshRef);
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