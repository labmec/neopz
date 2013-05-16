//
//  toolstransienttime.cpp
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>

#include "toolstransienttime.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSpStructMatrix.h"
#include "pzelastpressure.h"
#include "pznlfluidstructure2d.h"
#include "pzfstrmatrix.h"
#include "../HydraulicFracturePropagation/PlaneFracture/TPZJIntegral.h"
#include "pzreducedspace.h"
#include "tpzcompmeshreferred.h"
#include "pzbndcond.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.toolstransienttime"));
#endif


const double timeScale = 1.;

ToolsTransient::ToolsTransient(){
    
}

ToolsTransient::~ToolsTransient(){
    
}

//TPZFMatrix<REAL> ToolsTransient::InitialSolution(TPZGeoMesh *gmesh, TPZCompMesh * cmesh, int matId, int porder, REAL valsol){
//    
//    TPZAnalysis an(cmesh);
//	int nrs = an.Solution().Rows();
//    TPZVec<REAL> initsol(nrs,valsol);
//    int dim = cmesh->Dimension();
//    
//    TPZCompMesh  * cmesh_projL2 = CMeshProjectionL2(gmesh, dim, matId, porder, initsol);
//    TPZAnalysis anL2(cmesh_projL2);
//	TPZSkylineStructMatrix full(cmesh_projL2);
//	anL2.SetStructuralMatrix(full);
//	TPZStepSolver<REAL> step;
//	step.SetDirect(ELDLt);
//	anL2.SetSolver(step);
//	anL2.Run();
//    
//    TPZFMatrix<REAL> InitialSolution=anL2.Solution();
//    cmesh->LoadSolution(InitialSolution);
//    
//    return InitialSolution;
//}

//TPZCompMesh * ToolsTransient::CMeshProjectionL2(TPZGeoMesh *gmesh, int dim, int matId, int pOrder, TPZVec<STATE> &solini)
//{
//    /// criar materiais
//	TPZL2Projection *material;
//	material = new TPZL2Projection(matId, dim, 1, solini, pOrder);
//    
//    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//    cmesh->SetDimModel(dim);
//    TPZMaterial * mat(material);
//    cmesh->InsertMaterialObject(mat);
//    
//	cmesh->SetAllCreateFunctionsContinuous();
//    
//	cmesh->SetDefaultOrder(pOrder);
//    cmesh->SetDimModel(dim);
//	
//	//Ajuste da estrutura de dados computacional
//	cmesh->AutoBuild();
//    
//	return cmesh;
//}

void ToolsTransient::StiffMatrixLoadVec(TPZNLFluidStructure2d *mymaterial, TPZCompMesh* mphysics,
                                        TPZAnalysis *an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec)
{
	mymaterial->SetCurrentState();
    TPZFStructMatrix matsk(mphysics);

	an->SetStructuralMatrix(matsk);
	TPZStepSolver<REAL> step;

	step.SetDirect(ELU);
	an->SetSolver(step);
    
    an->Assemble();
	
	matK1 = an->StructMatrix(); //<<< essa bodega está retornando matriz 0x0 !!!

	fvec = an->Rhs();
}

void ToolsTransient::MassMatrix(TPZNLFluidStructure2d *mymaterial, TPZCompMesh *mphysics, TPZFMatrix<REAL> & Un)
{
    mymaterial->SetLastState();
	TPZSpStructMatrix matsp(mphysics);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
    matsp.CreateAssemble(Un,guiInterface);
}


void ToolsTransient::SolveSistTransient(REAL deltaT,REAL maxTime, TPZFMatrix<REAL> InitialSolution , TPZAnalysis *an, TPZNLFluidStructure2d * &mymaterial, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics)
{
    std::ofstream outPW("OUTPUT.txt");
    
    std::string outputfile;
	outputfile = "TransientSolution";
    
    std::stringstream outP, outW, outJ;
    
    outP << "Saida" << 0 << "={";
    SaidaMathPressao(meshvec, mphysics, outP);
    outP << "};\n";
    
    meshvec[0]->LoadReferences();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    PlotWIntegral(meshvec[0], outW, 0.);
    
	int nrows;
	nrows = an->Solution().Rows();
	TPZFMatrix<REAL> res_total(nrows,1,0.0);
    TPZFMatrix<REAL> chutenewton(meshvec[0]->Solution().Rows(),1,1.);
    meshvec[0]->LoadSolution(chutenewton);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    TPZFMatrix<REAL> SolIterK = mphysics->Solution();
	TPZFMatrix<REAL> matK;
	TPZFMatrix<REAL> fres(mphysics->NEquations(),1);
    TPZFMatrix<REAL> fmat(mphysics->NEquations(),1);
    fres.Zero();
    fmat.Zero();
    
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT;
    
    ////////// Calculo da Integral-J ///
    REAL radius = 0.5;
    ComputeKI(meshvec[0], radius, outJ);
    ////////////////////////////////////
    
    std::map<REAL,REAL> notUsedHere;
    REAL pfracMedio = 0.;
    
    std::ofstream saida("s.txt");
	while (TimeValue <= maxTime) //passo de tempo
	{
        outP << "Saida" << (int)(TimeValue/timeScale) << "={";

        fres.Zero();
        StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fres);
        int nr = matK.Rows();
        
        res_total = fres + fmat;
        REAL res = Norm(res_total);
        REAL tol = 1.e-8;
        int maxit = 15;
        int nit = 0;
        
        while(res > tol && nit < maxit) //itercao de Newton
        {
            an->Rhs() = res_total;
            an->Solve();
            an->LoadSolution(SolIterK + an->Solution());
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            
            SolIterK = an->Solution();
            
            fres.Zero();
            StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fres);
            res_total = fres + fmat;

            res = Norm(res_total);
            std::cout << res << std::endl;
            nit++;
        }
        
        fmat.Zero();
        MassMatrix(mymaterial, mphysics,fmat);
        
        std::map<REAL,REAL> notUsedHere;
        pfracMedio = PressaoMedia(meshvec, mphysics, notUsedHere);
        mymaterial->UpdateLeakoff(pfracMedio);
        
        ComputeKI(meshvec[0], radius, outJ, cent, TimeValue, false);
        
        if(res >= tol)
        {
            std::cout << cent << " , normRes = " << res << std::endl;
            //DebugStop();
        }
        if(nit >= maxit)
        {
            std::cout << cent << " , nitTot = " << nit << std::endl;
            //DebugStop();
        }
        
        {
            SaidaMathPressao(meshvec, mphysics, outP);
            outP << "};\n";
        }
        
        {
            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            PosProcessMult(meshvec,mphysics,an,plotfile);
        }

        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
        PlotWIntegral(meshvec[0], outW, TimeValue);
        
        InitialSolution = mphysics->Solution();
        cent++;
        TimeValue = cent*deltaT;
        
        //CheckConv(an->Solution() , an, mymaterial, meshvec, mphysics);
    }
    
    outW << "\nTrapArea[displ_]:=Block[{displSize,area,baseMin,baseMax,h},\n";
    outW << "displSize = Length[displ];\n";
    outW << "area = 0;\n";
    outW << "For[i = 1, i < displSize,\n";
    outW << "baseMin = displ[[i, 2]];\n";
    outW << "baseMax = displ[[i + 1, 2]];\n";
    outW << "h = displ[[i + 1, 1]] - displ[[i, 1]];\n";
    outW << "area += (baseMin+baseMax)*h/2;\n";
    outW << "i++;];\n";
    outW << "area];\n\n";
    outW << "Areas = {";
    int nsteps = maxTime/deltaT;
    for(int ig = 0; ig<=nsteps; ig++)
    {
        outW << "{" << ig*deltaT << ",TrapArea[displ" << (int)(ig*deltaT/timeScale) << "]}";
        if(ig != nsteps)
        {
            outW << ",";
        }
    }
    outW << "};\n\n";
    
    outW << "Qinj=" << fabs(mymaterial->Qinj()) << ";\n";
    outW << "Ttot = " << maxTime << ";\n";
    outW << "nsteps = " << cent-1 << ";\n";
    outW << "dt = Ttot/nsteps;\n\n";
    outW << "clnum = " << mymaterial->Cl()*100000 << "*10^-5;\n";
    outW << "pfracnum[t_] := press[t + dt/20];\n";
    outW << "SigConfnum = " << mymaterial->SigmaConf() << ";\n";
    outW << "Penum = " << mymaterial->Pe() << ";\n";
    outW << "Prefnum = " << mymaterial->Pref() << ";\n";
    outW << "vspnum = " << mymaterial->vsp()*100000 << "*10^-5;\n";
    outW << "Lfnum = " << mymaterial->Lf() << ";\n\n";
    outW << "vlini = 0;\n";
    outW << "(* FIM DOS INPUTS *)\n\n";
    
    outW << "(* KERNEL *)\n";
    
    outW << "vlF\\[Tau][cl_, pfrac_, SigConf_, Pe_, Pref_, vsp_, \\[Tau]_] :=Block[{clcorr, vlcomputed},\n";
    outW << "clcorr = N[cl Sqrt[(pfrac + SigConf - Pe)/Pref]];\n";
    outW << "vlcomputed = N[2 clcorr Sqrt[\\[Tau]] + vsp];\n";
    outW << "vlcomputed\n";
    outW << "];\n\n";
    
    outW << "qlFvlacumANDvlacumnext[vlacum_, dt_, cl_, pfrac_, SigConf_, Pe_,Pref_, vsp_] := Block[{clcorr, tstar, vlacumNext, qlcomputed},\n";
    outW << "clcorr = N[cl Sqrt[(pfrac + SigConf - Pe)/Pref]];\n";
    outW << "tstar = If[vlacum < vsp, 0, N[(vlacum - vsp)^2/(4 clcorr^2)]];\n";
    outW << "vlacumNext =vlF\\[Tau][cl, pfrac, SigConf, Pe, Pref, vsp, tstar + dt];\n";
    outW << "qlcomputed = N[(vlacumNext - vlacum)/dt];\n";
    outW << "{qlcomputed, vlacumNext, tstar, tstar + dt}\n";
    outW << "];\n\n";
    
    outW << "qlVec = {};\n";
    outW << "vlVec = {{0, 0}};\n";
    outW << "tstarsss = {};\n\n";
    
    outW << "pass = False;\n";
    outW << "vlnext = vlini;\n";
    outW << "For[step = 0, step < nsteps, step++,\n";
    outW << "QlactVlnext = qlFvlacumANDvlacumnext[vlnext, dt, clnum, pfracnum[step*dt],SigConfnum, Penum, Prefnum, vspnum];\n";
    outW << "qlact = QlactVlnext[[1]];\n";
    outW << "vlnext = QlactVlnext[[2]];\n";
    outW << "AppendTo[qlVec, {(step)*dt, qlact}];\n";
    outW << "AppendTo[vlVec, {(step + 1)*dt, vlnext}];\n";
    outW << "AppendTo[tstarsss, {QlactVlnext[[3]], QlactVlnext[[4]]}];\n";
    outW << "];\n";
    outW << "(* FIM DO KERNEL *)\n\n";
    
    outW << "(* OUTPUTS *)\n";
    outW << "plVLnum =ListPlot[vlVec, AxesOrigin -> {0,0}, Joined -> True,PlotRange -> All, Filling -> Axis, AxesLabel -> {\"T\", \"VL\"},PlotStyle -> {Red, Thickness[0.012]}];\n";
    outW << "plQLnum =ListPlot[qlVec, AxesOrigin -> {0,0}, Joined -> True,PlotRange -> All, PlotStyle -> {Red, Thickness[0.012]},Filling -> Axis, AxesLabel -> {\"T\", \"QL\"}];\n\n";
    
    outW << "Show[plQLnum]\n";
    outW << "Show[plVLnum]\n\n";
    
    outW << "WintegralPlot =ListPlot[Areas, AxesOrigin -> {0,0},PlotStyle -> {PointSize[0.015]},AxesLabel->{\"t\",\"V\"},Filling->Axis];\n";
    outW << "vInjTable = Table[{t,Qinj*t},{t,0,Ttot, dt}];\n";
    outW << "vInj[t_] = Interpolation[vInjTable, InterpolationOrder -> 0][t];\n";
    outW << "vfiltrado[t_] =Interpolation[vlVec, InterpolationOrder -> 0][t];\n";
    outW << "vf[t_]=vInj[t]-vfiltrado[t];\n";
    outW << "QfinalPlot = Plot[vf[t], {t, 0, Ttot}, PlotStyle -> Red,AxesOrigin->{0,0}];\n";
    outW << "Show[WintegralPlot,QfinalPlot]\n";
    outW << "(* FIM DOS OUTPUTS *)\n\n";
    
    //saida para mathematica
    outP << "SAIDAS={";
    for(int i = 1; i < cent; i++)
    {
        outP << "Saida" << (int)(i*deltaT/timeScale);
        if(i < cent-1)
        {
            outP << ",";
        }
    }
    outP << "};\n\n";
    
    outP << "minx=Min[Transpose[Flatten[SAIDAS,1]][[1]]];\n";
    outP << "maxx=Max[Transpose[Flatten[SAIDAS,1]][[1]]];\n";
    outP << "miny=Min[Transpose[Flatten[SAIDAS,1]][[2]]];\n";
    outP << "maxy=Max[Transpose[Flatten[SAIDAS,1]][[2]]];\n\n";
    
    outP << "Manipulate[ListPlot[SAIDAS[[n]],Joined->True,AxesOrigin->{0,0},PlotRange->{{0,maxx},{0,maxy}},AxesLabel->{\"pos\",\"p\"}],{n,1,Length[SAIDAS],1}]\n\n";
    outP << "TimePressVec = {};\n";
    outP << "For[pos = 1, pos <= Length[SAIDAS], pos++,\n";
    outP << "AppendTo[TimePressVec, {pos*" << deltaT << ",(SAIDAS[[pos]])[[Round[Length[SAIDAS[[pos]]]/2], 2]]}];\n];\n";
    outP << "ListPlot[TimePressVec, Joined -> True,AxesLabel->{\"t\",\"p\"},AxesOrigin->{0,0}]\n";
    outP << "press = Interpolation[TimePressVec,InterpolationOrder->0];\n";
    
    
    outJ << "};\n";
    outJ << "ListPlot[J, Joined -> True,AxesLabel->{\"t\",\"J\"},AxesOrigin->{0,0},Filling->Axis]\n";
    
    outPW << "(*** PRESSAO ***)\n" << outP.str() << "\n\n\n\n\n(*** W ***)\n" << outW.str() << "\n\n\n\n\n(*** J ***)\n" << outJ.str();
    outPW.close();
}

void ToolsTransient::ComputeKI(TPZCompMesh * elastMesh, REAL radius, std::stringstream & outFile, int cent, REAL TimeValue, bool firstCall)
{
    if(firstCall == true)
    {
        outFile << "J={";
    }
    else
    {
        REAL XcrackTip = -1.;
        TPZGeoMesh * gm = elastMesh->Reference();
        int bcfluxOutMat = -20;//ver materialId no arquivo main_elast.cpp
        for(int ell = 0; ell < gm->NElements(); ell++)
        {
            if(gm->ElementVec()[ell] && gm->ElementVec()[ell]->MaterialId() == bcfluxOutMat)
            {
                int nodeIndex = gm->ElementVec()[ell]->NodeIndex(0);
                XcrackTip = gm->NodeVec()[nodeIndex].Coord(0);
                break;
            }
        }
        if(XcrackTip < 1.E-3)
        {
            DebugStop();
        }
        TPZVec<REAL> Origin(3,0.);
        Origin[0] = XcrackTip;
        TPZVec<REAL> normalDirection(3,0.);
        normalDirection[2] = 1.;
        
        //////////// Computing pressure at middle point of crack length /////
        TPZVec<REAL> xx(3,0.), qsii(2,0.);
        xx[0] = XcrackTip/2.;
        int initialEl = 0;
        TPZGeoEl * geoEl = elastMesh->Reference()->FindElement(xx, qsii, initialEl, 2);
        if(!geoEl)
        {
            DebugStop();
        }
        TPZCompEl * compEl = geoEl->Reference();
        if(!compEl)
        {
            DebugStop();
        }
        TPZInterpolationSpace * intpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
        TPZMaterialData data;
        intpEl->InitMaterialData(data);
        intpEl->ComputeShape(qsii, data);
        intpEl->ComputeSolution(qsii, data);
        TPZElasticityMaterial * elast2D = dynamic_cast<TPZElasticityMaterial *>(compEl->Material());
        TPZVec<REAL> Solout(3);
        int var = 10;//Stress Tensor
        elast2D->Solution(data, var, Solout); 
        REAL pressure = -Solout[1];
        /////////////////////////////////////////////////////////////////////
        
        Path2D * Jpath = new Path2D(elastMesh, Origin, normalDirection, radius, pressure);
        JIntegral2D integralJ;
        integralJ.PushBackPath2D(Jpath);
        TPZVec<REAL> KI(2,0.);
        KI = integralJ.IntegratePath2D(0);
        
        if(cent != 1)
        {
            outFile << ",";
        }
        outFile << "{" << TimeValue << "," << KI[0] << "}";
    }
}

void ToolsTransient::PlotWIntegral(TPZCompMesh *cmesh, std::stringstream & outW, int solNum)
{
    TPZCompMeshReferred * cmeshref = dynamic_cast<TPZCompMeshReferred*>(cmesh);
    outW << "displ" << (int)(solNum/timeScale) << "={";
    int npts = 1;
    
    bool isFirstTime = true;
    for(int el = 0; el < cmeshref->NElements(); el++)
    {
        TPZCompEl *cel = cmeshref->ElementVec()[el];
        if(!cel) continue;
        TPZGeoEl * gel1D = cel->Reference();
        int crak1DMatId = 2;
        if(!gel1D || gel1D->HasSubElement() || gel1D->Dimension() != 1 || gel1D->MaterialId() != crak1DMatId)
        {
            continue;
        }
        
        {
            TPZVec<REAL> qsi1D(1,0.), qsi2D(2,0.), XX(3,0.);
            
            for(int p = -npts; p <= +npts; p++)
            {
                if(isFirstTime == false && p == -npts)
                {
                    continue;
                }
                qsi1D[0] = double(p)/npts;
                gel1D->X(qsi1D, XX);
                int inilIndex = 0;
                TPZGeoEl * gel2D = cmesh->Reference()->FindElement(XX, qsi2D, inilIndex, 2);
                
                TPZElasticityMaterial * elast2D = dynamic_cast<TPZElasticityMaterial *>(gel2D->Reference()->Material());
                
                TPZVec<REAL> Solout(3);
                
                int var = 0;
                TPZReducedSpace * intpEl = dynamic_cast<TPZReducedSpace *>(gel2D->Reference());
                TPZMaterialData data;
                intpEl->InitMaterialData(data);
                
                intpEl->ComputeShape(qsi2D, data);
                intpEl->ComputeSolution(qsi2D, data);
                elast2D->Solution(data, var, Solout);
                
                REAL posX = XX[0];
                REAL posY = 2.*(Solout[1]);
                if(fabs(posX) < 1.E-5)
                {
                    posX = 0;
                }
                if(fabs(posY) < 1.E-5)
                {
                    posY = 0;
                }
                if(isFirstTime)
                {
                    isFirstTime = false;
                }
                else
                {
                    outW << ",";
                }
                outW << "{" << posX << "," << posY << "}";
            }
        }//<<<
    }
    outW << "};\n";
}

void ToolsTransient::SaidaMathPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, std::stringstream & outP)
{
    std::map<REAL,REAL> pos_pressure;
    PressaoMedia(meshvec,mphysics,pos_pressure);
    
    int posCount = 0;
    std::map<REAL,REAL>::iterator it, itaux;
    for(it = pos_pressure.begin(); it != pos_pressure.end(); it++)
    {
        itaux = it;
        itaux++;
        REAL pos = it->first;
        REAL press = it->second;
        outP << "{" << pos << "," << press << "}";
        if(itaux != pos_pressure.end())
        {
            outP << ",";
        }
        posCount++;
    }
}

REAL ToolsTransient::PressaoMedia(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh * mphysics, std::map<REAL,REAL> & pos_pressure)
{
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);

    int pressureMatId = 2;//ver arquivo main_elast (pois contem os materialIds)
    
    for(int i = 0;  i< meshvec[1]->ElementVec().NElements(); i++)
    {
        TPZCompEl * cel = meshvec[1]->ElementVec()[i];
        if(cel->Reference()->MaterialId() != pressureMatId)
        {
            continue;
        }
        TPZInterpolatedElement * sp = dynamic_cast <TPZInterpolatedElement*> (cel);
        if(!sp) continue;
        TPZVec<REAL> qsi(1,0.), Xqsi(3,0.);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        
        qsi[0] = -1.;
        sp->ComputeShape(qsi, data);
        sp->ComputeSolution(qsi, data);
        TPZVec<REAL> SolP = data.sol[0];
        cel->Reference()->X(qsi,Xqsi);
        REAL pos = Xqsi[0];
        REAL press = data.sol[0][0];
        pos_pressure[pos] = press;
    }
    
    REAL middlePress = -1.;
    int sz = pos_pressure.size();
    int middle = sz/2, posCount = 0;
    
    if(middle == 0)
    {
        DebugStop();
    }

    std::map<REAL,REAL>::iterator it, itaux;
    for(it = pos_pressure.begin(); it != pos_pressure.end(); it++)
    {
        itaux = it;
        itaux++;
        REAL press = it->second;
        if(posCount == middle)
        {
            middlePress = press;
        }
        posCount++;
    }
    
    return middlePress;
}

void ToolsTransient::PosProcessMult(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis *an, std::string plotfile)
{
    //TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(5), vecnames(1);
	
	scalnames[0] = "DisplacementX";
	scalnames[1] = "DisplacementY";
    scalnames[2] = "SigmaX";
    scalnames[3] = "SigmaY";
    scalnames[4] = "Pressure";
    vecnames[0] = "Displacement";
	
	const int dim = 2;
	int div =0;
	an->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an->PostProcess(div);
}

//TPZFMatrix<REAL> ToolsTransient::SetSolution(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder, int matId, REAL valIni){
//    
//    TPZAnalysis an(cmesh);
//    int dim = cmesh->Dimension();
//	int nrs = an.Solution().Rows();
//    TPZVec<REAL> loadvec(nrs,valIni);
//    TPZCompMesh  * cmesh_projL2 = ToolsTransient::CMeshProjectionL2(gmesh,dim, matId, pOrder, loadvec);
//    TPZAnalysis anL2(cmesh_projL2);
//    
//    //Solve
//	TPZSkylineStructMatrix full(cmesh_projL2);
//	anL2.SetStructuralMatrix(full);
//	TPZStepSolver<REAL> step;
//	step.SetDirect(ELDLt);
//	anL2.SetSolver(step);
//	anL2.Run();
//    
//    TPZFMatrix<REAL> InitSol;
//    InitSol = anL2.Solution();
//    
//    return InitSol;
//}


void ToolsTransient::CheckConv(TPZFMatrix<REAL> actQsi , TPZAnalysis *an, TPZNLFluidStructure2d * &mymaterial,
                               TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics)
{
	int neq = mphysics->NEquations();
    int nsteps = 5;
    
    TPZFMatrix<REAL> fxIni(neq,1);
    TPZFMatrix<REAL> fxAprox(neq,1);
    TPZFMatrix<REAL> fxExato(neq,1);
    
    TPZFMatrix<REAL> errorVec(neq,1,0.);
	TPZFMatrix<REAL> errorNorm(nsteps,1,0.);
    
    TPZFMatrix<REAL> fLIni;
    TPZFMatrix<REAL> fLtemp;
    
    TPZFMatrix<REAL> dFx(neq,1);
    
    TPZFMatrix<REAL> qsiIni = actQsi;

    fxIni.Zero();
    fLIni.Zero();
    StiffMatrixLoadVec(mymaterial, mphysics, an, fLIni, fxIni);
    
    if(fLIni.Rows() != neq || fLIni.Cols() != neq)
    {
        DebugStop();
    }
    
    TPZVec<REAL> deltaQsi(neq,0.1);
    double alpha;
    
    for(int i = 0; i < nsteps; i++)
    {
        alpha = i/10.;
        
        ///Fx aproximado
        dFx.Zero();
        for(int r = 0; r < neq; r++)
        {
            for(int c = 0; c < neq; c++)
            {
                dFx(r,0) +=  fLIni.GetVal(r,c) * (alpha * deltaQsi[c]);
            }
        }
        fxAprox = fxIni + dFx;
        
        ///Fx exato
        for(int r = 0; r < neq; r++)
        {
            actQsi(r,0) = qsiIni(r,0) + (alpha * deltaQsi[r]);
        }
        an->LoadSolution(actQsi);
        fxExato.Zero();
        fLtemp.Zero();
        StiffMatrixLoadVec(mymaterial, mphysics, an, fLtemp, fxExato);
        
        ///Erro
        errorVec.Zero();
        for(int r = 0; r < neq; r++)
        {
            errorVec(r,0) = fxExato(r,0) - fxAprox(r,0);
        }
        
        ///Norma do erro
        double XDiffNorm = Norm(errorVec);
        errorNorm(i,0) = XDiffNorm;
    }
    
    std::cout << "Convergence Order:\n";
    for(int j = 2; j < nsteps; j++)
    {
        std::cout << ( log(errorNorm(1,0)) - log(errorNorm(j,0)) )/( errorNorm(0.1) - errorNorm(j/10.) ) << "\n";
    }
}