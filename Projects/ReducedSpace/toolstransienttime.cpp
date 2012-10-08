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

#include <fstream>
#include <sstream>

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.toolstransienttime"));
#endif


ToolsTransient::ToolsTransient(){
    
}

ToolsTransient::~ToolsTransient(){
    
}

#include "pzfstrmatrix.h"
void ToolsTransient::StiffMatrixLoadVec(TPZNLFluidStructure2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec){
    
	mymaterial->SetCurrentState();
   // TPZSpStructMatrix matsk(mphysics);
    TPZFStructMatrix matsk(mphysics);
   // TPZSkylineStructMatrix matsk(mphysics);
	an.SetStructuralMatrix(matsk); 
	TPZStepSolver<REAL> step; 
	//step.SetDirect(ELDLt);
	step.SetDirect(ELU);
	an.SetSolver(step); 
	an.Run(); 
	
	matK1 = an.StructMatrix();
	fvec = an.Rhs();
    
}

TPZAutoPointer <TPZMatrix<REAL> > ToolsTransient::MassMatrix(TPZNLFluidStructure2d *mymaterial, TPZCompMesh *mphysics){
    
    mymaterial->SetLastState();
    //TPZSkylineStructMatrix matsp(mphysics);
	TPZSpStructMatrix matsp(mphysics);
    
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<REAL> Un;
    TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}


TPZCompMesh * ToolsTransient::CMeshProjectionL2(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 1;
    int matId2 = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(matId2, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
	cmesh->SetAllCreateFunctionsContinuous();
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
	return cmesh;
}

std::ofstream outfile("SaidaPressao.nb");
void ToolsTransient::SolveSistTransient(REAL deltaT,REAL maxTime, TPZNLFluidStructure2d * &mymaterial, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics)
{
    TPZAnalysis an(mphysics);
	TPZFMatrix<REAL> Initialsolution = an.Solution();
    an.LoadSolution(Initialsolution);
   // mphysics->LoadSolution(Initialsolution);
    Initialsolution.Print();
    
    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<REAL> > matM = MassMatrix(mymaterial, mphysics);
    
    outfile << "Saida" << 0 << "={";
    SaidaMathPressao(meshvec, mphysics);
    outfile << "};\n";
    
    
    //Criando matriz de rigidez (tangente) matK e vetor de carga (residuo) 
	TPZFMatrix<REAL> matK;	
	TPZFMatrix<REAL> fres;
    //an.Solution().Zero();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fres);
    
    
	int nrows;
	nrows = matK.Rows();
	TPZFMatrix<REAL> res_total(nrows,1,0.0);
	TPZFMatrix<REAL> Mass_X_SolTimeN(nrows,1,0.0);
	TPZFMatrix<REAL> SolTimeN = Initialsolution;
    TPZFMatrix<REAL> SolIterK = SolTimeN;
	
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT; 
	while (TimeValue <= maxTime)
	{	
        outfile << "Saida" << cent << "={";
        
		// This time solution i for Transient Analytic Solution
        
        SolTimeN.Print();
		matM->Multiply(SolTimeN,Mass_X_SolTimeN);
        Mass_X_SolTimeN.Print();
		
        for(int i=0; i<5; i++){
            
            res_total = fres + Mass_X_SolTimeN;
            res_total.Print();
            
            an.Rhs() = res_total;
            an.Solve();
            an.Solution().Print();
            an.LoadSolution(an.Solution() + SolIterK);
            an.Solution().Print();
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            
            SolIterK = an.Solution();
            StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fres);
        }
        
        SaidaMathPressao(meshvec, mphysics);
        outfile << "};\n";
        
        cent++;
        TimeValue = cent*deltaT;
        SolTimeN = SolIterK;

	}
    
    //saida para mathematica
    outfile << "SAIDAS={";
    for(int i = 1; i < cent; i++)
    {
        outfile << "Saida" << i;
        if(i < cent-1)
        {
            outfile << ",";
        }
    }
    outfile << "};\n\n";
    
    outfile << "minx=Min[Transpose[Flatten[SAIDAS,1]][[1]]];\n";
    outfile << "maxx=Max[Transpose[Flatten[SAIDAS,1]][[1]]];\n";
    outfile << "miny=Min[Transpose[Flatten[SAIDAS,1]][[2]]];\n";
    outfile << "maxy=Max[Transpose[Flatten[SAIDAS,1]][[2]]];\n\n";

    outfile << "Manipulate[ListPlot[SAIDAS[[n]],Joined->True,AxesOrigin->{0,0},PlotRange->{{minx,maxx},{miny,maxy}}],{n,1,Length[SAIDAS],1}]\n";
    
    outfile.close();
}

void ToolsTransient::SaidaMathPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics)
{
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    
    std::map<REAL,REAL> time_pressure;
    
    for(int i = 0;  i< meshvec[1]->ElementVec().NElements(); i++)
    {
        TPZCompEl * cel =meshvec[1]->ElementVec()[i];
        TPZInterpolatedElement * sp = dynamic_cast <TPZInterpolatedElement*>(cel);
        if(!sp) continue;
        TPZVec<REAL> qsi(1,0.),out(3,0.);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        
        for(int j = 0; j < 1; j++)
        {
            qsi[0] = -1.;//+2.*i/10.;
            sp->ComputeShape(qsi, data);
            sp->ComputeSolution(qsi, data);
            TPZVec<REAL> SolP = data.sol[0]; 
            cel->Reference()->X(qsi,out);
            REAL pos = out[0];
            REAL press = data.sol[0][0];
            time_pressure[pos] = press;
        }
        if(out[0] > 50.) continue;
    }
    
    
    std::map<REAL,REAL>::iterator it, itaux;
    for(it = time_pressure.begin(); it != time_pressure.end(); it++)
    {
        itaux = it;
        itaux++;
        REAL pos = it->first;
        REAL press = it->second;
        outfile << "{" << pos << "," << press << "}";
        if(itaux != time_pressure.end())
        {
            outfile << ",";
        }
    }
}
