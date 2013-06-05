//
//  File.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/31/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//


#include "pzfstrmatrix.h"
#include "toolstransienttime.h"

#include "TPZVTKGeoMesh.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
//#endif



int main(int argc, char *argv[])
{
    //Propagation criterion
    REAL Lx = 400.;
    REAL Ly = 400.;
    REAL Lf = 50.;
    REAL Hf = 1.;
    REAL Young = 3.9E4;
    REAL Poiss = 0.25;
    REAL Fx = 0.;
    REAL Fy = 0.;
    REAL Visc = 0.001E-6;
    REAL SigN = 61.5;
    REAL Qinj  = -0.1/Hf;
    REAL Ttot = 50.;
    REAL Nsteps = 20.;
    REAL Cl = 0.005;
    REAL Pe = 10.;
    REAL SigmaConf = 11.;
    REAL Pref = 60000.;
    REAL vsp = 0.001;
    REAL KIc = 300.;
    int p = 2;
    
    ToolsTransient ToolTrans(p, Lx, Ly, Lf, Hf, Young, Poiss, Fx, Fy, Visc, SigN, Qinj, Ttot, Nsteps, Cl, Pe, SigmaConf, Pref, vsp, KIc);
	//primeira malha
	
	// geometric mesh (initial)
    TPZGeoMesh * gmesh = ToolTrans.Mesh2D(10.);
    

    /** >>> Resolvendo um problema modelo de elastica linear para utilizar a
    solucao como espaco de aproximacao do problema nao linear (acoplado) */
	TPZCompMesh * cmesh_elast = ToolTrans.CMeshElastic(gmesh);
    TPZAnalysis an1(cmesh_elast);
    ToolTrans.MySolve(an1, cmesh_elast);
    
    /** >>> Passando a malha computacional jah processada do problema modelo elastico para a malha do problema elastico que serah acoplado.
    Esta malha que foi passada servirah como espaco de aproximacao da malha do problema elastico que serah acoplado */
    TPZCompMeshReferred * cmesh_referred = ToolTrans.CMeshReduced(gmesh, cmesh_elast);
    cmesh_referred->ComputeNodElCon();
    TPZFStructMatrix fstr(cmesh_referred);
    TPZFMatrix<STATE> rhs;//(1);
    TPZAutoPointer<TPZMatrix<STATE> > strmat = fstr.CreateAssemble(rhs,NULL);
	
    /** >>> Criando a malha computacional para o problema de fluxo de fluido */
    TPZCompMesh * cmesh_pressure = ToolTrans.CMeshPressure(gmesh);
    
    /** >>> Problema multifisico : acoplamento das 2 malhas computacionais anteriores (elastica e fluxo) */
    TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh_referred;
	meshvec[1] = cmesh_pressure;
    
    gmesh->ResetReference();
	TPZNLFluidStructure2d * mymaterial = NULL;
    
    /** >>> Serah utilizado o material criado (pznlfluidstructure2D do Agnaldo) */
    TPZCompMesh * mphysics = ToolTrans.MalhaCompMultphysics(gmesh, meshvec, mymaterial);
    mphysics->SetDefaultOrder(p);
    
    REAL deltaT = Ttot/Nsteps;
    mymaterial->SetTimeStep(deltaT);
    REAL maxTime = Ttot;
    
    /** >>> Metodo de resolucao de problema transient */
    TPZAnalysis *an = new TPZAnalysis(mphysics);
    TPZFMatrix<REAL> InitialSolution = an->Solution();
    ToolTrans.SolveSistTransient(deltaT, maxTime, InitialSolution, an, mymaterial, meshvec, mphysics);
    
    return 0;
}


