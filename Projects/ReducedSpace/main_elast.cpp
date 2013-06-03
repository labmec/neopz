//
//  File.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/31/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//


#include "pzfstrmatrix.h"
#include "toolstransienttime.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
#endif


    //Dimensions:
    const REAL Lx = 400.;
    const REAL Ly = 400.;
    const REAL Lf = 50.;
    const REAL Hf = 1.;
    //Elastic properties:
    const REAL ED = 3.9E4;//MPa
    const REAL nu = 0.25;
    const REAL fx = 0.;
    const REAL fy = 0.;
    //Fluid property:
    const REAL visc = 0.001E-6;//MPa.s
    //BCs:
    const REAL sigN = 61.5;/// <<< sigma.n no problema elastico que servira de espaco de aproximacao para o elastico multifisico
    const REAL Qinj  = -0.1/Hf;///vazao de 1 asa de fratura dividido pela altura da fratura
    //time:
    const REAL Ttot = 50.;//sem bug -> Ttot = 50.
    const REAL nsteps = 20.;//sem bug -> nsteps = 20.
    const REAL deltaT = Ttot/nsteps;
    //Leakoff:
    const REAL Cl = 0.005;
    const REAL Pe = 10.;//MPa
    const REAL SigmaConf = 11.;//MPa
    const REAL Pref = 60000.;//MPa
    const REAL vsp = 0.001;
    //Propagation criterion
    const REAL KIc = 300.;
    //



int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
#ifdef LOG4CXX
//    std::stringstream str;
//    str << "TESTANDO!!!";
//    LOGPZ_DEBUG(logger,str.str());
#endif
    
    int p = 2;
	//primeira malha
	
	// geometric mesh (initial)
    TPZGeoMesh * gmesh = ToolsTransient::Mesh2D(Lf, Lx, Ly, 10.);
    
    //computational mesh elastic
/** >>> Resolvendo um problema modelo de elastica linear para utilizar a solucao 
    como espaco de aproximacao do problema nao linear (acoplado) */
	TPZCompMesh * cmesh_elast = ToolsTransient::CMeshElastic(gmesh, p, ED, nu, sigN);
    TPZAnalysis an1(cmesh_elast);
    ToolsTransient::MySolve(an1, cmesh_elast);
    
/** >>> Passando a malha computacional jah processada do problema modelo elastico para a malha do problema elastico que serah acoplado.
    Esta malha que foi passada servirah como espaco de aproximacao da malha do problema elastico que serah acoplado */
    TPZCompMeshReferred * cmesh_referred = ToolsTransient::CMeshReduced(gmesh, cmesh_elast, p, ED, nu);
    cmesh_referred->ComputeNodElCon();
    TPZFStructMatrix fstr(cmesh_referred);
    TPZFMatrix<STATE> rhs;//(1);
    TPZAutoPointer<TPZMatrix<STATE> > strmat = fstr.CreateAssemble(rhs,NULL);
	
/** >>> Criando a malha computacional para o problema de fluxo de fluido */
    TPZCompMesh * cmesh_pressure = ToolsTransient::CMeshPressure(gmesh, p, Qinj);
    
/** >>> Problema multifisico : acoplamento das 2 malhas computacionais anteriores (elastica e fluxo) */
    //multiphysic mesh
    TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh_referred;
	meshvec[1] = cmesh_pressure;
    
    gmesh->ResetReference();
	TPZNLFluidStructure2d * mymaterial = NULL;
/** >>> Serah utilizado o material criado (pznlfluidstructure2D do Agnaldo) */
    TPZCompMesh * mphysics = ToolsTransient::MalhaCompMultphysics(gmesh, meshvec, mymaterial, ED, nu, fx, fy, Hf, Lf, visc,
                                                                  Qinj, Cl, Pe, SigmaConf, sigN, Pref, vsp, KIc, Lx);
    mphysics->SetDefaultOrder(p);
    
    mymaterial->SetTimeStep(deltaT);
    REAL maxTime = Ttot;
    
/** >>> Metodo de resolucao de problema transient */
    TPZAnalysis *an = new TPZAnalysis(mphysics);
    TPZFMatrix<REAL> InitialSolution = an->Solution();
    ToolsTransient::SolveSistTransient(deltaT, maxTime, InitialSolution, an, mymaterial, meshvec, mphysics);
    
    return 0;
}


