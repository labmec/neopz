//Classes utilit�ias
#include "pzvec.h"
#include "pzstack.h"
#include "TPZTimer.h"

//Classes Geom�ricas
#include "pzgmesh.h"
#include "pzgeoel.h"

//Classes Computacionais
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzintel.h"

#include "pzpoisson3d.h"
#include "pzbndcond.h"

//Matrizes de armazenamento e estruturais
#include "pzfmatrix.h"
//#include "pzskylmat.h"
#include "pzysmp.h"
#include "pzskylstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"

//Solver
#include "pzstepsolver.h"

//An�ise
#include "pzanalysis.h"


//Pos-processamento
#include "pzvisualmatrix.h"
#include "TPZVTKGeoMesh.h"

// logging
#include "pzlog.h"

//Bibliotecas std, math etc
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>


#include "pzgengrid.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.conductivity"));
#endif

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

using namespace std;

///Funcao para a criacao da malha
TPZAutoPointer<TPZCompMesh> CriarMalha(REAL L, REAL deltat, REAL height);

///Funcao que calcula a condutividade da configuracao
REAL Conductivity(REAL L, REAL delta, REAL height);

///Refine a malha em torno da singularidade
void DivideAround(TPZAutoPointer<TPZGeoMesh> gmesh, int bc, int nrefine);

///Funcao para ajustar a ordem p de aproximacao
void PRefineAround(TPZAutoPointer<TPZCompMesh> cmesh, int bc, int nrefine, int pstart);

///Funcao que calcula o fluxo de saida
REAL FluxoSaida(TPZAutoPointer<TPZCompMesh> cmesh, TPZFMatrix<REAL> &solution, int bc);

///Parametro para a funcao
static REAL g_Length = 0.;
///Insercao de uma condicao de contorno dada por uma funcao
void forcingfunction(const TPZVec<REAL> &ponto, TPZVec<REAL> &force);

int main()
{
    InitializePZLOG();
    
    REAL minL = 1.;
    REAL maxL = 4.;
    REAL mindelt = 0.001;
    REAL maxdelt = 0.01;
    REAL minheight = 40;
    REAL maxheight = 60;
    int numL = 4;
    int numdelta = 6;
    int numheight = 2;
    int iL, id, ih;
    std::ofstream cond("condutividade.txt",ios::app);
    //cond << "L delta delta/L height conductivity\n";
    for(iL=0; iL<numL; iL++)
    {
        REAL L = minL + iL*(maxL-minL)/numL;
        for(id = 0; id<numdelta; id++)
        {
            REAL delta = mindelt + id*(maxdelt-mindelt)/numdelta;
            for(ih = 0; ih<numheight; ih++)
            {
                REAL height = minheight + ih*(maxheight-minheight)/numheight;
                REAL conductivity = Conductivity(L,delta,height);
                cond << L << " " << delta << " " << delta/L << " " << height << " " << conductivity << std::endl;
                std::cout << L << " " << delta << " " << delta/L << " " << height << " " << conductivity << std::endl;
            }
        }
    }
    return 0;
}

TPZAutoPointer<TPZCompMesh> CriarMalha(REAL L, REAL delta, REAL height)
{
//    TPZGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numl = 1, REAL rot = 0.5);
    int nelem = 32;
    int nrefine = 5;
    int defaultporder = 2;
    TPZVec<int> nx(2,nelem);
    TPZVec<REAL> x0(3,0.),x1(3,L);
    x1[1] = height;
    x1[2] = 0.;
    TPZGenGrid grid(nx,x0,x1);
    TPZManVector<REAL> minsizes(2,delta);
    TPZManVector<REAL,2> progression(2,1.);
    grid.ComputeGeometricProgression(minsizes, progression);
    progression[0] = 1./progression[0];
    grid.SetGeometricProgression(progression);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    grid.Read(gmesh);
    grid.SetBC(gmesh,1,-4);
    grid.SetBC(gmesh,3,-1);
    TPZVec<REAL> xnext(x0);
    xnext[0] = x1[0]-delta;
    grid.SetBC(gmesh,x0,xnext,-2);
    grid.SetPointBC(gmesh,xnext,-3);
    DivideAround(gmesh, -3, nrefine);
    std::ofstream output("conductivitygmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(),output,true);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    TPZMatPoisson3d *poiss = new TPZMatPoisson3d(1,2);
    TPZAutoPointer<TPZMaterial> autopoiss(poiss);
    TPZFMatrix<REAL> val1(1,1,0.),val2(1,1,0.);
    cmesh->InsertMaterialObject(autopoiss);
    TPZBndCond *bc = new TPZBndCond(autopoiss, -1, 0, val1, val2);
    TPZAutoPointer<TPZMaterial> bcauto(bc);
    cmesh->InsertMaterialObject(bcauto);
    TPZBndCond *bc2 = new TPZBndCond(autopoiss, -2, 0, val1, val2);
    g_Length = L;
    TPZAutoPointer<TPZFunction> force = new TPZDummyFunction(forcingfunction);
    bc2->SetForcingFunction(force);
    TPZAutoPointer<TPZMaterial> bcauto2(bc2);
    bcauto2->SetForcingFunction(force);
    cmesh->InsertMaterialObject(bcauto2);

    TPZBndCond *bc4 = new TPZBndCond(autopoiss, -4, 0, val1, val2);
    TPZAutoPointer<TPZMaterial> bcauto4(bc4);
    cmesh->InsertMaterialObject(bcauto4);
    
    cmesh->SetDefaultOrder(defaultporder);
    cmesh->AutoBuild();
    PRefineAround(cmesh, -3, nrefine, defaultporder);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    return cmesh;
}

///Refine a malha em torno da singularidade
void DivideAround(TPZAutoPointer<TPZGeoMesh> gmesh, int bc, int nrefine)
{
    int iel;
    int nelem = gmesh->NElements();
    TPZGeoEl *gel;
    for(iel = 0; iel<nelem; iel++)
    {
        gel = gmesh->ElementVec()[iel];
        if(!gel) continue;
        if(gel->MaterialId() == bc)
        {
            break;
        }
    }
    int iref = 0;
    for( iref = 0; iref < nrefine; iref++)
    {
        TPZStack<TPZGeoElSide> connected;
        TPZGeoElSide gels(gel,0);
        TPZGeoElSide neigh = gels.Neighbour();
        while(neigh != gels)
        {
            connected.Push(neigh);
            neigh = neigh.Neighbour();
        }
        int nelcon = connected.NElements();
        int ic;
        for(ic = 0; ic<nelcon; ic++)
        {
            if(! connected[ic].Element()->HasSubElement())
            {
                TPZStack<TPZGeoEl *> subels;
                connected[ic].Element()->Divide(subels);
            }
        }
    }
}

///Funcao para ajustar a ordem p de aproximacao
void PRefineAround(TPZAutoPointer<TPZCompMesh> cmesh, int bc, int nrefine, int pstart)
{
    int iel;
    TPZGeoMesh * gmesh = cmesh->Reference();
    int nelem = gmesh->NElements();
    TPZGeoEl *gel;
    for(iel = 0; iel<nelem; iel++)
    {
        gel = gmesh->ElementVec()[iel];
        if(!gel) continue;
        if(gel->MaterialId() == bc)
        {
            break;
        }
    }
    int iref = 0;
    for( iref = 0; iref < nrefine; iref++)
    {
        int nel = cmesh->NElements();
        int el;
        for( el = 0; el<nel; el++)
        {
            TPZCompEl *cel = cmesh->ElementVec()[el];
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
            if(!intel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Level() == iref+1)
            {
                intel->PRefine(pstart+iref+1);
            }
        }
    }
    TPZStack<TPZGeoElSide> connected;
    TPZGeoElSide gels(gel,0);
    TPZGeoElSide neigh = gels.Neighbour();
    while(neigh != gels)
    {
        connected.Push(neigh);
        neigh = neigh.Neighbour();
    }
    int nelcon = connected.NElements();
    int ic;
    for(ic = 0; ic<nelcon; ic++)
    {
        if(connected[ic].Element()->Reference())
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (connected[ic].Element()->Reference());
            if(!intel) continue;

            intel->PRefine(pstart);
        }
    }
    cmesh->ExpandSolution();
}


///Funcao que calcula a condutividade da configuracao
REAL Conductivity(REAL L, REAL delta, REAL height)
{
    TPZAutoPointer<TPZCompMesh> cmesh = CriarMalha(L,delta,height);
    TPZAnalysis an(cmesh);
    TPZSkylineStructMatrix strskyl(cmesh);
    strskyl.SetNumThreads(2);
    an.SetStructuralMatrix(strskyl);
    TPZStepSolver<REAL> solve;
    solve.SetDirect(ECholesky);
    an.SetSolver(solve);
    std::cout << "before running\n";
    std::cout << "neq " << cmesh->NEquations() << std::endl;
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    an.Run();
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
    std::cout << "t1 " << t1 << " t2 " << t2 << " elapse " << t2-t1 << "finished\n";
#endif
    
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("state");
    std::stringstream sout;
    sout << "conductivity" << L << "_" << delta << "_" << height << ".vtk";
    an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
    an.PostProcess(2);
    
    REAL fluxo = FluxoSaida(cmesh,an.Solution(),-4);
    
    // mostrar a estrutura da matriz graficamente
    TPZFMatrix<REAL> vismat;
    cmesh->ComputeFillIn(100,vismat);
    VisualMatrixVTK(vismat,"matrixstruct.vtk");
    
    return fluxo;
    
}

///Funcao que calcula o fluxo de saida
REAL FluxoSaida(TPZAutoPointer<TPZCompMesh> cmesh, TPZFMatrix<REAL> &solution, int bc)
{
    TPZSkylineStructMatrix strskyl(cmesh);
    strskyl.SetNumThreads(2);
    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;

    std::set<int> matids;
    matids.insert(1);
    strskyl.SetMaterialIds(matids);
    TPZFMatrix<REAL> rhs;
    TPZAutoPointer<TPZMatrix<REAL> > mat = strskyl.CreateAssemble(rhs, gui);

    // identify equations associated with linear shape functions of elements with materialid bc
    std::set<int> cornereqs;
    int nel = cmesh->NElements();
    int el;
    for( el = 0; el<nel; el++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if(!intel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if(gel->MaterialId() == bc)
        {
            int nc = intel->NCornerConnects();
            int ic;
            for (ic = 0; ic<nc; ic++)
            {
                TPZConnect &c = intel->Connect(ic);
                if(c.HasDependency()) continue;
                
                int seqnum = c.SequenceNumber();
                int eq = cmesh->Block().Position(seqnum);
                cornereqs.insert(eq);
            }
        }
    }
    mat->Multiply(solution,rhs);
    
    std::cout << "number of boundary nodes " << cornereqs.size() << std::endl;
    
    REAL fluxo=0.;
    std::set<int>::iterator it;
    for(it=cornereqs.begin(); it!= cornereqs.end(); it++)
    {
        fluxo += rhs(*it,0);
    }
    return fluxo;

}



void forcingfunction(const TPZVec<REAL> &p, TPZVec<REAL> &f) 
{
    f[0] = p[0]/g_Length;
}
