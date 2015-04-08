/**
 * @file
 * @brief Tests for hdiv pyramid
 * @author Nathan Shauer
 * @since 2016
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>

#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeotetrahedra.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZTimer.h"

#include <sys/time.h>

#include "run_stats_table.h"
#ifdef USING_TBB
#include <tbb/tbb.h>
#endif



#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.pyramtests"));
#endif


using namespace std;

TPZGeoMesh *MalhaCubo(string &projectpath, const int &nref);
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);
void InsertElasticityCubo(TPZCompMesh *mesh);

struct TTimer {
    struct timeval fini, fend;
    REAL fsec;
    
    TTimer(){
        fsec = 0.;
    }
    
    ~TTimer(){}
    
    void start(){
        gettimeofday(&fini, NULL);
    }
    
    void stop(){
        gettimeofday(&fend, NULL);
        fsec = fend.tv_sec - fini.tv_sec;
        fsec += (fend.tv_usec - fini.tv_usec)/1000000.;
    }
    
    REAL seconds(){
        return fsec;
    }
    
};


int main(int argc, char *argv[])
{
    
    TTimer tref;
    tref.start();
    gRefDBase.InitializeUniformRefPattern(ETetraedro);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EPiramide);
    tref.stop();
    cout << "\ntempo de refpattern = " << tref.seconds() << endl;
    
    string projectpath = "/Projects/PyramidHdivTests/";

#ifdef LOG4CXX
        std::string dirname = PZSOURCEDIR;
        std::string FileName = dirname;
        FileName = dirname + projectpath;
        FileName += "pyramlogfile.cfg";
        InitializePZLOG(FileName);
#endif
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream str;
        str << "\nRodando testes de piramede Hdiv" << std::endl;
        LOGPZ_DEBUG(logger,str.str())
    }
#endif
    
    // Parametros
    const int dim = 3;
    int plevel = 3;
    int numthreads = 8;
    int nref = 2;
    if (argc == 1) {
        cout << "\nATENCAO: voce nao passou argumentos, rodando c/ parametros hardcode!";
    }
    else if (argc == 4) {
        nref = atoi(argv[1]);
        plevel = atoi(argv[2]);
        numthreads = atoi(argv[3]);
        cout << "\nRodando com:" << endl;
        cout << "nref = " << nref << endl;
        cout << "plevel = " << plevel << endl;
        cout << "numthreads = " << numthreads << endl;
    }
    else{
        cout << "\nERRO - Num de argumento nao especificado" << endl;
        DebugStop();
    }
    
#ifdef USING_TBB
    if(numthreads == 0)
    {
        numthreads = 1;
    }
    tbb::task_scheduler_init init(numthreads);
#endif
    
    // Malha Geometrica
    cout << "\nCriando a gmesh... "; cout.flush();
    TTimer timer;
    timer.start();
    TPZGeoMesh *gmesh = MalhaCubo(projectpath,nref);
    if (0) {
        std::ofstream out("Cube.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    }
    timer.stop();
    cout << timer.seconds() << " s" << endl;

    // Malha computacional
    cout << "\nCriando a cmesh... "; cout.flush();
    timer.start();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(plevel);
    InsertElasticityCubo(cmesh);
    cmesh->AutoBuild();
    timer.stop();
    cout << timer.seconds() << " s" << endl;
    std::cout << "Numero de equacoes = " << cmesh->NEquations() << std::endl;
    
    // Analysis
    cout << "\nCriando o analysis... "; cout.flush();
    timer.start();
#define NOORDER
#ifdef NOORDER
    TPZAnalysis an(cmesh,false);
#else
    TPZAnalysis an(cmesh);
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    TPZSkylineStructMatrix skyl(cmesh);
    skyl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(skyl);
    an.SetSolver(step);
    timer.stop();
    cout << timer.seconds() << " s" << endl;
    
    // Resolvendo
    cout << "\nComecando o assemble... "; cout.flush();
    timer.start();
    an.Assemble();
    timer.stop();
    cout << timer.seconds() << " s" << endl;
    
    
    std::string exePath(argv[0]);
    stringstream ss;
    ss << numthreads;
    string strnthreads = ss.str();
    string filename = exePath + ".txt";// + "_nthreads_" + strnthreads + ".txt";
    ofstream out;
    std::cout << "Output filename " << filename.c_str() << std::endl;
    if (numthreads == 0) {
        out.open(filename.c_str());
    }else{
        out.open(filename.c_str(), ofstream::out | ofstream::app);
    }

    out << "\nRodando com:" << endl;
    out << "nref = " << nref << endl;
    out << "plevel = " << plevel << endl;
    out << "numthreads = " << numthreads << endl;
    out << "T assemble = " << timer.seconds() << endl;
    
    return 0;
    
    an.Solve();
    // Pos Processamento
    TPZStack<string> scalnames, vecnames;
    scalnames.Push("StressX");
    vecnames.Push("state");
    string plotfile = "cubovtk.vtk";
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    an.PostProcess(0);
    
    std::cout << "FINISHED" << std::endl;
	return 0;
}



// ------------------------ Para testes do assemble -----------------------------

TPZGeoMesh *MalhaCubo(string &projectpath, const int &nref)
{
    long numnodes=-1;
    long numelements=-1;
    
    string FileName, dirname = PZSOURCEDIR;
    FileName = dirname + projectpath;
    FileName += "cube1.txt";
    
    {
        bool countnodes = false;
        bool countelements = false;
        
        ifstream read (FileName.c_str());
        
        while(read)
        {
            char buf[1024];
            read.getline(buf, 1024);
            std::string str(buf);
            if(str == "Coordinates") countnodes = true;
            if(str == "end coordinates") countnodes = false;
            if(countnodes) numnodes++;
            
            if(str == "Elements") countelements = true;
            if(str == "end elements") countelements = false;
            if(countelements) numelements++;
        }
    }
    
    TPZGeoMesh * gMesh = new TPZGeoMesh;
    
    gMesh -> NodeVec().Resize(numnodes);
    
    TPZManVector <long> TopolTetra(4);
    
    const long Qnodes = numnodes;
    TPZVec <TPZGeoNode> Node(Qnodes);
    
    //setting nodes coords
    long nodeId = 0, elementId = 0, matElId = 1;
    
    ifstream read;
    read.open(FileName.c_str());
    
    double nodecoordX , nodecoordY , nodecoordZ ;
    
    char buf[1024];
    read.getline(buf, 1024);
    read.getline(buf, 1024);
    std::string str(buf);
    long in;
    for(in=0; in<numnodes; in++)
    {
        read >> nodeId;
        read >> nodecoordX;
        read >> nodecoordY;
        read >> nodecoordZ;
        Node[nodeId-1].SetNodeId(nodeId);
        Node[nodeId-1].SetCoord(0,nodecoordX);
        Node[nodeId-1].SetCoord(1,nodecoordY);
        Node[nodeId-1].SetCoord(2,nodecoordZ);
        gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
    }
    
    {
        read.close();
        read.open(FileName.c_str());
        
        long l , m = numnodes+5;
        for(l=0; l<m; l++)
        {
            read.getline(buf, 1024);
        }
        
        
        long el;
        int neumann1 = -4, neumann2 = -5;
        //std::set<int> ncoordz; //jeitoCaju
        for(el=0; el<numelements; el++)
        {
            read >> elementId;
            read >> TopolTetra[0]; //node 1
            read >> TopolTetra[1]; //node 2
            read >> TopolTetra[2]; //node 3
            read >> TopolTetra[3]; //node 4
            
            // O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
            TopolTetra[0]--;
            TopolTetra[1]--;
            TopolTetra[2]--;
            TopolTetra[3]--;
            
            long index = el;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
        }
        
        gMesh->BuildConnectivity();
        
        // Colocando as condicoes de contorno
        for(el=0; el<numelements; el++)
        {
            TPZManVector <TPZGeoNode,4> Nodefinder(4);
            TPZManVector <REAL,3> nodecoord(3);
            TPZGeoEl *tetra = gMesh->ElementVec()[el];
            
            // na face x = 1
            TPZVec<long> ncoordzVec(0); long sizeOfVec = 0;
            for (int i = 0; i < 4; i++)
            {
                long pos = tetra->NodeIndex(i);
                Nodefinder[i] = gMesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[0] == 1.)
                {
                    sizeOfVec++;
                    ncoordzVec.Resize(sizeOfVec);
                    ncoordzVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordzVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordzVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,neumann1);
            }
            
            // Na face x = -1
            ncoordzVec.Resize(0);
            sizeOfVec = 0;
            for (int i = 0; i < 4; i++) 
            {
                long pos = tetra->NodeIndex(i);
                Nodefinder[i] = gMesh->NodeVec()[pos];
                
                Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[0] == -1.)
                {
                    sizeOfVec++;
                    ncoordzVec.Resize(sizeOfVec);
                    ncoordzVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordzVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordzVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,neumann2);	
            }
            
        }
        
        TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
        yz[0] = 1.;
        z[2] = -1;
        int bcidxyz = -1, bcidyz = -2, bcidz = -3;
        SetPointBC(gMesh, xyz, bcidxyz);
        SetPointBC(gMesh, yz, bcidyz);
        SetPointBC(gMesh, z, bcidz);
    }
    
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        const int nel = gMesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gMesh->Element(iel);
            if (gel || !gel->HasSubElement() || gel->Dimension() > 0) {
                gel->Divide(sons);
            }
        }
    }
    
    return gMesh;
}

/// Generate a boundary geometric element at the indicated node
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
{
    // look for an element/corner node whose distance is close to start
    TPZGeoNode *gn1 = gr->FindNode(x);
    long iel;
    long nelem = gr->ElementVec().NElements();
    TPZGeoEl *gel;
    for (iel = 0; iel<nelem; iel++) {
        gel = gr->ElementVec()[iel];
        if(!gel) continue;
        int nc = gel->NCornerNodes();
        int c;
        for (c=0; c<nc; c++) {
            TPZGeoNode *gn = gel->NodePtr(c);
            if (gn == gn1) {
                break;
            }
        }
        if (c<nc) {
            TPZGeoElBC(gel, c, bc);
            return;
        }
    }
}

void InsertElasticityCubo(TPZCompMesh *mesh)
{
    mesh->SetDimModel(3);
    int nummat = 1, neumann = 1, mixed = 2;
    //	int dirichlet = 0;
    int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5;   //, dirp2 = -6;
    TPZManVector<STATE> force(3,0.);
    //force[1] = 0.;
    
    STATE ElaE = 1000., poissonE = 0.2;   //, poissonV = 0.1, ElaV = 100.;
    
    STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
    lambdaV = 11.3636;
    muV = 45.4545;
    alpha = 1.;
    deltaT = 0.01;
    
    //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
    //viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
    //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, ElaE, poissonE, lambdaV, muV, alphaT, force);
    TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, ElaE, poissonE, force);
    
    TPZFNMatrix<6> qsi(6,1,0.);
    //viscoelast->SetDefaultMem(qsi); //elast
    //int index = viscoelast->PushMemItem(); //elast
    TPZMaterial * viscoelastauto(viscoelast);
    mesh->InsertMaterialObject(viscoelastauto);
    
    // Neumann em x = 1;
    TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
    val2(0,0) = 1.;
    TPZBndCond *bc4 = viscoelast->CreateBC(viscoelastauto, neumann1, neumann, val1, val2);
    TPZMaterial * bcauto4(bc4);
    mesh->InsertMaterialObject(bcauto4);
    
    // Neumann em x = -1;
    val2(0,0) = -1.;
    TPZBndCond *bc5 = viscoelast->CreateBC(viscoelastauto, neumann2, neumann, val1, val2);
    TPZMaterial * bcauto5(bc5);
    mesh->InsertMaterialObject(bcauto5);
    
    val2.Zero();
    // Dirichlet em -1 -1 -1 xyz;
    val1(0,0) = 1e4;
    val1(1,1) = 1e4;
    val1(2,2) = 1e4;
    TPZBndCond *bc1 = viscoelast->CreateBC(viscoelastauto, dir1, mixed, val1, val2);
    TPZMaterial * bcauto1(bc1);
    mesh->InsertMaterialObject(bcauto1);
    
    // Dirichlet em 1 -1 -1 yz;
    val1(0,0) = 0.;
    val1(1,1) = 1e4;
    val1(2,2) = 1e4;
    TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
    TPZMaterial * bcauto2(bc2);
    mesh->InsertMaterialObject(bcauto2);
    
    // Dirichlet em 1 1 -1 z;
    val1(0,0) = 0.;
    val1(1,1) = 0.;
    val1(2,2) = 1e4;
    TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
    TPZMaterial * bcauto3(bc3);
    mesh->InsertMaterialObject(bcauto3);
    
}
