/**
 * @file
 * @brief Tests for Decompose_Cholesky, Decompose_LDLt and MultAdd with Skyline Matrices
 * @author Julia Ramos Beltrao
 * @since 2013
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "TPZFileStream.h"
#include "TPZBFileStream.h" // TPZBFileStream, TPZFileStream
#include "pzmd5stream.h"

#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrprecond.h"
#include "pzdohrstructmatrix.h"
#include "pzstepsolver.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"

#include "pzelast3d.h"
#include "pzbndcond.h"

#include "tpzdohrassembly.h"

#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "tpzpairstructmatrix.h"
#include "pzviscoelastic.h"
#include "TPZTimer.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeotetrahedra.h"
#include "pzskylstrmatrix.h"

#include "pzskylmat.h"

#include "tpzarc3d.h"
#include "tpzdohrmatrix.h"

#include "pzvtkmesh.h"

#include "pzlog.h"

#include <fstream>
#include <string>

#ifdef PZ_LOG
static TPZLogger loggerconverge("pz.converge");
static TPZLogger logger("main");
#endif

#include "pzskylmat.h"

//#include "timing_analysis.h"
#include "arglib.h"
#include "run_stats_table.h"

#ifdef HAS_GETRUSAGE
#include <sys/resource.h> // getrusage
#endif

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
using namespace tbb;
// If you have issues with: dyld: Library not loaded: libtbb.dylib
// try setting the LD path. Ex:
//   export DYLD_FALLBACK_LIBRARY_PATH=/Users/borin/Desktop/neopz/tbb40_297oss/lib/
#endif

using namespace std;

void InsertElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh);
TPZGeoMesh *MalhaCubo(string FileName);
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);

void CopyTo(TPZSkylMatrix <REAL> * skylmat, TPZSkylMatrix<TPZFlopCounter> &fp_counter);
void PrintInfo (string funct, TPZCounter info);

void help(const char* prg)
{
	cout << "Compute the Decompose_Cholesky, Decompose_LDLt, MultAdd methods for a matrix generated of a mesh" << endl;
	cout << endl;
    cout << "Usage: " << prg << "-if file [-v verbose_level] "
    << "[-clk_rdt rdt_file1] [-ldlt_rdt rdt_file2] [-mult_rdt rdt_file3] [-h]" << endl << endl;
    clarg::arguments_descriptions(cout, "  ", "\n");
}


clarg::argBool h("-h", "help message", false);
clarg::argBool print_flops("-pf", "print floating point operations", false);

clarg::argInt verb_level("-v", "verbosity level", 0);
int verbose = 0;
/* Verbose macro. */
#define VERBOSE(level,...) if (level <= verbose) cout << __VA_ARGS__

clarg::argInt porder("-porder", "polinomial order", 1);
clarg::argString input("-if", "input file", "cube1.txt");

/** Uncoment for parsing matrix from input **/
//clarg::argString ifn("-ifn", "input matrix file name (use -bi to read from binary files)", "matrix.txt");
//clarg::argBool br("-br", "binary reference. Reference decomposed matrix file format == binary.", false);
//clarg::argBool bi("-bi", "binary input. Input file format == binary.", false);
//clarg::argBool bd("-bd", "binary dump. Dump file format == binary.", false);

/** Uncoment for md5 checking of resulting matrix **/
//clarg::argString gen_dm_sig("-gen_dm_md5", "generates MD5 signature for decomposed matrix into file.", "decomposed_matrix.md5");
//clarg::argString chk_dm_sig("-chk_dm_md5", "compute MD5 signature for decomposed matrix and check against MD5 at file.", "decomposed_matrix.md5");
//clarg::argString chk_dm_error("-chk_dm_error", "check the decomposed matrix error against a reference matrix. (use -br to read from binary files)", "ref_decomposed_matrix.txt");
//clarg::argDouble error_tol("-error_tol", "error tolerance.", 1.e-12);

/** Uncoment for dumping resulting matrix **/
//clarg::argString dump_dm("-dump_dm", "dump decomposed matrix. (use -bd for binary format)", "dump_matrix.txt");


/* Run statistics. */
RunStatsTable clk_rst("-clk_rdt", "Decompose_Cholesky() statistics raw data table");
RunStatsTable ldlt_rst("-ldlt_rdt", "Decompose_LDLt() statistics raw data table");
RunStatsTable sfwd_rst("-sfwd_rdt", "Subst_Forward(...) statistics raw data table");
RunStatsTable sbck_rst("-sbck_rdt", "Subst_Backward(...) statistics raw data table");
RunStatsTable mult_rst("-mult_rdt", "MultAdd(...) statistics raw data table");
RunStatsTable sor_rst("-sor_rdt", "SolveSOR(...) statistics raw data table");


class FileStreamWrapper
{
public:
    FileStreamWrapper(bool b) : binary(b) {}
    ~FileStreamWrapper() {}
    
    void OpenWrite(const std::string& fn)
    {
        if (binary)
            bfs.OpenWrite(fn);
        else
            fs.OpenWrite(fn);
    }
    
    void OpenRead(const std::string& fn)
    {
        if (binary)
            bfs.OpenRead(fn);
        else
            fs.OpenRead(fn);
    }
    
    operator TPZStream&()
    {
        if (binary)
            return bfs;
        else
            return fs;
    }
    
protected:
    
    bool binary;
    TPZFileStream  fs;
    TPZBFileStream bfs;
};


#include <sched.h>     //sched_getcpu

//template class TPZSkylMatrix<TPZFlopCounter>;
//template class TPZMatrix<TPZFlopCounter>;
//template class TPZAutoPointer<TPZMatrix<TPZFlopCounter> >;

#ifdef USING_PAPI
#include <papi.h>
#endif

int main(int argc, char *argv[])
{
	int dim = 3;
    
	/* Parse the arguments */
	if (clarg::parse_arguments(argc, argv)) {
		cerr << "Error when parsing the arguments!" << endl;
        return 1;
    }
    
    verbose = verb_level.get_value();
    
    if (h.get_value() == true) {
        help(argv[0]);
        return 1;
    }
    
    if (verbose >= 1) {
        std::cout << "- Arguments -----------------------" << std::endl;
        clarg::values(std::cout, false);
        std::cout << "-----------------------------------" << std::endl;
    }
    
    /* Generating matrix. */
    
    
	TPZGeoMesh *gmesh = 0;
    
    TPZAutoPointer<TPZCompMesh> cmesh;
    
    gmesh = MalhaCubo(input.get_value());
    cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    InsertElasticityCubo(cmesh);
    cmesh->SetDefaultOrder(porder.get_value());
    cmesh->AutoBuild();
	int64_t neq = cmesh->NEquations();
	cout << "numero de equacoes = " << neq << endl;
    
    // Gerando a matriz
    TPZLinearAnalysis an(cmesh);
    TPZSkylineStructMatrix skyl(cmesh);
    an.SetStructuralMatrix(skyl);
    TPZStepSolver<REAL> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Assemble();
    TPZAutoPointer<TPZMatrix<REAL> > skylmat1 = new TPZSkylMatrix<REAL>();
    skylmat1 = an.Solver().Matrix();
	TPZSkylMatrix<REAL> *orig = dynamic_cast<TPZSkylMatrix<REAL> *> (skylmat1.operator->());
	TPZAutoPointer<TPZMatrix<REAL> > skylmat2 = skylmat1->Clone();
	TPZAutoPointer<TPZMatrix<REAL> > skylmat3 = skylmat1->Clone();
	TPZAutoPointer<TPZMatrix<REAL> > skylmat4 = skylmat1->Clone();
	TPZFMatrix<REAL> f(neq,1,M_PI);
	TPZFMatrix<TPZFlopCounter > f2(neq,1,M_PI);
	TPZSkylMatrix <TPZFlopCounter> fp1;
	TPZSkylMatrix <TPZFlopCounter> fp2;
	TPZSkylMatrix <TPZFlopCounter> fp3;
	TPZSkylMatrix <TPZFlopCounter> fp4;
	
	CopyTo(orig,fp1);
	CopyTo(orig,fp2);
	CopyTo(orig,fp3);
	CopyTo(orig,fp4);
    
	if (clk_rst.was_set()) {
        clk_rst.start();
        skylmat1->Decompose_Cholesky();
        clk_rst.stop();
		if (print_flops.was_set()) {
            TPZCounter c = TPZFlopCounter::gCount;
            fp1.Decompose_Cholesky();
            c = TPZFlopCounter::gCount - c;
            PrintInfo("Decompose_Cholesky()", c);
		}
	}
    
	if (ldlt_rst.was_set()) {
        ldlt_rst.start();
        skylmat2->Decompose_LDLt();
        ldlt_rst.stop();
		if (print_flops.was_set()) {
            TPZCounter c = TPZFlopCounter::gCount;
            fp2.Decompose_LDLt();
            c = TPZFlopCounter::gCount - c;
            PrintInfo("Decompose_LDLt()", c);
		}
	}
    
	if (sfwd_rst.was_set()) {
		if (!clk_rst.was_set()) skylmat1->Decompose_Cholesky();
        sfwd_rst.start();
#ifdef USING_PAPI
        float rtime1, ptime1, mflops1;
        int64_t flpops1;
        PAPI_flops ( &rtime1, &ptime1, &flpops1, &mflops1 );
#endif
        skylmat1->Subst_Forward(&f);
#ifdef USING_PAPI
        float rtime2, ptime2, mflops2;
        int64_t flpops2;
        PAPI_flops ( &rtime2, &ptime2, &flpops2, &mflops2 );
#endif
        sfwd_rst.stop();
#ifdef USING_PAPI
		cout << "skylmat1->Subst_Forward(&f)" << endl;
		cout << " rtime:  " << rtime2-rtime1 << endl;
		cout << " ptime:  " << ptime2-ptime1 << endl;
		cout << " flpops: " << flpops2-flpops1 << endl;
		cout << " mflops: " << mflops2-mflops1 << endl;
#endif
		if (print_flops.was_set()) {
            TPZCounter c = TPZFlopCounter::gCount;
            fp1.Subst_Forward(&f2);
            c = TPZFlopCounter::gCount - c;
            PrintInfo("Subst_Forward(...)", c);
		}
	}
    
	if (sbck_rst.was_set()){
		if (!clk_rst.was_set()) skylmat1->Decompose_Cholesky();
        sbck_rst.start();
        skylmat1->Subst_Backward(&f);
        sbck_rst.stop();
		if (print_flops.was_set()) {
            TPZCounter c = TPZFlopCounter::gCount;
            fp1.Subst_Backward(&f2);
            c = TPZFlopCounter::gCount - c;
            PrintInfo("Subst_Backward(...)", c);
		}
	}
    
	if (sor_rst.was_set()) {
        TPZFMatrix<REAL> res(neq,1);
        TPZFMatrix<REAL> residual(neq,1);
        TPZFMatrix<REAL> scratch(neq,neq);
		int64_t niter = 10;
		REAL overrelax = 1.1;
		REAL tol = 10e-7;
        sor_rst.start();
		//skylmat4->SolveSOR(niter, *f, res, &residual, scratch, overrelax, tol);
        sor_rst.stop();
		if (print_flops.was_set()) {
            TPZFMatrix<TPZFlopCounter> f3(neq, 1);
            TPZFMatrix<TPZFlopCounter> res2(neq,1);
            TPZFMatrix<TPZFlopCounter> residual2(neq,1);
            TPZFMatrix<TPZFlopCounter> scratch2(neq,neq);
            TPZCounter c = TPZFlopCounter::gCount;
            fp3.SolveSOR(niter, f3, res2, &residual2, scratch2, overrelax, tol);
            c = TPZFlopCounter::gCount - c;
            PrintInfo("SolveSOR(...)", c);
		}
	}
    
	if (mult_rst.was_set()) {
        TPZFMatrix<REAL> result(neq,neq);
        mult_rst.start();
        skylmat3->MultAdd(result, result, result, 1., 0);
        mult_rst.stop();
		if (print_flops.was_set()) {
            TPZFMatrix<TPZFlopCounter> result2(neq,neq);
            TPZCounter c = TPZFlopCounter::gCount;
            fp4.MultAdd(result2, result2, result2, 1., 0);
            c = TPZFlopCounter::gCount - c;
            PrintInfo("MultAdd(...)", c);
		}
	}
    
	return 0;
    
}

void InsertElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
    int nummat = 1, neumann = 1, mixed = 2;
    //      int dirichlet = 0;
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

TPZGeoMesh *MalhaCubo(string FileName)
{
    int numnodes=-1;
    int numelements=-1;
    
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
    
    TPZManVector <int64_t> TopolTetra(4);
    
    const int Qnodes = numnodes;
    TPZVec <TPZGeoNode> Node(Qnodes);
    
    //setting nodes coords
    int64_t nodeId = 0, elementId = 0, matElId = 1;
    
    ifstream read;
    read.open(FileName.c_str());
    
	double nodecoordX , nodecoordY , nodecoordZ ;
    
    char buf[1024];
    read.getline(buf, 1024);
    read.getline(buf, 1024);
    std::string str(buf);
    int in;
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
        
        int l , m = numnodes+5;
        for(l=0; l<m; l++)
        {
            read.getline(buf, 1024);
        }
        
        
        int el;
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
            
            int64_t index = el;
            
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
            TPZVec<int64_t> ncoordzVec(0); int64_t sizeOfVec = 0;
            for (int i = 0; i < 4; i++)
            {
                int pos = tetra->NodeIndex(i);
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
                int pos = tetra->NodeIndex(i);
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
    
	int iref, nref=1;
	TPZVec <TPZGeoEl*> sons;
	for (iref=0; iref<nref; iref++)
	{
		int nelements = gMesh->NElements();
		for (int iel=0; iel<nelements; iel++)
        {
			TPZGeoEl * gel = gMesh->ElementVec()[iel];
			if (!gel) DebugStop();
			if(!gel->HasSubElement()) gel->Divide(sons);
		}
	}
    
    ofstream arg("malhaPZ1BC.txt");
    gMesh->Print(arg);
    
    std::ofstream out("Cube.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gMesh, out, true);
    
    return gMesh;
}

/// Generate a boundary geometric element at the indicated node
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
{
    // look for an element/corner node whose distance is close to start
    TPZGeoNode *gn1 = gr->FindNode(x);
    int64_t iel;
    int64_t nelem = gr->ElementVec().NElements();
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

void CopyTo(TPZSkylMatrix <REAL> * skylmat1, TPZSkylMatrix<TPZFlopCounter> &fp_counter)
{
	fp_counter.Redim(skylmat1->Rows(),skylmat1->Cols());
	TPZVec <int64_t> skylvec(skylmat1->Rows());
	for (int64_t i = 0; i < skylmat1->Rows(); i++) {
		skylvec[i] = i - skylmat1->Size(i) + 1 ;
	}
    
	fp_counter.SetSkyline(skylvec);
	for (int i = 0; i < skylmat1->Rows(); i++) {
		for (int j = 0 ; j < skylmat1->Cols(); j++) {
			fp_counter.PutVal(i, j, skylmat1->GetVal(i,j));
		}
	}
}

void PrintInfo (string funct, TPZCounter info)
{
	cout << "Printing flop count info of " << funct << endl;
	info.Print();
	cout << "****************************" << endl;
	return;
}
