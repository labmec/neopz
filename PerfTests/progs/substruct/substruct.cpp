/**
 * @file
 * @brief Tests for sub structuration
 * @author Philippe Devloo
 * @since 2006
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>

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
#include "TPZVTKGeoMesh.h"
#include "pzgeotetrahedra.h"
#include "pzskylstrmatrix.h"

#include "tpzarc3d.h"
#include "tpzdohrmatrix.h"

#include "pzvtkmesh.h"
#include "pzfstrmatrix.h"

#include "pzlog.h"

#include "pzbfilestream.h" // TPZBFileStream, TPZFileStream
#include "pzmd5stream.h"

#include <fstream>
#include <string>

#ifdef PZ_LOG
static TPZLogger loggerconverge("pz.converge");
static TPZLogger logger("main");
#endif

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

#ifdef USING_MKL
#include "mkl_service.h"
#endif

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh);
void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh);
void InsertElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh);
TPZGeoMesh *MalhaPredio();
TPZGeoMesh *MalhaCubo();
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);
REAL Height(TPZGeoMesh *gmesh);
int SubStructure(TPZAutoPointer<TPZCompMesh> cmesh, REAL height);

using namespace std;

void help(const char* prg)
{
    cout << "Compute ...." << endl;
    cout << "The application is divided in four main steps: s1, s2, s3 and s4" << endl;
    cout << endl;
    cout << "  Step 1 (mp/mc  -> ckpt 1): " << endl;
    cout << "  Step 2 (ckpt 1 -> ckpt 2): dohrstruct->Create()" << endl;
    cout << "  Step 3 (ckpt 2 -> ckpt 3): dohrstruct->Assemble(), dohrstruct->Preconditioner() " << endl;
    cout << "  Step 4 (ckpt 3 -> end   ): cg.solve() ..." << endl << endl;
    cout << "Usage: " << prg << "starting_point stop_point file [-of output_file]" << endl << endl;
    cout << " starting_point = {-cf1|-cf2|-cf3|-mc|-mp} input_file" << endl;
    cout << " starting_point = {-st1|-st2|-st3}" << endl;
    
    clarg::arguments_descriptions(cout, "  ", "\n");
}

clarg::argString cf1("-cf1", "starts execution from checkpoint 1 (read checkpoint file)", "ckpt1.ckpt");
clarg::argString cf2("-cf2", "starts execution from checkpoint 2 (read checkpoint file)", "ckpt2.ckpt");
clarg::argString cf3("-cf3", "starts execution from checkpoint 3 (read checkpoint file)", "ckpt3.ckpt");
clarg::argBool   st1("-st1", "stop at checkpoint 1 (after dump)", false);
clarg::argBool   st2("-st2", "stop at checkpoint 2 (after dump)", false);
clarg::argBool   st3("-st3", "stop at checkpoint 3 (after dump)", false);

clarg::argString gen_sig_ckpt1("-gen_c1_md5", "generates MD5 signature for checkpoint 1 and dump into file.", "ckpt1.md5");
clarg::argString chk_sig_ckpt1("-chk_c1_md5", "compute MD5 signature for checkpoint 1 and check against MD5 at file.", "ckpt1.md5");
clarg::argString gen_sig_ckpt2("-gen_c2_md5", "generates MD5 signature for checkpoint 2 and dump into file.", "ckpt2.md5");
clarg::argString chk_sig_ckpt2("-chk_c2_md5", "compute MD5 signature for checkpoint 2 and check against MD5 at file.", "ckpt2.md5");
clarg::argString gen_sig_ckpt3("-gen_c3_md5", "generates MD5 signature for checkpoint 3 and dump into file.", "ckpt3.md5");
clarg::argString chk_sig_ckpt3("-chk_c3_md5", "compute MD5 signature for checkpoint 3 and check against MD5 at file.", "ckpt3.md5");
clarg::argString gen_sig_ckpt4("-gen_c4_md5", "generates MD5 signature for checkpoint 4 and dump into file.", "ckpt4.md5");
clarg::argString chk_sig_ckpt4("-chk_c4_md5", "compute MD5 signature for checkpoint 4 and check against MD5 at file.", "ckpt4.md5");

clarg::argString dc1("-dc1", "dump checkpoint 1 to file", "ckpt1.ckpt");
clarg::argString dc2("-dc2", "dump checkpoint 2 to file", "ckpt2.ckpt");
clarg::argString dc3("-dc3", "dump checkpoint 3 to file", "ckpt3.ckpt");
clarg::argString mc("-mc", "starts execution from beginning - read a \"malha_cubo\" input file",
                    "../cube1.txt");
clarg::argString mp("-mp", "starts execution from beginning - read a \"malha_predio\" input file",
                    "../8andares02.txt");

clarg::argInt plevel     ("-p", "plevel", 1);
clarg::argInt num_it ("-num_it", "number limit of iterations to the CG solver", 1000);
clarg::argInt nt_sm("-nt_sm", "Dohr (l1): number of threads to process the submeshes", 0);
clarg::argInt nt_d("-nt_d", "Dohr (l1): number of threads to decompose each submesh", 0);
clarg::argInt nt_a("-nt_a", "Pair (l2): number of threads to assemble each submesh (multiply by nt_sm)", 0);
clarg::argBool dohr_tbb("-dohr_tbb", "assemble TPZDohrStructMatrix (level 1) using TBB", false);
clarg::argBool pair_tbb("-pair_tbb", "assemble TPZPairStructMatrix (level 2) using TBB", false);
clarg::argInt nt_multiply("-nt_m", "number of threads to multiply", 0);
clarg::argInt nsub("-nsub", "number of substructs", 4);

clarg::argInt dim_arg("-dim", "dim???", 3);
clarg::argInt maxlevel("-maxlevel", "maxlevel???", 5);
clarg::argInt sublevel("-sublevel", "sublevel???", 3);

clarg::argInt verb_level("-v", "verbosity level", 0);

clarg::argBool bc("-bc", "binary checkpoints", false);

clarg::argBool h("-h", "help message", false);

clarg::argInt refin("-ref", "refine mesh", 1);

/* Run statistics. */
RunStatsTable total_rst   ("-tot_rdt", "Whole program (total) statistics raw data table");
RunStatsTable create_rst  ("-cre_rdt", "Create statistics raw data table (step 2)");
RunStatsTable assemble_rst("-ass_rdt", "Assemble statistics raw data table (step 3)");
RunStatsTable precond_rst ("-pre_rdt", "Precond statistics raw data table (step 3)");
RunStatsTable solve_rst   ("-sol_rdt", "Solver statistics raw data table (step 4)");

#ifdef USING_LIKWID
#define PERF_START(obj)				\
likwid_markerStartRegion(#obj);		\
obj.start()
#define PERF_STOP(obj)				\
obj.stop();					\
likwid_markerStopRegion(#obj)
#else
#define PERF_START(obj) obj.start()
#define PERF_STOP(obj) obj.stop()
#endif


class FileStreamWrapper
{
public:
    FileStreamWrapper() {}
    ~FileStreamWrapper() {}
    
    void OpenWrite(const std::string& fn)
    {
        if (bc.was_set())
            bfs.OpenWrite(fn);
        else
            fs.OpenWrite(fn);
    }
    
    void OpenRead(const std::string& fn)
    {
        if (bc.was_set())
            bfs.OpenRead(fn);
        else
            fs.OpenRead(fn);
    }
    
    operator TPZStream&()
    {
        if (bc.was_set())
            return bfs;
        else
            return fs;
    }
    
protected:
    TPZFileStream  fs;
    TPZBFileStream bfs;
};

#ifdef USING_LIKWID
#include<likwid.h>

struct likwid_manager_t {
    likwid_manager_t() {
        std::cout << "Calling likwid_markerInit()" << std::endl;
        likwid_markerInit();
    }
    ~likwid_manager_t() {
        likwid_markerClose();
        std::cout << "Calling likwid_markerClose()" << std::endl;
    }
};

#endif

int main(int argc, char *argv[])
{
    
#ifdef USING_LIKWID
    likwid_manager_t likwid_manager;
#endif
    
#ifdef USING_LAPACK
    setenv("VECLIB_MAXIMUM_THREADS", "1", true);
#endif
    
#ifdef USING_MKL
    mkl_set_num_threads(1);
#endif
    
    int main_ret = EXIT_SUCCESS;
    
    /* Parse the arguments */
    if (clarg::parse_arguments(argc, argv)) {
        cerr << "Error when parsing the arguments!" << endl;
        return 1;
    }
    
    if (h.get_value() == true) {
        help(argv[0]);
        return 1;
    }
    
#ifdef USING_TBB
    int number_tbb=nt_a.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif
    
    /* Verbose macro. */
    unsigned verbose = verb_level.get_value();
#   define VERBOSE(level,...) if (level <= verbose) cout << __VA_ARGS__
    
    if (verbose >= 1) {
        std::cout << "- Arguments -----------------------" << std::endl;
        clarg::values(std::cout, false);
        std::cout << "-----------------------------------" << std::endl;
    }
    
    if (!mp.was_set() && !mc.was_set() && !cf1.was_set() &&
        !cf2.was_set() && !cf3.was_set())
    {
        cerr << "A \"starting_point\" must be provided!" << endl;
        help(argv[0]);
        return 1;
    }
    
    PERF_START(total_rst);
    
    if (pair_tbb.was_set())
        TPZPairStructMatrix::gNumThreads = -1;
    else
        TPZPairStructMatrix::gNumThreads = nt_a.get_value();
    
    TPZGeoMesh  *gmesh = 0;
    TPZAutoPointer<TPZCompMesh> cmeshauto = 0;
    TPZDohrStructMatrix* dohrstruct = 0;
    TPZFMatrix<STATE> *rhs = NULL;
    TPZMatrix<STATE> *matptr = 0;
    int dim = dim_arg.get_value();
    TPZCompEl::SetgOrder(plevel.get_value());
    
    bool running = false;
    
    /* Start from malha_cubo or malha_predio? */
    if (mp.was_set() || mc.was_set())
    {
        if (mp.was_set()) // Predio Elastisco
        {
            if (running) {
                cerr << "ERROR: you must select only one of the start modes: "
                << "mp, mc, cf1, cf2 or cf3" << endl;
                exit(1);
            }
            else
                running = true;
            
            gmesh = MalhaPredio();
            cmeshauto = new TPZCompMesh(gmesh);
            cmeshauto->SetDimModel(dim);
            InsertElasticity(cmeshauto);
            cmeshauto->AutoBuild();
        }
        if (mc.was_set()) // Cubo Elastico
        {
            if (running) {
                cerr << "ERROR: you must select only one of the start modes: "
                << "mp, mc, cf1, cf2 or cf3" << endl;
                exit(1);
            }
            else  running = true;
            
            VERBOSE(1, "Reading MalhaCubo from file: " << mc.get_value() << endl);
            gmesh = MalhaCubo();
            cmeshauto = new TPZCompMesh(gmesh);
            cmeshauto->SetDimModel(dim);
            cmeshauto->SetDefaultOrder(plevel.get_value());
            //cmeshauto->SetAllCreateFunctionsContinuousWithMem();
            //cmeshauto->SetAllCreateFunctionsContinuous();
            InsertElasticityCubo(cmeshauto);
            cmeshauto->AutoBuild();
        }
        
        VERBOSE(1, "Number of equations " << cmeshauto->NEquations() << endl);

//#define ASSEMBLE_PERF
#ifdef ASSEMBLE_PERF
		TPZFStructMatrix fullstruct(cmeshauto);
		fullstruct.SetNumThreads(nt_a.get_value());

		int64_t sz = cmeshauto->NEquations();

		/*
		// ************** PARALELO **************
		TPZFMatrix<STATE> rhs_t(sz, 1, 0.);
		fullstruct.Assemble(rhs_t, 0);

		// ************** SERIAL *************
		fullstruct.SetNumThreads(0);
		TPZFMatrix<STATE> rhs_b(sz, 1, 0.);
		fullstruct.Assemble(rhs_b, 0);

		for (int i=0; i<sz; i++) {
			if(fabs(rhs_b(i,0)-rhs_t(i,0))>1.e-9) {
				printf("%d - %.5f %.5f\n",i, rhs_b(i,0),rhs_t(i,0));
				exit(101);
			}
		}
		return 0;
*/
		TPZFMatrix<STATE> rhs_t(sz, 1, 0.);
		PERF_START(assemble_rst);
		fullstruct.Assemble(rhs_t, 0);
		PERF_STOP(assemble_rst);

		return 0;
#endif
        
        dohrstruct = new TPZDohrStructMatrix(cmeshauto);
        dohrstruct->IdentifyExternalConnectIndexes();
        
        VERBOSE(1, "Substructuring the mesh" << endl);
        
        dohrstruct->SubStructure(nsub.get_value());
        
#ifdef PZ_LOG
        {
            std::stringstream str;
            cmeshauto->Print(str);
            LOGPZ_DEBUG(logger,str.str());
        }
#endif
        
        /* Dump checkpoint 1? */
        if (dc1.was_set() && running)
        {
            VERBOSE(1, "Dumping checkpoint 1 into: " << dc1.get_value() << endl);
            FileStreamWrapper CheckPoint1;
            CheckPoint1.OpenWrite(dc1.get_value());
            cmeshauto->Reference()->Write(CheckPoint1, 0);
            cmeshauto->Write(CheckPoint1, 0);
            dohrstruct->Write(CheckPoint1);
        }
        /* Gen/Check checkpoint 1 MD5 signature? */
        if ((gen_sig_ckpt1.was_set() || chk_sig_ckpt1.was_set()) && running)
        {
            TPZMD5Stream sig;
            cmeshauto->Reference()->Write(sig, 0);
            cmeshauto->Write(sig, 0);
            dohrstruct->Write(sig);
            if (chk_sig_ckpt1.was_set()) {
                int ret;
                if ((ret=sig.CheckMD5(chk_sig_ckpt1.get_value()))) {
                    cerr << "ERROR: MD5 Signature for checkpoint 1 does not match. (ret = " << ret << ")" << endl;
                    return 1;
                }
            }
            if (gen_sig_ckpt1.was_set()) {
                int ret;
                if ((ret = sig.WriteMD5(gen_sig_ckpt1.get_value()))) {
                    cerr << "ERROR when writing ckpt 1 MD5 Signature to file (ret = " << ret << "): "
                    << gen_sig_ckpt1.get_value() << endl;
                    return 1;
                }
            }
        }
    }
    
    if(st1.was_set()) running = false;
    
    // Start from Checkpoint 1
    if (cf1.was_set())
    {
        if (running) {
            cerr << "ERROR: you must select only one of the start modes: mp, mc, cf1, cf2 or cf3" << endl;
            exit(1);
        }
        else
            running = true;
        
        gmesh = new TPZGeoMesh;
        cmeshauto = new TPZCompMesh(gmesh);
        dohrstruct = new TPZDohrStructMatrix(cmeshauto);
        /* Read the checkpoint. */
        {
            FileStreamWrapper CheckPoint1;
            CheckPoint1.OpenRead(cf1.get_value());
            gmesh->Read(CheckPoint1, 0);
            cmeshauto->Read(CheckPoint1, gmesh);
            dohrstruct->Read(CheckPoint1);
        }
        
        dim = cmeshauto->Dimension();
        VERBOSE(1, "Reading dim from file. new dim = " << dim << ", old dim = " << dim_arg.get_value() << endl);
        
    }
    
    /* Work between checkpoint 1 and checkpoint 2 */
    if (running) {
        
        PERF_START(create_rst);
        matptr = dohrstruct->Create();
        PERF_STOP(create_rst);
        
        if (dc2.was_set())
        {
            VERBOSE(1, "Dumping checkpoint 2 into: " << dc2.get_value() << endl);
            FileStreamWrapper CheckPoint2;
            CheckPoint2.OpenWrite(dc2.get_value());
            SAVEABLE_STR_NOTE(CheckPoint2,"cmeshauto->Reference()->Write()");
            cmeshauto->Reference()->Write(CheckPoint2, 0);
            SAVEABLE_STR_NOTE(CheckPoint2,"cmeshauto->Write()");
            cmeshauto->Write(CheckPoint2, 0);
            SAVEABLE_STR_NOTE(CheckPoint2,"matptr->Write()");
            matptr->Write(CheckPoint2, 1);
            SAVEABLE_STR_NOTE(CheckPoint2,"dohrstruct->Write()");
            dohrstruct->Write(CheckPoint2);
        }
        
        /* Gen/Check checkpoint 2 MD5 signature? */
        if (gen_sig_ckpt2.was_set() || chk_sig_ckpt2.was_set())
        {
            TPZMD5Stream sig;
            cmeshauto->Reference()->Write(sig, 0);
            cmeshauto->Write(sig, 0);
            matptr->Write(sig, 1);
            dohrstruct->Write(sig);
            
            if (chk_sig_ckpt2.was_set()) {
                if (sig.CheckMD5(chk_sig_ckpt2.get_value())) {
                    cerr << "ERROR: MD5 Signature for checkpoint 2 does not match." << endl;
                    return 1;
                }
            }
            if (gen_sig_ckpt2.was_set()) {
                if (sig.WriteMD5(gen_sig_ckpt2.get_value())) {
                    cerr << "ERROR when writing ckpt 2 MD5 Signature to file: "
                    << gen_sig_ckpt2.get_value() << endl;
                    return 1;
                }
            }
        }
    }
    
    if(st2.was_set()) running = false;
    
    // Start from Checkpoint 2
    if (cf2.was_set())
    {
        if (running) {
            cerr << "ERROR: you must select only one of the start modes: mp, mc, cf1, cf2 or cf3" << endl;
            exit(1);
        }
        else
            running = true;
        
        FileStreamWrapper CheckPoint2;
        CheckPoint2.OpenRead(cf2.get_value());
        gmesh = new TPZGeoMesh;
        SAVEABLE_SKIP_NOTE(CheckPoint2);
        gmesh->Read(CheckPoint2,0);
        cmeshauto = new TPZCompMesh(gmesh);
        SAVEABLE_SKIP_NOTE(CheckPoint2);
        cmeshauto->Read(CheckPoint2, &gmesh);
        SAVEABLE_SKIP_NOTE(CheckPoint2);
        matptr = dynamic_cast<TPZMatrix<STATE> *>(TPZSavable::Restore(CheckPoint2, 0));
        dohrstruct = new TPZDohrStructMatrix(cmeshauto);
        SAVEABLE_SKIP_NOTE(CheckPoint2);
        dohrstruct->Read(CheckPoint2);
    }
    
    TPZAutoPointer<TPZMatrix<STATE> > precond = NULL;
    /* Work between checkpoint 2 and checkpoint 3 */
    if (running)
    {
        if (nt_multiply.was_set() && nt_multiply.get_value() != 1)
        {
            dohrstruct->SetNumThreads(nt_multiply.get_value());
        }
        
        PERF_START(assemble_rst);
        TPZAutoPointer<TPZGuiInterface> gui;
        rhs = new TPZFMatrix<STATE>(cmeshauto->NEquations(),1,0.);
        VERBOSE(1,"dohrstruct->Assemble()" << endl);
        if (dohr_tbb.was_set())
            dohrstruct->AssembleTBB(*matptr,*rhs, gui);
        else
            dohrstruct->Assemble(*matptr,*rhs, gui, nt_sm.get_value(), nt_d.get_value());
        PERF_STOP(assemble_rst);
        
        PERF_START(precond_rst);
        precond = dohrstruct->Preconditioner();
        PERF_STOP(precond_rst);
        
        if (dc3.was_set())
        {
            VERBOSE(1, "Dumping checkpoint 3 into: " << dc3.get_value() << endl);
            FileStreamWrapper CheckPoint3;
            CheckPoint3.OpenWrite(dc3.get_value());
            cmeshauto->Reference()->Write(CheckPoint3, 0);
            cmeshauto->Write(CheckPoint3, 0);
            matptr->Write(CheckPoint3, 1);
            precond->Write(CheckPoint3, 1);
            rhs->Write(CheckPoint3, 0);
        }
        
        /* Gen/Check checkpoint 3 MD5 signature? */
        if (gen_sig_ckpt3.was_set() || chk_sig_ckpt3.was_set())
        {
            TPZMD5Stream sig;
            cmeshauto->Reference()->Write(sig, 0);
            cmeshauto->Write(sig, 0);
            matptr->Write(sig, 1);
            precond->Write(sig, 1);
            rhs->Write(sig, 0);
            int ret;
            if (chk_sig_ckpt3.was_set()) {
                if ((ret=sig.CheckMD5(chk_sig_ckpt3.get_value()))) {
                    cerr << "ERROR(ret=" << ret << ") : MD5 Signature for checkpoint 3 does not match." << endl;
                    return 1;
                }
            }
            if (gen_sig_ckpt3.was_set()) {
                if ((ret=sig.WriteMD5(gen_sig_ckpt3.get_value()))) {
                    cerr << "ERROR (ret=" << ret << ") when writing ckpt 3 MD5 Signature to file: "
                    << gen_sig_ckpt3.get_value() << endl;
                    return 1;
                }
            }
        }
    }
    
    if(st3.was_set()) running = false;
    
    // Start from Checkpoint 3
    if (cf3.was_set())
    {
        if (running) {
            cerr << "ERROR: you must select only one of the start modes: mp, mc, cf1, cf2 or cf3" << endl;
            exit(1);
        }
        else
            running = true;
        gmesh = new TPZGeoMesh;
        cmeshauto = new TPZCompMesh(gmesh);
        dohrstruct = new TPZDohrStructMatrix(cmeshauto);
        
        dim = cmeshauto->Dimension();
        VERBOSE(1, "Reading dim from file. new dim = " << dim << ", old dim = " << dim_arg.get_value() << endl);
        
        FileStreamWrapper CheckPoint3;
        CheckPoint3.OpenRead(cf3.get_value());
        gmesh->Read(CheckPoint3, 0);
        cmeshauto->Read(CheckPoint3, gmesh);
        matptr = dynamic_cast<TPZMatrix<STATE> *>(TPZSavable::Restore(CheckPoint3, 0));
        precond = dynamic_cast<TPZMatrix<STATE> *>(TPZSavable::Restore(CheckPoint3, matptr));
        rhs = new TPZFMatrix<STATE>(cmeshauto->NEquations(),1,0.);
        rhs->Read(CheckPoint3, 0);
    }
    
    int neq, iterations;
    neq = iterations = 0;
    
    
    
    if (running) {
        
        /* Work after checkpoint 3 */
        TPZAutoPointer<TPZMatrix<STATE> > dohr = matptr;
        
        neq = dohr->Rows();
        TPZFMatrix<STATE> diag(neq,1,0.), produto(neq,1);
        
        VERBOSE(1, "Number of equations " << neq << endl);
        
        TPZStepSolver<STATE> pre(precond);
        pre.SetMultiply();
        TPZStepSolver<STATE> cg(dohr);
        
        /* Configure the CG solver to iterate:
         - until it converges (residual <= 1.e-8), or
         - until it reaches 500 itearations.
         */
        cg.SetCG(num_it.get_value(),pre,1.e-8,0);
        
        PERF_START(solve_rst);
        cg.Solve(*rhs,diag);
        PERF_STOP(solve_rst);
        
        
        iterations = cg.NumIterations();
        
        /* checking if the solver converged */
        if (cg.GetTolerance() > 1.e-8)
        {
            cerr << "ERROR: solver do not converged with the limit of iterations."  << endl;
            exit(1);
        }
        
        TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohrptr =
        dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohr.operator->());
        
        if (!dohrptr) {
            DebugStop();
        }
        
        dohrptr->AddInternalSolution(diag);
        
        
        /* Gen/Check checkpoint 4 MD5 signature? */
        if (gen_sig_ckpt4.was_set() || chk_sig_ckpt4.was_set())
        {
            TPZMD5Stream sig;
            
            dohrptr->Write(sig, 0);
            diag.Write(sig, 0);
            cg.Write(sig, 0);
            //cmeshauto->Reference()->Write(sig, 0);
            //cmeshauto->Write(sig, 0);
            //matptr->Write(sig, 1);
            //precond->Write(sig, 1);
            //rhs->Write(sig, 0);
            
            int ret;
            if (chk_sig_ckpt4.was_set()) {
                if ((ret=sig.CheckMD5(chk_sig_ckpt4.get_value()))) {
                    cerr << "ERROR(ret=" << ret << ") : MD5 Signature for checkpoint 4 does not match." << endl;
                    main_ret = 1;
                }
            }
            if (gen_sig_ckpt4.was_set()) {
                if ((ret=sig.WriteMD5(gen_sig_ckpt4.get_value()))) {
                    cerr << "ERROR (ret=" << ret << ") when writting ckpt 4 MD5 Signature to file: "
                    << gen_sig_ckpt4.get_value() << endl;
                    main_ret = 2;
                }
            }
        }
        
        typedef std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > > subtype;
        const subtype &sublist = dohrptr->SubStructures();
        subtype::const_iterator it = sublist.begin();
        int subcount=0;
        while (it != sublist.end()) {
            
            TPZFMatrix<STATE> subext,subu;
            dohrptr->fAssembly->Extract(subcount,diag,subext);
            (*it)->UGlobal(subext,subu);
            TPZCompMesh *submesh = SubMesh(cmeshauto, subcount);
            submesh->LoadSolution(subu);
            
            //ViscoElastico
            //Atualizando memoria do material
            std::map<int ,TPZMaterial * > materialmap(submesh->MaterialVec());
            std::map<int ,TPZMaterial * >::iterator itmat;
            for (itmat = materialmap.begin(); itmat != materialmap.end() ; itmat++)
            {
                TPZMaterial * mat = itmat->second;
                TPZViscoelastic *vmat = dynamic_cast< TPZViscoelastic *> (mat);
                if(vmat)
                {
                    DebugStop(); // Should never enter because it is using elasticity
                    vmat->SetUpdateMem(true);
                }
            }
            subcount++;
            it++;
        }
        
#ifdef PZ_LOG
        {
            std::stringstream sout;
            diag.Print("Resultado do processo iterativo",sout);
            LOGPZ_INFO(loggerconverge,sout.str())
        }
#endif
        
        TPZMaterial * mat = cmeshauto->FindMaterial(1);
        int nstate = mat->NStateVariables();
        int nscal = 0, nvec = 0;
        if(nstate ==1)
        {
            nscal = 1;
        }
        else
        {
            nvec = 1;
        }
        TPZManVector<std::string> scalnames(nscal),vecnames(nvec);
        if(nscal == 1)
        {
            scalnames[0]="state";
        }
        else
        {
            vecnames[0] = "state";
        }
        std::string postprocessname("dohrmann_visco.vtk"); // Remember it is not viscoelastic, just elastic!
        TPZVTKGraphMesh vtkmesh(cmeshauto.operator->(),dim,mat,scalnames,vecnames);
        vtkmesh.SetFileName(postprocessname);
        vtkmesh.SetResolution(1);
        int numcases = 1;
        
        
        // Iteracoes de tempo
        int istep = 0;
        vtkmesh.DrawMesh(numcases);
        vtkmesh.DrawSolution(istep, 1.);
    }
    
    PERF_STOP(total_rst);
    
    cout << " -- Execution Data -- " << endl;
    cout << "Input                              :   ";
    if (mc.was_set())
    {
        cout << mc.get_value() << endl;
    } else {
        cout << mp.get_value() << endl;
    }
    cout << "Nsub                               :   " << nsub.get_value() << endl;
    cout << "Equations (cmeshauto->NEquations)  :   " << cmeshauto->NEquations() << endl;
    cout << "Corner Equations                   :   " << dohrstruct->NumberCornerEqs() << endl;
    cout << "Equations (dohr->Rows)             :   " << neq << endl;
    cout << "Iterations (CG)                    :   " << iterations << endl;
    
    if (gmesh != NULL) delete gmesh;
    
    return main_ret;
}

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
    mesh->SetDimModel(3);
    int nummat = 1;
    STATE E = 1.e6;
    STATE poisson = 0.3;
    TPZManVector<STATE> force(3,0.);
    force[1] = 20.;
    TPZElasticity3D *elast = new TPZElasticity3D(nummat,E,poisson,force);
    TPZMaterial * elastauto(elast);
    TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
    TPZBndCond *bc = elast->CreateBC(elastauto, -1, 0, val1, val2);
    TPZMaterial * bcauto(bc);
    mesh->InsertMaterialObject(elastauto);
    mesh->InsertMaterialObject(bcauto);
}

void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
    mesh->SetDimModel(3);
    int nummat = 1;
    TPZManVector<STATE> force(3,0.);
    force[1] = 20.;
    STATE ElaE = 1000., poissonE = 0.2, ElaV = 100., poissonV = 0.1;
    
    STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
    lambdaV = 11.3636;
    muV = 45.4545;
    alpha = 1.;
    deltaT = 0.01;
    
    TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
    viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
    
    TPZMaterial * viscoelastauto(viscoelast);
    TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
    TPZBndCond *bc = viscoelast->CreateBC(viscoelastauto, -1, 0, val1, val2);
    TPZFNMatrix<6,STATE> qsi(6,1,0.);
    viscoelast->SetDefaultMem(qsi); //elast
    viscoelast->PushMemItem(); //elast
    TPZMaterial * bcauto(bc);
    mesh->InsertMaterialObject(viscoelastauto);
    mesh->InsertMaterialObject(bcauto);
}

void InsertElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh)
{
    mesh->SetDimModel(3);
    int nummat = 1, neumann = 1, mixed = 2;
    //      int dirichlet = 0;
    int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5;
    TPZManVector<STATE> force(3,0.);
    //force[1] = 0.;
    REAL Ela = 1000., poisson = 0.;
    REAL lambdaV = 0, muV = 0, alphaT = 0;
    lambdaV = 11.3636;
    muV = 45.4545;
    alphaT = 0.01;
    
    
    //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, Ela, poisson, lambdaV, muV, alphaT, force);
    TPZElasticity3D *elast = new TPZElasticity3D(nummat, Ela, poisson, force);
    
    TPZFNMatrix<6> qsi(6,1,0.);
    //viscoelast->SetDefaultMem(qsi); //elast
    //int index = viscoelast->PushMemItem(); //elast
    TPZMaterial * elastauto(elast);
    mesh->InsertMaterialObject(elastauto);
    
    // Neumann em x = 1;
    TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
    val2(0,0) = 1.;
    TPZBndCond *bc4 = elast->CreateBC(elastauto, neumann1, neumann, val1, val2);
    TPZMaterial * bcauto4(bc4);
    mesh->InsertMaterialObject(bcauto4);
    
    // Neumann em x = -1;
    val2(0,0) = -1.;
    TPZBndCond *bc5 = elast->CreateBC(elastauto, neumann2, neumann, val1, val2);
    TPZMaterial * bcauto5(bc5);
    mesh->InsertMaterialObject(bcauto5);
    
    val2.Zero();
    // Dirichlet em -1 -1 -1 xyz;
    val1(0,0) = 1.e4;
    val1(1,1) = 1.e4;
    val1(2,2) = 1.e4;
    TPZBndCond *bc1 = elast->CreateBC(elastauto, dir1, mixed, val1, val2);
    TPZMaterial * bcauto1(bc1);
    mesh->InsertMaterialObject(bcauto1);
    
    // Dirichlet em 1 -1 -1 yz;
    val1(0,0) = 0.;
    val1(1,1) = 1.e4;
    val1(2,2) = 1.e4;
    TPZBndCond *bc2 = elast->CreateBC(elastauto, dir2, mixed, val1, val2);
    TPZMaterial * bcauto2(bc2);
    mesh->InsertMaterialObject(bcauto2);
    
    // Dirichlet em 1 1 -1 z;
    val1(0,0) = 0.;
    val1(1,1) = 0.;
    val1(2,2) = 1.e4;
    TPZBndCond *bc3 = elast->CreateBC(elastauto, dir3, mixed, val1, val2);
    TPZMaterial * bcauto3(bc3);
    mesh->InsertMaterialObject(bcauto3);
}

TPZGeoMesh *MalhaPredio()
{
    //int nBCs = 1;
    int numnodes=-1;
    int numelements=-1;
    
    string FileName = mp.get_value();
    {
        bool countnodes = false;
        bool countelements = false;
        
        ifstream read (FileName.c_str());
        if (!read.is_open()) {
            cerr << "Could not open file: " << FileName << endl;
            exit(1);
        }
        
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
    
    TPZVec <int64_t> TopolTetra(4);
    
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
        int matBCid = -1;
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
            
            int index = el;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
        }
        
        gMesh->BuildConnectivity();
        // Colocando as condicoes de contorno
        
        for(el=0; el<numelements; el++)
        {
            TPZManVector <TPZGeoNode,4> Nodefinder(4);
            TPZManVector <REAL,3> nodecoord(3);
            TPZGeoEl *tetra = gMesh->ElementVec()[el];
            // na face z = 0
            TPZVec<int64_t> ncoordzVec(0); int64_t sizeOfVec = 0;
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gMesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[2] == 0.)
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
                TPZGeoElBC(tetraSide,matBCid);
            }
        }
    }
    
    
    ofstream arg("malhaPZ.txt");
    gMesh->Print(arg);
    ofstream predio("GeoPredio.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true);
    
    return gMesh;
}

TPZGeoMesh *MalhaCubo()
{
    int numnodes=-1;
    int numelements=-1;
    
    string FileName = mc.get_value();
    {
        bool countnodes = false;
        bool countelements = false;
        
        ifstream read (FileName.c_str());
        if (!read.is_open()) {
            cerr << "Could not open file: " << FileName << endl;
            exit(1);
        }
        
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
                int64_t pos = tetra->NodeIndex(i);
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
                int64_t pos = tetra->NodeIndex(i);
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
    
    /* refine mesh */
    if (refin.was_set()) {
        
        int nh = refin.get_value();
        
        for ( int ref = 0; ref < nh; ref++ ){
            TPZVec<TPZGeoEl *> filhos;
            int n = gMesh->NElements();
            for ( int i = 0; i < n; i++ ){
                TPZGeoEl * gel = gMesh->ElementVec() [i];
                
                if(!gel) continue;
                if(gel->Dimension() < 1) continue;
                if(gel->HasSubElement()) continue;
                
                gel->Divide (filhos);
            }
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

REAL Height(TPZGeoMesh *gmesh)
{
    TPZAdmChunkVector<TPZGeoNode> &nodevec = gmesh->NodeVec();
    int64_t nnodes = nodevec.NElements();
    int64_t in;
    REAL maxz = 0.;
    for (in=0; in<nnodes; in++) {
        REAL z = nodevec[in].Coord(2);
        maxz = (maxz < z) ? z : maxz;
    }
    return maxz;
}

int SubStructure(TPZAutoPointer<TPZCompMesh> cmesh, REAL height)
{
    int nelem = cmesh->NElements();
    TPZManVector<int> subindex(nelem,-1);
    int iel;
    int nsub = 0;
    for (iel=0; iel<nelem; iel++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        int nsides = gel->NSides();
        TPZManVector<REAL> center(gel->Dimension(),0.), xco(3,0.);
        gel->CenterPoint(nsides-1,center);
        gel->X(center,xco);
        REAL z = xco[2];
        int floor = (int) z/height;
        nsub = (floor+1) > nsub ? (floor+1) : nsub;
        subindex[iel] = floor;
    }
    
#ifdef PZDEBUG
    {
        TPZGeoMesh *gmesh = cmesh->Reference();
        int nelgeo = gmesh->NElements();
        TPZVec<int> domaincolor(nelgeo,-999);
        int cel;
        int nel = cmesh->NElements();
        for (cel=0; cel<nel; cel++) {
            TPZCompEl *compel = cmesh->ElementVec()[cel];
            if(!compel) continue;
            TPZGeoEl *gel = compel->Reference();
            if (!gel) {
                continue;
            }
            domaincolor[gel->Index()] = subindex[cel];
        }
        ofstream vtkfile("partition.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
    }
#endif
    
    int isub;
    TPZManVector<TPZSubCompMesh *> submeshes(nsub,0);
    for (isub=0; isub<nsub; isub++) 
    {
        int64_t index;
        std::cout << '^'; std::cout.flush();
        submeshes[isub] = new TPZSubCompMesh(cmesh,index);
        
        if (index < subindex.NElements()) 
        {
            subindex[index] = -1;
        }
    }
    for (iel=0; iel<nelem; iel++) 
    {
        int domindex = subindex[iel];
        if (domindex >= 0) 
        {
            TPZCompEl *cel = cmesh->ElementVec()[iel];
            if (!cel) 
            {
                continue;
            }
            submeshes[domindex]->TransferElement(cmesh.operator->(),iel);
        }
    }
    cmesh->ComputeNodElCon();
    for (isub=0; isub<nsub; isub++) 
    {
        submeshes[isub]->MakeAllInternal();
        std::cout << '*'; std::cout.flush();
    }
    
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    return nsub;
}
