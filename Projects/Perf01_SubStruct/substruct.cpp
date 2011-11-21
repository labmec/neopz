/***************************************************************************
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrprecond.h"
#include "pzdohrstructmatrix.h"
#include "pzstepsolver.h"
#include "pzcompel.h"

#include "pzelast3d.h"
#include "pzbndcond.h"

#include "tpzdohrassembly.h"

#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "tpzpairstructmatrix.h"

#include "pzvtkmesh.h"

#include "pzlog.h"

#include <fstream>
#include <string>

#include "tpzdohrmatrix.h"

#include "building_mesh.h"

#ifdef LOG4CXX
// #define LOG4CXX_TRACE(logger, expression) ...   
// #define LOG4CXX_DEBUG(logger, expression) ...   
// #define LOG4CXX_INFO(logger, expression) ...   
// #define LOG4CXX_WARN(logger, expression) ...   
// #define LOG4CXX_ERROR(logger, expression) ...   
// #define LOG4CXX_FATAL(logger, expression) ...   
static LoggerPtr perflog(Logger::getLogger("perf"));
#else
#define LOG4CXX_TRACE(logger, expression)
#define LOG4CXX_DEBUG(logger, expression)
#define LOG4CXX_INFO(logger, expression)
#define LOG4CXX_WARN(logger, expression)
#define LOG4CXX_ERROR(logger, expression)
#define LOG4CXX_FATAL(logger, expression)
#endif

#include "clock_timer.h"
#include "timing_analysis.h"
#include <sys/resource.h> // getrusage

using namespace std;

//EB: notes
// Subst_Backward e o dobro do Subst_Forward, mas possuem o mesmo numero de operacoes. pq?
//LOG4CXX

//EBORIN: dim2_2threas
//main: 28.1% (Create Assemble -- 12.9%)
//Threaded execution
// ThreadDohrmanAssemblyList::ThreadWork(void*): 39% 
// TPZPairStructMatrix::ThreadData::ThreadWork(void*): 12.9%
// TPZPairStructMatrix::ThreadData::ThreadAssembly1(void*): 3%

void printrusage(ostream& out, const char* msg, const struct rusage& ru)
{
  double elapsed;
  out << "============================================" << std::endl;
  out << setw(12) << " " << " " << msg << std::endl;
  elapsed = (ru.ru_utime.tv_sec * 1000.0) + (ru.ru_utime.tv_usec / 1000.0);
  out << setw(12) << setprecision(2) << 
    fixed << elapsed << " ms. user time used" << std::endl;
  elapsed = (ru.ru_stime.tv_sec * 1000.0) + (ru.ru_stime.tv_usec / 1000.0);
  out << setw(12) << setprecision(2) << 
    fixed << elapsed << " ms. system time used" << std::endl;
  out << setw(12) << ru.ru_maxrss <<   " integral max resident set size" << std::endl;
  out << setw(12) << ru.ru_ixrss <<    " integral shared text memory size" << std::endl;
  out << setw(12) << ru.ru_idrss <<    " integral unshared data size" << std::endl;
  out << setw(12) << ru.ru_isrss <<    " integral unshared stack size" << std::endl;
  out << setw(12) << ru.ru_minflt <<   " page reclaims" << std::endl;
  out << setw(12) << ru.ru_majflt <<   " page faults" << std::endl;
  out << setw(12) << ru.ru_nswap <<    " swaps" << std::endl;
  out << setw(12) << ru.ru_inblock <<  " block input operations" << std::endl;
  out << setw(12) << ru.ru_oublock <<  " block output operations" << std::endl;
  out << setw(12) << ru.ru_msgsnd <<   " messages sent" << std::endl;
  out << setw(12) << ru.ru_msgrcv <<   " messages received" << std::endl;
  out << setw(12) << ru.ru_nsignals << " signals received" << std::endl;
  out << setw(12) << ru.ru_nvcsw <<    " voluntary context switches" << std::endl;
  out << setw(12) << ru.ru_nivcsw <<   " involuntary context switches" << std::endl;
}

void printrusage2(ostream& out, const char* msg1, const struct rusage& ru1, 
		  const char* msg2, const struct rusage& ru2)
{
  double elapsed1, elapsed2;
  out << "==============================================================" << std::endl;
  out << setw(12) << msg1 << " | " << setw(12) << msg2 << " | " << std::endl;

  elapsed1 = (ru1.ru_utime.tv_sec * 1000.0) + (ru1.ru_utime.tv_usec / 1000.0);
  elapsed2 = (ru2.ru_utime.tv_sec * 1000.0) + (ru2.ru_utime.tv_usec / 1000.0);
  out << setw(12) << setprecision(2) << fixed << elapsed1 << " | "
      << setw(12) << setprecision(2) << fixed << elapsed2 
      << " | ms. user time used" << std::endl;

  elapsed1 = (ru1.ru_stime.tv_sec * 1000.0) + (ru1.ru_stime.tv_usec / 1000.0);
  elapsed2 = (ru2.ru_stime.tv_sec * 1000.0) + (ru2.ru_stime.tv_usec / 1000.0);
  out << setw(12) << setprecision(2) << fixed << elapsed1  << " | "
      << setw(12) << setprecision(2) << fixed << elapsed2 
      << " | ms. system time used" << std::endl;
  out << setw(12) << ru1.ru_maxrss <<   " | " << setw(12) << ru2.ru_maxrss <<   " | integral max resident set size" << std::endl;
  out << setw(12) << ru1.ru_ixrss <<    " | " << setw(12) << ru2.ru_ixrss <<    " | integral shared text memory size" << std::endl;
  out << setw(12) << ru1.ru_idrss <<    " | " << setw(12) << ru2.ru_idrss <<    " | integral unshared data size" << std::endl;
  out << setw(12) << ru1.ru_isrss <<    " | " << setw(12) << ru2.ru_isrss <<    " | integral unshared stack size" << std::endl;
  out << setw(12) << ru1.ru_minflt <<   " | " << setw(12) << ru2.ru_minflt <<   " | page reclaims" << std::endl;
  out << setw(12) << ru1.ru_majflt <<   " | " << setw(12) << ru2.ru_majflt <<   " | page faults" << std::endl;
  out << setw(12) << ru1.ru_nswap <<    " | " << setw(12) << ru2.ru_nswap <<    " | swaps" << std::endl;
  out << setw(12) << ru1.ru_inblock <<  " | " << setw(12) << ru2.ru_inblock <<  " | block input operations" << std::endl;
  out << setw(12) << ru1.ru_oublock <<  " | " << setw(12) << ru2.ru_oublock <<  " | block output operations" << std::endl;
  out << setw(12) << ru1.ru_msgsnd <<   " | " << setw(12) << ru2.ru_msgsnd <<   " | messages sent" << std::endl;
  out << setw(12) << ru1.ru_msgrcv <<   " | " << setw(12) << ru2.ru_msgrcv <<   " | messages received" << std::endl;
  out << setw(12) << ru1.ru_nsignals << " | " << setw(12) << ru2.ru_nsignals << " | signals received" << std::endl;
  out << setw(12) << ru1.ru_nvcsw <<    " | " << setw(12) << ru2.ru_nvcsw <<    " | voluntary context switches" << std::endl;
  out << setw(12) << ru1.ru_nivcsw <<   " | " << setw(12) << ru2.ru_nivcsw <<   " | involuntary context switches" << std::endl;
  out << "==============================================================" << std::endl;
}

int main(int argc, char *argv[])
{
  ClockTimer timer;
  ClockTimer total_timer;
  TimingAnalysis ta;

  total_timer.start();

#ifdef LOG4CXX
  log4cxx::PropertyConfigurator::configure("log4cxx_perf.cfg");
#endif

  LOG4CXX_INFO(perflog, "LOG4CXX: perflog is working with level " << 
	       perflog->getEffectiveLevel()->toString());
	
  //Options
  bool dump_mesh = false; // Dump building mesh after reading.
  int dim        = 3;     // # of dimensions
  int plevel     = 1;     // Interpolation order (gOrder)
  int maxlevel   = 5;	  // ?
  int sublevel   = 3;     // ?
#ifdef ASSEMBLE_THREADS
  int TPZStructMatrix_assemble_numthreads = ASSEMBLE_THREADS;
#warning "assemble_numthreads set to DECOMPOSE_THREADS"
#else
  int TPZStructMatrix_assemble_numthreads = 1;
#endif
#ifdef DECOMPOSE_THREADS
  int TPZDohrStructMatrix_decompose_numthreads = DECOMPOSE_THREADS;
#warning "decompose_numthreads set to DECOMPOSE_THREADS"
#else
  int TPZDohrStructMatrix_decompose_numthreads = 1;
#endif
#ifdef GNUM_THREADS
  //# of therads to build the matrix?
  TPZPairStructMatrix::gNumThreads = GNUM_THREADS;
#warning "TPZPairStructMatrix::gNumThreads set to GNUM_THREADS"
#else
  //# of therads to build the matrix?
  TPZPairStructMatrix::gNumThreads = 0;
#endif
#ifdef SUBSTR_SZ
#warning "Setting dohrstruct substructure_sz with SUBSTR_SZ"
  int substructure_sz = SUBSTR_SZ;
#else
  int substructure_sz = 8;
#endif

  LOG4CXX_INFO(perflog, "dim    = " << dim);
  LOG4CXX_INFO(perflog, "plevel = " << plevel);
  LOG4CXX_INFO(perflog, "dim    = " << dim);
  LOG4CXX_INFO(perflog, "TPZPairStructMatrix::gNumThreads         = " << TPZPairStructMatrix::gNumThreads);
  LOG4CXX_INFO(perflog, "TPZStructMatrix_assemble_numthreads      = " << TPZStructMatrix_assemble_numthreads);
  LOG4CXX_INFO(perflog, "TPZDohrStructMatrix_decompose_numthreads = " << TPZDohrStructMatrix_decompose_numthreads);
  LOG4CXX_INFO(perflog, "dohrstruct.SubStructure size             = " << substructure_sz);

  //Sets the interpolation order (gOrder).
  TPZCompEl::SetgOrder(plevel);
  
//  Construir um cubo (em vez do Predio). Dim=3 gera um cubo, dim=2 gera um quadrado.
//  maxlevel: quantidade de vezes que aplicamos o refinamento. Quanto maior mais sub-cubos
//  sublevel: nivel a ser utilizado para dividir o cubo em subdominios
//  (ex: cubo de 64x64x64 e sublevel = 8 dividiremos o problema em
//  cubos de 8x8x8)
// 			TPZGenSubStruct sub(dim,maxlevel,sublevel);
// 			cmesh = sub.GenerateMesh();
// 			cmesh->SetDimModel(dim);
// 			gmesh = cmesh->Reference();

  string input_file("../8andares02.txt");
  TIME_SEC_BEG(timer,"Reading building mesh from " << input_file);
  TPZGeoMesh* gmesh = BuildBuildingMesh(input_file.c_str());
  TIME_SEC_END(ta,timer,"Reading building mesh from " << input_file);
  
  if (dump_mesh) {
    ofstream output_f("malhaPZ.txt");
    gmesh->Print(output_f);
  }

  //goto skip_all;
  
  {
    TIME_SEC_BEG(timer,"Build cmesh");
    TPZAutoPointer<TPZCompMesh> cmesh;
    cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    InsertElasticity(cmesh);
    cmesh->AutoBuild();
    TIME_SEC_END(ta,timer,"Build cmesh");
		
    LOG4CXX_INFO(perflog, "# of equations " << cmesh->NEquations());

    TPZAutoPointer<TPZCompMesh> cmeshauto(cmesh);
    TPZDohrStructMatrix dohrstruct(cmeshauto,
				   TPZStructMatrix_assemble_numthreads,
				   TPZDohrStructMatrix_decompose_numthreads);
    dohrstruct.IdentifyExternalConnectIndexes();
  
    // FIX: remove TPZfTime?
    //		TPZfTime timetosub; // init of timer
  
    TIME_SEC_BEG(timer,"Substructuring the mesh");

    dohrstruct.SubStructure(substructure_sz);

    //EBORIN: # threads? It looks like it is already set at the constructor. (REMOVE)
    // dohrstruct.SetNumThreads(numthreads);
    TPZAutoPointer<TPZGuiInterface> gui;
    TPZFMatrix rhs(cmesh->NEquations(),1,0.);
    TIME_SEC_END(ta,timer,"Substructuring the mesh");

    //EBORIN: CreateAssemble -- dim2_2threads: 12.9%
    TIME_SEC_BEG(timer,"CreateAssemble");
    TPZAutoPointer<TPZMatrix> dohr = dohrstruct.CreateAssemble(rhs, gui);
    TIME_SEC_END(ta,timer,"CreateAssemble");

    TIME_SEC_BEG(timer,"Preconditioner");
    TPZAutoPointer<TPZMatrix> precond = dohrstruct.Preconditioner();
    TIME_SEC_END(ta,timer,"Preconditioner");
  
    TPZFMatrix diag(dohr->Rows(),1,5.);
    TPZFMatrix produto(dohr->Rows(),1);
    LOG4CXX_INFO(perflog, "# of equations " << dohr->Rows());

    TIME_SEC_BEG(timer,"Multiply started");
    dohr->Multiply(diag,produto);
    TIME_SEC_END(ta,timer,"Multiply started");
		
    TPZDohrMatrix<TPZDohrSubstructCondense> *dohrptr = 
      dynamic_cast<TPZDohrMatrix<TPZDohrSubstructCondense> *> (dohr.operator->());

    if (!dohrptr) {
      DebugStop();
    }
  
    TIME_SEC_BEG(timer,"AdjustResidual");
    dohrptr->AdjustResidual(produto);
    TIME_SEC_END(ta,timer,"AdjustResidual");
		
    TIME_SEC_BEG(timer,"Solver setup");
    diag.Zero();
    TPZStepSolver pre(precond);
    pre.SetMultiply();
    TPZStepSolver cg(dohr);
    cg.SetCG(500,pre,1.e-8,0);
    TIME_SEC_END(ta,timer,"Solver setup");

    TIME_SEC_BEG(timer,"cg.Solve");
    cg.Solve(rhs,diag);
    TIME_SEC_END(ta,timer,"cg.Solve");

    TIME_SEC_BEG(timer,"AddInternalSolution");
    dohrptr->AddInternalSolution(diag);
    TIME_SEC_END(ta,timer,"AddInternalSolution");

    TIME_SEC_BEG(timer,"Final steps");
    typedef std::list<TPZAutoPointer<TPZDohrSubstructCondense> > subtype;
    const subtype &sublist = dohrptr->SubStructures(); 
    subtype::const_iterator it = sublist.begin();
    int subcount=0;
    while (it != sublist.end()) {
      TPZFMatrix subext,subu;
      dohrptr->fAssembly->Extract(subcount,diag,subext);
      (*it)->UGlobal(subext,subu);
      TPZCompMesh *submesh = SubMesh(cmeshauto, subcount);
      submesh->LoadSolution(subu);
      subcount++;
      it++;
    }
    
    TPZAutoPointer<TPZMaterial> mat = cmeshauto->FindMaterial(1);
    int nstate = mat->NStateVariables();
    int nscal = 0, nvec = 0;
    if(nstate ==1) {
      nscal = 1;
    }
    else {
      nvec = 1;
    }
    TPZManVector<std::string> scalnames(nscal),vecnames(nvec);
    if(nscal == 1) {
      scalnames[0]="state";            
    }
    else {
      vecnames[0] = "state";
    }
    std::string postprocessname("dohrmann.vtk");
    TPZVTKGraphMesh vtkmesh(cmesh.operator->(),dim,mat,scalnames,vecnames);
    vtkmesh.SetFileName(postprocessname);
    vtkmesh.SetResolution(1);
    int numcases = 1;
    vtkmesh.DrawMesh(numcases);
    int step = 0;
    vtkmesh.DrawSolution(step, 1.);
  }

  // TODO: gerar um resultado para comparar a corretude.

  // Must delete out of scope so that Mesh autopointer does not break... :-(
  delete gmesh;

  TIME_SEC_END(ta,timer,"Final steps");

 skip_all:

  total_timer.stop();

  struct rusage self, children;
  getrusage(RUSAGE_SELF, &self);
  getrusage(RUSAGE_CHILDREN, &children);

#ifdef LOG4CXX
  std::ostringstream ostr;
  ta.share_report(ostr, total_timer.getUnits());
  LOG4CXX_INFO(perflog, ostr.str());
  printrusage2(ostr, "Self", self, "Children", children);
  LOG4CXX_INFO(perflog, ostr.str());
#else
  ta.share_report(std::cout, total_timer.getUnits());
  printrusage2(std::cout, "Self", self, "Children", children);
#endif;

  return EXIT_SUCCESS;
}
