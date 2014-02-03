/*******************************************************************************
 *   Copyright (C) 2014 by:                                                    *
 *   Gilvan Vieira (gilvandsv@gmail.com)                                       *
 *                                                                             *
 *   This program is free software; you can redistribute it and/or modify      *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation; either version 2 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program; if not, write to the                             *
 *   Free Software Foundation, Inc.,                                           *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                 *
 ******************************************************************************/

/**
 * Standard Input/Output from C
 */
#include <cstdio>
/**
 * Functions to Read Data
 */
#include "input.h"
/**
 * Command Line Arguments - Borin
 */
#include "arglib.h"
/**
 * STL Vector
 */
#include <vector>
/**
 * Threading Building Blocks Headers
 */
#ifdef USING_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#endif

#ifdef USING_HWLOC
#include <hwloc.h>
hwloc_topology_t topology;
#endif

/**
 * Time
 */
#include "run_stats_table.h"
/**
 * Commmand Line Options
 */
clarg::argInt       plevel("-p","Polinomial Order", 1);
clarg::argInt       num_matrices("-nmat", "Number of Matrices", 64);
clarg::argInt       num_threads("-nt", "Number of Threads", 1);
clarg::argString    input_file("-mc", "Cubo Input File", "cube1.txt");
/**
 * Types of Executions:
 *  1 - Serial
 *  2 - Parallel with Main Memory Allocation
 *  3 - Parallel with On Thread Memory Allocation
 *  4 - Parallel with ReallocForNuma (Borin)
 *  5 - Parallel with HWloc Realloc
 */
clarg::argInt       exec_type("-type", "Execution Type", 1);

RunStatsTable dec_rst   ("-perf", "Decompose Cholesky Statistics");

#ifdef USING_TBB
class TBBWork {
protected:
    std::vector <TPZFMatrix<REAL>* > *matrices;
    TPZFMatrix<REAL> *orig;
public:
    TBBWork(std::vector <TPZFMatrix<REAL>* > *m, TPZFMatrix<REAL> *o): matrices(m), orig(o) {};
    virtual void operator()(const tbb::blocked_range<int>& range) const = 0;
    virtual void execute() = 0;
};

class MainTBB : public TBBWork {
public:
    MainTBB(std::vector <TPZFMatrix<REAL>* > *m, TPZFMatrix<REAL> *o) : TBBWork(m, o) {};
    void operator()(const tbb::blocked_range<int>& range) const {
        for( int i=range.begin(); i!=range.end(); ++i ) {
            (*matrices)[i]->Decompose_Cholesky();
        }
    }
    void execute() {
        for (int i=0; i<num_matrices.get_value(); i++) {
            (*matrices)[i] = new TPZFMatrix<REAL>(*orig);
        }
        tbb::parallel_for(tbb::blocked_range<int>(0, matrices->size()), *this);
    }
};

class OnThreadTBB : public TBBWork {
public:
    OnThreadTBB(std::vector <TPZFMatrix<REAL>* > *m, TPZFMatrix<REAL> *o) : TBBWork(m, o) {};
    void operator()(const tbb::blocked_range<int>& range) const {
        for( int i=range.begin(); i!=range.end(); ++i ) {
            (*matrices)[i] = new TPZFMatrix<REAL>(*orig);
            (*matrices)[i]->Decompose_Cholesky();
        }
    }
    void execute() {
        tbb::parallel_for(tbb::blocked_range<int>(0, matrices->size()), *this);
    }
};

class ReallocNumaTBB : public TBBWork {
public:
    ReallocNumaTBB(std::vector <TPZFMatrix<REAL>* > *m, TPZFMatrix<REAL> *o) : TBBWork(m, o) {};
    void operator()(const tbb::blocked_range<int>& range) const {
        for( int i=range.begin(); i!=range.end(); ++i ) {
            /* Realloc - Using First Touch */
            TPZFMatrix<REAL> *copy = new TPZFMatrix<REAL>(*(*matrices)[i]);
            TPZFMatrix<REAL> *old = (*matrices)[i];
            (*matrices)[i] = copy;
            delete old;
            /* Decompose */
            (*matrices)[i]->Decompose_Cholesky();
        }
    }
    void execute() {
        for (int i=0; i<num_matrices.get_value(); i++) {
            (*matrices)[i] = new TPZFMatrix<REAL>(*orig);
        }
        tbb::parallel_for(tbb::blocked_range<int>(0, matrices->size()), *this);
    }
};

#ifdef USING_HWLOC
class HWlocTBB : public TBBWork {
public:
    HWlocTBB(std::vector <TPZFMatrix<REAL>* > *m, TPZFMatrix<REAL> *o) : TBBWork(m, o) {};
    
    void migrate(TPZFMatrix<REAL>* m) const {
        hwloc_bitmap_t set = hwloc_bitmap_alloc();
		hwloc_get_cpubind(topology, set, HWLOC_CPUBIND_THREAD);
		hwloc_get_last_cpu_location(topology, set, HWLOC_CPUBIND_THREAD);
		hwloc_bitmap_singlify(set);
        hwloc_set_area_membind ( topology, (const void*)m->Adress(), m->MemoryFootprint(), (hwloc_const_cpuset_t)set, HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_MIGRATE );
    }
    
    void operator()(const tbb::blocked_range<int>& range) const {
        for( int i=range.begin(); i!=range.end(); ++i ) {
            /* Migrate Data Using HWLoc */
            migrate((*matrices)[i]);
            /* Decompose */
            (*matrices)[i]->Decompose_Cholesky();
        }
    }
    void execute() {
        for (int i=0; i<num_matrices.get_value(); i++) {
            (*matrices)[i] = new TPZFMatrix<REAL>(*orig);
        }
        tbb::parallel_for(tbb::blocked_range<int>(0, matrices->size()), *this);
    }
};
#endif /* USING_HWLOC */
#endif /* USING_TBB */

int main (int argc, char **argv)
{
#ifdef USING_HWLOC
    hwloc_topology_init(&topology);
	hwloc_topology_load(topology);
#endif
    
    /**
     * Read and Parse the Command Line Arguments
     */
    if (clarg::parse_arguments(argc, argv)) {
        cerr << "Error when parsing the arguments!" << endl;
        return 1;
    }
    
#ifdef USING_TBB
    tbb::task_scheduler_init init(num_threads.get_value());
    printf ("-- Number of Threads: %d\n", num_threads.get_value());
#endif
    /**
     * Read and Create a Skyline Matrix
     */
    TPZAutoPointer<TPZMatrix<REAL> > ret = Input::CreateCuboSkyMatrix(input_file.get_value(),plevel.get_value());
    /**
     * Cast from AutoPointer to TPZFMatrix
     */
    TPZFMatrix<REAL> *orig = new TPZFMatrix<REAL>(*ret.operator->());
    /**
     * Printing the Memory Footprint of the Original Skyline Matrix
     */
    printf ("Memory Footprint TPZFMatrix: %.2f Mbytes\n.", ((double)orig->MemoryFootprint()/(1024.0*1024.0)));
    /**
     * Vector to Store Copies of the Original Matrix
     */
    std::vector <TPZFMatrix<REAL>* > matrices;
    /**
     * Serial Copy
     */
    matrices.resize(num_matrices.get_value(), 0);
    

    dec_rst.start();
    
    TBBWork *work = 0;
    switch (exec_type.get_value()) {
            /* - Serial */
        case 1:
            for (int i=0; i<num_matrices.get_value(); i++) {
                matrices[i] = new TPZFMatrix<REAL>(*orig);
                matrices[i]->Decompose_Cholesky();
            }
            break;
        case 2:
            /* - Parallel with Main Thread Memory Allocation */
            work = new MainTBB(&matrices, orig);
            break;
        case 3:
            /* - Parallel with On Thread Memory Allocation */
            work = new OnThreadTBB(&matrices, orig);
            break;
        case 4:
            /* - Parallel with ReallocForNuma (Borin) */
            work = new ReallocNumaTBB(&matrices, orig);
            break;
        case 5:
            /* - Parallel with HWloc Realloc */
            work = new HWlocTBB(&matrices, orig);
            break;
        default:
            printf("Option not available!\n");
            return 1;
    }
    
    if(work) work->execute();
    
    delete orig;
    for (int i=0; i<num_matrices.get_value(); i++)
        if(matrices[i]) delete matrices[i];
    
    dec_rst.stop();
    
} /* main */
