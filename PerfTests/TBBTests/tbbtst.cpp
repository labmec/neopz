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

/**
 * Commmand Line Options
 */
clarg::argInt       plevel("-p","Polinomial Order", 2);
clarg::argInt       num_matrices("-nmat", "Number of Matrices", 64);
clarg::argInt       num_threads("-nt", "Number of Threads", 0);
clarg::argString    input_file("-mc", "Cubo Input File", "cube1.txt");

#ifdef USING_TBB
class tbb_work {
public:
    tbb_work(std::vector <TPZSkylMatrix<REAL>* > *m): matrices(m) {};
    std::vector <TPZSkylMatrix<REAL>* > *matrices;
    
    void serial() {
        for (int i=0; i<matrices->size(); i++) {
            (*matrices)[i]->Decompose_Cholesky();
            printf("-");
        }
    }
    void operator()(const tbb::blocked_range<int>& range) const {
        for( int i=range.begin(); i!=range.end(); ++i ) {
            (*matrices)[i]->Decompose_Cholesky();
            printf("-");
        }
    }
    
};
#endif
int main (int argc, char **argv)
{
    /**
     * Read and Parse the Command Line Arguments
     */
	if (clarg::parse_arguments(argc, argv)) {
		cerr << "Error when parsing the arguments!" << endl;
        return 1;
    }
    /**
     * Read and Create a Skyline Matrix
     */
    TPZAutoPointer<TPZMatrix<REAL> > ret = Input::CreateCuboSkyMatrix(input_file.get_value(),plevel.get_value());
    /**
     * Cast from AutoPointer to TPZSkylMatrix
     */
    TPZSkylMatrix<REAL> *orig = dynamic_cast<TPZSkylMatrix<REAL> *> (ret.operator->());
    /**
     * Printing the Memory Footprint of the Original Skyline Matrix
     */
    printf ("Memory Footprint SkylMatrix: %ld Kbytes\n.", (orig->MemoryFootprint()/1024));
    /**
     * Vector to Store Copies of the Original Matrix
     */
    std::vector <TPZSkylMatrix<REAL>* > matrices;
    /**
     * Serial Copy
     */
    matrices.resize(num_matrices.get_value());
    
    for (int i=0; i<num_matrices.get_value(); i++)
        matrices[i] = new TPZSkylMatrix<REAL>(*orig);
    
    /**
     * Choose Between Serial and Parallel Processing */
    if (num_threads.get_value() == 0) {
        /**
         * Serial Decomposition
         */
        for (int i=0; i<num_matrices.get_value(); i++) {
            //matrices[i]->SetIsDecomposed(0);
            matrices[i]->Decompose_Cholesky();
            printf("-");
        }
    } else {
        
#ifdef USING_TBB
        tbb::task_scheduler_init init(num_threads.get_value());
        
        tbb_work work(&matrices);
        
        //work.serial();
        tbb::parallel_for(tbb::blocked_range<int>(0, matrices.size()), work);
#endif
    }
    
} /* main */
