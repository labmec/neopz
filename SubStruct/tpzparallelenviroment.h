/**
 * @file
 * @brief Contains the TPZParallelEnviroment class which store the 
 * parallel enviroment variables
 * @since 2/21/14.
 */

#ifndef TPZPARALLELENVIROMENT_H
#define TPZPARALLELENVIROMENT_H

#ifdef USING_TBB
#include <tbb/partitioner.h>
#endif

class TPZParallelEnviroment 
{
	public:
		/** @brief */
		TPZParallelEnviroment() {};

		/** @brief */
		~TPZParallelEnviroment() {};

#ifdef USING_TBB
		tbb::affinity_partitioner fSubstructurePartitioner;
#endif

};

// Extern Variable
extern TPZParallelEnviroment pzenviroment;

#endif // TPZPARALLELENVIROMENT_H
