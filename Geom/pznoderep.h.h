/**
 * @file
 * @brief Contains the implementation of the TPZNodeRep methods.
 */

#ifndef PZNODEREPHH
#define PZNODEREPHH

#include "pznoderep.h"
#include "pzlog.h"
#include <sstream>

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpzpyramid.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

static TPZLogger loggernoderep("pz.geom.noderep");


namespace pzgeom {
	
	/** @brief Constructor with node map */
	template<int N, class Topology>
	TPZNodeRep<N,Topology>::TPZNodeRep(const TPZNodeRep<N,Topology> &cp, std::map<int64_t,int64_t> & gl2lcNdMap)
    : TPZRegisterClassId(&TPZNodeRep::ClassId)
	{
		int64_t i;
		for(i = 0; i < N; i++)
		{
			if (gl2lcNdMap.find(cp.fNodeIndexes[i]) == gl2lcNdMap.end())
			{
				std::stringstream sout;
				sout << "ERROR in - " << __PRETTY_FUNCTION__
				<< " trying to clone a node " << i << " index " << cp.fNodeIndexes[i]
				<< " wich is not mapped";
				LOGPZ_ERROR(loggernoderep,sout.str().c_str());
				fNodeIndexes[i] = -1;
				continue;
			}
			fNodeIndexes[i] = gl2lcNdMap [ cp.fNodeIndexes[i] ];
		}
	}
	
	
}

#endif
