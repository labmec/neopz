
#include "pzgmesh.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "pzbndcond.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "TPZVTKGeoMesh.h"
#include "tpzautopointer.h"
#include "pzgeoquad.h"
#include "tpzgeoelrefpattern.h"
#include "TPZFrontStructMatrix.h"

#include "ReservoirData.h"
#include "SimulationData.h"


#include <time.h>

#include "fad.h"

REAL mypow(const REAL &a, const int &n)
{
	if(n == 0) return 1.;
	return a*mypow(a,n-1);
}

void FillProblemData();


int main()
{
	FillProblemData();
    
    {
        
    }
    
	return 0;
}

void FillProblemData()
{
    
    
}

