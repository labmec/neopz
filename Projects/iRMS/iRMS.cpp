#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZIntQuadQuarterPoint.h"
#include <time.h>
#include "pzgmesh.h"
#include "TPZRefPatternTools.h"
#include "pzmaterial.h"

#include "TRMRawData.h"
#include "TRMSimworxMeshGenerator.h"
#include "TRMOrchestra.h"
#include "TRMSpaceOdissey.h"

#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"


// iMRS (i) innovatory Multiscale Reservoir Simulator

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.iRMS"));
#endif

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void MultiScaleSimulation();

int main()
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    // This code use normalized piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
    TPZMaterial::gBigNumber = 1.0e12;
    
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    // Running whole process
    MultiScaleSimulation();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    std::cout  << "iMRS:: Overal execution time = " << (t2-t1) << std::endl;
#endif
    
    std::cout << "iMRS:: Process complete normally." << std::endl;
    return 0;
}

void MultiScaleSimulation()
{
    // Materials ids and boundary settings
    TPZAutoPointer<TRMRawData> RawData  = new TRMRawData;
    
    //  Dimension on gmsh reservoir    
    bool Is3DGeometry = true;
    
    bool IsSinglePhaseQ = true;
    if(IsSinglePhaseQ){
//        RawData->SinglePhaseReservoirHMM(Is3DGeometry); // FEM and HMM chapter
        RawData->SinglePhaseReservoir(Is3DGeometry); // Single-phase flow
    }
    else{
        RawData->CaseTracerTransport(Is3DGeometry); // Case 1 and 2 Tracer transport
//        RawData->TwoPhaseWaterOilReservoir(Is3DGeometry); // Two-phase flow
    }
    
    TRMSimulationData * SimData = new TRMSimulationData;
    SimData->SetRawData(RawData);
    TRMOrchestra  * SymphonyX           = new TRMOrchestra;
    SymphonyX->SetSimulationData(SimData);
    
    SymphonyX->SetSegregatedQ(true);
    SymphonyX->CreateSegregatedAnalysis(true); //  Static Solution
    SymphonyX->RunStaticProblem();
    SymphonyX->CreateSegregatedAnalysis(false);  // Evolutionary Solution
    SymphonyX->RunEvolutionaryProblem();

    std::cout << "Dual complete normally." << std::endl;
    
}

