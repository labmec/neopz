#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZRefPatternTools.h"

#include "TRMRawData.h"
#include "TRMSimworxMeshGenerator.h"


// iRMS (i rise for / innovatory Reservoir Muli-scale Simulator)

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.iRMS"));
#endif

void LinearTracer();

void CreateExampleRawData(TRMRawData &data)
{
    data.fLw = 500.;
    data.fHasLiner = true;
    data.fHasCasing = true;
    
    data.fReservoirWidth = 500.;
    data.fReservoirLength = 1000.;
    data.fReservoirHeight = 50.;
    data.fProdVertPosition = 25;
}

int main()
{
    
    // This code use normalized piola contravariant mapping for nonlinear mappings
    // HDivPiola = 0;
    
    gRefDBase.ReadRefPatternDBase("../RefPatterns.rpt");
    
    TRMRawData rawdata;
    
    CreateExampleRawData(rawdata);
    
    TRMSimworxMeshGenerator meshGen;
    TPZGeoMesh *gmesh = meshGen.CreateSimworxGeoMesh(rawdata);
    
    std::cout << " Process complete normally." << std::endl;
    return 0;
}


void LinearTracer()
{
    // Simulation Data in SI units
    
    
    
    // Rock petrophysical description Data in SI units
    
    
    
    // Fluid description Data in SI units
    
}