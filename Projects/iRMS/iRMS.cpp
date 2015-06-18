#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZIntQuadQuarterPoint.h"
#include <time.h>
#include "pzgmesh.h"
#include "TPZRefPatternTools.h"

#include "TRMRawData.h"
#include "TRMSimworxMeshGenerator.h"
#include "TRMOrchestra.h"
#include "TRMSpaceOdissey.h"


// iRMS (i rise for / innovatory Reservoir Muli-scale Simulator)

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.iRMS"));
#endif

void LinearTracer();
void CheckQuarterPoint();
void CreateExampleRawData(TRMRawData &data);

int main()
{
    
    // This code use normalized piola contravariant mapping for nonlinear mappings
    // HDivPiola = 0;
    
    gRefDBase.ReadRefPatternDBase("../RefPatterns.rpt");
    TRMRawData rawdata;
    CreateExampleRawData(rawdata);
    TRMSimworxMeshGenerator meshGen;
    
    TRMSpaceOdissey * spacegenerator = new TRMSpaceOdissey;
    
    TPZManVector<int,2> dx(2,10), dy(2,1), dz(2,1);
    dx[0] = 100;
    dy[0] = 100;
    dz[0] = 100;

    spacegenerator->CreateGeometricBoxMesh(dx, dy, dz);
    spacegenerator->PrintGeometry();
    
    std::cout << "Process complete normally." << std::endl;
    return 0;
}

void LinearTracer()
{
    // Simulation Data in SI units
    
    
    
    // Rock petrophysical description Data in SI units
    
    
    
    // Fluid description Data in SI units
    
}

void CheckQuarterPoint()
{
    TPZIntQuadQuarterPoint qt(10);
    qt.SetCorner(1);
//    TPZIntQuad qt(60);
    int np = qt.NPoints();
    TPZManVector<REAL,2> pt(2,0.);
    std::cout << "Numpoints = " << np << std::endl;
    REAL weight;
    REAL integral = 0.;
    for (int ip = 0; ip<np; ip++) {
        qt.Point(ip, pt, weight);
        std::cout << "ip " << ip << " pt " << pt << " weight " << weight << std::endl;
        REAL r = sqrt((pt[0]-1)*(pt[0]-1)+(pt[1]+1)*(pt[1]+1));
        integral += weight/r;
    }
    std::cout << "Integral " << integral << std::endl;
    
}

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