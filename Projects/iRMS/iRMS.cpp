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

void LinearTracerPrimal();
void LinearTracerDual();
void CheckQuarterPoint();
void CreateExampleRawData(TRMRawData &data);

int main()
{
    
    // This code use normalized piola contravariant mapping for nonlinear mappings
    HDivPiola = 0;
    


    // Running primal problem
    LinearTracerPrimal();
    
//    // Running dual problem
//    LinearTracerDual();
    
    
    std::cout << "Process complete normally." << std::endl;
    return 0;
}

void LinearTracerPrimal()
{

    TRMSpaceOdissey SpaceGenerator;
    TRMOrchestra  * SymphonyX = new TRMOrchestra;
    SymphonyX->CreateAnalPrimal(SpaceGenerator);
    std::cout << "Primal complete normally." << std::endl;
    
}

void LinearTracerDual()
{
    
    TRMSpaceOdissey SpaceGenerator;
    TRMOrchestra  * SymphonyX = new TRMOrchestra;
    SymphonyX->CreateAnalDual(SpaceGenerator);
    std::cout << "Dual complete normally." << std::endl;       
    
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

