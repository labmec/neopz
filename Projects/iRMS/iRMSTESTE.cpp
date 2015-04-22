#include "pzlog.h"
#include "tpzautopointer.h"
#include <time.h>

// iRMS (i rise for / innovatory Reservoir Muli-scale Simulator)

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.iRMS"));
#endif

int main()
{
    
    // This code use normalized piola contravariant mapping for nonlinear mappings
    // HDivPiola = 0;
    
    // Simulation Data SI units
    
    // Rock petrophysical description Data SI units
    
    // Fluid description Data SI units
    
    
    std::cout << " Process complete normally." << std::endl;
    return 0;
}


