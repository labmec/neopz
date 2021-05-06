#include "TPZMatTransientSingleSpace.h"
#include "TPZLagrangeMultiplier.h"
//just to test if the interfaces for discontinuous materials have not changed
template class TPZMatTransientSingleSpace<TPZLagrangeMultiplier<STATE>>;