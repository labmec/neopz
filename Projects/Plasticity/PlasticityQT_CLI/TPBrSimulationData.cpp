#include "TPBrSimulationData.h"
#include "TPBrStrainStressDataBase.h"
#include "TPBrDataControl.h"

TPBrSimulationData::TPBrSimulationData()
{
    fstart_idx = -1;
    fend_idx = -1;
    fmed_idx = -1;
}


TPBrSimulationData::TPBrSimulationData(int startidx, int endidx, int medidx)
{
    fstart_idx = startidx;
    fend_idx = endidx;
    fmed_idx = medidx;
}

TPBrSimulationData::~TPBrSimulationData()
{
    
}