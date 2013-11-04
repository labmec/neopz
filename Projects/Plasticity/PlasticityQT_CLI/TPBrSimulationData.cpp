#include "TPBrSimulationData.h"
#include "TPBrStrainStressDataBase.h"
#include "TPBrDataControl.h"

TPBrSimulationData::TPBrSimulationData()
{
    start_idx = -1;
    end_idx = -1;
    med_idx = -1;
}


TPBrSimulationData::TPBrSimulationData(int startidx, int endidx, int medidx)
{
    start_idx = startidx;
    end_idx = endidx;
    med_idx = medidx;
}
