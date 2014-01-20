#include "TPBrSimulationData.h"
#include "TPBrStrainStressDataBase.h"
#include "TPBrDataControl.h"

TPBrSimulationData::TPBrSimulationData()
{
    fstart_idx = -1;
    fend_idx = -1;
    fMedicao_idx = -1;
    felastic_idx = -1;
}


TPBrSimulationData::TPBrSimulationData(int startidx, int endidx, int medidx)
{
    fstart_idx = startidx;
    fend_idx = endidx;
    fMedicao_idx = medidx;
    felastic_idx = -1;
}

TPBrSimulationData::~TPBrSimulationData()
{
    
}
