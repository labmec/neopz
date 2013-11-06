#include "TPBrStrainStressDataBase.h"

TPBrStrainStressDataBase::TPBrStrainStressDataBase()
{
}

/// gera os dados que serao mostrados na tela
void TPBrStrainStressDataBase::GeneratePlot(ECurveType curvetype, std::vector<REAL> &X, std::vector<REAL> &Y)
{
    X.resize(fEps_Ax.size());
    Y.resize(fEps_Ax.size());
    for (int i = 0; i<fEps_Ax.size(); i++) {
        GetXY(i, curvetype, X[i], Y[i]);
    }
}

/// calcula os valores de X e Y correspondente a um indice
void TPBrStrainStressDataBase::GetXY(int index, ECurveType curvetype, REAL &X, REAL &Y)
{
    
    switch (curvetype) {
        case EEpsrSigr:
            X = this->fEps_Lat[index];
            Y = this->fSig_Lat[index];
            break;
        case EI1SqJ2:
            X = I1(index);
            Y = SqJ2(index);
        default:
            DebugStop();
            break;
    }
}

/// retorno o valor da primeira invariante para o index
REAL TPBrStrainStressDataBase::I1(int index)
{
    REAL I1 = this->fSig_Ax[index]+2.*this->fSig_Lat[index];
    return I1;
}

/// retorno o valor de Sq(J2) para o index
REAL TPBrStrainStressDataBase::SqJ2(int index)
{
    REAL I1 = this->I1(index);
    REAL sigdesv[3] ={fSig_Ax[index]-I1/3.,fSig_Lat[index]-I1/3.,fSig_Lat[index]-I1/3.};
    REAL J2 = (sigdesv[0]*sigdesv[0]+sigdesv[1]*sigdesv[1]+sigdesv[2]*sigdesv[2])/2.;
    return sqrt(J2);
}


