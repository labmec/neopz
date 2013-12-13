#include "TPBrStrainStressDataBase.h"
#include "TPBrDataControl.h"

TPBrStrainStressDataBase::TPBrStrainStressDataBase() : fGlobalId(-1)
{
}

TPBrStrainStressDataBase::~TPBrStrainStressDataBase()
{
    if(fGlobalId != -1)
    {
        DADOS.DeleteGlobalId(fGlobalId);
    }
}

// BEGIN ENVOLTORIA
void TPBrStrainStressDataBase::GenerateEnvelope(std::vector<REAL> &X, std::vector<REAL> &Y)
{
    X.resize(fEps_Ax.size());
    Y.resize(fEps_Ax.size());
    for (int i = 0; i<fEps_Ax.size(); i++) {
        GetEnvelope(i, X[i], Y[i]);
    }
}

void TPBrStrainStressDataBase::SetEnvelope(REAL A, REAL B, REAL C)
{
    fA = A;
    fB = B;
    fC = C;
    qDebug() <<"------------------------->"<< A << B << C;
}

void TPBrStrainStressDataBase::GetEnvelope(int index, REAL &X, REAL &Y)
{
    X = I1(index);
    Y = F1(index);
}

// END ENVOLTORIA

/// gera os dados que serao mostrados na tela
void TPBrStrainStressDataBase::GeneratePlot(ECurveType curvetype, std::vector<REAL> &X, std::vector<REAL> &Y)
{
    X.resize(fEps_Ax.size());
    Y.resize(fEps_Ax.size());
    for (int i = 0; i<fEps_Ax.size(); i++) {
        GetXY(i, curvetype, X[i], Y[i]);
    }
}
/// gera os dados que serao mostrados na tela
void TPBrStrainStressDataBase::GeneratePlot(ECurveType curvetype, std::vector<REAL> &X, std::vector<REAL> &Y, std::vector<REAL> &X2, std::vector<REAL> &X3)
{
    X.resize(fEps_Ax.size());
    Y.resize(fEps_Ax.size());
    X2.resize(fEps_Ax.size());
    X3.resize(fEps_Ax.size());
    for (int i = 0; i<fEps_Ax.size(); i++) {
        GetXY(i, curvetype, X[i], Y[i], X2[i], X3[i]);
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
	    break;
	case EEpsaxSigax:
        X = this->fEps_Ax[index];
        Y = this->fSig_Ax[index];
            break;	    
	case EEpsvSigv:
	    X = Epsv(index);
	    Y = Sigv(index);
	    break;
    default:
            DebugStop();
            break;
    }
}

/// calcula os valores de X, X2, X3 e Y correspondente a um indice
void TPBrStrainStressDataBase::GetXY(int index, ECurveType curvetype, REAL &X, REAL &Y, REAL &X2, REAL &X3)
{
    switch (curvetype) {
	case EEpsaxEpsrEpsvSigax:
	    X = this->fEps_Ax[index];
	    Y = this->fSig_Ax[index];
	    X2 = this->fEps_Lat[index];
	    X3 = Epsv(index);
	    break;
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

/// retorna o valor de Eps Volumetrico para o index
REAL TPBrStrainStressDataBase::Epsv(int index)
{
    REAL epsvol = fEps_Ax[index] + 2 * fEps_Lat[index];
    return epsvol;
}

/// retorna o valor de Sig Volumetrico para o index
REAL TPBrStrainStressDataBase::Sigv(int index)
{
    REAL sigvol = fSig_Ax[index] + 2 * fSig_Ax[index];
    return sigvol;
}

// retorna o valor da envoltoria para o index
REAL TPBrStrainStressDataBase::F1(int index)
{
    REAL I1 = this->I1(index);
    REAL F1 = fA - fC*exp(fB*I1);
    return F1;
}

