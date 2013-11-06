#ifndef TPBrStrainStressDataBase_H
#define TPBrStrainStressDataBase_H

#include "pzreal.h"
#include "pzvec.h"

#include <vector>

class TPBrStrainStressDataBase
{

public:
    
    /// tipos de curvas permitidas
    enum ECurveType {EEpsrSigr, EI1SqJ2, EEpsaxSigax, EEpsvSigv, EEpsaxEpsrEpsvSigax};
    
    TPBrStrainStressDataBase();
    
    virtual ~TPBrStrainStressDataBase()
    {
        
    }
    
    TPBrStrainStressDataBase(const TPBrStrainStressDataBase &copy) : fGlobalId(copy.fGlobalId),
        fSig_Ax(copy.fSig_Ax), fSig_Lat(copy.fSig_Lat), fEps_Ax(copy.fEps_Ax), fEps_Lat(copy.fEps_Lat)
    {
        
    }
    
    TPBrStrainStressDataBase &operator=(const TPBrStrainStressDataBase &copy)
    {
        fGlobalId = copy.fGlobalId;
        fSig_Ax = copy.fSig_Ax;
        fSig_Lat = copy.fSig_Lat;
        fEps_Ax = copy.fEps_Ax;
        fEps_Lat = copy.fEps_Lat;
        return *this;
    }
    
    /// estabelece o identificador global
    void SetGlobalId(int globalid)
    {
        fGlobalId = globalid;
    }
    
    /// retorna o identificador global do objeto
    int GlobalId()
    {
        return fGlobalId;
    }

    /// gera os dados que serao mostrados na tela
    virtual void GeneratePlot(ECurveType curvetype, std::vector<REAL> &X, std::vector<REAL> &Y);
    
    /// gera os dados que serao mostrados na tela
    virtual void GeneratePlot(ECurveType curvetype, std::vector<REAL> &X, std::vector<REAL> &Y, std::vector<REAL> &X2, std::vector<REAL> &X3);
    
    /// calcula os valores de X e Y correspondente a um indice
    void GetXY(int index, ECurveType curvetype, REAL &X, REAL &Y);
    
    /// calcula os valores de X, X2, X3 e Y correspondente a um indice
    void GetXY(int index, ECurveType curvetype, REAL &X, REAL &Y, REAL &X2, REAL &X3);
    
    /// retorno o valor da primeira invariante para o index
    REAL I1(int index);
    
    /// retorno o valor de Sq(J2) para o index
    REAL SqJ2(int index);
    
    /// retorna o valor de Eps Volumetrico para o index
    REAL Epsv(int index);

    /// retorna o valor de Sig Volumetrico para o index
    REAL Sigv(int index);
    
public:
    
    int fGlobalId;
    TPZVec<REAL> fSig_Ax;
    TPZVec<REAL> fSig_Lat;
    TPZVec<REAL> fEps_Ax;
    TPZVec<REAL> fEps_Lat;

};

#endif // TPBrStrainStressDataBase_H
