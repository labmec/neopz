#ifndef TPBrStrainStressDataBase_H
#define TPBrStrainStressDataBase_H

#include "pzreal.h"
#include "pzvec.h"

#include <vector>

class TPBrStrainStressDataBase
{

public:
    
    /// tipos de curvas permitidas
    enum ECurveType {EEpsrSigr = 0, EI1SqJ2 = 1, EEpsaxSigax = 2, EEpsvSigv = 3, EEpsaxEpsrEpsvSigax = 4};
    
    TPBrStrainStressDataBase();
    
    virtual ~TPBrStrainStressDataBase();

    TPBrStrainStressDataBase(const TPBrStrainStressDataBase &copy) : fGlobalId(-1),
        fSig_Ax(copy.fSig_Ax), fSig_Lat(copy.fSig_Lat), fEps_Ax(copy.fEps_Ax), fEps_Lat(copy.fEps_Lat)
    {
        
    }
    
    TPBrStrainStressDataBase &operator=(const TPBrStrainStressDataBase &copy)
    {
        fGlobalId = -1;
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
    
    /// set the strain and stress vectors
    void SetStrainStress(const TPZVec<REAL> &sigax, const TPZVec<REAL> &epsax, const TPZVec<REAL> &sigr, const TPZVec<REAL> &epsr)
    {
        fSig_Ax = sigax;
        fEps_Ax = epsax;
        fSig_Lat = sigr;
        fEps_Lat = epsr;
    }

    /// gera os dados da envoltoria que sera mostrada na tela
    virtual void GenerateEnvelope(std::vector<REAL> &X, std::vector<REAL> &Y);

    /// gera os dados da envoltoria que sera mostrada na tela
    virtual void GenerateEnvelope(std::vector<REAL> &X, std::vector<REAL> &Y, REAL I1min, REAL I1max);

    /// seta os valores de A, B e C que serao usados no calculo da envoltoria
    void SetEnvelope(REAL A, REAL B, REAL C);


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

    // retorna o valor da envoltoria para o index
    REAL F1(REAL I1);

    virtual int Get_start_idx() const {
        DebugStop();
        return 0;
    }
    virtual int Get_end_idx() const {
        DebugStop();
        return fSig_Ax.size()-1;
    }
    virtual int Get_elastic_trans_idx() const {
        DebugStop();
        return 0;
    }


    virtual void Print() const
    {
        std::cout << "fGlobalId " << fGlobalId << std::endl;
    }
    
public:
    
    int fGlobalId;
    TPZVec<REAL> fSig_Ax;
    TPZVec<REAL> fSig_Lat;
    TPZVec<REAL> fEps_Ax;
    TPZVec<REAL> fEps_Lat;
    TPZVec<REAL> fF1;
    REAL fA;
    REAL fB;
    REAL fC;

};

#endif // TPBrStrainStressDataBase_H
