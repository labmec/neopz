/**

 * @file TPZWaveguideModalAnalysisPML.h
 * @brief Header file for class TPZWaveguideModalAnalysisPML.\n
 */

#ifndef TPZWAVEGUIDEMODALANALYSISPML_H
#define TPZWAVEGUIDEMODALANALYSISPML_H


#include "TPZWaveguideModalAnalysis.h"
/**
 * @ingroup material
 * @brief This class implements a rectangular PML for the TPZWaveguideModalAnalysis class.
 */
class  TPZWaveguideModalAnalysisPML :
    public TPZWaveguideModalAnalysis
{
protected:
    bool fAttX{false};
    REAL fPmlBeginX{-1};
    bool fAttY{false};
    REAL fPmlBeginY{-1};
    
    STATE fAlphaMaxX{-1};
    STATE fAlphaMaxY{-1};
    STATE fDX{-1};
    STATE fDY{-1};
    TPZWaveguideModalAnalysisPML() = default;

    void ComputeSParameters(const TPZVec<REAL> &x, CSTATE&sx, CSTATE&sy) const;
public:
    //! Creates PML based on another domain region
    TPZWaveguideModalAnalysisPML(const int id,
                                 const TPZWaveguideModalAnalysis &mat);
    //! Sets information regarding the attenuation of the PML in the x-direction
    void SetAttX(const REAL pmlBegin, const STATE alpha, const REAL d);
    //! Sets information regarding the attenuation of the PML in the y-direction
    void SetAttY(const REAL pmlBegin, const STATE alpha, const REAL d);

    //! Gets the permeability of the material
    void GetPermeability(const TPZVec<REAL> &x,TPZVec<CSTATE> &ur) const override;
    //! Gets the permittivity of the material
    void GetPermittivity(const TPZVec<REAL> &x,TPZVec<CSTATE> &er) const override;
    
    TPZWaveguideModalAnalysisPML * NewMaterial() const override;
    
    std::string Name() const override { return "TPZWaveguideModalAnalysisPML";}

    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;
};

#endif
