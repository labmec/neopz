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
    bool fAttX;
    STATE fPmlBeginX;
    bool fAttY;
    STATE fPmlBeginY;
    
    STATE fAlphaMax;
    STATE fD;
    TPZWaveguideModalAnalysisPML() = default;

    void ComputeSParameters(const TPZVec<REAL> &x, CSTATE&sx, CSTATE&sy);
public:
    TPZWaveguideModalAnalysisPML(const int id,const TPZWaveguideModalAnalysis &mat,
                    const bool &attX, REAL &pmlBeginX,
                    const bool &attY, REAL &pmlBeginY,
                    const REAL &alphaMax, const REAL &d);

    TPZWaveguideModalAnalysisPML * NewMaterial() const override;
    
    std::string Name() const override { return "TPZWaveguideModalAnalysisPML";}
    /**
       @name ContributeMethods
       @{
    */
    void Contribute(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;
    /**@}*/
    /**
       @name SolutionMethods
       @{*/
    //! Computes the solution at an integration point
    void Solution(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
                  int var, TPZVec<CSTATE> &solout) override;
    /**@}*/

    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;
};

#endif
