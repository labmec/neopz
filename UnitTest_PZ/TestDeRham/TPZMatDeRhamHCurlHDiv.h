/**

 * @file TPZMatDeRhamHCurlHDiv.h
 * @brief Header file for class TPZMatDeRhamHCurlHDiv.\n
 */

#ifndef TPZMATDERHAMHCURLHDIV_H
#define TPZMATDERHAMHCURLHDIV_H


#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
/**
 * @ingroup material
 * @brief This class checks if the gradients of H1 are contained in the HCurl 
 * approximation space.
 */
class  TPZMatDeRhamHCurlHDiv :
    public TPZMatBase<STATE,TPZMatCombinedSpacesT<STATE>>
{
    using TBase = TPZMatBase<STATE,TPZMatCombinedSpacesT<STATE>>;
public:
    //! Constructor taking dimension and material id
    inline TPZMatDeRhamHCurlHDiv(int id, int dim);
    
    TPZMatDeRhamHCurlHDiv * NewMaterial() const override;
    
    std::string Name() const override { return "TPZMatDeRhamHCurlHDiv"; }
    
    /** @brief Returns the integrable dimension of the material */
    inline int Dimension() const override {return fDim;}

    [[nodiscard]] int NStateVariables() const override{return 1;}
    
    /** @brief Index of the HCurl approximation space*/
    [[nodiscard]] inline static constexpr int HCurlIndex() { return fHCurlMeshIndex;}
    /** @brief Index of the H1 approximation space*/
    [[nodiscard]] inline static constexpr int HDivIndex() { return fHDivMeshIndex;}
    /**
       @name ContributeMethods
       @{
    */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                      TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override{}
    /**@}*/
    /**
       @name SolutionMethods
       @{*/
    //! Variable index of a given solution
    int VariableIndex(const std::string &name) const override;
    //! Number of variables associated with a given solution
    int NSolutionVariables(int var) const override;
    //! Computes the solution at an integration point
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                  int var, TPZVec<STATE> &solout) override;
    /**@}*/

protected:
    TPZMatDeRhamHCurlHDiv() = default;
    int fDim{-1};
    static constexpr int fHDivMeshIndex{1};
    static constexpr int fHCurlMeshIndex{0};
};

inline TPZMatDeRhamHCurlHDiv::TPZMatDeRhamHCurlHDiv(int id, int dim) :
    TBase(id), fDim(dim)
{
    if(dim != 3) {
        DebugStop();
    }
}

#endif
