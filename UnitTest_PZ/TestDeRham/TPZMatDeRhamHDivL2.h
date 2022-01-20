/**

 * @file TPZMatDeRhamH1HCurl.h
 * @brief Header file for class TPZMatDeRhamH1HCurl.\n
 */

#ifndef TPZMATDERHAMHDIVL2_H
#define TPZMATDERHAMHDIVL2_H


#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
/**
 * @ingroup material
 * @brief This class checks if the gradients of H1 are contained in the HCurl 
 * approximation space.
 */
class  TPZMatDeRhamHDivL2 :
    public TPZMatBase<STATE,TPZMatCombinedSpacesT<STATE>>
{
    using TBase = TPZMatBase<STATE,TPZMatCombinedSpacesT<STATE>>;
public:
    //! Constructor taking dimension and material id
    inline TPZMatDeRhamHDivL2(int id, int dim);
    
    TPZMatDeRhamHDivL2 * NewMaterial() const override;
    
    std::string Name() const override { return "TPZMatDeRhamHDivL2"; }
    
    /** @brief Returns the integrable dimension of the material */
    inline int Dimension() const override {return fDim;}

    [[nodiscard]] int NStateVariables() const override{return 1;}
    
    /** @brief Index of the HCurl approximation space*/
    [[nodiscard]] inline static constexpr int HDivIndex() { return fHDivMeshIndex;}
    /** @brief Index of the H1 approximation space*/
    [[nodiscard]] inline static constexpr int L2Index() { return fL2MeshIndex;}
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
    TPZMatDeRhamHDivL2() = default;
    int fDim{-1};
    static constexpr int fHDivMeshIndex{0};
    static constexpr int fL2MeshIndex{1};
};

inline TPZMatDeRhamHDivL2::TPZMatDeRhamHDivL2(int id, int dim) :
    TBase(id), fDim(dim)
{}

#endif
