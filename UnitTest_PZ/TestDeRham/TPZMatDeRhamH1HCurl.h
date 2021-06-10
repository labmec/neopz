/**

 * @file TPZMatDeRhamH1HCurl.h
 * @brief Header file for class TPZMatDeRhamH1HCurl.\n
 */

#ifndef TPZMATDERHAMH1HCURL_H
#define TPZMATDERHAMH1HCURL_H


#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
/**
 * @ingroup material
 * @brief This class checks if the gradients of H1 are contained in the HCurl 
 * approximation space.
 */
class  TPZMatDeRhamH1HCurl :
    public TPZMatBase<STATE,TPZMatCombinedSpacesT<STATE>>
{
    using TBase = TPZMatBase<STATE,TPZMatCombinedSpacesT<STATE>>;
public:
    //! Constructor taking dimension and material id
    inline TPZMatDeRhamH1HCurl(int id, int dim);
    
    TPZMatDeRhamH1HCurl * NewMaterial() const override;
    
    std::string Name() const override { return "TPZMatDeRhamH1HCurl"; }
    
    /** @brief Returns the integrable dimension of the material */
    inline int Dimension() const override {return fDim;}

    [[nodiscard]] int NStateVariables() const override{return 1;}
    
    /** @brief Index of the HCurl approximation space*/
    [[nodiscard]] inline static constexpr int HCurlIndex() { return fHCurlMeshIndex;}
    /** @brief Index of the H1 approximation space*/
    [[nodiscard]] inline static constexpr int H1Index() { return fH1MeshIndex;}
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
    TPZMatDeRhamH1HCurl() = default;
    int fDim{-1};
    static constexpr int fH1MeshIndex{0};
    static constexpr int fHCurlMeshIndex{1};
};

inline TPZMatDeRhamH1HCurl::TPZMatDeRhamH1HCurl(int id, int dim) :
    TBase(id), fDim(dim)
{}

#endif
