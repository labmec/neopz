//
// Created by natalia on 22/05/19.
//
#include "pzreal.h"
#include "pzmatrix.h"
#ifndef INTPOINTSFEM_TPZINTPOINTSDATA_H
#define INTPOINTSFEM_TPZINTPOINTSDATA_H


class TPZIntPointsData {
public:
    /** @brief Default constructor */
    TPZIntPointsData();

    /** @brief Default destructor */
    ~TPZIntPointsData();

    /** @brief Copy constructor
     * @param copy : original object
     */
    TPZIntPointsData(const TPZIntPointsData &copy);

    /** @brief operator= */
    TPZIntPointsData &operator=(const TPZIntPointsData &copy);

    void GatherSolution(TPZFMatrix<REAL> &solution, TPZFMatrix<REAL> &gather_solution);

    void ColoredAssemble(TPZFMatrix<REAL>  &nodal_forces, TPZFMatrix<REAL>  &residual);

    void SetNColor(int64_t ncolor) {
        fNColor = ncolor;
    }

    int64_t NColor() {
        return fNColor;
    }

    void SetIndexes(TPZVec<int> indexes) {
        fIndexes = indexes;
    }

    TPZVec<int> Indexes() {
        return fIndexes;
    }

    void SetIndexesColor(TPZVec<int> indexescolor) {
        fIndexesColor = indexescolor;
    }

    TPZVec<int> IndexesColor() {
        return fIndexesColor;
    }

    void SetWeight(TPZVec<REAL> weight) {
        fWeight = weight;
    }

    TPZVec<REAL> Weight() {
        return fWeight;
    }

private:
    /** @brief Number of colors of colored mesh */
    int64_t fNColor;

    /** @brief DOF indexes vector ordered by element */
    TPZVec<int> fIndexes;

    /** @brief Colored DOF indexes vector ordered by element */
    TPZVec<int> fIndexesColor;

    /** @brief Weight Vector */
    TPZVec<REAL> fWeight;


};


#endif //INTPOINTSFEM_TPZINTPOINTSDATA_H
