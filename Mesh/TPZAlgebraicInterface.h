#ifndef TPZALGEBRAICINTERFACE_H
#define TPZALGEBRAICINTERFACE_H

#include "pzcompel.h"
#include <map>
#include <vector>

/// @brief A computational element that implements an algebraic interface, which can be used to apply constraints to the solution of a problem
class TPZAlgebraicInterface : public TPZCompEl
{
public:
    TPZAlgebraicInterface();
    TPZAlgebraicInterface(TPZCompMesh &mesh, TPZGeoEl *reference);
    TPZAlgebraicInterface(TPZCompMesh &mesh, const TPZCompEl &copy);
    TPZAlgebraicInterface(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<int64_t,int64_t> &gl2lcElMap);
    virtual ~TPZAlgebraicInterface();

    TPZCompEl *Clone(TPZCompMesh &mesh) const override;
    TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                            std::map<int64_t,int64_t> &gl2lcConMap,
                            std::map<int64_t,int64_t> &gl2lcElMap) const override;

    virtual int NConnects() const override;
    virtual int64_t ConnectIndex(int i) const override;
    virtual int Dimension() const override;
    virtual void SetConnectIndex(int inode, int64_t index) override;
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override {

    }

    void SetMultiplier(REAL multiplier) {
        fMultiplier = multiplier;
    }

    REAL Multiplier() const {
        return fMultiplier;
    }

    virtual void CalcStiff(TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef) override;
    virtual void CalcStiff(TPZElementMatrixT<CSTATE> &ek, TPZElementMatrixT<CSTATE> &ef) override;

    int ClassId() const override;

protected:

  // connect indexes of the element
  std::array<int64_t,2> fConnectIndexes;

  /// @brief multiplier of the Lagrange multiplier, usually plus or minus one
  REAL fMultiplier = 1.0;

  // algebraic interface does not have a geometric element
  // the algebraic interface will establish a one to one relation between the connects of the first connect and the second connect
    template<class TVar>
    void CalcStiffInternal(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef);
};

#endif // TPZALGEBRAICINTERFACE_H
