/**
 * \file
 * @brief Contains the TPZBndCond class which implements a boundary condition for TPZMaterial objects.
 */

#ifndef TPZBNDCOND_H
#define TPZBNDCOND_H

#include "TPZMaterial.h"

#include <iostream>

class TPZMaterial;
/**
 * @ingroup material
 * @brief This class defines a type-agnostic interface for boundary condition in TPZMaterial objects.
 * @note Type-related methods are declared and defined in TPZBndCondT
 */
class TPZBndCond : public virtual TPZSavable {
	
	friend class TPZMaterial;
protected:
    //! Default constructor. 
	TPZBndCond() :
        TPZRegisterClassId(&TPZBndCond::ClassId) {}
    //! Constructor taking type
    explicit TPZBndCond(int type) : fType(type) {}
public :
    //! Set associated material. It also gets its attributes.
    virtual void SetMaterial(TPZMaterial * mat) = 0;
    //! Get associated material.
    [[nodiscard]] virtual TPZMaterial * Material();
    //! Get type of BC.
    [[nodiscard]] int Type() const { return fType; }
	//! Set type of BC.
	void SetType(int type){ this->fType = type; }
    //! Material identifier
    [[nodiscard]] virtual int Id() const = 0;

    [[nodiscard]] int ClassId() const override;
    //! Whether the boundary condition has a forcing function
    [[nodiscard]] virtual bool HasForcingFunctionBC() const = 0;
    void Print(std::ostream &out = std::cout) const;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
protected:
	/** @brief Boundary condition type. */
	int fType{-1};
    /** @brief Pointer to associated TPZMaterial instance*/
    TPZMaterial *fMaterial{nullptr};
};

#endif
