/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2014  <copyright holder> <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MAT2DPOSPRO_H
#define MAT2DPOSPRO_H
#include "TPZMaterial.h"

class Mat2Dpospro : public TPZMaterial {

protected:
  
    /** @brief Problem dimension */
    int fDim;  

     /** @brief Delta */
    REAL fDelta;   

     /** @brief Alpha */
    REAL fAlpha;  
  
public:

    /** @brief Default Constructor */
    Mat2Dpospro();
    
    /** @brief Constructor with matid and dimension */    
    Mat2Dpospro(int matid, int dim);
    
    /** @brief Destructor */      
    virtual ~Mat2Dpospro();
    
    /** @brief copy constructor */
    Mat2Dpospro(const Mat2Dpospro &other);
    
    /** @brief operator equal */    
    Mat2Dpospro &operator=(const Mat2Dpospro &copy);
    
    /** @brief Print Method */
    virtual void Print(std::ostream & out);
    
    /** @brief Name of the material */
    virtual std::string Name() { return "Mat2Dpospro"; }
    
    /** @brief Returns the integrable dimension */    
    virtual int Dimension() const;
    
    /** @brief Return the number of state variables */
    virtual int NStateVariables(); 
    
    /** @brief Not used contribute methods */
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
      return;
    }
    /** @brief Used contribute methods */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

    /** @brief Fill material data parameter with necessary requirements for the Contribute method*/    
    virtual void FillDataRequirements(TPZMaterialData &data);
    
    /** @brief Fill material data parameter with necessary requirements for the ContributeBC method*/    
    virtual void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data);
    
    /** @brief Returns the variable index associated with the name */    
    virtual int VariableIndex(const std::string &name);
    
    /** @brief Returns the number of variables associated with the variable indexed by var */
    virtual int NSolutionVariables(int var);
    
    /** @brief Calculates a solution given datavec*/    
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    void SetParameters(REAL alpha, REAL delta);

    void GetParameters(REAL &alpha, REAL &delta);   
  
};

#endif // MAT2DPOSPRO_H

