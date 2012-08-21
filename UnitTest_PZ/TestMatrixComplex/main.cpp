/* 
 * File:   main.cpp
 * Author: adriano
 *
 * Created on July 25, 2012, 1:03 AM
 */

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzblockdiag.h"
#include "pzbndmat.h"
#include "pzespmat.h"
#include "pzsbndmat.h"
#include "pzsfulmat.h"


using namespace std;

/*
 * 
 */

template<class T>
class ComplexMatrixTester {

public:
    ComplexMatrixTester() {
    
        rightMatrix.AutoFill();
        leftMatrix.AutoFill();
        resultMatrix.AutoFill();
    }
    
    ComplexMatrixTester(T left, T right) {
        
        SetLeftMatrix(left);
        SetRightMatrix(right);  
        //resultMatrix = new T;
    }
    
    ComplexMatrixTester(const ComplexMatrixTester &orig) {}
    
    virtual ~ComplexMatrixTester() {}
    
    void SetResultMatrix(T resultMatrix) {
    
        this->resultMatrix = resultMatrix;
    }
    
    T GetResultMatrix() const {
    
        return resultMatrix;
    }
    void SetRightMatrix(T rightMatrix) {
    
        this->rightMatrix = rightMatrix;
    }
    
    T GetRightMatrix() const {
    
        return rightMatrix;
    }
    
    void SetLeftMatrix(T leftMatrix) {
    
        this->leftMatrix = leftMatrix;
    }
    
    T GetLeftMatrix() const {
    
        return leftMatrix;
    }
    
    void Add() {
    
         leftMatrix.Add(rightMatrix, resultMatrix);
    }
    
    void Substract() {
    
        leftMatrix.Substract(rightMatrix, resultMatrix);
    }
    
    void Multiply() {
    
        leftMatrix.Multiply(rightMatrix, resultMatrix, 0, 1);
    }
    
    void DecomposeLU() {
    
        leftMatrix.Decompose_LU();
    }
    
    void PrintMatlabFormat(const char *matlabFileName, std::ostream &out, T matrix) const {
    
        out << matlabFileName << " = [";
    
        for (int x = 0; x < matrix.Rows(); x++) {
            for (int y = 0; y < matrix.Cols(); y++) {
                out << matrix(x, y).real() << " + i*" << matrix(x, y).imag() << " ";
            }
        
            if (x != matrix.Rows() - 1)
                out << ";";  
        }
   
        out<<"]; \n";
        //matrix.Print(matlabFileName, out, EMathematicaInput);
    }

private:
    T leftMatrix;
    T rightMatrix;
    T resultMatrix;
};

int main() {
    
    // dimension: square matrix
    const int dim = 3;
    
    // full matrix
    TPZFNMatrix<dim, STATE> leftMatrix(dim, dim, 0.0);
    TPZFNMatrix<dim, STATE> rightMatrix(dim, dim, 0.0);
    // fill matrices with random values
    leftMatrix.AutoFill();
    rightMatrix.AutoFill();
    
    // band matrix
    //TPZFBMatrix<STATE> leftMatrix(dim, dim);
    //TPZFBMatrix<STATE> rightMatrix(dim, dim);
    //leftMatrix.AutoFill();
    //rightMatrix.AutoFill();
    
    // sparse matrix
    //TPZSpMatrix<STATE> leftMatrix;
    //TPZSpMatrix<STATE> rightMatrix;
    //leftMatrix.AutoFill();
    //rightMatrix.AutoFill();
    
    // Full Matrix testing
    ComplexMatrixTester<TPZFNMatrix<dim, STATE> > operationTester(leftMatrix, rightMatrix);
    // Band Matrix testing
    //ComplexMatrixTester<TPZFBMatrix<STATE> > operationTester(leftMatrix, rightMatrix);
    //ComplexMatrixTester<TPZFBMatrix<STATE> > operationTester;
    // Sparse Matrix testing
    //ComplexMatrixTester<TPZSpMatrix<STATE> > operationTester(leftMatrix, rightMatrix);
    
    // printing left and right matrices original values
    ofstream outLeftMatrix("ComplexLeftMatrix.dat");
    ofstream outRightMatrix("ComplexRightMatrix.dat");
    operationTester.PrintMatlabFormat("leftMatrix", outLeftMatrix, operationTester.GetLeftMatrix());
    operationTester.PrintMatlabFormat("rightMatrix", outRightMatrix, operationTester.GetRightMatrix());
    
    // Add operation
    operationTester.Add();
    ofstream outSumMatrix("SumResult.dat");
    operationTester.PrintMatlabFormat("sumMatrix", outSumMatrix, operationTester.GetResultMatrix());
    
     // Sub operation
    operationTester.Substract();
    ofstream outSubMatrix("SubResult.dat");
    operationTester.PrintMatlabFormat("subMatrix", outSubMatrix, operationTester.GetResultMatrix());
    
     // Multiplication operation
    operationTester.Multiply();
    ofstream outMultiMatrix("MultiplicationResult.dat");
    operationTester.PrintMatlabFormat("multiMatrix", outMultiMatrix, operationTester.GetResultMatrix());
    
    // Lu decompose operation
    // left Matrix
    operationTester.DecomposeLU();
    ofstream outLUMatrix("LUDecompeResult.dat");
    operationTester.PrintMatlabFormat("LUMatrix", outLUMatrix, operationTester.GetLeftMatrix());

   return 0;
}








