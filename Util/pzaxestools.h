/**
 * @file
 * @brief Contains declaration of the TPZAxesTools class which implements verifications over axes.
 */

#ifndef AXESTOOLS
#define AXESTOOLS

#include "pzfmatrix.h"
#include "pzerror.h"

/**
 * @ingroup util
 * @brief Implements method to verify whether axes is an orthogonal normalizes matrix and to transformation from given axes
 * to euclidian basis and viceverse. \ref util "Utility"
 */
template<class TVar>
class TPZAxesTools {
public:
    /// simple constructor
    TPZAxesTools() {}
    /// destructor
    ~TPZAxesTools() {}
    
    /**
     * @brief Verify whether parameter axes is an orthogonal normalized matrix.
     * @param axes Object to check if it is a orthogonal and normalized matrix
     */
    static void VerifyAxes(const TPZFMatrix<TVar> &axes) {
#ifdef PZDEBUG
//        const double tol = 1.e-8;
        bool check = true;
        for(int i = 0; i < axes.Rows(); i++){
            for(int j = 0; j < axes.Rows(); j++) {
                if(i == j) continue;
                TVar sum = 0.;
                for(int k = 0; k < axes.Cols(); k++){
                    sum += axes.GetVal(i,k)*axes.GetVal(j,k);
                }//k
                if(!IsZero(sum)){
                    check = false;
                }//if
            }//for j
        }//for i
        
        if(check == false){
            PZError << "\nError at " << __PRETTY_FUNCTION__;
            axes.Print(", axes = ", PZError);
            PZError <<  "\n";
        }
#endif
    }
    
    /**
     * @brief Makes the basis transformation from axes basis to euclidian basis.
     * @param dudx Output matrix
     * @param dudaxes Input matrix
     * @param axesv Must be an orthogonal normalized matrix. Axes vectors are written in rows.
     * @param colMajor If true, dudaxes is stored as a column-major matrix. Else, row-major.
     */
    static void Axes2XYZ(const TPZFMatrix<TVar> &dudaxes, TPZFMatrix<TVar> &dudx, const TPZFMatrix<REAL> &axesv, bool colMajor = true){
        TPZFNMatrix<9,TVar> axes(axesv.Rows(),axesv.Cols());
        for (int r=0; r<axes.Rows(); r++) {
            for (int c=0; c<axes.Cols(); c++) {
                axes(r,c) = axesv.GetVal(r,c);
            }
        }
#ifdef PZDEBUG
        TPZAxesTools<REAL>::VerifyAxes(axesv);
#endif
        if( colMajor ){
            TPZFNMatrix<9,TVar> axesT;
            axes.Transpose(&axesT);
            if(dudx.Rows() != axesT.Rows() || dudx.Cols() != dudaxes.Cols())
            {
                dudx.Redim(axesT.Rows(), dudaxes.Cols());
            }
            dudx.Zero();
            axesT.Multiply(dudaxes,dudx);
        }
        else{
            dudx.Resize(dudaxes.Rows(), axes.Cols());
            dudx.Zero();
            dudaxes.Multiply(axes,dudx);
        }
    }
    
    /**
     * @brief Makes the basis transformation from euclidian basis to axes basis
     * @param dudaxes Input matrix
     * @param dudx Output matrix
     * @param axes must be an orthogonal normalized matrix. Axes vectors are written in rows.
     * @param colMajor If true, dudaxes is stored as a column-major matrix. Else, row-major.
     */
    static void XYZ2Axes(TPZFMatrix<TVar> &dudaxes, const TPZFMatrix<TVar> &dudx, const TPZFMatrix<REAL> &axes, bool colMajor = true){
        TPZAxesTools::VerifyAxes(axes);
        if( colMajor ){
            dudaxes.Resize(axes.Rows(), dudx.Cols());
            dudaxes.Zero();
            axes.Multiply(dudx,dudaxes);
        }
        else{
            TPZFNMatrix<9,TVar> axesT;
            axes.Transpose(&axesT);
            dudaxes.Resize(dudx.Rows(), axesT.Cols());
            dudaxes.Zero();
            dudx.Multiply(axesT,dudaxes);
        }
    }
    
    static TVar ComputeDetjac(TPZFMatrix<TVar> &gradx){
        
        int dim = gradx.Cols();
        TVar detjac =0.;
        
        if(dim==1){
            detjac = sqrt(gradx(0,0)*gradx(0,0)+gradx(1,0)*gradx(1,0)+gradx(2,0)*gradx(2,0));
        }else if(dim==2){
            TVar vec[3] = {gradx(1,0)*gradx(2,1)-gradx(2,0)*gradx(1,1),
                gradx(2,0)*gradx(0,1)-gradx(0,0)*gradx(2,1),gradx(0,0)*gradx(1,1) - gradx(0,1)*gradx(1,0)};
            detjac = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
            
        }else if(dim==3){
            TVar val = gradx(0,0)*gradx(1,1)*gradx(2,2) + gradx(0,1)*gradx(1,2)*gradx(2,0) + gradx(0,2)*gradx(1,0)*gradx(2,1) - gradx(0,2)*gradx(1,1)*gradx(2,0) - gradx(0,0)*gradx(1,2)*gradx(2,1) - gradx(0,1)*gradx(1,0)*gradx(2,2);
            detjac = sqrt(val*val);
        }else{
            DebugStop();
        }
        return detjac;
    }
    
    
    /** @brief Compute GradX as a function of jac and axes */
    static void ComputeGradX(TPZFMatrix<TVar> &jac, TPZFMatrix<TVar> &axes, TPZFMatrix<TVar> &gradx)
    {
        int nc = jac.Rows();
        gradx.Redim(3, jac.Rows());
        for (int i=0; i<3; i++) {
            for (int j=0; j<nc; j++) {
                for (int k=0; k<nc; k++) {
                    gradx(i,j) += axes(k,i)*jac(k,j);
                }
            }
        }
    }
    
    /** Test code */
    static int main(){
        TPZFMatrix<> axes(2,3,0.);
        axes(0,0) = sqrt(2.)/2.;
        axes(0,1) = 0.;
        axes(0,2) = sqrt(2.)/2.;
        axes(1,0) = 0.;
        axes(1,1) = 1.;
        axes(1,2) = 0.;
        TPZFMatrix<> dudaxes(2,1);
        dudaxes(0,0) = 0.3;
        dudaxes(1,0) = 0.7;
        TPZFMatrix<> dudx;
        dudaxes.Print("dudaxes=",std::cout);
        TPZAxesTools::Axes2XYZ(dudaxes,dudx,axes);
        dudx.Print("dudx=",std::cout);
        TPZAxesTools::XYZ2Axes(dudaxes,dudx,axes);
        dudaxes.Print("dudaxes=",std::cout);
        return 1;
    }
    
};

#endif

