#include "pzfunction.h"

template<class TVar>
class TPZVec;
template<class TVar>
class TPZFMatrix;

/*****************************************************
USAGE EXAMPLE #1


void mySol(const TPZVec<REAL> &loc, TPZVec<STATE> &res,
           TPZFMatrix<STATE>&deriv){
   //fancy stuff
}


//more fancy stuff
const int porder = 4;
TPZAutoPointer<TPZFunction<STATE>> solPtr(new TPZExactFunction<STATE>(mySol,
                                        pOrder));
for(auto imat : cmesh->MaterialVec()){
  imat.second->SetExactSol(solPtr);
}

USAGE EXAMPLE #2

//more fancy stuff
auto mySol = [](const TPZVec<REAL> &loc, TPZVec<STATE> &res,
                TPZFMatrix<STATE>&deriv){
//fancy stuff
};
const int porder = 4;
TPZAutoPointer<TPZFunction<STATE>> solPtr(new TPZExactFunction<STATE>(mySol,
                                        pOrder));
for(auto imat : cmesh->MaterialVec()){
  imat.second->SetExactSol(solPtr);
}
*****************************************************/


//! Class meant for setting exact solutions to TPZMaterial instances
template<class TVar>
class TPZExactFunction : public virtual TPZFunction<TVar> {
private:
    //!Non-managed pointer to actual function
    std::function<void(const TPZVec<REAL> &loc, TPZVec<STATE> &result,
                     TPZFMatrix<STATE> &deriv)> fExact;
    int fPolyOrder = -1;
public:
	
	//! Default constructor.
	TPZExactFunction() = default;
    //! Constructor setting exact function.
    TPZExactFunction(
        std::function<void(const TPZVec<REAL> &loc,
                           TPZVec<STATE> &result,
                           TPZFMatrix<STATE> &deriv)>
        f, int polyOrder) {fExact = f; fPolyOrder = polyOrder;}
    //! Set exact function.
    void SetExact(
        std::function<void(const TPZVec<REAL> &loc,
                           TPZVec<STATE> &result,
                           TPZFMatrix<STATE> &deriv)>
        f, int polyOrder) {fExact = f; fPolyOrder = polyOrder;}

        /**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df) override
    {
        if(fExact){fExact(x,f,df);}
        else{
            PZError<<__PRETTY_FUNCTION__;
            PZError<<" no exact solution was set.\n";
            PZError<<"Aborting...\n";
            DebugStop();
        }
    }

    //! Polynomial order to be used in the integration rule
    int PolynomialOrder() const override{
        return fPolyOrder;
    }

    //! Used when saving to file. Not really useful in this case
    int ClassId() const override;
};


template<class TVar>
int TPZExactFunction<TVar>::ClassId() const{
    return Hash("TPZExactFunction")
        ^ TPZFunction<TVar>::ClassId() << 1;
}
