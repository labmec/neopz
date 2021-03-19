#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzsysmp.h"
#include "pzlog.h"

#include "pzstepsolver.h"

#include "PZMatrixMarket.h"


#include <cmath>
#include <set>

#ifdef PZ_LOG
static PZLogger logger("pz.elasticity");
#endif

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

template<class TFloat, class TDouble>
void SolvePreconditioned(TDouble *objdouble, TFloat *objfloat, DecomposeType EDecompose);

int main()
{
    
    DecomposeType dec = ECholesky;
    
    
    
    {
//        TPZSYsmpMatrix<double> fmat;
        TPZSkylMatrix<double> fmat;
//        TPZSBMatrix<double> fmat;
//        TPZMatrixMarket::Read("../bcsstk15.mtx",fmat);
        TPZMatrixMarket::Read("../Andrade.mtx",fmat);
#ifdef USING_BOOST
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
        
        fmat.Decompose_Cholesky();
        
#ifdef USING_BOOST
        boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
        std::cout << "Time for decompose double " << t2-t1 << std::endl;
        std::cout.flush();
#endif
    }
    {
//        TPZSYsmpMatrix<float> fmat;
        TPZSkylMatrix<float> fmat;
//        TPZSBMatrix<float> fmat;
//        TPZMatrixMarket::Read("../bcsstk15.mtx",fmat);
        TPZMatrixMarket::Read("../Andrade.mtx",fmat);
#ifdef USING_BOOST
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
        
        fmat.Decompose_Cholesky();
        
#ifdef USING_BOOST
        boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
        std::cout << "Time for decompose float " << t2-t1 << std::endl;
        std::cout.flush();
#endif
    }
    
//    DecomposeType dec = ECholesky;
    
    {
        TPZSkylMatrix<double> *fmatDouble = new TPZSkylMatrix<double>();
        TPZSkylMatrix<float> *fmatFloat = new TPZSkylMatrix<float>();
        TPZMatrixMarket::Read("../Andrade.mtx",*fmatDouble);
        fmatFloat->CopyFrom(*fmatDouble);
        
#ifdef USING_BOOST
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
        
        SolvePreconditioned(fmatDouble, fmatFloat, dec);
        
#ifdef USING_BOOST
        boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
        std::cout << "Time for preconditioned solution " << t2-t1 << std::endl;
        std::cout.flush();
#endif
        
    }
    
    {
        TPZSkylMatrix<double> *fmatDouble = new TPZSkylMatrix<double>();
        TPZMatrixMarket::Read("../Andrade.mtx",*fmatDouble);
        TPZAutoPointer<TPZMatrix<double> > matrix = fmatDouble;
        
#ifdef USING_BOOST
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
        
        TPZStepSolver<double> step(matrix);
        step.SetDirect(dec);
        TPZFMatrix<double> rhs(fmatDouble->Rows(),1,1.), sol(fmatDouble->Rows(),1.,0.);
        step.Solve(rhs, sol);
        
#ifdef USING_BOOST
        boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
        std::cout << "Time for direct solution " << t2-t1 << std::endl;
        std::cout.flush();
#endif
    }
    
    
    return 0;
    
}

template<class TFloat, class TDouble>
void SolvePreconditioned(TDouble *objdouble, TFloat *objfloat, DecomposeType EDecompose)
{
    TPZAutoPointer<TPZMatrix<double> > matrix = objdouble;
    TPZAutoPointer<TPZMatrix<float> > floatmat = objfloat;
    TPZStepSolver<double> step(matrix);
    
    objfloat->CopyFrom(*objdouble);
    std::list<int64_t> singular;
    objfloat->Decompose(EDecompose,singular);
    TDouble *floatdec = new TDouble();
    floatdec->CopyFrom(*objfloat);
    TPZStepSolver<double> stepfloat;
    stepfloat.SetMatrix(floatdec);
    stepfloat.SetDirect(EDecompose);
    step.SetCG(10, stepfloat, 1.e-6, 0);
    
    TPZFMatrix<double> rhs(objdouble->Rows(),1,1.), sol(objdouble->Rows(),1.,0.);
    step.Solve(rhs, sol);

}

