#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "PZMatrixMarket.h"


#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

int main()
{
    {
        TPZSBMatrix<double> fmat;
//        TPZMatrixMarket::Read("../bcsstk15.mtx",fmat);
        TPZMatrixMarket::Read("../Andrade.mtx",fmat);
#ifdef USING_BOOST
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
        
        fmat.Decompose_Cholesky();
        
#ifdef USING_BOOST
        boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
        std::cout << "Time for decompose " << t2-t1 << std::endl;
#endif
    }
    {
        TPZSBMatrix<float> fmat;
//        TPZMatrixMarket::Read("../bcsstk15.mtx",fmat);
        TPZMatrixMarket::Read("../Andrade.mtx",fmat);
#ifdef USING_BOOST
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
        
        fmat.Decompose_Cholesky();
        
#ifdef USING_BOOST
        boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
        std::cout << "Time for decompose " << t2-t1 << std::endl;
#endif
    }
    return 0;
    
}