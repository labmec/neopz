/**
 * @file
 * @brief Contains the implementation of the TPZFStructMatrix methods. 
 */

#include "pzfstrmatrix.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include <sstream>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.tpzfstructmatrix"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
#endif


using namespace std;

TPZMatrix<STATE> * TPZFStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix<STATE> *stiff = Create();
	long neq = this->fMesh->NEquations();
	rhs.Redim(neq,1);
	Assemble(*stiff,rhs,guiInterface);
	
#ifdef LOG4CXX
	if(loggerel->isDebugEnabled())
	{
		std::stringstream sout;
		stiff->Print("Stiffness matrix",sout);
		rhs.Print("Right hand side", sout);
		LOGPZ_DEBUG(loggerel,sout.str())
	}
#endif
    return stiff;
}

TPZMatrix<STATE> * TPZFStructMatrix::Create(){
	long neq = fEquationFilter.NActiveEquations();
    
	return new TPZFMatrix<STATE>(neq,neq,0.);
}

TPZFStructMatrix::TPZFStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

TPZFStructMatrix::TPZFStructMatrix(TPZAutoPointer<TPZCompMesh> mesh) : TPZStructMatrix(mesh)
{
}

TPZStructMatrix * TPZFStructMatrix::Clone(){
    return new TPZFStructMatrix(*this);
}
