//$Id: pzcompelpostproc.h,v 1.13 2010-11-24 17:48:10 diogo Exp $

#ifndef PZCOMPELPOSTPROC_H
#define PZCOMPELPOSTPROC_H

class TPZMaterialData;

#include "pzreferredcompel.h"
#include "pzinterpolationspace.h"
#include "tpzautopointer.h"
#include "pzmaterial.h"
#include "pzmaterialdata.h"
#include "pzelmat.h"
#include "pzstack.h"
#include "pzcmesh.h"
#include "pzquad.h" //TPZIntPoints
#include <cmath>
#include "pzlog.h"
#include "pzpostprocmat.h"
#include "pzvec.h"
#include "pzreal.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

#include "TPZShapeDisc.h" 

#ifdef LOG4CXX
static LoggerPtr CompElPostProclogger(Logger::getLogger("mesh.TPZCompElPostProc"));
#endif

using namespace std;

/**
This class implements the TPZCompEl structure to enable copying the solution
of the referred compEl at the integration points to itself and interpolating it
inside the element

@since May, 1 2009
*/




template <class TCOMPEL >

class TPZCompElPostProc : public TPZReferredCompEl<TCOMPEL>
{
public:

  TPZCompElPostProc();

  virtual ~TPZCompElPostProc();
		
  TPZCompElPostProc(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

  TPZCompElPostProc(TPZCompMesh &mesh, const TPZCompElPostProc<TCOMPEL> &copy);

  /**
   * used to generate patch mesh... generates a map of connect index from
   * global mesh to clone mesh
   */
  TPZCompElPostProc(TPZCompMesh &mesh,
              const TPZCompElPostProc<TCOMPEL> &copy,
              std::map<int,int> & gl2lcConMap,
              std::map<int,int> & gl2lcElMap);

  /**
   * Initializes the shape function type in order to allow non ill-conditioned L2 Transfer matrix
   */
  void InitializeShapeFunctions();
	
	
  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const;

  /**
   * Create a copy of the given element. The clone copy have the connect indexes
   * mapped to the local clone connects by the given map
   * @param mesh Patch clone mesh
   * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
   * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
   */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> & gl2lcConMap,std::map<int,int>&gl2lcElMap) const;
	
  void ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi);
  
  /**
   * The CalcResidual reimplementation is in charge of extracting the data from the
   * referred compEl at the integration points and building a minimum square residual
   * method to extrapolate this information throughout the element subdomain.
   * The final ef vector shall be copied onto the solution vector, as it represents
   * the shape functions multipliers of the extrapolation functions.
   */
  virtual void CalcResidual(TPZElementMatrix &ef);
  
  /**
   * Null implementation of the CalcStiff in order to ensure it wouldn't produce
   * any valid system of equations.
   */
  virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * Compare some fields of 2 TPZMaterialData and return true if these do match.
   */
  bool dataequal(TPZMaterialData &d1,TPZMaterialData &d2);

 /**
  * Avoids the calling of the TPZCompElReferred::ComputeSolution wich would attempt to
  * append solution of the referred element.
  * Computes solution and its derivatives in local coordinate qsi
  * @param qsi master element coordinate
  * @param phi matrix containing shape functions compute in qsi point
  * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
  * @param axes direction of the derivatives
  * @param sol finite element solution
  * @param dsol solution derivatives
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                               const TPZFMatrix &axes, TPZSolVec &sol, TPZGradSolVec &dsol);

 /**
  * Avoids the calling of the TPZCompElReferred::ComputeSolution wich would attempt to
  * append solution of the referred element.
  * Computes solution and its derivatives in the local coordinate qsi.
  * This method will function for both volumetric and interface elements
  * @param qsi master element coordinate of the interface element
  * @param sol finite element solution
  * @param dsol solution derivatives
  * @param axes axes associated with the derivative of the solution
  * @param leftsol finite element solution
  * @param dleftsol solution derivatives
  * @param leftaxes axes associated with the left solution
  * @param rightsol finite element solution
  * @param drightsol solution derivatives
  * @param rightaxes axes associated with the right solution
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi,
                             TPZVec<REAL> &normal,
                             TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix &leftaxes,
                             TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix &rightaxes);

	
  /**
   * Reimplemented in order to ensure the shape functions are computed with local support
   * and thus benefits from the orthogonality properties of the legendre polynomials defined
   * by the method InitializeShapeFunctios().
   */
  virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                 TPZFMatrix &jacobian, TPZFMatrix &axes,
                                 REAL &detjac, TPZFMatrix &jacinv,
                                 TPZFMatrix &phi, TPZFMatrix &dphix);
	
public:

  /**
   * Save the element data to a stream
   */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
   * Read the element data from a stream
   */
  virtual void Read(TPZStream &buf, void *context);
	
};

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::TPZCompElPostProc() : TPZReferredCompEl<TCOMPEL>() {
	TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions();
}

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::~TPZCompElPostProc() {

}

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::TPZCompElPostProc(TPZCompMesh &mesh, TPZGeoEl *gel, int &index) :
  TPZReferredCompEl<TCOMPEL>(mesh, gel, index){
	TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions();
}

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::TPZCompElPostProc(TPZCompMesh &mesh, const TPZCompElPostProc<TCOMPEL> &copy) :
  TPZReferredCompEl<TCOMPEL>(mesh, copy) {
	TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions();
}

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::TPZCompElPostProc(TPZCompMesh &mesh,
                                      const TPZCompElPostProc<TCOMPEL> &copy,
                                      std::map<int,int> & gl2lcConMap,
                                      std::map<int,int> & gl2lcElMap):
    TPZReferredCompEl<TCOMPEL>(mesh,copy,gl2lcConMap,gl2lcElMap)
{
	TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions();
}

template <class TCOMPEL>
inline TPZCompEl * TPZCompElPostProc<TCOMPEL>::Clone(TPZCompMesh &mesh) const{
    return new TPZCompElPostProc<TCOMPEL> (mesh, *this);
}

template <class TCOMPEL>
inline TPZCompEl * TPZCompElPostProc<TCOMPEL>::ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> & gl2lcConMap,std::map<int,int>&gl2lcElMap)const{
    return new TPZCompElPostProc<TCOMPEL> (mesh, *this, gl2lcConMap, gl2lcElMap);
}

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions(){
	//TPZReferredCompEl<TCOMPEL>::fShapefunctionType = pzshape::TPZShapeDisc::ETensorial;//pzshape::TPZShapeDisc::EOrdemTotal;
	// an orthogonal (or one closest possible to) polynomial function is very important
	// here to ensure the L2 solution transfer matrix isn't ill-conditioned
	//TPZReferredCompEl<TCOMPEL>::SetOrthogonalFunction(pzshape::TPZShapeDisc::ChebyshevWithoutScale);
 //   TPZReferredCompEl<TCOMPEL>::SetOrthogonalFunction(pzshape::TPZShapeDisc::LegendreWithoutScale);
}

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
  TCOMPEL::ComputeRequiredData(data, qsi); 
}

/**
 * write the element data to a stream
 */
template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::Write(TPZStream &buf, int withclassid)
{
  TCOMPEL::Write(buf,withclassid);
}

/**
 * Read the element data from a stream
 */
template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::Read(TPZStream &buf, void *context)
{
  TCOMPEL::Read(buf,context);
}




template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::CalcResidual(TPZElementMatrix &ef){
  ef.Reset();


	
  this->InitializeElementMatrix(ef);// the inintialization of the ef matrix
	//preceeds the verifications below because it is advisable to esit
	// with proper ef size, no matter the return reason
	
  TPZCompEl * pCompElRef = TPZReferredCompEl<TCOMPEL>::ReferredElement();
	
  TPZInterpolationSpace * pIntSpRef = dynamic_cast<TPZInterpolationSpace *>(pCompElRef);
  if(!pCompElRef){
      PZError << "Error at " << __PRETTY_FUNCTION__ << " Referred CompEl == NULL\n";
      return;
  }
	
  if(pCompElRef == this){
      PZError << "Error at " << __PRETTY_FUNCTION__ << " Referred CompEl == itself\n";
      return;
  }
	
  TPZPostProcMat * pPostProcMat = dynamic_cast<TPZPostProcMat *>(this->Material().operator->());
  if(!pPostProcMat){
      return;
      //skipping unhandled element
  }
  
  TPZMaterial * pMaterialRef = pIntSpRef->Material().operator->();
  if(!pMaterialRef){
      PZError << "Error at " << __PRETTY_FUNCTION__ << " Referred CompEl->Material() == NULL\n";
      return;
  }
	
  if (this->NConnects() == 0) return;///boundary discontinuous elements have this characteristic

  int numeq = ef.fMat.Rows();
  TPZFMatrix efTemp(numeq,1,0.);

  TPZMaterialData data, dataRef;
  this      ->InitMaterialData(data);
  pIntSpRef ->InitMaterialData(dataRef);
  data.p	= this      ->MaxOrder();
  dataRef.p = pIntSpRef ->MaxOrder();
  int dim = pIntSpRef ->Dimension();
  TPZManVector<REAL,3> intpoint(dim,0.);
  TPZManVector<REAL,3> intpointRef(dim,0);
  REAL weight = 0.;
  REAL weightRef = 0;
  const TPZIntPoints &intrule    = this      ->GetIntegrationRule();
  const TPZIntPoints &intruleRef = pIntSpRef ->GetIntegrationRule();
/*
  TPZManVector<int,3> p2(dim,data.p*2);
  TPZManVector<int,3> p2Ref(dim,dataRef.p*2);
  intrule.   SetOrder(p2);
  intruleRef.SetOrder(p2Ref);
  if(pMaterialRef->HasForcingFunction() || !this->Reference()->IsLinearMapping())
 {
      TPZManVector<int,3> order(dim,intrule.GetMaxOrder());
	  TPZManVector<int,3> orderRef(dim,intruleRef.GetMaxOrder());
      intrule   .SetOrder(order);
      intruleRef.SetOrder(orderRef);
  }
 */

  int intrulepoints    = intrule.   NPoints();
  int intrulepointsRef = intruleRef.NPoints();
  if(intrulepoints != intrulepointsRef)
  {
      PZError << "Error at " << __PRETTY_FUNCTION__ << " Referred CompEl with different number of integration points\n";
      return;
  }
	
  int nshape = this->NShapeF();
  TPZFMatrix ekTemp(nshape, nshape, 0.);
  
  TPZVec<int> varIndex;
  int stackedVarSize = pPostProcMat->NStateVariables();
  pPostProcMat->GetPostProcessVarIndexList(varIndex);
  TPZVec<REAL> Sol;
  
  for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
  {
      intrule.   Point(int_ind,intpoint,   weight);
      intruleRef.Point(int_ind,intpointRef,weightRef);
      this->      ComputeShape(intpoint, data.x, data.jacobian, 
                               data.axes, data.detjac, data.jacinv, 
                               data.phi, data.dphix);
//cout << "\n data.phi = " << data.phi;
      pIntSpRef ->ComputeShape(intpointRef, dataRef.x, dataRef.jacobian, 
                               dataRef.axes, dataRef.detjac, dataRef.jacinv, 
                               dataRef.phi, dataRef.dphix); 
//cout << "\n dataRef.phi = " << dataRef.phi;
	  

	  
      weight    *= fabs(data.detjac);
      weightRef *= fabs(dataRef.detjac);
      data   .intPtIndex = int_ind;
      dataRef.intPtIndex = int_ind;
      this      ->ComputeRequiredData(data,    intpoint);
      pIntSpRef ->ComputeRequiredData(dataRef, intpointRef);

//#ifdef LOG4CXX
//	  {
//		  std::stringstream sout;
//		//  sout << "\n data.phi = " << data.phi;
//		  //sout << "\n intpoint = " << data.;
//		  LOGPZ_INFO(CompElPostProclogger,sout.str().c_str());
//	  }
//#endif
	  
	 
//cout << "\tEvaluated data.sol=" << data.sol << "\n";
	  
      if(!dataequal(data,dataRef)){
	      PZError << "Error at " << __PRETTY_FUNCTION__ << " this and Referred CompEl TPZMaterialData(s) do not match\n";
	      ef.Reset();
	      return;
      }
      data.sol[0].Resize(stackedVarSize,0.);
      int index = 0;
      // stacking the solutions to post process.
      for(int var_ind = 0; var_ind < varIndex.NElements(); var_ind++)
      {
		  int nsolvars = pMaterialRef->NSolutionVariables(varIndex[var_ind]);
		  Sol.Resize(nsolvars);
          pMaterialRef->Solution(dataRef, varIndex[var_ind], Sol);
          for(int i = 0; i <nsolvars; i++)data.sol[0][index+i] = Sol[i];
		  index += nsolvars;		
      }
      
//cout << "\timposed data.sol=" << data.sol << "\n";
	  

	  
	  pPostProcMat->Contribute(data,weight,ekTemp,efTemp);
	  
  }//loop over integration points
	
  TPZFMatrix ekCopy(ekTemp);

	//cout << "\nekTemp = \n " << ekTemp;
	//cout << endl << endl;

  TPZFMatrix rhsTemp(nshape, 1, 0.);
  for(int i_st = 0; i_st < stackedVarSize; i_st++)
  {
	  
		efTemp.GetSub(i_st*nshape, 0, nshape, 1, rhsTemp);
	  
	    TPZFMatrix rhsCopy(rhsTemp), result;
  	//	int status = ekTemp.Solve_Cholesky(&(rhsTemp));
	    int status = ekTemp.Solve_LU(&(rhsTemp));
	  
	    ekCopy.MultAdd(rhsTemp, rhsCopy, result, 1., -1.);
	    REAL invRes = Norm(result);
 		if(!status ){
	  		PZError << "Error at " << __PRETTY_FUNCTION__ << " Unable to solve the transference linear system\n";
	  		ef.Reset();
	  		return;
		}
 		
	    if(invRes > 1.e-8)PZError << "Error at " << __PRETTY_FUNCTION__ 
		  << " Transference linear system solved with residual norm = " 
		  << invRes << " at " << i_st << " export variable\n";
	  
	    for(int i_sh = 0; i_sh < nshape; i_sh++)
		  	ef.fMat(i_sh * stackedVarSize + i_st, 0) = rhsTemp(i_sh); 
  }
	
	cout << "*";
	cout.flush();
	
//cout << "\n solved solution:\n" << ef.fMat << endl;

}//CalcResidual

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){
	PZError << "\nTPZCompElPostProc<TCOMPEL>::CalcStiff() Should never be called!!!\n";
	return;
}

template <class TCOMPEL>
inline bool TPZCompElPostProc<TCOMPEL>::dataequal(TPZMaterialData &d1,TPZMaterialData &d2)
{
	const REAL SMALLNUMBER = 1.e-8;
	int i;
	if(d1.p!=d2.p)
	{
		DebugStop();
		return 0;
	}
	REAL res = 0;
	int dim = d1.x.NElements();
	int nshape = d1.phi.Rows();
	int nshape2 = d2.phi.Rows();
	if(dim != d2.x.NElements() || nshape!= nshape2) 
	{
		DebugStop();
		return 0; // dimensions and number of integration points shall match
	}
	for(i = 0; i < dim; i++)res += pow(d1.x[i]-d2.x[i],2.); // integration points must be at the same locations
	/*for(i = 0; i < nshape; i++)
	{
		int j;
	    res += pow(d1.phi(i,0)-d2.phi(i,0),2.);
	    for(j = 0; j < dim; j++) res += pow(d1.dphix(j,i)-d2.dphix(j,i),2.);
	}
	res += pow(d1.detjac - d2.detjac, 2.);*/
    
//res += sqr(d1.weight - d2.weight);
    if(res > SMALLNUMBER)
	{
		DebugStop();
		return 0;
	}
    return 1;
}

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::ComputeSolution(TPZVec<REAL> &qsi,
                                                   TPZFMatrix &phi,
                                                   TPZFMatrix &dphix,
                                                   const TPZFMatrix &axes,
                                                   TPZSolVec &sol,
                                                   TPZGradSolVec &dsol){
  TCOMPEL::ComputeSolution(qsi, phi, dphix, axes, sol, dsol);
}//method

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::ComputeSolution(TPZVec<REAL> &qsi,
                                                   TPZVec<REAL> &normal,
                                                   TPZSolVec &leftsol, TPZGradSolVec &dleftsol, TPZFMatrix &leftaxes,
                                                   TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix &rightaxes){
  TCOMPEL::ComputeSolution(qsi, normal, leftsol, dleftsol, leftaxes, rightsol, drightsol, rightaxes);
}

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                 TPZFMatrix &jacobian, TPZFMatrix &axes,
                                 REAL &detjac, TPZFMatrix &jacinv,
                                 TPZFMatrix &phi, TPZFMatrix &dphix){
  TPZGeoEl * ref = this->Reference();
  if (!ref){
    PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
    return;
  }//if
  ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);

  ref->X(intpoint, X);
//  this->Shape(intpoint,intpoint,phi,dphix);
	this->Shape(intpoint,phi,dphix);
  //this->Shape(intpoint,X,phi,dphix);

  ///axes is identity in discontinuous elements
  axes.Resize(dphix.Rows(), dphix.Rows());
  axes.Identity();
}

#endif
