/**
 * \file
 * @brief Contains implementations of the TPZMaterialTest3D methods.
 */
#include "pzmattest3d.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzadmchunk.h"
#include <math.h>
#include <cmath>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.tpzmattest3d"));
#endif
using namespace std;
int TPZMaterialTest3D::geq3 = -1;//Cedric : para testes no programa main 3dmaterial.c

TPZMaterialTest3D::TPZMaterialTest3D() : TPZMaterial(), fXf()
{
}

TPZMaterialTest3D::TPZMaterialTest3D(int nummat) : TPZMaterial(nummat), fXf()
{
}

TPZMaterialTest3D::~TPZMaterialTest3D()
{
}

int TPZMaterialTest3D::NStateVariables()
{
	return 1;
}

void TPZMaterialTest3D::Print(std::ostream &out)
{
	std::stringstream sout;
	sout << " Name of material : " << Name() << "properties: " ;
	fXf.Print( sout.str().c_str(),out,EFormatted );
	sout << " Common properties: " ;
	TPZMaterial::Print(out);
}

void TPZMaterialTest3D::Contribute( TPZMaterialData &data,REAL weight,
								   TPZFMatrix &ek,TPZFMatrix &ef )
{
	TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	int phr = phi.Rows();
	if(fForcingFunction) {            // phi(in, 0) = phi_in
		TPZManVector<REAL> res(1);
		fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
		fXf(0,0) = res[0];
	}
	if(geq3==0) {
		//proje�o L2 da carga fXf : fForcingFunction
		/*for( int jn = 0; jn < phr; jn++ ) {
		 ef(jn, 0) += weight * fXf(0,0) * phi(jn, 0);
		 for( int in = 0; in < phr; in++ ) {
		 ek(in,jn) += weight * phi(in,0) * phi(jn,0);
		 }
		 }*/
		for( int in = 0; in < phr; in++ )
		{
			ef(in, 0) += weight * fXf(0,0) * phi(in, 0);
			for( int jn = 0; jn < phr; jn++ )
			{
				ek(in,jn) += weight * phi(in,0) * phi(jn,0);
			}
		}
	}else if( geq3==1 )
	{
		//Equa�o de Poisson
		for( int in = 0; in < phr; in++ )
		{
			ef(in, 0) += weight * fXf(0,0) * phi(in, 0);
			for( int jn = 0; jn < phr; jn++ )
			{
				/*REAL dphix = axes(0,0)*dphi(0,in)+axes(1,0)*dphi(1,in)+axes(2,0)*dphi(2,in);
				 REAL dphiy = axes(0,1)*dphi(0,in)+axes(1,1)*dphi(1,in)+axes(2,1)*dphi(2,in);
				 REAL dphiz = axes(0,2)*dphi(0,in)+axes(1,2)*dphi(1,in)+axes(2,2)*dphi(2,in);
				 REAL dphjx = axes(0,0)*dphi(0,jn)+axes(1,0)*dphi(1,jn)+axes(2,0)*dphi(2,jn);
				 REAL dphjy = axes(0,1)*dphi(0,jn)+axes(1,1)*dphi(1,jn)+axes(2,1)*dphi(2,jn);
				 REAL dphjz = axes(0,2)*dphi(0,jn)+axes(1,2)*dphi(1,jn)+axes(2,2)*dphi(2,jn);
				 ek(in,jn) += weight*(dphix*dphjx+dphiy*dphjy+dphiz*dphjz);*/
				ek(in,jn) += weight * ( dphi(0,in) * dphi(0,jn) + dphi(1,in) * dphi(1,jn) );//
				if(dphi.Rows() == 3) ek(in,jn) += weight *  dphi(2,in) * dphi(2,jn);// );
			}
		}
	}
}
/*  //2D
 for(int jn=0 ; jn<phi.Rows() ; ++jn){
 for(ic=0; ic<r; ic++) for(jc=0; jc<r; jc++) {
 REAL dphix = dphi(0,in)*axes(0,0)+axes(1,0)*dphi(1,in);
 REAL dphiy = dphi(0,in)*axes(0,1)+axes(1,1)*dphi(1,in);
 REAL dphjx = dphi(0,jn)*axes(0,0)+axes(1,0)*dphi(1,jn);
 REAL dphjy = dphi(0,jn)*axes(0,1)+axes(1,1)*dphi(1,jn);
 ek(in*r+ic,jn*r+jc) += weight*(dphix*dphjx+dphiy*dphjy);
 }
 }
 */
void TPZMaterialTest3D::ContributeBC( TPZMaterialData &data,REAL weight,
									 TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc)
{
	// TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	int phr = phi.Rows();
	short in,jn;
	REAL v2[1];
	v2[0] = bc.Val2()(0,0);
	
	switch (bc.Type())
	{
		case 0 : 			// Dirichlet condition
		{
			for(in = 0 ; in < phr; in++)
			{
				ef(in,0) += gBigNumber * v2[0] * phi(in,0) * weight;
				for (jn = 0 ; jn < phr; jn++)
				{
					ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
				}
			}
			break;
		}
		case 1 :			// Neumann condition
		{
			for(in = 0 ; in < phi.Rows(); in++)
			{
				ef(in,0) += v2[0] * phi(in,0) * weight;
			}
			break;
		}
		case 2 :		// condicao mista
		{
			for(in = 0 ; in < phi.Rows(); in++)
			{
				ef(in, 0) += v2[0] * phi(in, 0) * weight;
				for (jn = 0 ; jn < phi.Rows(); jn++)
				{
					// peso de contorno => integral de contorno
					ek(in,jn) += bc.Val1()(0,0) * phi(in,0) * phi(jn,0) * weight;
				}
			}
		}
	}
}


/** returns the variable index associated with the name*/
int TPZMaterialTest3D::VariableIndex(const std::string &name)
{
	if(!strcmp("Displacement6",name.c_str()))return  0;
	if(!strcmp("Solution",name.c_str()))     return  1;
	if(!strcmp("Derivate",name.c_str()))     return  2;
	return TPZMaterial::VariableIndex(name);
}


int TPZMaterialTest3D::NSolutionVariables(int var)
{
	if(var == 0 || var == 1 || var == 2 || var == 10) return 1;
	return TPZMaterial::NSolutionVariables(var);
}


void TPZMaterialTest3D::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,
                                 TPZFMatrix &axes,int var,TPZVec<REAL> &Solout)
{
	if(var == 0 || var == 1) Solout[0] = Sol[0];//function
	else if(var == 2)
	{
		Solout[0] = DSol(0,0);//derivate
		Solout[1] = DSol(1,0);//derivate
		if(DSol.Rows()>2) Solout[2] = DSol(2,0);//derivate
	}
	else TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}


void TPZMaterialTest3D::Flux( TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/,
							 TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/)
{
	//Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
	LOGPZ_WARN( logger,"ERROR - Not Implemented yet!");
}


void TPZMaterialTest3D::Errors( TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,TPZFMatrix &dudx,
							   TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,TPZVec<REAL> & u_exact,
							   TPZFMatrix &du_exact,TPZVec<REAL> &values)
{
	//TPZVec<REAL> sol(1),dsol(3);
	TPZManVector<REAL> sol(1),dsol(3);
	Solution(u,dudx,axes,1,sol);
	Solution(u,dudx,axes,2,dsol);
	if(dudx.Rows()<3)
	{
		REAL dx = du_exact(0,0)*axes(0,0)+du_exact(1,0)*axes(0,1);
		REAL dy = du_exact(0,0)*axes(1,0)+du_exact(1,0)*axes(1,1);
		REAL parc1 = fabs(dx-dudx(0,0));
		REAL parc2 = fabs(dy-dudx(1,0));
		//Norma L2
		values[1] = pow(fabs(u[0] - u_exact[0]),(REAL)2.0);
		//seminorma
		values[2] = pow(parc1,(REAL)2.)+pow(parc2,(REAL)2.);
		//Norma Energia
		values[0] = values[1]+values[2];
		return;
	}
	//values[1] : eror em norma L2
	values[1]  = pow(sol[0] - u_exact[0],(REAL)2.0);
	//values[2] : erro em semi norma H1
	values[2]  = pow(dsol[0] - du_exact(0,0),(REAL)2.0);
	if(dudx.Rows()>1) values[2] += pow(dsol[1] - du_exact(1,0),(REAL)2.0);
	if(dudx.Rows()>2) values[2] += pow(dsol[2] - du_exact(2,0),(REAL)2.0);
	//values[0] : erro em norma H1 <=> norma Energia
	values[0]  = values[1]+values[2];
}

/*
 REAL dphix = axes(0,0)*dsol[0]+axes(1,0)*dsol[1]+axes(2,0)*dsol[2];
 REAL dphiy = axes(0,1)*dsol[0]+axes(1,1)*dsol[1]+axes(2,1)*dsol[2];
 REAL dphiz = axes(0,2)*dsol[0]+axes(1,2)*dsol[1]+axes(2,2)*dsol[2];
 
 //values[2]  = pow(dphix - du_exact(0,0),2.0);
 //values[2] += pow(dphiy - du_exact(1,0),2.0);
 //values[2] += pow(dphiz - du_exact(2,0),2.0);
 if(dudx.Rows()<3) {
 REAL dx = du_exact(0,0)*axes(0,0);//+du_exact(1,0)*axes(0,1);
 REAL dy;
 if(dudx.Rows()>1)
 dx += du_exact(1,0)*axes(0,1);
 dy = du_exact(0,0)*axes(1,0)+du_exact(1,0)*axes(1,1);
 //values[1] : eror em norma L2
 values[1]  = pow(sol[0] - u_exact[0],2.0);
 //values[2] : erro em semi norma H1
 values[2]  = pow(dx - dudx(0,0),2.0);
 if(dudx.Rows()>1) values[2] += pow(dy - dudx(1,0),2.0);
 //values[0] : erro em norma H1 <=> norma Energia
 values[0]  = values[1]+values[2];
 return;
 }
 */

TPZAutoPointer<TPZMaterial>  TPZMaterialTest3D::NewMaterial()
{
	int matid = Id();
	TPZMaterialTest3D *mat = new TPZMaterialTest3D(matid);
	mat->fXf = fXf;
	return mat;
}


void TPZMaterialTest3D::SetMaterial(TPZFMatrix &xfin)
{
	fXf = xfin;
}


int TPZMaterialTest3D::Dimension()
{
	return 3;
}


void TPZMaterialTest3D::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
	fXf.Read( buf,0 );
	buf.Read( &geq3, 1);
#ifdef DEBUG2
	int classid = -1;
	buf.Read( &classid,1 );
	if( classid != ClassId() )
	{
		std::stringstream sout;
		sout << "Error restoring object " << __PRETTY_FUNCTION__
		<< " waiting for ClassId()= " << ClassId()
		<< " restored : " << classid;
		LOGPZ_ERROR( logger,sout.str().c_str() );
	}
#endif
}


void TPZMaterialTest3D::Write(TPZStream &buf, int withclassid)
{
	TPZMaterial::Write(buf,withclassid);
	fXf.Write( buf,0);
	buf.Write( &geq3,1 );
#ifdef DEBUG2
	int classid = ClassId();
	buf.Write( &classid,1 );
#endif
}

int TPZMaterialTest3D::ClassId() const
{
	return TPZMATTEST3DID;
}

template class
TPZRestoreClass < TPZMaterialTest3D,TPZMATTEST3DID > ;

