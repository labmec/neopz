/**
 * \file
 * @brief Contains implementations of the TPZNLMat1d methods.
 */
//$Id: pznlmat1d.cpp,v 1.7 2008-10-20 11:56:21 longhin Exp $

#include "pznlmat1d.h"
#include "pzbndcond.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.tpznlmat1d"));
#endif

TPZNLMat1d::TPZNLMat1d(int id)
: TPZMaterial(id)
{
	fArea = 0.;
	fE = 0.;
}


TPZNLMat1d::~TPZNLMat1d()
{
	//Nothing to be done here!
}

void TPZNLMat1d::Print(std::ostream &out)
{
	out << __PRETTY_FUNCTION__ << std::endl;
	out << "Cross section initial area = " << fArea << std::endl;
	out << "Young module = " << fE << std::endl;
}

int TPZNLMat1d::VariableIndex(const std::string &name)
{
	if (!strcmp(name.c_str(),"Tens�")) return 0;
	else return -1;
}

int TPZNLMat1d::NSolutionVariables(int var)
{
	return (var==0) ? 1 : -1;
}

int TPZNLMat1d::ClassId() const
{
	return TPZNLMAT1D;
}

void TPZNLMat1d::SetData(std::istream &data)
{
	data >> fArea >> fE;
}

void TPZNLMat1d::Write(TPZStream &buf, int withclassid)
{
	TPZMaterial::Write(buf,withclassid);
	buf.Write(&fArea);
	buf.Write(&fE);
	int classid = ClassId();
	buf.Write(&classid);
}

void TPZNLMat1d::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
	buf.Read(&fArea);
	buf.Read(&fE);
	int classid;
	buf.Read (&classid);
	
	if (classid != ClassId())
	{
		LOGPZ_ERROR (logger, " wrong classid");
	}
}

void TPZNLMat1d::Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef)
{
	TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	// TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	TPZVec<REAL> &sol=data.sol[0];
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	TPZFMatrix &axes=data.axes;
	
	
	//   //sol -> (u1,w1,u2,w2)
	//   //dphi = 2/l => l0 = 2/dphi
	//   //x21 = delta x = l cos theta
	//   //z21 = delta y = l sen theta
	//   //theta = arctg (axes (0,1) / axes (0,0) )
	double theta = atan(axes (0,1) / axes (0,0) );
	double l = 2./ dphi(0,0);
	double x21 = l * cos(theta);
	double z21 = l * sin(theta);
	double u1,u2,w1,w2;
	u1 = sol [0];
	u2 = sol [2];
	w1 = sol [1];
	w2 = sol [3];
	double u21 = u2 - u1;
	double w21 = w2 - w1;
	double alpha0 = l/2.;
	
	double eps = Eps(sol,axes,dphi);
	double sigma = eps * fE;
	
	int i,j;
	double fqi = 0.;
	fqi += 2.*alpha0*fArea*sigma;
	
	//b1
	TPZVec<REAL> b1 (4,0.), b2(4,0.), b(4,0.);
	b1[0] = -x21;  b2[0] = -u21;
	b1[1] =  x21;  b2[1] =  u21;
	b1[2] = -z21;  b2[2] = -w21;
	b1[3] =  z21;  b2[3] =  w21;
	
	double fac = (1. / 4.) / (alpha0 * alpha0);
	for (i=0;i<4;i++)
	{
		b1[i] *= fac;
		b2[i] *= fac;
		b [i] = b1[i] + b2[i];
		ef(i,0) = fqi * b[i] / l;
	}
	
	TPZFMatrix Kt1(4,4,0.);
	TPZFMatrix Kt2(4,4,0.);
	TPZFMatrix Kts(4,4,0.);
	
	fac = 2. * fE * fArea * alpha0;
	Kts(1,0) = -1.;
	Kts(0,1) = -1.;
	Kts(2,3) = -1.;
	Kts(3,2) = -1.;
	for (i=0;i<4;i++)
	{
		Kts(i,i) = 1.;
		for (j=0;j<4;j++)
		{
			Kt1(i,j) = fac * b1[i] * b1[j];
			Kt2(i,j) = fac * (b1[i] * b2[j] + b2[i] * b1[j] + b2[i] * b2[j]);
			Kts(i,j) *= fArea * sigma / (2. * alpha0);
			ek(i,j) = ( Kt1(i,j) + Kt2(i,j) + Kts(i,j) ) / l;
		}
	}
}

void TPZNLMat1d::ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ek,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
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
	
	switch (bc.Type()) {
		case 0 :      // Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				ef(in,0) += gBigNumber * v2[0] * phi(in,0) * weight;
				for (jn = 0 ; jn < phr; jn++) {
					ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
				}
			}
			break;
		case 1 :      // Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(in,0) += v2[0] * phi(in,0) * weight;
			}
			break;
		default :
			PZError << __PRETTY_FUNCTION__ << " boundary condition type not implemented\n";
	}
	
}


void TPZNLMat1d::Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ef)
{
	TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	// TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	TPZVec<REAL> &sol=data.sol[0];
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	TPZFMatrix &axes=data.axes;
	
	
	//   //sol -> (u1,w1,u2,w2)
	//   //dphi = 2/l => l0 = 2/dphi
	//   //x21 = delta x = l cos theta
	//   //z21 = delta y = l sen theta
	//   //theta = arctg (axes (0,1) / axes (0,0) )
	double theta = atan(axes (0,1) / axes (0,0) );
	double l = 2./ dphi(0,0);
	double x21 = l * cos(theta);
	double z21 = l * sin(theta);
	double u1,u2,w1,w2;
	u1 = sol [0];
	u2 = sol [2];
	w1 = sol [1];
	w2 = sol [3];
	double u21 = u2 - u1;
	double w21 = w2 - w1;
	double alpha0 = l/2.;
	
	double eps = Eps(sol,axes,dphi);
	double sigma = eps * fE;
	
	int i;
	double fqi = 0.;
	fqi += 2.*alpha0*fArea*sigma;
	
	//b1
	TPZVec<REAL> b1 (4,0.), b2(4,0.), b(4,0.);
	b1[0] = -x21;  b2[0] = -u21;
	b1[1] =  x21;  b2[1] =  u21;
	b1[2] = -z21;  b2[2] = -w21;
	b1[3] =  z21;  b2[3] =  w21;
	
	double fac = (1. / 4.) / (alpha0 * alpha0);
	for (i=0;i<4;i++)
	{
		b1[i] *= fac;
		b2[i] *= fac;
		b [i] = b1[i] + b2[i];
		ef(i,0) = fqi * b[i] / l;
	}
}

