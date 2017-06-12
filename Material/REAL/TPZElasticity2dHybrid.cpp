/**
 * @file
 * @brief Contains implementations of the TPZElasticityMaterial methods.
 */

#include "TPZElasticity2DHybrid.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.elasticity.data"));
#endif

#include <fstream>
using namespace std;

TPZElasticity2DHybrid::TPZElasticity2DHybrid() : TPZElasticityMaterial(0) {
}

TPZElasticity2DHybrid::TPZElasticity2DHybrid(int id) : TPZElasticityMaterial(id) {
}

TPZElasticity2DHybrid::TPZElasticity2DHybrid(int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress) : TPZElasticityMaterial(num,E,nu,fx,fy,plainstress) {
}

TPZElasticity2DHybrid::~TPZElasticity2DHybrid() {
}




void TPZElasticity2DHybrid::ContributeBC(TPZMaterialData &data,REAL weight,
										 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
    
    
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype==data.EVecShape){
        ContributeVecShapeBC(data,weight,ek, ef,bc);
        return;
    }
    
	TPZFMatrix<REAL> &phi = data.phi;
     int dim = Dimension();

	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    
	int phr = phi.Rows();
	short in,jn;
    
    if (ef.Cols() != bc.NumLoadCases()) {
        DebugStop();
    }
	
//		In general when the problem is  needed to stablish any convention for ContributeBC implementations

//     REAL v2[2];
// 	v2[0] = bc.Val2()(0,0);
// 	v2[1] = bc.Val2()(1,0);
    int nstate = NStateVariables();

    TPZFMatrix<STATE> &v1 = bc.Val1();


    switch (bc.Type()) {
        case 1 :			// Neumann condition
        {
            for(in = 0 ; in < phr; in++) {
                for (int il = 0; il<NumLoadCases(); il++)
                {
                    REAL v2[2];
                    v2[0] = bc.Val2(il)(0,0);
                    v2[1] = bc.Val2(il)(1,0);
                    ef(2*in,il)   += BIGNUMBER * v2[0] * phi(in,0) * weight;        // forced v2 displacement
                    ef(2*in+1,il) += BIGNUMBER * v2[1] * phi(in,0) * weight;        // forced v2 displacement
                }
                for (jn = 0 ; jn < phi.Rows(); jn++)
                {
                    ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
                }
            }
        }
            break;
            
        case 0 :		// Dirichlet condition
        {
            for (in = 0; in < phr; in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(2*in,il) += -v2(0,0) * phi(in,0) * weight;        // force in x direction
                    ef(2*in+1,il) +=  -v2(1,0) * phi(in,0) * weight;      // force in y direction
                }
            }
        }
            break;
            
        case 2 :		// Mixed Condition
        {
            for(in = 0 ; in < phi.Rows(); in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(2*in,il) += v2(0,0) * phi(in,0) * weight;        // force in x direction
                    ef(2*in+1,il) += v2(1,0) * phi(in,0) * weight;      // forced in y direction
                }
                
                for (jn = 0 ; jn < phi.Rows(); jn++) {
                    ek(2*in,2*jn) += bc.Val1()(0,0) * phi(in,0) * phi(jn,0) * weight;         // peso de contorno => integral de contorno
                    ek(2*in+1,2*jn) += bc.Val1()(1,0) * phi(in,0) * phi(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += bc.Val1()(1,1) * phi(in,0) * phi(jn,0) * weight;
                    ek(2*in,2*jn+1) += bc.Val1()(0,1) * phi(in,0) * phi(jn,0) * weight;
                }
            }   // este caso pode reproduzir o caso 0 quando o deslocamento
            
            break;
            
        }
    }
}


void TPZElasticity2DHybrid::ContributeVecShapeBC(TPZMaterialData &data,REAL weight,
										 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
    
    TPZFMatrix<REAL> &phi = data.phi;
    
	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    
	int phc = phi.Cols();
	short in,jn;
	
	switch (bc.Type()) {
		case 1 :			// Dirichlet condition
			for(in = 0 ; in < phc; in++) {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    
                    ef(in,il) += weight*BIGNUMBER*(v2(0,il)*phi(0,in) + v2(1,il) * phi(1,in));
                }
				for (jn = 0 ; jn < phc; jn++) {
                    
                    ek(in,jn) += weight*BIGNUMBER*(phi(0,in)*phi(0,jn) + phi(1,in)*phi(1,jn));
				}
			}
			break;
			
		case 0 :			// Neumann condition
            for (in = 0; in < phc; in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(in,il)+= -weight*(v2(0,il)*phi(0,in) + v2(1,il)*phi(1,in));
                }
            }
			break;
			
		case 2 :		// condicao mista
			for(in = 0 ; in < phc; in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                     ef(in,il) += weight * (v2(0,il)*phi(0,in) + v2(1,il)*phi(1,in));
                }
				
				for (jn = 0; jn <phc; jn++) {
                    
                    ek(in,jn) += bc.Val1()(0,0)*phi(0,in)*phi(0,jn)*weight 
                    
                                + bc.Val1()(1,0)*phi(1,in)*phi(0,jn)*weight
                    
                                + bc.Val1()(0,1)*phi(0,in)*phi(1,jn)*weight
                    
                                + bc.Val1()(1,1)*phi(1,in)*phi(1,jn)*weight;
				}
			}// este caso pode reproduzir o caso 0 quando o deslocamento
	}      //  eh nulo introduzindo o BIGNUMBER pelos valores da condicao
}


TPZElasticity2DHybrid::TPZElasticity2DHybrid(const TPZElasticity2DHybrid &copy) :
TPZElasticityMaterial(copy)
{
}


int TPZElasticity2DHybrid::ClassId() const
{
    return /** @brief Id of Elasticity material */
    TPZELASTICITY2DHYBRIDMATERIALID;

}

template class TPZRestoreClass<TPZElasticity2DHybrid,TPZELASTICITY2DHYBRIDMATERIALID>;

void TPZElasticity2DHybrid::Read(TPZStream &buf, void *context)
{
	TPZElasticityMaterial::Read(buf,context);
	
}

void TPZElasticity2DHybrid::Write(TPZStream &buf, int withclassid)
{
	TPZElasticityMaterial::Write(buf,withclassid);
	
}

