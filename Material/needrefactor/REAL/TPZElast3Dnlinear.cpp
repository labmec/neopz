//
//  TPZElast3Dnlinear.cpp
//  PZ
//
//  Created by Cesar Lucci on 22/10/13.
//
//

#include "TPZElast3Dnlinear.h"
#include "pzbndcond.h"

TPZElast3Dnlinear::TPZElast3Dnlinear() : TPZRegisterClassId(&TPZElast3Dnlinear::ClassId),
TPZElasticity3D()
{
    
}

TPZElast3Dnlinear::TPZElast3Dnlinear(int nummat, STATE E, STATE poisson, TPZVec<STATE> &force,
                                     STATE preStressXX, STATE preStressYY, STATE preStressZZ) :
TPZRegisterClassId(&TPZElast3Dnlinear::ClassId),
TPZElasticity3D(nummat, E, poisson, force, preStressXX, preStressYY, preStressZZ)
{
    
}

TPZElast3Dnlinear::~TPZElast3Dnlinear()
{
    
}

void TPZElast3Dnlinear::Contribute(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ef)
{
    DebugStop();//Not implemented!!!
}

void TPZElast3Dnlinear::Contribute(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef)
{
#ifdef PZDEBUG
    if(ek.Rows() != ef.Rows())
    {
        std::cout << "\n\n" << "ek and ef should have same number of rows!!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
    
    if(ef.Cols() != this->fNumLoadCases)
    {
        std::cout << "\n\n" << "ef should have fNumLoadCases equals to its NCols()!!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
#endif
    
    ContributeVecShapeAux(data, weight, ek, ef);
}

void TPZElast3Dnlinear::ContributeBC(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef,
                                     TPZBndCond &bc)
{
#ifdef PZDEBUG
    if(ek.Rows() != ef.Rows())
    {
        std::cout << "\n\n" << "ek and ef should have same number of rows!!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
    
    if(ef.Cols() != this->fNumLoadCases)
    {
        std::cout << "\n\n" << "ef should have fNumLoadCases equals to its NCols()!!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
#endif
    
    ContributeVecShapeBCAux(data, weight, ek, ef, bc);
}

int TPZElast3Dnlinear::ClassId() const{
    return Hash("TPZElast3Dnlinear") ^ TPZElasticity3D::ClassId() << 1;
}

void TPZElast3Dnlinear::FillDataRequirements(TPZMaterialData &data)
{
	TPZMaterial::FillDataRequirements(data);
	data.fNeedsSol = true;
}


//--------------------------------------------------------------------------------


void TPZElast3Dnlinear::ContributeVecShapeAux(TPZMaterialData &data,
                                              REAL weight,
                                              TPZFMatrix<STATE> &ek,
                                              TPZFMatrix<STATE> &ef)
{
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype != data.EVecShape)
    {
        DebugStop();
    }
    
	TPZFMatrix<REAL> & phi = data.phi;
    TPZFMatrix<REAL> & dphi = data.dphix;
    TPZFMatrix<STATE> & dsol = data.dsol[0];
	
	int phc = phi.Cols();
	int efc = ef.Cols();
	
	if(fForcingFunction)
    {
		TPZManVector<STATE> res(3);
		fForcingFunction->Execute(data.x,res);
		fForce[0] = res[0];
		fForce[1] = res[1];
		fForce[2] = res[2];
	}
	
	REAL dvxdx, dvxdy, dvxdz;
    REAL dvydx, dvydy, dvydz;
    REAL dvzdx, dvzdy, dvzdz;
    
    REAL duxdx, duxdy, duxdz;
    REAL duydx, duydy, duydz;
    REAL duzdx, duzdy, duzdz;
    
	/*
	 * Plain strain materials values
	 */
    REAL lambda = fE*fPoisson/((1.+fPoisson)*(1.-2.*fPoisson));
    REAL mu = fE/(2.*(1.+fPoisson));
    
	for( int in = 0; in < phc; in++ )
    {
        //x
		dvxdx = dphi(0,in);
		dvxdy = dphi(1,in);
        dvxdz = dphi(2,in);
        
        //y
		dvydx = dphi(3,in);
		dvydy = dphi(4,in);
		dvydz = dphi(5,in);
        
        //z
        dvzdx = dphi(6,in);
		dvzdy = dphi(7,in);
		dvzdz = dphi(8,in);
		
        for(int col = 0; col < efc; col++)
        {
            ef(in,col) += weight*(  fForce[0] * phi(0, in) +
                                    fForce[1] * phi(1, in) +
                                    fForce[2] * phi(2, in) -
                                    dvxdx * fPreStress[0]  -
                                    dvydy * fPreStress[1]  -
                                    dvzdz * fPreStress[2]    );
            
            //x
            duxdx = dsol(0,0);
            duxdy = dsol(1,0);
            duxdz = dsol(2,0);
            
            //y
            duydx = dsol(0,1);
            duydy = dsol(1,1);
            duydz = dsol(2,1);
            
            //z
            duzdx = dsol(0,2);
            duzdy = dsol(1,2);
            duzdz = dsol(2,2);
            
            REAL eq1 =  duydy*dvxdx*lambda + duzdz*dvxdx*lambda + duxdy*dvydx*mu +
                        duydx*dvydx*mu + duxdz*dvzdx*mu + duzdx*dvzdx*mu +
                        duxdx*dvxdx*(lambda + 2.*mu);
            
            REAL eq2 =  duxdx*dvydy*lambda + duzdz*dvydy*lambda + duxdy*dvxdy*mu +
                        duydx*dvxdy*mu + duydz*dvzdy*mu + duzdy*dvzdy*mu +
                        duydy*dvydy*(lambda + 2.*mu);
            
            REAL eq3 =  duxdx*dvzdz*lambda + duydy*dvzdz*lambda + duxdz*dvxdz*mu +
                        duzdx*dvxdz*mu + duydz*dvydz*mu + duzdy*dvydz*mu +
                        duzdz*dvzdz*(lambda + 2.*mu);
            
            ef(in,col) -= weight * (eq1 + eq2 + eq3);
        }
		for( int jn = 0; jn < phc; jn++ )
        {
            //x
            duxdx = dphi(0,jn);
            duxdy = dphi(1,jn);
            duxdz = dphi(2,jn);
            
            //y
            duydx = dphi(3,jn);
            duydy = dphi(4,jn);
            duydz = dphi(5,jn);
            
            //z
            duzdx = dphi(6,jn);
            duzdy = dphi(7,jn);
            duzdz = dphi(8,jn);
            
            REAL eq1 =  duydy*dvxdx*lambda + duzdz*dvxdx*lambda + duxdy*dvydx*mu +
                        duydx*dvydx*mu + duxdz*dvzdx*mu + duzdx*dvzdx*mu +
                        duxdx*dvxdx*(lambda + 2.*mu);
            
            REAL eq2 =  duxdx*dvydy*lambda + duzdz*dvydy*lambda + duxdy*dvxdy*mu +
                        duydx*dvxdy*mu + duydz*dvzdy*mu + duzdy*dvzdy*mu +
                        duydy*dvydy*(lambda + 2.*mu);
            
            REAL eq3 =  duxdx*dvzdz*lambda + duydy*dvzdz*lambda + duxdz*dvxdz*mu +
                        duzdx*dvxdz*mu + duydz*dvydz*mu + duzdy*dvydz*mu +
                        duzdz*dvzdz*(lambda + 2.*mu);
            
            ek(in,jn) += weight * (eq1 + eq2 + eq3);
		}
	}
}


void TPZElast3Dnlinear::ContributeVecShapeBCAux(TPZMaterialData &data,
                                                REAL weight,
                                                TPZFMatrix<STATE> &ek,
                                                TPZFMatrix<STATE> &ef,
                                                TPZBndCond &bc)
{
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype != data.EVecShape)
    {
        DebugStop();
    }
    
    TPZFMatrix<REAL> & phi = data.phi;
    TPZManVector<STATE,3> sol = data.sol[0];
    
	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    
	int phc = phi.Cols();
	short in,jn;
	
	switch (bc.Type())
    {
		case 0:// Dirichlet condition
        {
			for(in = 0 ; in < phc; in++)
            {
                for(int il = 0; il < fNumLoadCases; il++)
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(in,il) += weight * BIGNUMBER * ( v2(0,il)*phi(0,in) + v2(1,il)*phi(1,in) + v2(2,il)*phi(2,in) );
                    
                    ef(in,il) -= weight * BIGNUMBER * ( sol[0]*phi(0,in) + sol[1]*phi(1,in) + sol[2]*phi(2,in) );
                }
				for(jn = 0 ; jn < phc; jn++)
                {
                    ek(in,jn) += weight * BIGNUMBER * ( phi(0,jn)*phi(0,in) + phi(1,jn)*phi(1,in) + phi(2,jn)*phi(2,in) );
				}
			}
			break;
		}
		case 1:// Neumann condition
        {
            for(in = 0; in < phc; in++)
            {
                for(int il = 0; il < fNumLoadCases; il++)
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(in,il) += weight * ( v2(0,il)*phi(0,in) + v2(1,il)*phi(1,in) + v2(2,il)*phi(2,in) );
                }
                //... continua no TPZPlaneFractCouplingMat
            }
			break;
		}
		case 2:// condicao mista
        {
            DebugStop();
			for(in = 0 ; in < phc; in++)
            {
                for(int il = 0; il < fNumLoadCases; il++)
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(in,il) += weight * ( v2(0,il)*phi(0,in) + v2(1,il)*phi(1,in) + v2(2,il)*phi(2,in) );
                    
                    ef(in,il)  -= (   bc.Val1()(0,0)*phi(0,in)*sol[0]*weight
                    
                                    + bc.Val1()(1,0)*phi(1,in)*sol[0]*weight
                                    
                                    + bc.Val1()(2,0)*phi(2,in)*sol[0]*weight
                                    
                                    
                                    + bc.Val1()(0,1)*phi(0,in)*sol[1]*weight
                                    
                                    + bc.Val1()(1,1)*phi(1,in)*sol[1]*weight
                                    
                                    + bc.Val1()(2,1)*phi(2,in)*sol[1]*weight
                                    
                                    
                                    + bc.Val1()(0,2)*phi(0,in)*sol[2]*weight
                                    
                                    + bc.Val1()(1,2)*phi(1,in)*sol[2]*weight
                                    
                                    + bc.Val1()(2,2)*phi(2,in)*sol[2]*weight  );
                }
				
				for(jn = 0; jn <phc; jn++)
                {
                    
                    ek(in,jn)  += bc.Val1()(0,0)*phi(0,in)*phi(0,jn)*weight
                    
                                + bc.Val1()(1,0)*phi(1,in)*phi(0,jn)*weight
                                
                                + bc.Val1()(2,0)*phi(2,in)*phi(0,jn)*weight
                                
                                
                                + bc.Val1()(0,1)*phi(0,in)*phi(1,jn)*weight
                                
                                + bc.Val1()(1,1)*phi(1,in)*phi(1,jn)*weight
                                
                                + bc.Val1()(2,1)*phi(2,in)*phi(1,jn)*weight
                                
                                
                                + bc.Val1()(0,2)*phi(0,in)*phi(2,jn)*weight
                                
                                + bc.Val1()(1,2)*phi(1,in)*phi(2,jn)*weight
                                
                                + bc.Val1()(2,2)*phi(2,in)*phi(2,jn)*weight;
                }
			}
            break;
        }
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
        {
            for(in = 0 ; in < phc; in++)
            {
                for(int il = 0; il < fNumLoadCases; il++)
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(in,il) -= weight * BIGNUMBER * ( v2(0,il)*phi(0,in)*sol[0] + v2(1,il)*phi(1,in)*sol[1] + v2(2,il)*phi(2,in)*sol[2] );
                    for(jn = 0 ; jn < phc; jn++)
                    {
                        ek(in,jn) += weight * BIGNUMBER * ( v2(0,il)*phi(0,in)*phi(0,jn) + v2(1,il)*phi(1,in)*phi(1,jn) + v2(2,il)*phi(2,in)*phi(2,jn) );
                    }
                }
			}
            break;
        }
        case 4: // stressField Neumann condition
        {
            DebugStop();//Nao implementado!!!
            break;
        }
        default:
            PZError << "TPZElastitity3D::ContributeBC error - Wrong boundary condition type" << std::endl;
	}
}
