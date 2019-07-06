//
//  TPZMDPMaterial.cpp
//  PZ
//
//  Created by Agnaldo Farias on 19/09/14.
//
//

#include "TPZMDPMaterial.h"
#include "TPZMatLaplacianLagrange.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"
#include "pzaxestools.h"
#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("TPZMDPMaterial.h"));
#endif


using namespace std;

TPZMDPMaterial::TPZMDPMaterial(int nummat, int dim) : TPZMatLaplacian(nummat,dim)
{
    fC[0]=0.;
    fC[1]=0.;
    fC[2]=0.;
}

TPZMDPMaterial::TPZMDPMaterial() : TPZMatLaplacian()
{
    fC[0]=0.;
    fC[1]=0.;
    fC[2]=0.;
}


TPZMDPMaterial & TPZMDPMaterial::operator=(const TPZMDPMaterial &copy)
{
	TPZMatLaplacian::operator = (copy);
    fC[0]=copy.fC[0];
    fC[1]=copy.fC[1];
    fC[2]=copy.fC[2];
	return *this;
}

TPZMDPMaterial::TPZMDPMaterial(const TPZMDPMaterial &copy) : TPZMatLaplacian(copy)
{
    this->operator=(copy);
}

TPZMDPMaterial::~TPZMDPMaterial() {
}


void TPZMDPMaterial::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
    TPZMatLaplacian::Print(out);
}

void TPZMDPMaterial::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    //contribuicoes referente a malha fina:
    TPZFMatrix<REAL>  &phif = datavec[0].phi;
    TPZFMatrix<REAL> &dphif = datavec[0].dphix;
    TPZVec<REAL>  &x = datavec[0].x;
    TPZFMatrix<REAL> &axes = datavec[0].axes;
    
    int phr_f = phif.Rows();
    
    STATE XfLoc = fXf;
    
    if(fForcingFunction) {
        TPZManVector<STATE,1> res(1);
        fForcingFunction->Execute(x,res);
        XfLoc = res[0];
    }
    
    //-------------------------------
    
    STATE ConvDirAx[3];

    switch(fDim) {
        case 1:
            ConvDirAx[0] = axes(0,0)*fC[0]+axes(0,1)*fC[1]+axes(0,2)*fC[2];
            break;
        case 2:
            ConvDirAx[0] = axes(0,0)*fC[0]+axes(0,1)*fC[1]+axes(0,2)*fC[2];
            ConvDirAx[1] = axes(1,0)*fC[0]+axes(1,1)*fC[1]+axes(1,2)*fC[2];
            break;
        case 3:
            ConvDirAx[0] = axes(0,0)*fC[0]+axes(0,1)*fC[1]+axes(0,2)*fC[2];
            ConvDirAx[1] = axes(1,0)*fC[0]+axes(1,1)*fC[1]+axes(1,2)*fC[2];
            ConvDirAx[2] = axes(2,0)*fC[0]+axes(2,1)*fC[1]+axes(2,2)*fC[2];
            break;
        default:
            PZError << "TPZMDPMaterial::Contribute dimension error " << fDim << endl;
    }

    //-------------------------------
    
    /**
     *@ingroup Obter espaco teste otimo: obter e in Vr
     *@brief Find ur in Vr such that: \f$ a(ur,vr) = f(vr), for all vr in Vr \f$
     *@brief being \f$ a(ur,vr) = grad(ur)grad(vr) \f$
     */
    for(int in = 0; in < phr_f; in++)
    {
        int kd;
        ef(in, 0) += (STATE)weight*XfLoc*(STATE)phif(in,0);
        
        
        for(int jn = 0; jn < phr_f; jn++)
        {
            for(kd=0; kd<fDim; kd++)
            {
                ek(in,jn) += (STATE)weight*(fK*(STATE)(dphif(kd,in)*dphif(kd,jn))
                                + (STATE)(dphif(kd,jn)*ConvDirAx[kd]*phif(in)));
            }
        }
    }
    
    
    /**
     *@ingroup orthogonal projection: Vr in Xh
     *@brief Find uh in Xh such that: \f$ b(ur-uh, vh) = 0, for all vh in Xh \f$
     *@brief being \f$ b(ur, vh) - b(uh, vh) = grad(ur)grad(vh) + ur*vh - grad(uh)grad(vh) - uh*vh \f$
     */
    TPZFMatrix<REAL>  &phic = datavec[1].phi;
    TPZFMatrix<REAL> &dphic = datavec[1].dphix;
    int phr_c = phic.Rows();
    
    //b(ur,vh)
    for( int in = 0; in < phr_c; in++ ) {
        int kd;
        for(int jn = 0; jn < phr_f; jn++)
        {
            //ur*vh
            //ek(phr_f + in, jn) += (STATE)weight*((STATE)(phic(in,0)*phif(jn,0)));
            
            for(kd=0; kd<fDim; kd++)
            {
                //grad(ur)*grad(vh)
                ek(phr_f + in, jn) += (STATE)weight*(fC[0]/fK*(STATE)(dphic(kd,in)*dphif(kd,jn)));
            }
        }
    }
    
    //-b(uh,vh)
    for( int in = 0; in < phr_c; in++ ) {
        int kd;
        for(int jn = 0; jn < phr_c; jn++)
        {
            //-(u,v)
           // ek(phr_f + in, phr_f + jn) += (-1.)*(STATE)weight*((STATE)(phic(in,0)*phic(jn,0)));
            
            for(kd=0; kd<fDim; kd++)
            {
                //-grad(uh)grad(vh)
                ek(phr_f + in, phr_f + jn) += (-1.)*(STATE)weight*(fC[0]/fK*(STATE)(dphic(kd,in)*dphic(kd,jn)));
            }
        }
    }
    
#ifdef PZDEBUG
//    if (ek.VerifySymmetry()){
//        cout << __PRETTY_FUNCTION__ << "\nMATRIZ SIMETRICA" << endl;
//        //DebugStop()
//    }
    
#endif
    
}


void TPZMDPMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	
    if(datavec.size()!=2) DebugStop();
    
	TPZFMatrix<REAL>  &phif = datavec[0].phi;
    TPZFMatrix<REAL>  &phic = datavec[1].phi;
     TPZFMatrix<REAL> &axes = datavec[0].axes;
    
	int phrf = phif.Rows();
    int phrc = phic.Rows();
	short in,jn;
	STATE v2[1];
    
	v2[0] = bc.Val2()(0,0);
	if(bc.HasForcingFunction()) {
		TPZManVector<STATE> res(1);
		bc.ForcingFunction()->Execute(datavec[1].x,res);
		v2[0] = res[0];
	}
    
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
            
            for(in = 0 ; in < phrf; in++)
            {
                ef(in,0) += (STATE)(gBigNumber* phif(in,0)*weight)*v2[0];
                for (jn = 0 ; jn < phrf; jn++)
                {
                    ek(in,jn) += gBigNumber * phif(in,0)*phif(jn,0)*weight;
                }
            }

			for(in = 0 ; in < phrc; in++)
            {
				ef(phrf+in,0) += (STATE)(gBigNumber* phic(in,0)*weight)*v2[0];
				for (jn = 0 ; jn < phrc; jn++)
                {
					ek(phrf+in,phrf+jn) += gBigNumber * phic(in,0)*phic(jn,0)*weight;
				}
			}
			break;
            
		case 1 :			// Neumann condition
			for(in = 0 ; in < phrf; in++)
            {
				ef(in,0) += v2[0]*(STATE)(phif(in,0)*weight);
			}
            
            for(in = 0 ; in < phrc; in++)
            {
				ef(phrf+in,0) += v2[0]*(STATE)(phic(in,0)*weight);
			}
			break;
            
		case 2 :		// mixed condition
			for(in = 0 ; in < phrf; in++)
            {
				ef(in, 0) += v2[0]*(STATE)(phif(in, 0)*weight);
				for (jn = 0 ; jn < phrf; jn++)
                {
                    ek(in,jn) += bc.Val1()(0,0)*(STATE)(phif(in,0)*phif(jn,0)*weight);
				}
			}
            
//            for(in = 0 ; in < phrc; in++)
//            {
//				ef(phrf+in, 0) += v2[0]*(STATE)(phic(in, 0)*weight);
//				for (jn = 0 ; jn < phrc; jn++)
//                {
//                    ek(phrf+in,phrf+jn) += bc.Val1()(0,0)*(STATE)(phic(in,0)*phic(jn,0)*weight);
//				}
//			}
			break;
            
        case 3: // outflow condition
            int id, il, jl;
            REAL normal[3];
            if (fDim == 1) PZError << __PRETTY_FUNCTION__ << " - ERROR! The normal vector is not available for 1D TPZInterpolatedElement\n";
            if (fDim == 2){
                normal[0] = axes(0,1);
                //normal[1] = axes(1,1);
            }
            if (fDim == 3){
                normal[0] = axes(0,2);
                normal[1] = axes(1,2);
                normal[2] = axes(2,2);
            }
            REAL ConvNormal = 0.;
            for(id=0; id<fDim; id++) ConvNormal += fC[id]*normal[id];
            if(ConvNormal > 0.) {
                for(il=0; il<phrf; il++) {
                    for(jl=0; jl<phrf; jl++) {
                        ek(il,jl) += weight * ConvNormal * phif(il)*phif(jl);
                    }
                }
            }
            else{
                if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
            }
            break;

	}
}


int TPZMDPMaterial::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str()))        return  1;
	if(!strcmp("ExactSolution",name.c_str()))   return  2;
    if(!strcmp("OptimalTestFunction",name.c_str())) return 3;
    if(!strcmp("Derivative",name.c_str()))      return  4;
	return TPZMaterial::VariableIndex(name);
}

int TPZMDPMaterial::NSolutionVariables(int var){
	if(var==1 || var==2 || var==3) return 1;
	if(var == 4) return fDim;
    
	return TPZMaterial::NSolutionVariables(var);
}



void TPZMDPMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    if (var == 1) {
        Solout[0] = datavec[1].sol[0][0];
        return;
    }
    
    if (var == 3) {
        Solout[0] = datavec[0].sol[0][0];
        return;
    }
    
    TPZVec<STATE> solExata(1);
	TPZFMatrix<STATE> flux(2,1);
    //Exact soluion
	if(var == 2){
		fForcingFunctionExact->Execute(datavec[1].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}//var6
    
    if (var == 4) {
        int id;
        TPZFNMatrix<9,STATE> dsoldx;
        TPZFMatrix<STATE> DSol;
        TPZFMatrix<REAL> axes;
        
        DSol=datavec[1].dsol[0];
        axes=datavec[1].axes;

        TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		for(id=0 ; id<fDim; id++)
        {
			Solout[id] = dsoldx(id,0);
		}
        return;
    }
}

void TPZMDPMaterial::Write(TPZStream &buf, int withclassid) const{
	TPZMatLaplacian::Write(buf, withclassid);
}

void TPZMDPMaterial::Read(TPZStream &buf, void *context){
	TPZMatLaplacian::Read(buf, context);
}

template class TPZRestoreClass<TPZMDPMaterial>;
