/**
 * @file
 * @brief Contains the implementation of the TPZShapeLinear methods.
 */

#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzerror.h"
#include "pzreal.h"
using namespace std;

namespace pzshape {
	
	
	void TPZShapeLinear::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
		// Quadratic or higher shape functions
		if(num <= 0) return;
        phi.Zero();
        dphi.Zero();
		phi.Put(0,0,1.0);
		dphi.Put(0,0, 0.0);
		if(num == 1) return;
		phi.Put(1,0, x);
		dphi.Put(0,1, 1.0);
		int ord;
       // dphi.Print("DphisAntes = ",std::cout,EMathematicaInput);
    
		for(ord = 2;ord<num;ord++) {
			phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
			dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
		}
       // dphi.Print("DphisDepois = ",std::cout,EMathematicaInput);
	}
	
    void TPZShapeLinear::Expo(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
        // Quadratic or higher shape functions
        if(num <= 0) return;
        phi.Put(0,0,1.0);
        dphi.Put(0,0, 0.0);
        if(num == 1) return;
        phi.Put(1,0, x);
        dphi.Put(0,1, 1.0);
        int ord;
        for(ord = 2;ord<num;ord++) {
            phi.Put(ord,0, x*phi(ord-1,0));
            dphi.Put(0,ord, ord*phi(ord-1,0));
        }
    }
    
	void TPZShapeLinear::Legendre(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
		
		// Quadratic or higher shape functions
		if (num <= 0) return;
		phi.Put(0, 0, 1.0);
		dphi.Put(0, 0, 0.0);
		
		if (num == 1) return;
		
		phi.Put(1, 0, x);
		dphi.Put(0, 1, 1.0);
		
		int ord;
		//Aqui fica diferente do Chebyshev
		REAL ord_real, value;
		for (ord = 2; ord < num; ord++)
		{
			//casting int ord to REAL ord_real
			ord_real = (REAL)ord;
			//computing the ord_th function
			value    = ( (2.0 * (ord_real - 1.0) + 1.0) * x * phi(ord - 1, 0) - (ord_real - 1.0) * phi(ord - 2 , 0) ) / (ord_real);
			phi.Put(ord, 0, value);
			
			//computing the ord_th function's derivative
			value    = (2.0 * (ord_real - 1.0) + 1.0) * phi(ord - 1, 0) + dphi(0, ord - 2);
			dphi.Put(0, ord, value);
		}
		
#ifdef PZDEBUG
		int printing = 0;
		if (printing){
			cout << "Legendre" << endl;
			for(ord = 0; ord < num; ord++)
			{
				cout << "x = " << x << endl;
				cout << "phi(" << ord << ", 0) = " << phi(ord, 0) << endl;
				cout << "dphi(0, " << ord << " = " << dphi(0, ord) << endl;
				cout << endl;
			}
		}
#endif
		
	} //end of method
	
	void TPZShapeLinear::Legendre(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, int nderiv){
		
		// Quadratic or higher shape functions
		if (num <= 0) return;
		phi.Put(0, 0, 1.0);
		dphi.Put(0, 0, 0.0);
		
		int ideriv;
		for (ideriv = 1; ideriv < nderiv; ideriv++)
			dphi.Put(ideriv, 0, 0.0);
		
		
		if (num == 1) return;
		
		phi.Put(1, 0, x);
		dphi.Put(0, 1, 1.0);
		
		for (ideriv = 1; ideriv < nderiv; ideriv++)
			dphi.Put(ideriv, 1, 0.0);
		
		int ord;
		//Aqui fica diferente do Chebyshev
		REAL ord_real, value;
		for (ord = 2; ord < num; ord++)
		{
			//casting int ord to REAL ord_real
			ord_real = (REAL)ord;
			//computing the ord_th function
			value    = ( (2.0 * (ord_real - 1.0) + 1.0) * x * phi(ord - 1, 0) - (ord_real - 1.0) * phi(ord - 2 , 0) ) / (ord_real);
			phi.Put(ord, 0, value);
			
			//computing the ord_th function's derivative
			value    = (2.0 * (ord_real - 1.0) + 1.0) * phi(ord - 1, 0) + dphi(0, ord - 2);
			dphi.Put(0, ord, value);
			
			for (ideriv = 1; ideriv < nderiv; ideriv++){
				value = (2.0 * (ord_real - 1.0) + 1.0) * dphi(ideriv - 1, ord - 1) + dphi(ideriv, ord - 2);
				dphi.Put(ideriv, ord, value);	    	 
			}
			
		}
		
#ifdef PZDEBUG
		int printing = 0;
		if (printing){
			cout << "Legendre" << endl;
			for(ord = 0; ord < num; ord++)
			{
				cout << "x = " << x << endl;
				cout << "phi(" << ord << ", 0) = " << phi(ord, 0) << endl;
				cout << "dphi(0, " << ord << " = " << dphi(0, ord) << endl;
				cout << endl;
			}
		}
#endif
		
	} //end of method
	
    /**
     * @brief Jacobi parameters orthogonal polynomials
     */
    REAL TPZShapeLinear::JacobiA = 0.5;
    REAL TPZShapeLinear::JacobiB = 0.5;
    
	void TPZShapeLinear::Jacobi(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		// Quadratic or higher shape functions
        REAL a = JacobiA, b = JacobiB;
		if(num <= 0) return;
		phi.Put(0,0,1.0);
		dphi.Put(0,0, 0.0);
		if(num == 1) return;
		phi.Put(1,0, 0.5 * (a - b + (2.0 + a + b) * x));
		dphi.Put(0,1, 0.5*(num + a + b + 1)*phi(num-1,0) );
        // Following http://mathworld.wolfram.com/JacobiPolynomial.html
        REAL A_ab=0.0, B_ab=0.0, C_ab=0.0;
		for(int ord = 2;ord<num;ord++) {
            //Computing Coefficients for three terms recursion relation
            A_ab = ((2.0*ord+a+b+1)*(2.0*ord+a+b+2))/(2.0*(ord+1)*(ord+a+b+1));
            B_ab = ((b*b-a*a)*(2.0*ord+a+b+1))/(2.0*(ord+1)*(ord+a+b+1)*(2.0*ord+a+b));
            C_ab = ((ord+a)*(ord+b)*(2.0*ord+a+b+2))/((ord+1)*(ord+a+b+1)*(2.0*ord+a+b));
			phi.Put(ord,0, (A_ab*x - B_ab)*phi(ord-1,0) - C_ab*phi(ord-2,0));
			dphi.Put(0,ord, 0.5*(ord + a + b + 1)*phi(ord-1,0));
		}
	} //end of method
    
	void TPZShapeLinear::Hermite(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		// Quadratic or higher shape functions (for physicists)
		if(num <= 0) return;
		phi.Put(0,0,1.0);
		dphi.Put(0,0, 0.0);
		if(num == 1) return;
		phi.Put(1,0, 2.0*x);
		dphi.Put(0,1, 2.0);
        // Following http://mathworld.wolfram.com/HermitePolynomial.html
		for(int ord = 2;ord<num;ord++) {
			phi.Put(ord,0, 2.0*x*phi(ord-1,0) - (2.0*(REAL)ord-1.0)*phi(ord-2,0));
			dphi.Put(0,ord, (2.0*(REAL)ord-1.0)*phi(ord-2,0));
		}
	} //end of method
	
	// Setting Chebyshev polynomials as orthogonal sequence generating shape functions
	void (*TPZShapeLinear::fOrthogonal)(REAL, int, TPZFMatrix<REAL> &, TPZFMatrix<REAL> &) = TPZShapeLinear::Chebyshev;
	
	

	
    void TPZShapeLinear::ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
    {
        int nshape = 2+(order[0]-1);
        if (shapeorders.Rows() != nshape) {
            DebugStop();
        }
        shapeorders(0,0) = 1;
        shapeorders(1,0) = 1;
        for (int i=2; i<nshape; i++) {
            shapeorders(i,0) = i;
        }
    }
    
    
    void TPZShapeLinear::SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders)
    {
        DebugStop();
    }
    
	
	int TPZShapeLinear::NConnectShapeF(int side, int order) {
		if(side<2) return 1;//0 a 4
		if(side<3) return (order-1);//6 a 14
		PZError << "TPZShapeLinear::NConnectShapeF, bad parameter side " << side << endl;
		return 0;
	}
	
	int TPZShapeLinear::NShapeF(const TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}
	
	
	void TPZShapeLinear::Chebyshev(const FADREAL & x,int num,TPZFMatrix<FADREAL> &phi,
                                   TPZFMatrix<FADREAL> &dphi){
		// Quadratic or higher shape functions
		if(num <= 0) return;
        FADREAL zero(x.size(),0.);
		dphi.Put(0,0, zero*0.0);
		phi(0,0) = 1.0+zero; // <!> Remark: the derivatives other than the 0th are set to null
		if(num == 1) return;
		phi.Put(1,0, x);
		dphi.Put(0,1, 1.0+zero);
		//phi[1].fastAccessDx(0)=1.0; // <!> Remark: the derivatives other than the 0th aren't set to null
		//just ensuring the derivatives are properly initialized and that FAD objects of more than
		// one derivative are used
		int ord;
		for(ord = 2;ord<num;ord++) {
			phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
			dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
		}
	}
	
	void (*TPZShapeLinear::FADfOrthogonal)(const FADREAL&,int ,TPZFMatrix<FADREAL> &, TPZFMatrix<FADREAL> &dphi) =  TPZShapeLinear::Chebyshev/*(FADREAL&, int, TPZVec<FADREAL>&)*/;//Chebyshev;


};
