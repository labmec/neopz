/**
 * @file
 * @brief Contains the implementation of the TPZShapePrism methods.
 */

#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"

using namespace std;

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	

    void TPZShapePrism::ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
    {
        //DebugStop();
        
        int64_t nsides = TPZShapePrism::NSides;
        int nshape;
        
        int linha = 0;
        for (int side = 0; side < nsides; side++)
        {
            
            nshape = 1;
            if(side >= NCornerNodes) nshape = NConnectShapeF(side,order[side-NCornerNodes]);
            int sideorder = 1;
            if(side >= NCornerNodes) sideorder = order[side-NCornerNodes];
            
            TPZGenMatrix<int> locshapeorders(nshape,3);
            SideShapeOrder(side, id, sideorder, locshapeorders);
            
            int nlin = locshapeorders.Rows();
            int ncol = locshapeorders.Cols();
            
            for (int il = 0; il<nlin; il++)
            {
                for (int jc = 0; jc<ncol; jc++)
                {
                    shapeorders(linha, jc) = locshapeorders(il, jc);
                }
                linha++;
            }
        }

    }
    
    
    void TPZShapePrism::SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders)
    {
        //DebugStop();
        if (side<=5)
        {
            if (shapeorders.Rows() != 1)
            {
                DebugStop();
            }
            shapeorders(0,0) = 1;
            shapeorders(0,1) = 0;
            shapeorders(0,2) = 0;
        }
        else if (side>5 && side<15)
        {
            int nshape = order-1;
            if (shapeorders.Rows() != nshape)
            {
                DebugStop();
            }
            for (int ioy = 0; ioy < order-1; ioy++)
            {
                shapeorders(ioy,0) = ioy+2;
            }
        }
        else if (side == 15||side == 19)
        {
            int nshape = (order-2)*(order-1)/2;
            if (shapeorders.Rows() != nshape)
            {
                DebugStop();
            }
            TPZStack<int> lowersides;
            LowerDimensionSides(side, lowersides);
            lowersides.Push(side);
            
            //TPZVec<int> locsideorder(lowersides.size(),order);
            
            int nnodes = NSideNodes(side);
            
            TPZManVector<int64_t, 4> locid(nnodes);
            for (int node=0; node<locid.size(); node++) {
                locid[node] = id[ContainedSideLocId(side, node)];// SideNodeLocId( side, node);
            }// sera que esta pegando os ids corretos mesmo?
            
            TPZGenMatrix<int> locshapeorders(nshape,3);
            
            
            TPZShapeTriang::SideShapeOrder(6,locid, order, locshapeorders);
            
            // temos que arrumar a saida de locshapeorders para adequar a orientacao dos vetores que geram
            // a face do lado side
            
            // aqui o locshapeorders esta armazenado so para x e y
            for (int il = 0; il<nshape; il++)
            {
                shapeorders(il, 0) = locshapeorders(il, 0);
                shapeorders(il, 1) = locshapeorders(il, 1);
                shapeorders(il, 2) = locshapeorders(il, 2);
            }


        }
        else if (side >= 16 && side <=18)
        {
            if (shapeorders.Rows() != (order-1)*(order-1))
            {
                DebugStop();
            }
            TPZStack<int> lowersides;
            LowerDimensionSides(side, lowersides);
            lowersides.Push(side);
            
            //TPZVec<int> locsideorder(lowersides.size(),order);
            
            int nnodes = NSideNodes(side);
            
            TPZManVector<int64_t, 4> locid(nnodes);
            for (int node=0; node<locid.size(); node++) {
                locid[node] = id[ContainedSideLocId(side, node)];// SideNodeLocId( side, node);
            }// sera que esta pegando os ids corretos mesmo?
            
            int nshape = (order-1)*(order-1);
            TPZGenMatrix<int> locshapeorders(nshape,3);
            
            
            TPZShapeQuad::SideShapeOrder(8,locid, order, locshapeorders);
            
            // temos que arrumar a saida de locshapeorders para adequar a orientacao dos vetores que geram
            // a face do lado side
            
            // aqui o locshapeorders esta armazenado so para x e y
            for (int il = 0; il<nshape; il++)
            {
                shapeorders(il, 0) = locshapeorders(il, 0);
                shapeorders(il, 1) = locshapeorders(il, 1);
                shapeorders(il, 2) = locshapeorders(il, 2);
            }

        }
        else
        {   // interno
            int nshape = (order-2)*(order-1)*(order-1)/2;
            if (shapeorders.Rows() != nshape) {
                DebugStop();
            }
            int count = 0;
            int ord1 = order - 2;
            int ord2 = order - 1;
            for (int i=0; i<ord1; i++) {
                for (int j=0; j<ord1; j++) {
                    for (int k=0; k<ord2; k++) {
                        int a = i;
                        int b = j;
                        int c = k;
                        int maxAB = a+b;//a>b? a : b;
                        if (   ( (a+b)<ord1 )  && (c < ord2)   ) // Duvida
                        {
                            shapeorders(count,0) = 3 + maxAB;
                            shapeorders(count,1) = 3 + maxAB;
                            shapeorders(count,2) = 2 + c;
                            count++;
                        }
                        
                    }
                    
                }
            }

//            int orderplus1 = order + 2;
//            int orderplus2 = order + 1;
//            for (int i=3; i<orderplus1; i++) {
//                for (int j=3; j<orderplus1; j++) {
//                    for (int k=3; k<orderplus2; k++) {
//                        int a = i;
//                        int b = j;
//                        int c = k;
//                        if (   ( (a+b)<orderplus1 )  && (c < orderplus2)   ) // Duvida
//                        {
//                            shapeorders(count,0) = a;
//                            shapeorders(count,1) = b;
//                            shapeorders(count,2) = c;
//                            count++;
//                        }
//                        
//                    }
//                    
//                }
//            }

        }

    }
	
	
	
	int TPZShapePrism::NConnectShapeF(int side, int order) {
		if(side<6) return 1;//0 a 4
		if(side<15) return (order-1);//6 a 14
		if(side==15 || side==19) {
			return ((order-2)*(order-1)/2);
		}
		if(side>15 && side<19) {//16,17,18
			return ((order-1)*(order-1));
		}
		if(side==20) {
			return ((order-2)*(order-1)*(order-1)/2);
		}
		PZError << "TPZShapePrism::NConnectShapeF, bad parameter side " << side << endl;
		return 0;
	}
	
	int TPZShapePrism::NShapeF(const TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}

};
