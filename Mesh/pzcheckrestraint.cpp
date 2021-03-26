/**
 * @file
 * @brief Contains the implementation of the TPZCheckRestraint methods.
 */

#include "pzcheckrestraint.h"
#include "pzintel.h"
#include "pzcmesh.h"
#include "pztrnsform.h"
#include "pzgeoelside.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzcheckgeom.h"
#include "TPZMaterial.h"

using namespace std;

TPZCheckRestraint::TPZCheckRestraint(TPZCompElSide smalll, TPZCompElSide large) {
	// Stores the small computational element with commom side
	fSmall = smalll;
	// Gets one element in lower level than small
	fLarge = smalll.LowerLevelElementList(1);
	// Checking whether element in lower level is neighbour than large element with commom side
	if (!large.Reference().NeighbourExists(fLarge.Reference())) {
		cout << "TPZCheckRestraint created for a wrong large  element\n";
	}
	// Stores large element with commom side
	fLarge = large;
    
	// The small element as interpolate element
	TPZInterpolatedElement *smallel = dynamic_cast<TPZInterpolatedElement *> (smalll.Element());
	// Gets the number of shape functions for small element over the commom side
	int smallsize = smallel->NSideShapeF(fSmall.Side());
	// Gets the number of connects over the commom side for small element
	int nsmallconnect = smallel->NSideConnects(fSmall.Side());
	int smallside = fSmall.Side();

	// The large element as interpolate element
	TPZInterpolatedElement *largel = dynamic_cast<TPZInterpolatedElement *> (fLarge.Element());
	if(!largel) {
		cout << "TPZCheckRestraint created for an element/side without restraint\n";
		largel = smallel;
		fLarge = fSmall;
	}
	// Gets the number of shape functions for large element over the commom side
	int largesize = largel->NSideShapeF(fLarge.Side());
	// Gets the number of connects over the commom side for large element
	int nlargeconnect = largel->NSideConnects(fLarge.Side());
	int largeside = fLarge.Side();
	
	// Redimensioning of restraint matrix
	fRestraint.Redim(smallsize,largesize);
	fMesh = smallel->Mesh();
	int nmat = fMesh->MaterialVec().size();  
	int nstate = 1;
	if(nmat) {
		std::map<int, TPZMaterial * >::iterator mit = fMesh->MaterialVec().begin();
		nstate = mit->second->NStateVariables();
	}

	fSmallSize.Resize(nsmallconnect);
	fSmallPos.Resize(nsmallconnect);
	fSmallConnect.Resize(nsmallconnect);
	fLargeSize.Resize(nlargeconnect);
	fLargePos.Resize(nlargeconnect);
	fLargeConnect.Resize(nlargeconnect);
	int ic;
	int nc = nsmallconnect;
	if(nc) fSmallPos[0] = 0;
	for(ic=0; ic<nc; ic++) {
		int connectindex = smallel->SideConnectLocId(ic,smallside);
        TPZConnect &connect = smallel->Connect(connectindex);
		fSmallSize[ic] = smallel->NConnectShapeF(connectindex,connect.Order());
		if(smallel->Connect(connectindex).CheckDependency(fSmallSize[ic], fMesh, nstate) == -1) {
			cout << "TPZCheckRestraint incompatible small element connect dimensions\n";
		}
		if(ic) fSmallPos[ic] = fSmallPos[ic-1] + fSmallSize[ic-1];
		fSmallConnect[ic] = smallel->SideConnectIndex(ic,smallside);
	}
	nc = nlargeconnect;
	if(nc) fLargePos[0] = 0;
	for(ic=0; ic<nc; ic++) {
		int connectindex = largel->SideConnectLocId(ic,largeside);
        TPZConnect &c = largel->Connect(connectindex);
		fLargeSize[ic] = largel->NConnectShapeF(connectindex,c.Order());
		if(largel->Connect(connectindex).CheckDependency(fLargeSize[ic], fMesh, nstate) == -1 ){
			cout << "TPZCheckRestraint incompatible large element connect dimensions\n";
		}
		if(ic) fLargePos[ic] = fLargePos[ic-1] + fLargeSize[ic-1];
		fLargeConnect[ic] = largel->SideConnectIndex(ic,largeside);
	}
	nc = nsmallconnect;
	for(ic=0; ic<nc; ic++) {
		AddConnect(fSmallConnect[ic]);
	}
}

int TPZCheckRestraint::SmallConnect(int connectid) {
	int nc = fSmallConnect.NElements();
	int ic;
	for(ic=0; ic<nc; ic++) if(fSmallConnect[ic] == connectid) return ic;
	return -1;
}

int TPZCheckRestraint::LargeConnect(int connectid) {
	int nc = fLargeConnect.NElements();
	int ic;
	for(ic=0; ic<nc; ic++) if(fLargeConnect[ic] == connectid) return ic;
	return -1;
}

void TPZCheckRestraint::AddConnect(int connectindex) {
	int smalll = SmallConnect(connectindex);
	TPZConnect &smallc = fMesh->ConnectVec()[connectindex];
	int large = LargeConnect(connectindex);
	
	if(large != -1) {
		if(smalll == -1) {
			cout << "TPZCheckRestraint::AddConnect data structure error 0\n";
			return;
		}
		int firstl = fSmallPos[smalll];
		//    int lastl = firstl+fSmallSize[small];
		int firstc = fLargePos[large];
		int lastc = firstc+fLargeSize[large];
		int ic,il;
		for(ic=firstc,il=firstl; ic<lastc; ic++,il++) {
			fRestraint(il,ic) = 1.;
		}
	} else {
		TPZConnect::TPZDepend *depend = smallc.FirstDepend();
		if(!depend) {
			cout << "TPZCheckRestraint::AddConnect data structure error 1\n";
		}
		while(depend) {
			AddDependency(connectindex,depend->fDepConnectIndex,depend->fDepMatrix);
			depend = depend->fNext;
		}
	}
}

void TPZCheckRestraint::AddDependency(int smallconnectindex, int largeconnectindex, TPZFMatrix<REAL> &dependmatrix) {
	
	int smalll = SmallConnect(smallconnectindex);
	//  TPZConnect &smallc = fMesh->ConnectVec()[smallconnectindex];
	int large = LargeConnect(largeconnectindex);
	if(large != -1) {
		int firstl = fSmallPos[smalll];
		int lastl = firstl+fSmallSize[smalll];
		int firstc = fLargePos[large];
		int lastc = firstc+fLargeSize[large];
		int ic,il;
		if (firstc != lastc && firstl != lastl && (firstc < 0 || lastc > fRestraint.Cols() || firstl < 0 || lastl > fRestraint.Rows())){
			cout << "TPZCheckRestraint::AddDependency : indexing error\n";
			int a;
			cin >> a;
			return;
		}
		if ((lastc - firstc) != dependmatrix.Cols() ||  (lastl-firstl) != dependmatrix.Rows()){
			cout << "TPZCheckRestraint::AddDependency : incompatible dimensions\n";
			int a;
			cin >> a;
			return;
		}
		for(il=firstl; il<lastl; il++) {
			int line = il - firstl;
			for (ic=firstc; ic<lastc; ic ++){
				int column = ic - firstc;
				if (ic < 0 || ic >= fRestraint.Cols() || il < 0 || il >= fRestraint.Rows()){
					cout << "TPZCheckRestraint::AddDependency : Should never pass here, was already checked 1\n";
					int a;
					cin >> a;
					return;
				}
				if (line >= dependmatrix.Rows() || column >= dependmatrix.Cols()){
					cout << "TPZCheckRestraint::AddDependency : Should never pass here, was already checked 2\n";
					int a;
					cin >> a;
					return;
				}
				fRestraint(il,ic) += (REAL)dependmatrix(line,column);
			}
		}
	} else {
		TPZConnect &largec = fMesh->ConnectVec()[largeconnectindex];
		TPZConnect::TPZDepend *depend = largec.FirstDepend();
		if(!depend) {
			// recado de erro
			cout << "TPZCheckRestraint::Error:\tLarge connect without dependency\n";
			cout.flush();
			return;
		}
		while(depend) {
			// comparar dimensao das matrizes
			int deprows,depcols;
			deprows = depend->fDepMatrix.Rows();
			depcols = dependmatrix.Cols();
			
			if (deprows != depcols){
				cout << "TPZCheckRestraint::Error:\tLarge connect without dependency\n";
				cout.flush();
				//return;
			}
			
			TPZFMatrix<REAL> depmat = dependmatrix * depend->fDepMatrix;
			AddDependency(smallconnectindex,depend->fDepConnectIndex,depmat);
			depend = depend->fNext;
		}
	}
}

TPZFMatrix<REAL>& TPZCheckRestraint::RestraintMatrix(){
	
	return fRestraint;
	
}

int TPZCheckRestraint::CheckRestraint() {
	
	TPZInterpolatedElement *smallel = dynamic_cast<TPZInterpolatedElement *> (fSmall.Element());
	TPZInterpolatedElement *largel = dynamic_cast<TPZInterpolatedElement *> (fLarge.Element());
	TPZIntPoints *intrule = smallel->Reference()->CreateSideIntegrationRule(fSmall.Side(),2);
	int dims = smallel->Reference()->SideDimension(fSmall.Side());
	int diml = largel->Reference()->SideDimension(fLarge.Side());
	TPZVec<int> ord(dims,5);
	intrule->SetOrder(ord);
	TPZGeoElSide smallside = fSmall.Reference();
	TPZGeoElSide largeside = fLarge.Reference();
	TPZTransform<> t(smallside.Dimension());
	smallside.SideTransform3(largeside,t);
	
	TPZTransform<> T = smallel->Reference()->ComputeParamTrans(largel->Reference(),fLarge.Side(), fSmall.Side());//transforma��o direta, sem acumulo
	if(T.CompareTransform(t))//caso erro � maior que tol=1.e-6 retorna 1
		PZError << "TPZCheckRestraint::CheckRestraint transformation error!\n";
	
	int numint = intrule->NPoints();
	int numshapes = fRestraint.Rows();
	int numshapel = fRestraint.Cols();
	TPZFMatrix<REAL> phis(numshapes,1),dphis(dims,numshapes),phil(numshapel,1),dphil(diml,numshapel);
	TPZFMatrix<REAL> philcheck(numshapel,1);
	TPZVec<REAL> points(dims),pointl(diml),point(3);
	int in,check = 0;
	REAL weight,error;
	for(int it=0; it<numint; it++) {
		intrule->Point(it,points,weight);
		smallel->SideShapeFunction(fSmall.Side(),points,phis,dphis);
		t.Apply(points,pointl);
		largel->SideShapeFunction(fLarge.Side(),pointl,phil,dphil);
		fRestraint.Multiply(phis,philcheck,1);
		error = 0.;
		for(in=0; in<numshapel; in++) {
			for(auto x = 0; x < philcheck.Cols(); x++){
				error += (philcheck(in,x)-phil(in,x))*(philcheck(in,x)-phil(in,x));
			}
		}
		error = sqrt(error/numshapel);
		if(error > 1.e-6) {
		    if(check == 0){
                std::cout<<std::setw(10)<<"fRestraint:"<<std::endl;
                for(auto x = 0; x < fRestraint.Rows(); x++) {
                    for (auto y = 0; y < fRestraint.Cols(); y++) {
                        std::cout << std::setw(10) << fRestraint(x, y) << "\t";
                    }
                    std::cout<<std::endl;
                }
                std::cout<<std::endl;
		    }
			check = 1;
            for(in=0; in<numshapel; in++) {
                std::cout<<std::setw(10)<<"philcheck:";
                for(auto x = 0; x < philcheck.Cols(); x++){
                    std::cout<<std::setw(10)<<philcheck(in,x)<<"\t";
                }
                std::cout<<std::endl;
                std::cout<<std::setw(10)<<"phil:";
                for(auto x = 0; x < phil.Cols(); x++){
                    std::cout<<std::setw(10)<<phil(in,x)<<"\t";
                }
                std::cout<<std::endl;
            }

            for(in=0; in<numshapes; in++) {
                std::cout<<std::setw(10)<<"phis:";
                for(auto x = 0; x < phis.Cols(); x++){
                    std::cout<<std::setw(10)<<phis(in,x)<<"\t";
                }
                std::cout<<std::endl;
            }
			cout << "TPZCheckRestraint failed error = " << error << endl;
		}
	}
	delete intrule;
	return check;
}

void TPZCheckRestraint::Print(ostream &out){
	
	fSmall.Element()->Print(out);
	out <<"Small side: " <<fSmall.Side() << endl;
	fLarge.Element()->Print(out);
	out << "Large side: " <<fLarge.Side() << endl;
	fRestraint.Print("Restraint Matrix",out);
	
	out <<  "hierarqui of elements ";
	fSmall.Reference().Element()->Print(out);
	TPZGeoElSide smallgeo = fSmall.Reference();
	TPZGeoElSide largegeo = fLarge.Reference();
	while(smallgeo.Exists() && !smallgeo.NeighbourExists(largegeo))
    {
		out <<  "Small element side =  " << smallgeo.Side() << endl;
		out <<  "Small geometric element printout \n";
		smallgeo.Element()->Print();
		smallgeo = smallgeo.Element()->Father2(smallgeo.Side());
    }
	if(!smallgeo.Exists()) 
    {
		out <<  "TPZCheckRestraint::Print inconsistent datastructure\n";
    } else {
		out <<  "Small element side =  " << smallgeo.Side() << endl;
		out <<  "Small geometric element printout \n";
		smallgeo.Element()->Print();
    }
	out <<  "Large element side =  " << largegeo.Side() << endl;
	out <<  "Large geometric element printout \n";
	largegeo.Element()->Print();
	
	int nsmal = fSmallConnect.NElements();
	int nlarge = fLargeConnect.NElements();
	int is, il;
	
	out << "fSmallConnect" << endl;
	for (is=0;is<nsmal;is++){
		out << "fSmallConnect [ " << is << "] = \t " << fSmallConnect[is] << endl;
	}
	
	
	out << "fLargeConnect" << endl;
	for (il=0;il<nlarge;il++){
		out << "fLargeConnect [ " << il << "] = \t " << fLargeConnect[il] << endl;
	}
	out << endl;
	
}

void TPZCheckRestraint::Diagnose() {
	
	TPZCheckGeom chkgeo;
	chkgeo.CheckSubFatherTransform(fSmall.Element()->Reference(),fSmall.Side());
	chkgeo.CheckRefinement(fSmall.Element()->Reference()->Father());
	
	
	TPZStack<int> smalldim;
	fSmall.Reference().Element()->LowerDimensionSides(fSmall.Side(),smalldim);
	int cap = smalldim.NElements();
	int s;
	for(s=0; s<cap; s++) {
		TPZCompElSide smalll(fSmall.Element(),smalldim[s]);
		TPZCompElSide largedim = smalll.LowerLevelElementList(1);
		if(!largedim.Exists()) continue;
		TPZCheckRestraint diag(smalll,largedim);
		if(diag.CheckRestraint()) {
			(cout) << "CheckRestraint failed for :" << endl;
		} else {
			(cout) << "CheckRestraint is consistent for : " << endl;
		}
		diag.Print(cout);
		int cind = largedim.Element()->ConnectIndex(largedim.Side());
		TPZCompElSide largelargedim = largedim.LowerLevelElementList(1);
		while(largelargedim.Exists() && LargeConnect(cind) == -1) {
			TPZCheckRestraint diag2(largedim,largelargedim);
			if(diag2.CheckRestraint()) {
				(cout) << "CheckRestraint failed for :" << endl;
			} else {
				(cout) << "CheckRestraint is consistent for : " << endl;
			}
			diag2.Print(cout);
			largedim = largelargedim;
			largelargedim = largelargedim.LowerLevelElementList(1);
			if(largelargedim.Exists()) {
				cind = largelargedim.Element()->ConnectIndex(largedim.Side());
			} else {
				(cout) << "TPZCheckRestraint::Diagnose inconsistent data structure 2\n";
			}
		}
	}
}

