//
//  pzgradientreconstruction.cpp
//  PZ
//
//  Created by Agnaldo Farias on 4/10/13.
//
//

#include "pzgradientreconstruction.h"
#include "tpzintpoints.h"
#include "pzintel.h"
#include "pzgnode.h"
#include <cmath>
using namespace std;

TPZGradientReconstruction::TPZGradientReconstruction(int var){
    
    fvar = var;
    fDistortedMesh = false;
    fSlopeLimiter = false;
}

TPZGradientReconstruction::~TPZGradientReconstruction(){
    
}

TPZGradientReconstruction::TPZGradientReconstruction(const TPZGradientReconstruction &copy){
    
    fvar = copy.fvar;
    fDistortedMesh = copy.fDistortedMesh;
    fSlopeLimiter = copy.fSlopeLimiter;
}

TPZGradientReconstruction & TPZGradientReconstruction::operator=(const TPZGradientReconstruction &copy){
    
	fvar  = copy.fvar;
    fDistortedMesh = copy.fDistortedMesh;
     fSlopeLimiter = copy.fSlopeLimiter;
    return *this;
}

void TPZGradientReconstruction::UseWeightCoefficients(){
    
    fDistortedMesh=true;
}

void TPZGradientReconstruction::UseSlopeLimiter(){
    
    fSlopeLimiter=true;
}

void TPZGradientReconstruction::NodeCloserCenterX(TPZCompEl *cel, TPZManVector<REAL> &nodecelX){
    
    if(!cel) DebugStop();
    TPZGeoEl* gel = cel->Reference();
   
    TPZVec<REAL> centerpsi(3,0.), centercelX(3,0.);
    TPZVec<REAL> coord(3,0.);
    TPZVec<REAL> coordX(3,0.);
    nodecelX.Resize(3, 0.);
    
    gel->CenterPoint(gel->NSides()-1,centerpsi);
    centercelX.Fill(0.);
    gel->X(centerpsi,centercelX);
    
    int nnodes = gel->NCornerNodes();
    TPZManVector<REAL> dist(nnodes,0.);
    
    for (int in = 0; in<nnodes; in++) {

        TPZGeoNode geonode = gel->Node(in);
        geonode.GetCoordinates(coord);
        gel->X(coord,coordX);
        
        dist[in]=gel->Distance(centercelX, coordX);
    }
    
    int id = 0;
    for(int j=1; j<dist.size(); j++){
        
        if(dist[id] > dist[j]){
            id = j;
        }
    }
    
    gel->Node(id).GetCoordinates(coord);
    gel->X(coord,nodecelX);
}

void TPZGradientReconstruction::CalcWeights(TPZCompEl *cel, std::set<TPZCompEl *> neighscel, int paramk, TPZManVector<REAL> &weights){
    
    if(paramk<1 || paramk>2) DebugStop();
    
    TPZManVector<REAL> nodecelX;
    NodeCloserCenterX(cel, nodecelX);
    
    int ineighs = -1;
    int nneighs = neighscel.size();
    TPZManVector<REAL> dist(nneighs,0.);
    TPZManVector<REAL> centerpsi(3,0.0), centerneigh(3,0.0);
    
    REAL sum=0.;
    std::set<TPZCompEl *>::iterator it;
    for(it=neighscel.begin(); it!=neighscel.end(); ++it)
    {
        ineighs++;
        TPZGeoEl * gel = (*it)->Reference();
        
        if(!gel) DebugStop();
        
        centerneigh.Fill(0.0);
        gel->CenterPoint(gel->NSides()-1,centerpsi);
        gel->X(centerpsi,centerneigh);
        
        dist[ineighs]=gel->Distance(centerneigh, nodecelX);
        
        sum += 1./pow((REAL)dist[ineighs],paramk);
    }
    
    weights.Resize(dist.size(), 0.);
    REAL temp;
    
    for (int i = 0; i<dist.size(); i++)
    {
        temp = 1./pow((REAL)dist[ineighs],paramk);
        weights[i] = temp/sum;
        if(weights[i]<0 || weights[i]>1) DebugStop();
    }
}

void TPZGradientReconstruction::InsertWeights(TPZVec<REAL> weights, TPZFMatrix<REAL> &DeltaH, TPZFMatrix<REAL> &DifSol){
    
    if(DeltaH.Rows()!=weights.size()) DebugStop();
    if(DifSol.Rows()!=weights.size()) DebugStop();
    
    long ncH = DeltaH.Cols();
    long ncD = DifSol.Cols();
    for (long i = 0; i<weights.size(); i++) {
        
        for (long j = 0; j<ncH; j++){
            DeltaH(i,j) = weights[i]*DeltaH(i,j);
        }
        
        for (long k = 0; k<ncD; k++){
            DifSol(i,k) = weights[i]*DifSol(i,k);
        }
    }
}


void  TPZGradientReconstruction::CalcSlopeLimiter(TPZCompEl *cel, REAL solcel, TPZStack<TPZCompEl *> neighscell, TPZStack<REAL> solneighscell, REAL &alpha){
    
    
    long nels = neighscell.size();
    if(nels<1) DebugStop();
    
    TPZCompEl *celneigh_max;
    TPZCompEl *celneigh_min;
    REAL solKmax, solKmin;//K is the cell 
    
    celneigh_min = neighscell[0];
    solKmin = solneighscell[0];
    celneigh_max = neighscell[0];
    solKmax = solneighscell[0];
    for(long i = 1; i<nels; i++)
    {
        if(solKmin > solneighscell[i])
        {
            solKmin = solneighscell[i];
            celneigh_min = neighscell[i];
            continue;
        }
        
        if(solKmax < solneighscell[i])
        {
            solKmax = solneighscell[i];
            celneigh_max = neighscell[i];
        }
    }
    
    if(solKmax < solKmin) DebugStop();
    
    int IntOrder =  cel->GetgOrder();
    REAL solKside;
    REAL temp;
    TPZStack<REAL> alphavec;
    
    int nsides = cel->Reference()->NSides();
    for (int is=0; is<nsides-1; is++)
    {
        TPZCompElSide celside(cel,is);
        solKside = MeanCell(celside.Element(),IntOrder);
        
        if(solKside > solKmax)
        {
            temp = (solKmax - solcel)/(solKside-solcel);
        }
        else if(solKside < solKmin)
        {
            temp = (solKmin - solcel)/(solKside-solcel);
        }
        else{
            temp = 1.;
        }
        
        if(temp < 0. || temp > 1.) DebugStop();
        alphavec.Push(temp);
    }
    
    alpha = alphavec[0];
    for (long j=1; j<alphavec.size(); j++){
        
        if(alpha > alphavec[j])
        {
            alpha = alphavec[j];
        }
    }
}

REAL TPZGradientReconstruction::MeanCell(TPZCompEl *cel,int IntOrder)
{
	TPZIntPoints *pointIntRule = ((TPZInterpolatedElement*)cel)->Reference()->CreateSideIntegrationRule((cel->Reference()->NSides())-1,IntOrder);
	int it, npoints = pointIntRule->NPoints();
	REAL integral = 0.0;
	TPZManVector<REAL> point(3,0.);
	TPZManVector<REAL> xpoint(3,0.);
	REAL weight;
	for(it=0;it<npoints;it++) {
		pointIntRule->Point(it,point,weight);
		weight /= cel->Reference()->RefElVolume();
		cel->Reference()->X(point,xpoint);
        
        TPZVec<STATE> sol;
        cel->Solution(xpoint, fvar, sol);
        if(sol.size()!=1) {
            PZError << "TPZGradientReconstruction::The number of solutions variable can not be other than 1.\n";
            DebugStop();
        }
#ifdef STATE_COMPLEX
        integral += weight*sol[0].real();
#else
        integral += weight*sol[0];
#endif
    }
	return integral;
}

void TPZGradientReconstruction::GradientReconstructionByLeastSquares(TPZCompEl *cel,TPZManVector<REAL,3> &center,TPZVec<REAL> &Grad) {
   
    TPZFMatrix<REAL> grad;
    int dim;
    dim = cel->Mesh()->Dimension();
    
    if(!cel || cel->Dimension()!=dim){
        
        PZError << "TPZGradientReconstruction:: Element has dimension equal to the dimension of the problem.\n";
        DebugStop();
    }
	REAL solalfa;
	REAL solbeta;
    
    int k, side;
	// Integration order to calculate cell mean of solution
	int intOrder = cel->GetgOrder();
    
	TPZStack<TPZCompElSide> neighs;
	int nneighs;
	
    center.Resize(3, 0.);
	TPZManVector<REAL> centerpsi(3,0.0), centerbeta(3,0.0);
    
    TPZManVector<REAL> nodecell(3,0.0);
	
	TPZFMatrix<REAL> A(dim,dim);  // Linear System matrix
	grad.Redim(dim,1);
	
	//matrices to apply the method of least squares
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	
    // Finding the center of the current element
    TPZGeoEl* gelalfa = cel->Reference();
    gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
    center.Fill(0.);
    gelalfa->X(centerpsi,center);
	
	solalfa = MeanCell(cel,intOrder);
    
    neighs.Resize(0);
    
    // Looking for all the neighbors of the element.
    for(side = 0; side <cel->Reference()->NSides()-1; side++)
    {
        TPZCompElSide celside(cel,side);
        celside.ConnectedElementList(neighs,0,0);
    }
    
    nneighs = neighs.NElements();
    if(!nneighs){
        
        PZError << "TPZGradientReconstruction::Element has no neighbors.\n";
        DebugStop();
    }
    
    //Set of neighbors
    std::set<TPZCompEl *> neighscel;
    for(int i =0; i<nneighs; i++)
    {
        TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(neighs[i].Element());
        if(!InterpEl || InterpEl->Dimension()!=dim) continue;
        neighscel.insert(neighs[i].Element());
    }
    
    nneighs=neighscel.size();
    
    // If there are neighbors, we make the process of least squares to calculate an approximation of the gradient.
    // For each neighbor calculate the DeltaH (from center to center of the current element)
    // the value of the solution at its center solbeta
    DeltaH.Redim(nneighs,dim);
    DeltaHTranspose.Redim(dim,nneighs);
    DifSol.Redim(nneighs,1);
    
    //cel neighs max and min
    TPZStack<TPZCompEl *> neighscell;
    TPZStack<REAL> solneighscell;
    
    // Mounting deltaH matrix and difference of the solutions DifSol 
    int ineighs=-1;
    long counter=0;
    std::set<TPZCompEl *>::iterator it;
    for(it=neighscel.begin(); it!=neighscel.end(); ++it)
    {
        ineighs++;
        TPZGeoEl * gelbeta = (*it)->Reference();
        
        if(!gelbeta) DebugStop();
        
        centerpsi.Fill(0.0);
        centerbeta.Fill(0.0);
        gelbeta->CenterPoint(gelbeta->NSides()-1,centerpsi);
        gelbeta->X(centerpsi,centerbeta);
        solbeta = MeanCell(gelbeta->Reference(),intOrder);
        
        for(k=0;k<dim;k++)
        {
            DeltaH(ineighs,k) = centerbeta[k] - center[k];
        }
        DifSol(ineighs,0) = solbeta - solalfa;
        counter ++;
        
        if(fSlopeLimiter==true)
        {
            neighscell.Push((*it));
            solneighscell.Push(solbeta);
        }
    }
    
    //insert weight
    if(fDistortedMesh==true)
    {
        TPZManVector<REAL> weights;
        CalcWeights(cel, neighscel, 1, weights);
        InsertWeights(weights, DeltaH, DifSol);
    }
    
    // Solving the system by least squares: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    grad = DeltaHTranspose*DifSol;
    A = DeltaHTranspose*DeltaH;
    if(counter > 0){
        A.SolveDirect(grad,ELU);
    }
	
	// Return gradient vector
	Grad.Resize(dim);
    REAL alphaK = 1.;
    if(fSlopeLimiter==true)
    {
        CalcSlopeLimiter(cel,solalfa, neighscell, solneighscell, alphaK);
    }
	for(int j=0;j<dim;j++)
		Grad[j] = (1./alphaK)*grad(j,0);
}


void TPZGradientReconstruction::GetDataGradient(TPZCompMesh *cmesh,TPZFMatrix<REAL> &datagradients){
    
    // Resizing the array data reconstruction gradients.
    int dim  = cmesh->Dimension();
    long nelem = cmesh->NElements();
    datagradients.Redim(nelem,2*dim+2);
	
    TPZManVector<REAL,3> center;
    TPZManVector<REAL> Grad(dim);
    
	long i, k;
	long counter = 0;
	
	
    TPZCompEl *cel;
    
    // Calculating the gradient in the computational element.
    for(i=0;i<nelem;i++) {
        
        cel = cmesh->ElementVec()[i];
        
        if(!cel || cel->Dimension()!=dim) continue;
		center.Fill(0.0);
		Grad.Fill(0.0);
        
        GradientReconstructionByLeastSquares(cel, center, Grad);
		
        //data of the vector gradiente
        for(k=0;k<dim;k++){
            datagradients(counter,k) = center[k];//center of the element
            datagradients(counter,dim+k) = Grad[k];//gradient value
        }
		
		datagradients(counter,2*dim) = cel->Index(); //element index
		datagradients(counter,2*dim+1) = cel->VolumeOfEl();//// Volume of the element
        
		counter++;
	}
    
    // Resizing the matrix of gradients
	k = datagradients.Cols();
    datagradients.Resize(counter,k);
}