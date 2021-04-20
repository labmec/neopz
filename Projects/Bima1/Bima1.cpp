#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "TPZGenGrid2D.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "pzmatrix.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "pzinterpolationspace.h"
#include "pzsubcmesh.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzfstrmatrix.h"
#include "TPZGenGrid2D.h"
#include "pzbndcond.h"
#include "TPZMaterial.h"
#include "pzelmat.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include "pzhdivpressure.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "TPZInterfaceEl.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

#include "TPZRefPattern.h"
#include "TPZRefPatternDataBase.h"

#include "TPZMatDualHybridPoisson.h"

#ifdef PZ_LOG

static TPZLogger logger("Bima.main");

#endif
static void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL n = 0;
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = pow((REAL)2,0.25 + n/2.)*pow(r,0.5 + n)*cos((0.5 + n)*t);
    u[0] = sol;
    
    du(0,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    
}

static void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFNMatrix<10,STATE> fake(2,1);
    SolExataSteklov(loc,result,fake);
}
static void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFNMatrix<10,REAL> fake(2,1);
    result[0] = loc[0]*loc[0];
}

static void SolExataSteklovSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    u.Resize(1, 0.);
    du.Resize(2, 1);
    du(0,0)=0.;
    du(1,0)=0.;
    
    REAL x = loc[0];
    REAL y = loc[1];
    //const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = 3.363585661014858*pow(pow(x,2) + pow(y,2),1.75)*cos(3.5*t);
    u[0] = sol;
    
    //flux = -k*grad(u), k=1 nesse problema
    du(0,0) = pow(pow(x,2) + pow(y,2),0.75)*(11.772549813552002*x*cos(3.5*t) + 11.772549813552002*y*sin(3.5*t));
    du(1,0) = pow(pow(x,2) + pow(y,2),0.75)*(11.772549813552002*y*cos(3.5*t) - 11.772549813552002*x*sin(3.5*t));
}

static void DirichletSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFNMatrix<10,STATE> fake(2,1);
    SolExataSteklovSuave(loc,result,fake);
}


TPZGeoMesh * MalhaGeo(const int h);
void GroupElements(TPZCompMesh *cmesh);

/** Resolver o problema do tipo
 * -Laplac(u) = 0
 * du/dn = lambda u em todo contorno
 */

using namespace std;

TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder);

TPZInterpolationSpace * FindInterpolationSpace(TPZVec<REAL> &xVec, TPZCompMesh *cmesh, TPZVec<REAL> &qsi);
bool EstaNoQuadrado(const TPZVec<REAL> &xVec, REAL r, REAL tol);
REAL Compute_dudnQuadradoError(int ndiv, TPZCompMesh *cmesh);
REAL Compute_dudn(TPZInterpolationSpace * sp, TPZVec<REAL> &intpoint, TPZVec<REAL> &normal);

#include "pzfmatrix.h"
#include "pzaxestools.h"

bool erronoquadrado = true;
int main()
{
    
    std::ofstream myerrorfile("erros.txt");
	
	for (int porder= 4; porder<5; porder++) {
		
		for(int h=1;h<6;h++){
			
			
			TPZGeoMesh *gmesh = MalhaGeo(h);//malha geometrica
			
			
			TPZCompMesh *cmesh = CreateHybridCompMesh(*gmesh,porder);//malha computacional
            
            std::ofstream out2("gmeshfinal.txt");
            gmesh->Print(out2);
            
            GroupElements(cmesh);
            
			cmesh->LoadReferences();//mapeia para a malha geometrica lo
            
            std::ofstream out("cmesh.txt");
            cmesh->Print(out);
			
			TPZAnalysis analysis(cmesh);
            
            TPZSkylineNSymStructMatrix str(cmesh);
            str.SetNumThreads(8);
            //TPZFStructMatrix str(cmesh);
            
            TPZAutoPointer<TPZMatrix<STATE> > mat = str.Create();
            str.EquationFilter().Reset();
            TPZAutoPointer<TPZMatrix<STATE> > mat2 = mat->Clone();
            
            analysis.SetStructuralMatrix(str);
            TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
            TPZStepSolver<STATE> *gmrs = new TPZStepSolver<STATE>(mat2);
            step->SetReferenceMatrix(mat2);
            step->SetDirect(ELU);
            gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
            TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
            TPZAutoPointer<TPZMatrixSolver<STATE> > autogmres = gmrs;
            analysis.SetSolver(autogmres);
           
            analysis.Run();
            
            TPZVec<std::string> scalnames(1),vecnames(0);
            scalnames[0] = "Solution";
            analysis.DefineGraphMesh(2,scalnames,vecnames,"bima.vtk");
            
            myerrorfile << "\nh = "<< h << " p = " << porder << "\n";
            myerrorfile << "neq = " << cmesh->NEquations() << "\n";
            
            //erro global
            if(erronoquadrado==false){
                analysis.PostProcess(0);
                analysis.SetExact(SolExataSteklovSuave);
                TPZVec<REAL> erros(3);
                analysis.PostProcessError(erros);
                myerrorfile << "H1 = " << erros[0];
                myerrorfile << " L2 = " << erros[1];
                myerrorfile << " semi H1 = " << erros[2] << "\n"<<std::endl;
            }
            else
            {
                REAL errofluxo =  Compute_dudnQuadradoError(h, cmesh);
                myerrorfile<<"Erro semi H1 para fluxo = " << errofluxo<<"\n\n";
                myerrorfile.flush();
            }

				
            cmesh->SetName("Malha depois de Analysis-----");
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
			{
				std::stringstream sout;
				cmesh->Print(sout);
			//	submesh->Block().Print("Block",sout);
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
            delete cmesh;
            delete gmesh;
	
        }
        myerrorfile << "\n"<<std::endl;
    }
	return 0;
}

TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder){
	//TPZCompEl::SetgOrder(porder);
	TPZCompMesh *comp = new TPZCompMesh(&gmesh);
    comp->SetDimModel(gmesh.Dimension());
	
    comp->ApproxSpace().CreateDisconnectedElements(true);
    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
	
	
	// Criar e inserir os materiais na malha
    REAL beta = 6;
	TPZMatDualHybridPoisson *mat = new TPZMatDualHybridPoisson(1,0.,beta);
	TPZMaterial * automat(mat);
	comp->InsertMaterialObject(automat);
	
	
	// Condicoes de contorno
	TPZFMatrix<STATE> val1(1,1,1.),val2(1,1,0.);
	
	TPZMaterial *bnd = automat->CreateBC (automat,-1,2,val1,val2);//misto tbem
        bnd->SetForcingFunction(Dirichlet2,porder);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-2,0,val1,val2);
	TPZBndCond *bndcond = dynamic_cast<TPZBndCond *> (bnd);
        bnd->SetForcingFunction(DirichletSuave,porder);
//	bndcond->SetValFunction(ValFunction);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-3,0,val1,val2);
        bnd->SetForcingFunction(DirichletSuave,porder);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-4,0,val1,val2);
        bnd->SetForcingFunction(DirichletSuave,porder);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-5,0,val1,val2);
        bnd->SetForcingFunction(DirichletSuave,porder);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-6,0,val1,val2);
        bnd->SetForcingFunction(DirichletSuave,porder);
	comp->InsertMaterialObject(bnd);
	
	// Ajuste da estrutura de dados computacional
    comp->ApproxSpace().CreateDisconnectedElements(true);
    comp->ApproxSpace().SetAllCreateFunctionsContinuous();

    std::set<int> matids;
    matids.insert(1);
	comp->AutoBuild(matids);
    
    comp->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    matids.clear();
    for (int i=2; i<=6; i++) {
        matids.insert(-i);
    }
    
    comp->SetDefaultOrder(porder);
	comp->AutoBuild(matids);
    comp->SetDimModel(2);
    
	comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
	comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
    comp->LoadReferences();
    comp->ApproxSpace().CreateInterfaceElements(comp,true);
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		comp->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	

    matids.insert(1);
    for (int i=2; i<=6; i++) {
        matids.insert(-i);
    }

    
    
    comp->ApproxSpace().Hybridize(*comp, matids);
    
    
	
	comp->SetName("Malha Computacional Original");
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		comp->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
    return comp;
}


void ValFunction(TPZVec<REAL> &loc, TPZFMatrix<STATE> &Val1, TPZVec<STATE> &Val2, int &BCType)
{
	BCType = 2;
	Val1.Redim(1, 1);
	Val1(0,0) = 1.;
	Val2[0] = loc[0];
}


TPZGeoMesh * MalhaGeo(const int ndiv){//malha quadrilatera
	
	///malha geometrica
    TPZGeoMesh * gmesh = new TPZGeoMesh();
    
	
    
    
    ///Criando nós
    const int nnodes = 6;
    double coord[nnodes][2] = {{-0.5,0},{0,0},{0.,0.5},{-0.5,0.5},{0.5,0},{0.5,0.5}};
    for(int i = 0; i < nnodes; i++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL,3> nodeCoord(3);
        nodeCoord[0] = coord[i][0];
        nodeCoord[1] = coord[i][1];
        nodeCoord[2] = 0.;
        gmesh->NodeVec()[nodind].Initialize(i,nodeCoord,*gmesh);
    }
    
    ///Criando elementos
    const int nel = 2;
    int els[nel][4] = {{0,1,2,3},{1,4,5,2}};
    for(int iel = 0; iel < nel; iel++){
        TPZManVector<int64_t,4> nodind(4);
        int64_t index;
        nodind[0] = els[iel][0];
        nodind[1] = els[iel][1];
        nodind[2] = els[iel][2];
        nodind[3] = els[iel][3];
        gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    ///Criando elementos de contorno
    const int nelbc = 6;
    int bcels[nelbc][3] = {{0,1,-3},{1,4,-2},{4,5,-4},{5,2,-6},{2,3,-6},{3,0,-5}};
    for(int iel = 0; iel < nelbc; iel++){
        TPZManVector<int64_t,4> nodind(2);
        int64_t index;
        nodind[0] = bcels[iel][0];
        nodind[1] = bcels[iel][1];
        int matid = bcels[iel][2];
        gmesh->CreateGeoElement(EOned,nodind,matid,index);
    }
    
    ///Construindo conectividade da malha
	gmesh->BuildConnectivity();
    
	{
		std::ofstream myfile("geoMinima.txt");
		gmesh->Print(myfile);
	}
	
    ///Refinamento uniforme da malha
	
	///Inicializando padrões de refinamento uniforme
    gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
	for (int i = 0; i < ndiv; i++){
        int nel = gmesh->NElements();
		for (int iel = 0; iel < nel; iel++){
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
			if (!gel) continue;
			if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
            TPZVec<TPZGeoEl*> filhos;
			gel->Divide(filhos);
        }///iel
    }///i
    
    {
		std::ofstream myfile("geoRefinado.txt");
		gmesh->Print(myfile);
	}
    
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif 		 		 
	return gmesh;
}

void GroupElements(TPZCompMesh *cmesh)
{
    cmesh->LoadReferences();
    int nel = cmesh->NElements();
    std::set<TPZCompEl *> celset;
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();
        if (dim ==2) {
            celset.insert(cel);
        }
    }
    std::set<int> elgroupindices;

    for (std::set<TPZCompEl *>::iterator it = celset.begin(); it != celset.end(); it++) {
        
        std::list<TPZCompEl *> group;
        group.push_back(*it);
        TPZCompEl *cel = *it;
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZStack<TPZCompElSide> connected;
            TPZCompElSide celside(cel,is);
            celside.EqualLevelElementList(connected, false, false);
            int neq = connected.NElements();
            for (int eq=0; eq<neq; eq++) {
                TPZCompElSide eqside = connected[eq];
                TPZCompEl *celeq = eqside.Element();
                TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(celeq);
                if (!intface) {
                    continue;
                }
                TPZCompEl *left = intface->LeftElement();
                if (left == cel) {
                    //put in the group
                    group.push_back(intface);
                }
            }
        }
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
                (*it)->Print(sout);
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        int64_t index;
        TPZElementGroup *celgroup = new TPZElementGroup(*cmesh,index);
        elgroupindices.insert(index);
        for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
            celgroup->AddElement(*it);
        }
    }
    cmesh->ComputeNodElCon();
    
    for (std::set<int>::iterator it = elgroupindices.begin(); it!=elgroupindices.end(); it++) {
        TPZCompEl *cel = cmesh->ElementVec()[*it];
        TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel);
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

TPZInterpolationSpace * FindInterpolationSpace(TPZVec<REAL> &xVec, TPZCompMesh *cmesh, TPZVec<REAL> &qsi){
    const int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZCondensedCompEl * cel = dynamic_cast<TPZCondensedCompEl *>(cmesh->ElementVec()[iel]);
        if(!cel) continue;
        if(cel->Dimension()==1) continue;
    
        TPZElementGroup *eg = dynamic_cast<TPZElementGroup *>(cel->ReferenceCompEl());

        TPZInterpolationSpace *sp;
        bool IsInDomain;
        for(int i=0; i<eg->GetElGroup().size(); i++)
        {
            sp = dynamic_cast<TPZInterpolationSpace *>(eg->GetElGroup()[i]);
            if(!sp) continue;
            IsInDomain = sp->Reference()->ComputeXInverse(xVec, qsi, 1e-8);
            if(IsInDomain) return sp;
        }
    }
    
    DebugStop();//nao achei ninguem
    return NULL;
    
}

///Gamma eh a uniao das regioes: 1) x=-0.25 e 0<=y<=0.25; 2) -0.25<=x<=0.25 e y=0.25; 3) x=0.25 e 0<=y<=0.25
//considera-se apenas malha quadrilateral
REAL Compute_dudnQuadradoError(int ndiv, TPZCompMesh *cmesh)
{
    
    REAL error = 0.;
    const REAL rsize = 0.25;
    const int nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    int facesFound = 0;
    for(int iel = 0; iel < nel; iel++)
    {
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[iel]);
        if(!sp) continue;
        
        int dimcel = sp->Dimension();
        if(dimcel != dim-1) continue;
        int matid = sp->Material()->Id();
        if(matid < 0) continue;
      
        TPZGeoEl * gel = sp->Reference();
        TPZManVector<REAL,3> qsi(gel->Dimension()), xCenter(3,0.);
        gel->CenterPoint(gel->NSides()-1, qsi);
        gel->X(qsi,xCenter);
        const REAL tol = 1e-3;
        if(EstaNoQuadrado(xCenter,rsize,tol) == false) continue;
        facesFound++;
        
        TPZManVector<REAL> faceNormal(3,0.);
        if(xCenter[0]==-0.25) faceNormal[0]=-1.;
        if(xCenter[0]== 0.25) faceNormal[0]=1.;
        if(xCenter[1]== 0.25) faceNormal[1]=1.;
        
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 2);
        TPZManVector<int,3> order(gel->Dimension(),intrule->GetMaxOrder());
        intrule->SetOrder(order);
        
        const int npoints = intrule->NPoints();
        //LOOP OVER INTEGRATION POINTS
        for(int ip = 0; ip < npoints; ip++){
            REAL weight, detjac;
            intrule->Point(ip,qsi,weight);
            TPZFNMatrix<9> jacobian(3,3), axes(3,3), jacinv(3,3);
            gel->Jacobian( qsi, jacobian, axes, detjac, jacinv);
            weight *= fabs(detjac);
            const REAL dudnval = Compute_dudn(sp,qsi,faceNormal);
            TPZManVector<REAL> xVec(3,0.);
            gel->X(qsi,xVec);
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<2,STATE> duExato(2,1);
            SolExataSteklovSuave(xVec, uExato, duExato);
            const REAL dudnExato = duExato(0,0)*faceNormal[0]+duExato(1,0)*faceNormal[1];
            error += weight*(dudnval - dudnExato)*(dudnval - dudnExato);
        }///for i
    }///iel
    
    
    int ninterf = 4*pow(2.,ndiv-1);
    //for(int i = 1; i < ndiv-1; i++) n *= 2;
    if(facesFound != ninterf){
        std::cout<<"Opa, numero errado de interfaces"<< std::endl;
        DebugStop();
    }
    
    error = sqrt(error);
    return error;
    
}///method


bool EstaNoQuadrado(const TPZVec<REAL> &xVec, REAL r, REAL tol){
    
    const REAL x = xVec[0];
    const REAL y = xVec[1];
    
    //esta aresta nunca seria bem integrada por causa da singularidade  ///aresta y = 0
    /*tamarindo  if( fabs(y-0.) <= tol && fabs(x) <= r + tol){
     return true;
     }        */
    
    ///aresta y = r
    if( fabs(y-r) <= tol && fabs(x) <= r +tol ){
        return true;
    }
    
    ///aresta x = -r ou x = +r
    if( fabs( fabs(x) -r ) <= tol && y <= r + tol){
        return true;
    }
    
    return false;
    
}///bool

REAL Compute_dudn(TPZInterpolationSpace * sp, TPZVec<REAL> &intpoint, TPZVec<REAL> &normal){
    
    TPZGeoEl *gel = sp->Reference();
    if(gel->Dimension()!=1) DebugStop();
    
    //procurando elemento de interface vizinho
    TPZStack<TPZCompElSide> neighequal;
    neighequal.Resize(0);
    TPZCompElSide celside(sp,gel->NSides()-1);
    celside.EqualLevelElementList(neighequal, 0, 0);
    
    int i, nneighs = neighequal.size();
    if(nneighs!=4) DebugStop();
//    for(i=0; i<nneighs; i++){
//        neighequal[i].Element()->Print();
//    }
    
    TPZFNMatrix<5,STATE> dsoldxL, dsoldxR;
    TPZStack<TPZFNMatrix<5,STATE>,2> dsolVec;
    for(i=0; i<nneighs; i++){
    
        TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement *>(neighequal[i].Element());
        if(!face) continue;
        
        //debug
        TPZManVector<REAL> faceNormal(3);
        face->CenterNormal(faceNormal);
        REAL inner = faceNormal[0]*normal[0] + faceNormal[1]*normal[1];
        if(fabs(inner) < 0.95){
            DebugStop();
        }
        
        //Transformacao do elemento de Lagrange para o elemento de interface
        TPZGeoElSide gels(gel,gel->NSides()-1);
        TPZGeoElSide gelintef = neighequal[i].Reference();
        TPZGeoElSide gelsideinterf(gelintef.Element(),gelintef.Element()->NSides()-1);
        TPZTransform<> tr1 = gels.NeighbourSideTransform(gelintef);
        TPZTransform<> tr2 = gelintef.SideToSideTransform(gelsideinterf);
        TPZTransform<> transfLagranToInterf = tr2.Multiply(tr1);
        TPZManVector<REAL,3> IntPointInterf(2);
        transfLagranToInterf.Apply(intpoint,IntPointInterf);
        
        
        //Transformacao do elemento de interface para o elemento 2D
        TPZCompElSide LeftSide = face->LeftElementSide();
        TPZTransform<> transfLeft;
        face->ComputeSideTransform(LeftSide, transfLeft);
        TPZInterpolationSpace * LeftEl = dynamic_cast<TPZInterpolationSpace*>(LeftSide.Element());
        
        TPZCompElSide RightSide = face->RightElementSide();
        TPZTransform<> transfRight;
        face->ComputeSideTransform(RightSide, transfRight);
        TPZInterpolationSpace * RightEl = dynamic_cast<TPZInterpolationSpace*>(RightSide.Element());
        
        TPZManVector<REAL,3> LeftIntPoint(2), RightIntPoint(2);
        TPZVec<STATE> solL, solR;
        TPZFMatrix<STATE> dsoldaxes(2,1);
        
        //Calculo da solucao
        int count = 0;
        if(LeftSide.Element()->Dimension() > sp->Dimension())
        {
            transfLeft.Apply(IntPointInterf, LeftIntPoint);
            TPZMaterialData dataL;
            LeftEl->InitMaterialData(dataL);
            LeftEl->ComputeShape(LeftIntPoint, dataL);
            LeftEl->ComputeSolution(LeftIntPoint, dataL);
            dsoldaxes(0,0) = dataL.dsol[0][0];
            dsoldaxes(1,0) = dataL.dsol[0][1];
            TPZAxesTools<STATE>::Axes2XYZ(dsoldaxes, dsoldxL, dataL.axes);
            dsolVec.push_back(dsoldxL);
            count++;
            continue;
        }
        
        if(RightSide.Element()->Dimension() > sp->Dimension())
        {
            transfRight.Apply(IntPointInterf, RightIntPoint);
            TPZMaterialData dataR;
            RightEl->InitMaterialData(dataR);
            RightEl->ComputeShape(LeftIntPoint, dataR);
            RightEl->ComputeSolution(LeftIntPoint, dataR);
            dsoldaxes(0,0) = dataR.dsol[0][0];
            dsoldaxes(1,0) = dataR.dsol[0][1];
            TPZAxesTools<STATE>::Axes2XYZ(dsoldaxes, dsoldxR, dataR.axes);
            dsolVec.push_back(dsoldxR);
            count++;
            continue;
        }
        
        if(count==0 || dsolVec.size()!=2){
            DebugStop();
        }
    }
    
    REAL result;
    result = (0.5*(dsolVec[0](0,0)*normal[0]+dsolVec[0](1,0)*normal[1]) + 0.5*(dsolVec[1](0,0)*normal[0]+dsolVec[1](1,0)*normal[1]));
    return result;
}


