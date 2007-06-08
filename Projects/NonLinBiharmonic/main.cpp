/**
 * @param type: 1 = quadrados; 2 = triangulos
 * @param resolution: nmero de refinamentos uniformes
 * @param porder: ordem de interpolacao
 * @since  Fev 01, 2005
 * @author Paulo Rafael B�ing
 */
#include <iostream>
#include <time.h>
#include <pzgmesh.h>
#include <pzcmesh.h>
#include <pzcompel.h>
#include <pzgeoel.h>
#include <pzvec.h>
#include <pzgnode.h>
#include <pzgeoelbc.h>
//#include <pzpoisson3d.h>
#include <pzmatrix.h>
#include <TPZCompElDisc.h>
#include <pzfstrmatrix.h>
#include <pzstepsolver.h>
#include <fstream>
#include <pzbstrmatrix.h>
#include <pzblockdiag.h>
#include <pzbdstrmatrix.h>
#include <pzbndmat.h>
#include <pznonlinbiharmonic.h>
#include "TPZGeoElement.h"
#include "TPZInterfaceEl.h"
#include <pzanalysiserror.h>
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "TPZShapeDisc.h"
#include "pzvec.h"
#include "pzquad.h"
#include "pzbndcond.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pznonlinanalysis.h"
#include "checkconv.h"
using namespace pzshape;
using namespace pzrefine;
using namespace pzgeom;
using namespace std;
static REAL Pi, lbd, Re;

TPZCompMesh * CreateMesh(int , int , int );

void SolveLU(TPZNonLinearAnalysis &, TPZCompMesh *, TPZGeoMesh *, ofstream &);

void SolveIterative(TPZNonLinearAnalysis &, TPZCompMesh *, TPZGeoMesh *);

void SolveIterative_2(TPZNonLinearAnalysis &, TPZCompMesh *, TPZGeoMesh *);

void SolveIterative_GMRES(TPZNonLinearAnalysis &, TPZCompMesh *, TPZGeoMesh *, ofstream &);

void SolveFrontal(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha);

void AnalysisErrorandOrder(void);
void AnalyseElementSolutionClass6(TPZCompElDisc *intel,TPZVec<REAL> &pont_ksi);
TPZCompElDisc * FindElement(TPZCompMesh *malhacomp, TPZVec<REAL> &pont,TPZVec<REAL> &pont_ksi);

void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] =  0.5*exp(lbd*x[0])*sin(2.*Pi*x[1])*( lbd - 2.*Pi)*(lbd + 2.*Pi)*
           (-lbd*lbd + lbd*Re + 4.*Pi*Pi)/(Re*Pi);
}

void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] =        x[1]-(exp(lbd*x[0])*sin(2.0* Pi*x[1]))/(Pi*2.0);
  deriv(0,0) = -lbd*exp(lbd*x[0])*sin(2.0* Pi*x[1])/(Pi*2.0); //dx
  deriv(1,0) =  1.0 - exp(lbd*x[0])*cos(2.0* Pi*x[1]);   // dy
  deriv(2,0) = -lbd*lbd *exp(lbd*x[0])*sin(2.0* Pi*x[1])/(Pi*2.0)
               + 2.0*exp(lbd*x[0])*sin(2.0* Pi*x[1])*Pi;   // laplaciano
  deriv(3,0) = -lbd*lbd*lbd *exp(lbd*x[0])*sin(2.0* Pi*x[1])/(Pi*2.0)
               + lbd*2.0*exp(lbd*x[0])*sin(2.0* Pi*x[1])*Pi;//D/dx do laplaciano
  deriv(4,0) = -lbd*lbd *exp(lbd*x[0])*cos(2.0* Pi*x[1])
               + 2.0*exp(lbd*x[0])*cos(2.0* Pi*x[1])*Pi*2.0*Pi;//D/dy do laplaciano
  deriv(5,0) = -lbd*lbd *exp(lbd*x[0])*sin(2.0* Pi*x[1])/(Pi*2.0);// dxx
  deriv(6,0) = 2.0*exp(lbd*x[0])*sin(2.0* Pi*x[1])*Pi;// dyy
  deriv(7,0) = -lbd *exp(lbd*x[0])*cos(2.0* Pi*x[1]);// dxy

}

void CC1(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = 0.0;  // Cond. de Cont. de Dirichlet
  f[1] = -1. + exp(lbd*x[0]);  // Cond. de Cont. de Neumann
}


void CC2(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] =  x[1]-exp(lbd*1.5)*sin(2.0* Pi*x[1])/(Pi*2.0);
  f[1] = -lbd*exp(lbd*1.5)*sin(2.0* Pi*x[1])/(Pi*2.0);
}
void CC3(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = 2.0;  // Cond. de Cont. de Dirichlet
  f[1] = 1. - exp(lbd*x[0]);  // Cond. de Cont. de Neumann
}


void CC4(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = x[1]-exp(-lbd*0.5)*sin(2.0* Pi*x[1])/(Pi*2.0);
  f[1] = lbd*exp(-lbd*0.5)*sin(2.0*Pi*x[1])/(Pi*2.0);
}

int main(){

//  TPZInterfaceElement::SetCalcStiffPenalty();

  TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);

  Pi = 4.0*atan(1.0);

  TPZNonLinBiharmonic::NorP = 1;

  ofstream out("teste.txt");
  ofstream nonlin("AnalysisNL.out");

  time_t tm1 = time(0);
  int p_inicial=3, p_final=3, n_inicial=2,  n_final=2;

  int Reynolds[]={1};  // Reynolds e NReynolds precisam ter a mesma dimensao.
  TPZVec<int> NReynolds(1,0);                       // Isso e manual
  for(int ip=0; ip<NReynolds.NElements(); ip++) NReynolds[ip] = Reynolds[ip];

  nonlin<<" Kovasznay Flow "<<endl;
  nonlin<<" The current time is :"<<ctime(&tm1)<<endl;

  out<<n_final-n_inicial<<"  "<<p_final-p_inicial<<"  "<<NReynolds.NElements()<<endl;

  for(int iR=0; iR<NReynolds.NElements(); iR++){

    out<<NReynolds[iR]<<endl;

    Re = TPZNonLinBiharmonic::Re = NReynolds[iR];
    cout<<" Reynolds = "<<Re<<endl;
    lbd = Re/2. - sqrt(Re*Re/4. + 4.*Pi*Pi);
    cout<<" lbd = "<<lbd<<endl;
    time_t tm3=0, tm2=0 ;
    clock_t ck3, ck2;

    for(int p_ordem=p_inicial; p_ordem<=p_final; p_ordem++){
      out<<p_ordem;
      for(int n_refin=n_inicial; n_refin<=n_final; n_refin++){

      tm2 = time(0);
      ck2 = clock();
      if (tm3){
        nonlin<<" Wall time = "<<difftime(tm2,tm3)/60.<<"minutes"<<endl;
	nonlin<<" CPU  time = "<<(double(ck2-ck3)/CLOCKS_PER_SEC)/60.<<"minutes"<<endl;
      }
      tm3 = time(0);
      ck3 = clock();

      nonlin<<" _________________________________________________________________________"<<endl;
      nonlin<<"Dados: "<<endl;
      nonlin<<"    Metodo               = "<<TPZNonLinBiharmonic::NorP<<" (1=Newton, 0=Picard)"<<endl;
      nonlin<<"    Reynolds             = "<<NReynolds[iR]<<endl;
      nonlin<<"    Ordem p              = "<<p_ordem<<endl;
      nonlin<<"    Nivel de Refinamento = "<<n_refin<<endl<<endl;



      TPZGeoMesh *geomalha = NULL;
      TPZCompMesh *malha = NULL;

      malha = CreateMesh(1, n_refin, p_ordem);

      geomalha = malha->Reference();
      TPZNonLinearAnalysis an(malha, nonlin);

      //SolveIterative(an, malha, geomalha);
      //SolveIterative_2(an, malha, geomalha);
      //SolveIterative_GMRES(an, malha, geomalha, out);


      SolveLU(an, malha, geomalha, nonlin);
      an.PostProcess(1);

      ofstream analy("Analysis.out");
      an.Print("Analysis.out", analy);

      ofstream geofile("GeoMalha.out");
      geomalha->Print(geofile);

      ofstream cfile("CompMalha.out");
      malha->Print(cfile);

      ofstream file("Solution.out");
      TPZFMatrix toprint = an.Solution();
      toprint.Print("solution", file);

      TPZVec<char *> scalar_names(1);
      TPZVec<char *> vec_names(1);
      scalar_names[0] = "Solution";
      vec_names[0] = "Derivate";
      an.DefineGraphMesh(2, scalar_names, vec_names, "Solution.dx");

      an.PostProcess(2);

      an.SetExact(ExactSolution);
      TPZVec<REAL> pos;
      out<<endl;
      an.PostProcess(pos,out);

      TPZVec<REAL> pt_x(3,0.0), pt_ksi(2,0);
      pt_x[0] = 0.5;
      pt_x[1] = 0.5;
      pt_x[2] = 0.0;
      TPZCompElDisc *myElem = NULL;
      myElem = FindElement(malha, pt_x, pt_ksi);
      if(myElem) {
        cout<<"myElem != 0"<<endl;
        AnalyseElementSolutionClass6(myElem, pt_ksi);
       }

      delete malha;
      delete geomalha;


   }//nivel
  }//pordem
 }// Reynolds

 AnalysisErrorandOrder();

}//end of main


TPZCompElDisc * FindElement(TPZCompMesh *malhacomp, TPZVec<REAL> &pont,TPZVec<REAL> &pont_ksi){
 cout<<endl<<"   ---------  "<<pont_ksi<<endl;
 cout<<endl<<"   =========  "<<pont<<endl;
 for(int i = 0; i< malhacomp->ElementVec().NElements() ;i++) {
     TPZCompEl* mycompel = malhacomp->ElementVec()[i];//Divide(0,subelindex,0);
     if(!mycompel) continue;
     TPZGeoEl * mygeoele = mycompel->Reference();
     if((mygeoele->WhichSide(pont)!=-1)) { //
         mygeoele->ComputeXInverse(pont,pont_ksi);
	 cout<<endl<<"   ---------  "<<pont_ksi<<endl;
	 cout<<endl<<"   =========  "<<pont<<endl;
	 //cout<<"ACHEI"<<endl;
         return dynamic_cast<TPZCompElDisc*> (mycompel);
	 //return dynamic_cast<TPZCompElDiscTPZInterpolatedElement *> (mycompel);
     }

 }
   return 0;
  }

void AnalyseElementSolutionClass6(TPZCompElDisc *intel,TPZVec<REAL> &pont_ksi)
{
  cout<<endl<<"   ---------"<<pont_ksi<<endl;
  //int dim = intel->Reference()->Dimension();

   int variable;
  TPZAutoPointer<TPZMaterial> mat = intel->Material();
  if(!mat) return;
//  char varname[] = "Solution";
//  char varname[] = "Derivate";
  variable = mat->VariableIndex("Solution");
  if(variable == -1)
  {
    cout << __PRETTY_FUNCTION__ << " Variable index for Solution not found" << endl;
    return;
  }
  int numvar = mat->NSolutionVariables(variable);
  TPZManVector<REAL> solution(numvar,0.);

  intel->Solution(pont_ksi,variable,solution);
  cout << "The solution is " << solution << endl;

}



void AnalysisErrorandOrder() {

 ifstream infile("teste.txt", ios_base::in);
 if(!infile)
   cout<<" Can not open file: teste.txt, for read"<<endl;

 ofstream outf("AnalysisErrorandOrder.tex", ios_base::out);
 outf.setf(ios_base::scientific, ios_base::floatfield);
 outf.precision(4);
 outf.width(6);

 int NNiveis, Np, NReynolds;
 infile >> NNiveis;
 infile >> Np;
 infile >> NReynolds;
 int tmp_p, tmp_Re;

 TPZVec<REAL> error_atual(5,0.), error_anter(5,0.), ordem(5,0.);


 for(int iR=0; iR<NReynolds; iR++){
   infile >> tmp_Re;
     for(int ip=0; ip<=Np; ip++){
     infile >> tmp_p;
     for(int iN=0; iN<=NNiveis; iN++){
       if(iN==0){
         for(int j=0; j<5; j++){
           infile >> error_atual[j];
         }
	 outf<<"\\begin{table}[ht]"<<endl;
         outf<<"\\centering \\caption{Re="<<tmp_Re<<"  p = "<<tmp_p<<"} \\label{tab:Re"<<tmp_Re<<"p"<<tmp_p<<"}"<<endl;
         outf<<"\\renewcommand\\arraystretch{1.5}"<<endl;
         outf<<"\\vspace{20pt}"<<endl;
         outf<<"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}"<<endl;
         outf<<"\\hline"<<endl;
         outf<<"  & \\multicolumn{2}{c|}{L2} & \\multicolumn{2}{c|}{H1} &\\multicolumn{2}{c|}{H2}&\\multicolumn{2}{c|}{U}  \\\\"<<endl;
         outf<<"\\hline"<<endl;
         outf<<" L & Error & Order & Error & Order & Error & Order & Error & Order \\\\"<<endl;
         outf<<"\\hline"<<endl;
         outf.setf(ios_base::scientific, ios_base::floatfield);
         outf<<iN<<" & "<<error_atual[0]<<" &        & "<<error_atual[2]<<" &        &"<<error_atual[4]<<" &        &"<<error_atual[1]<<" &         \\\\"<<endl;
	 outf<<"\\hline"<<endl;
       }
       else{
         for(int j=0; j<5; j++){
	   error_anter[j] = error_atual[j];
	   infile >> error_atual[j];
	   ordem[j] = log10(error_atual[j]/error_anter[j])/log10(0.5);
         }
	 outf.setf(ios_base::scientific, ios_base::floatfield);
	 outf<<iN<<" & "<<error_atual[0]<<" & ";
	 outf.setf(ios_base::floatfield, ios_base::floatfield);
	 //outf.width(4);
 	 outf<<ordem[0]<< " & ";
	 outf.setf(ios_base::scientific, ios_base::floatfield);
	 outf<<error_atual[2]<<" & ";
	 outf.setf(ios_base::floatfield, ios_base::floatfield);
	 outf<<ordem[2]<< " & ";
	 outf.setf(ios_base::scientific, ios_base::floatfield);
	 outf<<error_atual[4]<<" & ";
	 outf.setf(ios_base::floatfield, ios_base::floatfield);
	 outf<<ordem[4]<< " & ";
	 outf.setf(ios_base::scientific, ios_base::floatfield);
	 outf<<error_atual[1]<<" & ";
	 outf.setf(ios_base::floatfield, ios_base::floatfield);
	 outf<<ordem[1]<< "     \\\\"<<endl;
	 outf<<"\\hline"<<endl;
       }
       if(iN==NNiveis){
         outf<<"\\end{tabular}"<<endl;
	 outf<<"\\end{table}"<<endl<<endl<<endl;
       }
     }
   }
 }


}


void SolveLU(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha, ofstream &iterproc){

/*
  cout << " TPZSpStructMatrix " << endl;
  an.SetStructuralMatrix(mat);
  cout << endl;

  TPZMatrix *sec_mat = NULL;


  TPZStepSolver solv;
  solv.SetGMRES(3000,200,Pre,1e-9,0);
  cout << "GMRES " << endl;
  solv.SetMatrix(sec_mat);
  an.SetSolver(solv);

  an.Solution().Redim(0,0);
  cout << "Assemble " << endl;

 */


  TPZBandStructMatrix mat(malha);
  //TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);
  //TPZParFrontStructMatrix mat(malha);
  //cout << "BandStructMatrix " << endl;
  an.SetStructuralMatrix(mat);
  TPZMatrix *sec_mat = NULL;

  TPZStepSolver solv;
  solv.SetMatrix(sec_mat);
  solv.SetDirect(ELU);

  cout << "ELU " << endl;
  an.SetSolver(solv);

  cout << endl;
  an.Solution().Redim(0,0);
  cout << "Assemble " << endl;


  an.IterativeProcess(iterproc,1e-8,60);
  //an.Assemble();   an.Run() faz a Assemble e Solve
  //an.Solve();
  cout << endl;
  //TPZStepSolver solv2;
  // TPZSpStructMatrix sparse(malha);
  //an.SetStructuralMatrix(sparse);
  //solv2.SetGMRES( 2000, 10, solv, 1.e-16, 0);
  //an.SetSolver(solv2);
  //an.Assemble();
  //cout << "Solve " << endl;
  //an.Solve();

  TPZFBMatrix * bandmat = dynamic_cast<TPZFBMatrix*>(solv.Matrix().operator ->());
  if (bandmat){
    cout << endl << "Banda = " << bandmat->GetBand() << endl;
    cout << "No equacoes = " << malha->NEquations() << endl;
    iterproc << "No equacoes = " << malha->NEquations() << endl;

//    cout << "No equacoes 2  TPZBandStructMatrix mat(malha)";
  //TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);
  //TPZParFrontStructMatrix mat(malha);

    //     bandmat->Print();

  }

}

void SolveIterative(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha){

  TPZBandStructMatrix mat(malha);
  //  TPZFStructMatrix mat(malha);
  an.SetStructuralMatrix(mat);
  cout << "BandStructMatrix " << endl;

  TPZBlockDiagonalStructMatrix strblock(malha);
  TPZBlockDiagonal * block = new TPZBlockDiagonal();
  strblock.AssembleBlockDiagonal(*block);

  TPZStepSolver solv, precond(block);

  precond.SetDirect(ELU);
  //TPZBandStructMatrix mat(malha);
  //TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);
  //TPZParFrontStructMatrix mat(malha);
  cout << "BandStructMatrix " << endl;

  cout << " GMRES " << endl;

  solv.SetGMRES( 2000, 10, precond, 1e-10, 0);

  an.SetSolver(solv);
  an.SetStructuralMatrix(mat);
  cout << endl;
  cout << "Assemble " << endl;

  //an.Assemble();  // an.Run() faz a Assemble e Solve
  cout << endl;
  cout << "Solve " << endl;

  //an.Solve();

  TPZFBMatrix * bandmat = dynamic_cast<TPZFBMatrix*>(solv.Matrix().operator ->());
  if (bandmat){
    cout << endl << "Banda = " << bandmat->GetBand() << endl;
    cout << "No equacoes = " << malha->NEquations() << endl;
    cout << "No equacoes 2  = " << bandmat->Rows() << endl;

    //     bandmat->Print();
  }

}

void SolveIterative_2(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha){

  TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);

  cout << " TPZParFrontStructMatrix " << endl;
  an.SetStructuralMatrix(mat);
  cout << endl;

  TPZMatrix *sec_mat = NULL;

  TPZStepSolver solv;
  solv.SetDirect(ELU);
  cout << "ELU " << endl;
  solv.SetMatrix(sec_mat);
  an.SetSolver(solv);

  an.Solution().Redim(0,0);
  cout << "Assemble " << endl;

  ofstream iterproc("IteraProc.dat");


  an.IterativeProcess(iterproc,1e-10,50);
  cout << endl;

  TPZFBMatrix * bandmat = dynamic_cast<TPZFBMatrix*>(solv.Matrix().operator ->());
  if (bandmat){
    cout << endl << "Banda = " << bandmat->GetBand() << endl;
    cout << "No equacoes = " << malha->NEquations() << endl;
    cout << "No equacoes 2  TPZBandStructMatrix mat(malha)";
  }

}

void SolveIterative_GMRES(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha, ofstream &iterproc){


  TPZSpStructMatrix mat(malha);

  cout << " TPZSpStructMatrix " << endl;
  an.SetStructuralMatrix(mat);
  cout << endl;

  TPZMatrix *sec_mat = NULL;

  TPZStepSolver Pre;

  TPZBlockDiagonalStructMatrix strBlockDiag(malha);

  TPZBlockDiagonal * block = new TPZBlockDiagonal();
  strBlockDiag.AssembleBlockDiagonal(*block);

  Pre.SetMatrix(block);
  Pre.SetDirect(ELU);

  TPZStepSolver solv;
  solv.SetGMRES(3000,200,Pre,1e-8,0);
  cout << "GMRES " << endl;
  solv.SetMatrix(sec_mat);
  an.SetSolver(solv);

  an.Solution().Redim(0,0);
  cout << "Assemble " << endl;

//  ofstream iterproc("IteraProc.dat");
//  for(int j1=1; j1<6; j1++){
//    TPZNonLinBiharmonic::Re = j1*400;
    an.IterativeProcess(iterproc,1e-10,50);
//  }
   cout << endl;

    cout << "No equacoes = " << malha->NEquations() << endl;
    iterproc<< "No equacoes = " << malha->NEquations() << endl;





}


TPZCompMesh * CreateMesh(int type, int n_refin, int p_ordem){

  TPZCompEl::SetgOrder(p_ordem);

  TPZGeoMesh *geomalha = new TPZGeoMesh;

  TPZVec<REAL> coord(2);

 int nodind;

  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = -0.5;
  coord[1] = 0.;
  geomalha->NodeVec()[nodind].Initialize(0, coord, *geomalha);

  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = 1.5;
  coord[1] = 0.;
  geomalha->NodeVec()[nodind].Initialize(1, coord, *geomalha);


  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = -0.5;
  coord[1] = 2.0;
  geomalha->NodeVec()[nodind].Initialize(2, coord, *geomalha);


  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = 1.5;
  coord[1] = 2.0;
  geomalha->NodeVec()[nodind].Initialize(3, coord, *geomalha);

  TPZVec<TPZGeoEl *> geoel(1);
  int index;
  TPZVec<int> incid(4);
  incid[0] = 0;
  incid[1] = 1;
  incid[2] = 3;
  incid[3] = 2;
  geoel[0] = geomalha->CreateGeoElement(EQuadrilateral, incid, 1, index);


  geomalha->BuildConnectivity();

  //   Refinar malha
  TPZVec<TPZGeoEl *> children, netos, bisnetos,bisbisnetos,bisbisbisnetos, trinetos;


  if(n_refin==1)
    for(int i=0;i<1;i++){
       geoel[i]->Divide(children);
    }
  else if(n_refin==2)
    for(int i=0;i<1;i++){
       geoel[i]->Divide(children);
         for(int j = 0; j < children.NElements(); j++){
            children[j]->Divide(netos);
        }
   }
  else if(n_refin==3)
   for(int i=0;i<1;i++){
       geoel[i]->Divide(children);
         for(int j = 0; j < children.NElements(); j++){
            children[j]->Divide(netos);
                  for(int j2 = 0; j2< netos.NElements(); j2++) {
                	 netos[j2]->Divide(bisnetos);
         	  }
	}

  }
  else if(n_refin==4)
    for(int i=0;i<1;i++){
       geoel[i]->Divide(children);
         for(int j = 0; j < children.NElements(); j++){
            children[j]->Divide(netos);
                  for(int j2 = 0; j2< netos.NElements(); j2++) {
                	 netos[j2]->Divide(bisnetos);
                          for(int j3 = 0; j3< bisnetos.NElements(); j3++) {
               	             bisnetos[j3]->Divide(bisbisnetos);
         		  }
		 }

       }
   }
  else if(n_refin==5)
   for(int i=0;i<1;i++){
       geoel[i]->Divide(children);
         for(int j = 0; j < children.NElements(); j++){
            children[j]->Divide(netos);
                  for(int j2 = 0; j2< netos.NElements(); j2++) {
                	 netos[j2]->Divide(bisnetos);
                          for(int j3 = 0; j3< bisnetos.NElements(); j3++) {
               	            bisnetos[j3]->Divide(bisbisnetos);
    			     for(int j4 = 0; j4< bisbisnetos.NElements(); j4++) {
               	               bisbisnetos[j4]->Divide(bisbisbisnetos);
    			     }
         		 }
		  }
      }
   }
 else if(n_refin==6)
   for(int i=0;i<1;i++){
       geoel[i]->Divide(children);
         for(int j = 0; j < children.NElements(); j++){
            children[j]->Divide(netos);
                  for(int j2 = 0; j2< netos.NElements(); j2++) {
                	 netos[j2]->Divide(bisnetos);
                          for(int j3 = 0; j3< bisnetos.NElements(); j3++) {
               	            bisnetos[j3]->Divide(bisbisnetos);
    			     for(int j4 = 0; j4< bisbisnetos.NElements(); j4++) {
               	               bisbisnetos[j4]->Divide(bisbisbisnetos);
       			         for(int j5 = 0; j5< bisbisbisnetos.NElements(); j5++) {
               	                    bisbisbisnetos[j5]->Divide(trinetos);
    			        }
    			     }
         		 }
		  }
      }
   }
   else {

      cout<<" 1 <= n_refin <= 5 "<<endl;
   }

  TPZGeoElBC elbc1(geoel[0],4,-1,*geomalha);
  TPZGeoElBC elbc2(geoel[0],5,-2,*geomalha);

  TPZGeoElBC elbc3(geoel[0],6,-3,*geomalha);
  TPZGeoElBC elbc4(geoel[0],7,-4,*geomalha);

//   TPZGeoElBC elbc5(geoel[2],5,-2,*geomalha);
//   TPZGeoElBC elbc6(geoel[2],6,-3,*geomalha);

//   TPZGeoElBC elbc7(geoel[3],6,-3,*geomalha);
//   TPZGeoElBC elbc8(geoel[3],7,-2,*geomalha);

  TPZCompMesh * malha = new TPZCompMesh(geomalha);
  malha->SetDimModel(2);

  TPZAutoPointer<TPZMaterial> mater;
  mater = new TPZNonLinBiharmonic(1,0.);  // segundo par. �a f(x)
                                   // primeiro par. �o material

  mater->SetForcingFunction(Forcing1);
  malha->InsertMaterialObject(mater);

 TPZFMatrix val1(1,1,0.), val2(2,1,0.);

 TPZAutoPointer<TPZMaterial> bnd1 = mater->CreateBC(mater,-1,0, val1, val2);
 TPZAutoPointer<TPZMaterial> bnd2 = mater->CreateBC(mater,-2,0, val1, val2);
 TPZAutoPointer<TPZMaterial> bnd3 = mater->CreateBC(mater,-3,0, val1, val2);
 TPZAutoPointer<TPZMaterial> bnd4 = mater->CreateBC(mater,-4,0, val1, val2);

 bnd1->SetForcingFunction(CC1);
 malha->InsertMaterialObject(bnd1);
 bnd2->SetForcingFunction(CC2);
 malha->InsertMaterialObject(bnd2);
 bnd3->SetForcingFunction(CC3);
 malha->InsertMaterialObject(bnd3);
 bnd4->SetForcingFunction(CC4);
 malha->InsertMaterialObject(bnd4);

  /*
  TPZBndCond *cond_front[4];

  cond_front[0] = mater->CreateBC(-1, 0, val1, val2);
  cond_front[1] = mater->CreateBC(-2, 0, val1, val2);
  cond_front[2] = mater->CreateBC(-3, 0, val1, val2);
  cond_front[3] = mater->CreateBC(-4, 0, val1, val2);

  cond_front[0]->SetForcingFunction(CC1);
  cond_front[1]->SetForcingFunction(CC2);
  cond_front[2]->SetForcingFunction(CC3);
  cond_front[3]->SetForcingFunction(CC4);


//                                             // 1 cond. fronteira (negativo)
//                                             // 2 Dirichlet
//                                             // 3 Caso as condicoes s� mistas
//                                             // 4 Dirichlet ou Newmann
//   // Se a cond. de front. for uma funcao , ver TPZMaterial::fForcingFunction para BndCond::Contribute
  // cond_front[0]->SetForcingFunction(Dirichlet1);             Tem algun problema.
  //  cond_front[1]->SetForcingFunction(Dirichlet2);

  malha->InsertMaterialObject(mater);
  // malha->InsertMaterialObject(cond_front[0]);
  malha->InsertMaterialObject(cond_front[0]);
  malha->InsertMaterialObject(cond_front[1]);
  malha->InsertMaterialObject(cond_front[2]);
  malha->InsertMaterialObject(cond_front[3]);
*/


  // Para usar elementos descontinuos. Caso for subtraido, ser�elementos continuos.
   TPZGeoElement<TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);


  malha->AutoBuild();

  return malha;

}//end of CreateMesh
