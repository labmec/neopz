// Obs:
// O pzintel.cc deste diretorio tem as inclusoes de
//
//  TPZVec<int> order(3,1);
//  intrule.SetOrder(order);
//
// depois das linhas de
//
//  TPZIntPoints &intrule = GetIntegrationRule();
//
// nos metodos calcStiff e CalcEnergy, para validar
// o codigo para 1 pto de integracao.
// Fazer o mesmo para demais pontos!!!!!


#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
//#include "pzerror.h"
#include "pzgeoel.h"
//#include "pzcosys.h"
#include "pzmatrix.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "pzpoisson3d.h"
#include "pzcheckgeom.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "pzmathyperelastic.h"
#include "pzskylmat.h"
#include "pzelmat.h"
#include "pzbndcond.h"

#ifdef _AUTODIFF
#include "fadType.h"
void FADToMatrix(FADFADREAL &U, TPZFMatrix & ek, TPZFMatrix & ef);
#endif

TPZCompMesh *CreateMesh();
int mainFull();
void Assemble(TPZMatrix & stiffness, TPZFMatrix & rhs, int method, TPZCompMesh & Mesh);
TPZCompMesh *CreateMesh();
TPZMatrix * CreateAssemble(TPZFMatrix &rhs, int method, TPZCompMesh & Mesh);

/*
void error(char * err)
{
  PZError << "FADERROR: " << err << endl;
};
*/

int main()
{

 const int numShape  = 1, ndof = 3;
 const int dim = 3;


//void TPZMatHyperElastic::Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &/*sol*/,TPZFMatrix &dsol,REAL weight,
//			  TPZFMatrix &/*axes*/,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {
 TPZFMatrix phi(numShape,1), dphi(dim,numShape), dsol(dim, ndof, 0.);
 TPZFMatrix ek(numShape * ndof, numShape * ndof, 0.), ef(numShape * ndof, 1, 0.),ekFAD(numShape * ndof, numShape * ndof, 0.), efFAD(numShape * ndof, 1, 0.);
 TPZVec<REAL> x(3,0.), sol(3,0.);
 TPZFMatrix axes(3,3,0.),jacinv(3,3,0.);
 REAL weight = 1.;
 int id, idf, ishape, i;
	double xdsol[]= {1.e-3,1.e-2,1.e-4,1.e-1,3.e-2,2.e-2,4.e-5,25.e-4,33.e-3};
	int xdsolcount=0;


 for(ishape = 0; ishape < numShape; ishape++) {
 	phi(ishape,0) = (random()%100)/41.;
	for(id = 0; id<dim; id++) {
	 	dphi(id,ishape) = xdsol[xdsolcount++];
	}
 }

 TPZVec<REAL> u(numShape * ndof);
 for(i = 0; i<numShape*ndof; i++) {
 	u[i] = (random()%100)/41.;
 }

 for(id = 0; id < dim; id++)
 {
     for(idf = 0; idf< ndof; idf++)
     {
          for(ishape = 0; ishape < numShape; ishape ++)
	  {
	       dsol(id, idf) += dphi(id, ishape) * u[idf + ishape* ndof];
	  }
     }
 }
 TPZMatHyperElastic hyp(1, 1.e4, 0.2);
	TPZMaterialData data;
	data.x = x;
	data.sol[0] = sol;
	data.dsol[0] = dsol;
	data.axes = axes;
	data.phi = phi;
	data.dphix = dphi;
    hyp.Contribute(data,weight,ek,ef);

////////////////

#ifdef _AUTODIFF
    
  TPZVec<FADFADREAL> solFAD(ndof);
  TPZVec<FADFADREAL> dsolFAD(ndof * dim);// x, y and z data aligned

  FADREAL defaultFAD(ndof*numShape, (REAL)0., (REAL)0.);
  if(defaultFAD.dx(0)==1.)PZError << "\nError: FAD doesn't have default constructor for parameters: (number of derivatives, default value, default derivative value) !";
  FADFADREAL defaultFADFAD(ndof*numShape, defaultFAD, defaultFAD);

  FADFADREAL U(defaultFADFAD); // Zeroed Energy Value -> ready for contribution

  solFAD.Fill(defaultFADFAD);
  dsolFAD.Fill(defaultFADFAD);

TPZVec<REAL> in;

    for(ishape=0; ishape<numShape; ishape++) {

      for(idf=0; idf<ndof; idf++) {

	solFAD[idf].val().val() += u[idf + ndof * ishape] * phi(ishape,0);
	solFAD[idf].val().fastAccessDx(idf+ishape*ndof) += phi(ishape,0);
	solFAD[idf].fastAccessDx(idf+ishape*ndof).val() += phi(ishape,0);
	for(id=0; id<dim; id++)
	{

	   dsolFAD[id+(idf)*dim].val().val() += u[idf + ishape*ndof] * dphi(id, ishape);
	   dsolFAD[id+(idf)*dim].val().fastAccessDx(idf+ishape*ndof) += dphi(id, ishape);
	   dsolFAD[id+(idf)*dim].fastAccessDx(idf+ishape*ndof).val() += dphi(id, ishape);
	}
      }
    }
/*
    cout << "\nCalcEnergy sol\n" << sol;
    cout << "\nCalcEnergy dsol\n" << dsol;
    cout << "\nCalcEnergy phi\n" << phi;
    cout << "\nCalcEnergy dphix\n" << dphix;
*/

    hyp.ContributeEnergy(x,solFAD,dsolFAD,U,weight);
				   
FADToMatrix(U, ekFAD, efFAD);
///////////////

 ek.Print("Stiffness matrix");
 ef.Print("Right hand side");
 ekFAD.Print("FAD Stiffness matrix");
 efFAD.Print("FAD Right hand side");
 ek -= ekFAD;
 ef -= efFAD;
 ek.Print("Stiffness matrix");
 ef.Print("Right hand side");
 REAL dif = Norm(ek);
 cout << "Difference in norm " << dif << endl;
    
#endif
    
  return 0;
}

/*
int mainFull(){

TPZCompEl::SetgOrder(2);
  TPZCompMesh *cmesh = CreateMesh();
	int nsol = cmesh->Solution().Rows();
	int is;
	for(is=0; is<nsol; is++) {
		cmesh->Solution()(is,0) = (random()%100)/41.;
	}
      TPZAnalysis an (cmesh);
      TPZSkylineStructMatrix strskyl(cmesh);
      TPZFMatrix rhs(24,1,0.);
      TPZMatrix *stiff1 = CreateAssemble(rhs,1, *cmesh);
      TPZMatrix *stiff2 = CreateAssemble(rhs,0, *cmesh);
      TPZFMatrix one(*stiff1);
      TPZFMatrix two(*stiff2);
      one -= two;
      one.Print("difference");
      REAL norm = Norm(one);
      cout << "Norm of difference " << norm << endl;
//      an.SetStructuralMatrix(strskyl);

//      TPZStepSolver direct;
//      direct.SetDirect(ECholesky);
//      an.SetSolver(direct);

//      an.Run();
//      an.Rhs().Print();
//      an.Solution().Print();
  //TPZMatrixSolver::Diagnose();
  return 0;
}



// *************************************
// ************Option 0*****************
// *******Shape Quadrilateral***********
// *************************************
TPZCompMesh *CreateMesh(){

  //malha quadrada de nr x nc
  const	int numrel = 1;
  const	int numcel = 1;
  const int numzel = 1;
  //  int numel = numrel*numcel;
  TPZVec<REAL> coord(3,0.);

  // criar um objeto tipo malha geometrica
  TPZGeoMesh *geomesh = new TPZGeoMesh();

  // criar nos
  int i,j,k;
  for(k=0;k<numzel+1;k++)
  for(i=0; i<(numrel+1); i++) {
    for (j=0; j<(numcel+1); j++) {
      int nodind = geomesh->NodeVec().AllocateNewElement();
      TPZVec<REAL> coord(3);
      coord[0] = j;
      coord[1] = i;
      coord[2] = k;
      geomesh->NodeVec()[nodind] = TPZGeoNode(k*(numzel+1)*(numcel+1)+i*(numcel+1)+j,coord,*geomesh);
    }
  }

  TPZVec<int> indices(8);
  // criação dos elementos
  TPZGeoEl *gel[1];

      indices[0] = 0;
      indices[1] = 1;
      indices[3] = 2;
      indices[2] = 3;
      indices[4] = 4;
      indices[5] = 5;
      indices[7] = 6;
      indices[6] = 7;
      // O proprio construtor vai inserir o elemento na malha
      gel[0] = new TPZGeoElC3d(0,indices,1,*geomesh);

  //Divisão dos elementos
  TPZVec<TPZGeoEl *> sub;

  geomesh->BuildConnectivity();
  //geomesh->Print(cout);

//  gel[0]->Divide(sub);

  // Criação das condições de contorno geométricas
  TPZGeoElBC t3(gel[0],4,-1,*geomesh);
  TPZGeoElBC t4(gel[0],6,-2,*geomesh);


  // Criação da malha computacional
  TPZCompMesh *comp = new TPZCompMesh(geomesh);

  // Criar e inserir os materiais na malha
  TPZMatHyperElastic *mat = new
  TPZMatHyperElastic(1,1.e5,.25);

  comp->InsertMaterialObject(mat);

  TPZMaterial *meumat = mat;

  // Condições de contorno
  // Dirichlet
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  TPZMaterial *bnd = meumat->CreateBC (-1,0,val1,val2);
  comp->InsertMaterialObject(bnd);

  // Neumann
  val2(0,0)=7.;
  val2(1,0)=3.;
  val2(2,0)=5.;
  bnd = meumat->CreateBC (-2,1,val1,val2);
  comp->InsertMaterialObject(bnd);

  //comp->Reference()->Print();
  //comp->Print(cout);

  // Ajuste da estrutura de dados computacional
  comp->AutoBuild();


  comp->AdjustBoundaryElements();
  comp->CleanUpUnconnectedNodes();
    //comp->Print(cout);
    comp->SetName("Malha Computacional Original");
    return comp;
}

TPZMatrix * CreateAssemble(TPZFMatrix &rhs, int method, TPZCompMesh & Mesh){
    int neq = Mesh.NEquations();
    TPZVec<int> skyline;
    Mesh.Skyline(skyline);
    TPZSkylMatrix *stiff = new TPZSkylMatrix(neq,skyline);
    rhs.Redim(neq,1);
    Assemble(*stiff,rhs,method, Mesh);
    return stiff;
}


void Assemble(TPZMatrix & stiffness, TPZFMatrix & rhs, int method, TPZCompMesh & Mesh){

  int iel;
  //int numel = 0;
  int nelem = Mesh.NElements();
  TPZElementMatrix ek,ef;
  TPZManVector<int> destinationindex(0);
  TPZManVector<int> sourceindex(0);
  REAL stor1[1000],stor2[1000],stor3[100],stor4[100];
//  ek.fMat = new TPZFMatrix(0,0,stor1,1000);
//  ek.fConstrMat = new TPZFMatrix(0,0,stor2,1000);
//  ef.fMat = new TPZFMatrix(0,0,stor3,100);
//  ef.fConstrMat = new TPZFMatrix(0,0,stor4,100);

  TPZAdmChunkVector<TPZCompEl *> &elementvec = Mesh.ElementVec();

  for(iel=0; iel < nelem; iel++) {
    TPZCompEl *el = elementvec[iel];
    if(!el) continue;
    //	  int dim = el->NumNodes();
//  if(method == 0) {
    el->CalcStiff(ek,ef);
//  } else {

//    TPZInterpolatedElement * pIntel = NULL;
//    pIntel = dynamic_cast<TPZInterpolatedElement *>(el);
//    if(pIntel)
//    {
//       pIntel->CalcEnergy(ek,ef);
//    }else
//    {
//       el->CalcStiff(ek,ef);
//    }
//   }

    if(!el->HasDependency()) {
      //ek.fMat->Print("stiff has no constraint",test);
      //ef.fMat->Print("rhs has no constraint",test);
      //test.flush();
      destinationindex.Resize(ek.fMat.Rows());
      int destindex = 0;
      int numnod = ek.NConnects();
      for(int in=0; in<numnod; in++) {
         int npindex = ek.ConnectIndex(in);
         TPZConnect &np = Mesh.ConnectVec()[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = Mesh.Block().Position(blocknumber);
         int ndf = Mesh.Block().Size(blocknumber);
	 //	 if (numnod == 27){
	 //   cout << "First equation " << firsteq <<"\t ndf " << ndf << endl;
	 //	 }
         for(int idf=0; idf<ndf; idf++) {
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      stiffness.AddKel(ek.fMat,destinationindex);
      rhs.AddFel(ef.fMat,destinationindex);
    } else {
      // the element has dependent nodes
      el->ApplyConstraints(ek,ef);
      //ek.fMat->Print("stif no constraint",test);
      //ek.fConstrMat->Print("stif constrained",test);
      //ef.fMat->Print("rhs no constraint",test);
      //ef.fConstrMat->Print("rhs constrained",test);
      //test.flush();
      //test << "sum of columns\n";
      int destindex = 0;
      int fullmatindex = 0;
      destinationindex.Resize(ek.fConstrMat.Rows());
      sourceindex.Resize(ek.fConstrMat.Rows());
      int numnod = ek.fConstrConnect.NElements();
      for(int in=0; in<numnod; in++) {
         int npindex = ek.fConstrConnect[in];
         TPZConnect &np = Mesh.ConnectVec()[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = Mesh.Block().Position(blocknumber);
         int ndf = Mesh.Block().Size(blocknumber);
         if(np.HasDependency()) {
           fullmatindex += ndf;
           continue;
         }
         for(int idf=0; idf<ndf; idf++) {
           sourceindex[destindex] = fullmatindex++;
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      sourceindex.Resize(destindex);
      destinationindex.Resize(destindex);
      stiffness.AddKel(ek.fConstrMat,sourceindex,destinationindex);
      rhs.AddFel(ef.fConstrMat,sourceindex,destinationindex);
//if(ek.fConstrMat->Decompose_LU() != -1) {
//    el->ApplyConstraints(ek,ef);
//    ek.Print(*this,check);
//    check.flush();
//}
    }
  }//fim for iel
//
  int neq = rhs.Rows();
//	if(nelem < 34 && neq < 100){
//    stiffness.Print("TPZStructMatrix::Assemble GLOBAL MATRIX (after Assemble)",out);
//    rhs.Print("TPZStructMatrix::Assemble GLOBAL LOAD (after Assemble)",out);
//  }
}

*/

#ifdef _AUTODIFF

void FADToMatrix(FADFADREAL &U, TPZFMatrix & ek, TPZFMatrix & ef)
{
//  int efsz = ef.Rows();
//  int ekrows = ek.Rows();
//  int ekcols = ek.Cols();

  int Ucols = U.size();
  int Urows = U.val().size();


  FADREAL * pBufferFAD;
  int i,j;
  for(j = 0; j < Urows; j++)
  {
     pBufferFAD = & U.fastAccessDx(j);
     ef(j,0) = - pBufferFAD->val();
     // U.val().fastAccessDx(i); must be the same as U.fastAccessDx(i).val();
     for(i = 0; i < Ucols; i++)
     {
        ek(i,j) = pBufferFAD->fastAccessDx(i);
     }
  }
}
#endif