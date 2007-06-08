#include "includes.h"

#include "c0-simplequad.cc"
#include "c01-lshape-quad.cc"
#include "c02-lshape-triang.cc"
#include "c03-lshape-quadtri.cc"
#include "c04-simplehexa.cc"
#include "c05-cornercube.cc"
#include "c06-Tetra.cc"
#include "c07-Prism.cc"
#include "c08-Mixed.cc"
#include "c09-Plate-Cedric.cc"
#include "c10-Prism-Pref.cc"
#include "c11-Piram-Tetra.cc"
#include "c15-Exp-Hexa.cc"

int gPrintLevel = 0;
void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp);
static REAL angle = 0.2;

//int gDebug = 0;
TPZCompMesh *ReadCase(int &nref, int &dim, int &opt);

TPZCompMesh *ReadKumar(char *filename);
int MaxLevel(TPZCompMesh *mesh);


void Exact(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol);
void ExactSimple3D(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol);
void Exact3D(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol);


static ofstream MALHAG("malhageometrica");//CEDRIC
static int mygorder = 1;
void CompareNeighbours(TPZGeoMesh *mesh);


int main(){

  int nref = 0;
  int dim = 0;
  int opt = 0;

  TPZCompEl::SetgOrder(mygorder);
  gDebug = 0;

  TPZCompMesh *cmesh = ReadCase(nref,dim,opt);
  TPZStack<TPZGeoEl *> gelstack;
  TPZStack<int> porders;

  //	TPZCheckGeom::main();
  //	return 0;
  //	return (TPZGeoCloneMesh::main());

  ofstream convergence("conv3d.txt");

  cmesh->Reference()->SetName("Malha Geométrica original");
  // cmesh->Reference()->Print(cout);

  cmesh->SetName("Malha Computacional Original");
  //  cmesh->Print(cout);

  cmesh->CleanUpUnconnectedNodes();
  TPZStack<char *> scalnames, vecnames;
  scalnames.Push("POrder");
  scalnames.Push("Error");

  if(nstate == 1) {
    //   an.SetExact(Exact);
    scalnames.Push("TrueError");
    scalnames.Push("EffectivityIndex");
  }else if(nstate == 2) {
    scalnames.Push("sig_x");
    scalnames.Push("sig_y");
    scalnames.Push("tau_xy");
  }

  vecnames.Push("state");
  if(dim < 3){
    vecnames.Push("displacement");
  }
  ofstream out("output.txt");

  //Multigrid======================

//   TPZMGAnalysis mgan (cmesh);
//   mgan.SetStructuralMatrix(strskyl);
//   TPZStepSolver *direct = new TPZStepSolver;
//   direct->SetDirect(ELDLt);
//   mgan.SetSolver(*direct);
//   delete direct;
//   direct = 0;
//   mgan.Run();
//   TPZCompMesh *finemesh = cmesh;
  // ===================================
  REAL valerror =0.;
  REAL valtruerror=0.;
  TPZManVector<REAL> ervec,truervec,effect;

  {
    int r;
    for(r=0; r<nref; r++) {
      //      if (r == 7) gDebug = 1;
      {
        TPZAnalysis an (cmesh);
        if (opt == 4 || opt == 7 || opt == 8){
          an.SetExact(ExactSimple3D);
          // an.SetExact(Exact3DExp);
        }
        else if(opt==1 || opt==2){
          an.SetExact(Exact);
        }
        else if(opt==12){
          an.SetExact(Exact3D);
        }
        char buf [256];
        sprintf(buf,"hptest%d.dx",r);
        an.DefineGraphMesh(dim,scalnames,vecnames,buf);
        cmesh->SetName("Malha computacional adaptada");

        if (gDebug == 1){
          cmesh->Reference()->Print(cout);
          cmesh->Print(cout);
        }

        TPZSkylineStructMatrix strskyl(cmesh);
        an.SetStructuralMatrix(strskyl);

        TPZStepSolver direct;
        direct.SetDirect(ECholesky);
        an.SetSolver(direct);

        an.Run();
        //an.Rhs().Print();
        //an.Solution().Print();

        if (r==nref -1)
          an.PostProcess(0,dim);

        ofstream outpira("/compile/pira.out",ios::app);
        TPZVec<REAL> pospira;
        an.PostProcess(pospira, outpira);

      }
      time_t sttime;

      {
        TPZAdaptMesh adapt(cmesh);
        //      adapt.SetCompMesh (cmesh);

        cout << "\n\n\n\nEntering Auto Adaptive Methods... step " << r << "\n\n\n\n";

        //	if(r==4) gDebug = 1;
        time (& sttime);
        TPZCompMesh *adptmesh;

        switch (opt){
          case (1) :{
            adapt.GetAdaptedMesh(gelstack,porders,valerror,valtruerror,ervec,Exact,truervec,effect,0);
            break;
          }
          case (2) :{
            adapt.GetAdaptedMesh(gelstack,porders,valerror,valtruerror, ervec,Exact,truervec,effect,0);
            break;
          }
          case (4) :{
            adapt.GetAdaptedMesh(gelstack,porders,valerror,valtruerror, ervec,ExactSimple3D,truervec,effect,0);
            break;
          }
          case (12) : {
            adapt.GetAdaptedMesh(gelstack,porders,valerror,valtruerror, ervec,Exact3D,truervec,effect,0);
            break;
          }
          case (15):{
            adapt.GetAdaptedMesh(gelstack,porders,valerror,valtruerror,ervec,Exact3DExp,truervec,effect,0);
            break;
          }
          default:
          adapt.GetAdaptedMesh(gelstack,porders,valerror,valtruerror, ervec,0,truervec,effect,0);
        }
      }
      time_t endtime;
      time (& endtime);

      int time_elapsed = endtime - sttime;
      cout << "\n\n\n\nExiting Auto Adaptive Methods....step " << r
           << "time elapsed " << time_elapsed << "\n\n\n\n";

      int prt;
      cout << "neq = " << cmesh->NEquations() << " erestimate = " << valerror
           << " true " << valtruerror <<  " effect " << valerror/valtruerror << endl;

      convergence   << cmesh->NEquations() << "\t"
                    << valerror << "\t" << valtruerror << "\t"
                    << ( valtruerror / valerror ) <<  "\t" << sttime <<endl;

      for (prt=0;prt<ervec.NElements();prt++){
        cout <<"error " << ervec[prt] << "  truerror = " << truervec[prt] << "  Effect " << effect[prt] << endl;
        // convergence << '\t' << ervec[prt] << '\t' << truervec[prt] << "  Effect " << effect[prt] <<  endl;
        //  adptmesh->Print(cout);
      }

      //      if (r==6){
      //adaptmesh->Print(cout);
      //      }
      //     adptmesh->Print(cout);
      cout.flush();
      cmesh->Reference()->ResetReference();
      TPZCompMesh *adaptmesh;
      adaptmesh = TPZAdaptMesh::CreateCompMesh(cmesh,gelstack,porders);
      gelstack.Resize(0);
      porders.Resize(0);
      delete cmesh;
      cmesh = adaptmesh;
      adaptmesh->Print(cout);
      cout.flush();
      //adptmesh->Print(MALHAG);//CEDRIC

      /*   if (r == (nref-1)){ */
      /*        an.PostProcess(2,2); */
      /*        cout << "The maximum level = " << MaxLevel(cmesh) << endl; */
      /*      } */
    }
  }
  TPZMatrixSolver::Diagnose();
  CompareNeighbours(cmesh->Reference());
  delete cmesh;
  return 0;
}


TPZCompMesh *ReadCase(int &nref, int &dim, int &opt){

  cout << "**************************************" << endl;
  cout << "******Auto adaptive test code********" << endl;
  cout << "**************************************" << endl;

  cout << "Select the analysis type: \n0 - Simple quadrilateral 2D \n1 - L Shape Quadrilateral\n"
       << "2 - Triangular Simples \n3 - Plane mesh (quadrilateral and triangular elements)"
       << "\n4 - 3D Simples \n5 - 3D Canto\n" <<"6 - Tetraedro\n7 - Prisma\n8 - All elements\n9 - All topologies\n10 Aleatorio\n"
       << "11 Pyramid and Tetrahedre\n12Exact 3d Poisson\n"
       << "15 Cube Exp\n";

  cin >> opt;

  TPZCompMesh *cmesh;

  switch (opt){
  case (0) :{
    cmesh = CreateSillyMesh();
    break;
  }
  case (1) :{
    cmesh = CreateMesh();
    break;
  }
  case (2) :{
    cmesh = CreateTriangularMesh();
    break;
  }
  case (3):{
    cmesh = CreatePlanMesh();
    break;
  }
  case (4) :{
    cmesh = CreateSimple3DMesh();
    break;
  }
  case (5) :{
    cmesh = Create3DMesh();
    break;
  }
  case (6) :{
    cmesh = Create3DTetraMesh();
    break;
  }
  case (7) :{
    cmesh = Create3DPrismMesh();
    break;
  }
  case (8) :{
    cmesh = CreateTestMesh();
    break;
  }
  case (10) :{
    cmesh = CreateAleatorioMesh();
    break;
  }
  case (11) :{
    cmesh = CreatePyramTetraMesh();
    break;
  }
  case (12) :{
    cmesh = Create3DDiscMesh();
    break;
  }
  case (15) :{
    cmesh = Create3DExpMesh();
    break;
  }

  default:
    cmesh = CreateMesh();
  }

  dim = 2;
  opt > 3 ? dim=3 : dim = 2;

  cout << "number of refinement steps : ";
  cin >> nref;

  cout << "Maximum p order:    ";
  int p;
  cin >> p;
  cout << endl;

  TPZOneDRef::gMaxP = p;

  return cmesh;
}


void Exact(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol) {
  	REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
  	REAL theta = atan2(x[1],x[0]);
  	REAL rexp = pow(r,onethird);
  	sol[0] = rexp*sin(onethird*(theta+PI/2));
  	dsol(0,0) = onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
  	dsol(1,0) = onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}

void Exact3D(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol) {
  	REAL r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  	REAL theta = atan2(x[1],x[0]);
  	REAL rexp = pow(r,onethird);
  	sol[0] = rexp*sin(onethird*(theta+PI/2));
  	dsol(0,0) = onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
  	dsol(1,0) = onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}


void ExactSimple3D(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol) {
	sol[0] = x[2];
  	dsol(0,0) = 0.;
  	dsol(1,0) = 0.;
	dsol(2,0) = 1.;
}


void NeumannExp(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = pow(exp(1.),1/(0.1 + pow(-0.5 + x[0],2) + pow(-0.5 + x[1],2)))*(1 - x[0])*x[0]*(1 - x[1])*x[1];
}


int MaxLevel(TPZCompMesh *mesh) {
  	int nel = mesh->NElements();
  	int el;
  	int level = 0;
  	for(el=0; el<nel; el++) {
    		TPZCompEl *cel = mesh->ElementVec()[el];
    		if(!cel) continue;
    		TPZGeoEl *gel = cel->Reference();
    		if(!gel) continue;
    		int gellev = gel->Level();
    		level = (level <gellev) ? gellev : level;
  	}
  	return level;
}

void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp){
        disp[0] = -(x[1]-0.5)*sin(angle)+(x[0]-0.5)*cos(angle)-(x[0]-0.5);
        disp[1] = (x[1]-0.5)*cos(angle)+(x[0]-0.5)*sin(angle)-(x[1]-0.5);
        disp[2] = 0.;
}

TPZCompMesh *ReadKumar(char *filename) {

  	int nnodes,nelem,nmat,nbcd,nbc;
  	TPZGeoMesh *gmesh = new TPZGeoMesh();
  	ifstream input(filename);
  	if(!input) {
    		cout << "file could not be opened " << filename << endl;
    		return 0;
  	}
  	char buf[256];
  	input.getline(buf,256);
  	input >> nnodes >> nelem >> nmat >> nbcd >> nbc;
  	gmesh->NodeVec().Resize(nnodes);
  	input.getline(buf,256);
  	while(buf[0] != '#') input.getline(buf,256);
  	int id;
  	TPZVec<REAL> coord(2);
  	int nod;
  	for(nod=0; nod< nnodes; nod++) {
    		input >> id >> coord[0] >> coord[1];
    		gmesh->NodeVec()[id] = TPZGeoNode(id,coord,*gmesh);
  	}

  	input.getline(buf,256);
  	while(buf[0] != '#') input.getline(buf,256);
  	char c = input.peek();
  	while(c == '#') {
    		input.getline(buf,256);
    		c = input.peek();
  	}

  	int nel;
  	TPZVec<int> elvertices(4);
  	int elnodes[9],matindex;
  	for(nel=0; nel<nelem; nel++) {
    		for(nod=0; nod<9; nod++) input>>elnodes[nod];
    			input >> matindex;
    		for(nod=0; nod<4; nod++) elvertices[nod] = elnodes[nod];
		//    		new TPZGeoElQ2d(nel,elvertices,matindex,*gmesh);
		int index;
		gmesh->CreateGeoElement(EQuadrilateral,elvertices,matindex,index);
  	}

  	gmesh->BuildConnectivity2();
	input.getline(buf,256);
	char *compare = strstr(buf,"# BC records");

  	while(!compare) {
    		input.getline(buf,256);
    		compare = strstr(buf,"# BC records");
  	}

  	int elnum, side, bcnum;
  	int bc;
  	for(bc=0; bc<nbc; bc++) {
    		input >> elnum >> side >> bcnum;
    		TPZGeoElBC(gmesh->ElementVec()[elnum-1],side+4-1,-bcnum-1,*gmesh);
  	}

 	while(!strstr(buf,"real value)")) input.getline(buf,256);
  	REAL e1111,e1122,e2222,e1212;
  	input >> buf >> e1111 >> buf >> e1122 >> buf >> e2222 >> buf >> e1212;

  	REAL E,nu;
  	nu = e1122/e1111;
  	E = e1212*(1+nu);

  	TPZElasticityMaterial *mat = new TPZElasticityMaterial(3,E,nu,0.,0.);
  	TPZFMatrix val1(2,2,0.),val2(2,1,0.);
  	TPZBndCond *bc1 = mat->CreateBC(-1,0,val1,val2);
  	val2(1,0) = -1.;
  	TPZBndCond *bc2 = mat->CreateBC(-2,1,val1,val2);
  	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  	cmesh->InsertMaterialObject(mat);
  	cmesh->InsertMaterialObject(bc1);
  	cmesh->InsertMaterialObject(bc2);

  	cmesh->AutoBuild();
  	cmesh->AdjustBoundaryElements();
  	cmesh->CleanUpUnconnectedNodes();
  	//  TPZCompMesh *fine = TPZMGAnalysis::UniformlyRefineMesh(cmesh);
  	//  delete cmesh;
  	return cmesh;
}

void CompareNeighbours(TPZGeoMesh *mesh) {

  TPZAdmChunkVector<TPZGeoEl *> &geovec = mesh->ElementVec();
  int nel = geovec.NElements();
  int iel;
  for(iel=0; iel<nel; iel++) {
    TPZGeoEl *gel = geovec[iel];
    if(!gel) continue;
    int nsides = gel->NSides();
    int is;
    for(is=0; is<nsides; is++) {
      TPZStack<TPZGeoElSide> st1,st2;
      TPZGeoElSide gelside = TPZGeoElSide(gel,is);
      gelside.AllNeighbours(st1);
      gelside.ComputeNeighbours(st2);
      Sort<TPZGeoElSide>(st1);
      Sort<TPZGeoElSide>(st2);
      int nlist1 = st1.NElements();
      int nlist2 = st2.NElements();
      if(nlist1 != nlist2) {
        cout << "AllNeighbours is different form ComputeNeighbours\n";
        continue;
      }
      int il;
      for(il=0; il<nlist1; il++) {
        if(st1[il] != st2[il]) {
          cout << "Different neighbours\n";
        }
      }
    }
  }
}


