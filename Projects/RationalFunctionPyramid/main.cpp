/**
 * @file
 * @brief Projeto elaborado para validar as definições feitas pelo Cedric para construir funções racionais para elementos piramidais.
 */

#include "pzshapelinear.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

#include "pzbstrmatrix.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "pzcheckmesh.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzpoisson3d.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "CedricTest.h"

#include "TPZRefPatternTools.h"

#include "TPZParSkylineStructMatrix.h"

#include <stdio.h>

#include "pzlog.h"

#include "pzgeoelbc.h"

/** Initialiazing file for Log4CXX for this project. */
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.Cedric"));
#endif

// output files  -> Because it has many energy faults
std::ofstream out("output.txt");

/** Using standard library for C++ */
using namespace std;
/** Using shape library from NeoPZ code */
using namespace pzshape;

static bool gDebug = false;
int POrder = 1;

int MaterialId = 1;
int MaterialBC1 = -1;

/** To print geometric elements with respective dimension */
void PrintGeoMeshAsCompMeshInVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=false, const int matidtodivided=1);

// MAIN FUNCTION FOR NUMERICAL TESTS TO CEDRIC CLASS
int main(int argc, char *argv[]) {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
    
    // launch arguments
    // porderinit porderend  gcaseinit(meshtype)  gcaseend nsubdivisoninit nsubdivisionend nsubdivisioninterval
    
    int gcaseinit = 1, gcaseend = 3;
    int nsubdivisionsinit = 3, nsubdivisionsend = 16, nsubdivisionsinterval = 3;
    int porderinit = 2, porderend = 4;
    switch(argc) {
        case 8:
            nsubdivisionsinterval = atoi(argv[7]);
        case 7:
            nsubdivisionsend = atoi(argv[6]);
        case 6:
            nsubdivisionsinit = atoi(argv[5]);
        case 5:
            gcaseend = atoi(argv[4]);
        case 4:
			gcaseinit = atoi(argv[3]);
        case 3:
            porderend = atoi(argv[2]);
        case 2:
            porderinit = atoi(argv[1]);
            break;
        case 1:
            std::cout << "\nCommand line: Program p_order geomesh_initial geomesh_last n_initial n_final dn_interval (h = 1/n)" << std::endl;
            std::cout << "\nRunning with: Executavel p " << porderinit << " - " << porderend << " " << gcaseinit << " " << gcaseend << " " << nsubdivisionsinit << " " << nsubdivisionsend << " " << nsubdivisionsinterval << "\n";
            out << "\nRunning with: Executavel " << porderinit << " - " << porderend << " " << gcaseinit << " " << gcaseend << " " << nsubdivisionsinit << " " << nsubdivisionsend << " " << nsubdivisionsinterval << "\n";
            break;
        default:
            std::cout << "\nBad number of arguments. Finishing." << std::endl;
            return 1;
    }
    if(porderinit < 1 || porderinit > porderend || gcaseinit < 0 || gcaseend < 0 || nsubdivisionsinit < 1 || nsubdivisionsend <= nsubdivisionsinit || nsubdivisionsinterval < 1) {
        std::cout << "\nBad parameter.";
        out << "\nBad parameter.";
    }
    if(porderend > 10 || gcaseinit > 4 || gcaseend > 4 || nsubdivisionsinit > 100 || nsubdivisionsend > 100) {
        std::cout << "\nParameter out of avaliable limit.";
        out << "\nParameter out of avaliable limit.";
    }

	// Initializing a ref patterns
	gRefDBase.InitializeAllUniformRefPatterns();
    
    // To store errors
    std::ofstream arq("Errors.txt",ios::app);
    std::cout << "\nRUNNING CedricTest.\n\n";
    out << "\nRUNNING CedricTest.\n\n";
    arq << "\nRUNNING CedricTest:\n\n";
    
    // setting p order
    /** Set polynomial order */
    for(POrder=porderinit;POrder<porderend+1;POrder++) {
        std::cout << "\nInterpolation order " << POrder << std::endl;
        out << "\nInterpolation order " << POrder << std::endl;
        TPZCompEl::SetgOrder(POrder);
        
        TCedricTest cedric;
        // Loop over type of element: geocase = 1(hexahedra), 2(Pyramid+Tetrahedra), 3(Tetrahedra)
        for(int gcase=gcaseinit;gcase<gcaseend;gcase++) {
            std::cout << "\n\tCase " << gcase;
            out << "\n\tCase " << gcase;
            for(int nsubdivisions=nsubdivisionsinit;nsubdivisions<nsubdivisionsend;nsubdivisions+=nsubdivisionsinterval) {
                std::cout << "\n\t\tNumber of sub-divisions " << nsubdivisions << " " << nsubdivisionsinit << " " << nsubdivisionsend << " " << nsubdivisionsinterval << " ";
                out << "\n\t\tNumber of sub-divisions " << nsubdivisions << " " << nsubdivisionsinit << " " << nsubdivisionsend << " " << nsubdivisionsinterval << " ";
                out << "\nApproximation error ";
                arq << "\nApproximation error ";
                cedric.Run(nsubdivisions,gcase,POrder,MaterialId,arq);
                if(POrder == 1)
                {
                    out << "\nInterpolation errors";
                    arq << "\nInterpolation errors";
                    cedric.InterpolationError(nsubdivisions,gcase,MaterialId,arq);
                }
            }
        }
    }
    return 0;
}


void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial, const int matidtodivided) {
  TPZManVector<TPZGeoEl*> filhos;
  for(int D=0; D<nDiv; D++)
  {
      int nels = gmesh->NElements();
      for(int elem = 0; elem < nels; elem++)
      {    
          TPZGeoEl * gel = gmesh->ElementVec()[elem];
          if(!gel || gel->HasSubElement())
              continue;
          if(dim > 0 && gel->Dimension() != dim) continue;
          if(!allmaterial) {
              if(gel->MaterialId() == matidtodivided){
                  gel->Divide(filhos);
              }
          }
          else{
              gel->Divide(filhos);
          }
          if(gDebug) {
              REAL volgel = fabs(gel->Volume());
              REAL sumvol = 0.;
              for(int nsubs=0;nsubs<gel->NSubElements();nsubs++)
                  sumvol += fabs(filhos[nsubs]->Volume());
              if(!IsZero(volgel-sumvol)) {
                  std::cout << "Division of geometric element " << elem << " is wrong.\n";
                  DebugStop();
              }
          }
      }
  }
  gmesh->ResetConnectivities();
  gmesh->BuildConnectivity();
}


////////////////////////////////////////////////////////////////////////////////////////
//////////   FICHERA CORNER - Problem as Anders Solin Presentation   ///////////////////
////////////////////////////////////////////////////////////////////////////////////////

TPZGeoMesh *ConstructingFicheraCorner() {
    REAL co[8][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {1.,1.,1.},
        {0.,1.,1.}
    };
    int indices[1][8] = {{0,1,2,3,4,5,6,7}};
    
    const int nelem = 1;
    int nnode = 8;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(8);
        for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
    
 //   TPZGeoElBC gbc1(elvec[0],20,-1);
 //   TPZGeoElBC gbc2(elvec[0],25,-2);
    return gmesh;
}

void ExactSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	sol[0] = sqrt( sqrt (quad_r) );
	if(!IsZero(sol[0])) {
		REAL den = sol[0]*sol[0]*sol[0];
		dsol(0,0) = .5*x[0]/den;
		dsol(1,0) = .5*x[1]/den;
		dsol(1,0) = .5*x[2]/den;
	}
	else {
		dsol(0,0) = dsol(1,0) = dsol(2,0) = 0.;
	}
}

void BCSolin(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	bcsol[0] = sqrt( sqrt (quad_r) );	
}

void Ff(const TPZVec<REAL> &x, TPZVec<STATE> &f) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	REAL raiz = sqrt( sqrt(quad_r));
	f[0] = -3./(4.*(raiz*raiz*raiz));
}

TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction) {
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
	// Creating Poisson material
	TPZMaterial *mat = new TPZMatPoisson3d(MaterialId,dim);
	TPZVec<REAL> convd(3,0.);
	((TPZMatPoisson3d *)mat)->SetParameters(1.,0.,convd);
	if(hasforcingfunction) {
		mat->SetForcingFunction(new TPZDummyFunction<STATE>(Ff, 5));
	}
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	cmesh->SetAllCreateFunctionsContinuous();
    
	// Boundary conditions
	// Dirichlet
	TPZAutoPointer<TPZFunction<STATE> > FunctionBC = new TPZDummyFunction<STATE>(BCSolin, 5);
	TPZFMatrix<STATE> val1(dim,dim,0.),val2(dim,1,0.);
	TPZMaterial *bnd = mat->CreateBC(mat,MaterialBC1,0,val1,val2);
	bnd->SetForcingFunction(FunctionBC);
	cmesh->InsertMaterialObject(bnd);
	
	cmesh->AutoBuild();
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        TPZCheckMesh tst(cmesh,&sout);
        tst.VerifyAllConnects();
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif    
    
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
	cmesh->CleanUpUnconnectedNodes();
	return cmesh;
}

TPZCompMesh *CreateMeshToLaplace(TPZGeoMesh *gmesh,int dim,int hasforcingfunction) {
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
	// Creating Poisson material
	TPZMaterial *mat = new TPZMatPoisson3d(MaterialId,dim);
//	TPZVec<REAL> convd(3,0.);
//	((TPZMatPoisson3d *)mat)->SetParameters(1.,0.,convd);
	if(hasforcingfunction) {
		mat->SetForcingFunction(new TPZDummyFunction<STATE>(Ff, 5));
	}
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	cmesh->SetAllCreateFunctionsContinuous();
    
	// Boundary conditions
	// Dirichlet
	TPZAutoPointer<TPZFunction<STATE> > FunctionBC = new TPZDummyFunction<STATE>(BCSolin, 5);
	TPZFMatrix<STATE> val1(dim,dim,0.),val2(dim,1,0.);
	TPZMaterial *bnd = mat->CreateBC(mat,MaterialBC1,0,val1,val2);
	bnd->SetForcingFunction(FunctionBC);
	cmesh->InsertMaterialObject(bnd);
	
	cmesh->AutoBuild();
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        TPZCheckMesh tst(cmesh,&sout);
        tst.VerifyAllConnects();
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
	cmesh->CleanUpUnconnectedNodes();
	return cmesh;
}

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int64_t nelem = 0;
	int64_t ngelem=gmesh->NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			REAL distancia = TPZGeoEl::Distance(center,FirstNode);
			if(distancia > distance) distance = distancia;
			isdefined = true;
		}
		REAL centerdist = TPZGeoEl::Distance(center,point);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
}

/** To print geometric elements with respective dimension */
void PrintGeoMeshAsCompMeshInVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename) {
	int64_t i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		TPZGeoEl *gel = gmesh->ElementVec()[i];
		if(gel && gel->Reference())
			DataElement[i] = gel->Dimension();
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}
