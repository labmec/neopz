/**
 * @file
 * @brief Implements a tutorial example using geometric NeoPZ module
 */
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzlog.h"

#include <fstream>

using std::ifstream;
using std::ofstream;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.extend"));
#endif

/// Read the geometric data from a file
void LerMalhaGeom(const std::string &arquivo,TPZGeoMesh &mesh);
/// Test the capability of generating a large geometric mesh
void LargeMesh(int nrefloop);
/// Test the prismatic extension of the topology
void TestTopology();

/// Program to exemplify the reading a a geometric mesh
int main_back() {
    
#ifdef LOG4CXX
	InitializePZLOG();
#endif
    
    /// test the prismatic extension of the topology
    TestTopology();
    
	TPZGeoMesh malha2;
	LerMalhaGeom("../quad_st_800_1R_16X.gri",malha2);
	
	ofstream output("output.dat");
	malha2.Print(output);
    
    std::ofstream graphfile("output.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&malha2, graphfile);
	return 0;
}

void LerMalhaGeom(const std::string &nome, TPZGeoMesh &grid) {
	ifstream infile(nome.c_str());
    
	int linestoskip;
	char buf[256];
	infile >> linestoskip;
	int i,j;
	for(i=0; i<linestoskip;i++) infile.getline(buf,255);
	infile.getline (buf,255);
	infile.getline (buf,255);
    /**
     * ntri  : number of triangles
     * npoin : number of pints
     * nbouf : number of boundary faces
     * nquad : number of quadrilaterals
     */
	int ntri,npoin,nbouf,nquad,nsidif;
	infile >> ntri >> npoin >> nbouf >> nquad >> nsidif;
	infile.getline (buf,255);
	infile.getline(buf,255);
    
	grid.NodeVec ().Resize(npoin+1);
	TPZVec<int> nodeindices(4);
	int mat, elid;
	for(i=0;i<nquad;i++) {
		infile >> elid;
		for(j=0; j<4;j++) infile >> nodeindices[j];
		infile >> mat;
        int index;
		grid.CreateGeoElement(EQuadrilateral,nodeindices,mat,index,1);
	}
	infile.getline(buf,255);
	infile.getline(buf,255);
    
	int nodeid,dum;
	char c;
	TPZVec<REAL> coord(3,0.);
	for(i=0; i<npoin; i++) {
		infile >> nodeid >> coord[0] >> coord[1] >> c >> dum;
		grid.NodeVec ()[nodeid].Initialize (nodeid,coord,grid);
	}
	infile.getline (buf,255);
	infile.getline (buf,255);
    
	TPZVec<int> sideid(2,0);
	for(i=0; i<nbouf; i++) {
		infile >> sideid[0] >> sideid[1] >> elid >> dum >> mat;
		TPZGeoEl *el = grid.ElementVec ()[elid-1];
		int side = el->WhichSide (sideid);
		TPZGeoElBC(el,side,-mat);
	}
	grid.BuildConnectivity();
    
	return;
}

#include "PrismExtend.h"
#include "TPZGeoExtend.h"
#include "pzshapepoint.h"
#include "pzshapeextend.h"
#include "tpzpoint.h"
#include "pzgeopoint.h"


void TestTopology()
{
    // A line is a prismatic extension of a point
//    pztopology::Pr<pztopology::TPZPoint> line;
    typedef pztopology::Pr<pztopology::TPZPoint> tline;
    // A quadrilateral is a prismatic extension of a line
    pztopology::Pr<pztopology::Pr<pztopology::TPZPoint> > quad;
    // A prismatic extension of a point shape function is a line
    pzshape::SPr<pzshape::TPZShapePoint> shline;
    // the topologic extension of a point is a line (geometry)
    pzgeom::GPr<pzgeom::TPZGeoPoint,pztopology::Pr<pztopology::TPZPoint> > geoline;
    TPZFMatrix<REAL> coordd(3,2,0.);
    int ic;
    for(ic=0;ic<3;ic++) coordd(ic,1) = 1.;
    
    tline::Diagnostic();
    quad.Diagnostic();
	geoline.Diagnostic(coordd);
}

void LargeMesh(int nrefloop)
{
    
	double coordstore[][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}
        , {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}};
	// criar um objeto tipo malha geometrica
	TPZGeoMesh malha;
    
	// criar quatro nos
	int i,j;
	TPZVec<REAL> coord(3,0.);
	malha.NodeVec().Resize(8);
	for(i=0; i<8; i++) {
		// initializar as coordenadas do no em um vetor
		for (j=0; j<3; j++) coord[j] = coordstore[i][j];
        
		// identificar um espaco no vetor onde podemos armazenar
		// este vetor
		//		int nodeindex = malha.NodeVec ().AllocateNewElement ();
        
		// initializar os dados do no
		malha.NodeVec ()[i].Initialize (i,coord,malha);
	}
    
	// criar um elemento
    
	// initializar os indices dos n�s
	TPZVec<int> indices(8);
	for(i=0; i<8; i++) indices[i] = i;
    
	// O proprio construtor vai inserir o elemento na malha
    int index;
	malha.CreateGeoElement(ECube,indices,1,index,0);
    
    
	malha.BuildConnectivity ();
	
	int iref;
	TPZManVector<TPZGeoEl *> subs;
	for (iref=0; iref<nrefloop; iref++) 
	{
		int nel = malha.NElements();
		int iel;
		for (iel=0; iel<nel; iel++) {
			malha.ElementVec()[iel]->Divide(subs);
		}
		std::cout << "Refinement step " << iref << " Number of elements " << malha.NElements() << std::endl;
		std::cout.flush();
	}
	//DebugStop();
}