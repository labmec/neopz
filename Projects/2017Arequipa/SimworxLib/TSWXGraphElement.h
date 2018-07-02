//---------------------------------------------------------------------------

#ifndef TSWXGraphElementH
#define TSWXGraphElementH
//---------------------------------------------------------------------------

#include <map>
#include <set>
#include <fstream>

#include "pzvec.h"
#include "pzeltype.h"
#include "pzgnode.h"
#include "pzgmesh.h"
//#include "TSWXVector.h"
#include "pzmanvector.h"
#include "pzvec.h"
#include "TSWXGraphMesh.h"

class TSWXGraphMesh;
class TPZInterpolationSpace;

class TSWXGraphElement
{
protected:

void CollectElements(TPZCompMesh *cmesh, std::set< TPZInterpolationSpace* > &AllEls);

public:
	// construtor com o número de refinamentos desejado
	TSWXGraphElement( int NRef );

	~TSWXGraphElement();

  void GenerateVTKData(TPZCompMesh * cmesh,
                       int dim,
                       double time,
                       const TPZVec<std::string> & nodalVarIndex,
                       const TPZVec<std::string> & cellVarIndex,
                       TSWXGraphMesh & graphMesh);

  void GenerateVTKData(TPZCompMesh * cmesh,
                       int dim,
                       double time,
                       const TPZVec<std::string> & nodalVarIndex,
                       const TPZVec<std::string> & cellVarIndex,
                       TSWXGraphMesh & graphMesh,
                       const std::set< int > &mySetMatIds);

  /** Define solucoes constantes no elemento
   * Apenas solucoes escalares, por enquanto
   */
  void SetConstScalarSolution(const std::map< std::string, TPZVec<REAL > > &sol);

	void Write (std::ostream &Out);

	bool Read (std::istream & In);

private:

  /** Define solucoes constantes no elemento
   * Apenas solucoes escalares, por enquanto
   */
  std::map< std::string, TPZVec<REAL > > fConstScalarSolution;

  /** Gera malha vtk.
   * Parametro mySetMatIds contem material ids que devem ser exportados
   * Se mySetMatIds.size = 0, nenhum material é excluido
   */
  void GenerateVTKData(TPZCompMesh * cmesh,
                       int dim,
                       const TPZVec<std::string> & nodalVarIndex,
                       const TPZVec<std::string> & cellVarIndex,
                       TPZStack< TSWXGraphNode > &nodes,
                       TPZStack< TSWXGraphEl > &elem,
                       TPZVec< TSWXGraphSol > &sol,
                       const std::set< int > &mySetMatIds,
                       bool updateMesh);

  void GetSolution(TPZInterpolationSpace *sp,
                         std::string varName,
                         TPZVec< TPZVec< REAL > > &ParametricCoord,
                         TPZVec< TSWXGraphSingleSol > &locSol);

	//Para cada tipo de elemento retorna as coordenadas de nós dos filhos
	// e a conectividade gerada
	void GetData ( MElementType Type, TPZVec<TPZGeoNode> & NodeVec,
									TPZVec < TPZVec < int > > & ConnectivityVec,
                  TPZVec< TPZGeoNode> &CenterOfSubEls );

	//Dada uma lista de conectividades PZ ( nós relativos a cada lado do elemento)
	//retorna o vetor de conectividades vtk
	EGraphElType PZToVtkConnecitivities ( TPZVec< int > & ConnectivityVec);

	std::map < MElementType, TPZVec<TPZGeoNode> > fNodeData;
	std::map < MElementType, TPZVec < TPZVec<int> > > fConnectivityData;

	void ProcessAllTopologies( int nRef );
	void ProcessElement ( MElementType Type, int NRef, 
												TPZVec <TPZGeoNode> & AllNodes, 
												TPZVec < TPZVec<int> > & AllElements );
	
	void CreateElement ( TPZGeoMesh & gMesh, MElementType & type );
	
	void Divide ( TPZGeoMesh &gMesh, int nRef );
	
};

#endif
