//---------------------------------------------------------------------------


#include "TSWXGraphElement.h"
#include "pzinterpolationspace.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzsubcmesh.h"

//---------------------------------------------------------------------------

TSWXGraphElement::TSWXGraphElement( int NRef ) : fNodeData(), fConnectivityData(),
                                                 fConstScalarSolution()
{
	ProcessAllTopologies (NRef);
}

TSWXGraphElement::~TSWXGraphElement()
{
	//Nothing to be done here
}

void TSWXGraphElement::Write (std::ostream &Out)
{
	std::map< MElementType, TPZVec < TPZGeoNode > >::iterator nodeIt;
	std::map< MElementType, TPZVec < TPZVec < int > > >::iterator connIt;

	int nelements = fNodeData.size();
	Out << nelements << std::endl;      //número de tipos de elementos
	
	for ( nodeIt = fNodeData.begin(), connIt = fConnectivityData.begin();
				nodeIt != fNodeData.end() && connIt != fConnectivityData.end();
				nodeIt++, connIt++ )
	{
		Out << nodeIt->first << "\n"; //type
		TPZVec < TPZGeoNode > & nodeVec = (nodeIt->second);
		int nnodes = nodeVec.NElements();
		Out << nnodes << "\t";    //número de nós
		for (int i = 0; i < nnodes; i++) 
		{
			Out << nodeVec[i].Id() << "\t";         //node id
			for (int j = 0; j < 3; j++) 
			{
				Out << nodeVec[i].Coord(j) << "\t";           //node coordinates
			}
		}
		Out << std::endl;
		
		TPZVec < TPZVec < int > > & connectVec = (connIt->second);
		int nelements = connectVec.NElements();
		Out << nelements << "\n";       //número de elementos
		for (int i = 0; i < nelements; i++) 
		{
			TPZVec <int> & connectIdx =  connectVec[i];
			Out << i << "\t";                          //elindex
			int ncon = connectIdx.NElements();
			Out << ncon << "\t";                     //número de connectividades
			for (int j = 0; j < ncon; j++) 
			{
				Out << connectIdx[ j ] << "\t";       //conectividade i
			}
		}
	}
  Out << "\n";

  const int size = fConstScalarSolution.size();
  Out << size << "\n";
  std::map< std::string, TPZVec<REAL > >::const_iterator w;
  for(w = fConstScalarSolution.begin(); w != fConstScalarSolution.end(); w++){
    Out << w->first << "\n";
    Out << w->second.NElements() << "\t";
    for(int i = 0; i < w->second.NElements(); i++) Out << w->second[i] << "\t";
    Out << "\n";
  }///for w

}

bool TSWXGraphElement::Read (std::istream & In)
{
	int neltypes;
	int nnodes;
	int nelements;
	int nconnects;
	MElementType type;
	int nodeid;
	int elindex;
	TPZVec<REAL> coord(3,0.);
	In >> neltypes;
	TPZGeoMesh gMesh;

	for (int i=0; i < neltypes; i++)
	{
		int auxtype;
		In >> auxtype;
		type = (MElementType) auxtype;
		In >> nnodes;
		TPZVec < TPZGeoNode > nodeVec ( nnodes );
		for (int j=0; i < nnodes; j++)
		{
			In >> nodeid;
			for (int k = 0; k < 3; k++)
			{
				In >> coord[k];
			}
			TPZGeoNode node (nodeid,coord,gMesh);
			nodeVec[j] = node;
		}
		fNodeData[type] = nodeVec;

		In >> nelements;
		TPZVec < TPZVec <int> > connVec (nelements);
		for (int j = 0; j < nelements; j++)
		{
			In >> elindex >> nconnects;
			connVec[j].Resize(nconnects);
			for (int k = 0; k < nconnects; k++)
			{
				In >> connVec[j][k];
			}
		}
		fConnectivityData[type] = connVec;
	}

  fConstScalarSolution.clear();
  int size;
  In >> size;
  for(int imap = 0; imap < size; imap++){
    std::string name;
    In >> name;
    int nvec;
    In >> nvec;
    TPZVec<REAL> vec(nvec);
    for(int i = 0; i < nvec; i++) In >> vec[i];
    fConstScalarSolution[name] = vec;
  }///for w


  return true;
}

void TSWXGraphElement::GetData ( MElementType Type,
																 TPZVec<TPZGeoNode> & NodeVec,
																 TPZVec < TPZVec < int > > & ConnectivityVec,
                                 TPZVec< TPZGeoNode> &CenterOfSubEls )
{
	if ( fNodeData.find(Type) == fNodeData.end() ||
			 fConnectivityData.find(Type) == fConnectivityData.end() )
	{
		DebugStop();
		NodeVec.Resize(0);
		ConnectivityVec.Resize(0);
		return;
	}
	NodeVec = fNodeData [ Type ];
	ConnectivityVec = fConnectivityData [ Type ];

  int nsubEl = ConnectivityVec.NElements();
  CenterOfSubEls.Resize(nsubEl);
  for(int isubel = 0; isubel < nsubEl; isubel++){
    int size = ConnectivityVec[isubel].NElements();
    if(size){
      int node = ConnectivityVec[isubel][ size -1 ];
      CenterOfSubEls[isubel] = NodeVec[node];
    }
    else{
      DebugStop();
    }
  }

}



EGraphElType TSWXGraphElement::PZToVtkConnecitivities ( TPZVec < int > & ConnectivityVec)
{
  EGraphElType elType = ENoneGraphElType;

	int nsides = ConnectivityVec.NElements();
	//Acerto das connectividades para compatibilidade com o padrão do vtk
	switch (nsides)
	{
		case ( 1 ) :
		{
			//EPoint
			break;
		}
		case ( 3 ) :
		{
			//EOned -> vtkQuadraticEdge
      elType = EvtkQuadraticEdge;
			break;
		}
		case ( 7 ) :
		{
			//ETriangle -> vtkBiQuadraticTriangle
      elType = EvtkBiQuadraticTriangle;
			break;
		}
		case ( 9 ) :
		{
			//EQuadrilateral: -> vtkBiQuadraticQuad
      elType = EvtkBiQuadraticQuad;
			break;
		}
		case ( 15 ) :
		{ //ETetraedro -> vtkQuadraticTetra
      elType = EvtkQuadraticTetra;
			// não tem as faces nem o centro do volume
			ConnectivityVec.Resize ( 10 );
			break;
		}
		case ( 19 ) :
		{ //EPiramide -> vtkQuadraticPyramid
      elType = EvtkQuadraticPyramid;
			// não tem as faces nem o centro do volume
			ConnectivityVec.Resize ( 13 );
			break;
		}
		case ( 21 ) :
		{ //EPrisma -> vtkQuadraticWedge
      elType = EvtkQuadraticWedge;
			// não tem o centro das faces triangulares
			ConnectivityVec [ 16 ] = ConnectivityVec [ 17 ];
			ConnectivityVec [ 17 ] = ConnectivityVec [ 18 ];
			ConnectivityVec [ 18 ] = ConnectivityVec [ 19 ];
			ConnectivityVec.Resize ( 19 );

      ConnectivityVec.Resize(15);///se comentar essa linha, tem-se o EvtkBiQuadraticQuadraticWedge


      TPZVec < int > PZ(ConnectivityVec), VTK(15);
      for(int i = 0; i < 9; i++) VTK[i] = PZ[i];
      VTK[12]=PZ[9];
      VTK[13]=PZ[10];
      VTK[14]=PZ[11];
      VTK[9]=PZ[12];
      VTK[10]=PZ[13];
      VTK[11]=PZ[14];

      ConnectivityVec = VTK;


			break;
		}
		case ( 27 ) :
		{ // ECube -> EvtkQuadraticHexahedron
      elType = EvtkQuadraticHexahedron;
			double aux = ConnectivityVec [ 12 ];
			ConnectivityVec [ 12 ] = ConnectivityVec [ 16 ];
			ConnectivityVec [ 16 ] = aux;

			aux = ConnectivityVec [ 13 ];
			ConnectivityVec [ 13 ] = ConnectivityVec [ 17 ];
			ConnectivityVec [ 17 ] = aux;

			aux = ConnectivityVec [ 14 ];
			ConnectivityVec [ 14 ] = ConnectivityVec [ 18 ];
			ConnectivityVec [ 18 ] = aux;

			aux = ConnectivityVec [ 15 ];
			ConnectivityVec [ 15 ] = ConnectivityVec [ 19 ];
			ConnectivityVec [ 19 ] = aux;

			aux = ConnectivityVec [ 20 ];
			ConnectivityVec [ 20 ] = ConnectivityVec [ 24 ];
			ConnectivityVec [ 24 ] = aux;
			aux = ConnectivityVec [ 21 ];
			ConnectivityVec [ 21 ] = ConnectivityVec [ 22 ];
			ConnectivityVec [ 22 ] = aux;
      ConnectivityVec.Resize(20);///se comentar essa linha, tem-se o TriQuadraticHexahedron
			break;
		}
		default:
		{
			DebugStop();
			break;
		}
	}
	return elType;
}


void TSWXGraphElement::ProcessAllTopologies( int nRef )
{
	{
		MElementType type = EOned;
		TPZVec <TPZGeoNode> AllNodes;
		TPZVec < TPZVec< int > > Connectivities;
		ProcessElement ( type, nRef, AllNodes, Connectivities);
		fNodeData[ type ] = AllNodes;
		fConnectivityData[ type ] = Connectivities;
	}
	{
		MElementType type = ETriangle;
		TPZVec <TPZGeoNode> AllNodes;
		TPZVec < TPZVec< int > > Connectivities;
		ProcessElement ( type, nRef, AllNodes, Connectivities);
		fNodeData[ type ] = AllNodes;
		fConnectivityData[ type ] = Connectivities;
	}
	{
		MElementType type = EQuadrilateral;
		TPZVec <TPZGeoNode> AllNodes;
		TPZVec < TPZVec< int > > Connectivities;
		ProcessElement ( type, nRef, AllNodes, Connectivities);
		fNodeData[ type ] = AllNodes;
		fConnectivityData[ type ] = Connectivities;
	}
	{
		MElementType type = ETetraedro;
		TPZVec <TPZGeoNode> AllNodes;
		TPZVec < TPZVec< int > > Connectivities;
		ProcessElement ( type, nRef, AllNodes, Connectivities);
		fNodeData[ type ] = AllNodes;
		fConnectivityData[ type ] = Connectivities;
	}
	{
		MElementType type = EPiramide;
		TPZVec <TPZGeoNode> AllNodes;
		TPZVec < TPZVec< int > > Connectivities;
		ProcessElement ( type, nRef, AllNodes, Connectivities);
		fNodeData[ type ] = AllNodes;
		fConnectivityData[ type ] = Connectivities;
	}
	{
		MElementType type = EPrisma;
		TPZVec <TPZGeoNode> AllNodes;
		TPZVec < TPZVec< int > > Connectivities;
		ProcessElement ( type, nRef, AllNodes, Connectivities);
		fNodeData[ type ] = AllNodes;
		fConnectivityData[ type ] = Connectivities;
	}
	{
		MElementType type = ECube;
		TPZVec <TPZGeoNode> AllNodes;
		TPZVec < TPZVec< int > > Connectivities;
		ProcessElement ( type, nRef, AllNodes, Connectivities);
		fNodeData[ type ] = AllNodes;
		fConnectivityData[ type ] = Connectivities;
	}
}

void TSWXGraphElement::ProcessElement ( MElementType Type, int nRef,
																				TPZVec <TPZGeoNode> & AllNodes,
																				TPZVec < TPZVec< int > > & Connectivities )
{
	//Malha auxiliar
	TPZGeoMesh gMesh;

	//criando um elemento "pai" com coordenas espaciais iguais às coordenadas do
	//seu element mestre
	CreateElement(gMesh, Type);
	
	if (gMesh.NElements() != 1) 
	{
		DebugStop();
		return;
	}

	//Refinar o elemento pelo número de vezes especificado
	Divide ( gMesh, nRef );

	int i;

	//Extrair apenas os elementos "mais refinados" (sem filhos)
	int nel = gMesh.NElements();
	TPZStack <TPZGeoEl *> dividedStack;
	for (i = 0; i < nel; i++)
	{
		TPZGeoEl * gel = gMesh.ElementVec()[i];
		if ( !gel || gel->HasSubElement() ) 
		{
			continue;
		}
		dividedStack.Push(gel);
	}
	int nsubel = dividedStack.NElements();
	Connectivities.Resize(nsubel);

	//Criar a lista de conectividades para elementos triquadráticos do vtk
	for (i=0;i<nsubel;i++)
	{
		TPZGeoEl *gel = dividedStack.Pop();
		if (!gel) continue;
		int nsides = gel->NSides();
		
		Connectivities[i].Resize(nsides);
		int nnodes = gMesh.NodeVec().NElements();
		int ncorner = gel->NCornerNodes();
		//inserindo conectividades dos nós de canto
		for (int j = 0; j < ncorner; j++)
		{
			Connectivities[i][j] = gel->NodeIndex(j);
		}
		
		gMesh.NodeVec().Resize( nnodes + nsides - ncorner );
		//criando nós no centro dos lados e inserindo os seus indexes na lista de conectividades
		for (int j = ncorner; j < nsides; j++)
		{
			TPZVec <REAL> centerPtMaster(3,0.);
			TPZVec <REAL> centerPtSpace(3,0.);
			gel->CenterPoint(j,centerPtMaster);
			gel->X(centerPtMaster,centerPtSpace);
			int nodeIdx = nnodes + j -ncorner;
			gMesh.NodeVec()[nodeIdx].Initialize ( centerPtSpace, gMesh );
			Connectivities[i][j] = nodeIdx;
		}
	}

	//copiando o vetor de nós (id versus coordenadas)
	int nnodes = gMesh.NodeVec().NElements();
	AllNodes.Resize ( nnodes );
	for (i = 0; i < nnodes; i++) 
	{
		AllNodes[i] = gMesh.NodeVec()[i];
	}
}

void TSWXGraphElement::Divide ( TPZGeoMesh &gMesh, int nRef )
{
	TPZGeoEl *gel = gMesh.ElementVec()[0];
	TPZStack < TPZGeoEl * > gelToDivide;
	gelToDivide.Push(gel);

	TPZVec<TPZGeoEl *> subgelVec;

	while( gelToDivide.NElements() )
	{
		TPZGeoEl *gelref = gelToDivide.Pop();
		if (!gelref|| gelref->Level() >= (nRef) )
		{
			continue;
		}
		gelref->Divide(subgelVec);
		for (int iel = 0; iel < subgelVec.NElements(); iel++) 
		{
			TPZGeoEl *son = subgelVec[iel];
			if (!son) 
			{
				DebugStop();
				continue;
			}
			gelToDivide.Push(son);
		}
		subgelVec.Resize(0);
	}
}

void TSWXGraphElement::CreateElement ( TPZGeoMesh & gMesh, MElementType & type )
{
	int matid = 37466155;
	TPZVec<REAL> coord(3,0.);
	TPZVec<int64_t> cornerIdxVec;
	
	switch ( type ) 
	{
		case ( EPoint ):
		{
			int nnodes = 1;
			gMesh.NodeVec().Resize(nnodes);
			cornerIdxVec.Resize(nnodes);
			
			gMesh.NodeVec()[0].Initialize(coord,gMesh);
			cornerIdxVec[0] = 0;
			break;
		}
		case ( EOned ):
		{
			//dois nos, -1, +1
			int nnodes = 2;
			gMesh.NodeVec().Resize(nnodes);
			cornerIdxVec.Resize(nnodes);

			coord[0] = -1.;
			cornerIdxVec[0] = 0;
			gMesh.NodeVec()[0].Initialize(coord,gMesh);
			coord[0] = 1.;
			cornerIdxVec[1] = 1;
			gMesh.NodeVec()[1].Initialize(coord,gMesh);
			break;
		}
		case ( ETriangle ) :
		{
		//três nós, (0,0), (1,0), (0,1)
			int nnodes = 3;
			gMesh.NodeVec().Resize(nnodes);
			cornerIdxVec.Resize(nnodes);

			coord[0] = 0.;
			coord[1] = 0.;
			cornerIdxVec[0] = 0;
			gMesh.NodeVec()[0].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = 0.;
			cornerIdxVec[1] = 1;
			gMesh.NodeVec()[1].Initialize(coord,gMesh);
			coord[0] = 0.;
			coord[1] = 1.;
			cornerIdxVec[2] = 2;
			gMesh.NodeVec()[2].Initialize(coord,gMesh);
			break;
		}
		case (EQuadrilateral):
		{
		//quatro nós, (-1,-1), (1,-1), (1,1), (-1,1)
			int nnodes = 4;
			gMesh.NodeVec().Resize(nnodes);
			cornerIdxVec.Resize(nnodes);

			coord[0] = -1.;
			coord[1] = -1.;
			cornerIdxVec[0] = 0;
			gMesh.NodeVec()[0].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = -1.;
			cornerIdxVec[1] = 1;
			gMesh.NodeVec()[1].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = 1.;
			cornerIdxVec[2] = 2;
			gMesh.NodeVec()[2].Initialize(coord,gMesh);
			coord[0] = -1.;
			coord[1] = 1.;
			cornerIdxVec[3] = 3;
			gMesh.NodeVec()[3].Initialize(coord,gMesh);
			break;
		}

		case (ETetraedro) :
		{
		//quatro nós, (0,0,0), (1,0,0), (0,0,0), (0,0,1)
			int nnodes = 4;
			gMesh.NodeVec().Resize(nnodes);
			cornerIdxVec.Resize(nnodes);

			coord[0] = 0.;
			coord[1] = 0.;
			coord[2] = 0.;
			cornerIdxVec[0] = 0;
			gMesh.NodeVec()[0].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = 0.;
			coord[2] = 0.;
			cornerIdxVec[1] = 1;
			gMesh.NodeVec()[1].Initialize(coord,gMesh);
			coord[0] = 0.;
			coord[1] = 1.;
			coord[2] = 0.;
			cornerIdxVec[2] = 2;
			gMesh.NodeVec()[2].Initialize(coord,gMesh);
			coord[0] = 0.;
			coord[1] = 0.;
			coord[2] = 1.;
			cornerIdxVec[3] = 3;
			gMesh.NodeVec()[3].Initialize(coord,gMesh);
			break;
		}

		case (EPiramide) :
		{
		//cinco nós, (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,1)
			int nnodes = 5;
			gMesh.NodeVec().Resize(nnodes);
			cornerIdxVec.Resize(nnodes);

			coord[0] = -1.;
			coord[1] = -1.;
			coord[2] = 0.;
			cornerIdxVec[0] = 0;
			gMesh.NodeVec()[0].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = -1.;
			coord[2] = 0.;
			cornerIdxVec[1] = 1;
			gMesh.NodeVec()[1].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = 1.;
			coord[2] = 0.;
			cornerIdxVec[2] = 2;
			gMesh.NodeVec()[2].Initialize(coord,gMesh);
			coord[0] = -1.;
			coord[1] = 1.;
			coord[2] = 0.;
			cornerIdxVec[3] = 3;
			gMesh.NodeVec()[3].Initialize(coord,gMesh);
			coord[0] = 0.;
			coord[1] = 0.;
			coord[2] = 1.;
			cornerIdxVec[4] = 4;
			gMesh.NodeVec()[4].Initialize(coord,gMesh);
			break;
		}

		case ( EPrisma ):
		{
		//seis nós, (0,0,-1), (1,0,-1), (0,1,-1),(0,0,1), (1,0,1), (0,1,1)
			int nnodes = 6;
			gMesh.NodeVec().Resize(nnodes);
			cornerIdxVec.Resize(nnodes);

			coord[0] = 0.;
			coord[1] = 0.;
			coord[2] = -1.;
			cornerIdxVec[0] = 0;
			gMesh.NodeVec()[0].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = 0.;
			coord[2] = -1.;
			cornerIdxVec[1] = 1;
			gMesh.NodeVec()[1].Initialize(coord,gMesh);
			coord[0] = 0.;
			coord[1] = 1.;
			coord[2] = -1.;
			cornerIdxVec[2] = 2;
			gMesh.NodeVec()[2].Initialize(coord,gMesh);
			coord[0] = 0.;
			coord[1] = 0.;
			coord[2] = 1.;
			cornerIdxVec[3] = 3;
			gMesh.NodeVec()[3].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = 0.;
			coord[2] = 1.;
			cornerIdxVec[4] = 4;
			gMesh.NodeVec()[4].Initialize(coord,gMesh);
			coord[0] = 0.;
			coord[1] = 1.;
			coord[2] = 1.;
			cornerIdxVec[5] = 5;
			gMesh.NodeVec()[5].Initialize(coord,gMesh);
			break;
		}

		case (ECube):
		{
		//oito nós, (-1,-1,-1), (1,-1,-1), (1,1,-1),(-1,1,-1),
		//          (-1,-1,1), (1,-1,1), (1,1,1),(-1,1,1)
			int nnodes = 8;
			gMesh.NodeVec().Resize(nnodes);
			cornerIdxVec.Resize(nnodes);

			coord[0] = -1.;
			coord[1] = -1.;
			coord[2] = -1.;
			cornerIdxVec[0] = 0;
			gMesh.NodeVec()[0].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = -1.;
			coord[2] = -1.;
			cornerIdxVec[1] = 1;
			gMesh.NodeVec()[1].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = 1.;
			coord[2] = -1.;
			cornerIdxVec[2] = 2;
			gMesh.NodeVec()[2].Initialize(coord,gMesh);
			coord[0] = -1.;
			coord[1] = 1.;
			coord[2] = -1.;
			cornerIdxVec[3] = 3;
			gMesh.NodeVec()[3].Initialize(coord,gMesh);
			coord[0] = -1.;
			coord[1] = -1.;
			coord[2] = 1.;
			cornerIdxVec[4] = 4;
			gMesh.NodeVec()[4].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = -1.;
			coord[2] = 1.;
			cornerIdxVec[5] = 5;
			gMesh.NodeVec()[5].Initialize(coord,gMesh);
			coord[0] = 1.;
			coord[1] = 1.;
			coord[2] = 1.;
			cornerIdxVec[6] = 6;
			gMesh.NodeVec()[6].Initialize(coord,gMesh);
			coord[0] = -1.;
			coord[1] = 1.;
			coord[2] = 1.;
			cornerIdxVec[7] = 7;
			gMesh.NodeVec()[7].Initialize(coord,gMesh);
			break;
		}
		default:
		{
			DebugStop();
			break;
		}
	}
	int64_t index = 0;
	gMesh.CreateGeoElement(type,cornerIdxVec,matid,index,0);
}

void TSWXGraphElement::GetSolution(TPZInterpolationSpace *sp,
                                   std::string varName,
                                   TPZVec< TPZVec< REAL > > &ParametricCoord,
                                   TPZVec< TSWXGraphSingleSol > &locSol){

  if(fConstScalarSolution.find(varName) != fConstScalarSolution.end()){
    const int npts = ParametricCoord.NElements();
    const int dimSol = 1;///scalarSolution apenas
    locSol.resize(npts);
    for(int i = 0; i < npts; i++){
      locSol[i].resize(dimSol);
      locSol[i][0] = fConstScalarSolution[varName].operator[](sp->Index());
    }
    return;
  }


  if(!sp->Material()) DebugStop();
  const int varIndex = sp->Material()->VariableIndex(varName);
  const int dimSol = sp->Material()->NSolutionVariables(varIndex);

  const int npts = ParametricCoord.NElements();
  locSol.resize(npts);
  for(int i = 0; i < npts; i++){
    TPZManVector<REAL> solution(dimSol);
    sp->Solution(ParametricCoord[i], varIndex, solution);
    locSol[i].resize(solution.NElements());
    for(int j = 0; j < solution.NElements(); j++){
      locSol[i][j] = solution[j];
    }
  }
}

void TSWXGraphElement::GenerateVTKData(TPZCompMesh * cmesh,
                                       int dim,
                                       double time,
                                       const TPZVec<std::string> & nodalVarIndex,
                                       const TPZVec<std::string> & cellVarIndex,
                                       TSWXGraphMesh & graphMesh){
  std::set< int > mySetMatIds;
  this->GenerateVTKData(cmesh,dim,time,nodalVarIndex,cellVarIndex,graphMesh,mySetMatIds);
}

void TSWXGraphElement::GenerateVTKData(TPZCompMesh * cmesh,
                                       int dim,
                                       double time,
                                       const TPZVec<std::string> & nodalVarIndex,
                                       const TPZVec<std::string> & cellVarIndex,
                                       TSWXGraphMesh & graphMesh,
                                       const std::set< int > &mySetMatIds){
  if(!cmesh) return;
  TPZVec< TSWXGraphSol > sol;
  bool updateMesh = graphMesh.fNodes.size() == 0 ? true : false;   //tamarindo
  this->GenerateVTKData(cmesh,dim,nodalVarIndex,cellVarIndex,graphMesh.fNodes,
                        graphMesh.fElem,sol,
                        mySetMatIds,updateMesh);
  graphMesh.fSol.AddSol(time,sol);
}


void TSWXGraphElement::CollectElements(TPZCompMesh *cmesh, std::set< TPZInterpolationSpace* > &AllEls){
  if(!cmesh) return;
  const int nel = cmesh->NElements();
  for(int iel = 0; iel < nel; iel++){
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    if(!cel) continue;

    TPZSubCompMesh * subcmesh = dynamic_cast< TPZSubCompMesh* >(cel);
    if(subcmesh){
      this->CollectElements(subcmesh,AllEls);
      continue;
    }

    TPZGeoEl * gel =cel->Reference();
    if(!gel) continue;
    TPZInterpolationSpace * sp = dynamic_cast<TPZInterpolationSpace*>(cel);
    if(!sp) continue;
    if(sp->NConnects() == 0) continue;
    if(!sp->Material()) continue;

    AllEls.insert(sp);
  }
}

void TSWXGraphElement::GenerateVTKData(TPZCompMesh * cmesh,
                                       int dim,
                                       const TPZVec<std::string> & nodalVarIndex,
                                       const TPZVec<std::string> & cellVarIndex,
                                       TPZStack< TSWXGraphNode > &nodes,
                                       TPZStack< TSWXGraphEl > &elem,
                                       TPZVec< TSWXGraphSol > &sol,
                                       const std::set< int > &mySetMatIds,
                                       bool updateMesh){
  if(!cmesh) return;
  if(updateMesh){
    nodes.resize(0);
    elem.resize(0);
  }
  sol.resize(0);
  const int nnodalBasedSol = nodalVarIndex.NElements();
  const int ncellBasedSol = cellVarIndex.NElements();
  sol.resize(nnodalBasedSol+ncellBasedSol);
  for(int i = 0; i < nnodalBasedSol; i++){
    sol[i].fSolType = ENodeSolution;
    sol[i].fTitle = nodalVarIndex[i];
  }
  for(int i = 0; i < ncellBasedSol; i++){
    sol[nnodalBasedSol+i].fSolType = ECellSolution;
    sol[nnodalBasedSol+i].fTitle = cellVarIndex[i];
  }

  TPZVec< TPZVec<REAL> > ParametricNodes;
  TPZVec< TPZVec<int> > ParametricIncid;
  TPZVec< TPZGeoNode> CenterOfSubEls;
  TPZVec< TPZVec<REAL> > ParametricCenterNodes;

  std::set< TPZInterpolationSpace* > AllEls;
  std::set< TPZInterpolationSpace* >::iterator w;

  this->CollectElements(cmesh, AllEls);

  for(w = AllEls.begin(); w != AllEls.end(); w++){

    TPZInterpolationSpace *sp = *w;
    if(!sp) continue;

    TPZGeoEl * gel = sp->Reference();
    if(gel->Dimension() != dim) continue;

    if(mySetMatIds.size()){
        if(mySetMatIds.find( sp->Material()->Id() ) == mySetMatIds.end()) continue;
    }

    {
      ///nodal
      TPZManVector<TPZGeoNode,1000> NodeVec;
      this->GetData(gel->Type(), NodeVec, ParametricIncid, CenterOfSubEls);
      int nCreatedNodes = NodeVec.NElements();
      ParametricNodes.Resize(nCreatedNodes);
      for(int in = 0; in < nCreatedNodes; in++){
        ParametricNodes[in].Resize(3);
        NodeVec[in].GetCoordinates( ParametricNodes[in] );
      }

      ///centers
      nCreatedNodes = CenterOfSubEls.NElements();
      ParametricCenterNodes.Resize(nCreatedNodes);
      for(int in = 0; in < nCreatedNodes; in++){
        ParametricCenterNodes[in].Resize(3);
        CenterOfSubEls[in].GetCoordinates( ParametricCenterNodes[in] );
      }
    }

    for( int inodalBaseSol = 0; inodalBaseSol < nnodalBasedSol; inodalBaseSol++){
      TPZVec< TSWXGraphSingleSol > locSol;
      ///solucao nodal
      GetSolution(sp, nodalVarIndex[inodalBaseSol], ParametricNodes, locSol);
      for(unsigned int isol = 0; isol < locSol.size(); isol++){
        sol[inodalBaseSol].fData.push_back(locSol[isol]);
      }
    }

    for( int icellBaseSol = 0; icellBaseSol < ncellBasedSol; icellBaseSol++){
      TPZVec< TSWXGraphSingleSol > locSol;
      ///solucao cell
      GetSolution(sp, cellVarIndex[icellBaseSol], ParametricCenterNodes, locSol);
      for(unsigned int isol = 0; isol < locSol.size(); isol++){
        sol[icellBaseSol+nnodalBasedSol].fData.push_back(locSol[isol]);
      }
    }

    if(updateMesh){
      ///nos
      int Node0 = nodes.size();
      for(int inode = 0; inode < ParametricNodes.size(); inode++){
        TPZManVector< REAL, 3> XVal(3);
        gel->X(ParametricNodes[inode], XVal);
        TSWXGraphNode newNode(3);
        newNode[0] = XVal[0]; newNode[1] = XVal[1]; newNode[2] = XVal[2];
        nodes.push_back(newNode);
      }///inode

      ///incidencia
      for(int i = 0; i < ParametricIncid.NElements(); i++){

        EGraphElType elType = this->PZToVtkConnecitivities( ParametricIncid[i] );

        const int nnodesEl = ParametricIncid[i].NElements();
        TSWXGraphEl locEl;
        locEl.fElType = elType;
        locEl.fMatId = sp->Material()->Id();
        locEl.fIncid.resize(nnodesEl);
        for(int j = 0; j < nnodesEl; j++){
          locEl.fIncid[j] = ParametricIncid[i][j] + Node0;
        }///j

        elem.push_back(locEl);

      }///i
    }

  }///for iel

}///void

void TSWXGraphElement::SetConstScalarSolution(const std::map< std::string, TPZVec<REAL > > &sol){
  fConstScalarSolution = sol;
}


