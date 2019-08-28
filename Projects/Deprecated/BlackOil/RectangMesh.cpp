TPZGeoMesh * QuarterFiveSpotReg(double lxy, double lz)
{	
  TPZGeoMesh * Mesh = new TPZGeoMesh;

  int Qnodes = 484;
  TPZVec < TPZVec <REAL> > NodesCoords(Qnodes);
  for(int i = 0; i < Qnodes; i++) NodesCoords[i].Resize(3);

  int QTDnodesByLayer = 121;

  for(int lay = 0; lay < 4; lay++) // Begin NodesGeneration by Layer
  {
    int firstNode = QTDnodesByLayer * lay; //First node ID on this layer

    for(int nx = 0; nx < 11; nx++)
    {
        for(int ny = 0; ny < 11; ny++)
        {
            NodesCoords[11*nx+ny + firstNode][0] = -lxy/2. + lxy/10.*nx;
            NodesCoords[11*nx+ny + firstNode][1] = -lxy/2. + lxy/10.*ny;
            NodesCoords[11*nx+ny + firstNode][2] = -lz/2. + lz/3.*lay;
        }
    }
  }

  Mesh->NodeVec().Resize(Qnodes);
  TPZVec <TPZGeoNode> Node(Qnodes);
  for(int n = 0; n < Qnodes; n++)
  {
      Node[n].SetNodeId(n);
      Node[n].SetCoord(&NodesCoords[n][0]);
      Mesh->NodeVec()[n] = Node[n];
  }

  TPZVec <int> Topol;

  // Materials IDs
  int reservMat = 1;
  int InjectionBC = -1;
  int ProductionBC = -2;
  int arcMat = -3;

  int id = 0;

  for(int lay = 0; lay < 3; lay++) // Begin NodesGeneration by Layer
  {
    int firstNode = QTDnodesByLayer * lay; //First node ID on this layer
    Topol.Resize(8);
    for(int ny = 0; ny < 10; ny++)
    {
        for(int nx = 0; nx < 10; nx++)
        {
            if((nx != 0 || ny != 0) && (nx != 9 || ny != 9))
            {
                ///Reservoir Cubes
                Topol[0] = 11*ny+nx + firstNode;
                Topol[1] = 11*ny+nx+1 + firstNode;
                Topol[2] = 11*(ny+1)+nx+1 + firstNode;
                Topol[3] = 11*(ny+1)+nx + firstNode;
                Topol[4] = 11*ny+nx + firstNode + QTDnodesByLayer;
                Topol[5] = 11*ny+nx+1 + firstNode + QTDnodesByLayer ;
                Topol[6] = 11*(ny+1)+nx+1 + firstNode + QTDnodesByLayer ;
                Topol[7] = 11*(ny+1)+nx + firstNode + QTDnodesByLayer;
                new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,reservMat,*Mesh);
                id++;
            }
        }
    }

    ///Quadrilaterals InjectionBC
    Topol.Resize(4);
    Topol[0] = firstNode + 12;
    Topol[1] = firstNode + 1;
    Topol[2] = firstNode + 1  + QTDnodesByLayer;
    Topol[3] = firstNode + 12 + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,InjectionBC,*Mesh);
    id++;

    Topol[0] = firstNode + 11;
    Topol[1] = firstNode + 12;
    Topol[2] = firstNode + 12  + QTDnodesByLayer;
    Topol[3] = firstNode + 11 + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,InjectionBC,*Mesh);
    id++;

    ///Quadrilaterals ProductionBC
    Topol[0] = firstNode + 109;
    Topol[1] = firstNode + 108;
    Topol[2] = firstNode + 108 + QTDnodesByLayer;
    Topol[3] = firstNode + 109 + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ProductionBC,*Mesh);
    id++;

    Topol[0] = firstNode + 108;
    Topol[1] = firstNode + 119;
    Topol[2] = firstNode + 119 + QTDnodesByLayer;
    Topol[3] = firstNode + 108 + QTDnodesByLayer;
    new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ProductionBC,*Mesh);
    id++;
  }


  Mesh->BuildConnectivity();

#ifdef LOG4CXX
  {
    std::stringstream sout;
    Mesh->Print(sout);
    LOGPZ_DEBUG(logger,sout.str());
  }
#endif

  return Mesh;
}
