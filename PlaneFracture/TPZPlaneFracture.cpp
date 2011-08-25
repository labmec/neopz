/**
 * \file
 * @brief Contains implementations of the TPZPlaneFracture methods.
 */
/*
 *  TPZPlaneFracture.cpp
 *  Crack
 *
 *  Created by Cesar Lucci on 09/08/10.
 *  Copyright 2010 LabMeC. All rights reserved.
 *
 */

#include "TPZPlaneFracture.h"
#include "TPZGeoLinear.h"
#include "pznoderep.h"
#include "pzgeoel.h"
#include "TPZVTKGeoMesh.h"
#include "TPZPoligonalChain.h"
#include "tpzgeoelrefpattern.h"

using namespace pzgeom;

const int __fractureRefinedEl_Mat = -999;
const int __fractureLine_Mat = 999;

const double __smallNum = 1.E-10;

TPZPlaneFracture::TPZPlaneFracture(TPZGeoMesh * planeMesh, int nodeOrigin, int nodeX, int nodeY, int TrimQTD)
{
	fplaneMesh = planeMesh;
	
	TPZFMatrix localAxes(3,2,0.);
	REAL valFrom, valTo;
	for(int c = 0; c < 3; c++)
	{
		valFrom = planeMesh->NodeVec()[nodeOrigin].Coord(c);
		
		//first column
		valTo = planeMesh->NodeVec()[nodeX].Coord(c);
		localAxes.PutVal(c, 0, valTo-valFrom);
		
		//second column
		valTo = planeMesh->NodeVec()[nodeY].Coord(c);
		localAxes.PutVal(c, 1, valTo-valFrom);
	}
	
	TPZFMatrix notUsedHere;
	localAxes.GramSchmidt(fFromR3toR2, notUsedHere);
	fFromR3toR2.Transpose();
	
	fTrimQTD = (TrimQTD > __minTrimQTD) ? TrimQTD : __minTrimQTD;
}

TPZPlaneFracture::~TPZPlaneFracture()
{
	
}

TPZGeoMesh * TPZPlaneFracture::GetFractureMesh(TPZVec<REAL> &poligonalChain)
{
	#ifdef DEBUG
	int ncoord = poligonalChain.NElements();
	if(ncoord%3 != 0)
	{
		std::cout << "poligonalChain boundary dont have groups of 3 coordinates (x,y,z)!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	#endif
	
	TPZGeoMesh * fractMesh = new TPZGeoMesh(*fplaneMesh);
	int nelem = fractMesh->NElements();
	
	std::map< int, std::set<double> > elId_TrimCoords;
	std::list< std::pair<int,double> > elIdSequence;
	
	DetectEdgesCrossed(poligonalChain, fractMesh, elId_TrimCoords, elIdSequence);
	
//	{	// Monotone Chains ... (nao sei se farei!!!)
//		TPZVec<REAL> poligonalChainUpdated(0);
//		UpdatePoligonalChain(fractMesh, elIdSequence, poligonalChainUpdated);
//		
//		TPZPoligonalChain thisPChain;
//		thisPChain.SplitInMonotoneChains(poligonalChainUpdated, this->fFromR3toR2);	
//	}
	
	//Refining 1D elements
	TPZVec<TPZGeoEl*> sons;
	std::map< int, std::set<double> >::iterator it;
	for(it = elId_TrimCoords.begin(); it != elId_TrimCoords.end(); it++)
	{
		int el1DId = it->first;
		TPZAutoPointer<TPZRefPattern> linRefp = Generate1DRefPatt(it->second);
		fractMesh->ElementVec()[el1DId]->SetRefPattern(linRefp);
		fractMesh->ElementVec()[el1DId]->Divide(sons);
	}
	
	//Refining 2D elements
	for(int el = 0; el < nelem; el++)
	{
		TPZGeoEl * gel = fractMesh->ElementVec()[el];
		TPZAutoPointer<TPZRefPattern> elRefp = TPZRefPatternTools::PerfectMatchRefPattern(gel);
		if(elRefp)
		{
			gel->SetRefPattern(elRefp);
			gel->Divide(sons);
		}
	}
	
	GenerateCrackBoundary(fractMesh, elIdSequence);
	
	fractMesh->BuildConnectivity();
	
	return fractMesh;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::DetectEdgesCrossed(TPZVec<REAL> &poligonalChain, TPZGeoMesh * fractMesh, std::map< int,
										  std::set<double> > &elId_TrimCoords, std::list< std::pair<int,double> > &elIdSequence)
{
	int npoints = (poligonalChain.NElements())/3;
	int nelem = fractMesh->NElements();
	
	TPZGeoEl * gel = PointElement(0, fractMesh, poligonalChain);
	TPZGeoEl * nextGel = NULL;
	
	if(!gel)
	{
		std::cout << "first point of crack tip boundary does NOT belong to any 2D element" << std::endl;
		std::cout << "(or was given a mesh with zero quantity of 2D elements)!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	
	double alphaMin;
	bool reachNextPoint;
	int nElsCrossed, thispoint, nextpoint;
	std::map< int, std::set<double> >::iterator it;
	std::set<double> trim;
	TPZVec< TPZVec<REAL> > intersectionPoint;
	TPZFMatrix x(3,1), dx(3,1);
	TPZVec<REAL> xNext(3), qsi2D(2);
	TPZVec <int> Topol(2), edgeVec;
	
	int p;
	for(p = 0; p < (npoints-1); p++)
	{
		nElsCrossed = 0;
		alphaMin = 0.;
		thispoint = 3*p;
		nextpoint = 3*(p+1);
		for(int c = 0; c < 3; c++)
		{
			x(c,0) = poligonalChain[thispoint+c];
			dx(c,0) = poligonalChain[nextpoint+c] - poligonalChain[thispoint+c];
			xNext[c] = poligonalChain[nextpoint+c];
		}
		
		double norm = 0.;
		for(int c = 0; c < 3; c++)
		{
			norm += dx(c,0)*dx(c,0);
		}
		norm = sqrt(norm);		
		for(int c = 0; c < 3; c++)
		{
			dx(c,0) = dx(c,0)/norm;
		}
		
		reachNextPoint = gel->ComputeXInverse(xNext, qsi2D);
		while(reachNextPoint == false && nElsCrossed < nelem)
		{
			nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, true);
			
			alphaMin = __smallNum;//qdo vai de vizinho a vizinho ateh chegar no proximo ponto, nao deve-se incluir os (alphaX=0)
			if(nextGel->ComputeXInverse(xNext, qsi2D))
			{
				reachNextPoint = true;
			}
			gel = nextGel;
			nElsCrossed++;
		}
		if(nElsCrossed == nelem)
		{
			//Deve ter alternado entre vizinhos!!!
			DebugStop();
		}
	}
	
	for(int c = 0; c < 3; c++)
	{
		dx(c,0) = -fFromR3toR2(0,c);//direcao oposta ao eixo x do sistema de coordenadas
	}
	
	p = 0;//Fechando o inicio da fratura
	thispoint = 3*p;
	nextpoint = 3*(p+1);
	for(int c = 0; c < 3; c++)
	{
		x(c,0) = poligonalChain[thispoint+c];
	}
	nextGel = PointElement(p, fractMesh, poligonalChain);
	gel = NULL;
	alphaMin = 0.;
	while(gel != nextGel)
	{
		gel = nextGel;
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, false);
		alphaMin = __smallNum;
	}
	
	p = npoints-1;//Fechando o final da fratura
	thispoint = 3*p;
	nextpoint = 3*(p+1);
	for(int c = 0; c < 3; c++)
	{
		x(c,0) = poligonalChain[thispoint+c];
	}
	nextGel = PointElement(p, fractMesh, poligonalChain);
	gel = NULL;
	alphaMin = 0.;
	while(gel != nextGel)
	{
		gel = nextGel;
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, true);
		alphaMin = __smallNum;
	}
}

TPZGeoEl * TPZPlaneFracture::CrossToNextNeighbour(TPZGeoEl * gel, TPZFMatrix &x, TPZFMatrix dx, double alphaMin,
												  std::map< int, std::set<double> > &elId_TrimCoords, std::list< std::pair<int,double> > &elIdSequence, bool pushback)
{
	bool thereIsAn1DElemAlready;
	int edge;
	std::map< int, std::set<double> >::iterator it;
	std::set<double> trim;
	TPZVec< TPZFMatrix > ExactIntersectionPoint, ModulatedIntersectionPoint;
	TPZVec<REAL> qsi1Dvec(1), xCrackBoundary(3);
	TPZVec<int> Topol(2), edgeVec;
	
	TPZGeoMesh * fractMesh = gel->Mesh();
	
	bool haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint, alphaMin);
	
	if(haveIntersection == false)
	{
		haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint, 0.);
		if(haveIntersection == false)
		{
			std::cout << "The point inside element does NOT intersect its edges! EdgeIntersection method face an exeption!" << std::endl;
			std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
			
			std::cout << "Element " << gel->Id() << std::endl;
			for(int n = 0; n < gel->NNodes(); n++)
			{
				std::cout << "Node " << n << std::endl;
				gel->NodePtr(n)->Print(std::cout);
				std::cout << std::endl;
			}
			std::cout << "x: " << x(0,0) << " , " << x(1,0) << " , " << x(2,0) << std::endl;
			std::cout << "dx: " << dx(0,0) << " , " << dx(1,0) << " , " << dx(2,0) << std::endl;
			std::cout << "alphaMin = " << alphaMin << std::endl << std::endl;
			
			DebugStop();
		}
	}
	TPZVec<REAL> xLin(3);
	for(int edg = 0; edg < edgeVec.NElements(); edg++)
	{
		edge = edgeVec[edg];
		TPZVec<REAL> n0(3), n1(3);
		fractMesh->NodeVec()[gel->SideNodeIndex(edge, 0)].GetCoordinates(n0);
		fractMesh->NodeVec()[gel->SideNodeIndex(edge, 1)].GetCoordinates(n1);
		for(int c = 0; c < 3; c++)
		{
			xLin[c] = ModulatedIntersectionPoint[edg](c,0);
		}
		double qsi1D = LinearComputeXInverse(xLin, n0, n1);
		
		TPZGeoElSide gelEdge(gel, edge);
		TPZGeoElSide neighEdge = gelEdge.Neighbour();
		thereIsAn1DElemAlready = false;
		while(neighEdge != gelEdge)
		{
			if(neighEdge.Element()->Dimension() == 1)//jah existe um elemento 1D inserido nesta aresta!
			{
				thereIsAn1DElemAlready = true;
				int neighEdgeId = neighEdge.Element()->Id();
				it = elId_TrimCoords.find(neighEdgeId);
				TPZTransform transBetweenNeigh = neighEdge.NeighbourSideTransform(gelEdge);
				qsi1D *= transBetweenNeigh.Mult()(0,0);
				it->second.insert(qsi1D);
				if(pushback) // push_BACK
				{
					elIdSequence.push_back(std::make_pair(neighEdgeId, qsi1D));
				}
				else // push_FRONT
				{
					elIdSequence.push_front(std::make_pair(neighEdgeId, qsi1D));
				}
				//
				break;
			}
			neighEdge = neighEdge.Neighbour();
		}
		if(thereIsAn1DElemAlready == false)//nao existe um elemento 1D nesta aresta!
		{
			trim.clear();
			trim.insert(qsi1D);
			Topol[0] = gel->SideNodeIndex(edge, 0);
			Topol[1] = gel->SideNodeIndex(edge, 1);
			TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * linGeo = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, __fractureRefinedEl_Mat, *fractMesh);
			fractMesh->BuildConnectivity();
			
			int linGeoId = linGeo->Id();
			elId_TrimCoords[linGeoId] = trim;
			
			if(pushback) // push_BACK
			{
				elIdSequence.push_back(std::make_pair(linGeoId, qsi1D));
			}
			else // push_FRONT
			{
				elIdSequence.push_front(std::make_pair(linGeoId, qsi1D));
			}
		}
	}
	x = ExactIntersectionPoint[ExactIntersectionPoint.NElements() - 1];
	TPZGeoElSide gelEdge(gel, edge);
	TPZGeoElSide neighEdge = gelEdge.Neighbour();
	while(neighEdge != gelEdge)
	{
		if(neighEdge.Element()->Dimension() == 2)
		{
			break;
		}
		neighEdge = neighEdge.Neighbour();
	}
	gel = neighEdge.Element();
	
	return gel;
}

bool TPZPlaneFracture::EdgeIntersection(TPZGeoEl * gel, TPZFMatrix &x, TPZFMatrix &dx, TPZVec<int> &edge,
										TPZVec< TPZFMatrix > &ExactIntersect, TPZVec< TPZFMatrix > &ModulatedIntersect, double alphaMin)
{
	int nearNode;
	bool IsNearNode = TPZPlaneFracture::NearestNode(gel, x, nearNode, __smallNum);
	
	edge.Resize(0);
	ExactIntersect.Resize(0);
	ModulatedIntersect.Resize(0);
	
	int ncnodes = gel->NCornerNodes();
	
	TPZVec< TPZFMatrix > node(ncnodes), dnode(ncnodes);
	TPZVec<REAL> nodeCoord(3);
	for(int n = 0; n < ncnodes; n++)
	{
		node[n].Resize(3,1);
		dnode[n].Resize(3,1);
		gel->NodePtr(n)->GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			node[n].PutVal(c,0,nodeCoord[c]);
		}
	}//n
	std::list<double> edgeNorm;
	std::map<double, TPZVec<REAL> > alpha; //<double: alpha , <double: alphaNodemod, double: alphaNodesmooth, double: norm, int: edge that intersect> >
	double alphaX, alphaNodemod, alphaNodesmooth, norm;
	for(int n = 0; n < ncnodes; n++)
	{	
		norm = 0.;
		for(int c = 0; c < 3; c++)//computing the direction from node N to node N+1
		{
			if(n != (ncnodes-1) )
			{
				dnode[n](c,0) = node[n+1](c,0) - node[n](c,0);				
			}
			else
			{
				dnode[n](c,0) = node[0](c,0) - node[n](c,0);	
			}
			norm += dnode[n](c,0)*dnode[n](c,0);
		}
		norm = sqrt(norm);
		for(int c = 0; c < 3; c++)//normalizing computed direction
		{
			dnode[n](c,0) /= norm;
		}
		
		alphaX = ComputeAlphaX(x, dx, node[n], dnode[n]);
		alphaNodemod = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, true, false);
		alphaNodesmooth = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, IsNearNode, true);
		if(alphaX >= alphaMin && alphaNodesmooth >= 0. && alphaNodesmooth <= norm)
		{
			int thisEdge = n+ncnodes;
			TPZVec<REAL> someData(4);
			someData[0] = alphaNodemod; someData[1] = alphaNodesmooth; someData[2] = norm; someData[3] = thisEdge;
			alpha[alphaX] = someData;
		}
	}//n
	
	// a aresta que serah interseccionada serah a que tiver menor alphaX não negativo, i.e.: o primeiro par do mapa alpha!
	std::map<double, TPZVec<REAL> >::iterator it = alpha.begin();
	if(it != alpha.end())
	{		
		alphaX = it->first;
		alphaNodemod = (it->second)[0];
		alphaNodesmooth = (it->second)[1];
		norm = (it->second)[2];
		edge.Resize(1); edge[0] = int((it->second)[3]);	
		
		ModulatedIntersect.Resize(1); ModulatedIntersect[0].Resize(3,1);
		for(int c = 0; c < 3; c++)
		{
			ModulatedIntersect[0](c,0) = node[edge[0]-ncnodes](c,0) + alphaNodemod*dnode[edge[0]-ncnodes](c,0);
		}
		ExactIntersect.Resize(1); ExactIntersect[0].Resize(3,1);
		for(int c = 0; c < 3; c++)
		{
			ExactIntersect[0](c,0) = node[edge[0]-ncnodes](c,0) + alphaNodesmooth*dnode[edge[0]-ncnodes](c,0);
		}
		
		if(alphaX <= alphaMin && alpha.size() > 1)
		{
			it++;
			alphaX = it->first;
			alphaNodemod = (it->second)[0];
			alphaNodesmooth = (it->second)[1];
			norm = (it->second)[2];
			edge.Resize(2); edge[1] = int((it->second)[3]);	
			
			ModulatedIntersect.Resize(2); ModulatedIntersect[1].Resize(3,1);
			for(int c = 0; c < 3; c++)
			{
				ModulatedIntersect[1](c,0) = node[edge[1]-ncnodes](c,0) + alphaNodemod*dnode[edge[1]-ncnodes](c,0);
			}
			ExactIntersect.Resize(2); ExactIntersect[1].Resize(3,1);
			for(int c = 0; c < 3; c++)
			{
				ExactIntersect[1](c,0) = node[edge[1]-ncnodes](c,0) + alphaNodesmooth*dnode[edge[1]-ncnodes](c,0);
			}
		}
		return true;
	}
	
	//se o ponto p estah sobre um noh e a direcao dp aponta para fora do elemento...
	if(IsNearNode)
	{
		std::cout << "Estah no noh " << nearNode << " e nao foi encontrada interseccao!" << std::endl;
		std::cout << "Tratar este caso!" << std::endl << std::endl;
	}
	
	return false;
}

TPZGeoEl * TPZPlaneFracture::PointElement(int p, TPZGeoMesh * fractMesh, TPZVec<REAL> &poligonalChain)
{
	TPZVec<REAL> x(3), qsi2D(2);
	for(int c = 0; c < 3; c++)
	{
		x[c] = poligonalChain[3*p+c];
	}
	TPZGeoEl * gel = NULL;
	int nelem = fractMesh->NElements();
	for(int el = 0; el < nelem; el++)//hunting the first element (that it contains the first point)
	{
		TPZGeoEl * firstGel = fractMesh->ElementVec()[el];
		if(firstGel->Dimension() != 2)
		{
			continue;
		}
		if(firstGel->ComputeXInverse(x, qsi2D))
		{
			gel = firstGel;
			break;
		}
	}
	
	return gel;
}

double TPZPlaneFracture::ComputeAlphaNode(TPZFMatrix &x, TPZFMatrix &dx, TPZFMatrix &node, TPZFMatrix &dnode, double norm, bool modulate, bool smooth)
{
	TPZFMatrix xR2(2,1), dxR2(2,1), nodeR2(2,1), dnodeR2(2,1);
	
	//basis change from R3 to R2
	fFromR3toR2.Multiply(x, xR2);
	fFromR3toR2.Multiply(dx, dxR2);
	fFromR3toR2.Multiply(node, nodeR2);
	fFromR3toR2.Multiply(dnode, dnodeR2);
	
	double fractionNumQ =	dxR2(1,0)*nodeR2(0,0) - dxR2(0,0)*nodeR2(1,0) -
	dxR2(1,0)*xR2(0,0) + 	dxR2(0,0)*xR2(1,0);
	
	double fractionDenomQ =	dnodeR2(1,0)*dxR2(0,0) - dnodeR2(0,0)*dxR2(1,0);
	
	double alphaNode = -1.;
	if(fabs(fractionDenomQ) > __smallNum)
	{
		alphaNode = fractionNumQ/fractionDenomQ;
		if(fabs(alphaNode) < __smallNum) alphaNode = 0.;
		if(fabs(alphaNode - norm) < __smallNum) alphaNode = norm;
	}
	
	if(modulate)
	{
		int TrimQTD = fTrimQTD;
		if(smooth)
		{
			TrimQTD *= __TrimQTDmultiplier;
		}
		int nsegm = int(alphaNode/(norm/TrimQTD) + 0.5);
		if(nsegm == 0)//Tirando a interseccao do noh inicial da aresta
		{
			nsegm = 1;
		}
		else if(nsegm == TrimQTD)//Tirando a interseccao do noh final da aresta
		{
			nsegm = TrimQTD - 1;
		}
		alphaNode = nsegm*(norm/TrimQTD);//modulando o ponto para multiplo de (norm/fTrimQTD)
	}
	
	return alphaNode;
}

double TPZPlaneFracture::ComputeAlphaX(TPZFMatrix &x, TPZFMatrix &dx, TPZFMatrix &node, TPZFMatrix &dnode)
{
	TPZFMatrix xR2(2,1), dxR2(2,1), nodeR2(2,1), dnodeR2(2,1);
	
	//basis change from R3 to R2
	fFromR3toR2.Multiply(x, xR2);
	fFromR3toR2.Multiply(dx, dxR2);
	fFromR3toR2.Multiply(node, nodeR2);
	fFromR3toR2.Multiply(dnode, dnodeR2);
	
	//computing alpha (dx vector multiplier to intersect element edges)
	double fractionNumP   =	dnodeR2(1,0)*nodeR2(0,0) - dnodeR2(0,0)*nodeR2(1,0) -
	dnodeR2(1,0)*xR2(0,0) + dnodeR2(0,0)*xR2(1,0);
	
	double fractionDenomP =	dnodeR2(1,0)*dxR2(0,0) - dnodeR2(0,0)*dxR2(1,0);
	
	double alphaX = -1.;
	if(fabs(fractionDenomP) > __smallNum)
	{
		alphaX = fractionNumP/fractionDenomP;
		if(fabs(alphaX) < __smallNum)
		{
			alphaX = 0.;
		}
	}
	
	return alphaX;
}

bool TPZPlaneFracture::NearestNode(TPZGeoEl * gel, TPZFMatrix &x, int &node, double tol)
{
	node = -1;
	bool IsNear = false;
	
	TPZVec<REAL> nodeCoord(3);
	int nnodes = gel->NNodes();
	
	for(int n = 0; n < nnodes; n++)
	{
		double dist = 0.;
		gel->NodePtr(n)->GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			dist += (x(c,0) - nodeCoord[c])*(x(c,0) - nodeCoord[c]);
		}
		dist = sqrt(dist);
		
		if(dist <= tol)
		{
			node = n;
			IsNear = true;
			break;
		}
	}
	
	return IsNear;
}

int TPZPlaneFracture::NearestNode(TPZGeoMesh * gmesh, TPZVec<REAL> &x, double tol)
{
	int node = -1;
	
	TPZVec<REAL> nodeCoord(3);
	int nnodes = gmesh->NNodes();
	
	for(int n = 0; n < nnodes; n++)
	{
		double dist = 0.;
		gmesh->NodeVec()[n].GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			dist += (x[c] - nodeCoord[c])*(x[c] - nodeCoord[c]);
		}
		dist = sqrt(dist);
		
		if(dist <= tol)
		{
			node = n;
			break;
		}
	}
	
	if(node == -1)
	{
		std::cout << "Node not found for coordinates ( " << x[0] << " , " << x[1] << " , " << x[2] << " )" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		
		DebugStop();
	}
	
	return node;
}

double TPZPlaneFracture::LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1)
{
	double dL = 0., L = 0.;
	for(int c = 0; c < 3; c++)
	{
		L += (n1[c]-n0[c])*(n1[c]-n0[c]);
		dL += (x[c]-n0[c])*(x[c]-n0[c]);
	}
	L = sqrt(L);
	dL = sqrt(dL);
	
	#ifdef DEBUG
	if(fabs(L) < __smallNum || fabs(dL) < __smallNum || fabs(L - dL) < __smallNum)
	{
		std::cout << "n0 and n1 are coincident nodes!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	#endif
	
	double qsi = -1. + dL/L*2.;
	
	return qsi;
}

TPZAutoPointer<TPZRefPattern> TPZPlaneFracture::Generate1DRefPatt(std::set<double> &TrimCoord)
{
	if(TrimCoord.size() == 0)
	{
		return NULL;
	}
	
	int Qnodes = TrimCoord.size() + 2;
	int Qelements = TrimCoord.size() + 2;
	int QsubElements = Qelements - 1;
	TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
	
	//1. setting initial and final nodes coordinates of 1D element mesh
	NodeCoord[0][0] = -1.;
	NodeCoord[1][0] =  1.;
	//.
	
	//2. setting intermediate nodes coordinates of 1D element mesh
	int c = 2;
	std::set<double>::iterator it;
	for(it = TrimCoord.begin(); it != TrimCoord.end(); it++)
	{
		double coord = *it;
		
		(NodeCoord[c])[0] = coord;
		c++;
	}
	//.
	
	//3. initializing internal mesh of refPattern
	TPZGeoMesh internalMesh;
	internalMesh.SetMaxNodeId(Qnodes-1);
	internalMesh.SetMaxElementId(Qelements-1);
	internalMesh.NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(NodeCoord[n]);
		internalMesh.NodeVec()[n] = Node[n]; 
	}
	//.
	
	//4. inserting 1D elements on internal mesh of refPattern
	int elId = 0;
	TPZVec <int> Topol(2);
	
	//4.1 inserting father element
	Topol[0] = 0; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__fractureRefinedEl_Mat,internalMesh);
	elId++;
	//.
	
	//4.2 inserting subelements
	//first subelement
	Topol[0] = 0; Topol[1] = 2;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son1 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__fractureRefinedEl_Mat,internalMesh);
	son1->SetFather(father);
	son1->SetFather(father->Index());
	elId++;
	//
	
	//last subelement
	Topol[0] = TrimCoord.size() + 1; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son2 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__fractureRefinedEl_Mat,internalMesh);
	son2->SetFather(father);
	son2->SetFather(father->Index());
	elId++;
	//
	
	for(int el = 2; el < QsubElements; el++)
	{
		Topol[0] = el; Topol[1] = el+1;
		TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__fractureRefinedEl_Mat,internalMesh);
		son->SetFather(father);
		son->SetFather(father->Index());
		elId++;
	}
	//.
	//.
	
	internalMesh.BuildConnectivity();
	
	TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(internalMesh);
	TPZAutoPointer<TPZRefPattern> Found = gRefDBase.FindRefPattern(refPattern);
	if(!Found)
	{
		gRefDBase.InsertRefPattern(refPattern);
		refPattern->InsertPermuted();
		
		return refPattern;
	}
	else
	{
		return Found;
	}
}

void TPZPlaneFracture::UpdatePoligonalChain(TPZGeoMesh * gmesh, std::list< std::pair<int,double> > &elIdSequence,
						  TPZVec<REAL> &poligonalChainUpdated)
{
	int nptos = elIdSequence.size();
	poligonalChainUpdated.Resize(3*nptos);

	TPZVec<REAL> qsi1Dvec(1), ptoCoord(3);
	int el1Did, posX, posY, posZ, p = 0;
	double qsi1D;
	
	std::list< std::pair<int,double> >::iterator it;
	for(it = elIdSequence.begin(); it != elIdSequence.end(); it++)
	{
		el1Did = it->first;
		TPZGeoEl * el1D = gmesh->ElementVec()[el1Did];
		
		#ifdef DEBUG
		int elDim = el1D->Dimension();
		if(elDim != 1)
		{
			std::cout << "The elIdSequence supposedly would contains ids of elements exclusively 1D!" << std::endl;
			std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
			DebugStop();
		}
		#endif
		
		qsi1D = it->second;
		qsi1Dvec[0] = qsi1D;
		
		el1D->X(qsi1Dvec, ptoCoord);
		
		posX = 3*p;
		posY = 3*p+1;
		posZ = 3*p+2;
		
		poligonalChainUpdated[posX] = ptoCoord[0];
		poligonalChainUpdated[posY] = ptoCoord[1];
		poligonalChainUpdated[posZ] = ptoCoord[2];
		
		p++;
	}
}

void TPZPlaneFracture::GenerateCrackBoundary(TPZGeoMesh * gmesh, std::list< std::pair<int,double> > &elIdSequence)
{
	TPZVec<REAL> qsi0vec(1), qsi1vec(1), node0coord(3), node1coord(3);
	TPZVec<int> Topol(2);
	int el0id, el1id, n0, n1;
	double qsi0, qsi1;
	std::list< std::pair<int,double> >::iterator crackit0, crackit1, crackitEnd;
	crackitEnd = elIdSequence.end(); crackitEnd--;
	for(crackit0 = elIdSequence.begin(); crackit0 != crackitEnd; crackit0++)
	{
		crackit1 = crackit0; crackit1++;
		
		el0id = crackit0->first;
		TPZGeoEl * el0 = gmesh->ElementVec()[el0id];
		qsi0 = crackit0->second;
		qsi0vec[0] = qsi0;
		el0->X(qsi0vec, node0coord);
		n0 = NearestNode(gmesh, node0coord, __smallNum);
		Topol[0] = n0;
		
		el1id = crackit1->first;
		TPZGeoEl * el1 = gmesh->ElementVec()[el1id];
		qsi1 = crackit1->second;
		qsi1vec[0] = qsi1;
		el1->X(qsi1vec, node1coord);
		n1 = NearestNode(gmesh, node1coord, __smallNum);
		Topol[1] = n1;
		
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, __fractureLine_Mat, *gmesh);
	}
}
