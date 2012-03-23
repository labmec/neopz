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

using namespace pzgeom;

TPZPlaneFracture::TPZPlaneFracture(TPZGeoMesh * planeMesh, int nodeOrigin, int nodeX, int nodeY)
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
}

TPZPlaneFracture::~TPZPlaneFracture()
{
	
}

TPZAutoPointer<TPZRefPattern> TPZPlaneFracture::Generate1DRefPatt(set<double> &TrimCoord)
{
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
	set<double>::iterator it;
	int c = 2;
	for(it = TrimCoord.begin(); it != TrimCoord.end(); it++)
	{
		(NodeCoord[c])[0] = *it;
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
		Node[n].SetCoord(&NodeCoord[n][0]);
		internalMesh.NodeVec()[n] = Node[n]; 
	}
	//.
	
	//4. inserting 1D elements on internal mesh of refPattern
	int elId = 0;
	int matElId = -100;
	TPZVec <int> Topol(2);
	
	//4.1 inserting father element
	Topol[0] = 0; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,matElId,internalMesh);
	elId++;
	//.
	
	//4.2 inserting subelements
	if(TrimCoord.size() > 0)
	{
		//first subelement
		Topol[0] = 0; Topol[1] = 2;
		TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son1 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,matElId,internalMesh);
		son1->SetFather(father);
		son1->SetFather(father->Index());
		elId++;
		//
		
		//last subelement
		Topol[0] = TrimCoord.size() + 1; Topol[1] = 1;
		TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son2 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,matElId,internalMesh);
		son2->SetFather(father);
		son2->SetFather(father->Index());
		elId++;
		//
		
		for(int el = 2; el < QsubElements; el++)
		{
			Topol[0] = el; Topol[1] = el+1;
			TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,matElId,internalMesh);
			son->SetFather(father);
			son->SetFather(father->Index());
			elId++;
		}
	}
	//.
	//.
	
	internalMesh.BuildConnectivity();
	
	TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(internalMesh);
	
	if(!gRefDBase.FindRefPattern(refPattern))
	{
		gRefDBase.InsertRefPattern(refPattern);
	}
	refPattern->InsertPermuted();
	
	return refPattern;
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
	
	double qsi = -1. + dL/L*2.;
	
	return qsi;
}

TPZVec<REAL> TPZPlaneFracture::EdgeIntersection(TPZGeoEl * gel, const TPZVec<REAL> &p, const TPZVec<REAL> &dp, int &edge)
{
	#ifdef DEBUG
	if(!gel)
	{
		std::cout << "Invalid element on " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	if(p.NElements() != 3 || dp.NElements() != 3)
	{
		std::cout << "Invalid point p and/or direction dp on " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	#endif

	double norm = 0.;
	for(int c = 0; c < 3; c++)
	{
		norm += dp[c]*dp[c];
	}
	norm = sqrt(norm);

	#ifdef DEBUG
	if(fabs(norm) < 1.E-18)
	{
		std::cout << "Null direction dp on " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	#endif

	TPZFMatrix x(3,1), xR2(2,1), dx(3,1), dxR2(2,1);
	for(int c = 0; c < 3; c++)
	{
		x.PutVal(c, 0, p[c]);
		dx.PutVal(c, 0, dp[c]/norm);
	}
	fFromR3toR2.Multiply(x, xR2);
	fFromR3toR2.Multiply(dx, dxR2);
	
	int ncnodes = gel->NCornerNodes();

	TPZVec< TPZFMatrix > node(ncnodes), dnode(ncnodes), nodeR2(ncnodes), dnodeR2(ncnodes);
	TPZVec<REAL> nodeCoord(3);
	for(int n = 0; n < ncnodes; n++)
	{
		node[n].Resize(3,1);	nodeR2[n].Resize(2,1);
		dnode[n].Resize(3,1);	dnodeR2[n].Resize(2,1);
		gel->Mesh()->NodeVec()[gel->NodeIndex(n)].GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			node[n].PutVal(c,0,nodeCoord[c]);
		}
	}
	std::map<double,int> alphaP; //<double: alpha , int: edge that intersect>
	double fractionNum, fractionDenom, alpha;
	for(int n = 0; n < ncnodes; n++)
	{	
		norm = 0.;
		for(int c = 0; c < 3; c++)//computing the directions from node N to node N+1
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
		for(int c = 0; c < 3; c++)//normalizing computed directions
		{
			dnode[n](c,0) /= norm;
		}
		//basis change from R3 to R2
		fFromR3toR2.Multiply(node[n], nodeR2[n]);
		fFromR3toR2.Multiply(dnode[n], dnodeR2[n]);
		
		//computing alphaP (dx vector multiplier to intersect element edges)
		fractionNum =	dnodeR2[n](1,0)*nodeR2[n](0,0) -
						dnodeR2[n](0,0)*nodeR2[n](1,0) -
						dnodeR2[n](1,0)*xR2(0,0) +
						dnodeR2[n](0,0)*xR2(1,0);
		
		fractionDenom =	dnodeR2[n](1,0)*dxR2(0,0) -
						dnodeR2[n](0,0)*dxR2(1,0);
		alpha = -1.;
		if(fabs(fractionDenom) > 1.E-15)
		{
			// alphaP eh a solucao do sistema: {p + alphaP.dp == q + alphaQ.dq}, ou seja,
			// a norma que multiplica o vetor dp e cruza a reta (q+alphaQ.dq)
			alpha = fractionNum/fractionDenom;
		}
		if(alpha > 0.)//entra somente positivos
		{
			alphaP[alpha] = (n+ncnodes);
		}
	}
	
	#ifdef DEBUG
	if(alphaP.size() == 0)
	{
		// significa que o ponto e a direcao dada intersecciona as arestas do elemento no sentido negativo da direcao,
		// portanto p nao estah dentro do dominio do elemento!
		std::cout << "Given point p does NOT belongs to element domain!!!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	#endif
	
	// a aresta que serah interseccionada serah a que tiver menor alpha não negativo, i.e.: o primeiro par do mapa alphaP!
	std::map<double,int>::iterator it = alphaP.begin();
	alpha = it->first;
	edge = it->second;
	
	TPZVec<REAL> intersect(3);
	for(int c = 0; c < 3; c++)
	{
		intersect[c] = x(c,0) + alpha*dx(c,0);
	}
	return intersect;
}

TPZGeoMesh * TPZPlaneFracture::GetFractureMesh(TPZVec< TPZVec<REAL> > &crackTip)
{
	#ifdef DEBUG
	for(int bound = 0; bound < crackTip.NElements(); bound++)
	{
		int ncoord = crackTip[bound].NElements();
		if(ncoord%3 != 0)
		{
			std::cout << "crackTip boundary #" << bound << " dont have groups of 3 coordinates (x,y,z)!" << std::endl;
			std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
			DebugStop();
		}
	}
	#endif
	
	TPZGeoMesh * fractMesh;
	*fractMesh = *fplaneMesh;
	int nelem = fractMesh->NElements();
	
	TPZVec<REAL> x(3), dx(3), qsi(2);
	int edge;
	for(int actualBoundary = 0; actualBoundary < crackTip.NElements(); actualBoundary++)
	{
		TPZVec<REAL> boundaryPoint = crackTip[actualBoundary];
		int npoints = (boundaryPoint.NElements())/3;
		for(int p = 0; p < (npoints-1); p++)
		{
			int thispoint = 3*p;
			int nextpoint = 3*(p+1);
			for(int c = 0; c < 3; c++)
			{
				x[c] = boundaryPoint[thispoint+c];
				dx[c] = boundaryPoint[nextpoint+c] - boundaryPoint[thispoint+c];
			}
			TPZGeoEl * gel = NULL;
			if(p == 0)
			{
				for(int el = 0; el < nelem; el++)
				{
					TPZGeoEl * firstGel = fractMesh->ElementVec()[el];
					if(firstGel->Dimension() != 2)
					{
						continue;
					}
					firstGel->ComputeXInverse(x, qsi);
					if(firstGel->IsInParametricDomain(qsi))
					{
						gel = firstGel;
						break;
					}
				}
				if(!gel)
				{
					std::cout << "first point of crack tip boundary does NOT belong to any element!" << std::endl;
					std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
					DebugStop();
				}
			}
			else
			{
					//TODO: pegar vizinho pelo lado "edge" e verificar se computexinverse estah nele!
					//se nao, continua a projetar na mesma direcao ateh achar o elemento que contem o proximo pto
					//se sim, muda a direcao conforme prox ponto etc...
			}

			TPZVec<REAL> intersectionPoint = EdgeIntersection(gel, x, dx, edge);
			TPZVec<REAL> n0(3), n1(3);
			fractMesh->NodeVec()[gel->SideNodeIndex(edge, 0)].GetCoordinates(n0);
			fractMesh->NodeVec()[gel->SideNodeIndex(edge, 1)].GetCoordinates(n1);
			double qsi = LinearComputeXInverse(intersectionPoint, n0, n1);//////// modular o qsi para minimizar opcoes de refpattern!!!
			std::set<double> trim;
			trim.insert(qsi);
			TPZAutoPointer<TPZRefPattern> linRefp = Generate1DRefPatt(trim);
			
			//Adding 1D elements on gMesh
			int lineId, lineMat = 9876;
			TPZVec <int> Topol(2);
			Topol[0] = gel->SideNodeIndex(edge, 0); Topol[1] = gel->SideNodeIndex(edge, 1);
			TPZGeoEl * linGel = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (lineId,Topol, lineMat,*fractMesh);
			fractMesh->BuildConnectivity();
			TPZVec<TPZGeoEl*> sons;
			linGel->SetRefPattern(linRefp);
			linGel->Divide(sons);	
		}
	}
	
	return fractMesh;
}
