#include <iostream>


#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "tpzquadrilateral.h"
#include "TPZQuadSphere.h"
#include "tpzgeoblend.h"
#include "TPZWavyLine.h"

#include "TPZVTKGeoMesh.h"

void QuadWavy();
void QuadSphere();
void ExemploIni();

int main(int argc, char *argv[])
{
	QuadSphere();
	return 0;
}

void QuadWavy()
{
	///INSTANCIACAO DA MALHA GEOMETRICA
	TPZGeoMesh * geomesh = new TPZGeoMesh;
	
	int nnodes = 4;//quantidade de nos da malha geometrica
	geomesh->NodeVec().Resize(nnodes);
	
	///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
	TPZGeoNode node0, node1, node2, node3;
	
	TPZVec<REAL> coord(3,0.);
	
	//no 0
	coord[0] = 0.;
	coord[1] = 0.;
	coord[2] = 0.;
	node0.SetNodeId(0);
	node0.SetCoord(coord);
	geomesh->NodeVec()[0] = node0;
	
	//no 1
	coord[0] = 1.;
	coord[1] = 0.;
	coord[2] = 0.;
	node1.SetNodeId(1);
	node1.SetCoord(coord);
	geomesh->NodeVec()[1] = node1;
	
	//no 2
	coord[0] = 1.;
	coord[1] = 2.;
	coord[2] = 0.;
	node2.SetNodeId(2);
	node2.SetCoord(coord);
	geomesh->NodeVec()[2] = node2;
	
	//no 3
	coord[0] = 0.;
	coord[1] = 1.;
	coord[2] = 0.;
	node3.SetNodeId(3);
	node3.SetCoord(coord);
	geomesh->NodeVec()[3] = node3;
		
	//INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
	TPZVec<int64_t> topology(4);//Quadrilatero ira utilizar 4 nos
	topology[0] = 0;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 2;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	int materialId = 1;
	
	TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > *quadrilatero =
	new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (topology,materialId,*geomesh);
	
	
	//INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
	TPZVec<int64_t> topologyWavy(2);//Quadrilatero ira utilizar 4 nos
	topologyWavy[0] = 1;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topologyWavy[1] = 2;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	
	int materialIdW = -1;
	
	TPZGeoElRefPattern< pzgeom::TPZWavyLine > * wavy =
	new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy,materialIdW,*geomesh);
	
	TPZVec<REAL> wavedir(3,0.);
	wavedir[0] = 0.;
	wavedir[2] = 0.1;
	int numwaves = 5;
	wavy->Geom().SetData(wavedir,numwaves);
		
	//CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
	geomesh->BuildConnectivity();
	
	///EXEMPLOS DE MAPEAMENTO GEOMETRICO ATRAVES DO METODO X
	TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
	
	std::cout << "\n\n";
	
	///Refinando o elemento quadrilateral em subelementos (chamados de "filhos")
	TPZVec<TPZGeoEl *> sons;
	
	const int nref = 7;
	for (int iref = 0; iref < nref; iref++) {
		int nel = geomesh->NElements();
		for (int iel = 0; iel < nel; iel++) {
			TPZGeoEl *gel = geomesh->ElementVec()[iel];
			if(!gel->HasSubElement()) gel->Divide(sons);
		}
	}

	std::ofstream outfile("malha.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile,true);
	
}

void QuadSphere()
{
	///INSTANCIACAO DA MALHA GEOMETRICA
	TPZGeoMesh * geomesh = new TPZGeoMesh;
	
	int nnodes = 6;//quantidade de nos da malha geometrica
	geomesh->NodeVec().Resize(nnodes);
	
	///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
	TPZGeoNode node0, node1, node2, node3, node4, node5;	
	TPZVec<REAL> coord(3,0.);
	const REAL r = 2.;
	TPZManVector<REAL,3> xc(3,0.);
	xc[0] = 1.;
	xc[1] = 2.;
	xc[2] = 3.;
	
	//no 0
	coord[0] = xc[0] + r;
	coord[1] = xc[1] + 0.;
	coord[2] = xc[2] + 0.;
	node0.SetNodeId(0);
	node0.SetCoord(coord);
	geomesh->NodeVec()[0] = node0;
	
	//no 1
	coord[0] = xc[0] + 0.;
	coord[1] = xc[1] + 0.;
	coord[2] = xc[2] + -r;
	node1.SetNodeId(1);
	node1.SetCoord(coord);
	geomesh->NodeVec()[1] = node1;
	
	//no 2
	coord[0] = xc[0] + 0.;
	coord[1] = xc[1] + r;
	coord[2] = xc[2] + 0.;
	node2.SetNodeId(2);
	node2.SetCoord(coord);
	geomesh->NodeVec()[2] = node2;
	
	//no 3
	coord[0] = xc[0] + 0.;
	coord[1] = xc[1] + 0.;
	coord[2] = xc[2] + r;
	node3.SetNodeId(3);
	node3.SetCoord(coord);
	geomesh->NodeVec()[3] = node3;
	
	//no 4
	coord[0] = xc[0] + 0.;
	coord[1] = xc[1] + -r;
	coord[2] = xc[2] + 0.;
	node4.SetNodeId(4);
	node4.SetCoord(coord);
	geomesh->NodeVec()[4] = node4;
	
	//no 5
	coord[0] = xc[0] + -r;
	coord[1] = xc[1] + 0.;
	coord[2] = xc[2] + 0.;
	node5.SetNodeId(5);
	node5.SetCoord(coord);
	geomesh->NodeVec()[5] = node5;
	
	int materialId = 1;
	
	// El 0
	//INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
	TPZVec<int64_t> topology(4);//Quadrilatero ira utilizar 4 nos
	topology[0] = 0;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 2;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > * quadrilatero =
	new TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > (topology,materialId,*geomesh);
	
	quadrilatero->Geom().SetData(r,xc);
	
	// El 1
	topology[0] = 2;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 5;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > * quadrilatero2 =
	new TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > (topology,materialId,*geomesh);
	
	quadrilatero2->Geom().SetData(r,xc);
	
	// El 2	
	topology[0] = 5;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 4;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > * quadrilatero3 =
	new TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > (topology,materialId,*geomesh);
	
	quadrilatero3->Geom().SetData(r,xc);
	
	// El 3
	topology[0] = 4;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 0;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > * quadrilatero4 =
	new TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > (topology,materialId,*geomesh);
	
	quadrilatero4->Geom().SetData(r,xc);
	
	//CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
	geomesh->BuildConnectivity();
	
	
	///EXEMPLOS DE MAPEAMENTO GEOMETRICO ATRAVES DO METODO X
	TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
	TPZFNMatrix<9,REAL> jac(2,2),axes(2,3),jacinv(2,2);
	REAL detjac;
	
	std::cout << "\n\n";
	
	qsi[0] = +0.23;
	qsi[1] = -0.17;
	quadrilatero->X(qsi,xqsi);
	std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
	quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
	jac.Print("Jac:");
	axes.Print("axes");
	std::cout << "detjac = " << detjac << std::endl;

	qsi[0] = -0.89;
	qsi[1] = +0.31;
	quadrilatero->X(qsi,xqsi);
	std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
	quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
	jac.Print("Jac:");
	axes.Print("axes");
	std::cout << "detjac = " << detjac << std::endl;		
	
	qsi[0] = +0.05;
	qsi[1] = +1.0;
	quadrilatero->X(qsi,xqsi);
	std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
	quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
	jac.Print("Jac:");
	axes.Print("axes");
	std::cout << "detjac = " << detjac << std::endl;
	
	qsi[0] = -0.73;
	qsi[1] = -0.35;
	quadrilatero->X(qsi,xqsi);
	std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
	quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
	jac.Print("Jac:");
	axes.Print("axes");
	std::cout << "detjac = " << detjac << std::endl;
	
	std::cout << "\n\n";
	
	
	TPZVec<TPZGeoEl *> sons;
	const int nref = 6;
	for (int iref = 0; iref < nref; iref++) {
		int nel = geomesh->NElements();
		for (int iel = 0; iel < nel; iel++) {
			TPZGeoEl *gel = geomesh->ElementVec()[iel];
			if(!gel->HasSubElement()) gel->Divide(sons);
		}
	}
	
	std::ofstream outfile("malhaesfera.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile);

	
}

/** Geracao de uma malha geometrica de 1 elemento quadrilateral apenas */
void ExemploIni()
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
	TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nnodes = 6;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node0, node1, node2, node3, node4, node5;
    
    TPZVec<REAL> coord(3,0.);
    
    //no 0
    coord[0] = 3.;
    coord[1] = 2.;
    coord[2] = 3.;
    node0.SetNodeId(0);
    node0.SetCoord(coord);
    geomesh->NodeVec()[0] = node0;
    
    //no 1
    coord[0] = 7.;
    coord[1] = 1.;
    coord[2] = 4.;
    node1.SetNodeId(1);
    node1.SetCoord(coord);
    geomesh->NodeVec()[1] = node1;
    
    //no 2
    coord[0] = 8.;
    coord[1] = 6.;
    coord[2] = 2.;
    node2.SetNodeId(2);
    node2.SetCoord(coord);
    geomesh->NodeVec()[2] = node2;
    
    //no 3
    coord[0] = 4.;
    coord[1] = 4.;
    coord[2] = 6.;
    node3.SetNodeId(3);
    node3.SetCoord(coord);
    geomesh->NodeVec()[3] = node3;
    
    //no 4
    coord[0] = 9.;
    coord[1] = 3.;
    coord[2] = 6.;
    node4.SetNodeId(4);
    node4.SetCoord(coord);
    geomesh->NodeVec()[4] = node4;
    
    //no 5
    coord[0] = 8.5;
    coord[1] = 7.;
    coord[2] = 6.;
    node5.SetNodeId(5);
    node5.SetCoord(coord);
    geomesh->NodeVec()[5] = node5;
    
    //INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
    TPZVec<int64_t> topology(4);//Quadrilatero ira utilizar 4 nos
    topology[0] = 0;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 2;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
    
    int materialId = 1;
    
    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * quadrilatero =
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (topology,materialId,*geomesh);
    
    
    topology[0] = 1;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 4;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 5;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    topology[3] = 2;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (topology,materialId,*geomesh);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    ///EXEMPLOS DE MAPEAMENTO GEOMETRICO ATRAVES DO METODO X
    TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
    
    std::cout << "\n\n";
    
    qsi[0] = +0.23;
    qsi[1] = -0.17;
    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";

    TPZFNMatrix<9,REAL> jac(2,2),axes(2,3),jacinv(2,2);
    REAL detjac;
    quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
    qsi[0] = -0.89;
    qsi[1] = +0.31;
    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    
    
    qsi[0] = +0.05;
    qsi[1] = +1.0;
    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    
    qsi[0] = -0.73;
    qsi[1] = -0.35;
    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    
    std::cout << "\n\n";
    
    
    ///Refinando o elemento quadrilateral em subelementos (chamados de "filhos")
    TPZVec<TPZGeoEl *> sons0, sons1;
    quadrilatero->Divide(sons0);
    
    ///Refinando os filhos do elemento quadrilateral
    for(int s = 0; s < sons0.NElements(); s++)
    {
        sons1.Resize(0);
        TPZGeoEl * actSon = sons0[s];
        actSon->Divide(sons1);
    }
    
    std::ofstream outfile("malha.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile);
    
    geomesh->Print();
    
}
