
#include "TPZExtendGridDimension.h"
#include "pzgmesh.h"
#include "pzgeoel.h"

using namespace std;

TPZExtendGridDimension::TPZExtendGridDimension(char *geofile,REAL thickness) : fFineFileMesh(geofile){
	
	fThickness = thickness;
}

TPZExtendGridDimension::TPZExtendGridDimension(TPZGeoMesh *finegeomesh,REAL thickness){
	
	fFineGeoMesh = finegeomesh;
	fThickness = thickness;
}

TPZGeoMesh * TPZExtendGridDimension::ExtendedMesh(){
	/**
	 * a malha 2D sera extendida para uma malha 3D: logo ela eh plana e conforme
	 * as incidencias devem estar dadas em sentido antihorario - vista superior do plano XY
	 * e as coordenadas sao da forma (x,y,0)
	 * a terceira componente devera ser thickness: altura da malha
	 * os elementos 2D podem ser triangulos ou quadrilateros
	 * si os elementos sao triangulos os elementos 3D serao prismas retos
	 * si os elementos sao quadrilateros os elementos 3D serao hexaedros retos
	 */
	TPZGeoMesh *extendedmesh = new TPZGeoMesh;
	int maxid = fFineGeoMesh->CreateUniqueNodeId();
	int nelem = fFineGeoMesh->ElementVec().NElements(),i,j;
	TPZGeoNode gnode;
	int nnodes = fFineGeoMesh->NodeVec().NElements();
	//o numero de nos sera duplicado
	extendedmesh->NodeVec().Resize(2*nnodes);
	TPZVec<REAL> coord(3);
	int index;
	//criacao dos nos da malha 3D
	for(i=0;i<nnodes;i++){
		gnode = fFineGeoMesh->NodeVec()[i];
		coord[0] = gnode.Coord(0);
		coord[1] = gnode.Coord(1);
		coord[2] = gnode.Coord(2);// = 0.0
		extendedmesh->NodeVec()[i].Initialize(coord,*extendedmesh);
		coord[2] = fThickness;
		extendedmesh->NodeVec()[i+maxid].Initialize(coord,*extendedmesh);
	}
	//criacao de elementos da malha 3D
	TPZGeoEl *gel;
	TPZVec<int> incidel;
	for(i=0;i<nelem;i++){
		gel = fFineGeoMesh->ElementVec()[i];
		if(!gel) continue;
		int type = gel->Type();
		if(type==2){//triangle
			incidel.Resize(6);
			if(fThickness > 0){
				for(j=0;j<3;j++) incidel[j] = gel->NodeIndex(j);
				for(j=3;j<6;j++) incidel[j] = incidel[j-3]+maxid;
			} else if(fThickness < 0){
				for(j=0;j<3;j++) incidel[j] = gel->NodeIndex(j)+maxid;
				for(j=3;j<6;j++) incidel[j] = gel->NodeIndex(j-3);
			}
			int matind = gel->MaterialId();
			extendedmesh->CreateGeoElement(EPrisma,incidel,matind,index);
		}
		if(type==3){//quadrilateral
			incidel.Resize(8);
			if(fThickness > 0){
				for(j=0;j<4;j++) incidel[j] = gel->NodeIndex(j);
				for(j=4;j<8;j++) incidel[j] = incidel[j-4]+maxid;
			} else if(fThickness < 0){
				for(j=0;j<4;j++) incidel[j] = gel->NodeIndex(j)+maxid;
				for(j=4;j<8;j++) incidel[j] = gel->NodeIndex(j-4);
			}
			int matind = gel->MaterialId();
			extendedmesh->CreateGeoElement(ECube,incidel,matind,index);     
		}
	}
	extendedmesh->BuildConnectivity();
	return extendedmesh;
}

void TPZExtendGridDimension::PrintGeneratedMesh(ostream &out){
	
	fFineGeoMesh->Print(out);
}
