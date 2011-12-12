/**
 * @file
 * @brief Contains the implementation of the TPZExtendGridDimension methods. 
 */

#include "TPZExtendGridDimension.h"
#include "pzgmesh.h"
#include "pzgeoel.h"

using namespace std;

TPZExtendGridDimension::TPZExtendGridDimension(char *geofile,REAL thickness) : fFineFileMesh(geofile){
	
	fThickness = thickness;
}

TPZExtendGridDimension::TPZExtendGridDimension(TPZAutoPointer<TPZGeoMesh> &finegeomesh,REAL thickness){
	
	fFineGeoMesh = finegeomesh;
	fThickness = thickness;
}

TPZGeoMesh *TPZExtendGridDimension::ExtendedMesh() {
	// a malha 2D sera extendida para uma malha 3D: logo ela eh plana e conforme
	// as incidencias devem estar dadas em sentido antihorario - vista superior do plano XY
	// e as coordenadas sao da forma (x,y,0)
	// a terceira componente devera ser thickness: altura da malha
	// os elementos 2D podem ser triangulos ou quadrilateros
	// si os elementos sao triangulos os elementos 3D serao prismas retos
	// si os elementos sao quadrilateros os elementos 3D serao hexaedros retos
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
		if(type==2) {             //triangle
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
		if(type==3) {             //quadrilateral
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

TPZGeoMesh *TPZExtendGridDimension::ExtendedMesh(int naumentedlayers,int matidbottom,int matidtop){
	if(naumentedlayers < 1)   // returns the same geometric mesh
		return fFineGeoMesh.operator->();
	
	// a malha 2D sera extendida para uma malha 3D: logo ela eh plana e conforme
	// as incidencias devem estar dadas em sentido antihorario - vista superior do plano XY
	// e as coordenadas sao da forma (x,y,0)
	// a terceira componente devera ser thickness: altura da malha
	// os elementos 2D podem ser triangulos ou quadrilateros
	// si os elementos sao triangulos os elementos 3D serao prismas retos
	// si os elementos sao quadrilateros os elementos 3D serao hexaedros retos
	
	TPZGeoMesh *extendedmesh = new TPZGeoMesh;
	int maxid = fFineGeoMesh->CreateUniqueNodeId();
	int nelem = fFineGeoMesh->ElementVec().NElements();
	int i,j,k;
	TPZGeoNode gnode;
	int nnodes = fFineGeoMesh->NodeVec().NElements();
	
	//o numero de nos sera duplicado
	extendedmesh->NodeVec().Resize((naumentedlayers+1)*nnodes);
	TPZVec<REAL> coord(3);
	int index;
	
	//criacao dos nos da malha 3D
	for(i=0;i<nnodes;i++) {
		gnode = fFineGeoMesh->NodeVec()[i];
		coord[0] = gnode.Coord(0);
		coord[1] = gnode.Coord(1);
		coord[2] = gnode.Coord(2);// = 0.0
		extendedmesh->NodeVec()[i].Initialize(coord,*extendedmesh);
		for(j=0;j<naumentedlayers;j++) {
			coord[2] += fThickness;
			extendedmesh->NodeVec()[i+(j+1)*maxid].Initialize(coord,*extendedmesh);
		}
	}
	//criacao de elementos da malha 3D
	TPZGeoEl *gel;
	TPZVec<int> incidelorig;
	TPZVec<int> incidel;
	int matind;
	for(i=0;i<nelem;i++) {
		gel = fFineGeoMesh->ElementVec()[i];
		if(!gel) continue;
		int type = gel->Type();
		if(type==2 || type==3) {//triangle
			nnodes = gel->NNodes();
			incidelorig.Resize(nnodes);
			for(j=0;j<nnodes;j++)
				incidelorig[j] = gel->NodeIndex(j);
			incidel.Resize(2*nnodes);
			for(k=0;k<naumentedlayers;k++) {
				if(fThickness > 0) {
					for(j=0;j<nnodes;j++) incidel[j] = incidelorig[j];
					for(j=nnodes;j<2*nnodes;j++) incidel[j] = incidel[j-nnodes]+maxid;
					// initial indexes of the nodes must to be update, upper triangle
					for(j=0;j<nnodes;j++) incidelorig[j] = incidel[j+nnodes];
				} else if(fThickness < 0) {
					for(j=0;j<nnodes;j++) incidel[j] = incidelorig[j]+maxid;
					for(j=nnodes;j<2*nnodes;j++) incidel[j] = incidelorig[j-nnodes];
					// initial indexes of the nodes must to be update, lower triangle
					for(j=0;j<nnodes;j++) incidelorig[j] = incidel[j];
				}
				matind = gel->MaterialId();
				if(type==2) gel = extendedmesh->CreateGeoElement(EPrisma,incidel,matind,index);
				else gel = extendedmesh->CreateGeoElement(ECube,incidel,matind,index);
			}
		}
		// When the geometric element has boundary condition
		else {
			if(gel->MaterialId() < 0) {
				nnodes = gel->NNodes();
				incidelorig.Resize(nnodes);
				for(j=0;j<nnodes;j++) incidelorig[j] = gel->NodeIndex(j);
				incidel.Resize(2*nnodes);
				for(k=0;k<naumentedlayers;k++) {
					if(fThickness > 0) {
						for(j=0;j<nnodes;j++) incidel[j] = incidelorig[j];
						for(j=nnodes;j<2*nnodes;j++) incidel[j] = incidel[j-nnodes]+maxid;
						// initial indexes of the nodes must to be update, upper triangle
						for(j=0;j<nnodes;j++) incidelorig[j] = incidel[j+nnodes];
					} else if(fThickness < 0) {
						for(j=0;j<nnodes;j++) incidel[j] = incidelorig[j]+maxid;
						for(j=nnodes;j<2*nnodes;j++) incidel[j] = incidelorig[j-nnodes];
						// initial indexes of the nodes must to be update, lower triangle
						for(j=0;j<nnodes;j++) incidelorig[j] = incidel[j];
					}
					matind = gel->MaterialId();
					if(nnodes==1) gel = extendedmesh->CreateGeoElement(EOned,incidel,matind,index);
					else gel = extendedmesh->CreateGeoElement(EQuadrilateral,incidel,matind,index);
				}
			}
		}
	}
	// Inserting all the elements of the original two-dimensional mesh as bc elements (bottom bc)
	if(matidbottom < 0) {
		for(i=0;i<nelem;i++) {
			gel = fFineGeoMesh->ElementVec()[i];
			if(!gel || gel->MaterialId() < 0) continue;
			nnodes = gel->NNodes();
			incidelorig.Resize(nnodes);
			for(j=0;j<nnodes;j++)
				incidelorig[j] = gel->NodeIndex(j);
			if(nnodes==3) gel = extendedmesh->CreateGeoElement(ETriangle,incidelorig,matidbottom,index);
			else gel = extendedmesh->CreateGeoElement(EQuadrilateral,incidelorig,matidbottom,index);
		}
	}
	// Inserting all the elements on last layer inserted as bc elements (top bc)
	if(matidtop < 0) {
		for(i=0;i<nelem;i++) {
			gel = fFineGeoMesh->ElementVec()[i];
			if(!gel || gel->MaterialId() < 0) continue;
			nnodes = gel->NNodes();
			incidelorig.Resize(nnodes);
			for(j=0;j<nnodes;j++)
				incidelorig[j] = gel->NodeIndex(j)+(naumentedlayers*maxid);
			if(nnodes==3) gel = extendedmesh->CreateGeoElement(ETriangle,incidelorig,matidtop,index);
			else gel = extendedmesh->CreateGeoElement(EQuadrilateral,incidelorig,matidtop,index);
		}
	}
	extendedmesh->BuildConnectivity();
	return extendedmesh;
}

void TPZExtendGridDimension::PrintGeneratedMesh(ostream &out){
	
	fFineGeoMesh->Print(out);
}
