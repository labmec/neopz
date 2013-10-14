/**
 * @file
 * @brief Contains the implementation of the TPZCheckGeom methods.
 */

#include <fstream>

#include "pzcheckgeom.h"
#include "pzquad.h"
#include "pztrnsform.h"
#include "pzstack.h"
#include "pzgeoelside.h"

using namespace std;

TPZCheckGeom::TPZCheckGeom() {
	fMesh = 0;
}

int TPZCheckGeom::CheckElement(TPZGeoEl *gel) {
	
	int check = 0;
	int nsides = gel->NSides();
	int geldim = gel->Dimension();
	int is;
	for(is=nsides-1; is>=0; is--) {
		TPZStack<TPZGeoElSide> highdim;
		int dim;
		int sidedim = gel->SideDimension(is);
		for(dim = sidedim+1; dim<= geldim; dim++) {
			gel->AllHigherDimensionSides(is,dim,highdim);
		}
		int nhighdim = highdim.NElements();
		int idim;
		for(idim=0; idim<nhighdim; idim++) {
			check = (CheckSideTransform(gel,is,highdim[idim].Side()) || check);
		}
	}
	return check;
}

int TPZCheckGeom::PerformCheck() {
	long nel = fMesh->NElements();
	long iel;
	int check = 0;
	for(iel = 0; iel<nel; iel++) {
		TPZGeoEl *gel = fMesh->ElementVec()[iel];
		if(!gel) continue;
		TPZStack<TPZGeoEl *> subel;
		gel->Divide(subel);		
		if(iel==3){//TESTE
			cout << "\nElemento Piramide\n";
		}//TESTE
		if(iel == 4)
			cout << "\nElemento de Linha";
		if(iel == 5)
			cout << "\nElemento Quadrilatero\n";
		check = (CheckElement(gel) || check);
		check = (CheckRefinement(gel) || check);
	}
	return check;
}

int TPZCheckGeom::CheckRefinement(TPZGeoEl *gel){
	
	int check = 0;
	if(!gel || !gel->HasSubElement()) return check;
	int nsides = gel->NSides();
	int is;
	for(is=0; is<nsides; is++) {
		TPZStack<TPZGeoElSide> subel;
		gel->GetSubElements2(is,subel);
		int nsub = subel.NElements();
		int isub;
		for(isub=0; isub<nsub; isub++) {
			TPZGeoElSide fath = subel[isub].Father2();
			int son = subel[isub].Element()->WhichSubel();
			if(fath.Side() != is) {
				PZError << "TPZCheckGeom::CheckRefinement non corresponding subelement/sides son "
				<< son << " sonside " << subel[isub].Side() << " fathside " << is <<
				" fath2side " << fath.Side() << endl;
				gel->Print();
				check = 1;
			}
		}
	}
	int nsub = gel->NSubElements();
	for(is=0; is<nsub; is++) {
		TPZGeoEl *sub = gel->SubElement(is);
		int nsubsides = sub->NSides();
		int iss;
		for(iss=0; iss<nsubsides; iss++) {
			check = (CheckSubFatherTransform(sub,iss) || check);
		}
	}
	
	return check;
}

int TPZCheckGeom::CheckSideTransform(TPZGeoEl *gel, int sidefrom, int sideto){
	
	int check = 0;
	int nsides = gel->NSides();
	TPZIntPoints *integ = gel->CreateSideIntegrationRule(sidefrom,2);
	TPZTransform trans = gel->SideToSideTransform(sidefrom,sideto);
	TPZTransform trans1 = gel->SideToSideTransform(sidefrom,nsides-1);
	TPZTransform trans2 = gel->SideToSideTransform(sideto,nsides-1);
	int sidefromdim = gel->SideDimension(sidefrom);
	int sidetodim = gel->SideDimension(sideto);
	int geldim = gel->Dimension();
	TPZVec<REAL> intpoint(sidefromdim);
	TPZVec<REAL> sidetopoint(sidetodim);
	TPZVec<REAL> elpoint1(geldim),elpoint2(geldim);
	TPZVec<REAL> x1(3),x2(3);
	int nintpoints = integ->NPoints();
	int ip;
	REAL w;
	for(ip=0; ip<nintpoints; ip++) {
		integ->Point(ip,intpoint,w);
		trans.Apply(intpoint,sidetopoint);
		trans1.Apply(intpoint,elpoint1);
		trans2.Apply(sidetopoint,elpoint2);
		gel->X(elpoint1,x1);
		gel->X(elpoint2,x2);
		REAL dif = 0;
		int nx = x1.NElements();
		int ix;
		for(ix=0; ix<nx; ix++) dif += (x1[ix]-x2[ix])*(x1[ix]-x2[ix]);
		if(dif > 1.e-6) {
			PZError << "TPZCheckGeom::CheckSideTransform sidefrom = "<< sidefrom
			<< " sideto = " << sideto << " dif = " << dif << endl;
			gel->Print();
			check = 1;
		}
	}
	delete integ;
	return check;
}

int TPZCheckGeom::CheckSubFatherTransform(TPZGeoEl *subel, int sidesub) {
	int check = 0;
	TPZGeoElSide father = subel->Father2(sidesub);
	if(!father.Exists()) return check;
	TPZIntPoints *integ = subel->CreateSideIntegrationRule(sidesub,2);
	int subsidedim = subel->SideDimension(sidesub);
	int subdim = subel->Dimension();
	TPZTransform trans(subsidedim);
	trans = subel->BuildTransform2(sidesub,father.Element(),trans);
	int fathsidedim = father.Dimension();
	int fathdim = father.Element()->Dimension();
	int nsubsides = subel->NSides();
	int nfathsides = father.Element()->NSides();
	TPZTransform trans1 = subel->SideToSideTransform(sidesub,nsubsides-1);
	TPZTransform trans2 = father.Element()->SideToSideTransform(father.Side(),nfathsides-1);
	TPZVec<REAL> intpoint(subsidedim);
	TPZVec<REAL> sidetopoint(fathsidedim);
	TPZVec<REAL> elpoint1(subdim),elpoint2(fathdim);
	TPZVec<REAL> x1(3),x2(3);
	int nintpoints = integ->NPoints();
	int ip;
	REAL w;
	for(ip=0; ip<nintpoints; ip++) {
		integ->Point(ip,intpoint,w);
		trans.Apply(intpoint,sidetopoint);
		trans1.Apply(intpoint,elpoint1);
		trans2.Apply(sidetopoint,elpoint2);
		subel->X(elpoint1,x1);
		father.Element()->X(elpoint2,x2);
		int otherfatherside = father.Element()->WhichSide(elpoint2);
		if(otherfatherside != father.Side()) {
			int son = subel->WhichSubel();
			PZError << "TPZCheckGeom::CheckSubFatherTransform son " << son << " sidesub = "<< sidesub
			<< " fathside = " << father.Side() << " otherfatherside = " << otherfatherside << endl;
			check=1;
		}
		REAL dif = 0;
		long nx = x1.NElements();
		long ix;
		for(ix=0; ix<nx; ix++) dif += (x1[ix]-x2[ix])*(x1[ix]-x2[ix]);
		if(dif > 1.e-6) {
			int son = subel->WhichSubel();
			PZError << "TPZCheckGeom::CheckSubFatherTransform son " << son << " sidesub = "<< sidesub
			<< " fathside = " << father.Side() << " dif = " << dif << endl;
			//			subel->Print();
			check = 1;
			TPZTransform t = subel->ComputeParamTrans(father.Element(),father.Side(),sidesub);
			t.PrintInputForm(cout);
			cout << endl;
			trans.PrintInputForm(cout);
			cout << endl;
			check = 1;
		}
		
	}
	if(check == 0) {
		TPZTransform t = subel->ComputeParamTrans(father.Element(),father.Side(),sidesub);
		check = t.Compare(trans);
		if(check == 1){
			int son = subel->WhichSubel();
			PZError << "TPZCheckGeom::CheckSubFatherTransform son " << son << " sidesub = "<< sidesub
			<< " fathside = " << father.Side()  << endl;
			t.PrintInputForm(cout);
			cout << endl;
			trans.PrintInputForm(cout);
			cout << endl;
		}
		// compare t with trans
	}
	
	delete integ;
	return check;
	
}

#include "pzshapelinear.h"
#include "pzshapecube.h"
#include "pzshapepiram.h"
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"

int TPZCheckGeom::main() {
	TPZCheckGeom local;
	local.CreateMesh();
	ofstream meshfile("mesh.txt");
	local.fMesh->Print(meshfile);
	local.PerformCheck();
	return 1;
}

static REAL nodeco[12][3] = {
	{0,0,0},
	{1,0,0},
	{2,0,0},
	{0,1,0},
	{1,1,0},
	{2,1,0},
	{0,0,1},
	{1,0,1},
	{2,0,1},
	{0,1,1},
	{1,1,1},
	{2,1,1}
};

static int nodind[7][8] = {
	{0,1,4,3,6,7,10,9},
	{1,4,10,7,2},
	{8,7,10,2},
	{2,5,4,8,11,10},
	{0,1},
	{0,1,7,6},
	{1,2,7}
};

static int numnos[7] = {8,5,4,6,2,4,3};

void TPZCheckGeom::CreateMesh() {
	
	if(fMesh) delete fMesh;
	fMesh = new TPZGeoMesh();
	long noind[12];
	int no;
	for(no=0; no<12; no++) {
		noind[no] = fMesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord(3);
		coord[0] = nodeco[no][0];
		coord[1] = nodeco[no][1];
		coord[2] = nodeco[no][2];
		fMesh->NodeVec()[noind[no]].Initialize(coord,*fMesh);
	}
	int matid = 1;
	TPZVec<long> nodeindex;
	int nel;
	for(nel=0; nel<7; nel++) {
		int in;
		nodeindex.Resize(numnos[nel]);
		for(in=0; in<numnos[nel]; in++) {
			nodeindex[in] = nodind[nel][in];
		}
		long index;  
		switch(nel) {
			case 0:
				fMesh->CreateGeoElement(ECube, nodeindex, matid, index);
				break;
			case 1:
				fMesh->CreateGeoElement(EPiramide, nodeindex,matid, index);
				break;
			case 2:
				fMesh->CreateGeoElement(ETetraedro, nodeindex,matid, index);
				break;
			case 3:
				fMesh->CreateGeoElement(EPrisma, nodeindex,matid, index);
				break;
			case 4:
				fMesh->CreateGeoElement(EOned, nodeindex,matid, index);
				break;
			case 5:
				fMesh->CreateGeoElement(EQuadrilateral, nodeindex,matid, index);
				break;
			case 6:
				fMesh->CreateGeoElement(ETriangle, nodeindex,matid, index);
				break;
			default:
				break;
		}
	}
	fMesh->BuildConnectivity();
}
