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

TPZCheckGeom::TPZCheckGeom(TPZGeoMesh *gmesh) : fMesh(gmesh) {
}

int TPZCheckGeom::CheckElement(TPZGeoEl *gel) {
	
	int check = 0;
    check = check || CheckInternalTransforms(gel);
    check = check || CheckRefinement(gel);
    check = check || CheckNeighbourMap(gel);
	return check;
}

/// check the internal side transformations
int TPZCheckGeom::CheckInternalTransforms(TPZGeoEl *gel)
{
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

/// divide all elements and call PerformCheck
int TPZCheckGeom::DivideandCheck()
{
    int64_t nel = fMesh->NElements();
    for(int64_t iel = 0; iel<nel; iel++) {
        TPZGeoEl *gel = fMesh->ElementVec()[iel];
        if(!gel) continue;
        TPZStack<TPZGeoEl *> subel;
        gel->Divide(subel);		
    }
    return PerformCheck();
}

/*** @brief  Check if all node and elements ids are unique */
int TPZCheckGeom::CheckIds()
{
    int64_t numnodes   =   fMesh->NNodes();
    int64_t numels     =   fMesh->NElements();
    
    std::set<int64_t> nodeids;
    std::set<int64_t> elsids;
    
    for (int64_t inode = 0; inode < numnodes; inode++) {
        nodeids.insert(fMesh->NodeVec()[inode].Id());
    }
    if (nodeids.size() != numnodes) {
        std::cout << "The nodes have duplicate ids - EXPECT TROUBLE!\n";
        return 1;
    }
    
    
    for (int64_t iel = 0; iel < numels; iel++) {
        elsids.insert(fMesh->ElementVec()[iel]->Id());
    }
    if (elsids.size() != numels) {
        std::cout << "The elements have duplicate ids - EXPECT TROUBLE!\n";
        return 1;
    }
    
    return 0;
}

int TPZCheckGeom::PerformCheck() {
	int64_t nel = fMesh->NElements();
	int check = 0;
	for(int64_t iel = 0; iel<nel; iel++) {
		TPZGeoEl *gel = fMesh->ElementVec()[iel];
		if(!gel) continue;
		check = (CheckElement(gel) || check);
	}
    check = (CheckIds() || check);
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
	TPZTransform<> trans = gel->SideToSideTransform(sidefrom,sideto);
	TPZTransform<> trans1 = gel->SideToSideTransform(sidefrom,nsides-1);
	TPZTransform<> trans2 = gel->SideToSideTransform(sideto,nsides-1);
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
	TPZTransform<> trans(subsidedim);
	trans = subel->BuildTransform2(sidesub,father.Element(),trans);
	int fathsidedim = father.Dimension();
	int fathdim = father.Element()->Dimension();
	int nsubsides = subel->NSides();
	int nfathsides = father.Element()->NSides();
	TPZTransform<> trans1 = subel->SideToSideTransform(sidesub,nsubsides-1);
	TPZTransform<> trans2 = father.Element()->SideToSideTransform(father.Side(),nfathsides-1);
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
		int64_t nx = x1.NElements();
		int64_t ix;
		for(ix=0; ix<nx; ix++) dif += (x1[ix]-x2[ix])*(x1[ix]-x2[ix]);
		if(dif > 1.e-6) {
			int son = subel->WhichSubel();
			PZError << "TPZCheckGeom::CheckSubFatherTransform son " << son << " sidesub = "<< sidesub
			<< " fathside = " << father.Side() << " dif = " << dif << endl;
			//			subel->Print();
			check = 1;
			TPZTransform<> t = subel->ComputeParamTrans(father.Element(),father.Side(),sidesub);
			t.PrintInputForm(cout);
			cout << endl;
			trans.PrintInputForm(cout);
			cout << endl;
			check = 1;
		}
		
	}
	if(check == 0) {
		TPZTransform<> t = subel->ComputeParamTrans(father.Element(),father.Side(),sidesub);
		check = t.CompareTransform(trans);
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
	local.DivideandCheck();
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
	
	if(fMesh) DebugStop();
	fMesh = new TPZGeoMesh();
	int64_t noind[12];
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
	TPZVec<int64_t> nodeindex;
	int nel;
	for(nel=0; nel<7; nel++) {
		int in;
		nodeindex.Resize(numnos[nel]);
		for(in=0; in<numnos[nel]; in++) {
			nodeindex[in] = nodind[nel][in];
		}
		int64_t index;  
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

/// verify if the mapping between neighbouring elements is conforming
int TPZCheckGeom::CheckNeighbourMap(TPZGeoEl *gel)
{
    int check = 0;
    int nsides = gel->NSides();
    for (int is=0; is<nsides; is++) {
        int sidedim = gel->SideDimension(is);
        if (sidedim == 0) {
            continue;
        }
        int order = 4;
        TPZIntPoints *integ = gel->CreateSideIntegrationRule(is, order);
        int npoints = integ->NPoints();
        REAL w;
        TPZManVector<REAL,3> X1(3),X2(3);
        TPZGeoElSide gelside(gel,is);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside)
        {
            TPZTransform<> tr(sidedim);
            gelside.SideTransform3(neighbour, tr);
            TPZManVector<REAL,3> pt1(sidedim),pt2(sidedim);
            for (int ip = 0; ip < npoints; ip++) {
                integ->Point(ip, pt1, w);
                tr.Apply(pt1, pt2);
                gelside.X(pt1, X1);
                neighbour.X(pt2, X2);
                REAL norm = 0;
                for (int i=0; i<3; i++) {
                    norm += (X1[i]-X2[i])*(X1[i]-X2[i]);
                }
                if (norm > 1.e-12) {
                    std::cout << "Incompatible geometry between neighbours pt1 = " << pt1 << " pt2 = " << pt2 << " X1 = " << X1 << " X2 = " << X2 << "\n";
                    gelside.Print(std::cout);
                    neighbour.Print(std::cout);
                    check = 1;
                }
            }
            neighbour = neighbour.Neighbour();
        }
    }
    return check;
}

/// Verify is the ids of the elements and nodes are unique
void TPZCheckGeom::CheckUniqueId()
{
    std::set<int64_t> elementid, nodeid;
    int64_t nnode = fMesh->NNodes();
    for (int64_t in=0; in<nnode; in++) {
        int64_t id = fMesh->NodeVec()[in].Id();
        if (id != -1) {
            if (nodeid.find(id) == nodeid.end()) {
                nodeid.insert(id);
            }
            else
            {
                std::cout << "Node id repetido id = " << id << std::endl;
                DebugStop();
            }
        }
    }
    int64_t nelem = fMesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZGeoEl *gel = fMesh->Element(el);
        if (!gel) {
            continue;
        }
        int64_t id = gel->Id();
        if (elementid.find(id) != elementid.end()) {
            std::cout << "Duplicate element id = " << id << std::endl;
            DebugStop();
        }
        else
        {
            elementid.insert(id);
        }
    }
}

void TPZCheckGeom::UniformRefine(int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = fMesh->ElementVec();
        int nels = fMesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZManVector< TPZGeoEl *,20 > filhos;
            TPZGeoEl * gel = gelvec[elem];
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
        }
    }
}
