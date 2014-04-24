#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "pzelasmat.h" 
#include "pzelasthybrid.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

#include "pzgraphmesh.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;

const int matInterno = 1;
const int matCoarse = 2;

const int dirichlet = 0;

const int bc1 = -1;
const int bc2 = -2;
const int bc3 = -3;
const int bc4 = -4;


TPZGeoMesh *MalhaGeom2(REAL Lx, REAL Ly);
TPZCompMesh *MalhaCompTemporaria(TPZGeoMesh * gmesh);
TPZCompMesh *MalhaComp2(TPZGeoMesh * gmesh,int pOrder,std::set<long> coarseindex);

void RefinamentoUniforme(TPZGeoMesh *gmesh, int nref,TPZVec<int> dims);
void RefinamentoAdaptado(TPZGeoMesh *gmesh, TPZStack<TPZManVector<REAL,3> > coordcentro);

TPZCompMesh *SkeletonCoarseCompMesh (TPZCompMesh *cmesh, int matId);

void ChangeIndex(TPZGeoMesh *gmesh, int matcoarse1D);

void GetElIndexCoarseMesh(TPZGeoMesh * gmesh, std::set<long> &coarseindex);

int main(int argc, char *argv[])
{
    TPZGeoMesh * gmesh = MalhaGeom2(1.,1.);
	ofstream arg0("gmesh0.txt");
	gmesh->Print(arg0);
    
//-------- construindo malha coarse ----------
    
    //1 refinamento uniforme
    TPZVec<int> dims(2,0);
    dims[0]=1; dims[1]=2;
    RefinamentoUniforme(gmesh, 1, dims);
    ofstream arg1("gmesh1.txt");
	gmesh->Print(arg1);
    
    //refinamento uniforme em alguns elementos
    TPZStack<TPZManVector<REAL,3> > coordcenter;
    TPZManVector<REAL,3> xcoord(3,0.);
    
    xcoord[0]=0.25; xcoord[1]=0.75;
    coordcenter.Push(xcoord);
    
    xcoord[0]=0.25; xcoord[1]=1.0;
    coordcenter.Push(xcoord);
    
    xcoord[0]=0.0; xcoord[1]=0.75;
    coordcenter.Push(xcoord);
    
    xcoord[0]=0.75; xcoord[1]=0.25;
    coordcenter.Push(xcoord);
    
    xcoord[0]=1.0; xcoord[1]=0.25;
    coordcenter.Push(xcoord);
    
    xcoord[0]=0.75; xcoord[1]= 0.0;
    coordcenter.Push(xcoord);
    
    RefinamentoAdaptado(gmesh,coordcenter);
    ofstream arg2("gmesh2.txt");
	gmesh->Print(arg2);
    
    //construir elementos 1D de interface
    TPZCompMesh * cmesh1 = MalhaCompTemporaria(gmesh);
    gmesh->ResetReference();
    
    //mudar matId dos elementos 1D de interface
    ChangeIndex(gmesh, matCoarse);
    ofstream arg3("gmesh3.txt");
	gmesh->Print(arg3);
    
    ofstream file1("malhageometricaCoarse.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file1, true);
    
    //index dos elementos da malha coarse
    std::set<long> coarseindex;
    GetElIndexCoarseMesh(gmesh, coarseindex);
    
    if(coarseindex.find(6) != coarseindex.end())
    {
        std::cout << "\n\n\nNAO ACHEI O NUMERO DESEJADO\n\n\n";
    }
    
    
    std::set<long>::iterator it;
    std::cout << "coarse index: \n";
    for (it=coarseindex.begin(); it!=coarseindex.end(); ++it)
        std::cout << ' ' << *it;
    std::cout << '\n';

//-------- malha mais fina -------
    //refinamento uniforme dos elementos 2D
    dims.Resize(1, 0);
    dims[0]=2;
    RefinamentoUniforme(gmesh, 1, dims);
    ofstream arg4("gmesh4.txt");
	gmesh->Print(arg4);
    
    ofstream file2("malhageometricaFina.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file2, true);
    
//malha computacional
    TPZCompMesh * cmesh = MalhaComp2(gmesh,1,coarseindex);
    ofstream arg5("cmesh.txt");
	cmesh->Print(arg5);
    
    return EXIT_SUCCESS;
}

TPZGeoMesh *MalhaGeom2(REAL Lx, REAL Ly)
{
    int Qnodes = 4;
	long dim = 2;
    
	TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <long> TopolQuad(4);
	TPZVec <long> TopolLine(2);
	
	//indice dos nos
	long id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	//indice dos elementos
	id = 0;
    
    //elementos internos
    TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 2;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matInterno,*gmesh);
	id++;
    
    //elementos de contorno
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
	id++;
    
    //construir a malha
	gmesh->BuildConnectivity();
	
	return gmesh;
}

TPZCompMesh* MalhaCompTemporaria(TPZGeoMesh * gmesh)
{
	/// criar materiais
	int dim = 2;
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matInterno,dim);
	TPZMaterial * mat1(material);
	material->NStateVariables();
    
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(mat1);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	
    TPZMaterial * BCondD1 = material->CreateBC(mat1, bc2,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondD1);
	
	TPZMaterial * BCondD2 = material->CreateBC(mat1, bc4,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondD2);
    
	TPZMaterial * BCondN1 = material->CreateBC(mat1, bc1,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondN1);
    
    TPZMaterial * BCondN2 = material->CreateBC(mat1, bc3,dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCondN2);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->AutoBuild();
    
    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
    
    return cmesh;
}

TPZCompMesh* MalhaComp2(TPZGeoMesh * gmesh, int pOrder,std::set<long> coarseindex)
{
	/// criar materiais
	int dim = 2;
    TPZMatPoisson3d *material1 = new TPZMatPoisson3d(matInterno,dim);
    TPZMatPoisson3d *material2 = new TPZMatPoisson3d(matCoarse,dim);
    
	TPZMaterial * mat1(material1);
    TPZMaterial * mat2(material2);
    
	material1->NStateVariables();
    material2->NStateVariables();
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
    TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(mat1);
    cmesh->InsertMaterialObject(mat2);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	
    TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc2,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondD1);
	
	TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondD2);

	TPZMaterial * BCondN1 = material1->CreateBC(mat1, bc1,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondN1);
    
    TPZMaterial * BCondN2 = material1->CreateBC(mat1, bc3,dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCondN2);
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Criar elementos computacionais malha MHM
    int nel = gmesh->NElements();
    int matid, eldim;
    long index;
    int hassubel, nsubels;
    int iel, is;

    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    
    for(iel = 0; iel<nel; iel++)
    {
        gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        eldim = gel->Dimension();
        matid = gel->MaterialId();
        
        //elementos de dimensao = dim (malha fina)
        if(eldim==dim)
        {
            index = gel->Index();
            if(coarseindex.find(index) != coarseindex.end())
            {
                nsubels = gel->NSubElements();
                for (is=0; is<nsubels; is++)
                {
                    gsubel = gel->SubElement(is);
                    hassubel = gsubel->HasSubElement();
                    if(!hassubel){
                        cmesh->CreateCompEl(gsubel,index);
                    }
                }
                for (is=0; is<nsubels; is++)
                {
                    gsubel = gel->SubElement(is);
                    hassubel = gsubel->HasSubElement();
                    if(!hassubel){
                        gsubel->ResetReference();
                    }
                }
            }
            continue;
        }
        
        //elementos de dimensao = dim-1
        
        //malha coarse
        if(matid==matCoarse)
        {
            cmesh->CreateCompEl(gel, index);
            gel->ResetReference();
            
            continue;
        }
        
        //elementos de contorno
        hassubel =gel->HasSubElement();
        if(!hassubel)
        {
            cmesh->CreateCompEl(gel, index);
            gel->ResetReference();
        }
    }
    
    cmesh->LoadReferences();
    cmesh->ExpandSolution();
    
    return cmesh;
}


void RefinamentoUniforme(TPZGeoMesh *gmesh, int nref,TPZVec<int> dims)
{

    int ir, iel, k;
    int nel=0, dim=0;
    int ndims = dims.size();
	for(ir = 0; ir < nref; ir++ )
    {
		TPZVec<TPZGeoEl *> filhos;
        nel = gmesh->NElements();
        
		for (iel = 0; iel < nel; iel++ )
        {
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            
            dim = gel->Dimension();
            
            for(k = 0; k<ndims; k++)
            {
                if(dim == dims[k])
                {
                    gel->Divide (filhos);
                    break;
                }
            }
		}
	}
    
}

void RefinamentoAdaptado(TPZGeoMesh *gmesh, TPZStack<TPZManVector<REAL,3> > coordcentro)
{
    int size = coordcentro.NElements();
    std::set<TPZGeoEl *> setgeo;
    
    TPZVec<REAL> qsi(3,0.);
    long iniEl = 0;
    
    TPZStack<TPZGeoEl *> vecgel;
    vecgel.Resize(0);
    int eldim = 2;
    for(int ix = 0; ix < size; ix++)
    {
        TPZGeoEl * gel = NULL;
        
        if(coordcentro[ix][0]== 0. || coordcentro[ix][1]==0. || coordcentro[ix][0]== 1. || coordcentro[ix][1]==1.)
        {
            eldim = 1;
            gel = gmesh->FindElement(coordcentro[ix], qsi, iniEl, eldim);
            eldim = 2;
        }else gel = gmesh->FindElement(coordcentro[ix], qsi, iniEl, eldim);
        if(!gel) DebugStop();
        setgeo.insert(gel);
    }
    
    
    int nel = gmesh->NElements();
    TPZVec<TPZGeoEl *> filhos;
    for(int iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        std::set<TPZGeoEl *>::iterator found = setgeo.find(gel);
        
        if(gel==(*found)){
            continue;
        }
        gel->Divide (filhos);
    }
}

void GetElIndexCoarseMesh(TPZGeoMesh * gmesh, std::set<long> &coarseindex)
{
    int nel = gmesh->NElements();
    int iel;
    int hassubel=0;
    int dim = gmesh->Dimension();
    int eldim;
    for(iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        hassubel = gel->HasSubElement();
        eldim = gel->Dimension();
        if(!hassubel && eldim ==dim)
        {
            coarseindex.insert(gel->Index());
        }
    }
    
}

void ChangeIndex(TPZGeoMesh *gmesh, int matcoarse1D)
{
    int nel = gmesh->NElements();
    int iel;
    int hassubel, ninterf;
    
    for(iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        hassubel = gel->HasSubElement();
        if(!hassubel)
        {
            ninterf = gel->NumInterfaces();
            if(ninterf > 1) DebugStop();
            if (ninterf==1)
            {
                gel->SetMaterialId(matcoarse1D);
                gel->DecrementNumInterfaces();
            }
        }
    }

}

//TPZCompMesh *SkeletonCoarseCompMesh (TPZCompMesh *cmesh)
//{
//    
//    //ponteiro para a malha geometrica de mesh
//    TPZGeoMesh *gmesh = cmesh->Reference();
//    if(!gmesh)
//    {
//        DebugStop();
//    }
//    
//    //Reseta as referencias do ponteiro para a malha geometrica criada
//    //e criar uma nova malha computacional baseada nesta malha geometrica
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    //TPZCompMesh *newcmesh = new TPZCompMesh(gmesh);
//    
//    int dim = cmesh->Dimension();
//    //newcmesh->SetDimModel(dim);
//    
//    TPZStack<TPZCompElSide> neighequal,neighsmaller;
//    TPZCompElSide neighbigger;
//    int nneighs=0;
//    
//    int nel = cmesh->ElementVec().NElements();
//    for(long el = 0; el < nel; el++)
//	{
//		TPZCompEl * compEl = cmesh->ElementVec()[el];
//		if(!compEl) continue;
//        
//        int eldim = compEl->Reference()->Dimension();
//        int elmatId = compEl->Reference()->MaterialId();
//        
//        //convencao PZ: elemento de contorno tem sempre Id negativo
//		if(elmatId < 0)
//		{
//			compEl->Reference()->ResetReference();
//		}
//        else if(eldim == dim)
//        {
//            compEl->Reference()->ResetReference();
//            
//            int nsides = compEl->Reference()->NSides();
//            int ncorn = compEl->Reference()->NCornerNodes();
//            for(int side = ncorn; side < nsides; side++)
//            {
//                neighequal.Resize(0);
//                neighsmaller.Resize(0);
//                
//                TPZCompElSide celside(compEl,side);
//                celside.EqualLevelElementList(neighequal, 0, 0);
//                neighbigger = celside.LowerLevelElementList(0);
//                
//                if(neighbigger){
//                    neighequal.Push(neighbigger);
//                }
//                nneighs = neighequal.NElements();
//                
//                //                celside.HigherLevelElementList(neighsmaller, 1, 1);
//                //                int nneighsmaller = neighsmaller.size();
//                //                if(nneighs && nneighsmaller)
//                //                {
//                //                    DebugStop();
//                //                }
//                
//                if(nneighs != 0)
//                {
//                    //Loop on neighboring elements greater or equal level
//                    for(int i =0; i<nneighs; i++)
//                    {
//                        TPZGeoEl *geoside = neighequal[i].Reference().Element();
//                        
//                        //Do not assume neighbors by node
//                        if(neighequal[i].Side() < geoside->NCornerNodes()) continue;
//                        
//                        //verificando se eh elemento 1d
//                        if(geoside->Dimension() == dim-1 && geoside->MaterialId() > 0) continue;
//                        
//                        long index;
//                        TPZInterpolatedElement *newcompel;
//                        newcompel = dynamic_cast<TPZInterpolatedElement *> (cmesh->CreateCompEl(geoside, index));
//                        newcompel->LoadElementReference();
//                    }
//                }
//            }
//        }
//        else continue;
//	}
//
//    return cmesh;
//}

TPZCompMesh *SkeletonCoarseCompMesh (TPZCompMesh *cmesh, int matId)
{
    
    //ponteiro para a malha geometrica de mesh
    TPZGeoMesh *gmesh = cmesh->Reference();
    if(!gmesh)
    {
        DebugStop();
    }
    
    //Resetar as referencias do ponteiro para a malha geometrica criada
    //e criar uma nova malha computacional baseada nesta malha geometrica
    //gmesh->ResetReference();
    
//    TPZCompMesh *newcmesh = new TPZCompMesh(gmesh);
//    gmesh->ResetReference();
//    newcmesh->LoadReferences();
    
    int dim = cmesh->Dimension();
    //newcmesh->SetDimModel(dim);
    
    TPZStack<TPZCompElSide> neighequal,neighsmaller;
    TPZCompElSide neighbigger;
    int nneighs=0;
    
    int nel = gmesh->ElementVec().NElements();
    for(long el = 0; el < nel; el++)
	{
		TPZGeoEl * gel = gmesh->ElementVec()[el];
		if(!gel) continue;
        if(gel->HasSubElement()) continue;
        
        int eldim = gel->Dimension();
        //int elmatId = gel->MaterialId();
        
        //        //convencao PZ: elemento de contorno tem sempre Id negativo
        //		if(elmatId < 0)
        //		{
        //			geo->Reference()->ResetReference();
        //		}
        if(eldim == dim)
        {
            //compEl->Reference()->ResetReference();
            
            int nsides = gel->NSides();
            int ncorn = gel->NCornerNodes();
            for(int side = ncorn; side < nsides; side++)
            {
                neighequal.Resize(0);
                neighsmaller.Resize(0);
                
                TPZGeoElSide gelside(gel,side);
                gelside.EqualLevelCompElementList(neighequal, 0, 0);
                neighbigger = gelside.LowerLevelCompElementList2(0);
                
                if(neighbigger){
                    neighequal.Push(neighbigger);
                }
                nneighs = neighequal.NElements();
                
//                gelside.HigherLevelCompElementList2(neighsmaller, 1, 1);
//                int nneighsmaller = neighsmaller.size();
//                if(nneighs && nneighsmaller)
//                {
//                    DebugStop();
//                }
                
                if(nneighs != 0)
                {
                    //Loop on neighboring elements greater or equal level
                    for(int i =0; i<nneighs; i++)
                    {
                        TPZGeoEl *geoside = neighequal[i].Reference().Element();
                        
                        //Do not assume neighbors by node
                        if(neighequal[i].Side() < geoside->NCornerNodes()) continue;
                        
                        //verificando se eh elemento 1d
                        if(geoside->Dimension() == dim-1 && geoside->MaterialId() > 0) continue;
                        
                        long index;
                        TPZInterpolatedElement *newcompel;
                        geoside->ResetReference();
                        newcompel = dynamic_cast<TPZInterpolatedElement *> (cmesh->CreateCompEl(geoside, index));
                        newcompel->Reference()->SetMaterialId(matId);
                        newcompel->LoadElementReference();
                        newcompel->Print();
                        cout<<"\n";
                        newcompel->Reference()->Print();
                    }
                }
            }
            gel->ResetReference();
        }
        else continue;
	}
    //cmesh->LoadReferences();
    cmesh->InitializeBlock();
    return cmesh;
}