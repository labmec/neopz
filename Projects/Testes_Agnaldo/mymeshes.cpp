//
//  mymeshes.cpp
//  PZ
//
//  Created by Agnaldo Farias on 29/04/13.
//
//

#include "mymeshes.h"

#include "pzgeoelbc.h"
#include "TPZCompElDisc.h"
#include "pzpoisson3d.h"
#include "pzelasmat.h"
#include "pzgeoelbc.h"
#include "TPZSkylineNSymStructMatrix.h"


DadosMalhas::DadosMalhas(){
    
    //dados do material
    fEyoung = 0.;
    fpoisson= 0.;
    falpha = 0.;
    fSe = 0.;
    fperm = 0.;
    fvisc = 0.;
    ffx = 0.;
    ffy = 0.;
    fsign = 0.;
    
    fpref= 0.;
    fLref = 0.;
    fkovervisc = 0.;
    
    fvalsourceterm = 0.;
    
    //dados da malha geometrica
    fmatId =1;
    fbcBottom = -1;
    fbcRight = -2;
    fbcTop = -3;
    fbcLeft = -4;
    fbcTopStripLoad = -5;
    fbcSourceTerm = 2;
    
    
    fdirichlet =0;
    fneumann = 11;
    fneumdir=10;
    fdirfreey_neum=300;
    fdirneum = 1;
    fmixedneum = 21;
    fmixeddirich = 20;
    
    fmixedFreeYXdirich = 400;
    fmixedFreeXYdirich = 500;
}

DadosMalhas::~DadosMalhas(){
    
}

DadosMalhas::DadosMalhas(const DadosMalhas &copy){
    
    this->operator=(copy);
}

DadosMalhas & DadosMalhas::operator=(const DadosMalhas &copy){
    
    fEyoung  = copy. fEyoung;
	fpoisson = copy.fpoisson;
	falpha  = copy.falpha;
    fSe  = copy.fSe;
    fperm = copy.fperm ;
    fvisc = copy.fvisc;
    ffx = copy.ffx;
    ffy = copy.ffy;
    fsign = copy.fsign;
    
    
    fmatId = copy.fmatId;
    fbcBottom = copy.fbcBottom;
    fbcRight = copy.fbcRight;
    fbcTop = copy.fbcTop;
    fbcLeft = copy.fbcLeft;
    fbcTopStripLoad = copy.fbcTopStripLoad;
    fbcSourceTerm = copy.fbcSourceTerm;
    
    fdirichlet = copy.fdirichlet;
    fneumann = copy.fneumann;
    fneumdir = copy.fneumdir;
    fdirfreey_neum = copy.fdirfreey_neum;
    fdirneum = copy.fdirneum;
    fmixedneum = copy.fmixedneum;
    fmixeddirich = copy.fmixeddirich;
    
    fmixedFreeXYdirich = copy.fmixedFreeXYdirich;
    fmixedFreeYXdirich = copy.fmixedFreeYXdirich;
    
    fpref= copy.fpref;
    fLref = copy.fLref;
    fkovervisc = copy.fkovervisc;
    
    fvalsourceterm = copy.fvalsourceterm;
    
	return *this;
}


TPZGeoMesh *DadosMalhas::GMesh(bool triang_elements, REAL L, REAL w){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*L;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = L - xi*L;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,w);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    if(triang_elements==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,fmatId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,fmatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcBottom,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcRight,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcTop,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcLeft,*gmesh);
//        
//        TopolPoint[0] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolLine,fmatId+1,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,fmatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcBottom,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcRight,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcTop,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcLeft,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Geometrica Inicial\n ";
    //        gmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
	return gmesh;
}

TPZGeoMesh * DadosMalhas::GMesh2(REAL L, REAL w, REAL La){
    
    int nrefx = L/La;
    int nrefy = w/La;
    int nnodesx = (nrefx+1);
    int nnodesy = (nrefy+1);
    int Qnodes = nnodesx*nnodesy;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx;
    REAL valy;
    
    for(int yi=0; yi< nnodesy; yi++)
    {
        valy = La*yi;
        for(int xi = 0; xi< nnodesx; xi++)
        {
            valx = La*xi;
            
            Node[id].SetNodeId(id);
            Node[id].SetCoord(0 ,valx);//coord X
            Node[id].SetCoord(1 ,valy);//coord Y
            gmesh->NodeVec()[id] = Node[id];
            id++;
        }
    }
    
   	
	//indice dos elementos
	id = 0;
    int elx, ely;
    for(ely=0; ely<nrefy; ely++)
    {
        for(elx = 0; elx<nrefx; elx++)
        {
            TopolQuad[0] = elx + ely*nnodesx;
            TopolQuad[1] = TopolQuad[0]+1;
            TopolQuad[2] = TopolQuad[1] + nnodesx;
            TopolQuad[3] = TopolQuad[2]-1;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,fmatId,*gmesh);
            id++;
        }
    }
   
    //boundary elements
    for(elx = 0; elx<nrefx; elx++)
    {
        TopolLine[0] = elx;
        TopolLine[1] = TopolLine[0]+1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcBottom,*gmesh);
        id++;
    }
    
    for(ely = 0; ely<nrefy; ely++)
    {
        TopolLine[0] = nrefx + ely*nnodesx;
        TopolLine[1] = TopolLine[0] + nnodesx;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcRight,*gmesh);
        id++;
    }
    
    for(elx = 0; elx<nrefx-1; elx++)
    {
        TopolLine[0] = (Qnodes-1) - elx;
        TopolLine[1] = TopolLine[0]-1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcTop,*gmesh);
        id++;
    }
    
    {
        TopolLine[0] = (Qnodes-1) - elx;
        TopolLine[1] = TopolLine[0]-1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcTopStripLoad,*gmesh);
        id++;
    }
    
    for(ely = 0; ely<nrefy; ely++)
    {
        TopolLine[0] = ely*nnodesx;
        TopolLine[1] = TopolLine[0] + nnodesx;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbcLeft,*gmesh);
        id++;
    }
    
    
	gmesh->BuildConnectivity();
    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Geometrica Inicial\n ";
    //        gmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    return gmesh;
}


TPZGeoMesh *DadosMalhas::GMesh3(REAL L, REAL w){
    
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
	REAL co[9][2] = {{0.,0.},{L/2,0.},{L,0.},{0.,w/2},{L/2,w/2},{L,w/2},{0.,w},{L/2,w},{L,w}};
//	int64_t indices[1][9] = {{0,1,2,3,4,5,6,7,8}};
	
	int nnode = 9;
	const int nelem = 5;
	TPZGeoEl *elvec[nelem];
	int64_t nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}

    int64_t el;
	for ( el=0; el<nelem; el++ )
	{
        if(el!=1 && el!=2){
            TPZVec<int64_t> nodind(4);
            nodind[0] = el;
            nodind[1] = el + 1;
            nodind[2] = el + 4;
            nodind[3] = el + 3;
            elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,fmatId,el);
        }
        else if (el==1){
            TPZVec<int64_t> nodind(3);
            nodind[0] = 1;
            nodind[1] = 2;
            nodind[2] = 4;
            elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,fmatId,el);
        }else{
            TPZVec<int64_t> nodind(3);
            nodind[0] = 5;
            nodind[1] = 2;
            nodind[2] = 4;
            elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,fmatId,el);
        }
	}
    
    gmesh->BuildConnectivity();
    
    //Cricao das condicoes de contorno
    TPZGeoElBC gbc1(elvec[0], 4,fbcBottom);
    TPZGeoElBC gbc2(elvec[1], 3,fbcBottom);
    TPZGeoElBC gbc3(elvec[2], 3,fbcRight);
    TPZGeoElBC gbc4(elvec[4], 5,fbcRight);
    TPZGeoElBC gbc5(elvec[4], 6,fbcTop);
    TPZGeoElBC gbc6(elvec[3], 6,fbcTop);
    TPZGeoElBC gbc7(elvec[3], 7,fbcLeft);
    TPZGeoElBC gbc8(elvec[0], 7,fbcLeft);
    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Geometrica Inicial\n ";
    //        gmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
	return gmesh;
}

TPZGeoMesh *DadosMalhas::GMesh4(REAL L, REAL w, int h, int nrefdir){
    
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
	REAL co[4][2] = {{0.,0.},{L,0.},{L,w},{0.,w}};
//	int indices[1][4] = {{0,1,2,3}};
	
	int nnode = 4;
	const int nelem = 1;
	TPZGeoEl *elvec[nelem];
	int nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}
    
    int64_t el;
	for ( el=0; el<nelem; el++ )
	{
        TPZVec<int64_t> nodind(4);
        nodind[0] = 0;
        nodind[1] = 1;
        nodind[2] = 2;
        nodind[3] = 3;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,fmatId,el);
    }
    
    gmesh->BuildConnectivity();
    
    //Cricao das condicoes de contorno
    TPZGeoElBC gbc1(elvec[0], 4,fbcBottom);
    TPZGeoElBC gbc2(elvec[0], 5,fbcRight);
    TPZGeoElBC gbc5(elvec[0], 6,fbcTop);
    TPZGeoElBC gbc8(elvec[0], 7,fbcLeft);
    
    for(int D = 0; D < h; D++)
    {
        int64_t nels = gmesh->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    //	gmesh->BuildConnectivity();
    
    int nrefUnif = nrefdir;
    for(int ref = 0; ref < nrefUnif; ref++)
    {
        int nelem = gmesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            if(gmesh->ElementVec()[el]->Dimension() < 1) continue;
            if(gmesh->ElementVec()[el]->HasSubElement()) continue;
            if(gmesh->ElementVec()[el]->MaterialId() == fbcBottom)
            {
                TPZVec<TPZGeoEl*> sons;
                gmesh->ElementVec()[el]->Divide(sons);
                continue;
            }
            //else...
            for(int s = 0; s < gmesh->ElementVec()[el]->NSides(); s++)
            {
                TPZGeoElSide gelside(gmesh->ElementVec()[el],s);
                TPZGeoElSide neighside(gelside.Neighbour());
                bool refinedAlready = false;
                while(neighside != gelside)
                {
                    if(neighside.Element()->MaterialId() == fbcBottom)
                    {
                        TPZVec<TPZGeoEl*> sons;
                        gmesh->ElementVec()[el]->Divide(sons);
                        refinedAlready = true;
                        break;
                    }
                    neighside = neighside.Neighbour();
                }
                if(refinedAlready == true)
                {
                    break;
                }
            }
        }
    }

    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Geometrica Inicial\n ";
    //        gmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
	return gmesh;
}


void DadosMalhas::UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int64_t nels = gmesh->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    //	gmesh->BuildConnectivity();
}


TPZCompMesh* DadosMalhas:: MalhaCompElast(TPZGeoMesh * gmesh,int pOrder, bool twomaterial, bool stripload)
{
    /// criar material
	int dim = 2;
	TPZVec<REAL> force(dim,0.);
    //REAL E = 0;
	//REAL poisson = 0;
    int planestress = -1;
	TPZElasticityMaterial *material;
    
    
	material = new TPZElasticityMaterial(fmatId, fEyoung, fpoisson, ffx, ffy, planestress);
	TPZMaterial * mat(material);
        
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    if(twomaterial==true){
        TPZElasticityMaterial *material2 = new TPZElasticityMaterial(fmatId+1, fEyoung, fpoisson, ffx, ffy, planestress);
        TPZMaterial * mat2(material2);
        cmesh->InsertMaterialObject(mat2);
    }

	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1 = material->CreateBC(mat, fbcBottom,fdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, fbcLeft,fdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, fbcTop,fdirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, fbcRight,fdirichlet, val1, val2);
    TPZMaterial * BCond5 = material->CreateBC(mat, fbcTopStripLoad,fdirichlet, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);
	cmesh->InsertMaterialObject(BCond4);
    if(stripload) cmesh->InsertMaterialObject(BCond5);
	
	
    if(twomaterial==true){
        
        //Inserir materiais
        std::set<int> MaterialIDs;
        MaterialIDs.insert(fmatId);
        MaterialIDs.insert(fmatId+1);
        MaterialIDs.insert(fbcBottom);
        MaterialIDs.insert(fbcRight);
        MaterialIDs.insert(fbcTop);
        MaterialIDs.insert(fbcLeft);
        
        cmesh->AutoBuild(MaterialIDs);
        
        //AQUI: AutoBuild(mat1)
        set<int> SETmat1;
//        SETmat1.insert(fmatId);
//        SETmat1.insert(fbcBottom);
//        SETmat1.insert(fbcRight);
//        SETmat1.insert(fbcTop);
//        SETmat1.insert(fbcLeft);
//        
//        cmesh->AutoBuild(SETmat1);
//        gmesh->ResetReference();
//        
//        //AQUI: AutoBuild(mat2)
//        set<int> SETmat2;
//        SETmat2.insert(fmatId+1);
//        cmesh->AutoBuild(SETmat2);
//        cmesh->LoadReferences();
        
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
        return cmesh;
    }

    
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}

TPZCompMesh *DadosMalhas::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, bool twomaterial, bool stripload)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(fmatId,dim);
	TPZMaterial * mat(material);
	material->NStateVariables();
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    if(twomaterial==true){
        TPZMatPoisson3d *material2 = new TPZMatPoisson3d(fmatId+1,dim);
        TPZMaterial * mat2(material2);
        cmesh->InsertMaterialObject(mat2);
    }
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond1 = material->CreateBC(mat, fbcBottom,fdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, fbcRight,fdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, fbcTop,fdirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, fbcLeft,fdirichlet, val1, val2);
    TPZMaterial * BCond5 = material->CreateBC(mat, fbcTopStripLoad,fdirichlet, val1, val2);
    
	cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    if(stripload) cmesh->InsertMaterialObject(BCond5);
    
    
    if(twomaterial==true) {
        
        //Inserir materiais
        std::set<int> MaterialIDs;
        MaterialIDs.insert(fmatId);
        MaterialIDs.insert(fmatId+1);
        MaterialIDs.insert(fbcBottom);
        MaterialIDs.insert(fbcRight);
        MaterialIDs.insert(fbcTop);
        MaterialIDs.insert(fbcLeft);
        
        cmesh->AutoBuild(MaterialIDs);
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
        return cmesh;
    }
    
    //Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}

TPZCompMesh *DadosMalhas::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, bool triang, bool twomaterial)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(fmatId,dim);
	material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    if(twomaterial==true) {
        TPZMatPoisson3d *material2 = new TPZMatPoisson3d(fmatId+1,dim);
        TPZMaterial * mat2(material2);
        cmesh->InsertMaterialObject(mat2);
    }
    
    ///Inserir condicao de contorno
//	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//	TPZMaterial * BCond1 = material->CreateBC(mat, fbcBottom,fdirichlet, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(mat, fbcRight,fdirichlet, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(mat, fbcTop,fdirichlet, val1, val2);
//    TPZMaterial * BCond4 = material->CreateBC(mat, fbcLeft,fdirichlet, val1, val2);
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
//    cmesh->SetAllCreateFunctionsContinuous();
//    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
//    cmesh->InsertMaterialObject(BCond1);
//    cmesh->InsertMaterialObject(BCond2);
//    cmesh->InsertMaterialObject(BCond3);
//    cmesh->InsertMaterialObject(BCond4);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
    ///inserir connect da pressao
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
	    newnod.SetLagrangeMultiplier(1);
    }
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) continue;
        celdisc->SetConstC(1.);
        celdisc->SetCenterPoint(0, 0.);
        celdisc->SetCenterPoint(1, 0.);
        celdisc->SetCenterPoint(2, 0.);
        celdisc->SetTrueUseQsiEta();
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(triang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }
    
    
#ifdef PZDEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    if(twomaterial==true){
        
        //Inserir materiais
        std::set<int> MaterialIDs;
        MaterialIDs.insert(fmatId);
        MaterialIDs.insert(fmatId+1);
        MaterialIDs.insert(fbcBottom);
        MaterialIDs.insert(fbcRight);
        MaterialIDs.insert(fbcTop);
        MaterialIDs.insert(fbcLeft);
        
        cmesh->AutoBuild(MaterialIDs);
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
        return cmesh;
    }

    
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	return cmesh;
}

TPZCompMesh *DadosMalhas::MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZPoroElasticMF2d * &mymaterial,TPZAutoPointer<TPZFunction<STATE> > solExata){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=fmatId;
    int dim =2;
    
    //	criar material
    int planestress = 0; // This is a Plain strain problem
    mymaterial = new TPZPoroElasticMF2d(MatId,dim);
    mymaterial->SetfPlaneProblem(planestress);
    
    mymaterial->SetParameters(fperm, fvisc);
    mymaterial->SetElasticityParameters(fEyoung, fpoisson, falpha, fSe, ffx, ffy);
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
	mymaterial->SetExactSol(solExata);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    REAL sig0 = fsign;
    REAL ptop = 0.;
    val2(0,0)= 0.;
    val2(1,0)= sig0;
    val2(2,0)= ptop;
	TPZMaterial * BCond1 = mymaterial->CreateBC(mat, fbcTop,fneumdir, val1, val2);
    
    val2.Redim(3,1);
    REAL big = mymaterial->gBigNumber;
    val1(0,0) = big;
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat,fbcRight, fmixedneum, val1, val2);
    TPZMaterial * BCond4 = mymaterial->CreateBC(mat,fbcLeft, fmixedneum, val1, val2);
    
    val1.Redim(3,2);
    val2.Redim(3,1);
    //val2(2,0)= piniD;
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat,fbcBottom,fdirneum, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

TPZCompMesh *DadosMalhas::MalhaCompTerzaghi(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZPoroElasticMF2d * &mymaterial,TPZAutoPointer<TPZFunction<STATE> > solExata){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=fmatId;
    int dim =2;
    
    int planestress = 0; // This is a Plain strain problem
    
    
    mymaterial = new TPZPoroElasticMF2d(MatId,dim);
    
    mymaterial->SetfPlaneProblem(planestress);
    mymaterial->SetParameters(fperm, fvisc);
    mymaterial->SetElasticityParameters(fEyoung, fpoisson, falpha, fSe, ffx, ffy);
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
	mymaterial->SetExactSol(solExata);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    REAL sig0 = fsign;
    REAL ptop = 0.;
    val2(0,0)= 0.;
    val2(1,0)= sig0;
    val2(2,0)= ptop;
	TPZMaterial * BCond1 = mymaterial->CreateBC(mat, fbcTop,fneumdir, val1, val2);
    
    val2.Redim(3,1);
    REAL big = mymaterial->gBigNumber;
    val1(0,0) = big;
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat,fbcRight, fmixedneum, val1, val2);
    TPZMaterial * BCond4 = mymaterial->CreateBC(mat,fbcLeft, fmixedneum, val1, val2);
    
    val1.Redim(3,2);
    val2.Redim(3,1);
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat,fbcBottom,fdirneum, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

TPZCompMesh *DadosMalhas::MalhaCompBarryMercer(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,TPZAutoPointer<TPZFunction<STATE> > sourceterm, TPZAutoPointer<TPZFunction<STATE> > solExata){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=fmatId;
    int dim =2;
    
    int planestress = 0; // This is a Plain strain problem
    
    TPZPoroElasticMF2d *mymaterial1 = new TPZPoroElasticMF2d(MatId,dim);
    TPZPoroElasticMF2d * mymaterial2 = new TPZPoroElasticMF2d(MatId+1,dim);
    
    mymaterial1->SetfPlaneProblem(planestress);
    mymaterial1->SetParameters(fperm, fvisc);
    mymaterial1->SetElasticityParameters(fEyoung, fpoisson, falpha, fSe, ffx, ffy);
    mymaterial1->SetExactSol(solExata);
    
    mymaterial2->SetfPlaneProblem(planestress);
    mymaterial2->SetParameters(fperm, fvisc);
    mymaterial2->SetElasticityParameters(fEyoung, fpoisson, falpha, fSe, ffx, ffy);
    mymaterial2->SetForcingFunction(sourceterm);
    mymaterial2->SetExactSol(solExata);
   
    
    TPZMaterial *mat1(mymaterial1);
    mphysics->InsertMaterialObject(mat1);
    
    TPZMaterial *mat2(mymaterial2);
    mphysics->InsertMaterialObject(mat2);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    REAL big = mymaterial1->gBigNumber;
    
    //ux=0 (Dirichlet), p = 0(Dirichlet), duy/dy=0 (Neumann, nao preciso aplicar pois sao iguais a zero)
    val1(0,0) = big;
    TPZMaterial * BCond1 = mymaterial1->CreateBC(mat1,fbcBottom, fmixeddirich, val1, val2);
    TPZMaterial * BCond3 = mymaterial1->CreateBC(mat1,fbcTop, fmixeddirich, val1, val2);
    
    //uy=0 (Dirichlet),  p = 0(Dirichlet), dux/dx=0, (Neumann, nao preciso aplicar pois sao iguais a zero)
    val1(0,0) = 0.;
    val1(1,1) = big;
    TPZMaterial * BCond2 = mymaterial1->CreateBC(mat1,fbcRight, fmixeddirich, val1, val2);
    TPZMaterial * BCond4 = mymaterial1->CreateBC(mat1,fbcLeft, fmixeddirich, val1, val2);
    
//     TPZMaterial * BCond1 = mymaterial1->CreateBC(mat1,fbcBottom, fmixedFreeYXdirich, val1, val2);
//     TPZMaterial * BCond3 = mymaterial1->CreateBC(mat1,fbcTop, fmixedFreeYXdirich, val1, val2);
//     TPZMaterial * BCond2 = mymaterial1->CreateBC(mat1,fbcRight, fmixedFreeXYdirich, val1, val2);
//     TPZMaterial * BCond4 = mymaterial1->CreateBC(mat1,fbcLeft, fmixedFreeXYdirich, val1, val2);
    
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    
    
    //Inserir materiais
    std::set<int> MaterialIDs;
    MaterialIDs.insert(fmatId);
    MaterialIDs.insert(fmatId+1);
    MaterialIDs.insert(fbcBottom);
    MaterialIDs.insert(fbcRight);
    MaterialIDs.insert(fbcTop);
    MaterialIDs.insert(fbcLeft);
    
    mphysics->AutoBuild(MaterialIDs);
  	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}


TPZCompMesh *DadosMalhas::MalhaCompBarryMercerPressureSolution(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZPoroElasticMF2d * &mymaterial, TPZAutoPointer<TPZFunction<STATE> > BCterm, TPZAutoPointer<TPZFunction<STATE> > solExata){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=fmatId;
    int dim =2;
    
    int planestress = 0; // This is a Plain strain problem
    mymaterial = new TPZPoroElasticMF2d(MatId,dim);
    
    mymaterial->SetfPlaneProblem(planestress);
    mymaterial->SetParameters(fperm, fvisc);
    mymaterial->SetElasticityParameters(fEyoung, fpoisson, falpha, fSe, ffx, ffy);
    mymaterial->SetExactSol(solExata);
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    REAL big = mymaterial->gBigNumber;
    
    //ux=0 (Dirichlet), p = 0(Dirichlet), duy/dy=0 (Neumann, nao preciso aplicar pois sao iguais a zero)
    
    TPZMaterial * BCond1 = mymaterial->CreateBC(mat,fbcBottom, fdirichlet, val1, val2);
    BCond1->SetForcingFunction(BCterm);
    val1(0,0) = big;
//    TPZMaterial * BCond3 = mymaterial->CreateBC(mat,fbcTop, fmixedneum, val1, val2);
//    BCond3->SetForcingFunction(BCterm);
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat,fbcTop, fmixeddirich, val1, val2);
    
    //uy=0 (Dirichlet),  p = 0(Dirichlet), dux/dx=0, (Neumann, nao preciso aplicar pois sao iguais a zero)
    val1(0,0) = 0.;
    val1(1,1) = big;
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat,fbcRight, fmixeddirich, val1, val2);
    TPZMaterial * BCond4 = mymaterial->CreateBC(mat,fbcLeft, fmixeddirich, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    
    mphysics->AutoBuild();
  	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

#include "pzbstrmatrix.h"

void DadosMalhas::SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	TPZBandStructMatrix full(fCmesh);
	//TPZSkylineStructMatrix<STATE> full(fCmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
    //	step.SetDirect(ELDLt); //caso simetrico
    step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
}

TPZAutoPointer <TPZMatrix<STATE> > DadosMalhas::MassMatrix(TPZPoroElasticMF2d * mymaterial, TPZCompMesh* mphysics, int nthreads){
    
    mymaterial->SetLastState();
    //TPZSkylineStructMatrix<STATE> matsp(mphysics);
    //TPZSkylineNSymStructMatrix matsp(mphysics);
	TPZSpStructMatrix matsp(mphysics);
    matsp.SetNumThreads(nthreads);
    
	std::set< int > materialid;
	int matid = mymaterial->MatId();
	materialid.insert(matid);
	matsp.SetMaterialIds (materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<STATE> Un;
    
    
    //TPZMatrix<REAL> *matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}

TPZAutoPointer <TPZMatrix<STATE> > DadosMalhas::MassMatrix(TPZCompMesh* mphysics, int nthreads){
    
    TPZMaterial * mat1 = mphysics->FindMaterial(fmatId);
    TPZMaterial * mat2 = mphysics->FindMaterial(fmatId+1);
    
    TPZPoroElasticMF2d * mat12 = dynamic_cast<TPZPoroElasticMF2d *>(mat1);
    TPZPoroElasticMF2d * mat22 = dynamic_cast<TPZPoroElasticMF2d *>(mat2);
    
    mat12->SetLastState();
    mat22->SetLastState();
    
    //TPZSkylineStructMatrix<STATE> matsp(mphysics);
	TPZSpStructMatrix matsp(mphysics);
    matsp.SetNumThreads(nthreads);
	std::set< int > materialid;
	materialid.insert(fmatId);
    materialid.insert(fmatId+1);
	matsp.SetMaterialIds (materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<STATE> Un;
    //TPZMatrix<REAL> *matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}


void DadosMalhas::StiffMatrixLoadVec(TPZPoroElasticMF2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<STATE> &matK1, TPZFMatrix<STATE> &fvec, int nthreads) {
    
    mymaterial->SetCurrentState();
    //TPZFStructMatrix<STATE> matsk(mphysics);
    TPZSkylineStructMatrix<STATE> matsk(mphysics);
    matsk.SetNumThreads(nthreads);
	an.SetStructuralMatrix(matsk);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	matK1 = an.StructMatrix();
	fvec = an.Rhs();
    
    
    
    //    TPZBandStructMatrix full(cmesh);
    //    an.SetStructuralMatrix(full);
    //    TPZStepSolver step;
    //    step.SetDirect(ELU);
    //    an.SetSolver(step);
}

void DadosMalhas::StiffMatrixLoadVec(TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<STATE> &matK1, TPZFMatrix<STATE> &fvec, int nthreads){
    
	TPZMaterial * mat1 = mphysics->FindMaterial(fmatId);
    TPZMaterial * mat2 = mphysics->FindMaterial(fmatId+1);
    
    TPZPoroElasticMF2d * mat12 = dynamic_cast<TPZPoroElasticMF2d *>(mat1);
    TPZPoroElasticMF2d * mat22 = dynamic_cast<TPZPoroElasticMF2d *>(mat2);
    
	mat12->SetCurrentState();
    mat22->SetCurrentState();

    //TPZFStructMatrix<STATE> matsk(mphysics);
    TPZSkylineStructMatrix<STATE> matsk(mphysics);
	an.SetStructuralMatrix(matsk);
    matsk.SetNumThreads(nthreads);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	matK1 = an.StructMatrix();
	fvec = an.Rhs();
    
    //    TPZBandStructMatrix full(cmesh);
    //    an.SetStructuralMatrix(full);
    //    TPZStepSolver step;
    //    step.SetDirect(ELU);
    //    an.SetSolver(step);
}


void DadosMalhas::PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
{
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    
	TPZManVector<std::string,10> scalnames(3), vecnames(2);
	//scalnames[0] = "DisplacementX";
	//scalnames[1] = "DisplacementY";
    //vecnames[0] = "Displacement";
    //scalnames[2] = "SigmaX";
	//scalnames[3] = "SigmaY";
	scalnames[0] = "PorePressure";
    scalnames[2] = "FluxoY";
	vecnames[0] = "Fluxo";
	//vecnames[1] = "MinusKMuGradP";
    
    scalnames[1] = "ExactPressure";
    //scalnames[6] = "FluxoX";
    
    //scalnames[4] = "ExactDisplacementX";
    //scalnames[5] = "ExactDisplacementY";
   // scalnames[8] = "ExactSigmaX";
    //scalnames[6] = "ExactSigmaY";
    vecnames[1]  = "ExactFluxo";
    //vecnames[1]  = "ExactDisplacement";
    //vecnames[4] = "MinusKMuGradP";
	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
//	std::ofstream out("malha.txt");
//	an.Print("nothing",out);
}

void DadosMalhas::SolveSistTransient(REAL deltaT,REAL maxTime, TPZPoroElasticMF2d * &mymaterial ,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, int ntimestep, REAL &timeatual){
	
    
    TPZAnalysis an(mphysics);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    std::string outputfile;
	outputfile = "TransientSolution";
    
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    PosProcessMultphysics(meshvec,mphysics,an,plotfile);
    
    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(mymaterial, mphysics, 8);
    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //            std::stringstream sout;
    //        	matM->Print("matM = ", sout,EMathematicaInput);
    //        	LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    //Criando matriz de rigidez (matK) e vetor de carga
	TPZFMatrix<STATE> matK;
	TPZFMatrix<STATE> fvec;
    StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fvec,8);
    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //
    //        std::stringstream sout;
    //        matK.Print("matK = ", sout,EMathematicaInput);
    //        fvec.Print("fvec = ", sout,EMathematicaInput);
    //        //Print the temporal solution
    //        Initialsolution.Print("Intial conditions = ", sout,EMathematicaInput);
    //        TPZFMatrix<REAL> Temp;
    //        TPZFMatrix<REAL> Temp2;
    //        matM->Multiply(Initialsolution,Temp);
    //        Temp.Print("Temp matM = ", sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    
	int nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
	
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT;
	while (TimeValue <= maxTime)
	{
        timeatual  = TimeValue;
		// This time solution i for Transient Analytic Solution
		mymaterial->SetTimeValue(TimeValue);
		matM->Multiply(Lastsolution,TotalRhstemp);
        
        //#ifdef PZ_LOG
        //        if(logdata.isDebugEnabled())
        //        {
        //            std::stringstream sout;
        //            sout<< " tempo = " << cent;
        //            Lastsolution.Print("\nIntial conditions = ", sout,EMathematicaInput);
        //            TotalRhstemp.Print("Mat Mass x Last solution = ", sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logdata,sout.str())
        //        }
        //#endif
        
		TotalRhs = fvec + TotalRhstemp;
		an.Rhs() = TotalRhs;
		an.Solve();
		Lastsolution = an.Solution();
		
        if(cent%ntimestep==0){
            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            PosProcessMultphysics(meshvec,mphysics,an,plotfile);
        }
		
        cent++;
		TimeValue = cent*deltaT;
//        
//        if(cent == 100){
//            deltaT = 1000*deltaT;
//            mymaterial->SetTimeStep(deltaT);
//            mphysics->Solution().Zero();
//            matM = MassMatrix(mymaterial, mphysics);
//            StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fvec);
//            
//        }
	}
}


void DadosMalhas::SolveSistWithError(REAL deltaT,REAL maxTime,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, int ntimestep, REAL &timeatual,TPZAutoPointer<TPZFunction<STATE> > solExata1, TPZAutoPointer<TPZFunction<STATE> > solExata2, int h,  ofstream saidaerro){
    
    TPZMaterial *mat = mphysics->FindMaterial(fmatId);
    TPZPoroElasticMF2d *matporoelast = dynamic_cast<TPZPoroElasticMF2d *>(mat);
    matporoelast->SetTimeStep(deltaT);
    
    TPZAnalysis an(mphysics);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    std::string outputfile;
	outputfile = "TransientSolution";
    
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    PosProcessMultphysics(meshvec,mphysics,an,plotfile);
    
    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(matporoelast, mphysics,8);
    
    //Criando matriz de rigidez (matK) e vetor de carga
	TPZFMatrix<STATE> matK;
	TPZFMatrix<STATE> fvec;
    StiffMatrixLoadVec(matporoelast, mphysics, an, matK, fvec,8);
    
    int nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
	
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT;
	while (TimeValue <= maxTime)
	{
        timeatual  = TimeValue;
		// This time solution i for Transient Analytic Solution
		matporoelast->SetTimeValue(TimeValue);
		matM->Multiply(Lastsolution,TotalRhstemp);
        
		TotalRhs = fvec + TotalRhstemp;
		an.Rhs() = TotalRhs;
		an.Solve();
		Lastsolution = an.Solution();
		
        if(cent%ntimestep==0){
            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            PosProcessMultphysics(meshvec,mphysics,an,plotfile);
        }
		
        TPZVec<REAL> erros;
        
//        saidaerro<<" Erro da simulacao multifisica do deslocamento (u)" <<endl;
//        TPZAnalysis an12(meshvec[0]);
//        an12.SetExact(*solExata1);
//        an12.PostProcessError(erros, saidaerro);
//        
//        
//        saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
//        TPZAnalysis an22(meshvec[1]);
//        an22.SetExact(*solExata2);
//        an22.PostProcessError(erros, saidaerro);
//        
//        saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
//        TPZAnalysis an32(meshvec[2]);
//        an32.SetExact(*solExata2);
//        an32.PostProcessError(erros, saidaerro);

        
        cent++;
		TimeValue = cent*deltaT;
	}
}

void DadosMalhas::SolveSistBarryMercert(REAL deltaT,REAL maxTime,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, int ntimestep, REAL &timeatual){
	
    
    
    TPZMaterial *mat1 = mphysics->FindMaterial(fmatId);
    TPZMaterial *mat2 = mphysics->FindMaterial(fmatId+1);
    TPZPoroElasticMF2d * mat12 = dynamic_cast<TPZPoroElasticMF2d *>(mat1);
    TPZPoroElasticMF2d * mat22 = dynamic_cast<TPZPoroElasticMF2d *>(mat2);
    
    mat12->SetTimeStep(deltaT);
    mat22->SetTimeStep(deltaT);
    
    TPZAnalysis an(mphysics);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    std::string outputfile;
	outputfile = "TransientSolution";
    
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    PosProcessMultphysics(meshvec,mphysics,an,plotfile);
    
    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(mphysics,8);
    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //            std::stringstream sout;
    //        	matM->Print("matM = ", sout,EMathematicaInput);
    //        	LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    //Criando matriz de rigidez (matK) e vetor de carga
	TPZFMatrix<STATE> matK;
	TPZFMatrix<STATE> fvec;
    StiffMatrixLoadVec(mphysics, an, matK, fvec,8);
    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //
    //        std::stringstream sout;
    //        matK.Print("matK = ", sout,EMathematicaInput);
    //        fvec.Print("fvec = ", sout,EMathematicaInput);
    //        //Print the temporal solution
    //        Initialsolution.Print("Intial conditions = ", sout,EMathematicaInput);
    //        TPZFMatrix<REAL> Temp;
    //        TPZFMatrix<REAL> Temp2;
    //        matM->Multiply(Initialsolution,Temp);
    //        Temp.Print("Temp matM = ", sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    
	int nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
	
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT;
	while (TimeValue <= maxTime)
	{
        timeatual  = TimeValue;
		// This time solution i for Transient Analytic Solution
		matM->Multiply(Lastsolution,TotalRhstemp);
        StiffMatrixLoadVec(mphysics, an, matK, fvec,8);
        
        //#ifdef PZ_LOG
        //        if(logdata.isDebugEnabled())
        //        {
        //            std::stringstream sout;
        //            sout<< " tempo = " << cent;
        //            Lastsolution.Print("\nIntial conditions = ", sout,EMathematicaInput);
        //            TotalRhstemp.Print("Mat Mass x Last solution = ", sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logdata,sout.str())
        //        }
        //#endif
        
		TotalRhs = fvec + TotalRhstemp;
		an.Rhs() = TotalRhs;
		an.Solve();
		Lastsolution = an.Solution();
		
        if(cent%ntimestep==0){
            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            PosProcessMultphysics(meshvec,mphysics,an,plotfile);
        }
		
        cent++;
		TimeValue = cent*deltaT;
        
//        if(cent == 100){
//            deltaT = 1000*deltaT;
//            mat12->SetTimeStep(deltaT);
//            mat22->SetTimeStep(deltaT);
            //mphysics->Solution().Zero();
            //matM = MassMatrix(mphysics);
            //StiffMatrixLoadVec(mphysics, an, matK, fvec);
            
        //}
	}
}


#include "pzl2projection.h"
TPZCompMesh *DadosMalhas::CMeshPressureL2(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini, bool triang)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(fmatId, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
    ///inserir connect da pressao
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
	    newnod.SetLagrangeMultiplier(1);
    }
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        celdisc->SetConstC(1.);
//        celdisc->SetCenterPoint(0, 0.);
//        celdisc->SetCenterPoint(1, 0.);
//        celdisc->SetCenterPoint(2, 0.);
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(triang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }
    
    
#ifdef PZDEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
	return cmesh;
}

void DadosMalhas::RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
	
	int i;
	
	// Refinando no local desejado
	TPZVec<REAL> point(3);
    //	point[0] = point[1] = point[2] = -0.25;
    //	if(dim==1) point[1] = point[2] = 0.0;
    //	else if(dim==2) point[2] = 0.0;
    
    point[0] = 0.;
    point[1] = 1.;
    point[2] = 0.;
	REAL r = sqrt(1.);
	
	if(ntyperefs==2) {
		REAL radius = 0.2;
		for(i=0;i<nref;i+=2) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius);
			RefineGeoElements(dim,gmesh,point,r,radius);
			if(nref < 5) radius *= 0.35;
			else if(nref < 7) radius *= 0.2;
			else radius *= 0.1;
		}
		if(i==nref) {
			RefineGeoElements(dim,gmesh,point,r,radius);
		}
	}
	else {
		REAL radius = 0.5;
		for(i=0;i<nref+1;i++) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius);
			radius *= 0.2;
		}
	}
	// Constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}

void DadosMalhas::RefiningNearLine(int dim,TPZGeoMesh *gmesh,int nref) {
	
	int i;
	
	// Refinando no local desejado
	TPZManVector<REAL> point(3);
	point[0] = 0; point[1] =  5.0; point[2] = 0.0;
	REAL r = 0.0;
	
	REAL radius = 6.;
	for(i=0;i<nref;i++) {
		// To refine elements with center near to points than radius
		RefineGeoElements(dim,gmesh,point,r,radius);
		radius *= 0.6;
        //radius /= 0.8
        
        // Constructing connectivities
        gmesh->ResetConnectivities();
        gmesh->BuildConnectivity();
        AjustarContorno(gmesh);
	}
	
}


void DadosMalhas::RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
    
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		REAL centerdist = TPZGeoEl::Distance(center,point);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
}

void DadosMalhas::AjustarContorno(TPZGeoMesh *gmesh)
{
    gmesh->BuildConnectivity();
    int nel = gmesh->NElements();
    for(int ie = 0; ie<nel; ie++){
        TPZGeoEl *gel = gmesh->ElementVec()[ie];
        if(!gel || gel->HasSubElement()) continue;
        
        int matid = gel->MaterialId();
        if(matid>0) continue;
        
        int nside = gel->NSides();
        TPZGeoElSide thisside(gel,nside-1);
        TPZGeoElSide neighbour = thisside.Neighbour();
        
        int dimel = gel->Dimension();
        int dimn = neighbour.Element()->Dimension();
        if(dimn==dimel) continue;
        
        int nsubel = neighbour.Element()->NSubElements();
        
        for(int is=0; is < nsubel; is++){
            TPZGeoEl *neigel = neighbour.Element()->SubElement(is);
            if(neigel->HasSubElement()) continue;
            
            int mylevel  = gel->Level();
            int neiglevel = neigel->Level();
            
            if(mylevel!=neiglevel){
                TPZVec< TPZGeoEl * > filhos;
                gel->Divide(filhos);
                break;
            }
        }
    }
}
