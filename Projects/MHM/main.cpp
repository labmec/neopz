#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzmatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "TPZInterfaceEl.h"


#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzbndcond.h"


#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "tpzsparseblockdiagonalstructmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzbndmat.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"


#include "pzlog.h"

#include "pzskylstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzpoisson3d.h"
#include "DarcyMHM.h"
#include <time.h>
#include <stdio.h>
#include <sstream>
#include <fstream>

#include "TPZVTKGeoMesh.h"

// Using Log4cXX as logging tool
//
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.poisson3d"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.poisson3d.data"));
#endif
//
// End Using Log4cXX as logging tool

TPZGeoMesh * MalhaGeo ( const int h );
TPZCompMesh *CreateCompMeshHybrid ( TPZGeoMesh &gmesh, int porder );
void BuildByParts(TPZCompMesh & result);
long FindMidNode(TPZCompMesh &cmesh);
void Solve ( TPZAnalysis &an );

int main()
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
	int h = 1;
    
	TPZGeoMesh *gmesh = MalhaGeo(h);
    
    TPZCompMesh *cmesh =CreateCompMeshHybrid(*gmesh, 1);
    
    std::ofstream arcg2 ( "gmesh2.txt" );
    cmesh->Reference()->Print ( arcg2 );

    
    std::ofstream arc1 ( "cmesh.vtk" );
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh->Reference(), arc1, true);
    
	TPZAnalysis an(cmesh,true);
    
	Solve( an );
    {
        std::ofstream arc ( "cmesh.txt" );
        cmesh->Print ( arc );
    }
    
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("state");
    an.DefineGraphMesh(2, scalnames, vecnames, "MHM.vtk");
    an.PostProcess(0, 2);

//	PosProcess(an, AE, h, p);


	
	return 0;
}

TPZGeoMesh * MalhaGeo( const int h )
{
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	REAL co[9][2] = {{0.,0.},{0.5,0.},{0.5,0.5},{0.,0.5},{1.,0.},{1.,0.5},{1.,1.},{0.5,1.},{0.,1.}};
	int indices[4][4] = {{0,1,2,3},{1,4,5,2},{2,5,6,7},{3,2,7,8}};
	TPZGeoEl *elvec[4];
	int nnode = 9;
	int nelem = 4;
	int nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}
	
	int el;
	for ( el=0; el<nelem; el++ )
	{
		TPZVec<long> nodind ( 4 );
		for ( nod=0; nod<4; nod++ ) nodind[nod]=indices[el][nod];
		long index;
		elvec[el] = gmesh->CreateGeoElement ( EQuadrilateral,nodind,1,index );
	}
	
	gmesh->BuildConnectivity();
	
	
	TPZGeoElBC gbc1 ( elvec[0],4,-1 );
	TPZGeoElBC gbc2 ( elvec[1],4,-1 );
	TPZGeoElBC gbc3 ( elvec[1],5,-1 );
	TPZGeoElBC gbc4 ( elvec[2],5,-1 );
	TPZGeoElBC gbc5 ( elvec[2],6,-1 );
	TPZGeoElBC gbc6 ( elvec[3],6,-1 );
	TPZGeoElBC gbc7 ( elvec[3],7,-1 );
	TPZGeoElBC gbc8 ( elvec[0],7,-1 );//	TPZGeoElBC gbc8 ( elvec[1],7,-8,*gmesh );

	TPZGeoElBC gbc9  ( elvec[0],5,2 );
	TPZGeoElBC gbc10 ( elvec[1],6,2 );
	TPZGeoElBC gbc11 ( elvec[2],7,2 );
	TPZGeoElBC gbc12 ( elvec[3],4,2 );
    
    
    gmesh->ResetReference();
    
	for ( int ref = 0; ref < h; ref++ )
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
   		TPZVec<TPZGeoEl *> filhos1;
		int n = gmesh->NElements();
        for ( int i = 0; i < n; i++ )
		{
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
//            if (gel->Dimension() == 1 && gel->MaterialId()<0) gel->Divide ( filhos1 );
		}//for i
	}//ref


    gmesh->AddInterfaceMaterial(1, -1, 10);
    gmesh->AddInterfaceMaterial(1, -2, 11);
    gmesh->AddInterfaceMaterial(1, -3, 12);
    gmesh->AddInterfaceMaterial(1, -4, 13);
    gmesh->AddInterfaceMaterial(1, -5, 14);
    gmesh->AddInterfaceMaterial(1, -6, 15);
    gmesh->AddInterfaceMaterial(1, -7, 16);
    gmesh->AddInterfaceMaterial(1, -8, 17);
    gmesh->AddInterfaceMaterial(1, 2, 18);

    gmesh->AddInterfaceMaterial(-1, 1, 10);
    gmesh->AddInterfaceMaterial(-2, 1, 11);
    gmesh->AddInterfaceMaterial(-3, 1, 12);
    gmesh->AddInterfaceMaterial(-4, 1, 13);
    gmesh->AddInterfaceMaterial(-5, 1, 14);
    gmesh->AddInterfaceMaterial(-6, 1, 15);
    gmesh->AddInterfaceMaterial(-7, 1, 16);
    gmesh->AddInterfaceMaterial(-8, 1, 17);
    gmesh->AddInterfaceMaterial( 2, 1, 18);
    
    
    std::ofstream arcg ( "gmesh.txt" );
    gmesh->Print ( arcg );
    //
    
    std::ofstream vtkfile("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

	return gmesh;
}


TPZCompMesh *CreateCompMeshHybrid ( TPZGeoMesh &gmesh, int porder ){
	TPZCompEl::SetgOrder ( porder );
	TPZCompMesh *result = new TPZCompMesh( &gmesh );
	result->SetDimModel ( 2 );
	
	//result->SetAllCreateFunctionsDiscontinuous();
    result->SetAllCreateFunctionsContinuous();//default, não precisaria ser setado novamente
	
//	TPZMatPoisson3d *material ;
//	material = new TPZMatPoisson3d ( 1,2 );

    TPZMatDarcyMHM *materialinterface1 = new TPZMatDarcyMHM( 10,2 );
  	result->InsertMaterialObject ( materialinterface1 );
    TPZMatDarcyMHM *materialinterface2 = new TPZMatDarcyMHM( 11,2 );
    materialinterface2->SetMultiplier(-1.);
  	result->InsertMaterialObject ( materialinterface2 );
    

    
    TPZMatDarcyMHM *material ;
	material = new TPZMatDarcyMHM( 1,2 );
	   
//	TPZVec<REAL> convdir(3);//direção da convecção
//	convdir[0]=0.0;//sqrt(2.)/2.;
//	convdir[1]=0.0;//sqrt(2.)/2.;
//	REAL diff= 1.;
//	REAL conv=0.;
//	material-> SetParameters(diff, conv, convdir);
//    TPZAutoPointer<TPZFunction<STATE> > ExactSol = new TPZDummyFunction<STATE> (ExactSolution);
//	material->SetForcingFunctionExact(ExactSol);
    //material->SetNonSymmetric();
	//material->SetSymmetric();
	TPZMaterial *mat ( material );
    
    
	
//	TPZAutoPointer<TPZFunction<STATE> > BC1 = new TPZDummyFunction<STATE> (CC1);
//	TPZAutoPointer<TPZFunction<STATE> > BC2 = new TPZDummyFunction<STATE> (CC2);
//	TPZAutoPointer<TPZFunction<STATE> > BC3 = new TPZDummyFunction<STATE> (CC3);
//	TPZAutoPointer<TPZFunction<STATE> > BC4 = new TPZDummyFunction<STATE> (CC4);
    
//	TPZAutoPointer<TPZFunction<STATE> > LoadVector = new TPZDummyFunction<STATE> (Forcing1);
//	material->SetForcingFunction ( LoadVector );
	result->InsertMaterialObject ( mat );
	
	TPZFMatrix<STATE> val1 ( 1,1,0. ), val2 ( 1,1,0. );// 0 é Dirichlet, 1 é Neumann, 2 é Robin(implementada apenas no Contínuo)
    
	TPZMaterial *bnd1 = material->CreateBC ( mat,-1,1, val1, val2 );
//	TPZMaterial *bnd2 = material->CreateBC ( mat,-2,0, val1, val2 );
//	TPZMaterial *bnd3 = material->CreateBC ( mat,-3,0, val1, val2 );
//	TPZMaterial *bnd4 = material->CreateBC ( mat,-4,0, val1, val2 );
//	TPZMaterial *bnd5 = material->CreateBC ( mat,-5,0, val1, val2 );
//	TPZMaterial *bnd6 = material->CreateBC ( mat,-6,0, val1, val2 );
//	TPZMaterial *bnd7 = material->CreateBC ( mat,-7,0, val1, val2 );
//	TPZMaterial *bnd8 = material->CreateBC ( mat,-8,0, val1, val2 );
    
   	TPZMaterial *bnd9 = material->CreateBC ( mat,2,1, val1, val2 );
    
//	bnd1->SetForcingFunction ( BC1 );
//	bnd2->SetForcingFunction ( BC2 );
//	bnd3->SetForcingFunction ( BC3 );
//	bnd4->SetForcingFunction ( BC4 );
	
	result->InsertMaterialObject ( bnd1 );
//	result->InsertMaterialObject ( bnd2 );
//	result->InsertMaterialObject ( bnd3 );
//	result->InsertMaterialObject ( bnd4 );
//	result->InsertMaterialObject ( bnd5 );
//	result->InsertMaterialObject ( bnd6 );
//	result->InsertMaterialObject ( bnd7 );
//	result->InsertMaterialObject ( bnd8 );
	result->InsertMaterialObject ( bnd9 );
	
    BuildByParts(*result);

	result->SetName("CMesh1");
    


	return result;
}

void BuildByParts(TPZCompMesh & cmesh)
{

    TPZGeoMesh *gmesh=cmesh.Reference();
    for (int igel=0; igel<gmesh->NElements(); igel++)
    {
        TPZGeoEl *gelroot=gmesh->ElementVec()[igel];
        if (gelroot && gelroot->MaterialId()== 1 && !gelroot->Father()) //se eh um macroelemento
        {
        
        // gelroot é um "macroelemento"
            for (int igel2=0; igel2<gmesh->NElements(); igel2++)
            {
                TPZGeoEl *gel2=gmesh->ElementVec()[igel2];
                if (!gel2|| gel2->LowestFather()!=gelroot||gel2->HasSubElement())
                {
                    continue;
                }
                long index;
                cmesh.CreateCompEl(gel2, index);// cria os elementos 2D
            }
            long midconnectindex = FindMidNode(cmesh);
            TPZConnect &c = cmesh.ConnectVec()[midconnectindex];
            c.SetLagrangeMultiplier(2);
            gmesh->ResetReference();
            gmesh->SetReference(&cmesh);
        }

        if (gelroot && gelroot->Dimension()==1  && !gelroot->Father()) //se eh um elemento 1D (Lambda)
        {
            long index3;
            cmesh.CreateCompEl(gelroot, index3);
            TPZCompEl *cel = cmesh.ElementVec()[index3];
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.SetLagrangeMultiplier(1);
            }
            cel->Reference()->ResetReference();
        }
    }
    cmesh.LoadReferences();
    cmesh.ExpandSolution();
    
//    // criacao dos elementos de interface
//    for (int icel=0; icel<cmesh.NElements(); icel++)
//    {
//        TPZCompEl * cel = cmesh.ElementVec()[icel];
//        TPZGeoEl *gelRef = cel->Reference();
//        for (int iSide=gelRef->NCornerNodes(); iSide < gelRef->NSides()-1; iSide++)
//        {
//            TPZGeoElSide gelRefSide = gelRef->Neighbour(iSide);
//            if (gelRefSide.LowerLevelCompElementList2(true)  )
//            {
//                TPZInterpolationSpace *sp= dynamic_cast<TPZInterpolationSpace *>(cel);
//                sp->CreateInterface(iSide,true);
//            }
//        }
//        
//    }
    
    //4- Create the interface elements
    int meshdim = cmesh.Dimension();
	for (int igel3=0; igel3<gmesh->NElements(); ++igel3)
    {
		TPZGeoEl *gel = gmesh->ElementVec()[igel3];// percorro todos os elementos de dimensao 1 da gmesh
		if (!gel || gel->Dimension() != meshdim-1 || !gel->Reference())
        {
            continue;
		}
		int matid = gel->MaterialId();
        // escolho dentre esses elementos os que tem o material id que eu quero
		if(matid != 2 && matid > 0)//2 sao os elementos 1D do interior e < 0 sao os 1D no contorno
        {
			continue;
		}
		//over the dimension-1 sides
		int nsides = gel->NSides();
        if(nsides!=3)
        {
            DebugStop();// todos os elementos 1D tem 3 lados
        }
		int is;
		for (is=0; is<nsides; ++is)// percorro os lados do elemento, embora esteja interessado apenas no lado 2
        {
			int sidedim = gel->SideDimension(is);
			if (sidedim != meshdim-1)// apenas o lado 2 tem dimensao 1
            {
				continue;
			}
			TPZStack<TPZCompElSide> celsides;
			TPZGeoElSide gelside(gel,is);
            TPZGeoElSide gelneigh(gelside.Neighbour());
//            gelneigh.Element()->Print(std::cout);
            int reorient = 0;
            while (gelneigh != gelside)
            {
                TPZStack<TPZGeoElSide> subgeoelements;
                gelneigh.GetSubElements2(subgeoelements);
                if (subgeoelements.NElements() == 0 && gelneigh.Element()->MaterialId() == 1) {
                    subgeoelements.Push(gelneigh);
                }
                for (int i = 0; i < subgeoelements.NElements(); i++)
                {
                    if (reorient > 1) {
                        DebugStop();
                    }
                    if (subgeoelements[i].Dimension() == gelside.Dimension())
                    {
                        TPZGeoElSide geosideneigh = subgeoelements[i];
//                      geosideneigh.Element()->Print();
                        TPZCompElSide right=gelside.Reference();
                        TPZCompElSide left =geosideneigh.Reference();
                        int matid = 10+reorient;
                        TPZGeoEl *interfaceEl = geosideneigh.Element()->CreateBCGeoEl(geosideneigh.Side(), matid);
                        long index;
                        TPZInterfaceElement *newel = new TPZInterfaceElement(cmesh,interfaceEl,index,left,right);
                    }
                    
                }
                if (subgeoelements.NElements() != 0) {
                    reorient++;
                }
                gelneigh = gelneigh.Neighbour();
            }
		}
	}
}


void Solve ( TPZAnalysis &an )
{
	TPZCompMesh *malha = an.Mesh();
    
	//TPZBandStructMatrix mat(malha);
	TPZSkylineStructMatrix mat(malha);// requer decomposição simétrica, não pode ser LU!
	//TPZBlockDiagonalStructMatrix mat(malha);//ok!
	//TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
	//TPZFStructMatrix mat( malha );// ok! matriz estrutural cheia
	//TPZSpStructMatrix mat( malha );//matriz estrutural esparsa (???? NÃO FUNCIONOU !!!!!!!!!!)
    mat.SetNumThreads(0);
	TPZStepSolver<STATE> solv;
	solv.SetDirect (  ELDLt );//ECholesky);// ELU , ELDLt ,
    
	
    //	cout << "ELDLt " << endl;
	an.SetSolver ( solv );
	an.SetStructuralMatrix ( mat );
    std::cout << std::endl;
	an.Solution().Redim ( 0,0 );
    std::cout << "Assemble " << std::endl;
	an.Assemble();
    {
        std::ofstream fileout("rigidez.nb");
        an.Solver().Matrix()->Print("Rigidez = ", fileout, EMathematicaInput);
        an.Rhs().Print("Rhs = ", fileout, EMathematicaInput);
    }
	an.Solve();
    std::cout << std::endl;
    std::cout << "No equacoes = " << malha->NEquations() << std::endl;
}

long FindMidNode(TPZCompMesh &cmesh)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    int nel = gmesh->NElements();
    TPZManVector<REAL,3> xmid(3);
    int ncont = 0;
    for (int el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        if(!gel->Reference()) continue;
        int side = gel->NSides()-1;
        TPZManVector<REAL,3> xi(gel->Dimension()),xco(3);
        gel->CenterPoint(side, xi);
        gel->X(xi, xco);
        for(int i=0; i<3; i++) xmid[i] += xco[i];
        ncont++;
    }
    if (ncont == 0) {
        DebugStop();
    }
    for(int i=0; i<3; i++) xmid[i] /= ncont;
    TPZManVector<REAL,2> xi(2);
    long zero = 0;
    TPZGeoEl *centerel = gmesh->FindElement(xmid, xi, zero , 2);
    if (!centerel || ! centerel->Reference()) {
        DebugStop();
    }
    return centerel->Reference()->ConnectIndex(0);
}
