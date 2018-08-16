/**
 * @file
 * @brief Project to save and load information of meshes into the disk.
 */
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzcheckgeom.h"

#include "pzmatrix.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZMaterial.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"

#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "TPZRefPatternTools.h"

#include "../Poisson_ArcTan/pzgclonemesh.h"
#include "../Poisson_ArcTan/pzcclonemesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

#include "pzcheckmesh.h"

#include "pzlog.h"
#include "TPZPersistenceManager.h"

int ExtractingCommandRegistered(std::ifstream &file,std::string &cmeshname,TPZStack<std::string> &commands);
void ApplyCommand(TPZCompMesh *cmesh,TPZVec<std::string> &command);

/** Ideia principal: Recuperar uma malha computacional do disco, realizar mudancas registradas em um arquivo (lendo apenas as linhas com o referencia ao nome da malha, e a malha final deve ser a mesma obtida no processo inicial.
 */
bool TestingLoadingSavedMeshes();
int gDebug = 0;

void MakeCompatibles(TPZGeoMesh *gmesh,TPZCompMesh *cmesh);

std::ofstream *out;


/** Ideia principal: Create a computational meshes on the order on side is different on Max Order preferred on this side
 */
bool TestingOrderIncompatibilityOnRestrainedSides();
TPZGeoMesh *CreateQuadrilateralMesh();
TPZGeoMesh *CreateQuadrilateralMesh2();
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);



int main() {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    if(TestingOrderIncompatibilityOnRestrainedSides())
        return 10;
    
    if(!TestingLoadingSavedMeshes())
        return 1;

    std::cout << std::endl << "\tThe END." << std::endl;
    return 0;
}

bool TestingOrderIncompatibilityOnRestrainedSides() {
    
    out = new (std::ofstream)("saida.txt");
	// Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);

    TPZGeoMesh *gmesh;
    gmesh = CreateQuadrilateralMesh();
//    gmesh = CreateQuadrilateralMesh2();
    if(!gmesh) return false;
//    gmesh->Print();
    
    int p = 1;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMesh(gmesh);
    if(!cmesh) return false;

    // Dividing overlapped elements
    std::cout << "\nDividing element 0 and 2. THESE ELEMENTS ARE OVERLAPPED!!\n";
	TPZVec<int64_t> subels;
    cmesh->ElementVec()[0]->Divide(cmesh->ElementVec()[0]->Index(),subels);
//    gmesh->Print();
    cmesh->ElementVec()[2]->Divide(cmesh->ElementVec()[2]->Index(),subels);
//    gmesh->Print();
    // Dividing upper left subelement of the element 2.
    cmesh->ElementVec()[10]->Divide(cmesh->ElementVec()[10]->Index(),subels);
    cmesh->ExpandSolution();
//    gmesh->Print();
//    cmesh->Print();

    // Setting order p=2 for element 1
    std::cout << "\nPRefine with order 2 for element 1.\n";
    ((TPZInterpolatedElement *)(cmesh->ElementVec()[1]))->PRefine(2);
    cmesh->ExpandSolution();
    gmesh->Print();
    cmesh->Print();

    // Dividing left sub elements of the element 10 (upper left sub element of the element 2)
 //   cmesh->ElementVec()[11]->Divide(cmesh->ElementVec()[11]->Index(),subels);
 //   cmesh->ElementVec()[14]->Divide(cmesh->ElementVec()[14]->Index(),subels);
 //   gmesh->Print();

    
    std::stringstream sout("incompatible.txt");
    TPZCheckMesh checker(cmesh,&sout);
    
    // If no incompatilibities found return true (1).
    if(!checker.VerifyAllConnects())
        sout << " Exists connect with incompatibility.\n";
    // If no order incompatibilities was found return -1.
    if(checker.CheckConnectOrderConsistency()>-1)
        sout << " Exists order inconsistency.\n";
    
    std::cout << "\n";
    return true;
}
TPZGeoMesh *CreateQuadrilateralMesh() {
    REAL co[6][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.},
        {-1.,1.,0.},
        {-1.,0.,0.},
    };
    int64_t indices[3][4] = {{0,1,2,3},{0,3,4,5},{0,1,2,3}};
    
    const int nelem = 3;
    int nnode = 6;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int64_t nod;
    for(nod=0; nod<nnode; nod++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int64_t el;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}
TPZGeoMesh *CreateQuadrilateralMesh2() {
    REAL co[8][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.},
        {-1.,1.,0.},
        {-1.,0.,0.},
        {2.,0.,0.},
        {2,1,0}
    };
    int64_t indices[3][4] = {{0,1,2,3},{0,3,4,5},{0,6,7,3}};
    
    const int nelem = 3;
    int nnode = 8;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int64_t nod;
    for(nod=0; nod<nnode; nod++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int64_t el;
    for(el=0; el<nelem; el++) {
        TPZManVector<int64_t> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;

}
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
	// Creating Poisson material
    int dim = 2;
    int MaterialId = 1;
	TPZMaterial *mat = new TPZMatPoisson3d(MaterialId,dim);
	TPZVec<REAL> convd(3,0.);
	((TPZMatPoisson3d *)mat)->SetParameters(1.,0.,convd);
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	cmesh->SetAllCreateFunctionsContinuous();
    
	// Boundary conditions
	// Dirichlet
	TPZFMatrix<STATE> val1(dim,dim,0.),val2(dim,1,0.);
	TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
	cmesh->InsertMaterialObject(bnd);
	
	cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
	cmesh->CleanUpUnconnectedNodes();
	return cmesh;
}

bool TestingLoadingSavedMeshes() {
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();
    // gRefDBase.InitializeRefPatterns();
    
    std::string filename, cmeshname;
    std::cout << std::endl << "INPUT - Name of file to load mesh ";
    std::cin >> filename;
	//char *q = strchr(((char *)cmeshname.c_str()),'.');
    //	if(!q)
	//    filename = cmeshname + ".txt";
	// Verifying if file exist and it is open
	std::ifstream in(filename.c_str());
	if(!in.is_open())
		return false;
	in.close();
    
    TPZPersistenceManager::OpenRead(filename);
	for(int i=0;i<filename.size();i++) {
		char p = filename[i];
		if(p=='_') break;
		cmeshname += ((char)toupper(p));
	}
    // Files with information meshes
    std::ofstream outfirstmesh("InitialCMesh.txt");
    std::ofstream outmesh("FinalCMesh.txt");
    
    // Creating geometric mesh
	TPZGeoMesh* gmesh;
    gmesh = dynamic_cast<TPZGeoMesh* >(TPZPersistenceManager::ReadFromFile());
    //    gmesh.Read(fstr,0);
    TPZCompMesh* cmesh;
    cmesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::ReadFromFile());
	MakeCompatibles(gmesh,cmesh);
    //    cmesh.Read(fstr,gmesh);
    //    cmesh->AutoBuild();
    int dim = cmesh->Dimension();
    cmesh->Print(outfirstmesh);
    outfirstmesh.close();
    
    /** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("Solution");
    
    // Resolver
    // Solve using symmetric matrix then using Cholesky (direct method)
    TPZAnalysis an(cmesh);
    TPZSkylineStructMatrix strskyl(cmesh);
    an.SetStructuralMatrix(strskyl);
    {
        std::stringstream sout;
        sout << "LoadMesh_" << dim << "D.vtk";
        //      an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
    }
    
    // Post processing
    an.PostProcess(0,dim);
    
    // Cleaning LOG directory
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Reading file "pzinterpolateddivide.txt"
    TPZStack<std::string> commands;
    std::ifstream log("pzinterpolatedrefine.txt");
    if(ExtractingCommandRegistered(log,cmeshname,commands)) {
        ApplyCommand(cmesh,commands);
        cmesh->AutoBuild();
        cmesh->ExpandSolution();
    }
    
    cmesh->Print(outmesh);
    outmesh.close();
    TPZAnalysis analysis(cmesh);
    {
        std::stringstream sout;
        sout << "LoadAndRefinedMesh_" << dim << "D.vtk";
        analysis.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
    }
    
    // Post processing
    analysis.PostProcess(0,dim);

    return true;
}

// Save information of the current mesh in disk
void SaveCompMesh(TPZCompMesh *cmesh, int timessave,TPZCompMesh *cmeshmodified,bool check) {
    if(!cmesh || timessave < 0) {
        std::cout << "SaveCompMesh - Bad argument: " << (void *)cmesh << " " << timessave << std::endl;
        return;
    }
#ifdef LOG4CXX
    {
        std::stringstream soutthis;
        if(cmeshmodified) soutthis << (void*)cmeshmodified;
        else soutthis << (void*)cmesh;
        // Rename the computational mesh
        cmesh->SetName(soutthis.str());
        soutthis << "_" << timessave;
        std::string filenamethis("LOG/");
        filenamethis.append(soutthis.str());
        filenamethis.append(".txt");
        TPZPersistenceManager::OpenWrite(filenamethis);
        
        // Renaming the geometric mesh
        std::stringstream gout;
        gout << (void*)cmesh->Reference();
        cmesh->Reference()->SetName(gout.str());
        
        // Save geometric mesh data
        TPZPersistenceManager::WriteToFile(cmesh->Reference());
        // Save computational mesh data
        TPZPersistenceManager::WriteToFile(cmesh);
        TPZPersistenceManager::CloseWrite();
        // To check printing computational mesh data in file
        if(check) {
            std::string filename("Mesh_");
            filename.append(soutthis.str());
            filename.append(".txt");
            std::ofstream arq(filename.c_str());
            cmesh->Print(arq);
        }
    }
#endif
}

void MakeCompatibles(TPZGeoMesh *gmesh,TPZCompMesh *cmesh) {
	int64_t ig, ngels = gmesh->NElements();
	int64_t ic, ncels = cmesh->NElements();
	TPZGeoEl *gel;
	TPZCompEl *cel;
	for(ig=0;ig<ngels;ig++) {
		gel = gmesh->ElementVec()[ig];
		if(!gel) continue;
		gel->ResetReference();
	}
	for(ic=0;ic<ncels;ic++) {
		cel = cmesh->ElementVec()[ic];
		if(!cel) continue;
		int64_t index = cel->GetRefElPatch()->Index();
		gel = gmesh->ElementVec()[index];
		if(!gel)
			DebugStop();
		gel->SetReference(cel);
	}
}
/* read a log line to extract command after a identifier (pointer of computational object)
int GetCommand(std::string &command,int nargs,TPZManVector<int,5> &argindex) {
    std::string commandknowed[32];
    commandknowed[0] = "Divide";
    commandknowed[1] = "PRefine";
    std::stringstream commandline(command);
    std::string commandname;
    commandline >> commandname >> nargs;
    argindex = new 
    
    for(int i=0;i<2;i++) {
        if(command.compare(commandknowed[i].c_str()))
            return i;
    }
    return -1;
}*/
void ApplyCommand(TPZCompMesh *cmesh,TPZVec<std::string> &commands) {
    int i;
    int64_t index, indexcel, indexgel;
    std::string commandname;
    TPZGeoMesh *gmesh = cmesh->Reference();
	gmesh->ResetReference();
	gmesh->SetReference(cmesh);
	cmesh->SetReference(gmesh);
	TPZGeoEl *gel;
	TPZCompEl *cel;

    // Making all divide first
    for(i=0;i<commands.NElements();i++) {
        std::stringstream commandline(commands[i]);
        commandline >> commandname;
        if(!commandname.compare("Divide")) {
            TPZVec<int64_t> subs;
            commandline >> index;
            commandline >> indexgel;
			cel = cmesh->ElementVec()[index];
			gel = gmesh->ElementVec()[indexgel];
			if(!cel) {
				std::cout << "\nComputational element was not exist. Exiting! \n";
			}
            indexcel = cel->Index();
			if(index != indexcel)
				std::cout << "\nDifferent index obtained from geometric element.\n";
            cmesh->Divide(index, subs, 1);
			// Verifying the indexes of new sub elements are the same as original mesh
			for(int j=0;j<subs.NElements();j++) {
				commandline >> indexgel;
				index = cmesh->ElementVec()[subs[j]]->Reference()->Index();
				if(indexgel != index) {
					cel = gmesh->ElementVec()[indexgel]->Reference();
					if(!cel) 
						cel = cmesh->ElementVec()[subs[j]];
				}
			}
        }
//    }
    // Making all refines
//    for(i=0;i<commands.NElements();i++) {
//        std::stringstream commandline(commands[i]);
//        commandline >> commandname;
        else if(!commandname.compare("PRefine")) {
            int order;
            commandline >> index;
			cel = cmesh->ElementVec()[index];
            commandline >> index;
            commandline >> order;
            gel = gmesh->ElementVec()[index];
			if(!gel)
				DebugStop();
			if(!gel->Reference()) {
				if(cel->Reference()->Index() == index)
					gel->SetReference(cel);
				else
					DebugStop();
			}
			indexcel = cel->Index();
            TPZInterpolatedElement *cint = dynamic_cast <TPZInterpolatedElement *> (cmesh->ElementVec()[indexcel]);
            cint->PRefine(order);
        }
    }
}
// Read tpzinterpolateddivide.txt and apply all the refinements on cmesh
// Returns the number of commands extracted
int ExtractingCommandRegistered(std::ifstream &file,std::string &cmeshname,TPZStack<std::string> &commands) {
    std::string meshname;
    std::string commandline;

    while(!file.eof()) {
        std::getline(file, commandline);
        std::stringstream strcommand(commandline);
        strcommand >> meshname;
        if(!meshname.compare(cmeshname.c_str())) {
            char command[256];
            char *p;
            strcommand.getline(command,256);
            p = command;
            while(*p == ' ') p++;
            if(*p) commands.Push(p);
        }
    }
    return commands.NElements();
}
