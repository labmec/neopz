/**
 * @file
 * @brief Project to save and load information of meshes into the disk.
 */
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzmatrix.h"
#include "pzsave.h"

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

#include "pzmaterial.h"
#include "pzbndcond.h"

#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "TPZRefPatternTools.h"

#include "../Poisson_ArcTan/pzgclonemesh.h"
#include "../Poisson_ArcTan/pzcclonemesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

#include "pzlog.h"

int ExtractingCommandRegistered(std::ifstream &file,std::string &cmeshname,TPZStack<std::string> &commands);
void ApplyCommand(TPZCompMesh *cmesh,TPZVec<std::string> &command);


int gDebug = 0;

void MakeCompatibles(TPZGeoMesh *gmesh,TPZCompMesh *cmesh);

/** Ideia principal: Recuperar uma malha computacional do disco, realizar mudancas registradas em um arquivo (lendo apenas as linhas com o referencia ao nome da malha, e a malha final deve ser a mesma obtida no processo inicial.
 */
int main() {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();
   // gRefDBase.InitializeRefPatterns();

    TPZFileStream fstr;
    std::string filename, cmeshname;
    std::cout << std::endl << "INPUT - Name of file to load mesh ";
    std::cin >> filename;
	//char *q = strchr(((char *)cmeshname.c_str()),'.');
//	if(!q)
	//    filename = cmeshname + ".txt";
	// Verifying if file exist and it is open
	std::ifstream in(filename);
	if(!in.is_open())
		return 1;
	in.close();

    fstr.OpenRead(filename);
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
    gmesh = dynamic_cast<TPZGeoMesh* >(TPZSaveable::Restore(fstr,0));
//    gmesh.Read(fstr,0);
    TPZCompMesh* cmesh;
    cmesh = dynamic_cast<TPZCompMesh *>(TPZSaveable::Restore(fstr,gmesh));
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


    std::cout << std::endl << "\tThe END." << std::endl;
    return 0;
}

// Save information of the current mesh in disk
void SaveCompMesh(TPZCompMesh *cmesh, int timessave,TPZCompMesh *cmeshmodified,bool check) {
    if(!cmesh || timessave < 0) {
        std::cout << "SaveCompMesh - Bad argument: " << (void *)cmesh << " " << timessave << std::endl;
        return;
    }
#ifdef LOG4CXX
    {
        TPZFileStream fstrthis;
        std::stringstream soutthis;
        if(cmeshmodified) soutthis << (void*)cmeshmodified;
        else soutthis << (void*)cmesh;
        // Rename the computational mesh
        cmesh->SetName(soutthis.str());
        soutthis << "_" << timessave;
        std::string filenamethis("LOG/");
        filenamethis.append(soutthis.str());
        filenamethis.append(".txt");
        fstrthis.OpenWrite(filenamethis);
        
        // Renaming the geometric mesh
        std::stringstream gout;
        gout << (void*)cmesh->Reference();
        cmesh->Reference()->SetName(gout.str());
        
        // Save geometric mesh data
        int classid = cmesh->Reference()->ClassId();
        fstrthis.Write(&classid,1);   // this first data is necessary to use TPZSaveable::Restore
        cmesh->Reference()->Write(fstrthis,0);
        // Save computational mesh data
        classid = cmesh->ClassId();
        fstrthis.Write(&classid,1);   // this first data is necessary to use TPZSaveable::Restore
        cmesh->Write(fstrthis,0);
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
	long ig, ngels = gmesh->NElements();
	long ic, ncels = cmesh->NElements();
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
		long index = cel->GetRefElPatch()->Index();
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
    long index, indexcel, indexgel;
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
            TPZVec<long> subs;
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
