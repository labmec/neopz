
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
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

#include "pzmaterial.h"
#include "pzbndcond.h"

#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

#include "pzlog.h"

void ExtractingCommandRegistered(std::ifstream &file,std::string &cmeshname,TPZStack<std::string> &commands);
void ApplyCommand(TPZCompMesh *cmesh,TPZVec<std::string> &command);

/** Ideia principal: Recuperar uma malha computacional do disco, realizar mudancas registradas em um arquivo (lendo apenas as linhas com o referencia ao nome da malha, e a malha final deve ser a mesma obtida no processo inicial.
 */
int main() {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    TPZFileStream fstr;
    std::string filename, cmeshname;
    std::cout << "Name of file to load mesh ";
    std::cin >> cmeshname;
    filename = cmeshname + ".txt";
    cmeshname.insert(11,1,0);
    fstr.OpenRead(filename);
	
    // Creating geometric mesh
	TPZGeoMesh gmesh;
    gmesh.Read(fstr,0);
    TPZCompMesh cmesh;
    cmesh.Read(fstr,&gmesh);
    int dim = cmesh.Dimension();
    std::ofstream outmesh("FinalCMesh.txt");
    
    // Reading file "tpzinterpolateddivide.log"
    TPZStack<std::string> commands;
    std::ifstream log("tpzinterpolateddivide.txt");
    ExtractingCommandRegistered(log,cmeshname,commands);
    cmesh.Print(outmesh);
    
    /** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("Solution");

    // Resolver
    // Solve using symmetric matrix then using Cholesky (direct method)
    TPZAnalysis an(&cmesh);
    TPZSkylineStructMatrix strskyl(&cmesh);
    an.SetStructuralMatrix(strskyl);
    {
        std::stringstream sout;
        sout << "Poisson" << dim << "D_MESH.vtk";
        an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
    }

    // Post processing
    an.PostProcess(0,dim);

    return 0;
}

// read a log line to extract command after a identifier (pointer of computational object)
int GetCommand(std::string &command,int &argindex) {
    std::string commandknowed[32];
    commandknowed[0] = "Divide";
    commandknowed[1] = "PRefine";
    
    for(int i=0;i<2;i++) {
        if(command.compare(commandknowed[i])>0)
            return i;
    }
    return -1;
}
void ApplyCommand(TPZCompMesh *cmesh,TPZVec<std::string> &command) {
    int i, index;
    TPZVec<int> subs;
    for(i=0;i<command.NElements();i++) {
        int icommand = GetCommand(command[i],index);
        switch(icommand) {
            case 0:
                cmesh->Divide(index,subs);
            default:
                break;
        }
    }
}
// Read tpzinterpolateddivide.log and apply all the refinements on cmesh
void ExtractingCommandRegistered(std::ifstream &file,std::string &cmeshname,TPZStack<std::string> &commands) {
    std::string meshname;
    std::string commandline;

    while(!file.eof()) {
        std::getline(file, commandline);
        std::stringstream strcommand(commandline);
        strcommand >> meshname;
        bool c = meshname==cmeshname;
        if(!c) {
            char command[256];
            strcommand.getline(command,256);
            commands.Push(command);
        }
    }
}
