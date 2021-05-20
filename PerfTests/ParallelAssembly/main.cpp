#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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
#include "TPZLinearAnalysis.h"

#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPZBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "mixedpoisson.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacian.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "pzintel.h"
#include <iostream>
#include <string>

#include <math.h>
#include <set>



using namespace std;


/// malha geometrica de grande porte
TPZAutoPointer<TPZGeoMesh> MalhaGeom(int dimension, TPZVec<int> &nblocks, int nref);


void RefinamentoUniforme(TPZAutoPointer<TPZGeoMesh> gmesh, int nref,TPZVec<int> dims);

TPZAutoPointer<TPZCompMesh> CreatePressureMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder);

#ifdef PZ_LOG
static TPZLogger logger("pz.mainskeleton");
#endif

enum MRegular {EUniform,EUnbalanced};

TPZAutoPointer<TPZCompMesh> CreateCompMesh(int dimension, int porder, int64_t nelem, MRegular regular);

int main(int argc, char *argv[])
{
    
    int porder = 1;
    int dimension = 2;
    int nelem = 1000000;
    TPZAutoPointer<TPZCompMesh> cmesh = CreateCompMesh(dimension, porder, nelem, EUniform);
    

    std::cout << "Computational mesh created\n";
    std::cout.flush();
//#ifdef PZDEBUG
//    {
//        std::ofstream out("../Pressure.txt");
//        cmesh->Print(out);
//    }
//#endif
    
    std::cout << "Number of equations " << cmesh->NEquations() << std::endl;
    std::cout << "Number of elements " << cmesh->NElements() << std::endl;
    
    //calculo solution
    TPZLinearAnalysis *an = new TPZLinearAnalysis(cmesh,false);

    TPZSkylineStructMatrix strmat(cmesh);
    
//#ifndef PZDEBUG
//    skyl.SetNumThreads(1);
//#endif
    an->SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    std::cout << "Assembly\n";
    an->AssembleResidual();
    std::cout << "Once\n";
    an->AssembleResidual();
    std::cout << "Finished\n";
    return 0;
}
    


void RefinamentoUniforme(TPZAutoPointer<TPZGeoMesh> gmesh, int nref,TPZVec<int> dims)
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



/// malha geometrica de grande porte
TPZAutoPointer<TPZGeoMesh> MalhaGeom(int dim, TPZVec<int> &nblocks, int nref)
{
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    REAL scale = pow(2, nref);
    REAL lx = nblocks[0]*scale;
    REAL ly = nblocks[1]*scale;
    REAL lz = scale;
    x1[0] = lx;
    x1[1] = ly;
    x1[2] = 0.;
    TPZManVector<int,3> nx(nblocks);
    TPZGenGrid2D gengrid(nx,x0,x1);
    TPZAutoPointer<TPZGeoMesh> meshresult2d = new TPZGeoMesh;
    gengrid.SetRefpatternElements(false);
    gengrid.Read(meshresult2d);
    
    gengrid.SetBC(meshresult2d.operator->(), 4, -1);
    gengrid.SetBC(meshresult2d.operator->(), 5, -1);
    gengrid.SetBC(meshresult2d.operator->(), 6, -1);
    gengrid.SetBC(meshresult2d.operator->(), 7, -1);
    if (dim ==2) {
        return meshresult2d;
    }
    if (dim != 3) {
        DebugStop();
    }
    TPZExtendGridDimension extend(meshresult2d,lz);
    // create uniform refinement elements
    extend.SetElType(0);
    TPZGeoMesh *res3d = extend.ExtendedMesh(nblocks[2],-2,-2);
    TPZAutoPointer<TPZGeoMesh> meshresult3d(res3d);

    TPZCheckGeom check(res3d);
    check.UniformRefine(nref);
    
    
#ifdef PZDEBUG
    std::ofstream vtkfile("../gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(res3d, vtkfile);
#endif
    return meshresult3d;
}



TPZAutoPointer<TPZCompMesh> CreatePressureMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder)
{
    TPZAutoPointer<TPZCompMesh> cmeshPressure = new TPZCompMesh(gmesh);
    cmeshPressure->SetDimModel(gmesh->Dimension());
    cmeshPressure->ApproxSpace().SetAllCreateFunctionsContinuous();
//    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshPressure->SetDefaultOrder(porder);
    TPZMatLaplacian *matl2 = new TPZMatLaplacian(1);
    matl2->SetDimension(cmeshPressure->Dimension());
    matl2->SetSymmetric();
    cmeshPressure->InsertMaterialObject(matl2);
    cmeshPressure->AutoBuild();
    std::cout << "Computational mesh created\n";
    return cmeshPressure;
}

TPZAutoPointer<TPZCompMesh> CreateCompMesh(int dimension, int porder, int64_t nelem, MRegular regular)
{
    TPZAutoPointer<TPZCompMesh> result;
    int pincr = 4;
    if (dimension == 3) {
        int blocksize = ((int)pow(nelem, 0.33))+1;
        TPZAutoPointer<TPZGeoMesh> gmesh;
        int nref = 0;
        TPZManVector<int> nblocks(3,blocksize);
        gmesh = MalhaGeom(dimension, nblocks, nref);
        gmesh->SetDimension(3);
        std::cout << "Geometric mesh created elements " << gmesh->NElements() << std::endl;
        result = CreatePressureMesh(gmesh, porder);
    }
    else if (dimension == 2)
    {
        pincr = 6;
        int blocksize = ((int)sqrt(nelem))+1;
        TPZAutoPointer<TPZGeoMesh> gmesh;
        int nref = 0;
        TPZManVector<int> nblocks(2,blocksize);
        gmesh = MalhaGeom(dimension, nblocks, nref);
        gmesh->SetDimension(2);
        std::cout << "Geometric mesh created elements " << gmesh->NElements() << std::endl;
        result = CreatePressureMesh(gmesh, porder);
    }
    else
    {
        DebugStop();
    }
    if (regular == EUnbalanced) {
        for (int64_t el =0; el < result->NElements(); el += result->NElements()/10) {
            TPZCompEl *cel = result->Element(el);
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if (!intel) {
                DebugStop();
            }
            intel->PRefine(porder+pincr);
        }
    }
    result->ExpandSolution();
    return result;
}


