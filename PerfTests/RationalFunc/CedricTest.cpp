#include "pzshapelinear.h"

#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

#include "pzbstrmatrix.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "pzcheckmesh.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzpoisson3d.h"

#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"

#include "CedricTest.h"

#include "TPZRefPatternTools.h"

#include "TPZParFrontStructMatrix.h"

#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "arglib.h"

#include <stdio.h>

#include "pzlog.h"

#include "pzgeoelbc.h"

/** Initialiazing file for Log4CXX for this project */
#ifdef PZ_LOG
static TPZLogger logger("pz.Cedric");
#endif

#define DEFORM

TPZManVector<REAL,3> TCedricTest::fX0(3,0.5), TCedricTest::fEps(3,0.05);

TCedricTest::TCedricTest()
{
    REAL coord[8][3] = {
//        {0,0,0},
//        {1,0,0},
//        {1,1,0},
//        {0,1,0},
//        {0,0,1},
//        {1,0,1},
//        {1,1,1},
//        {0,1,1}
        {-0.5,-0.5,-0.5},
        {2,-1,-1},
        {1.1,1.1,-0.1},
        {-1,2,-1},
        {-1,-1,2},
        {1.2,-0.2,1.2},
        {2,2,2},
        {-0.5,1.5,1.5}
    };
    fDeformed.NodeVec().Resize(8);
    TPZManVector<int64_t,8> indices(8);
    for (int i=0; i<8; i++) {
        indices[i] = i;
        for (int c=0; c<3; c++) {
            fDeformed.NodeVec()[i].SetCoord(c, coord[i][c]);
        }
    }
    int64_t index;
    fDeformed.CreateGeoElement(ECube, indices, 1, index);
}

/// Deform the geometric mesh according to the coordinates of fDeformed
void TCedricTest::DeformGMesh(TPZGeoMesh &gmesh)
{
    int64_t nnodes = gmesh.NodeVec().NElements();
    TPZManVector<REAL,3> xbefore(3),xafter(3);
    for (int64_t nod=0; nod<nnodes; nod++) {
        gmesh.NodeVec()[nod].GetCoordinates(xbefore);
        for (int i=0; i<3; i++) {
            xbefore[i] = 2.*xbefore[i]-1.;
        }
        fDeformed.ElementVec()[0]->X(xbefore, xafter);
        gmesh.NodeVec()[nod].SetCoord(xafter);
    }
}
void TCedricTest::InterpolationError(int nsubdivisions,int geocase, int MaterialId,std::ostream &out)
{
    
    TPZGeoMesh *gmesh;
    switch(geocase) {
        case 1:
            gmesh = HexahedralMesh(nsubdivisions,MaterialId);
            break;
        case 2:
            gmesh = PyramidalAndTetrahedralMesh(nsubdivisions,MaterialId);
            break;
        case 3:
            gmesh = TetrahedralMesh(nsubdivisions,MaterialId);
            break;
        case 4:
            gmesh = TetrahedralMeshUsingRefinement(nsubdivisions,MaterialId);
            break;
    }

#ifdef DEFORM
    DeformGMesh(*gmesh);
    out << "Deformed ";
#else
    CheckConsistency(gmesh);
    out << "Regular ";
#endif
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int nelembc = AddBoundaryElements(gmesh);

    TPZCompEl::SetgOrder(1);
    
    /** Generating computational mesh */
    
    TPZCompMesh *cmesh = GenerateCompMesh(gmesh);
    if(!cmesh) {
        if(gmesh) {
            delete gmesh;
            gmesh = NULL;
        }
        out << "Interp_err nsubdivision " << nsubdivisions << " eltype " << geocase << " aborted\n";;
        return;
    }

    //int dim = cmesh->Dimension();

    TPZLinearAnalysis analysis(cmesh);

    out << "Interp_err nsubdivision " << nsubdivisions << " nelem " << (cmesh->NElements()-nelembc) << " eltype ";
    
    LoadInterpolation(cmesh);
    
    switch(geocase) {
        case 1:
            out << "Hexahedra ";
            break;
        case 2:
            out << " Pyramid ";
            break;
        case 4:
            out << " Tetraedra ";
            break;
        case 3:
            out << " TetraedraRef ";
            break;
        default:
            out << "Undefined ";
            break;
    }
    
    out << "POrder " << 1 << " Neq " << cmesh->NEquations() ;
    analysis.SetExact(Exact);
    TPZManVector<REAL> errvec;
    
    analysis.Solution() = cmesh->Solution();
    
    analysis.PostProcessError(errvec);
    
	// printing error
    out << " ErrH1 " << errvec[0] << " ErrL2 " << errvec[1] << " ErrH1Semi " << errvec[2];
    
	// printing error
	
    out << std::endl;
    
    
    int64_t nelem = cmesh->NElements();
    TPZFMatrix<STATE> &sol = cmesh->Solution();
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZMaterial *mat = cel->Material();
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        if (!bc) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int ncorner = gel->NCornerNodes();
        for (int ic=0; ic<ncorner; ic++) {
            TPZConnect &c = cmesh->ConnectVec()[ic];
            int64_t seqnum = c.SequenceNumber();
            STATE val = sol(cmesh->Block().Index(seqnum,0));
            if (fabs(val) > 1.e-6) {
                TPZManVector<REAL,3> x(3);
                gel->NodePtr(ic)->GetCoordinates(x);
                std::cout << "Coordinate " << x << " solution " << val << std::endl;
            }
        }
    }
	
    
    /** Cleaning allocated meshes */
    if(cmesh) {
        delete cmesh;
        cmesh = NULL;
    }
    if(gmesh) {
        delete gmesh;
        gmesh = NULL;
    }


}
void TCedricTest::LoadInterpolation(TPZCompMesh *cmesh)
{
    cmesh->Solution().Zero();
    TPZGeoMesh *gmesh = cmesh->Reference();
    TPZManVector<REAL,3> value(1);
    TPZFNMatrix<3,REAL> deriv(3,1);
    int64_t nel = gmesh->NElements();
    TPZFMatrix<STATE> &sol = cmesh->Solution();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        TPZCompEl *cel = gel->Reference();
        TPZManVector<REAL,3> x(3);
        int ncorner = gel->NCornerNodes();
        for (int ic=0; ic<ncorner; ic++) {
            gel->NodePtr(ic)->GetCoordinates(x);
            Exact(x, value, deriv);
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            sol(cmesh->Block().Index(seqnum,0)) = value[0];
        }
    }
}

#define CONDENSE

clarg::argInt  nthreads ("-nt", "Number of threads.", 8);

void TCedricTest::Run(int nsubdivisions,int geocase,int POrder,int MaterialId,std::ostream &out) {
    TPZGeoMesh *gmesh;
    switch(geocase) {
        case 1:
            gmesh = HexahedralMesh(nsubdivisions,MaterialId);
            break;
        case 2:
            gmesh = PyramidalAndTetrahedralMesh(nsubdivisions,MaterialId);
            break;
        case 3:
            gmesh = TetrahedralMesh(nsubdivisions,MaterialId);
            break;
        case 4:
            gmesh = TetrahedralMeshUsingRefinement(nsubdivisions,MaterialId);
            break;
    }
#ifdef DEFORM
    out << "Deformed ";
    DeformGMesh(*gmesh);
#else
    out << "Regular ";
    CheckConsistency(gmesh);
#endif
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int nelembc = AddBoundaryElements(gmesh);
    
    /** Generating computational mesh */

    TPZCompMesh *cmesh = GenerateCompMesh(gmesh);
    if(!cmesh) {
        if(gmesh) {
            delete gmesh;
            gmesh = NULL;
        }
        out << "Approx_err nsubdivision " << nsubdivisions << " eltype " << geocase << " aborted\n";;
        return;
    }
    
#ifdef CONDENSE
    CreateCondensedElements(cmesh);
#endif
    
    int dim = cmesh->Dimension();
    
    TPZLinearAnalysis analysis(cmesh);
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    out << "Approx_err nsubdivision " << nsubdivisions << " nelem " << (cmesh->NElements()-nelembc) << " eltype ";
    switch(geocase) {
        case 1:
            out << "Hexahedra ";
            break;
        case 2:
            out << " Pyramid ";
            break;
        case 4:
            out << " Tetraedra ";
            break;
        case 3:
            out << " TetraedraRef ";
            break;
        default:
            out << "Undefined ";
            break;
    }
    
    out << "POrder " << POrder << " Neq " << cmesh->NEquations() ;
    analysis.SetExact(Exact);
    TPZManVector<REAL> errvec;
    
//    TPZSkylineStructMatrix skylstr(cmesh);
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > skylstr(cmesh);
    skylstr.SetNumThreads(nthreads.get_value());
    analysis.SetStructuralMatrix(skylstr);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    analysis.SetSolver(step);
    
    // To post process
    /** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Solution");
    
    std::stringstream sout;
    sout << std::setprecision(2) << "Laplace_P" << POrder << "_MESH" << geocase <<  "_Div" << nsubdivisions << ".vtk";
    analysis.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
    
    analysis.Run();
        
	static int printsol = 0;
//	if(printsol == 4)
	//	analysis.PostProcess(3,dim);
//	else
	    analysis.PostProcess(0,dim);
	printsol++;
    
#ifdef CONDENSE
    UnwrapElements(cmesh);
#endif
    
    analysis.PostProcessError(errvec);
    
    out << " ErrH1 " << errvec[0] << " ErrL2 " << errvec[1] << " ErrH1Semi " << errvec[2];

	// printing error
	
    out << std::endl;
    
    /** Cleaning allocated meshes */
    if(cmesh) {
        delete cmesh;
        cmesh = NULL;
    }
    if(gmesh) {
        delete gmesh;
        gmesh = NULL;
    }
}


class ForceFunction : public TPZFunction<STATE>
{
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &df)
    {
        val[0] = 0.;
        for (int i=0; i<3; i++) {
            REAL vx = TCedricTest::fx(x[i], TCedricTest::fX0[i], TCedricTest::fEps[i]);
            int j = (i+1)%3;
            REAL vy = TCedricTest::fx(x[j], TCedricTest::fX0[j], TCedricTest::fEps[j]);
            int k = (j+1)%3;
            REAL vz = TCedricTest::d2fx(x[k], TCedricTest::fX0[k], TCedricTest::fEps[k]);
            val[0] += vx*vy*vz;
        }
    }
    
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &val)
    {
        int i = 0;
        REAL vx = TCedricTest::fx(x[i], TCedricTest::fX0[i], TCedricTest::fEps[i]);
        int j = (i+1)%3;
        REAL vy = TCedricTest::fx(x[j], TCedricTest::fX0[j], TCedricTest::fEps[j]);
        int k = (j+1)%3;
        REAL vz = TCedricTest::fx(x[k], TCedricTest::fX0[k], TCedricTest::fEps[k]);
        val[0] = vx*vy*vz;
        
    }
    virtual int NFunctions()
    {
        return 1;
    }
    
    virtual int PolynomialOrder()
    {
        return 5;
    }
    
};

static int pyramid[2][5]=
{
    {0,1,2,3,4},
    {4,5,6,7,2}
};
static int tetraedra[2][4]=
{
    {1,2,5,4},
    {4,7,3,2}
};
static int tetraedra_2[6][4]=
{
    {1,2,5,4},
    {4,7,3,2},
    {0,1,2,4},
    {0,2,3,4},
    {4,5,6,2},
    {4,6,7,2}
};


void TCedricTest::GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem)
{
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (int64_t i=0; i<=nelem; i++) {
        for (int64_t j=0; j<=nelem; j++) {
            for (int64_t k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
}

TPZGeoMesh *TCedricTest::PyramidalAndTetrahedralMesh(int64_t nelem,int MaterialId)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh, nelem);
    
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef PZ_LOG
                {
                    std::stringstream sout;
                    sout << "Pyramid and tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el=0; el<2; el++)
                {
                    TPZManVector<int64_t,5> elnodes(5);
                    for (int il=0; il<5; il++) {
                        elnodes[il] = nodes[pyramid[el][il]];
                    }
                    int64_t index;
                    gmesh->CreateGeoElement(EPiramide, elnodes, MaterialId, index);
                    elnodes.resize(4);
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

TPZGeoMesh *TCedricTest::TetrahedralMesh(int64_t nelem,int MaterialId)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh,nelem);
    
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef PZ_LOG
                {
                    std::stringstream sout;
                    sout << "Tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el=0; el<6; el++)
                {
                    TPZManVector<int64_t,4> elnodes(4);
                    int64_t index;
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra_2[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=false, const int matidtodivided=1);

TPZGeoMesh *TCedricTest::TetrahedralMeshUsingRefinement(int64_t nelemdata,int MaterialId)
{
    // CONSIDERING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // And dividing into five tetrahedras
    int nrefs = nelemdata/5;
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    REAL InitialL = 1.0;
    
    const int nelem = 5;
    const int nnode = 8;
    REAL co[nnode][3] = {
        {0.,0.,0.},
        {InitialL,0.,0.},
        {InitialL,InitialL,0.},
        {0.,InitialL,0.},
        {0.,0.,InitialL},
        {InitialL,0.,InitialL},
        {InitialL,InitialL,InitialL},
        {0.,InitialL,InitialL},
    };
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    TPZVec<TPZVec<int64_t> > indices(nelem);
    int nnodebyelement = 4;
    int el;
    for(el=0;el<nelem;el++)
        indices[el].Resize(nnodebyelement);
    // nodes to first element
    indices[0][0] = 0;
    indices[0][1] = 1;
    indices[0][2] = 3;
    indices[0][3] = 4;
    // nodes to second element
    indices[1][0] = 1;
    indices[1][1] = 2;
    indices[1][2] = 3;
    indices[1][3] = 6;
    // nodes to third element
    indices[2][0] = 4;
    indices[2][1] = 5;
    indices[2][2] = 6;
    indices[2][3] = 1;
    // nodes to fourth element
    indices[3][0] = 6;
    indices[3][1] = 7;
    indices[3][2] = 4;
    indices[3][3] = 3;
    // nodes to fifth element
    indices[4][0] = 1;
    indices[4][1] = 4;
    indices[4][2] = 6;
    indices[4][3] = 3;
    
    TPZGeoEl *elvec[nelem];
    for(el=0; el<nelem; el++) {
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement(ETetraedro,indices[el],MaterialId,index);
    }
    gmesh->BuildConnectivity();
    
    UniformRefinement(nrefs,gmesh,3);
    
    return gmesh;
}

TPZGeoMesh *TCedricTest::HexahedralMesh(int64_t nelem,int MaterialId)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh, nelem);
    
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Cube nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                int64_t index;
                gmesh->CreateGeoElement(ECube, nodes, MaterialId, index);
            }
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

/// verify if the faces without neighbour should be orthogonal to the main planes
void TCedricTest::CheckConsistency(TPZGeoMesh *mesh)
{
    int64_t nel = mesh->NElements();
    for(int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = mesh->ElementVec()[el];
        int nsides = gel->NSides();
        for(int is=0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            if(gelside.Dimension() != 2) {
                continue;
            }
            if(gelside.Neighbour() != gelside) {
                continue;
            }
            TPZManVector<REAL,2> xi(2,0.);
            gelside.CenterPoint(xi);
            TPZFNMatrix<6,REAL> axes(2,3);
            TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
            REAL detjac;
            gelside.Jacobian(xi, jac, axes, detjac, jacinv);
            TPZManVector<REAL,3> x(3,0.);
            gelside.X(xi, x);
            TPZManVector<REAL,3> normal(3);
            normal[0] = fabs(axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1));
            normal[1] = fabs(-axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0));
            normal[2] = fabs(axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0));
            REAL tol = 1.e-6;
            REAL xmin = 1., xmax = 0.;
            int numtol = 0;
            for(int i=0; i<3; i++) {
                if(xmin > x[i]) xmin = x[i];
                if(xmax < x[i]) {
                    xmax = x[i];
                }
                if(normal[i] > tol) {
                    numtol++;
                }
            }
            if(numtol != 1) {
                DebugStop();
            }
            if(xmin > tol && xmax < 1.-tol) {
                DebugStop();
            }
        }
    }
}

TPZCompMesh *TCedricTest::GenerateCompMesh(TPZGeoMesh *gmesh)
{
    int dim=3;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    /** Inserting material */
    TPZMatPoisson3d *poiss = new TPZMatPoisson3d(1,dim);
    TPZAutoPointer<TPZFunction<STATE> > force = new ForceFunction;
    poiss->SetForcingFunction(force );
    cmesh->InsertMaterialObject(poiss);
    
    /** Inserting boundary condition - Dirichlet with value 0.0 */
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    TPZBndCond *bc = new TPZBndCond(poiss, -1, 0, val1, val2);
    bc->TPZMaterial::SetForcingFunction(force);
    cmesh->InsertMaterialObject(bc);
    
    /** Constructing computational mesh */
    cmesh->AutoBuild();
    
    return cmesh;
}

int TCedricTest::AddBoundaryElements(TPZGeoMesh *gmesh)
{
    int nelembc = 0;
    int64_t nelem = gmesh->NElements();
    for (int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            TPZGeoElSide gelside(gel,is);
            if (gelside.Dimension() != 2) {
                continue;
            }
            if (gelside.Neighbour() == gelside) {
                TPZGeoElBC(gelside, -1);
                nelembc++;
            }
        }
    }
    return nelembc;
}

void TCedricTest::CreateCondensedElements(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != 3) continue;
        int64_t index;
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh,index);
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            int sidedim = gel->SideDimension(is);
            TPZCompElSide celside(cel,is);
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> elsidevec;
            gelside.ConnectedCompElementList(elsidevec, 0, 0);
            int nelside = elsidevec.size();
            for (int neigh=0; neigh<nelside; neigh++) {
                TPZCompElSide celneigh = elsidevec[neigh];
                TPZGeoElSide gelneigh = celneigh.Reference();
                TPZGeoEl *geln = gelneigh.Element();
                if (geln->Dimension() == sidedim) {
                    elgr->AddElement(celneigh.Element());
                }
            }
        }
        elgr->AddElement(intel);
    }
    cmesh->ComputeNodElCon();
    
    nel = cmesh->NElements();
    for (int el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        TPZElementGroup *celgr = dynamic_cast<TPZElementGroup *>(cel);
        if (!celgr) {
            DebugStop();
        }
        //TPZCondensedCompEl *cond = new TPZCondensedCompEl(celgr);
    }
}

void TCedricTest::UnwrapElements(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->ElementVec().NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
        if(cond)
        {
            cond->Unwrap();
        }
    }
    nel = cmesh->ElementVec().NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
        if (elgr) {
            elgr->Unwrap();
        }
    }
}
