/**
 * @file
 * @brief Projeto elaborado para resolver el problema de Poisson 2D sobre una placa plana con circunferencia de radio 1/2 centrada en (0.5,0.5)
 */

#include "pzlog.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzsbstrmatrix.h"
#include "pzfstrmatrix.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzpoisson3d.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"

#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzshapelinear.h"

#include "TPZRefPatternTools.h"

#include "pzgradient.h"
#include "pzl2projection.h"


#include <time.h>
#include <stdio.h>
#include <fstream>
#include <cmath>

using namespace std;
using namespace pzshape;

int materialId = 1;
int const matIdL2Proj = 2;
//int anothertests = 1;
char saida[512];
ofstream out("ConsolePoisson2D.txt");             // To store output of the console

STATE ValueK = 10000;

std::string Archivo = PZSOURCEDIR;

TPZGeoMesh *CreateGeoMesh();
TPZGeoMesh *CreateGeoMesh(std::string &nome);
// Crea malla computacional sem forcingfunction quando hasforcingfunction = 0, ou toma diferentes forcingfuncition para diferentes
// valores de hasforcingfunction
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int porder,int dim,bool isdiscontinuous);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<TPZVec<REAL> > &points,REAL &distance,bool &isdefined);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &points,REAL r,REAL &distance,bool &isdefined);
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);

void PrintGeoMeshVTKWithGradientAsData(TPZCompMesh *gmesh,char *filename);

void RightTermCircle(const TPZVec<REAL> &x, TPZVec<STATE> &force);

void ExactSolCircle(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);

void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points);

void formatTimeInSec(char *strtime,int timeinsec);

void SaidaSolVTK(TPZAnalysis &an, TPZCompMesh *Cmesh, std::string plotfile);
void SystemSolve(TPZAnalysis &an, TPZCompMesh *Cmesh);

/*
 *@brief Method to replace the solution by finite element method from the reconstruction of the gradient.
 *Using the L2 projection
 *@param cmesh in:  Computational mesh
 *@param var in: Index of the Solution Variable
 *@param matid_l2proj in: Id of the l2 projection material
 */
void ProjectionGradientReconstructedInFESpace(TPZCompMesh *cmesh,int var, int matid_l2proj);

/*
 *@brief Method to reconstruction gradient by Least Squares
 *@param cel in:  Computational element
 *@param center out:  Center point of the element
 *@param solalfa out: Value of the approximate solution (uh) at the center point
 *@param grad out: Value of the gradient reconstructed
 *@param var in: Index of the Solution Variable
 */
void GradientReconstructionByLeastSquares(TPZCompEl *cel,TPZManVector<REAL,3> &center, TPZManVector<STATE> &solalfa,  TPZFMatrix<STATE> &grad, int var);

//Trocar todos os elementos do cmesh apontando para o material TPZL2ProjectionFromGradient
void ChangeMaterialIdIntoCompElement(TPZCompEl *cel, int oldmatid, int newmatid);

void AssembleGlobalMatrix(TPZCompEl *el, TPZElementMatrix &ek, TPZElementMatrix &ef,TPZMatrix<STATE> & stiffmatrix, TPZFMatrix<STATE> &rhs);

#ifdef LOG4CXX
static PZLogger logger("pz.Poisson2D_ArcTan_Disc");

#endif

// MAIN FUNCTION TO NUMERICAL SOLVE
/** Laplace equation on square - Volker John article 2000 */

int main() {
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeAllUniformRefPatterns();
    // To compute processing times
    time_t sttime;
    time_t endtime;
    int time_elapsed;
    char tempo[256];
    
    ofstream fileerrors("ErrorsHP2D_ArcTan.txt");   // To store all errors calculated by TPZAnalysis (PosProcess)
    time (& sttime);
    
    int dim = 2;
    int p = 2;
    
    //Sem adaptatividade
    TPZGeoMesh *gmesh = CreateGeoMesh();
    UniformRefine(gmesh, 2);
    RefiningNearCircunference(dim,gmesh,3,1);
    
    TPZCompMesh *cmesh = CreateMesh(gmesh, p, dim, true);
    cmesh->ComputeNodElCon();
    
#ifdef LOG4CXX
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
    
#endif
    
    TPZAnalysis an(cmesh);
    SystemSolve(an, cmesh);
    
    // Calculando o tempo que demorou para calcular em cada cenario
    time (& endtime);
    time_elapsed = endtime - sttime;
    formatTimeInSec(tempo, time_elapsed);
    // Printing information
    fileerrors << "Approximation Error: Solving time = " << time_elapsed << std::endl;
    
    //char filename[256];
    //sprintf(filename,"Poisson2DGradient.vtk");
    //PrintGeoMeshVTKWithGradientAsData(cmesh, filename);
    std::string filename = "SolutionSemRecGrad.vtk";
    SaidaSolVTK(an, cmesh, filename);
    
    
    //------ RESOLVENDO COM RECONST. GRADIENT ------
    
    //l2 projection of the gradient into finite element space
    ProjectionGradientReconstructedInFESpace(cmesh,1, matIdL2Proj);
    an.LoadSolution(cmesh->Solution());
    std::string filename2 = "L2PROJSolution.vtk";
    SaidaSolVTK(an,cmesh,filename2);
    
    cmesh->CleanUp();
    delete cmesh;
    delete gmesh;
    
    return 0;
}
void PrintGeoMeshVTKWithGradientAsData(TPZCompMesh *cmesh,char *filename) {
    int i, size = cmesh->NElements();
    int dim = cmesh->Dimension();
    TPZChunkVector<REAL> DataElement;
    DataElement.Resize(dim*size);
    // Making dimension of the elements as data element
    for(i=0;i<size;i++) {
        if(cmesh->ElementVec()[i])
            DataElement[i] = (cmesh->ElementVec()[i])->Dimension();
        else
            DataElement[i] = -10.0;
    }
    // Printing geometric mesh to visualization in Paraview
    //TPZVTKGeoMesh::PrintGMeshVTK(cmesh, filename, DataElement);
}

void ProjectionGradientReconstructedInFESpace(TPZCompMesh *cmesh,int var, int matid_l2proj){
    
    // Redimensionando a matriz dos dados da reconstruca de gradientes
    int dim  = cmesh->Dimension();
    int nelem = cmesh->NElements();
    
    TPZManVector<REAL,3> center;
    TPZManVector<STATE> solalfa;
    TPZFMatrix<STATE> Grad;
    
    //criar ponteiro para TPZFunction
    TPZGradient *pGrad = new TPZGradient;
    TPZAutoPointer<TPZFunction<STATE> > fp(pGrad);
    
    //Criar matrix de rigidez e vetor de carga
    int numloadcases=0;
    unsigned int im;
    for(im=0; im<cmesh->MaterialVec().size(); im++){
        if(!cmesh->MaterialVec()[im]) continue;
        numloadcases = cmesh->MaterialVec()[im]->NumLoadCases();
        break;
    }
    
    int neq = cmesh->NEquations();
    TPZFMatrix<STATE> rhs;
    rhs.Redim(neq,numloadcases);
    
    //TPZBandStructMatrix stmatrix(cmesh);
    TPZSkylineStructMatrix stmatrix(cmesh);
    TPZMatrix<STATE> *stiffmatrix = stmatrix.Create();
    
    int matid;
    
    for(int i=0; i<nelem; i++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZElementMatrix ek(cel->Mesh(), TPZElementMatrix::EK);
        TPZElementMatrix ef(cel->Mesh(), TPZElementMatrix::EF);
        
        // Nada sera realizado para elementos con dimension diferente de la dimension del problema
        if(cel->Dimension()!=dim) continue;
        
        matid = cel->Material()->Id();
        
        //gradient reconstruction
        GradientReconstructionByLeastSquares(cel, center, solalfa, Grad,var);
        
        TPZManVector<STATE,3> gradient;
        
        gradient.Resize(Grad.Rows());
        for (int i = 0; i < Grad.Rows(); i++) {
            gradient[i] = Grad(i);
        }
        
        //set data of the gradient reconstructed
        pGrad->SetData(center, gradient, (STATE)1., solalfa[0]);
        
        //change material id current to material id of the L2ProjectionMaterial
        ChangeMaterialIdIntoCompElement(cel, matid, matid_l2proj);
        
        //set forcing function of l2 projection material
        TPZMaterial *mat = cel->Material();
        mat->SetForcingFunction(fp);
        
        //load the matrix ek and vector ef of the element
        cel->CalcStiff(ek,ef);
        
        //assemble pos l2 projection
        AssembleGlobalMatrix(cel, ek, ef, *stiffmatrix, rhs);
        
        //Return for original material and current solution of the mesh
        ChangeMaterialIdIntoCompElement(cel, matid_l2proj, matid);
    }
    
    
    //Solve linear system and transfer the solution to computational mesh
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    step.SetMatrix(stiffmatrix);
    TPZFMatrix<STATE> result;
    step.Solve(rhs, result);
    cmesh->LoadSolution(result);
}

void GradientReconstructionByLeastSquares(TPZCompEl *cel,TPZManVector<REAL,3> &center, TPZManVector<STATE> &solalfa,  TPZFMatrix<STATE> &grad, int var) {
    
    int dim;
    dim = cel->Mesh()->Dimension();
    
    // Nada sera realizado para elementos com dimensao diferente da dimensao do problema
    if(!cel || cel->Dimension()!=dim) DebugStop();
    
    int nstates=0;
    nstates = cel->Material()->NSolutionVariables(var);
    
    int k, side;
    
    TPZStack<TPZCompElSide> neighs;
    int nneighs;
    
    center.Resize(3, 0.);
    solalfa.Resize(nstates,0.0);
    TPZManVector<REAL> centerpsi(3,0.0), centerbeta(3,0.0);
    TPZManVector<STATE> solbeta(nstates,0.0);
    
    TPZFMatrix<STATE> A(dim,dim);  // Linear System matrix
    grad.Redim(dim,1);
    
    // matrizes para aplicar o metodo dos minimos quadrados
    TPZFMatrix<STATE> DeltaH;
    TPZFMatrix<STATE> DeltaHTranspose;
    TPZFMatrix<STATE> DifSol;
    
    // Encontrando o centro do elemento atual (cel)
    TPZGeoEl* gelalfa = cel->Reference();
    gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
    center.Fill(0.);
    gelalfa->X(centerpsi,center);
    cel->Solution(centerpsi,var,solalfa);
    
    neighs.Resize(0);
    
    // Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
    for(side = 0; side <cel->Reference()->NSides()-1; side++)
    {
        TPZCompElSide celside(cel,side);
        celside.ConnectedElementList(neighs,0,0);
    }
    
    // si no hay vecinos continuamos con el siguiente elemento
    nneighs = neighs.NElements();
    if(!nneighs) DebugStop();
    
    std::set<TPZCompEl *> neighscel;
    for(int i =0; i<nneighs; i++)
    {
        TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(neighs[i].Element());
        if(!InterpEl || InterpEl->Dimension()!=dim) continue;
        neighscel.insert(neighs[i].Element());
    }
    
    nneighs=neighscel.size();
    
    // si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente
    // Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
    // y el valor de la solucion en su centro solbeta
    DeltaH.Redim(nneighs,dim);
    DeltaHTranspose.Redim(dim,nneighs);
    DifSol.Redim(nneighs,1);
    
    // Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
    int ineighs=-1;
    int counter=0;
    std::set<TPZCompEl *>::iterator it;
    for(it=neighscel.begin(); it!=neighscel.end(); ++it)
    {
        //(*it)->Print();
        ineighs++;
        TPZGeoEl * gelbeta = (*it)->Reference();
        
        if(!gelbeta) DebugStop();
        
        centerpsi.Fill(0.0);
        centerbeta.Fill(0.0);
        gelbeta->CenterPoint(gelbeta->NSides()-1,centerpsi);
        gelbeta->X(centerpsi,centerbeta);
        gelbeta->Reference()->Solution(centerpsi,var,solbeta);
        
        for(k=0;k<dim;k++)
        {
            DeltaH(ineighs,k) = centerbeta[k] - center[k];
        }
        DifSol(ineighs,0) = solbeta[nstates-1] - solalfa[nstates-1];
        
        counter ++;
    }
    
    // Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    grad = DeltaHTranspose*DifSol;
    A = DeltaHTranspose*DeltaH;
    if(counter > 0){
        A.SolveDirect(grad,ELU);
    }
}

void ChangeMaterialIdIntoCompElement(TPZCompEl *cel, int oldmatid, int newmatid) {
    
    // Changes material Id only elements with required id (matid)
    if(cel->Material()->Id() != oldmatid) return;
    
    //mudar o material id
    TPZGeoEl *gel;
    gel = cel->Reference();
    gel->SetMaterialId(newmatid);
}

void AssembleGlobalMatrix(TPZCompEl *el, TPZElementMatrix &ek, TPZElementMatrix &ef,TPZMatrix<STATE> & stiffmatrix, TPZFMatrix<STATE> &rhs){
    
    if(!el->HasDependency()) {
        ek.ComputeDestinationIndices();
        
        stiffmatrix.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        
    } else {
        // the element has dependent nodes
        ek.ApplyConstraints();
        ef.ApplyConstraints();
        ek.ComputeDestinationIndices();
        
        stiffmatrix.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
    }
}

void SystemSolve(TPZAnalysis &an, TPZCompMesh *Cmesh){
    
    TPZSkylineNSymStructMatrix strskyl(Cmesh);
    //    TPZBandStructMatrix strskyl(Cmesh);
    an.SetStructuralMatrix(strskyl);
    
    TPZStepSolver<STATE> direct;
    direct.SetDirect(ELU);
    an.SetSolver(direct);
    
    // Solving
    an.Run();
}

void SaidaSolVTK(TPZAnalysis &an, TPZCompMesh *Cmesh, std::string plotfile){
    
    TPZStack<std::string> scalarnames, vecnames;
    const int dim = Cmesh->Dimension();
    int div = 0;
    
    scalarnames.Push("Solution");
    scalarnames.Push("POrder");
    scalarnames.Push("KDuDx");
    scalarnames.Push("KDuDy");
    
    vecnames.Push("Derivative");
    vecnames.Push("MinusKGradU");
    
    an.DefineGraphMesh(dim,scalarnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    
}

void formatTimeInSec(char *strtime,int timeinsec) {
    if(!strtime) return;
    memset(strtime,0,strlen(strtime));
    //	strtime[0] = '\0';
    int anos=0, meses=0, dias=0, horas=0, minutos=0, segundos=0;
    while(1) {
        if(timeinsec < 60) {
            segundos = timeinsec;
            break;
        }
        else {
            timeinsec -= 60;
            minutos++;
            if(minutos > 59) {
                minutos -= 60;
                horas++;
                if(horas > 23) {
                    horas -= 24;
                    dias++;
                    if(dias > 29) {
                        dias -= 30;
                        meses++;
                        if(meses > 11) {
                            meses -= 12;
                            anos++;
                        }
                    }
                }
            }
        }
    }
    // Formating
    if(anos)
        sprintf(strtime,"%d a, %d m, %d d, %02d:%02d:%02d",anos,meses,dias,horas,minutos,segundos);
    else {
        if(meses)
            sprintf(strtime,"%d m, %d d, %02d:%02d:%02d",meses,dias,horas,minutos,segundos);
        else {
            if(dias)
                sprintf(strtime,"%d d, %02d:%02d:%02d",dias,horas,minutos,segundos);
            else
                sprintf(strtime,"%02d:%02d:%02d",horas,minutos,segundos);
        }
    }
}

void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points) {
    Points.Resize(npoints);
    TPZManVector<REAL> point(3,0.);
    REAL angle = (2*M_PI)/npoints;
    for(int i=0;i<npoints;i++) {
        point[0] = center[0]+radius*cos(i*angle);
        point[1] = center[1]+radius*sin(i*angle);
        Points[i] = point;
    }
}

void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
    
    int i;
    bool isdefined = false;
    
    // Refinando no local desejado
    int npoints = 1000;
    TPZVec<REAL> point(3);
    point[0] = point[1] = 0.5; point[2] = 0.0;
    REAL r = 0.25;
    TPZVec<TPZManVector<REAL> > Points(npoints);
    GetPointsOnCircunference(npoints,point,r,Points);
    
    if(ntyperefs==2) {
        REAL radius = 0.19;
        for(i=0;i<nref;i+=2) {
            // To refine elements with center near to points than radius
            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
            if(nref < 5) radius *= 0.35;
            else if(nref < 7) radius *= 0.2;
            else radius *= 0.1;
        }
        if(i==nref) {
            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
        }
    }
    else {
        REAL radius = 0.2;
        for(i=0;i<nref+1;i++) {
            // To refine elements with center near to points than radius
            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
            if(nref < 5) radius *= 0.6;
            else if(nref < 7) radius *= 0.3;
            else radius *= 0.15;
        }
    }
    // Constructing connectivities
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
    TPZManVector<REAL> centerpsi(3), center(3);
    // Refinamento de elementos selecionados
    TPZGeoEl *gel;
    TPZVec<TPZGeoEl *> sub;
    
    int nelem = 0;
    int ngelem=gmesh->NElements();
    // na esquina inferior esquerda Nó = (0,-1,0)
    while(nelem<ngelem) {
        gel = gmesh->ElementVec()[nelem++];
        if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
        gel->CenterPoint(gel->NSides()-1,centerpsi);
        gel->X(centerpsi,center);
        if(!isdefined) {
            TPZVec<REAL> FirstNode(3,0.);
            gel->CenterPoint(0,centerpsi);
            gel->X(centerpsi,FirstNode);
            REAL distancia = TPZGeoEl::Distance(center,FirstNode);
            if(distancia > distance) distance = distancia;
            isdefined = true;
        }
        REAL centerdist = TPZGeoEl::Distance(center,point);
        if(fabs(r-centerdist) < distance) {
            gel->Divide(sub);
        }
    }
}

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<TPZVec<REAL> > &points,REAL &distance,bool &isdefined) {
    TPZManVector<REAL> centerpsi(3), center(3);
    // Refinamento de elementos selecionados
    TPZGeoEl *gel;
    TPZVec<TPZGeoEl *> sub;
    
    int nelem = 0;
    int ngelem=gmesh->NElements();
    int i, npoints = points.NElements();
    // na esquina inferior esquerda Nó = (0,-1,0)
    while(nelem<ngelem) {
        gel = gmesh->ElementVec()[nelem++];
        if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
        gel->CenterPoint(gel->NSides()-1,centerpsi);
        gel->X(centerpsi,center);
        if(!isdefined) {
            TPZVec<REAL> FirstNode(3,0.);
            gel->CenterPoint(0,centerpsi);
            gel->X(centerpsi,FirstNode);
            distance = 1.1*TPZGeoEl::Distance(center,FirstNode);
            isdefined = true;
        }
        for(i=0;i<npoints;i++) {
            REAL semidiag = TPZGeoEl::Distance(center,points[i]);
            if(semidiag < distance) {
                gel->Divide(sub);
                break;
            }
        }
    }
}


#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

//**** Creating Geometric Mesh as square */
TPZGeoMesh *CreateGeoMesh() {
    
    int Qnodes = 4;
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolLine(2);
    
    //indice dos nos
    int64_t id;
    id = 0;
    REAL dx = 1.;
    for (int i=0; i<Qnodes/2;i++) {
        
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0, i*dx);//coord x
        Node[id].SetCoord(1, 0.);//coord y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    for (int i=0; i<Qnodes/2;i++) {
        
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0, 1. - i*dx);//coord x
        Node[id].SetCoord(1, 1.);//coord y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    //indice dos elementos
    id = 0;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 2;
    TopolQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,materialId,*gmesh);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
    
    gmesh->BuildConnectivity();
    return gmesh;
}

//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh(std::string &archivo) {
    
    // Ejemplo uni-dimensional para la generacion de una malla para un reservatorio
    TPZReadGIDGrid grid;
    TPZGeoMesh *meshgrid = grid.GeometricGIDMesh(archivo);
    if(!meshgrid->NElements())
        return 0;
    
    return meshgrid;
}

//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int porder, int dim, bool isdiscontinuous) {
    
    // Creating Poisson material
    TPZMatPoisson3d *material = new TPZMatPoisson3d(materialId,dim);
    TPZVec<REAL> convd(3,0.);
    material->SetParameters(ValueK,0.,convd);
    material->SetNoPenalty();
    material->SetNonSymmetric();
    
    TPZAutoPointer<TPZFunction<STATE> > myforce = new TPZDummyFunction<STATE>(RightTermCircle,5);
    material->SetForcingFunction(myforce);
    
    
    TPZMaterial * mat(material);
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetDimModel(material->Dimension());
    
    // Creating four boundary condition
    TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);
    
    // Condicion de Dirichlet fijando la posicion de la placa
    TPZMaterial * BCond = material->CreateBC(mat,-1, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond);
    
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,material->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
    
    
    if (isdiscontinuous==true) {
        //Set discontinuous functions
        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->AutoBuild();
        cmesh->ExpandSolution();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
        //Adding interface element
        for(int64_t el = 0; el < cmesh->ElementVec().NElements(); el++)
        {
            TPZCompEl * compEl = cmesh->ElementVec()[el];
            if(!compEl) continue;
            int64_t index = compEl ->Index();
            if(compEl->Dimension() == cmesh->Dimension())
            {
                TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[index]);
                if(!InterpEl) continue;
                InterpEl->CreateInterfaces(false);
            }
        }
    }
    else{
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();
        cmesh->ExpandSolution();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    }
    return cmesh;
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
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
    // Re-constructing connectivities
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

/** We are considering - f, because is as TPZMatPoisson3d was implemented in Contribute method */
void RightTermCircle(const TPZVec<REAL> &x, TPZVec<STATE> &force) {
    //	REAL Epsilon = 1000000;
    REAL B = (16.*ValueK)/M_PI;
    REAL F = 2*sqrt(ValueK);
    REAL G = -0.4375;
    
    REAL sum = x[0]*(x[0]-1) + x[1]*(x[1]-1);
    REAL prod = x[0]*(x[0]-1)*x[1]*(x[1]-1);
    
    REAL temp = F*(G-sum);
    REAL arctan = atan(temp);
    REAL den = (1+temp*temp)*(1+temp*temp);
    REAL num = 2*F*(sum*(2*F*F*prod*(8*G+1)-(1+F*F*G*G)+F*F*sum*(2*G-6*prod-sum))-2*prod*(F*F*G+5*F*F*G*G+5));
    
    force[0] = B*(sum*(M_PI+2*arctan)+(num/den));
    /*	REAL B = (-16.0*ValueK)/M_PI;
     // Computing Q(x,y) = 2*Sqrt[ValueK]*(.25^2-(x-0.5)^2-(y-0.5)^2)  Doing F = -0.5*Sqrt[ValueK]  Then Q=F*(7/4 + 4 Sum)
     REAL F = (-0.5)*sqrt(ValueK);
     REAL sum = x[0]*(x[0]-1.) + x[1]*(x[1]-1.);
     REAL temp = F*((7./4.)+4*sum);
     REAL arctan = atan(temp);
     
     REAL prod = x[0]*(x[0]-1.)*x[1]*(x[1]-1.);
     REAL den = (1+temp*temp);
     
     force[0] = B*(sum*(M_PI+2.*arctan)+((8*F*(2.+5*sum))/den)-((32*F*prod*temp)/(den*den)));*/
}

void ExactSolCircle(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
    /*	REAL B = 8./M_PI;
     REAL F = (-0.5)*sqrt(ValueK);
     REAL prodx = x[0]*(x[0]-1.);
     REAL prody = x[1]*(x[1]-1.);
     REAL sum = prodx + prody;
     REAL temp = F*((7./4.)+4*sum);
     REAL arctan = atan(temp);
     
     REAL prod = prodx*prody;
     // Solution
     sol[0] = B*prod*(M_PI+ 2*arctan);
     // Partial derivaties
     REAL den = 1.+temp*temp;
     dsol(0,0) = B*prody*(2*x[0] - 1.)*(M_PI+(2*arctan)+((8*F*prodx)/den));
     dsol(1,0) = B*prodx*(2*x[1] - 1.)*(M_PI+(2*arctan)+((8*F*prody)/den));*/
    REAL F = 2*sqrt(ValueK);
    REAL arc = F*((0.25*0.25) - (x[0] - 0.5)*(x[0] - 0.5) - (x[1] - 0.5)*(x[1] - 0.5));
    REAL prodx = x[0]*(x[0]-1.);
    REAL prody = x[1]*(x[1]-1.);
    REAL prod = prodx*prody;
    sol[0] = 8*prod*(1+(2./M_PI)*(atan(arc)));
    REAL temp = prody*(2*x[0]-1.)*(M_PI + 2*atan(arc));
    REAL frac = 2*prod*F*(1.-2*x[0]);
    frac = frac/(1+arc*arc);
    dsol(0,0) = (8./M_PI)*(temp + frac);
    temp = prodx*(2*x[1]-1.)*(M_PI + 2*atan(arc));
    frac = 2*prod*F*(1.-2*x[1]);
    frac = frac/(1+arc*arc);
    dsol(1,0) = (8./ M_PI)*(temp + frac);    
}
REAL PartialDerivateX(const TPZVec<REAL> &x) {
    /*	REAL F = 2*sqrt(ValueK);
     REAL arc = F*((0.25*0.25) - (x[0] - 0.5)*(x[0] - 0.5) - (x[1] - 0.5)*(x[1] - 0.5));
     REAL prodx = x[0]*(x[0]-1.);
     REAL prody = x[1]*(x[1]-1.);
     REAL result = (8./M_PI)*prody*(2*x[0]-1);
     REAL temp = M_PI + 2*atan(arc);
     REAL frac = 2*F*prodx;
     frac = frac/(1+arc*arc);
     temp -= frac;
     return (result*temp);*/
    REAL B = 8./M_PI;
    REAL F = (-0.5)*sqrt(ValueK);
    
    REAL prodx = x[0]*(x[0]-1.);
    REAL prody = x[1]*(x[1]-1.);
    REAL sum = prodx + prody;
    
    REAL temp = F*((7./4.)+4*sum);
    REAL arctan = atan(temp);
    REAL den = 1.+temp*temp;
    return ( B*prody*(2*x[0] - 1.)*(M_PI+(2*arctan)+((8*F*prodx)/den)));
}

REAL PartialDerivateY(const TPZVec<REAL> &x) {
    /*	REAL F = 2*sqrt(ValueK);
     REAL arc = F*((0.25*0.25) - (x[0] - 0.5)*(x[0] - 0.5) - (x[1] - 0.5)*(x[1] - 0.5));
     REAL prodx = x[0]*(x[0]-1.);
     REAL prody = x[1]*(x[1]-1.);
     REAL result = (8./M_PI)*prodx*(2*x[1]-1);
     REAL temp = M_PI + 2*atan(arc);
     REAL frac = 2*F*prody;
     frac = frac/(1+arc*arc);
     temp -= frac;
     return (result*temp);*/
    REAL B = 8./M_PI;
    REAL F = (-0.5)*sqrt(ValueK);
    
    REAL prodx = x[0]*(x[0]-1.);
    REAL prody = x[1]*(x[1]-1.);
    REAL sum = prodx + prody;
    
    REAL temp = F*((7./4.)+4*sum);
    REAL arctan = atan(temp);
    REAL den = 1.+temp*temp;
    
    return (B*prodx*(2*x[1] - 1.)*(M_PI+(2*arctan)+((8*F*prody)/den)));
}
