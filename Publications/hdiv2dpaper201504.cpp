//
//  Hdiv2dPaper201504.cpp
//  PZ
//
//  Created by Douglas Castro on 18/03/15.
//
//

#include "hdiv2dpaper201504.h"

//#define DEFORMED
#define SENOSENO

Hdiv2dPaper201504::Hdiv2dPaper201504()
{
    
    fDim = 2;
    
    fmatId = 1;
    
    fdirichlet = 0;
    fneumann = 1;
    
    fbc1 = -2;
    fbc2 = -3;
    fbc3 = -4;
    fbc4 = -5;
    
    fmatskeleton = -7;
    
    ftriang = false;
}

Hdiv2dPaper201504::~Hdiv2dPaper201504()
{
    
}

void Hdiv2dPaper201504::Run(ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZFMatrix< REAL > &errors)
{
    if (POrderBeginAndEnd.size()!=2)
    {
        std::cout << " POrderBeginAndEnd Vector must have size == 2 and contain just two values, the first and last approximation orders. " << std::endl;
        std::cout << " Some like POrderBeginAndEnd[0] = 1, POrderBeginAndEnd[1] = 3 " << std::endl;
        return;
    }
    
    if (ndivinterval.size()!=2)
    {
        std::cout << " ndivinterval Vector must have size == 2 and contain just two values, the first and last number of uniform refinements. " << std::endl;
        std::cout << " Some like ndivinterval[0] = 0, ndivinterval[1] = 3 for simulations with 0, 1, 2 and 3 refinements " << std::endl;
        std::cout << " Some like ndivinterval[0] = 1, ndivinterval[1] = 1 for simulations with 1 refinement " << std::endl;
        return;
    }
    
    int sizeErrorMatrix = (ndivinterval[1]-ndivinterval[0])*(POrderBeginAndEnd[1] - POrderBeginAndEnd[0]) + 1;
    int errorposition = 0;
    
    errors.Resize(sizeErrorMatrix, 4);
    
    for(int p = POrderBeginAndEnd[0]; p <= POrderBeginAndEnd[1]; p++)
    {
        
        for (int ndiv = ndivinterval[0]; ndiv<=ndivinterval[1]; ndiv++)
        {
            int ordemP = p;
            
            std::cout<< " BEGIN (Quadrilateral Domain) - Polinomial degree: " << ordemP << " rerfinement size (h): " << ndiv << std::endl;
        
            TPZGeoMesh *gmesh;
            switch (element) {
                case EQuad:
                {
                    gmesh = this->GMesh(fDim, false, ndiv);
                }
                    break;
                case ETriangle:
                {
                    setTriangTrue();
                    gmesh = this->GMesh(fDim, true, ndiv);
                }
                    break;
                default:
                    break;
            }
            
            
            switch (problem) {
                case EH1:
                {
                    
                    gmesh->SetDimension(fDim);
                    TPZCompMesh *cmeshH1 = this->CMeshH1(gmesh, ordemP, fDim);
                    //condensar
                    for (int64_t iel=0; iel<cmeshH1->NElements(); iel++) {
                        TPZCompEl *cel = cmeshH1->Element(iel);
                        if(!cel) continue;
                        new TPZCondensedCompEl(cel);
                    }
                    
                    cmeshH1->ExpandSolution();
                    cmeshH1->CleanUpUnconnectedNodes();
                    
                    TPZAnalysis anh1(cmeshH1, true);
                    
                    SolveSyst(anh1, cmeshH1);
                    
                    ErrorH1(cmeshH1, ordemP, ndiv, errorposition, errors);

                }
                    break;
                case EHDiv:
                {
                    DebugStop(); //Mudanca na topologia
                }
                    break;
                case EHDivStar:
                {

                    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
                    
                    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
                    
                    // P**
                    ChangeExternalOrderConnects(cmesh1);
                    
                    
                    //malha multifisica
                    TPZVec<TPZCompMesh *> meshvec(2);
                    meshvec[0] = cmesh1;
                    meshvec[1] = cmesh2;
                    
                    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
                    
                    TPZAnalysis an(mphysics, true);
                    
                    SolveSyst(an, mphysics);
                    
                    //Calculo do erro
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    std::cout << "Computing Errors\n";
                    
                    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, errorposition, errors);
                }
                    break;
                    
                case EHDivStarStar:
                {

                    ordemP  =  ordemP + 1;
                    
                    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
                    
                    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
                    //P**
                    ChangeExternalOrderConnects(cmesh1);
                    
                    //malha multifisica
                    TPZVec<TPZCompMesh *> meshvec(2);
                    meshvec[0] = cmesh1;
                    meshvec[1] = cmesh2;
                    
                    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
                    
                    TPZAnalysis an(mphysics, true);
                    
                    SolveSyst(an, mphysics);
                    
                    //Calculo do erro
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    std::cout << "Computing Errors\n";
                    
                    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, errorposition, errors);
                }
                    break;
                    
                default:
                    break;
            }
            
            std::cout<< " END (QUAD) - polinomial degree: " << ordemP << " refinement size (h): " << ndiv << std::endl;
            
            errorposition++;

        }
    }
    
    std::cout<< " The End " << std::endl;
    
    
}

void Hdiv2dPaper201504::PrintErrors(ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZVec<REAL> &errors, std::ostream &output)
{
    if (POrderBeginAndEnd.size()!=2)
    {
        std::cout << " POrderBeginAndEnd Vector must have size == 2 and contain just two values, the first and last approximation orders. " << std::endl;
        std::cout << " Some like POrderBeginAndEnd[0] = 1, POrderBeginAndEnd[1] = 3 " << std::endl;
        return;
    }
    
    if (ndivinterval.size()!=2)
    {
        std::cout << " ndivinterval Vector must have size == 2 and contain just two values, the first and last number of uniform refinements. " << std::endl;
        std::cout << " Some like ndivinterval[0] = 0, ndivinterval[1] = 3 for simulations with 0, 1, 2 and 3 refinements " << std::endl;
        std::cout << " Some like ndivinterval[0] = 1, ndivinterval[1] = 1 for simulations with 1 refinement " << std::endl;
        return;
    }
    
    for(int p = POrderBeginAndEnd[0]; p <= POrderBeginAndEnd[1]; p++)
    {
                output << "\n WHEN p = " << p << " \n " << endl;
                output << "ndiv " << setw(6) << "DoFTot" << setw(20) << "DofCond" << setw(28) << "PrimalL2Error" << setw(35) << "L2DualError"  << endl;
        
        for (int ndiv = ndivinterval[0]; ndiv<=ndivinterval[1]; ndiv++)
        {
            int ordemP = p;
            
            std::cout<< " BEGIN (Quadrilateral Domain) - Polinomial degree: " << ordemP << " rerfinement size (h): " << ndiv << std::endl;
            
            TPZGeoMesh *gmesh;
            switch (element) {
                case EQuad:
                {
                    gmesh = this->GMesh(fDim, false, ndiv);
                }
                    break;
                case ETriangle:
                {
                    setTriangTrue();
                    gmesh = this->GMesh(fDim, true, ndiv);
                }
                    break;
                default:
                    break;
            }
            
            gmesh->SetDimension(fDim);
            
            switch (problem) {
                case EH1:
                {
                    
                    TPZCompMesh *cmeshH1 = this->CMeshH1(gmesh, ordemP, fDim);
                    
                    int dofTotal = cmeshH1->NEquations();
                    
                    //condensar
                    for (int64_t iel=0; iel<cmeshH1->NElements(); iel++) {
                        TPZCompEl *cel = cmeshH1->Element(iel);
                        if(!cel) continue;
                        new TPZCondensedCompEl(cel);
                    }
                    
                    cmeshH1->ExpandSolution();
                    cmeshH1->CleanUpUnconnectedNodes();
                    
                    int dofCondensed = cmeshH1->NEquations();
                    
                    TPZAnalysis anh1(cmeshH1, true);
                    
                    SolveSyst(anh1, cmeshH1);
                    
                    ErrorH1(cmeshH1, ordemP, ndiv, output, dofTotal, dofCondensed);
                    
                }
                    break;
                case EHDiv:
                {
                    DebugStop(); //Mudanca na topologia
                }
                    break;
                case EHDivStar:
                {

                    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
                    
                    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
                    
                    // para rodar P**
                    ChangeExternalOrderConnects(cmesh1);
                    
                    int DofCond, DoFT;
                    DoFT = cmesh1->NEquations() + cmesh2->NEquations();
                    
                    //malha multifisica
                    TPZVec<TPZCompMesh *> meshvec(2);
                    meshvec[0] = cmesh1;
                    meshvec[1] = cmesh2;
                    
                    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
                    
                    DofCond = mphysics->NEquations();
                    
                    TPZAnalysis an(mphysics, true);
                    
                    SolveSyst(an, mphysics);
                    
                    //Calculo do erro
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    std::cout << "Computing Errors\n";
                    
                    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, output, DoFT, DofCond);
                }
                    break;
                    
                case EHDivStarStar:
                {

                    
                    ordemP  =  ordemP + 1;
                    
                    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
                    
                    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
                    
                    // para rodar P**
                    ChangeExternalOrderConnects(cmesh1);
                    
                    int DofCond, DoFT;
                    DoFT = cmesh1->NEquations() + cmesh2->NEquations();
                    
                    //malha multifisica
                    TPZVec<TPZCompMesh *> meshvec(2);
                    meshvec[0] = cmesh1;
                    meshvec[1] = cmesh2;
                    
                    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
                    
                    DofCond = mphysics->NEquations();
                    
                    TPZAnalysis an(mphysics, true);
                    
                    SolveSyst(an, mphysics);
                    
                    //Calculo do erro
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    std::cout << "Computing Errors\n";
                    
                    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, output, DoFT, DofCond);
                }
                    break;
                    
                default:
                    break;
            }

            
        }
        
        output << "\n ----------------------------------------------------------------------------- " << endl;
        std::cout<< " END (Quadrilateral Domain) - Polinomial degree: " << p << std::endl;
        
    }
    
    std::cout<< " The End " << std::endl;
    
    
}

TPZGeoMesh *Hdiv2dPaper201504::GMesh(int dim, bool ftriang, int ndiv)
{
    
    if(dim!=2)
    {
        std::cout << "Wrong dimension" << std::endl;
        dim = 2;
        DebugStop();
    }
    
    int Qnodes =  4;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    gmesh->SetDimension(dim);
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
    TPZVec <int64_t> TopolLine(2);
    
    //indice dos nos
    int64_t id = 0;

    TPZManVector<REAL,3> coord(2,0.);
    int in = 0;
    //c0
    coord[0] = 0.0;
    coord[1] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c1
    coord[0] =  1.0;
    coord[1] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c2
    coord[0] =  0.0;
    coord[1] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c3
    coord[0] = 1.0;
    coord[1] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //indice dos elementos
    id = 0;
    
    if(ftriang) // triangulo
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,fmatId,*gmesh);
        id++;
        
        TopolTriang[0] = 0;
        TopolTriang[1] = 3;
        TopolTriang[2] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,fmatId,*gmesh);
        id++;
        
    }
    else{
        
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 3;
        TopolQuad[3] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,fmatId,*gmesh);
        id++;
    }
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc1,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc2,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc3,*gmesh);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc4,*gmesh);
    
    
    
    gmesh->BuildConnectivity();
    
    /// gmesh para aqui
    
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < ndiv; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    return gmesh;
}

void Hdiv2dPaper201504::SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    solp[0]=0.;
    
    int dim = 3; //getDimension();
    
    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    flux.Resize(dim, 1);
    
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = sin(M_PI*x)*sin(M_PI*y);
    flux(0,0) = -M_PI*cos(M_PI*x)*sin(M_PI*y);
    flux(1,0) = -M_PI*cos(M_PI*y)*sin(M_PI*x);
    
}

void Hdiv2dPaper201504::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){
    
    int dim = 3; //getDimension();
    
    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    double x = pt[0];
    double y = pt[1];
    
    ff[0] = 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    

    
    
}

void Hdiv2dPaper201504::SolExataH1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    solp[0]=0.;
    
    int dim = 3; //getDimension();
    
    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    flux.Resize(dim, 1);
    
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = sin(M_PI*x)*sin(M_PI*y);
    flux(0,0) = -M_PI*cos(M_PI*x)*sin(M_PI*y);
    flux(1,0) = -M_PI*cos(M_PI*y)*sin(M_PI*x);
    
    
        flux(0,0) *= -1.;
        flux(1,0) *= -1.;
    
}

void Hdiv2dPaper201504::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
{
    int dim = 2; //getDimension();
    
    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    flux.Resize(dim, 1);
    
    double x = pt[0];
    double y = pt[1];
    
    ff[0] = -2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    
}


void Hdiv2dPaper201504::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0] = 0.0;
    
}

void Hdiv2dPaper201504::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    solp[0] = 0.0;
    
}

void Hdiv2dPaper201504::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0] = 0.0;
    
}

void Hdiv2dPaper201504::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    
    solp[0] = 0.0;

}



TPZCompMesh *Hdiv2dPaper201504::CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    //    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExataH1, 5);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingH1,5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(fDim);
    cmesh->SetDefaultOrder(pOrder);
    
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;

    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;

    
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}


TPZCompMesh *Hdiv2dPaper201504::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,fDim);
    //TPZMatPoissonD3 *material = new TPZMatPoissonD3(fmatId,fDim);
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;

    
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(fDim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);

    
    cmesh->SetDefaultOrder(pOrder);
    
    
    TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(fmatskeleton, fDim-1, 1);
    TPZMaterial * mat2(matskelet);
    cmesh->InsertMaterialObject(mat2);
    
    

    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}

TPZCompMesh *Hdiv2dPaper201504::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
    //TPZMatPoissonD3 *material = new TPZMatPoissonD3(fmatId,fDim);
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(fDim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    cmesh->SetDefaultOrder(pOrder);
    
    bool h1function = true;// com esqueleto precisa disso
    if(pOrder>0/*h1function*/){
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    else{
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
        newnod.SetLagrangeMultiplier(1);
    }
    
    if(!h1function)
    {
        
        int nel = cmesh->NElements();
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            celdisc->SetConstC(1.);
            celdisc->SetCenterPoint(0, 0.);
            celdisc->SetCenterPoint(1, 0.);
            celdisc->SetCenterPoint(2, 0.);
            celdisc->SetTrueUseQsiEta();
            
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                if(ftriang==true) celdisc->SetTotalOrderShape();
                else celdisc->SetTensorialShape();
            }
            
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

TPZCompMesh *Hdiv2dPaper201504::CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
{
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim = gmesh->Dimension();
    
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,dim);
    
    //incluindo os dados do problema
    TPZFNMatrix<9,REAL> PermTensor(3,3,0.);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3,0.);
    
    
    // tensor de permutacao
    TPZFNMatrix<9,REAL> TP(3,3,0.0);
    TPZFNMatrix<9,REAL> InvTP(3,3,0.0);
    
    // Hard coded
    for (int id = 0; id < 3; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    PermTensor = TP;
    InvPermTensor = InvTP;
    
    material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    
    solexata = new TPZDummyFunction<STATE>(SolExata,5);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing,5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(10);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;

    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1D, 5);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2D, 5);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond2->SetForcingFunction(FBCond2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3D, 5);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond3->SetForcingFunction(FBCond3);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4D, 5);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    BCond4->SetForcingFunction(FBCond4);


    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();

    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    //Creating multiphysic elements containing skeletal elements.
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    
    int64_t nel = mphysics->ElementVec().NElements();
    
    std::map<int64_t, int64_t> bctoel, eltowrap;
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = mphysics->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if (matid < 0) {
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == dim && neighbour.Element()->Reference()) {
                    // got you!!
                    bctoel[el] = neighbour.Element()->Reference()->Index();
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                DebugStop();
            }
        }
    }
    
    TPZStack< TPZStack< TPZMultiphysicsElement *,7> > wrapEl;
    for(int64_t el = 0; el < nel; el++)
    {
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(mphysics->Element(el));
        if(mfcel->Dimension()==dim) TPZBuildMultiphysicsMesh::AddWrap(mfcel, fmatId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
    }
    
    for (int64_t el =0; el < wrapEl.size(); el++) {
        TPZCompEl *cel = wrapEl[el][0];
        int64_t index = cel->Index();
        eltowrap[index] = el;
    }
    
    meshvec[0]->CleanUpUnconnectedNodes();
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    std::map<int64_t, int64_t>::iterator it;
    for (it = bctoel.begin(); it != bctoel.end(); it++) {
        int64_t bcindex = it->first;
        int64_t elindex = it->second;
        if (eltowrap.find(elindex) == eltowrap.end()) {
            DebugStop();
        }
        int64_t wrapindex = eltowrap[elindex];
        TPZCompEl *bcel = mphysics->Element(bcindex);
        TPZMultiphysicsElement *bcmf = dynamic_cast<TPZMultiphysicsElement *>(bcel);
        if (!bcmf) {
            DebugStop();
        }
        wrapEl[wrapindex].Push(bcmf);
        
    }
    
    //------- Create and add group elements -------
    int64_t index, nenvel;
    nenvel = wrapEl.NElements();
    TPZStack<TPZElementGroup *> elgroups;
    for(int64_t ienv=0; ienv<nenvel; ienv++){
        TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
        elgroups.Push(elgr);
        nel = wrapEl[ienv].NElements();
        for(int jel=0; jel<nel; jel++){
            elgr->AddElement(wrapEl[ienv][jel]);
        }
    }
    
    mphysics->ComputeNodElCon();
    // create condensed elements
    // increase the NumElConnected of one pressure connects in order to prevent condensation
    for (int64_t ienv=0; ienv<nenvel; ienv++) {
        TPZElementGroup *elgr = elgroups[ienv];
        int nc = elgr->NConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = elgr->Connect(ic);
            if (c.LagrangeMultiplier() > 0) {
                c.IncrementElConnected();
                break;
            }
        }
        new TPZCondensedCompEl(elgr);
    }

    mphysics->CleanUpUnconnectedNodes();
    mphysics->ExpandSolution();

    return mphysics;
    
}


void Hdiv2dPaper201504::ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, int pos,  TPZFMatrix< REAL > &errors)
{
    
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExataH1, elerror, 0);
        
        int nerr = elerror.size();
        globalerrors.resize(nerr);

        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
    }
    // sequence of storage
    // for each
    errors.PutVal(pos, 0, p);      // polinomial order
    errors.PutVal(pos, 1, ndiv);   // refinement
    errors.PutVal(pos, 2, sqrt(globalerrors[1]));  // pressure
    errors.PutVal(pos, 3, sqrt(globalerrors[2]));  // flux
    
}

void Hdiv2dPaper201504::ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
{
    
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExataH1, elerror, 0);
        
        int nerr = elerror.size();
        globalerrors.resize(nerr);
        
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
    }
    
    out << ndiv << setw(10) << DoFT << setw(20) << DofCond << setw(28) << sqrt(globalerrors[1]) << setw(35)  << sqrt(globalerrors[2])  << endl;
    
}


void Hdiv2dPaper201504::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, int pos, TPZFMatrix< REAL > &errors)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globalerrorsDual(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrorsDual[i] += elerror[i]*elerror[i];
        }
        
        
    }
    
    
    nel = l2mesh->NElements();
    
    TPZManVector<REAL,10> globalerrorsPrimal(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        globalerrorsPrimal.resize(nerr);
        
        for (int i=0; i<nerr; i++) {
            globalerrorsPrimal[i] += elerror[i]*elerror[i];
        }
        
    }
    
    // sequence of storage
    // for each
    errors.PutVal(pos, 0, p);      // polinomial order
    errors.PutVal(pos, 1, ndiv);   // refinement
    errors.PutVal(pos, 2, sqrt(globalerrorsPrimal[1]));  // pressure
    errors.PutVal(pos, 3, sqrt(globalerrorsDual[1]));  // flux
    
}

void Hdiv2dPaper201504::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globalerrorsDual(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrorsDual[i] += elerror[i]*elerror[i];
        }
        
        
    }
    
    
    nel = l2mesh->NElements();

    TPZManVector<STATE,10> globalerrorsPrimal(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        globalerrorsPrimal.resize(nerr);
        
        for (int i=0; i<nerr; i++) {
            globalerrorsPrimal[i] += elerror[i]*elerror[i];
        }
        
    }
    
    out << ndiv << setw(10) << DoFT << setw(20) << DofCond << setw(28) << sqrt(globalerrorsPrimal[1]) << setw(35)  << sqrt(globalerrorsDual[1])  << endl;
    
}

void Hdiv2dPaper201504::ChangeExternalOrderConnects(TPZCompMesh *mesh){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    int cordermin = -1;
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        
        if(cel->Dimension()== dim){
            for (int icon=0; icon<ncon-1; icon++){
                TPZConnect &co  = cel->Connect(icon);
                corder = co.Order();
                nshape = co.NShape();
                if(corder!=cordermin){
                    cordermin = corder-1;
                    co.SetOrder(cordermin,cel->ConnectIndex(icon));
                    co.SetNShape(nshape-1);
                    mesh->Block().Set(co.SequenceNumber(),nshape-1);
                }
            }
        }
    }
    mesh->ExpandSolution();
    mesh->CleanUpUnconnectedNodes();
}



void Hdiv2dPaper201504::SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    std::cout <<"Numero de equacoes "<< fCmesh->NEquations()<< std::endl;
    
    bool isdirect = true;
    bool simetrico = true;
    bool isfrontal = false;
    if (isdirect)
    {
        if (simetrico)
        {
            //TPZSkylineStructMatrix strmat(fCmesh);
            if (isfrontal) {
                TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
                strmat.SetDecomposeType(ELDLt);
                
                int numthreads = 1;
                
                strmat.SetNumThreads(numthreads);
                
                an.SetStructuralMatrix(strmat);
            }
            else
            {
                //TPZBandStructMatrix full(fCmesh);
                TPZSkylineStructMatrix skylstr(fCmesh); //caso simetrico
                //    TPZSkylineNSymStructMatrix full(fCmesh);
                an.SetStructuralMatrix(skylstr);
            }
            
            
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt); //caso simetrico
            an.SetSolver(step);
            an.Run();
        }
        else
        {
            TPZBandStructMatrix full(fCmesh);
            an.SetStructuralMatrix(full);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELU);
            an.SetSolver(step);
            an.Run();
        }
        
    }
    else
    {
        TPZSkylineStructMatrix skylstr(fCmesh); //caso simetrico
        skylstr.SetNumThreads(10);
        an.SetStructuralMatrix(skylstr);
        
        TPZAutoPointer<TPZMatrix<STATE> > matbeingcopied = skylstr.Create();
        TPZAutoPointer<TPZMatrix<STATE> > matClone = matbeingcopied->Clone();
        
        TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>(matClone);
        TPZStepSolver<STATE> *Solver = new TPZStepSolver<STATE>(matbeingcopied);
        precond->SetReferenceMatrix(matbeingcopied);
        precond->SetDirect(ELDLt);
        Solver->SetGMRES(20, 20, *precond, 1.e-18, 0);
        //        Solver->SetCG(10, *precond, 1.0e-10, 0);
        an.SetSolver(*Solver);
        an.Run();
    }
    
    
    
}






