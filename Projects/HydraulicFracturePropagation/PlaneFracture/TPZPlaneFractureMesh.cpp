/**
 * \file
 * @brief Contains implementations of the TPZPlaneFractureMesh methods.
 * @author Cesar Lucci
 * @since 09/08/2010
 */

#include "TPZPlaneFractureMesh.h"

#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "TPZVTKGeoMesh.h"
#include "tpzgeoelrefpattern.h"
#include "tpzchangeel.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzreducedspace.h"
#include "TPZElast3Dnlinear.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZInterfaceEl.h"

#include "pzmat2dlin.h"
#include <algorithm>

using namespace pztopology;

#define just1IntersectionByEdge

TPZPlaneFractureMesh::TPZPlaneFractureMesh()
{
    std::cout << "Default constructor would not be used in this class!\n";
    DebugStop();
}
//------------------------------------------------------------------------------------------------------------


TPZPlaneFractureMesh::TPZPlaneFractureMesh(REAL bulletTVDIni, REAL bulletTVDFin,
                                           REAL xLength, REAL yLength, REAL Lmax,
                                           bool just1stripe, bool layersStripesToo)
{
    this->fInitialElIndex = 0;
    
    this->fPreservedMesh = new TPZGeoMesh;
    this->fRefinedMesh = NULL;
    this->fLastRefinedMesh = this->fRefinedMesh;
    
    this->fLmax = Lmax;
    this->fLfrac = 0.;
    
    this->fnstripes = 0;
    
    this->fCouplingMatVec.Resize(0);
    
    this->fjust1stripe = just1stripe;
    this->flayerStripesToo = layersStripesToo;
    
    std::set<REAL> espacamentoVerticalTVD;
    
    espacamentoVerticalTVD.insert(bulletTVDIni);
    espacamentoVerticalTVD.insert(bulletTVDFin);
    
    //Inserindo TVDs impostos (TVDs das camadas fornecidas)
    // Obs.: Eh considerado que o TVDfin da camada (s) corresponde ao TVDini da camada (s+1)
    int nstretches = globLayerStruct.NLayers();
    for(int s = 0; s < nstretches; s++)
    {
        REAL pos = globLayerStruct.GetLayer(s).fTVDini;
        espacamentoVerticalTVD.insert(pos);
        if(s == (nstretches-1))
        {
            pos = globLayerStruct.GetLayer(s).fTVDfin;
            espacamentoVerticalTVD.insert(pos);
        }
    }
    
    //Inserindo TVDs intermediarios para geracao da malha geometrica (baseados no Lmax)
    std::set<REAL>::iterator itS = espacamentoVerticalTVD.begin();
    REAL pos0 = *itS;
    itS++;
    for(; itS != espacamentoVerticalTVD.end(); itS++)
    {
        REAL pos1 = *itS;
        REAL deltaZ = fabs(pos1 - pos0);
        
        int nrows = 1;
        REAL deltaZused = deltaZ/nrows;
        while(deltaZused > Lmax)
        {
            nrows++;
            deltaZused = deltaZ/nrows;
        }
        
        for(int r = 1; r <= nrows; r++)
        {
            REAL z = pos0 + r*deltaZused;
            espacamentoVerticalTVD.insert(z);
        }
        pos0 = pos1;
    }
    
    
    //Conversao de TVD para profundidade (MD)
    std::list<REAL> espacamentoVerticalDEPTH;
    for(itS = espacamentoVerticalTVD.begin(); itS != espacamentoVerticalTVD.end(); itS++)
    {
        REAL posDepth = *itS;
        espacamentoVerticalDEPTH.push_back(-posDepth);
    }
    
    this->GeneratePreservedMesh(espacamentoVerticalDEPTH, bulletTVDIni, bulletTVDFin, xLength, yLength);
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureMesh::~TPZPlaneFractureMesh()
{
    delete this->fPreservedMesh;
    this->fPreservedMesh = NULL;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::InitializeFractureGeoMesh(TPZVec<std::pair<REAL,REAL> > &poligonalChain)
{
    std::cout << "\n************** GERANDO MALHA GEOMETRICA DA FRATURA\n";
    
    this->fLastRefinedMesh = this->fRefinedMesh;//ainda nao podemos deletar a ultima pois haverah transferencia de leakoff e matriz de massa
    this->fRefinedMesh = new TPZGeoMesh(*this->fPreservedMesh);
    
    int64_t nelem = fRefinedMesh->NElements();
    
    //Local structure data
    std::map< int64_t, std::set<REAL> > auxElIndex_TrimCoords;
    std::list< std::pair<int64_t,REAL> > auxElIndexSequence;
    
    //Searching for the edges that are intersected by poligonal chain
    this->DetectEdgesCrossed(poligonalChain, fRefinedMesh, auxElIndex_TrimCoords, auxElIndexSequence);
    
    //Refining auxiliar 1D elements
    TPZVec<TPZGeoEl*> sons;
    std::map< int64_t, std::set<REAL> >::iterator it;
    for(it = auxElIndex_TrimCoords.begin(); it != auxElIndex_TrimCoords.end(); it++)
    {
        int el1DIndex = it->first;
        TPZGeoEl * el1D = fRefinedMesh->ElementVec()[el1DIndex];
        
        TPZAutoPointer<TPZRefPattern> linRefp = Generate1DRefPatt(it->second);
        
        el1D->SetRefPattern(linRefp);
        el1D->Divide(sons);
    }
    
    //Refining 2D and 3D elements with the intention to match the geometry of the crack boundary
    for(int64_t el = 0; el < nelem; el++)
    {
        TPZGeoEl * gel = fRefinedMesh->ElementVec()[el];//2D element in 2D mesh
        if(globMaterialIdGen.IsOutsideFractMat(gel->MaterialId()) == false ||
           gel->HasSubElement())
        {
            continue;
        }
        
#ifdef PZDEBUG
        if(gel->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
        TPZAutoPointer<TPZRefPattern> elRefp = TPZRefPatternTools::PerfectMatchRefPattern(gel);
        if(elRefp)
        {
            gel->SetRefPattern(elRefp);
            gel->Divide(sons);
            
            int innerSide = gel->NSides() - 1;
            gel = gel->Neighbour(innerSide).Element();//3D element in 3D mesh
            
#ifdef PZDEBUG
            if(gel->Dimension() != 3)
            {
                DebugStop();
            }
#endif
            
            elRefp = TPZRefPatternTools::PerfectMatchRefPattern(gel);
            if(elRefp)
            {
                gel->SetRefPattern(elRefp);
                gel->Divide(sons);
            }
            int firstFace = 20;
            int lastFace = 25;
            for(int f = firstFace; f < lastFace; f++)//2D BC element in 3D mesh
            {
                TPZGeoElSide hexaFace(gel,f);
                TPZGeoElSide quadriFace = hexaFace.Neighbour();
                if( quadriFace != hexaFace && globMaterialIdGen.IsBoundaryMaterial(quadriFace.Element()->MaterialId()) )
                {
                    elRefp = TPZRefPatternTools::PerfectMatchRefPattern(quadriFace.Element());
                    if(elRefp)
                    {
                        quadriFace.Element()->SetRefPattern(elRefp);
                        quadriFace.Element()->Divide(sons);
                    }
                }
            }
        }
    }
    
    this->GenerateCrackBoundary(this->fRefinedMesh, auxElIndexSequence);
    this->SeparateElementsInMaterialSets(this->fRefinedMesh);
    this->UpdatePoligonalChain(this->fRefinedMesh, auxElIndexSequence, poligonalChain);
    
    this->RefineDirectionalToCrackTip(1);
    
//        {
//            std::ofstream outRefinedMesh("RefinedMesh.vtk");
//            TPZVTKGeoMesh::PrintGMeshVTK(fRefinedMesh, outRefinedMesh, true);
//            DebugStop();
//        }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::DeleteLastRefinedMesh()
{
    if(this->fLastRefinedMesh)
    {
        delete this->fLastRefinedMesh;
        this->fLastRefinedMesh = NULL;
    }
}
//------------------------------------------------------------------------------------------------------------

TPZCompMesh * TPZPlaneFractureMesh::GetFractureCompMesh()
{
    fRefinedMesh->ResetReference();
    
    TPZCompMesh * cmesh = new TPZCompMesh(this->fRefinedMesh);
    
    cmesh->SetDimModel(3);
    cmesh->SetDefaultOrder(1);
    
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int neumann = 1, dirich = 0, dirichDir = 3;
    
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        STATE young = globLayerStruct.GetLayer(lay).fYoung;
        STATE poisson = globLayerStruct.GetLayer(lay).fPoisson;
        TPZVec<STATE> force(3,0.);
        
        STATE prestressYY = globLayerStruct.GetLayer(lay).fSigYY;
        
        ////Rock
        TPZMaterial * materialLin = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force,
                                                        0., prestressYY, 0.);
        
        cmesh->InsertMaterialObject(materialLin);
        
        ////BCs
        k.Zero();
        f.Zero();
        f(0,0) = 1.;
        f(2,0) = 1.;
        TPZBndCond * dirichDirLeft = new TPZBndCond(materialLin,globMaterialIdGen.LeftMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(dirichDirLeft);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichRight = new TPZBndCond(materialLin,globMaterialIdGen.RightMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichRight);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichOutFracture = new TPZBndCond(materialLin,globMaterialIdGen.OutSideFractMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichOutFracture);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichFarfield = new TPZBndCond(materialLin,globMaterialIdGen.FarfieldMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichFarfield);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichTop = new TPZBndCond(materialLin,globMaterialIdGen.TopMatId(), dirich, k, f);
        cmesh->InsertMaterialObject(dirichTop);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichBottom = new TPZBndCond(materialLin,globMaterialIdGen.BottomMatId(), dirich, k, f);
        cmesh->InsertMaterialObject(dirichBottom);
        
        ///////////insideFract
        k.Zero();
        f.Zero();
        if(this->fjust1stripe)
        {
            int fractMatId = globMaterialIdGen.InsideFractMatId(lay, 0);
            if(this->fInsideFractureMatId.find(fractMatId) != this->fInsideFractureMatId.end())
            {
                TPZBndCond * neumannInsideFract = new TPZBndCond(materialLin,fractMatId, neumann, k, f);
                cmesh->InsertMaterialObject(neumannInsideFract);
            }
        }
        else
        {
            for(int stripe = 0; stripe < this->fnstripes; stripe++)
            {
                int fractMatId = globMaterialIdGen.InsideFractMatId(lay, stripe);
                if(this->fInsideFractureMatId.find(fractMatId) != this->fInsideFractureMatId.end())
                {
                    TPZBndCond * neumannInsideFract = new TPZBndCond(materialLin,fractMatId, neumann, k, f);
                    cmesh->InsertMaterialObject(neumannInsideFract);
                }
            }
        }
    }
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    this->IncreasePOrderOnFracture(cmesh,1);
    
    return cmesh;
}
//------------------------------------------------------------------------------------------------------------

//#define apenasPrimeiraLeva
void TPZPlaneFractureMesh::IncreasePOrderOnFracture(TPZCompMesh * cmesh, int oldPorder)
{
    cmesh->LoadReferences();
    
    //Capturando os elementos envolvidos diretamente
    std::set<int> porder2Indexes;
    for(int el = 0; el < cmesh->Reference()->NElements(); el++)
    {
        TPZGeoEl * elastGel = cmesh->Reference()->ElementVec()[el];
        
        if(globMaterialIdGen.IsInsideFractMat(elastGel->MaterialId()))
        {
            TPZGeoEl * lowestFather = elastGel->LowestFather();
            porder2Indexes.insert(lowestFather->Index());
            
            for(int sd = 0; sd < lowestFather->NSides(); sd++)
            {
                TPZGeoElSide fractSide(lowestFather,sd);
                TPZGeoElSide neighSide(fractSide.Neighbour());
                
                while(neighSide != fractSide)
                {
                    TPZGeoEl * lowestFather2 = neighSide.Element()->LowestFather();
                    porder2Indexes.insert(lowestFather2->Index());
                    
                    neighSide = neighSide.Neighbour();
                }
            }
        }
    }
    
    //Capturando os elementos envolvidos indiretamente
    std::set<int> porder2IndexesTemp;
    
#ifdef apenasPrimeiraLeva
    //Aqui serao tornados p2 somente os que faltaram na primeira camada em y que encosta na fratura
    //(nesta etapa, serao pegos os vizinhos ao cracktip que ficaram de fora na caçada acima)
    //DO NOTHING
#else
    //Aqui serao tornados p2 os da segunda camada em y que encosta na fratura
    for(int el = 0; el < cmesh->Reference()->NElements(); el++)
    {
        if(porder2Indexes.find(el) == porder2Indexes.end())
        {
            continue;
        }
        
        TPZGeoEl * elastGel = cmesh->Reference()->ElementVec()[el];
        
        for(int sd = 0; sd < elastGel->NSides(); sd++)
        {
            TPZGeoElSide fractSide(elastGel,sd);
            TPZGeoElSide neighSide(fractSide.Neighbour());
            
            while(neighSide != fractSide)
            {
                TPZGeoEl * lowestFather = neighSide.Element()->LowestFather();
                porder2IndexesTemp.insert(lowestFather->Index());
                
                neighSide = neighSide.Neighbour();
            }
        }
    }
#endif
    
    porder2Indexes.insert(porder2IndexesTemp.begin(),porder2IndexesTemp.end());
    
    //Subindo pOrder na regiao da fratura
    int greatherPOrder = oldPorder + 1;
    
    //    TPZVec<REAL> elData(cmesh->NElements(),1);//VTK
    std::set<int>::iterator it;
    for(it = porder2Indexes.begin(); it != porder2Indexes.end(); it++)
    {
        TPZGeoEl * elastGel = cmesh->Reference()->ElementVec()[*it];
        if(elastGel->HasSubElement() == false && elastGel->Reference())
        {
            TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace*>(elastGel->Reference());
            intel->PRefine(greatherPOrder);
            //            elData[elastGel->Reference()->Index()] = greatherPOrder;//VTK
        }
        else
        {
            TPZVec<TPZGeoEl *> sons(0);
            elastGel->GetHigherSubElements(sons);
            
            for(int ss = 0; ss < sons.NElements(); ss++)
            {
                TPZGeoEl * sonGel = sons[ss];
                if(sonGel->Reference())
                {
                    TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace*>(sonGel->Reference());
                    intel->PRefine(greatherPOrder);
                    //                    elData[sonGel->Reference()->Index()] = greatherPOrder;//VTK
                }
            }
        }
    }
    //    std::ofstream outPOrder("CMeshPOrder.vtk");
    //    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, outPOrder, elData, "porder");
}

TPZCompMeshReferred * TPZPlaneFractureMesh::GetFractureCompMeshReferred(TPZCompMesh * cmeshRef)
{
    fRefinedMesh->ResetReference();
    
    TPZCompMeshReferred * cmesh = new TPZCompMeshReferred(fRefinedMesh);
    
    cmesh->SetDimModel(3);
    cmesh->SetDefaultOrder(2);
    
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int neumann = 1, dirich = 0, dirichDir = 3;
    
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        STATE young = globLayerStruct.GetLayer(lay).fYoung;
        STATE poisson = globLayerStruct.GetLayer(lay).fPoisson;
        TPZVec<STATE> force(3,0.);
        
        STATE prestressYY = 0.;//A integral-J nao inclui a translacao do pre-stress (ver LayerStruct::GetStressAppliedJustForJIntegral)
        
        ////Rock
        TPZMaterial * materialLin = new TPZElasticity3D(globMaterialIdGen.RockMatId(lay), young, poisson, force, 0., prestressYY, 0.);
        cmesh->InsertMaterialObject(materialLin);
        
        ////BCs
        k.Zero();
        f.Zero();
        f(0,0) = 1.;
        f(2,0) = 1.;
        TPZBndCond * dirichDirLeft = new TPZBndCond(materialLin,globMaterialIdGen.LeftMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(dirichDirLeft);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichRight = new TPZBndCond(materialLin,globMaterialIdGen.RightMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichRight);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichOutFracture = new TPZBndCond(materialLin,globMaterialIdGen.OutSideFractMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichOutFracture);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichFarfield = new TPZBndCond(materialLin,globMaterialIdGen.FarfieldMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichFarfield);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichTop = new TPZBndCond(materialLin,globMaterialIdGen.TopMatId(), dirich, k, f);
        cmesh->InsertMaterialObject(dirichTop);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichBottom = new TPZBndCond(materialLin,globMaterialIdGen.BottomMatId(), dirich, k, f);
        cmesh->InsertMaterialObject(dirichBottom);
        
        ///////////insideFract
        k.Zero();
        f.Zero();
        if(this->fjust1stripe)
        {
            int fractMatId = globMaterialIdGen.InsideFractMatId(lay, 0);
            if(this->fInsideFractureMatId.find(fractMatId) != this->fInsideFractureMatId.end())
            {
                TPZBndCond * neumannInsideFract = new TPZBndCond(materialLin,fractMatId, neumann, k, f);
                cmesh->InsertMaterialObject(neumannInsideFract);
            }
        }
        else
        {
            for(int stripe = 0; stripe < this->fnstripes; stripe++)
            {
                int fractMatId = globMaterialIdGen.InsideFractMatId(lay, stripe);
                if(this->fInsideFractureMatId.find(fractMatId) != this->fInsideFractureMatId.end())
                {
                    TPZBndCond * neumannInsideFract = new TPZBndCond(materialLin,fractMatId, neumann, k, f);
                    cmesh->InsertMaterialObject(neumannInsideFract);
                }
            }
        }
    }
    
    int numsol = cmeshRef->Solution().Cols();
    cmesh->AllocateNewConnect(numsol, 1, 1);
    
    TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmesh);
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->LoadReferred(cmeshRef);
    
    return cmesh;
}
//------------------------------------------------------------------------------------------------------------

TPZCompMesh * TPZPlaneFractureMesh::GetPressureCompMesh()
{
    fRefinedMesh->ResetReference();
    
    TPZCompMesh * cmesh = new TPZCompMesh(fRefinedMesh);
    
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(2);
    
    TPZFMatrix<STATE> xk(1,1,1.);
    TPZFMatrix<STATE> xc(1,1,0.);
    TPZFMatrix<REAL> xb(1,1,0.);
    TPZFMatrix<STATE> xf(1,1,-2.);
    int neumannFluxIn = 4;
    
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        ///////////insideFract
        if(this->fjust1stripe)
        {
            int fractMatId = globMaterialIdGen.InsideFractMatId(lay, 0);
            if(this->fInsideFractureMatId.find(fractMatId) != this->fInsideFractureMatId.end())
            {
                TPZMat2dLin * mat = new TPZMat2dLin(fractMatId);
                mat->SetMaterial(xk,xc,xf);
                cmesh->InsertMaterialObject(mat);
                
                ///////////bullet
                TPZFMatrix<STATE> k(2,2,0.), f(2,1,0.);
                TPZBndCond * fluxInBC = new TPZBndCond(mat, globMaterialIdGen.BulletMatId(lay), neumannFluxIn, k, f);
                cmesh->InsertMaterialObject(fluxInBC);
            }
        }
        else
        {
            for(int stripe = 0; stripe < fnstripes; stripe++)
            {
                ///////////insideFract
                int fractMatId = globMaterialIdGen.InsideFractMatId(lay, stripe);
                if(this->fInsideFractureMatId.find(fractMatId) != this->fInsideFractureMatId.end())
                {
                    TPZMat2dLin * mat = new TPZMat2dLin(fractMatId);
                    mat->SetMaterial(xk,xc,xf);
                    cmesh->InsertMaterialObject(mat);
                    
                    ///////////bullet
                    if(stripe == 0)
                    {
                        TPZFMatrix<STATE> k(2,2,0.), f(2,1,0.);
                        TPZBndCond * fluxInBC = new TPZBndCond(mat, globMaterialIdGen.BulletMatId(lay), neumannFluxIn, k, f);
                        cmesh->InsertMaterialObject(fluxInBC);
                    }
                }
            }
        }
    }
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    /*{//Elementos de Interface : PIRA
     this->fRefinedMesh->ResetReference();
     
     int64_t index;
     for(int el = 0; el < this->fRefinedMesh->NElements(); el++)
     {
     TPZGeoEl * gel = this->fRefinedMesh->ElementVec()[el];
     if(!gel || gel->HasSubElement() || globMaterialIdGen.IsBulletMaterial(gel->MaterialId()) == false)
     {
     continue;
     }
     if(gel->Dimension() != 1)
     {
     DebugStop();
     }
     {
     REAL Diameter = 0.3;//por enquanto hardcode
     REAL visc = 200.02E-10;
     REAL Qinj = -0.05333333333333/2/20.;//
     TPZPlaneFractBulletMat * interfMat = new TPZPlaneFractBulletMat(gel->MaterialId(),
     Diameter,
     visc,
     Qinj);
     cmesh->InsertMaterialObject(interfMat);
     }
     cmesh->CreateCompEl(gel, index);
     }
     
     cmesh->LoadReferences();
     
     int nel = this->fRefinedMesh->NElements();
     for(int el = 0; el < nel; el++)
     {
     TPZGeoEl * gel = this->fRefinedMesh->ElementVec()[el];
     if(!gel || gel->HasSubElement() || globMaterialIdGen.IsBulletMaterial(gel->MaterialId()) == false)
     {
     continue;
     }
     
     TPZGeoElSide gelRight1Dside(gel,2);//right
     TPZCompElSide celRight1Dside = gelRight1Dside.Reference();
     
     bool found = false;
     TPZGeoElSide gelLeft2Dside = gelRight1Dside.Neighbour();//left
     while(gelLeft2Dside != gelRight1Dside)
     {
     int matId2D = gelLeft2Dside.Element()->MaterialId();
     if(globMaterialIdGen.IsInsideFractMat(matId2D) &&
     gelLeft2Dside.Element()->Reference())
     {
     found = true;
     }
     else
     {
     gelLeft2Dside = gelLeft2Dside.Neighbour();
     }
     if(found)
     {
     break;
     }
     }
     if(gelLeft2Dside.Element()->Dimension() != 2 || !gelLeft2Dside.Element()->Reference())
     {
     std::cout << "\n\n\nDimension != 2 or No reference (no compel)!!!\n\n\n";
     DebugStop();
     }
     
     TPZCompElSide celLeft2Dside = gelLeft2Dside.Reference();
     
     TPZGeoEl * face = gel->CreateBCGeoEl(2, gel->MaterialId());
     
     new TPZInterfaceElement(*cmesh, face, index, celLeft2Dside, celRight1Dside);
     }
     }*/
    
    return cmesh;
}
//------------------------------------------------------------------------------------------------------------

TPZCompMesh * TPZPlaneFractureMesh::GetMultiPhysicsCompMesh(TPZVec<TPZCompMesh *> & meshvec, TPZCompMesh * lastElastCMesh,
                                                            REAL Qinj1wing_Hbullet, REAL visc)
{
    fRefinedMesh->ResetReference();
    
    TPZCompMesh * cmesh = new TPZCompMesh(this->fRefinedMesh);
    
    cmesh->SetDimModel(3);
    cmesh->SetDefaultOrder(2);
    
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    
    int neumann = 1, dirich = 0, dirichDir = 3, neumannFluxIn = 4;
    
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        STATE young = globLayerStruct.GetLayer(lay).fYoung;
        STATE poisson = globLayerStruct.GetLayer(lay).fPoisson;
        TPZVec<STATE> force(3,0.);
        
        STATE prestressYY = globLayerStruct.GetLayer(lay).fSigYY;
        
        REAL Cl = globLayerStruct.GetLayer(lay).fCl;
        REAL Pe = globLayerStruct.GetLayer(lay).fPe;
        REAL gradPref = globLayerStruct.GetLayer(lay).fgradPref;
        REAL vsp = globLayerStruct.GetLayer(lay).fvsp;
        
        ////Rock
        TPZPlaneFractCouplingMat * couplingMat = new TPZPlaneFractCouplingMat(globMaterialIdGen.RockMatId(lay), young, poisson, force,
                                                                              0., prestressYY, 0.,
                                                                              visc,
                                                                              Cl,
                                                                              Pe,
                                                                              gradPref,
                                                                              vsp);
        couplingMat->SetLastElastCMesh(lastElastCMesh);
        
        //Inserindo no vetor de materiais de acoplamento
        int oldSz = fCouplingMatVec.NElements();
        fCouplingMatVec.Resize(oldSz+1);
        fCouplingMatVec[oldSz] = couplingMat;
        //
        
        cmesh->InsertMaterialObject(couplingMat);
        
        ////BCs
        k.Zero();
        f.Zero();
        f(0,0) = 1.;
        f(2,0) = 1.;
        TPZBndCond * dirichDirLeft = new TPZBndCond(couplingMat,globMaterialIdGen.LeftMatId(lay), dirichDir, k, f);
        cmesh->InsertMaterialObject(dirichDirLeft);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichRight = new TPZBndCond(couplingMat,globMaterialIdGen.RightMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichRight);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichOutFracture = new TPZBndCond(couplingMat,globMaterialIdGen.OutSideFractMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichOutFracture);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichFarfield = new TPZBndCond(couplingMat,globMaterialIdGen.FarfieldMatId(lay), dirich, k, f);
        cmesh->InsertMaterialObject(dirichFarfield);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichTop = new TPZBndCond(couplingMat,globMaterialIdGen.TopMatId(), dirich, k, f);
        cmesh->InsertMaterialObject(dirichTop);
        
        k.Zero();
        f.Zero();
        TPZBndCond * dirichBottom = new TPZBndCond(couplingMat,globMaterialIdGen.BottomMatId(), dirich, k, f);
        cmesh->InsertMaterialObject(dirichBottom);
        
        ///////////insideFract
        if(this->fjust1stripe)
        {
            ///////////insideFract
            int fractMatId = globMaterialIdGen.InsideFractMatId(lay, 0);
            if(this->fInsideFractureMatId.find(fractMatId) != this->fInsideFractureMatId.end())
            {
                k.Zero();
                f.Zero();
                TPZBndCond * neumannInsideFract = new TPZBndCond(couplingMat,fractMatId, neumann, k, f);
                cmesh->InsertMaterialObject(neumannInsideFract);
                
                ///////////bullet
                k.Zero();
                f.Zero();
                f(0,0) = Qinj1wing_Hbullet;
                TPZBndCond * fluxInBC = new TPZBndCond(couplingMat, globMaterialIdGen.BulletMatId(lay), neumannFluxIn, k, f);
                cmesh->InsertMaterialObject(fluxInBC);
            }
        }
        else
        {
            for(int stripe = 0; stripe < this->fnstripes; stripe++)
            {
                ///////////insideFract
                int fractMatId = globMaterialIdGen.InsideFractMatId(lay, stripe);
                if(this->fInsideFractureMatId.find(fractMatId) != this->fInsideFractureMatId.end())
                {
                    k.Zero();
                    f.Zero();
                    TPZBndCond * neumannInsideFract = new TPZBndCond(couplingMat,fractMatId, neumann, k, f);
                    cmesh->InsertMaterialObject(neumannInsideFract);
                    
                    ///////////bullet
                    if(stripe == 0)
                    {
                        k.Zero();
                        f.Zero();
                        f(0,0) = Qinj1wing_Hbullet;
                        TPZBndCond * fluxInBC = new TPZBndCond(couplingMat, globMaterialIdGen.BulletMatId(lay), neumannFluxIn, k, f);
                        cmesh->InsertMaterialObject(fluxInBC);
                    }
                }
            }
        }
    }
    
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec,cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec,cmesh);
    
    return cmesh;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::SetNeumannOnEntireFracture(TPZCompMesh * cmeshref)
{
    bool neumanWasApplied = false;
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        if(this->fjust1stripe)
        {
            int matId = globMaterialIdGen.InsideFractMatId(lay, 0);
            if(cmeshref->MaterialVec().find(matId) == cmeshref->MaterialVec().end())
            {//Fratura ainda nao entrou nesta camada!
                continue;
            }
            
            TPZMaterial * mat = cmeshref->MaterialVec().find(matId)->second;
            if(!mat)
            {
                DebugStop();
            }
            TPZBndCond * bcmat = dynamic_cast<TPZBndCond *>(mat);
            if(!bcmat)
            {
                DebugStop();
            }
            
            neumanWasApplied = true;
            bcmat->Val2()(1,0) = globLayerStruct.StressAppliedOnFractureStripe();
        }
        else
        {
            for(int stripe = 0; stripe < this->fnstripes; stripe++)
            {
                int matId = globMaterialIdGen.InsideFractMatId(lay, stripe);
                if(cmeshref->MaterialVec().find(matId) == cmeshref->MaterialVec().end())
                {//Fratura ainda nao entrou nesta camada!
                    continue;
                }
                
                TPZMaterial * mat = cmeshref->MaterialVec().find(matId)->second;
                if(!mat)
                {
                    DebugStop();
                }
                TPZBndCond * bcmat = dynamic_cast<TPZBndCond *>(mat);
                if(!bcmat)
                {
                    DebugStop();
                }
                
                neumanWasApplied = true;
                bcmat->Val2()(1,0) = globLayerStruct.StressAppliedOnFractureStripe();
            }
        }
    }
    
    if(neumanWasApplied == false)
    {
        std::cout << "\n\n\nNeumann was not applied on entire fracture\n";
        DebugStop();
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::SetNeumannOnThisStripe(TPZCompMesh * cmeshref, int actStripe)
{
    if(this->fjust1stripe)
    {
        DebugStop();//Who is calling this method if is just 1 stripe???
    }
    
    bool neumanWasApplied = false;
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        for(int stripe = 0; stripe < this->fnstripes; stripe++)
        {
            int matId = globMaterialIdGen.InsideFractMatId(lay, stripe);
            if(cmeshref->MaterialVec().find(matId) == cmeshref->MaterialVec().end())
            {//Fratura ainda nao entrou nesta camada!
                continue;
            }
            
            TPZMaterial * mat = cmeshref->MaterialVec().find(matId)->second;
            if(!mat)
            {
                DebugStop();
            }
            TPZBndCond * bcmat = dynamic_cast<TPZBndCond *>(mat);
            if(!bcmat)
            {
                DebugStop();
            }
            if(stripe == actStripe)
            {
                neumanWasApplied = true;
                
                bcmat->Val2()(1,0) = globLayerStruct.StressAppliedOnFractureStripe();
            }
            else
            {
                bcmat->Val2()(1,0) = 0.;
            }
        }
    }
    
    if(neumanWasApplied == false)
    {
        std::cout << "\n\n\nNeumann was not applied on stripe " << actStripe << "\n";
        DebugStop();
    }
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureMesh::SetNeumannOnThisLayerAndStripe(TPZCompMesh * cmeshref, int actLayer, int actStripe)
{
    if(this->flayerStripesToo == false)
    {
        DebugStop();//Who is calling this method if is no LayerStripes???
    }
    
    bool neumanWasApplied = false;
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        for(int stripe = 0; stripe < this->fnstripes; stripe++)
        {
            int matId = globMaterialIdGen.InsideFractMatId(lay, stripe);
            if(cmeshref->MaterialVec().find(matId) == cmeshref->MaterialVec().end())
            {//Fratura ainda nao entrou nesta camada!
                continue;
            }
            
            TPZMaterial * mat = cmeshref->MaterialVec().find(matId)->second;
            if(!mat)
            {
                DebugStop();
            }
            TPZBndCond * bcmat = dynamic_cast<TPZBndCond *>(mat);
            if(!bcmat)
            {
                DebugStop();
            }
            if(lay == actLayer && stripe == actStripe)
            {
                neumanWasApplied = true;
                
                bcmat->Val2()(1,0) = globLayerStruct.StressAppliedOnFractureStripe();
            }
            else
            {
                bcmat->Val2()(1,0) = 0.;
            }
        }
    }
    
    return neumanWasApplied;
}
//------------------------------------------------------------------------------------------------------------

int TPZPlaneFractureMesh::NStripes()
{
    return this->fnstripes;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::SetActualState()
{
    for(int coupmat = 0; coupmat < this->fCouplingMatVec.NElements(); coupmat++)
    {
        this->fCouplingMatVec[coupmat]->SetActualState();
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::SetPastState()
{
    for(int coupmat = 0; coupmat < this->fCouplingMatVec.NElements(); coupmat++)
    {
        this->fCouplingMatVec[coupmat]->SetPastState();
    }
}
//------------------------------------------------------------------------------------------------------------

int TPZPlaneFractureMesh::NCrackTipElements()
{
    return this->fcrackBoundaryElementsIndexes.NElements();
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl * TPZPlaneFractureMesh::GetCrackTipGeoElement(int pos)
{
    int elIndex = this->fcrackBoundaryElementsIndexes[pos];
    TPZGeoEl * gel1D = this->fRefinedMesh->ElementVec()[elIndex];
    
    return gel1D;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::GetKIc_FromLayerOfThisZcoord(REAL zCoord, REAL &KIc)
{
    KIc = globLayerStruct.GetLayerFromZcoord(zCoord).fKIc;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureMesh::GeoElementIsOnPreservedMesh(TPZGeoEl * gel)
{
    if(gel->Index() >= this->fPreservedMesh->NElements())
    {
        return false;
    }
    TPZGeoEl * preservedGel = this->fPreservedMesh->ElementVec()[gel->Index()];
    bool gelIsOnPreservedMesh = (preservedGel->Id() == gel->Id());
    
    return gelIsOnPreservedMesh;
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl * TPZPlaneFractureMesh::Find2DElementNearCrackTip(int posCrackTip, TPZVec<REAL> & x)
{
    int64_t crackElIndex = this->fcrackBoundaryElementsIndexes[posCrackTip];
    
    TPZVec<REAL> qsi(2,0.);
    TPZGeoEl * gel = this->fRefinedMesh->FindElementCaju(x, qsi, crackElIndex, 2);
    
    return gel;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::GeneratePreservedMesh(std::list<REAL> & espacamentoVerticalDEPTH,
                                                 REAL bulletTVDIni, REAL bulletTVDFin,
                                                 REAL xLength, REAL yLength)
{
    std::cout << "\n************** GERANDO MALHA GEOMETRICA DE REFERENCIA (PRESERVED)\n\n\n";
    
    REAL bulletDEPTHIni = -bulletTVDIni;
    REAL bulletDEPTHFin = -bulletTVDFin;
    
    TPZVec< TPZVec<REAL> > NodeCoord(0);
    int64_t nrows, ncols;
    
    int nDirRef = int(log(fabs(yLength)/this->fLmax)/log(2.));
    int nLayersY = nDirRef + 2;
    
    REAL Y = 0.;
    GenerateNodesAtPlaneY(espacamentoVerticalDEPTH, xLength, NodeCoord, nrows, ncols, Y);
    
    Y = yLength/pow(2.,nDirRef);
    for(int lay = 1; lay < nLayersY; lay++)
    {
        GenerateNodesAtPlaneY(espacamentoVerticalDEPTH, xLength, NodeCoord, nrows, ncols, Y);
        Y *= 2.;
    }
    int64_t nNodesByLayer = nrows*ncols;
    int64_t Qnodes = nNodesByLayer * nLayersY;
    
    //initializing gmesh
    this->fPreservedMesh->NodeVec().Resize(Qnodes);
    
    int64_t pos = 0;
    TPZGeoNode Node;
    for(int64_t l = 0; l < nLayersY; l++)
    {
        for(int64_t n = 0; n < nNodesByLayer; n++)
        {
            Node.SetNodeId(pos);
            Node.SetCoord(NodeCoord[pos]);
            this->fPreservedMesh->NodeVec()[pos] = Node;
            pos++;
        }
    }
    
    //inserting quadrilaterals and 1D bullet region
    TPZVec<int64_t> Topol2(2), Topol4(4), Topol8(8);
    for(int64_t r = 0; r < (nrows-1); r++)
    {
        for(int64_t c = 0; c < (ncols-1); c++)
        {
            Topol4[0] = ncols*r+c; Topol4[1] = ncols*r+c+1; Topol4[2] = ncols*(r+1)+c+1; Topol4[3] = ncols*(r+1)+c;
            
            REAL z0 = this->fPreservedMesh->NodeVec()[Topol4[0]].Coord(2);
            REAL z3 = this->fPreservedMesh->NodeVec()[Topol4[3]].Coord(2);
            REAL zMed = (z0 + z3)/2.;
            
            int whatLayer = globLayerStruct.WhatLayer(zMed);
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.OutSideFractMatId(whatLayer),*this->fPreservedMesh);
            
            //detecting bullet region
            if(c == 0)
            {
                REAL tol = 1.E-3;
                if(zMed < (bulletDEPTHIni + tol) && zMed > (bulletDEPTHFin - tol))
                {
                    Topol2[0] = ncols*r;
                    Topol2[1] = ncols*(r+1);
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol2,globMaterialIdGen.BulletMatId(whatLayer),*this->fPreservedMesh);
                }
            }
        }
    }
    
    //inserting hexaedrons
    for(int64_t ly = 0; ly < (nLayersY-1); ly++)
    {
        for(int64_t r = 0; r < (nrows-1); r++)
        {
            for(int64_t c = 0; c < (ncols-1); c++)
            {
                Topol8[0] = ncols*r+c + ly*nNodesByLayer;
                Topol8[1] = ncols*r+c+1 + ly*nNodesByLayer;
                Topol8[2] = ncols*(r+1)+c+1 + ly*nNodesByLayer;
                Topol8[3] = ncols*(r+1)+c + ly*nNodesByLayer;
                //
                Topol8[4] = ncols*r+c + (ly+1)*nNodesByLayer;
                Topol8[5] = ncols*r+c+1 + (ly+1)*nNodesByLayer;
                Topol8[6] = ncols*(r+1)+c+1 + (ly+1)*nNodesByLayer;
                Topol8[7] = ncols*(r+1)+c + (ly+1)*nNodesByLayer;
                
                REAL z0 = this->fPreservedMesh->NodeVec()[Topol8[0]].Coord(2);
                REAL z3 = this->fPreservedMesh->NodeVec()[Topol8[3]].Coord(2);
                REAL zMed = (z0+z3)/2.;
                
                int whatLayer = globLayerStruct.WhatLayer(zMed);
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (Topol8,globMaterialIdGen.RockMatId(whatLayer),*this->fPreservedMesh);
                
                if(ly == (nLayersY-2))//farfield cc
                {
                    Topol4[0] = ncols*r+c + (ly+1)*nNodesByLayer;
                    Topol4[1] = ncols*r+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[2] = ncols*(r+1)+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*(r+1)+c + (ly+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.FarfieldMatId(whatLayer),*this->fPreservedMesh);
                }
                if(c == 0)//left cc
                {
                    Topol4[0] = ncols*r+c + ly*nNodesByLayer;
                    Topol4[1] = ncols*r+c + (ly+1)*nNodesByLayer;
                    Topol4[2] = ncols*(r+1)+c + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*(r+1)+c + ly*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.LeftMatId(whatLayer),*this->fPreservedMesh);
                }
                else if(c == (ncols - 2))//right cc
                {
                    Topol4[0] = ncols*r+c+1 + ly*nNodesByLayer;
                    Topol4[1] = ncols*(r+1)+c+1 + ly*nNodesByLayer;
                    Topol4[2] = ncols*(r+1)+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*r+c+1 + (ly+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.RightMatId(whatLayer),*this->fPreservedMesh);
                }
                if(r == 0)//top cc
                {
                    Topol4[0] = ncols*r+c + ly*nNodesByLayer;
                    Topol4[1] = ncols*r+c+1 + ly*nNodesByLayer;
                    Topol4[2] = ncols*r+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*r+c + (ly+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.TopMatId(),*this->fPreservedMesh);
                }
                else if(r == (nrows - 2))//bottom cc
                {
                    Topol4[0] = ncols*(r+1)+c + ly*nNodesByLayer;
                    Topol4[1] = ncols*(r+1)+c + (ly+1)*nNodesByLayer;
                    Topol4[2] = ncols*(r+1)+c+1 + (ly+1)*nNodesByLayer;
                    Topol4[3] = ncols*(r+1)+c+1 + ly*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (Topol4,globMaterialIdGen.BottomMatId(),*this->fPreservedMesh);
                }
            }
        }
    }
    
    this->fPreservedMesh->BuildConnectivity();
    this->fPreservedMesh->SetMaxElementId(this->fPreservedMesh->NElements()-1);
    this->fPreservedMesh->SetMaxNodeId(this->fPreservedMesh->NNodes()-1);
    
    //    this->RefineUniformAllFracturePlane(1);//AQUICAJU
    
    //    std::ofstream outPreservedMesh("PreservedMesh.vtk");
    //    TPZVTKGeoMesh::PrintGMeshVTK(this->fPreservedMesh, outPreservedMesh, true);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::RefineUniformAllFracturePlane(int ndiv)
{
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ECube);
    
    for(int div = 0; div < ndiv; div++)
    {
        int nelem = this->fPreservedMesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            TPZGeoEl * gel = this->fPreservedMesh->ElementVec()[el];
            if(gel->HasSubElement())
            {
                continue;
            }
            if(globMaterialIdGen.IsOutsideFractMat(gel->MaterialId()))
            {
                gel->SetRefPattern(NULL);//garantia de que serah utilizado refinamento uniforme!!!
                TPZVec<TPZGeoEl*> sons(0);
                gel->Divide(sons);
                
                TPZGeoEl * rock3D = gel->Neighbour(gel->NSides()-1).Element();
                sons.Resize(0);
                rock3D->SetRefPattern(NULL);
                rock3D->Divide(sons);
                
                {//Procurando por faces de contorno para refinar tambem!
                    for(int edge = 4; edge < 8; edge++)//exclusivo para quadrilatero!
                    {
                        TPZGeoElSide gelSide(gel, edge);
                        TPZGeoElSide boundarySide = gelSide.Neighbour();
                        while(boundarySide != gelSide)
                        {
                            if(globMaterialIdGen.IsBoundaryMaterial(boundarySide.Element()->MaterialId()) ||
                               globMaterialIdGen.IsBulletMaterial(boundarySide.Element()->MaterialId()))
                            {
                                sons.Resize(0);
                                boundarySide.Element()->SetRefPattern(NULL);
                                boundarySide.Element()->Divide(sons);
                            }
                            boundarySide = boundarySide.Neighbour();
                        }
                    }
                }
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::RefineDirectionalToCrackTip(int ndiv)
{
    std::set<int> crackTipMat;
    crackTipMat.insert(globMaterialIdGen.CrackTipMatId());
    for(int div = 0; div < ndiv; div++)
    {
        int nelem = this->fRefinedMesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            TPZGeoEl * gel = this->fRefinedMesh->ElementVec()[el];
            if(gel->HasSubElement() == false)
            {
                TPZRefPatternTools::RefineDirectional(gel, crackTipMat);
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::TurnIntoQuarterPoint(TPZGeoMesh * refinedMesh)
{
    for(int i = 0; i < fcrackBoundaryElementsIndexes.NElements(); i++)
    {
        TPZGeoEl * gel1D = refinedMesh->ElementVec()[fcrackBoundaryElementsIndexes[i]];
        
#ifdef PZDEBUG
        if(gel1D->Dimension() != 1)
        {
            DebugStop();
        }
#endif
        
        for(int s = 0; s < 2; s++)
        {
            TPZGeoElSide edge(gel1D,s);
            TPZGeoElSide neigh(edge.Neighbour());
            
            while(edge != neigh)
            {
                if(neigh.Element()->HasSubElement() == false && neigh.Element()->IsLinearMapping())
                {
                    int neighSide = neigh.Side();
                    TPZGeoEl * neighEl = TPZChangeEl::ChangeToQuarterPoint(refinedMesh, neigh.Element()->Index(), neighSide);
                    neigh = neighEl->Neighbour(neighSide);
                }
                else
                {
                    neigh = neigh.Neighbour();
                }
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::GenerateNodesAtPlaneY(std::list<REAL> & espacamentoVerticalDEPTH, REAL xLength,
                                                 TPZVec< TPZVec<REAL> > & NodeCoord, int64_t & nrows, int64_t & ncols,
                                                 REAL Y)
{
    nrows = espacamentoVerticalDEPTH.size();
    
    int nStretches = 1;
    REAL deltaX = xLength/nStretches;
    while(deltaX > this->fLmax)
    {
        nStretches++;
        deltaX = xLength/nStretches;
    }
    ncols = nStretches+1;
    
    const int64_t nNodesByLayer = nrows*ncols;
    int64_t oldNNodes = NodeCoord.NElements();
    int64_t newNNodes = oldNNodes + nNodesByLayer;
    
    NodeCoord.Resize(newNNodes);
    
    int64_t nodeId = oldNNodes;
    
    //setting nodes coords
    std::list<REAL>::iterator it;
    for(it = espacamentoVerticalDEPTH.begin(); it != espacamentoVerticalDEPTH.end(); it++)
    {
        for(int64_t c = 0; c < ncols; c++)
        {
            NodeCoord[nodeId].Resize(3,0.);
            NodeCoord[nodeId][0] = c*deltaX;
            NodeCoord[nodeId][1] = Y;
            NodeCoord[nodeId][2] = *it;
            nodeId++;
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::DetectEdgesCrossed(TPZVec<std::pair<REAL,REAL> > &poligonalChain, TPZGeoMesh * planeMesh,
                                              std::map< int64_t, std::set<REAL> > &auxElIndex_TrimCoords,
                                              std::list< std::pair<int64_t,REAL> > &auxElIndexSequence)
{
    int64_t npoints = (poligonalChain.NElements());
    int64_t nelem = planeMesh->NElements();
    
    TPZVec<REAL> coord(3,0.);
    int firstPt = 0;
    coord[0] = poligonalChain[firstPt].first;
    coord[2] = poligonalChain[firstPt].second;
    
    int axe0 = 0;//axe X
    int axe1 = 2;//axe Z
    int axeNormal = 1;//axe Y
    TPZVec<REAL> qsi(2,0.);
    TPZGeoEl * firstGel = planeMesh->FindElementCaju(coord, qsi, this->fInitialElIndex, 2);
    
    TPZGeoEl * gel = firstGel;
    TPZGeoEl * nextGel = NULL;
    
#ifdef PZDEBUG
    if(!gel)
    {
        std::cout << "first point of crack tip boundary does NOT belong to any 2D element" << std::endl;
        std::cout << "(or was given a mesh with zero quantity of 2D elements)!" << std::endl;
        std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
#endif
    
    REAL alphaMin;
    bool reachNextPoint;
    int64_t nElsCrossed, thispoint, nextpoint;
    std::map< int, std::set<REAL> >::iterator it;
    std::set<REAL> trim;
    TPZVec< TPZVec<REAL> > intersectionPoint;
    TPZManVector<REAL,3> x(3,0.), dx(3,0.);
    TPZManVector<REAL,3> xNext(3,0.), qsi2D(2,0.);
    TPZVec <int64_t> Topol(2), edgeVec;
    
    int64_t p;
    for(p = 0; p < (npoints-1); p++)
    {
        nElsCrossed = 0;
        alphaMin = 0.;
        thispoint = p;
        nextpoint = p+1;
        
        x[0] = poligonalChain[thispoint].first;
        x[2] = poligonalChain[thispoint].second;
        
        dx[0] = poligonalChain[nextpoint].first - poligonalChain[thispoint].first;
        dx[2] = poligonalChain[nextpoint].second - poligonalChain[thispoint].second;
        
        xNext[0] = poligonalChain[nextpoint].first;
        xNext[2] = poligonalChain[nextpoint].second;
        
        this->fLfrac = Max(this->fLfrac,xNext[0]);
        
        REAL norm = 0.;
        for(int c = 0; c < 3; c++)
        {
            norm += dx[c]*dx[c];
        }
        norm = sqrt(norm);
        for(int c = 0; c < 3; c++)
        {
            dx[c] = dx[c]/norm;
        }
        
        REAL Tol;
        ZeroTolerance(Tol);
        reachNextPoint = gel->ComputeXInverse(xNext, qsi2D,Tol);
        int axe0 = 0;//axe X
        int axe1 = 2;//axe Z
        int axeNormal = 1;//axe Y
        while(reachNextPoint == false && nElsCrossed < nelem)
        {
            nextGel = this->CrossToNextNeighbour(gel, x, dx, alphaMin, auxElIndex_TrimCoords, auxElIndexSequence, true, axe0, axe1, axeNormal, false);
            
            alphaMin = __smallNum;//qdo vai de vizinho a vizinho ateh chegar no proximo ponto, nao deve-se incluir os (alphaX=0)
            if(nextGel->ComputeXInverse(xNext, qsi2D,Tol))
            {
                reachNextPoint = true;
            }
            gel = nextGel;
            nElsCrossed++;
        }
        if(nElsCrossed == nelem)
        {
            //Deve ter alternado entre vizinhos!!!
            DebugStop();
        }
    }
    
    TPZGeoEl * lastGel = gel;
    
    dx[0] = -1.;//direcao oposta ao eixo x do sistema de coordenadas
    dx[1] =  0.;//direcao oposta ao eixo x do sistema de coordenadas
    dx[2] =  0.;//direcao oposta ao eixo x do sistema de coordenadas
    
    p = 0;//Fechando o inicio da fratura
    thispoint = p;
    nextpoint = p+1;
    
    x[0] = poligonalChain[thispoint].first;
    x[1] = 0.;
    x[2] = poligonalChain[thispoint].second;
    
    nextGel = firstGel;
    gel = NULL;
    alphaMin = 0.;
    while(gel != nextGel)
    {
        gel = nextGel;
        nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, auxElIndex_TrimCoords, auxElIndexSequence, false, axe0, axe1, axeNormal, true);
        if(!nextGel)
        {
            DebugStop();
        }
        alphaMin = __smallNum;
    }
    
    //*********** TESTE ENGLAND-GREEN
    //    {
    //        std::cout << "\n\n\n\n\nTESTE DO PHIL DE ENGLAND-GREEN (POLIGONAL CHAIN EH HORIZONTAL)\n\n\n\n\n";
    //        dx[0] = +1.;//direcao oposta ao eixo x do sistema de coordenadas
    //        dx[1] =  0.;//direcao oposta ao eixo x do sistema de coordenadas
    //        dx[2] =  0.;//direcao oposta ao eixo x do sistema de coordenadas
    //    }
    
    p = npoints;//Fechando o final da fratura
    thispoint = p-1;
    
    x[0] = poligonalChain[thispoint].first;
    x[1] = 0.;
    x[2] = poligonalChain[thispoint].second;
    
    nextGel = lastGel;
    gel = NULL;
    alphaMin = 0.;
    while(gel != nextGel)
    {
        gel = nextGel;
        nextGel = this->CrossToNextNeighbour(gel, x, dx, alphaMin, auxElIndex_TrimCoords, auxElIndexSequence, true, axe0, axe1, axeNormal, true);
        if(!nextGel)
        {
            DebugStop();
        }
        alphaMin = __smallNum;
    }
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl * TPZPlaneFractureMesh::CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, REAL alphaMin,
                                                      std::map< int64_t, std::set<REAL> > &auxElIndex_TrimCoords,
                                                      std::list< std::pair<int64_t,REAL> > &auxElIndexSequence,
                                                      bool pushback, int planeAxe0, int planeAxe1, int planeNormal,
                                                      bool closingFracture)
{
    bool thereIsAn1DElemAlready;
    int edge;
    std::map< int64_t, std::set<REAL> >::iterator it;
    TPZVec< TPZVec<REAL> > ExactIntersectionPoint, ModulatedIntersectionPoint;
    TPZVec<int> edgeVec;
    REAL Tol;
    ZeroTolerance(Tol);
    
    TPZGeoMesh * planeMesh = gel->Mesh();
    
    bool haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint,
                                             alphaMin, planeAxe0, planeAxe1, planeNormal);
    
    if(haveIntersection == false)
    {
        haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint,
                                            0., planeAxe0, planeAxe1, planeNormal);
        
#ifdef PZDEBUG
        if(haveIntersection == false)
        {
            TPZManVector<REAL,3> qsi2D(2,0.);
            if(gel->ComputeXInverse(x, qsi2D,Tol))
            {
                std::cout << "The point inside element does NOT intersect its edges! EdgeIntersection method face an exeption!" << std::endl;
                std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
                
                std::cout << "Element " << gel->Index() << std::endl;
                for(int64_t n = 0; n < gel->NNodes(); n++)
                {
                    std::cout << "Node " << n << std::endl;
                    gel->NodePtr(n)->Print(std::cout);
                    std::cout << std::endl;
                }
                std::cout << "x: " << x[0] << " , " << x[1] << " , " << x[2] << std::endl;
                std::cout << "dx: " << dx[0] << " , " << dx[1] << " , " << dx[2] << std::endl;
                std::cout << "alphaMin = " << alphaMin << std::endl << std::endl;
                
                //Pontos da corrente poligonal coincidem?
                DebugStop();
            }
            else
            {
                std::cout << "Trying to find edge intersection from an point that doesnt belong to given element!" << std::endl;
                std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
                
                DebugStop();
            }
        }
#endif
    }
    TPZVec<REAL> xLin(3);
    for(int64_t edg = 0; edg < edgeVec.NElements(); edg++)
    {
        edge = edgeVec[edg];
        TPZVec<REAL> n0(3), n1(3);
        planeMesh->NodeVec()[gel->SideNodeIndex(edge, 0)].GetCoordinates(n0);
        planeMesh->NodeVec()[gel->SideNodeIndex(edge, 1)].GetCoordinates(n1);
        for(int c = 0; c < 3; c++)
        {
            REAL coordM = ModulatedIntersectionPoint[edg][c];
            xLin[c] = coordM;
        }
        
        REAL qsi1D = LinearComputeXInverse(xLin, n0, n1);
        
        TPZGeoElSide gelEdge(gel, edge);
        TPZGeoElSide neighEdge = gelEdge.Neighbour();
        thereIsAn1DElemAlready = false;
        while(neighEdge != gelEdge)
        {
            if(neighEdge.Element()->Dimension() == 1)//jah existe um elemento 1D inserido nesta aresta!
            {
                thereIsAn1DElemAlready = true;
                TPZTransform<> transBetweenNeigh = neighEdge.NeighbourSideTransform(gelEdge);
                
                qsi1D *= transBetweenNeigh.Mult()(0,0);
                qsi1D = ((int)(qsi1D*100.))/100.;//arredondamento para evitar problemas de comparacao entre doubles na insercao do set!
                
                int64_t neighEdgeIndex = neighEdge.Element()->Index();
                it = auxElIndex_TrimCoords.find(neighEdgeIndex);
                if(it == auxElIndex_TrimCoords.end())//Devia ser entao elemento 1D de materialId == globMaterialIdGen.__1DbulletMat
                {
                    std::set<REAL> trim;
                    trim.insert(qsi1D);
                    auxElIndex_TrimCoords[neighEdgeIndex] = trim;
                    
#ifdef just1IntersectionByEdge
                    if(pushback)
                    {
                        auxElIndexSequence.push_back(std::make_pair(neighEdgeIndex, qsi1D));
                    }
                    else // push_FRONT
                    {
                        auxElIndexSequence.push_front(std::make_pair(neighEdgeIndex, qsi1D));
                    }
#endif
                }
#ifndef just1IntersectionByEdge
                else
                {
                    it->second.insert(qsi1D);
                }
                if(pushback)
                {
                    auxElIndexSequence.push_back(std::make_pair(neighEdgeIndex, qsi1D));
                }
                else // push_FRONT
                {
                    auxElIndexSequence.push_front(std::make_pair(neighEdgeIndex, qsi1D));
                }
#endif
                break;
            }
            neighEdge = neighEdge.Neighbour();
        }
        if(thereIsAn1DElemAlready == false)//nao existe um elemento 1D nesta aresta!
        {
            qsi1D = ((int)(qsi1D*100.))/100.;//arredondamento para evitar problemas de comparacao entre doubles na insercao do set!
            
            TPZVec<int64_t> Topol(2);
            
            std::set<REAL> trim;
            trim.insert(qsi1D);
            Topol[0] = gel->SideNodeIndex(edge, 0);
            Topol[1] = gel->SideNodeIndex(edge, 1);
            TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * linGeo =
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, globMaterialIdGen.Aux1DMatId(), *planeMesh);
            
            planeMesh->BuildConnectivity();
            
            int64_t linGeoIndex = linGeo->Index();
            auxElIndex_TrimCoords[linGeoIndex] = trim;
            
            if(pushback) // push_BACK
            {
                auxElIndexSequence.push_back(std::make_pair(linGeoIndex, qsi1D));
            }
            else // push_FRONT
            {
                auxElIndexSequence.push_front(std::make_pair(linGeoIndex, qsi1D));
            }
        }
    }
    x = ExactIntersectionPoint[ExactIntersectionPoint.NElements() - 1];
    TPZGeoElSide gelEdge(gel, edge);
    TPZGeoElSide neighEdge = gelEdge.Neighbour();
    bool neighFound = false;
    while(neighEdge != gelEdge)
    {
        if(neighEdge.Element()->Dimension() == 2 &&
           globMaterialIdGen.IsOutsideFractMat(neighEdge.Element()->MaterialId()) &&
           neighEdge.Element()->HasSubElement() == false)//substitui por essa verificacao!!!
        {
            neighFound = true;
            break;
        }
        neighEdge = neighEdge.Neighbour();
    }
    if(neighFound)
    {
        return neighEdge.Element();
    }
    else if(closingFracture == false)
    {
        DebugStop();//fratura saiu para fora do dominio
    }
    else if(gelEdge.Element()->Dimension() == 2 &&
            globMaterialIdGen.IsOutsideFractMat(gelEdge.Element()->MaterialId()) &&
            gelEdge.Element()->HasSubElement() == false)//substitui por essa verificacao!!!
    {
        //fechando fratura no inicio ou no final, portanto o
        //ultimo elemento nao encontra vizinho do tipo 2d com
        //materialId == globMaterialIdGen.__2DfractureMat_outside, retornando a si mesmo
        return gelEdge.Element();
    }
    
    return NULL;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureMesh::EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
                                            TPZVec< TPZVec<REAL> > &ExactIntersect, REAL alphaMin,
                                            int planeAxe0, int planeAxe1, int planeNormal)
{
    int64_t nearNode;
    bool IsNearNode = TPZChangeEl::NearestNode(gel, x, nearNode, __smallNum);
    
    edge.Resize(0);
    ExactIntersect.Resize(0);
    
    int ncnodes = gel->NCornerNodes();
    
    TPZVec< TPZVec<REAL> > node(ncnodes), dnode(ncnodes);
    TPZVec<REAL> nodeCoord(3);
    for(int n = 0; n < ncnodes; n++)
    {
        node[n].Resize(3,1);
        dnode[n].Resize(3,1);
        gel->NodePtr(n)->GetCoordinates(nodeCoord);
        for(int c = 0; c < 3; c++)
        {
            (node[n])[c] = nodeCoord[c];
        }
    }
    std::list<REAL> edgeNorm;
    /** <REAL: alpha , <REAL: alphaNodemod, REAL: alphaNodesmooth, REAL: norm, int: edge that intersect> > */
    std::map<REAL, TPZVec<REAL> > alpha;
    
    REAL alphaX, alphaNodemod, alphaNodesmooth, norm;
    for(int n = 0; n < ncnodes; n++)
    {
        norm = 0.;
        for(int c = 0; c < 3; c++)//computing the direction from node N to node N+1
        {
            if(n != (ncnodes-1) )
            {
                (dnode[n])[c] = node[n+1][c] - node[n][c];
            }
            else
            {
                dnode[n][c] = node[0][c] - node[n][c];
            }
            norm += dnode[n][c]*dnode[n][c];
        }
        norm = sqrt(norm);
        for(int c = 0; c < 3; c++)//normalizing computed direction
        {
            dnode[n][c] /= norm;
        }
        
        alphaX = ComputeAlphaX(x, dx, node[n], dnode[n], planeAxe0, planeAxe1, planeNormal);
        alphaNodemod = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, true, false, planeAxe0, planeAxe1, planeNormal);
        alphaNodesmooth = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, IsNearNode, true, planeAxe0, planeAxe1, planeNormal);
        if(alphaX >= alphaMin && alphaNodesmooth >= 0. && alphaNodesmooth <= norm)
        {
            int thisEdge = n+ncnodes;
            TPZVec<REAL> someData(4);
            someData[0] = alphaNodemod; someData[1] = alphaNodesmooth; someData[2] = norm; someData[3] = thisEdge;
            alpha[alphaX] = someData;
        }
    }
    
    // a aresta que serah interseccionada serah a que tiver menor alphaX não negativo, i.e.: o primeiro par do mapa alpha!
    std::map<REAL, TPZVec<REAL> >::iterator it = alpha.begin();
    if(it != alpha.end())
    {
        alphaX = it->first;
        alphaNodemod = (it->second)[0];
        alphaNodesmooth = (it->second)[1];
        norm = (it->second)[2];
        edge.Resize(1); edge[0] = int((it->second)[3]);
        
        ExactIntersect.Resize(1); ExactIntersect[0].Resize(3,1);
        for(int c = 0; c < 3; c++)
        {
            ExactIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodesmooth*dnode[edge[0]-ncnodes][c];
        }
        
        if(alphaX <= alphaMin && alpha.size() > 1)
        {
            it++;
            alphaX = it->first;
            alphaNodemod = (it->second)[0];
            alphaNodesmooth = (it->second)[1];
            norm = (it->second)[2];
            edge.Resize(2); edge[1] = int((it->second)[3]);
            
            ExactIntersect.Resize(2); ExactIntersect[1].Resize(3,1);
            for(int c = 0; c < 3; c++)
            {
                ExactIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodesmooth*dnode[edge[1]-ncnodes][c];
            }
        }
        return true;
    }
    
    //se o ponto p estah sobre um noh e a direcao dp aponta para fora do elemento...
    if(IsNearNode)
    {
        std::cout << "Estah no noh " << nearNode << " e nao foi encontrada interseccao!" << std::endl;
        std::cout << "Tratar este caso!" << std::endl << std::endl;
    }
    
    return false;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureMesh::EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
                                            TPZVec< TPZVec<REAL> > &ExactIntersect,
                                            TPZVec< TPZVec<REAL> > &ModulatedIntersect, REAL alphaMin,
                                            int planeAxe0, int planeAxe1, int planeNormal)
{
    int64_t nearNode;
    bool IsNearNode = TPZChangeEl::NearestNode(gel, x, nearNode, __smallNum);
    
    edge.Resize(0);
    ExactIntersect.Resize(0);
    ModulatedIntersect.Resize(0);
    
    int ncnodes = gel->NCornerNodes();
    
    TPZVec< TPZVec<REAL> > node(ncnodes), dnode(ncnodes);
    TPZVec<REAL> nodeCoord(3);
    for(int n = 0; n < ncnodes; n++)
    {
        node[n].Resize(3,1);
        dnode[n].Resize(3,1);
        gel->NodePtr(n)->GetCoordinates(nodeCoord);
        for(int c = 0; c < 3; c++)
        {
            (node[n])[c] = nodeCoord[c];
        }
    }
    std::list<REAL> edgeNorm;
    /** <REAL: alpha , <REAL: alphaNodemod, REAL: alphaNodesmooth, REAL: norm, int: edge that intersect> > */
    std::map<REAL, TPZVec<REAL> > alpha;
    
    REAL alphaX, alphaNodemod, alphaNodesmooth, norm;
    for(int n = 0; n < ncnodes; n++)
    {
        norm = 0.;
        for(int c = 0; c < 3; c++)//computing the direction from node N to node N+1
        {
            if(n != (ncnodes-1) )
            {
                (dnode[n])[c] = node[n+1][c] - node[n][c];
            }
            else
            {
                dnode[n][c] = node[0][c] - node[n][c];
            }
            norm += dnode[n][c]*dnode[n][c];
        }
        norm = sqrt(norm);
        for(int c = 0; c < 3; c++)//normalizing computed direction
        {
            dnode[n][c] /= norm;
        }
        
        alphaX = ComputeAlphaX(x, dx, node[n], dnode[n], planeAxe0, planeAxe1, planeNormal);
        alphaNodemod = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, true, false, planeAxe0, planeAxe1, planeNormal);
        alphaNodesmooth = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, IsNearNode, true, planeAxe0, planeAxe1, planeNormal);
        if(alphaX >= alphaMin && alphaNodesmooth >= 0. && alphaNodesmooth <= norm)
        {
            int thisEdge = n+ncnodes;
            TPZVec<REAL> someData(4);
            someData[0] = alphaNodemod; someData[1] = alphaNodesmooth; someData[2] = norm; someData[3] = thisEdge;
            alpha[alphaX] = someData;
        }
    }
    
    // a aresta que serah interseccionada serah a que tiver menor alphaX não negativo, i.e.: o primeiro par do mapa alpha!
    std::map<REAL, TPZVec<REAL> >::iterator it = alpha.begin();
    if(it != alpha.end())
    {
        alphaX = it->first;
        alphaNodemod = (it->second)[0];
        alphaNodesmooth = (it->second)[1];
        norm = (it->second)[2];
        edge.Resize(1); edge[0] = int((it->second)[3]);
        
        ModulatedIntersect.Resize(1); ModulatedIntersect[0].Resize(3,1);
        for(int c = 0; c < 3; c++)
        {
            ModulatedIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodemod*dnode[edge[0]-ncnodes][c];
        }
        ExactIntersect.Resize(1); ExactIntersect[0].Resize(3,1);
        for(int c = 0; c < 3; c++)
        {
            ExactIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodesmooth*dnode[edge[0]-ncnodes][c];
        }
        
        if(alphaX <= alphaMin && alpha.size() > 1)
        {
            it++;
            alphaX = it->first;
            alphaNodemod = (it->second)[0];
            alphaNodesmooth = (it->second)[1];
            norm = (it->second)[2];
            edge.Resize(2); edge[1] = int((it->second)[3]);
            
            ModulatedIntersect.Resize(2); ModulatedIntersect[1].Resize(3,1);
            for(int c = 0; c < 3; c++)
            {
                ModulatedIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodemod*dnode[edge[1]-ncnodes][c];
            }
            ExactIntersect.Resize(2); ExactIntersect[1].Resize(3,1);
            for(int c = 0; c < 3; c++)
            {
                ExactIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodesmooth*dnode[edge[1]-ncnodes][c];
            }
        }
        return true;
    }
    
    //se o ponto p estah sobre um noh e a direcao dp aponta para fora do elemento...
    if(IsNearNode)
    {
        std::cout << "Estah no noh " << nearNode << " e nao foi encontrada interseccao!" << std::endl;
        std::cout << "Tratar este caso!" << std::endl << std::endl;
    }
    
    return false;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureMesh::ComputeAlphaNode(TPZVec<REAL> &x, TPZVec<REAL> &dx,
                                            TPZVec<REAL> &node, TPZVec<REAL> &dnode,
                                            REAL norm, bool modulate, bool smooth,
                                            int planeAxe0, int planeAxe1, int planeNormal)
{
    REAL fractionNumQ =	dx[planeAxe1]*node[planeAxe0] - dx[planeAxe0]*node[planeAxe1] - dx[planeAxe1]*x[planeAxe0] + dx[planeAxe0]*x[planeAxe1];
    
    REAL fractionDenomQ =	dnode[planeAxe1]*dx[planeAxe0] - dnode[planeAxe0]*dx[planeAxe1];
    
    REAL alphaNode = -1.;
    if(fabs(fractionDenomQ) > __smallNum)
    {
        alphaNode = fractionNumQ/fractionDenomQ;
        if(fabs(alphaNode) < __smallNum) alphaNode = 0.;
        if(fabs(alphaNode - norm) < __smallNum) alphaNode = norm;
    }
    
    if(modulate)
    {
        int stretchesQTD = __EdgeStretchesQTD;
        if(smooth)
        {
            stretchesQTD *= __TrimQTDmultiplier;
        }
        int nsegm = int(alphaNode/(norm/stretchesQTD) + 0.5);
        if(nsegm == 0)//Tirando a interseccao do noh inicial da aresta
        {
            nsegm = 1;
        }
        else if(nsegm == stretchesQTD)//Tirando a interseccao do noh final da aresta
        {
            nsegm = stretchesQTD - 1;
        }
        alphaNode = nsegm*(norm/stretchesQTD);//modulando o ponto para multiplo de (norm/__EdgeStretchesQTD)
    }
    
    return alphaNode;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureMesh::ComputeAlphaX(TPZVec<REAL> &x, TPZVec<REAL> &dx,
                                         TPZVec<REAL> &node, TPZVec<REAL> &dnode,
                                         int planeAxe0, int planeAxe1, int planeNormal)
{
    //computing alpha (dx vector multiplier to intersect element edges)
    REAL fractionNumP   =	dnode[planeAxe1]*node[planeAxe0] - dnode[planeAxe0]*node[planeAxe1] -
    dnode[planeAxe1]*x[planeAxe0] + dnode[planeAxe0]*x[planeAxe1];
    
    REAL fractionDenomP =	dnode[planeAxe1]*dx[planeAxe0] - dnode[planeAxe0]*dx[planeAxe1];
    
    REAL alphaX = -1.;
    if(fabs(fractionDenomP) > __smallNum)
    {
        alphaX = fractionNumP/fractionDenomP;
        if(fabs(alphaX) < __smallNum)
        {
            alphaX = 0.;
        }
    }
    
    return alphaX;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureMesh::LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1)
{
    REAL dL = 0., L = 0.;
    for(int c = 0; c < 3; c++)
    {
        L += (n1[c]-n0[c])*(n1[c]-n0[c]);
        dL += (x[c]-n0[c])*(x[c]-n0[c]);
    }
    L = sqrt(L);
    dL = sqrt(dL);
    
#ifdef PZDEBUG
    if(fabs(L) < __smallNum || fabs(dL) < __smallNum || fabs(L - dL) < __smallNum)
    {
        std::cout << "n0 and n1 are coincident nodes!" << std::endl;
        std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
#endif
    
    REAL qsi = -1. + dL/L*2.;
    
    return qsi;
}
//------------------------------------------------------------------------------------------------------------

TPZAutoPointer<TPZRefPattern> TPZPlaneFractureMesh::Generate1DRefPatt(std::set<REAL> &TrimCoord)
{
    if(TrimCoord.size() == 0)
    {
        return NULL;
    }
    
    int64_t Qnodes = TrimCoord.size() + 2;
    int64_t Qelements = TrimCoord.size() + 2;
    int64_t QsubElements = Qelements - 1;
    TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
    for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
    
    //1. setting initial and final nodes coordinates of 1D element mesh
    NodeCoord[0][0] = -1.;
    NodeCoord[1][0] =  1.;
    //.
    
    //2. setting intermediate nodes coordinates of 1D element mesh
    int c = 2;
    std::set<REAL>::iterator it, itnext;
    for(it = TrimCoord.begin(); it != TrimCoord.end(); it++)
    {
        itnext = it;
        itnext++;
        
        REAL coord = *it;
        if(itnext != TrimCoord.end())
        {
            REAL coordNext = *itnext;
            if(fabs(coordNext - coord) < 1.E-3)
            {
                DebugStop();
            }
        }
        
        (NodeCoord[c])[0] = coord;
        c++;
    }
    //.
    
    //3. initializing internal mesh of refPattern
    TPZGeoMesh internalMesh;
    internalMesh.NodeVec().Resize(Qnodes);
    TPZVec <TPZGeoNode> Node(Qnodes);
    for(int n = 0; n < Qnodes; n++)
    {
        Node[n].SetNodeId(n);
        Node[n].SetCoord(NodeCoord[n]);
        internalMesh.NodeVec()[n] = Node[n];
    }
    //.
    
    //4. inserting 1D elements on internal mesh of refPattern
    TPZVec<int64_t> Topol(2);
    
    //4.1 inserting father element
    Topol[0] = 0; Topol[1] = 1;
    TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol,globMaterialIdGen.Aux1DMatId(),internalMesh);
    //.
    
    //4.2 inserting subelements
    //first subelement
    Topol[0] = 0; Topol[1] = 2;
    TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son1 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol,globMaterialIdGen.Aux1DMatId(),internalMesh);
    son1->SetFather(father);
    son1->SetFather(father->Index());
    //
    
    //last subelement
    Topol[0] = TrimCoord.size() + 1; Topol[1] = 1;
    TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son2 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol,globMaterialIdGen.Aux1DMatId(),internalMesh);
    son2->SetFather(father);
    son2->SetFather(father->Index());
    //
    
    for(int64_t el = 2; el < QsubElements; el++)
    {
        Topol[0] = el; Topol[1] = el+1;
        TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (Topol,globMaterialIdGen.Aux1DMatId(),internalMesh);
        son->SetFather(father);
        son->SetFather(father->Index());
    }
    //.
    //.
    
    internalMesh.BuildConnectivity();
    internalMesh.SetMaxNodeId(internalMesh.NNodes()-1);
    internalMesh.SetMaxElementId(internalMesh.NElements()-1);
    
    TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(internalMesh);
    TPZAutoPointer<TPZRefPattern> Found = gRefDBase.FindRefPattern(refPattern);
    if(!Found)
    {
        gRefDBase.InsertRefPattern(refPattern);
        refPattern->InsertPermuted();
        
        return refPattern;
    }
    else
    {
        return Found;
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::UpdatePoligonalChain(TPZGeoMesh * gmesh,
                                                std::list< std::pair<int64_t,REAL> > &auxElIndexSequence,
                                                TPZVec<std::pair<REAL,REAL> > &poligonalChainUpdated)
{
    int nptos = auxElIndexSequence.size();
    poligonalChainUpdated.Resize(nptos);
    
    TPZVec<REAL> qsi1Dvec(1), ptoCoord(3);
    int el1Did, p = 0;
    REAL qsi1D;
    
    std::list< std::pair<int64_t,REAL> >::iterator it;
    for(it = auxElIndexSequence.begin(); it != auxElIndexSequence.end(); it++)
    {
        el1Did = it->first;
        TPZGeoEl * el1D = gmesh->ElementVec()[el1Did];
        
#ifdef PZDEBUG
        int elDim = el1D->Dimension();
        if(elDim != 1)
        {
            std::cout << "The auxElIndexSequence supposedly would contains ids of elements exclusively 1D!" << std::endl;
            std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
#endif
        
        qsi1D = it->second;
        qsi1Dvec[0] = qsi1D;
        
        el1D->X(qsi1Dvec, ptoCoord);
        
        poligonalChainUpdated[p] = std::make_pair(ptoCoord[0],ptoCoord[2]);
        
        p++;
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureMesh::GenerateCrackBoundary(TPZGeoMesh * refinedMesh,
                                                 std::list< std::pair<int64_t,REAL> > &auxElIndexSequence)
{
    this->fcrackBoundaryElementsIndexes.Resize(0);
    TPZVec<REAL> qsi0vec(1), qsi1vec(1), node0coord(3), node1coord(3);
    TPZVec<int64_t> Topol(2);
    
    REAL qsi0, qsi1;
    std::list< std::pair<int64_t,REAL> >::iterator crackit0, crackit1, crackitEnd;
    crackitEnd = auxElIndexSequence.end(); crackitEnd--;
    for(crackit0 = auxElIndexSequence.begin(); crackit0 != crackitEnd; crackit0++)
    {
        crackit1 = crackit0; crackit1++;
        
        int64_t el0index = crackit0->first;
        TPZGeoEl * el0 = refinedMesh->ElementVec()[el0index];
        qsi0 = crackit0->second;
        qsi0vec[0] = qsi0;
        el0->X(qsi0vec, node0coord);
        int64_t n0 = TPZChangeEl::NearestNode(refinedMesh, node0coord, __smallNum);
        Topol[0] = n0;
        
        int64_t el1index = crackit1->first;
        TPZGeoEl * el1 = refinedMesh->ElementVec()[el1index];
        qsi1 = crackit1->second;
        qsi1vec[0] = qsi1;
        el1->X(qsi1vec, node1coord);
        int64_t n1 = TPZChangeEl::NearestNode(refinedMesh, node1coord, __smallNum);
        Topol[1] = n1;
        
        TPZGeoEl * crack1D = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, globMaterialIdGen.CrackTipMatId(), *refinedMesh);
        int oldSize = this->fcrackBoundaryElementsIndexes.NElements();
        this->fcrackBoundaryElementsIndexes.Resize(oldSize+1);
        this->fcrackBoundaryElementsIndexes[oldSize] = crack1D->Index();
    }
    
    refinedMesh->BuildConnectivity();
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFractureMesh::SeparateElementsInMaterialSets(TPZGeoMesh * refinedMesh)
{
    this->fInsideFractureMatId.clear();
    
    int n1Dels = this->fcrackBoundaryElementsIndexes.NElements();
    std::map<int,locFracture2DEl> fracturedElems;
    
    
    int maxNStripes = 9;//limitacao do globMaterialIdGen
    int mult = 1; // (stripe width) = (mult * this->fLmax)
    int predictedNStripes = 0;
    
    if(this->fjust1stripe == false)
    {
        predictedNStripes = this->fLfrac / (mult * this->fLmax) + 1;
        while(predictedNStripes > maxNStripes)
        {
            mult++;
            predictedNStripes = this->fLfrac / (mult * this->fLmax) + 1;
        }
    }
    
    this->fnstripes = 0;
    
    //Capturando subelementos que encostam no contorno da fratura
    for(int el = 0; el < n1Dels; el++)
    {
        int cracktipIndex = this->fcrackBoundaryElementsIndexes[el];
        TPZGeoEl * gel = refinedMesh->ElementVec()[cracktipIndex];//1D element of crach boundary
        
#ifdef PZDEBUG
        if(gel->Dimension() != 1)
        {
            DebugStop();
        }
#endif
        int sd = 2;
        TPZGeoElSide side1D(gel,sd);
        TPZGeoElSide sideNeigh = side1D.Neighbour();
        while(sideNeigh != side1D)
        {
            int sideNeighbyside = sideNeigh.Side();
            TPZGeoEl * neigh = sideNeigh.Element();
            
            if(neigh->HasSubElement() || neigh->Dimension() != 2)
            {
                sideNeigh = sideNeigh.Neighbour();
                continue;
            }
            
            //Primeiro estagio da identificacao dos elementos 2D que estao no interior da fratura.
            //brief: primeiramente sao identificados os elementos no interior da fratura que encostam no crack tip.
            //Obs.: A condicao (IsBoundaryMaterial(neigh) == false) eh para excluir os quadrilateros de condicao de contorno dirichlet.
            if(globMaterialIdGen.IsBoundaryMaterial(neigh->MaterialId()) == false)
            {
                TPZVec<REAL> neighCenterQSI(neigh->Dimension()), neighCenterX(3);
                neigh->CenterPoint(neigh->NSides()-1, neighCenterQSI);
                neigh->X(neighCenterQSI, neighCenterX);
                
                TPZVec<REAL> n0(3), n1(3);
                refinedMesh->NodeVec()[gel->SideNodeIndex(sd, 0)].GetCoordinates(n0);
                refinedMesh->NodeVec()[gel->SideNodeIndex(sd, 1)].GetCoordinates(n1);
                
                //Como o contorno da fratura foi construido no sentido antihorario no plano x,z (normal Y > 0),
                //interessam os elementos aa direita do elemento 1D. Portanto eh feito produto vetorial entre os vetores
                //frac=(n1-n0) e cg_neigh=(cg-n0). O vizinho aa direita apresentarah componente em Y positiva.
                REAL crossYcomp = n0[2]*n1[0] - n0[0]*n1[2] -
                n0[2]*neighCenterX[0] + n1[2]*neighCenterX[0] +
                n0[0]*neighCenterX[2] - n1[0]*neighCenterX[2];
                
                if(crossYcomp > 0.)
                {
                    if(neigh->LowestFather() != neigh)
                    {//O elemento pai define que stripe seus filhos estao pois uma stripe tem largura do elemento pai.
                        neigh->LowestFather()->CenterPoint(neigh->LowestFather()->NSides()-1, neighCenterQSI);
                        neigh->LowestFather()->X(neighCenterQSI, neighCenterX);
                    }
                    REAL Xc = neighCenterX[0];
                    REAL Zc = neighCenterX[2];
                    
                    int stripe = 0;
                    if(this->fjust1stripe == false)
                    {
                        stripe = (int)(Xc/(mult*this->fLmax));
                        this->fnstripes = MAX(this->fnstripes,stripe+1);
                    }
                    int layer = globLayerStruct.WhatLayer(Zc);
                    
                    int fractMatId = globMaterialIdGen.InsideFractMatId(layer, stripe);
                    neigh->SetMaterialId(fractMatId);
                    this->fInsideFractureMatId.insert(fractMatId);
                    
                    std::map<int,locFracture2DEl>::iterator edgIt_temp = fracturedElems.find(neigh->Index());
                    if(edgIt_temp != fracturedElems.end())
                    {
                        locFracture2DEl neighFractEl = edgIt_temp->second;
                        neighFractEl.RemoveThisEdge(sideNeighbyside);
                        edgIt_temp->second = neighFractEl;
                    }
                    else
                    {
                        locFracture2DEl fractEl(neigh);
                        fractEl.RemoveThisEdge(sideNeighbyside);
                        fracturedElems[fractEl.Index()] = fractEl;
                    }
                }
            }
            
            sideNeigh = sideNeigh.Neighbour();
        }
    }
    
    //capturanto demais elementos no interior da fratura
    std::set<int> finishedFracturedElems;
    while(fracturedElems.size() > 0)
    {
        std::map<int,locFracture2DEl>::iterator edgIt = fracturedElems.begin();
        
        locFracture2DEl actEl = edgIt->second;
        
        std::set<int>::iterator sideIt;
        for(sideIt = actEl.fEdge.begin(); sideIt != actEl.fEdge.end(); sideIt++)
        {
            int side = *sideIt;
            TPZGeoElSide actElEdge(actEl.fElem2D,side);
            TPZGeoElSide neighElSide = actElEdge.Neighbour();
            
            bool wellDone = false;
            while(actElEdge != neighElSide && wellDone == false)
            {
                int sideNeighbyside = neighElSide.Side();
                TPZGeoEl * neighEl = neighElSide.Element();
                if(neighEl->Dimension() == 2 && globMaterialIdGen.IsOutsideFractMat(neighEl->MaterialId()) && !neighEl->HasSubElement())
                {
                    int neighElIndex = neighEl->Index();
                    if(finishedFracturedElems.find(neighElIndex) != finishedFracturedElems.end())
                    {
                        wellDone = true;
                        continue;
                    }
                    std::map<int,locFracture2DEl>::iterator edgIt_temp = fracturedElems.find(neighElIndex);
                    if(edgIt_temp != fracturedElems.end())
                    {
                        locFracture2DEl neighFractEl = edgIt_temp->second;
                        neighFractEl.RemoveThisEdge(sideNeighbyside);
                        
                        if(neighFractEl.IsOver())
                        {
                            finishedFracturedElems.insert(neighFractEl.Index());
                            fracturedElems.erase(edgIt_temp);
                        }
                        else
                        {
                            edgIt_temp->second = neighFractEl;
                        }
                    }
                    else
                    {
                        TPZVec<REAL> neighCenterQSI(neighEl->Dimension()), neighCenterX(3);
                        neighEl->LowestFather()->CenterPoint(neighEl->LowestFather()->NSides()-1, neighCenterQSI);
                        neighEl->LowestFather()->X(neighCenterQSI, neighCenterX);
                        
                        REAL Xc = neighCenterX[0];
                        REAL Zc = neighCenterX[2];
                        
                        int stripe = 0;
                        if(this->fjust1stripe == false)
                        {
                            stripe = (int)(Xc/(mult*this->fLmax));
                            this->fnstripes = MAX(this->fnstripes,stripe+1);
                        }
                        int layer = globLayerStruct.WhatLayer(Zc);
                        
                        int fractMatId = globMaterialIdGen.InsideFractMatId(layer, stripe);
                        neighEl->SetMaterialId(fractMatId);
                        this->fInsideFractureMatId.insert(fractMatId);
                        
                        locFracture2DEl fractEl(neighEl);
                        fractEl.RemoveThisEdge(sideNeighbyside);
                        fracturedElems[fractEl.Index()] = fractEl;
                    }
                    
                    wellDone = true;
                }
                neighElSide = neighElSide.Neighbour();
            }
        }
        
        finishedFracturedElems.insert(actEl.Index());
        fracturedElems.erase(edgIt);
    }
}
//------------------------------------------------------------------------------------------------------------
