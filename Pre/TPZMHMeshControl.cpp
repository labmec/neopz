//
//  TPZMHMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/3/14.
//
//

#include "TPZMHMeshControl.h"
#include "TPZRefPattern.h"
#include "Projection/TPZL2Projection.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include "TPZCompElLagrange.h"
#include "TPZBndCond.h"
#include "TPZInterfaceEl.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzsubcmesh.h"
#include "TPZRefPatternTools.h"
#include "pzlog.h"
#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#ifdef PZ_LOG
static TPZLogger logger("pz.mhmeshcontrol");
static TPZLogger loggerWRAP("pz.mhmeshwrap");
#endif

// toto


TPZMHMeshControl::TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &geotomhm) : fGMesh(gmesh), fProblemType(EScalar), fNState(1), fGeoToMHMDomain(geotomhm), fMHMtoSubCMesh(), fGlobalSystemSize(-1), fGlobalSystemWithLocalCondensationSize(-1), fNumeq(-1)
{
    if(!gmesh) DebugStop();
#ifdef PZDEBUG
    if (geotomhm.size() != fGMesh->NElements()) {
        DebugStop();
    }
#endif
    fpOrderInternal = 2;
    fpOrderSkeleton = 1;

    int64_t nc = geotomhm.size();
    for (int64_t c=0; c<nc; c++) {
        fMHMtoSubCMesh[geotomhm[c] ] = -1;
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coarse element indexes ";
        for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fCMesh = new TPZMultiphysicsCompMesh(fGMesh);
    fCMesh->SetDimModel(fGMesh->Dimension());
    fPressureFineMesh = new TPZCompMesh(fGMesh);
    fFluxMesh = new TPZCompMesh(fGMesh);
}

/// Define the partitioning information of the MHM mesh
void TPZMHMeshControl::DefinePartition(TPZVec<int64_t> &partitionindex, std::map<int64_t,std::pair<int64_t,int64_t> > &skeleton)
{
#ifdef PZDEBUG
    if (partitionindex.size() != partitionindex.size()) {
        DebugStop();
    }
#endif
    fGeoToMHMDomain = partitionindex;
    fInterfaces = skeleton;
    fMHMtoSubCMesh.clear();
    int64_t nc = partitionindex.size();
    for (int64_t c=0; c<nc; c++) {
        fMHMtoSubCMesh[partitionindex[c] ] = -1;
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coarse element indexes ";
        for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// Define the partitioning information of the MHM mesh
// This method calculates the skeleton indexes
void TPZMHMeshControl::DefinePartition(TPZVec<int64_t> &partitionindex) {

    TPZGeoMesh *gmesh = this->fGMesh.operator->();

    // Creates skeleton elements
    std::map<int64_t, std::pair<int64_t, int64_t>> skeleton;

    int dim = gmesh->Dimension();

    const int64_t nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {

        TPZGeoEl *gel = gmesh->Element(iel);

        if (gel->HasSubElement()) continue;
        if (gel->Dimension() == dim-1)
        {
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if(neighbour.Element()->Dimension() != dim) DebugStop();
            int partition = partitionindex[neighbour.Element()->Index()];
            int64_t neighindex = neighbour.Element()->Index();
            auto val = std::pair<int64_t,int64_t>(partition,gel->Index());
            skeleton[gel->Index()] = val;
            continue;
        }

        // Iterates through the sides of the element
        int nsides = gel->NSides();

        for (int iside = 0; iside < nsides; iside++) {

            TPZGeoElSide gelside(gel, iside);

            // Filters boundary sides
            if (gelside.Dimension() != dim - 1) continue;

            TPZGeoElSide neighbour = gelside.Neighbour();

            // Filters neighbour sides that belong to volume elements
            if (neighbour.Element()->Dimension() != dim) continue;

            int64_t thisId = gel->Index();
            int64_t neighId = neighbour.Element()->Index();

            // This lambda checks if a skeleton has already been created over gelside
            auto hasSkeletonNeighbour = [&]() -> bool {
                neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    int neighbourMatId = neighbour.Element()->MaterialId();
                    if (neighbourMatId == fSkeletonMatId) {
                        return true;
                    }
                   neighbour = neighbour.Neighbour();
                }
                return false;
            };

            int64_t thisGroup = partitionindex[thisId];
            int64_t neighGroup = partitionindex[neighId];

            if (thisGroup != neighGroup) {
                if (!hasSkeletonNeighbour()) {
                    TPZGeoElBC bc(gelside, fSkeletonMatId);
                    std::pair<int64_t, int64_t> adjacentGroupIndexes(thisGroup, neighGroup);
                    int64_t bc_index = bc.CreatedElement()->Index();
                    if(bc_index == neighGroup)
                    {
                        // you are unwillingly creating a boundary skeleton
                        DebugStop();
                    }
                    skeleton.insert({bc.CreatedElement()->Index(), adjacentGroupIndexes});
                }
            }
        }
    }

    DefinePartition(partitionindex, skeleton);
}

TPZMHMeshControl::TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : fGMesh(gmesh), fGeoToMHMDomain(),
                                                                       fMHMtoSubCMesh(), fGlobalSystemSize(-1),
                                                                       fGlobalSystemWithLocalCondensationSize(-1),
                                                                       fNumeq(-1) {


#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coarse element indexes ";
        for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fCMesh = new TPZMultiphysicsCompMesh(fGMesh);
    fCMesh->SetDimModel(fGMesh->Dimension());
    fPressureFineMesh = new TPZCompMesh(fGMesh);
    fFluxMesh = new TPZCompMesh(fGMesh);

}

TPZMHMeshControl::TPZMHMeshControl(const TPZMHMeshControl &copy){

    this->operator=(copy);
}

TPZMHMeshControl &TPZMHMeshControl::operator=(const TPZMHMeshControl &cp){

    fGMesh = cp.fGMesh;
    fProblemType = cp.fProblemType;
    fNState = cp.fNState;
    fSkeletonMatId = cp.fSkeletonMatId;
    fLagrangeMatIdLeft = cp.fLagrangeMatIdLeft;
    fLagrangeMatIdRight = cp.fLagrangeMatIdRight;
    fGeoToMHMDomain = cp.fGeoToMHMDomain;
    fMHMtoSubCMesh = cp.fMHMtoSubCMesh;
    fLagrangeAveragePressure = cp.fLagrangeAveragePressure;
    fCMesh = cp.fCMesh;
    fCMesh->SetReference(fGMesh);
    fPressureFineMesh = cp.fPressureFineMesh;
    fFluxMesh = cp.fFluxMesh;
    fpOrderInternal = cp.fpOrderInternal;
    fSkeletonWrapMatId = cp.fSkeletonWrapMatId;
    fBoundaryWrapMatId = cp.fBoundaryWrapMatId;
    fInternalWrapMatId = cp.fInternalWrapMatId;
    fHybridize = cp.fHybridize;
    fSwitchLagrangeSign = cp.fSwitchLagrangeSign;
    fGlobalSystemSize = cp.fGlobalSystemSize;
    fGlobalSystemWithLocalCondensationSize = cp.fGlobalSystemWithLocalCondensationSize;
    fNumeq = cp.fNumeq;
    return *this;
}

TPZMHMeshControl::~TPZMHMeshControl()
{
    fGMesh->ResetReference();
    int64_t nel = fCMesh->NElements();
    for (int64_t el = 0; el<nel ; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        if(!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            continue;
        }
        delete cel;
    }
    for (int64_t el = 0; el<nel ; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        cel->LoadElementReference();
        delete cel;
    }
}


/// Define the MHM partition by the coarse element indices
void TPZMHMeshControl::DefinePartitionbyCoarseIndices(TPZVec<int64_t> &coarseindices)
{
    int64_t ncoarse = coarseindices.size();
    int64_t nel = fGMesh->NElements();
    int dim = fGMesh->Dimension();
    fGeoToMHMDomain.Resize(nel);
    for (int64_t el=0; el<nel; el++) {
        fGeoToMHMDomain[el] = -1;
    }
    for (int64_t el=0; el<ncoarse; el++) {
        fGeoToMHMDomain[coarseindices[el]] = coarseindices[el];
        TPZGeoEl *gel = fGMesh->Element(coarseindices[el]);
        // initialize the index of the subdomain associated with the geometric element
        if (gel) {
            fMHMtoSubCMesh[gel->Index()] = -1;
        }
    }
    // set the MHM domain index for all elements
    // look for a father element that has a MHM domain index
    // adopt the domain index of the ancester
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || fGeoToMHMDomain[el] != -1) {
            continue;
        }
        int domain = -1;
        TPZGeoEl *father = gel->Father();
        while (father) {
            int64_t index = father->Index();
            if (fGeoToMHMDomain[index] != -1) {
                fGeoToMHMDomain[el] = fGeoToMHMDomain[index];
                break;
            }
            father = father->Father();
        }
    }
    // if all neighbours of a lower dimensional element belong to the same subdomain, then
    // the element is put in this domain
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || fGeoToMHMDomain[el] != -1) {
            continue;
        }
        if(gel->Dimension() == dim) continue;
        TPZGeoElSide gelside(gel);
        int domain = -1;
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            auto index = neighbour.Element()->Index();
            int neighdomain = fGeoToMHMDomain[index];
            if(domain == -1 && neighdomain != -1) domain = neighdomain;
            // we found neighbours of different domain
            if(neighdomain != -1 && neighdomain != domain) break;
            neighbour = neighbour.Neighbour();
        }
        if(domain != -1 && gelside == neighbour)
        {
            fGeoToMHMDomain[el] = domain;
        }
    }
    CreateSkeletonElements();
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Subdomain indices " << fGeoToMHMDomain << std::endl;
        sout << "Recognized domains ";
        for (auto it = fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << it->first << ' ';
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}


/// will create 1D elements on the interfaces between the coarse element indices
void TPZMHMeshControl::CreateSkeletonElements()
{
    if (fInterfaces.size()) {
        DebugStop();
    }

    if(fSkeletonMatId < 0) DebugStop();

    TPZCompMesh *cmesh = CriaMalhaTemporaria();

    int64_t nel = fGMesh->NElements();
    int dimension = fGMesh->Dimension();
    int ninterf;

    for(int64_t iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = fGMesh->ElementVec()[iel];
        if(!gel) continue;

        ninterf = gel->NumInterfaces();
        if(ninterf > 1) DebugStop();
        if (ninterf==1)
        {
            TPZCompEl *cel = gel->Reference();
            TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
            if (!intface) {
                DebugStop();
            }
            int interfacematid = fSkeletonMatId;
            TPZCompEl *left = intface->LeftElement();
            TPZCompEl *right = intface->RightElement();
            int64_t leftind = left->Reference()->Index();
            int64_t rightind = right->Reference()->Index();
            if (left->Reference()->Dimension() == dimension-1) {
                // the boundary element is the skeleton element
                continue;
            }
            if (right->Reference()->Dimension() == dimension-1) {
                // the boundary element is the skeleton element
                continue;
            }
            fInterfaces[iel] = std::make_pair(leftind, rightind);
            gel->SetMaterialId(interfacematid);
            // in order to prevent the element from being deleted
            gel->DecrementNumInterfaces();
        }
    }

    fGMesh->ResetReference();
    delete cmesh;

    BuildWrapMesh();
//    BuildWrapMesh(fGMesh->Dimension()-1);

    fGeoToMHMDomain.Resize(fGMesh->NElements(), -1);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it = fInterfaces.begin();
        while (it != fInterfaces.end()) {
            int leftdim = fGMesh->ElementVec()[it->second.first]->Dimension();
            int rightdim = fGMesh->ElementVec()[it->second.second]->Dimension();
            sout << "Interface index " << it->first << " Left Element " << it->second.first << "/" << leftdim
            << " Right Element " << it->second.second << "/" << rightdim << std::endl;
            it++;
        }
        sout << "Geometric mesh after creating the skeleton\n";
        fGMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// divide the skeleton elements
void TPZMHMeshControl::DivideSkeletonElements(int ndivide)
{
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it;
    for (int divide=0; divide<ndivide; divide++)
    {
        std::map<int64_t, std::pair<int64_t,int64_t> > mapdivided;
        for (it=fInterfaces.begin(); it!=fInterfaces.end(); it++) {
            int64_t elindex = it->first;
//            if (elindex == it->second.second) {
//                mapdivided[elindex] = it->second;
//                continue;
//            }
            TPZGeoEl *gel = fGMesh->Element(elindex);
            TPZAutoPointer<TPZRefPattern> refpat = TPZRefPatternTools::PerfectMatchRefPattern(gel);
            gel->SetRefPattern(refpat);
            TPZManVector<TPZGeoEl *,10> subels;
            gel->Divide(subels);
            int64_t nsub = subels.size();
            for (int is=0; is<nsub; is++) {
                if (subels[is]->Index() >= fGeoToMHMDomain.size()) {
                    fGeoToMHMDomain.Resize(subels[is]->Index()+1000, -1);
                }
                fGeoToMHMDomain[subels[is]->Index()] = fGeoToMHMDomain[elindex];
                mapdivided[subels[is]->Index()] = it->second;
                // for boundary elements, the second element is the interface element
                if(elindex == it->second.second)
                {
                    mapdivided[subels[is]->Index()].second = subels[is]->Index();
                }
            }
        }
        fInterfaces = mapdivided;
    }
//    BuildWrapMesh(fGMesh->Dimension());
//    BuildWrapMesh(fGMesh->Dimension()-1);
    fGeoToMHMDomain.Resize(fGMesh->NElements(), -1);

    std::cout<<"WrapMatId created \n";
    std::cout << "fSkeletonWrapMatId "<<fSkeletonWrapMatId<<std::endl;
    std::cout << "fBoundaryWrapMatId "<<fBoundaryWrapMatId<<std::endl;
    //std::cout << "fHdivWrapMatId "<<fHDivWrapMatid<<std::endl;
}


TPZCompMesh* TPZMHMeshControl::CriaMalhaTemporaria()
{
    fGMesh->ResetReference();
	/// criar materiais
	int dim = fGMesh->Dimension();
    std::set<int> matids, bcids;
    TPZGeoMesh &gmesh = fGMesh;
    int64_t nel = gmesh.NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fMHMtoSubCMesh.find(el) == fMHMtoSubCMesh.end()) {
            continue;
        }
        int materialid = gel->MaterialId();
        matids.insert(materialid);
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == dim-1) {
                    bcids.insert(neighbour.Element()->MaterialId());
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }

    TPZCompMesh *cmesh = new TPZCompMesh(fGMesh);
	cmesh->SetDimModel(dim);
    TPZManVector<STATE,1> sol(1,0.);
    int nstate = 1;
    std::set<int>::iterator it = matids.begin();
    TPZMaterialT<STATE> *meshmat = nullptr;
    while (it != matids.end()) {
        TPZL2Projection<STATE> *material = new TPZL2Projection(*it,dim,nstate,sol);
        cmesh->InsertMaterialObject(material);
        if (!meshmat) {
            meshmat = material;
        }
        it++;

    }

    if (!meshmat) {
        DebugStop();
    }

	it = bcids.begin();
    while (it != bcids.end()) {
        ///Inserir condicao de contorno
        TPZFMatrix<STATE> val1(2,2,0.);
        TPZVec<STATE> val2(2,0.);

        TPZBndCond * BCondD1 = meshmat->CreateBC(meshmat, *it,0, val1, val2);
        cmesh->InsertMaterialObject(BCondD1);

        it++;
    }

    cmesh->SetAllCreateFunctionsDiscontinuous();
    TPZCreateApproximationSpace &create = cmesh->ApproxSpace();

    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fMHMtoSubCMesh.find(el) == fMHMtoSubCMesh.end()) {
            continue;
        }
        create.CreateCompEl(gel, *cmesh);
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == dim-1) {
                    int matid = neighbour.Element()->MaterialId();
                    if (bcids.find(matid) != bcids.end()) {
                        create.CreateCompEl(neighbour.Element(), *cmesh);
                    }
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }

    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
    cmesh->ExpandSolution();
    return cmesh;
}

/// Create all data structures for the computational mesh
void TPZMHMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    int nstate = fNState;
    int dim = fGMesh->Dimension();
    fCMesh->SetDimModel(dim);
    InsertPeriferalMaterialObjects();
    // create the H1 elements
    CreateInternalElements();
    // create the HDiv elements
    CreateSkeleton();
    if(fLagrangeAveragePressure)
    {
        CreateLagrangeMultiplierMesh();
    }
    // creating the multiphysics space
    {
        auto meshvecauto = GetMeshes();
        TPZVec<TPZCompMesh *> meshvec(meshvecauto.size());
        for(int m=0; m<meshvec.size(); m++)
        {
            meshvec[m] = meshvecauto[m].operator->();
        }
        fCMesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
        fCMesh->BuildMultiphysicsSpace(meshvec);
        int64_t nel = fCMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCMesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            int mhmdomain = fGeoToMHMDomain[gel->Index()];
            SetSubdomain(cel, mhmdomain);
        }
    }
    CreateInterfaceElements();
    
    
//    AddBoundaryInterfaceElements();
    fCMesh->ExpandSolution();
    fCMesh->CleanUpUnconnectedNodes();


#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "*********** BEFORE SUBSTRUCTURING *************\n";
        fCMesh->Print(sout);
        int64_t nel = fCMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMesh->Element(el);
            if(!cel) continue;
            sout << "el index " << el << " subdomain " << WhichSubdomain(cel) << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();
    if(usersubstructure==true){
        this->SubStructure();
    }
    fNumeq = fCMesh->NEquations();
}

/// will create the internal elements, one coarse element at a time
void TPZMHMeshControl::CreateInternalElements()
{
    TPZCompEl::SetgOrder(fpOrderInternal);

    if(!fPressureFineMesh) fPressureFineMesh = new TPZCompMesh(fGMesh);
    fPressureFineMesh->ApproxSpace().SetCreateLagrange(false);
    fPressureFineMesh->SetAllCreateFunctionsContinuous();

    //Criar elementos computacionais malha MHM

    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    TPZCompMesh *cmesh = fPressureFineMesh.operator->();
    fConnectToSubDomainIdentifier[cmesh].Expand(10000);

    int64_t nel = fGMesh->NElements();
    for (auto itMHM = fMHMtoSubCMesh.begin(); itMHM != fMHMtoSubCMesh.end(); itMHM++)
    {
        fGMesh->ResetReference();
        bool LagrangeCreated = false;
        for (int64_t el=0; el< nel; el++)
        {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (!gel || gel->HasSubElement() || fGeoToMHMDomain[el] != itMHM->first) {
                continue;
            }
            
            // create the flux element
            TPZCompEl *cel = fPressureFineMesh->CreateCompEl(gel);
            const int64_t index = cel->Index();
        
            /// associate the connects with the subdomain
            SetSubdomain(cel, itMHM->first);
            // we need to create a lagrange multiplier element in order to delay decomposition of an equation
            if (!LagrangeCreated)
            {
                LagrangeCreated = true;
                TPZCompEl *cel = fPressureFineMesh->ElementVec()[index];
                int64_t cindex = cel->ConnectIndex(0);
                int nshape(1), nvar(1), order(1);
                int lagrangelevel = 0;
                // in case we have a space representing the average pressure, then
                // one connect will be ordered after the distributed flux which has lagrange level 0
                // in case we dont have an average pressure, one connect needs to be ordered after
                // the skeleton fluxes. The skeleton fluxes have lagrange level 2
                if (this->fLagrangeAveragePressure) {
                    lagrangelevel = 1;
                }
                else
                {
                    lagrangelevel = 3;
                }
                fPressureFineMesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                if (fProblemType == EElasticity2D) {
                    cindex = cel->ConnectIndex(2);
                    fPressureFineMesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                }
                if (fProblemType == EElasticity3D) {
                    cindex = cel->ConnectIndex(6);
                    fPressureFineMesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                }
            }
        }
    }
    fGMesh->ResetReference();
    fPressureFineMesh->ExpandSolution();
}

/// will create the elements on the skeleton
void TPZMHMeshControl::CreateSkeleton()
{
    if(!fFluxMesh) fFluxMesh = new TPZCompMesh(fGMesh);
    // comment this line or not to switch the type of skeleton elements
    int meshdim = fFluxMesh->Dimension();
    fFluxMesh->SetDimModel(meshdim);
//    fFluxMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    fFluxMesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    int order = fpOrderSkeleton;
    if (order < 0) {
        order = 0;
    }
    fFluxMesh->SetDefaultOrder(order);
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        int64_t elindex = it->first;
        // skip the boundary elements
        if (elindex == it->second.second) {
            it++;
            DebugStop();
            continue;
        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        // create a discontinuous element to model the flux
        TPZCompEl *cel = fFluxMesh->CreateCompEl(gel);
        int Side = gel->NSides()-1;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            int lagrangelevel = 0;
            if (this->fLagrangeAveragePressure) {
                lagrangelevel = 2;
            }
            else
            {
                lagrangelevel = 2;
            }
            cel->Connect(ic).SetLagrangeMultiplier(lagrangelevel);
        }
        SetSubdomain(cel, -1);

        if (it->second.first < it->second.second) {
            // set the flux orientation depending on the relative value of the element ids
            intel->SetSideOrient(Side, 1);
        }
        else
        {
            intel->SetSideOrient(Side, -1);
        }
        SetSubdomain(cel, -1);
        gel->ResetReference();
        it++;
    }
    fFluxMesh->SetDimModel(meshdim);
    fFluxMesh->ExpandSolution();
    fGMesh->ResetReference();
}

/// will create the interface elements between the internal elements and the skeleton
void TPZMHMeshControl::CreateInterfaceElements()
{
    if(!fCMesh)
    {
        std::cout << "You should create the multiphysics mesh first\n" << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    fCMesh->LoadReferences();
    int dim = fGMesh->Dimension();
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it;
    /// loop over the skeleton elements
    for (it=fInterfaces.begin(); it != fInterfaces.end(); it++) {
        // index of the skeleton element
        int64_t elindex = it->first;
        // left and right indexes in the coarse mesh
        int64_t leftelindex = it->second.first;
        int64_t rightelindex = it->second.second;
        // skip boundary elements
        if(rightelindex == elindex) continue;
        int matid = 0, matidleft = 0, matidright = 0;

        // second condition indicates a boundary element
        if (leftelindex < rightelindex || rightelindex == elindex)
        {
            matidleft = fLagrangeMatIdLeft;
            matidright = fLagrangeMatIdRight;
        }
        else
        {
            matidleft = fLagrangeMatIdRight;
            matidright = fLagrangeMatIdLeft;
        }
        int numlado = 2;
        if (rightelindex == elindex) {
            numlado = 1;
        }

        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZCompElSide celskeleton = gelside.Reference();
        // the skeleton element must exist
        if(!celskeleton.Element()) DebugStop();
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            int neighmatid = neighbour.Element()->MaterialId();
            if(neighmatid == fSkeletonWrapMatId || neighmatid == fBoundaryWrapMatId)
            {
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        if (neighbour == gelside) {
            DebugStop();
        }
        TPZStack<TPZGeoElSide> gelstack;
        gelstack.Push(neighbour);
        while (gelstack.size())
        {
            TPZGeoElSide smallGeoElSide = gelstack.Pop();
            // the smaller elements returned by GetSubElements include element/sides of lower dimension
            if (smallGeoElSide.Dimension() != gel->Dimension()) {
                continue;
            }
            // look for the neighbours of smallGeoElSide for elements which are of dimension of the fGMesh
            // and which beint64_t to the subdomains
            TPZGeoElSide neighbour = smallGeoElSide.Neighbour();
            while (neighbour != smallGeoElSide)
            {
                if (neighbour.Element()->Dimension() != fGMesh->Dimension()) {
                    neighbour = neighbour.Neighbour();
                    continue;
                }
                TPZCompElSide csmall = neighbour.Reference();
                if (csmall)
                {
                    int matid = -1;
                    if(WhichSubdomain(csmall.Element()) == leftelindex) {
                        matid = matidleft;
                    }
                    else if(WhichSubdomain(csmall.Element()) == rightelindex)
                    {
                        matid = matidright;
                    }
                    if (matid == -1) {
                        DebugStop();
                    }
                    // create an interface between the finer element and the MHM flux
                    TPZGeoEl *gelnew = smallGeoElSide.Element()->CreateBCGeoEl(smallGeoElSide.Side(), matid);
                    new TPZMultiphysicsInterfaceElement(fCMesh, gelnew, csmall, celskeleton);
#ifdef PZ_LOG
                    if (logger.isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "New interface left " << smallGeoElSide.Element()->Index() << " right " << gel->Index() << " matid " << matid;
                        sout << " interface index " << gelnew->Reference()->Index() << " beint64_ts to subdomain " << WhichSubdomain(gelnew->Reference());
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
                neighbour = neighbour.Neighbour();
            }

            smallGeoElSide.GetSubElements2(gelstack);
        }
    }
}


/// verify if the element is a sibling of
bool TPZMHMeshControl::IsSibling(int64_t son, int64_t father)
{
    TPZGeoEl *gel = fGMesh->ElementVec()[son];
    while (gel && gel->Index() != father) {
        gel = gel->Father();
    }
    if (gel) {
        return true;
    }
    else
    {
        return false;
    }
}

/// Print diagnostics
void TPZMHMeshControl::PrintDiagnostics(std::ostream &out)
{
    int nelgroups = this->fMHMtoSubCMesh.size();
    out << __PRETTY_FUNCTION__ << " Number of coarse elements " << nelgroups << std::endl;

    for (std::map<int64_t,int64_t>::iterator it = fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
        PrintSubdomain(it->first, out);
    }
    PrintBoundaryInfo(out);
}

/// print the diagnostics for a subdomain
void TPZMHMeshControl::PrintSubdomain(int64_t elindex, std::ostream &out)
{
    int64_t nel = fCMesh->NElements();
    std::set<int64_t> celindices;
    std::multimap<int64_t,int64_t> interfaces;
    TPZStack<TPZCompElSide> boundaries;
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        int64_t gelindex = gel->Index();
        if(fGeoToMHMDomain[gelindex] == elindex)
        {
            celindices.insert(el);
            // verify whether their neighbours are in the subdomain
            AddElementBoundaries(elindex, el, boundaries);
        }
        TPZMultiphysicsInterfaceElement *iface = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (iface) {
            TPZCompEl *left = iface->LeftElement();
            TPZGeoEl *gleft = left->Reference();
            if (IsSibling(gleft->Index(), elindex)) {
                interfaces.insert(std::make_pair(left->Index(),el));
            }
        }
    }
    out << "Diagnostics for subdomain formed by element seed " << elindex << std::endl;
    out << "Number of computational elements contained in the group " << celindices.size() << std::endl;
    out << "Comp Element indices ";
    for (std::set<int64_t>::iterator i=celindices.begin(); i != celindices.end(); i++) {
        out << *i << " ";
    }
    out << std::endl;
    out << "Number of boundary interfaces " << interfaces.size() << std::endl;
    for (std::multimap<int64_t, int64_t>::iterator i = interfaces.begin(); i != interfaces.end(); i++) {
        int64_t el = i->second;
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        TPZMultiphysicsInterfaceElement *face = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (!face) {
            DebugStop();
        }
        out << "Comp Element index " << face->LeftElement()->Index() << " side " << face->Left().Side() << " skeleton element " << face->RightElement()->Index() << " interface index "
            << face->Index() << " face material id " << face->Reference()->MaterialId() << std::endl;
    }
    out << "Number of elements on the sides of elseed " << boundaries.size() << std::endl;
    for (int64_t i=0; i<boundaries.size(); i++) {
        out << "Comp Element index " << boundaries[i].Element()->Index() << " side " << boundaries[i].Side() << std::endl;
    }
    out << "Geom Element indices ";
    for (std::set<int64_t>::iterator i=celindices.begin(); i != celindices.end(); i++) {
        out << fCMesh->ElementVec()[*i]->Reference()->Index() << " ";
    }
    out << std::endl;
    out << "Number of boundary interfaces " << interfaces.size() << std::endl;
    for (std::multimap<int64_t, int64_t>::iterator i = interfaces.begin(); i != interfaces.end(); i++) {
        int64_t el = i->second;
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        TPZMultiphysicsInterfaceElement *face = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (!face) {
            DebugStop();
        }
        out << "Geom Element index " << face->LeftElement()->Reference()->Index() << " side " << face->Left().Side() << " skeleton element " << face->RightElement()->Reference()->Index() << " interface index "
        << face->Index() << " face material id " << face->Reference()->MaterialId() << std::endl;
    }
    out << "Number of elements on the sides of elseed " << boundaries.size() << std::endl;
    for (int64_t i=0; i<boundaries.size(); i++) {
        out << "Geom Element index " << boundaries[i].Element()->Reference()->Index() << " side " << boundaries[i].Side() << std::endl;
    }
}

/// put the element side which face the boundary on the stack
void TPZMHMeshControl::AddElementBoundaries(int64_t elseed, int64_t compelindex, TPZStack<TPZCompElSide> &result)
{
    TPZCompEl *cel = fCMesh->ElementVec()[compelindex];
    TPZGeoEl *gel = cel->Reference();
    int ns = gel->NSides();
    int dim = gel->Dimension();
    for (int is=0; is<ns; is++) {
        if (gel->SideDimension(is) != dim-1) {
            continue;
        }
        TPZGeoElSide gelside(gel,is);
        TPZGeoElSide father(gelside);
        while (father && father.Element()->Index() != elseed) {
            father = father.Father2();
        }
        // if the elseed is not a father aint64_t the boundary, the element is not on the boundary
        if (!father || father.Dimension() != dim-1) {
            continue;
        }
        TPZCompElSide celside(cel,is);

        result.Push(celside);
    }
}


/// print the indices of the boundary elements and interfaces
void TPZMHMeshControl::PrintBoundaryInfo(std::ostream &out)
{
    out << "Output for the structure of the boundary elements\n";
    std::set<int> MHMmatids;
    MHMmatids.insert(fLagrangeMatIdLeft);
    MHMmatids.insert(fLagrangeMatIdRight);
    MHMmatids.insert(fSkeletonMatId);
    int dim = fGMesh->Dimension();

    out << "\nCOMPUTATIONAL ELEMENT INFORMATION\n";

    int64_t nelem = fCMesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *celf = dynamic_cast<TPZInterfaceElement *>(cel);
        if (celf) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() > dim-1) {
            continue;
        }
        int matid = gel->MaterialId();
        if (MHMmatids.find(matid) != MHMmatids.end()) {
            continue;
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        out << "Comp Boundary Element Index " << el << " matid " << matid << std::endl;
        for (int64_t i=0; i< celstack.NElements(); i++) {
            TPZCompElSide cels = celstack[i];
            TPZGeoElSide gels = cels.Reference();
            TPZGeoEl *g = gels.Element();
            int matid = g->MaterialId();
            if (MHMmatids.find(matid) != MHMmatids.end()) {
                continue;
            }
            TPZCompEl *c = cels.Element();
            TPZInterfaceElement *f = dynamic_cast<TPZInterfaceElement *>(c);
            if (f) {
                TPZCompEl *left = f->LeftElement();
                TPZCompEl *right = f->RightElement();
                out << "Comp Interface leftel index " << left->Index() << " rightel index " << right->Index() << std::endl;
            }
        }

    }

    out << "\nGEOMETRIC ELEMENT INFORMATION\n";
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *celf = dynamic_cast<TPZInterfaceElement *>(cel);
        if (celf) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() > dim-1) {
            continue;
        }
        int matid = gel->MaterialId();
        if (MHMmatids.find(matid) != MHMmatids.end()) {
            continue;
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        out << "Geom Boundary Element Index " << gel->Index() << " matid " << matid << std::endl;
        for (int64_t i=0; i< celstack.NElements(); i++) {
            TPZCompElSide cels = celstack[i];
            TPZGeoElSide gels = cels.Reference();
            TPZGeoEl *g = gels.Element();
            int matid = g->MaterialId();
            if (MHMmatids.find(matid) != MHMmatids.end()) {
                continue;
            }
            TPZCompEl *c = cels.Element();
            TPZInterfaceElement *f = dynamic_cast<TPZInterfaceElement *>(c);
            if (f) {
                TPZCompEl *left = f->LeftElement();
                TPZCompEl *right = f->RightElement();
                out << "Geom Interface leftel index " << left->Reference()->Index() << " rightel index " << right->Reference()->Index() << std::endl;
            }
        }

    }
}

#include "pzvec_extras.h"

/// create the lagrange multiplier mesh, one element for each subdomain
void TPZMHMeshControl::CreateLagrangeMultiplierMesh()
{
    if(fCMeshDomainFlux) fCMeshDomainFlux = new TPZCompMesh(fGMesh);
    int dim = fGMesh->Dimension();
    fCMeshDomainFlux->SetDimModel(dim);
    fCMeshDomainFlux->SetAllCreateFunctionsDiscontinuous();
    fCMeshDomainFlux->SetDefaultOrder(0);
    fGMesh->ResetReference();
    int nstate = 1;
    switch (fProblemType) {
        case EScalar:
            nstate = 1;
            break;
        case EElasticity3D:
            nstate = 6;
            break;
        case EElasticity2D:
            nstate = 3;
            break;
        default:
            DebugStop();
            break;
    }
    int64_t connectcounter = fCMesh->NConnects();
	/// criar materiais
    std::set<int> matids;
    TPZGeoMesh &gmesh = fGMesh;
    int64_t nel = gmesh.NElements();
    // this code needs to be modified to create lagrange computational elements which share a connect
    // between each other
    //DebugStop();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fMHMtoSubCMesh.find(el) == fMHMtoSubCMesh.end()) {
            continue;
        }
        int materialid = gel->MaterialId();
        matids.insert(materialid);
    }

    TPZManVector<STATE,1> sol(1,0.);
    std::set<int>::iterator it = matids.begin();
    TPZMaterial *meshmat = 0;
    while (it != matids.end()) {
        auto *material = new TPZNullMaterial(*it,dim,nstate);
        fCMeshDomainFlux->InsertMaterialObject(material);
        if (!meshmat) {
            meshmat = material;
        }
        it++;

    }
    if (!meshmat) {
        DebugStop();
    }
    
    // data structure to compute the center and area of each subdomain
    std::map<int64_t,TPZManVector<REAL,3> > DomainCenters;
    std::map<int64_t,REAL> DomainAreas;
    
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fMHMtoSubCMesh.find(el) == fMHMtoSubCMesh.end()) {
            continue;
        }
        TPZManVector<REAL,3> centerksi(dim),centerx(3);
        gel->CenterPoint(gel->NSides()-1, centerksi);
        gel->X(centerksi, centerx);
        REAL area = gel->SideArea(gel->NSides()-1);
        sscal(centerx,area);
        int64_t subdomain = fMHMtoSubCMesh[el];
        auto it = DomainCenters.find(subdomain);
        if(it ==  DomainCenters.end())
        {
            DomainCenters[subdomain] = centerx;
            DomainAreas[subdomain] = area;
        }
        else
        {
            saxpy(it->second,centerx,1.);
            DomainAreas[subdomain] += area;
        }

        TPZCompElDisc *disc = new TPZCompElDisc(fCMeshDomainFlux,gel);
        disc->SetTotalOrderShape();
        disc->SetFalseUseQsiEta();
        int64_t cindex = disc->ConnectIndex(0);
#ifdef PZDEBUG
        static int count = 0;
        if (count == 0)
        {
            TPZConnect &c = disc->Connect(0);
            std::cout << "Number of shape functions of discontinuous element " << c.NShape() << std::endl;
            count++;
        }
#endif
        SetSubdomain(disc, el);
    }
    {
        int64_t nel = fCMeshDomainFlux->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCMeshDomainFlux->Element(el);
            TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
            auto gelindex = cel->Reference()->Index();
            auto domain = fMHMtoSubCMesh[gelindex];
            TPZManVector<REAL,3> domaincenter = DomainCenters[domain];
            REAL area = DomainAreas[domain];
            for (int i=0; i<3; i++) {
                disc->SetCenterPoint(i, domaincenter[i]/area);
            }
        }
    }
    fCMeshDomainFlux->ExpandSolution();
    fGMesh->ResetReference();

    fCMeshDomainPressure = new TPZCompMesh(fCMeshDomainFlux);
    {
        int64_t nel = fCMeshDomainPressure->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMeshDomainPressure->Element(el);
            TPZCompEl *cel2 = fCMeshDomainFlux->Element(el);
            int subdomain = WhichSubdomain(cel2);
            SetSubdomain(cel, subdomain);
        }
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        fCMeshDomainFlux->Print(sout);
        fCMeshDomainPressure->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fGMesh->ResetReference();
}

/// transform the computational mesh into a multiphysics mesh
void TPZMHMeshControl::TransferToMultiphysics()
{
    fGMesh->ResetReference();
    this->fCMesh = new TPZMultiphysicsCompMesh(fGMesh);
    this->fCMesh->SetDimModel(fGMesh->Dimension());
    fCMesh->SetAllCreateFunctionsMultiphysicElem();

    // copy the material objects
    std::map<int,TPZMaterial *>::iterator it = fPressureFineMesh->MaterialVec().begin();
    while (it != fPressureFineMesh->MaterialVec().end()) {
        it->second->Clone(fCMesh->MaterialVec());
        it++;
    }

    int64_t nel = fPressureFineMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
        if (intface) {
            continue;
        }
        TPZCompElLagrange *lagr = dynamic_cast<TPZCompElLagrange *>(cel);
        if (lagr) {
            new TPZCompElLagrange(fCMesh,*cel);
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }

        TPZCompEl *celnew = fCMesh->CreateCompEl(gel);
        TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(celnew);
        if (!mult) {
            DebugStop();
        }
        mult->AddElement(cel, 0);
    }
    fCMesh->LoadReferences();
    nel = fCMeshDomainFlux->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMeshDomainFlux->ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.ConnectedCompElementList(celstack, 0, 0);
        if (gel->Reference()) {
            celstack.Push(gelside.Reference());
        }
        int nstack = celstack.size();
        for (int ist=0; ist<nstack; ist++) {
            TPZCompEl *celmult = celstack[ist].Element();
            TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(celmult);
            mult->AddElement(cel, 1);
        }
    }
    fGMesh->ResetReference();
    fCMesh->LoadReferences();
    nel = fCMeshDomainPressure->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMeshDomainPressure->ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.ConnectedCompElementList(celstack, 0, 0);
        if (gel->Reference()) {
            celstack.Push(gelside.Reference());
        }
        int nstack = celstack.size();
        for (int ist=0; ist<nstack; ist++) {
            TPZCompEl *celmult = celstack[ist].Element();
            TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(celmult);
            mult->AddElement(cel, 2);
        }
    }

    nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZMultiphysicsElement *multel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!multel) {
            continue;
        }
        int nmeshes = multel->NMeshes();
        for (int im=nmeshes; im<3; im++) {
            multel->AddElement(0, im);
        }
    }

    //void TPZBuildMultiphysicsMesh::AddConnects(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh)
    TPZManVector<TPZCompMesh *,3> cmeshvec(3,0);
    cmeshvec[0] = fPressureFineMesh.operator->();
    cmeshvec[1] = fCMeshDomainFlux.operator->();
    cmeshvec[2] = fCMeshDomainPressure.operator->();
    TPZCompMesh *cmesh = fCMesh.operator->();
    TPZBuildMultiphysicsMesh::AddConnects(cmeshvec,cmesh);

    nel = fPressureFineMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!intface) {
            continue;
        }
        //        TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right);
        TPZCompElSide pressleft = intface->LeftElementSide();
        TPZCompElSide pressright = intface->RightElementSide();
        TPZGeoElSide gleft = pressleft.Reference();
        TPZGeoElSide gright = pressright.Reference();
        TPZCompElSide multleft = gleft.Reference();
        TPZCompElSide multright = gright.Reference();
        new TPZMultiphysicsInterfaceElement(fCMesh,intface->Reference(),multleft,multright);

    }


    nel = fCMeshDomainPressure->NElements();
    int64_t npressconnect = fPressureFineMesh->NConnects();
    int64_t nlagrangeconnect = fCMeshDomainFlux->NConnects();
    // nel numero de dominios MHM, tem um connect associado a cada um e os mesmos estao no final
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = this->fCMeshDomainPressure->Element(el);
        int64_t pressureconnect = cel->ConnectIndex(0);
        int64_t cindex = npressconnect+nlagrangeconnect+pressureconnect;
        fCMesh->ConnectVec()[cindex].SetLagrangeMultiplier(3);
    }
    fCMesh->ExpandSolution();
    //fCMesh->SaddlePermute();
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        fCMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}


/// substructure the mesh
void TPZMHMeshControl::SubStructure()
{
    // for each connect index, the submesh index
    std::map<int64_t, int64_t > connectdest;
    // for each coarse geometric index, a subcompmesh
    std::map<int64_t, TPZSubCompMesh *> submeshes;
    std::map<int64_t,int64_t>::iterator it = fMHMtoSubCMesh.begin();

    // create the submeshes
    while (it != fMHMtoSubCMesh.end()) {
        TPZSubCompMesh *submesh = new TPZSubCompMesh(fCMesh);
        submeshes[it->first] = submesh;
        it++;
    }
    for (std::map<int64_t, TPZSubCompMesh *>::iterator it = submeshes.begin(); it != submeshes.end(); it++) {
        fMHMtoSubCMesh[it->first] = it->second->Index();
    }

    fGMesh->ResetReference();
    fCMesh->LoadReferences();

    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = fCMesh->Element(el);
        if(!cel) continue;
        if (dynamic_cast<TPZSubCompMesh *>(cel)) {
            continue;
        }
        int64_t domain = WhichSubdomain(cel);

        if (domain == -1) {
            continue;
        }
        if (submeshes.find(domain) == submeshes.end()) {
            DebugStop();
        }
        submeshes[domain]->TransferElement(fCMesh.operator->(), cel->Index());
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "Transferring element index " << cel->Index() << " geometric index ";
            TPZGeoEl *gel = cel->Reference();
            if (gel) {
                sout << gel->Index();
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
    }
    fCMesh->ComputeNodElCon();


    std::map<int64_t, TPZSubCompMesh *>::iterator itsub = submeshes.begin();
    while (itsub != submeshes.end()) {
        TPZSubCompMesh *submesh = itsub->second;
        int nc = submesh->NConnects();
        std::set<int64_t> internals;
        // put all connects with one element connection internal in the submesh
        for (int ic=0; ic<nc; ic++) {
            int64_t connectindex = submesh->ConnectIndex(ic);
            TPZConnect &c = submesh->Connect(ic);
            int lagrange = c.LagrangeMultiplier();
            if (c.NElConnected() >1) {
                continue;
            }
            bool makeinternal = false;
            // if hybridizing all internal connects can be condensed
            if (fHybridize) {
                makeinternal = true;
            }
            else if (lagrange < 3) {
                makeinternal = true;
            }
            int64_t internal = submesh->InternalIndex(connectindex);
            if (makeinternal)
            {
                internals.insert(internal);
            }
            else
            {
                c.IncrementElConnected();
#ifdef PZDEBUG
                std::cout << "For subdomain " << itsub->first << " connect index " << connectindex << " left external as lagrange multiplier\n";
#endif
            }
        }
        submesh->MakeAllInternal();
//        for (std::set<int64_t>::iterator it = internals.begin(); it != internals.end(); it++) {
//            submesh->MakeInternal(*it);
//        }
        submesh->InitializeBlock();
        itsub++;
    }
    fCMesh->CleanUpUnconnectedNodes();
    itsub = submeshes.begin();
    while (itsub != submeshes.end()) {
        TPZSubCompMesh *submesh = itsub->second;
        int numthreads = 0;
        int preconditioned = 0;
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "Newly created submesh for element " << *it << "\n";
            submesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZAutoPointer<TPZGuiInterface> guiInterface;
        submesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
        itsub++;
    }

    fCMesh->SaddlePermute();
}

/// print the data structure
void TPZMHMeshControl::Print(std::ostream &out)
{

    /// geometric mesh used to create the computational mesh
    if (fGMesh)
    {
        out << "******************* GEOMETRIC MESH *****************\n";
        fGMesh->Print(out);
    }

    /// computational mesh to contain the pressure elements
    // this mesh is the same as fCMesh if there are no lagrange multipliers assocated with the average pressure
    if (fPressureFineMesh)
    {
        out << "******************* PRESSURE MESH *****************\n";
        fPressureFineMesh->Print(out);
    }


    /// computational MHM mesh being built by this class
    if (fCMesh)
    {
        out << "******************* COMPUTATIONAL MESH *****************\n";
        fCMesh->Print(out);
    }

    /// computational mesh to represent the constant states
    if (fCMeshDomainFlux)
    {
        out << "******************* LAGRANGE MULTIPLIER MESH *****************\n";
        fCMeshDomainFlux->Print(out);
    }

    /// computational mesh to represent the constant states
    if (fCMeshDomainPressure)
    {
        out << "******************* CONSTANTE PRESSURE MESH *****************\n";
        fCMeshDomainPressure->Print(out);
    }

    /// material id associated with the skeleton elements
    out << "Skeleton Mat Id " <<  fSkeletonMatId << std::endl;

    /// material id associated with the lagrange multiplier elements
    out << "Lagrange mat id left - right " <<  fLagrangeMatIdLeft << " - " << fLagrangeMatIdRight << std::endl;

    /// interpolation order of the internal elements
    out << "Internal polynomial order " <<  fpOrderInternal << std::endl;

    /// interpolation order of the skeleton elements
    out << "Skeleton polynomial order " << fpOrderSkeleton << std::endl;

    /// indices of the geometric elements which define the skeleton mesh
    {
        out << "Geometric element indices of the coarse mesh ";
        for (std::map<int64_t,int64_t>::iterator it= fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            out << it->first << " " << it->second << " ";
        }
        out << std::endl;
    }
    /// indices of the skeleton elements and their left/right elements of the skeleton mesh
    out << "Skeleton element indices with associated left and right coarse element indices\n";
    {
        std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it;
        for (it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
            out << "skel index " << it->first << " Left Right indices " << it->second.first << " " << it->second.second << std::endl;
        }
    }

    /// flag to determine whether a lagrange multiplier is included to force zero average pressures in the subdomains
    /**
     * when imposing average pressure to be zero, a multiphysics mesh is created
     */
    out << "Will generate a constant pressure mesh " <<  fLagrangeAveragePressure << std::endl;
    /// subdonain indices of the connects
    out << "Subdomain indices of the connects\n";
    for (auto it = fConnectToSubDomainIdentifier.begin(); it != fConnectToSubDomainIdentifier.end(); it++)
    {
        out << "Mesh address " << (void *) it->first;
        TPZManVector<int64_t> &cvec = it->second;
        for (int64_t i=0; i<cvec.size(); i++) {
            if (i && !(i%20)) {
                out << std::endl;
            }
            out << "(" << i << "->" << cvec[i] << ") ";
        }
        out << std::endl;
    }

}

/// Insert Boundary condition objects that do not perform any actual computation
void TPZMHMeshControl::InsertPeriferalMaterialObjects()
{
    // The significant material objects have been inserted already
    int matid = *fMaterialIds.begin();
    TPZMaterial *mat = fCMesh->FindMaterial(matid);
    if (!mat) {
        DebugStop();
    }
    int nstate = fNState;
    int dim = fGMesh->Dimension();

    TPZFNMatrix<4,STATE> val1(fNState,fNState,0.), val2Flux(fNState,1,0.);
    TPZMaterial *matPerif = nullptr;

    if (fCMesh->FindMaterial(fSkeletonMatId)) {
        DebugStop();
    }
    matPerif = new TPZNullMaterialCS<STATE>(fSkeletonMatId,dim-1,nstate);
    fCMesh->InsertMaterialObject(matPerif);

    if (1) {

        int nstate = fNState;
        int dim = fGMesh->Dimension();

        auto *matwrap = new TPZNullMaterialCS<>(fSkeletonWrapMatId,dim-1,nstate);
        fCMesh->InsertMaterialObject(matwrap);
        
        if (fCMesh->FindMaterial(fLagrangeMatIdLeft)) {
            DebugStop();
        }
        if (fCMesh->FindMaterial(fLagrangeMatIdRight)) {
            DebugStop();
        }
        auto *matleft = new TPZLagrangeMultiplierCS<STATE>(fLagrangeMatIdLeft,dim,nstate);
        auto *matright = new TPZLagrangeMultiplierCS<STATE>(fLagrangeMatIdRight,dim,nstate);
        if (fSwitchLagrangeSign) {
            matleft->SetMultiplier(-1.);
            matright->SetMultiplier(1.);
        }
        else
        {
            matleft->SetMultiplier(1.);
            matright->SetMultiplier(-1.);
        }
        fCMesh->InsertMaterialObject(matleft);
        fCMesh->InsertMaterialObject(matright);
    }
    // insert material objects in the pressure mesh
    {
        auto *mat = new TPZNullMaterial<>(fSkeletonWrapMatId,this->fGMesh->Dimension()-1,fNState);
        fPressureFineMesh->InsertMaterialObject(mat);
        for(auto matid : fMaterialIds)
        {
            auto *mat = new TPZNullMaterial<>(matid,this->fGMesh->Dimension(),fNState);
            fPressureFineMesh->InsertMaterialObject(mat);
        }
        for(auto matid : fMaterialBCIds)
        {
            auto *mat = new TPZNullMaterial<>(matid,this->fGMesh->Dimension()-1,fNState);
            fPressureFineMesh->InsertMaterialObject(mat);
        }
    }
    // insert material objects in the flux mesh
    {
        {
            auto *mat = new TPZNullMaterial<>(this->fSkeletonMatId,this->fGMesh->Dimension(),fNState);
            fFluxMesh->InsertMaterialObject(mat);
        }
    }
    if(fLagrangeAveragePressure)
    {
        if(!fCMeshDomainFlux) fCMeshDomainFlux = new TPZCompMesh(fGMesh);
        if(!fCMeshDomainPressure) fCMeshDomainPressure = new TPZCompMesh(fGMesh);
        // insert material objects in the domain flux and average meshes
        for(auto matid : fMaterialIds)
        {
            int nstate = 1;
            if(fProblemType == EElasticity2D) nstate = 3;
            if(fProblemType == EElasticity3D) nstate = 6;
            TPZNullMaterial<> *mat = new TPZNullMaterial<>(matid,dim,nstate);
            fCMeshDomainFlux->InsertMaterialObject(mat);
            TPZNullMaterial<> *mat2 = new TPZNullMaterial<>(matid,dim,nstate);
            fCMeshDomainPressure->InsertMaterialObject(mat2);
        }
    }

}

/// associates the connects of an element with a subdomain
/// adds the association of geometric element with subdomain as well
void TPZMHMeshControl::SetSubdomain(TPZCompEl *cel, int64_t subdomain)
{
    int ncon = cel->NConnects();
    for (int ic=0; ic<ncon; ic++) {
        int64_t cindex = cel->ConnectIndex(ic);
        SetSubdomain(cel->Mesh(), cindex, subdomain);
    }
    TPZGeoEl *gel = cel->Reference();
    int64_t index = gel->Index();

    if (index >= fGeoToMHMDomain.size()) {
        fGeoToMHMDomain.Resize(index+1, -1);
        fGeoToMHMDomain[index] = subdomain;
    }
}


/// associates the connects index with a subdomain
void TPZMHMeshControl::SetSubdomain(TPZCompMesh *cmesh, int64_t cindex, int64_t subdomain)
{
    // cmesh indicates the atomic mesh. It can be a flux mesh or pressure mesh
    if (cindex >= fConnectToSubDomainIdentifier[cmesh].size()) {
        fConnectToSubDomainIdentifier[cmesh].Resize(cindex+1, -1);
    }
    fConnectToSubDomainIdentifier[cmesh][cindex] = subdomain;

}

/// returns to which subdomain a given element beint64_ts
// this method calls debugstop if the element beint64_ts to two subdomains
int64_t TPZMHMeshControl::WhichSubdomain(TPZCompEl *cel)
{
    int ncon = cel->NConnects();
    std::set<int64_t> domains;
    TPZCompMesh *cmesh = cel->Mesh();
    TPZManVector<int64_t> &cvec = fConnectToSubDomainIdentifier[cmesh];
    for (int ic=0; ic<ncon; ic++)
    {
        int64_t cindex = cel->ConnectIndex(ic);
        if (cvec[cindex] != -1) {
            domains.insert(cvec[cindex]);
        }
    }
    // if the element has connects in two different subdomains then something is wrong
    if (domains.size() > 1) {
        for (int ic=0; ic<ncon; ic++) {
            int64_t cindex = cel->ConnectIndex(ic);
            std::cout << cindex << "|" << cvec[cindex] << " ";
        }
        std::cout << std::endl;
        DebugStop();
    }
    if (domains.size() ==0) {
        return -1;
    }
    int64_t domain = *domains.begin();
    return domain;
}

/// Subdomains are identified by computational mesh, this method will join
// the subdomain indices of the connects to the multi physics mesh
void TPZMHMeshControl::JoinSubdomains(TPZVec<TPZCompMesh *> &meshvec, TPZCompMesh *multiphysicsmesh)
{
    int nmeshes = meshvec.size();
    int64_t totalconnects = 0;
    for (int im=0; im<nmeshes; im++) {
        if (fConnectToSubDomainIdentifier.find(meshvec[im]) == fConnectToSubDomainIdentifier.end()) {
            DebugStop();
        }
        totalconnects += fConnectToSubDomainIdentifier[meshvec[im]].size();
    }
    fConnectToSubDomainIdentifier[multiphysicsmesh].Resize(totalconnects);
    TPZManVector<int64_t> &mfvec = fConnectToSubDomainIdentifier[multiphysicsmesh];
    int64_t count = 0;
    for (int im=0; im<nmeshes; im++) {
        TPZManVector<int64_t> &cvec = fConnectToSubDomainIdentifier[meshvec[im]];
        int64_t nc = cvec.size();
        for (int64_t ic=0; ic<nc; ic++) {
            mfvec[count++] = cvec[ic];
        }
    }
}



bool IsAncestor(TPZGeoEl *son, TPZGeoEl *father)
{
    TPZGeoEl *check = son;
    while (check && check != father) {
        check = check->Father();
    }
    if (check==father) {
        return true;
    }

    return false;
}

void TPZMHMeshControl::BuildConnectedElementList(TPZGeoElSide gelside, std::list<TPZCompElSide> &ellist)
{
    int interfacedim = gelside.Dimension();
    TPZGeoEl *gel = gelside.Element();
    if(gel->Dimension() == interfacedim+1 && gel->Reference())
    {
        TPZCompElSide celside = gelside.Reference();
        ellist.push_back(celside);

    }
    // look for neighbours of gelside
    TPZGeoElSide neighbour = gelside.Neighbour();
    TPZGeoElSide neighdivided;
    while(neighbour != gelside)
    {
        TPZGeoEl *neighgel = neighbour.Element();
        if(neighgel->Dimension() == interfacedim+1 && neighgel->Reference())
        {
            TPZCompElSide neighcel = neighbour.Reference();
            ellist.push_back(neighcel);
        }
        if(neighbour.HasSubElement() && neighbour.NSubElements() > 1)
        {
            neighdivided = neighbour;
        }
        neighbour = neighbour.Neighbour();
    }
    if(neighdivided)
    {
        TPZStack<TPZGeoElSide> subelements;
        neighdivided.GetSubElements2(subelements);
        int64_t nsub = subelements.size();
        for(int64_t isub = 0; isub < nsub; isub++)
        {
            TPZGeoElSide gelsub = subelements[isub];
            if(gelsub.Dimension() == interfacedim)
            {
                BuildConnectedElementList(gelsub, ellist);
            }
        }
    }
}

/// identify connected elements to the skeleton elements
// the computational mesh is determined by the element pointed to by the geometric element
// skeleton index of the geometric element which defines the skeleton
// leftright indices of the left and right geometric elements which define the domains
// ellist : two elements : list of elements connected to the left and right indices respectively
// ellist : one element : skeleton is a boundary list contains the compelsides linked to the skeleton on the interior of the mech
void TPZMHMeshControl::ConnectedElements(int64_t skeleton, std::pair<int64_t,int64_t> &leftright, std::map<int64_t, std::list<TPZCompElSide> > &ellist)
{
    TPZGeoEl *gelskeleton = fGMesh->Element(skeleton);
    int meshdim = fGMesh->Dimension();
    if (gelskeleton->Dimension() != meshdim-1) {
        DebugStop();
    }
    TPZGeoElSide skelside(gelskeleton,gelskeleton->NSides()-1);

    std::list<TPZCompElSide> allconnected;
    BuildConnectedElementList(skelside, allconnected);
    
#ifdef PZDEBUG
    // verify if the acumulated area of allconnected is the area of the element
    REAL area = gelskeleton->SideArea(gelskeleton->NSides()-1);
    REAL listarea = 0.;
    for(auto it : allconnected)
    {
        TPZGeoElSide gside = it.Reference();
        listarea += gside.Area();
    }
    REAL diff = fabs(area-listarea);
    REAL diff2 = fabs(2*area - listarea);
    if(diff > 1.e-3*area && diff2 >= 1.e-3*area)
    {
        std::cout << "The area of the connected elements does not correspond to the area " <<
        "of the skeleton element area = " << area << " listarea = " << listarea <<
        " diff/area = " << diff/area << " diff2/area " << diff2/area << std::endl;
    }
#endif
    
    for(auto celside : allconnected)
    {
        TPZGeoElSide gelside = celside.Reference();
        TPZGeoEl *gel = gelside.Element();
        int64_t gelindex = gel->Index();
        int mhm_domain = fGeoToMHMDomain[gelindex];
        if(mhm_domain == leftright.first)
        {
            ellist[mhm_domain].push_back(celside);
        }
        else if(mhm_domain == leftright.second)
        {
            ellist[mhm_domain].push_back(celside);
        }
        else
        {
            std::cout << "I dont understand leftright : " << leftright.first << " " <<
            leftright.second << " mhm_domain " << mhm_domain << std::endl;
        }
    }
    // a skeleton element with no computational elements is a bug
    if (ellist.size() == 0) {
        DebugStop();
    }
    // if the skeleton is boundary and if there are more than one subdomain, it is a bug
    if (skeleton == leftright.second && ellist.size() != 1) {
        DebugStop();
    }
    // if the skeleton is not boundary then there must be exactly two subdomains connected
    if(skeleton != leftright.second && ellist.size() != 2)
    {
        DebugStop();
    }
}




/// identify connected elements to the skeleton elements
// the computational mesh is determined by the element pointed to by the geometric element
// skeleton index of the geometric element which defines the skeleton
// leftright indices of the left and right geometric elements which define the domains
// ellist : two elements : list of elements connected to the left and right indices respectively
// ellist : one element : skeleton is a boundary list contains the compelsides linked to the skeleton on the interior of the mech
void TPZMHMeshControl::ConnectedElements2(int64_t skeleton, std::pair<int64_t,int64_t> &leftright, std::map<int64_t, std::list<TPZCompElSide> > &ellist)
{
    TPZGeoEl *gelskeleton = fGMesh->Element(skeleton);
    int meshdim = fGMesh->Dimension();
    if (gelskeleton->Dimension() != meshdim-1) {
        DebugStop();
    }
    TPZGeoElSide gelside(gelskeleton,gelskeleton->NSides()-1);
    TPZGeoElSide neighbour(gelside.Neighbour());
    while (neighbour != gelside) {
        if (neighbour.Element()->MaterialId() == fSkeletonWrapMatId || neighbour.Element()->MaterialId() == fBoundaryWrapMatId) {
            break;
        }
        neighbour = neighbour.Neighbour();
    }
    if(neighbour == gelside)
    {
        std::cout << "For element " << gelskeleton->Index() << " with mat id " << gelskeleton->MaterialId() << " no skeleton found\n";
        neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            std::cout << "Neighbour matid " << neighbour.Element()->MaterialId() << std::endl;
            neighbour = neighbour.Neighbour();
        }
        DebugStop();
    }
    TPZGeoElSide skeletonwrap = neighbour;
    TPZStack<TPZGeoElSide> wrapstack;
    wrapstack.Push(skeletonwrap);
    while(wrapstack.size())
    {
        skeletonwrap = wrapstack.Pop();
        neighbour = skeletonwrap.Neighbour();
        while (neighbour != skeletonwrap)
        {
            if (neighbour.Element()->Dimension() == meshdim && neighbour.Element()->Reference()) {
                if (fGeoToMHMDomain[neighbour.Element()->Index()] == leftright.first)
                {
                    ellist[leftright.first].push_back(neighbour.Reference());
                }
                if (skeleton != leftright.second && fGeoToMHMDomain[neighbour.Element()->Index()] == leftright.second)
                {
                    ellist[leftright.second].push_back(neighbour.Reference());
                }
            }
            neighbour = neighbour.Neighbour();
        }
        if (skeletonwrap.HasSubElement()) {
            int nsub = skeletonwrap.Element()->NSubElements();
            for (int isub = 0; isub < nsub; isub++) {
                TPZGeoEl *subel = skeletonwrap.Element()->SubElement(isub);
                TPZGeoElSide subelside(subel,subel->NSides()-1);
                wrapstack.Push(subelside);
            }
        }
    }
    // a skeleton element with no computational elements is a bug
    if (ellist.size() == 0) {
        DebugStop();
    }
    // if the skeleton is boundary and if there are more than one subdomain, it is a bug
    if (skeleton == leftright.second && ellist.size() != 1) {
        DebugStop();
    }
    // if the skeleton is not boundary then there must be exactly two subdomains connected
    if(skeleton != leftright.second && ellist.size() != 2)
    {
        DebugStop();
    }
}

/// identify interface elements connected to the skeleton elements
// the computational mesh is determined by the element pointed to by the geometric element
void TPZMHMeshControl::ConnectedInterfaceElements(int64_t skeleton, std::pair<int64_t,int64_t> &leftright, std::map<int64_t, std::list<TPZInterfaceElement *> > &intfaces)
{
    std::map<int64_t, std::list<TPZCompElSide> > ellist;
    ConnectedElements(skeleton, leftright, ellist);
    // neighbour to each element there is an interface element
    for (std::map<int64_t, std::list<TPZCompElSide> >::iterator it = ellist.begin(); it != ellist.end(); it++) {
        std::list<TPZCompElSide> &loclist = it->second;
        for (std::list<TPZCompElSide>::iterator it2 = loclist.begin(); it2 != loclist.end(); it2++) {
            TPZCompElSide celside = *it2;
            TPZGeoElSide gelside = celside.Reference();
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                TPZCompElSide celneigh = neighbour.Reference();
                if (celneigh) {
                    TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(celneigh.Element());
                    if (intface && (intface->LeftElementSide() == celside || intface->RightElementSide() == celside)) {
                        intfaces[it->first].push_back(intface);
                        break;
                    }
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                DebugStop();
            }
        }
    }
}


/// Create the wrap elements
void TPZMHMeshControl::BuildWrapMesh()
{
    // all the elements should be neighbour of a wrap element
    int64_t nel = fGMesh->NElements();
    int dim = fGMesh->Dimension();
    fGeoToMHMDomain.Resize(nel+1000,-1);
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = fGMesh->Element(el);
        if(!gel || gel->Dimension() != dim) continue;
        int nsides = gel->NSides();
        int nfaces = gel->NSides(dim-1);
        for(int is = nsides-nfaces-1; is<nsides-1; is++)
        {
            TPZGeoElSide gelside(gel,is);
            if(gelside.HasNeighbour(fSkeletonMatId))
            {
                TPZGeoElBC gbc(gelside,fSkeletonWrapMatId);
                TPZGeoEl *gelcreate = gbc.CreatedElement();
                int64_t index = gelcreate->Index();
                if(index >= fGeoToMHMDomain.size()) fGeoToMHMDomain.Resize(index+1000, -1);
                fGeoToMHMDomain[index] = fGeoToMHMDomain[el];
            }
        }
    }
    fGeoToMHMDomain.Resize(fGMesh->NElements(),-1);
}
