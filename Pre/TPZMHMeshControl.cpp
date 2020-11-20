//
//  TPZMHMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/3/14.
//
//

#include "TPZMHMeshControl.h"
#include "TPZRefPattern.h"
#include "pzl2projection.h"
#include "TPZNullMaterial.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZCompElLagrange.h"
#include "pzbndcond.h"
#include "pzmat1dlin.h"
#include "TPZInterfaceEl.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzsubcmesh.h"
#include "TPZRefPatternTools.h"
#include "pzlog.h"
#include "TPZVTKGeoMesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "TPZGeoElSideAncestors.h"

#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mhmeshcontrol"));
static LoggerPtr loggerWRAP(Logger::getLogger("pz.mhmeshwrap"));
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
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
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
    fFluxMesh->SetName("SkeletonFluxMesh");
    fCMeshConstantPressure = new TPZCompMesh(fGMesh);
    fCMeshLagrange = new TPZCompMesh(fGMesh);
}

/// Define the partitioning information of the MHM mesh
void TPZMHMeshControl::DefinePartition(TPZVec<int64_t> &partitionindex, std::map<int64_t,std::pair<int64_t,int64_t> > &skeleton)
{
#ifdef PZDEBUG
    if (partitionindex.size() != fGMesh->NElements()) {
        DebugStop();
    }
#endif
    fGeoToMHMDomain = partitionindex;
    fInterfaces = skeleton;
    fMHMtoSubCMesh.clear();
    int64_t nc = partitionindex.size();
    for (int64_t c=0; c<nc; c++) {
        if(partitionindex[c] != -1)
        {
            fMHMtoSubCMesh[partitionindex[c] ] = -1;
        }
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
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
        // all volumetric elements need to belong to a partition
        if(gel->Dimension()==dim && partitionindex[iel] == -1)
        {
            DebugStop();
        }

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

    partitionindex.Resize(fGMesh->NElements(), -1);
    DefinePartition(partitionindex, skeleton);
}

TPZMHMeshControl::TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : fGMesh(gmesh), fGeoToMHMDomain(),
    fMHMtoSubCMesh(), fGlobalSystemSize(-1),
    fGlobalSystemWithLocalCondensationSize(-1),
    fNumeq(-1)
{

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
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
    fFluxMesh->SetName("SkeletonFluxMesh");
    fCMeshConstantPressure = new TPZCompMesh(fGMesh);
    fCMeshLagrange = new TPZCompMesh(fGMesh);
}

TPZMHMeshControl::TPZMHMeshControl(const TPZMHMeshControl &copy){

    this->operator=(copy);
}

TPZMHMeshControl &TPZMHMeshControl::operator=(const TPZMHMeshControl &cp){

    fGMesh = new TPZGeoMesh(*cp.fGMesh.operator->());
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
    // is this really desirable??? or should we make a "deep copy" (Philippe)
    // look at the inconsistency of the geometric mesh... the pressure mesh
    // will not point to the geometric mesh!!
    DebugStop();
    fPressureFineMesh = cp.fPressureFineMesh;
    fFluxMesh = cp.fFluxMesh;
    fCMeshConstantPressure = cp.fCMeshConstantPressure;
    fCMeshLagrange = fCMeshLagrange;
    fpOrderSkeleton = cp.fpOrderSkeleton;
    fpOrderInternal = cp.fpOrderInternal;
    fWrapMatId = cp.fWrapMatId;
    fBoundaryWrapMatId = cp.fBoundaryWrapMatId;
    fInternalWrapMatId = cp.fInternalWrapMatId;
    fHybridizeSkeleton = cp.fHybridizeSkeleton;
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
    fGeoToMHMDomain.Resize(nel);
    for (int64_t el=0; el<nel; el++) {
        fGeoToMHMDomain[el] = -1;
    }
    for (int64_t el=0; el<ncoarse; el++) {
        fGeoToMHMDomain[coarseindices[el]] = coarseindices[el];
        TPZGeoEl *gel = fGMesh->Element(coarseindices[el]);
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
    CreateSkeletonElements();
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
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


/// will create dim-1 elements on the interfaces between the coarse element indices
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
                fInterfaces[leftind]=std::make_pair(rightind, leftind);
                continue;
            }
            if (right->Reference()->Dimension() == dimension-1) {
                // the boundary element is the skeleton element
                fInterfaces[rightind]=std::make_pair(leftind, rightind);
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

//    BuildWrapMesh(fGMesh->Dimension());
//    BuildWrapMesh(fGMesh->Dimension()-1);

    fGeoToMHMDomain.Resize(fGMesh->NElements(), -1);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
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
    std::cout << "fWrapMatId "<<fWrapMatId<<std::endl;
}

/// divide the skeleton elements
void TPZMHMeshControl::DivideBoundarySkeletonElements()
{
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it;
    bool hasdivided = true;
    while (hasdivided)
    {
        hasdivided = false;
        // mapdivided will contain the new fInterface structure
        std::map<int64_t, std::pair<int64_t,int64_t> > mapdivided;
        for (it=fInterfaces.begin(); it!=fInterfaces.end(); it++) {
            int64_t elindex = it->first;
            // if the following condition is not satisfied then the interface is not a boundary
//            if (elindex != it->second.second) {
//                mapdivided[it->first] = it->second;
//                continue;
//            }
            TPZGeoEl *gel = fGMesh->Element(elindex);
            // if the geometric element associated with the interface was not divided
            // then the interface will not be divided either
            if(!gel->HasSubElement())
            {
                mapdivided[it->first] = it->second;
                continue;
            }
            hasdivided = true;
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
                // update the interface data structure
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
    TPZMaterial *meshmat = 0;
    while (it != matids.end()) {
        TPZL2Projection *material = new TPZL2Projection(*it,dim,nstate,sol);
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
        TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

        TPZMaterial * BCondD1 = meshmat->CreateBC(meshmat, *it,0, val1, val2);
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
        int64_t index;
        create.CreateCompEl(gel, *cmesh, index);
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
                        create.CreateCompEl(neighbour.Element(), *cmesh, index);
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
    CreateInternalElements();
    //    AddBoundaryElements();
    CreateSkeleton();
    if (fLagrangeAveragePressure) {
        // I want to change this mesh to consist of discontinuous elements with common connect
        this->CreateLagrangeMultiplierMesh();
    }
    
    #ifdef PZDEBUG
        if(1)
        {
            std::ofstream out1("HDivMesh.txt");
            fFluxMesh->Print(out1);
            std::ofstream out2("PressureMesh.txt");
            fPressureFineMesh->Print(out2);
            std::ofstream out5("DistributedFluxMHMesh.txt");
            fCMeshLagrange->Print(out5);
            std::ofstream out6("AverageSolutionMHMesh.txt");
            fCMeshConstantPressure->Print(out6);
            std::ofstream out7("MultiPhysicsBeforeCondensing.txt");
            fCMesh->Print(out7);

        }
    #endif

    BuildMultiPhysicsMesh();
    
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream out7("MultiPhysicsBeforeCondensing.txt");
        fCMesh->Print(out7);
    }
#endif

    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh, out);
    }
    
    CreateInterfaceElements();
//    AddBoundaryInterfaceElements();
    fCMesh->ExpandSolution();
    fCMesh->CleanUpUnconnectedNodes();


#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
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

    SetLagrangeMultiplierLevels();
    
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

    fPressureFineMesh->SetName("Pressure mesh");

    TPZCompMesh *cmesh = fPressureFineMesh.operator->();
    cmesh->ApproxSpace().SetCreateLagrange(false);
    cmesh->SetAllCreateFunctionsContinuous();

    //Criar elementos computacionais malha MHM

    TPZGeoEl *gsubel = NULL;
    fConnectToSubDomainIdentifier[cmesh].Expand(10000);

    int64_t nel = fGMesh->NElements();
    int meshdim = fGMesh->Dimension();
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
            if(gel->Dimension() == meshdim && fGeoToMHMDomain[el] == -1)
            {
                DebugStop();
            }
            int64_t index;
            // create the flux element
            cmesh->CreateCompEl(gel, index);
            TPZCompEl *cel = cmesh->Element(index);
            /// associate the connects with the subdomain
            SetSubdomain(cel, itMHM->first);
            // we need to create a lagrange multiplier element in order to delay decomposition of an equation
            if (!LagrangeCreated)
            {
                LagrangeCreated = true;
                TPZCompEl *cel = cmesh->ElementVec()[index];
                int64_t cindex = cel->ConnectIndex(0);
                int nshape(1), nvar(1), order(1);
                int lagrangelevel = 0;
                if (this->fLagrangeAveragePressure) {
                    lagrangelevel = 2;
                }
                else
                {
                    // this connect takes the job of the average pressure
                    lagrangelevel = 4;
                }
                cmesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                if (fProblemType == EElasticity2D) {
                    cindex = cel->ConnectIndex(1);
                    cmesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                }
                if (fProblemType == EElasticity3D) {
                    cindex = cel->ConnectIndex(2);
                    cmesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                }
            }
        }
        cmesh->ComputeNodElCon();
        // create wrap computational elements on the element/sides which have only one connectivty
        int64_t nelcreated = 0;
        int64_t nel_adjust = fGMesh->NElements();
        for (int64_t el=0; el< nel_adjust; el++)
        {
            TPZGeoEl *gel = fGMesh->Element(el);
            if(!gel) continue;
            // for all elements that have a referred element
            TPZCompEl *cel = gel->Reference();
            if(!cel) continue;
            if(gel->Dimension() != meshdim)
            {
                continue;
            }
            int nsides = gel->NSides();
            int firstside = nsides-gel->NSides(meshdim-1)-1;
            for(int side=firstside; side<nsides-1; side++)
            {
                TPZConnect &c = cel->Connect(side);
                // if the connect has dependency, then the side must be internal
                if(c.HasDependency()) continue;
                if(c.NElConnected() == 1)
                {
                    TPZGeoElSide gelside(gel,side);
                    TPZGeoElBC gbc(gelside,fWrapMatId);
                    int64_t index;
                    TPZCompEl *cel = cmesh->CreateCompEl(gbc.CreatedElement(), index);
                    SetSubdomain(cel, itMHM->first);
                    nelcreated++;
#ifdef PZDEBUG
                    TPZGeoElSideAncestors ancestor(gelside);
                    if(!ancestor.HasLargerorEqual(fSkeletonMatId) && !ancestor.HasLargerorEqual(fMaterialBCIds))
                    {
                        gel->Print(std::cout);
                        TPZGeoElSide neighbour(gelside.Neighbour());
                        while(neighbour != gelside)
                        {
                            std::cout << "matid " << neighbour.Element()->MaterialId();
                            neighbour = neighbour.Neighbour();
                        }
                        DebugStop();
                    }
#endif
                }
            }
        }
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "For subdomain " << itMHM->first << "nBoundary created " << nelcreated;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
    }
    fPressureFineMesh->ExpandSolution();
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fPressureFineMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fGMesh->ResetReference();
}

/// will create the elements on the skeleton
void TPZMHMeshControl::CreateSkeleton()
{
    // comment this line or not to switch the type of skeleton elements
    TPZCompMesh *cmesh = fFluxMesh.operator->();
    int meshdim = cmesh->Dimension();
    cmesh->SetDimModel(meshdim);
//    fCMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    int order = fpOrderSkeleton;
    if (order < 0) {
        DebugStop();
        order = 0;
    }
    cmesh->SetDefaultOrder(order);
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        int64_t elindex = it->first;
        // skip the boundary elements
//        if (elindex == it->second.second) {
//            it++;
//            continue;
//        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        int64_t index;
        // create a boundary flux element to model the flux
        
        TPZCompEl *cel = cmesh->CreateCompEl(gel, index);
        int Side = gel->NSides()-1;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            int lagrangelevel = 0;
            if (this->fLagrangeAveragePressure) {
                lagrangelevel = 3;
            }
            else
            {
                lagrangelevel = 3;
            }
            cel->Connect(ic).SetLagrangeMultiplier(lagrangelevel);
        }
        SetSubdomain(cel, -1);

        if (elindex == it->second.second) {
            // set the side orientation of the boundary elements
            intel->SetSideOrient(Side, 1);
            SetSubdomain(cel, it->second.first);
        }
        else
        {
            if (it->second.first < it->second.second) {
                // set the flux orientation depending on the relative value of the element ids
                intel->SetSideOrient(Side, 1);
            }
            else
            {
                intel->SetSideOrient(Side, -1);
            }
            SetSubdomain(cel, -1);
        }
        gel->ResetReference();
        it++;
    }
    fFluxMesh->ExpandSolution();
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fFluxMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// will create the interface elements between the internal elements and the skeleton
void TPZMHMeshControl::CreateInterfaceElements()
{
    if(fHybridizeSkeleton)
    {
        // include some more combinations in matmap
        // fSkeletonMatId -> fPressureSkeletonMatId
        // fSecondSkeletonMatId -> fPressureSkeletonMatId
        DebugStop();
    }
    fGMesh->ResetReference();
    fCMesh->LoadReferences();
#ifdef LOG4CXX
    {
        std::stringstream sout;
        for(auto gel: fGMesh->ElementVec())
        {
            int matid = gel->MaterialId();
            if(matid == fWrapMatId)
            {
                TPZManVector<REAL> xcenter(3);
                TPZGeoElSide gelside(gel);
                gelside.CenterX(xcenter);
                sout << "gel index " << gel->Index() << " matid = " << matid
                << " xcenter " << xcenter << " compel " << (void*) gel->Reference() << std::endl;
            }
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    std::set<int> AllMatids(fMaterialBCIds);
    AllMatids.insert(fSkeletonMatId);
    int64_t nel = fCMesh->NElements();
    for(int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = fCMesh->Element(el);
        if(!cel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        int matid = gel->MaterialId();
        if(matid != fWrapMatId) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSideAncestors ancestor(gelside);
        TPZGeoElSide largeside = ancestor.HasLargerorEqual(AllMatids);
        if(largeside && !largeside.Reference()) {
            largeside = ancestor.HasLarger(AllMatids);
        }
        // if the large element has no associated computational element, look for a larger element that has
        while(largeside && !largeside.Reference())
        {
            ancestor.SetCurrent(largeside);
        }
        int targetmatid = largeside.Element()->MaterialId();
        int lagrangematid = fLagrangeMatIdLeft;
        if(targetmatid == fSkeletonMatId)
        {
            // detect whether the geometric element is on the left or right side of the skeleton
            int MHMDomain = fGeoToMHMDomain[gel->Index()];
            if(MHMDomain == -1) DebugStop();
            int64_t skelindex = largeside.Element()->Index();
            // the indexes of the skeleton elements should be contained in the interfaces data structure
            auto skelit = fInterfaces.find(skelindex);
            if(skelit == fInterfaces.end()) DebugStop();
            // the geometric element should be belong to either the left of right MHM domain
            auto leftMHM = skelit->second.first;
            auto rightMHM = skelit->second.second;
            if(MHMDomain != leftMHM && MHMDomain != rightMHM) DebugStop();
            // if it belongs to the right domain, change the lagrange multiplier matid
            if(MHMDomain == rightMHM) lagrangematid = fLagrangeMatIdRight;
            if(fHybridizeSkeleton)
            {
                // point the large element to the second skeleton conditionally
                DebugStop();
            }
        }
        #ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    TPZManVector<REAL> xcenter(3);
                    gelside.CenterX(xcenter);
                    sout << "gel index " << gel->Index() << " matid = " << matid
                    << " xcenter " << xcenter << std::endl;
                    sout << "large matid " << targetmatid << " lagrange matid " << lagrangematid;
                    LOGPZ_DEBUG(logger, sout.str())
                }
        #endif

        // create the lagrange multiplier element
        TPZCompElSide celside(gelside.Reference());
        TPZCompElSide largecelside(largeside.Reference());
        if(!celside || !largecelside)
        {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelside,lagrangematid);
        int64_t index;
//        /** @brief Constructor */
//        TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right);
        
        new TPZMultiphysicsInterfaceElement(fCMesh,gelbc.CreatedElement(),index,celside,largecelside);
        
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
        TPZInterfaceElement *iface = dynamic_cast<TPZInterfaceElement *>(cel);
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
        TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!face) {
            DebugStop();
        }
        out << "Comp Element index " << face->LeftElement()->Index() << " side " << face->LeftElementSide().Side() << " skeleton element " << face->RightElement()->Index() << " interface index "
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
        TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!face) {
            DebugStop();
        }
        out << "Geom Element index " << face->LeftElement()->Reference()->Index() << " side " << face->LeftElementSide().Side() << " skeleton element " << face->RightElement()->Reference()->Index() << " interface index "
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

/// Add the boundary elements to the computational mesh
void TPZMHMeshControl::AddBoundaryElements()
{
    fGMesh->ResetReference();
    int dim = fGMesh->Dimension();
    std::set<int> notincluded;
    notincluded.insert(fSkeletonMatId);
    notincluded.insert(fLagrangeMatIdLeft);
    notincluded.insert(fLagrangeMatIdRight);
    int64_t nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->ElementVec()[el];
        if (!gel || gel->Dimension() == dim || gel->HasSubElement()) {
            continue;
        }
        int matid = gel->MaterialId();
        if (notincluded.find(matid) != notincluded.end()) {
            continue;
        }
        int64_t index;
        fCMesh->CreateCompEl(gel, index);
        gel->ResetReference();
    }
}

/// Add the boundary interface elements to the computational mesh
void TPZMHMeshControl::AddBoundaryInterfaceElements()
{
    fCMesh->LoadReferences();
    int dim = fGMesh->Dimension();
    std::set<int> notincluded;
    notincluded.insert(fSkeletonMatId);
    notincluded.insert(fLagrangeMatIdLeft);
    notincluded.insert(fLagrangeMatIdRight);
    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() == dim || gel->HasSubElement()) {
            continue;
        }
        int matid = gel->MaterialId();
        if (notincluded.find(matid) != notincluded.end()) {
            continue;
        }
        if (gel->Dimension() != dim-1) {
            DebugStop();
        }
        TPZStack<TPZCompElSide> celstack;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZCompElSide celside = gelside.Reference();
        gelside.EqualorHigherCompElementList2(celstack, 1, 0);
        int64_t nst = celstack.size();
        for (int64_t i=0; i<nst; i++) {
            TPZCompElSide cs = celstack[i];
            TPZGeoElSide gs = cs.Reference();
            if (gs == gelside) {
                continue;
            }
            int64_t index;
            TPZGeoEl *gelnew = gs.Element()->CreateBCGeoEl(gs.Side(), matid);
            new TPZInterfaceElement(fCMesh, gelnew , index, cs, celside);

        }
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

/// create the lagrange multiplier mesh, one element for each subdomain
void TPZMHMeshControl::CreateLagrangeMultiplierMesh()
{
    int dim = fGMesh->Dimension();
    fCMeshLagrange->SetAllCreateFunctionsDiscontinuous();
    fCMeshLagrange->SetDefaultOrder(0);
    fCMeshConstantPressure->SetAllCreateFunctionsDiscontinuous();
    fCMeshConstantPressure->SetDefaultOrder(0);
    fGMesh->ResetReference();
    int64_t connectcounter = 0;
    TPZCompElDisc::SetTotalOrderShape(fCMeshLagrange.operator->());
    int64_t nel = fPressureFineMesh->NElements();
    struct subdomain_data
    {
        int64_t fConnectIndex;
        TPZManVector<REAL,3> fCenter;
        REAL fTotalarea = 0.;
    };
    std::map<int,subdomain_data> MHMtoConnectMap;
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        int64_t gelindex = gel->Index();
        int subdomain = WhichSubdomain(cel);
        if (gel->Dimension() != dim || fMHMtoSubCMesh.find(subdomain) == fMHMtoSubCMesh.end()) {
            continue;
        }
        auto it = MHMtoConnectMap.find(subdomain);
        int64_t connectindex = -1;
        if(it == MHMtoConnectMap.end())
        {
            MHMtoConnectMap[subdomain].fConnectIndex = connectcounter;
            it = MHMtoConnectMap.find(subdomain);
            it->second.fCenter.Resize(3, 0.);
            it->second.fTotalarea = 0.;
            connectindex = connectcounter;
            connectcounter++;
        }
        else
        {
            connectindex = it->second.fConnectIndex;
        }
        // we will need to set the center point of all discontinuous element to the gravity
        // center of the subdomain
        TPZManVector<REAL,3> xcenter(3);
        TPZGeoElSide gelside(gel);
        gelside.CenterX(xcenter);
        REAL area = gelside.Area();
        for(int i=0; i<3; i++) it->second.fCenter[i] += xcenter[i]*area;
        it->second.fTotalarea += area;
        {
            int64_t index;
            TPZCompElDisc *disc = new TPZCompElDisc(fCMeshLagrange,gel,index);
            disc->SetTotalOrderShape();
            disc->SetFalseUseQsiEta();
            disc->SetConnectIndex(0, connectindex);
            SetSubdomain(disc, subdomain);
        }
        {
            int64_t index;
            TPZCompElDisc *disc = new TPZCompElDisc(fCMeshConstantPressure,gel,index);
            disc->SetTotalOrderShape();
            disc->SetFalseUseQsiEta();
            disc->SetConnectIndex(0, connectindex);
            SetSubdomain(disc, subdomain);
        }

//        fCMeshConstantStates->CreateCompEl(gel, index);
    }
    fCMeshLagrange->CleanUpUnconnectedNodes();
    fCMeshConstantPressure->CleanUpUnconnectedNodes();
    int64_t nc = fCMeshLagrange->NConnects();
    for(int64_t ic=0; ic<nc; ic++)
    {
        TPZConnect &c = fCMeshLagrange->ConnectVec()[ic];
        if(c.NElConnected() == 0)
        {
            SetSubdomain(fCMeshLagrange.operator->(), ic, -1);
            SetSubdomain(fCMeshConstantPressure.operator->(), ic,-1);
        }
    }

    for(auto it : MHMtoConnectMap)
    {
        for(int i=0;i<3;i++) it.second.fCenter[i] /= it.second.fTotalarea;
        int64_t cindex = it.second.fConnectIndex;
        fCMeshLagrange->ConnectVec()[cindex].SetLagrangeMultiplier(1);
        fCMeshConstantPressure->ConnectVec()[cindex].SetLagrangeMultiplier(4);
    }
    int64_t nel_constant = fCMeshConstantPressure->NElements();
    for (int64_t el=0; el<nel_constant; el++) {
        TPZCompEl *cel = fCMeshConstantPressure->ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        int subdomain = WhichSubdomain(cel);
        int64_t gelindex = gel->Index();
        if (gel->Dimension() != dim || fMHMtoSubCMesh.find(subdomain) == fMHMtoSubCMesh.end()) {
            DebugStop();
        }
        auto it = MHMtoConnectMap.find(subdomain);
        if(it == MHMtoConnectMap.end()) DebugStop();
        // impose the center point on the discontinuous elements
        {
            TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
            if(!disc) DebugStop();
            disc->SetCenterPoint(it->second.fCenter);
        }
        {
            TPZCompEl *cel = fCMeshLagrange->Element(el);
            TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
            if(!disc) DebugStop();
            disc->SetCenterPoint(it->second.fCenter);
        }
    }
    fCMeshLagrange->ExpandSolution();
    fCMeshConstantPressure->ExpandSolution();
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        fCMeshLagrange->Print(sout);
        fCMeshConstantPressure->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#ifdef PZDEBUG
    {
        int64_t numc = fConnectToSubDomainIdentifier[fCMeshConstantPressure.operator->()].size();
        int64_t numc2 = fCMeshConstantPressure->NConnects();
        if(numc != numc2) DebugStop();
    }
    {
        int64_t numc = fConnectToSubDomainIdentifier[fCMeshLagrange.operator->()].size();
        int64_t numc2 = fCMeshLagrange->NConnects();
        if(numc != numc2) DebugStop();
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
            int64_t index;
            new TPZCompElLagrange(fCMesh,*cel,index);
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        int64_t index;
        fCMesh->CreateCompEl(gel, index);
        TPZCompEl *celnew = fCMesh->ElementVec()[index];
        TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(celnew);
        if (!mult) {
            DebugStop();
        }
        mult->AddElement(cel, 0);
    }
    fCMesh->LoadReferences();
    nel = fCMeshLagrange->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMeshLagrange->ElementVec()[el];
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
    nel = fCMeshConstantPressure->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMeshConstantPressure->ElementVec()[el];
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
    cmeshvec[1] = fCMeshLagrange.operator->();
    cmeshvec[2] = fCMeshConstantPressure.operator->();
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
        int64_t index;
        new TPZMultiphysicsInterfaceElement(fCMesh,intface->Reference(),index,multleft,multright);

    }


    nel = fCMeshConstantPressure->NElements();
    int64_t npressconnect = fPressureFineMesh->NConnects();
    int64_t nlagrangeconnect = fCMeshLagrange->NConnects();
    // nel numero de dominios MHM, tem um connect associado a cada um e os mesmos estao no final
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = this->fCMeshConstantPressure->Element(el);
        int64_t pressureconnect = cel->ConnectIndex(0);
        int64_t cindex = npressconnect+nlagrangeconnect+pressureconnect;
        fCMesh->ConnectVec()[cindex].SetLagrangeMultiplier(3);
    }
    fCMesh->ExpandSolution();
    //fCMesh->SaddlePermute();
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
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

    // create the submeshes
    for (auto it : fMHMtoSubCMesh) {
        int64_t index;
        TPZSubCompMesh *submesh = new TPZSubCompMesh(fCMesh,index);
        submeshes[it.first] = submesh;
    }
    for (auto it : submeshes) {
        fMHMtoSubCMesh[it.first] = it.second->Index();
    }

    fGMesh->ResetReference();
    TPZCompMesh *multiphys = fCMesh.operator->();
    multiphys->LoadReferences();

    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = multiphys->Element(el);
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
        submeshes[domain]->TransferElement(multiphys, cel->Index());
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
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
    multiphys->ComputeNodElCon();


    for (auto itsub : submeshes) {
        TPZSubCompMesh *submesh = itsub.second;
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
            if (fHybridizeSkeleton) {
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
                std::cout << "For subdomain " << itsub.first << " connect index " << connectindex << " left external as lagrange multiplier\n";
                submesh->ExpandSolution();
                c.Print(*multiphys);
#endif
            }
        }
        submesh->MakeAllInternal();
//        for (std::set<int64_t>::iterator it = internals.begin(); it != internals.end(); it++) {
//            submesh->MakeInternal(*it);
//        }
        submesh->InitializeBlock();
    }
    fCMesh->CleanUpUnconnectedNodes();
    for (auto itsub : submeshes) {
        TPZSubCompMesh *submesh = itsub.second;
        int numthreads = 0;
        int preconditioned = 0;
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "Submesh for MHM domain " << itsub.first << "\n";
            submesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZAutoPointer<TPZGuiInterface> guiInterface;
        submesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
    }

    multiphys->SaddlePermute();
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
    if (fCMeshLagrange)
    {
        out << "******************* LAGRANGE MULTIPLIER MESH *****************\n";
        fCMeshLagrange->Print(out);
    }

    /// computational mesh to represent the constant states
    if (fCMeshConstantPressure)
    {
        out << "******************* CONSTANTE PRESSURE MESH *****************\n";
        fCMeshConstantPressure->Print(out);
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
#ifdef PZDEBUG
    if(fMaterialIds.size() == 0)
    {
        std::cout << __PRETTY_FUNCTION__ << " you forgot to initialize the fMaterialIds data structure\n";
        DebugStop();
    }
#endif
    for(auto matid : fMaterialIds)
    {
        TPZMaterial *mat = fCMesh->FindMaterial(matid);
        if (!mat) {
            DebugStop();
        }
    }
    for(auto matid : fMaterialBCIds)
    {
        TPZMaterial *mat = fCMesh->FindMaterial(matid);
        if (!mat) {
            DebugStop();
        }
    }

    int dim = fCMesh->Dimension();
    TPZNullMaterial *matPerif = NULL;

    if (fCMesh->FindMaterial(fSkeletonMatId)) {
        DebugStop();
    }
    matPerif = new TPZNullMaterial(fSkeletonMatId);
    matPerif->SetDimension(dim-1);
    matPerif->SetNStateVariables(fNState);
    fCMesh->InsertMaterialObject(matPerif);

    if (fHybridizeSkeleton) {
        if (fCMesh->FindMaterial(fPressureSkeletonMatId)) {
            DebugStop();
        }
        matPerif = new TPZNullMaterial(fPressureSkeletonMatId);
        matPerif->SetDimension(dim-1);
        matPerif->SetNStateVariables(fNState);
        fCMesh->InsertMaterialObject(matPerif);

        if (fCMesh->FindMaterial(fSecondSkeletonMatId)) {
            DebugStop();
        }
        matPerif = new TPZNullMaterial(fSecondSkeletonMatId);
        matPerif->SetDimension(dim-1);
        matPerif->SetNStateVariables(fNState);
        fCMesh->InsertMaterialObject(matPerif);
    }
    // insert the lagrange material ids
    {
        int LagrangeMatIdLeft = 50;
        int LagrangeMatIdRight = 51;
        int nstate = fNState;
        int dim = fGMesh->Dimension();

        if (fCMesh->FindMaterial(fLagrangeMatIdLeft)) {
            DebugStop();
        }
        if (fCMesh->FindMaterial(fLagrangeMatIdRight)) {
            DebugStop();
        }
        TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(fLagrangeMatIdLeft,dim,nstate);
        TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(fLagrangeMatIdRight,dim,nstate);
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
    // insert the wrap material ids
    {
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(fWrapMatId);
            nullmat->SetNStateVariables(fNState);
            nullmat->SetDimension(dim-1);
            fCMesh->InsertMaterialObject(nullmat);
        }
    }

    // insert material objects in the pressure mesh
    {
        TPZCompMesh *cmesh = fPressureFineMesh.operator->();
        for(auto it : fMaterialIds)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(it);
            nullmat->SetDimension(dim);
            nullmat->SetNStateVariables(fNState);
            cmesh->InsertMaterialObject(nullmat);
        }
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(fWrapMatId);
            nullmat->SetNStateVariables(fNState);
            nullmat->SetDimension(dim-1);
            cmesh->InsertMaterialObject(nullmat);
        }
    }
    // insert material objects in the flux mesh
    {
        TPZCompMesh *cmesh = fFluxMesh.operator->();
        for(auto it : fMaterialBCIds)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(it);
            nullmat->SetDimension(dim);
            nullmat->SetNStateVariables(fNState);
            cmesh->InsertMaterialObject(nullmat);
        }
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(fSkeletonMatId);
            nullmat->SetNStateVariables(fNState);
            nullmat->SetDimension(dim-1);
            cmesh->InsertMaterialObject(nullmat);
        }
    }
    {
        // insert the material in the constant pressure mesh and lagrange mesh
        TPZCompMesh *cmesh = fCMeshConstantPressure.operator->();
        int nStateVariables = 1;
        if (fProblemType == EScalar) {
            nStateVariables = 1;
        }
        else if(fProblemType == EElasticity2D)
        {
            nStateVariables = 3;
        }
        else if(fProblemType == EElasticity3D)
        {
            nStateVariables = 6;
        }
        for(auto it : fMaterialIds)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(it);
            nullmat->SetDimension(dim);
            nullmat->SetNStateVariables(nStateVariables);
            cmesh->InsertMaterialObject(nullmat);
        }
        cmesh = fCMeshLagrange.operator->();
        for(auto it : fMaterialIds)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(it);
            nullmat->SetDimension(dim);
            nullmat->SetNStateVariables(nStateVariables);
            cmesh->InsertMaterialObject(nullmat);
        }
    }
}

void TPZMHMeshControl::HybridizeSkeleton(int skeletonmatid, int pressurematid)
{
    // comment this line or not to switch the type of skeleton elements
    int meshdim = fCMesh->Dimension();
//    fCMesh->SetDimModel(meshdim-1);
//    fCMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    fCMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    fCMesh->ApproxSpace().CreateDisconnectedElements(true);
    int order = fpOrderSkeleton;
    if (order < 1) {
        DebugStop();
        order = 1;
    }
    TPZCompMesh *cmesh = fCMesh.operator->();
    // expand the vector
    fConnectToSubDomainIdentifier[cmesh].Resize(cmesh->NConnects()+10000, -1);
    cmesh->SetDefaultOrder(order);
    int64_t nelem = fGMesh->NElements();
    fCMesh->LoadReferences();
    // keep the pointers of the skeleton elements we will work
    TPZManVector<TPZCompEl *> skeletoncomp(nelem,0);
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it;
    // loop over the skeleton elements
    for (it=fInterfaces.begin(); it != fInterfaces.end(); it++) {
        int64_t elindex = it->first;
        // skip the boundary interfaces
        if (it->first == it->second.second) {
            continue;
        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        // skip the boundary elements
        if (gel->MaterialId() < 0) {
            continue;
        }
        TPZCompEl *cel = gel->Reference();
        if (!cel) {
            DebugStop();
        }
        skeletoncomp[elindex] = cel;
    }
    fGMesh->ResetReference();
    for (it=fInterfaces.begin(); it != fInterfaces.end(); it++) {
        int64_t elindex = it->first;
        // skip the boundary elements
        if (elindex == it->second.second) {
            continue;
        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        TPZCompEl *celskeleton = skeletoncomp[elindex];
        // skip the boundary elements
        if (gel->MaterialId() < 0) {
            DebugStop();
            continue;
        }
        int64_t index1,index2;
        TPZManVector<TPZGeoEl *,4> GelVec(3);
        GelVec[0] = gel;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        // create geometric elements aint64_t the skeleton element
        TPZGeoElBC gbc1(gelside,fSecondSkeletonMatId);
        TPZGeoElBC gbc2(gelside,fPressureSkeletonMatId);
        GelVec[1] = gbc1.CreatedElement();
        GelVec[2] = gbc2.CreatedElement();

        fCMesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
        fCMesh->ApproxSpace().CreateDisconnectedElements(true);
        // create a discontinuous element to model the flux
        // index1 is the new flux element
        fCMesh->CreateCompEl(GelVec[1], index1);
        GelVec[1]->ResetReference();
        // index 2 is the new pressure element
        fCMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        fCMesh->ApproxSpace().CreateDisconnectedElements(true);
        fCMesh->CreateCompEl(GelVec[2], index2);
        GelVec[2]->ResetReference();

        SetSubdomain(fCMesh->Element(index1), -1);
        SetSubdomain(fCMesh->Element(index2), -1);

        // swap the elements. The skeleton element is now a pressure element
        // the former skeleton element is now a newly created geometric element
        celskeleton->SetReference(GelVec[2]->Index());
        fCMesh->Element(index2)->SetReference(gel->Index());
        GelVec[2]->SetMaterialId(fSkeletonMatId);
        gel->SetMaterialId(fPressureSkeletonMatId);


    }
    // the flux elements and pressure elements have been created, now generate the interface elements
    // the existing interface elements will be adjusted
    // two new interfaces need to be created
    fCMesh->LoadReferences();
    for(it=fInterfaces.begin(); it != fInterfaces.end(); it++) {
        int64_t elindex = it->first;
        // skip the boundary elements
        if (elindex == it->second.second) {
            continue;
        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        // skip the boundary elements
        if (gel->MaterialId() < 0) {
            DebugStop();
            continue;
        }
        // the element should be a pressure element
        if (!gel->Reference() || gel->MaterialId() != fPressureSkeletonMatId) {
            DebugStop();
        }
        // identify the pressure element
        TPZCompEl *celpressure = gel->Reference();
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        // look for the flux elements
        TPZCompEl *celskeleton = 0, *celsecondskeleton = 0;

        // find all elements connected to the pressure element
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            if (neighbour.Element()->MaterialId() == fSkeletonMatId) {
                celskeleton = neighbour.Element()->Reference();
            }
            if (neighbour.Element()->MaterialId() == fSecondSkeletonMatId) {
                celsecondskeleton = neighbour.Element()->Reference();
            }
            neighbour = neighbour.Neighbour();
        }

        if (!celpressure || !celskeleton || !celsecondskeleton) {
            DebugStop();
        }

        std::map<int64_t, std::list<TPZInterfaceElement *> > interfaces;
        ConnectedInterfaceElements(it->first, it->second, interfaces);

#ifdef PZDEBUG
        {
            std::stringstream sout;
            sout << "For skeleton element gelindex " << gel->Index() << " area " << gel->SideArea(gel->NSides()-1) << " we found interfaces\n";
            for (std::map<int64_t, std::list<TPZInterfaceElement *> >::iterator it = interfaces.begin(); it != interfaces.end(); it++) {
                sout << "domain " << it->first << std::endl;
                for (std::list<TPZInterfaceElement *>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                    TPZGeoEl *g = (*it2)->Reference();
                    int ns = g->NSides();
                    sout << "intface gelindex " << (*it2)->Reference()->Index() << " celindex " << (*it2)->Index() << " area " << g->SideArea(ns-1) << std::endl;
                }
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

        // find all interface elements (and others) connected to the pressure element
        for (std::map<int64_t, std::list<TPZInterfaceElement *> >::iterator it = interfaces.begin(); it != interfaces.end(); it++) {
            for (std::list<TPZInterfaceElement *>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                TPZInterfaceElement *intface = *it2;
                int matid = intface->Reference()->MaterialId();
                if (matid != fLagrangeMatIdLeft && matid != fLagrangeMatIdRight) {
                    DebugStop();
                }
                TPZCompElSide right = intface->RightElementSide();
                TPZCompElSide left = intface->LeftElementSide();
                if (right.Element() != celskeleton) {
                    DebugStop();
                }
                // switch the interface element with this material id to the newly created flux element
                if (matid == fLagrangeMatIdRight) {
                    TPZCompElSide cside(celsecondskeleton,gel->NSides()-1);
                    intface->SetLeftRightElements(left, cside);
                    right = cside;
                }

                int64_t leftcindex = left.Element()->ConnectIndex(0);
                int subdomain = fConnectToSubDomainIdentifier[cmesh][leftcindex];
                if (subdomain == -1) {
                    DebugStop();
                }
                fGeoToMHMDomain[right.Element()->Reference()->Index()] = subdomain;
                SetSubdomain(right.Element(), subdomain);
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Interface right flux element " << right.Element()->Index() << " has been adjusted subdomain is " << WhichSubdomain(right.Element());
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
            }
        }

#ifdef PZDEBUG
        {
            std::stringstream sout;
            sout << "For skeleton element gelindex " << gel->Index() << " area " << gel->SideArea(gel->NSides()-1) << " we found interfaces\n";
            for (std::map<int64_t, std::list<TPZInterfaceElement *> >::iterator it = interfaces.begin(); it != interfaces.end(); it++) {
                sout << "domain " << it->first << std::endl;
                for (std::list<TPZInterfaceElement *>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                    sout << "intface gelindex " << (*it2)->Reference()->Index() << " celindex " << (*it2)->Index() << " subdomain " << WhichSubdomain(*it2) << std::endl;
                }
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif


        // create the interfaces between the flux elements and the newly created pressure element
        TPZCompElSide leftflux(celskeleton,gel->NSides()-1);
        TPZCompElSide rightflux(celsecondskeleton,gel->NSides()-1);
        TPZCompElSide pressure(celpressure,gel->NSides()-1);
        TPZGeoElBC gbc3(gelside,fLagrangeMatIdRight);
        TPZGeoElBC gbc4(gelside,fLagrangeMatIdLeft);
        int64_t index3, index4;
        TPZInterfaceElement *intfaceleft = new TPZInterfaceElement(fCMesh,gbc3.CreatedElement(),index3);
        TPZInterfaceElement *intfaceright = new TPZInterfaceElement(fCMesh,gbc4.CreatedElement(),index4);
        intfaceleft->SetLeftRightElements(leftflux, pressure);
        intfaceright->SetLeftRightElements(rightflux, pressure);
#ifdef PZDEBUG
        {
            std::stringstream sout;
            sout << "Interfaces created by hybridization interface left subdomain el index " << intfaceleft->Index() << " dom " << WhichSubdomain(intfaceleft) <<
            " el index " << intfaceright->Index() << " interface right subdomain " << WhichSubdomain(intfaceright);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // adjust the lagrange level of the flux and pressure connects
        // pressure element
        int nc = celpressure->NConnects();
        for (int ic=0; ic<nc; ic++) {
            int lagrangelevel = 4;
            celpressure->Connect(ic).SetLagrangeMultiplier(lagrangelevel);
        }
        // skeleton element
        nc = celskeleton->NConnects();
        for (int ic=0; ic<nc; ic++) {
            int lagrangelevel = 2;
            celskeleton->Connect(ic).SetLagrangeMultiplier(lagrangelevel);
        }
        nc = celsecondskeleton->NConnects();
        for (int ic=0; ic<nc; ic++) {
            int lagrangelevel = 2;
            celsecondskeleton->Connect(ic).SetLagrangeMultiplier(lagrangelevel);
        }
    }
    fCMesh->ExpandSolution();
    fCMesh->SetDimModel(meshdim);
    fConnectToSubDomainIdentifier[fCMesh.operator->()].Resize(fCMesh->NConnects(), -1);

    // the new geometric elements are not internal
    fGeoToMHMDomain.Resize(fGMesh->NElements(), -1);
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
    }
#ifdef PZDEBUG
    if(fGeoToMHMDomain[index] != -1 && subdomain != fGeoToMHMDomain[index])
    {
        DebugStop();
    }
#endif
    fGeoToMHMDomain[index] = subdomain;
}

/// associates a geometric element with a subdomain
void TPZMHMeshControl::SetSubdomain(TPZGeoEl *gel, int64_t subdomain)
{
    int64_t index = gel->Index();

    if (index >= fGeoToMHMDomain.size()) {
        fGeoToMHMDomain.Resize(index+1, -1);
    }
    fGeoToMHMDomain[index] = subdomain;
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
    if(fConnectToSubDomainIdentifier.find(cmesh) == fConnectToSubDomainIdentifier.end())
    {
        DebugStop();
    }
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
        int64_t numc = fConnectToSubDomainIdentifier[meshvec[im]].size();
        int64_t numc2 = meshvec[im]->NConnects();
        if(numc != numc2)
        {
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
        if (neighbour.Element()->MaterialId() == fWrapMatId) {
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
void TPZMHMeshControl::BuildWrapMesh(int dim)
{
    // all the elements should be neighbour of a wrap element
    int64_t nel = fGMesh->NElements();
  
//    int meshdim = fGMesh->Dimension();
// first create the neighbours of boundary elements (works for dim = fGMesh->Dimension()-2)
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->MaterialId() == fBoundaryWrapMatId) {
                    int wrapmat = HasWrapNeighbour(gelside);
                    if (wrapmat && wrapmat != fBoundaryWrapMatId) {
                        DebugStop();
                    }
                    else
                    {
                        CreateWrap(gelside,fBoundaryWrapMatId);
                    }
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
    // create the wrap elements around the skeleton elements
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                int neighmatid = neighbour.Element()->MaterialId();
                if (IsSkeletonMatid(neighmatid)) {
                    int wrapmat = HasWrapNeighbour(gelside);
                    if (wrapmat && wrapmat != fWrapMatId) {
                        std::cout << "neighbour is skeleton with matid " << neighmatid <<
                        " but the element has a wrap neighbour with matid " << wrapmat << std::endl;
                        gelside.Print(std::cout);
                        std::cout << std::endl;
                        gelside.Element()->Print(std::cout);
                        std::ofstream out("gmesh.txt");
                        fGMesh->Print(out);
                        std::ofstream outvtk("gmesh.vtk");
                        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh,outvtk);
                        DebugStop();
                    }
                    else if(!wrapmat)
                    {
                        CreateWrap(neighbour,fWrapMatId);
                    }
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            if (!HasWrapNeighbour(gelside)) {
                CreateWrap(gelside);
            }
        }
    }
}

/// Verify if the element side contains a wrap neighbour
int TPZMHMeshControl::HasWrapNeighbour(TPZGeoElSide gelside)
{

    int element_dimension = gelside.Dimension();
    std::set<int> wrapmat;
    wrapmat.insert(fWrapMatId);
    wrapmat.insert(fInternalWrapMatId);
#ifdef PZDEBUG
    {
        // gelside cannot have a material of type wrap
        if(gelside.Element()->Dimension() == element_dimension && wrapmat.find(gelside.Element()->MaterialId()) != wrapmat.end())
        {
            DebugStop();
        }
    }
#endif
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        if (neighbour.Element()->Dimension() == element_dimension && wrapmat.find(neighbour.Element()->MaterialId()) != wrapmat.end()) {
            return neighbour.Element()->MaterialId();
        }
        neighbour = neighbour.Neighbour();
    }
    return 0;
}

/// Verify if the mesh datastructure is consistent
/// if the current gelside has no father, then none of its neighbours should either
void TPZMHMeshControl::CheckDivisionConsistency(TPZGeoElSide gelside)
{
#ifdef PZDEBUG
    TPZGeoElSide father = gelside.StrictFather();
    if (father && father.Dimension() == gelside.Dimension()) {
        DebugStop();
    }
#endif
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        TPZGeoElSide neighfather = neighbour.StrictFather();
        if (neighfather && neighfather.Element()->Dimension() == gelside.Dimension()) {
            std::cout << "neighfather dimension " << neighfather.Element()->Dimension() << " my dimension " << gelside.Dimension() << std::endl;
            DebugStop();
        }
        neighbour = neighbour.Neighbour();
    }
}


/// Return the wrap material id (depends on being boundary, neighbour of skeleton or interior
int TPZMHMeshControl::WrapMaterialId(TPZGeoElSide gelside)
{
    TPZGeoElSide father = gelside.StrictFather();
    while(father && father.Dimension() == gelside.Dimension())
    {
        gelside = gelside.StrictFather();
        father = gelside.StrictFather();
    }
    int meshdim = gelside.Element()->Mesh()->Dimension();
//#ifdef PZDEBUG
//    if (gelside.Dimension() != meshdim-1) {
//        DebugStop();
//    }
//#endif
    int wrapmat = HasWrapNeighbour(gelside);
    if(wrapmat)
    {
        return wrapmat;
    }
    if (gelside.Dimension() == meshdim-1)
    {
        int hasDimNeighbour = 0;
        TPZGeoElSide neighbour = gelside.Neighbour();
        bool isBoundary = false;
        while (neighbour != gelside) {
            if (neighbour.Element()->MaterialId() <0) {
                isBoundary = true;
            }
            if (neighbour.Element()->Dimension() == meshdim && !isBoundary) {
                hasDimNeighbour = 1;
            }
            if (IsSkeletonMatid(neighbour.Element()->MaterialId())) {
                return fWrapMatId;
            }
            neighbour = neighbour.Neighbour();
        }
        if (hasDimNeighbour) {
            return fInternalWrapMatId;
        }
        return fBoundaryWrapMatId;
    }
    else if(gelside.Dimension() == meshdim-2)
    {
        TPZGeoElSide neighbour = gelside.Neighbour();
        int hasBoundaryNeighbour = 0;
        while (neighbour != gelside) {
            if (IsSkeletonMatid(neighbour.Element()->MaterialId())) {
                return fWrapMatId;
            }
            if (neighbour.Element()->MaterialId() == fBoundaryWrapMatId)
            {
                hasBoundaryNeighbour = 1;
            }
            neighbour = neighbour.Neighbour();
        }
        if (hasBoundaryNeighbour) {
            return fBoundaryWrapMatId;
        }
        return fInternalWrapMatId;
    }
    else
    {
        DebugStop();
    }
    return wrapmat;
}

/// CreateWrapMesh of a given material id
void TPZMHMeshControl::CreateWrap(TPZGeoElSide gelside)
{
    TPZGeoElSide father = gelside.StrictFather();
    while(father && father.Dimension() == gelside.Dimension())
    {
        gelside = gelside.StrictFather();
        father = gelside.StrictFather();
    }
#ifdef PZDEBUG
//    CheckDivisionConsistency(gelside);
#endif
    int wrapmat = HasWrapNeighbour(gelside);
    if(!wrapmat)
    {
        wrapmat = WrapMaterialId(gelside);
        TPZGeoElBC gbc(gelside, wrapmat);
#ifdef LOG4CXX
        if(loggerWRAP->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Creating a wrap element on element/side " << gelside.Element()->Index() <<
            " " << gelside.Side() << " with wrapmat " << wrapmat << " el index " << gbc.CreatedElement()->Index();
            LOGPZ_DEBUG(loggerWRAP, sout.str())
        }
#endif
        DivideWrap(gbc.CreatedElement());
    }
    else
    {
        TPZGeoElSide neighbour = gelside;
        while (neighbour != gelside) {
            if (neighbour.Element()->MaterialId() == wrapmat) {
                DivideWrap(neighbour.Element());
                break;
            }
            neighbour = neighbour.Neighbour();
        }
    }
}

/// CreateWrapMesh of a given material id
void TPZMHMeshControl::CreateWrap(TPZGeoElSide gelside, int wrapmaterial)
{
    TPZGeoElSide father = gelside.StrictFather();
    while(father && father.Dimension() == gelside.Dimension())
    {
        gelside = gelside.StrictFather();
        father = gelside.StrictFather();
    }
#ifdef PZDEBUG
    if(gelside.Element()->Mesh()->Dimension() == 2)
    {
        CheckDivisionConsistency(gelside);
    }
#endif
    int wrapmat = HasWrapNeighbour(gelside);
    if(wrapmat == 0)
    {
        wrapmat = wrapmaterial;
        TPZGeoElBC gbc(gelside, wrapmat);
        DivideWrap(gbc.CreatedElement());
    }
    else if(wrapmaterial != wrapmat)
    {
        DebugStop();
    }
    else
    {
        TPZGeoElSide neighbour = gelside;
        while (neighbour != gelside) {
            if (neighbour.Element()->MaterialId() == wrapmat) {
                DivideWrap(neighbour.Element());
                break;
            }
            neighbour = neighbour.Neighbour();
        }
    }

}


/// Recursively divide the wrap element while it has divided neighbours
void TPZMHMeshControl::DivideWrap(TPZGeoEl *wrapelement)
{
    TPZGeoElSide gelside(wrapelement, wrapelement->NSides()-1);
    int dim = gelside.Dimension();
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        if(neighbour.Element()->HasSubElement() && neighbour.NSubElements() > 1)
        {
            // verify if the subelements are connected to a wrap element
            TPZStack<TPZGeoElSide> subelsides;
            neighbour.GetSubElements2(subelsides);
            int nsub = subelsides.size();
            int nwrap = 0;
            int nsidesdim = 0;
            for(int isub = 0; isub<nsub; isub++)
            {
                if(subelsides[isub].Dimension() == dim)
                {
                    nsidesdim++;
                }
                else continue;
                if(HasWrapNeighbour(subelsides[isub])) nwrap++;
            }
            if(nwrap != 0 && nwrap != nsidesdim)
            {
                std::cout << "I don't understand\n";
                DebugStop();
            }
            if(nwrap != 0) return;
            if (neighbour.Side() == neighbour.NSides()-1) {
                TPZAutoPointer<TPZRefPattern> siderefpattern = neighbour.Element()->GetRefPattern();
                if(!siderefpattern)
                {
                    std::cout << __PRETTY_FUNCTION__ << " We expect elements to have refinement patterns\n";
//                    DebugStop();
                }
                else
                {
                    wrapelement->SetRefPattern(siderefpattern);
                }
            }
            else
            {
                TPZAutoPointer<TPZRefPattern> elrefpattern = neighbour.Element()->GetRefPattern();
                if(!elrefpattern)
                {
                    std::cout << __PRETTY_FUNCTION__ << " We expect elements to have refinement patterns\n";
//                    DebugStop();
                }
                else
                {
                    TPZAutoPointer<TPZRefPattern> siderefpattern = elrefpattern->SideRefPattern(neighbour.Side());
                    if(!siderefpattern) DebugStop();
                    wrapelement->SetRefPattern(siderefpattern);
                }
            }
            TPZStack<TPZGeoEl *> subels;

            wrapelement->Divide(subels);
#ifdef PZDEBUG
            if(subels.size() == 1)
            {
                DebugStop();
            }
#endif
#ifdef LOG4CXX
            if(loggerWRAP->isDebugEnabled())
            {
                std::stringstream sout;
                sout << std::endl;
                for (int is=0; is<subels.size(); is++) {
                    sout << "Created a son wrap element of " << wrapelement->Index() <<
                    " matid " << wrapelement->MaterialId() << " with index " <<
                    subels[is]->Index() << std::endl;
                }
                LOGPZ_DEBUG(loggerWRAP, sout.str())
            }
#endif
            for (int is=0; is<subels.size(); is++) {
                DivideWrap(subels[is]);
            }
            break;
        }
        neighbour = neighbour.Neighbour();
    }
}

static void PrintConnect(TPZConnect &c, std::ostream &out)
{
    out << "Order " << (int)c.Order() << " NShape " << c.NShape() << " Lagrange " << (int)c.LagrangeMultiplier() << " condensed " << c.IsCondensed();
}

static void PrintGeo(TPZGeoEl *gel, std::ostream &out)
{
    TPZManVector<REAL,3> xi(gel->Dimension()),x(3);
    gel->CenterPoint(gel->NSides()-1, xi);
    gel->X(xi, x);
    out << " center " << x << " matid " << gel->MaterialId();
}

static void PrintElementInfo(TPZCompEl *cel, std::ostream &out);

/// print info on a condensed computational element
static void PrintCondensedElementInfo(TPZCondensedCompEl *condense, std::ostream &out)
{
    TPZCompEl *compel = condense->ReferenceCompEl();
    out << "&&&&& Condensed element index " << condense->Index() << std::endl;
    PrintElementInfo(compel,out);
    out << "End Condensed element &&&&&\n";
}


/// print subdomain info, focussing on condensed elements, geometric location of Lagrange level of connects
static void PrintSubDomainInfo(TPZSubCompMesh *submesh, std::ostream &out)
{
    out << "***** Sub domain info for domain index " << submesh->Index() <<
    "*****" <<std::endl;
    out << "External connects : \n";
    int nc = submesh->NConnects();
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = submesh->Connect(ic);
        if(c.NElConnected()) PrintConnect(c,out); out << std::endl;
    }
    out << "Elements\n";
    int64_t nel = submesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = submesh->Element(el);
        if(cel)
        {
            PrintElementInfo(cel, out);
        }
    }
}

/// print info on an element group
static void PrintElementGroupInfo(TPZElementGroup *group, std::ostream &out)
{
    out << "+++++ Group index = " << group->Index() << std::endl;
    TPZVec<TPZCompEl *> elvec = group->GetElGroup();
    int numel = elvec.size();
    for (int el=0; el<numel; el++) {
        PrintElementInfo(elvec[el],out);
    }
    out << "End Group Index +++++\n";
}

static void PrintElementInfo(TPZCompEl *cel, std::ostream &out)
{
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
    TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(cel);
    TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
    if(elgr)
    {
        PrintElementGroupInfo(elgr, out);
        return;
    }
    if(submesh)
    {
        PrintSubDomainInfo(submesh, out);
        return;
    }
    if(condense)
    {
        PrintCondensedElementInfo(condense, out);
        return;
    }
    TPZGeoEl *gel = cel->Reference();
    out << "Compel index " << cel->Index() << " ";
    if(gel) PrintGeo(gel,out);
    out << std::endl;
    int nc = cel->NConnects();
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = cel->Connect(ic);
        if(c.NElConnected())
        {
            out << "Cindex " << cel->ConnectIndex(ic) << " ";
            PrintConnect(c,out); out << std::endl;
        }
    }
}

/// print info about the element
void TPZMHMeshControl::PrintMeshInfo(TPZCompMesh *cmesh, std::ostream &out)
{
    out << "mesh name " << cmesh->Name() << std::endl;
    int64_t ncon = cmesh->NConnects();
    for (int64_t ic=0; ic<ncon; ic++) {
        TPZConnect &c = cmesh->ConnectVec()[ic];
        if(c.NElConnected())
        {
            out << "con index " << ic << " ";
            PrintConnect(c, out);
            out << std::endl;
        }
    }
    int64_t nel = cmesh->NElements();
    for(int64_t el = 0; el < nel; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        if(cel) PrintElementInfo(cel, out);
    }
}

/// get the element indices for which a multiphysics element needs to be created
void TPZMHMeshControl::GetGeometricElementPartition(TPZVec<int64_t> &gelindices)
{
    TPZManVector<TPZAutoPointer<TPZCompMesh> > meshvec;
    GetMeshVec(meshvec);
    int64_t nelgeo = fGMesh->NElements();
    TPZVec<int64_t> shouldcreate(nelgeo,0);
    int nmeshes = meshvec.size();
    for (int imesh = 0; imesh < nmeshes; imesh++) {
        TPZCompMesh *cmesh = meshvec[imesh].operator->();
        int64_t nel = cmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            shouldcreate[cel->ReferenceIndex()] = 1;
        }
    }
    int64_t numgeoel = 0;
    for (int64_t el = 0; el<nelgeo; el++) {
        if(shouldcreate[el] == 1)
        {
            numgeoel++;
            TPZGeoEl *gel = fGMesh->Element(el);
            TPZGeoEl *father = gel->Father();
            // make shure only the leaf element will create a multiphysics element
            while (father) {
                shouldcreate[father->Index()] = 0;
                father = father->Father();
            }
        }
    }
    gelindices.resize(numgeoel);
    int64_t counter = 0;
    for (int64_t el = 0; el<nelgeo; el++) {
        if(shouldcreate[el] == 1)
        {
            gelindices[counter] = el;
            counter++;
        }
    }
    gelindices.resize(counter);
}


/// build the multi physics mesh (not at the finest geometric mesh level
void TPZMHMeshControl::BuildMultiPhysicsMesh()
{
    if (fCMesh->NElements() != 0) {
        DebugStop();
    }
    fCMesh->SetAllCreateFunctionsMultiphysicElem();
    TPZMultiphysicsCompMesh *mphysics = fCMesh.operator->();
    int vecsize = 4;
    TPZManVector<TPZCompMesh *> meshvec(vecsize);
    GetMeshVec(meshvec);
    TPZStack<int64_t> gelindexes;
    GetGeometricElementPartition(gelindexes);
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Geometric indices for which we will create multiphysics elements" << std::endl;
        sout << gelindexes;
//        std::cout << sout.str() << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fGMesh->ResetReference();
    mphysics->BuildMultiphysicsSpace(meshvec,gelindexes);
    mphysics->LoadReferences();
    // associate the connects with subdomains
    for(auto index : gelindexes)
    {
        TPZGeoEl *gel = fGMesh->Element(index);
        if(!gel) DebugStop();
        TPZCompEl *cel = gel->Reference();
        if(!cel || cel->Mesh() != mphysics) DebugStop();
        int MHMDomain = fGeoToMHMDomain[index];
        if(MHMDomain >= 0)
        {
            SetSubdomain(cel, MHMDomain);
        }
    }

}

/// Set the Lagrange levels of the connects of the different meshes
void TPZMHMeshControl::SetLagrangeMultiplierLevels()
{
    // NO HYBRIDIZED SKELETON
    //    internal pressure 0
    //    MHM distributed flux 1
    //    delayed pressures (1, 2 or 3 connects) 2
    //    mhm flux 3
    //    MHM average pressure 4
    
    // HYBRIDIZED SKELETON
    //    internal pressure 0
    //    distributed flux 1
    //    delayed pressures (1, 2 or 3 connects) 2
    //    mhm flux 3
    //    average pressure 4
    //    skeleton pressure  5
    
    int dim = fGMesh->Dimension();
    {
        int npress_delay = 0;
        switch(fProblemType)
        {
            case EScalar:
                npress_delay = 1;
                break;
            case EElasticity2D:
                npress_delay = 2;
                break;
            case EElasticity3D:
                npress_delay = 3;
                break;
            default:
                DebugStop();
        }
        std::map<int,int> count_delay;
        int64_t nel = fPressureFineMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fPressureFineMesh->Element(el);
            int subdomain = WhichSubdomain(cel);
            TPZGeoEl *gel = cel->Reference();
            int matid = gel->MaterialId();
            if(fMaterialIds.find(matid) != fMaterialIds.end())
            {
                int ncon = cel->NConnects();
                for (int ic=0; ic<ncon; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    if(count_delay[subdomain] < npress_delay){
                        c.SetLagrangeMultiplier(2);
                        count_delay[subdomain]++;
                    } else
                    {
                        c.SetLagrangeMultiplier(0);
                    }
                }
            } else if(matid == fPressureSkeletonMatId)
            {
                int ncon = cel->NConnects();
                for (int ic=0; ic<ncon; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    c.SetLagrangeMultiplier(5);
                }
            } else if (matid == fWrapMatId)
            {
                
            } else {
                DebugStop();
            }
        }
    }
    // set the lagrange levels of the flux mesh
    {
        int64_t nel = fFluxMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fFluxMesh->Element(el);
            int subdomain = WhichSubdomain(cel);
            TPZGeoEl *gel = cel->Reference();
            int matid = gel->MaterialId();
            if(fMaterialBCIds.find(matid) != fMaterialBCIds.end() || matid == fSkeletonMatId)
            {
                int ncon = cel->NConnects();
                for (int ic=0; ic<ncon; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    c.SetLagrangeMultiplier(3);
                }
            }
            else
            {
                DebugStop();
            }
        }
    }
    if(fLagrangeAveragePressure)
    {
        int64_t nel = fCMeshConstantPressure->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fCMeshConstantPressure->Element(el);
            int ncon = cel->NConnects();
            for (int ic=0; ic<ncon; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.SetLagrangeMultiplier(4);
            }
            cel = fCMeshLagrange->Element(el);
            for (int ic=0; ic<ncon; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.SetLagrangeMultiplier(1);
            }
        }
    }
}

