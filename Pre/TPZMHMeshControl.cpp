//
//  TPZMHMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/3/14.
//
//

#include "TPZMHMeshControl.h"

#include "pzl2projection.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZCompElLagrange.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzsubcmesh.h"
#include "TPZRefPatternTools.h"
#include "pzlog.h"

#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mhmeshcontrol"));
#endif

// toto

TPZMHMeshControl::TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices) : fGMesh(gmesh),
fSkeletonMatId(0), fSecondSkeletonMatId(0), fPressureSkeletonMatId(0), fLagrangeMatIdLeft(50), fLagrangeMatIdRight(51), fLagrangeAveragePressure(false), fHybridize(false)
{
    for (std::set<long>::iterator it=coarseindices.begin(); it != coarseindices.end(); it++) {
        fCoarseIndices[*it]=-1;
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coarse element indexes ";
        for (std::map<long,long>::iterator it=fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fpOrderInternal = 2;
    fpOrderSkeleton = 1;
    fCMesh = new TPZCompMesh(fGMesh);
    fPressureFineMesh = fCMesh;
    fCMesh->SetDimModel(fGMesh->Dimension());
}

TPZMHMeshControl::TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<long> &coarseindices) : fGMesh(gmesh),
fSkeletonMatId(0), fSecondSkeletonMatId(0), fPressureSkeletonMatId(0), fLagrangeMatIdLeft(50), fLagrangeMatIdRight(51), fCoarseIndices(), fLagrangeAveragePressure(false), fHybridize(false)
{
    
    fpOrderInternal = 2;
    fpOrderSkeleton = 1;

    long nc = coarseindices.size();
    for (long c=0; c<nc; c++) {
        fCoarseIndices[coarseindices[c] ] = -1;
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coarse element indexes ";
        for (std::map<long,long>::iterator it=fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fCMesh = new TPZCompMesh(fGMesh);
    fPressureFineMesh = fCMesh;
    fCMesh->SetDimModel(fGMesh->Dimension());
}

TPZMHMeshControl::TPZMHMeshControl(const TPZMHMeshControl &copy){
    
    this->operator=(copy);
}

TPZMHMeshControl &TPZMHMeshControl::operator=(const TPZMHMeshControl &cp){
    
    fGMesh = new TPZGeoMesh(*cp.fGMesh.operator->());
    fSkeletonMatId = cp.fSkeletonMatId;
    fLagrangeMatIdLeft = cp.fLagrangeMatIdLeft;
    fLagrangeMatIdRight = cp.fLagrangeMatIdRight;
    fCoarseIndices = cp.fCoarseIndices;
    fLagrangeAveragePressure = cp.fLagrangeAveragePressure;
    fCMesh = cp.fCMesh;
    fCMesh->SetReference(fGMesh);
    fPressureFineMesh = cp.fPressureFineMesh;
    fPressureFineMesh = cp.fPressureFineMesh;
    fpOrderSkeleton = cp.fpOrderSkeleton;
    fpOrderInternal = cp.fpOrderInternal;
    fHybridize = cp.fHybridize;
    return *this;
}

/// will create 1D elements on the interfaces between the coarse element indices
void TPZMHMeshControl::CreateSkeletonElements(int matid)
{
    if (fInterfaces.size()) {
        DebugStop();
    }
    
    if(matid < 0) DebugStop();
    
    fSkeletonMatId = matid;
    
    TPZCompMesh *cmesh = CriaMalhaTemporaria();
    
    long nel = fGMesh->NElements();
    int dimension = fGMesh->Dimension();
    int ninterf;
    
    for(long iel = 0; iel<nel; iel++)
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
            int interfacematid = matid;
            TPZCompEl *left = intface->LeftElement();
            TPZCompEl *right = intface->RightElement();
            long leftind = left->Reference()->Index();
            long rightind = right->Reference()->Index();
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
    
    
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        std::map<long, std::pair<long,long> >::iterator it = fInterfaces.begin();
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
    std::map<long, std::pair<long,long> >::iterator it;
    for (int divide=0; divide<ndivide; divide++)
    {
        std::map<long, std::pair<long,long> > mapdivided;
        for (it=fInterfaces.begin(); it!=fInterfaces.end(); it++) {
            long elindex = it->first;
//            if (elindex == it->second.second) {
//                mapdivided[elindex] = it->second;
//                continue;
//            }
            TPZGeoEl *gel = fGMesh->Element(elindex);
            TPZAutoPointer<TPZRefPattern> refpat = TPZRefPatternTools::PerfectMatchRefPattern(gel);
            gel->SetRefPattern(refpat);
            TPZManVector<TPZGeoEl *,10> subels;
            gel->Divide(subels);
            long nsub = subels.size();
            for (int is=0; is<nsub; is++) {
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
}


TPZCompMesh* TPZMHMeshControl::CriaMalhaTemporaria()
{
    fGMesh->ResetReference();
	/// criar materiais
	int dim = fGMesh->Dimension();
    std::set<int> matids, bcids;
    TPZGeoMesh &gmesh = fGMesh;
    long nel = gmesh.NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fCoarseIndices.find(el) == fCoarseIndices.end()) {
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
    
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fCoarseIndices.find(el) == fCoarseIndices.end()) {
            continue;
        }
        long index;
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
    int nstate = 1;
    int dim = fGMesh->Dimension();
    fCMesh->SetDimModel(dim);
    TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(fLagrangeMatIdLeft,dim,nstate);
    TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(fLagrangeMatIdRight,dim,nstate);
    matleft->SetMultiplier(1.);
    matright->SetMultiplier(-1.);
    fCMesh->InsertMaterialObject(matleft);
    fCMesh->InsertMaterialObject(matright);
    CreateInternalElements();
    //    AddBoundaryElements();
    CreateSkeleton();
    CreateInterfaceElements();
//    AddBoundaryInterfaceElements();
    fCMesh->ExpandSolution();
    fCMesh->CleanUpUnconnectedNodes();
    if (fHybridize) {
        Hybridize();
    }
    if (fLagrangeAveragePressure) {
        this->CreateLagrangeMultiplierMesh();
        this->TransferToMultiphysics();
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "*********** BEFORE SUBSTRUCTURING *************\n";
        fCMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    if(usersubstructure==true){
        this->SubStructure();
    }
}

/// will create the internal elements, one coarse element at a time
void TPZMHMeshControl::CreateInternalElements()
{
    TPZCompEl::SetgOrder(fpOrderInternal);
    
    fCMesh->ApproxSpace().SetCreateLagrange(false);
    fCMesh->SetAllCreateFunctionsContinuous();
    
    //Criar elementos computacionais malha MHM
    
    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    fConnectToSubDomainIdentifier.Expand(10000);
    
    for (std::map<long,long>::iterator it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++)
    {
        std::set<long> elset;
        bool LagrangeCreated = false;
        long iel = it->first;
        elset.insert(iel);
//        std::map<long, std::pair<long,long> >::iterator it2;
//        for (it2 = fInterfaces.begin(); it2 != fInterfaces.end(); it2++) {
//            if (it2->first == it2->second.second && it2->second.first == *it) {
//                elset.insert(it2->first);
//            }
//        }

        while (elset.size()) {
            std::set<long>::iterator itel = elset.begin();
            long elfirst = *itel;
            elset.erase(itel);
            gel = fGMesh->ElementVec()[elfirst];
            if(!gel) DebugStop();
            // we work with only the leaf elements
            if (! gel->HasSubElement()) {
                long index;
                // create the flux element
                fCMesh->CreateCompEl(gel, index);
                TPZCompEl *cel = fCMesh->Element(index);
                /// associate the connects with the subdomain
                SetSubdomain(cel, it->first);
                // we need to create a lagrange multiplier element in order to delay decomposition of an equation
                if (!LagrangeCreated)
                {
                    LagrangeCreated = true;
                    TPZCompEl *cel = fCMesh->ElementVec()[index];
                    long cindex = cel->ConnectIndex(0);
                    int nshape(1), nvar(1), order(1);
                    int lagrangelevel = 0;
                    if (this->fLagrangeAveragePressure) {
                        lagrangelevel = 1;
                    }
                    else
                    {
                        lagrangelevel = 3;
                    }
                    fCMesh->ConnectVec()[cindex].SetLagrangeMultiplier(lagrangelevel);
                    
                }
                continue;
            }
            int nsubels = gel->NSubElements();
            for (int is=0; is<nsubels; is++)
            {
                gsubel = gel->SubElement(is);
                elset.insert(gsubel->Index());
            }
        }
        
        fGMesh->ResetReference();
        
    }
}

/// will create the elements on the skeleton
void TPZMHMeshControl::CreateSkeleton()
{
    // comment this line or not to switch the type of skeleton elements
    int meshdim = fCMesh->Dimension();
    fCMesh->SetDimModel(meshdim);
//    fCMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    fCMesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    int order = fpOrderSkeleton;
    if (order < 0) {
        order = 0;
    }
    fCMesh->SetDefaultOrder(order);
    std::map<long, std::pair<long,long> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        long elindex = it->first;
        // skip the boundary elements
//        if (elindex == it->second.second) {
//            it++;
//            continue;
//        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        long index;
        // create a discontinuous element to model the flux
        fCMesh->CreateCompEl(gel, index);
        TPZCompEl *cel = fCMesh->ElementVec()[index];
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
    fCMesh->SetDimModel(meshdim);
}

/// will create the interface elements between the internal elements and the skeleton
void TPZMHMeshControl::CreateInterfaceElements()
{
    fCMesh->LoadReferences();
    int dim = fGMesh->Dimension();
    std::map<long, std::pair<long,long> >::iterator it = fInterfaces.begin();
    /// loop over the skeleton elements
    while (it != fInterfaces.end()) {
        // index of the skeleton element
        long elindex = it->first;
        // left and right indexes in the coarse mesh
        long leftelindex = it->second.first;
        long rightelindex = it->second.second;
        // skip boundary elements
//        if (elindex == rightelindex) {
//            it++;
//            continue;
//        }
        
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
//        if (matid != fSkeletonMatId) {
//            DebugStop();
//        }
        std::set<long> celindices;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZCompElSide celskeleton = gelside.Reference();
        // the skeleton element must exist
        if(!celskeleton.Element()) DebugStop();
        TPZGeoElSide neighbour = gelside.Neighbour();
        // the geometric neighbour of the skeleton elements are the macro elements
        while (neighbour != gelside) {
            TPZStack<TPZCompElSide> celstack;
            bool onlyinterpolated = true;
            gelside.HigherLevelCompElementList2(celstack, onlyinterpolated, 0);
            // I dont remember why we need to call this procedure twice??
            long nst = celstack.size();
            for (long ist = 0; ist<nst; ist++) {
                TPZGeoElSide sideloc = celstack[ist].Reference();
                sideloc.HigherLevelCompElementList2(celstack, onlyinterpolated, 0);
            }
#ifdef PZDEBUG
            if (celstack.size() != nst) {
                std::cout << __PRETTY_FUNCTION__ << " check it out to improve documentation\n";
            }
#endif
            // for the case where the macro elements havent been refined
            gelside.EqualLevelCompElementList(celstack, onlyinterpolated, 0);
            // check if there is a neighbour larger than the skeleton element???
            TPZCompElSide tmp = gelside.LowerLevelCompElementList2(0);
            if (tmp) {
                celstack.Push(tmp);
            }
            long nstack = celstack.NElements();
            for (long i=0; i<nstack; i++) {
                TPZCompElSide celside = celstack[i];
                TPZGeoElSide gelneigh = celside.Reference();
                // the sons of an element include entities of lower dimension
                // we are only interested in entities of the size of boundary conditions
                if (gelneigh.Dimension() != dim-1) {
                    continue;
                }
                long celindex = celside.Element()->Index();
                // avoid creating two interfaces for a single element
                if (celindices.find(celindex) != celindices.end()) {
                    continue;
                }
                long gelneighindex = gelneigh.Element()->Index();
                // determine if the element is son of the left or right side of the MHM flux element
                bool issiblingleft = IsSibling(gelneighindex, leftelindex);
                bool issiblingright = IsSibling(gelneighindex, rightelindex);
#ifdef PZDEBUG
                if (!issiblingleft && !issiblingright) {
                    DebugStop();
                }
#endif
                int matid = 0;
                if (issiblingleft && leftelindex < rightelindex)
                {
                    matid = fLagrangeMatIdLeft;
                }
                else
                {
                    matid = fLagrangeMatIdRight;
                }
                // this means the element is a boundary. Retain its material id
                if (gelneigh.Element()->Dimension() < dim || matid == 0) {
                    DebugStop();
//                    matid = gelneigh.Element()->MaterialId();
                }
                celindices.insert(celindex);
                if (!tmp || i < nstack-1)
                {
                    // create an interface between the finer element and the MHM flux
                    long index;
                    TPZGeoEl *gelnew = gelneigh.Element()->CreateBCGeoEl(gelneigh.Side(), matid);
                    new TPZInterfaceElement(fCMesh, gelnew , index, celside, celskeleton);
    #ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "New interface left " << gelneigh.Element()->Index() << " right " << gel->Index();
                        LOGPZ_DEBUG(logger, sout.str())
                    }
    #endif
                }
                else
                {
                    // create an interface between the larger element and the MHM flux
                    long index;
                    TPZGeoEl *gelnew = gelside.Element()->CreateBCGeoEl(gelside.Side(), matid);
                    new TPZInterfaceElement(fCMesh, gelnew , index, celside, celskeleton);
#ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "New interface left " << gelneigh.Element()->Index() << " right " << gel->Index();
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
                
            }
            neighbour = neighbour.Neighbour();
        }
        it++;
    }
    // loop over the boundary elements and create the interfaces

}

/// verify if the element is a sibling of
bool TPZMHMeshControl::IsSibling(long son, long father)
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
    int nelgroups = this->fCoarseIndices.size();
    out << __PRETTY_FUNCTION__ << " Number of coarse elements " << nelgroups << std::endl;
    
    for (std::map<long,long>::iterator it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
        PrintSubdomain(it->first, out);
    }
    PrintBoundaryInfo(out);
}

/// print the diagnostics for a subdomain
void TPZMHMeshControl::PrintSubdomain(long elindex, std::ostream &out)
{
    long nel = fCMesh->NElements();
    std::set<long> celindices;
    std::multimap<long,long> interfaces;
    TPZStack<TPZCompElSide> boundaries;
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        long gelindex = gel->Index();
        if(IsSibling(gelindex, elindex))
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
    for (std::set<long>::iterator i=celindices.begin(); i != celindices.end(); i++) {
        out << *i << " ";
    }
    out << std::endl;
    out << "Number of boundary interfaces " << interfaces.size() << std::endl;
    for (std::multimap<long, long>::iterator i = interfaces.begin(); i != interfaces.end(); i++) {
        long el = i->second;
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!face) {
            DebugStop();
        }
        out << "Comp Element index " << face->LeftElement()->Index() << " side " << face->LeftElementSide().Side() << " skeleton element " << face->RightElement()->Index() << " interface index "
            << face->Index() << " face material id " << face->Reference()->MaterialId() << std::endl;
    }
    out << "Number of elements on the sides of elseed " << boundaries.size() << std::endl;
    for (long i=0; i<boundaries.size(); i++) {
        out << "Comp Element index " << boundaries[i].Element()->Index() << " side " << boundaries[i].Side() << std::endl;
    }
    out << "Geom Element indices ";
    for (std::set<long>::iterator i=celindices.begin(); i != celindices.end(); i++) {
        out << fCMesh->ElementVec()[*i]->Reference()->Index() << " ";
    }
    out << std::endl;
    out << "Number of boundary interfaces " << interfaces.size() << std::endl;
    for (std::multimap<long, long>::iterator i = interfaces.begin(); i != interfaces.end(); i++) {
        long el = i->second;
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!face) {
            DebugStop();
        }
        out << "Geom Element index " << face->LeftElement()->Reference()->Index() << " side " << face->LeftElementSide().Side() << " skeleton element " << face->RightElement()->Reference()->Index() << " interface index "
        << face->Index() << " face material id " << face->Reference()->MaterialId() << std::endl;
    }
    out << "Number of elements on the sides of elseed " << boundaries.size() << std::endl;
    for (long i=0; i<boundaries.size(); i++) {
        out << "Geom Element index " << boundaries[i].Element()->Reference()->Index() << " side " << boundaries[i].Side() << std::endl;
    }
}

/// put the element side which face the boundary on the stack
void TPZMHMeshControl::AddElementBoundaries(long elseed, long compelindex, TPZStack<TPZCompElSide> &result)
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
        // if the elseed is not a father along the boundary, the element is not on the boundary
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
    long nel = fGMesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->ElementVec()[el];
        if (!gel || gel->Dimension() == dim || gel->HasSubElement()) {
            continue;
        }
        int matid = gel->MaterialId();
        if (notincluded.find(matid) != notincluded.end()) {
            continue;
        }
        long index;
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
    long nel = fCMesh->NElements();
    for (long el=0; el<nel; el++) {
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
        long nst = celstack.size();
        for (long i=0; i<nst; i++) {
            TPZCompElSide cs = celstack[i];
            TPZGeoElSide gs = cs.Reference();
            if (gs == gelside) {
                continue;
            }
            long index;
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

    long nelem = fCMesh->NElements();
    for (long el=0; el<nelem; el++) {
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
        for (long i=0; i< celstack.NElements(); i++) {
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
    for (long el=0; el<nelem; el++) {
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
        for (long i=0; i< celstack.NElements(); i++) {
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
    fCMeshLagrange = new TPZCompMesh(fGMesh);
    int dim = fGMesh->Dimension();
    fCMeshLagrange->SetDimModel(dim);
    fCMeshLagrange->SetAllCreateFunctionsDiscontinuous();
    fCMeshLagrange->SetDefaultOrder(0);
    fGMesh->ResetReference();
    long connectcounter = fCMesh->NConnects();
	/// criar materiais
    std::set<int> matids;
    TPZGeoMesh &gmesh = fGMesh;
    long nel = gmesh.NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fCoarseIndices.find(el) == fCoarseIndices.end()) {
            continue;
        }
        int materialid = gel->MaterialId();
        matids.insert(materialid);
    }
    
    TPZManVector<STATE,1> sol(1,0.);
    int nstate = 1;
    std::set<int>::iterator it = matids.begin();
    TPZMaterial *meshmat = 0;
    while (it != matids.end()) {
        TPZL2Projection *material = new TPZL2Projection(*it,dim,nstate,sol);
        fCMeshLagrange->InsertMaterialObject(material);
        if (!meshmat) {
            meshmat = material;
        }
        it++;
        
    }
    if (!meshmat) {
        DebugStop();
    }
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim || fCoarseIndices.find(el) == fCoarseIndices.end()) {
            continue;
        }
        long index;
        TPZCompEl *disc = new TPZCompElDisc(fCMeshLagrange,gel,index);
        long cindex = disc->ConnectIndex(0);
        SetSubdomain(cindex, el,connectcounter);
        
//        fCMeshConstantStates->CreateCompEl(gel, index);
    }
    fCMeshLagrange->ExpandSolution();
    fGMesh->ResetReference();
    
    long connectcounter2 = fCMeshLagrange->NConnects();
    for (long i=connectcounter+connectcounter2; i<connectcounter+2*connectcounter2; i++) {
        SetSubdomain(i, fConnectToSubDomainIdentifier[i-connectcounter2]);
    }
    fCMeshConstantPressure = new TPZCompMesh(fCMeshLagrange);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        fCMeshLagrange->Print(sout);
        fCMeshConstantPressure->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fGMesh->ResetReference();
}

/// transform the computational mesh into a multiphysics mesh
void TPZMHMeshControl::TransferToMultiphysics()
{
    fGMesh->ResetReference();
    this->fCMesh = new TPZCompMesh(fGMesh);
    this->fCMesh->SetDimModel(fGMesh->Dimension());
    fCMesh->SetAllCreateFunctionsMultiphysicElem();
    
    // copy the material objects
    std::map<int,TPZMaterial *>::iterator it = fPressureFineMesh->MaterialVec().begin();
    while (it != fPressureFineMesh->MaterialVec().end()) {
        it->second->Clone(fCMesh->MaterialVec());
        it++;
    }
    
    long nel = fPressureFineMesh->NElements();
    for (long el=0; el<nel; el++) {
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
            long index;
            new TPZCompElLagrange(fCMesh,*cel,index);
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        long index;
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
    for (long el=0; el<nel; el++) {
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
    for (long el=0; el<nel; el++) {
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
    for (long el=0; el<nel; el++) {
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
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!intface) {
            continue;
        }
        //        TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, long &index, TPZCompElSide left, TPZCompElSide right);
        TPZCompElSide pressleft = intface->LeftElementSide();
        TPZCompElSide pressright = intface->RightElementSide();
        TPZGeoElSide gleft = pressleft.Reference();
        TPZGeoElSide gright = pressright.Reference();
        TPZCompElSide multleft = gleft.Reference();
        TPZCompElSide multright = gright.Reference();
        long index;
        new TPZMultiphysicsInterfaceElement(fCMesh,intface->Reference(),index,multleft,multright);
        
    }

    
    nel = fCMeshConstantPressure->NElements();
    long npressconnect = fPressureFineMesh->NConnects();
    long nlagrangeconnect = fCMeshLagrange->NConnects();
    // nel numero de dominios MHM, tem um connect associado a cada um e os mesmos estao no final
    for (long el=0; el<nel; el++)
    {
        TPZCompEl *cel = this->fCMeshConstantPressure->Element(el);
        long pressureconnect = cel->ConnectIndex(0);
        long cindex = npressconnect+nlagrangeconnect+pressureconnect;
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
void TPZMHMeshControl::SubStructure2()
{
    // for each connect index, the submesh index
    std::map<long, long > connectdest;
    // for each coarse geometric index, a subcompmesh
    std::map<long, TPZSubCompMesh *> submeshes;
    std::map<long,long>::iterator it = fCoarseIndices.begin();
    
    // create the submeshes
    while (it != fCoarseIndices.end()) {
        long index;
        TPZSubCompMesh *submesh = new TPZSubCompMesh(fCMesh,index);
        submeshes[it->first] = submesh;
        it++;
    }
    for (std::map<long, TPZSubCompMesh *>::iterator it = submeshes.begin(); it != submeshes.end(); it++) {
        fCoarseIndices[it->first] = it->second->Index();
    }
    fGMesh->ResetReference();
    fCMesh->LoadReferences();
    
    // put the internal elements in the submeshes
    // build the list of connects which are internal
    it = fCoarseIndices.begin();
    while (it != fCoarseIndices.end()) {
        long index = it->first;
        
        // put all the sons of gel in the submesh
        TPZGeoEl *gel = fGMesh->ElementVec()[index];
        if (!gel) {
            DebugStop();
        }
        TPZStack<TPZCompElSide> celstack;
        int ns = gel->NSides();
        TPZGeoElSide gelside(gel,ns-1);
        gelside.ConnectedCompElementList(celstack, 0, 0);
        if (gel->Reference()) {
            celstack.Push(gelside.Reference());
        }
        // find boundary elements which neighbour the coarse element
        std::map<long,std::pair<long, long> >::iterator it2;
        for (it2 = fInterfaces.begin(); it2 != fInterfaces.end(); it2++) {
            long elindex = it2->first;
            // if we have a boundary element
            // put the element or its siblings in the celstack
            // all elements in celstack will be transferred to the subdomain
            if (elindex == it2->second.second && it2->second.first == index) {
                TPZGeoEl *bound = fGMesh->Element(elindex);
                TPZGeoElSide boundside(bound,bound->NSides()-1);
                // if the skeleton has the same resolution as the internal mesh
                if (bound->Reference()) {
                    celstack.Push(boundside.Reference());
                }
                // put all subelements on the stack - THIS ONLY WORKS FOR DISCONTINUOUS ELEMENTS
                else
                {
                    boundside.ConnectedCompElementList(celstack, 0, 0);
                }
            }
        }
        int ncel = celstack.size();
        for (int icel=0; icel<ncel; icel++) {
            TPZCompElSide celside = celstack[icel];
            TPZCompEl *cel = celside.Element();
            if (cel->Mesh() != fCMesh.operator->()) {
                continue;
            }
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
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                long connectindex = cel->ConnectIndex(ic);
//                if (connectdest.find(connectindex) != connectdest.end()) {
//                    DebugStop();
//                }
                long submeshindex = submeshes[it->first]->Index();
                connectdest[connectindex] = submeshindex;
            }
            submeshes[it->first]->TransferElement(fCMesh.operator->(), cel->Index());
        }
        it++;
    }
    // transfer all other elements (interface and lagrange)
    // if the element contains a connect which belongs to the subdomain make it part of the subdomain
    long nel = fCMesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (submesh) {
            continue;
        }
        int subdomain = -1;
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            long connectindex = cel->ConnectIndex(ic);
            // connectdest : for each connect index the subdomain destination
            if (connectdest.find(connectindex) != connectdest.end()) {
                long domain = connectdest[connectindex];
                if (subdomain != -1 && subdomain != domain) {
                    DebugStop();
                }
                else
                {
                    subdomain = domain;
                }
            }
        }
        if (subdomain != -1)
        {
            TPZCompEl *subcel = fCMesh->ElementVec()[subdomain];
            submesh = dynamic_cast<TPZSubCompMesh *>(subcel);
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
            submesh->TransferElement(fCMesh.operator->(), el);
        }
    }
    fCMesh->ComputeNodElCon();
    
    std::map<long, TPZSubCompMesh *>::iterator itsub = submeshes.begin();
    while (itsub != submeshes.end()) {
        TPZSubCompMesh *submesh = itsub->second;
        int nc = submesh->NConnects();
        std::set<long> internals;
        // put all connects with one element connection internal in the submesh
        for (int ic=0; ic<nc; ic++) {
            long connectindex = submesh->ConnectIndex(ic);
            TPZConnect &c = submesh->Connect(ic);
            int lagrange = c.LagrangeMultiplier();
            if (c.NElConnected() >1) {
                continue;
            }
            if ((this->fLagrangeAveragePressure && lagrange < 3) || lagrange < 2) {
                long internal = submesh->InternalIndex(connectindex);
                internals.insert(internal);
            }
        }
        for (std::set<long>::iterator it = internals.begin(); it != internals.end(); it++) {
            submesh->MakeInternal(*it);
        }
        submesh->ExpandSolution();
        itsub++;
    }
    fCMesh->CleanUpUnconnectedNodes();
    itsub = submeshes.begin();
    while (itsub != submeshes.end()) {
        TPZSubCompMesh *submesh = itsub->second;
        int numthreads = 0;
        int preconditioned = 0;
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            submesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZGuiInterface *gui = 0;
        submesh->SetAnalysisSkyline(numthreads, preconditioned, gui);
        itsub++;
    }
    
    fCMesh->SaddlePermute();
}

/// substructure the mesh
void TPZMHMeshControl::SubStructure()
{
    // for each connect index, the submesh index
    std::map<long, long > connectdest;
    // for each coarse geometric index, a subcompmesh
    std::map<long, TPZSubCompMesh *> submeshes;
    std::map<long,long>::iterator it = fCoarseIndices.begin();
    
    // create the submeshes
    while (it != fCoarseIndices.end()) {
        long index;
        TPZSubCompMesh *submesh = new TPZSubCompMesh(fCMesh,index);
        submeshes[it->first] = submesh;
        it++;
    }
    for (std::map<long, TPZSubCompMesh *>::iterator it = submeshes.begin(); it != submeshes.end(); it++) {
        fCoarseIndices[it->first] = it->second->Index();
    }

    fGMesh->ResetReference();
    fCMesh->LoadReferences();

    long nel = fCMesh->NElements();
    for (long el=0; el<nel; el++)
    {
        TPZCompEl *cel = fCMesh->Element(el);
        if(!cel) continue;
        if (dynamic_cast<TPZSubCompMesh *>(cel)) {
            continue;
        }
        long domain = WhichSubdomain(cel);

        if (domain == -1) {
            continue;
        }
        if (submeshes.find(domain) == submeshes.end()) {
            DebugStop();
        }
        submeshes[domain]->TransferElement(fCMesh.operator->(), cel->Index());
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
    fCMesh->ComputeNodElCon();
    
        
    std::map<long, TPZSubCompMesh *>::iterator itsub = submeshes.begin();
    while (itsub != submeshes.end()) {
        TPZSubCompMesh *submesh = itsub->second;
        int nc = submesh->NConnects();
        std::set<long> internals;
        // put all connects with one element connection internal in the submesh
        for (int ic=0; ic<nc; ic++) {
            long connectindex = submesh->ConnectIndex(ic);
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
            else if ((this->fLagrangeAveragePressure && lagrange < 3) || lagrange < 3) {
                makeinternal = true;
            }
            if (makeinternal)
            {
                long internal = submesh->InternalIndex(connectindex);
                internals.insert(internal);
            }
        }
        for (std::set<long>::iterator it = internals.begin(); it != internals.end(); it++) {
            submesh->MakeInternal(*it);
        }
        submesh->ExpandSolution();
        itsub++;
    }
    fCMesh->CleanUpUnconnectedNodes();
    itsub = submeshes.begin();
    while (itsub != submeshes.end()) {
        TPZSubCompMesh *submesh = itsub->second;
        int numthreads = 0;
        int preconditioned = 0;
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
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
        for (std::map<long,long>::iterator it= fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
            out << it->first << " " << it->second << " ";
        }
        out << std::endl;
    }
    /// indices of the skeleton elements and their left/right elements of the skeleton mesh
    out << "Skeleton element indices with associated left and right coarse element indices\n";
    {
        std::map<long, std::pair<long,long> >::iterator it;
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
    for (long i=0; i<fConnectToSubDomainIdentifier.size(); i++) {
        if (i && !(i%20)) {
            out << std::endl;
        }
        out << "(" << i << "->" << fConnectToSubDomainIdentifier[i] << ") ";
    }
    out << std::endl;

}

void TPZMHMeshControl::Hybridize()
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
    // expand the vector
    fConnectToSubDomainIdentifier.Resize(fCMesh->NConnects()+10000, -1);
    fCMesh->SetDefaultOrder(order);
    long nelem = fGMesh->NElements();
    fCMesh->LoadReferences();
    // keep the pointers of the skeleton elements we will work
    TPZManVector<TPZCompEl *> skeletoncomp(nelem,0);
    std::map<long, std::pair<long,long> >::iterator it;
    // loop over the skeleton elements
    for (it=fInterfaces.begin(); it != fInterfaces.end(); it++) {
        long elindex = it->first;
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
        long elindex = it->first;
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
        long index1,index2;
        TPZManVector<TPZGeoEl *,4> GelVec(3);
        GelVec[0] = gel;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        // create geometric elements along the skeleton element
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
        long elindex = it->first;
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

        // find all interface elements (and others) connected to the pressure element
        TPZStack<TPZCompElSide> celstack;
        neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            if (neighbour.Reference()) {
                celstack.Push(neighbour.Reference());
            }
            else
            {
                neighbour.HigherLevelCompElementList2(celstack, 0, 0);
            }
            neighbour = neighbour.Neighbour();
        }
        long nelst = celstack.size();
        // keep track of the elements that have been adjusted
        std::set<TPZCompEl *> checked;
        for (long elst = 0; elst < nelst; elst++) {
            TPZCompEl *celst = celstack[elst].Element();
            TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(celst);
            if (!intface) {
                continue;
            }
            // adjust each interface element only once
            if (checked.find(intface) != checked.end()) {
                continue;
            }
            checked.insert(intface);
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
            
            long leftcindex = left.Element()->ConnectIndex(0);
            int subdomain = fConnectToSubDomainIdentifier[leftcindex];
            if (subdomain == -1) {
                DebugStop();
            }
            SetSubdomain(right.Element(), subdomain);
        }
        
        // create the interfaces between the flux elements and the newly created pressure element
        TPZCompElSide leftflux(celskeleton,gel->NSides()-1);
        TPZCompElSide rightflux(celsecondskeleton,gel->NSides()-1);
        TPZCompElSide pressure(celpressure,gel->NSides()-1);
        TPZGeoElBC gbc3(gelside,fLagrangeMatIdRight);
        TPZGeoElBC gbc4(gelside,fLagrangeMatIdLeft);
        long index3, index4;
        TPZInterfaceElement *intfaceleft = new TPZInterfaceElement(fCMesh,gbc3.CreatedElement(),index3);
        TPZInterfaceElement *intfaceright = new TPZInterfaceElement(fCMesh,gbc4.CreatedElement(),index4);
        intfaceleft->SetLeftRightElements(leftflux, pressure);
        intfaceright->SetLeftRightElements(rightflux, pressure);
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
    fConnectToSubDomainIdentifier.Resize(fCMesh->NConnects(), -1);

}

/// associates the connects of an element with a subdomain
void TPZMHMeshControl::SetSubdomain(TPZCompEl *cel, long subdomain, long offset)
{
    int ncon = cel->NConnects();
    for (int ic=0; ic<ncon; ic++) {
        long cindex = cel->ConnectIndex(ic);
        SetSubdomain(cindex, subdomain, offset);
    }
}

/// associates the connects index with a subdomain
void TPZMHMeshControl::SetSubdomain(long cindex, long subdomain, long offset)
{
    if (cindex+offset >= fConnectToSubDomainIdentifier.size()) {
        fConnectToSubDomainIdentifier.Resize(cindex+offset+1, -1);
    }
    fConnectToSubDomainIdentifier[cindex+offset] = subdomain;

}

/// returns to which subdomain a given element belongs
// this method calls debugstop if the element belongs to two subdomains
long TPZMHMeshControl::WhichSubdomain(TPZCompEl *cel, long offset)
{
    int ncon = cel->NConnects();
    std::set<long> domains;
    for (int ic=0; ic<ncon; ic++)
    {
        long cindex = cel->ConnectIndex(ic);
        if (fConnectToSubDomainIdentifier[cindex] != -1) {
            domains.insert(fConnectToSubDomainIdentifier[cindex]);
        }
    }
    if (domains.size() > 1) {
        for (int ic=0; ic<ncon; ic++) {
            long cindex = cel->ConnectIndex(ic);
            std::cout << cindex << "|" << fConnectToSubDomainIdentifier[cindex] << " ";
        }
        std::cout << std::endl;
        DebugStop();
    }
    if (domains.size() ==0) {
        return -1;
    }
    long domain = *domains.begin();
    return domain;
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
    else
    {
        return false;
    }
}
/// identify connected elements to the skeleton elements
// the computational mesh is determined by the element pointed to by the geometric element
void TPZMHMeshControl::ConnectedElements(long skeleton, std::pair<long,long> &leftright, std::map<long, std::list<TPZCompElSide> > &ellist)
{
    TPZGeoEl *gelskeleton = fGMesh->Element(skeleton);
    int meshdim = fGMesh->Dimension();
    if (gelskeleton->Dimension() != meshdim-1) {
        DebugStop();
    }
    TPZGeoElSide gelside(gelskeleton,gelskeleton->NSides()-1);
    TPZStack<TPZGeoElSide> connected;
    TPZGeoElSide neighbour(gelside.Neighbour());
    bool leftfound = false;
    bool rightfound = false;
    std::pair<TPZGeoEl *, TPZGeoEl *> leftrightgeo(fGMesh->Element(leftright.first),fGMesh->Element(leftright.second));
    while (neighbour != gelside) {
        if (neighbour.Element()->Dimension() == meshdim) {
            if (!leftfound && IsAncestor(neighbour.Element(),leftrightgeo.first))
            {
                leftfound = true;
                connected.Push(neighbour);
            }
            if (!rightfound && IsAncestor(neighbour.Element(),leftrightgeo.second))
            {
                rightfound = true;
                connected.Push(neighbour);
            }
        }
        neighbour = neighbour.Neighbour();
    }
    if (skeleton != leftright.second && connected.size() != 2) {
        DebugStop();
    }
    if (skeleton == leftright.second && connected.size() != 1)
    {
        DebugStop();
    }
    // vou ter que pegar qualquer connected filho de left ou right e verificar que o tamanho eh pelo menos 1 ou 2
    for (int i=0; i<connected.size(); i++) {
        long rootel = -1;
        if (IsAncestor(connected[i].Element(),leftrightgeo.first))
        {
            rootel = leftrightgeo.first->Index();
        }
        if (IsAncestor(connected[i].Element(),leftrightgeo.second))
        {
            rootel = leftrightgeo.second->Index();
        }
        if (rootel == -1) {
            DebugStop();
        }
        TPZStack<TPZGeoElSide> tocheck;
        tocheck.Push(connected[i]);
        while (tocheck.size())
        {
            TPZGeoElSide work = tocheck.Pop();
            if (work.Element()->Reference()) {
                ellist[rootel].push_back(work.Reference());
            }
            else
            {
                if (!work.HasSubElement()) {
                    DebugStop();
                }
                TPZStack<TPZGeoElSide> subels;
                work.GetSubElements2(subels);
                for (int i=0; i<subels.size(); i++) {
                    TPZGeoElSide son = subels[i];
                    if (son.Dimension() == meshdim - 1) {
                        tocheck.Push(son);
                    }
                }
            }
        }
    }
}
