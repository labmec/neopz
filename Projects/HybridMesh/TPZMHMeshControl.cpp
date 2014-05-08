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
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mhmeshcontrol"));
#endif

TPZMHMeshControl::TPZMHMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices) : fGMesh(gmesh),
fSkeletonMatId(0), fLagrangeMatIdLeft(50), fLagrangeMatIdRight(51), fCoarseIndices(coarseindices)
{
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coarse element indexes ";
        for (std::set<long>::iterator it=fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    fCMesh = new TPZCompMesh(fGMesh);
    fCMesh->SetDimModel(fGMesh->Dimension());
}

/// will create 1D elements on the interfaces between the coarse element indices
void TPZMHMeshControl::CreateCoarseInterfaces(int matid)
{
    if (fInterfaces.size()) {
        DebugStop();
    }
    
    fSkeletonMatId = matid;
    
    TPZCompMesh *cmesh = CriaMalhaTemporaria();
    
    long nel = fGMesh->NElements();
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
            TPZCompEl *left = intface->LeftElement();
            TPZCompEl *right = intface->RightElement();
            long leftind = left->Reference()->Index();
            long rightind = right->Reference()->Index();
            fInterfaces[iel] = std::make_pair(leftind, rightind);
            gel->SetMaterialId(matid);
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
                    //bcids.insert(neighbour.Element()->MaterialId());
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
    
    return cmesh;
}

/// Create all data structures for the computational mesh
void TPZMHMeshControl::BuildComputationalMesh()
{
    int nstate = 1;
    int dim = fGMesh->Dimension();
    TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(fLagrangeMatIdLeft,dim,nstate);
    TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(fLagrangeMatIdRight,dim,nstate);
    matright->SetMultiplier(-1.);
    fCMesh->InsertMaterialObject(matleft);
    fCMesh->InsertMaterialObject(matright);
    CreateInternalElements();
    AddBoundaryElements();
    CreateSkeleton();
    CreateInterfaceElements();
    AddBoundaryInterfaceElements();
    fCMesh->ExpandSolution();
}

/// will create the internal elements, one coarse element at a time
void TPZMHMeshControl::CreateInternalElements()
{
    TPZCompEl::SetgOrder(fpOrderInternal);
    
    
    fCMesh->SetAllCreateFunctionsContinuous();
    
    //Criar elementos computacionais malha MHM
    
    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    
    for (std::set<long>::iterator it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++)
    {
        std::set<long> elset;
        bool LagrangeCreated = false;
        long iel = *it;
        elset.insert(iel);
        while (elset.size()) {
            std::set<long>::iterator itel = elset.begin();
            long elfirst = *itel;
            elset.erase(itel);
            gel = fGMesh->ElementVec()[elfirst];
            if(!gel) DebugStop();
            if (! gel->HasSubElement()) {
                long index;
                fCMesh->CreateCompEl(gel, index);
                if (!LagrangeCreated) {
                    LagrangeCreated = true;
                    TPZCompEl *cel = fCMesh->ElementVec()[index];
                    long cindex = cel->ConnectIndex(0);
                    int nshape(1), nvar(1), order(1);
                    long newconnect = fCMesh->AllocateNewConnect(nshape,nvar,order);
                    fCMesh->ConnectVec()[newconnect].SetLagrangeMultiplier(2);
                    long index;
                    /*TPZCompElLagrange *lagrange = */new TPZCompElLagrange(fCMesh,cindex,0,newconnect,0,index);
                    
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
//    fCMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    std::map<long, std::pair<long,long> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        long elindex = it->first;
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        long index;
        fCMesh->CreateCompEl(gel, index);
        TPZCompEl *cel = fCMesh->ElementVec()[index];
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            cel->Connect(ic).SetLagrangeMultiplier(1);
        }
        fCMesh->ElementVec()[index]->Reference()->ResetReference();
        it++;
    }
}

/// will create the interface elements between the internal elements and the skeleton
void TPZMHMeshControl::CreateInterfaceElements()
{
    fCMesh->LoadReferences();
    int dim = fGMesh->Dimension();
    std::map<long, std::pair<long,long> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        long elindex = it->first;
        long leftelindex = it->second.first;
        long rightelindex = it->second.second;
        
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        int matid = gel->MaterialId();
        if (matid != fSkeletonMatId) {
            DebugStop();
        }
        std::set<long> celindices;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZCompElSide celskeleton = gelside.Reference();
        if(!celskeleton.Element()) DebugStop();
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            TPZStack<TPZCompElSide> celstack;
            gelside.HigherLevelCompElementList2(celstack, 1, 0);
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            TPZCompElSide tmp = gelside.LowerLevelCompElementList2(0);
            if (tmp) {
                celstack.Push(tmp);
            }
            long nstack = celstack.NElements();
            for (long i=0; i<nstack; i++) {
                TPZCompElSide celside = celstack[i];
                TPZGeoElSide gelneigh = celside.Reference();
                if (gelneigh.Dimension() != dim-1) {
                    continue;
                }
                long celindex = celside.Element()->Index();
                if (celindices.find(celindex) != celindices.end()) {
                    continue;
                }
                long gelneighindex = gelneigh.Element()->Index();
                int matid = 0;
                if (IsSibling(gelneighindex, leftelindex)) {
                    matid = fLagrangeMatIdLeft;
                }
                if (IsSibling(gelneighindex, rightelindex)) {
                    matid = fLagrangeMatIdRight;
                }
                // this means the element is a boundary. Retain its material id
                if (gelneigh.Element()->Dimension() < dim) {
                    matid = gelneigh.Element()->MaterialId();
                }
                celindices.insert(celindex);
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
    
    for (std::set<long>::iterator it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
        PrintSubdomain(*it, out);
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
