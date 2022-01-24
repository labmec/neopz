/**
 * @file pzcondensedcompel.cpp
 * @brief Contains the implementations of the TPZElementGroup methods.
 */

#include "pzelementgroup.h"
#include "TPZElementMatrixT.h"
#include "pzlog.h"
#include "pzstepsolver.h"
#include "pzcmesh.h"
#include <algorithm>

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzelementgroup");
#endif

TPZElementGroup::TPZElementGroup() : 
TPZRegisterClassId(&TPZElementGroup::ClassId),
TPZCompEl(), fElGroup(), fConnectIndexes()
{
}


TPZElementGroup::~TPZElementGroup()
{
    if (fElGroup.size()) {
        DebugStop();
    }
}

/** @brief create a copy of the condensed computational element in the other mesh */
TPZElementGroup::TPZElementGroup(TPZCompMesh &mesh, const TPZElementGroup &copy) : TPZRegisterClassId(&TPZElementGroup::ClassId),
TPZCompEl(mesh, copy)
{
    TPZStack<TPZCompEl *> newel;
    int nel = copy.fElGroup.size();
    for (int el=0; el<nel; el++) {
        newel.Push(copy.fElGroup[el]->Clone(mesh) );
    }
    for (int el=0; el<nel; el++) {
        AddElement(newel[el]);
    }
}

/** @brief add an element to the element group
 */
void TPZElementGroup::AddElement(TPZCompEl *cel)
{
    fElGroup.Push(cel);
    std::set<int64_t> connects;
    int nc = fConnectIndexes.size();
    for (int ic=0; ic<nc; ic++) {
        connects.insert(fConnectIndexes[ic]);
    }
    nc = cel->NConnects();
    for (int ic=0; ic<nc; ic++) {
        connects.insert(cel->ConnectIndex(ic));
    }
    nc = connects.size();
    if (nc != fConnectIndexes.size()) {
        fConnectIndexes.Resize(nc, 0);
        std::set<int64_t>::iterator it = connects.begin();
        for (int ic = 0; it != connects.end(); it++,ic++) {
            fConnectIndexes[ic] = *it;
        }
    }
    int64_t elindex = cel->Index();
    Mesh()->ElementVec()[elindex] = 0;
    std::list<TPZOneShapeRestraint> ellist = cel->GetShapeRestraints();
    for (std::list<TPZOneShapeRestraint>::iterator it=ellist.begin(); it != ellist.end(); it++) {
        int64_t cindex = it->fFaces[0].first;
        fRestraints[cindex] = *it;
    }
#ifdef PZ_LOG2
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Hiding element index " << elindex << " from the mesh data structure";
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
}

/** @brief unwrap the condensed element from the computational element and delete the condensed element */
void TPZElementGroup::Unwrap()
{
    int nel = fElGroup.size();
    for (int el=0; el<nel; el++) {
        int64_t elindex = fElGroup[el]->Index();
        Mesh()->ElementVec()[elindex] = fElGroup[el];
    }
    fElGroup.Resize(0);
    fConnectIndexes.Resize(0);
    delete this;
}

/**
 * @brief Set the index i to node inode
 * @param inode node to set index
 * @param index index to be seted
 */
void TPZElementGroup::SetConnectIndex(int inode, int64_t index)
{
    LOGPZ_ERROR(logger,"SetConnectIndex should never be called")
    DebugStop();
}


/// Reorder the connects in increasing number of elements connected
void TPZElementGroup::ReorderConnects()
{
    TPZManVector<std::pair<int,int64_t>, 100 > orderedindexes(fConnectIndexes.size());
    int nc = fConnectIndexes.size();
    std::set<int64_t> depreceive;
    // if a connect of the condensed element receives a contribution through connect dependency
    // this will be tracked by depreceive
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = Connect(ic);
        TPZConnect::TPZDepend * dep = c.FirstDepend();
        while (dep) {
            depreceive.insert(dep->fDepConnectIndex);
            dep = dep->fNext;
        }
    }
    for (int i=0; i<nc ; ++i) {
        TPZConnect &c = Connect(i);
        int64_t cindex = ConnectIndex(i);
        if(c.NElConnected() == 1 && c.HasDependency() == 0 && depreceive.find(cindex) == depreceive.end())
        {
            orderedindexes[i] = {0,cindex};
        }
        else
        {
            orderedindexes[i] = {1,cindex};
        }
    }
    std::sort(&orderedindexes[0],&orderedindexes[0]+nc);
    for(int ic=0; ic<nc; ic++) fConnectIndexes[ic] = orderedindexes[ic].second;
}

void TPZElementGroup::ReorderConnects(TPZManVector<int64_t> &connects)
{
    int64_t nc = connects.size();

    // TPZManVector<std::pair<int,int64_t>, 100 > orderedindexes(connects.size());
    for (int ic=0; ic<nc; ic++)
    {
        fConnectIndexes[ic] = connects[ic];
    }
}


/**
 * @brief Method for creating a copy of the element in a patch mesh
 * @param mesh Patch clone mesh
 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
 * @param gl2lcElMap map the computational elements
 */
/**
 * Otherwise of the previous clone function, this method don't
 * copy entire mesh. Therefore it needs to map the connect index
 * from the both meshes - original and patch
 */
TPZCompEl *TPZElementGroup::ClonePatchEl(TPZCompMesh &mesh,
                                std::map<int64_t,int64_t> & gl2lcConMap,
                                std::map<int64_t,int64_t> & gl2lcElMap) const
{    
    TPZElementGroup *result = new TPZElementGroup(mesh);
    int nel = fElGroup.size();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el]->ClonePatchEl(mesh,gl2lcConMap,gl2lcElMap);
        result->AddElement(cel);
    }
    return result;
}

void TPZElementGroup::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) const {
    InitializeElementMatrix(ek);
    int64_t rows = ek.Matrix().Rows();
    ek.Matrix().Redim(rows, rows);
    ek.fType = TPZElementMatrix::EK;
    InitializeElementMatrix(ef);
    std::map<int64_t,TPZOneShapeRestraint>::const_iterator it;
    for (it = fRestraints.begin(); it != fRestraints.end(); it++) {
        ek.fOneRestraints.push_back(it->second);
        ef.fOneRestraints.push_back(it->second);
    }
}//void

void TPZElementGroup::InitializeElementMatrix(TPZElementMatrix &ef) const {
	const int ncon = this->NConnects();
	int numeq = 0;
	ef.Block().SetNBlocks(ncon);
    for (int ic=0; ic<ncon; ic++) {
        TPZConnect &c = Connect(ic);
        int blsize = c.NShape()*c.NState();
        numeq += blsize;
        ef.Block().Set(ic, blsize);
    }
    const int numloadcases = 1;
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
	ef.Matrix().Redim(numeq,numloadcases);
	ef.fConnect.Resize(ncon);
	for(int i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
	}
    std::map<int64_t,TPZOneShapeRestraint>::const_iterator it;
    for (it = fRestraints.begin(); it != fRestraints.end(); it++) {
        ef.fOneRestraints.push_back(it->second);
    }
}//void


/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
template<class TVar>
void TPZElementGroup::CalcStiffInternal(TPZElementMatrixT<TVar> &ek,TPZElementMatrixT<TVar> &ef)
{
    std::map<int64_t,int64_t> locindex;
    int64_t ncon = fConnectIndexes.size();
    for (int64_t ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Calcstiff Element Group Index " << Index();
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    InitializeElementMatrix(ek, ef);
    int64_t nel = fElGroup.size();
    TPZElementMatrixT<TVar> ekloc,efloc;
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        
#ifdef PZ_LOG
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Assembling element " << el << " out of " << nel;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#ifdef PZDEBUG
        if(!cel){
            DebugStop();
        }
#endif
        
        cel->CalcStiff(ekloc, efloc);
#ifdef PZ_LOG
        if (logger.isDebugEnabled() ) {
            TPZGeoEl *gel = cel->Reference();
            
            int matid = 0;
            if(gel) matid = gel->MaterialId();
            std::stringstream sout;
            if (gel) {
                sout << "Material id " << matid <<std::endl;
            }
            else
            {
                sout << "No associated geometry\n";
            }
            sout << "Connect indexes ";
            for (int i=0; i<cel->NConnects(); i++) {
                sout << cel->ConnectIndex(i) << " ";
            }
            sout << std::endl;
            sout << "Local indexes ";
            for (int i=0; i<cel->NConnects(); i++) {
                sout << locindex[cel->ConnectIndex(i)] << " ";
            }
            sout << std::endl;
            ekloc.fMat.Print("EKElement =",sout,EMathematicaInput);
//            ekloc.fBlock.Print("EKBlock =",sout,&ekloc.fMat);
            efloc.fMat.Print("Vetor de carga",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
        
#endif
        int nelcon = ekloc.NConnects();
        for (int ic=0; ic<nelcon; ic++) {
            int iblsize = ekloc.fBlock.Size(ic);
            int icindex = ekloc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int idf = 0; idf<iblsize; idf++) {
                ef.fMat.at(ef.fBlock.at(ibldest,0,idf,0)) += efloc.fMat.at(efloc.fBlock.at(ic,0,idf,0));
            }
            for (int jc = 0; jc<nelcon; jc++) {
                int jblsize = ekloc.fBlock.Size(jc);
                int jcindex = ekloc.fConnect[jc];
                int jbldest = locindex[jcindex];
                for (int idf = 0; idf<iblsize; idf++) {
                    for (int jdf=0; jdf<jblsize; jdf++) {
                        ek.fMat.at(ek.fBlock.at(ibldest,jbldest,idf,jdf)) += ekloc.fMat.at(ekloc.fBlock.at(ic,jc,idf,jdf));
                    }
                }
            }

        }
    }
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "Connect indices " << fConnectIndexes << std::endl;
            //ek.fBlock.Print("EKBlockAssembled = ",sout,&ek.fMat);
            ek.fMat.Print("EKAssembled = ",sout,EMathematicaInput);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
}

/** @brief Verifies if the material associated with the element is contained in the set */
bool TPZElementGroup::HasMaterial(const std::set<int> &materialids) const {
    
    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
#ifdef PZDEBUG
        if(!cel){
            DebugStop();
        }
#endif
        bool has_material_Q = cel->HasMaterial(materialids);
        if (has_material_Q) {
            return true;
        }
    }
    return false;
}


/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
template<class TVar>
void TPZElementGroup::CalcResidualInternal(TPZElementMatrixT<TVar> &ef)
{
    std::map<int64_t,int64_t> locindex;
    int64_t ncon = fConnectIndexes.size();
    for (int64_t ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
    InitializeElementMatrix(ef);
    int64_t nel = fElGroup.size();
    TPZElementMatrixT<TVar> efloc;
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
#ifdef PZDEBUG
        if(!cel){
            DebugStop();
        }
#endif
        cel->CalcResidual(efloc);
        const int nelcon = efloc.NConnects();

        for (int ic=0; ic<nelcon; ic++) {
            int iblsize = efloc.fBlock.Size(ic);
            int icindex = efloc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int idf = 0; idf<iblsize; idf++) {
                ef.fMat.at(ef.fBlock.at(ibldest,0,idf,0)) += efloc.fMat.at(efloc.fBlock.at(ic,0,idf,0));
            }
        }        
    }
}

/**
 * @brief Performs an error computation on the element
 * @param errors [out] the L2 norm of the error of the solution
 * @param flux [in] value of the interpolated flux values
 */
void TPZElementGroup::EvaluateError(TPZVec<REAL> &errors, bool store_error)
{
    int nerr = errors.size();
    errors.Fill(0.);
    int meshdim = Mesh()->Dimension();
    for (auto cel : fElGroup)
    {
#ifdef PZDEBUG
        if(!cel){
            DebugStop();
        }
#endif
        TPZManVector<REAL,10> errloc(nerr,0.);
        cel->EvaluateError(errloc, store_error);
        if (errloc.size() != nerr) {
            nerr = errloc.size();
            errors.Resize(nerr, 0.);
        }
        for (int i=0; i<errloc.size(); i++) {
            errors[i] += errloc[i]*errloc[i];
        }
    }
    for (int i=0; i<errors.size(); i++) {
        errors[i] = sqrt(errors[i]);
    }
}

int TPZElementGroup::ClassId() const{
    return Hash("TPZElementGroup") ^ TPZCompEl::ClassId() << 1;
}

/** @brief Verifies if any element needs to be computed corresponding to the material ids */
bool TPZElementGroup::NeedsComputing(const std::set<int> &matids)
{
    bool result = false;
    for (int el=0; el<fElGroup.size(); el++) {
        result = fElGroup[el]->NeedsComputing(matids);
        if (result == true) {
            return result;
        }
    }
    return result;
}

