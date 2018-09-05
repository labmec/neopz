/**
 * @file pzcondensedcompel.cpp
 * @brief Contains the implementations of the TPZElementGroup methods.
 */

#include "pzelementgroup.h"
#include "pzlog.h"
#include "pzstepsolver.h"
#include "pzcmesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzelementgroup"));
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
#ifdef LOG4CXX2
    if (logger->isDebugEnabled())
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
    int64_t index;
    TPZElementGroup *result = new TPZElementGroup(mesh,index);
    int nel = fElGroup.size();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el]->ClonePatchEl(mesh,gl2lcConMap,gl2lcElMap);
        result->AddElement(cel);
    }
    return result;
}

void TPZElementGroup::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) const {
    InitializeElementMatrix(ek);
    int64_t rows = ek.fMat.Rows();
    ek.fMat.Redim(rows, rows);
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
	ef.fBlock.SetNBlocks(ncon);
    for (int ic=0; ic<ncon; ic++) {
        TPZConnect &c = Connect(ic);
        int blsize = c.NShape()*c.NState();
        numeq += blsize;
        ef.fBlock.Set(ic, blsize);
    }
    const int numloadcases = 1;
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
	ef.fMat.Redim(numeq,numloadcases);
	ef.fNumStateVars = 1;
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
void TPZElementGroup::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    std::map<int64_t,int64_t> locindex;
    int64_t ncon = fConnectIndexes.size();
    for (int64_t ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
    InitializeElementMatrix(ek, ef);
    int64_t nel = fElGroup.size();
    TPZElementMatrix ekloc,efloc;
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        cel->CalcStiff(ekloc, efloc);
#ifdef LOG4CXX
        if (logger->isDebugEnabled() ) {
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
            efloc.Print(sout);
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
                ef.fBlock(ibldest,0,idf,0) += efloc.fBlock(ic,0,idf,0);
            }
            for (int jc = 0; jc<nelcon; jc++) {
                int jblsize = ekloc.fBlock.Size(jc);
                int jcindex = ekloc.fConnect[jc];
                int jbldest = locindex[jcindex];
                for (int idf = 0; idf<iblsize; idf++) {
                    for (int jdf=0; jdf<jblsize; jdf++) {
                        ek.fBlock(ibldest,jbldest,idf,jdf) += ekloc.fBlock(ic,jc,idf,jdf);
                    }
                }
            }

        }
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "Connect indices " << fConnectIndexes << std::endl;
            ek.fBlock.Print("EKBlockAssembled = ",sout,&ek.fMat);
            ek.fMat.Print("EKAssembled = ",sout,EMathematicaInput);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
    }
}


/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
void TPZElementGroup::CalcResidual(TPZElementMatrix &ef)
{
    std::map<int64_t,int64_t> locindex;
    int64_t ncon = fConnectIndexes.size();
    for (int64_t ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
    InitializeElementMatrix(ef);
    int64_t nel = fElGroup.size();
    TPZElementMatrix efloc;
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        cel->CalcResidual(efloc);
        
        int nelcon = efloc.NConnects();
        for (int ic=0; ic<nelcon; ic++) {
            int iblsize = efloc.fBlock.Size(ic);
            int icindex = efloc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int idf = 0; idf<iblsize; idf++) {
                ef.fBlock(ibldest,0,idf,0) += efloc.fBlock(ic,0,idf,0);
            }
        }        
    }
}

/**
 * @brief Performs an error estimate on the elemen
 * @param fp function pointer which computes the exact solution
 * @param errors [out] the L2 norm of the error of the solution
 * @param flux [in] value of the interpolated flux values
 */
void TPZElementGroup::EvaluateError(std::function<void(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> fp,
                           TPZVec<REAL> &errors, bool store_error)
{
    int nerr = errors.size();
    errors.Fill(0.);
    int meshdim = Mesh()->Dimension();
    int nel = fElGroup.size();
    for (int el=0; el<nel; el++) {
        TPZManVector<REAL,10> errloc(nerr,0.);
        TPZGeoEl *elref = fElGroup[el]->Reference();
        if (elref && elref->Dimension() != meshdim) {
            continue;
        }
        fElGroup[el]->EvaluateError(fp, errloc, store_error);
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

