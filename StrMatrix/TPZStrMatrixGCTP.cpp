/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixGCTP methods.
 */

#include "TPZStrMatrixGCTP.h"

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"
#include "pzsfulmat.h"

#include "pzgnode.h"
#include "TPZTimer.h"
#include "TPZThreadTools.h"


#include "pzcheckconsistency.h"
#include "TPZMaterial.h"
#include <functional>

using namespace std;

#include "pzlog.h"

#include "pz_pthread.h"
#include "run_stats_table.h"
#include "TPZRenumbering.h"
#include "TPZThreadPool.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrixGCTP"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

RunStatsTable stat_ass_graph_gctp("-ass_graph_gctp", "Run statistics table for the graph creation and coloring TPZStructMatrixGCTP.");

TPZStructMatrixGCTP::TPZStructMatrixGCTP(TPZCompMesh *mesh) : TPZStructMatrixBase(mesh) {
    stat_ass_graph_gctp.start();
    TPZStructMatrixGCTP::OrderElement(this->Mesh(), fElementOrder);
    fNColors = TPZRenumbering::ColorElements(this->Mesh(), fElementOrder, fElementColors);
    stat_ass_graph_gctp.stop();
}

TPZStructMatrixGCTP::TPZStructMatrixGCTP(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrixBase(cmesh) {
    stat_ass_graph_gctp.start();
    TPZStructMatrixGCTP::OrderElement(this->Mesh(), fElementOrder);
    fNColors = TPZRenumbering::ColorElements(this->Mesh(), fElementOrder, fElementColors);
    stat_ass_graph_gctp.stop();
}

TPZStructMatrixGCTP::TPZStructMatrixGCTP(const TPZStructMatrixGCTP &copy) : TPZStructMatrixBase(copy),
fElementOrder(copy.fElementOrder),
fElementColors(copy.fElementColors),
fNColors(copy.fNColors) {
}

TPZMatrix<STATE> *TPZStructMatrixGCTP::Create() {
    cout << "TPZStructMatrixGCTP::Create should never be called\n";
    return 0;
}

TPZStructMatrixGCTP *TPZStructMatrixGCTP::Clone() {
    cout << "TPZStructMatrixGCTP::Clone should never be called\n";
    return 0;
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixGCTP::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    ass_stiff.start();
    if (fEquationFilter.IsActive()) {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense, rhs.Cols(), 0.);
        if (this->fNumThreads) {
            this->MultiThread_Assemble(stiffness, rhsloc, guiInterface);
        } else {
            this->Serial_Assemble(stiffness, rhsloc, guiInterface);
        }

        fEquationFilter.Scatter(rhsloc, rhs);
    } else {
        if (this->fNumThreads) {
            this->MultiThread_Assemble(stiffness, rhs, guiInterface);
        } else {
            this->Serial_Assemble(stiffness, rhs, guiInterface);
        }
    }
    ass_stiff.stop();
}

void TPZStructMatrixGCTP::Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    ass_rhs.start();
    if (fEquationFilter.IsActive()) {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
        int64_t neqexpand = fEquationFilter.NEqExpand();
        if (rhs.Rows() != neqexpand || Norm(rhs) != 0.) {
            DebugStop();
        }
        TPZFMatrix<STATE> rhsloc(neqcondense, rhs.Cols(), 0.);
        if (this->fNumThreads) {
            this->MultiThread_Assemble(rhsloc, guiInterface);
        } else {
            this->Serial_Assemble(rhsloc, guiInterface);
        }
        fEquationFilter.Scatter(rhsloc, rhs);
    } else {
        if (this->fNumThreads) {
            this->MultiThread_Assemble(rhs, guiInterface);
        } else {
            this->Serial_Assemble(rhs, guiInterface);
        }
    }
    ass_rhs.stop();
}

void TPZStructMatrixGCTP::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    if (!fMesh) {
        LOGPZ_ERROR(logger, "Serial_Assemble called without mesh")
        DebugStop();
    }
#ifdef LOG4CXX
    if (loggerelmat->isDebugEnabled()) {
        if (dynamic_cast<TPZSubCompMesh *> (fMesh)) {
            std::stringstream sout;
            sout << "AllEig = {};";
            LOGPZ_DEBUG(loggerelmat, sout.str())
        }
    }
#endif
#ifdef PZDEBUG
    if (rhs.Rows() != fEquationFilter.NActiveEquations()) {
        DebugStop();
    }
#endif

    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
#ifdef LOG4CXX
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();

    int64_t nelem = fMesh->NElements();
    int64_t count = 0;
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el) continue;
        int matidsize = fMaterialIds.size();
        if (matidsize) {
            TPZMaterial * mat = el->Material();
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
            if (!mat) {
                if (!submesh) {
                    continue;
                } else if (!submesh->NeedsComputing(fMaterialIds)) continue;
            } else {
                int matid = mat->Id();
                if (!this->ShouldCompute(matid)) continue;
            }
        }

        count++;
        if (!(count % 1000)) {
            std::cout << '*';
            std::cout.flush();
        }
        if (!(count % 20000)) {
            std::cout << std::endl;
        }
        calcstiff.start();

        ek.Reset(fMesh, TPZElementMatrix::EK);
        ef.Reset(fMesh, TPZElementMatrix::EF);
        el->CalcStiff(ek, ef);

        if (guiInterface && guiInterface->AmIKilled()) {
            return;
        }

#ifdef LOG4CXX
        if (loggerelmat->isDebugEnabled()) {
            if (dynamic_cast<TPZSubCompMesh *> (fMesh)) {
                std::stringstream objname;
                objname << "Element" << iel;
                std::string name = objname.str();
                objname << " = ";
                std::stringstream sout;
                ek.fMat.Print(objname.str().c_str(), sout, EMathematicaInput);
                sout << "AppendTo[AllEig,Eigenvalues[" << name << "]];";

                LOGPZ_DEBUG(loggerelmat, sout.str())
            }
        }

#endif

#ifdef CHECKCONSISTENCY
        //extern TPZCheckConsistency stiffconsist("ElementStiff");
        stiffconsist.SetOverWrite(true);
        bool result;
        result = stiffconsist.CheckObject(ek.fMat);
        if (!result) {
            globalresult = false;
            std::stringstream sout;
            sout << "element " << iel << " computed differently";
            LOGPZ_ERROR(loggerCheck, sout.str())
        }

#endif

        calcstiff.stop();
        assemble.start();

        if (!ek.HasDependency()) {
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
            rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);

#ifdef LOG4CXX
            if (loggerel->isDebugEnabled()) {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                if (gel) {
                    TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                    gel->CenterPoint(gel->NSides() - 1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for computational element index " << el->Index() << std::endl;
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                } else {
                    sout << "Stiffness for computational element without associated geometric element" << std::endl;
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel, sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
#ifdef LOG4CXX
            if (loggerel->isDebugEnabled() && !dynamic_cast<TPZSubCompMesh *> (fMesh)) {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                if (gel) {
                    TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                    gel->CenterPoint(gel->NSides() - 1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                } else {
                    sout << "Stiffness for computational element index " << iel << std::endl;
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel, sout.str())
            }
#endif
        }
        assemble.stop();
    }

#ifdef LOG4CXX
    if (loggerCheck->isDebugEnabled()) {
        std::stringstream sout;
        sout << "The comparison results are : consistency check " << globalresult << " write read check " << writereadresult;
        stiffness.Print("Matriz de Rigidez: ", sout, EMathematicaInput);
        rhs.Print("Right Hand side", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }

#endif

}

void TPZStructMatrixGCTP::Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> /*guiInterface*/) {
    TPZTimer calcresidual("Computing the residual vector");
    TPZTimer assemble("Assembling the residual vector");

    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();

    TPZElementMatrix ef(fMesh, TPZElementMatrix::EF);
    
    for (auto iel : fElementOrder) {
        TPZCompEl *el = elementvec[iel];
        if (!el) continue;

        TPZMaterial * mat = el->Material();
        if (!mat) continue;
        int matid = mat->Id();
        if (!this->ShouldCompute(matid)) continue;


        calcresidual.start();
        ef.Reset(fMesh, TPZElementMatrix::EF);
        el->CalcResidual(ef);

        calcresidual.stop();

        assemble.start();

        if (!ef.HasDependency()) {
            ef.ComputeDestinationIndices();
            fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ef.ComputeDestinationIndices();
            fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat, ef.fSourceIndex, ef.fDestinationIndex);
        }

        assemble.stop();

    }
#ifdef LOG4CXX
    {
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << calcresidual.processName() << " " << calcresidual << std::endl;
            sout << assemble.processName() << " " << assemble;
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
    }
#endif
}

TPZMatrix<STATE> * TPZStructMatrixGCTP::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    TPZMatrix<STATE> *stiff = Create();

    //int64_t neq = stiff->Rows();
    int64_t cols = MAX(1, rhs.Cols());
    rhs.Redim(fEquationFilter.NEqExpand(), cols);
    Assemble(*stiff, rhs, guiInterface);

#ifdef LOG4CXX2
    if (loggerel->isDebugEnabled()) {
        std::stringstream sout;
        stiff->Print("Stiffness matrix", sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel, sout.str())
    }
#endif
    return stiff;

}

void TPZStructMatrixGCTP::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    if (guiInterface && guiInterface->AmIKilled()) {
        return;
    }

    auto cmesh = this->Mesh();

    std::function<void(TPZStack<int64_t>)> job = [this, cmesh, &mat, &rhs, &guiInterface] (TPZStack<int64_t> indices) {
        for (auto index : indices){
            int64_t iel = fElementOrder[index];
            if (iel >= 0) {
                TPZCompEl *el = cmesh->ElementVec()[iel];
                if (!el) {
                    return;
                }
                TPZElementMatrix ek(cmesh, TPZElementMatrix::EK);
                TPZElementMatrix ef(cmesh, TPZElementMatrix::EF);

                el->CalcStiff(ek, ef);
                if (guiInterface && guiInterface->AmIKilled()) {
                    return;
                }

                if (!el->HasDependency()) {
                    ek.ComputeDestinationIndices();
                    if (this->EquationFilter().IsActive()){
                        this->FilterEquations(ek.fSourceIndex, ek.fDestinationIndex);
                    }
                } else {
                    // the element has dependent nodes
                    ek.ApplyConstraints();
                    ef.ApplyConstraints();
                    ek.ComputeDestinationIndices();
                    if (this->EquationFilter().IsActive()){
                        this->FilterEquations(ek.fSourceIndex, ek.fDestinationIndex);
                    }
                }

                // Assemble matrix and rhs
                if (!ek.HasDependency()) {
                    mat.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
                    rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
                } else {
                    mat.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
                    rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
                }

            } else {
                std::cout << "the element in ElColorSequence is negative???\n";
                DebugStop();
            }
        }
    };

    int64_t n_elements = fElementOrder.size();
    for (int64_t color = 0; color < fNColors; ++color) {
        TPZTaskGroup colorGroup;
        TPZManVector<TPZStack<int64_t>,10> indices(fNumThreads);
        int64_t elem_count = 0;
        for (int64_t i = 0; i < n_elements; ++i) {
            if (fElementColors[i] == color) {
                indices[(elem_count++)%fNumThreads].push_back(i);
            }
        }
        for (auto &ind_vec:indices) {
            TPZThreadPool::globalInstance().run(1, &colorGroup, job, ind_vec);
        }
        colorGroup.Wait();
    }

#ifdef LOG4CXX
    if (loggerCheck->isDebugEnabled()) {
        std::stringstream sout;
        mat.Print("Matriz de Rigidez: ", sout, EMathematicaInput);
        rhs.Print("Right hand side", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
#endif
}

void TPZStructMatrixGCTP::MultiThread_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    if (guiInterface && guiInterface->AmIKilled()) {
        return;
    }

    auto cmesh = this->Mesh();

    std::function<void(TPZStack<int64_t>) > job = [this, cmesh, &rhs, &guiInterface] (TPZStack<int64_t> indices) {
        for (auto index: indices){
            int64_t iel = fElementOrder[index];
            if (iel >= 0) {
                TPZCompEl *el = cmesh->ElementVec()[iel];
                if (!el) {
                    return;
                }
                TPZElementMatrix ef(cmesh, TPZElementMatrix::EF);

                el->CalcResidual(ef);
                if (guiInterface && guiInterface->AmIKilled()) {
                    return;
                }

                if (!el->HasDependency()) {
                    ef.ComputeDestinationIndices();
                    if (this->EquationFilter().IsActive()){
                        this->FilterEquations(ef.fSourceIndex, ef.fDestinationIndex);
                    }
                } else {
                    // the element has dependent nodes
                    ef.ApplyConstraints();
                    ef.ComputeDestinationIndices();
                    if (this->EquationFilter().IsActive()){
                        this->FilterEquations(ef.fSourceIndex, ef.fDestinationIndex);
                    }
                }

                // Assemble rhs
                if (!ef.HasDependency()) {
                    rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
                } else {
                    rhs.AddFel(ef.fConstrMat, ef.fSourceIndex, ef.fDestinationIndex);
                }

            } else {
                std::cout << "the element in ElColorSequence is negative???\n";
                DebugStop();
            }
        }
    };

    int64_t n_elements = fElementOrder.size();
    for (int64_t color = 0; color < fNColors; ++color) {
        TPZTaskGroup colorGroup;
        TPZManVector<TPZStack<int64_t>,10> indices(fNumThreads);
        int64_t elem_count = 0;
        for (int64_t i = 0; i < n_elements; ++i) {
            if (fElementColors[i] == color) {
                indices[(elem_count++)%fNumThreads].push_back(i);
            }
        }
        for (auto &ind_vec:indices) {
            TPZThreadPool::globalInstance().run(1, &colorGroup, job, ind_vec);
        }
        colorGroup.Wait();
    }

#ifdef LOG4CXX
    if (loggerCheck->isDebugEnabled()) {
        std::stringstream sout;
        rhs.Print("Right hand side", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
#endif
}

void TPZStructMatrixGCTP::OrderElement(TPZCompMesh *cmesh, TPZVec<int64_t> &ElementOrder) {
    int64_t numelconnected = 0;
    int64_t nconnect = cmesh->NConnects();
    //firstelconnect contains the first element index in the elconnect vector
    TPZVec<int64_t> firstelconnect(nconnect + 1);
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        numelconnected += cmesh->ConnectVec()[ic].NElConnected();
        firstelconnect[ic + 1] = firstelconnect[ic] + cmesh->ConnectVec()[ic].NElConnected();
    }
    TPZVec<int64_t> elconnect(numelconnected, -1);
    int64_t el;
    TPZCompEl *cel;

#ifdef NOORDER
    int64_t count = 0;
    std::cout << __PRETTY_FUNCTION__ << " no element order\n";
    ElementOrder.Resize(cmesh->NElements(), -1);
    for (el = 0; el < cmesh->ElementVec().NElements(); el++) {
        cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        ElementOrder[count] = el;
        count++;
    }
    ElementOrder.Resize(count);
    return;
#endif
    for (el = 0; el < cmesh->ElementVec().NElements(); el++) {
        cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int64_t nc = connectlist.NElements();
        for (int64_t ic = 0; ic < nc; ic++) {
            int64_t cindex = connectlist[ic];
            elconnect[firstelconnect[cindex]] = el;
            firstelconnect[cindex]++;
        }
    }
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        firstelconnect[ic + 1] = firstelconnect[ic] + cmesh->ConnectVec()[ic].NElConnected();
    }

    ElementOrder.Resize(cmesh->ElementVec().NElements(), -1);
    ElementOrder.Fill(-1);
    TPZVec<int64_t> nodeorder(cmesh->ConnectVec().NElements(), -1);
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        int64_t seqnum = cmesh->ConnectVec()[ic].SequenceNumber();
        if (seqnum >= 0) nodeorder[seqnum] = ic;
    }
    int64_t elsequence = 0;
    TPZVec<int> elorderinv(cmesh->ElementVec().NElements(), -1);
    for (int64_t seq = 0; seq < nconnect; seq++) {
        int64_t ic = nodeorder[seq];
        if (ic == -1) continue;
        int64_t firstind = firstelconnect[ic];
        int64_t lastind = firstelconnect[ic + 1];
        for (int64_t ind = firstind; ind < lastind; ind++) {
            el = elconnect[ind];
            if (el == -1) {
                continue;
            }
            if (elorderinv[el] == -1) elorderinv[el] = elsequence++;
        }
    }
    elsequence = 0;
    for (int64_t seq = 0; seq < cmesh->ElementVec().NElements(); seq++) {
        if (elorderinv[seq] == -1) continue;
        ElementOrder[elorderinv[seq]] = seq;
    }
    int64_t seq;
    for (seq = 0; seq < cmesh->ElementVec().NElements(); seq++) {
        if (ElementOrder[seq] == -1) break;
    }

    ElementOrder.Resize(seq);
}

int TPZStructMatrixGCTP::ClassId() const {
    return Hash("TPZStructMatrixGCTP") ^ TPZStructMatrixBase::ClassId() << 1;
}

void TPZStructMatrixGCTP::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf, context);

    buf.Read(fElementOrder);
    buf.Read(fElementColors);
    buf.Read(&fNColors);
}

void TPZStructMatrixGCTP::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf, withclassid);

    buf.Write(fElementOrder);
    buf.Write(fElementColors);
    buf.Write(&fNColors);
}

template class TPZRestoreClass<TPZStructMatrixGCTP>;