#include "TPZStructMatrixBase.h"

#include "pzcmesh.h"
#include "pzcheckmesh.h"
#include "pzerror.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"
#include "TPZThreadPool.h"
#include "TPZTimer.h"


#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrixOR"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
static LoggerPtr loggerGlobStiff(Logger::getLogger("pz.strmatrix.globalstiffness"));
#endif

TPZStructMatrixBase::TPZStructMatrixBase() : fMesh(NULL), fEquationFilter(0) {
    this->SetNumThreads(TPZThreadPool::globalInstance().threadCount());
}

TPZStructMatrixBase::TPZStructMatrixBase(TPZCompMesh *mesh)
    : fEquationFilter(0) {
    SetMesh(mesh);
    this->SetNumThreads(TPZThreadPool::globalInstance().threadCount());
}

TPZStructMatrixBase::TPZStructMatrixBase(TPZAutoPointer<TPZCompMesh> mesh)
    : fEquationFilter(0) {
    SetMesh(mesh);
    this->SetNumThreads(TPZThreadPool::globalInstance().threadCount());
}

TPZStructMatrixBase::TPZStructMatrixBase(const TPZStructMatrixBase &copy)
    : fMesh(copy.fMesh), fCompMesh(copy.fCompMesh),
      fEquationFilter(copy.fEquationFilter), fMaterialIds(copy.fMaterialIds),
      fNumThreads(copy.fNumThreads) {
}

void TPZStructMatrixBase::SetMesh(TPZCompMesh *mesh) {
    fMesh = mesh;
    fEquationFilter.SetNumEq(mesh ? mesh->NEquations() : 0);
#ifdef PZDEBUG
    if (mesh) {
        TPZCheckMesh checkmesh(fMesh, &std::cout);
        if (checkmesh.CheckConnectSeqNumberConsistency() != 0) {
            DebugStop();
        }
    }
#endif
}

void TPZStructMatrixBase::SetMesh(TPZAutoPointer<TPZCompMesh> mesh) {
    fCompMesh = mesh;
    SetMesh(mesh.operator->());
}

TPZMatrix<STATE> *TPZStructMatrixBase::CreateAssemble(
    TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    TPZMatrix<STATE> *stiff = Create();
    
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

void TPZStructMatrixBase::FilterEquations(TPZVec<int64_t> &origindex, TPZVec<int64_t> &destindex) const
{
    //destindex = origindex;
    fEquationFilter.Filter(origindex, destindex);
    
}

void TPZStructMatrixBase::SetMaterialIds(const std::set<int> &materialids)
{
    fMaterialIds = materialids;
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::set<int>::const_iterator it;
        std::stringstream sout;
        sout << "setting input material ids ";
        for(it=materialids.begin(); it!= materialids.end(); it++)
        {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    if(!fMesh)
    {
        LOGPZ_WARN(logger,"SetMaterialIds called without mesh")
        return;
    }
    int64_t iel;
    TPZAdmChunkVector<TPZCompEl*> &elvec = fMesh->ElementVec();
    int64_t nel = elvec.NElements();
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = elvec[iel];
        if(!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (cel);
        if(!subcmesh) continue;
        TPZAutoPointer<TPZAnalysis> analysis = subcmesh->Analysis();
        if(!analysis)
        {
            LOGPZ_ERROR(logger,"SetMaterialIds called for substructure without analysis object")
            DebugStop();
        }
        TPZStructMatrixBase *str = analysis->StructMatrix().operator->();
        if(!str)
        {
            LOGPZ_WARN(logger,"SetMaterialIds called for substructure without structural matrix")
            continue;
        }
        str->SetMaterialIds(materialids);
    }
}

int TPZStructMatrixBase::ClassId() const{
    return Hash("TPZStructMatrixBase");
}

void TPZStructMatrixBase::Read(TPZStream& buf, void* context) {
    fMesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::GetInstance(&buf));
    fCompMesh = TPZAutoPointerDynamicCast<TPZCompMesh>(TPZPersistenceManager::GetAutoPointer(&buf));
    fEquationFilter.Read(buf,context);
    buf.Read(fMaterialIds);
    buf.Read(&fNumThreads);
}

void TPZStructMatrixBase::Write(TPZStream& buf, int withclassid) const {
    TPZPersistenceManager::WritePointer(fMesh, &buf);
    TPZPersistenceManager::WritePointer(fCompMesh.operator ->(), &buf);
    fEquationFilter.Write(buf,withclassid);
    buf.Write(fMaterialIds);
    buf.Write(&fNumThreads);
}

void TPZStructMatrixBase::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
#ifdef PZDEBUG
    TExceptionManager activateExceptions;
#endif
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

    int64_t iel;
    int64_t nelem = fMesh->NElements();
    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
#ifdef LOG4CXX
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();

    int64_t count = 0;
    for (iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el) continue;
        int matid = 0;
        TPZGeoEl *gel = el->Reference();
        if (gel) {
            matid = gel->MaterialId();
        }
        int matidsize = fMaterialIds.size();
        if(matidsize){
            if(!el->NeedsComputing(fMaterialIds)) continue;
        }

        count++;
        if (!(count % 1000)) {
            std::cout << '*';
            std::cout.flush();
        }
        if (!(count % 20000)) {
            std::cout << "\n";
        }
        calcstiff.start();
        ek.Reset();
        ef.Reset();
        el->CalcStiff(ek, ef);
        if (guiInterface) if (guiInterface->AmIKilled()) {
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
                        /*          if(iel == 133)
                         {
                         std::stringstream sout2;
                         el->Reference()->Print(sout2);
                         el->Print(sout2);
                         LOGPZ_DEBUG(logger,sout2.str())
                         }
                         */
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
            //            TPZSFMatrix<STATE> test(stiffness);
            //            TPZFMatrix<STATE> test2(stiffness.Rows(),stiffness.Cols(),0.);
            //            stiffness.Print("before assembly",std::cout,EMathematicaInput);
            stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
#ifdef PZDEBUG
            REAL rhsnorm = Norm(ef.fMat);
            REAL eknorm = Norm(ek.fMat);
            if (std::isnan(rhsnorm) || std::isnan(eknorm)) {
                std::cout << "element " << iel << " has norm " << rhsnorm << std::endl;
                el->Print();
                ek.fMat.Print("ek",std::cout);
                ef.fMat.Print("ef",std::cout);
            }
#endif
            rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
            //            stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
            //            rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
            //            test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //            test -= stiffness;
            //            test.Print("matriz de rigidez diference",std::cout);
            //            test2.Print("matriz de rigidez interface",std::cout);
#ifdef LOG4CXX
            if (loggerel->isDebugEnabled()) {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                if (gel) {
                    TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                    gel->CenterPoint(gel->NSides() - 1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for computational element index " << el->Index() << " material id " << gel->MaterialId() << std::endl;
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                } else {
                    sout << "Stiffness for computational element without associated geometric element index " << el->Index() << "\n";
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
            if (loggerel->isDebugEnabled()) {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                //                el->Print();
                //                int nc = el->NConnects();
                //                for (int ic=0; ic<nc; ic++) {
                //                    std::cout << "Index " << el->ConnectIndex(ic) << " ";
                //                    el->Connect(ic).Print(*fMesh);
                //                    fMesh->ConnectVec()[ic].Print(*fMesh);
                //                }
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
        // tototototo
        //        GK.Multiply(Sol, GKSol);
        //        GKSol -= GF;
        //        GKSol.Transpose();
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "GKSol" << iel << " = ";
        //            GKSol.Print(str.str().c_str(),sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        stiffness.Multiply(Sol, GKSol);
        //        GKSol -= rhs;
        //        GKSol.Transpose();
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "StiffSol" << iel << " = ";
        //            GKSol.Print(str.str().c_str(),sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "GK" << iel << " = ";
        //            GK.Print(str.str().c_str(),sout,EMathematicaInput);
        //            std::stringstream str2;
        //            str2 << "ST" << iel << " = ";
        //            stiffness.Print(str2.str().c_str(),sout,EMathematicaInput);
        //            sout << "GK-ST\n";
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //
        //        stiffness.Zero();
        //        rhs.Zero();
        assemble.stop();
    }//fim for iel
    if (count > 1000) std::cout << std::endl;

#ifdef LOG4CXX
    if (loggerCheck->isDebugEnabled()) {
        std::stringstream sout;
        sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
    if (loggerGlobStiff->isDebugEnabled())
    {
        std::stringstream sout;
        stiffness.Print("GK = ",sout,EMathematicaInput);
        rhs.Print("GR = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerel,sout.str())
    }

#endif

}
