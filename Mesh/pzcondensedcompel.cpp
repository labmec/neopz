/**
 * @file pzcondensedcompel.cpp
 * @brief Contains the implementations of the TPZCondensedCompEl methods.
 */

#include "pzcondensedcompel.h"
#include "pzlog.h"
#include "pzstepsolver.h"
#include "pzelementgroup.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcondensedcompel"));
#endif

#ifdef USING_LAPACK
#define USING_DGER2
#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#endif
#endif

TPZCondensedCompEl::TPZCondensedCompEl(TPZCompEl *ref) :
TPZRegisterClassId(&TPZCondensedCompEl::ClassId)
{
    if(!ref)
    {
        DebugStop();
    }
    fReferenceCompEl = ref;
    fMesh = ref->Mesh();
    TPZGeoEl *gel = ref->Reference();
    if (gel) {
        SetReference(gel->Index());
    }
    SetIndex(ref->Index());
    fMesh->ElementVec()[fIndex] = this;
    int ncon = ref->NConnects();
    fIndexes.resize(ncon);
    for (int ic=0; ic< ncon ; ic++) {
        fIndexes[ic] = ic;
    }
    Resequence();
}


TPZCondensedCompEl::~TPZCondensedCompEl()
{
    if (fMesh->ElementVec()[fIndex] == this) {
        fMesh->ElementVec()[fIndex] = fReferenceCompEl;
        delete fReferenceCompEl;
    }
}

/** @brief create a copy of the condensed computational element in the other mesh */
TPZCondensedCompEl::TPZCondensedCompEl(const TPZCondensedCompEl &copy, TPZCompMesh &mesh) : TPZRegisterClassId(&TPZCondensedCompEl::ClassId)
{
    TPZCompEl *ref = fReferenceCompEl->Clone(mesh);
    if(!ref)
    {
        DebugStop();
    }
    fReferenceCompEl = ref;
    fMesh = ref->Mesh();
    SetReference(ref->Reference()->Index());
    SetIndex(ref->Index());
    fMesh->ElementVec()[fIndex] = this;
    int ncon = ref->NConnects();
    fIndexes.resize(ncon);
    for (int i=0; i<ncon ; ++i) {
        fIndexes[i] = i;
    }
}


/** @brief unwrap the condensed element from the computational element and delete the condensed element */
void TPZCondensedCompEl::Unwrap()
{
    long myindex = fIndex;
    fMesh->ElementVec()[myindex] = 0;
    TPZCompEl *ReferenceEl = fReferenceCompEl;
    int ncon = NConnects();
    for (int ic=0; ic<ncon ; ic++) {
        Connect(ic).SetCondensed(false);
    }
    TPZGeoEl *gel = Reference();
    if (gel) {
        gel->ResetReference();
    }
    delete this;
    ReferenceEl->Mesh()->ElementVec()[myindex] = ReferenceEl;
    if (gel) {
        gel->SetReference(ReferenceEl);
    }
}

/**
 * @brief Set the index i to node inode
 * @param inode node to set index
 * @param index index to be seted
 */
void TPZCondensedCompEl::SetConnectIndex(int inode, long index)
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
TPZCompEl *TPZCondensedCompEl::ClonePatchEl(TPZCompMesh &mesh,
                                std::map<long,long> & gl2lcConMap,
                                std::map<long,long> & gl2lcElMap) const
{
    TPZCompEl *cel = fReferenceCompEl->ClonePatchEl(mesh,gl2lcConMap,gl2lcElMap);
    TPZCondensedCompEl *result = new TPZCondensedCompEl(cel);
    result->fIndexes = fIndexes;
    result->fCondensed = fCondensed;
    return result;
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param sol finite element solution
 * @param dsol solution derivatives
 * @param axes axes associated with the derivative of the solution
 */
void TPZCondensedCompEl::ComputeSolution(TPZVec<REAL> &qsi,
                             TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes)
{
    fReferenceCompEl->ComputeSolution(qsi,sol,dsol,axes);
}

void TPZCondensedCompEl::ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data)
{
    fReferenceCompEl->ComputeSolution(qsi,data);
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
 * This method will function for both volumetric and interface elements
 * @param qsi master element coordinate of the interface element
 * @param normal vector
 * @param leftsol finite element solution
 * @param dleftsol solution derivatives
 * @param leftaxes axes associated with the left solution
 * @param rightsol finite element solution
 * @param drightsol solution derivatives
 * @param rightaxes axes associated with the right solution
 */
void TPZCondensedCompEl::ComputeSolution(TPZVec<REAL> &qsi,
                             TPZVec<REAL> &normal,
                             TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
                             TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes)
{
    fReferenceCompEl->ComputeSolution(qsi,normal,leftsol,dleftsol,leftaxes,rightsol,drightsol,rightaxes);
}

/**
 * @brief Computes solution and its derivatives in local coordinate qsi
 * @param qsi master element coordinate
 * @param phi matrix containing shape functions compute in qsi point
 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
 * @param axes [in] axes indicating the direction of the derivatives
 * @param sol finite element solution
 * @param dsol solution derivatives
 */
void TPZCondensedCompEl::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                             const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol)
{
    fReferenceCompEl->ComputeSolution(qsi,phi,dphix,axes,sol,dsol);
}

void TPZCondensedCompEl::Resequence()
{
    TPZStack<int> condensed;
    TPZStack<int> notcondensed;
    int nint=0,next=0;
    std::set<long> depreceive;
    int ncon = NConnects();
    for (int ic=0; ic<ncon; ic++) {
        TPZConnect &c = Connect(ic);
        TPZConnect::TPZDepend * dep = c.FirstDepend();
        while (dep) {
            depreceive.insert(dep->fDepConnectIndex);
            dep = dep->fNext;
        }
    }
    for (int i=0; i<ncon ; ++i) {
        TPZConnect &c = Connect(i);
        long cindex = ConnectIndex(i);
        if(c.NElConnected() == 1 && c.HasDependency() == 0 && depreceive.find(cindex) == depreceive.end())
        {
            c.SetCondensed(true);
        }
        if (c.IsCondensed()) {
            if(c.HasDependency())
            {
                std::cout << "Not Implemented yet\n";
                DebugStop();
            }
            condensed.Push(i);
            nint += c.NDof();
        }
        else
        {
            notcondensed.Push(i);
            next += c.NDof();
        }
    }
    int ncond = condensed.size();
    for (int i=0; i<ncond; ++i) {
        fIndexes[i] = condensed[i];
    }
    for (int i=0; i<notcondensed.size(); ++i) {
        fIndexes[i+ncond] = notcondensed[i];
    }
    //TPZAutoPointer<TPZMatrix<STATE> > k00 = new TPZFMatrix<STATE>(nint,nint,0.);
	TPZAutoPointer<TPZMatrix<STATE> > k00 = new TPZFMatrix<STATE>(nint, nint, 0.);
    //TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(k00);
    TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(k00);
    TPZAutoPointer<TPZMatrix<STATE> > mat2 = k00->Clone();
    
    TPZStepSolver<STATE> *gmrs = new TPZStepSolver<STATE>(mat2);
    step->SetReferenceMatrix(mat2);
    step->SetDirect(ELDLt);
    gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
    TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
    TPZAutoPointer<TPZMatrixSolver<STATE> > autogmres = gmrs;
    

    
    fCondensed.SetSolver(autostep);
    fCondensed.Redim(nint+next,nint);
}

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZCondensedCompEl::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    fReferenceCompEl->CalcStiff(ek,ef);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        Print(sout);
        sout << "Connect indices of element stiffness" << ek.fConnect << std::endl;
        //ek.fMat.Print("EKOrig = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    ek.PermuteGather(fIndexes);
    ef.PermuteGather(fIndexes);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Element connects\n";
        int nc = NConnects();
        for (int ic=0; ic<nc; ic++) {
            sout << "ic = " << ic << ' ' << " index " << ConnectIndex(ic) << ' ';
            Connect(ic).Print(*Mesh(),sout);
        }
        sout << "Permutations " << fIndexes << std::endl;
        sout << "Connect indices " << ek.fConnect << std::endl;

        ek.fMat.Print("EKPermute = ",sout,EMathematicaInput);
        ef.fMat.Print("EFPermute = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#ifdef USING_DGER
    
    fCondensed.Zero();
    
    long dim0 = fCondensed.Dim0();
    long dim1 = fCondensed.Dim1();
    long rows = ek.fMat.Rows();
    long cols = ek.fMat.Cols()+ef.fMat.Cols();
    
    
    TPZFMatrix<STATE> KF(rows,cols);
    
    for(long i=0; i<rows;i++) // Montando a matriz KF a partir de ek e ef
        for (long j=0; j<cols; j++)
        {
            if (j<rows)
                KF(i,j) = ek.fMat(i,j);
            else
                KF(i,j) = ef.fMat(i,j-rows);
        }
    
    for (long i=0; i<dim1; i++) // Aplicando os valores na matriz K10 de fCondensed
        for (long j=0; j<dim0; j++)
            fCondensed.K10().operator()(i,j)=KF(i+dim0,j);
    
    for (long i=0; i<dim1; i++) // Aplicando os valores na matriz K11 de fCondensed
        for (long j=0; j<dim1; j++)
            fCondensed(i+dim0,j+dim0)=KF(i+dim0,j+dim0);
    
    fCondensed.SetF(ef.fMat); // Definindo ek.fMat como vetor de forcas de fCondensed
    
    
    for (long i=0; i<rows-dim1; i++) // Realização da condensação estática em KF usando BLAS
    {
        for(long j=i+1;j<cols;j++) // DECOMPOSIÇÃO LDLt
        {
            if (j<rows)
            {
                KF(j,i)/=KF(i,i);
                KF(i,j)/=KF(i,i);
            }
            else
                KF(i,j)/=KF(i,i);
        }
        
#ifdef STATEdouble
        cblas_dger (CblasColMajor, rows-i-1, cols-i-1,
                    -KF(i,i), &KF(i+1,i), 1,
                    &KF(i,i+1), rows, &KF(i+1,i+1), rows);
#else
        DebugStop();
#endif
    }
    
    
    for (long i=dim0; i< rows; i++) // Aplicando matriz e vetor forças condensadas em ek e ef
    {
        ef.fMat(i,0) = KF.GetVal(i,ek.fMat.Rows());
        for (long j=dim0; j<ek.fMat.Rows(); j++)
        {
            ek.fMat(i,j) = KF.GetVal(i,j);
        }
    }
    
    TPZAutoPointer<TPZMatrix<STATE> > K00 = fCondensed.K00();
    for (long i=0; i<dim0; i++) // Substituindo novos valores de K00 e definindo-o como decomposto
        for (long j=0; j<dim0; j++)
            K00->operator()(i, j) = KF(i,j);
    
    fCondensed.K00()->SetIsDecomposed(ELDLt);
    
    
    for (long i=0; i<dim0; i++) // Substituindo valores obtidos para K01 usando o BLAS
        for (long j=0; j<dim1; j++)
            fCondensed.K01().operator()(i,j) = KF(i,j+dim0);
    
    fCondensed.K00()->Subst_LBackward(&fCondensed.K01()); //Com SubstL_Back chegamos ao K01 desejado
    
//    void cblas_dtrsm(const enum CBLAS_ORDER __Order, const enum CBLAS_SIDE __Side,
//                     const enum CBLAS_UPLO __Uplo, const enum CBLAS_TRANSPOSE __TransA,
//                     const enum CBLAS_DIAG __Diag, const int __M, const int __N,
//                     const double __alpha, const double *__A, const int __lda, double *__B,
//                     const int __ldb) __OSX_AVAILABLE_STARTING(__MAC_10_2,__IPHONE_4_0);

    cblas_dtrsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasUnit,)
    fCondensed.SetK01IsComputed(1);
    fCondensed.SetReduced();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        ek.fMat.Print("Rigidez",sout);
        ef.fMat.Print("Carga",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    /*******************USING MATRED***************************/
    
#else
    
    fCondensed.Zero();
    
    long dim = ek.fMat.Rows();
    for (long i=0; i<dim ; ++i) {
        for (long j=0; j<dim ; ++j) {
            fCondensed(i,j) = ek.fMat(i,j);
        }
    }
    
    fCondensed.SetF(ef.fMat);
    
    
    long dim1 = fCondensed.Dim1();
    TPZFNMatrix<200,STATE> K11(dim1,dim1),F1(dim1,ef.fMat.Cols());
    //const TPZFMatrix<REAL> &k11 = fCondensed.K11Red();
    
	fCondensed.K11Reduced(K11, F1);

#ifdef LOG4CXX
    if(logger->isDebugEnabled() && (Index() == 927 || Index() == 923))
    {
        std::stringstream sout;
        sout << "Index = " << Index() << std::endl;
        fCondensed.Print("Reduced = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    fCondensed.SetReduced();
    
    //const TPZFMatrix<REAL> &f1 = fCondensed.F1Red();
    long dim0 = dim-K11.Rows();
    for (long i=dim0; i<dim; i++) {
        ef.fMat(i,0) = F1.GetVal(i-dim0,0);
        for (long j=dim0; j<dim; j++) {
            ek.fMat(i,j) = K11.GetVal(i-dim0,j-dim0);
        }
    }
    
    
#endif
    
#ifdef USING_DGER
#ifdef USING_LAPACK
    TPZFMatrix<STATE> * K00_temp = dynamic_cast<TPZFMatrix<STATE> * >(fCondensed.K00().operator->());
    K00_temp->InitializePivot();
#endif
#endif
    
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        int nc = NConnects();
        for (int ic=0; ic<nc; ic++) {
            sout << "ic = " << ic << ' ';
            Connect(ic).Print(*Mesh(),sout);
        }
        ek.fMat.Print("EK11Reduced",sout,EMathematicaInput);
        ef.fMat.Print("EF11Reduced",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}


/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
void TPZCondensedCompEl::CalcResidual(TPZElementMatrix &ef)
{
    fReferenceCompEl->CalcResidual(ef);
    ef.PermuteGather(fIndexes);
    fCondensed.SetF(ef.fMat);
    //const TPZFMatrix<REAL> &f1 = fCondensed.F1Red();
    TPZFNMatrix<100,STATE> f1(fCondensed.Dim1(),ef.fMat.Cols());
	fCondensed.F1Red(f1);
    long dim1 = f1.Rows();
    long dim = ef.fMat.Rows();
    long dim0 = dim-dim1;
    for (long i= dim0; i<dim; i++) {
        ef.fMat(i,0) = f1.GetVal(i-dim0,0);
    }
}

/**
 * @brief Prints element data
 * @param out Indicates the device where the data will be printed
 */
void TPZCondensedCompEl::Print(std::ostream &out) const
{
    out << "Output for a condensed element\n";
    TPZCompEl::Print(out);
    
    out << "Index of grouped elements: ";
    TPZElementGroup *eg = dynamic_cast<TPZElementGroup *>(fReferenceCompEl);
    if(eg)
    {
        int nel = eg->GetElGroup().size();
        for(int i=0; i<nel-1; i++){
            out << eg->GetElGroup()[i]->Index() <<", ";
        }
        out << eg->GetElGroup()[nel-1]->Index() <<std::endl;
        out << "Connect indexes of the contained elements\n";
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = eg->GetElGroup()[i];
            TPZGeoEl *gel = cel->Reference();
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                out << cel->ConnectIndex(ic) << " ";
            }
            if (gel) {
                out << "matid " << gel->MaterialId();
                TPZManVector<REAL,3> xi(gel->Dimension()), xco(3);
                gel->CenterPoint(gel->NSides()-1, xi);
                gel->X(xi,xco);
                out << " center x " << xco;
            }
            out << std::endl;
        }
    }
    else
    {
        if (fReferenceCompEl) {
            fReferenceCompEl->Print(out);
        }
        else
        {
            DebugStop();
        }
    }
    out << "Internal index resequencing: " << fIndexes << std::endl;
    fCondensed.Print("Condensed matrix",out);
}


/** @brief Loads the solution within the internal data structure of the element */
/**
 * Is used to initialize the solution of connect objects with dependency \n
 * Is also used to load the solution within SuperElements
 */
void TPZCondensedCompEl::LoadSolution()
{
    // initialize the solution of the constrained connects
    TPZCompEl::LoadSolution();
    
    // if the matrix has not been condensed then nothing to do
    if (fCondensed.Rows() != fCondensed.Dim1())
    {
        return;
    }
    // compute the solution of the internal equations
    int dim0=0, dim1=0;
    int nc = NConnects(),nc0 = 0, nc1 = 0;
    int ic;
    for (ic=0; ic<nc ; ic++) {
        long connectindex = ConnectIndex(ic);
        TPZConnect &con = Connect(ic);
        int sz = con.NShape()*con.NState();
        if (con.IsCondensed()) {
#ifdef PZDEBUG
            if (dim1) {
                DebugStop();
            }
#endif
            dim0 += sz;
            nc0++;
        }
        else
        {
            dim1 += sz;
            nc1++;
        }
    }
    //TPZBlock<REAL> &bl = Mesh()->Block();
	TPZBlock<STATE> &bl = Mesh()->Block();
    long count = 0;
    //TPZFMatrix<REAL> u1(dim1,1,0.);
	TPZFMatrix<STATE> u1(dim1,1,0.);
    //TPZFMatrix<REAL> elsol(dim0+dim1,1,0.);
	TPZFMatrix<STATE> elsol(dim0+dim1,1,0.);
    for (ic=nc0; ic<nc ; ic++) {
        TPZConnect &c = Connect(ic);
        long seqnum = c.SequenceNumber();
        int blsize = bl.Size(seqnum);
        for (int ibl=0; ibl<blsize; ibl++) {
            u1(count++,0) = bl(seqnum,0,ibl,0);
        }
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        LOGPZ_DEBUG(logger, "Computing UGlobal")
    }
#endif
    fCondensed.UGlobal(u1, elsol);
    count = 0;
    for (ic=0; ic<nc0 ; ic++) {
        TPZConnect &c = Connect(ic);
        long seqnum = c.SequenceNumber();
        int blsize = bl.Size(seqnum);
        for (int ibl=0; ibl<blsize; ibl++) {
            bl(seqnum,0,ibl,0) = elsol(count++,0);
        }
    }
}

/** @brief adds the connect indexes associated with base shape functions to the set */
void TPZCondensedCompEl::BuildCornerConnectList(std::set<long> &connectindexes) const
{
    int nc = NConnects();
    std::set<long> refconn;
    fReferenceCompEl->BuildCornerConnectList(refconn);
    for (int ic=0; ic<nc ; ic++) {
        TPZConnect &c = Connect(ic);
        if (!c.IsCondensed() || !c.HasDependency()) {
            long index = ConnectIndex(ic);
            if (refconn.find(index) != refconn.end()) {
                connectindexes.insert(index);
            }
        }
    }
}

int TPZCondensedCompEl::ClassId(){
    //CLASSIDFRANreturn TPZCompEl::ClassId()^Hash("TPZCondensedCompEl");
}
