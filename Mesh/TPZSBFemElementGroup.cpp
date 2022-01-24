//
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#include "pzelementgroup.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzcmesh.h"
#include "tpzintpoints.h"
#include "pzintel.h"
#include "pzelctemp.h"
#include "TPZMaterial.h"

#include <numeric>

int TPZSBFemElementGroup::gDefaultPolynomialOrder = 0;

bool TPZSBFemElementGroup::gPolynomialShapeFunctions = false;

#ifdef USING_MKL
#include <mkl.h>
#elif MACOSX
#include <Accelerate/Accelerate.h>
#endif

#ifdef COMPUTE_CRC
#ifdef USING_BOOST
#include "boost/crc.hpp"
extern TPZVec<boost::crc_32_type::value_type> matglobcrc, eigveccrc, stiffcrc, matEcrc, matEInvcrc, matPhicrc;
#endif
#endif

#ifdef USING_BLAZE
//#define BLAZE_DEFAULT_ALIGMENT_FLAG = blaze::unaligned
//#define BLAZE_DEFAULT_PADDING_FLAG = blaze::unpadded
#define BLAZE_USE_VECTORIZATION 0
#define BLAZE_USE_SHARED_MEMORY_PARALLELIZATION 0
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/config/Thresholds.h>
#include <blaze/Math.h>
using blaze::columnMajor;
using blaze::DynamicMatrix;
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemelementgroup"));
static LoggerPtr loggercoefmatrices(Logger::getLogger("pz.mesh.sbfemcoefmatrices"));
static LoggerPtr loggercoefmatriceslocal(Logger::getLogger("pz.mesh.sbfemcoefmatriceslocal"));
static LoggerPtr loggerMT(Logger::getLogger("pz.mesh.sbfemelementgroupMT"));
static LoggerPtr loggersbfemcoef(Logger::getLogger("pz.mesh.sbfemelementgroupcoef"));
static LoggerPtr loggerBF(Logger::getLogger("pz.mesh.sbfemelementgroupBF"));
static LoggerPtr loggerbubble(Logger::getLogger("pz.mesh.sbfembubbleparam"));
static LoggerPtr loggerstiffnessbubble(Logger::getLogger("pz.mesh.sbfemstiffnessbubble"));
static LoggerPtr loggersbfemstifnessdata(Logger::getLogger("pz.mesh.sbfemstifnessdata"));
#endif

TPZSBFemElementGroup::TPZSBFemElementGroup(TPZCompMesh &mesh) : TPZElementGroup(mesh)
{
    fInternalPolynomialOrder = TPZSBFemElementGroup::gDefaultPolynomialOrder;
    fPolynomialShapeFunctions = TPZSBFemElementGroup::gPolynomialShapeFunctions;
    if (fInternalPolynomialOrder != 0) {
        int nshape = 0;
        int nvar = 1;
        int64_t newindex = Mesh()->AllocateNewConnect(nshape, nvar, fInternalPolynomialOrder);
        fInternalConnectIndex = newindex;
        Mesh()->ConnectVec()[fInternalConnectIndex].IncrementElConnected();
    }
}

/**
 * @brief Computes the SBFEM coefficient matrices and the mass matrix
 * @param E0 E0 coefficient matrix
 * @param E1 E1 coefficient matrix
 * @param E2 E2 coefficient matrix
 * @param M0 Mass matrix
 */
void TPZSBFemElementGroup::ComputeMatrices(TPZElementMatrixT<STATE> &E0, TPZElementMatrixT<STATE> &E1, TPZElementMatrixT<STATE> &E2, TPZElementMatrixT<STATE> &M0)
{
    std::map<int64_t,int64_t> locindex;
    int64_t ncon = fConnectIndexes.size();
    for (int64_t ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
    TPZElementMatrixT<STATE> ef(Mesh(),TPZElementMatrixT<STATE>::EF);
    int64_t nel = fElGroup.size();
    InitializeElementMatrix(E0, ef);
    InitializeElementMatrix(E1, ef);
    InitializeElementMatrix(E2, ef);
    InitializeElementMatrix(M0, ef);

    if(TPZSBFemElementGroup::gDefaultPolynomialOrder != 0){
        int ndof=0;
        for (int64_t ic=0; ic<ncon-1; ic++) {
            ndof += Mesh()->ConnectVec()[fConnectIndexes[ic]].NShape()* Mesh()->ConnectVec()[fConnectIndexes[ic]].NState();
        }
        E0.fMat.Resize(ndof, ndof);
        E1.fMat.Resize(ndof, ndof);
        E2.fMat.Resize(ndof, ndof);
    }

    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        TPZElementMatrixT<STATE> E0Loc(Mesh(),TPZElementMatrixT<STATE>::EK);
        TPZElementMatrixT<STATE> E1Loc(Mesh(),TPZElementMatrixT<STATE>::EK);
        TPZElementMatrixT<STATE> E2Loc(Mesh(),TPZElementMatrixT<STATE>::EK);
        TPZElementMatrixT<STATE> M0Loc(Mesh(),TPZElementMatrixT<STATE>::EK);
        sbfem->ComputeKMatrices(E0Loc, E1Loc, E2Loc,M0Loc);
        
#ifdef LOG4CXX
        if (loggercoefmatriceslocal->isDebugEnabled()) {
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
            E0Loc.fMat.Print("Matriz elementar E0",sout);
            E1Loc.fMat.Print("Matriz elementar E1",sout);
            E2Loc.fMat.Print("Matriz elementar E2",sout);
            M0Loc.fMat.Print("Matriz elementar M0",sout);
            LOGPZ_DEBUG(loggercoefmatriceslocal, sout.str())
        }
#endif
        int nelcon = E0Loc.NConnects();
        for (int ic=0; ic<nelcon; ic++) {
            int iblsize = E0Loc.fBlock.Size(ic);
            int icindex = E0Loc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int jc = 0; jc<nelcon; jc++) {
                int jblsize = E0Loc.fBlock.Size(jc);
                int jcindex = E0Loc.fConnect[jc];
                int jbldest = locindex[jcindex];
                for (int idf = 0; idf<iblsize; idf++) {
                    for (int jdf=0; jdf<jblsize; jdf++) {
                        E0.at(ibldest,jbldest,idf,jdf) += E0Loc.at(ic,jc,idf,jdf);
                        E1.at(ibldest,jbldest,idf,jdf) += E1Loc.at(ic,jc,idf,jdf);
                        E2.at(ibldest,jbldest,idf,jdf) += E2Loc.at(ic,jc,idf,jdf);
                        M0.at(ibldest,jbldest,idf,jdf) += M0Loc.at(ic,jc,idf,jdf);
                    }
                }
            }
        }
    }
}

/**
 * @brief Compute the eigenvalues and eigenvectors using blaze-lib
 */
void TPZSBFemElementGroup::CalcStiffBlaze(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef)
{
#ifdef USING_BLAZE
    
    TPZElementMatrixT<STATE> E0, E1, E2, M0;
    ComputeMatrices(E0, E1, E2, M0);

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "BLAZE VERSION\n";
        E0.fMat.Print("E0 = ",sout, EMathematicaInput);
        E1.fMat.Print("E1 = ",sout, EMathematicaInput);
        E2.fMat.Print("E2 = ",sout, EMathematicaInput);
        M0.fMat.Print("M0 = ",sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    int n = E0.fMat.Rows();
    
    int dim = Mesh()->Dimension();
    
    blaze::DynamicMatrix<STATE,blaze::columnMajor> E0blaze(n,n), E1blaze(n,n), E2blaze(n,n), E0Invblaze(n,n);

    memcpy(&E0blaze.data()[0], E0.fMat.Adress(), E0blaze.spacing()*n*sizeof(STATE));
    memcpy(&E1blaze.data()[0], E1.fMat.Adress(), E0blaze.spacing()*n*sizeof(STATE));
    memcpy(&E2blaze.data()[0], E2.fMat.Adress(), E0blaze.spacing()*n*sizeof(STATE));

#ifdef COMPUTE_CRC
    {
        boost::crc_32_type crc;
        if(E0blaze.spacing() != n) DebugStop();
        crc.process_bytes(&E0blaze(0,0), E0blaze.spacing()*n*sizeof(STATE));
        crc.process_bytes(&E1blaze(0,0), E1blaze.spacing()*n*sizeof(STATE));
        crc.process_bytes(&E2blaze(0,0), E2blaze.spacing()*n*sizeof(STATE));
        matEcrc[Index()] = crc.checksum();
    }
#endif
    
    E0Invblaze = inv( E0blaze );  // Compute the inverse of E0

#ifdef COMPUTE_CRC
    {
        boost::crc_32_type crc;
        if(E0Invblaze.spacing() != n) DebugStop();
        crc.process_bytes(&E0Invblaze(0,0), E0Invblaze.spacing()*n*sizeof(STATE));
        matEInvcrc[Index()] = crc.checksum();
    }
#endif

    blaze::DynamicMatrix<STATE,blaze::columnMajor> globmatblaze(2*n,2*n);
    blaze::DynamicMatrix<STATE,blaze::columnMajor> E0InvE1Tblaze = E0Invblaze*trans(E1blaze);
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
        	globmatblaze(i, j) = E0InvE1Tblaze(i,j);
            globmatblaze(i, j+n) = -E0Invblaze(i,j);
        }
    }

    blaze::DynamicMatrix<STATE,blaze::columnMajor> E1E0InvE1T = E1blaze*E0InvE1Tblaze;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmatblaze(i+n,j) = E1E0InvE1T(i,j)-E2blaze(i,j);
        }
    }

    blaze::DynamicMatrix<STATE,blaze::columnMajor> E1E0Invblaze = E1blaze*E0Invblaze;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++){
            globmatblaze(i+n,j+n) = -E1E0Invblaze(i,j);
        }
    }
    for(int i=0; i<n; i++){
    	globmatblaze(i,i) -= (dim-2)*0.5;
        globmatblaze(i+n,i+n) += (dim-2)*0.5;
    }
    
    blaze::DynamicVector<blaze::complex<double>,blaze::columnVector> eigvalblaze( 2*n );       // The vector for the real eigenvalues
	blaze::DynamicMatrix<blaze::complex<double>,blaze::columnMajor> eigvecblaze( 2*n, 2*n );
    
    // The matrix for the left eigenvectors
    static pthread_mutex_t mutex_serial = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex_serial);
	eigen(globmatblaze, eigvalblaze, eigvecblaze);
    pthread_mutex_unlock(&mutex_serial);

#ifdef COMPUTE_CRC2
    pthread_mutex_lock(&mutex);
    extern int gnumthreads;
    std::stringstream sout;
    sout << "eigval" << gnumthreads << ".nb";
    static int count_loc = 0;
    std::ofstream file;
    if (count_loc == 0) {
        file.open(sout.str());
    }
    else
    {
        file.open(sout.str(),std::ios::app);
    }
    std::stringstream eigv;
    eigv << "EigVec" << Index() << " = ";
    if(count_loc < 1)
    {
        eigvalblaze.Print(eigv.str().c_str(),file,EMathematicaInput);
    }
    count_loc++;
    pthread_mutex_unlock(&mutex);
#endif
    
    fQVectors.Resize(n, n);
    fPhi.Resize(n, n);
    TPZManVector<std::complex<double> > eigvalsel(n,0);
    TPZFMatrix<std::complex<double> > eigvecsel(2*n,n,0.),eigvalmat(1,n,0.);
    int count = 0;
    for (int i=0; i<2*n; i++) {
        if (eigvalblaze[i].real() < -1.e-6) {
            double maxvaleigenvec = 0;
            for (int j=0; j<n; j++) {
                fQVectors(j,count) = eigvecblaze(j+n,i);
                eigvecsel(j,count) = eigvecblaze(j,i);
                eigvecsel(j+n,count) = eigvecblaze(j+n,i);
                fPhi(j,count) = eigvecblaze(j,i);
                double realvalabs = fabs(fPhi(j,count).real());
                if (realvalabs > maxvaleigenvec) {
                    maxvaleigenvec = realvalabs;
                }
            }
            eigvalsel[count] = eigvalblaze[i];
            eigvalmat(0,count) = eigvalblaze[i];
            for (int j=0; j<n; j++) {
                fQVectors(j,count) /= maxvaleigenvec;
                eigvecsel(j,count) /= maxvaleigenvec;
                eigvecsel(j+n,count) /= maxvaleigenvec;
                fPhi(j,count) /= maxvaleigenvec;
            }
            count++;
        }
    }
    
#ifdef PZDEBUG
    if(loggercoefmatrices->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "eigval = {" << eigvalsel << "};\n";
        LOGPZ_DEBUG(loggercoefmatrices, sout.str())
    }
#endif

    if (dim == 2)
    {
        int nstate = Connect(0).NState();
        if (nstate != 2 && nstate != 1) {
            DebugStop();
        }
        if(count != n-nstate) {
            DebugStop();
        }
        int ncon = fConnectIndexes.size();
        int eq=0;
        std::set<int64_t> cornercon;
        BuildCornerConnectList(cornercon);
        for (int ic=0; ic<ncon; ic++) {
            int64_t conindex = ConnectIndex(ic);
            if (cornercon.find(conindex) != cornercon.end())
            {
                fPhi(eq,count) = 1;
                eigvecsel(eq,count) = 1;
                if (nstate == 2)
                {
                    fPhi(eq+1,count+1) = 1;
                    eigvecsel(eq+1,count+1) = 1;
                }
            }
            eq += Connect(ic).NShape()*Connect(ic).NState();
        }
    }
    if(dim==3 && count != n)
    {
        std::cout << __PRETTY_FUNCTION__ << __LINE__ << " count = " << count << " n = " << n << std::endl;
        for(int i=0; i< 2*n; i++) std::cout << eigvalblaze[i] << std::endl;
        DebugStop();
    }
    fEigenvalues = eigvalsel;
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "eigenvalues BLAZE " << eigvalsel << std::endl;
        fPhi.Print("Phivec =",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef COMPUTE_CRC
    {
        boost::crc_32_type crc;
        crc.process_bytes(&fPhi(0,0), fPhi.Rows()*fPhi.Cols()*sizeof(STATE));
        matPhicrc[Index()] = crc.checksum();
    }
#endif
        
    fPhiInverse.Redim(n, n);
    blaze::DynamicMatrix<std::complex<double>,blaze::columnMajor> phiblaze(n,n), PhiInverseblaze(n,n);
    memcpy(&phiblaze.data()[0], fPhi.Adress(),n*n*sizeof(std::complex<double>));
    PhiInverseblaze = inv( phiblaze );  // Compute the inverse of A
    if (PhiInverseblaze.spacing() != n) DebugStop();
    memcpy(fPhiInverse.Adress(), &PhiInverseblaze.data()[0], PhiInverseblaze.spacing()*n*sizeof(std::complex<double>));


#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        fPhiInverse.Print("fPhiInverse = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif

    if(1)
    {
        std::ofstream out("EigenProblem.nb");
        TPZFMatrix<STATE> globmatkeep(2*n,2*n,0);
        memcpy(globmatkeep.Adress(), &globmatblaze.data()[0], 2*n*2*n*sizeof(STATE));
        globmatkeep.Print("matrix = ",out,EMathematicaInput);
        eigvecsel.Print("eigvec =",out,EMathematicaInput);
        eigvalmat.Print("lambda =",out,EMathematicaInput);
        fPhi.Print("phi = ",out,EMathematicaInput);
        fPhiInverse.Print("phiinv = ",out,EMathematicaInput);
        fQVectors.Print("qvec = ",out,EMathematicaInput);
    }
    
    ComputeMassMatrix(M0);
    if (fInternalPolynomialOrder > 0)
    {
        ComputeBubbleParameters();
    }

    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->SetPhiEigVal(fPhi, fEigenvalues);
        if (fInternalPolynomialOrder > 0)
        {
            sbfem->SetCoefNonHomogeneous(fPhiBubble, fEigenvaluesBubble, fPhiInverse, fMatBubble);
        }
    }

    InitializeElementMatrix(ek, ef);

    if (fComputationMode == EOnlyMass) {
        ek.fMat = fMassMatrix;
        ek.fMat *= fMassDensity;
        ef.fMat.Zero();
        return;
    }

    if (fPolynomialShapeFunctions)
    {
        OverwritePhis(E0,E1,E2,ek,ef);
        return;
    }

    TPZFMatrix<std::complex<double> > ekloc;
    fQVectors.Multiply(fPhiInverse, ekloc);
    
    TPZFMatrix<double> ekimag(ekloc.Rows(),ekloc.Cols());
    for (int i=0; i<ekloc.Rows(); i++) {
        for (int j=0; j<ekloc.Cols(); j++) {
            ek.fMat(i,j) = ekloc(i,j).real();
            ekimag(i,j) = ekloc(i,j).imag();
        }
    }

    ComputeMassMatrix(M0);
    
    if(fComputationMode == EMass)
    {
        int nr = ek.fMat.Rows();
        for (int r=0; r<nr; r++) {
            for (int c=0; c<nr; c++) {
                ek.fMat(r,c) += fMassMatrix(r,c)/fDelt;
            }
        }
    }

    if (fInternalPolynomialOrder > 0)
    {
        // Computing the stiffness matrix related to the bubbles
        int64_t nbubbles = fEigenvaluesBubble.size();
        int n = fPhi.Rows();

        TPZFNMatrix<100, std::complex<double>> k01(n,n,0),  k02(n,n,0), k03(n,n,0);
        for (int64_t i = 0; i < n; i++)
        {
            for (int64_t j = 0; j < n; j++)
            {
                k01(i,j) = E0.fMat(i,j);
                k02(i,j) = E1.fMat(i,j);
                k03(i,j) = E2.fMat(i,j);
            }
        }

        TPZFNMatrix<100, std::complex<double>> temp1(n,nbubbles,0), temp2(n,nbubbles,0), temp3(n,nbubbles,0);

        k01.Multiply(fPhiBubble, temp1); 
        k02.Multiply(fPhiBubble, temp2);
        k03.Multiply(fPhiBubble, temp3);

        bool transpose = 1;
        fPhiBubble.Multiply(temp1, k01, transpose); // k01 = fPhi^T * E0 * fPhi
        fPhiBubble.Multiply(temp2, k02, transpose); // k02 = fPhi^T * E1 * fPhi
        fPhiBubble.Multiply(temp3, k03, transpose); // k02 = fPhi^T * E2 * fPhi
        
        int dim = Dimension();
        TPZFNMatrix<200,std::complex<double>> K0(nbubbles,nbubbles,0);
        // K0 = ( k01* (-eigval[i])*(-eigval[j]) + k02^T * (-eigval[i]) + k02 * (-eigval[j]) + k03 / (-eigval[i]-eigval[j])
        for (int i=0; i<nbubbles; i++)
        {
            for (int j=0; j<nbubbles; j++)
            {
                if( IsZero(double(dim) - 2. - fEigenvaluesBubble[i] - fEigenvaluesBubble[j]) )
                {
                    K0(i,j) += 0;
                } else
                {
                    K0(i,j) = (k01(i,j) * (-fEigenvaluesBubble[i]*-fEigenvaluesBubble[j]) + k02(j,i) * -fEigenvaluesBubble[i] +
                     k02(i,j) * -fEigenvaluesBubble[j] + k03(i,j))/(double(dim) - 2. - fEigenvaluesBubble[i] - fEigenvaluesBubble[j]);
                }
            }
        }
    
        K0.Multiply(fMatBubble, temp1);
        
        int ndofbubbles = fMatBubble.Cols();
        TPZFNMatrix<100, std::complex<double>> K(ndofbubbles,ndofbubbles,0);
        fMatBubble.Multiply(temp1, K, transpose);  // K = fMatBubble^T K fMatBubble;

        ek.fMat.Resize(n+ndofbubbles,n+ndofbubbles);
        for (int i=0; i<ndofbubbles; i++)
        {
            for (int j=0; j<ndofbubbles; j++)
            {
                ek.fMat(i+n,j+n) = K(i,j).real();
            }
        }

        if (fComputationMode == EStiffBubble)
        {
            TPZFNMatrix<100, std::complex<double>> k01(n,n,0),  k02(n,n,0), k03(n,n,0), k02t(n,n,0);
            for (int64_t i = 0; i < n; i++)
            {
                for (int64_t j = 0; j < n; j++)
                {
                    k01(i,j) = E0.fMat(i,j);
                    k02(i,j) = E1.fMat(i,j);
                    k02t(i,j) = E1.fMat(j,i);
                    k03(i,j) = E2.fMat(i,j);
                }
            }

            TPZFNMatrix<100, std::complex<double>> temp1(n,nbubbles,0), temp2(n,nbubbles,0), temp2t(n,nbubbles,0), temp3(n,nbubbles,0);
            k01.Multiply(fPhiBubble, temp1); 
            k02.Multiply(fPhiBubble, temp2);
            k02t.Multiply(fPhiBubble, temp2t);
            k03.Multiply(fPhiBubble, temp3);

            bool transpose = 1;
            fPhi.Multiply(temp1, k01, transpose); // k01 = fPhi^T * E0 * fPhi
            fPhi.Multiply(temp2, k02, transpose); // k02 = fPhi^T * E1 * fPhi
            fPhi.Multiply(temp2t, k02t, transpose); // k02 = fPhi^T * E1 * fPhi
            fPhi.Multiply(temp3, k03, transpose); // k03 = fPhi^T * E2 * fPhi
            K0.Resize(n,nbubbles);
            K0.Zero();
            for (int i=0; i<n; i++)
            {
                for (int j=0; j<nbubbles; j++)
                {
                    if( IsZero(double(dim) - 2. - fEigenvalues[i] - fEigenvaluesBubble[j]) )
                    {
                        K0(i,j) = 0;
                    } else
                    {
                        K0(i,j) = (k01(i,j) * (-fEigenvalues[i]*-fEigenvaluesBubble[j]) + k02t(i,j) * -fEigenvalues[i] +
                        k02(i,j) * -fEigenvaluesBubble[j] + k03(i,j))/(double(dim) - 2. - fEigenvalues[i] - fEigenvaluesBubble[j]);
                    }
                }
            }
    
            K0.Multiply(fMatBubble, temp1);
            
            int ndofbubbles = fMatBubble.Cols();
            fPhiInverse.Multiply(temp1, K, transpose);  // K = fMatBubble^T K fMatBubble;

            for (int i=0; i<n; i++)
            {
                for (int j=0; j<ndofbubbles; j++)
                {
                    ek.fMat(i,j+n) = K(i,j).real();
                    ek.fMat(j+n,i) = K(i,j).real();
                }
            }
        }

        // Computing the force vector
        int icon = this->ConnectIndex(NConnects()-1);

        TPZFNMatrix<100,std::complex<double>> f(n,1,0);
        TPZFNMatrix<100,std::complex<double>> fbubble(nbubbles,1,0);

        int64_t nel = fElGroup.size();
        
        for (int64_t j = 0; j<nel; j++) {
            TPZCompEl *cel = fElGroup[j];
            TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
            if (!sbfem) DebugStop();
#endif
            sbfem->LocalBodyForces(f, fbubble, fEigenvalues, fEigenvaluesBubble, icon);
        }
        ef.fMat.Zero();

        TPZFNMatrix<200, std::complex<REAL>> ef0;
        TPZFNMatrix<200, std::complex<REAL>> efbubbles;

        f.Transpose();
        fbubble.Transpose();
        f.Multiply(fPhiInverse, ef0);
        fbubble.Multiply(fMatBubble, efbubbles);
        
        ef.fMat.Resize(n+ndofbubbles,1);
        for (int i=0; i<n; i++) {
            ef.fMat(i,0) = ef0(0,i).real();
        }
        for (int i=0; i<ndofbubbles; i++) {
            ef.fMat(i+n,0) = efbubbles(0,i).real();
        }
    }
#endif
}

void TPZSBFemElementGroup::CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef)
{
#ifdef USING_BLAZE
    CalcStiffBlaze(ek,ef);
    return;
#endif

    TPZElementMatrixT<STATE> E0, E1, E2, M0;
    ComputeMatrices(E0, E1, E2, M0);

#ifdef LOG4CXX
    if (loggercoefmatrices->isDebugEnabled()) {
        std::stringstream sout;
        E0.fMat.Print("E0 = ",sout, EMathematicaInput);
        E1.fMat.Print("E1 = ",sout, EMathematicaInput);
        E2.fMat.Print("E2 = ",sout, EMathematicaInput);
        M0.fMat.Print("M0 = ",sout, EMathematicaInput);
        LOGPZ_DEBUG(loggercoefmatrices, sout.str())
    }
#endif
    
    int n = E0.fMat.Rows();
    int dim = Mesh()->Dimension();
    
    TPZFMatrix<STATE> E0Inv(E0.fMat);
    if(1)
    {
        try
        {
            TPZFMatrix<STATE> E0copy(E0.fMat);
            TPZVec<int> pivot;
            E0copy.Decompose_LU(pivot);
            E0Inv.Identity();
            E0copy.Substitution(&E0Inv, pivot);
        }
        catch(...)
        {
            exit(-1);
        }
    }
    else
    {
        TPZVec<int> pivot(E0Inv.Rows(),0);
        int nwork = 4*n*n + 2*n;
        TPZVec<STATE> work(2*nwork,0.);
        int info=0;
#ifdef STATEdouble
        dgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
#endif
#ifdef STATEfloat
        sgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
#endif
        if (info != 0) {
            DebugStop();
        }
#ifdef STATEdouble
        dgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
#endif
#ifdef STATEfloat
        sgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
#endif
        if (info != 0) {
            DebugStop();
        }
    }
    
    TPZFMatrix<STATE> globmat(2*n,2*n,0.);
    
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i,j+n) = -E0Inv(i,j);
        }
    }
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1.fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1.fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i+n,j) -= E2.fMat(i,j);
        }
    }

#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1.fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1.fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif

    for (int i=0; i<n; i++) {
        globmat(i,i) -= (dim-2)*0.5;
        globmat(i+n,i+n) += (dim-2)*0.5;
    }
    
    TPZFMatrix<STATE> globmatkeep(globmat);
    TPZFNMatrix<100,std::complex<double> > eigenVectors;
    TPZManVector<std::complex<double> > eigenvalues;
    SolveEigenProblemSBFEM(globmatkeep, eigenvalues, eigenVectors);
 
    {
#ifdef COMPUTE_CRC
        static pthread_mutex_t mutex =PTHREAD_MUTEX_INITIALIZER;
        pthread_mutex_lock(&mutex);
        extern int gnumthreads;
        std::stringstream sout;
        sout << "eigval" << gnumthreads << ".nb";
        static int count = 0;
        std::ofstream file;
        if (count == 0) {
            file.open(sout.str());
        }
        else
        {
            file.open(sout.str(),std::ios::app);
        }
        std::stringstream eigv;
        eigv << "EigVec" << Index() << " = ";
        if(count < 1)
        {
            eigenVectors.Print(eigv.str().c_str(),file,EMathematicaInput);
        }
        count++;
        pthread_mutex_unlock(&mutex);
#endif
    }

    if(0)
    {
        TPZManVector<STATE> eigvalreal(2*n,0.);
        TPZFMatrix<STATE> eigvecreal(2*n,2*n,0.);
        for (int i=0; i<2*n; i++) {
            eigvalreal[i] = eigenvalues[i].real();
            for (int j=0; j<2*n; j++) {
                eigvecreal(i,j) = eigenVectors(i,j).real();
            }
        }
    }
        
    fQVectors.Resize(n, n);
    fPhi.Resize(n, n);
    TPZManVector<std::complex<double> > eigvalsel(n,0);
    TPZFMatrix<std::complex<double> > eigvecsel(2*n,n,0.), eigvalmat(1,n,0.);
    int count = 0;
    for (int i=0; i<2*n; i++) {
        if (eigenvalues[i].real() < -1.e-6) {
            double maxvaleigenvec = 0;
            for (int j=0; j<n; j++) {
                fQVectors(j,count) = eigenVectors(j+n,i);
                eigvecsel(j,count) = eigenVectors(j,i);
                eigvecsel(j+n,count) = eigenVectors(j+n,i);
                fPhi(j,count) = eigenVectors(j,i);
                double realvalabs = fabs(fPhi(j,count));
                if (realvalabs > maxvaleigenvec) {
                    maxvaleigenvec = realvalabs;
                }
            }
            eigvalsel[count] = eigenvalues[i];
            eigvalmat(0,count) = eigenvalues[i];
            for (int j=0; j<n; j++) {
                fQVectors(j,count) /= maxvaleigenvec;
                eigvecsel(j,count) /= maxvaleigenvec;
                eigvecsel(j+n,count) /= maxvaleigenvec;
                fPhi(j,count) /= maxvaleigenvec;
            }
            count++;
        }
    }
    
#ifdef LOG4CXX
    if(loggercoefmatrices->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "eigval = {" << eigvalsel << "};\n";
        LOGPZ_DEBUG(loggercoefmatrices, sout.str())
    }
#endif

    if (dim == 2)
    {
        int nstate = Connect(0).NState();
        if (nstate != 2 && nstate != 1) {
            DebugStop();
        }
        if(count != n-nstate) {
            DebugStop();
        }
        int ncon = fConnectIndexes.size();
        int eq=0;
        std::set<int64_t> cornercon;
        BuildCornerConnectList(cornercon);
        for (int ic=0; ic<ncon; ic++) {
            int64_t conindex = ConnectIndex(ic);
            if (cornercon.find(conindex) != cornercon.end())
            {
                fPhi(eq,count) = 1;
                eigvecsel(eq,count) = 1;
                if (nstate == 2)
                {
                    fPhi(eq+1,count+1) = 1;
                    eigvecsel(eq+1,count+1) = 1;
                }
            }
            eq += Connect(ic).NShape()*Connect(ic).NState();
        }
    }
    if(dim==3 && count != n)
    {
        DebugStop();
    }

    fEigenvalues = eigvalsel;

#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "eigenvalues " << eigvalsel << std::endl;
        fPhi.Print("Phivec =",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
   double phinorm = Norm(fPhi);
   if (loggerMT->isDebugEnabled()) {
       std::stringstream sout;
       sout << "Element index " << Index() << " phinorm = " << phinorm;
       LOGPZ_DEBUG(loggerMT, sout.str())
   }
#endif
        
    TPZFMatrix<std::complex<double> > phicopy(fPhi);
    fPhiInverse.Redim(n, n);
    fPhiInverse.Identity();
    
    try
    {
        TPZVec<int> pivot;
        phicopy.Decompose_LU(pivot);
        phicopy.Substitution(&fPhiInverse, pivot);
    }
    catch(...)
    {
        exit(-1);
    }

#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Index = " << Index() << "\n";
        sout << "eigval = {" << EigenvaluesReal() << "};\n";
        // globmat.Print("Zp = ", sout,EMathematicaInput);
        fPhi.Print("fPhi = ", sout,EMathematicaInput);
        fPhiInverse.Print("fPhiInverse = ",sout,EMathematicaInput);
        fQVectors.Print("QVectors = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    if (fInternalPolynomialOrder > 0)
    {
        ComputeBubbleParameters();
    }

    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->SetPhiEigVal(fPhi, fEigenvalues);
        if (fInternalPolynomialOrder > 0)
        {
            sbfem->SetCoefNonHomogeneous(fPhiBubble, fEigenvaluesBubble, fPhiInverse, fMatBubble);
        }
    }

    InitializeElementMatrix(ek, ef);

    if (fComputationMode == EOnlyMass) {
        ek.fMat = fMassMatrix;
        ek.fMat *= fMassDensity;
        ef.fMat.Zero();
        return;
    }

    if (fPolynomialShapeFunctions)
    {
        OverwritePhis(E0,E1,E2,ek,ef);
        return;
    }

    TPZFMatrix<std::complex<double> > ekloc;
    fQVectors.Multiply(fPhiInverse, ekloc);
    
    TPZFMatrix<double> ekimag(ekloc.Rows(),ekloc.Cols());
    for (int i=0; i<ekloc.Rows(); i++) {
        for (int j=0; j<ekloc.Cols(); j++) {
            ek.fMat(i,j) = ekloc(i,j).real();
            ekimag(i,j) = ekloc(i,j).imag();
        }
    }

    ComputeMassMatrix(M0);
    
    if(fComputationMode == EMass)
    {
        int nr = ek.fMat.Rows();
        for (int r=0; r<nr; r++) {
            for (int c=0; c<nr; c++) {
                ek.fMat(r,c) += fMassMatrix(r,c)/fDelt;
            }
        }
    }

    if (fInternalPolynomialOrder > 0)
    {
        // Computing the stiffness matrix related to the bubbles
        int64_t nbubbles = fEigenvaluesBubble.size();
        int n = fPhi.Rows();

        TPZFNMatrix<100, std::complex<double>> k01(n,n,0),  k02(n,n,0), k03(n,n,0);
        for (int64_t i = 0; i < n; i++)
        {
            for (int64_t j = 0; j < n; j++)
            {
                k01(i,j) = E0.fMat(i,j);
                k02(i,j) = E1.fMat(i,j);
                k03(i,j) = E2.fMat(i,j);
            }
        }

        TPZFNMatrix<100, std::complex<double>> temp1(n,nbubbles,0), temp2(n,nbubbles,0), temp3(n,nbubbles,0);

        k01.Multiply(fPhiBubble, temp1); 
        k02.Multiply(fPhiBubble, temp2);
        k03.Multiply(fPhiBubble, temp3);

        bool transpose = 1;
        fPhiBubble.Multiply(temp1, k01, transpose); // k01 = fPhi^T * E0 * fPhi
        fPhiBubble.Multiply(temp2, k02, transpose); // k02 = fPhi^T * E1 * fPhi
        fPhiBubble.Multiply(temp3, k03, transpose); // k02 = fPhi^T * E2 * fPhi
        
        int dim = Dimension();
        TPZFNMatrix<200,std::complex<double>> K0(nbubbles,nbubbles,0);
        // K0 = ( k01* (-eigval[i])*(-eigval[j]) + k02^T * (-eigval[i]) + k02 * (-eigval[j]) + k03 / (-eigval[i]-eigval[j])
        for (int i=0; i<nbubbles; i++)
        {
            for (int j=0; j<nbubbles; j++)
            {
                if( IsZero(double(dim) - 2. - fEigenvaluesBubble[i] - fEigenvaluesBubble[j]) )
                {
                    K0(i,j) += 0;
                } else
                {
                    K0(i,j) = (k01(i,j) * (-fEigenvaluesBubble[i]*-fEigenvaluesBubble[j]) + k02(j,i) * -fEigenvaluesBubble[i] +
                     k02(i,j) * -fEigenvaluesBubble[j] + k03(i,j))/(double(dim) - 2. - fEigenvaluesBubble[i] - fEigenvaluesBubble[j]);
                }
            }
        }
    
        K0.Multiply(fMatBubble, temp1);
        
        int ndofbubbles = fMatBubble.Cols();
        TPZFNMatrix<100, std::complex<double>> K(ndofbubbles,ndofbubbles,0);
        fMatBubble.Multiply(temp1, K, transpose);  // K = fMatBubble^T K fMatBubble;

        ek.fMat.Resize(n+ndofbubbles,n+ndofbubbles);
        for (int i=0; i<ndofbubbles; i++)
        {
            for (int j=0; j<ndofbubbles; j++)
            {
                ek.fMat(i+n,j+n) = K(i,j).real();
            }
        }

        if (fComputationMode == EStiffBubble)
        {
            TPZFNMatrix<100, std::complex<double>> k01(n,n,0),  k02(n,n,0), k03(n,n,0), k02t(n,n,0);
            for (int64_t i = 0; i < n; i++)
            {
                for (int64_t j = 0; j < n; j++)
                {
                    k01(i,j) = E0.fMat(i,j);
                    k02(i,j) = E1.fMat(i,j);
                    k02t(i,j) = E1.fMat(j,i);
                    k03(i,j) = E2.fMat(i,j);
                }
            }

            TPZFNMatrix<100, std::complex<double>> temp1(n,nbubbles,0), temp2(n,nbubbles,0), temp2t(n,nbubbles,0), temp3(n,nbubbles,0);
            k01.Multiply(fPhiBubble, temp1); 
            k02.Multiply(fPhiBubble, temp2);
            k02t.Multiply(fPhiBubble, temp2t);
            k03.Multiply(fPhiBubble, temp3);

            bool transpose = 1;
            fPhi.Multiply(temp1, k01, transpose); // k01 = fPhi^T * E0 * fPhi
            fPhi.Multiply(temp2, k02, transpose); // k02 = fPhi^T * E1 * fPhi
            fPhi.Multiply(temp2t, k02t, transpose); // k02 = fPhi^T * E1 * fPhi
            fPhi.Multiply(temp3, k03, transpose); // k03 = fPhi^T * E2 * fPhi
            K0.Resize(n,nbubbles);
            K0.Zero();
            for (int i=0; i<n; i++)
            {
                for (int j=0; j<nbubbles; j++)
                {
                    if( IsZero(double(dim) - 2. - fEigenvalues[i] - fEigenvaluesBubble[j]) )
                    {
                        K0(i,j) = 0;
                    } else
                    {
                        K0(i,j) = (k01(i,j) * (-fEigenvalues[i]*-fEigenvaluesBubble[j]) + k02t(i,j) * -fEigenvalues[i] +
                        k02(i,j) * -fEigenvaluesBubble[j] + k03(i,j))/(double(dim) - 2. - fEigenvalues[i] - fEigenvaluesBubble[j]);
                    }
                }
            }
    
            K0.Multiply(fMatBubble, temp1);
            
            int ndofbubbles = fMatBubble.Cols();
            fPhiInverse.Multiply(temp1, K, transpose);  // K = fMatBubble^T K fMatBubble;

            for (int i=0; i<n; i++)
            {
                for (int j=0; j<ndofbubbles; j++)
                {
                    ek.fMat(i,j+n) = K(i,j).real();
                    ek.fMat(j+n,i) = K(i,j).real();
                }
            }
        }
        
        
#ifdef LOG4CXX
        if(loggerstiffnessbubble->isDebugEnabled())
        {
            std::stringstream sout;
            K.Print("KBubble = ",sout,EMathematicaInput);
            LOGPZ_DEBUG(loggerstiffnessbubble, sout.str())
        }
#endif

        // Computing the force vector
        int icon = this->ConnectIndex(NConnects()-1);

        TPZFNMatrix<100,std::complex<double>> f(n,1,0);
        TPZFNMatrix<100,std::complex<double>> fbubble(nbubbles,1,0);

        int64_t nel = fElGroup.size();
        
        for (int64_t j = 0; j<nel; j++) {
            TPZCompEl *cel = fElGroup[j];
            TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
            if (!sbfem) DebugStop();
#endif
            sbfem->LocalBodyForces(f, fbubble, fEigenvalues, fEigenvaluesBubble, icon);
        }
        ef.fMat.Zero();

        TPZFNMatrix<200, std::complex<REAL>> ef0;
        TPZFNMatrix<200, std::complex<REAL>> efbubbles;

        f.Transpose();
        fbubble.Transpose();
        f.Multiply(fPhiInverse, ef0);
        fbubble.Multiply(fMatBubble, efbubbles);
        
        ef.fMat.Resize(n+ndofbubbles,1);
        for (int i=0; i<n; i++) {
            ef.fMat(i,0) = ef0(0,i).real();
        }
        for (int i=0; i<ndofbubbles; i++) {
            ef.fMat(i+n,0) = efbubbles(0,i).real();
        }
        
#ifdef LOG4CXX
        if (loggerBF->isDebugEnabled()) {
            std::stringstream sout;

            sout << Index() << "\n";

            K0.Print("K0 = ", sout, EMathematicaInput);
            f.Print("f = ", sout, EMathematicaInput);
            fbubble.Print("fbubble = ", sout, EMathematicaInput);
            
            ek.fMat.Print("ek = ", sout, EMathematicaInput);
            ef.fMat.Print("ef = ", sout, EMathematicaInput);
            LOGPZ_DEBUG(loggerBF, sout.str())
        }
#endif
    }
}

void TPZSBFemElementGroup::LoadSolution()
{
    
    int nc = NConnects();
    int ndofs = fPhiInverse.Cols()+fMatBubble.Cols();
    int ncoef = fPhiInverse.Rows()+fMatBubble.Rows();
    
    TPZFNMatrix<100, std::complex<double> > uh_local(ncoef, fMesh->Solution().Cols(),0.);
    fCoef.Resize(ncoef,fMesh->Solution().Cols());

    TPZFMatrix<STATE> &meshSol = fMesh->Solution();
    
    int count = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = Connect(ic);
        int nshape = c.NShape();
        int nstate = c.NState();
        int blsize = nshape*nstate;
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = fMesh->Block().Position(seqnum);
        for (int seq=0; seq < blsize; seq++) {
            for (int j = 0; j < uh_local.Cols(); j++)
            {
                uh_local(count+seq,j) = meshSol(pos+seq,j);
            }
        }
        count += blsize;
    }
    if(fInternalPolynomialOrder == 0)
    {
        fPhiInverse.Multiply(uh_local, fCoef);
    } else
    {
        TPZFMatrix<std::complex<REAL> > coefbubble;
        TPZFMatrix<std::complex<REAL> > uh_localbubble(fMatBubble.Cols(),1);
        for (int64_t i = 0; i < fMatBubble.Cols(); i++)
        {
            uh_localbubble(i,0) = uh_local(i+fPhiInverse.Cols());
        }
        fMatBubble.Multiply(uh_localbubble, coefbubble);
        
        uh_local.Resize(fPhiInverse.Rows(),1);
        fPhiInverse.Multiply(uh_local, fCoef);

        fCoef.Resize(ncoef,1);
        for (int64_t i = 0; i < fMatBubble.Rows(); i++)
        {
            fCoef(i+fPhiInverse.Rows(),0) = coefbubble(i);
        }
    }

    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->LoadCoef(fCoef);
    }
    #ifdef LOG4CXX
        if (loggersbfemcoef->isDebugEnabled()) {
            std::stringstream sout;
            fCoef.Print("fCoef = ", sout, EMathematicaInput);
            LOGPZ_DEBUG(loggersbfemcoef, sout.str())
        }
#endif

}

/// Compute the mass matrix based on the value of M0 and the eigenvectors
void TPZSBFemElementGroup::ComputeMassMatrix(TPZElementMatrixT<STATE> &M0)
{
    //    M0.fMat.Print("Mass = ",std::cout,EMathematicaInput);
    TPZFMatrix<std::complex<REAL> > temp;
    REAL alpha = 1.;
    REAL beta = 0.;
    int transpose = 1;
    int nrow = fEigenvalues.size();
    TPZFNMatrix<100,std::complex<REAL> > M0Loc(nrow,nrow),MassLoc(nrow,nrow);
    for (int i=0; i<nrow; i++) {
        for(int j=0; j<nrow; j++)
        {
            M0Loc(i, j) = M0.fMat(i,j);
        }
    }
    fPhi.MultAdd(M0Loc, M0Loc, temp, alpha, beta, transpose);
    //    temp.Print("Temp = ",std::cout,EMathematicaInput);
    temp.Multiply(fPhi,MassLoc);
    for (int i=0; i<nrow; i++) {
        for (int j=0; j<nrow; j++) {
            MassLoc(i,j) /= (2.-fEigenvalues[i]-fEigenvalues[j]);
        }
    }
    fPhiInverse.MultAdd(MassLoc, M0Loc, temp,alpha, beta, transpose);
    temp.Multiply(fPhiInverse,MassLoc);
    
    fMassMatrix.Redim(nrow, nrow);
    for (int i=0; i<nrow; i++) {
        for(int j=0; j<nrow; j++)
        {
            fMassMatrix(i, j) = MassLoc(i,j).real();
        }
    }
}

void TPZSBFemElementGroup::LoadEigenVector(int64_t eig)
{
    fCoef.Zero();
    fCoef(eig,0) = 1.;
    if(eig>7 && eig < 11)
        fCoef(eig+1,0) = -1;
    
    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->LoadCoef(fCoef);
    }
    
}

/** @brief add an element to the element group
 */
void TPZSBFemElementGroup::AddElement(TPZCompEl *cel)
{
    std::set<int64_t> connects;
    int nc = fConnectIndexes.size();
    for (int ic=0; ic<nc; ic++) {
        connects.insert(fConnectIndexes[ic]);
    }
    TPZSBFemVolume *celvol = dynamic_cast<TPZSBFemVolume *>(cel);
    TPZCompEl *celskeleton = Mesh()->Element(celvol->SkeletonIndex());
    nc = celskeleton->NConnects();
    for (int ic=0; ic<nc; ic++) {
        connects.insert(celskeleton->ConnectIndex(ic));
    }
    nc = connects.size();
    if (nc != fConnectIndexes.size()) {
        fConnectIndexes.Resize(nc, 0);
        std::set<int64_t>::iterator it = connects.begin();
        for (int ic = 0; it != connects.end(); it++,ic++) {
            fConnectIndexes[ic] = *it;
        }
    }
    TPZElementGroup::AddElement(cel);
}

void TPZSBFemElementGroup::InitializeInternalConnect()
{
    if(fElGroup.size() == 0) return;
    int64_t ncon = NConnects();
    int nshapeboundary = 0;
    int64_t nshapeinternal = 0;

    std::vector<int64_t> connects(ncon);
    for (int64_t i=0; i<ncon; i++) {
        connects[i] = fConnectIndexes[i];
        nshapeboundary += Mesh()->ConnectVec()[fConnectIndexes[i]].NShape()*Mesh()->ConnectVec()[fConnectIndexes[i]].NState();
    }
    if(!(fElGroup[0]->Material()) ){
        std::cout << "TPZSBFemElementGroup::InitializeInternalConnect - No material associated\n";
        DebugStop();
    }
    std::vector<int64_t>::iterator it = std::find(connects.begin(), connects.end(), fInternalConnectIndex);
    if (it != connects.end()) {
        return;
    }

    // Adding the polynomial connects
    nshapeinternal += (fInternalPolynomialOrder-1)*nshapeboundary; // if p=1, only the hat function is added
    int nstate = fElGroup[0]->Material()->NStateVariables();
    nshapeinternal += nstate; // hat function

    // Adding the non-integer bubble functions
    // The number of non-integer exponents is the number eigenvalues - (dim*nstate*pOrder+1)
    // dim*nstate*pOrder+1 is the number of eigenvalues <= pOrder
    int integerexponents = this->Dimension()*nstate*Mesh()->GetDefaultOrder() + nstate;
    nshapeinternal += nshapeboundary - integerexponents;
    
    TPZConnect &c = Mesh()->ConnectVec()[fInternalConnectIndex];
    int64_t seq = c.SequenceNumber();

    fConnectIndexes.Resize(ncon+1);
    fConnectIndexes[ncon] = fInternalConnectIndex;
    
    c.SetNShape(nshapeinternal);
    Mesh()->Block().Set(seq, nshapeinternal);
    Mesh()->ExpandSolution();
}

void TPZSBFemElementGroup::ComputeBubbleParameters()
{
    int cont = 0;
    int nstate = fElGroup[0]->Material()->NStateVariables();
    int n = fPhi.Rows();

    // Finding the rational eigenvalues

    /*
    
    Number of exponents applied to compute the bubble functions based in the basis functions:
    - It is composed of the number of eigenvalues except for the \lambda = -1 and \lambda = 0,
    \lambda = -1 and \lambda = 0 will be the hat functions.
    - Moreover, the basic modes of displacement, characterized by the polynomial order adopted must be included only once. For
    instance, a mode for \lambda = 2 and nstate = 1 in a 2D problem will appear two times (two different eigenvectors), but one is
    the linear combination of the other one.

    */
    std::map<int,REAL> eigmap, eignull; // indexing the position of the eigenvalue
    int aproxeig = 0;

    if(fPolynomialShapeFunctions == false)
    {
        if (Dimension() == 2)
        {
            for (int i=0; i<n; i++)
            {
                // Including bubbles composed by the basis functions - except the hat function : it will be included in the polynomial space
                aproxeig = -nearbyint(fEigenvalues[i].real());
                if (aproxeig > fInternalPolynomialOrder)
                {
                    eigmap[i] = fEigenvalues[i].real();
                }
                if (IsZero(fEigenvalues[i].real()))
                {
                    eignull[i] = fEigenvalues[i].real();
                }
            }
        }
        else if(Dimension() == 3)
        {
            for (int i=0; i<n; i++)
            {
                // Including bubbles composed by the basis functions - except the hat function : it will be included in the polynomial space
                aproxeig = -nearbyint(fEigenvalues[i].real()+0.5);
                if (aproxeig > fInternalPolynomialOrder)
                {
                    eigmap[i] = fEigenvalues[i].real();
                }
                if (IsZero(fEigenvalues[i].real()+0.5))
                {
                    eignull[i] = fEigenvalues[i].real()+0.5;
                }
            }
        }
    }
    
    // Updating the connect for the bubble function
    /*
    Now I know how many bubble functions based on the basis functions I have:
    the number of different exponents of the eigenvalue problem (size of indeig) + a hat function + polynomial bubbles
    */
    // TPZConnect &c = Mesh()->ConnectVec()[fInternalConnectIndex];
    // int64_t seq = c.SequenceNumber();

    int nbubbleseig = eigmap.size();
    int neq = nbubbleseig + nstate + (fInternalPolynomialOrder-1)*n;
    
    // c.SetNShape(neq);
    // Mesh()->Block().Set(seq, neq);
    // Mesh()->ExpandSolution();

    // ##################################################################################
    // UPDATING EIGENVALUES

    // number of exponents
    
    int nexp = 2*nbubbleseig + n*(fInternalPolynomialOrder-1) + 2*nstate;
    if (fInternalPolynomialOrder > 1)
    {
        nexp += n;
    }

    fEigenvaluesBubble.resize(nexp);

    // \lambda_u & 1
    int pos = 0;
    for (auto & eigmap : eigmap)
    {
        if (Dimension() == 2)
        {
            fEigenvaluesBubble[2*pos] = fEigenvalues[eigmap.first];
        }
        else if(Dimension() == 3)
        {
            fEigenvaluesBubble[2*pos] = fEigenvalues[eigmap.first] + 0.5;
        }
        fEigenvaluesBubble[2*pos+1] = -1;//-double(fInternalPolynomialOrder);// + I*fEigenvaluesBubble[2*pos].imag();
        pos++;
    }
    
    // 0
    pos = 2*nbubbleseig;
    for (int i = 0; i < nstate; i++)
    {
        fEigenvaluesBubble[2*i + pos] = 0;
        fEigenvaluesBubble[2*i + 1 + pos] = -1;
    }
    // // Updating the eigenvalue vector for polynomial bubbles
    pos = 2*nbubbleseig + 2*nstate;

    if (fInternalPolynomialOrder > 1)
    {
        for (int p = 0; p < fInternalPolynomialOrder; p++)
        {
            for (int j = 0; j < n; j++)
            {
                fEigenvaluesBubble[j + n*p + pos] = -p-1;
            }
        }
    }
    // cout << fEigenvaluesBubble << endl;
    // ##################################################################################
    // UPDATING EIGENVECTORS
    
    fPhiBubble.Resize(n, nexp);
    fPhiBubble.Zero();

    // Rational exponents \xi^\lambda_u
    pos  = 0; // first equations are the bubbles composed of the sbfem basis functions
    for (auto & eigmap : eigmap)
    {
        for (int i = 0; i < n; i++)
        {
            fPhiBubble(i, pos*2) = fPhi(i, eigmap.first);
            fPhiBubble(i, pos*2 + 1) = fPhi(i, eigmap.first);
        }
        pos++;
    }

    // Constant function 1
    pos = 2*nbubbleseig;
    cont=0;
    for (auto & eign : eignull)
    {
        for (int i = 0; i < n; i++)
        {
            fPhiBubble(i, pos+ 2*cont) = fPhi(i, eign.first); // Eigenvector for the \xi^0
            fPhiBubble(i, pos+ 2*cont + 1) = fPhi(i, eign.first); // Eigenvector for the \xi^0
        }
        cont++;
    }

    // xi^p, p = 1, ..., k
    pos = 2*nbubbleseig + 2*nstate;
    if (fInternalPolynomialOrder > 1)
    {
        for (int p = 0; p < fInternalPolynomialOrder; p++)
        {
            for (int i = 0; i < n; i++)
            {
                fPhiBubble(i, i + n*p + pos) = 1.; // Identity matrices 
            }
        }
    }
    
    // ##################################################################################
    // UPDATING THE MATRIX THAR WILL COMPOSE THE LINEAR COMBINATIONS OF \xi^i
    
    fMatBubble.Resize(nexp, neq);
    // Line represents the number of exponents - related to the basis functions
    // Column the number of bubbles - equations
    
    // cont * rational function: \xi - \xi^\lambda
    for (int i = 0; i < nbubbleseig; i++)
    {
        fMatBubble(2*i, i) = -1;  // -\lambda
        fMatBubble(2*i+1, i) = 1; // 1
    }

    // nstate * hat function: 1 - \xi
    cont = 0;
    int posexp0 = 2*nbubbleseig; // 0 
    int ifunc = nbubbleseig;
    for (auto & eign : eignull)
    {
        fMatBubble(posexp0 + 2*cont, ifunc + cont) = 1.; // 0 
        fMatBubble(posexp0 + 2*cont + 1, ifunc + cont) = -1.; // 1
        cont++;
    }

    posexp0 = 2*nbubbleseig + 2*nstate;
    int posexp1 = 2*nbubbleseig + 2*nstate + n;
    ifunc = nbubbleseig + nstate;
    if (fInternalPolynomialOrder > 1)
    {
        for (int j = 0; j < fInternalPolynomialOrder-1; j++)
        {
            for (int i = 0; i < n; i++) 
            {
                fMatBubble(posexp0 + i + n*j, ifunc) = 1; // p-1
                fMatBubble(posexp1 + i + n*j, ifunc) = -1; // p
                ifunc++;
            }
        }
    }
    
#ifdef LOG4CXX
    if (loggerbubble->isDebugEnabled()) {
        std::stringstream sout;
        sout << Index() << "\n";
        
        fPhiBubble.Print("fPhiBubble = ", sout, EMathematicaInput);
        fMatBubble.Print("fMatBubble = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerbubble, sout.str())
    }
#endif

}

void TPZSBFemElementGroup::OverwritePhis(TPZElementMatrixT<STATE> &E0, TPZElementMatrixT<STATE> &E1, TPZElementMatrixT<STATE> &E2, TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef)
{
    auto n = fEigenvalues.size();
    fPhi.Zero();
    fPhiInverse.Zero();
    for (auto i = 0; i < n; i++)
    {
        fPhi(i,i) = 1;
        fPhiInverse(i,i) = 1;
        fEigenvalues[i] = -1;
    }

    int64_t nbubbles = fEigenvaluesBubble.size();
    int64_t ndofs = fMatBubble.Cols();
    TPZFNMatrix<100, std::complex<double>> phiu(n, n+nbubbles, 0), phiuinv(n+nbubbles, n+ndofs, 0);
    TPZManVector<std::complex<double>,10> eigval(n+nbubbles);
    for (auto i = 0; i < nbubbles; i++)
    {
        for (auto j=0; j < n; j++)
        {
            phiu(j,i+n) = fPhiBubble(j,i);
        }
        for (auto j = 0; j < ndofs; j++)
        {
            phiuinv(i+n,j+n) = fMatBubble(i,j);
        }
        eigval[i+n] = fEigenvaluesBubble[i];
    }
    for (auto i = 0; i < n; i++)
    {
        phiu(i,i) = 1;
        phiuinv(i,i) = 1;

        phiu(i,n+nbubbles-1) = 1;
        phiuinv(i,n) = 0;

        eigval[i] = -1;
    }
    TPZFNMatrix<100, std::complex<double>> k01(n,n,0),  k02(n,n,0), k03(n,n,0);
    for (int64_t i = 0; i < n; i++)
    {
        for (int64_t j = 0; j < n; j++)
        {
            k01(i,j) = E0.fMat(i,j);
            k02(i,j) = E1.fMat(i,j);
            k03(i,j) = E2.fMat(i,j);
        }
    }

    TPZFNMatrix<100, std::complex<double>> temp1(n,n+nbubbles,0), temp2(n,n+nbubbles,0), temp3(n,n+nbubbles,0);

    k01.Multiply(phiu, temp1); 
    k02.Multiply(phiu, temp2);
    k03.Multiply(phiu, temp3);

    bool transpose = 1;
    phiu.Multiply(temp1, k01, transpose); // k01 = fPhi^T * E0 * fPhi
    phiu.Multiply(temp2, k02, transpose); // k02 = fPhi^T * E1 * fPhi
    phiu.Multiply(temp3, k03, transpose); // k02 = fPhi^T * E2 * fPhi
        
    TPZFNMatrix<200,std::complex<double>> K0(n+nbubbles,n+nbubbles,0);
    // K0 = ( k01* (-eigval[i])*(-eigval[j]) + k02^T * (-eigval[i]) + k02 * (-eigval[j]) + k03 / (-eigval[i]-eigval[j])
    for (int i=0; i<n+nbubbles; i++) {
        for (int j=0; j<n+nbubbles; j++) {
            if(IsZero((-eigval[i]-eigval[j]).real())) {
                K0(i,j) += 0;
            } else {
                K0(i,j) = (k01(i,j) * (-eigval[i]*-eigval[j]) + k02(j,i) * -eigval[i] +
                    k02(i,j) * -eigval[j] + k03(i,j))/(-eigval[i]-eigval[j]);
            }
        }
    }
    
    K0.Multiply(phiuinv, temp1);
        
    TPZFNMatrix<100, std::complex<double>> K(n+ndofs,n+ndofs,0);
    phiuinv.Multiply(temp1, K, transpose);  // K = fMatBubble^T K fMatBubble;
        
    ek.fMat.Resize(n+ndofs,n+ndofs);
    for (int i=0; i<n+ndofs; i++) {
        for (int j=0; j<n+ndofs; j++) {
            ek.fMat(i,j) = K(i,j).real();
        }
    }
    int64_t nel = fElGroup.size();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->SetPhiEigVal(fPhi, fEigenvalues);
        if (fInternalPolynomialOrder > 0)
        {
            sbfem->SetCoefNonHomogeneous(fPhiBubble, fEigenvaluesBubble, fPhiInverse, fMatBubble);
        }
    }

    // Computing the force vector
    int icon = this->ConnectIndex(NConnects()-1);

    TPZFNMatrix<100,std::complex<double>> f(n,1,0);
    TPZFNMatrix<100,std::complex<double>> fbubble(nbubbles,1,0);

    for (int64_t j = 0; j<nel; j++) {
        TPZCompEl *cel = fElGroup[j];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->LocalBodyForces(f, fbubble, fEigenvalues, fEigenvaluesBubble, icon);
    }
    ef.fMat.Zero();

    TPZFNMatrix<200, std::complex<REAL>> ef0;
    TPZFNMatrix<200, std::complex<REAL>> efbubbles;

    fPhiInverse.Multiply(f, ef0, transpose);
    fMatBubble.Multiply(fbubble, efbubbles, transpose);
    
    ef.fMat.Resize(n+ndofs,1);
    ef.fMat.Zero();
    
}

void TPZSBFemElementGroup::SolveEigenProblemSBFEM(TPZFMatrix<STATE> &globmatkeep, TPZManVector<std::complex<double> > &eigenvalues, 
    TPZFNMatrix<100,std::complex<double> > &eigenvectors)
{
    if (globmatkeep.Rows() != globmatkeep.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< double > VL(globmatkeep.Rows(),globmatkeep.Cols(),0.),VR(globmatkeep.Rows(),globmatkeep.Cols(),0.);

    int dim = globmatkeep.Rows();
    int ldvr = dim;
    int lwork = 10+50*dim;
    int info = 0;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);
    
    TPZFMatrix<double> temp(globmatkeep);
    TPZVec<double> work(lwork,0.);
    dgeev_(jobvl, jobvr, &dim, &temp(0,0) , &dim, &realeigen[0], &imageigen[0], &VL(0,0), &dim, &VR(0,0), &ldvr, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if(IsZero(imageigen[i])){
            eigenvalues[i] = realeigen[i];
        } else {
            eigenvalues[i] = realeigen[i] + I*imageigen[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        if(IsZero(imageigen[i]))
        {
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i ) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }

#ifdef LOG4CXX
    for (int i = 0; i < dim; i++)
    {
        TPZFMatrix< std::complex<double> > res(dim, 1, 0.);
        TPZFMatrix< std::complex<double> > x(dim, 1, 0.);
        
        eigenvectors.GetSub(0, i, dim, 1, x);

        TPZFMatrix<std::complex<double>> cpma(dim,dim,0.);
        for(int i = 0; i < dim ; i++){
            for(int j = 0; j < dim; j++){
                cpma(i,j) = globmatkeep.GetVal(i,j);
            }
        }
        res = cpma * x - eigenvalues[i] * x;
        const auto norm = Norm(res);
        if(!IsZero(norm))
        {
            std::cout << "eigval " << eigenvalues[i]  << ", diff " << norm << std::endl;
            if (std::isnan(norm)) std::cout << "res" << res << std::endl;
        }
    }
#endif
    
}

//http://www.netlib.org/lapack/lug/node50.html
//https://software.intel.com/en-us/node/521079
