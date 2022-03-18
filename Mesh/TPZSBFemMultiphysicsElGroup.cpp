//
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/2021.
//
//
#include "TPZSBFemMultiphysicsElGroup.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSBFemVolumeMultiphysics.h"
#include "TPZMaterialDataT.h"
#include "TPZGeoLinear.h"

#ifdef USING_MKL
#include <mkl.h>
#elif MACOSX
#include <Accelerate/Accelerate.h>
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("sbfemmultiphysicselgroup"));
static LoggerPtr loggercoef(Logger::getLogger("sbfemmultiphysicscoefmatrices"));
static LoggerPtr loggereigensystem(Logger::getLogger("sbfemmultiphysicseigensystem"));
static LoggerPtr loggerlocalstiff(Logger::getLogger("loggerlocalstiff"));
static LoggerPtr loggerfulleigensys(Logger::getLogger("sbfemmultiphysicsfulleigensys"));
static LoggerPtr loggersbfemcondensed(Logger::getLogger("sbfemcondensed"));
static LoggerPtr loggerflux(Logger::getLogger("sbfemfluxhdiv"));
#endif

void TPZSBFemMultiphysicsElGroup::AddElement(TPZCompEl *cel)
{
    std::set<int64_t> connects;
    for (auto ic : fConnectIndexes)
    {
        connects.insert(ic);
    }
    
    auto celvol = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cel);
    TPZCompEl *celskeleton = Mesh()->Element(celvol->SkeletonIndex());

    auto nc = celskeleton->NConnects();
    for (int ic=0; ic<nc; ic++)
    {
        connects.insert(celskeleton->ConnectIndex(ic));
    }

    nc = connects.size();
    if (nc != fConnectIndexes.size())
    {
        fConnectIndexes.Resize(nc, 0);
        auto it = connects.begin();
        for (int ic = 0; it != connects.end(); it++,ic++)
        {
            fConnectIndexes[ic] = *it;
        }
    }
    TPZElementGroup::AddElement(cel);
}

void TPZSBFemMultiphysicsElGroup::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZElementGroup::Print(out);

    out << "Element indexes of the volume elements ";
    for (auto cel : fElGroup)
    {
        out << cel->Index() << " ";
    }
    out << std::endl;
    out << "Indices of the associated computational skeleton elements\n";
    for (auto cel : fElGroup)
    {
        auto vol = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cel);
        if(!vol) DebugStop();
        out << vol->SkeletonIndex() << " ";
    }
    out << std::endl;
    out << "Connect indexes \n";
    int nc = NConnects();
    for (int ic=0; ic<nc; ic++)
    {
        out << ConnectIndex(ic) << " ";
    }
    out << std::endl;
    out << "Connect indexes of the contained elements\n";
    for (auto cel : fElGroup)
    {
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++)
        {
            out << cel->ConnectIndex(ic) << " ";
        }
        out << std::endl;
    }
    out << "Connect indexes of the condensed element \n";
    nc = fCondEl->NConnects();
    for (int ic=0; ic<nc; ic++)
    {
        out << fCondEl->ConnectIndex(ic) << " ";
    }
    out << std::endl;
    
    out << "End of " << __PRETTY_FUNCTION__ << std::endl;
}


void TPZSBFemMultiphysicsElGroup::GroupandCondense(set<int> & condensedmatid)
{
#ifdef PZDEBUG
    if (fElGroup.size() == 0 || condensedmatid.size() == 0) DebugStop();
#endif

    fCondensedEls = new TPZElementGroup(*Mesh());
    const int64_t index = fCondensedEls->Index();

    for (auto celvol : fElGroup)
    {
        if (!celvol) continue;

        auto sbfemvol = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> * >(celvol);

#ifdef PZDEBUG
        if (!sbfemvol) DebugStop();
#endif

        for (auto cel : sbfemvol->SBFemElementVec())
        {
            if (!cel) DebugStop(); 
            auto matid = cel->Reference()->MaterialId();
            auto it = condensedmatid.find(matid);

            if (it == condensedmatid.end()) continue;
            
            fCondensedEls->AddElement(cel);
        }
    }
    
    // Updating NElConnected()
    // fDifPressure and fAverPressure must have NElConnected = 2, and the others = 1
    {
        for(auto cel : fElGroup)
        {
            auto sbfem = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cel);
            if(!sbfem) DebugStop();

            auto elvec = sbfem->SBFemElementVec();

            for (auto el : elvec)
            {
                if ( !el || !(el->Reference()) ) DebugStop();
            
                auto matid = el->Reference()->MaterialId();
                auto it = condensedmatid.find(matid);

                // If it is a condensed element, NElConnected = 1
                if (it != condensedmatid.end() && *(condensedmatid.begin()) != matid)
                {
                    auto nconcel = el->NConnects();
                    for (auto ic = 0; ic < nconcel; ic++)
                    {
                        el->Connect(ic).ResetElConnected();
                        el->Connect(ic).IncrementElConnected();
                    }
                } 
                // Not condensed elements are: fDifPressure and fAverPressure
            }
        }
    }

    bool keepmatrix = true;
    fCondEl = new TPZCondensedCompEl(fCondensedEls, keepmatrix);

    AdjustConnectivities();

    SetLocalIndices(index);
}

void TPZSBFemMultiphysicsElGroup::CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef)
{
    TPZElementMatrixT<STATE> E0, E1, E2;
    ComputeMatrices(E0, E1, E2);
    
    InitializeElementMatrix(ek, ef);

    int n = E0.fMat.Rows();
    auto dim = Mesh()->Dimension();
    
    TPZFMatrix<STATE> E0Inv(E0.fMat);

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
    
    TPZFMatrix<STATE> globmat(2*n,2*n,0.);
    
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);
#else
    cout << "SBFem does not execute for this configuration\n";
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
    cout << "SBFem does not execute for this configuration\n";
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
    cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif

    for (int i=0; i<n; i++) {
        globmat(i,i) -= (dim-2)*0.5;
        globmat(i+n,i+n) += (dim-2)*0.5;
    }
    
    
    static pthread_mutex_t mutex_serial = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex_serial);

    TPZFMatrix<STATE> globmatkeep(globmat);
    TPZFNMatrix<100,complex<double> > eigenVectors;
    TPZManVector<complex<double> > eigenvalues;
    // globmatkeep.SolveEigenProblem(eigenvalues, eigenVectors);
    this->SolveEigenProblemSBFEM(globmatkeep, eigenvalues, eigenVectors);
#ifdef LOG4CXX
    if (loggerfulleigensys->isDebugEnabled())
    {
        stringstream sout;
        EigenVectorsReal(eigenVectors).Print("eigvecfull = ", sout, EMathematicaInput);
        EigenvaluesReal(eigenvalues).Print(sout);
        LOGPZ_DEBUG(loggerfulleigensys, sout.str())
    }
#endif
    
    pthread_mutex_unlock(&mutex_serial);
        
    TPZFMatrix<std::complex<double>> QVectors(n,n,0.);
    fPhi.Resize(n, n);
    TPZManVector<std::complex<double>> eigvalsel(n,0);
    TPZFMatrix<std::complex<double> > eigvecsel(2*n,n,0.),eigvalmat(1,n,0.);
    int count = 0;
    for (int i=0; i<2*n; i++) {
        if (eigenvalues[i].real() < -1.e-6) {
            for (int j=0; j<n; j++) {
                QVectors(j,count) = eigenVectors(j+n,i);
                eigvecsel(j,count) = eigenVectors(j,i);
                eigvecsel(j+n,count) = eigenVectors(j+n,i);
                fPhi(j,count) = eigenVectors(j,i);
            }
            eigvalsel[count] = eigenvalues[i].real();
            eigvalmat(0,count) = eigenvalues[i];
            count++;
        }
    }
    fPhiQ = QVectors;

    if (dim == 2)
    {
        int nstate = Connect(0).NState();
        if (nstate != 2 && nstate != 1) {
            DebugStop();
        }
        if(count != n-nstate) {
            DebugStop();
        }
        int ncon = fCondEl->NConnects();
        int eq=0;
        std::set<int64_t> cornercon;
        BuildCornerConnectList(cornercon);
        for (int ic=ncon/2; ic<ncon; ic++) {
            int64_t conindex = fCondEl->ConnectIndex(ic);
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
            eq += fCondEl->Connect(ic).NShape()*fCondEl->Connect(ic).NState();
        }
    }
    if(dim==3 && count != n)
    {
        DebugStop();
    }
    fEigenvalues = eigvalsel;
        
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

    TPZFMatrix<std::complex<double>> ekloc;
    QVectors.Multiply(fPhiInverse, ekloc);
#ifdef PZDEBUG
    cout << eigenvalues << "\n";
#endif

    ek.fMat.Resize(ekloc.Rows(),ekloc.Cols());
    
    TPZFMatrix<double> ekimag(ekloc.Rows(),ekloc.Cols());
    for (int i=0; i<ekloc.Rows(); i++) {
        for (int j=0; j<ekloc.Cols(); j++) {
            ek.fMat(i,j) = ekloc(i,j).real();
            ekimag(i,j) = ekloc(i,j).imag();
        }
    }

#ifdef LOG4CXX
    if (loggereigensystem->isDebugEnabled())
    {
        std::stringstream sout;
        // sout << fLocalindices << endl;
        ek.fMat.Print("ekloc = ", sout, EMathematicaInput);
        // globmatkeep.Print("matrix = ", sout, EMathematicaInput);
        // eigvecsel.Print("eigvec =", sout, EMathematicaInput);
        eigvalmat.Print("lambda =", sout, EMathematicaInput);
        fPhi.Print("phi = ", sout, EMathematicaInput);
        fPhiInverse.Print("phiinv = ", sout, EMathematicaInput);
        QVectors.Print("qvec = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggereigensystem, sout.str())
    }
#endif
    
    for (auto cel : fElGroup)
    {
        auto sbfem = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad>*>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->SetPhiEigVal(fPhi, fEigenvalues);
    }

    // FLUX
    TPZFMatrix<double> k01 = fCondEl->Matrix().K01();
    n = fPhi.Rows();

    TPZFMatrix<std::complex<double>> fPhiD(2*n, 2*n,0.);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fPhiD(i,j) = fPhi(i,j);
            fPhiD(i+n,j+n) = fPhi(i,j);
        }
    }

    TPZFMatrix<std::complex<double>> tempcomplex(k01.Rows(),k01.Cols(),0.);
    for (auto i = 0; i < k01.Rows(); i++)
    {
        for (auto j = 0; j < k01.Cols(); j++)
        {
            tempcomplex(i,j) = k01(i,j);
        }   
    }
    
    tempcomplex.Multiply(fPhiD, fPhiFlux);
    
    for (auto cel : fElGroup)
    {
        auto sbfem = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->SetPhiFlux(fPhiFlux, fEigenvalues);
    }

#ifdef LOG4CXX  
    if (loggerflux->isDebugEnabled())
    {
        stringstream sout;
        k01.Print("k01 = ", sout, EMathematicaInput);
        fPhiD.Print("phid = ", sout, EMathematicaInput);
        fPhiFlux.Print("fPhiFlux = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerflux, sout.str())
    }
#endif

}

void TPZSBFemMultiphysicsElGroup::ComputeMatrices(TPZElementMatrixT<STATE> &E0, TPZElementMatrixT<STATE> &E1, TPZElementMatrixT<STATE> &E2)
{
    TPZElementMatrixT<STATE> ek(Mesh(),TPZElementMatrixT<STATE>::EK);
    TPZElementMatrixT<STATE> ef(Mesh(),TPZElementMatrixT<STATE>::EF);
    fCondEl->CalcStiff(ek,ef);

    int n = ek.fMat.Rows()/2;
    
    E0.fMat.Resize(n,n);
    E1.fMat.Resize(n,n);
    E2.fMat.Resize(n,n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            E0.fMat(i,j) = ek.fMat(i,j);
            E1.fMat(i,j) = ek.fMat(i+n,j);
            E2.fMat(i,j) = ek.fMat(i+n,j+n);
        }
    }
    {
        std::ofstream out("coefmatrices.txt");
        E0.fMat.Print("E0pz = ", out, EMathematicaInput);
        E1.fMat.Print("E1pz = ", out, EMathematicaInput);
        E2.fMat.Print("E2pz = ", out, EMathematicaInput);
    }

#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCondEl->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
    if (loggercoef->isDebugEnabled())
    {
        std::stringstream out;
        E0.fMat.Print("E0pz = ", out, EMathematicaInput);
        E1.fMat.Print("E1pz = ", out, EMathematicaInput);
        E2.fMat.Print("E2pz = ", out, EMathematicaInput);
        LOGPZ_DEBUG(loggercoef,out.str())
    }
#endif
}

void TPZSBFemMultiphysicsElGroup::InitializeElementMatrix(TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef) const
{
	const int ncon = fCondEl->NConnects();
	int numeq = 0;
	ef.fBlock.SetNBlocks(ncon);
	ek.fBlock.SetNBlocks(ncon);
    
    for (int ic=0; ic<ncon; ic++)
    {
        TPZConnect &c = Connect(ic);
        int blsize = c.NShape()*c.NState();
        numeq += blsize;
        ef.fBlock.Set(ic, blsize);
        ek.fBlock.Set(ic, blsize);
    }
    const int numloadcases = 1;
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
	ef.fMat.Redim(numeq,numloadcases);
	ef.fConnect.Resize(ncon);

    ek.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
	ek.fMat.Redim(numeq,numeq);
	ek.fConnect.Resize(ncon);

	for(int i=0; i<ncon; i++)
    {
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
    std::map<int64_t,TPZOneShapeRestraint>::const_iterator it;
    for (it = fRestraints.begin(); it != fRestraints.end(); it++) 
    {
        ef.fOneRestraints.push_back(it->second);
        ek.fOneRestraints.push_back(it->second);
    }
}//void

void TPZSBFemMultiphysicsElGroup::LoadSolution()
{
    int ndofs = fPhiInverse.Cols();
    int ncoef = fPhiInverse.Rows();
    
    TPZFNMatrix<100, std::complex<double>> uh_local(ncoef, fMesh->Solution().Cols(),0.);
    fCoef.Resize(ncoef,fMesh->Solution().Cols());
    
    TPZFMatrix<STATE> &meshSol = fMesh->Solution();

    int count = 0;
    int nc = fCondEl->NConnects();
    for (int ic=nc/2; ic<nc; ic++)
    {
        TPZConnect &c = fCondEl->Connect(ic);
        int nshape = c.NShape();
        int nstate = c.NState();
        int blsize = nshape*nstate;
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = fMesh->Block().Position(seqnum);
        for (int seq=0; seq < blsize; seq++) {
            for (int c=0; c<uh_local.Cols(); c++)
            {
                uh_local(count+seq,c) = meshSol(pos+seq,c);
            }
        }
        count += blsize;
    }
    fPhiInverse.Multiply(uh_local, fCoef);

    count = 0;
    TPZFMatrix<std::complex<double>> diageigval(ndofs,ndofs);
    for (int i = 0; i < ndofs; i++)
    {
        diageigval(i,i) = -fEigenvalues[i];
    }

    auto duh_local = fCoef;
    TPZFMatrix<std::complex<double>> temp;
    diageigval.Multiply(fCoef, temp);
    fPhi.Multiply(temp, duh_local);

    for (int ic=0; ic<nc/2; ic++)
    {
        TPZConnect &c = fCondEl->Connect(ic);
        int nshape = c.NShape();
        int nstate = c.NState();
        int blsize = nshape*nstate;
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = fMesh->Block().Position(seqnum);
        for (int seq=0; seq < blsize; seq++) {
            for (int c=0; c<uh_local.Cols(); c++)
            {
                meshSol(pos+seq,c) = duh_local(count+seq,c).real();
            }
        }
        count += blsize;
    }

    fPhiInverse.Multiply(duh_local, fCoefD);

    for (auto cel : fElGroup)
    {
        auto sbfem = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->LoadCoef(fCoef, fCoefD); // first coef is related to the average pressure, and the second is the dif pressure
    }

    fCondEl->LoadSolution();
    
#ifdef LOG4CXX
    if (loggereigensystem->isDebugEnabled())
    {
        std::stringstream sout;
        Mesh()->Solution().Print("uh = ", sout, EMathematicaInput);
        uh_local.Print("uh = ", sout, EMathematicaInput);
        fCoef.Print("fcoef = ", sout, EMathematicaInput);
        fCoefD.Print("fcoefd = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggereigensystem,sout.str())
    }
#endif
}

void TPZSBFemMultiphysicsElGroup::AdjustConnectivities()
{
    // Permuting the condensed element's connectivity
    auto ncon = fCondEl->NConnects();
    auto ncontot = this->NConnects();

    TPZManVector<int64_t> perm(ncon), indexes(ncontot);
    int posdif = 0;
    int posaver = ncon/2;
    
    for (auto cel : fElGroup)
    {
        auto sbfem  = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cel);
        auto elvec = sbfem->SBFemElementVec();
        // Adjusting connectivity - dif pressure        
        {
            auto celloc = elvec[0];

            auto nconlocal = celloc->NConnects();
            for (int ic = 0; ic < nconlocal; ic++)
            {
                perm[posdif+ic] = celloc->ConnectIndex(ic);
            }
            posdif += nconlocal;
        }
        // Adjusting connectivity - average pressure
        {
            auto celloc = elvec[6];

            auto nconlocal = celloc->NConnects();
            for (int ic = 0; ic < nconlocal; ic++)
            {
                perm[posaver+ic] = celloc->ConnectIndex(ic);
            }
            posaver += nconlocal;
        }
    }
    fCondEl->PermuteActiveConnects(perm);
    
    for (int64_t i = 0; i < fCondEl->NConnects()/2; i++)
    {
        auto &c = fCondEl->Connect(i);
        c.SetCondensed(true);
    }

    auto nconelgr = this->NConnects();
    auto notcondensed = nconelgr - ncon;
    TPZManVector<int64_t> connectsids(nconelgr);
    for (auto i = 0; i < notcondensed; i++)
    {
        connectsids[i] = this->ConnectIndex(i);
    }
    auto pos = nconelgr - ncon;
    for (auto i = 0; i < ncon; i++)
    {
        connectsids[i+pos] = fCondEl->ConnectIndex(i);
    }
    fCondensedEls->ReorderConnects(connectsids);


    for (auto i = 0; i < ncon/2; i++)
    {
        connectsids[i+ncon/2] = fCondEl->ConnectIndex(i);
    }
    for (auto i = ncon/2; i < ncon; i++)
    {
        connectsids[i-ncon/2] = fCondEl->ConnectIndex(i);
    }
    for (auto i = ncon; i < nconelgr; i++)
    {
        connectsids[i] = this->ConnectIndex(i-ncon);
    }
    
    this->ReorderConnects(connectsids);

#ifdef LOG4CXX
    if (loggersbfemcondensed->isDebugEnabled())
    {
        stringstream sout;
        fCondEl->Print(sout);
        LOGPZ_DEBUG(loggersbfemcondensed, sout.str());
    }
#endif
}

void TPZSBFemMultiphysicsElGroup::SetLocalIndices(int64_t index)
{
    int nc = fCondEl->NConnects();
    TPZManVector<int, 10> firsteqcond(nc + 1, 0);
    for (int ic = 0; ic < nc; ic++)
    {
        TPZConnect &c = fCondEl->Connect(ic);
        firsteqcond[ic + 1] = firsteqcond[ic] + c.NShape() * c.NState();
    }
    auto numeqcondensed = *max_element(firsteqcond.begin(), firsteqcond.end());

    std::map<int64_t, int> globtolocal;
    nc = this->NConnects();
    TPZManVector<int, 10> firsteq(nc + 1, 0);
    for (int ic = 0; ic < nc; ic++)
    {
        globtolocal[this->ConnectIndex(ic)] = ic;
        TPZConnect &c = this->Connect(ic);
        firsteq[ic + 1] = firsteq[ic] + c.NShape() * c.NState();
    }

    for (auto cel : fElGroup)
    {
        auto sbfem  = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cel);
        auto elvec = sbfem->SBFemElementVec();
        TPZCompEl *celskeleton = elvec[6];

        // Determining the number of equations for each element:
        int neq = 0;
        auto ncskeleton = celskeleton->NConnects();
        for (int ic = 0; ic < ncskeleton; ic++)
        {
            TPZConnect &c = celskeleton->Connect(ic);
            neq += c.NShape() * c.NState();
        }

        // LOCAL INDICES PRESSURE
        TPZManVector<int64_t> localindices(neq);
        int count = 0;
        for (int ic = 0; ic < celskeleton->NConnects(); ic++)
        {
            int64_t cindex = celskeleton->ConnectIndex(ic);
#ifdef PZDEBUG
            if (globtolocal.find(cindex) == globtolocal.end()) DebugStop();
#endif
            int locfirst = firsteq[globtolocal[cindex]];
            TPZConnect &c = celskeleton->Connect(ic);
            int neq = c.NShape() * c.NState();
            for (int eq = 0; eq < neq; eq++)
            {
                localindices[count++] = locfirst + eq;
            }
        }

        if (sbfem->Element(6)->Reference()->NodeIndex(0) != sbfem->Element(5)->Reference()->NodeIndex(0))
        {
// #ifdef PZDEBUG
//             cout << "Inverting local indices for subelement " << sbfem->Reference()->Index() << endl ;
// #endif
            TPZManVector<int64_t> localindicescopy(localindices);
            localindices[0] = localindicescopy[1];
            localindices[1] = localindicescopy[0];
        }

        // LOCAL INDICES FLUX
        int neqflux = 0;
        TPZCompEl *celflux = elvec[1];
        auto ncflux = celflux->NConnects();
        for (int ic = 0; ic < ncflux; ic++)
        {
            TPZConnect &c = celflux->Connect(ic);
            neqflux += c.NShape() * c.NState();
        }
        int neqint = 0;
        TPZCompEl *celinternal = elvec[3];
        auto ncint = celinternal->NConnects() - ncskeleton - 2*ncflux;
        for (int ic = 0; ic < ncint; ic++)
        {
            TPZConnect &c = celinternal->Connect(ic);
            neqint += c.NShape() * c.NState();
        }
        
        TPZManVector<int64_t> localindicesflux(neqflux);
        count = 0;
        for (int ic = 0; ic < ncflux; ic++)
        {
            int64_t cindex = celflux->ConnectIndex(ic);
#ifdef PZDEBUG
            if (globtolocal.find(cindex) == globtolocal.end()) DebugStop();
#endif
            TPZConnect &c = celflux->Connect(ic);
            int neqflux = c.NShape() * c.NState();
            int locfirst = firsteq[globtolocal[cindex]];
            for (int eq = 0; eq < neqflux; eq++)
            {
                localindicesflux[count++] = locfirst + eq;
            }
        }
        for (int id = 0; id < localindicesflux.size(); id++)
        {
            localindicesflux[id] -= numeqcondensed;
        }
        
        TPZManVector<int64_t> localindicesint(neqint);
        count = 0;
        for (int ic = 0; ic < ncint; ic++)
        {
            int64_t cindex = celinternal->ConnectIndex(ic);
#ifdef PZDEBUG
            if (globtolocal.find(cindex) == globtolocal.end()) DebugStop();
#endif
            TPZConnect &c = celinternal->Connect(ic);
            int neqint = c.NShape() * c.NState();
            int locfirst = firsteq[globtolocal[cindex]];
            for (int eq = 0; eq < neqint; eq++)
            {
                localindicesint[count++] = locfirst + eq;
            }
        }
        for (int id = 0; id < localindicesint.size(); id++)
        {
            localindicesint[id] -= numeqcondensed;
        }

        sbfem->SetLocalIndices(localindices,localindicesint,localindicesflux);
    }
}

/*

Helpful doc to understand connectivity in SBFEMMultiphysics approach:

ORDER CONNECTS INSIDE A TPZSBFemVolumeHdiv
fElGroup[0]: 36 37 38 -> fDifPressure
fElGroup[1]: 3 -> fFluxLeft
fElGroup[2]: 3 36 37 38 -> fInterface
fElGroup[3]: 0 1 2 3 4 | 24 25 26 -> fFlux | fInternalPressure
fElGroup[4]: 4 -> fFluxRight
fElGroup[5]: 4 48 49 50 -> fInterface
fElGroup[6]: 48 49 50 -> fSkeleton (External Pressure)

Material IDs:
fSkeleton = 6
fInterface = 7 
fRightFlux = 8
fInternal = 9
fLeftFlux = 10
fDifPressure = 11

ORDER CONNECTS INSIDE A TPZSBFemMultiphysicsElGroup
48 49 50 51 52 53 54 55 56 57 58 59 -> fSkeleton
36 37 38 39 40 41 42 43 44 45 46 47 -> fDifPressure
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 -> fFlux
24 25 26 27 28 29 30 31 32 33 34 35 -> fInternalPressure

ORDER CONNECTS INSIDE THE TPZCondensedCompEl (fCondEl)
36 37 38 39 40 41 42 43 44 45 46 47 -> fDifPressure
48 49 50 51 52 53 54 55 56 57 58 59 -> fSkeleton
Note: it is different from the order of TPZSBFemMultiphysicsElGroup. Here, 1st we have the fDifPressure 1st.

*/
