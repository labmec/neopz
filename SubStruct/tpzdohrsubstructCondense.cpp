/**
 * @file
 * @brief Contains the implementation of the TPZDohrSubstructCondense methods. 
 */

#include "tpzdohrsubstructCondense.h"
#include "tpzverysparsematrix.h"
#include <iostream>
#include "pzlog.h"
#include "TPZSimpleTimer.h"

#ifdef PZ_LOG
static TPZLogger logger("substruct.dohrsubstructcondense");
#endif

using namespace std;

template<class TVar>
typename TPZDohrSubstructCondense<TVar>::EWeightType TPZDohrSubstructCondense<TVar>::fWeightType = TPZDohrSubstructCondense<TVar>::CorrectWeight;


template<class TVar>
TPZDohrSubstructCondense<TVar>::TPZDohrSubstructCondense() : fNEquations(-1), fNumInternalEquations(-1), fNumExternalEquations(-1), wasRealloc(false)
{
	//Inicializacao
}

template<class TVar>
TPZDohrSubstructCondense<TVar>::~TPZDohrSubstructCondense()
{
	//Limpesa
}

/**
 * It computes the local contribution to r(c).
 * The method LoadWeightedResidual must be called before this one.
 */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::Contribute_rc_local(TPZFMatrix<TVar> &residual_local, TPZFMatrix<TVar> &rc_local) const
{
	fPhiC_Weighted_Condensed.Multiply(residual_local, rc_local, 1);
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::Contribute_Kc(TPZMatrix<TVar> &Kc, TPZVec<int> &coarseindex) {
	int i;
	int j;
	for (i=0;i<fCoarseNodes.NElements();i++) {
		for (j=0;j<fCoarseNodes.NElements();j++) {
			if ((Kc.IsSymmetric() && coarseindex[j] >= coarseindex[i]) || !Kc.IsSymmetric()) {
				Kc(coarseindex[i],coarseindex[j]) += fKCi(i,j);				
			}
		}
	}
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::Contribute_v1_local(TPZFMatrix<TVar> &v1_local, TPZFMatrix<TVar> &invKc_rc_local) const {
	int neqs = fNumExternalEquations;
	v1_local.Resize(neqs, 1);
	fPhiC_Weighted_Condensed.Multiply(invKc_rc_local,v1_local);
}

/**
 * It computes the local contribution to v2.
 */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::Contribute_v2_local(TPZFMatrix<TVar> &residual_local, TPZFMatrix<TVar> &v2_local)
{
	const TPZVec<int> &scatter = ScatterVec(ExternalFirst, Submesh);
	int ncoarse = fCoarseNodes.NElements();
    int ncols = residual_local.Cols();
	TPZFMatrix<TVar> LocalWeightedResidual(fNEquations+ncoarse,ncols,0.);
	int ninput = residual_local.Rows();
	int i;
    for (int ic=0; ic<ncols; ic++) 
    {
        for (i=0;i<ninput;i++) 
        {
            LocalWeightedResidual(scatter[i],ic) += fWeights[scatter[i]] * residual_local(i,ic);
        }
    }
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		LocalWeightedResidual.Print("LocalWeightedResidual ",sout);
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif
	fMatRedComplete->SetF(LocalWeightedResidual);
	TPZFMatrix<TVar> U1(ncoarse,ncols,0.), UGlobal(fNEquations+ncoarse,ncols,0.);
	fMatRedComplete->U1(U1);
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		U1.Print("U1 ",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fMatRedComplete->UGlobal(U1,UGlobal);
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		UGlobal.Print("UGlobal ",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	TPZVec<int> &gather2 = GatherVec(Submesh, ExternalFirst);
	v2_local.Resize(ninput, ncols);
    for (int ic=0; ic<ncols; ic++) 
    {
        for (i=0;i<ninput;i++) 
        {
            v2_local(i,ic) = fWeights[gather2[i]] * UGlobal(gather2[i],ic);
        }
    }
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::Contribute_v3_local(TPZFMatrix<TVar> &v1Plusv2, TPZFMatrix<TVar> &v3) const {
	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::Print(std::ostream &out) const
{
	out << __PRETTY_FUNCTION__ << std::endl;
	out << "fNEquations " << fNEquations << std::endl;
	out << "fNumInternalEquations " << fNumInternalEquations << std::endl;
	out << "fNumExternalEquations " << fNumExternalEquations << std::endl;
	const TPZVec<int> &internaleq = GatherVec(Submesh, InternalFirst);
	out << "Internal equations-first-Gather " << internaleq << std::endl;
	const TPZVec<int> &exteq = GatherVec(Submesh, ExternalFirst);
	out << "External equations-first-Gather " << exteq << std::endl;
	out << "Numbering of the Coarse nodes in the global mesh fCoarseNodes " << fCoarseNodes << std::endl;
	fMatRedComplete->Print("The matrix which computes phi and fKCi ",out,EMathematicaInput);
	fMatRed->Print("The matrix with the internal nodes condensed ", out,EMathematicaInput);
	fKCi.Print("Coarse Matrix fKCi",out);
	fPhiC.Print("fPhiC : ",out,EMathematicaInput);
	out << "fWeights = " << fWeights  << endl;
	fPhiC_Weighted_Condensed.Print("fPhiC_Weighted_Condensed", out);
	fLocalLoad.Print("fLocalLoad ", out);
	fLocalWeightedResidual.Print("fLocalWeightedResidual ", out);
	fAdjustSolution.Print("fAdjustSolution", out);
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::SolveSystemPhi() {
	int ncoarse = fCoarseNodes.NElements();
	int i,j;
	TPZFMatrix<TVar> rhs(this->fNEquations+ncoarse,ncoarse,0.);
	for(i=0; i<ncoarse; i++)
	{
		rhs(fNEquations+i,i) = 1.;
	}
	fMatRedComplete->SetF(rhs);
	fKCi.Resize(ncoarse,ncoarse);
	fMatRedComplete->U1(fKCi);
	fMatRedComplete->UGlobal(fKCi,rhs);
	for(i=0; i<ncoarse; i++) 
	{
		fKCi(i,i) += 10.;
		for(j=0; j<ncoarse; j++)
		{
			fKCi(i,j) = -fKCi(i,j);
		}
		fWeights[fCoarseNodes[i]] = fKCi(i,i);
	}
	fPhiC.Resize(fNEquations,ncoarse);
	for(i=0; i<fNEquations; i++)
	{
		for(j=0; j<ncoarse; j++)
		{
			fPhiC(i,j) = rhs(i,j);
		}
	}
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::ContributeDiagonalLocal(TPZFMatrix<TVar> &StiffnessDiagLocal) {
	int i;
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Weight used for assembly" << fWeights;
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
	
	TPZVec<int> &gather = GatherVec(Submesh, ExternalFirst);
	
	int neqs = fNumExternalEquations;
	StiffnessDiagLocal.Resize(neqs,1);
	for (i=0;i<neqs;i++) {
		StiffnessDiagLocal(i,0) = fWeights[gather[i]];
	}
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::ComputeWeightsLocal(TPZFMatrix<TVar> &StiffnessDiagLocal) {
	int i;
	//fWeights.Fill(1.);
	TPZVec<int> &scatter = ScatterVec(ExternalFirst, Submesh);
	int neqs = fNumExternalEquations;
	for (i=0;i<neqs;i++) {
		fWeights[scatter[i]] = fWeights[scatter[i]] / StiffnessDiagLocal(i,0);
	}
	for(; i<fNEquations; i++)
	{
		fWeights[scatter[i]] = 1.;
	}
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Weights = " <<  fWeights;
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
	TPZVec<int> &gather = GatherVec(Submesh, ExternalFirst);
	int c,nc = fPhiC.Cols();
	fPhiC_Weighted_Condensed.Resize(neqs,nc);
	for (i=0;i<neqs;i++) 
	{
		for(c=0; c<nc; c++)
		{
			fPhiC_Weighted_Condensed(i,c) = fPhiC(gather[i],c)*fWeights[gather[i]];
		}
	}
}

/**
 * Computes the condensed right hand side for the substructure
 * @param rhs the right hand side ordered external first
 */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::ContributeRhs(TPZFMatrix<TVar> &rhs)
{
	int nglob = fNEquations;
	int ncols = rhs.Cols();
	typedef std::pair<ENumbering,ENumbering> Numbering;
	Numbering relat(ExternalFirst,Submesh);
	Numbering relat2(InternalFirst,Submesh);
	typename std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> >::const_iterator itrelat, itrelat2, itrelat3;
	itrelat = fPermutationsScatter.find(relat);
	itrelat2 = fPermutationsScatter.find(relat2);
	if(itrelat == fPermutationsScatter.end() || itrelat2 == fPermutationsScatter.end() )
	{
		DebugStop();
		return;
	}
	TPZFMatrix<TVar> resglobal(nglob,ncols,0.),resloc(fNumExternalEquations,ncols,0.);
	PermuteGather (itrelat2->second, fLocalLoad, resglobal, 0, nglob);
	fMatRed->SetF(resglobal);
	fMatRed->F1Red(resloc);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        resglobal.Print("resglobal ",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	resglobal.Zero();
	PermuteScatter(itrelat->second, resloc, resglobal, 0, fNumExternalEquations);
	PermuteGather(itrelat->second, resglobal, rhs, 0, fNumExternalEquations);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		resloc.Print("Condensed F", sout);
		resglobal.Print("Scattered F", sout);
		rhs.Print("External first", sout);
		sout << "vector for scatter " << itrelat->second << std::endl;
		sout << "vector for gather " << itrelat->second << std::endl;
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
#ifdef PZDEBUG 
	TPZFMatrix<TVar> test(resloc);
	test -= rhs;
	//	TVar err = Norm(test);
#endif
}

/**
 * Computes the global solution based on the interface solution
 */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::UGlobal(TPZFMatrix<TVar> &UGlob, TPZFMatrix<TVar> &USub)
{
	int nglob = fNEquations;
	int ncols = UGlob.Cols();
	USub.Redim (nglob, ncols);
	typedef std::pair<ENumbering,ENumbering> Numbering;
	Numbering relat(Submesh,ExternalFirst);
	Numbering relat2(InternalFirst,Submesh);
	typename std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> >::const_iterator itrelat, itrelat2, itrelat3;
	itrelat2 = fPermutationsScatter.find(relat2);
	itrelat = fPermutationsScatter.find(relat);
	if(itrelat2 == fPermutationsScatter.end() || itrelat == fPermutationsScatter.end())
	{
		DebugStop();
		return;
	}
    TPZFMatrix<TVar> uloc(nglob,ncols,0.);
    {
        TPZFNMatrix<200> uext(fNumExternalEquations,ncols);
        fMatRed->UGlobal2(UGlob,uloc);
//			PermuteGather(itrelat->second, UGlob, uext, 0, fNumExternalEquations);
        
    }
    PermuteScatter(itrelat2->second, uloc, USub , 0, nglob);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		//		uext.Print("Boundary node solution", sout);
        uloc.Print("Complete solution internal first", sout);
		UGlob.Print("submesh solution", sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::ContributeKULocal(const TVar alpha, const TPZFMatrix<TVar> &u, TPZFMatrix<TVar> &z) const
{
	int i,j;
	int nglob = fNEquations;
	int ncols = u.Cols();
	int neqs = fNEquations-fNumInternalEquations;
	TPZFMatrix<TVar> uloc(nglob,ncols,0.), uborder(fNEquations-fNumInternalEquations,ncols,0.),resborder(fNEquations-fNumInternalEquations,ncols,0.);
	typedef std::pair<ENumbering,ENumbering> Numbering;
	Numbering relat(ExternalFirst,Submesh);
	Numbering relat2(InternalFirst,Submesh);
	typename std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> >::const_iterator itrelat, itrelat2;
	itrelat = fPermutationsScatter.find(relat);
	itrelat2 = fPermutationsScatter.find(relat2);
	if(itrelat == fPermutationsScatter.end() || itrelat2 == fPermutationsScatter.end())
	{
		return;
	}
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Scatter from external do submesh" << itrelat->second << std::endl;
		sout << "Scatter from submesh to internal" << itrelat2->second << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	PermuteScatter(itrelat->second, u, uloc, 0, u.Rows());
	PermuteGather(itrelat2->second, uloc, uborder, fNumInternalEquations, fNEquations);
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		u.Print("Input matrix ", sout);
		uloc.Print("Natural ordering matrix",sout);
		uborder.Print("Matrix passed to the matred object", sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fMatRed->Multiply(uborder,resborder);
	
	TPZFMatrix<TVar> resglobal(nglob,ncols,0.),resloc(fNumExternalEquations,ncols,0.);
	PermuteScatter(itrelat2->second, resborder, resglobal, fNumInternalEquations, fNEquations);
	PermuteGather(itrelat->second, resglobal, resloc, 0, u.Rows());
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Value of the local multiply = ";
		resborder.Print("resborder " ,sout);
		resglobal.Print("resglobal ",sout);
		resloc.Print("resloc ",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
	int zcols = z.Cols();
	for (i=0;i<neqs;i++) {
		/* Sum row "i" of temp1 with row "fGlobalNodes[i]" of z */
		for (j=0;j<zcols;j++) {
			z(i,j) += alpha*resloc(i,j);
		}
	}
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::Initialize() {
	PrepareSystems();
	SolveSystemPhi();
	
}

template<class TVar>
void TPZDohrSubstructCondense<TVar>::PrepareSystems() {
}

/**
 * Adjust the residual to reflect a static condensation
 * The residual corresponding to the internal nodes will be zeroed
 */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::AdjustResidual(TPZFMatrix<TVar> &r_global)
{
	//	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
}

/**
 * Add the internal solution to the final result
 */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::AddInternalSolution(TPZFMatrix<TVar> &sol)
{	
}

/**
 * Apply a scatter permutation to the input vector using a scatter permutation
 * output[permute[i]] = input[i-first], first <= i < last
 * this method does not resize the elements
 */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::PermuteScatter(const TPZVec<int> &permute, const TPZFMatrix<TVar> &input, TPZFMatrix<TVar> &output, int first, int last) 
{
	int i,j,ncol = input.Cols();
	for(i=first; i<last; i++) for(j=0; j<ncol; j++)
	{
		output(permute[i],j) = input.GetVal(i-first,j);
	}
}

/**
 * Apply a gather permutation to the input vector using a scatter permutation
 * output[i-first] = input[permute[i]], first <= i < last
 * this method does not resize the elements
 */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::PermuteGather(const TPZVec<int> &permute, const TPZFMatrix<TVar> &input, TPZFMatrix<TVar> &output, int first, int last)
{
	int i,j,ncol = input.Cols();
	for(i=first; i<last; i++) for(j=0; j<ncol; j++)
	{
		output(i-first,j) = input.GetVal(permute[i],j);
	}
}

static TPZVec<int> dummyvec;
template<class TVar>
const TPZVec<int> &TPZDohrSubstructCondense<TVar>::GatherVec(ENumbering origin, ENumbering destination) const
{
	typename std::map<std::pair<ENumbering,ENumbering>, TPZVec<int> >::const_iterator it;
	it = fPermutationsScatter.find(std::pair<ENumbering,ENumbering>(destination,origin));
	if(it != fPermutationsScatter.end())
	{
		return it->second;
	}
	else
	{
		LOGPZ_ERROR(logger,"Gathervec not found");
		dummyvec.Resize(0);
		return dummyvec;
	}
}

template<class TVar>
const TPZVec<int> &TPZDohrSubstructCondense<TVar>::ScatterVec(ENumbering origin, ENumbering destination) const
{
	typename std::map<std::pair<ENumbering,ENumbering>, TPZVec<int> >::const_iterator it;
	it = fPermutationsScatter.find(std::pair<ENumbering,ENumbering>(destination,origin));
	if(it != fPermutationsScatter.end())
	{
		return it->second;
	}
	else
	{
		LOGPZ_ERROR(logger,"Gathervec not found");
		dummyvec.Resize(0);
		return dummyvec;
	}
}

/** @brief method for streaming the object to a stream */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::Write(TPZStream &out, int withclassid) const
{    
    SAVEABLE_STR_NOTE(out,"fMatRedComplete");
    if(fMatRedComplete)
    {
        int one = 1;
        out.Write(&one);
        fMatRedComplete->Write(out,0);
    }
    else {
        int zero = 0;
        out.Write(&zero);
    }
    SAVEABLE_STR_NOTE(out,"fNEquations");
    out.Write(&fNEquations);
    SAVEABLE_STR_NOTE(out,"fNumInternalEquations");
    out.Write(&fNumInternalEquations);
    SAVEABLE_STR_NOTE(out,"fNumExternalEquations");
    out.Write(&fNumExternalEquations);
    std::cout << fNEquations << " " << fNumInternalEquations << " " << fNumExternalEquations << std::endl;
    out.Write( fCoarseNodes);
    std::cout << fCoarseNodes << std::endl;
    int one(1),two(2),three(3),four(4);
    
    out.Write(&one);
    fPhiC.Write(out, 0);
    out.Write(&two);
    fPhiC_Weighted_Condensed.Write(out, 0);
    out.Write(&three);
    out.Write( fWeights);
    fKCi.Write(out, 0);
    out.Write(&four);
    typename std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> >::const_iterator it;
    int sc = fPermutationsScatter.size();
    out.Write(&sc);
    for (it=fPermutationsScatter.begin(); it != fPermutationsScatter.end(); it++) {
        int a = it->first.first;
        int b = it->first.second;
        out.Write(&a);
        out.Write(&b);
        out.Write( it->second);
    }
    if (fMatRed) {
        int one = 1;
        out.Write(&one);
        fMatRed->Write(out, 0);
    }
    else {
        int zero = 0;
        out.Write(&zero);
    }
    
    fLocalLoad.Write(out, 0);
    fLocalWeightedResidual.Write(out, 0);
    fAdjustSolution.Write(out, 0);

}

/** @brief method for reading the object for a stream */
template<class TVar>
void TPZDohrSubstructCondense<TVar>::Read(TPZStream &input, void *context)
{
    int a;
    SAVEABLE_SKIP_NOTE(input);
    input.Read(&a);
    if (a) {
        fMatRedComplete = new TPZMatRed<TVar, TPZFMatrix<TVar> >;
        fMatRedComplete->Read(input,0);
    }
    
    SAVEABLE_SKIP_NOTE(input);
    input.Read(&fNEquations);
    SAVEABLE_SKIP_NOTE(input);
    input.Read(&fNumInternalEquations);
    SAVEABLE_SKIP_NOTE(input);
    input.Read(&fNumExternalEquations);
    std::cout << fNEquations << " " << fNumInternalEquations << " " << fNumExternalEquations << std::endl;
    input.Read( fCoarseNodes);
    std::cout << fCoarseNodes << std::endl;
    int one(-1),two(-2),three(-3),four(-4);

    input.Read(&one);
    fPhiC.Read(input, 0);
    input.Read(&two);
    fPhiC_Weighted_Condensed.Read(input, 0);
    input.Read(&three);
    input.Read( fWeights);
    fKCi.Read(input, 0);
    input.Read(&four);
    int nc;
    input.Read(&nc);
    for (int ic=0; ic<nc; ic++) {
        int a;
        int b;
        input.Read(&a);
        input.Read(&b);
        ENumbering orig = (ENumbering)(a),dest = (ENumbering)(b);
        std::pair<ENumbering, ENumbering> p(orig,dest);
        input.Read( fPermutationsScatter[p]);
    }
    int control;
    input.Read(&control);
    if(control)
    {
        fMatRed = new TPZMatRed<TVar, TPZFMatrix<TVar> >;
        fMatRed->Read(input, 0);
    }
    fLocalLoad.Read(input, 0);
    fLocalWeightedResidual.Read(input, 0);
    fAdjustSolution.Read(input, 0);
    
    
}

template class TPZDohrSubstructCondense<float>;
template class TPZDohrSubstructCondense<double>;
template class TPZDohrSubstructCondense<long double>;

template class TPZDohrSubstructCondense<std::complex<float> >;
template class TPZDohrSubstructCondense<std::complex<double> >;
template class TPZDohrSubstructCondense<std::complex<long double> >;

