/**
 * @file
 * @brief Contains the implementation of the TPZDohrSubstructCondense methods. 
 */
/***************************************************************************
 *   Copyright (C) 2006 by Philippe Devloo   *
 *   phil@fec.unicamp.br   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "tpzdohrsubstructCondense.h"
#include <iostream>
#include "pzlog.h"
#include "TPZfTime.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("substruct.dohrsubstructcondense"));
#endif

using namespace std;

TPZDohrSubstructCondense::EWeightType TPZDohrSubstructCondense::fWeightType = TPZDohrSubstructCondense::CorrectWeight;


TPZDohrSubstructCondense::TPZDohrSubstructCondense() 
{
	//Inicializacao
}


TPZDohrSubstructCondense::~TPZDohrSubstructCondense()
{
	//Limpesa
}


/**
 * It computes the local contribution to r(c).
 * The method LoadWeightedResidual must be called before this one.
 */
void TPZDohrSubstructCondense::Contribute_rc_local(TPZFMatrix &residual_local, TPZFMatrix &rc_local)
{
	fPhiC_Weighted_Condensed.Multiply(residual_local, rc_local, 1, 1);
}


void TPZDohrSubstructCondense::Contribute_Kc(TPZMatrix &Kc, TPZVec<int> &coarseindex) {
	int i;
	int j;
	for (i=0;i<fCoarseNodes.NElements();i++) {
		for (j=0;j<fCoarseNodes.NElements();j++) {
			if ((Kc.IsSimetric() && coarseindex[j] >= coarseindex[i]) || !Kc.IsSimetric()) {
				Kc(coarseindex[i],coarseindex[j]) += fKCi(i,j);				
			}
		}
	}
}


void TPZDohrSubstructCondense::Contribute_v1_local(TPZFMatrix &v1_local, TPZFMatrix &invKc_rc_local) {
	int neqs = fNumExternalEquations;
	v1_local.Resize(neqs, 1);
	fPhiC_Weighted_Condensed.Multiply(invKc_rc_local,v1_local);
}


/**
 * It computes the local contribution to v2.
 */
void TPZDohrSubstructCondense::Contribute_v2_local(TPZFMatrix &residual_local, TPZFMatrix &v2_local)
{
	TPZVec<int> &scatter = ScatterVec(ExternalFirst, Submesh);
	int ncoarse = fCoarseNodes.NElements();
	TPZFMatrix LocalWeightedResidual(fNEquations+ncoarse,1,0.);
	int ninput = residual_local.Rows();
	int i;
	for (i=0;i<ninput;i++) 
	{
		LocalWeightedResidual(scatter[i],0) += fWeights[scatter[i]] * residual_local(i,0);
	}
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		LocalWeightedResidual.Print("LocalWeightedResidual ",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fMatRedComplete->SetF(LocalWeightedResidual);
	TPZFMatrix U1(ncoarse,1,0.), UGlobal(fNEquations+ncoarse,0.);
	fMatRedComplete->U1(U1);
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		U1.Print("U1 ",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fMatRedComplete->UGlobal(U1,UGlobal);
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		UGlobal.Print("UGlobal ",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	TPZVec<int> &gather2 = GatherVec(Submesh, ExternalFirst);
	v2_local.Resize(ninput, 1);
	for (i=0;i<ninput;i++) 
	{
		v2_local(i,0) = fWeights[gather2[i]] * UGlobal(gather2[i],0);
	}
}

void TPZDohrSubstructCondense::Contribute_v3_local(TPZFMatrix &v1Plusv2, TPZFMatrix &v3) const {
	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
}

void TPZDohrSubstructCondense::Print(std::ostream &out) const
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
	fMatRedComplete->Print("The matrix which computes phi and fKCi ",out);
	fMatRed->Print("The matrix with the internal nodes condensed ", out);
	fKCi.Print("Coarse Matrix fKCi",out);
	fPhiC.Print("fPhiC : ",out,EMathematicaInput);
	out << "fWeights = " << fWeights  << endl;
	fPhiC_Weighted_Condensed.Print("fPhiC_Weighted_Condensed", out);
	fLocalLoad.Print("fLocalLoad ", out);
	fLocalWeightedResidual.Print("fLocalWeightedResidual ", out);
	fAdjustSolution.Print("fAdjustSolution", out);
	
}

void TPZDohrSubstructCondense::SolveSystemPhi() {
	int ncoarse = fCoarseNodes.NElements();
	int i,j;
	TPZFMatrix rhs(this->fNEquations+ncoarse,ncoarse,0.);
	for(i=0; i<ncoarse; i++)
	{
		rhs(fNEquations+i,i) = 1.;
	}
	fMatRedComplete->SetF(rhs);
	fMatRedComplete->SetF0IsComputed(true);
	fMatRedComplete->SetF1IsReduced(true);
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
	//	fPhiC *= -1.;
}



void TPZDohrSubstructCondense::ContributeDiagonalLocal(TPZFMatrix &StiffnessDiagLocal) {
	int i;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Weight used for assembly" << fWeights;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	TPZVec<int> &gather = GatherVec(Submesh, ExternalFirst);
	
	int neqs = fNumExternalEquations;
	StiffnessDiagLocal.Resize(neqs,1);
	for (i=0;i<neqs;i++) {
		StiffnessDiagLocal(i,0) = fWeights[gather[i]];
	}
}

void TPZDohrSubstructCondense::ComputeWeightsLocal(TPZFMatrix &StiffnessDiagLocal) {
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
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Weights = " <<  fWeights;
		LOGPZ_DEBUG(logger,sout.str())
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
void TPZDohrSubstructCondense::ContributeRhs(TPZFMatrix &rhs)
{
	int nglob = fNEquations;
	int ncols = rhs.Cols();
	typedef std::pair<ENumbering,ENumbering> Numbering;
	Numbering relat(ExternalFirst,Submesh);
	Numbering relat2(InternalFirst,Submesh);
	std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> >::const_iterator itrelat, itrelat2, itrelat3;
	itrelat = fPermutationsScatter.find(relat);
	itrelat2 = fPermutationsScatter.find(relat2);
	if(itrelat == fPermutationsScatter.end() || itrelat2 == fPermutationsScatter.end() )
	{
		DebugStop();
		return;
	}
	TPZFMatrix resglobal(nglob,ncols,0.),resloc(fNumExternalEquations,ncols,0.);
	PermuteGather (itrelat2->second, fLocalLoad, resglobal, 0, nglob);
	fMatRed->SetF(resglobal);
	resloc = fMatRed->F1Red();
	resglobal.Zero();
	PermuteScatter(itrelat->second, resloc, resglobal, 0, fNumExternalEquations);
	PermuteGather(itrelat->second, resglobal, rhs, 0, fNumExternalEquations);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		resloc.Print("Condensed F", sout);
		resglobal.Print("Scattered F", sout);
		rhs.Print("External first", sout);
		sout << "vector for scatter " << itrelat->second << std::endl;
		sout << "vector for gather " << itrelat->second << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
#ifdef DEBUG 
	TPZFMatrix test(resloc);
	test -= rhs;
	//	REAL err = Norm(test);
#endif
}

/**
 * Computes the global solution based on the interface solution
 */
void TPZDohrSubstructCondense::UGlobal(TPZFMatrix &UGlob, TPZFMatrix &USub)
{
	int nglob = fNEquations;
	int ncols = UGlob.Cols();
	USub.Redim (nglob, ncols);
	typedef std::pair<ENumbering,ENumbering> Numbering;
	Numbering relat(Submesh,ExternalFirst);
	Numbering relat2(InternalFirst,Submesh);
	std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> >::const_iterator itrelat, itrelat2, itrelat3;
	itrelat2 = fPermutationsScatter.find(relat2);
	itrelat = fPermutationsScatter.find(relat);
	if(itrelat2 == fPermutationsScatter.end() || itrelat == fPermutationsScatter.end())
	{
		DebugStop();
		return;
	}
		TPZFMatrix uloc(nglob,ncols,0.);
		{
			TPZFNMatrix<200> uext(fNumExternalEquations,ncols);
			fMatRed->UGlobal2(UGlob,uloc);
//			PermuteGather(itrelat->second, UGlob, uext, 0, fNumExternalEquations);
			
		}
		PermuteScatter(itrelat2->second, uloc, USub , 0, nglob);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		//		uext.Print("Boundary node solution", sout);
        uloc.Print("Complete solution internal first", sout);
		UGlob.Print("submesh solution", sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}


void TPZDohrSubstructCondense::ContributeKULocal(const REAL alpha, const TPZFMatrix &u, TPZFMatrix &z) const
{
	int i,j;
	int nglob = fNEquations;
	int ncols = u.Cols();
	int neqs = fNEquations-fNumInternalEquations;
	TPZFMatrix uloc(nglob,ncols,0.), uborder(fNEquations-fNumInternalEquations,ncols,0.),resborder(fNEquations-fNumInternalEquations,ncols,0.);
	typedef std::pair<ENumbering,ENumbering> Numbering;
	Numbering relat(ExternalFirst,Submesh);
	Numbering relat2(InternalFirst,Submesh);
	std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> >::const_iterator itrelat, itrelat2;
	itrelat = fPermutationsScatter.find(relat);
	itrelat2 = fPermutationsScatter.find(relat2);
	if(itrelat == fPermutationsScatter.end() || itrelat2 == fPermutationsScatter.end())
	{
		return;
	}
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Scatter from external do submesh" << itrelat->second << std::endl;
		sout << "Scatter from submesh to internal" << itrelat2->second << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	PermuteScatter(itrelat->second, u, uloc, 0, u.Rows());
	PermuteGather(itrelat2->second, uloc, uborder, fNumInternalEquations, fNEquations);
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		u.Print("Input matrix ", sout);
		uloc.Print("Natural ordering matrix",sout);
		uborder.Print("Matrix passed to the matred object", sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fMatRed->Multiply(uborder,resborder);
	
	TPZFMatrix resglobal(nglob,ncols,0.),resloc(fNumExternalEquations,ncols,0.);
	PermuteScatter(itrelat2->second, resborder, resglobal, fNumInternalEquations, fNEquations);
	PermuteGather(itrelat->second, resglobal, resloc, 0, u.Rows());
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Value of the local multiply = ";
		resborder.Print("resborder " ,sout);
		resglobal.Print("resglobal ",sout);
		resloc.Print("resloc ",sout);
		LOGPZ_DEBUG(logger,sout.str())
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

void TPZDohrSubstructCondense::Initialize() {
	PrepareSystems();
	SolveSystemPhi();
	
}

void TPZDohrSubstructCondense::PrepareSystems() {
}

/**
 * Adjust the residual to reflect a static condensation
 * The residual corresponding to the internal nodes will be zeroed
 */
void TPZDohrSubstructCondense::AdjustResidual(TPZFMatrix &r_global)
{
	//	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
}

/**
 * Add the internal solution to the final result
 */
void TPZDohrSubstructCondense::AddInternalSolution(TPZFMatrix &sol)
{	
}

/**
 * Apply a scatter permutation to the input vector using a scatter permutation
 * output[permute[i]] = input[i-first], first <= i < last
 * this method does not resize the elements
 */
void TPZDohrSubstructCondense::PermuteScatter(const TPZVec<int> &permute, const TPZFMatrix &input, TPZFMatrix &output, int first, int last) 
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
void TPZDohrSubstructCondense::PermuteGather(const TPZVec<int> &permute, const TPZFMatrix &input, TPZFMatrix &output, int first, int last)
{
	int i,j,ncol = input.Cols();
	for(i=first; i<last; i++) for(j=0; j<ncol; j++)
	{
		output(i-first,j) = input.GetVal(permute[i],j);
	}
}

static TPZVec<int> dummyvec;
const TPZVec<int> &TPZDohrSubstructCondense::GatherVec(ENumbering origin, ENumbering destination) const
{
	std::map<std::pair<ENumbering,ENumbering>, TPZVec<int> >::const_iterator it;
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

const TPZVec<int> &TPZDohrSubstructCondense::ScatterVec(ENumbering origin, ENumbering destination) const
{
	std::map<std::pair<ENumbering,ENumbering>, TPZVec<int> >::const_iterator it;
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
