#include "TPZElementMatrixT.h"
#include "pzcmesh.h"
#include "TPZMatrixWindow.h"
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzelmat");
#else
static int logger;
#endif


template<class TVar>
TPZElementMatrixT<TVar>::TPZElementMatrixT(const TPZElementMatrixT<TVar> &cp) :
    TPZElementMatrix(cp), fMat(cp.fMat), fBlock(cp.fBlock),
    fConstrMat(cp.fConstrMat),fConstrBlock(cp.fConstrBlock)
{
    fBlock.SetMatrix(&fMat);
    fConstrBlock.SetMatrix(&fConstrMat);
}


template<class TVar>
TPZElementMatrixT<TVar> &TPZElementMatrixT<TVar>::operator=(const TPZElementMatrixT<TVar> &cp)
{
    TPZElementMatrix::operator=(cp);
    fMat = cp.fMat;
    fBlock = cp.fBlock;
    fConstrMat = cp.fConstrMat;
    fConstrBlock = cp.fConstrBlock;
    fBlock.SetMatrix(&fMat);
    fConstrBlock.SetMatrix(&fConstrMat);
    return *this;
}

template<class TVar>
void TPZElementMatrixT<TVar>::Print(std::ostream &out){
	if(fType == EK)
	{
        int ncon = fConnect.NElements();
        int64_t ic;
        out << "Connect vector\n";
        for(ic=0; ic<ncon; ic++) {
            //	out << "Connect index " << fConnect[ic] << endl;
            out << "ic = " << ic << " index " << fConnect[ic] << " ";
            this->fMesh->ConnectVec()[fConnect[ic]].Print(*fMesh,out);
        }
        out << "Constrained connect vector\n";
        //fConstrMat.Print("Constrained matrix",out);
        ncon = fConstrConnect.NElements();
        for(ic=0; ic<ncon; ic++) {
            //	out << "Connect index " << fConstrConnect[ic] << endl;
            out << "ic = " << ic << " index " << fConstrConnect[ic]  << ' ';
            this->fMesh->ConnectVec()[fConstrConnect[ic]].Print(*fMesh,out);
        }
		ComputeDestinationIndices();
		bool hasdepend = HasDependency();
		int64_t size = fSourceIndex.NElements();
		//TPZFMatrix<REAL> constrmatrix(size,size,0.);
		TPZFNMatrix<400,TVar> constrmatrix(size,size,0.);
		int64_t in,jn;
		for(in=0; in<size; in++)
		{
			for (jn=0; jn<size; jn++) {
				if(hasdepend)
				{
					constrmatrix(in,jn) = fConstrMat(fSourceIndex[in],fSourceIndex[jn]);
				}
				else {
					constrmatrix(in,jn) = fMat(fSourceIndex[in],fSourceIndex[jn]);
				}
				
			}
		}
        out << "SourceIndex = " << fSourceIndex << std::endl;
        fConstrMat.Print("EKOrig=",out,EMathematicaInput);
		std::stringstream sout;
        sout << "EK = ";
		//out << "Matrix size " << constrmatrix.Rows() << "\n";
		//sout << "ConstrainedMatrix = ";
		constrmatrix.Print(sout.str().c_str(), out, EMathematicaInput);
	}
    else if(fType == EF)
	{
		ComputeDestinationIndices();
		bool hasdepend = HasDependency();
		int64_t size = fSourceIndex.NElements();
		//TPZFMatrix<REAL> constrmatrix(size,size,0.);
		TPZFMatrix<TVar> constrmatrix(size,fMat.Cols(),0.);
		int64_t in,jn;
		for(in=0; in<size; in++)
		{
			for (jn=0; jn<fMat.Cols(); jn++) {
				if(hasdepend)
				{
					constrmatrix(in,jn) = fConstrMat(fSourceIndex[in],jn);
				}
				else {
					constrmatrix(in,jn) = fMat(fSourceIndex[in],jn);
				}
				
			}
		}
		std::stringstream sout;
        sout << "EF = ";
		//out << "Matrix size " << constrmatrix.Rows() << "\n";
		//sout << "ConstrainedMatrix = ";
		constrmatrix.Print(sout.str().c_str(), out, EMathematicaInput);
	}
    else
    {
        DebugStop();
    }
}

template<class TVar>
void TPZElementMatrixT<TVar>::SetMatrixSize(short NumBli, short NumBlj,
									 short BlSizei, short BlSizej) {
	
	
	if(fMat.Rows() != NumBli*BlSizei || fMat.Cols() != NumBlj*BlSizej) {
		fMat.Redim(NumBli*BlSizei,NumBlj*BlSizej);
	}
}
template<class TVar>
void TPZElementMatrixT<TVar>::SetMatrixMinSize(short NumBli, short NumBlj, 
										short BlSizei, short BlSizej) {
	
	
	if(fMat.Rows() < NumBli*BlSizei || fMat.Cols() < NumBlj*BlSizej) {
		fMat.Redim(NumBli*BlSizei,NumBlj*BlSizej);
	}
}

template<class TVar>
void TPZElementMatrixT<TVar>::ApplyConstraints(){
    
    if (this->fMat.Rows() == 0){
        LOGPZ_FATAL(logger, "this->fMat not initialized");
    }
    
    int64_t totalnodes= this->NConnects();
    this->fConstrConnect.Resize(totalnodes);
    if(totalnodes) this->fConstrConnect.Fill(0,0);

    std::set<int64_t> origlist,connectlist;
    for(int64_t in=0; in<totalnodes; in++) connectlist.insert(this->fConnect[in]);
    origlist = connectlist;
    // total number of nodes of the constrained element
    TPZConnect::BuildConnectList(connectlist, origlist, *this->fMesh);
    this->fConstrConnect.resize(connectlist.size());
    std::set<int64_t>::iterator it = connectlist.begin();
    for (int64_t i=0; i<connectlist.size(); i++) {
        fConstrConnect[i] = *it;
        it++;
    }
    std::map<int64_t,int> TargetConnectIndex;
    totalnodes = this->fConstrConnect.NElements();
    for(int64_t in = 0; in < totalnodes; in++) {
        TargetConnectIndex[fConstrConnect[in]] = in;
    }
    
    // compute the list of nodes and their proper order of processing
    TPZManVector<int,400> DependenceOrder;
    // this->fConstrNod, totalnodes and DependenceOrder
    // are initialized using codes documented above
    BuildDependencyOrder(this->fConstrConnect,DependenceOrder,*this->fMesh);
    
    // compute the number of statevariables
    // the number of state variables is the number of unknowns associated with
    // each shapefunction
    // numstate is best initialized during computation of the stiffness matrix
    //   TPZMaterial * mat = Material();
    //   int numstate = mat->NStateVariables();
    
    // initialize the block structure
    this->fConstrBlock.SetNBlocks(totalnodes);
    
    // toteq contains the total number of equations of the constrained matrix
    int64_t toteq = 0;
    for(int64_t in=0; in<totalnodes; in++) {
        int64_t dfnindex = this->fConstrConnect[in];
        TPZConnect &dfn = fMesh->ConnectVec()[dfnindex];
        int ndf = dfn.NDof(*fMesh);
        int ndfcheck = dfn.NState()*dfn.NShape();
        if(ndf != ndfcheck)
        {
            DebugStop();
        }
        this->fConstrBlock.Set(in,ndf);
        toteq += ndf;
    }
    
    this->fConstrBlock.Resequence();
    this->fConstrBlock.SetMatrix(&this->fConstrMat);
    
    int64_t nrhs = this->fMat.Cols();
    if (this->fType == TPZElementMatrix::EK){
        this->fConstrMat.Redim(toteq,toteq);
    }
    else{
        this->fConstrMat.Redim(toteq,nrhs);
    }
    
    // copy the original matrix to the constrained matrix
    int64_t numnod = this->fConnect.NElements();
    for(int64_t in=0; in<numnod; in++) {
        int irnode =0;
        int64_t idfn = this->fConnect[in];
        // find the index of the node in the destination (constrained) matrix
        // while(irnode < totalnodes && this->fConstrConnect[irnode] != idfn) irnode++;
#ifdef PZDEBUG
        if(TargetConnectIndex.find(idfn) == TargetConnectIndex.end()) {
            LOGPZ_ERROR(logger, "node not found in node list");
            DebugStop();
        }
#endif
        irnode = TargetConnectIndex[idfn];
        // first and last rows in the original matrix
        int64_t ifirst = this->fBlock.Position(in);
        int64_t ilast = ifirst+this->fBlock.Size(in);
        
        // first and last rows in the desination (reception) matrix
        int64_t irfirst = this->fConstrBlock.Position(irnode);
        //	   int irlast = irfirst+this->fConstrBlock->Size(irnode);
        
        int64_t i,ir,ieq;
        if (this->fType == TPZElementMatrix::EF){
            for(i=ifirst,ir=irfirst;i<ilast;i++,ir++) {
                for(ieq=0; ieq<nrhs; ieq++) {
                    (this->fConstrMat)(ir,ieq) += (this->fMat)(i,ieq);
                }
            }
        }
        else {
            int jn;
            for(jn=0; jn<numnod; jn++) {
                int jrnode = 0;
                int64_t jdfn = this->fConnect[jn];
                // find the index of the node in the destination (constrained) matrix
#ifdef PZDEBUG
                if(TargetConnectIndex.find(jdfn) == TargetConnectIndex.end()) {
                    LOGPZ_ERROR(logger, "node not found in node list");
                    DebugStop();
                }
#endif
                jrnode = TargetConnectIndex[jdfn];
                // first and last columns in the original matrix
                int64_t jfirst = this->fBlock.Position(jn);
                int64_t jlast = jfirst+this->fBlock.Size(jn);
                // first and last columns in the desination (reception) matrix
                int64_t jrfirst = this->fConstrBlock.Position(jrnode);
                //int jrlast = irfirst+this->fConstrBlock->Size(jrnode);
                int64_t j,jr;
                for(i=ifirst,ir=irfirst;i<ilast; i++,ir++) {
                    for(j=jfirst,jr=jrfirst;j<jlast; j++,jr++) {
                        (this->fConstrMat)(ir,jr) += (this->fMat)(i,j);
                    }
                }
            }
        }//else
    }
    
    int numnodes_processed = 0;
    int current_order = 0;
    TPZFNMatrix<150,TVar> localmat;
    
    TPZFMatrix<TVar> *fulldepmat = fUserAllocMat ? fUserAllocMat : &localmat;
    
    while(numnodes_processed < totalnodes) {
        int in;
        for(in=0; in<totalnodes; in++) {
            int64_t dfnindex = this->fConstrConnect[in];
            TPZConnect *dfn = &(fMesh->ConnectVec()[dfnindex]);
            if(DependenceOrder[in] != current_order) continue;
            
            // only nodes which have dependency order equal to the
            // current order are processed
            numnodes_processed++;
            
            int64_t inpos = this->fConstrBlock.Position(in);
            int64_t insize = this->fConstrBlock.Size(in);
            // inpos : position of the dependent equation
            // insize : number of equations processed
            
            // loop over the nodes from which dfn depends
            TPZConnect::TPZDependBase *dep = dfn->FirstDepend();
            
            
            while(dep) {
                auto dept = dynamic_cast<TPZConnect::TPZDepend<TVar>*>(dep);
                if(!dept){DebugStop();}
                const auto &depmat = dept->GetDepMatrix();
                
                int64_t depnodeindex = dep->fDepConnectIndex;
                // look for the index where depnode is found
                int depindex=0;
#ifdef PZDEBUG
                if(TargetConnectIndex.find(depnodeindex) == TargetConnectIndex.end()) {
                    LOGPZ_ERROR(logger, "node not found in node list");
                    DebugStop();
                }
#endif
                depindex = TargetConnectIndex[depnodeindex];
                
                int64_t deppos = this->fConstrBlock.Position(depindex);
                int64_t depsize = this->fConstrBlock.Size(depindex);
                // deppos : position of the receiving equation
                // depsize : number of receiving equations
                
                // process the rows of the constrained matrix
                int64_t send;
                int64_t receive;
                int ieq;
                TVar coef;
                int idf;
                const int numstate = dfn->NState();
                /*
                 now we fill fulldepmat, which takes into account numstate:
                 basically, each position of the original dep mat becomes
                 a diagonal block in the full dep mat
                 */

                    const int64_t orig_row=depmat.Rows();
                    const int64_t orig_col=depmat.Cols();
                    const int64_t full_row=numstate*orig_row;
                    const int64_t full_col=numstate*orig_col;
                    bool isShapeRestraint = (full_row == insize && full_col == depsize);
                    bool isAlgebraicRestraint = (orig_row == insize && orig_col == depsize);
                    if(!isShapeRestraint && ! isAlgebraicRestraint) DebugStop();
                    fulldepmat->Redim(full_row,full_col);
                    if(isShapeRestraint) {
                        for(int ic = 0; ic < orig_col; ic++){
                            const auto first_c = ic*numstate;
                            for(int ir = 0; ir < orig_row; ir++){
                                const auto first_r = ir*numstate;
                                const auto val = depmat.GetVal(ir,ic);
                                for(int istate = 0; istate < numstate; istate++){
                                    fulldepmat->PutVal(first_r+istate,first_c+istate,val);
                                }
                            }
                        }
                    } else {
                        *fulldepmat = depmat;
                    }
                    /*
                     now that we have dep mat, we need to compute the proper windows
                     such that we can compute
                     D^H K D,
                     where ^H stands for conjugate transpose.
                     
                     note: for the load vector, we compute only D^H F,
                     since there are no trial functions
                     */
                    
                    const int deprows = (int) fulldepmat->Rows();
                    const int depcols = (int) fulldepmat->Cols();
                    //the window is the full matrix
                    
                    {
                        TPZMatrixWindow<TVar> dep_window(*fulldepmat,0,0,deprows,depcols);
                        TPZMatrixWindow<TVar> send_window(this->fConstrMat,inpos,0,insize,(int)this->fConstrMat.Cols());
                        TPZMatrixWindow<TVar> receive_window(this->fConstrMat,deppos,0,depsize,(int)this->fConstrMat.Cols());
                        const TVar alpha{1};
                        const TVar beta{1};
                        //dep window is conjugate transpose
                        const int transp_a{2};
                        //we do not transpose send window
                        const int transp_x{0};
                        dep_window.MultAdd(send_window,receive_window,receive_window,alpha,beta,transp_a,transp_x);
                        
                    }
                    if (this->fType == TPZElementMatrix::EK){
                        //now we multiply it on the right side too
                        TPZMatrixWindow<TVar> dep_window(*fulldepmat,0,0,deprows,depcols);
                        TPZMatrixWindow<TVar> send_window(this->fConstrMat,0,inpos,this->fConstrMat.Rows(),insize);
                        TPZMatrixWindow<TVar> receive_window(this->fConstrMat,0,deppos,this->fConstrMat.Rows(),depsize);
                        const TVar alpha{1};
                        const TVar beta{1};
                        //we do not transpose send window
                        const int transp_a{0};
                        //we do not transpose dep window
                        const int transp_x{0};
                        send_window.MultAdd(dep_window,receive_window,receive_window,alpha,beta,transp_a,transp_x);
                    }
                    
                    dep = dep->fNext;
                } // end of while
                
                
            } // end of loop over all nodes
            current_order++;
        } // end of while loop
    }//void



template<class TVar>
void TPZElementMatrixT<TVar>::PermuteGather(TPZVec<int64_t> &permute)
{
    if (permute.size() != fConnect.size()) {
        DebugStop();
    }
    TPZElementMatrixT<TVar> cp(*this);
    for (int64_t i=0; i<fConnect.size(); ++i) {
        fConnect[i] = cp.fConnect[permute[i]];
        fBlock.Set(i, cp.fBlock.Size(permute[i]));
    }
    fBlock.Resequence();
#ifdef PZ_LOG2
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        cp.fBlock.Print("cp.fBlock ",sout);
        fBlock.Print("fBlock ",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    if (fType == EK) {
        int64_t ibl,jbl;
        for (ibl=0; ibl<fBlock.NBlocks(); ++ibl) {
            int64_t iblsize = fBlock.Size(ibl);
            for (jbl=0; jbl<fBlock.NBlocks(); ++jbl) {
                int64_t jblsize = fBlock.Size(jbl);
                for (int64_t idf=0; idf<iblsize; ++idf) {
                    for (int64_t jdf=0; jdf<jblsize; ++jdf) {
                      const auto [my_r,my_c] = fBlock.at(ibl,jbl,idf,jdf);
                      const auto [cp_r,cp_c] =
                        cp.fBlock.at(permute[ibl],permute[jbl],idf,jdf);
                      fMat.g(my_r,my_c) = cp.fMat.g(cp_r,cp_c);
                    }
                }
            }
        }
    }
    else if (fType == EF)
    {
        int64_t ibl;
        for (ibl=0; ibl<fBlock.NBlocks(); ++ibl) {
            int64_t iblsize = fBlock.Size(ibl);
            int64_t jblsize = fMat.Cols();
            for (int64_t idf=0; idf<iblsize; ++idf) {
                int64_t fromindex = fBlock.Index(ibl,idf);
                int64_t toindex = cp.fBlock.Index(permute[ibl],idf);
                for (int64_t jdf=0; jdf<jblsize; ++jdf) {
                    fMat(fromindex,jdf) = cp.fMat(toindex,jdf);
                }
            }
        }
    }
}

/// @brief Transfer the uncontrained equations from the constrained matrix to the uncontrained matrix
/// this method will change the configuration of the element matrix
template<class TVar>
void TPZElementMatrixT<TVar>::MakeUnconstrained() {
    if(!fMesh) DebugStop();
    if(!HasDependency()) return;
#ifdef PZDEBUG
    if(HasDependency() && fConstrMat.Rows() == 0) {
        // this means the method has been called without calling ApplyConstraints
        DebugStop();
    }
#endif
    auto &connectvec = fMesh->ConnectVec();
    std::set<int64_t> unconstrained;
    for(auto i : fConstrConnect) {
        TPZConnect &c = connectvec[i];
        if(!c.HasDependency()) {
            unconstrained.insert(i);
        }
    }
    TPZManVector<int64_t> connectindexes(unconstrained.size());
    int64_t count = 0;
    for(auto i : unconstrained) {
        connectindexes[count] = i;
        count++;
    }
    MakeUnconstrained(connectindexes);
}

/// @brief Transfer the uncontrained equations from the constrained matrix to the uncontrained matrix
/// this method will change the configuration of the element matrix
/// the resulting matrix will have the indicated order defined by connectindexes
template<class TVar>
void TPZElementMatrixT<TVar>::MakeUnconstrained(TPZVec<int64_t> &connectindexes) {
    if(!fMesh) DebugStop();
    if(!HasDependency()) {
        int64_t ncon = connectindexes.size();
        if(ncon != fConnect.size()) DebugStop();
        // reorder the connects defining a permutation
        std::map<int64_t,int64_t> connectposition;
        for(int64_t i=0; i<fConnect.size(); i++) connectposition[fConnect[i]] = i;
        TPZManVector<int64_t> permute(ncon);
        for(int64_t i = 0; i<ncon; i++) {
            int64_t cindex = connectindexes[i];
            if(connectposition.find(cindex) == connectposition.end()) DebugStop();
            permute[i] = connectposition[cindex];
        }
        PermuteGather(permute);
        return;
    }
    // at this point Hasdependency is true
#ifdef PZDEBUG
    if(fConstrMat.Rows() == 0) {
        // this means the method has been called without calling ApplyConstraints
        DebugStop();
    }
#endif
    auto &connectvec = fMesh->ConnectVec();
    std::set<int64_t> unconstrained;
    for(auto i : connectindexes) unconstrained.insert(i);
    for(auto i : fConstrConnect) {
        TPZConnect &c = connectvec[i];
        if(!c.HasDependency()) {
            if(unconstrained.find(i) == unconstrained.end()) {
                std::cout << "Inconsistent data structure\n";
                std::cout << "unconstrained connect index " << i << " not found\n";
                std::cout << "Constrained connects in the data structure ";
                for(auto z : fConstrConnect) std::cout << z << " ";
                std::cout << std::endl;
                std::cout << "Unconstrained connects ";
                for(auto z : unconstrained) std::cout << z << " ";
                std::cout << std::endl;
                std::cout << "connectindexes " << connectindexes << std::endl;
                std::cout << "fConnect " << this->fConnect << std::endl;
                std::cout << "fConstConnect " << this->fConstrConnect << std::endl;
                DebugStop();
            }
        }
    }
    // copy the data structure to the unconstrained matrix and block structure
    int64_t totalnodes = connectindexes.size();
    // initialize the block structure
    this->fBlock.SetNBlocks(totalnodes);
    this->fConnect.resize(totalnodes);
    
    // toteq contains the total number of equations of the constrained matrix
    int64_t toteq = 0;
    int count = 0;
    // map the connect index to the local vector
    std::map<int64_t,int> globaltolocal;
    for(auto cindex : connectindexes) {
        int64_t dfnindex = cindex;
        globaltolocal[cindex] = count;
        TPZConnect &dfn = connectvec[dfnindex];
        int ndf = dfn.NDof(*fMesh);
        int ndfcheck = dfn.NState()*dfn.NShape();
        if(ndf != ndfcheck)
        {
            DebugStop();
        }
        this->fBlock.Set(count,ndf);
        fConnect[count] = dfnindex;
        toteq += ndf;
        count++;
    }
    
    this->fBlock.Resequence();
    int64_t ncols = -1;
    if(fType == EK) {
        ncols = toteq;
    } else if(fType == EF) {
        ncols = fMat.Cols();
    } else {
        DebugStop();
    }
    fMat.Redim(toteq, ncols);
    for(int64_t i = 0; i<fConstrConnect.size(); i++) {
        int64_t constrindex = fConstrConnect[i];
        TPZConnect &c = connectvec[constrindex];
        if(c.HasDependency()) continue;
        if(globaltolocal.find(constrindex) == globaltolocal.end()) DebugStop();
        int ndof = fConstrBlock.Size(i);
        int64_t locindex = globaltolocal[constrindex];
        fBlock.Set(locindex,ndof);
    }
    fBlock.Resequence();
    if(fType == EK) {
        for(int64_t i = 0; i<fConstrConnect.size(); i++) {
            int64_t constrindex = fConstrConnect[i];
            TPZConnect &c = connectvec[constrindex];
            if(c.HasDependency()) continue;
            if(globaltolocal.find(constrindex) == globaltolocal.end()) DebugStop();
            int64_t ilocindex = globaltolocal[constrindex];
            int64_t ilocpos = fBlock.Position(ilocindex);
            int64_t ilocsize = fBlock.Size(ilocindex);
            int64_t iconstrpos = fConstrBlock.Position(i);
            if(ilocsize == 0) continue;
            for(int64_t j = 0; j<fConstrConnect.size(); j++) {
                int64_t constrindex = fConstrConnect[j];
                TPZConnect &c = connectvec[constrindex];
                if(c.HasDependency()) continue;
                if(globaltolocal.find(constrindex) == globaltolocal.end()) DebugStop();
                int64_t jlocindex = globaltolocal[constrindex];
                int64_t jlocpos = fBlock.Position(jlocindex);
                int64_t jlocsize = fBlock.Size(jlocindex);
                int64_t jconstrpos = fConstrBlock.Position(j);
                if(jlocsize == 0) continue;
                TPZMatrixWindow<TVar> MatW(fMat,ilocpos,jlocpos,ilocsize,jlocsize);
                TPZMatrixWindow<TVar> ConstrMatW(fConstrMat,iconstrpos,jconstrpos,ilocsize,jlocsize);
                MatW = ConstrMatW;

            }

        }
    } else if(fType == EF) {
        int64_t nrhs = fMat.Cols();
        for(int64_t i = 0; i<fConstrConnect.size(); i++) {
            int64_t constrindex = fConstrConnect[i];
            TPZConnect &c = connectvec[constrindex];
            if(c.HasDependency()) continue;
            if(globaltolocal.find(constrindex) == globaltolocal.end()) DebugStop();
            int64_t ilocindex = globaltolocal[constrindex];
            int64_t ilocpos = fBlock.Position(ilocindex);
            int64_t ilocsize = fBlock.Size(ilocindex);
            if(ilocsize == 0) continue;
            int64_t iconstrpos = fConstrBlock.Position(i);
            TPZMatrixWindow<TVar> MatW(fMat,ilocpos,0,ilocsize,nrhs);
            TPZMatrixWindow<TVar> ConstrMatW(fConstrMat,iconstrpos,0,ilocsize,nrhs);
            MatW = ConstrMatW;
        }
    }
    fConstrBlock.SetNBlocks(0);
    fConstrMat.Redim(0,0);
    fConstrConnect.resize(0);
#ifdef PZDEBUG
    if(HasDependency())
    {
        int in, nconnects = this->NConnects();
        int64_t index;
        for(in=0; in<nconnects; in++)
        {
            index = this->ConnectIndex(in);
           // bool val =this->fMesh->ConnectVec()[index].HasDependency();
            if(this->fMesh->ConnectVec()[index].HasDependency())
            {
                std::cout << "Connect index " << index << " has dependency\n";
            }
        }
        DebugStop();
    }
#endif

}




template class TPZElementMatrixT<STATE>;
template class TPZElementMatrixT<CSTATE>;
