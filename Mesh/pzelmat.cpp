/**
 * @file
 * @brief Contains the implementation of the TPZElementMatrix methods.
 */

#include "pzelmat.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"

#include "pzlog.h"

#ifdef LOG4CXX
static PZLogger logger("pz.mesh.tpzelmat");
#endif

using namespace std;

void TPZElementMatrix::SetMatrixSize(short NumBli, short NumBlj,
									 short BlSizei, short BlSizej) {
	
	
	if(fMat.Rows() != NumBli*BlSizei || fMat.Cols() != NumBlj*BlSizej) {
		fMat.Redim(NumBli*BlSizei,NumBlj*BlSizej);
	}
}
TPZElementMatrix &TPZElementMatrix::operator=(const TPZElementMatrix &cp)
{
    fType = cp.fType;
    fMesh = cp.fMesh;
    fConnect = cp.fConnect;
    fMat = cp.fMat;
    fBlock = cp.fBlock;
    fConstrConnect = cp.fConstrConnect;
    fConstrMat = cp.fConstrMat;
    fConstrBlock = cp.fConstrBlock;
    fDestinationIndex = cp.fDestinationIndex;
    fSourceIndex = cp.fSourceIndex;
    fBlock.SetMatrix(&fMat);
    fConstrBlock.SetMatrix(&fConstrMat);
    return *this;
}

void TPZElementMatrix::SetMatrixMinSize(short NumBli, short NumBlj, 
										short BlSizei, short BlSizej) {
	
	
	if(fMat.Rows() < NumBli*BlSizei || fMat.Cols() < NumBlj*BlSizej) {
		fMat.Redim(NumBli*BlSizei,NumBlj*BlSizej);
	}
}

TPZElementMatrix::TPZElementMatrix(const TPZElementMatrix &cp) : 
    fType(cp.fType), fMesh(cp.fMesh), fConnect(cp.fConnect), fMat(cp.fMat),
    fBlock(cp.fBlock), fConstrConnect(cp.fConstrConnect), fConstrMat(cp.fConstrMat),
    fConstrBlock(cp.fConstrBlock), fDestinationIndex(cp.fDestinationIndex),
    fSourceIndex(cp.fSourceIndex)
{
    fBlock.SetMatrix(&fMat);
    fConstrBlock.SetMatrix(&fConstrMat);
}

void TPZElementMatrix::Print(std::ostream &out){
	if(fType == EK)
	{
        int ncon = fConnect.NElements();
        int ic;
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
		TPZFNMatrix<400,STATE> constrmatrix(size,size,0.);
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
		TPZFMatrix<STATE> constrmatrix(size,fMat.Cols(),0.);
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

void TPZElementMatrix::ComputeDestinationIndices(){
    if (!this->HasDependency()){
		this->fSourceIndex.Resize(this->fMat.Rows());
		this->fDestinationIndex.Resize(this->fMat.Rows());
		int64_t destindex = 0L;
        int64_t fullmatindex = 0L;
		const int numnod = this->NConnects();
		for(int in = 0; in < numnod; in++){
			const int64_t npindex = this->ConnectIndex(in);
			TPZConnect &np = this->fMesh->ConnectVec()[npindex];
			int64_t blocknumber = np.SequenceNumber();
			int64_t firsteq = this->fMesh->Block().Position(blocknumber);
			int ndf = this->fMesh->Block().Size(blocknumber);
            if(np.HasDependency() || np.IsCondensed()) {
                fullmatindex += ndf;
                continue;
            }//for (np)
			for(int idf=0; idf<ndf; idf++){
				this->fSourceIndex[destindex] = fullmatindex++;
				this->fDestinationIndex[destindex++] = firsteq+idf;
			}//for idf
		}//for in
        this->fSourceIndex.Resize(destindex);
        this->fDestinationIndex.Resize(destindex);		
#ifdef LOG4CXX
        if (logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout<<"fSourceIndex " <<fSourceIndex<< "\nfDestinationIndex "<<fDestinationIndex<<std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	}//if
	else{
    
        int64_t destindex = 0L;
        int64_t fullmatindex = 0L;
        this->fDestinationIndex.Resize(this->fConstrMat.Rows());
        this->fSourceIndex.Resize(this->fConstrMat.Rows());
        int numnod = this->fConstrConnect.NElements();
        for(int in = 0; in < numnod; in++){
            const int64_t npindex = this->fConstrConnect[in];
            TPZConnect &np = this->fMesh->ConnectVec()[npindex];
            int64_t blocknumber = np.SequenceNumber();
            int64_t firsteq = this->fMesh->Block().Position(blocknumber);
            int ndf = this->fMesh->Block().Size(blocknumber);
            if(np.HasDependency() || np.IsCondensed()) {
                fullmatindex += ndf;
                continue;
            }//for (np)
            for(int idf=0; idf<ndf; idf++) {
                this->fSourceIndex[destindex] = fullmatindex++;
                this->fDestinationIndex[destindex++] = firsteq+idf;
            }//for idf
        }//for in
        this->fSourceIndex.Resize(destindex);
        this->fDestinationIndex.Resize(destindex);
    }
}//void

void TPZElementMatrix::ApplyConstraints(){
	
	if (this->fMat.Rows() == 0){
		LOGPZ_FATAL(logger, "this->fMat not initialized");
	}
	
	int totalnodes= this->NConnects();
	this->fConstrConnect.Resize(totalnodes);
	if(totalnodes) this->fConstrConnect.Fill(0,0);
	int in;
    std::set<int64_t> origlist,connectlist;
	for(in=0; in<totalnodes; in++) connectlist.insert(this->fConnect[in]);
    for (std::list<TPZOneShapeRestraint>::iterator it = fOneRestraints.begin(); it != fOneRestraints.end(); it++) {
        for (int c=0; c< it->fFaces.size(); c++) {
            connectlist.insert(it->fFaces[c].first);
        }
    }
    origlist = connectlist;
	// total number of nodes of the constrained element
	TPZConnect::BuildConnectList(connectlist, origlist, *this->fMesh);
    this->fConstrConnect.resize(connectlist.size());
    std::set<int64_t>::iterator it = connectlist.begin();
    for (int64_t i=0; i<connectlist.size(); i++) {
        fConstrConnect[i] = *it;
        it++;
    }
	totalnodes = this->fConstrConnect.NElements();
	
	// compute the list of nodes and their proper order of processing
	TPZVec<int> DependenceOrder;
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
	for(in=0; in<totalnodes; in++) {
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
	int numnod = this->fConnect.NElements();
	for(in=0; in<numnod; in++) {
		int irnode =0;
		int64_t idfn = this->fConnect[in];
		// find the index of the node in the destination (constrained) matrix
		while(irnode < totalnodes && this->fConstrConnect[irnode] != idfn) irnode++;
		
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
					(this->fConstrMat)(ir,ieq) = (this->fMat)(i,ieq);
				}
			}
		}
		else{
			int jn;
			for(jn=0; jn<numnod; jn++) {
				int jrnode = 0;
				int64_t jdfn = this->fConnect[jn];
				// find the index of the node in the destination (constrained) matrix
				while(jrnode < totalnodes && this->fConstrConnect[jrnode] != jdfn) jrnode++;
				if(jrnode == totalnodes) {
					LOGPZ_WARN(logger, "node not found in node list");
				}
				// first and last columns in the original matrix
				int64_t jfirst = this->fBlock.Position(jn);
				int64_t jlast = jfirst+this->fBlock.Size(jn);
				// first and last columns in the desination (reception) matrix
				int64_t jrfirst = this->fConstrBlock.Position(jrnode);
				//int jrlast = irfirst+this->fConstrBlock->Size(jrnode);
				int64_t j,jr;
				for(i=ifirst,ir=irfirst;i<ilast; i++,ir++) {
					for(j=jfirst,jr=jrfirst;j<jlast; j++,jr++) {
						(this->fConstrMat)(ir,jr) = (this->fMat)(i,j);
					}
				}
			}
		}//else
	}
	
	int numnodes_processed = 0;
	int current_order = 0;
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
			TPZConnect::TPZDepend *dep = dfn->FirstDepend();
			while(dep) {
				int64_t depnodeindex = dep->fDepConnectIndex;
				// look for the index where depnode is found
				int depindex=0;
				while(depindex < totalnodes && this->fConstrConnect[depindex] != depnodeindex) depindex++;
				if(depindex == totalnodes) {
					LOGPZ_WARN(logger,"node not found in node list");
				}
				
				int64_t deppos = this->fConstrBlock.Position(depindex);
				int64_t depsize = this->fConstrBlock.Size(depindex);
				// deppos : position of the receiving equation
				// depsize : number of receiving equations
				
				// process the rows of the constrained matrix
				int64_t send;
				int64_t receive;
				int ieq;
				STATE coef;
				int idf;
				int numstate = dfn->NState();
				for(send=inpos; send<inpos+insize; send += numstate) {
					for(receive=deppos; receive<deppos+depsize; receive += numstate) {
						coef = dep->fDepMatrix((send-inpos)/numstate,(receive-deppos)/numstate);
						if (this->fType == TPZElementMatrix::EK){
							for(ieq=0; ieq<toteq; ieq++) for(idf=0; idf<numstate; idf++)  {
								(this->fConstrMat)(receive+idf,ieq) += coef*(this->fConstrMat)(send+idf,ieq);
							}
						}//EK
						else{
							for(ieq=0; ieq<nrhs; ieq++) for(idf=0; idf<numstate; idf++) {
								(this->fConstrMat)(receive+idf,ieq) += coef*(this->fConstrMat)(send+idf,ieq);
							}
						}//EF
					}
				}
				
				if (this->fType == TPZElementMatrix::EK){
					for(send=inpos; send<inpos+insize; send += numstate) {
						for(receive=deppos; receive<deppos+depsize; receive += numstate) {
							coef = dep->fDepMatrix((send-inpos)/numstate,(receive-deppos)/numstate);
							for(ieq=0; ieq<toteq; ieq++) for(idf=0; idf<numstate; idf++) {
								(this->fConstrMat)(ieq,receive+idf) += coef*(this->fConstrMat)(ieq,send+idf);
							}
						}
					}
				}//EK
				
				dep = dep->fNext;
			} // end of while
            
            /// check whether the connect has a one shape restraint
            if (fOneRestraints.size())
            {
                ApplyOneShapeConstraints(in);
            }
            
            
		} // end of loop over all nodes
		current_order++;
	} // end of while loop
}//void

/// Apply the constraint of the one shape restraints
void TPZElementMatrix::ApplyOneShapeConstraints(int constraintindex)
{
    int64_t dfnindex = this->fConstrConnect[constraintindex];


#ifdef LOG4CXX
    int count = 0;
    for (std::list<TPZOneShapeRestraint>::iterator it = fOneRestraints.begin(); it != fOneRestraints.end(); it++) {
        if (it->fFaces[0].first != dfnindex) {
            continue;
        }
        count++;
    }
    if (count && logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Element matrix before ApplyOneShapeConstraint\n";
        fConstrMat.Print("EKBefore = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int64_t inpos = this->fConstrBlock.Position(constraintindex);
    int64_t toteq = this->fConstrMat.Rows();
    int64_t nrhs = this->fConstrMat.Cols();

    for (std::list<TPZOneShapeRestraint>::iterator it = fOneRestraints.begin(); it != fOneRestraints.end(); it++) {
        if (it->fFaces[0].first != dfnindex) {
            continue;
        }
        int64_t send = inpos+it->fFaces[0].second;
        for (int id=1; id<4; id++) {
            int64_t depindex = it->fFaces[id].first;
            int locdep = 0;
            for (locdep = 0; locdep < fConstrConnect.size(); locdep++) {
                if (fConstrConnect[locdep] == depindex) {
                    break;
                }
            }
            if (locdep == fConstrConnect.size()) {
                DebugStop();
            }
            int64_t deppos = this->fConstrBlock.Position(locdep);
            int64_t receive = deppos+it->fFaces[id].second;
            REAL coef = -it->fOrient[id]/it->fOrient[0];
            if (this->fType == TPZElementMatrix::EK){
                for(int ieq=0; ieq<toteq; ieq++) {
                    (this->fConstrMat)(receive,ieq) += coef*(this->fConstrMat)(send,ieq);
                }
            }//EK
            else
            {
                
                for(int ieq=0; ieq<nrhs; ieq++) {
                    (this->fConstrMat)(receive,ieq) += coef*(this->fConstrMat)(send,ieq);
                }
            }//EF

            if (this->fType == TPZElementMatrix::EK){
                for(int ieq=0; ieq<toteq; ieq++)
                {
                    (this->fConstrMat)(ieq,receive) += coef*(this->fConstrMat)(ieq,send);
                }
            }//EK

        }
        if (this->fType == TPZElementMatrix::EK){
            for(int ieq=0; ieq<toteq; ieq++)
            {
                (this->fConstrMat)(ieq,send) = 0.;
                (this->fConstrMat)(send,ieq) = 0.;
            }
            (this->fConstrMat)(send,send) = 1.;
        }//EK
        else
        {
            for(int ieq=0; ieq<nrhs; ieq++) {
                (this->fConstrMat)(send,ieq) = 0.;
            }
        }

    }
#ifdef LOG4CXX
    if (count && logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Element matrix after ApplyOneShapeConstraint\n";
        fConstrMat.Print("EKAfter = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}



bool TPZElementMatrix::HasDependency()
{
    if (fOneRestraints.size()) {
        return true;
    }
	int in, nconnects = this->NConnects();
	int64_t index;
	for(in=0; in<nconnects; in++)
    {
		index = this->ConnectIndex(in);
       // bool val =this->fMesh->ConnectVec()[index].HasDependency();
		if(this->fMesh->ConnectVec()[index].HasDependency())
        {
			return true;
		}
	}
	return false;
}

/** @brief permute the order of the connects */
void TPZElementMatrix::PermuteGather(TPZVec<int64_t> &permute)
{
    if (permute.size() != fConnect.size()) {
        DebugStop();
    }
    TPZElementMatrix cp(*this);
    for (int64_t i=0; i<fConnect.size(); ++i) {
        fConnect[i] = cp.fConnect[permute[i]];
        fBlock.Set(i, cp.fBlock.Size(permute[i]));
    }
    fBlock.Resequence();
#ifdef LOG4CXX2
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
                        fBlock(ibl,jbl,idf,jdf) = cp.fBlock(permute[ibl],permute[jbl],idf,jdf);
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
                for (int64_t jdf=0; jdf<jblsize; ++jdf) {
                    fBlock(ibl,0,idf,jdf) = cp.fBlock(permute[ibl],0,idf,jdf);
                }
            }
        }
    }
}


void TPZElementMatrix::BuildDependencyOrder(TPZVec<int64_t> &connectlist, TPZVec<int> &DependenceOrder, TPZCompMesh &mesh) {
    // nodelist (input) : vector which contains pointers to all nodes which
    // are in the dependency chain of the nodes of the element
    int64_t totalnodes = connectlist.NElements();
    DependenceOrder.Resize(totalnodes);
    DependenceOrder.Fill(0,0);
    // initialize the vector which contains the
    // dependency order to zero
    int CurrentOrder = 0;
    // order which is currently processed
    int64_t numnodes_processed = totalnodes;

    for (std::list<TPZOneShapeRestraint>::iterator it = fOneRestraints.begin(); it != fOneRestraints.end(); it++) {
        for (int i=1; i<4; i++)
        {
            int64_t index = it->fFaces[1].first;
            TPZConnect &dfn = mesh.ConnectVec()[index];
            dfn.SetDependenceOrder(index,mesh,1,connectlist,DependenceOrder);
        }
    }

    
    // number of nodes processed during the current cycle
    while(numnodes_processed) {
        
        numnodes_processed = 0;
        int64_t i;
        for(i=0; i<totalnodes; i++) {
            int64_t dfnindex = connectlist[i];
            TPZConnect &dfn = mesh.ConnectVec()[dfnindex];
            if(dfn.HasDependency() && DependenceOrder[i] == CurrentOrder) {
                dfn.SetDependenceOrder(dfnindex,mesh,CurrentOrder,connectlist,DependenceOrder);
                // this method will fill in the DependenceOrder vector by recursively
                // calling SetDependenceOrder for the nodes upon which dfn depends
                numnodes_processed++;
            }
        }
        // force the loop to process the order one connects
        if (fOneRestraints.size() && CurrentOrder == 0) {
            numnodes_processed++;
        }
        CurrentOrder++;
    }
}
