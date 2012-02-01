/**
 * @file
 * @brief Contains the implementation of the TPZElementMatrix methods.
 */
//$Id: pzelmat.cpp,v 1.13 2011-05-30 20:27:35 denise Exp $

#include "pzelmat.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzelmat"));
#endif

using namespace std;

void TPZElementMatrix::SetMatrixSize(short NumBli, short NumBlj,
									 short BlSizei, short BlSizej) {
	
	
	if(fMat.Rows() != NumBli*BlSizei || fMat.Cols() != NumBlj*BlSizej) {
		fMat.Redim(NumBli*BlSizei,NumBlj*BlSizej);
	}
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
    fSourceIndex(cp.fSourceIndex), fNumStateVars(cp.fNumStateVars)
{
    
}


void TPZElementMatrix::Print(std::ostream &out){
	fMat.Print("Unconstrained matrix",out,EMathematicaInput);
	int ncon = fConnect.NElements();
	int ic;
	for(ic=0; ic<ncon; ic++) {
		out << "Connect index " << fConnect[ic] << endl;
		this->fMesh->ConnectVec()[fConnect[ic]].Print(*fMesh,out);
	}
	fConstrMat.Print("Constrained matrix",out);
	ncon = fConstrConnect.NElements();
	for(ic=0; ic<ncon; ic++) {
		out << "Connect index " << fConstrConnect[ic] << endl;
		this->fMesh->ConnectVec()[fConstrConnect[ic]].Print(*fMesh,out);
	}
	if(fType == EK)
	{
		ComputeDestinationIndices();
		bool hasdepend = HasDependency();
		int size = fSourceIndex.NElements();
		TPZFMatrix constrmatrix(size,size,0.);
		int in,jn;
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
		std::stringstream sout;
		out << "Matrix size " << constrmatrix.Rows() << "\n";
		sout << "ConstrainedMatrix = ";
		constrmatrix.Print(sout.str().c_str(), out, EMathematicaInput);
	}
}

void TPZElementMatrix::ComputeDestinationIndices(){
    if (!this->HasDependency()){
		this->fSourceIndex.Resize(this->fMat.Rows());
		this->fDestinationIndex.Resize(this->fMat.Rows());
		int destindex = 0;
        int fullmatindex = 0;
		const int numnod = this->NConnects();
		for(int in = 0; in < numnod; in++){
			const int npindex = this->ConnectIndex(in);
			TPZConnect &np = this->fMesh->ConnectVec()[npindex];
			int blocknumber = np.SequenceNumber();
			int firsteq = this->fMesh->Block().Position(blocknumber);
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
		{
			std::stringstream sout;
			sout<<" fSourceIndex " <<fSourceIndex<< " fDestinationIndex "<<fDestinationIndex<<std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	}//if
	else{
    
        int destindex = 0;
        int fullmatindex = 0;
        this->fDestinationIndex.Resize(this->fConstrMat.Rows());
        this->fSourceIndex.Resize(this->fConstrMat.Rows());
        int numnod = this->fConstrConnect.NElements();
        for(int in = 0; in < numnod; in++){
            const int npindex = this->fConstrConnect[in];
            TPZConnect &np = this->fMesh->ConnectVec()[npindex];
            int blocknumber = np.SequenceNumber();
            int firsteq = this->fMesh->Block().Position(blocknumber);
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
	
	if (!this->fNumStateVars && this->fMat.Rows() != 0){
		LOGPZ_FATAL(logger, "this->fNumStateVars not initialized");
	}
	
	int totalnodes= this->NConnects();
	this->fConstrConnect.Resize(totalnodes);
	if(totalnodes) this->fConstrConnect.Fill(0,0);
	int in;
	for(in=0; in<totalnodes; in++) this->fConstrConnect[in] = this->fConnect[in];
	// total number of nodes of the constrained element
	TPZConnect::BuildConnectList(this->fConstrConnect, this->fConnect, *this->fMesh);
	totalnodes = this->fConstrConnect.NElements();
	
	// compute the list of nodes and their proper order of processing
	TPZVec<int> DependenceOrder(0);
	// this->fConstrNod, totalnodes and DependenceOrder
	// are initialized using codes documented above
	TPZConnect::BuildDependencyOrder(this->fConstrConnect,DependenceOrder,*this->fMesh);
	
	// compute the number of statevariables
	// the number of state variables is the number of unknowns associated with
	// each shapefunction
	// numstate is best initialized during computation of the stiffness matrix
	//   TPZAutoPointer<TPZMaterial> mat = Material();
	//   int numstate = mat->NStateVariables();
	
	// initialize the block structure
	this->fConstrBlock.SetNBlocks(totalnodes);
	
	// toteq contains the total number of equations of the constrained matrix
	int toteq = 0;
	for(in=0; in<totalnodes; in++) {
		int dfnindex = this->fConstrConnect[in];
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
	
	int nrhs = this->fMat.Cols();
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
		int idfn = this->fConnect[in];
		// find the index of the node in the destination (constrained) matrix
		while(irnode < totalnodes && this->fConstrConnect[irnode] != idfn) irnode++;
		
		// first and last rows in the original matrix
		int ifirst = this->fBlock.Position(in);
		int ilast = ifirst+this->fBlock.Size(in);
		
		// first and last rows in the desination (reception) matrix
		int irfirst = this->fConstrBlock.Position(irnode);
		//	   int irlast = irfirst+this->fConstrBlock->Size(irnode);
		
		int i,ir,ieq;
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
				int jdfn = this->fConnect[jn];
				// find the index of the node in the destination (constrained) matrix
				while(jrnode < totalnodes && this->fConstrConnect[jrnode] != jdfn) jrnode++;
				if(jrnode == totalnodes) {
					LOGPZ_WARN(logger, "node not found in node list");
				}
				// first and last columns in the original matrix
				int jfirst = this->fBlock.Position(jn);
				int jlast = jfirst+this->fBlock.Size(jn);
				// first and last columns in the desination (reception) matrix
				int jrfirst = this->fConstrBlock.Position(jrnode);
				//int jrlast = irfirst+this->fConstrBlock->Size(jrnode);
				int j,jr;
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
			int dfnindex = this->fConstrConnect[in];
			TPZConnect *dfn = &(fMesh->ConnectVec()[dfnindex]);
			if(DependenceOrder[in] != current_order) continue;
			
			// only nodes which have dependency order equal to the
			// current order are processed
			numnodes_processed++;
			
			int inpos = this->fConstrBlock.Position(in);
			int insize = this->fConstrBlock.Size(in);
			// inpos : position of the dependent equation
			// insize : number of equations processed
			
			// loop over the nodes from which dfn depends
			TPZConnect::TPZDepend *dep = dfn->FirstDepend();
			while(dep) {
				int depnodeindex = dep->fDepConnectIndex;
				// look for the index where depnode is found
				int depindex=0;
				while(depindex < totalnodes && this->fConstrConnect[depindex] != depnodeindex) depindex++;
				if(depindex == totalnodes) {
					LOGPZ_WARN(logger,"node not found in node list");
				}
				
				int deppos = this->fConstrBlock.Position(depindex);
				int depsize = this->fConstrBlock.Size(depindex);
				// deppos : position of the receiving equation
				// depsize : number of receiving equations
				
				// process the rows of the constrained matrix
				int send;
				int receive;
				int ieq;
				REAL coef;
				int idf;
				int numstate = dfn->NState();//this->fNumStateVars;
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
		} // end of loop over all nodes
		current_order++;
	} // end of while loop
}//void

bool TPZElementMatrix::HasDependency(){
	int nconnects = this->NConnects();
	int in, index;
	for(in=0; in<nconnects; in++){
		index = this->ConnectIndex(in);
		if(this->fMesh->ConnectVec()[index].HasDependency()){
			return true;
		}
	}
	return false;
}

/** @brief permute the order of the connects
 */
void TPZElementMatrix::PermuteGather(TPZVec<int> &permute)
{
    if (permute.size() != fConnect.size()) {
        DebugStop();
    }
    TPZElementMatrix cp(*this);
    for (int i=0; i<fConnect.size(); ++i) {
        fConnect[i] = cp.fConnect[permute[i]];
        fBlock.Set(i, cp.fBlock.Size(permute[i]));
    }
    fBlock.Resequence();
    if (fType == EK) {
        int ibl,jbl;
        for (ibl=0; ibl<fBlock.NBlocks(); ++ibl) {
            int iblsize = fBlock.Size(ibl);
            for (jbl=0; jbl<fBlock.NBlocks(); ++jbl) {
                int jblsize = fBlock.Size(jbl);
                for (int idf=0; idf<iblsize; ++idf) {
                    for (int jdf=0; jdf<jblsize; ++jdf) {
                        fBlock(ibl,jbl,idf,jdf) = cp.fBlock(permute[ibl],permute[jbl],idf,jdf);
                    }
                }
            }
        }
    }
    else if (fType == EF)
    {
        int ibl;
        for (ibl=0; ibl<fBlock.NBlocks(); ++ibl) {
            int iblsize = fBlock.Size(ibl);
            int jblsize = fMat.Cols();
            for (int idf=0; idf<iblsize; ++idf) {
                for (int jdf=0; jdf<jblsize; ++jdf) {
                    fBlock(ibl,0,idf,jdf) = cp.fBlock(permute[ibl],0,idf,jdf);
                }
            }
        }
    }
}

/*
 bandmat & bandmat::operator+=(elmat & ek) {
 
 for (int i = 0; i<ek.numnod(); i++) {
 for (int j = 0; j<ek.numnod(); j++) {
 addsub( (ek.node(i))->eq_number(), 
 (ek.node(j))->eq_number(), ek.mat->extract(i,j) );
 }
 }
 return *this;
 }
 
 void matrix::operator+=(elmat & ek) {
 
 if(ek.mat->colblocks() == ek.mat->rowblocks() ) {
 for (int i = 0; i<ek.numnod(); i++) {
 for (int j = 0; j<ek.numnod(); j++) {
 addsub( (ek.node(i))->eq_number(), 
 (ek.node(j))->eq_number(), ek.mat->extract(i,j) );
 }
 }
 } else if (ek.mat->colblocks() == 1) {
 for (int i = 0; i<ek.numnod(); i++) {
 addsub( (ek.node(i))->eq_number(), 
 0, ek.mat->extract(i,0) );
 }
 } else {
 pzerror << "matrix.+= doesn t know how to handle the\n"
 "assembly process\n";
 pzerror.show();
 }
 return;
 }
 */

