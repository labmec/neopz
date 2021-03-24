/**
 * @file
 * @brief Contains the implementation of the TPZCheckMesh methods.
 */

#include "pzcheckmesh.h"
#include "pzgeoelside.h"
#include "pzstack.h"
#include "pzintel.h"
#include "TPZMaterial.h"

using namespace std;

TPZCheckMesh::TPZCheckMesh(TPZCompMesh *mesh, std::ostream *out) {
	fMesh = mesh;
	fOut = out;
	fNState = 1;   /** WARNING */
}

void TPZCheckMesh::BuildDependList(int connect, TPZStack<int> &dependlist) {
	int nconnects = fMesh->ConnectVec().NElements();
	int ic;
	for(ic=0; ic<nconnects; ic++) {
		TPZConnect::TPZDepend *dep = fMesh->ConnectVec()[ic].FirstDepend();
		if(dep && dep->HasDepend(connect)) {
			dependlist.Push(ic);
			continue;
		}
	}
}

TPZCompElSide TPZCheckMesh::FindElement(int connect) {
	int nelem = fMesh->ElementVec().NElements();
	int el;
	for(el=0; el<nelem; el++) {
		TPZCompEl *cel = fMesh->ElementVec()[el];
		if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZGeoEl *gel = cel->Reference();
        int ns = gel->NSides();
		
		for(int side=0; side<ns; side++) {
            int nsc = intel->NSideConnects(side);
            for (int ic=0; ic<nsc; ic++)
            {
                int64_t index = intel->SideConnectIndex(ic, side);
                if(index == connect) {
                    return TPZCompElSide(cel,side);
                }
			}
		}
	}
	TPZCompElSide nullside;
	return nullside;
}

int TPZCheckMesh::VerifyCompatibilityBetweenNShapesAndBlockSize(int connect) {
	int check = 1;
	TPZConnect &df = fMesh->ConnectVec()[connect];
	if(df.HasDependency() || df.IsCondensed() || !df.NElConnected() || df.SequenceNumber() == -1) return check;
	int dofsize = df.NShape()*df.NState();
	// check the consistency between the block size and the data structure of the connect
	{
		int seqnum = df.SequenceNumber();
		int blsize = fMesh->Block().Size(seqnum);
		if (blsize != dofsize) {
			check = 0;
		}
	}
	return check;
}

int TPZCheckMesh::VerifyConnect(int connect) {

    int check = 1;

    {
        TPZCompElSide elcon = FindElement(connect);
        TPZConnect &c = fMesh->ConnectVec()[connect];
        TPZCompElSide large = elcon.LowerLevelElementList(1);
        if (c.HasDependency() && !large) {
            *fOut << "Connect " << connect << " has dependency but no large element\n";
            check = 0;
        }
        if(c.NShape() && !c.HasDependency() && large)
        {
            *fOut << "Connect " << connect << " has no dependency but has large element\n";
            check = 0;
        }
    }
	TPZStack<int> deplist;
	BuildDependList(connect,deplist);
	int ndep = deplist.NElements();
	int id;
	for(id=0; id<ndep; id++) {
		TPZCompElSide smalll = FindElement(deplist[id]);
		TPZCompElSide large = smalll.LowerLevelElementList(1);
		if(!large.Exists()) {
			check = 0;
            (*fOut) << "VerifyConnect of connect " << connect << std::endl;
            (*fOut) << "deplist = " << deplist << " deplist[id] = " << deplist[id] << std::endl;
            (*fOut) << "Connect index of connect which is restrained deplist[id] = " << deplist[id] << std::endl;
            TPZConnect &c = fMesh->ConnectVec()[deplist[id]];
            c.Print(*fMesh,(*fOut));
            (*fOut) << "Element/side which contains deplist[id] side = " << smalll.Side() << "\n";
            smalll.Element()->Print(*fOut);
			(*fOut) << "VerifyConnect of " << connect << " inconsistent\n";
			continue;
		}
		TPZCompEl *cel = large.Element();
		int ncl = cel->NConnects();
		int icl;
		for(icl=0; icl<ncl; icl++) {
			if(cel->ConnectIndex(icl) == connect) break;
		}
		if(icl == ncl) {
			check = 0;
			(*fOut) << "VerifyConnect of " << connect << " inconsistent\n";
		}
	}
	return check;
}

void TPZCheckMesh::DependencyReport(int connect) {
	
	(*fOut) << "DependencyReport for connect " << connect << endl;
	TPZCompElSide large = FindElement(connect);
	DependencyReport(connect, large);
}

void TPZCheckMesh::DependencyReport(int connect, TPZCompElSide & large) {
	
	TPZStack<TPZCompElSide> smalll;
	large.HigherLevelElementList(smalll,1,1);
	int nsmall = smalll.NElements();
	int is;
	for(is=0; is<nsmall; is++) {
		int conind = smalll[is].Element()->ConnectIndex(smalll[is].Side());
		(*fOut) << "Connect " << conind << " may depend on " << connect << endl;
	}
	if(large.Reference().Dimension() == 1) {
		TPZStack<TPZCompElSide> equal;
		large.EqualLevelElementList(equal,1,0);
		equal.Push(large);
		TPZStack<TPZCompElSide> highdim;
		int neq = equal.NElements();
		int eq;
		for(eq=0; eq<neq; eq++) equal[eq].HigherDimensionElementList(highdim,1,0);
		large.HigherDimensionElementList(highdim,1,0);
		int nhigh = highdim.NElements();
		int ih;
		for(ih=0; ih<nhigh; ih++) {
			DependencyReport(connect,highdim[ih]);
		}
	}
}

int TPZCheckMesh::VerifyAllConnects() {
	int ncon = fMesh->NConnects();
	int i;
	int check = 1;
	for (i=0; i<ncon; i++){
		//(*fOut) << "Startint to check consistency for connect " << i << endl;
		if (!VerifyConnect(i) || !VerifyCompatibilityBetweenNShapesAndBlockSize(i)) {
			check = 0;
			(*fOut) << "Check failed for connect: " << i << endl;
		}
	}
	return check;
}

int TPZCheckMesh::CheckDimensions() {
	int check = 0;
	check = CheckElementShapeDimension();
	int check2 = CheckConstraintDimension();
	return (check || check2);
}

int TPZCheckMesh::CheckElementShapeDimension() {
	int check = 0;
	int nelem = fMesh->ElementVec().NElements();
	int el;
	for(el=0; el<nelem; el++) {
		TPZCompEl *cel = fMesh->ElementVec()[el];
		if(!cel) continue;
		TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (cel);
		if(!cint) continue;
		int nc = cint->NConnects();
		int ic;
		for(ic=0; ic<nc; ic++) {
			int cind = cint->ConnectIndex(ic);
            TPZConnect &c = cint->Connect(ic);
			int nshape = cint->NConnectShapeF(ic,c.Order());
            if (c.NShape() != nshape) {
                cout << "TPZCheckMesh::CheckElement.. el = " << el << " connect = "
                << cind << " nshape = " << nshape << " c.NShape " << c.NShape() << endl;
                check =1 ;
            }
			fNState = 1;                                                      /// WARNING
			int seqnum = fMesh->ConnectVec()[cind].SequenceNumber();
			int blsize = fMesh->Block().Size(seqnum);
			if(nshape*fNState != blsize) {
				cout << "TPZCheckMesh::CheckElement.. el = " << el << " connect = "
				<< cind << " nshape = " << nshape << " blsize " << blsize << endl;
				cint->Print(cout);
				check = 1;
			}
		}
	}
	return check;
}

int TPZCheckMesh::CheckConstraintDimension()
{
	int check = 0;
	int ncon = fMesh->ConnectVec().NElements();
    std::set<int64_t> badcon;
	for(int64_t ic=0; ic<ncon; ic++) {
		TPZConnect &con = fMesh->ConnectVec()[ic];
		int iseqnum = con.SequenceNumber();
        if (iseqnum < 0) {
            continue;
        }
		int iblsize = fMesh->Block().Size(iseqnum);
		if(con.HasDependency()) {
			TPZConnect::TPZDepend *dep = con.FirstDepend();
			while(dep) {
				int r = dep->fDepMatrix.Rows();
				int c = dep->fDepMatrix.Cols();
				int jcind = dep->fDepConnectIndex;
				TPZConnect &jcon = fMesh->ConnectVec()[jcind];
				int jseqnum = jcon.SequenceNumber();
				int jblsize = fMesh->Block().Size(jseqnum);
				if(r*fNState != iblsize || c*fNState != jblsize) {
                    badcon.insert(ic);
					cout << "TPZCheckMesh::CheckConstraintDi.. ic = " << ic
					<< " depends on " << jcind << " r " << r << " iblsize "
					<< iblsize  << " c " << c << " jblsize " << jblsize << endl;
					check = 1;
				}
				dep = dep->fNext;
			}
		}
	}
    TPZVec<int64_t> connecttoel(fMesh->NConnects(),-1);
    int64_t nelem = fMesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = fMesh->Element(el);
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            int64_t cindex = cel->ConnectIndex(ic);
            if (connecttoel[cindex] == -1) {
                connecttoel[cindex] = el;
            }
        }
    }
    std::set<int64_t> badelements;
    for (std::set<int64_t>::iterator it = badcon.begin(); it != badcon.end(); it++) {
        /// find the element that contains badcon
        if (connecttoel[*it] == -1) {
            std::cout << "Connect " << *it << " has dependency but does not belong to an element\n";
        }
        TPZCompEl *cel = fMesh->Element(connecttoel[*it]);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            std::cout << "I dont understand\n";
            continue;
        }
        int nc = intel->NConnects();
        for (int ic = 0; ic<nc; ic++) {
            if (intel->ConnectIndex(ic) == *it)
            {
                TPZGeoEl *ref = intel->Reference();
                TPZGeoElSide gelside(ref,ic);
                TPZCompElSide large = gelside.LowerLevelCompElementList2(1);
                if (!large) {
                    std::cout << "I dont understand\n";
                    continue;
                }
                TPZGeoElSide gellarge = large.Reference();
                TPZStack<TPZCompElSide> equal;
                gellarge.EqualLevelCompElementList(equal, 1, 0);
                equal.Push(large);
                std::cout << "Connect index " << *it << " belongs to element " << connecttoel[*it] << " and is restrained by element ";
                for(int i=0; i< equal.size(); i++) std::cout << equal[i].Element()->Index() << " ";
                std::cout << std::endl;
                badelements.insert(large.Element()->Index());
                
                int nclarge = large.Element()->NConnects();
                std::set<int64_t> largecon;
                for (int icl=0; icl< nclarge; icl++) {
                    largecon.insert(large.Element()->ConnectIndex(icl));
                }
                TPZConnect &c = fMesh->ConnectVec()[*it];
                TPZConnect::TPZDepend *dep = c.FirstDepend();
                while(dep)
                {
                    if (largecon.find(dep->fDepConnectIndex) == largecon.end()) {
                        std::cout << "The connect " << *it << " belongs to element " << connecttoel[*it] << " as connect " << ic << std::endl;
                        std::cout << "The element is restrained along this side to element " << large.Element()->Index() << std::endl;
                        std::cout << "The connect " << *it << " depends on connect " << dep->fDepConnectIndex << " The connect of the large element are ";
                        for (int ic=0; ic<nclarge; ic++) {
                            std::cout << large.Element()->ConnectIndex(ic) << " ";
                        }
                        std::cout << std::endl;
                    }
                    dep = dep->fNext;
                }
            }
        }
    }
    for (std::set<int64_t>::iterator it=badelements.begin(); it != badelements.end(); it++) {
        fMesh->Element(*it)->Print();
    }
	return check;
}

int TPZCheckMesh::CheckConnectOrderConsistency() {
	
	int nel = fMesh->ElementVec().NElements();
	int nstate = (fMesh->MaterialVec().begin()->second)->NStateVariables();
	int iel;
	for(iel = 0; iel<nel; iel++) {
		TPZCompEl *cel = fMesh->ElementVec()[iel];
		if(!cel) continue;
		TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
		if(!intel->CheckElementConsistency()) {
			intel->Print();
			return iel;
		}
		int nc = intel->NConnects();
		int ic;
		for(ic = 0; ic<nc; ic++) {
			TPZConnect &c = intel->Connect(ic);
			int nshape = intel->NConnectShapeF(ic,c.Order());
			if(c.CheckDependency(nshape, fMesh, nstate) == -1) {
				cout << "TPZCheckMesh inconsistent dependency" << endl;
				intel->Print();
				return iel;
			}
		}
	}
	return -1;
}

/**
 * @brief This method verifies if the sequence numbers of dependent connects and/or condensed connect are ordered at the back of the sequence
 */
int TPZCheckMesh::CheckConnectSeqNumberConsistency()
{
    bool wrong = false;
    int64_t numindepconnect = 0;
    int64_t nconnect = fMesh->NConnects();
    for (int64_t ic=0; ic<nconnect; ic++) {
        TPZConnect &c = fMesh->ConnectVec()[ic];
        if (!c.HasDependency() && c.NElConnected() && !c.IsCondensed()) {
            numindepconnect++;
        }
    }
    if (numindepconnect == 0) {
        return 0;
    }
    for (int64_t ic=0; ic<nconnect; ic++) {
        TPZConnect &c = fMesh->ConnectVec()[ic];
        if (c.HasDependency() || !c.NElConnected() || c.IsCondensed()) {
            int64_t seqnum = c.SequenceNumber();
            if (seqnum < numindepconnect && seqnum != -1) { // @omar need more information about this case seqnum != -1
                *fOut << "Connect index " << ic << " has inconsistent sequence number " << seqnum << " nindependent connects " << numindepconnect << std::endl;
                wrong = true;
            }
        }
    }
    if (wrong == false) {
        return 0;
    }
    else
    {
        return 1;
    }
}


