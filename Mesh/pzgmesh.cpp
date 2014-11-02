/**
 * @file
 * @brief Contains the implementation of the TPZGeoMesh methods.
 */

#include "pzgmesh.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzvec_extras.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzmatrix.h"

#include "TPZRefPattern.h"
#include "TPZRefPatternDataBase.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"

#ifdef BORLAND
#include <io.h>
#include <fcntl.h>
#endif

#include <sstream>
#include <string>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzgeomesh"));
#endif


using namespace std;

TPZGeoMesh::TPZGeoMesh() :  fName(), fElementVec(0), fNodeVec(0)
{
	fReference = 0;
	fNodeMaxId = -1;
	fElementMaxId = -1;
    fDim = -1;
}

TPZGeoMesh::TPZGeoMesh(const TPZGeoMesh &cp) : TPZSaveable(cp)
{
	this->operator =(cp);
}

TPZGeoMesh & TPZGeoMesh::operator= (const TPZGeoMesh &cp )
{
	this->CleanUp();
	
	this->fName = cp.fName;
	long i, n = cp.fNodeVec.NElements();
	this->fNodeVec.Resize( n );
	for(i = 0; i < n; i++)
	{
		this->fNodeVec[i] = cp.fNodeVec[i];
	}//for
	n = cp.fElementVec.NElements();
	this->fElementVec.Resize( n );
	for(i = 0; i < n; i++)
	{
		this->fElementVec[i] = cp.fElementVec[i]->Clone(*this);
	}
	
	this->fNodeMaxId = cp.fNodeMaxId;
	this->fElementMaxId = cp.fElementMaxId;
	this->fInterfaceMaterials = cp.fInterfaceMaterials;
	
	this->fReference = NULL;
	
    this->fDim = cp.fDim;
	return *this;
}

TPZGeoMesh::~TPZGeoMesh()
{
	CleanUp();
}

/**Delete element, nodes, Cosys, boundary elements and boundary nodes in list*/
void TPZGeoMesh::CleanUp()
{
	long i, nel = fElementVec.NElements();
	for(i=0; i<nel; i++)
	{
		TPZGeoEl *el = fElementVec[i];
		if(el)
		{
			delete el;
			fElementVec[i] = 0;
		}
	}
	fElementVec.Resize(0);
	fElementVec.CompactDataStructure(1);
	fNodeVec.Resize(0);
	fNodeVec.CompactDataStructure(1);
}

void TPZGeoMesh::SetName (const std::string &nm)
{
	fName = nm;
}

void TPZGeoMesh::Print (std::ostream & out)
{
	out << "\n\t\t GEOMETRIC TPZGeoMesh INFORMATIONS:\n\n";
	out << "TITLE-> " << fName << "\n\n";
	out << "number of nodes               = " << fNodeVec.NElements() << "\n";
	out << "number of elements            = " << fElementVec.NElements()-fElementVec.NFreeElements() << "\n";
	
	out << "\n\tGeometric Node Information:\n\n";
	long i;
	long nnodes = fNodeVec.NElements();
	for(i=0; i<nnodes; i++)
	{
		out << "Index: " << i << " ";
		fNodeVec[i].Print(out);
	}
	out << "\n\tGeometric Element Information:\n\n";
	long nelem = fElementVec.NElements();
	for(i=0; i<nelem; i++)
	{
		if(fElementVec[i]) fElementVec[i]->Print(out);
		out << "\n";
	}
	
	out << "\nInterface materials : \n\n";
	InterfaceMaterialsMap::iterator w, e = this->fInterfaceMaterials.end();
	
	const int n = this->fInterfaceMaterials.size();
	int l, r, m;
	out << "number = " << n << "\n";
	out << "left material / right material -> interface material\n";
	for(w = this->fInterfaceMaterials.begin(); w != e; w++)
	{
		l = w->first.first;
		r = w->first.second;
		m = w->second;
		out << l << " / " << r << " -> " << m << "\n";
	}
}
void TPZGeoMesh::PrintTopologicalInfo(std::ostream & out)
{
	out << "TITLE-> " << fName << "\n";
	out << "Number of nodes       = " << fNodeVec.NElements() << "\n";
	out << "Number of elements    = " << fElementVec.NElements()-fElementVec.NFreeElements() << "\n";
	
	out << "\n\tGeometric Node Information:\n";
	long i;
	long nnodes = fNodeVec.NElements();
	for(i=0; i<nnodes; i++)
	{
		fNodeVec[i].Print(out);
	}
	out << "\n\tGeometric Element Information:\n";
	long nelem = fElementVec.NElements();
	for(i=0; i<nelem; i++)
	{
		if(fElementVec[i]) fElementVec[i]->PrintTopologicalInfo(out);
		out << "\n";
	}
}

void TPZGeoMesh::GetNodePtr(TPZVec<long> &nos,TPZVec<TPZGeoNode *> &nodep)
{
	long i,nnodes=nos.NElements();
	for(i = 0; i < nnodes; i++) nodep[i] = &fNodeVec[nos[i]];
}

void  TPZGeoMesh::ResetReference()
{
	TPZGeoEl *elp;
	long i,nelements=fElementVec.NElements();
	for(i=0;i<nelements;i++)
	{
		elp = fElementVec[i];
		if(elp) elp->ResetReference();
	}
	fReference = 0;
}

void TPZGeoMesh::RestoreReference(TPZCompMesh *cmesh)
{
	ResetReference();
	fReference = cmesh;
	TPZGeoEl *gel;
	TPZCompEl *cel;
	long i,nelem = cmesh->ElementVec().NElements();
	for(i=0;i<nelem;i++)
	{
		cel = cmesh->ElementVec()[i];
		if(cel)
		{
			gel = cel->Reference();
			if(!gel)
			{
				PZError << "RestoreReference incomplete. Exist computational element with geometrical\n";
				PZError << "element not belongs to the current geometrical mesh.\n";
				return;
			}
			gel->SetReference(cel);
		}
	}
}

// GetBoundaryElements returns all elements beweeen NodFrom and NodTo counterclock wise
//		this method uses the connectivity of the elements
//		BuildConnectivity should be called to initialize the connectivity information
// 	this method will only work for grid with 2-D topology
//		the current version will only work for a grid with only one level
void TPZGeoMesh::GetBoundaryElements(long NodFrom, long NodTo,TPZStack<TPZGeoEl *> &ElementVec,TPZStack<int> &Sides)
{
	//Find a first element whose first node on the side is NodFrom
	//TPZGeoEl *def = 0;
	//TPZAVLMap<int,TPZGeoEl *> elmap(def);
	map<long,TPZGeoEl *> elmap;
	long i,nelements=NElements();
	for(i=0;i<nelements;i++)
	{
		TPZGeoEl *el = fElementVec[i];
		if(el) elmap[el->Id()]=fElementVec[i];
	}
	
	long currentnode = NodFrom;
	TPZGeoEl *candidate = 0;
	int candidateside = 0;
	while(currentnode != NodTo)
	{
		// put all elements connected to currentnode in elmap, eliminate the elements
		//	from elmap which do not contain the node
		BuildElementsAroundNode(currentnode,elmap);
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
		{
			std::stringstream sout;
			std::map<long, TPZGeoEl *>::iterator it;
			sout << "Elements around node " << currentnode << " : ";
			for(it=elmap.begin(); it!=elmap.end(); it++)
			{
				sout << it->second->Index() << "|";
				int in;
				for(in=0; in<it->second->NNodes(); in++)
				{
					sout << it->second->NodeIndex(in) << ":";
				}
				sout << " -- ";
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif
		
		// find, within elmap the element which has currentnode as its first boundary side
		// node
		FindElement(elmap, currentnode, candidate, candidateside);
		//	if the element found is already contained in the list, we have a circular list
		// if no element was found, the topology may not be two dimensional
		if(!candidate)
		{
			LOGPZ_WARN(logger,"GetBoundaryElements no adjacent element found");
			break;
		}
		long index = 0;
		long nelvec = ElementVec.NElements();
		while(index<nelvec && ElementVec[index] != candidate) index++;
		if(index <nelvec && Sides[index]==candidateside) break;
		ElementVec.Push(candidate);
		Sides.Push(candidateside);
		elmap.erase(elmap.begin(), elmap.end());//CleanUp();
		elmap[candidate->Id()] = candidate;
		// initialize the list in which to look for connected elements
		currentnode = candidate->SideNodeIndex(candidateside,1);
	}
}

// Find all elements in elmap or neighbour of elements in elmap which contain a node
//void TPZGeoMesh::BuildElementsAroundNode(int currentnode,TPZAVLMap<int,TPZGeoEl*> &elmap){
void TPZGeoMesh::BuildElementsAroundNode(long currentnode,map<long,TPZGeoEl*> &elmap)
{
	// first eliminate all elements which do not contain currentnode
	//TPZPix iel = elmap.First();
	map<long, TPZGeoEl *>::iterator ielprev,iel=elmap.begin();
	TPZGeoEl *el;
	long i;
	while(iel!=elmap.end())
	{
		el = iel->second;
		ielprev=iel;
		iel++;//elmap.Next(iel);
		int numnode = el->NNodes();
		for(i=0; i< numnode; i++)
		{
			if(el->NodeIndex(i) == currentnode) break;
		}
		if(i == numnode)
		{
			elmap.erase(ielprev);
		}
	}
	iel = elmap.begin();
	while(iel!=elmap.end())
	{
		el = iel->second;//elmap.Contents(iel);
		iel++;//elmap.Next(iel);
		int nside = el->NSides();
		for(int is=0; is<nside; is++)
		{
			TPZGeoElSide neigh = el->Neighbour(is);
			if(!neigh.Exists()) continue;
			int numnode = neigh.Element()->NNodes();
			for(i=0; i< numnode; i++)
			{
				if(neigh.Element()->NodeIndex(i) == currentnode)
				{
					if(elmap.find(neigh.Element()->Id())==elmap.end())
					{
						// this should be implemented as a stack, so that we dont have to
						// go through the list again each time
						elmap[neigh.Element()->Id()] = neigh.Element();
						iel = elmap.begin();
					}
					break;	// get out of the loop over the nodes
				}
			}
		}
	}
}

// find, within elmap the element which has currentnode as its first boundary side
//  	node
//void TPZGeoMesh::FindElement(TPZAVLMap<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside) {
void TPZGeoMesh::FindElement(std::map<long,TPZGeoEl *> &elmap,long currentnode,TPZGeoEl* &candidate,int &candidateside)
{
    candidate = 0;
    //TPZPix iel = elmap.First();
    map<long , TPZGeoEl *>::iterator ielprev, iel = elmap.begin();
    while(iel!=elmap.end()) {
		TPZGeoEl *el = iel->second;//elmap.Contents(iel);
		ielprev=iel;
		iel++;//elmap.Next(iel);
		int ns = el->NSides();
		int is = el->NCornerNodes();
		for(; is < ns; is++) {
			TPZGeoElSide thisside(el,is);
			if (thisside.Dimension() > 1) continue;
			TPZGeoElSide neigh = el->Neighbour(is);
			TPZGeoElSide father = el->Father2(is);
			
#ifdef DEBUG2
			std::stringstream sout;
			sout << __PRETTY_FUNCTION__ << " for elidx " << el->Index() << std::endl;
			sout << "\tthiside " << el->Index() << "/" << is << std::endl;
			if (neigh.Element()) sout << "\tneighsd " << neigh.Element()->Index() << "/" << neigh.Side() << std::endl;
			if (father.Element())sout << "\tfathers " << father.Element()->Index() << "/" << father.Side() << std::endl;
#endif
			
			if(neigh == thisside && !father.Exists() && el->SideNodeIndex(is,0) == currentnode) {
				candidate = el;
				candidateside = is;
				
#ifdef DEBUG2
				sout << "\t\t\tNew Candidate found el/side = " << el << "/" << is << std::endl;
				LOGPZ_DEBUG (logger,sout.str().c_str());
#endif
				
				return;
			}
			
#ifdef DEBUG2
			else
			{
				sout << "candidate doesn't match...";
				LOGPZ_DEBUG(logger,sout.str().c_str());
			}
#endif
		}
	}
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "No neighbour found for node " << currentnode << " elmap \n";
		for(iel=elmap.begin(); iel!=elmap.end(); iel++)
		{
			iel->second->Print(sout);
		}
		LOGPZ_WARN(logger,sout.str());
	}
#endif
	
}

TPZGeoNode *TPZGeoMesh::FindNode(TPZVec<REAL> &co)
{
    long i=0, in, nnodes = fNodeVec.NElements();
    while(i<nnodes && fNodeVec[i].Id() == -1) i++;
    if(i == nnodes) return 0;
    TPZGeoNode *gnkeep = &fNodeVec[i];
    REAL distkeep = 0.;
    for(in=0;in<3;in++)
		distkeep += (co[in]-(gnkeep->Coord(in)))*(co[in]-(gnkeep->Coord(in)));
    while(i< nnodes) {
		TPZGeoNode *gn = &fNodeVec[i];
		REAL dist = 0.;
		for(in=0;in<3;in++) dist += (co[in]-gn->Coord(in))*(co[in]-gn->Coord(in));
		if(dist < distkeep)
		{
			gnkeep = gn;
			distkeep = dist;
		}
		i++;
		while(i<nnodes && fNodeVec[i].Id() == -1) i++;
	}
	return gnkeep;
}

/** by Philippe 2013 */
/** @brief Returns the element that is close to the given point x */
TPZGeoEl * TPZGeoMesh::FindCloseElement(TPZVec<REAL> &x, long & InitialElIndex, int targetDim) const
{
    // this method will not work if targetDim == 0 because it navigates through elements of dimension targetDim
    TPZManVector<REAL,3> xcenter(3);
    TPZManVector<TPZManVector<REAL, 3>, 8> cornercenter;
    // what happens if gelnext == 0 or if InitialIndex is out of scope??
    if (InitialElIndex >= fElementVec.NElements()) {
        DebugStop();
    }
    TPZGeoEl *gelnext = ElementVec()[InitialElIndex];
    
    if (gelnext == 0) {
        DebugStop();
    }
    
    TPZGeoEl *gel = 0;
    REAL geldist;
    std::map<REAL,int> cornerdist;
    while (gel != gelnext) {
        gel = gelnext;
        gelnext = 0;
        cornerdist.clear();
        // compute the corner coordinates
        for (int ic=0; ic<gel->NCornerNodes(); ic++) {
            TPZManVector<REAL,3> xcorner(3);
            gel->NodePtr(ic)->GetCoordinates(xcorner);
            REAL dist = 0.;
            for (int i=0; i<3; i++) {
                dist += (x[i]-xcorner[i])*(x[i]-xcorner[i]);
            }
            dist = sqrt(dist);
            cornerdist[dist]=ic;
        }
        // compute the distance of the center of the element
        {
            geldist = 0.;
            TPZManVector<REAL,3> xcenter(3);
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            gelside.CenterX(xcenter);
            geldist = dist(x,xcenter);
        }
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "FindClosestElement Tried element index " << gel->Index() << std::endl;
            sout << "Distance from the center " << geldist << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // find the closest corner node
        REAL closestcorner = cornerdist.begin()->first;
        // if the center node is closer than the cornernode, return the element
        if (geldist < closestcorner || closestcorner < 1.e-15) {
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Distance from the closest corner " << closestcorner << "bailing out " << std::endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            InitialElIndex = gel->Index();
            return gel;
        }
        // look for all neighbours of the corner side
        TPZGeoElSide gelside(gel,cornerdist.begin()->second);
        TPZGeoElSide neighbour = gelside.Neighbour();
        std::map<REAL,int> distneigh;
        while (neighbour != gelside)
        {
            if (neighbour.Element()->Dimension() == targetDim)
            {
                TPZManVector<REAL,3> center(3);
                TPZGeoElSide centerneigh(neighbour.Element(),neighbour.Element()->NSides()-1);
                centerneigh.CenterX(center);
                REAL distcenter = dist(center,x);
                distneigh[distcenter] = neighbour.Element()->Index();
            }
            neighbour = neighbour.Neighbour();
        }
        // choose the element whose center is closest to the coordinate
        REAL gelnextdist = 0.;
        if (distneigh.size() == 0) {
            gelnext = gel;
            gelnextdist = geldist;
        }
        else {
            gelnext = ElementVec()[distneigh.begin()->second];
            gelnextdist = distneigh.begin()->first;
        }
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "Closest element index " << gelnext->Index() << std::endl;
            sout << "Distance from the center " << gelnextdist << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // return if its center distance is larger than the distance of the current element
        if (geldist <= gelnextdist) {
            InitialElIndex = gel->Index();
            return gel;
        }
    }
    InitialElIndex = gel->Index();
    return gel;
}

TPZGeoEl * TPZGeoMesh::FindElement(TPZVec<REAL> &x, TPZVec<REAL> & qsi, long & InitialElIndex, int targetDim)
{
    TPZGeoEl *res = FindApproxElement(x, qsi, InitialElIndex, targetDim);
    TPZManVector<REAL,3> xaprox(3);
    res->X(qsi, xaprox);
    REAL dist = 0.;
    for (int i=0; i<3; i++) {
        dist += (xaprox[i]-x[i])*(xaprox[i]-x[i]);
    }
    dist = sqrt(dist);
    REAL zero;
    ZeroTolerance(zero);
    if (dist > zero*100)
    {
        std::stringstream sout;
        sout << "Element not found for coordinate x " << x << std::endl;
        sout << "Coordinate found " << xaprox << std::endl;
        sout << "Distance error " << dist << std::endl;
        sout << "Closest element index " << res->Index() << " El param " << qsi << std::endl;
        LOGPZ_ERROR(logger, sout.str())
    }
    return res;
}



/** @brief find an element/parameter close to the point */
TPZGeoEl *TPZGeoMesh::FindApproxElement(TPZVec<REAL> &x, TPZVec<REAL> & qsi, long & InitialElIndex, int targetDim)
{
    FindCloseElement(x, InitialElIndex, targetDim);
    TPZGeoEl * gel = this->ElementVec()[InitialElIndex]->LowestFather();
 
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "x coordinate " << x << std::endl;
        sout << "element index " << gel->Index() << std::endl;
        sout << "coordinates of the corner nodes of the element\n";
        for (int i=0; i<gel->NNodes(); i++) {
            TPZGeoNode *ptr = gel->NodePtr(i);
            for (int ic=0; ic<3; ic++) {
                sout << ptr->Coord(ic) << " ";
            }
            sout << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    qsi.Resize(targetDim, 0.);
    qsi.Fill(0.);
    
    if(gel->Dimension() != targetDim)
    {
        // the targetDim should be respected because it is a parameter of FindCloseElement
        DebugStop();
    }
    if(gel->Dimension() != targetDim)
    {
        DebugStop();
    }
    REAL zero;
    ZeroTolerance(zero);
    
    std::set<TPZGeoEl *> tested;
    // this method will call ComputeXInverse if the element dimension != 3
    if(gel->ComputeXInverse(x, qsi,zero*100.) == true)
    {
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            LOGPZ_DEBUG(logger, "Going into the FindSubElement alternative")
        }
#endif
        gel = FindSubElement(gel, x, qsi, InitialElIndex);
        return gel;
    }
    tested.insert(gel);
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Looking for x = " << x << std::endl;
        sout << "Tried out gel index " << gel->Index() << std::endl;
        sout << "found qsi " << qsi << std::endl;
        TPZManVector<REAL> locx(3);
        gel->X(qsi, locx);
        sout << "found x co " << locx << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZManVector<REAL,3> projection(gel->Dimension());
    int side = -1;
    side = gel->ProjectInParametricDomain(qsi, projection);

    
    TPZGeoElSide mySide(gel,side);
    TPZManVector<REAL,3> xclose(3);
    mySide.X(projection, xclose);
    REAL mindist = 0.;
    for (int i=0; i<3; i++) {
        mindist += (xclose[i]-x[i])*(xclose[i]-x[i]);
    }
    TPZGeoElSide bestside(mySide);
    TPZManVector<REAL,3> bestproj(projection);
    
    
    TPZStack<TPZGeoElSide> allneigh;
    TPZGeoElSide neighSide(mySide.Neighbour());
    while (neighSide != mySide) {
        if (neighSide.Element()->Dimension() == targetDim && tested.find(neighSide.Element()) == tested.end()) {
            allneigh.Push(neighSide);
        }
        neighSide = neighSide.Neighbour();
    }
    
    while (allneigh.size())
    {
        TPZGeoElSide  gelside = allneigh.Pop();
        TPZGeoEl *locgel = gelside.Element();
        if (tested.find(locgel) != tested.end()) {
            continue;
        }
        qsi.Fill(0.);
        if(locgel->ComputeXInverse(x, qsi, zero*100.) == true)
        {
#ifdef LOG4CXX
            if (logger->isDebugEnabled())LOGPZ_DEBUG(logger, "FOUND ! Going into the FindSubElement alternative")
#endif
            gel = FindSubElement(locgel, x, qsi, InitialElIndex);
            return gel;
        }
        
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Looking for x = " << x << std::endl;
            sout << "Tried out gel index " << locgel->Index() << std::endl;
            sout << "found qsi " << qsi << std::endl;
            TPZManVector<REAL> locx(3);
            gel->X(qsi, locx);
            sout << "found x co " << locx << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

        
        tested.insert(locgel);
        int side = -1;
        side = locgel->ProjectInParametricDomain(qsi, projection);
        TPZManVector<REAL,3> xclose(3);
        TPZGeoElSide locgelside(locgel,side);
        locgelside.X(projection, xclose);
        REAL dist = 0.;
        for (int i=0; i<3; i++) {
            dist += (xclose[i]-x[i])*(xclose[i]-x[i]);
        }
        if (dist < mindist) {
            mindist = dist;
            bestside = locgelside;
            bestproj = projection;
        }
        mySide = TPZGeoElSide(locgel,side);
        TPZGeoElSide neighSide(mySide.Neighbour());
        while (neighSide != mySide) {
            if (neighSide.Element()->Dimension() == targetDim && tested.find(neighSide.Element()) == tested.end()) {
                allneigh.Push(neighSide);
            }
            neighSide = neighSide.Neighbour();
        }
    }
    TPZGeoEl *bestgel = bestside.Element();
    qsi = bestproj;
    gel = FindSubElement(bestgel, x, qsi, InitialElIndex);
    return gel;

    return bestgel;
}

TPZGeoEl * TPZGeoMesh::FindSubElement(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> & qsi, long & InitialElIndex)
{
    REAL Tol;
    ZeroTolerance(Tol);
    
#ifdef DEBUG
    TPZManVector<REAL,3> locqsi(qsi);
    if(gel->ComputeXInverse(x,locqsi,Tol*100.) == false)
    {
        //The given gel does NOT contains the given x!!!
        std::cout << "FindSubElement called for non conforming point\n";
    }
#endif
    TPZManVector<REAL,3> locx(x);
    
    if(gel->HasSubElement() == false)
    {
        InitialElIndex = gel->Index();
        return gel;
    }
    else
    {
        TPZStack<TPZGeoEl*> subElements;
        gel->GetAllSiblings(subElements);
        
        int nsons = subElements.NElements();
        TPZGeoEl * son = NULL;
        TPZVec< TPZVec<REAL> > qsiSonVec(nsons);
        for(int s = 0; s < nsons; s++)
        {
            son = subElements[s];
            qsiSonVec[s].Resize(son->Dimension(),0.);
            if(son->ComputeXInverse(locx, qsiSonVec[s],Tol*100.))
            {
                qsi = qsiSonVec[s];
                InitialElIndex = son->Index();
                return FindSubElement(son, locx, qsi, InitialElIndex);
            }
        }
        
        //no son found (accuracy problems)
        {
            std::map<REAL,int> dist;
            for(int s = 0; s < nsons; s++)
            {
                son = subElements[s];
                TPZVec<REAL> qsiSonProj(son->Dimension());
                son->ProjectInParametricDomain(qsiSonVec[s], qsiSonProj);
                REAL distToProj = 0.;
                for(int c = 0; c < son->Dimension(); c++)
                {
                    distToProj += (qsiSonProj[c] - qsiSonVec[s][c]) * (qsiSonProj[c] - qsiSonVec[s][c]);
                }
                dist[sqrt(distToProj)] = s;
            }
            int sonPosition = dist.begin()->second;//smaller distance to parametric domain
            
            qsi = qsiSonVec[sonPosition];
            son = subElements[sonPosition];
            InitialElIndex = son->Index();
            
            return FindSubElement(son, locx, qsi, InitialElIndex);
        }
    }
    
    return NULL;
}

void TPZGeoMesh::BuildConnectivity()
{
	
	TPZVec<int> SideNum(NNodes(),-1);
	TPZVec<TPZGeoEl *> NeighNode(NNodes(),0);
	long nelem = NElements();
	long iel = 0;
	for(iel=0; iel<nelem; iel++)
    {
		TPZGeoEl *gel = fElementVec[iel];
		if(!gel) continue;
		int ncor = gel->NCornerNodes();
		int in;
		for(in=0; in<ncor; in++) {
			long nod = gel->NodeIndex(in);
			if(SideNum[nod] == -1)
			{
				NeighNode[nod] = gel;
				SideNum[nod] = in;
				if(gel->SideIsUndefined(in)) gel->SetSideDefined(in);
			}
			else
			{
				TPZGeoElSide neigh(NeighNode[nod],SideNum[nod]);
				TPZGeoElSide gelside(gel,in);
				if(!neigh.NeighbourExists(gelside))
				{
					neigh.SetConnectivity(gelside);
				}
			}
		}
    }
	for(iel=0; iel<nelem; iel++)
	{
		TPZGeoEl *gel = fElementVec[iel];
		if(!gel) continue;
		int ncor = gel->NCornerNodes();
		int nsides = gel->NSides();
		int is;
		for(is=ncor; is<nsides; is++)
		{
			if( gel->SideIsUndefined(is))
			{
				gel->SetSideDefined(is);
				TPZGeoElSide gelside(gel,is);
				TPZStack<TPZGeoElSide> neighbours;
				gelside.ComputeNeighbours(neighbours);
				long nneigh = neighbours.NElements();
				long in;
				for(in=0; in<nneigh; in++) {
					if(neighbours[in].Side() == -1)
					{
						std::cout << "TPZGeoMesh::BuildConnectivity : Inconsistent mesh detected analysing element/side:" ;
						std::cout << gelside ;
						std::cout << std::endl;
						continue;
					}
					if(neighbours[in].Element()->SideIsUndefined(neighbours[in].Side()))
					{
						neighbours[in].Element()->SetSideDefined(neighbours[in].Side());
					}
					gelside.SetConnectivity(neighbours[in]);
				}
			}
		}
    }
	
	//Verify node coordinates for curved elements
#ifdef DEBUG
	const long nel = this->NElements();
	for(long el = 0; el < nel; el++)
	{
		TPZGeoEl * gel = this->ElementVec()[el];
		if(!gel) continue;
		gel->VerifyNodeCoordinates();
	}//for el
#endif
	
	//Build the data structure of blend elements
	long Qelem = this->NElements();
	for(long el = 0; el < Qelem; el++)
	{
		TPZGeoEl * gel = this->ElementVec()[el];
		if(!gel) continue;
		gel->BuildBlendConnectivity();
	}//for el
	
}

void TPZGeoMesh::BuildConnectivityOld() {
	
	TPZVec<int> SideNum(NNodes(),-1);
	TPZVec<TPZGeoEl *> NeighNode(NNodes(),0);
	long nelem = NElements();
	long iel = 0;
	while(iel<nelem && fElementVec[iel] == 0) iel++;
	
	long numsearch =1;
	// if there are no elements, do nothing
	while(iel < nelem) {
		TPZGeoEl *el = fElementVec[iel];
		int numsides = el->NSides();
		int side;
		for(side = 0;side<numsides;side++)
		{
			// check whether all entries in NeighNode are equal
			
			int equalnode = 1;
			int numsidenodes = el->NSideNodes(side);
			long sidenode = el->SideNodeIndex(side,0);
			TPZGeoEl *neigh = NeighNode[sidenode];
			int sidenumber = SideNum[sidenode];
			for(long sn = 0;sn < numsidenodes; sn++) 
			{
				sidenode = el->SideNodeIndex(side,sn);
				if (neigh != NeighNode[sidenode])
				{
					equalnode=0;
					break;
				}
			}
			
			if(equalnode && neigh == 0)
			{
				if(el->SideIsUndefined(side))
				{
					int elloaded = 0;
					for(int in=0; in<el->NNodes(); in++)
					{
						if(NeighNode[el->NodeIndex(in)] == el) elloaded = 1;
					}
					// this element is not loaded and its side is undefined
					
					// load the element side in the NeighNode vector
					for(int sn=0;!elloaded && sn < numsidenodes; sn++) \
					{
						sidenode = el->SideNodeIndex(side,sn);
						NeighNode[sidenode] = el;
						SideNum[sidenode] = side;
					}
					numsearch++;
				}
			} else if(equalnode && side == sidenumber && neigh == el)
			{
				// unload the element side
				for(int sn=0;sn < numsidenodes; sn++)
				{
					sidenode = el->SideNodeIndex(side,sn);
					NeighNode[sidenode] = 0;
					SideNum[sidenode] = -1;
				}
				// if no neighbouring element was detected during the loop
				//    define the element side as undefined
				TPZGeoElSide neighbour = el->Neighbour(side);
				if(!neighbour.Exists()) el->SetSideDefined(side);
				numsearch++;
			}
			else if(equalnode && neigh != el)
			{
				// we found a neigbour
				TPZManVector<long> SideNodes(numsidenodes);
				// detect which side of the neigbour is loaded witin NeighNode
				for(int sn=0;sn < numsidenodes; sn++)
				{
					sidenode = el->SideNodeIndex(side,sn);
					SideNodes[sn] = sidenode;
				}
				int neighside = neigh->WhichSide(SideNodes);
				TPZGeoElSide neighbour(neigh,neighside);
				// WhichSide will tell the side number which contains the vector
				//    of node ids SideNodes
				if(neighbour.Side() != -1 && !el->NeighbourExists(side,neighbour))
				{
					TPZGeoElSide(el,side).SetConnectivity(neighbour);
					numsearch++;
				}
			}
		} // loop over the sides
		iel++;
		while(iel<nelem && fElementVec[iel] == 0) iel++;
		if(iel==nelem && numsearch)
		{
			numsearch = 0;
			iel = 0;
			while(iel<nelem && fElementVec[iel] == 0) iel++;
		}
	}
}

//Cedric : 03/03/99
TPZGeoEl *TPZGeoMesh::FindElement(long elid)
{
	long nel = fElementVec.NElements();
	TPZGeoEl *gel = 0;
	for(long i=0;i<nel;i++)
	{
    	gel = fElementVec[i];
		if(gel && gel->Id()==elid) break;
	}
	return gel;
}

long TPZGeoMesh::ElementIndex(TPZGeoEl *gel)
{
	long i=0;
	long index = gel->Index();
	if (ElementVec()[index] == gel) return index;
	
	long numel = ElementVec().NElements();
	while ( i < numel )
	{
		if (ElementVec()[i] == gel) break;
		i++;
	}
	if(i<numel) return i;
	else return -1;
}

long TPZGeoMesh::NodeIndex(TPZGeoNode *nod)
{
	long i=0;
	long numel = NodeVec().NElements();
	while ( i < numel )
	{
		if (&NodeVec()[i] == nod) break;
		i++;
	}
	return i;
}

#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzrefpoint.h"
#include "pzshapepoint.h"

using namespace pzgeom;
using namespace pzrefine;
using namespace pzshape;

TPZGeoEl *TPZGeoMesh::CreateGeoElement(MElementType type,
                                       TPZVec<long>& nodeindexes,
                                       int matid,
                                       long& index,
                                       int reftype){
	if (reftype == 0)  
    {
        switch( type )
        {
            case 0://point
                return new TPZGeoElement< TPZGeoPoint, TPZRefPoint>(nodeindexes, matid, *this, index );
            case 1://line
                return new TPZGeoElement< TPZGeoLinear, TPZRefLinear>(nodeindexes, matid, *this, index );
            case 2://triangle
                return new TPZGeoElement< TPZGeoTriangle, TPZRefTriangle >(nodeindexes, matid, *this, index );
            case 3://quadrilatera
                return  new TPZGeoElement< TPZGeoQuad, TPZRefQuad >(nodeindexes, matid, *this, index );
            case 4://tetraedra
                return new TPZGeoElement< TPZGeoTetrahedra, TPZRefTetrahedra >(nodeindexes, matid, *this, index );
            case 5:
                return new TPZGeoElement< TPZGeoPyramid, TPZRefPyramid >(nodeindexes, matid, *this, index );
            case 6:
                return new TPZGeoElement< TPZGeoPrism, TPZRefPrism >(nodeindexes, matid, *this, index );
            case 7:
                return new TPZGeoElement< TPZGeoCube, TPZRefCube >(nodeindexes, matid, *this, index );
            default:
                PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
                << " type = " << type << std::endl;
                return NULL;
        }
    }
	else
	{
		//TPZAutoPointer<TPZRefPattern> ref = GetUniformPattern(type);
		switch( type )
		{
			case 0://point
			{
				TPZGeoElRefPattern<TPZGeoPoint> * gel = new TPZGeoElRefPattern<TPZGeoPoint> (nodeindexes, matid, *this, index);
				return gel;
			}
			case 1://line
			{
				TPZGeoElRefPattern < TPZGeoLinear > *gel = new TPZGeoElRefPattern < TPZGeoLinear >(nodeindexes, matid, *this, index);
				//gel->SetRefPattern (ref);
				return gel;
			}
			case 2://triangle
			{
				TPZGeoElRefPattern < TPZGeoTriangle > *gel = new TPZGeoElRefPattern < TPZGeoTriangle >(nodeindexes, matid, *this, index);
				//gel->SetRefPattern (ref);
				return gel;
			}
			case 3://quadrilatera
			{
				TPZGeoElRefPattern < TPZGeoQuad > * gel = new TPZGeoElRefPattern < TPZGeoQuad >(nodeindexes, matid, *this, index);
				//gel->SetRefPattern (ref);
				return gel;
			}
			case 4://tetraedra
			{
				TPZGeoElRefPattern < TPZGeoTetrahedra > *gel = new TPZGeoElRefPattern < TPZGeoTetrahedra >(nodeindexes, matid, *this, index);
				//gel->SetRefPattern (ref);
				return gel;
			}
			case 5://pyramid
			{
				TPZGeoElRefPattern < TPZGeoPyramid > *gel = new TPZGeoElRefPattern < TPZGeoPyramid >(nodeindexes, matid, *this, index);
				//gel->SetRefPattern (ref);
				return gel;
			}
			case 6://prism
			{
				TPZGeoElRefPattern < TPZGeoPrism > *gel = new TPZGeoElRefPattern < TPZGeoPrism >(nodeindexes, matid, *this, index);
				//gel->SetRefPattern (ref);
				return gel;
			}
			case 7://cube
			{
				TPZGeoElRefPattern < TPZGeoCube > *gel = new TPZGeoElRefPattern < TPZGeoCube >(nodeindexes, matid, *this, index);
				//gel->SetRefPattern (ref);
				return gel;
			}
			default:
			{
				PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
                << " type = " << type << std::endl;
				return NULL;
			}
		}
	}
	//return NULL;
}

TPZGeoEl *TPZGeoMesh::CreateGeoBlendElement(MElementType type, TPZVec<long>& nodeindexes, int matid, long& index)
{ 
	switch( type ){
		case 0://point
		{
			TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPoint> > (nodeindexes,matid,*this);
			return gel;
		}
		case 1://line
		{
			TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoLinear> > (nodeindexes,matid,*this);
			return gel;
		}
		case 2://triangle
		{
			TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (nodeindexes,matid,*this);
			return gel;
		}
		case 3://quadrilateral
		{
			TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (nodeindexes,matid,*this);
			return gel;
		}
		case 4://tetraedron
		{
			TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTetrahedra> > (nodeindexes,matid,*this);
			return gel;
		}
		case 5://pyramid
		{
			TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPyramid> > (nodeindexes,matid,*this);
			return gel;
		}
		case 6://prism
		{
			TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> > (nodeindexes,matid,*this);
			return gel;
		}
		case 7://cube
		{
			TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (nodeindexes,matid,*this);
			return gel;
		}
		default:
		{
			PZError << "TPZGeoMesh::CreateGeoElementRefPattern type element not exists:"
			<< " type = " << type << std::endl;
			return NULL;
		}
	}
}

int TPZGeoMesh::ClassId() const
{
	return TPZGEOMESHID;
}

void TPZGeoMesh::DeleteElement(TPZGeoEl *gel,long index)
{
	if(index < 0 || gel != fElementVec[index])
	{
		index = ElementIndex(gel);
		if(index < 0)
		{
			PZError << "TPZGeoMesh::DeleteElement index error\n";
			return;
		}
	}
	if (gel->HasSubElement())
	{
		for (int i=0;i<gel->NSubElements();i++)
		{
			TPZGeoEl* son = gel->SubElement(i);
			DeleteElement(son,son->Index());
		}
	}
	gel->RemoveConnectivities();
	if(gel) delete gel;
	fElementVec[index] = NULL;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZGeoMesh,TPZGEOMESHID>;
#endif

void TPZGeoMesh::Read(TPZStream &buf, void *context)
{
	try
	{
		TPZSaveable::Read(buf,context);
		int classid;
		buf.Read(&classid,1);
		
		if (classid != ClassId() )
		{
			std::cout << "ERROR RESTORING GEOMETRIC MESH!!\n";
		}
		
		buf.Read(&fName,1);
		ReadObjects(buf,fNodeVec,this);
		ReadObjectPointers(buf,fElementVec,this);
		buf.Read(&fNodeMaxId,1);
		buf.Read(&fElementMaxId,1);
		long ninterfacemaps;
		buf.Read(&ninterfacemaps,1);
		long c;
		for(c=0; c< ninterfacemaps; c++)
		{
			int vals[3];
			buf.Read(vals,3);
			fInterfaceMaterials[pair<int,int>(vals[0],vals[1])]=vals[2];
		}
        buf.Read(&fDim);
	}
	catch(const exception& e)
	{
		cout << "Exception catched! " << e.what() << std::endl;
		cout.flush();
		DebugStop();
	}
}

void TPZGeoMesh::Write(TPZStream &buf, int withclassid)
{
	try
	{
		TPZSaveable::Write(buf,withclassid);
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            LOGPZ_DEBUG(logger,__PRETTY_FUNCTION__);
        }
#endif
		int classid = ClassId();
		buf.Write(&classid,1);
		buf.Write(&fName,1);
		WriteObjects(buf,fNodeVec);
		WriteObjectPointers(buf,fElementVec);
		buf.Write(&fNodeMaxId,1);
		buf.Write(&fElementMaxId,1);
		long ninterfacemaps = fInterfaceMaterials.size();
		buf.Write(&ninterfacemaps,1);
		InterfaceMaterialsMap::iterator it = fInterfaceMaterials.begin();
		for(; it != fInterfaceMaterials.end(); it++)
		{
			int vals[3];
			vals[0] = (it->first).first;
			vals[1] = (it->first).second;
			vals[2] = it->second;
			buf.Write(vals,3);
		}
        buf.Write(&fDim);
		
	}
	catch(const exception& e)
	{
		cout << "Exception catched! " << e.what() << std::endl;
		cout.flush();
		DebugStop();
	}
}//method

int TPZGeoMesh::AddInterfaceMaterial(int leftmaterial, int rightmaterial, int interfacematerial)
{
	std::pair<int, int> leftright(leftmaterial, rightmaterial);
	InterfaceMaterialsMap::iterator w, e;
	e = fInterfaceMaterials.end();
	w = fInterfaceMaterials.find(leftright);
	if (w == e)
	{ 
		fInterfaceMaterials[leftright] = interfacematerial;
		return 1;
	}
	return 0;
}

int TPZGeoMesh::InterfaceMaterial(int leftmaterial, int rightmaterial)
{
	std::pair<int, int> leftright(leftmaterial, rightmaterial);
	InterfaceMaterialsMap::iterator w, e;
	e = fInterfaceMaterials.end();
	
	//trying to find an interface material associated to left and right materials
	w = fInterfaceMaterials.find(leftright);
	if (w != e) return w->second;
	
	//when left and right are equal and no exception case was inserted in interfacematerialmap return the same material of left and right
	if (leftmaterial == rightmaterial) return leftmaterial;
	
	//error message
	std::stringstream mess;
	mess << "\nTPZGeoMesh::InterfaceMaterial - Interface material not found left " << leftmaterial << " right " << rightmaterial ;
	PZError << mess.str()  << std::endl;
	return GMESHNOMATERIAL;
}

void TPZGeoMesh::ClearInterfaceMaterialsMap()
{
	InterfaceMaterialsMap::iterator b, e;
	b = fInterfaceMaterials.begin();
	e = fInterfaceMaterials.end();
	fInterfaceMaterials.erase(b, e);
}

void TPZGeoMesh::ResetConnectivities()
{
	TPZGeoElSide side;
	long iel;
	const long nelem = this->NElements();
	for(iel = 0; iel < nelem; iel++)
	{
		TPZGeoEl * gel = this->ElementVec()[iel];
		if (!gel) continue;
		const int nsides = gel->NSides();
		int is;
		for(is = 0; is < nsides; is++)
		{
			gel->SetNeighbour(is, side);
		}
	}
	this->SetName("reset");
}

/// Compute the area of the domain
REAL TPZGeoMesh::Area()
{
    TPZVec<int> NeedsComputing(NElements(),1);
    int meshdim = Dimension();
    long nel = NElements();
    REAL result = 0.;
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = ElementVec()[el];
        if (!gel || !NeedsComputing[el]) {
            continue;
        }
        if (gel->Dimension() != meshdim) {
            NeedsComputing[el] = 0;
            continue;
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            long neighindex = neighbour.Element()->Index();
            NeedsComputing[neighindex] = 0;
            neighbour = neighbour.Neighbour();
        }
        result += gel->SideArea(gel->NSides()-1);
    }
    return result;
}
