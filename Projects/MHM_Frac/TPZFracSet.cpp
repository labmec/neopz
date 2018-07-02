//
//  FracPre.cpp
//  PZ
//
//  Created by Philippe Devloo on 07/09/17.
//
//

#include "TPZFracSet.h"
#include "pzgeoel.h"

/// will project the endnodes on the MHM mesh if they fall within the tolerance
int64_t TPZFracSet::InsertNode(TPZGeoNode &gnode)
{
    uint64_t nodekey = GetLoc(gnode);
    if (fPointMap.find(nodekey) == fPointMap.end()) {
        int64_t index = fNodeVec.AllocateNewElement();
        gnode.SetCoord(ConvertNodekey(nodekey));
        fNodeVec[index] = gnode;
        fPointMap[nodekey] = index;
        return index;
    }
    else
    {
        return fPointMap[nodekey];
    }
}


static bool IsBetween(const double& x0, const double& x, const double& x1){
    return (x >= x0) && (x <= x1);
}

static bool FindIntersection(const double& x0, const double& y0,
                                  const double& x1, const double& y1,
                                  const double& a0, const double& b0,
                                  const double& a1, const double& b1,
                                  double& xy, double& ab) {
    // four endpoints are x0, y0 & x1,y1 & a0,b0 & a1,b1
    // returned values xy and ab are the fractional distance along xy and ab
    // and are only defined when the result is true
    
    bool partial = false;
    double denom = (b0 - b1) * (x0 - x1) - (y0 - y1) * (a0 - a1);
    if (denom == 0) {
        xy = -1;
        ab = -1;
    } else {
        xy = (a0 * (y1 - b1) + a1 * (b0 - y1) + x1 * (b1 - b0)) / denom;
        partial = IsBetween(-0.02, xy, 1.02);
        if (partial) {
            // no point calculating this unless xy is between 0 & 1
            ab = (y1 * (x0 - a1) + b1 * (x1 - x0) + y0 * (a1 - x1)) / denom;
        }
    }
    if ( partial && IsBetween(-0.02, ab, 1.02)) {
        ab = 1-ab;
        xy = 1-xy;
        return true;
    }  else return false;
}

/// split the fractures if they intersect
void TPZFracSet::ComputeFractureIntersections()
{
    int64_t ifr=0;
    while (ifr < fFractureVec.NElements()) {
        int64_t jfr = ifr+1;
        int64_t nel = fFractureVec.NElements();
        double x0,x1,y0,y1;
        int64_t inode0 = fFractureVec[ifr].fNodes[0];
        int64_t inode1 = fFractureVec[ifr].fNodes[1];
        x0 = fNodeVec[inode0].Coord(0);
        y0 = fNodeVec[inode0].Coord(1);
        x1 = fNodeVec[inode1].Coord(0);
        y1 = fNodeVec[inode1].Coord(1);
        while (jfr < nel) {
            double a0,a1,b0,b1,xy,ab;
            int64_t jnode0 = fFractureVec[jfr].fNodes[0];
            int64_t jnode1 = fFractureVec[jfr].fNodes[1];
            a0 = fNodeVec[jnode0].Coord(0);
            b0 = fNodeVec[jnode0].Coord(1);
            a1 = fNodeVec[jnode1].Coord(0);
            b1 = fNodeVec[jnode1].Coord(1);
            if(FindIntersection(x0, y0, x1, y1, a0, b0, a1, b1, xy, ab) && ((xy > 0.01 && xy < 0.99) || (ab > 0.01 && ab < 0.99)))
            {
                SplitFractures(ifr, xy, jfr, ab);
                // the fracture ifr has changed!
                inode0 = fFractureVec[ifr].fNodes[0];
                inode1 = fFractureVec[ifr].fNodes[1];
                x0 = fNodeVec[inode0].Coord(0);
                y0 = fNodeVec[inode0].Coord(1);
                x1 = fNodeVec[inode1].Coord(0);
                y1 = fNodeVec[inode1].Coord(1);
            }
            jfr++;
        }
        ifr++;
    }
}


/// split the fractures
void TPZFracSet::SplitFractures(int64_t ifrac, double xy, int64_t jfrac, double ab)
{
    TPZManVector<REAL,3> x(3),y(3),xyv(3);
    int64_t inode0 = fFractureVec[ifrac].fNodes[0];
    int64_t inode1 = fFractureVec[ifrac].fNodes[1];
    fNodeVec[inode0].GetCoordinates(x);
    fNodeVec[inode1].GetCoordinates(y);
    if (jfrac>=0)
    {
        TPZManVector<REAL,3> a(3),b(3),abv(3);
        int64_t jnode0 = fFractureVec[jfrac].fNodes[0];
        int64_t jnode1 = fFractureVec[jfrac].fNodes[1];
        fNodeVec[jnode0].GetCoordinates(a);
        fNodeVec[jnode1].GetCoordinates(b);
        REAL diff = 0;
        for (int i=0; i<3; i++) {
            xyv[i] = x[i]+xy*(y[i]-x[i]);
            abv[i] = a[i]+ab*(b[i]-a[i]);
            diff += (xyv[i]-abv[i])*(xyv[i]-abv[i]);
        }
        if(xy < 0.) xy = 0.;
        if(xy > 1.) xy = 1.;
        if(ab < 0.) ab = 0.;
        if(ab > 1.) ab = 1.;
        for (int i=0; i<3; i++) {
            xyv[i] = x[i]+xy*(y[i]-x[i]);
            abv[i] = a[i]+ab*(b[i]-a[i]);
        }
        if (jfrac != -1 && diff > 1.e-10) {
            DebugStop();
        }
    }
    else
    {
        for (int i=0; i<3; i++) {
            xyv[i] = x[i]+xy*(y[i]-x[i]);
        }
    }
    TPZGeoNode node;
    node.SetCoord(xyv);
    int64_t nodeindex = InsertNode(node);
    if(nodeindex != inode0 && nodeindex != inode1)
    {
        TPZFracture frac1(fFractureVec[ifrac].fOrigId, fFractureVec[ifrac].fMatId, fFractureVec[ifrac].fNodes[0], nodeindex);
        frac1.fPhysicalName = fFractureVec[ifrac].fPhysicalName;
        frac1.fFracPerm = fFractureVec[ifrac].fFracPerm;
        TPZFracture frac2(fFractureVec[ifrac].fOrigId, fFractureVec[ifrac].fMatId, nodeindex, fFractureVec[ifrac].fNodes[1]);
        frac2.fPhysicalName = fFractureVec[ifrac].fPhysicalName;
        frac2.fFracPerm = fFractureVec[ifrac].fFracPerm;
        fFractureVec[ifrac] = frac1;
        int64_t nfrac2 = fFractureVec.AllocateNewElement();
        fFractureVec[nfrac2] = frac2;
    }
    if (jfrac != -1 && nodeindex != fFractureVec[jfrac].fNodes[0] && nodeindex != fFractureVec[jfrac].fNodes[1])
    {
        TPZFracture frac1(fFractureVec[jfrac].fOrigId, fFractureVec[jfrac].fMatId, fFractureVec[jfrac].fNodes[0], nodeindex);
        TPZFracture frac2(fFractureVec[jfrac].fOrigId, fFractureVec[jfrac].fMatId, nodeindex, fFractureVec[jfrac].fNodes[1]);
        frac1.fPhysicalName = fFractureVec[jfrac].fPhysicalName;
        frac1.fFracPerm = fFractureVec[jfrac].fFracPerm;
        frac2.fPhysicalName = fFractureVec[jfrac].fPhysicalName;
        frac2.fFracPerm = fFractureVec[jfrac].fFracPerm;

        fFractureVec[jfrac] = frac1;
        int64_t nfrac2 = fFractureVec.AllocateNewElement();
        fFractureVec[nfrac2] = frac2;
        
    }
    
}

/// Split fractures by the MHM grid
void TPZFracSet::SplitFracturesByMHM()
{
    int64_t nfrac = fFractureVec.NElements();
    for (int64_t ifr=0; ifr<nfrac; ifr++) {
        std::pair<uint32_t,uint32_t> p1,p2;
        TPZManVector<REAL,3> co1(3),co2(3);
        int64_t no1 = fFractureVec[ifr].fNodes[0];
        int64_t no2 = fFractureVec[ifr].fNodes[1];
        fNodeVec[no1].GetCoordinates(co1);
        fNodeVec[no2].GetCoordinates(co2);
        p1 = NodeKey(fFractureVec[ifr].fNodes[0]);
        p2 = NodeKey(fFractureVec[ifr].fNodes[1]);
        int64_t domain1 = p1.first/fMHMSpacingInt[0];
        int64_t domain2 = p2.first/fMHMSpacingInt[0];
        if(p2.first > p1.first && p2.first % fMHMSpacingInt[0] == 0) domain2--;
        if (domain2<domain1) {
            DebugStop();
        }
        while (domain2 > domain1) {
            REAL frac = ((domain2)*fMHMSpacing[0]-co1[0])/(co2[0]-co1[0]);
            SplitFractures(ifr, frac, -1, -1.);
            fNodeVec[fFractureVec[ifr].fNodes[1]].GetCoordinates(co2);
            std::pair<uint32_t, uint32_t> pc = NodeKey(co2);
            if (pc.first%fMHMSpacingInt[0] != 0) {
                DebugStop();
            }
            domain2--;
        }
    }
    nfrac = fFractureVec.NElements();
    for (int64_t ifr=0; ifr<nfrac; ifr++) {
        std::pair<uint32_t,uint32_t> p1,p2;
        TPZManVector<REAL,3> co1(3),co2(3);
        fNodeVec[fFractureVec[ifr].fNodes[0]].GetCoordinates(co1);
        fNodeVec[fFractureVec[ifr].fNodes[1]].GetCoordinates(co2);
        p1 = NodeKey(fFractureVec[ifr].fNodes[0]);
        p2 = NodeKey(fFractureVec[ifr].fNodes[1]);
        int64_t domain1 = p1.second/fMHMSpacingInt[1];
        int64_t domain2 = p2.second/fMHMSpacingInt[1];
        if (domain1 < domain2)
        {
            if(p2.second % fMHMSpacingInt[1] == 0) domain2--;
            if (domain2<domain1) {
                DebugStop();
            }
            while (domain2 > domain1) {
                fNodeVec[fFractureVec[ifr].fNodes[1]].GetCoordinates(co2);
                REAL frac = ((domain2)*fMHMSpacing[1]-co1[1])/(co2[1]-co1[1]);
                SplitFractures(ifr, frac, -1, -1.);
                fNodeVec[fFractureVec[ifr].fNodes[1]].GetCoordinates(co2);
                std::pair<uint32_t, uint32_t> pc = NodeKey(co2);
                if (pc.second%fMHMSpacingInt[1] != 0) {
                    DebugStop();
                }
                domain2--;
            }
        
        }
        else if(domain1 > domain2)
        {
            if(p1.second % fMHMSpacingInt[1] == 0) domain1--;
            while (domain1 > domain2) {
                REAL frac = ((domain2+1)*fMHMSpacing[0]-co1[1])/(co2[1]-co1[1]);
                SplitFractures(ifr, frac, -1, -1.);
                fNodeVec[fFractureVec[ifr].fNodes[1]].GetCoordinates(co2);
                std::pair<uint32_t, uint32_t> pc = NodeKey(co2);
                if (pc.second%fMHMSpacingInt[1] != 0) {
                    DebugStop();
                }                
                domain2++;

            }
        }
    }
}

/// add the cornernodes of the MHM mesh
void TPZFracSet::AddMHMNodes()
{
    int64_t nfrac = fFractureVec.NElements();
    int numx = (fTopRight[0]-fLowLeft[0])/fMHMSpacing[0];
    int numy = (fTopRight[1]-fLowLeft[1])/fMHMSpacing[1];
    for (int i=0; i<=numx; i++) {
        for (int j=0; j<=numy; j++) {
            TPZManVector<REAL,3> co(3,0.);
            co[0] = fLowLeft[0]+i*fMHMSpacing[0];
            co[1] = fLowLeft[1]+j*fMHMSpacing[1];
            TPZGeoNode node;
            node.SetCoord(co);
            int64_t index = InsertNode(node);
            fNodeVec[index].GetCoordinates(co);
            std::pair<uint32_t, uint32_t> p1 = NodeKey(co);
            if (p1.first % fMHMSpacingInt[0] != 0 || p1.second%fMHMSpacingInt[1] != 0) {
                DebugStop();
            }
        }
    }
    // build the datastructure for horizontal and vertical lines
    fHorizontalLines.Resize(numy+1);
    fVerticalLines.Resize(numx+1);
    int64_t nnodes = fNodeVec.NElements();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> co(3);
        fNodeVec[in].GetCoordinates(co);
        uint64_t locprev = GetLoc(fNodeVec[in]);
        std::pair<uint32_t, uint32_t> p1 = NodeKey(co);
        if (p1.first%fMHMSpacingInt[0] == 0) {
            int line = p1.first/fMHMSpacingInt[0];
            REAL coprev = fNodeVec[in].Coord(0);
            REAL conew = line*fTol*fMHMSpacingInt[0];
            fNodeVec[in].SetCoord(0, conew);
            uint64_t locafter = GetLoc(fNodeVec[in]);
            if (locprev != locafter) {
                DebugStop();
            }
            fVerticalLines[line][co[1]] = in;
        }
        if (p1.second%fMHMSpacingInt[1] == 0) {
            int line = p1.second/fMHMSpacingInt[1];
            REAL coprev = fNodeVec[in].Coord(1);
            REAL conew = line*fTol*fMHMSpacingInt[1];
            fNodeVec[in].SetCoord(1, conew);
            uint64_t locafter = GetLoc(fNodeVec[in]);
            if (locprev != locafter) {
                DebugStop();
            }
            fHorizontalLines[line][co[0]] = in;
        }
    }
    
    int64_t nhor = fHorizontalLines.NElements();
    for (int64_t hor = 0; hor < nhor; hor++) {
        for (auto it = fHorizontalLines[hor].begin(); it != fHorizontalLines[hor].end(); it++) {
            auto it2 = it;
            it2++;
            if (it2 == fHorizontalLines[hor].end()) {
                continue;
            }
            TPZManVector<int64_t,2> nodes(2);
            nodes[0] = it->second;
            nodes[1] = it2->second;
            int64_t index;
            int matid = matid_MHM_line;
            if (hor == 0 || hor == nhor-1) {
                matid = matid_BC;
            }
            TPZFracture frac = TPZFracture(nfrac++, matid, nodes[0], nodes[1]);
            frac.fFracPerm = 0;
            frac.fPhysicalName="MHMLine";
            if (hor == 0 || hor == nhor-1) {
                frac.fPhysicalName = "BC";
            }
            index  = fFractureVec.AllocateNewElement();
            fFractureVec[index] = frac;
        }
    }
    int64_t nver = fVerticalLines.NElements();
    for (int64_t ver = 0; ver < nver; ver++) {
        for (auto it = fVerticalLines[ver].begin(); it != fVerticalLines[ver].end(); it++) {
            auto it2 = it;
            it2++;
            if (it2 == fVerticalLines[ver].end()) {
                continue;
            }
            TPZManVector<int64_t,2> nodes(2);
            nodes[0] = it->second;
            nodes[1] = it2->second;
            int matid = matid_MHM_line;
            if (ver==0 || ver==nver-1) {
                matid = matid_BC;
            }
            int64_t index;
            TPZFracture frac = TPZFracture(nfrac++, matid, nodes[0], nodes[1]);
            frac.fFracPerm = 0.;
            frac.fPhysicalName = "MHMLine";
            if (ver == 0 || ver == nver-1) {
                frac.fPhysicalName = "BC";
            }
            index  = fFractureVec.AllocateNewElement();
            fFractureVec[index] = frac;
        }
    }

}

/// return the index of the MHM domain of a fracture
int TPZFracSet::MHMDomain(TPZFracture &frac)
{
    TPZManVector<REAL,3> x1(3), x2(3), xmid(3);
    fNodeVec[frac.fNodes[0]].GetCoordinates(x1);
    fNodeVec[frac.fNodes[1]].GetCoordinates(x2);
    std::pair<uint32_t,uint32_t> key0 = NodeKey(frac.fNodes[0]);
    std::pair<uint32_t,uint32_t> key1 = NodeKey(frac.fNodes[1]);
    if(key0.first == key1.first && key0.first%fMHMSpacingInt[0] == 0)
    {
        return -1;
    }
    if(key0.second == key1.second && key0.second%fMHMSpacingInt[1] == 0)
    {
        return -1;
    }
    for (int i=0; i<3; i++) {
        xmid[i] = (x1[i]+x2[i])*0.5-fLowLeft[i];
    }
    int numfacex = (fTopRight[0]-fLowLeft[0])/fMHMSpacing[0];

    int numx = (xmid[0])/fMHMSpacing[0];
    int numy = (xmid[1])/fMHMSpacing[1];
    return numy*numfacex+numx;
}

/// CleanUp FractureNetwork
void TPZFracSet::CleanupFractureNetwork()
{
    ComputeFractureOverlap();
    ToGeoMesh();
    MergeParallelLines();
    DeleteVeryShortFractures(6);
    FromGeoMesh();
}

/// Transfer to the geometric mesh
void TPZFracSet::ToGeoMesh()
{
    fgmesh.CleanUp();
    fgmesh.NodeVec() = fNodeVec;
    int64_t nfrac = fFractureVec.NElements();
    for (int64_t ifr = 0; ifr < nfrac; ifr++) {
        TPZManVector<int64_t,2> nodes(2);
        nodes = fFractureVec[ifr].fNodes;
        int64_t index;
        fgmesh.CreateGeoElement(EOned, nodes, fFractureVec[ifr].fMatId, index);
    }
    fgmesh.BuildConnectivity();
    
    int64_t nel = fgmesh.NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fgmesh.Element(el);
        if(!gel) continue;
        if (gel->Neighbour(2).Element() != gel) {
            std::cout << "gel index " << el << " is overlapping with " << gel->Neighbour(2).Element()->Index() << std::endl;
            TPZGeoEl *neigh = gel->Neighbour(2).Element();
            int matid = min(gel->MaterialId(),neigh->MaterialId());
            if(gel->MaterialId() == neigh->MaterialId())
            {
                matid = gel->MaterialId();
            }
            else if(gel->MaterialId() == matid_BC || neigh->MaterialId() == matid_BC)
            {
                matid = matid_BC;
            }
            else if((gel->MaterialId() == matid_internal_frac && neigh->MaterialId() == matid_MHM_line) ||
                (gel->MaterialId() == matid_MHM_line && neigh->MaterialId() == matid_internal_frac))
            {
                matid = matid_MHM_frac;
            }
            else
            {
                DebugStop();
            }
            gel->SetMaterialId(matid);
            neigh->SetMaterialId(matid);
            // remove the element with the lowest index
            if (gel->Index() > neigh->Index())
            {
                // delete neigh
                neigh->RemoveConnectivities();
                int64_t neighindex = neigh->Index();
                delete neigh;
                fgmesh.ElementVec()[neighindex] = 0;
                std::string matname = fFractureVec[neighindex].fPhysicalName;
                fFractureVec[gel->Index()].fPhysicalName = matname + "_MHM";
            }
            else
            {
                gel->RemoveConnectivities();
                std::string matname = fFractureVec[gel->Index()].fPhysicalName;
                fFractureVec[neigh->Index()].fPhysicalName = matname + "_MHM";
                delete gel;
                fgmesh.ElementVec()[el] = 0;
            }
        }
        
    }
}

/// Transfer from the geometric mesh
void TPZFracSet::FromGeoMesh()
{
    if (fFractureVec.NElements() != fgmesh.ElementVec().NElements()) {
        DebugStop();
    }
    int64_t nel = 0;
    for (int64_t el=0; el<fgmesh.NElements(); el++) {
        if (fgmesh.Element(el)) {
            nel++;
        }
    }
    int64_t count = 0;
    for (int64_t el=0; el<fgmesh.NElements(); el++) {
        TPZGeoEl *gel = fgmesh.Element(el);
        if (gel && el != count) {
            fFractureVec[count] = fFractureVec[el];
        }
        if(gel) count++;
    }
    fFractureVec.Resize(nel);
    
}

static void Direction(TPZGeoEl *gel, TPZVec<REAL> &dir)
{
    REAL dirnorm = 0;
    for (int i=0; i<3; i++) {
        dir[i] = gel->NodePtr(1)->Coord(i) - gel->NodePtr(0)->Coord(i);
        dirnorm += dir[i]*dir[i];
    }
    dirnorm = sqrt(dirnorm);
    for (int i=0; i<3; i++) {
        dir[i] /= dirnorm;
    }
}

static REAL Length(TPZGeoEl *gel)
{
    TPZManVector<REAL,3> dir(3);
    REAL dirnorm = 0;
    for (int i=0; i<3; i++) {
        dir[i] = gel->NodePtr(1)->Coord(i) - gel->NodePtr(0)->Coord(i);
        dirnorm += dir[i]*dir[i];
    }
    dirnorm = sqrt(dirnorm);
    return dirnorm;
}

/// Merge lines which are parallel
void TPZFracSet::MergeParallelLines()
{
    int64_t nel = fgmesh.NElements();
    REAL maxcos = 0.;
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fgmesh.Element(el);
        if(!gel) continue;
        TPZManVector<REAL,3> dir1(3);
        Direction(gel, dir1);
        int nnodes = gel->NCornerNodes();
        for(int is = 0; is<nnodes; is++)
        {
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                TPZManVector<REAL,3> dir2(3);
                Direction(neighbour.Element(), dir2);
                if (neighbour.Side() != is) {
                    for (int i=0; i<3; i++) {
                        dir2[i] *= -1.;
                    }
                }
                REAL cosangle = 0.;
                for (int i=0; i<3; i++) {
                    cosangle += dir1[i]*dir2[i];
                }
                if (cosangle> maxcos) {
                    maxcos = cosangle;
                }
                if (cosangle > 0.99) {
                    std::cout << "Fractures " << gel->Index() << " and " << neighbour.Element()->Index() << " are parallel " << cosangle << "\n";
                    std::cout << "Index " << gel->NodeIndex(0) << " ";
                    gel->Node(0).Print();
                    std::cout << "Index " << gel->NodeIndex(1) << " ";
                    gel->Node(1).Print();
                    std::cout << "Index " << neighbour.Element()->NodeIndex(0) << " ";
                    neighbour.Element()->Node(0).Print();
                    std::cout << "Index " << neighbour.Element()->NodeIndex(1) << " ";
                    neighbour.Element()->Node(1).Print();
                    REAL l1 = Length(gel);
                    REAL l2 = Length(neighbour.Element());
                    if (l1 < l2) {
                        gel->RemoveConnectivities();
                        delete gel;
                        fgmesh.ElementVec()[el] = 0;
                        gel = 0;
                        break;
                    }
                    else
                    {
                        neighbour.Element()->RemoveConnectivities();
                        int64_t neighindex = neighbour.Element()->Index();
                        delete neighbour.Element();
                        fgmesh.ElementVec()[neighindex] = 0;
                        neighbour = gelside;
                    }
                }
                neighbour = neighbour.Neighbour();
            }
            if(!gel) break;
        }
    }
    std::cout << "max cosine angle " << maxcos << std::endl;
}

/// Delete very short fractures
void TPZFracSet::DeleteVeryShortFractures(REAL length)
{
    int64_t nel = fgmesh.NElements();
    double lmin = 2000.;
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fgmesh.Element(el);
        if (!gel) {
            continue;
        }
        int64_t inode0 = gel->NodeIndex(0);
        int64_t inode1 = gel->NodeIndex(1);
        if (gel->Neighbour(0).Element() == gel || gel->Neighbour(1).Element() == gel)
        {
            REAL l = Length(gel);
            if (l<lmin) {
                lmin = l;
            }
            if(l < length)
            {
                std::cout << "Deleting Fracture " << gel->Index() << " length " << l << std::endl;
                std::cout << "Index " << gel->NodeIndex(0) << " ";
                gel->Node(0).Print();
                std::cout << "Index " << gel->NodeIndex(1) << " ";
                gel->Node(1).Print();

                gel->RemoveConnectivities();
                delete gel;
                fgmesh.ElementVec()[el] = 0;
                gel = 0;
                
            }
        }
    }
    std::cout << "shortest fracture length " << lmin << std::endl;
}

static bool ParallelOverlap(const double& x0, const double& y0,
                             const double& x1, const double& y1,
                             const double& a0, const double& b0,
                             const double& a1, const double& b1)
{
 
    // four endpoints are x0, y0 & x1,y1 & a0,b0 & a1,b1
    // returned values xy and ab are the fractional distance along xy and ab
    // and are only defined when the result is true
    
    double norm1 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
    double norm2 = sqrt((b1-b0)*(b1-b0)+(a1-a0)*(a1-a0));
    double denom = (b0 - b1) * (x0 - x1) - (y0 - y1) * (a0 - a1)/(norm1*norm2);
    int count=0;
    if (fabs(denom) < 1.e-3)
    {
        double alphax0,alphax1,alphaa0,alphaa1;
        double distx0, distx1, dista0, dista1;
        alphax0 = ((x0-a0)*(a1-a0)+(y0-b0)*(b1-b0))/(norm2*norm2);
        {
            double t = y0*(a0-a1)-x0*(b0-b1)+a1*b0-a0*b1;
            distx0 = t*t/(norm2*norm2);
        }
        if (alphax0 > 0 && alphax0 < 1. && distx0 < 1) {
            count++;
        }
        alphax1 = ((x1-a0)*(a1-a0)+(y1-b0)*(b1-b0))/(norm2*norm2);
        {
            double t = y1*(a0-a1)-x1*(b0-b1)+a1*b0-a0*b1;
            distx1 = t*t/(norm2*norm2);
        }
        if (alphax1 > 0. && alphax1 < 1. && distx1 < 1) {
            count++;
        }
        alphaa0 = ((a0-x0)*(x1-x0)+(b0-y0)*(y1-y0))/(norm1*norm1);
        {
            double t = b0*(x0-x1)-a0*(y0-y1)+x1*y0-x0*y1;
            dista0 = t*t/(norm1*norm1);
        }
        if (alphaa0 > 0. && alphaa0 < 1. && dista0 <1) {
            count++;
        }
        alphaa1 = ((a1-x0)*(x1-x0)+(b1-y0)*(y1-y0))/(norm1*norm1);
        {
            double t = b1*(x0-x1)-a1*(y0-y1)+x1*y0-x0*y1;
            dista1 = t*t/(norm1*norm1);
        }
        if (alphaa1 > 0 && alphaa1 < 1. && dista1 < 1) {
            count++;
        }
        if (count == 1 || count == 3 || count == 4) {
            DebugStop();
        }
        if (count == 2) {
            return true;
        }
        return false;
    }
    else
    {
        return false;
    }
}

static void NodeSequence(TPZFracSet &fracset, TPZVec<int64_t> &nodes, TPZVec<int64_t> &orderedNodes)
{
    TPZManVector<REAL, 3> a(3,0.);
    TPZManVector<TPZManVector<REAL,3>, 4> X(4,a);
    for (int i=0; i<4; i++) {
        for (int j=0; j<3; j++) {
            X[i][j] = fracset.fNodeVec[nodes[i]].Coord(j);
        }
    }
    std::map<REAL,int64_t> indices;
    indices[0.] = nodes[0];
    indices[1.] = nodes[1];
    double denom = sqrt((X[1][0]-X[0][0])*(X[1][0]-X[0][0])+(X[1][1]-X[0][1])*(X[1][1]-X[0][1]));
    double x0=X[0][0],x1=X[1][0],y0=X[0][1],y1=X[1][1],a0=X[2][0],a1=X[3][0],b0=X[2][1],b1=X[3][1];
    double alphaa0 = ((a0-x0)*(x1-x0)+(b0-y0)*(y1-y0))/(denom*denom);
    double alphaa1 = ((a1-x0)*(x1-x0)+(b1-y0)*(y1-y0))/(denom*denom);
    indices[alphaa0] = nodes[2];
    indices[alphaa1] = nodes[3];
    int count = 0;
    for (auto it = indices.begin(); it != indices.end(); it++) {
        orderedNodes[count++] = it->second;
    }
    if (count != 4) {
        DebugStop();
    }

}
/// split the fractures if they intersect
void TPZFracSet::ComputeFractureOverlap()
{
    int64_t ifr=0;
    while (ifr < fFractureVec.NElements()) {
        int64_t jfr = ifr+1;
        int64_t nel = fFractureVec.NElements();
        double x0,x1,y0,y1;
        int64_t inode0 = fFractureVec[ifr].fNodes[0];
        int64_t inode1 = fFractureVec[ifr].fNodes[1];
        x0 = fNodeVec[inode0].Coord(0);
        y0 = fNodeVec[inode0].Coord(1);
        x1 = fNodeVec[inode1].Coord(0);
        y1 = fNodeVec[inode1].Coord(1);
        while (jfr < nel) {
            double a0,a1,b0,b1;
            int64_t jnode0 = fFractureVec[jfr].fNodes[0];
            int64_t jnode1 = fFractureVec[jfr].fNodes[1];
            std::set<int64_t> nodes;
            nodes.insert(inode0);
            nodes.insert(inode1);
            nodes.insert(jnode1);
            nodes.insert(jnode0);
            a0 = fNodeVec[jnode0].Coord(0);
            b0 = fNodeVec[jnode0].Coord(1);
            a1 = fNodeVec[jnode1].Coord(0);
            b1 = fNodeVec[jnode1].Coord(1);
            if(nodes.size() == 4 && ParallelOverlap(x0, y0, x1, y1, a0, b0, a1, b1))
            {
                TPZManVector<int64_t,4> nodes(4), nodeseq(4);
                nodes[0] = inode0; nodes[1] = inode1; nodes[2] = jnode0; nodes[3] = jnode1;
                NodeSequence(*this, nodes, nodeseq);
                fFractureVec[ifr].fNodes[0] = nodeseq[0];
                fFractureVec[ifr].fNodes[1] = nodeseq[1];
                fFractureVec[jfr].fNodes[0] = nodeseq[1];
                fFractureVec[jfr].fNodes[1] = nodeseq[2];
                int64_t index = fFractureVec.AllocateNewElement();
                TPZFracture frac(fFractureVec[ifr]);
                frac.fNodes[0] = nodeseq[2];
                frac.fNodes[1] = nodeseq[3];
                fFractureVec[index] = frac;
            }
            jfr++;
        }
        ifr++;
    }
    
}

/// Compute the mesh size at the nodes
void TPZFracSet::ComputeMeshSizeAtNodes()
{
    REAL minsize = this->fMinElementSize;
    REAL normalsize = this->fElementSize;
    ToGeoMesh();
    fMeshSizeAtNodes.Resize(fNodeVec.NElements(), normalsize);
    fMeshSizeAtNodes.Fill(normalsize);
    int64_t nel = fgmesh.NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fgmesh.Element(el);
        if (!gel) {
            continue;
        }
        REAL length = Length(gel);
        if (length < normalsize) {
            REAL meshsize = length < minsize ? minsize : length;
            int64_t inode0 = gel->NodeIndex(0);
            int64_t inode1 = gel->NodeIndex(1);
            fMeshSizeAtNodes[inode0] = min(fMeshSizeAtNodes[inode0],meshsize);
            fMeshSizeAtNodes[inode1] = min(fMeshSizeAtNodes[inode1],meshsize);
        }
    }
    int64_t nnode = fgmesh.NodeVec().NElements();
    for (int64_t node = 0; node < nnode; node++) {
        TPZManVector<REAL,3> co(3);
        fgmesh.NodeVec()[node].GetCoordinates(co);
        for (int ico =0; ico < 2; ico++)
        {
            int mhm = co[ico]/fMHMSpacing[ico];
            REAL mhm_min = mhm*fMHMSpacing[ico];
            REAL mhm_max = (mhm+1)*fMHMSpacing[ico];
            REAL dist = min(co[ico]-mhm_min,mhm_max-co[ico]);
            if (dist > 1.e-6 && dist < normalsize) {
                dist = max(dist,minsize);
                fMeshSizeAtNodes[node] = min(fMeshSizeAtNodes[node],dist);
            }
        }
    }
    int64_t nhor = fHorizontalLines.NElements();
    for (int64_t hr=0; hr<nhor; hr++) {
        for (auto it = fHorizontalLines[hr].begin(); it != fHorizontalLines[hr].end(); it++) {
            auto it2 = it;
            it2++;
            if (it2 == fHorizontalLines[hr].end()) {
                continue;
            }
            int64_t node0 = it->second;
            int64_t node1 = it2->second;
            REAL length = it2->first - it->first;
            if (length > normalsize) {
                continue;
            }
            if (length < minsize) {
                length = minsize;
            }
            fMeshSizeAtNodes[node0] = min(fMeshSizeAtNodes[node0],length);
            fMeshSizeAtNodes[node1] = min(fMeshSizeAtNodes[node1],length);
        }
    }
    int64_t nver = fVerticalLines.NElements();
    for (int64_t vr=0; vr<nver; vr++) {
        for (auto it = fVerticalLines[vr].begin(); it != fVerticalLines[vr].end(); it++) {
            auto it2 = it;
            it2++;
            if (it2 == fVerticalLines[vr].end()) {
                continue;
            }
            int64_t node0 = it->second;
            int64_t node1 = it2->second;
            REAL length = it2->first - it->first;
            if (length > normalsize) {
                continue;
            }
            if (length < minsize) {
                length = minsize;
            }
            fMeshSizeAtNodes[node0] = min(fMeshSizeAtNodes[node0],length);
            fMeshSizeAtNodes[node1] = min(fMeshSizeAtNodes[node1],length);
        }
    }
    FromGeoMesh();

}


/// verify the data structure consistency
void TPZFracSet::CheckPointMapConsistency()
{
    int64_t nnodes = fNodeVec.NElements();
    for (int64_t in=0; in<nnodes; in++) {
        uint64_t locnode  = GetLoc(fNodeVec[in]);
        if (fPointMap.find(locnode) == fPointMap.end()) {
            DebugStop();
        }
        int64_t nodeindex = fPointMap[locnode];
        if (nodeindex != in) {
            DebugStop();
        }
    }
}

