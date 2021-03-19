/**
 * @file
 * @brief Contains the implementation of the TPZReadMeshHR methods. 
 */

#include "pzhyperplane.h"

#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzgeoelrefless.h"
#include "pzlog.h"
#include "pzfmatrix.h"
#include "pzvec_extras.h"
#include "pznumeric.h"

#include <iostream>
#include <sstream>

#ifdef PZ_LOG
static TPZLogger logger("pz.pre.pzhyperplane");
#endif

// verify if the element is not warped
static void CheckElement(TPZGeoEl *gel);

/// Compute the intersection of the geometric mesh with the plane (only leaf elements)
void TPZHyperPlaneIntersect::Intersect(TPZGeoMesh &inputmesh, const TPZHyperPlane &plane, TPZGeoMesh &targetmesh)
{
    // data structure for keeping the correspondence between the edges and the intersecting node
    typedef std::map<std::pair<int64_t, int64_t>, int64_t> midsidetype;
    midsidetype midsidenodes;

    int64_t nelem = inputmesh.NElements();
    for (int64_t iel = 0; iel<nelem ; iel++) 
    {
        TPZGeoEl *gel = inputmesh.ElementVec()[iel];
        // we only consider leaf elements
        if (gel->HasSubElement()) {
            continue;
        }
        // loop over the sides
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            // we only consider ribs
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZManVector<int64_t,2> sidenodeindexes(2);
            sidenodeindexes[0] = gel->SideNodeIndex(is, 0);
            sidenodeindexes[1] = gel->SideNodeIndex(is, 1);
            // order the side
            if (sidenodeindexes[0] > sidenodeindexes[1]) {
                int64_t temp = sidenodeindexes[0];
                sidenodeindexes[0] = sidenodeindexes[1];
                sidenodeindexes[1] = temp;
            }
            std::pair<int64_t, int64_t> nodepair(sidenodeindexes[0],sidenodeindexes[1]);
            // verify if this pair was already identified
            if (midsidenodes.find(nodepair) != midsidenodes.end()) {
                continue;
            }
            TPZManVector<REAL,3> x0(3),x1(3);
            inputmesh.NodeVec()[sidenodeindexes[0]].GetCoordinates(x0);
            inputmesh.NodeVec()[sidenodeindexes[1]].GetCoordinates(x1);
            bool left0 = plane.IsLeft(x0);
            bool left1 = plane.IsLeft(x1);
            if (left0 != left1) {
                // we found a new edge
                REAL param = EdgeIntersect(gelside, plane);
                TPZManVector<REAL,1> par(1,param);
                TPZManVector<REAL, 3> xnew(3);
                gelside.X(par, xnew);
                int64_t newnode = targetmesh.NodeVec().AllocateNewElement();
                targetmesh.NodeVec()[newnode].Initialize(xnew, targetmesh);
                midsidenodes[nodepair] = newnode;
            }
        }
    }
    // at this point we have create a "cloud" of point intersecting the plane
    
    // now we will create the elements
    for (int64_t iel = 0; iel<nelem ; iel++) 
    {
        TPZGeoEl *gel = inputmesh.ElementVec()[iel];
        // we only consider leaf elements
        if (gel->HasSubElement()) {
            continue;
        }
        if (gel->Dimension() != 3) {
            continue;
        }
        // loop over the sides
        int nsides = gel->NSides();
        TPZStack<std::pair<int64_t,int64_t> > sidenodepair;
        for (int is=0; is<nsides; is++) {
            // we only consider ribs
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZManVector<int64_t,2> sidenodeindexes(2);
            sidenodeindexes[0] = gel->SideNodeIndex(is, 0);
            sidenodeindexes[1] = gel->SideNodeIndex(is, 1);
            // order the side
            if (sidenodeindexes[0] > sidenodeindexes[1]) {
                int64_t temp = sidenodeindexes[0];
                sidenodeindexes[0] = sidenodeindexes[1];
                sidenodeindexes[1] = temp;
            }
            std::pair<int64_t, int64_t> nodepair(sidenodeindexes[0],sidenodeindexes[1]);
            // verify if this pair was already identified
            if (midsidenodes.find(nodepair) == midsidenodes.end()) {
                continue;
            }
            sidenodepair.Push(std::pair<int64_t,int64_t>(is,midsidenodes[nodepair]));
        }
        // if the element has no intersecting sides, there is nothing to do
        if (sidenodepair.size() == 0) {
            continue;
        }
        if (sidenodepair.size() == 3) {
            TPZManVector<int64_t,3> nodeindices(3);
            for (int i=0; i<3; i++) {
                nodeindices[i] = sidenodepair[i].second;
            }
            // we should create an element with 3 nodes
            new TPZGeoElRefLess<pzgeom::TPZGeoTriangle>(nodeindices,1,targetmesh);
        }
        else if (sidenodepair.size() == 4)
        {
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Corner coords\n";
                TPZManVector<int64_t,4> isleft(gel->NNodes());
                for (int i=0; i<gel->NNodes(); i++) {
                    TPZManVector<REAL,3> x(3,0.), jac(3,0.);
                    REAL distance;
                    gel->NodePtr(i)->GetCoordinates(x);
                    isleft[i] = plane.IsLeft(x);
                    distance = plane.Distance(x, jac);
                    sout << x << " distance " << distance << " isleft " << isleft[i] << std::endl;
                }
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            Reorder(gel, targetmesh, sidenodepair);
            
            
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                // identify the X coordinates in the new order
                TPZManVector<TPZManVector<REAL,3>,4> XNodes(4);
                for (int i=0; i<4; i++) {
                    XNodes[i].Resize(3, 0.);
                    targetmesh.NodeVec()[sidenodepair[i].second].GetCoordinates(XNodes[i]);
                    sout << "Node " << i << " coordinate " << XNodes[i] << std::endl;
                }
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            
            TPZManVector<int64_t,4> nodeindices(4);
            for (int i=0; i<4; i++) {
                nodeindices[i] = sidenodepair[i].second;
            }
            // create a quadrilateral or 2 triangles
            TPZGeoEl *gelnew = new TPZGeoElRefLess<pzgeom::TPZGeoQuad>(nodeindices,1,targetmesh);
            CheckElement(gelnew);
        }
        else {
            int centerindex = ReorderGeneral(targetmesh, sidenodepair);
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "cornercoordinates\n";
                int nc = gel->NCornerNodes();
                for (int ic=0; ic<nc; ic++) {
                    TPZManVector<REAL,3> x(3);
                    gel->NodePtr(ic)->GetCoordinates(x);
                    sout << "node " << ic << ' ' << x << std::endl;
                }
                sout << "Intersection coordinates\n";
                int64_t sz = sidenodepair.size();
                // identify the X coordinates in the new order
                TPZManVector<TPZManVector<REAL,3>,6> XNodes(sz);
                for (int64_t i=0; i<sz; i++) {
                    XNodes[i].Resize(3, 0.);
                    targetmesh.NodeVec()[sidenodepair[i].second].GetCoordinates(XNodes[i]);
                    sout << "Node " << i << " side " << sidenodepair[i].first << " coordinate " << XNodes[i] << std::endl;
                }
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            int64_t numnodes = sidenodepair.size();
            for (int64_t in=0; in<numnodes; in++) {
                TPZManVector<int64_t,3> nodeindices(3);
                nodeindices[0] = centerindex;
                nodeindices[1] = sidenodepair[in].second;
                nodeindices[2] = sidenodepair[(in+1)%numnodes].second;
                // we should create an element with 3 nodes
                new TPZGeoElRefLess<pzgeom::TPZGeoTriangle>(nodeindices,1,targetmesh);
            }
            //DebugStop();
            std::cout << "caso com " << sidenodepair.size() << " nos " << std::endl;
        }
    }
    targetmesh.BuildConnectivity();
}

/// Compute the intersection between the edge and the plane
REAL TPZHyperPlaneIntersect::EdgeIntersect(const TPZGeoElSide &gelside, const TPZHyperPlane &plane)
{
    TPZManVector<REAL, 1> edgeparam(1,0.);
    TPZManVector<REAL,3> planejac(3), X(3);
    gelside.X(edgeparam, X);
    REAL distance = plane.Distance(X,planejac);
    int niter = 0;
    while (std::abs(distance) > 1.e-6 && niter < 5) {
        TPZFNMatrix<1,REAL> edgejac(1,1), jacinv(1,1);
        TPZFNMatrix<3,REAL> axes(1,3);
        REAL detjac;
        gelside.Jacobian(edgeparam,edgejac, axes, detjac, jacinv);
        REAL inner = 0;
        for (int i=0; i<3; i++) {
            inner += axes(0,i)*planejac[i];
        }
        edgeparam[0] -= jacinv(0,0)*distance/inner;
        gelside.X(edgeparam, X);
        distance = plane.Distance(X, planejac);
        niter++;
    }
    if (niter == 5) {
        DebugStop();
    }
    return edgeparam[0];
}

void TPZHyperPlaneIntersect::Reorder(TPZGeoEl *gel, TPZGeoMesh &target, TPZVec<std::pair<int64_t,int64_t> > &sidenodepair)
{
    if (sidenodepair.size() != 4) {
        // method is only made for reordering 4 nodes
        // the methodology could be extended - ask me (philippe)
        DebugStop();
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        gel->Print(sout);
        // identify the X coordinates in the new order
        TPZManVector<TPZManVector<REAL,3>,4> XNodes(4);
        for (int i=0; i<4; i++) {
            XNodes[i].Resize(3, 0.);
            target.NodeVec()[sidenodepair[i].second].GetCoordinates(XNodes[i]);
            sout << "Node " << i << " coordinate " << XNodes[i] << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZManVector<std::pair<int64_t, int64_t> > locnodes(4);
    // look for the best permutation of the first three nodes to form a quadrilateral
    for (int perm=0; perm<3; perm++) {
        // permute the nodes
        for (int i=0; i<3; i++) {
            locnodes[(i+perm)%3] = sidenodepair[i];
        }
        locnodes[3] = sidenodepair[3];
        
        // identify the X coordinates in the new order
        TPZManVector<TPZManVector<REAL,3>,4> XNodes(4);
        for (int i=0; i<4; i++) {
            XNodes[i].Resize(3, 0.);
            target.NodeVec()[locnodes[i].second].GetCoordinates(XNodes[i]);
        }
        // identify the vector corresponding to the second and third edge
        TPZManVector<REAL,3> vecedge1(3),vecedge2(3,0.),vecedge3(3,0.), normaltri(3),vecpoint3(3),normaledge3(3);
        vecedge1 = XNodes[1]-XNodes[0];
        vecedge2 = XNodes[2]-XNodes[1];
        vecedge3 = XNodes[0]-XNodes[2];
        vecpoint3 = XNodes[3]-XNodes[2];
        TPZNumeric::ProdVetorial(vecedge2, vecedge3, normaltri);
        TPZNumeric::ProdVetorial(vecedge3, normaltri, normaledge3);
        REAL inner = 0;
        for (int i=0; i<3; i++) {
            inner += vecpoint3[i]*normaledge3[i];
        }
        REAL normvecpoint3 = TPZNumeric::Norm(vecpoint3);
        REAL normvecedge1 = TPZNumeric::Norm(vecedge1);
        if (inner > 1.e-10 || normvecpoint3 < 1.e-8 || normvecedge1 < 1.e-8) {
            // FOUND the orientation!!!
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Permutation = " << perm << std::endl;
                sout << "vecpoint3 " << vecpoint3 << std::endl;
                sout << "vecedge1 " << vecedge1 << std::endl;
                sout << "vecedge2 " << vecedge2 << std::endl;
                sout << "vecedge3 " << vecedge3 << std::endl;
                sout << "normaledge3 " << normaledge3 << std::endl;
                sout << "inner " << inner << " normvecpoint3 " << normvecpoint3 << " normvecedge1 " << normvecedge1 << std::endl;
                sout << "locnodes " << locnodes << std::endl;
                sout << "sidenodepair = " << sidenodepair << std::endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            sidenodepair = locnodes;
            return;
        }
    }
    // couldnt find any
    DebugStop();
}

// verify if the element is not warped
void CheckElement(TPZGeoEl *gel)
{
    TPZManVector<TPZManVector<REAL,3>, 4> normals(4);
    TPZFNMatrix<16,REAL> inner(4,4);
    for (int node=0; node<4; node++) {
        TPZFNMatrix<6,REAL> axes(2,3);
        TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
        TPZTransform<> tr(0,2);
        tr = gel->SideToSideTransform(node, 8);
        TPZManVector<REAL,2> pt0(0),ptnode(2);
        tr.Apply(pt0, ptnode);
        REAL detjac;
        gel->Jacobian(ptnode, jac, axes, detjac, jacinv);
        normals[node].Resize(3,0.);
        if (fabs(detjac) < 1.e-10) {
            continue;
        }
        TPZManVector<REAL,3> ax0(3),ax1(3);
        for (int i=0; i<3; i++) {
            ax0[i] = axes(0,i);
            ax1[i] = axes(1,i);
        }
        TPZNumeric::ProdVetorial(ax0, ax1, normals[node]);
    }
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            inner(i,j)= 0.;
            for (int k=0; k<3; k++) {
                inner(i,j) += normals[i][k]*normals[j][k];
            }
        }
    }
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            if (inner(i,j) < 0.) {
                inner.Print("inner");
                DebugStop();
            }
        }
    }
}

/// a version which allows for more than 4 nodes
int TPZHyperPlaneIntersect::ReorderGeneral(TPZGeoMesh &target, TPZVec<std::pair<int64_t,int64_t> > &sidenodepair)
{
    int64_t numnodes = sidenodepair.size();
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        // identify the X coordinates in the new order
        TPZManVector<TPZManVector<REAL,3>,8> XNodes(numnodes);
        for (int64_t i=0; i<numnodes; i++) {
            XNodes[i].Resize(3, 0.);
            target.NodeVec()[sidenodepair[i].second].GetCoordinates(XNodes[i]);
            sout << "Node " << i << " coordinate " << XNodes[i] << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    // compute the center coordinate
    TPZManVector<REAL,3> center(3,0.);
    TPZManVector<TPZManVector<REAL,3>,8> XNodes(numnodes);
    for (int64_t i=0; i<numnodes; i++) {
        XNodes[i].Resize(3, 0.);
        target.NodeVec()[sidenodepair[i].second].GetCoordinates(XNodes[i]);
        for (int j=0; j<3; j++) {
            center[j] += XNodes[i][j]/numnodes;
        }
    }
    int64_t centerindex = target.NodeVec().AllocateNewElement();
    target.NodeVec()[centerindex].Initialize(center, target);
    // compute the normal to the plane formed by the points
    TPZFNMatrix<100,REAL> normalnorm(numnodes,numnodes,0.);
    int imax=0, jmax = 0, maxnorm = 0.;
    
    for (int64_t in=0; in < numnodes; in++) {
        TPZManVector<REAL,3> veci(3);
        veci = XNodes[in]-center;
        TPZNumeric::NormalizeVetor3(veci);
        for (int64_t jn=0; jn<numnodes; jn++) {
            TPZManVector<REAL,3> vecj(3), vecprod(3);
            vecj = XNodes[jn]-center;
            TPZNumeric::NormalizeVetor3(vecj);
            TPZNumeric::ProdVetorial(veci, vecj, vecprod);
            REAL vecprodnorm = TPZNumeric::Norm(vecprod);
            normalnorm(in,jn) = vecprodnorm;
            if (vecprodnorm > maxnorm) {
                maxnorm = vecprodnorm;
                imax = in;
                jmax = jn;
            }
        }
    }
    // compute the normal with imax,jmax
    TPZManVector<REAL,3> veci(3),vecj(3), vecprod(3);
    veci = XNodes[imax]-center;
    vecj = XNodes[jmax]-center;
    TPZNumeric::NormalizeVetor3(veci);
    TPZNumeric::NormalizeVetor3(vecj);
    TPZNumeric::ProdVetorial(veci, vecj, vecprod);
    TPZNumeric::NormalizeVetor3(vecprod);
    TPZNumeric::ProdVetorial(vecj, vecprod, veci);
    // compute the theta angles for each node
    TPZManVector<REAL,20> angles(numnodes);
    for (int64_t in=0; in<numnodes; in++) {
        REAL x(0.),y(0.);
        TPZManVector<REAL,3> vecnod(3);
        vecnod = XNodes[in]-center;
        for (int i=0; i<3; i++) {
            x += veci[i]*vecnod[i];
            y += vecj[i]*vecnod[i];
        }
        angles[in] = atan2(x, y);
    }
    
    // ordena os nos de acordo com os angulos
    std::multimap<REAL, std::pair<int64_t,int64_t> > ordered;
    for (int64_t in=0; in<numnodes; in++) {
        ordered.insert(std::pair<REAL,std::pair<int64_t,int64_t> >(angles[in],sidenodepair[in]));
    }
    std::multimap<REAL, std::pair<int64_t,int64_t> >::iterator it;
    int64_t count = 0;
    for (it=ordered.begin(); it!= ordered.end(); it++) {
        sidenodepair[count] = it->second;
        count++;
    }
    return centerindex;
}
