//
//  TPZFracSimulation.cpp
//  PZ
//
//  Created by Philippe Devloo on 26/09/17.
//

#include "TPZFracSimulation.h"

#include "mixedpoisson.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZMixedPoissonParabolic.h"
#include "pzbndcond.h"

#include "TPZVTKGeoMesh.h"

#include "TPZRefPattern.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.fracsimulation"));
#endif


/// Build an MHM object with information given by the root name
void TPZFracSimulation::ReadDataFile(const std::string &rootname)
{
    std::string datafilename = rootname + ".data";
    std::ifstream datafile(datafilename);
    ReadDataFile(datafile);
    std::string meshfilename = rootname + ".msh";
    fGmsh.GeometricGmshMesh(meshfilename,fMHM->GMesh().operator->());
    AdjustGeometricMesh(rootname);
    // in order to build the wrap mesh
    fMHM->DivideSkeletonElements(0);
#ifdef PZDEBUG
    {
        std::string filename = rootname + "_gmesh_wrap.vtk";
        std::ofstream out(filename);
        TPZVTKGeoMesh::PrintGMeshVTK(fMHM->GMesh(), out);
    }
#endif
}

/// read the data in the data file and put the information in right places
// creates material objects for the computational meshes
void TPZFracSimulation::ReadDataFile(std::ifstream &input)
{
    ReadPreamble(input);
    ReadFractures(input);
}

static void ReadNextLine(std::ifstream &input, std::string &line)
{
    std::getline(input,line);
    while (input && line[0] == '#') {
        std::getline(input,line);
    }
}

/// reads the preamble of the data file
void TPZFracSimulation::ReadPreamble(std::ifstream &input)
{

    std::string line;
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fFracSet.fLowLeft[0] >> fFracSet.fLowLeft[1] >> fFracSet.fTopRight[0] >> fFracSet.fTopRight[1];
    }
    REAL delxmin = min(fFracSet.fTopRight[0]-fFracSet.fLowLeft[0],fFracSet.fTopRight[1]-fFracSet.fLowLeft[1]);
    REAL tol = delxmin/200.; // will generate a 200x200 raster for points
    fFracSet.SetTol(tol);
    ReadNextLine(input, line);
    int numMHM;
    {    std::istringstream stin(line);
        stin >> numMHM;
    }
    fFracSet.SetMHMSpacing(numMHM);
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fFracSet.fElementSize >> fFracSet.fMinElementSize;
    }
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fSimulationType;
        if(fSimulationType != 0 && fSimulationType != 1) DebugStop();
    }
    ReadNextLine(input, line);
    int nummat;
    int matidcounter = 20;
    {
        std::istringstream stin(line);
        stin >> nummat;
    }
    fGmsh.fPZMaterialId[1]["MHMLine"] = fMHM->fSkeletonMatId;
    std::string planemat;
    for (int i=0; i<nummat; i++) {
        ReadNextLine(input, line);
        std::string matname , mattype;
        int matid, dimension;
        REAL density, perm;
        {
            std::istringstream stin(line);
            stin >> matname >> dimension >> matid >> mattype >> density >> perm;
        }
        matname.erase(0,1);
        matname.erase(matname.end()-1,matname.end());
        if (mattype == "flow" && dimension == 2) {
            fFracSet.fPhysicalname = matname;
            fGmsh.fPZMaterialId[dimension][matname] = matidcounter;
            InsertDarcyMaterial(matidcounter, perm, density);
            fMaterialIds[matname] = matidcounter;
            matidcounter++;
        }
        if (mattype == "boundary" && dimension == 1)
        {
            fGmsh.fPZMaterialId[dimension][matname] = matidcounter;
            InsertDarcyBCMaterial(matidcounter, dimension, (int)(density+0.5), perm);
            fMaterialIds[matname] = matidcounter;
            matidcounter++;
        }
    }
    ReadNextLine(input, line);
    
    int numstepdefinitions = 0;
    {
        std::istringstream stin(line);
        stin >> numstepdefinitions;
    }
    fTimeSteps.Resize(numstepdefinitions);
    for (int ist=0; ist<numstepdefinitions; ist++) {
        ReadNextLine(input, line);
        std::istringstream stin(line);
        stin >> fTimeSteps[ist].first >> fTimeSteps[ist].second;
    }
    
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fPostProcessRootname;
    }
    
    int numpostprocess;
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> numpostprocess;
    }
    fPostProcnames.Resize(numpostprocess);
    for (int pp=0; pp<numpostprocess; pp++) {
        ReadNextLine(input, line);
        std::istringstream stin(line);
        stin >> fPostProcnames[pp];
    }

    ReadNextLine(input, line);
    while (line.find("end") == std::string::npos) {
        ReadNextLine(input, line);
        if (!input) {
            std::cout << "No end found in the preamble of the file\n";
            DebugStop();
        }
    }
}



#include "pzvec_extras.h"

TPZTransform GetTransform(std::set<int64_t> &nodes, TPZGeoMesh *gmesh)
{
    TPZTransform result(1,3);
    TPZManVector<int64_t,5> nodevec(nodes.size());
    int i=0;
    for (auto it = nodes.begin() ; it != nodes.end(); it++) {
        nodevec[i++] = *it;
    }
    TPZFNMatrix<25> dist(nodes.size(),nodes.size(),0.);
    for (int i=0; i<nodevec.size(); i++) {
        for (int j=0; j<nodevec.size(); j++) {
            dist(i,j) = 0;
            int64_t in = nodevec[i];
            int64_t jn = nodevec[j];
            for (int c=0; c<3; c++) {
                dist(i,j) += (gmesh->NodeVec()[in].Coord(c)-gmesh->NodeVec()[jn].Coord(c))*
                (gmesh->NodeVec()[in].Coord(c)-gmesh->NodeVec()[jn].Coord(c));
            }
            dist(i,j) = sqrt(dist(i,j));
        }
    }
    REAL maxdist = 0.;
    int jmax = 0;
    int imax = 0;
    for (int i=0; i<nodevec.size(); i++) {
        for (int j=0; j<nodevec.size(); j++) {
            if (dist(i,j) > maxdist) {
                maxdist = dist(i,j);
                imax = i;
                jmax = j;
            }
        }
    }
    TPZManVector<REAL,3> ivec(3),jvec(3),vector(3);
    gmesh->NodeVec()[nodevec[imax]].GetCoordinates(ivec);
    gmesh->NodeVec()[nodevec[jmax]].GetCoordinates(jvec);
    vector = jvec-ivec;
    REAL sum = 0.;
    for (int i=0; i<3; i++) {
        result.Mult()(0,i) = vector[i]*2./(maxdist*maxdist);
        sum += ivec[i]*result.Mult()(0,i);
    }
    result.Sum()(0,0) = -1.-sum;
#ifdef PZDEBUG
    TPZManVector<REAL,2> itr(1),jtr(1);
    result.Apply(ivec, itr);
    result.Apply(jvec, jtr);
    if (abs(itr[0]+1.) > 1.e-6 || abs(jtr[0]-1.) > 1.e-6) {
        DebugStop();
    }
#endif
    return result;
}

int64_t GroupElements(TPZStack<int64_t> &elems, TPZGeoMesh *gmesh)
{
    std::set<int64_t> nodes;
    for (int el=0; el<elems.size(); el++) {
        TPZGeoEl *gel = gmesh->Element(elems[el]);
        int nnodes = gel->NCornerNodes();
        for (int in=0; in<nnodes; in++) {
            nodes.insert(gel->NodeIndex(in));
        }
    }
    TPZTransform tr = GetTransform(nodes, gmesh);
    std::map<int64_t,int> nodemap;
    int nodemin = -1, nodemax = -1;
    {
        REAL minco = 0.;
        REAL maxco = 0.;
        int i=0;
        for (auto it=nodes.begin(); it != nodes.end(); it++) {
            TPZManVector<REAL,3> coord(3), cotr(1);
            gmesh->NodeVec()[*it].GetCoordinates(coord);
            tr.Apply(coord, cotr);
            if (cotr[0] < minco) {
                minco = cotr[0];
                nodemin = *it;
            }
            if (cotr[0] > maxco) {
                maxco = cotr[0];
                nodemax = *it;
            }
            nodemap[*it] = i++;
        }
    }
    TPZAutoPointer<TPZRefPattern> refpat;
    {
        TPZGeoMesh refpatmesh;
        refpatmesh.NodeVec().Resize(nodemap.size());
        for (auto it=nodemap.begin(); it != nodemap.end(); it++)
        {
            TPZManVector<REAL,3> co(3);
            gmesh->NodeVec()[it->first].GetCoordinates(co);
            refpatmesh.NodeVec()[it->second].Initialize(co, refpatmesh);
        }
        TPZManVector<int64_t,2> nodeindices(2);
        nodeindices[0] = nodemap[nodemin];
        nodeindices[1] = nodemap[nodemax];
        int64_t index;
        refpatmesh.CreateGeoElement(EOned, nodeindices, 1, index);
        for (int64_t el=0; el<elems.size(); el++) {
            TPZGeoEl *gel = gmesh->Element(elems[el]);
            nodeindices[0] = nodemap[gel->NodeIndex(0)];
            nodeindices[1] = nodemap[gel->NodeIndex(1)];
            refpatmesh.CreateGeoElement(EOned, nodeindices, 1, index);
            TPZGeoEl *gelsub = refpatmesh.Element(index);
            gelsub->SetFather((int64_t)0);
        }
        refpatmesh.BuildConnectivity();
        refpat = new TPZRefPattern(refpatmesh);
        for (int64_t el=0; el<refpatmesh.NElements(); el++) {
            refpatmesh.Element(el)->SetFather(-1);
        }
    }
    TPZManVector<int64_t,2> nodeindices(2);
    int64_t index;
    nodeindices[0] = nodemin;
    nodeindices[1] = nodemax;
    int matid = gmesh->Element(elems[0])->MaterialId();
    gmesh->CreateGeoElement(EOned, nodeindices, matid, index);
    gmesh->Element(index)->SetRefPattern(refpat);
    for (int el = 0; el<elems.size(); el++) {
        gmesh->Element(index)->SetSubElement(el, gmesh->Element(elems[el]));
        gmesh->Element(elems[el])->SetFather(gmesh->Element(index));
    }
    return index;
}

void CreateRefPatterns(TPZGeoMesh *gmesh, TPZVec<int64_t> &elemententity, std::set<int> &matids)
{
    int64_t nel = gmesh->NElements();
    int64_t el =0;
    while(el<nel)
    {
        while(el<nel && gmesh->Element(el)->Dimension() != gmesh->Dimension()-1) el++;
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if (matids.find(matid) == matids.end()) {
            el++;
            continue;
        }
        if(el == nel) break;
        TPZStack<int64_t> elstack;
        elstack.Push(el);
        int64_t ident = elemententity[el];
        el++;
        while(el < nel && elemententity[el] == ident)
        {
            elstack.Push(el);
            el++;
        }
        if (elstack.size() > 1) {
            int64_t index = GroupElements(elstack, gmesh);
            if (elemententity.size() <= index) {
                elemententity.Resize(index+1, -1);
            }
            elemententity[index] = ident;
        }
    }
    gmesh->BuildConnectivity();
}

/// verify is each macro element links only two subdomains
void VerifySubdomainConsistency(TPZGeoMesh *gmesh, TPZVec<int64_t> &elemententity, std::set<int> &matids)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if (gel->Dimension() != gmesh->Dimension()-1 || matids.find(matid) == matids.end()) {
            continue;
        }
        if (!gel->HasSubElement()) {
            continue;
        }
        std::set<int64_t> connected;
        int nsubel = gel->NSubElements();
        for (int isub = 0; isub<nsubel; isub++) {
            TPZGeoEl *subel = gel->SubElement(isub);
            TPZGeoElSide subelside(subel,subel->NSides()-1);
            TPZGeoElSide neighbour = subelside.Neighbour();
            while (neighbour != subelside) {
                if (neighbour.Element()->Dimension() == gmesh->Dimension()) {
                    connected.insert(elemententity[neighbour.Element()->Index()]);
                }
                neighbour = neighbour.Neighbour();
            }
        }
        if (connected.size() != 2) {
            std::cout << "Coordinates of the element\n";
            TPZManVector<REAL,3> co0(3), co1(3);
            gel->Node(0).GetCoordinates(co0);
            gel->Node(1).GetCoordinates(co1);
            std::cout << co0 << " " << co1 << std::endl;
            DebugStop();
        }
    }
}

/// Build skeleton data structure
void BuildSkeleton(TPZGeoMesh *gmesh, TPZVec<int64_t> &elemententity, std::set<int> &matids, std::map<int64_t,std::pair<int64_t,int64_t>> &skeleton)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if (gel->Dimension() != gmesh->Dimension()-1 || matids.find(matid) == matids.end()) {
            continue;
        }
        if (gel->Father()) {
            continue;
        }
        std::set<int64_t> connected;
        TPZStack<TPZGeoElSide> elstack;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        elstack.Push(gelside);
        while(elstack.size())
        {
            TPZGeoElSide search = elstack.Pop();
            TPZGeoElSide neighbour = search.Neighbour();
            while (neighbour != search)
            {
                if (neighbour.Element()->Dimension() == gmesh->Dimension())
                {
                    connected.insert(elemententity[neighbour.Element()->Index()]);
                }
                neighbour = neighbour.Neighbour();
            }
            if (search.Element()->HasSubElement())
            {
                int nsubel = search.Element()->NSubElements();
                for (int isub = 0; isub<nsubel; isub++)
                {
                    TPZGeoEl *subel = search.Element()->SubElement(isub);
                    TPZGeoElSide subelside(subel,subel->NSides()-1);
                    elstack.Push(subelside);
                }
            }
        }
#ifdef PZDEBUG
        if (connected.size() > 2)
        {
            std::cout << "Coordinates of the element\n";
            TPZManVector<REAL,3> co0(3), co1(3);
            gel->Node(0).GetCoordinates(co0);
            gel->Node(1).GetCoordinates(co1);
            std::cout << co0 << " " << co1 << std::endl;
            DebugStop();
        }
#endif
        if(connected.size() == 2)
        {
            auto it = connected.begin();
            int64_t domain1 = *it;
            it++;
            int64_t domain2 = *it;
            skeleton[el] = std::pair<int64_t,int64_t>(domain1,domain2);
        } else if (connected.size() == 1)
        {
            auto it = connected.begin();
            int64_t domain1 = *it;
            int64_t domain2 = el;
            skeleton[el] = std::pair<int64_t, int64_t>(domain1,domain2);
        }
        else
        {
            DebugStop();
        }
    }
}

/// Adjust the entity index of the fracture elemenents
void AdjustEntityOfFractures(TPZGeoMesh *gmesh,TPZVec<int64_t> &EntityIndex, std::set<int> fracmatid)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (fracmatid.find(gel->MaterialId()) != fracmatid.end()) {
            std::set<int64_t> entities;
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == gmesh->Dimension()) {
                    entities.insert(EntityIndex[neighbour.Element()->Index()]);
                }
                neighbour = neighbour.Neighbour();
            }
            if (entities.size() != 1) {
                DebugStop();
            }
            EntityIndex[gel->Index()] = *entities.begin();
        }
    }
}


/// creates and inserts a Darcy or ParabolicDarcy object in the mesh
void TPZFracSimulation::InsertDarcyMaterial(int matid, REAL permeability, REAL rho)
{
#ifdef PZDEBUG
    if (fMHM->fMaterialIds.find(matid) != fMHM->fMaterialIds.end()) {
        DebugStop();
    }
#endif
    int dimension = 2;
    TPZMixedPoisson * mat = new TPZMixedPoisson(matid,dimension);
    //    TPZMixedPoissonParabolic *mat = new TPZMixedPoissonParabolic(matid,dimension);
    //    mat->SetDeltaT(1000.);
    mat->SetSymmetric();
    mat->SetPermeability(permeability);
    
    fMHM->CMesh()->InsertMaterialObject(mat);
    
    TPZVecL2 *vecmat = new TPZVecL2(matid);
    fMHM->FluxMesh()->InsertMaterialObject(vecmat);
    TPZMatLaplacian *presmat = new TPZMatLaplacian(matid);
    presmat->SetDimension(dimension);
    fMHM->PressureMesh()->InsertMaterialObject(presmat);
    
    fMHM->fMaterialIds.insert(matid);
    
}

/// creates and inserts the boundary condition objects
void TPZFracSimulation::InsertDarcyBCMaterial(int matid, int dimension, int bctype, REAL val)
{
#ifdef PZDEBUG
    if (fMHM->fMaterialBCIds.find(matid) != fMHM->fMaterialBCIds.end()) {
        DebugStop();
    }
#endif
    int rootmat = *fMHM->fMaterialIds.begin();
    TPZMaterial *mat = fMHM->CMesh()->FindMaterial(rootmat);
    if (!mat) {
        DebugStop();
    }
    TPZFNMatrix<1,STATE> val1(1,1,0.),val2(1,1,val);
    TPZBndCond * bcmat = new TPZBndCond(mat,matid,bctype,val1,val2);
    //    TPZMixedPoissonParabolic *mat = new TPZMixedPoissonParabolic(matid,dimension);
    //    mat->SetDeltaT(1000.);
    
    fMHM->CMesh()->InsertMaterialObject(bcmat);
    
    TPZMaterial *fluxmat = fMHM->FluxMesh()->FindMaterial(rootmat);
    TPZBndCond *bnd = new TPZBndCond(fluxmat,matid,bctype,val1,val2);
    fMHM->FluxMesh()->InsertMaterialObject(bnd);
    
    fMHM->fMaterialBCIds.insert(matid);
    
}


void TPZFracSimulation::ReadFractures(std::ifstream &input)
{
    std::string line;
    int id = 0;
    while(input)
    {
        if(std::getline(input, line))
        {
            if (line[0] == '#') {
                continue;
            }
            if(line.find("CORNER") != std::string::npos)
            {
                TPZGeoNode first, second;
                TPZManVector<REAL,3> x0(3), x1(3);
                std::istringstream sin(line);
                std::string corner;
                sin >> corner >> x0[0] >> x0[1] >> x0[2] >> x1[0] >> x1[1] >> x1[2];
                first.SetCoord(x0);
                second.SetCoord(x1);
                int64_t index0 = fFracSet.InsertNode(first);
                int64_t index1 = fFracSet.InsertNode(second);
                first = fFracSet.fNodeVec[index0];
                second = fFracSet.fNodeVec[index1];
                uint64_t keyfirst = fFracSet.GetLoc(first);
                uint64_t keysecond = fFracSet.GetLoc(second);
                if (fFracSet.fPointMap.find(keyfirst) == fFracSet.fPointMap.end()) {
                    DebugStop();
                }
                if (fFracSet.GetLoc(first) != keyfirst) {
                    DebugStop();
                }
                if (fFracSet.fPointMap.find(keysecond) == fFracSet.fPointMap.end()) {
                    DebugStop();
                }
                if (fFracSet.GetLoc(second) != keysecond) {
                    DebugStop();
                }
                if (first.Coord(0) < second.Coord(0)) {
                    int matid = fFracSet.matid_internal_frac;
                    TPZFracture frac(id,matid,index0,index1);
                    int64_t index = fFracSet.fFractureVec.AllocateNewElement();
                    fFracSet.fFractureVec[index] = frac;
                }
                else
                {
                    int matid = fFracSet.matid_internal_frac;
                    TPZFracture frac(id,matid,index1,index0);
                    int64_t index = fFracSet.fFractureVec.AllocateNewElement();
                    fFracSet.fFractureVec[index] = frac;
                }
            }
            {
                uint64_t pos = line.find("PERM");
                if(pos != std::string::npos)
                {
                    int64_t lastfrac = fFracSet.fFractureVec.NElements()-1;
                    if(lastfrac < 0) DebugStop();
                    std::string sub = line.substr(pos+4,std::string::npos);
                    std::istringstream sin(sub);
                    sin >> fFracSet.fFractureVec[lastfrac].fFracPerm;
                }
            }
            {
                uint64_t pos = line.find("NAME");
                if(pos != std::string::npos)
                {
                    int64_t lastfrac = fFracSet.fFractureVec.NElements()-1;
                    if (lastfrac < 0) {
                        DebugStop();
                    }
                    std::string sub = line.substr(pos+4,std::string::npos);
                    std::istringstream sin(sub);
                    std::string fracname;
                    sin >> fracname;
                    fracname.erase(0,1);
                    fracname.erase(fracname.end()-1,fracname.end());
                    fFracSet.fFractureVec[lastfrac].fPhysicalName = fracname;
                }
            }
        }
        else
        {
            break;
        }
    }
    int64_t nfrac = fFracSet.fFractureVec.NElements();
    int matid = *(fMHM->fMaterialBCIds.rbegin())+1;
    matid = (matid+10)-matid%10;
    for (int ifr=0; ifr<nfrac; ifr++) {
        int meshdim = fMHM->GMesh()->Dimension();
        fMHM->InsertFractureFlowMaterial(matid);
        // Material medio poroso
        TPZMixedPoisson * mat = new TPZMixedPoisson(matid,meshdim-1);
        REAL perm = fFracSet.fFractureVec[ifr].fFracPerm;
        std::string matname = fFracSet.fFractureVec[ifr].fPhysicalName;
        fGmsh.fPZMaterialId[meshdim-1][matname] = matid;
        fMaterialIds[matname] = matid;
        mat->SetSymmetric();
        mat->SetPermeability(perm);
        TPZFNMatrix<9,REAL> K(3,3,0.),KInv(3,3,0.);
        K(0,0) = perm;
        K(1,1) = perm;
        KInv(0,0) = 1./K(0,0);
        KInv(1,1) = 1./K(1,1);
        mat->SetPermeabilityTensor(K, KInv);
        fMHM->CMesh()->InsertMaterialObject(mat);
        TPZVecL2 *vecl21 = new TPZVecL2(matid);
        fMHM->FluxMesh()->InsertMaterialObject(vecl21);
        TPZMat1dLin *mat1d = new TPZMat1dLin(matid);
        fMHM->PressureMesh()->InsertMaterialObject(mat1d);
        matname = matname + "_MHM";
        matid++;
        mat = new TPZMixedPoisson(*mat);
        mat->SetId(matid);
        fGmsh.fPZMaterialId[meshdim-1][matname] = matid;
        fMaterialIds[matname] = matid;
        fMHM->CMesh()->InsertMaterialObject(mat);
        TPZVecL2 *vecl2 = new TPZVecL2(matid);
        fMHM->FluxMesh()->InsertMaterialObject(vecl2);
        TPZMat1dLin *mat1d2 = new TPZMat1dLin(matid);
        fMHM->PressureMesh()->InsertMaterialObject(mat1d2);

        fMHM->fSkeletonWithFlowMatId.insert(matid);
        matid++;
    }
}

/// adjust the geometric element read from gmesh
void TPZFracSimulation::AdjustGeometricMesh(const std::string &rootname)
{
    {
        std::string filename = rootname + "_gmesh.vtk";
        std::ofstream out(filename);
        TPZVTKGeoMesh::PrintGMeshVTK(fMHM->GMesh().operator->(), out);
    }
    
    TPZGeoMesh *gmesh = fMHM->GMesh().operator->();
    std::set<int> matids = fMHM->fSkeletonWithFlowMatId;
    matids.insert(fMHM->fSkeletonMatId);
    CreateRefPatterns(gmesh, fGmsh.fEntityIndex, matids);
    VerifySubdomainConsistency(gmesh, fGmsh.fEntityIndex, matids);
    matids = fMHM->fFractureFlowDim1MatId;
    AdjustEntityOfFractures(gmesh,fGmsh.fEntityIndex,matids);
    std::map<int64_t,std::pair<int64_t,int64_t>> skeletonstruct;
    matids = fMHM->fSkeletonWithFlowMatId;
    matids.insert(fMHM->fSkeletonMatId);
    std::copy(fMHM->fMaterialBCIds.begin(),fMHM->fMaterialBCIds.end(), std::inserter(matids, matids.begin()));
    BuildSkeleton(gmesh, fGmsh.fEntityIndex, matids, skeletonstruct);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        gmesh->Print(sout);
        sout << "Entity index \n" << fGmsh.fEntityIndex << endl;
        for (auto it = skeletonstruct.begin(); it != skeletonstruct.end(); it++) {
            sout << "skeleton geo " << it->first << " linked to subdomains " << it->second.first << " " << it->second.second << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef PZDEBUG
    {
        std::string filename = rootname + ".gmesh_mhm.vtk";
        std::ofstream out(filename);
        TPZVTKGeoMesh::PrintGMeshVTK(fMHM->GMesh(), out);
    }
#endif
    fMHM->DefinePartition(fGmsh.fEntityIndex, skeletonstruct);
}



