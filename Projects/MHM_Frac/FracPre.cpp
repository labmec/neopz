//
//  FracPre.cpp
//  PZ
//
//  Created by Philippe Devloo on 07/09/17.
//
//

#include <stdio.h>
#include <list>
#include <fstream>
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzgnode.h"
#include "pzgeoel.h"
#include "TPZFracSet.h"


void ReadFracDefinition(const std::string &filename,TPZFracSet &fracset);

void ExportGMsh(TPZFracSet &fracset, const std::string &filename);

int main()
{
    TPZFracSet fracset;
    ReadFracDefinition("HorizontalFrac.data", fracset);
//    fracset.SetTol(5);
//    fracset.SetMHMSpacing(11);
    fracset.ComputeFractureIntersections();
    fracset.CleanupFractureNetwork();
    fracset.SplitFracturesByMHM();
    fracset.CleanupFractureNetwork();
    fracset.CheckPointMapConsistency();
    fracset.AddMHMNodes();
    fracset.CheckPointMapConsistency();
//    int nelem_MHMside = 2;
//    REAL mhm_size = fracset.fMHMSpacing[0];
//    REAL elem_size = mhm_size/nelem_MHMside;
    fracset.ComputeMeshSizeAtNodes();

    ExportGMsh(fracset,"HorizontalFrac.geo");
    
    return 0;
}

void ReadNextLine(std::ifstream &input, std::string &line)
{
    std::getline(input,line);
    while (input && line[0] == '#') {
        std::getline(input,line);
    }
}

void ReadPreamble(std::ifstream &input, TPZFracSet &fracset)
{
    std::string line;
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fracset.fLowLeft[0] >> fracset.fLowLeft[1] >> fracset.fTopRight[0] >> fracset.fTopRight[1];
    }
    REAL delxmin = min(fracset.fTopRight[0]-fracset.fLowLeft[0],fracset.fTopRight[1]-fracset.fLowLeft[1]);
    REAL tol = delxmin/200.; // will generate a 200x200 raster for points
    fracset.SetTol(tol);
    ReadNextLine(input, line);
    int numMHM;
    {    std::istringstream stin(line);
        stin >> numMHM;
    }
    fracset.SetMHMSpacing(numMHM);
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fracset.fElementSize >> fracset.fMinElementSize;
    }
    ReadNextLine(input, line);
    int simulationtype;
    {
        std::istringstream stin(line);
        stin >> simulationtype;
    }
    ReadNextLine(input, line);
    int numbermaterials;
    {
        std::istringstream stin(line);
        stin >> numbermaterials;
    }
    std::string planemat;
    for (int i=0; i<numbermaterials; i++) {
        ReadNextLine(input, line);
        std::string matname , mattype;
        int matid, dimension;
        REAL density, perm;
        {
            std::istringstream stin(line);
            stin >> matname >> dimension >> matid >> mattype >> density >> perm;
        }
        if (mattype == "flow" && dimension == 2) {
            matname.erase(0,1);
            matname.erase(matname.end()-1,matname.end());
            fracset.fPhysicalname = matname;
        }
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

void ReadFracDefinition(const std::string &filename,TPZFracSet &fracset)
{
    std::ifstream input(filename);
    ReadPreamble(input, fracset);
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
                int64_t index0 = fracset.InsertNode(first);
                int64_t index1 = fracset.InsertNode(second);
                first = fracset.fNodeVec[index0];
                second = fracset.fNodeVec[index1];
                uint64_t keyfirst = fracset.GetLoc(first);
                uint64_t keysecond = fracset.GetLoc(second);
                if (fracset.fPointMap.find(keyfirst) == fracset.fPointMap.end()) {
                    DebugStop();
                }
                if (fracset.GetLoc(first) != keyfirst) {
                    DebugStop();
                }
                if (fracset.fPointMap.find(keysecond) == fracset.fPointMap.end()) {
                    DebugStop();
                }
                if (fracset.GetLoc(second) != keysecond) {
                    DebugStop();
                }
                if (first.Coord(0) < second.Coord(0)) {
                    int matid = fracset.matid_internal_frac;
                    TPZFracture frac(id,matid,index0,index1);
                    int64_t index = fracset.fFractureVec.AllocateNewElement();
                    fracset.fFractureVec[index] = frac;
                }
                else
                {
                    int matid = fracset.matid_internal_frac;
                    TPZFracture frac(id,matid,index1,index0);
                    int64_t index = fracset.fFractureVec.AllocateNewElement();
                    fracset.fFractureVec[index] = frac;
                }
            }
            {
                uint64_t pos = line.find("PERM");
                if(pos != std::string::npos)
                {
                    int64_t lastfrac = fracset.fFractureVec.NElements()-1;
                    if(lastfrac < 0) DebugStop();
                    std::string sub = line.substr(pos+4,std::string::npos);
                    std::istringstream sin(sub);
                    sin >> fracset.fFractureVec[lastfrac].fFracPerm;
                }
            }
            {
                uint64_t pos = line.find("NAME");
                if(pos != std::string::npos)
                {
                    int64_t lastfrac = fracset.fFractureVec.NElements()-1;
                    if (lastfrac < 0) {
                        DebugStop();
                    }
                    std::string sub = line.substr(pos+4,std::string::npos);
                    std::istringstream sin(sub);
                    std::string fracname;
                    sin >> fracname;
                    fracname.erase(0,1);
                    fracname.erase(fracname.end()-1,fracname.end());
                    fracset.fFractureVec[lastfrac].fPhysicalName = fracname;
                }
            }
        }
        else
        {
            break;
        }
        
    }
}

void PreambleGMsh(std::ofstream &out)
{
    out << "SetFactory(\"OpenCASCADE\");\n\n\n";
}

void ExportPoints(TPZFracSet &fracset, std::ofstream &out)
{
    int64_t np = fracset.fNodeVec.NElements();
    for (int64_t in = 0; in<np; in++) {
        TPZManVector<REAL, 3> co(3,0.);
        fracset.fNodeVec[in].GetCoordinates(co);
        out << "Point(" << in << ") = {" << co[0] << ',' << co[1] << ',' << co[2] << ',' << fracset.fMeshSizeAtNodes[in] << "};\n";
    }
}

void ExportFractures(TPZFracSet &fracset, std::ofstream &out)
{
    int64_t nfrac = fracset.fFractureVec.NElements();
    for (int64_t ifr = 0; ifr <nfrac; ifr++) {
        out << "Line(" << ifr << ") = {" << fracset.fFractureVec[ifr].fNodes[0] << ", " << fracset.fFractureVec[ifr].fNodes[1] << "};\n";
    }
}

void ComputeLines(TPZFracSet &fracset, TPZManVector<int> &firsthor, TPZManVector<int> &firstver, int face, TPZVec<int> &lines)
{
    int numMHMx = (fracset.fTopRight[0]-fracset.fLowLeft[0])/fracset.fMHMSpacing[0];
    int numMHMy = (fracset.fTopRight[1]-fracset.fLowLeft[1])/fracset.fMHMSpacing[1];
    int facey = face/numMHMx;
    int facex = face%numMHMx;
    int ipx[2] = {0,0};
    for (auto it = fracset.fHorizontalLines[facey].begin(); it != fracset.fHorizontalLines[facey].end(); it++) {
        std::pair<uint32_t, uint32_t> key = fracset.NodeKey(it->second);
        if (key.first < facex*fracset.fMHMSpacingInt[0]) {
            ipx[0]++;
        }
        if (key.first < (facex+1)*fracset.fMHMSpacingInt[0]) {
            ipx[1]++;
        }
        if (key.first >= (facex+1)*fracset.fMHMSpacingInt[0]) {
            break;
        }
        
    }
    int ipxp1[2] = {0,0};
    for (auto it = fracset.fHorizontalLines[facey+1].begin(); it != fracset.fHorizontalLines[facey+1].end(); it++) {
        std::pair<uint32_t, uint32_t> key = fracset.NodeKey(it->second);
        if (key.first < facex*fracset.fMHMSpacingInt[0]) {
            ipxp1[0]++;
        }
        if (key.first < (facex+1)*fracset.fMHMSpacingInt[0]) {
            ipxp1[1]++;
        }
        if (key.first >= (facex+1)*fracset.fMHMSpacingInt[0]) {
            break;
        }
    }
    int ipy[2] = {0,0};
    for (auto it = fracset.fVerticalLines[facex].begin(); it != fracset.fVerticalLines[facex].end(); it++) {
        std::pair<uint32_t, uint32_t> key = fracset.NodeKey(it->second);
        if (key.second < facey*fracset.fMHMSpacingInt[1]) {
            ipy[0]++;
        }
        if (key.second < (facey+1)*fracset.fMHMSpacingInt[1]) {
            ipy[1]++;
        }
        if (key.second >= (facey+1)*fracset.fMHMSpacingInt[1]) {
            break;
        }
    }
    int ipyp1[2] = {0,0};
    for (auto it = fracset.fVerticalLines[facex+1].begin(); it != fracset.fVerticalLines[facex+1].end(); it++) {
        std::pair<uint32_t, uint32_t> key = fracset.NodeKey(it->second);
        if (key.second < facey*fracset.fMHMSpacingInt[1]) {
            ipyp1[0]++;
        }
        if (key.second < (facey+1)*fracset.fMHMSpacingInt[1]) {
            ipyp1[1]++;
        }
        if (key.second >= (facey+1)*fracset.fMHMSpacingInt[1]) {
            break;
        }
    }
    int nlines = ipx[1]+ipxp1[1]+ipy[1]+ipyp1[1]-(ipx[0]+ipxp1[0]+ipy[0]+ipyp1[0]);
    int linecount = 0;
    lines.Resize(nlines);
    int firstline;
    firstline = firsthor[facey];
    for (int i = ipx[0]; i<ipx[1]; i++) {
        lines[linecount++] = firstline+i;
    }
    firstline = firstver[facex+1];
    for (int i = ipyp1[0]; i < ipyp1[1]; i++) {
        lines[linecount++] = firstline+i;
    }
    firstline = firsthor[facey+1];
    for (int i = ipxp1[1]-1; i >= ipxp1[0]; i--) {
        lines[linecount++] = -(firstline+i);
    }
    firstline = firstver[facex];
    for (int i = ipy[1]-1; i >= ipy[0]; i--) {
        lines[linecount++] = -(firstline+i);
    }
}
void CreateGMshFaces(TPZFracSet &fracset, std::ofstream &out)
{
    int numMHMx = (fracset.fTopRight[0]-fracset.fLowLeft[0])/fracset.fMHMSpacing[0];
    int numMHMy = (fracset.fTopRight[1]-fracset.fLowLeft[1])/fracset.fMHMSpacing[1];
    int totalmhmlines = 0;
    int64_t nhor = fracset.fHorizontalLines.size();
    for (int64_t hor = 0; hor<nhor; hor++) {
        totalmhmlines += fracset.fHorizontalLines[hor].size()-1;
    }
    int64_t nver = fracset.fVerticalLines.size();
    for (int64_t ver = 0; ver<nver; ver++) {
        totalmhmlines += fracset.fVerticalLines[ver].size()-1;
    }
    int64_t firstMHM = fracset.fFractureVec.NElements()-totalmhmlines;
    TPZManVector<int> FirstHorizontal(numMHMy+2);
    TPZManVector<int> FirstVertical(numMHMx+2);
    
    int counter = firstMHM;
    for (int hor=0; hor<numMHMy+1; hor++) {
        FirstHorizontal[hor] = counter;
        counter += fracset.fHorizontalLines[hor].size()-1;
    }
    for (int ver=0; ver<numMHMx+1; ver++) {
        FirstVertical[ver] = counter;
        counter += fracset.fVerticalLines[ver].size()-1;
    }
    int nfaces = numMHMx*numMHMy;
    for (int iface = 0; iface<nfaces; iface++) {
        TPZManVector<int> lines;
        ComputeLines(fracset, FirstHorizontal, FirstVertical, iface, lines);
        out << "Line Loop(" << iface << ") = {";
        for (int i=0; i<lines.size(); i++) {
            out << lines[i];
            if(i<lines.size()-1) out << ", ";
        }
        out << "};\n";
        out << "Plane Surface(" << iface << ") = {" << iface << "};\n";
//        std::cout << "face " << iface << " lines " << lines << std::endl;
    }
}

void FracturesInFaces(TPZFracSet &fracset, std::ofstream &out)
{
    int64_t nfrac = fracset.fFractureVec.NElements();
    int maxfrac = 10000;
    nfrac = maxfrac < nfrac ? maxfrac : nfrac;
    for (int64_t ifr = 0; ifr < nfrac; ifr++)
    {
        int face = fracset.MHMDomain(fracset.fFractureVec[ifr]);
        if (face >= 0)
        {
            out << "Line{" << ifr << "} In Surface{" << face << "};\n";
        }
    }
}

void ExportPhysicalGroups(TPZFracSet &fracset, std::ofstream &out)
{
    std::set<int> matids;
    int64_t nfrac = 0;
    int64_t nel = fracset.fgmesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fracset.fgmesh.Element(el);
        if (!gel) {
            continue;
        }
        nfrac++;
        matids.insert(gel->MaterialId());
    }
    if (nfrac != fracset.fFractureVec.NElements()) {
        DebugStop();
    }
    out << "Physical Surface(\"" << fracset.fPhysicalname << "\") = {0};\n";
    for (int64_t in=1; in < fracset.fNumMHMDomains*fracset.fNumMHMDomains; in++) {
        out << "Physical Surface(\"" << fracset.fPhysicalname << "\") += {"<< in << "};\n";
    }
    std::map<std::string,int> matidcount;

    nel = fracset.fgmesh.NElements();
    int counter = 0;
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fracset.fgmesh.Element(el);
        if (!gel) {
            continue;
        }
        int matid = gel->MaterialId();
        std::string matname = fracset.fFractureVec[counter].fPhysicalName;
        out << "Physical Line(\"" << matname <<"\") ";
        if (matidcount[matname] != 0) {
            out << '+';
        }
        out << "= {" << counter << "};\n";
        matidcount[matname]++;
        counter++;
    }

}

void ExportGMsh(TPZFracSet &fracset, const std::string &filename)
{
    std::ofstream out(filename);
    PreambleGMsh(out);
    ExportPoints(fracset, out);
    ExportFractures(fracset, out);
    CreateGMshFaces(fracset, out);
    FracturesInFaces(fracset, out);
    ExportPhysicalGroups(fracset,out);
}
