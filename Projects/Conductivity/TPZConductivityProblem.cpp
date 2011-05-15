//
//  TPZConductivityProblem.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/14/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "TPZConductivityProblem.h"
#include "pzgengrid.h"
#include "pzvoidflux.h"
#include "pzbndcond.h"
#include "pzgeoel.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzanalysis.h"
#include "TPZInterfaceEl.h"

/// set up the finite element mesh and compute the flux
REAL TPZConductivityProblem::ComputeFlux()
{
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateCompMesh();
    TPZFStructMatrix fstrmatrix(cmesh);
    TPZAnalysis an(cmesh);
    an.SetStructuralMatrix(fstrmatrix);
    TPZStepSolver solve;
    solve.SetDirect(ELU);
    an.SetSolver(solve);
    an.Run();
    TPZVec<std::string> scalnames(1),vecnames;
    scalnames[0] = "state";
    an.DefineGraphMesh(2, scalnames, vecnames, this->fGraphicsFile );
    an.PostProcess(1);
//    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
//    TPZFMatrix rhs;
//    TPZMatrix *matrix = fstrmatrix.CreateAssemble(rhs, gui);
//    matrix->Print("Global matrix",std::cout);
//    rhs.Transpose();
//    rhs.Print("Rhs",std::cout);
    return MeshFlux(cmesh);
}

/// Method to compare the current object with a copy
bool TPZConductivityProblem::Compare (TPZSaveable *copy, bool override)
{
    TPZConductivityProblem *cp = dynamic_cast<TPZConductivityProblem *>(copy);
    if(!cp) return false;
    bool result = true;
    if(fabs(fDelx[0] - cp->fDelx[0]) > 1.e-6 || fabs(fDelx[1] - cp->fDelx[1] > 1.e-6)) result = false; 
    if(fNx[0] != cp->fNx[0] || fNx[1] != fNx[1]) result = false;
    if(fabs((fConductivity-cp->fConductivity)/fConductivity) > 1.e-10) result = false;
    if(fabs((fBridgeVoidRatio-cp->fBridgeVoidRatio)/fBridgeVoidRatio) > 1.e-10) result = false;
    if(fabs((fDelPressure-cp->fDelPressure)/fDelPressure) > 1.e-10) result = false;
    if(fabs((fFluidFlux-cp->fFluidFlux)/fFluidFlux) > 1.e-8) result = false;
    if(fGraphicsFile != cp->fGraphicsFile) result = false;
    if (result == false && override) {
        *this = *cp;
    }
    return result;
}

/// Method to compare the current object with a copy
bool TPZConductivityProblem::Compare (TPZSaveable *copy, bool override) const
{
    TPZConductivityProblem *cp = dynamic_cast<TPZConductivityProblem *>(copy);
    if(!cp) return false;
    bool result = true;
    if(fabs(fDelx[0] - cp->fDelx[0]) > 1.e-6 || fabs(fDelx[1] - cp->fDelx[1] > 1.e-6)) result = false; 
    if(fNx[0] != cp->fNx[0] || fNx[1] != fNx[1]) result = false;
    if(fabs((fConductivity-cp->fConductivity)/fConductivity) > 1.e-10) result = false;
    if(fabs((fBridgeVoidRatio-cp->fBridgeVoidRatio)/fBridgeVoidRatio) > 1.e-10) result = false;
    if(fabs((fDelPressure-cp->fDelPressure)/fDelPressure) > 1.e-10) result = false;
    if(fabs((fFluidFlux-cp->fFluidFlux)/fFluidFlux) > 1.e-8) result = false;
    if(fGraphicsFile != cp->fGraphicsFile) result = false;
    if (result == false && override) {
        DebugStop();
    }
    return result;
}

int TPZConductivityProblem::ClassId () const
{
    return TPZCONDUCTIVITYID;
}

template class TPZRestoreClass<TPZConductivityProblem, TPZCONDUCTIVITYID>;

/// write this object to the TPZStream buffer. Include the classid if withclassid = true
void TPZConductivityProblem::Write(TPZStream &buf, int withclassid)
{
    buf.Write(fDelx);
    buf.Write(fNx);
    buf.Write(&fBridgeVoidRatio,1);
    buf.Write(&fConductivity,1);
    buf.Write(&fDelPressure,1);
    buf.Write(&fGraphicsFile);
    buf.Write(&fFluidFlux,1);
    
}
/*
fDelx(cp.fDelx), fNx(cp.fNx),
fBridgeVoidRatio(cp.fBridgeVoidRatio), fConductivity(cp.fConductivity), 
fDelPressure(cp.fDelPressure), fGraphicsFile(cp.fGraphicsFile)
 */
/// read objects from the stream
void TPZConductivityProblem::Read(TPZStream &buf, void *context)
{
    buf.Read(fDelx);
    buf.Read(fNx);
    buf.Read(&fBridgeVoidRatio,1);
    buf.Read(&fConductivity,1);
    buf.Read(&fDelPressure,1);
    buf.Read(&fGraphicsFile);
    buf.Read(&fFluidFlux,1);
    
}

TPZAutoPointer<TPZCompMesh> TPZConductivityProblem::GenerateCompMesh()
{
    REAL solidvoidratio = fBridgeVoidRatio;
    REAL conductivity = fConductivity;
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    TPZVec<REAL> x0(2,0.);
    TPZGenGrid gengrid(fNx,x0,fDelx);
    gengrid.Read(gmesh);
    gengrid.SetBC(gmesh, 3, -1);
    gengrid.SetBC(gmesh, 1, -2);
    REAL perimeter,area;
    perimeter = Perimeter(gmesh);
    area = DomainArea(gmesh);
//    BOOST_CHECK_SMALL(area-1., 1.e-6);
//    BOOST_CHECK_SMALL(perimeter-20., 1.e-6);
    REAL bridgesize = area*solidvoidratio/perimeter; 
//    std::cout << "Bridge size " << bridgesize << std::endl;
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    TPZAutoPointer<TPZMaterial> vflux = new TPZVoidFlux(1,conductivity,bridgesize);
    cmesh->InsertMaterialObject(vflux);
    TPZFMatrix val1(1,1,0.),val2(1,1,fDelPressure);
    TPZAutoPointer<TPZMaterial> bc1 = vflux->CreateBC(vflux, -1, 0, val1, val2);
    cmesh->InsertMaterialObject(bc1);
    val2(0,0) = 0.;
    TPZAutoPointer<TPZMaterial> bc2 = vflux->CreateBC(vflux, -2, 0, val1, val2);
    cmesh->InsertMaterialObject(bc2);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(0);
    cmesh->AutoBuild();
    return cmesh;
}

/// compute the perimeter of all two dimensional elements
REAL TPZConductivityProblem::Perimeter(TPZGeoMesh &gmesh)
{
    REAL result = 0.;
    int el, nel;
    nel = gmesh.NElements();
    for (el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel || gel->Dimension() != 2) {
            continue;
        }
        int ns = gel->NSides();
        int is;
        for (is=0; is<ns; is++) {
            if(gel->SideDimension(is) == 1)
            {
                result += gel->SideArea(is);
            }
        }
    }
    return result;
}

/// compute the perimeter of all two dimensional elements
REAL TPZConductivityProblem::DomainArea(TPZGeoMesh &gmesh)
{
    REAL result = 0.;
    int el, nel;
    nel = gmesh.NElements();
    for (el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel || gel->Dimension() != 2) {
            continue;
        }
        result += gel->Volume();
    }
    return result;    
}

/// compute the flux over the right side
REAL TPZConductivityProblem::MeshFlux(TPZAutoPointer<TPZCompMesh> cmesh)
{
    TPZFStructMatrix fstr(cmesh);
    std::set<int> matids;
    matids.insert(1);
    fstr.SetMaterialIds(matids);
    TPZFMatrix rhs;
    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
    TPZMatrix *stiff = fstr.CreateAssemble(rhs,gui);
    stiff->Multiply(cmesh->Solution(),rhs);
//    rhs.Print("residual", std::cout);
    std::set<int> lefteq,righteq;
    int iel, nel = cmesh->NElements();
    for(iel = 0; iel<nel; iel++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        TPZInterfaceElement *iface = dynamic_cast<TPZInterfaceElement *> (cel);
        if(!iface) continue;
        TPZGeoEl *gel = iface->Reference();
        if(!gel)
        {
            DebugStop();
        }
        int materialid = gel->MaterialId();
        if(materialid < 0)
        {
            TPZConnect &c = iface->Connect(0);
            int seqnum = c.SequenceNumber();
            int eq = cmesh->Block().Position(seqnum);
            if(materialid == -1)
            {
                lefteq.insert(eq);
            }
            if(materialid == -2)
            {
                righteq.insert(eq);
            }
        }
    }
    REAL leftflux(0.),rightflux(0.);
    std::set<int>::iterator it;
    for(it=lefteq.begin(); it!=lefteq.end(); it++) 
    {
        leftflux += rhs(*it,0);
    }
    for(it=righteq.begin(); it!=righteq.end(); it++) 
    {
        rightflux += rhs(*it,0);
    }
    std::cout << "leftflux = " << leftflux << " rightflux = " << rightflux << std::endl;
    return rightflux;
}


