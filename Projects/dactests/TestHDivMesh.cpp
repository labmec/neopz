//
//  TestHDivMesh.cpp
//  PZ
//
//  Created by Philippe Devloo on 10/19/14.
//
//

#include "TestHDivMesh.h"
#include "pzmultiphysicselement.h"
#include "pzintel.h"
#include "pzlog.h"
#include "pzcmesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.dactest.testhdiv"));
#endif


int CompareShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB)
{
    TPZGeoElSide gelsideA = celsideA.Reference();
    TPZGeoElSide gelsideB = celsideB.Reference();
    int sideA = gelsideA.Side();
    int sideB = gelsideB.Side();
    TPZCompEl *celA = celsideA.Element();
    TPZCompEl *celB = celsideB.Element();    TPZMultiphysicsElement *MFcelA = dynamic_cast<TPZMultiphysicsElement *>(celA);
    TPZMultiphysicsElement *MFcelB = dynamic_cast<TPZMultiphysicsElement *>(celB);
    TPZInterpolatedElement *interA = dynamic_cast<TPZInterpolatedElement *>(MFcelA->Element(0));
    TPZInterpolatedElement *interB = dynamic_cast<TPZInterpolatedElement *>(MFcelB->Element(0));
    
    TPZMaterialData dataA;
    TPZMaterialData dataB;
    interA->InitMaterialData(dataA);
    interB->InitMaterialData(dataB);
    TPZTransform<> tr = gelsideA.NeighbourSideTransform(gelsideB);
    TPZGeoEl *gelA = gelsideA.Element();
    TPZTransform<> trA = gelA->SideToSideTransform(gelsideA.Side(), gelA->NSides()-1);
    TPZGeoEl *gelB = gelsideB.Element();
    TPZTransform<> trB = gelB->SideToSideTransform(gelsideB.Side(), gelB->NSides()-1);
    
    int dimensionA = gelA->Dimension();
    int dimensionB = gelB->Dimension();
    
    int nSideshapeA = interA->NSideShapeF(sideA);
    int nSideshapeB = interB->NSideShapeF(sideB);
    int is;
    int firstShapeA = 0;
    int firstShapeB = 0;
    for (is=0; is<sideA; is++) {
        firstShapeA += interA->NSideShapeF(is);
    }
    for (is=0; is<sideB; is++) {
        firstShapeB += interB->NSideShapeF(is);
    }
    
    TPZIntPoints *intrule = gelA->CreateSideIntegrationRule(gelsideA.Side(), 4);
    int nwrong = 0;
    int npoints = intrule->NPoints();
    int ip;
    for (ip=0; ip<npoints; ip++) {
        TPZManVector<REAL,3> pointA(gelsideA.Dimension()),pointB(gelsideB.Dimension()), pointElA(gelA->Dimension()),pointElB(gelB->Dimension());
        REAL weight;
        intrule->Point(ip, pointA, weight);
        int sidedim = gelsideA.Dimension();
        TPZFNMatrix<9> jacobian(sidedim,sidedim),jacinv(sidedim,sidedim),axes(sidedim,3);
        REAL detjac;
        gelsideA.Jacobian(pointA, jacobian, jacinv, detjac, jacinv);
        TPZManVector<REAL,3> normal(3,0.), xA(3),xB(3);
        normal[0] = axes(0,1);
        normal[1] = -axes(0,0);
        tr.Apply(pointA, pointB);
        trA.Apply(pointA, pointElA);
        trB.Apply(pointB, pointElB);
        gelsideA.Element()->X(pointElA, xA);
        gelsideB.Element()->X(pointElB, xB);
        for (int i=0; i<3; i++) {
            if(fabs(xA[i]- xB[i])> 1.e-6) DebugStop();
        }
        int nshapeA = 0, nshapeB = 0;
        interA->ComputeRequiredData(dataA, pointElA);
        interB->ComputeRequiredData(dataB, pointElB);
        nshapeA = dataA.phi.Rows();
        nshapeB = dataB.phi.Rows();
        if(nSideshapeA != nSideshapeB) DebugStop();
        
        TPZManVector<REAL> shapesA(nSideshapeA), shapesB(nSideshapeB);
        int nwrongkeep(nwrong);
        int i,j;
        for(i=firstShapeA,j=firstShapeB; i<firstShapeA+nSideshapeA; i++,j++)
        {
            int Ashapeind = i;
            int Bshapeind = j;
            int Avecind = -1;
            int Bvecind = -1;
            // if A or B are boundary elements, their shapefunctions come in the right order
            if (dimensionA != sidedim) {
                Ashapeind = dataA.fVecShapeIndex[i].second;
                Avecind = dataA.fVecShapeIndex[i].first;
            }
            if (dimensionB != sidedim) {
                Bshapeind = dataB.fVecShapeIndex[j].second;
                Bvecind = dataB.fVecShapeIndex[j].first;
            }
            if (dimensionA != sidedim && dimensionB != sidedim) {
                // vefify that the normal component of the normal vector corresponds
                Avecind = dataA.fVecShapeIndex[i].first;
                Bvecind = dataB.fVecShapeIndex[j].first;
                REAL vecnormalA = dataA.fNormalVec(0,Avecind)*normal[0]+dataA.fNormalVec(1,Avecind)*normal[1];
                REAL vecnormalB = dataB.fNormalVec(0,Bvecind)*normal[0]+dataB.fNormalVec(1,Bvecind)*normal[1];
                if(fabs(vecnormalA-vecnormalB) > 1.e-6)
                {
                    nwrong++;
                    LOGPZ_ERROR(logger, "normal vectors aren't equal")
                }
                
            }
            shapesA[i-firstShapeA] = dataA.phi(Ashapeind,0);
            shapesB[j-firstShapeB] = dataB.phi(Bshapeind,0);
            REAL valA = dataA.phi(Ashapeind,0);
            REAL valB = dataB.phi(Bshapeind,0);
            REAL diff = valA-valB;
            REAL decision = fabs(diff)-1.e-6;
            if(decision > 0.)
            {
                nwrong ++;
                std::cout << "valA = " << valA << " valB = " << valB << " Avecind " << Avecind << " Bvecind " << Bvecind <<
                " Ashapeind " << Ashapeind << " Bshapeind " << Bshapeind <<
                " sideA " << sideA << " sideB " << sideB << std::endl;
                LOGPZ_ERROR(logger, "shape function values are different")
            }
        }
        if(nwrong != nwrongkeep)
        {
            std::cout << "shapeA " << shapesA << std::endl;
            std::cout << "shapeB " << shapesB << std::endl;
        }
        //        if(nwrong)
        //        {
        //            std::cout << "\nNumber of different shape functions " << nwrong << std::endl;
        //        }
        //        BOOST_CHECK(nwrong == 0);
    }
    delete intrule;
    return nwrong;
}

int CompareSideShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB)
{
    TPZGeoElSide gelsideA = celsideA.Reference();
    TPZGeoElSide gelsideB = celsideB.Reference();
    int sideA = gelsideA.Side();
    int sideB = gelsideB.Side();
    const int sidedimension = gelsideA.Dimension();
    TPZCompEl *celA = celsideA.Element();
    TPZCompEl *celB = celsideB.Element();
    TPZMultiphysicsElement *MFcelA = dynamic_cast<TPZMultiphysicsElement *>(celA);
    TPZMultiphysicsElement *MFcelB = dynamic_cast<TPZMultiphysicsElement *>(celB);
    TPZInterpolatedElement *interA = dynamic_cast<TPZInterpolatedElement *>(MFcelA->Element(0));
    TPZInterpolatedElement *interB = dynamic_cast<TPZInterpolatedElement *>(MFcelB->Element(0));
    if (!interA || !interB) {
        DebugStop();
    }
    TPZGeoEl *gel = gelsideA.Element();
    int nshapeA = interA->NSideShapeF(sideA);
    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gelsideA.Side(), 4);
    int nwrong = 0;
    int npoints = intrule->NPoints();
    int ip;
    for (ip=0; ip<npoints; ip++) {
        TPZManVector<REAL,3> pointA(gelsideA.Dimension()),pointB(gelsideB.Dimension());
        REAL weight;
        intrule->Point(ip, pointA, weight);
        TPZTransform<> tr = gelsideA.NeighbourSideTransform(gelsideB);
        tr.Apply(pointA, pointB);
        TPZFNMatrix<200> phiA(nshapeA,1),dphiA(sidedimension,nshapeA),phiB(nshapeA,1),dphiB(sidedimension,nshapeA);
        interA->SideShapeFunction(sideA, pointA, phiA, dphiA);
        interB->SideShapeFunction(sideB, pointB, phiB, dphiB);
        int nshapeA = phiA.Rows();
        int nshapeB = phiB.Rows();
        if(nshapeA != nshapeB) DebugStop();
        int ish;
        for (ish=0; ish<nshapeA; ish++) {
            REAL Aval = phiA(ish,0);
            REAL Bval = phiB(ish,0);
            if(fabs(Aval-Bval) > 1.e-6)
            {
                std::cout << "i " << ish << " " << Aval << " " << Bval << std::endl;
                nwrong++;
            }
        }
        if(nwrong)
        {
            std::cout << "\nNumber of different shape functions " << nwrong << std::endl;
        }
        //        BOOST_CHECK(nwrong == 0);
    }
    return nwrong;
}

#include "TPZMultiphysicsInterfaceEl.h"

void TestMesh(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZMultiphysicsInterfaceElement *mfint = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (mfint) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        int nsides = gel->NSides();
        for (int side=0; side<nsides; side++) {
            if (gel->SideDimension(side) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                TPZCompEl *neighcel = neighbour.Element()->Reference();
                TPZMultiphysicsInterfaceElement *mfint = dynamic_cast<TPZMultiphysicsInterfaceElement *>(neighcel);
                if (!mfint) {
                    TPZCompElSide compelsideA = gelside.Reference();
                    TPZCompElSide compelsideB = neighbour.Reference();
                    CompareShapeFunctions(compelsideA, compelsideB);
                    CompareSideShapeFunctions(compelsideA, compelsideB);
                }
                
                neighbour = neighbour.Neighbour();
            }
        }
    }
}
