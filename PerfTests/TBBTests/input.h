/*******************************************************************************
 *   Copyright (C) 2014 by:                                                    *
 *   Gilvan Vieira (gilvandsv@gmail.com)                                       *
 *                                                                             *
 *   This program is free software; you can redistribute it and/or modify      *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation; either version 2 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program; if not, write to the                             *
 *   Free Software Foundation, Inc.,                                           *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                 *
 ******************************************************************************/

#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrprecond.h"
#include "pzdohrstructmatrix.h"
#include "pzstepsolver.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "tpzdohrassembly.h"
#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "tpzpairstructmatrix.h"
#include "pzviscoelastic.h"
#include "TPZTimer.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeotetrahedra.h"
#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "tpzarc3d.h"
#include "tpzdohrmatrix.h"
#include "pzvtkmesh.h"
#include "pzlog.h"

#include <string>

namespace Input {
    /// Generate a boundary geometric element at the indicated node
    void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
    {
        // look for an element/corner node whose distance is close to start
        TPZGeoNode *gn1 = gr->FindNode(x);
        int64_t iel;
        int64_t nelem = gr->ElementVec().NElements();
        TPZGeoEl *gel;
        for (iel = 0; iel<nelem; iel++) {
            gel = gr->ElementVec()[iel];
            if(!gel) continue;
            int nc = gel->NCornerNodes();
            int c;
            for (c=0; c<nc; c++) {
                TPZGeoNode *gn = gel->NodePtr(c);
                if (gn == gn1) {
                    break;
                }
            }
            if (c<nc) {
                TPZGeoElBC(gel, c, bc);
                return;
            }
        }
    }
    
    void InsertElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh)
    {
        mesh->SetDimModel(3);
        int nummat = 1, neumann = 1, mixed = 2;
        //      int dirichlet = 0;
        int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5;   //, dirp2 = -6;
        TPZManVector<STATE> force(3,0.);
        //force[1] = 0.;
        
        STATE ElaE = 1000., poissonE = 0.2;   //, poissonV = 0.1, ElaV = 100.;
        
        STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
        lambdaV = 11.3636;
        muV = 45.4545;
        alpha = 1.;
        deltaT = 0.01;
        
        //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
        //viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
        //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, ElaE, poissonE, lambdaV, muV, alphaT, force);
        TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, ElaE, poissonE, force);
        
        TPZFNMatrix<6> qsi(6,1,0.);
        //viscoelast->SetDefaultMem(qsi); //elast
        //int index = viscoelast->PushMemItem(); //elast
        TPZMaterial * viscoelastauto(viscoelast);
        mesh->InsertMaterialObject(viscoelastauto);
        
        // Neumann em x = 1;
        TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
        val2(0,0) = 1.;
        TPZBndCond *bc4 = viscoelast->CreateBC(viscoelastauto, neumann1, neumann, val1, val2);
        TPZMaterial * bcauto4(bc4);
        mesh->InsertMaterialObject(bcauto4);
        
        // Neumann em x = -1;
        val2(0,0) = -1.;
        TPZBndCond *bc5 = viscoelast->CreateBC(viscoelastauto, neumann2, neumann, val1, val2);
        TPZMaterial * bcauto5(bc5);
        mesh->InsertMaterialObject(bcauto5);
        
        val2.Zero();
        // Dirichlet em -1 -1 -1 xyz;
        val1(0,0) = 1e4;
        val1(1,1) = 1e4;
        val1(2,2) = 1e4;
        TPZBndCond *bc1 = viscoelast->CreateBC(viscoelastauto, dir1, mixed, val1, val2);
        TPZMaterial * bcauto1(bc1);
        mesh->InsertMaterialObject(bcauto1);
        
        // Dirichlet em 1 -1 -1 yz;
        val1(0,0) = 0.;
        val1(1,1) = 1e4;
        val1(2,2) = 1e4;
        TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
        TPZMaterial * bcauto2(bc2);
        mesh->InsertMaterialObject(bcauto2);
        
        // Dirichlet em 1 1 -1 z;
        val1(0,0) = 0.;
        val1(1,1) = 0.;
        val1(2,2) = 1e4;
        TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
        TPZMaterial * bcauto3(bc3);
        mesh->InsertMaterialObject(bcauto3);
        
    }
    
    TPZGeoMesh *MalhaCubo(string FileName)
    {
        int numnodes=-1;
        int numelements=-1;
        
        {
            bool countnodes = false;
            bool countelements = false;
            
            ifstream read (FileName.c_str());
            
            while(read)
            {
                char buf[1024];
                read.getline(buf, 1024);
                std::string str(buf);
                if(str == "Coordinates") countnodes = true;
                if(str == "end coordinates") countnodes = false;
                if(countnodes) numnodes++;
                
                if(str == "Elements") countelements = true;
                if(str == "end elements") countelements = false;
                if(countelements) numelements++;
            }
        }
        
        TPZGeoMesh * gMesh = new TPZGeoMesh;
        
        gMesh -> NodeVec().Resize(numnodes);
        
        TPZManVector <int64_t> TopolTetra(4);
        
        const int Qnodes = numnodes;
        TPZVec <TPZGeoNode> Node(Qnodes);
        
        //setting nodes coords
        int64_t nodeId = 0, elementId = 0, matElId = 1;
        
        ifstream read;
        read.open(FileName.c_str());
        
        double nodecoordX , nodecoordY , nodecoordZ ;
        
        char buf[1024];
        read.getline(buf, 1024);
        read.getline(buf, 1024);
        std::string str(buf);
        int in;
        for(in=0; in<numnodes; in++)
        {
            read >> nodeId;
            read >> nodecoordX;
            read >> nodecoordY;
            read >> nodecoordZ;
            Node[nodeId-1].SetNodeId(nodeId);
            Node[nodeId-1].SetCoord(0,nodecoordX);
            Node[nodeId-1].SetCoord(1,nodecoordY);
            Node[nodeId-1].SetCoord(2,nodecoordZ);
            gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
        }
        
        {
            read.close();
            read.open(FileName.c_str());
            
            int l , m = numnodes+5;
            for(l=0; l<m; l++)
            {
                read.getline(buf, 1024);
            }
            
            
            int el;
            int neumann1 = -4, neumann2 = -5;
            //std::set<int> ncoordz; //jeitoCaju
            for(el=0; el<numelements; el++)
            {
                read >> elementId;
                read >> TopolTetra[0]; //node 1
                read >> TopolTetra[1]; //node 2
                read >> TopolTetra[2]; //node 3
                read >> TopolTetra[3]; //node 4
                
                // O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
                TopolTetra[0]--;
                TopolTetra[1]--;
                TopolTetra[2]--;
                TopolTetra[3]--;
                
                int64_t index = el;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
            }
            
            gMesh->BuildConnectivity();
            
            // Colocando as condicoes de contorno
            for(el=0; el<numelements; el++)
            {
                TPZManVector <TPZGeoNode,4> Nodefinder(4);
                TPZManVector <REAL,3> nodecoord(3);
                TPZGeoEl *tetra = gMesh->ElementVec()[el];
                
                // na face x = 1
                TPZVec<int64_t> ncoordzVec(0); int64_t sizeOfVec = 0;
                for (int i = 0; i < 4; i++)
                {
                    int pos = tetra->NodeIndex(i);
                    Nodefinder[i] = gMesh->NodeVec()[pos];
                    Nodefinder[i].GetCoordinates(nodecoord);
                    if (nodecoord[0] == 1.)
                    {
                        sizeOfVec++;
                        ncoordzVec.Resize(sizeOfVec);
                        ncoordzVec[sizeOfVec-1] = pos;
                    }
                }
                if(ncoordzVec.NElements() == 3)
                {
                    int lado = tetra->WhichSide(ncoordzVec);
                    TPZGeoElSide tetraSide(tetra, lado);
                    TPZGeoElBC(tetraSide,neumann1);
                }
                
                // Na face x = -1
                ncoordzVec.Resize(0);
                sizeOfVec = 0;
                for (int i = 0; i < 4; i++)
                {
                    int pos = tetra->NodeIndex(i);
                    Nodefinder[i] = gMesh->NodeVec()[pos];
                    
                    Nodefinder[i].GetCoordinates(nodecoord);
                    if (nodecoord[0] == -1.)
                    {
                        sizeOfVec++;
                        ncoordzVec.Resize(sizeOfVec);
                        ncoordzVec[sizeOfVec-1] = pos;
                    }
                }
                if(ncoordzVec.NElements() == 3)
                {
                    int lado = tetra->WhichSide(ncoordzVec);
                    TPZGeoElSide tetraSide(tetra, lado);
                    TPZGeoElBC(tetraSide,neumann2);
                }
                
            }
            
            TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
            yz[0] = 1.;
            z[2] = -1;
            int bcidxyz = -1, bcidyz = -2, bcidz = -3;
            SetPointBC(gMesh, xyz, bcidxyz);
            SetPointBC(gMesh, yz, bcidyz);
            SetPointBC(gMesh, z, bcidz);
            
        }
        
        return gMesh;
    }
    
    TPZAutoPointer<TPZMatrix<REAL> >  CreateCuboSkyMatrix(string filename, int plevel)
    {
        TPZGeoMesh *gmesh = 0;
        
        TPZAutoPointer<TPZCompMesh> cmesh;
        
        gmesh = MalhaCubo(filename);
        cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDimModel(3);
        InsertElasticityCubo(cmesh);
        cmesh->SetDefaultOrder(plevel);
        cmesh->AutoBuild();
        
        // Gerando a matriz
        TPZLinearAnalysis an(cmesh);
        TPZSkylineStructMatrix skyl(cmesh);
        an.SetStructuralMatrix(skyl);
        TPZStepSolver<REAL> step;
        step.SetDirect(ECholesky);
        an.SetSolver(step);
        an.Assemble();
        TPZAutoPointer<TPZMatrix<REAL> > skylmat1 = new TPZSkylMatrix<REAL>();
        skylmat1 = an.Solver().Matrix();
        
        return skylmat1;
        //	TPZSkylMatrix<REAL> *orig = dynamic_cast<TPZSkylMatrix<REAL> *> (skylmat1.operator->());
        //
        //    return orig;
    }
    
}


