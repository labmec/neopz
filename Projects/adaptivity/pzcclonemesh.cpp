//METHODS DEFINITIONS FOR CLASS COMPUTATIONAL MESH

#include "pzerror.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzcclonemesh.h"
#include "pzgclonemesh.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzconnect.h"
#include "pzbndcond.h"
#include "pzmaterial.h"
#include "pzsolve.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzelmat.h"
#include "pztrnsform.h"
#include "pztransfer.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzmetis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzquad.h"
#include "pzonedref.h"
#include "pzvec.h"
#include "pzskylstrmatrix.h"
#include "pzdebug.h"
#include "pzcheckmesh.h"
#include "pzanalysis.h"

using namespace std;

template class TPZVec<TPZCompCloneMesh::TPZRefPattern>;

static ofstream gDeduce("deduce.txt");
static int zero = 0;

TPZCompCloneMesh::TPZCompCloneMesh (TPZGeoCloneMesh* gr, TPZCompMesh *cmesh) : TPZCompMesh(gr), fMapConnects(){
    fCloneReference = cmesh;
    //Cria um clone do vetor de materiais da malha mesh
    std::map<int, TPZAutoPointer<TPZMaterial> >::const_iterator it;
    for(it=cmesh->MaterialVec().begin(); it != cmesh->MaterialVec().end() ; it++) {
        //    mat->Print();
        it->second->Clone(MaterialVec());
    }
    //	Reference()->Print(cout);
}

TPZCompCloneMesh::~TPZCompCloneMesh() {
}

void TPZCompCloneMesh::AutoBuild() {
    TPZAdmChunkVector<TPZGeoEl *> &elvec = Reference()->ElementVec();
    int i,j, nelem = elvec.NElements();
    int index;
    TPZGeoCloneMesh *gclm =  dynamic_cast<TPZGeoCloneMesh *>(Reference());
    if (!gclm) {
        cout << "TPZCompCloneMesh::AutoBuild : clone mesh not initialised" <<endl;
    }
    
    if (gDebug) {
        gclm->SetName("Malha Clone Geometrica");
        gclm->Print(cout);
        //Reference()->Print(cout);
    }
    
    
    for(i=0; i<nelem; i++) {
        TPZGeoEl *gel = elvec[i];
        TPZGeoEl *clgel = gclm->ReferenceElement(i);
        if (clgel) {
            TPZCompEl *cel = clgel->Reference();
            if (cel){
                if(!gel){
                    cout << "TPZCompCloneMesh::AutoBuild: null geometric element detected" << endl;
                    continue;
                }
                if(gclm->IsPatchSon(gel)) {
                    
                    if (gDebug){
                        cout << "TPZCompCloneMesh::AutoBuild : Creating computational element \n Geometric Reference Element:\n"
                        << endl;
                        gel->Print();
                    }
                    
                    TPZCompEl *clcel = this->CreateCompEl(gel,index);
                    
                    if (gDebug){
                        cout << "TPZCompCloneMesh::AutoBuild : Computational element created:\n" << endl;
                        clcel->Print();
                    }
                    
                    int ncon = clcel->NConnects();
                    for (j=0; j<ncon; j++){
                        int refcon = cel->ConnectIndex(j);
                        int conid = clcel->ConnectIndex(j);
                        
                        if (gDebug){
                            cout << "Connects --- Reference :  " << refcon << "   Clone   " << conid << endl;
                        }
                        
                        if (! HasConnect(refcon) ){
                            fMapConnects [refcon] = conid;
                            if(conid == fOriginalConnects.NElements()) {
                                fOriginalConnects.Push(refcon);
                            } else if(conid < fOriginalConnects.NElements()) {
                                fOriginalConnects[conid] = refcon;
                            } else {
                                fOriginalConnects.Resize(conid+1);
                                fOriginalConnects[conid] = refcon;
                            }
                        }
                    }
                    TPZInterpolatedElement *orgintel = dynamic_cast<TPZInterpolatedElement *> (cel);
                    TPZInterpolatedElement *clintel = dynamic_cast<TPZInterpolatedElement *> (clcel);
                    for (j=0;j<cel->Reference()->NSides();j++){
                        
                        if (gDebug){
                            cout << "TPZCompCloneMesh::AutoBuild :Computational Element Before  PRefine:\n " << endl;
                            clcel->Print();
                        }
                        
                        int porder = orgintel->SideOrder(j);
                        // Check if everything is fine
                        //DebugStop();
                        TPZConnect &c = clintel->Connect(j);
                        c.SetOrder(porder);
                        //	    clintel->PRefine(j,porder);
                        
                        if (gDebug){
                            cout << "Porder " << porder << "    Side : " << j << endl;
                            cout << "TPZCompCloneMesh::AutoBuild :Computational Element After  PRefine:\n " << endl;
                            clcel->Print();
                        }
                    }
                }
            }
        }
    }
    //CleanUp
    ExpandSolution();
    CleanUpUnconnectedNodes();
    if (gDebug){
        cout << "TPZCompCloneMesh::AutoBuild :Computational Mesh Before BC\n " << endl;
        Print(cout);
    }
    
    CreateCloneBC();
    
/*    
    TPZCompEl *cel;
    TPZAdmChunkVector<TPZGeoElBC> &elbcvec = Reference()->BCElementVec();
    nelem = elbcvec.NElements();
    
    if (gDebug){
        cout << "TPZCompCloneMesh::AutoBuild : After Clone BCs generation\n " << endl;
    }
    
    for(i=0; i<nelem; i++) {
        if(!elbcvec[i].fBCElement) {
            cel = elbcvec[i].fElement->CreateBCCompEl(elbcvec[i].fSide,elbcvec[i].fId,*this);
            if(cel){
                elbcvec[i].fBCElement = cel->Reference();
                
                if (gDebug){
                    cout << "TPZCompCloneMesh::AutoBuild : BC Clone\n " << endl;
                    cel->Print(cout);
                }
                
            }
        }
    }
*/
    InitializeBlock();
    
    //	Print(cout);
    
    //Copiar Solução Bloco a Bloco
    int nc = fCloneReference->NConnects();
    for (i=0;i<nc;i++){
        if(! HasConnect(i)) continue;
        int clseqnum 	= ConnectVec()[fMapConnects[i]].SequenceNumber();
        int orgseqnum 	= fCloneReference->ConnectVec()[i].SequenceNumber();
        int ndoforg 	= fCloneReference->ConnectVec()[i].NDof(*fCloneReference);
        int ndofclone 	= ConnectVec()[fMapConnects[i]].NDof(*this);
        if( ndoforg != ndofclone) {
            //      Print(cout);
            cout << "Number of degree of freedom incompatible between clone and original mesh!\n";
            cout << "Clone connect id: " << i << "  Clone dof: " << ndofclone << "  Original dof: " << ndoforg << endl;
            DebugStop();
            continue;
        }
        for (j=0;j<ndoforg;j++){//clpos; j<(clpos + ndoforg); j++){
            //	cout << "Clone connect id: " << i << "  Clone dof: " << ndofclone << "  Original dof: " << ndoforg << "  Value:  " <<fCloneReference->Block()(orgseqnum,0,j,0) << endl;
            fSolutionBlock (clseqnum,0,j,0) = fCloneReference->Block()(orgseqnum,0,j,0);
        }
    }
    int printing=0;
    if(printing) {
        ofstream out("clonesolution.txt");
        fSolution.Print("Solution of the clone",out);
    }
}

/*
 void TPZCompCloneMesh::CreateCloneBC(){
 int i,j;//elementos e lados de elementos
 int nstate = MaterialVec()[0]->NStateVariables();
 TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
 TPZMaterial *bnd = MaterialVec()[0]->CreateBC (-1000,50,val1,val2);
 InsertMaterialObject(bnd);
 
 int ncon = ConnectVec().NElements();
 TPZVec<int> flagConn (ncon,0);
 
 int tmporder = TPZCompEl::gOrder;
 TPZCompEl::gOrder = 10;
 
 int printing = 0;
 if(printing) {
 ofstream test("test.txt",ios::app);
 Print(test);
 }
 
 for (i=0;i<NElements();i++){
 TPZInterpolatedElement *el = dynamic_cast<TPZInterpolatedElement *> ( ElementVec()[i]);
 int nsid = el->Reference()->NSides();
 for (j=nsid-1;j>=0;j--){
 TPZCompElSide side (el,j);
 TPZGeoElSide geoside = side.Reference();
 int clconid = side.ConnectIndex();
 int orgconid = fOriginalConnects[clconid];
 if (flagConn[clconid] == 1) continue;
 int orgelcon = fCloneReference->ConnectVec()[orgconid].NElConnected();
 int clelcon = ConnectVec()[clconid].NElConnected();
 TPZInterpolatedElement *origel = GetOriginalElement(el);
 int constraints =  fCloneReference->ConnectVec()[orgconid].HasDependency();
 int orgconstraints = ConnectVec()[clconid].HasDependency();
 if (orgelcon > clelcon || constraints != orgconstraints){
 TPZCompEl *celbc = el->Reference()->CreateBCCompEl(j,-1000,*this);
 TPZInterpolatedElement *celbcint = dynamic_cast<TPZInterpolatedElement *>(celbc);
 if(celbc){
 //elbcvec[i].fBCElement = cel->Reference();
 //TPZGeoElBC bc (geoside,-1000,*Reference());
 int ncreatedcon = celbc->NConnects();
 //celbcint->PRefine(ncreatedcon-1,origel->SideOrder(j));
 //el->PRefine(j,origel->SideOrder(j));
 int k;
 for (k=0;k<ncreatedcon;k++){
 flagConn[celbc->ConnectIndex(k)] = 1;
 }
 }
 }
 }
 }
 TPZCompEl::gOrder = tmporder;
 //  AdjustBoundaryElements();
 }
 */

void TPZCompCloneMesh::CreateCloneBC(){
    int i,j;//elementos e lados de elementos
    TPZAutoPointer<TPZMaterial> mat = MaterialVec().rbegin()->second;
    int nstate = mat->NStateVariables();
    int dim = mat->Dimension();
    TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
    TPZMaterial *bnd = mat->CreateBC (mat,-1000,50,val1,val2);
    InsertMaterialObject(bnd);
    
    Reference()->ResetReference();
    LoadReferences();
    fCloneReference->Reference()->ResetReference();
    this->fCloneReference->LoadReferences();
    TPZStack<TPZGeoElSide> bcelsides;
    int ncon = ConnectVec().NElements();
    TPZVec<int> flagConn (ncon,0);
    
    /*   int aux = fOriginalConnects.NElements(); */
    /*   cout << "Número de connects original: " << aux << endl; */
    /*   for (i=0;i<aux;i++){ */
    /*     cout << "i= " << fOriginalConnects[i] << endl; */
    /*   } */
    /*   cout << "----------------------------\n"; */
    
    /*  fReference->Print(cout); */
    /*  Print(cout); */
//    int tmporder = cmesh->GetDefaultOrder();
//    cmesh->SetDefaultOrder(10);
    int tmporder = GetDefaultOrder();
    SetDefaultOrder(10);
    
    int printing = 0;
    if(printing) {
        ofstream test("test.txt",ios::app);
        Print(test);
    }
    
    for (i=0;i<NElements();i++){
        TPZInterpolatedElement *el = dynamic_cast<TPZInterpolatedElement *> ( ElementVec()[i]);
        int nsid = el->Reference()->NSides();
        for (j=nsid-1;j>=0;j--){
            TPZCompElSide side (el,j);
            TPZGeoElSide geoside = side.Reference();
            if(geoside.Dimension() != dim-1) continue;
            int clconid = side.ConnectIndex();
            int orgconid = fOriginalConnects[clconid];
            if (flagConn[clconid] == 1) continue;
            int orgelcon = fCloneReference->ConnectVec()[orgconid].NElConnected();
            int clelcon = ConnectVec()[clconid].NElConnected();
            TPZInterpolatedElement *origel = GetOriginalElement(el);
            TPZCompElSide orgside(origel,j);
            TPZCompElSide cllarge = side.LowerLevelElementList(1);
            TPZCompElSide large = orgside.LowerLevelElementList(1);
            //      int constraints =  fCloneReference->ConnectVec()[orgconid].HasDependency();
            //     int orgconstraints = ConnectVec()[clconid].HasDependency();
            if (orgelcon > clelcon || cllarge.Exists() != large.Exists()){
                bcelsides.Push(geoside);
            }
        }
    }
    Reference()->ResetReference();
    this->LoadReferences();
    int nbc = bcelsides.NElements();
    int ibc;
    for(ibc = 0; ibc<nbc; ibc++) {
        TPZCompEl *celbc = bcelsides[ibc].Element()->CreateBCCompEl(bcelsides[ibc].Side(),-1000,*this);
    }
//    cmesh->SetDefaultOrder(tmporder);
    SetDefaultOrder(tmporder);
    //  AdjustBoundaryElements();
}

int TPZCompCloneMesh::HasConnect(int cnid){
    return (fMapConnects.find(cnid) != fMapConnects.end());
}

TPZCompMesh * TPZCompCloneMesh::UniformlyRefineMesh() {
    ExpandSolution();
    if(gPrintLevel ==2) {
        TPZCheckMesh chck(this,&cout);
        chck.VerifyAllConnects();
        Reference()->Print(cout);
        this->Print(cout);
    }
    TPZStack <TPZGeoEl *> bcgelstack;
    TPZStack <int> bcporderstack;
    TPZStack<int> elementindex;
    TPZStack<TPZCompEl *> elementpointers;
    int el, nelem;
    nelem = ElementVec().NElements();
    // take out the boundary elements with inconsistent father structure
    for(el=0; el<nelem; el++) {
        TPZInterpolatedElement *bcel =  dynamic_cast<TPZInterpolatedElement *>(ElementVec()[el]);
        if(bcel) {
            if (bcel->Material()->Id() == -1000) {
                bcgelstack.Push(bcel->Reference());
                int nsides = bcel->NConnects();
                int sideorder = bcel->SideOrder(nsides-1);
                bcporderstack.Push(sideorder);
                elementindex.Push(el);
                elementpointers.Push(bcel);
                ElementVec()[el] = 0;
            }
        }
    }
    TPZCompMesh *cmesh = Clone();
    // put the elements back in
    nelem = elementpointers.NElements();
    for(el=0; el<nelem; el++) {
        ElementVec()[elementindex[el]] = elementpointers[el];
    }
    
    if(gPrintLevel == 2) {
        TPZCheckMesh chk(cmesh,&cout);
        chk.VerifyAllConnects();
    }
    Reference()->ResetReference();
    cmesh->LoadReferences();
    
    Reference()->Print(std::cout);
    
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    nelem = elementvec.NElements();
    
    TPZStack<TPZCompEl *> copyel;
    for(el=0; el<nelem; el++) {
        TPZCompEl *cel = elementvec[el];
        if(!cel) continue;
        copyel.Push(cel);
    }
    
    nelem = copyel.NElements();
    for(el=0; el<nelem; el++) {
        
        TPZCompEl *cel = copyel[el];
        TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (cel);
        if(!cint) continue;
        int ncon = cint->NConnects();
        int porder = cint->PreferredSideOrder(ncon-1);
        
        if (gDebug){
            cout << "TPZCompCloneMesh::UniformlyRefineMesh Element Data Before PRefine\n";
            cout << "NConnects() :  " << ncon << "   POrder : " << porder;
        }
        
        TPZVec<int> subelindex;
        cint->Divide(el,subelindex,1);
        
        if (gDebug){
            cout << "TPZCompCloneMesh::UniformlyRefineMesh Element Data After Divide\n";
            int idbg, indbg, neldbg = subelindex.NElements();
            for (idbg = 0; idbg<neldbg; idbg++){
                cout << "\nElement : " << idbg << endl;
                indbg = subelindex[idbg];
                TPZCompEl *celdbg = elementvec[indbg];
                if (celdbg) celdbg->Print(cout);
                else cout << "SubElement not initialized!!\n";
            }
        }
        
        int isub;
        for (isub=0; isub<subelindex.NElements();isub++){
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelindex[isub]]);
            intel->PRefine(porder+1);
            
            if(gPrintLevel == 2) {
                TPZCheckMesh chk(cmesh,&cout);
                chk.VerifyAllConnects();
            }
            if (gDebug){
                cout << "TPZCompCloneMesh::UniformlyRefineMesh Element Data After PRefine\n";
                intel->Print(cout);
            }
            
        }
        cmesh->ExpandSolution();
    }
    
    int tempgorder = cmesh->GetDefaultOrder();
    
    int nbc = bcgelstack.NElements();
    int i;
    for (i=0;i<nbc;i++){
        TPZGeoEl *bcgel =  bcgelstack[i];
        
        TPZVec <TPZGeoEl *> bcsubgel;
        bcgel->Divide(bcsubgel);
        
        int nsubel = bcsubgel.NElements();
        int is;
        for (is=0;is<nsubel;is++){
            int elord = bcporderstack[i]+1;
            int maxp = TPZOneDRef::gMaxP;
            elord = elord > maxp ? maxp : elord;
            cmesh->SetDefaultOrder(elord);
            int indexsub;
            cmesh->CreateCompEl(bcsubgel[is],indexsub);
        }
    }
    cmesh->SetDefaultOrder(tempgorder);
    cmesh->ExpandSolution();
    cmesh->InitializeBlock();
    
    return cmesh;
}

void TPZCompCloneMesh::MeshError(TPZCompMesh *fine,
                                 TPZVec<REAL> &ervec,
                                 void(*f)(TPZVec<REAL> &loc,
                                          TPZVec<REAL> &val,
                                          TPZFMatrix &deriv),
                                 TPZVec<REAL> &truervec){
    //Evaluates the solution f
    
    if (fine->Reference()->Reference() != fine){
        fine->Reference()->ResetReference();
        fine->LoadReferences();
    }
    
    TPZGeoCloneMesh *gclmesh = dynamic_cast<TPZGeoCloneMesh *> (fine->Reference());
    if (gclmesh->GetMeshReferenceElement()->MaterialId() < 0) return;
    
    int diagnostic = 0;
    if(diagnostic) {
        ofstream test("test.txt",ios::app);
        Print(test);
        Solution().Print("coarse mesh solution",test);
        fine->Reference()->Print(test);
        fine->Print(test);
        fine->Solution().Print("fine mesh solution", test);
    }
    
    
    int computesolution = 1;
    //  fine->Solution().Print();
    if(computesolution) {
        TPZSkylineStructMatrix clfstr(fine);
        TPZAnalysis clfan(fine);
        TPZStepSolver cldirect;
        cldirect.SetDirect(ELDLt);
        clfan.SetStructuralMatrix(clfstr);
        clfan.SetSolver(cldirect);
        clfan.Run();
        
        int printing=0;
        if(printing) {
            clfan.Solution().Print();
        }
    }
    /*   if(computesolution) { */
    /*     TPZSkylineStructMatrix clfstr(fine); */
    /*     TPZAnalysis clfan(fine); */
    /*     TPZStepSolver jacobi; */
    /*     TPZStepSolver cg; */
    /*     clfan.SetStructuralMatrix(clfstr); */
    /*     clfan.SetSolver(cg); */
    /*     clfan.Assemble(); */
    /*     jacobi.ShareMatrix(clfan.Solver()); */
    /*     jacobi.SetJacobi(1,0,0); */
    /*     cg.ShareMatrix(clfan.Solver()); */
    /*     cg.SetCG(5,jacobi,1E-6,1); */
    /*     clfan.SetSolver(cg); */
    /*     clfan.Solve(); */
    /*     int printing=0; */
    /*     if(printing) { */
    /*       clfan.Solution().Print(); */
    /*     } */
    /*   } */
    //  clfan.Rhs().Print();
    //  cout << "Solution of the fine mesh\n";
    //   clfan.Solution().Print();
    
    Reference()->ResetReference();
    LoadReferences();
    TPZAutoPointer<TPZMaterial> mat = fine->MaterialVec().rbegin()->second;
    int dim = mat->Dimension();
    TPZCompEl *cel;
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fine->ElementVec();
    int numel = elementvec.NElements();
    int el;
    //  TPZGeoCloneMesh *gclmesh = dynamic_cast<TPZGeoCloneMesh *> (Reference());
    /*   fine->SetName("Malha Fina"); */
    /*   fine->Print(); */
    /*   fine->Reference()->Print(); */
    for(el=0; el<numel; el++) {
        cel = elementvec[el];
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!IsFather(gel)) continue;
        
        //convergencia...
        //    TPZGeoMesh *orgeomesh = fCloneReference->Reference();
        //    TPZCompMesh *aux =  orgeomesh->Reference();
        //    orgeomesh->ResetReference();
        //    fCloneReference->LoadReferences();
        //   int orgcelind = GetOriginalElementIndex(el);
        //    if (orgcelind >= 0) cel = fCloneReference->ElementVec()[orgcelind];
        //    else cel = 0;
        //    orgeomesh->ResetReference();
        //    aux->LoadReferences();
        
        
        if(!cel) continue;
        TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!cint) continue;
        int ncon = cint->NConnects();
        TPZGeoElSide gelside(cint->Reference(),ncon-1);
        if(gelside.Dimension() != dim) continue;
        TPZGeoElSide gellarge(gelside);
        while(!gellarge.Reference().Exists() && gellarge.Father2().Exists()) gellarge = gellarge.Father2();
        if(!gellarge.Reference().Exists()) {
            cout << "TPZCompCloneMesh::BuildTransferMatrix element " << el << " found no corresponding element\n";
            continue;
        }
        TPZCompElSide cellargeside = gellarge.Reference();
        TPZCompEl *cellarge = cellargeside.Element();
        TPZInterpolatedElement *cintlarge = (TPZInterpolatedElement *) cellarge;
        TPZTransform transform(gelside.Dimension(),gellarge.Dimension());
        gelside.SideTransform3(gellarge,transform);
        
        int anelindex = cellarge->Index();
        if (anelindex < 0 || anelindex >=  NElements()) {
            anelindex = cellarge->Reference()->Father()->Reference()->Index();
            cout << "TPZCompCloneMesh::ERROR\nFine clone reference called!!\n";
        }
        int index = GetOriginalElementIndex(anelindex);
        REAL truerror = 0.;
        REAL erro =  ElementError(cint,cintlarge,transform,f,truerror);
        ervec[index] += erro;
        if(erro > 1.0e-8 || truerror > 1.e-8) {//CEDRIC
            cout << "index: " << index << "  Erro somado: " << erro;
        }
        if(f){
            if (truerror > 0)  truervec[index]  += truerror;
            if(erro > 1.0e-8 || truerror > 1.e-8) {//CEDRIC
                cout << " erro real " << truerror;
            }
        }
        if(erro > 1.0e-8 || truerror > 1.e-8) cout << endl;
    }
}

int TPZCompCloneMesh::GetOriginalElementIndex(int elindex){
    TPZGeoCloneMesh* gclmesh = (TPZGeoCloneMesh *) fReference;
    TPZGeoEl *gel = fElementVec[elindex]->Reference();//gclmesh->ReferenceElement(elindex);
    if (!gel) return -1;
    TPZGeoEl * gelref = gclmesh->ReferenceElement(gclmesh->Index(gel));
    
    if (!gelref) return -1;
    TPZCompEl *cel = gelref->Reference();
    if (!cel) return -1;
    int orgindex = cel->Index();
    return orgindex;
}

int TPZCompCloneMesh::IsFather(TPZGeoEl *el){
    TPZGeoEl *father = el;
    TPZGeoCloneMesh *clgmesh = (TPZGeoCloneMesh *) fReference;
    TPZGeoEl *georef = clgmesh->GetMeshReferenceElement();
    while (father){
        if (father == georef) return 1;
        father = father->Father();
    }
    return 0;
}

REAL TPZCompCloneMesh::ElementError(TPZInterpolatedElement *fine, TPZInterpolatedElement *coarse, TPZTransform &tr,
                                    void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix &deriv),REAL &truerror){
    // accumulates the transfer coefficients between the current element and the
    // coarse element into the transfer matrix, using the transformation t
    int locnod = fine->NConnects();
    int cornod = coarse->NConnects();
    int locmatsize = fine->NShapeF();
    int cormatsize = coarse->NShapeF();
    REAL error = 0.;
    truerror = 0.;
    
    REAL loclocmatstore[500] = {0.},loccormatstore[500] = {0.};
    TPZFMatrix loclocmat(locmatsize,locmatsize,loclocmatstore,500);
    TPZFMatrix loccormat(locmatsize,cormatsize,loccormatstore,500);
    
    TPZAutoPointer<TPZIntPoints> intrule = fine->GetIntegrationRule().Clone();
    int dimension = fine->Dimension();
    int numdof = fine->Material()->NStateVariables();
    TPZBlock &locblock = fine->Mesh()->Block();
    TPZFMatrix &locsolmesh = fine->Mesh()->Solution();
    
    TPZBlock &corblock = coarse->Mesh()->Block();
    TPZFMatrix &corsolmesh = coarse->Mesh()->Solution();
    
    TPZVec<REAL> locsol(numdof);
    TPZFMatrix locdsol(dimension,numdof);
    
    TPZVec<REAL> corsol(numdof);
    TPZFMatrix cordsol(dimension,numdof);
    
    TPZManVector<int> prevorder(dimension),order(dimension);
    intrule->GetOrder(prevorder);
    
    TPZManVector<int> interpolation(dimension);
    fine->GetInterpolationOrder(interpolation);
    
    // compute the interpolation order of the shapefunctions squared
    int dim;
    int maxorder = interpolation[0];
    for(dim=0; dim<interpolation.NElements(); dim++) {
        maxorder = interpolation[dim] < maxorder ? maxorder : interpolation[dim];
    }
    for(dim=0; dim<dimension; dim++) {
        order[dim] = 20;//Cedric
        //order[dim] = maxorder;
    }
    intrule->SetOrder(order);
    
    REAL locphistore[50]={0.},locdphistore[150]={0.};
    TPZFMatrix locphi(locmatsize,1,locphistore,50);
    TPZFMatrix locdphi(dimension,locmatsize,locdphistore,150),locdphix(dimension,locmatsize);
    // derivative of the shape function
    // in the master domain
    
    REAL corphistore[50]={0.},cordphistore[150]={0.};
    TPZFMatrix corphi(cormatsize,1,corphistore,50);
    TPZFMatrix cordphi(dimension,cormatsize,cordphistore,150), cordphix(dimension,cormatsize);
    // derivative of the shape function
    // in the master domain
    REAL jacobianstore[9],axesstore[9];
    TPZManVector<REAL> int_point(dimension),coarse_int_point(dimension);
    TPZFMatrix jacfine(dimension,dimension,jacobianstore,9),jacinvfine(dimension,dimension);
    TPZFMatrix axesfine(3,3,axesstore,9);
    TPZManVector<REAL> xfine(3);
    TPZFMatrix jaccoarse(dimension,dimension),jacinvcoarse(dimension,dimension);
    TPZFMatrix axescoarse(3,3), axesinner(3,3);
    TPZManVector<REAL> xcoarse(3);
    
    REAL jacdetcoarse;
    int numintpoints = intrule->NPoints();
    REAL weight;
    //  int lin,ljn,cjn;
    int i,j,k;
    
    TPZVec<REAL> truesol(numdof);
    TPZFMatrix truedsol(dimension,numdof);
    for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {
        intrule->Point(int_ind,int_point,weight);
        REAL jacdetfine;
        fine->Reference()->Jacobian( int_point, jacfine , axesfine, jacdetfine, jacinvfine);
        fine->Reference()->X(int_point, xfine);
        if(f) f(xfine,truesol,truedsol);
        fine->Shape(int_point,locphi,locdphi);
        tr.Apply(int_point,coarse_int_point);
        coarse->Shape(coarse_int_point,corphi,cordphi);
        coarse->Reference()->Jacobian( coarse_int_point,jaccoarse,axescoarse, jacdetcoarse, jacinvcoarse);
        coarse->Reference()->X(coarse_int_point,xcoarse);
        REAL dist = (xfine[0]-xcoarse[0])*(xfine[0]-xcoarse[0])+(xfine[1]-xcoarse[1])*(xfine[1]-xcoarse[1])+(xfine[2]-xcoarse[2])*(xfine[2]-xcoarse[2]);
        if(dist > 1.e-6) cout << "TPZCompCloneMesh::ElementError transformation between fine and coarse is wrong\t dist = " << dist << endl;
        int dim = axesfine.Rows();
        if (dim != dimension) {
            DebugStop();
        }
        for(i=0; i<dim; i++) {
            for(j=0; j<dim; j++) {
                axesinner(i,j) = 0.;
                for(k=0; k<3; k++) axesinner(i,j) += axesfine(i,k)*axescoarse(j,k);
            }
        }
        
        // Cesar 2003-01-07
        //if(fabs(axesinner(0,0) - 1.) > 1.e-6 || fabs(axesinner(1,1) - 1.) > 1.e-6 || fabs(axesinner(0,1)) > 1.e-6 || fabs(axesinner(1,0)) > 1.e-6) {
        // if(fabs(fabs(axesinner(0,0)) - 1.) > 1.e-6 || fabs(fabs(axesinner(1,1)) - 1.) > 1.e-6 || fabs(fabs(axesinner(0,1))) > 1.e-6 || fabs(fabs(axesinner(1,0))) > 1.e-6) {
        //      cout << "TPZCompCloneMesh::Axesinner is not identify?\n";
        //    }
        weight *= fabs(jacdetfine);
        locdphix.Zero();
        int ieq,d;
        switch(dim) {
            case 0:
                //dphix.Redim(1,1);
                //dphix(0,0) = dphi(0,0);
                break;
            case 1:
                for(d=0; d<dimension; d++) {
                    for(ieq=0; ieq<locmatsize; ieq++) locdphix(d,ieq) = locdphi(d,ieq)*(1./jacdetfine);
                    for(ieq=0; ieq<cormatsize; ieq++) cordphix(d,ieq) = cordphi(d,ieq)*(axesinner(0,0)/jacdetcoarse);
                }
                break;
            case 2:
                for(ieq = 0; ieq < locmatsize; ieq++) {
                    locdphix(0,ieq) = jacinvfine(0,0)*locdphi(0,ieq) + jacinvfine(1,0)*locdphi(1,ieq);
                    locdphix(1,ieq) = jacinvfine(0,1)*locdphi(0,ieq) + jacinvfine(1,1)*locdphi(1,ieq);
                    REAL tmp[2];
                    tmp[0] = locdphix(0,ieq)*axesfine(0,0)+locdphix(1,ieq)*axesfine(1,0);
                    tmp[1] = locdphix(0,ieq)*axesfine(0,1)+locdphix(1,ieq)*axesfine(1,1);
                    locdphix(0,ieq) = tmp[0];
                    locdphix(1,ieq) = tmp[1];
                }
                for(ieq = 0; ieq < cormatsize; ieq++) {
                    cordphix(0,ieq) = jacinvcoarse(0,0)*cordphi(0,ieq) + jacinvcoarse(1,0)*cordphi(1,ieq);
                    cordphix(1,ieq) = jacinvcoarse(0,1)*cordphi(0,ieq) + jacinvcoarse(1,1)*cordphi(1,ieq);
                    REAL tmp[2];
                    tmp[0] = cordphix(0,ieq)*axescoarse(0,0)+cordphix(1,ieq)*axescoarse(1,0);
                    tmp[1] = cordphix(0,ieq)*axescoarse(0,1)+cordphix(1,ieq)*axescoarse(1,1);
                    cordphix(0,ieq) = tmp[0];
                    cordphix(1,ieq) = tmp[1];
                }
                break;
            case 3:
                for(ieq = 0; ieq < locmatsize; ieq++) {
                    locdphix(0,ieq) = jacinvfine(0,0)*locdphi(0,ieq) + jacinvfine(1,0)*locdphi(1,ieq) + jacinvfine(2,0)*locdphi(2,ieq);
                    locdphix(1,ieq) = jacinvfine(0,1)*locdphi(0,ieq) + jacinvfine(1,1)*locdphi(1,ieq) + jacinvfine(2,1)*locdphi(2,ieq);
                    locdphix(2,ieq) = jacinvfine(0,2)*locdphi(0,ieq) + jacinvfine(1,2)*locdphi(1,ieq) + jacinvfine(2,2)*locdphi(2,ieq);
                }
                for(ieq = 0; ieq < cormatsize; ieq++) {
                    cordphix(0,ieq) = jacinvcoarse(0,0)*cordphi(0,ieq) + jacinvcoarse(1,0)*cordphi(1,ieq) + jacinvcoarse(2,0)*cordphi(2,ieq);
                    cordphix(1,ieq) = jacinvcoarse(0,1)*cordphi(0,ieq) + jacinvcoarse(1,1)*cordphi(1,ieq) + jacinvcoarse(2,1)*cordphi(2,ieq);
                    cordphix(2,ieq) = jacinvcoarse(0,2)*cordphi(0,ieq) + jacinvcoarse(1,2)*cordphi(1,ieq) + jacinvcoarse(2,2)*cordphi(2,ieq);
                    REAL tmp[3];
                    tmp[0] = cordphix(0,ieq)*axesinner(0,0)+cordphix(1,ieq)*axesinner(0,1)+cordphix(2,ieq)*axesinner(0,2);
                    tmp[1] = cordphix(0,ieq)*axesinner(1,0)+cordphix(1,ieq)*axesinner(1,1)+cordphix(2,ieq)*axesinner(1,2);
                    tmp[2] = cordphix(0,ieq)*axesinner(2,0)+cordphix(1,ieq)*axesinner(2,1)+cordphix(2,ieq)*axesinner(2,2);
                    cordphix(0,ieq) = tmp[0];
                    cordphix(1,ieq) = tmp[1];
                    cordphix(2,ieq) = tmp[2];
                }
                break;
            default:
                PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
                PZError.flush();
        }
        int iv=0;
        locsol.Fill(0.);
        locdsol.Zero();
        int in;
        for(in=0; in<locnod; in++) {
            TPZConnect *df = &fine->Connect(in);
            int dfseq = df->SequenceNumber();
            int dfvar = locblock.Size(dfseq);
            int pos = locblock.Position(dfseq);
            
            for(int jn=0; jn<dfvar; jn++) {
                locsol[iv%numdof] += locphi(iv/numdof,0)*locsolmesh(pos+jn,0);
                for(d=0; d<dim; d++)
                    locdsol(d,iv%numdof) += locdphix(d,iv/numdof)*locsolmesh(pos+jn,0);
                iv++;
            }
        }
        corsol.Fill(0.);
        cordsol.Zero();
        iv=0;
        for(in=0; in<cornod; in++) {
            TPZConnect *df = &coarse->Connect(in);
            int dfseq = df->SequenceNumber();
            int dfvar = corblock.Size(dfseq);
            int pos = corblock.Position(dfseq);
            for(int jn=0; jn<dfvar; jn++) {
                corsol[iv%numdof] += corphi(iv/numdof,0)*corsolmesh(pos+jn,0);
                for(d=0; d<dim; d++)
                    cordsol(d,iv%numdof) += cordphix(d,iv/numdof)*corsolmesh(pos+jn,0);
                iv++;
            }
        }
        int jn;
        for(jn=0; jn<numdof; jn++) {
            //       error += (locsol[jn]-corsol[jn])*(locsol[jn]-corsol[jn])*weight;
            //       if(f) truerror += (corsol[jn]-truesol[jn])*(corsol[jn]-truesol[jn])*weight;
            for(d=0; d<dim; d++) {
                error += (locdsol(d,jn)-cordsol(d,jn))*(locdsol(d,jn)-cordsol(d,jn))*weight;
                if(f) truerror += (cordsol(d,jn)-truedsol(d,jn))*(cordsol(d,jn)-truedsol(d,jn))*weight;
            }
        }
    }
    intrule->SetOrder(prevorder);
    return error;
}

void TPZCompCloneMesh::ApplyRefPattern(REAL minerror, TPZVec<REAL> &ervec, TPZCompMesh *fine,
                                       TPZStack<TPZGeoEl *> &gelstack, TPZStack<int> &porder){
    int i;
    TPZGeoCloneMesh *gclmesh = (TPZGeoCloneMesh *) Reference();
    
    gclmesh->ResetReference();
    //  LoadReferences();
    fine->LoadReferences();
    
    TPZAutoPointer<TPZMaterial> mat = MaterialVec().rbegin()->second;
    int nstate = mat->NStateVariables();
    TPZOneDRef fn (nstate);
    
    /*   TPZStack<TPZGeoEl *> gelstack; */
    /*   TPZStack<int> porders; */
    
    //  int test = gclmesh->NReference();
    
    for (i=0;i<NElements();i++){
        /*     TPZGeoEl *clgel =  gclmesh->ReferenceElement(i); */
        /*     if (!clgel){ */
        /*       cout << "TPZCompCloneMesh::ApplyRefPattern ERROR!\nFine mesh is loaded, pattern analysis needs coarse mesh!"; */
        /*       return; */
        /*     } */
        /*     TPZCompEl *orgel = clgel->Reference(); */
        TPZCompEl *cel = ElementVec()[i];
        TPZGeoEl *gel = cel->Reference();
        
        //if the element is a reference element ...
        if (!IsFather(gel)) continue;
        int anelindex = cel->Index();
        if (anelindex < 0 || anelindex >=  NElements()) {
            anelindex = cel->Reference()->Father()->Reference()->Index();
        }
        int orgelindex = GetOriginalElementIndex(anelindex);
        //    int orgelindex = orgel->Index();
        REAL curerror = ervec[orgelindex];
        if (curerror >= minerror){
            TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (ElementVec()[i]);
            //  TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (fCloneReference->ElementVec()[orgelindex]);
            if(!cint) {
                if(ElementVec()[i]) {
                    //if (fCloneReference->ElementVec()[orgelindex]){
                    cout << "TPZAdaptMesh RefinementPattern only for interpolated meshes - I don't understand\n";
                    //	  cout << " i = " << i  << " nelem = " << nelem << endl;
                }
                continue;
            }
            //      cout << "analysing element " << i << endl;
            AnalyseElement(fn,cint,gelstack,porder);//gelstack para a malha original....
        }
        else {
            TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (fCloneReference->ElementVec()[orgelindex]);
            int ncon = cint->NConnects();
            int po = cint->PreferredSideOrder(ncon-1);
            //    cout << "Elemento nao refinado " << endl;
            //
            
            //cint->Reference()->Print();
            gelstack.Push(cint->Reference());
            porder.Push(po);
        }
    }
    //Com os parâmetros de analyseElement cria uma malha refinada
    //  AdaptElements(gelstack,porder);//CreateRefineMesh;
}

void TPZCompCloneMesh::AnalyseElement( TPZOneDRef &f, TPZInterpolatedElement *cint,
                                      TPZStack<TPZGeoEl *> &subels, TPZStack<int> &porders) {
    
    //obtencão do elemento geométrico de cint
    TPZGeoEl *gel = cint->Reference();
    //número de conects de cint
    int ncon = cint->NConnects();
    //ordem do elemento
    int intorder = cint->SideOrder(ncon-1);
    gDeduce << "Internal order = " << intorder << endl;
    //Identifica a malha geométrica de geo
    TPZGeoMesh *gmesh = gel->Mesh();
    //número de nós
    int ncorners = gel->NCornerNodes();
    TPZVec<int> cornerid(ncorners);
    TPZVec<int> cornerindexes(ncorners);
    //ordem p para cada nó
    TPZVec<int> localporders(ncorners);
    int ids;
    for(ids=0; ids<ncorners; ids++) {
        cornerid[ids] = gel->NodePtr(ids)->Id();
        cornerindexes[ids] = gel->NodeIndex(ids);
    }
    int n1dsides = 0;
    int nsides = gel->NSides();
    int side;
    
    int maxp = TPZOneDRef::gMaxP;
    
    //Calcula o número de arestas
    for(side=0; side<nsides; side++) if(gel->SideDimension(side) == 1) n1dsides++;
    
    //Vetor de padrões de refinamento com dimensão igual ao
    //número de arestas refinamento unidimensional
    TPZVec<TPZRefPattern> refpattern(n1dsides);
    n1dsides = 0;
    for(side=0; side<nsides; side++) {
        //só considera as arestas
        int sidedimension = gel->SideDimension(side);
        if(sidedimension != 1) continue;
        //obtém a ordem do elemento
        //    int level = gel->Level();
        TPZStack<TPZCompElSide> subelsides;
        TPZStack<TPZCompElSide> auxsubelsides;
        TPZGeoElSide gelside(gel,side);
        //     if(gelside.Neighbour().Exists()) gelside = gelside.Neighbour();
        //     else {
        //       cout << "TPZAnalyseElement coarse element without neighbour\n";
        //       continue;
        //     }
        //obtém a lista de elementos computacionais que são
        //derivados pelos lados subelsides de gelside
        //	gel->GetSubElements2(side,subelsides,1);
        //gelside.HigherLevelCompElementList2(subelsides,1,1);
        gelside.HigherLevelCompElementList2(auxsubelsides,1,1);
        //    gelside.SmallConnect(level+1,subelsides,1);
        int nsubel = auxsubelsides.NElements();
        int i;
        for (i=0;i<nsubel;i++){
            if (auxsubelsides[i].Reference().Dimension() == 1){
                subelsides.Push(auxsubelsides[i]);
            }
        }
        if(subelsides.NElements() != 2 ||
           subelsides[0].Reference().Element()->Father() != gel ||
           subelsides[1].Reference().Element()->Father() != gel ) {
            cout << "A one dimensional side with more than one subelement or inconsistent mesh\n";
            subelsides.Resize(0);
            cint->Reference()->Print(cout);
            int is;
            int nsub = cint->Reference()->NSubElements();
            for (is=0;is<nsub;is++){
                cint->Reference()->SubElement(is)->Print(cout);
            }
            gelside.HigherLevelCompElementList2(subelsides,1,1);
            //      gelside.SmallConnect(level+1,subelsides,1);
            continue;
        }
        
        TPZGeoElSide gels1 = subelsides[0].Reference();
        TPZGeoElSide gels2 = subelsides[1].Reference();
        //verifica a ordem dos nós do centro da  aresta
        if(gels1.SideNodeIndex(1) != gels2.SideNodeIndex(0)) {
            TPZCompElSide temp = subelsides[0];
            subelsides[0]=subelsides[1];
            subelsides[1] = temp;
            gels1 = subelsides[0].Reference();
            gels2 = subelsides[1].Reference();
        }
        
        if(gels1.SideNodeIndex(1) != gels2.SideNodeIndex(0)) {
            cout << "Unexpected situation\n";
            cout << gels1.SideNodeIndex(0) << ' ' << gels1.SideNodeIndex(1) << ' '
            << gels2.SideNodeIndex(0) << ' ' << gels2.SideNodeIndex(1) <<endl;
            cout << "Element Analysed: " << cint->Index() << endl;
            cint->Reference()->Print();
            int is;
            int nsub = cint->Reference()->NSubElements();
            cout <<"SubElements\n";
            for (is=0;is<nsub;is++){
                cout << "Subelement: " << is << endl;
                cint->Reference()->SubElement(is)->Print(cout);
            }
            cout << "Subelement 1:" << endl;
            gels1.Element()->Print();
            cout << "Subelement 2:" << endl;
            gels2.Element()->Print();
            cint->Mesh()->SetName("Unconsistent mesh");
            cint->Mesh()->Print(cout);
            continue;
        }
        
        TPZVec<int> index(3),id(3);
        index[0] = gels1.SideNodeIndex(0);
        id[0] = gmesh->NodeVec()[index[0]].Id();
        index[1] = gels1.SideNodeIndex(1);
        id[1] = gmesh->NodeVec()[index[1]].Id();
        index[2] = gels2.SideNodeIndex(1);
        id[2] = gmesh->NodeVec()[index[2]].Id();
        //     cout << "side = " << side << " indexes " << id[0] << ' ' << id[1] << ' ' << id[2] << endl;
        REAL del[3];
        //calcula a distancia entre os nos da aresta
        for(i=0; i<3; i++) del[i] = gmesh->NodeVec()[index[1]].Coord(i)-gmesh->NodeVec()[index[0]].Coord(i);
        REAL delx = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);
        
        //obtém a lista de connects da aresta
        TPZCompMesh *cmesh = subelsides[0].Element()->Mesh();
        TPZConnect *connects[5];
        TPZInterpolatedElement *c1 = (TPZInterpolatedElement *) subelsides[0].Element();
        int s1 = subelsides[0].Side();
        TPZInterpolatedElement *c2 = (TPZInterpolatedElement *) subelsides[1].Element();
        int s2 = subelsides[1].Side();
        connects[0] = c1->SideConnect(0,s1);
        connects[2] = c1->SideConnect(1,s1);
        connects[1] = c1->SideConnect(2,s1);
        connects[4] = c2->SideConnect(1,s2);
        connects[3] = c2->SideConnect(2,s2);
        
        //calcula o número de graus de liberdade da aresta
        int dof = 0;
        for(i=0; i<5; i++) dof += connects[i]->NDof(*cmesh);         //??
        TPZFMatrix U(dof,1);
        int p1 = c1->SideOrder(s1);
        int p2 = c2->SideOrder(s2);
        //    p1 = p1 > maxp ? maxp : p1;
        //   p2 = p2 > maxp ? maxp : p2;
        
        gDeduce << "side = " << side << " order = " << p1 << ' ' << p2 << endl;
        int c,counter=0;
        for(c=0; c<5; c++) {
            int nd = connects[c]->NDof(*cmesh);
            int seq = connects[c]->SequenceNumber();
            int d;
            //obtem o bloco de cada connect e o coloca na matriz U
            for(d=0; d<nd; d++) U(counter++,0) = (*cmesh).Block()(seq,0,d,0);
        }
        int hp1, hp2;
        REAL hperror;
        
        //Para cada aresta é calculado o menor erro através do cálculo do refinamento unidimensional
        REAL error = f.BestPattern(U,id,p1,p2,hp1, hp2, hperror,delx);
        //define o refinamento para o elemento??
        TPZRefPattern optimal = {id[0],id[1],id[2],p1,p2,hp1,hp2,hperror,error};
        refpattern[n1dsides] = optimal;
        n1dsides++;
    }
    refpattern.Resize(n1dsides);
    //gDeduceRefPattern(refpattern,cornerids,localporders,intorder);
    DeduceRefPattern(refpattern,cornerid,localporders,intorder);
    
    TPZGeoCloneMesh *gclmesh = dynamic_cast<TPZGeoCloneMesh *> (Reference());
    int gelindex =  gclmesh->Index(gel);
    TPZGeoEl *orgel = gclmesh->ReferenceElement(gelindex);
    
    if(localporders[1] == -1) {
        //    cout << "Inserindo elemento com refinamento p\n";
        //    orgel->Print();
        subels.Push(orgel);
        porders.Push(localporders[0]);
        //     cout << "prefine order " << localporders[0] << endl;
        return;
    }
    TPZStack<TPZGeoElSide> gelsides;
    
    TPZVec<TPZGeoEl *> orgelsub;
    orgel->Divide(orgelsub);
    
    for (ids=0;ids<ncorners;ids++) cornerindexes[ids]=orgel->NodeIndex(ids);
    
    // gel->GetSubElement(nsides-1,cornerids,gelsides);
    int dim = orgel->Dimension();
    orgel->GetSubElements2(nsides-1,gelsides,dim);
    //  orgel->Print(cout);
    //  for (ids=0;ids<ncorners;ids++){
    //    cout << cornerindexes[ids] << "\t";
    //  }
    //  cout << endl;
    for(ids=0; ids<ncorners; ids++) {
        //    cout << "Subelemento inserido " << endl;
        //    gelsides[ids].Element()->Print();
        subels.Push(gelsides[ids].Element());
        porders.Push(localporders[ids]);
    }
    // hook for triangular, tetrahedral, prism and pyramid elements
    if(gelsides.NElements() > ncorners) {
        // compute the maximum interpolation order of all "corner" elements"
        int maxorder = 0;
        for(ids=0; ids<ncorners; ids++) maxorder = (maxorder<localporders[ids]) ? localporders[ids] : maxorder;
        // Assign the maximum order to all remaining elements
        for(ids=ncorners; ids<gelsides.NElements(); ids++) {
            subels.Push(gelsides[ids].Element());
            porders.Push(maxorder);
        }
    }
}

void TPZCompCloneMesh::DeduceRefPattern(TPZVec<TPZRefPattern> &refpat,
                                        TPZVec<int> &cornerids,
                                        TPZVec<int> &porders,
                                        int originalp) {
    
    // Eliminate the refinement pattern suggestion if
    // the error is smaller than 10% of the total error
    int nref = refpat.NElements();
    REAL totalerror = 0.;
    int ir;
    
    // cálcula o erro total do elemento -
    // somatório dos erros nas arestas
    for(ir=0; ir<nref; ir++) {
        totalerror += refpat[ir].fError;
    }
    
    //verifica o número de nós de canto
    int ncorners = cornerids.NElements();
    
    //Print the incoming refpattern to the log file
    for(ir=0; ir<ncorners; ir++) gDeduce << cornerids[ir] << ' ';
    gDeduce << endl;
    for(ir=0; ir<nref; ir++) {
        gDeduce << refpat[ir].fId[0] << ' ' << refpat[ir].fId[1] << ' ' << refpat[ir].fId[2] << ' ' << refpat[ir].fp[0] <<' ' << refpat[ir].fp[1] <<" error " << refpat[ir].fError;
        gDeduce << ' ' << refpat[ir].fh[0] << ' ' << refpat[ir].fh[1] << ' ' << refpat[ir].fhError << endl;
    }
    
    //Desconsidera erros da ordem de 10^-3 do erro total
    for(ir=0; ir<nref; ir++) {
        if(refpat[ir].fError < totalerror*1.e-3) {
            refpat[ir].fp[1] = -1;
            refpat[ir].fp[0] = originalp+1;
        }
    }
    
    //Imprime novamente a malha sem os refinamentos desconsiderados
    gDeduce << "originalp " << endl;
    for(ir=0; ir<nref; ir++) {
        gDeduce << refpat[ir].fId[0] << ' ' << refpat[ir].fId[1] << ' ' << refpat[ir].fId[2] << ' ' << refpat[ir].fp[0] <<' ' << refpat[ir].fp[1] <<" error " << refpat[ir].fError;
        gDeduce << ' ' << refpat[ir].fh[0] << ' ' << refpat[ir].fh[1] << ' ' << refpat[ir].fhError << endl;
    }
    // for each corner id, identify the edges which are connected to it
    // if any edge suggests an h-refinement, use the h-refinement
    // if no edge suggests an h-refinement return a unique parameter p in porders
    // lembrar que quando o não refinamento p fornece o menor erro, p2 = -1!1
    int pref = 1;
    for(ir=0; ir<nref; ir++) if(refpat[ir].fp[1] != -1) pref = 0;
    // if pref is still equal 1, we will use the prefinement
    // determine the maximum p-order suggested and return it.
    if(pref) {
        int maxp = 0;
        //busca o maior p
        for(ir=0; ir<nref; ir++) {
            if(refpat[ir].fp[0] > maxp) maxp = refpat[ir].fp[0];
            if(refpat[ir].fp[1] > maxp) maxp = refpat[ir].fp[1];
        }
        porders.Fill(-1);
        porders[0] = maxp; //apenas um p é retornado
        gDeduce << "prefinement order " << maxp << endl;
        return;
    }
    
    
    TPZVec<int> perm(refpat.NElements());
    TPZVec<REAL> error(refpat.NElements());
    totalerror = 0.;
    for(ir=0; ir<nref; ir++) {
        error[ir] = refpat[ir].fhError;
        totalerror += refpat[ir].fhError;
        perm[ir] = ir;
    }
    Sort(error,perm);
    
    // h-refinement will be used
    // determine the order of interpolation of the sub elements
    int ic;
    gDeduce << "errors in their order ";
    for(ir=0; ir<nref; ir++) gDeduce << perm[ir] << ' ' << refpat[perm[ir]].fhError << ' ';
    gDeduce << endl;
    gDeduce << "h-refinement order of the elements ";
    for(ic=0; ic<ncorners; ic++) porders[ic]=0;
    for(ir=0; ir<nref; ir++) {
        //Despreza os erros pequenos
        if(refpat[perm[ir]].fhError < 1.e-3*totalerror) continue;
        for(ic=0; ic<ncorners; ic++) {
            int id = cornerids[ic];
            //identifica o corner 0 ou 2
            if(refpat[perm[ir]].fId[0] == id) {
                //porders[ic] não está sendo zerado no início do for???
				if(porders[ic] == 0) porders[ic] = refpat[perm[ir]].fh[0];
            }
            if(refpat[perm[ir]].fId[2] == id) {
                if(porders[ic] == 0) porders[ic] = refpat[perm[ir]].fh[1];
            }
        }
    }
    for(ic=0; ic<ncorners; ic++) {
        //para os lados os corners que não precisavam de
        //refinamento, sua ordem é aumentada
        if(porders[ic] == 0) porders[ic] = originalp/2+1;
        gDeduce << porders[ic] << ' ';
    }
    gDeduce << endl;
}

void TPZCompCloneMesh::AdaptElements(TPZVec<TPZGeoEl *> &gelstack,TPZVec<int> &porders) {
    
    //Idenifica o vetor de elementos computacionais de mesh
    //  TPZAdmChunkVector<TPZCompEl *> &elementvec = ElementVec();
    
    int el,nelem = gelstack.NElements();
    
    Reference()->ResetReference();
    LoadReferences();
    TPZGeoCloneMesh *gclmesh = (TPZGeoCloneMesh *)Reference();
    
    for(el=0; el<nelem; el++) {
        //identifica os elementos geométricos passados em gelstack
        //    TPZGeoEl *gel = gelstack[el];
        TPZGeoEl *clgel = gelstack[el];
        int clelindex = gclmesh->Index(clgel);
        TPZGeoEl *gel = gclmesh->ReferenceElement(clelindex);
        
        if(!gel) {
            cout << "TPZCompCloneMesh::CreateCompMesh encountered an null element\n";
            continue;
        }
        //int celindex;
        
        //Cria um TPZIntel baseado no gel identificado
        //TPZInterpolatedElement *csint;
        //    csint = dynamic_cast<TPZInterpolatedElement *> (gel->CreateCompEl(*cmesh,celindex));
        //    csint = dynamic_cast<TPZInterpolatedElement *> (gel->CreateCompEl(*fCloneReference,celindex));
        //if(!csint) continue;
        
        //Refina em p o elemento criado
        //csint->PRefine(porders[el]);
    }
    
    //Mais einh!!
    //  cmesh->AdjustBoundaryElements();
    fCloneReference->AdjustBoundaryElements();
    //  return cmesh;
}

void TPZCompCloneMesh::Sort(TPZVec<REAL> &vec, TPZVec<int> &perm) {
    int i,j;
    int imin = 0;
    int imax = vec.NElements();
    for(i=imin; i<imax; i++) {
        for(j=i+1; j<imax; j++) {
            if(vec[perm[i]] < vec[perm[j]]) {
                int kp = perm[i];
                perm[i] = perm[j];
                perm[j] = kp;
            }
        }
    }
}


void TPZCompCloneMesh::Print (ostream & out) {
    
    out <<  "\n\t\tCOMPUTABLE CLONE GRID INFORMATIONS:\n\n";
    out << "\tREFERENCE MESH:\t" << fCloneReference->Name() << endl;
    ComputeNodElCon();
    out << "\n\t\tCOMPUTABLE GRID INFORMATIONS:\n\n";
    out << "TITLE-> " << fName << "\n\n";
    
    out << "number of connects            = " << NConnects() << endl;
    out << "number of elements            = " << NElements() << endl;
    out << "number of materials           = " << NMaterials() << endl;
    
    out << "\n\t Cloned Connect Information:\n\n";
    int i, nelem = NConnects();
    for(i=0; i<nelem; i++) {
        if(fConnectVec[i].SequenceNumber() == -1) {
            if(fConnectVec[i].HasDependency()) {
                out << "TPZCompMesh::Print inconsistency of connect\n";
                out << "Index " << i << ' ';
                fConnectVec[i].Print(*this,out);
            }
            continue;
        }
        out << "Clone connection - ";
        out << "Index " << i << ' ';
        fConnectVec[i].Print(*this,out);
        int orgcid = fOriginalConnects[i];
        out << "Original connection - " << "Index " << orgcid << " ";
        fCloneReference->ConnectVec()[orgcid].Print(*fCloneReference,out);
    }
    out << "\n\t Computable Element Information:\n\n";
    nelem = NElements();
    for(i=0; i<nelem; i++) {
        if(!fElementVec[i]) continue;
        TPZCompEl *el = fElementVec[i];
        out << "Index " << i << ' ';
        el->Print(out);
        if(!el->Reference()) continue;
        out << "\tReference Id = " << el->Reference()->Id() << endl;
    }
    out << "\n\tMaterial Information:\n\n";
    std::map<int, TPZAutoPointer<TPZMaterial> >::const_iterator it;
    for (it=MaterialVec().begin(); it!=MaterialVec().end(); it++) {
        it->second->Print(out);
    }
}

TPZInterpolatedElement *TPZCompCloneMesh::GetOriginalElement(TPZCompEl *el) {
    int index = el->Index();
    int orgind = GetOriginalElementIndex(index);
    return dynamic_cast<TPZInterpolatedElement *> (fCloneReference->ElementVec()[orgind]);
}
