#include "TPZElastoPlasticIntPointsStructMatrix.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "pzmetis.h"
#include "TPZMyLambdaExpression.h"
#include "TPZElasticCriterion.h"

#ifdef USING_MKL
#include <mkl.h>
#endif


TPZElastoPlasticIntPointsStructMatrix::TPZElastoPlasticIntPointsStructMatrix(TPZCompMesh *cmesh) : TPZSymetricSpStructMatrix(cmesh), fLambdaExp(), fSparseMatrixLinear(), fRhsLinear(), fCoefToGradSol() {

}

TPZElastoPlasticIntPointsStructMatrix::~TPZElastoPlasticIntPointsStructMatrix() {
}

TPZStructMatrix * TPZElastoPlasticIntPointsStructMatrix::Clone(){
    return new TPZElastoPlasticIntPointsStructMatrix(*this);
}

TPZMatrix<STATE> * TPZElastoPlasticIntPointsStructMatrix::Create(){

    if(!isBuilt()) {
        this->SetUpDataStructure();
    }
    
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    fMesh->ComputeElGraph(elgraph,elgraphindex,fMaterialIds);
    TPZMatrix<STATE> * mat = SetupMatrixData(elgraph, elgraphindex);
    return mat;
}

TPZMatrix<STATE> *TPZElastoPlasticIntPointsStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {

    int64_t neq = fMesh->NEquations();
    TPZMatrix<STATE> *stiff = Create(); ///
    TPZSYsmpMatrix<STATE> *mat = dynamic_cast<TPZSYsmpMatrix<STATE> *> (stiff);
    rhs.Redim(neq,1);
    
    
    
    int n_state = 2; /// or dim
    TPZVec<STATE> & K_g = mat->A();
    TPZVec<int> & indexes = fCoefToGradSol.Indexes();
    int64_t n_vols = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fNumBlocks;
    int64_t n_cols = fCoefToGradSol.IrregularBlocksMatrix().Cols();
    TPZVec<int> el_n_dofs = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColSizes; // Making copy here!
    TPZVec<int> mat_indexes = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColColPosition; // Making copy here!
    TPZVec<int> cols_first_index = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColFirstIndex;

    TPZVec<REAL> depxx;
    TPZVec<REAL> depyy;
    TPZVec<REAL> depxy;
    Dep(depxx, depyy, depxy);
    
#ifdef USING_CUDA
    int nblocks = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fNumBlocks;
    TPZVecGPU<REAL> Kxx(fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColColPosition[nblocks]);
    TPZVecGPU<REAL> Kyy(fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColColPosition[nblocks]);
    TPZVecGPU<REAL> Kxy(fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColColPosition[nblocks]);
    
    TPZVecGPU<REAL> d_dep(depxx.size());
    
    d_dep.set(&depxx[0], depxx.size());
    fCoefToGradSol.IrregularBlocksMatrix().KMatrix(d_dep.getData(), Kxx.getData());
    
    d_dep.set(&depyy[0], depyy.size());
    fCoefToGradSol.IrregularBlocksMatrix().KMatrix(d_dep.getData(), Kyy.getData());
    
    d_dep.set(&depxy[0], depxy.size());
    fCoefToGradSol.IrregularBlocksMatrix().KMatrix(d_dep.getData(), Kxy.getData());
#else
    int nblocks = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fNumBlocks;
    TPZVec<REAL> Kxx(fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColColPosition[nblocks]);
    TPZVec<REAL> Kyy(fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColColPosition[nblocks]);
    TPZVec<REAL> Kxy(fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColColPosition[nblocks]);
    
    fCoefToGradSol.IrregularBlocksMatrix().KMatrix(&depxx[0], &Kxx[0]); /// Can
    fCoefToGradSol.IrregularBlocksMatrix().KMatrix(&depyy[0], &Kyy[0]);
    fCoefToGradSol.IrregularBlocksMatrix().KMatrix(&depxy[0], &Kxy[0]);

    std::cout << Kxx << std::endl;

#endif
    
    /// implement OptV1
    if(1){
        /// Serial
        for (int iel = 0; iel < n_vols; iel++) {
            int el_dof = el_n_dofs[iel];
            int pos = cols_first_index[iel];
            int mat_pos = mat_indexes[iel];

            TPZFMatrix<REAL> Kel(n_state*el_dof, n_state*el_dof, 0.);
//            TPZManVector<REAL,64> kg_el(el_dof*el_dof*n_state*n_state);
            
            for (int i_dof = 0; i_dof < el_dof; i_dof++){
                int i_dest_1 = indexes[pos+i_dof];
                int j_dest_1 = indexes[pos+i_dof+n_cols];
                
                for (int j_dof = 0; j_dof < el_dof; j_dof++){
                    Kel.PutVal(2*i_dof, 2*j_dof, Kxx[mat_pos+i_dof*el_dof+j_dof]);
                    Kel.PutVal(2*i_dof + 1, 2*j_dof + 1, Kyy[mat_pos+i_dof*el_dof+j_dof]);
                    Kel.PutVal(2*i_dof, 2*j_dof + 1, Kxy[mat_pos+i_dof*el_dof+j_dof]);
                    Kel.PutVal(2*i_dof + 1, 2*j_dof, Kxy[mat_pos+i_dof+j_dof*el_dof]);
                    
                    STATE val_xx = Kxx[mat_pos+i_dof*el_dof+j_dof];
                    STATE val_yy = Kyy[mat_pos+i_dof*el_dof+j_dof];
                    STATE val_xy = Kxy[mat_pos+i_dof*el_dof+j_dof];
                    STATE val_yx = Kxy[mat_pos+i_dof+j_dof*el_dof];
//                    std::cout << "val_xx = " << val_xx << std::endl;
                    int i_dest_2 = indexes[pos+j_dof];
                    int j_dest_2 = indexes[pos+j_dof+n_cols];
                    
                    val_xx += stiff->GetVal(i_dest_1, i_dest_2);
                    val_yy += stiff->GetVal(j_dest_1, j_dest_2);
                    val_xy += stiff->GetVal(i_dest_1, j_dest_2);
                    val_yx += stiff->GetVal(j_dest_2, i_dest_1);

                    stiff->PutVal(i_dest_1, i_dest_2, val_xx);
                    stiff->PutVal(j_dest_1, j_dest_2, val_yy);
                    stiff->PutVal(i_dest_1, j_dest_2, val_xy);
                    stiff->PutVal(j_dest_2, i_dest_1, val_yx);
                }
            }
//            TPZCompEl *cel = fMesh->Element(iel);
//            TPZElementMatrix ek(fMesh,TPZElementMatrix::EK);
//            TPZElementMatrix ef(fMesh,TPZElementMatrix::EF);
//            cel->CalcStiff(ek, ef);
//
//            TPZFMatrix<REAL> res(n_state*el_dof, n_state*el_dof);
//            ek.fMat.Print(std::cout);
//            Kel.Print(std::cout);
//            res = ek.fMat - Kel;
//            res.Print(std::cout);
        }
    
        std::ofstream out("kg.txt");
        stiff->Print("kg = ",out,EMathematicaInput);
        out.flush();
        
    }
    

    
    Assemble(*stiff,rhs,guiInterface);
//    std::ofstream out("Kref.txt");
//    stiff->Print("Kref = ",out,EMathematicaInput);
//    out.flush();
    mat->ComputeDiagonal();
    return stiff;
}

void TPZElastoPlasticIntPointsStructMatrix::SetUpDataStructure() {

    if(isBuilt()) {
        std::cout << __PRETTY_FUNCTION__ << " Data structure has been setup." << std::endl;
        return;
    }
    
    TPZStack<int> elindex_domain;
    std::set<int> boundary_matids;
    this->GetDomainElements(elindex_domain, boundary_matids); // Candidate to tbb or openmp

    TPZIrregularBlocksMatrix::IrregularBlocks blocksData;
    this->SetUpIrregularBlocksData(elindex_domain, blocksData);

    int64_t rows = blocksData.fRowFirstIndex[blocksData.fNumBlocks];
    int64_t cols = blocksData.fColFirstIndex[blocksData.fNumBlocks];
    TPZIrregularBlocksMatrix blocksMatrix(rows, cols);
    blocksMatrix.SetBlocks(blocksData);
    fCoefToGradSol.SetIrregularBlocksMatrix(blocksMatrix);

    TPZVec<int> indexes;
    this->SetUpIndexes(elindex_domain, indexes);
    fCoefToGradSol.SetIndexes(indexes);

    TPZVec<int> coloredindexes;
    int ncolor;
    this->ColoredIndexes(elindex_domain, indexes, coloredindexes, ncolor);
    fCoefToGradSol.SetIndexesColor(coloredindexes);
    fCoefToGradSol.SetNColors(ncolor);

    AssembleBoundaryData(boundary_matids);

#ifdef USING_CUDA
    std::cout << "Transfering data to GPU..." << std::endl;
    fCoefToGradSol.TransferDataToGPU();
    std::cout << "Done!" << std::endl;
#endif

}

void TPZElastoPlasticIntPointsStructMatrix::Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    TPZSymetricSpStructMatrix::Assemble(mat,rhs, guiInterface);

//    auto it_end = fSparseMatrixLinear.MapEnd();
//    
//    for (auto it = fSparseMatrixLinear.MapBegin(); it!=it_end; it++) {
//        int64_t row = it->first.first;
//        int64_t col = it->first.second;
//        STATE val = it->second;
//        STATE vol_val = mat.GetVal(row, col);
//        vol_val += val; /// TODO:: Add val
//        mat.PutVal(row, col, vol_val);
//    }

    rhs+=fRhsLinear;
}

void TPZElastoPlasticIntPointsStructMatrix::Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){

    if(!isBuilt()) {
        this->SetUpDataStructure();
    }

    int neq = fMesh->NEquations();

#ifdef USING_CUDA
    TPZVecGPU<REAL> solution(neq);
    solution.set(&fMesh->Solution()(0,0), neq);

    TPZVecGPU<REAL> dgrad_u;
    TPZVecGPU<REAL> drhs(neq);
    rhs.Resize(neq, 1);
    rhs.Zero();

    fCoefToGradSol.Multiply(solution, dgrad_u);

    TPZFMatrix<REAL> grad_u(dgrad_u.getSize(),1);
    TPZFMatrix<REAL> sigma;
    dgrad_u.get(&grad_u(0,0), dgrad_u.getSize());

    fLambdaExp.ComputeSigma(grad_u, sigma);

    TPZVecGPU<REAL> dsigma(sigma.Rows());
    dsigma.set(&sigma(0,0), sigma.Rows());

    fCoefToGradSol.MultiplyTranspose(dsigma, drhs);
    drhs.get(&rhs(0,0), neq);
#else
    TPZFMatrix<REAL> grad_u;
    TPZFMatrix<REAL> sigma;
    rhs.Resize(neq, 1);
    rhs.Zero();

    fCoefToGradSol.Multiply(fMesh->Solution(), grad_u);
    fLambdaExp.ComputeSigma(grad_u, sigma);
    fCoefToGradSol.MultiplyTranspose(sigma, rhs);
#endif
    rhs += fRhsLinear;



}

void TPZElastoPlasticIntPointsStructMatrix::AssembleBoundaryData(std::set<int> &boundary_matids) {
    int64_t neq = fMesh->NEquations();

    TPZStructMatrix str(fMesh);
    str.SetMaterialIds(boundary_matids);
    TPZAutoPointer<TPZGuiInterface> guiInterface;
    fRhsLinear.Resize(neq, 1);
    fRhsLinear.Zero();
    fSparseMatrixLinear.Resize(neq, neq);
    str.Assemble(fSparseMatrixLinear, fRhsLinear, guiInterface);
}

void TPZElastoPlasticIntPointsStructMatrix::GetDomainElements(TPZStack<int> &elindex_domain, std::set<int> &boundary_matids) {
    boundary_matids.clear();
    int dim = fMesh->Dimension();
    for (int64_t i = 0; i < fMesh->NElements(); i++) {
        TPZCompEl *cel = fMesh->Element(i);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;
        int mat_id = gel->MaterialId();
        if(gel->Dimension() == dim){
            fMaterialIds.insert(mat_id);
            elindex_domain.Push(cel->Index());
        } else {
            boundary_matids.insert(mat_id);
        }
    }
}

void TPZElastoPlasticIntPointsStructMatrix::SetUpIrregularBlocksData(TPZStack<int> &elindex_domain, TPZIrregularBlocksMatrix::IrregularBlocks &blocksData) {
    int nblocks = elindex_domain.size();

    blocksData.fNumBlocks = nblocks;
    blocksData.fRowSizes.resize(nblocks);
    blocksData.fColSizes.resize(nblocks);
    blocksData.fMatrixPosition.resize(nblocks + 1);
    blocksData.fRowFirstIndex.resize(nblocks + 1);
    blocksData.fColFirstIndex.resize(nblocks + 1);

    blocksData.fMatrixPosition[0] = 0;
    blocksData.fRowFirstIndex[0] = 0;
    blocksData.fColFirstIndex[0] = 0;

    blocksData.fRowRowPosition.resize(nblocks + 1);
    blocksData.fColColPosition.resize(nblocks + 1);
    blocksData.fRowRowPosition[0] = 0;
    blocksData.fColColPosition[0] = 0;

    int64_t rows = 0;
    int64_t cols = 0;
    for(int iel = 0; iel < nblocks; iel++) {
        TPZCompEl *cel = fMesh->Element(elindex_domain[iel]);
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element
        int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

        blocksData.fRowSizes[iel] = dim * npts;
        blocksData.fColSizes[iel] = nf;
        blocksData.fMatrixPosition[iel + 1] = blocksData.fMatrixPosition[iel] + blocksData.fRowSizes[iel] * blocksData.fColSizes[iel];
        blocksData.fRowFirstIndex[iel + 1] =  blocksData.fRowFirstIndex[iel] + blocksData.fRowSizes[iel];
        blocksData.fColFirstIndex[iel + 1] = blocksData.fColFirstIndex[iel] + blocksData.fColSizes[iel];

        blocksData.fRowRowPosition[iel + 1] = blocksData.fRowRowPosition[iel] + blocksData.fRowSizes[iel] * blocksData.fRowSizes[iel];
        blocksData.fColColPosition[iel + 1] = blocksData.fColColPosition[iel] + blocksData.fColSizes[iel] * blocksData.fColSizes[iel];

        rows += blocksData.fRowSizes[iel];
        cols += blocksData.fColSizes[iel];
    }

    blocksData.fStorage.resize(blocksData.fMatrixPosition[nblocks]);

    for (int iel = 0; iel < nblocks; ++iel) {
        TPZCompEl *cel = fMesh->Element(elindex_domain[iel]);
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int row_el = blocksData.fRowSizes[iel];
        int col_el = blocksData.fColSizes[iel];
        int pos_el = blocksData.fMatrixPosition[iel];
        int dim = cel->Dimension();

        TPZFMatrix<REAL> elmatrix;
        elmatrix.Resize(row_el, col_el);

        TPZMaterialData data;
        cel_inter->InitMaterialData(data);

        for (int64_t ipts = 0; ipts < row_el / dim; ipts++) {
            TPZVec<REAL> qsi(dim);
            REAL w;
            int_rule->Point(ipts, qsi, w);
            cel_inter->ComputeRequiredData(data, qsi);

            TPZFMatrix<REAL> dphiXY;
            data.axes.Transpose();
            data.axes.Multiply(data.dphix, dphiXY);

            for (int inf = 0; inf < col_el; inf++) {
                for (int idim = 0; idim < dim; idim++)
                    elmatrix(ipts * dim + idim, inf) = dphiXY(idim, inf);
            }
        }
        elmatrix.Transpose(); // Using CSR format
        TPZFMatrix<REAL> elmatloc(row_el, col_el, &blocksData.fStorage[pos_el], row_el * col_el);
        elmatloc = elmatrix;
    }
}

void TPZElastoPlasticIntPointsStructMatrix::SetUpIndexes(TPZStack<int> &elindex_domain, TPZVec<int> & dof_indexes) {
    int64_t nblocks = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fNumBlocks;
    int64_t rows = fCoefToGradSol.IrregularBlocksMatrix().Rows();
    int64_t cols = fCoefToGradSol.IrregularBlocksMatrix().Cols();

    dof_indexes.resize(fMesh->Dimension() * cols);
    TPZVec<REAL> weight(rows / fMesh->Dimension());

    int64_t cont1 = 0;
    int64_t cont2 = 0;
    int64_t wit = 0;
    for (int iel = 0; iel < nblocks; ++iel) {
        TPZCompEl *cel = fMesh->Element(elindex_domain[iel]);
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element

        TPZMaterialData data;
        cel_inter->InitMaterialData(data);

        for (int64_t ipts = 0; ipts < npts; ipts++) {
            TPZVec<REAL> qsi(dim);
            REAL w;
            int_rule->Point(ipts, qsi, w);
            cel_inter->ComputeRequiredData(data, qsi);
            weight[wit] = w * std::abs(data.detjac);
            wit++;
        }

        int64_t ncon = cel->NConnects();
        for (int64_t icon = 0; icon < ncon; icon++) {
            int64_t id = cel->ConnectIndex(icon);
            TPZConnect &df = fMesh->ConnectVec()[id];
            int64_t conid = df.SequenceNumber();
            if (df.NElConnected() == 0 || conid < 0 || fMesh->Block().Size(conid) == 0) continue;
            else {
                int64_t pos = fMesh->Block().Position(conid);
                int64_t nsize = fMesh->Block().Size(conid);
                for (int64_t isize = 0; isize < nsize; isize++) {
                    if (isize % 2 == 0) {
                        dof_indexes[cont1] = pos + isize;
                        cont1++;
                    } else {
                        dof_indexes[cont2 + cols] = pos + isize;
                        cont2++;
                    }
                }
            }
        }
    }

    TPZMaterial *material = fMesh->FindMaterial(1);
    fLambdaExp.SetMaterial(material);
    fLambdaExp.SetIntPoints(rows / fMesh->Dimension());
    fLambdaExp.SetWeightVector(weight);
}

void TPZElastoPlasticIntPointsStructMatrix::ColoredIndexes(TPZStack<int> &elindex_domain, TPZVec<int> &indexes, TPZVec<int> &coloredindexes, int &ncolor) {
    int64_t nblocks = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fNumBlocks;
    int64_t cols = fCoefToGradSol.IrregularBlocksMatrix().Cols();

    TPZVec<int64_t> connects_vec(fMesh->NConnects(),0);
    TPZVec<int64_t> elemcolor(nblocks,-1);

    int64_t contcolor = 0;
    bool needstocontinue = true;

    while (needstocontinue)
    {
        int it = 0;
        needstocontinue = false;
        for (auto iel : elindex_domain) {
            TPZCompEl *cel = fMesh->Element(iel);
            if (!cel || cel->Dimension() != fMesh->Dimension()) continue;

            it++;
            if (elemcolor[it-1] != -1) continue;

            TPZStack<int64_t> connectlist;
            fMesh->Element(iel)->BuildConnectList(connectlist);
            int64_t ncon = connectlist.size();

            int64_t icon;
            for (icon = 0; icon < ncon; icon++) {
                if (connects_vec[connectlist[icon]] != 0) break;
            }
            if (icon != ncon) {
                needstocontinue = true;
                continue;
            }
            elemcolor[it-1] = contcolor;
            for (icon = 0; icon < ncon; icon++) {
                connects_vec[connectlist[icon]] = 1;
            }
        }
        contcolor++;
        connects_vec.Fill(0);
    }

    ncolor = contcolor;
    coloredindexes.resize(fMesh->Dimension() * cols);
    int64_t neq = fMesh->NEquations();
    for (int64_t iel = 0; iel < nblocks; iel++) {
        int64_t elem_col = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColSizes[iel];
        int64_t cont_cols = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fColFirstIndex[iel];

        for (int64_t icols = 0; icols < elem_col; icols++) {
            coloredindexes[cont_cols + icols] = indexes[cont_cols + icols] + elemcolor[iel]*neq;
            coloredindexes[cont_cols + cols + icols] = indexes[cont_cols + cols + icols] + elemcolor[iel]*neq;
        }
    }
}

void TPZElastoPlasticIntPointsStructMatrix::Dep(TPZVec<REAL> &depxx, TPZVec<REAL> &depyy, TPZVec<REAL> &depxy) {
    int nblocks = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fNumBlocks;
    int sizedep = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fRowRowPosition[nblocks];
    // for (int i = 0; i < nblocks; ++i) {
    //     int rows = fCoefToGradSol.IrregularBlocksMatrix().Blocks().fRowSizes[i];
    //     sizedep += rows * rows;
    // }

    depxx.resize(sizedep);
    depyy.resize(sizedep);
    depxy.resize(sizedep);
    depxx.Fill(0.);
    depyy.Fill(0.);
    depxy.Fill(0.);

    int depel_pos = 0;
    for (int iel = 0; iel < fMesh->NElements(); ++iel) {
        TPZCompEl *cel = fMesh->Element(iel);
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!cel_inter) DebugStop();
        if (cel->Reference()->Dimension() != fMesh->Dimension()) continue;
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());
        TPZMaterialData data;
        cel_inter->InitMaterialData(data);

        TPZMaterial *cel_mat = cel->Material();
        TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem> *mat = dynamic_cast<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem> *>(cel_mat);

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element

        for (int64_t ipts = 0; ipts < npts; ipts++) {
            TPZVec<REAL> qsi(dim);
            REAL w;
            int_rule->Point(ipts, qsi, w);
            cel_inter->ComputeRequiredData(data, qsi);
            REAL weight = w * fabs(data.detjac);

            TPZTensor<REAL> deltastrain;
            TPZTensor<REAL> stress;
            TPZFMatrix<REAL> dep(6,6);
            deltastrain.Zero();
            stress.Zero();
            dep.Zero();
            mat->GetPlasticModel().ApplyStrainComputeSigma(deltastrain,stress,&dep);

            int pos1 = depel_pos + ipts * (dim * dim * npts + dim);
            int pos2 = depel_pos + ipts * (dim * dim * npts + dim) + 1;
            int pos3 = depel_pos + ipts * (dim * dim * npts + dim) + dim * npts;
            int pos4 = depel_pos + ipts * (dim * dim * npts + dim) + dim * npts + 1;

            depxx[pos1] = weight * dep.GetVal(_XX_, _XX_);
            depxx[pos2] = dep.GetVal(_XX_, _XY_) * 0.5;
            depxx[pos3] = dep.GetVal(_XY_, _XX_);
            depxx[pos4] = weight * dep.GetVal(_XY_, _XY_) * 0.5;

            depyy[pos1] = weight * dep.GetVal(_XY_, _XY_) * 0.5;
            depyy[pos2] = dep.GetVal(_XY_, _YY_);
            depyy[pos3] = dep.GetVal(_YY_, _XY_) * 0.5;
            depyy[pos4] = weight * dep.GetVal(_YY_, _YY_);

            depxy[pos1] = dep.GetVal(_XX_, _XY_) * 0.5;
            depxy[pos2] = weight * dep.GetVal(_XX_, _YY_);
            depxy[pos3] = weight * dep.GetVal(_XY_, _XY_) * 0.5;
            depxy[pos4] = dep.GetVal(_XY_, _YY_);
        }
        depel_pos = depel_pos + (npts * dim) * (npts * dim);
    }
}
