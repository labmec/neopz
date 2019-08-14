#include "TPZElastoPlasticIntPointsStructMatrix.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "pzmetis.h"
#include "TPZConstitutiveLawProcessor.h"
#include "TPZElasticCriterion.h"
#include "pzbndcond.h"
#include "Timer.h"

#ifdef USING_MKL
#include <mkl.h>
#endif


TPZElastoPlasticIntPointsStructMatrix::TPZElastoPlasticIntPointsStructMatrix(TPZCompMesh *cmesh) : TPZSymetricSpStructMatrix(cmesh), fSparseMatrixLinear(), fRhsLinear(), fIntegrator(), fBCMaterialIds() {

    if (!cmesh->Reference()->Dimension()) {
        DebugStop();
    }
    fDimension = cmesh->Reference()->Dimension();
    m_IA_to_sequence.resize(0);
    m_JA_to_sequence.resize(0);

    #ifdef USING_CUDA
    d_IA_to_sequence.resize(0);    
    d_JA_to_sequence.resize(0);
    d_el_color_indexes.resize(0); 
    #endif
}

TPZElastoPlasticIntPointsStructMatrix::~TPZElastoPlasticIntPointsStructMatrix() {
}

TPZStructMatrix * TPZElastoPlasticIntPointsStructMatrix::Clone(){
    return new TPZElastoPlasticIntPointsStructMatrix(*this);
}

TPZMatrix<STATE> * TPZElastoPlasticIntPointsStructMatrix::Create(){

    if(!isBuilt()) {
        this->SetUpDataStructure(); // When basis functions are computed and storaged
    }
    
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    fMesh->ComputeElGraph(elgraph,elgraphindex,fMaterialIds); // This method seems to be efficient.
    TPZMatrix<STATE> * mat = SetupMatrixData(elgraph, elgraphindex);
    
    /// Sparsify global indexes
    // Filling local std::map
    TPZSYsmpMatrix<STATE> *stiff = dynamic_cast<TPZSYsmpMatrix<STATE> *> (mat);

    m_IA_to_sequence.resize(stiff->IA().size());
    for (int i = 0; i < stiff->IA().size(); ++i) {
        m_IA_to_sequence[i] = stiff->IA()[i];
    }

    m_JA_to_sequence.resize(stiff->JA().size());
    for (int i = 0; i < stiff->JA().size(); ++i) {
        m_JA_to_sequence[i] = stiff->JA()[i];
    }

    FillLIndexes();
    
#ifdef USING_CUDA
    Timer timer;   
    timer.TimeUnit(Timer::ESeconds);
    timer.TimerOption(Timer::ECudaEvent);
    timer.Start();
    std::cout << "Transfering data to GPU..." << std::endl;
    fIntegrator.TransferDataToGPU();
    this->TransferDataToGPU();
    timer.Stop();
    std::cout << "Done! It took " <<  timer.ElapsedTime() << timer.Unit() << std::endl;
#endif

    return mat;
}

void TPZElastoPlasticIntPointsStructMatrix::FillLIndexes(){

    TPZVec<int> & indexes = fIntegrator.DoFIndexes();
    TPZVec<int> & el_n_dofs = fIntegrator.IrregularBlocksMatrix().Blocks().fColSizes;
    TPZVec<int> & cols_first_index = fIntegrator.IrregularBlocksMatrix().Blocks().fColFirstIndex;
    
    int n_colors = m_first_color_index.size()-1;
    m_first_color_l_index.resize(n_colors+1);
    m_first_color_l_index[0]=0;
    for (int ic = 0; ic < n_colors; ic++) {
        
        int first = m_first_color_index[ic];
        int last = m_first_color_index[ic + 1];
        int nel_per_color = last - first;
        int64_t c = 0;
        for (int i = 0; i < nel_per_color; i++) {
            int iel = m_el_color_indexes[first + i];
            int el_dof = el_n_dofs[iel];
            int n_entries = (el_dof*el_dof + el_dof)/2;
            c += n_entries;
        }
        m_first_color_l_index[ic+1] = c + m_first_color_l_index[ic];
    }
    m_color_l_sequence.resize(m_first_color_l_index[n_colors]);
    
    for (int ic = 0; ic < n_colors; ic++) {
        int first = m_first_color_index[ic];
        int last = m_first_color_index[ic + 1];
        int nel_per_color = last - first;
        int64_t c = m_first_color_l_index[ic];
        for (int i = 0; i < nel_per_color; i++) {
            int iel = m_el_color_indexes[first + i];
            int el_dof = el_n_dofs[iel];
            int pos = cols_first_index[iel];
            for (int i_dof = 0; i_dof < el_dof; i_dof++) {
                int64_t i_dest = indexes[pos + i_dof];
                for (int j_dof = i_dof; j_dof < el_dof; j_dof++) {
                    int64_t j_dest = indexes[pos + j_dof];
                    int64_t l_index = me(m_IA_to_sequence, m_JA_to_sequence, i_dest, j_dest);
                    m_color_l_sequence[c] = l_index;
                    c++;
                }
            }
        }
    }
    
    
}

#ifdef USING_CUDA
void TPZElastoPlasticIntPointsStructMatrix::TransferDataToGPU() {
    d_IA_to_sequence.resize(m_IA_to_sequence.size());
    d_IA_to_sequence.set(&m_IA_to_sequence[0], m_IA_to_sequence.size());

    d_JA_to_sequence.resize(m_JA_to_sequence.size());
    d_JA_to_sequence.set(&m_JA_to_sequence[0], m_JA_to_sequence.size());

    d_el_color_indexes.resize(m_el_color_indexes.size());
    d_el_color_indexes.set(&m_el_color_indexes[0], m_el_color_indexes.size());

    d_RhsLinear.resize(fRhsLinear.Rows());
    d_RhsLinear.set(&fRhsLinear(0,0), fRhsLinear.Rows());

    d_color_l_sequence.resize(m_color_l_sequence.size());
    d_color_l_sequence.set(&m_color_l_sequence[0], m_color_l_sequence.size());

}
#endif

int64_t TPZElastoPlasticIntPointsStructMatrix::me(TPZVec<int> &IA, TPZVec<int> &JA, int64_t & i_dest, int64_t & j_dest) {
    int64_t row(i_dest),col(j_dest);
    if (i_dest > j_dest) {
        int64_t temp = i_dest;
        row = col;
        col = temp;
    }
    for(int ic=IA[row] ; ic < IA[row+1]; ic++ ) {
        if ( JA[ic] == col )
        {
            return ic;
        }
    }
    return 0;
}

TPZMatrix<STATE> *TPZElastoPlasticIntPointsStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {

    int64_t neq = fMesh->NEquations();
    TPZMatrix<STATE> *stiff = Create(); // @TODO:: Requires optimization.
    rhs.Redim(neq,1);
    Assemble(*stiff,rhs,guiInterface);
    return stiff;
}

void TPZElastoPlasticIntPointsStructMatrix::SetUpDataStructure() {

    if(isBuilt()) {
        std::cout << __PRETTY_FUNCTION__ << " Data structure has been setup." << std::endl;
        return;
    }
    
    TPZVec<int> element_indexes;
    ComputeDomainElementIndexes(element_indexes);
    
    ClassifyMaterialsByDimension();

    TPZIrregularBlocksMatrix::IrregularBlocks blocksData;
    this->SetUpIrregularBlocksData(element_indexes, blocksData);

    int64_t rows = blocksData.fRowFirstIndex[blocksData.fNumBlocks];
    int64_t cols = blocksData.fColFirstIndex[blocksData.fNumBlocks];
    TPZIrregularBlocksMatrix blocksMatrix(rows, cols);
    blocksMatrix.SetBlocks(blocksData);
    fIntegrator.SetIrregularBlocksMatrix(blocksMatrix);

    TPZVec<int> dof_indexes;
    this->SetUpIndexes(element_indexes, dof_indexes);
    fIntegrator.SetDoFIndexes(dof_indexes);

    TPZVec<int> colored_element_indexes;
    int ncolor;
    
    this->ColoredIndexes(element_indexes, dof_indexes, colored_element_indexes, ncolor);
    
    fIntegrator.SetColorIndexes(colored_element_indexes);
    fIntegrator.SetNColors(ncolor);

    AssembleBoundaryData();
}

 void TPZElastoPlasticIntPointsStructMatrix::Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    TPZFMatrix<REAL> Dep;
     fIntegrator.ResidualIntegration(fMesh->Solution(),rhs, Dep);
     rhs += fRhsLinear;
//     Dep.Print(std::cout);

    TPZSYsmpMatrix<STATE> &stiff = dynamic_cast<TPZSYsmpMatrix<STATE> &> (mat);
    TPZVec<STATE> &Kg = stiff.A();
    int64_t nnz = Kg.size();
    
    TPZVec<int> & indexes = fIntegrator.DoFIndexes();
    TPZVec<int> & el_n_dofs = fIntegrator.IrregularBlocksMatrix().Blocks().fColSizes;
    TPZVec<int> & cols_first_index = fIntegrator.IrregularBlocksMatrix().Blocks().fColFirstIndex;

    int n_colors = m_first_color_index.size()-1;

#ifdef USING_CUDA      
    TPZVecGPU<REAL> d_Kg(nnz);
    d_Kg.Zero();
    fIntegrator.ConstitutiveLawProcessor().De();
#endif

    Timer timer;   
    timer.TimeUnit(Timer::ESeconds);
    timer.TimerOption(Timer::EChrono);
    timer.Start();
#ifdef USING_CUDA
    #ifdef COMPUTE_K_GS
    std::cout << "COMPUTE_K_GPU_GS" << std::endl;          
        int first = m_first_color_index[0];
        int last = m_first_color_index[n_colors];
        int el_dofs = el_n_dofs[0];
        int nel_per_color = last - first;

        TPZVecGPU<REAL> d_Kc(m_color_l_sequence.size());
        d_Kc.Zero();
        fCudaCalls.MatrixAssembleGS(d_Kc.getData(), first, last, d_el_color_indexes.getData(), fIntegrator.ConstitutiveLawProcessor().WeightVectorDev().getData(),
            fIntegrator.DoFIndexesDev().getData(), fIntegrator.IrregularBlocksMatrix().BlocksDev().dStorage.getData(),
            fIntegrator.IrregularBlocksMatrix().BlocksDev().dRowSizes.getData(), fIntegrator.IrregularBlocksMatrix().BlocksDev().dColSizes.getData(),
            fIntegrator.IrregularBlocksMatrix().BlocksDev().dRowFirstIndex.getData(), fIntegrator.IrregularBlocksMatrix().BlocksDev().dColFirstIndex.getData(),
            fIntegrator.IrregularBlocksMatrix().BlocksDev().dMatrixPosition.getData());

    for (int ic = 0; ic < n_colors; ic++) {
        int first_l = m_first_color_l_index[ic];
        int last_l = m_first_color_l_index[ic + 1];
        int n_l_indexes = last_l - first_l;
        TPZVecGPU<REAL> aux(n_l_indexes);
        fCudaCalls.GatherOperation(n_l_indexes, d_Kg.getData(), aux.getData(), &d_color_l_sequence.getData()[first_l]);
        fCudaCalls.DaxpyOperation(n_l_indexes, 1., &d_Kc.getData()[first_l], aux.getData());
        fCudaCalls.ScatterOperation(n_l_indexes, aux.getData(), d_Kg.getData(), &d_color_l_sequence.getData()[first_l]);
    }
    d_Kg.get(&Kg[0], d_Kg.getSize()); // back to CPU
    #else
    std::cout << "COMPUTE_K_GPU" << std::endl;
    for (int ic = 0; ic < n_colors; ic++) /// Serial by color
    {     
        int first = m_first_color_index[ic];
        int last = m_first_color_index[ic + 1];
        fCudaCalls.MatrixAssemble(d_Kg.getData(), first, last, &d_el_color_indexes.getData()[first], fIntegrator.ConstitutiveLawProcessor().WeightVectorDev().getData(),
            fIntegrator.DoFIndexesDev().getData(), fIntegrator.IrregularBlocksMatrix().BlocksDev().dStorage.getData(),
            fIntegrator.IrregularBlocksMatrix().BlocksDev().dRowSizes.getData(), fIntegrator.IrregularBlocksMatrix().BlocksDev().dColSizes.getData(),
            fIntegrator.IrregularBlocksMatrix().BlocksDev().dRowFirstIndex.getData(), fIntegrator.IrregularBlocksMatrix().BlocksDev().dColFirstIndex.getData(),
            fIntegrator.IrregularBlocksMatrix().BlocksDev().dMatrixPosition.getData(), d_IA_to_sequence.getData(), d_JA_to_sequence.getData());
    }
    d_Kg.get(&Kg[0], d_Kg.getSize()); // back to CPU    
    #endif
#else
#ifdef COMPUTE_K_GS
     std::cout << "COMPUTE_K_CPU_GS" << std::endl;
     
     int first = m_first_color_index[0];
     int last = m_first_color_index[n_colors];
     int el_dofs = el_n_dofs[0];
     int nel_per_color = last - first;
     TPZVec<REAL> Kc(m_color_l_sequence.size(),0.0);
     
#ifdef USING_TBB
         tbb::parallel_for(size_t(0),size_t(nel_per_color),size_t(1),[&](size_t i)
#else
        for (int i = 0; i < nel_per_color; i++)
#endif
         {
             int iel = m_el_color_indexes[first + i];
             
             // Compute Elementary Matrix.
             TPZFMatrix<STATE> K;
             fIntegrator.ComputeTangentMatrix(iel,K);
             int stride = i*(el_dofs * el_dofs + el_dofs)/2;
             int c = stride;
             for(int i_dof = 0 ; i_dof < el_dofs; i_dof++){
                 for(int j_dof = i_dof; j_dof < el_dofs; j_dof++){
                     Kc[c] += K(i_dof,j_dof);
                     c++;
                 }
             }
         }
#ifdef USING_TBB
);
#endif
         
       for (int ic = 0; ic < n_colors; ic++) {
           int first_l = m_first_color_l_index[ic];
           int last_l = m_first_color_l_index[ic + 1];
           int n_l_indexes = last_l - first_l;
           /// Gather from Kg
           TPZVec<REAL> aux(n_l_indexes);
           cblas_dgthr(n_l_indexes, &Kg[0], &aux[0], &m_color_l_sequence[first_l]);
           
           // Contributing
           cblas_daxpy(n_l_indexes, 1., &Kc[first_l], 1., &aux[0],1);
           
           /// Scatter to Kg
           cblas_dsctr(n_l_indexes, &aux[0], &m_color_l_sequence[first_l], &Kg[0]);
       }

#else
     std::cout << "COMPUTE_K_CPU" << std::endl;
     
      for (int ic = 0; ic < n_colors; ic++) {
             
          
#ifdef USING_TBB
       tbb::parallel_for(size_t(m_first_color_index[ic]),size_t(m_first_color_index[ic+1]),size_t(1),[&](size_t i)
                                                 
#else
       for (int i = m_first_color_index[ic]; i < m_first_color_index[ic+1]; i++)
#endif
       {
             int iel = m_el_color_indexes[i];
           
             // Compute Elementary Matrix.
             TPZFMatrix<STATE> K;
//             fIntegrator.ComputeTangentMatrix(iel,K);
             fIntegrator.ComputeTangentMatrix(iel,Dep,K);

             int el_dof = el_n_dofs[iel];
             int pos = cols_first_index[iel];
           
             for (int i_dof = 0; i_dof < el_dof; i_dof++) {
                 
                 int64_t i_dest = indexes[pos + i_dof];
                 
                 for (int j_dof = i_dof; j_dof < el_dof; j_dof++) {
                     
                     int64_t j_dest = indexes[pos + j_dof];
                     STATE val = K(i_dof,j_dof);
                     int64_t  index = me(m_IA_to_sequence, m_JA_to_sequence, i_dest, j_dest);
                     Kg[index] += val;
                     
                 }
             }
           
         }
#ifdef USING_TBB
 );
#endif
          
}
         
#endif

#endif

    timer.Stop();
    std::cout << "K Assemble: Elasped time [sec] = " << timer.ElapsedTime() << std::endl;
     
    auto it_end = fSparseMatrixLinear.MapEnd();
    for (auto it = fSparseMatrixLinear.MapBegin(); it!=it_end; it++) {
      int64_t row = it->first.first;
      int64_t col = it->first.second;
      STATE val = it->second + stiff.GetVal(row, col);
      stiff.PutVal(row, col, val);
    }

//    timer.Start();
//    Assemble(rhs,guiInterface);
//    timer.Stop();
//    std::cout << "R Assemble: Elasped time [sec] = " << timer.ElapsedTime() << std::endl;
}

int TPZElastoPlasticIntPointsStructMatrix::StressRateVectorSize(){

    switch (fDimension) {
        case 1:{
            return 1;
        }
        break;
        case 2:{
            return 3;
        }
        break;
        case 3:{
            return 6;
        }
        break;
        default:
        {
            DebugStop();
            return 0;
        }
        break;
    }
}

void TPZElastoPlasticIntPointsStructMatrix::Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    int neq = fMesh->NEquations();

    rhs.Resize(neq, 1);
    rhs.Zero();  
#ifdef USING_CUDA
    TPZVecGPU<REAL> d_rhs(rhs.Rows());
    d_rhs.Zero();
    fIntegrator.ResidualIntegration(fMesh->Solution(),d_rhs);
    fCudaCalls.DaxpyOperation(neq, 1., d_RhsLinear.getData(), d_rhs.getData());
    d_rhs.get(&rhs(0,0), neq); //back to CPU
#else
    fIntegrator.ResidualIntegration(fMesh->Solution(),rhs);
    rhs += fRhsLinear;
#endif

}

void TPZElastoPlasticIntPointsStructMatrix::AssembleBoundaryData() {

    int64_t neq = fMesh->NEquations();
    TPZStructMatrix str(fMesh);
    str.SetMaterialIds(fBCMaterialIds);
    str.SetNumThreads(fNumThreads);
    TPZAutoPointer<TPZGuiInterface> guiInterface;
    fRhsLinear.Resize(neq, 1);
    fRhsLinear.Zero();
    fSparseMatrixLinear.Resize(neq, neq);
    str.Assemble(fSparseMatrixLinear, fRhsLinear, guiInterface);
}

void TPZElastoPlasticIntPointsStructMatrix::ComputeDomainElementIndexes(TPZVec<int> &element_indexes) {

    TPZStack<int> el_indexes_loc;
    for (int64_t i = 0; i < fMesh->NElements(); i++) {
        TPZCompEl *cel = fMesh->Element(i);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;
        if(gel->Dimension() == fDimension){
            el_indexes_loc.Push(cel->Index());
        }
    }
    element_indexes = el_indexes_loc;
}

void TPZElastoPlasticIntPointsStructMatrix::ClassifyMaterialsByDimension() {

    for (auto material : fMesh->MaterialVec()) {
        TPZBndCond * bc_mat = dynamic_cast<TPZBndCond *>(material.second);
        bool domain_material_Q = !bc_mat;
        if (domain_material_Q) {
            fMaterialIds.insert(material.first);
        }else{
            fBCMaterialIds.insert(material.first);
        }
    }
}

void TPZElastoPlasticIntPointsStructMatrix::SetUpIrregularBlocksData(TPZVec<int> &element_indexes, TPZIrregularBlocksMatrix::IrregularBlocks &blocksData) {

    /// Number of elements
    int nblocks = element_indexes.size();

    blocksData.fNumBlocks = nblocks;
    blocksData.fRowSizes.resize(nblocks);
    blocksData.fColSizes.resize(nblocks);
    blocksData.fMatrixPosition.resize(nblocks + 1);
    blocksData.fRowFirstIndex.resize(nblocks + 1);
    blocksData.fColFirstIndex.resize(nblocks + 1);

    blocksData.fMatrixPosition[0] = 0;
    blocksData.fRowFirstIndex[0] = 0;
    blocksData.fColFirstIndex[0] = 0;

    // @TODO Candidate for TBB ParallelScan
    // Example with lambda expression https://www.threadingbuildingblocks.org/docs/help/reference/algorithms/parallel_scan_func.html
    for(int iel = 0; iel < nblocks; iel++) {

        TPZCompEl *cel = fMesh->Element(element_indexes[iel]);
        
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!cel_inter) DebugStop();
#ifdef PZDEBUG
        int dim = cel->Reference()->Dimension();
        if (dim !=fDimension) {
            DebugStop();
        }
#endif
        
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());
        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t n_el_dof = cel_inter->NShapeF(); // number of shape functions of the element (Element DoF)

        blocksData.fRowSizes[iel] = StressRateVectorSize() * npts;
        blocksData.fColSizes[iel] = n_el_dof * fDimension;
        
        blocksData.fMatrixPosition[iel + 1] = blocksData.fMatrixPosition[iel] + blocksData.fRowSizes[iel] * blocksData.fColSizes[iel];
        blocksData.fRowFirstIndex[iel + 1] =  blocksData.fRowFirstIndex[iel] + blocksData.fRowSizes[iel];
        blocksData.fColFirstIndex[iel + 1] = blocksData.fColFirstIndex[iel] + blocksData.fColSizes[iel];
    }

    blocksData.fStorage.resize(blocksData.fMatrixPosition[nblocks]);
    
    
    
    TPZFNMatrix<8,REAL> dphiXY;
    TPZMaterialData data;
    TPZVec<REAL> qsi(fDimension);
    
    // @TODO Candidate for TBB ParallelFor
    for (int iel = 0; iel < nblocks; ++iel) {

        TPZCompEl *cel = fMesh->Element(element_indexes[iel]);
        
        int sigma_entries_times_npts      = blocksData.fRowSizes[iel];
        int n_el_dof            = blocksData.fColSizes[iel]/fDimension;
        int pos_el              = blocksData.fMatrixPosition[iel];
        
        
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!cel_inter) DebugStop();
        
#ifdef PZDEBUG
        int dim = cel->Reference()->Dimension();
        if (dim !=fDimension) {
            DebugStop();
        }
#endif
        
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());
        
        int npts = int_rule->NPoints();
        int n_sigma_entries = StressRateVectorSize();
        
#ifdef PZDEBUG
        if (npts != sigma_entries_times_npts / n_sigma_entries) {
            DebugStop();
        }
#endif

        cel_inter->InitMaterialData(data);
        TPZFMatrix<REAL> B_el(npts*n_sigma_entries,n_el_dof * fDimension,0.0);
        for (int64_t i_pts = 0; i_pts < npts; i_pts++) {

            REAL w;
            int_rule->Point(i_pts, qsi, w);
            cel_inter->ComputeRequiredData(data, qsi);

            data.axes.Transpose();
            data.axes.Multiply(data.dphix, dphiXY);

            for (int j_dim = 0; j_dim < 2; j_dim++){

                for (int i_dof = 0; i_dof < n_el_dof; i_dof++) {

                    REAL val = dphiXY(j_dim, i_dof);
                    for (int i_sigma_comp = 0; i_sigma_comp < fDimension; i_sigma_comp++){

                        B_el(i_pts * n_sigma_entries + i_sigma_comp + j_dim, i_dof*2 + i_sigma_comp) = val;
                    }
                }
            }

        }
        B_el.Transpose();
        TPZFMatrix<REAL> B_el_loc(npts * n_sigma_entries, n_el_dof * fDimension, &blocksData.fStorage[pos_el], npts * n_sigma_entries * n_el_dof * fDimension);
        B_el_loc = B_el;
    }
}

void TPZElastoPlasticIntPointsStructMatrix::SetUpIndexes(TPZVec<int> &element_indexes, TPZVec<int> & dof_indexes) {

    int64_t nblocks = fIntegrator.IrregularBlocksMatrix().Blocks().fNumBlocks;
    int64_t rows = fIntegrator.IrregularBlocksMatrix().Rows();
    int64_t cols = fIntegrator.IrregularBlocksMatrix().Cols();
    TPZVec<int> & dof_positions = fIntegrator.IrregularBlocksMatrix().Blocks().fColFirstIndex;

    dof_indexes.resize(cols);
    int64_t npts = rows / StressRateVectorSize();
    TPZVec<REAL> weight(npts);

    int64_t wit = 0;
    for (int iel = 0; iel < nblocks; ++iel) {

        int dof_pos = dof_positions[iel];
        TPZCompEl *cel = fMesh->Element(element_indexes[iel]);
        
        
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!cel_inter) DebugStop();
        TPZGeoEl * gel = cel->Reference();
        if (!gel) DebugStop();
        
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t el_npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element

        TPZFMatrix<REAL> jac,axes, jacinv;
        REAL detjac;
        for (int64_t ipts = 0; ipts < el_npts; ipts++) {
            TPZVec<REAL> qsi(dim);
            REAL w;
            int_rule->Point(ipts, qsi, w);
            gel->Jacobian(qsi, jac, axes, detjac, jacinv);
            weight[wit] = w * std::abs(detjac);
            wit++;
        }

        int64_t ncon = cel->NConnects();
        int i_dof = 0;
        for (int64_t icon = 0; icon < ncon; icon++) {
            int64_t id = cel->ConnectIndex(icon);
            TPZConnect &df = fMesh->ConnectVec()[id];
            int64_t conid = df.SequenceNumber();
            if (df.NElConnected() == 0 || conid < 0 || fMesh->Block().Size(conid) == 0) continue;
            int64_t pos = fMesh->Block().Position(conid);
            int64_t nsize = fMesh->Block().Size(conid);
            for (int64_t isize = 0; isize < nsize; isize++) {
                dof_indexes[dof_pos+i_dof] = pos + isize;
                i_dof++;
            }
            
        }
    }

    TPZMaterial *material = fMesh->FindMaterial(1);
    fIntegrator.ConstitutiveLawProcessor().SetMaterial(material);
    fIntegrator.ConstitutiveLawProcessor().SetUpDataByIntPoints(npts);
    fIntegrator.ConstitutiveLawProcessor().SetWeightVector(weight);
}

void TPZElastoPlasticIntPointsStructMatrix::ColoredIndexes(TPZVec<int> &element_indexes, TPZVec<int> &indexes, TPZVec<int> &coloredindexes, int &ncolor) {

    int64_t nblocks = fIntegrator.IrregularBlocksMatrix().Blocks().fNumBlocks;
    int64_t cols = fIntegrator.IrregularBlocksMatrix().Cols();

    TPZVec<int64_t> connects_vec(fMesh->NConnects(),0);
    TPZVec<int64_t> elemcolor(nblocks,-1);

    int64_t contcolor = 0;
    bool needstocontinue = true;

    while (needstocontinue)
    {
        int it = 0;
        needstocontinue = false;
        for (auto iel : element_indexes) {
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
    coloredindexes.resize(cols);
    int64_t neq = fMesh->NEquations();
    for (int64_t iel = 0; iel < nblocks; iel++) {
        int64_t elem_col = fIntegrator.IrregularBlocksMatrix().Blocks().fColSizes[iel];
        int64_t cont_cols = fIntegrator.IrregularBlocksMatrix().Blocks().fColFirstIndex[iel];

        for (int64_t icols = 0; icols < elem_col; icols++) {
            coloredindexes[cont_cols + icols] = indexes[cont_cols + icols] + elemcolor[iel]*neq;
        }
    }
    
    std::map<int64_t, std::vector<int64_t> > color_map;
    for (int64_t iel = 0; iel < nblocks; iel++) {
        color_map[elemcolor[iel]].push_back(iel);
    }
    
    if (contcolor != color_map.size()) {
        DebugStop();
    }
    
    m_el_color_indexes.resize(nblocks);
    m_first_color_index.resize(color_map.size()+1);

    int c_color = 0;

    m_first_color_index[c_color] = 0;
    for (auto color_data : color_map) {
        int n_el_per_color = color_data.second.size();
        int iel = m_first_color_index[c_color];
        for (int i = 0; i < n_el_per_color ; i++) {
            int el_index = color_data.second[i];
            m_el_color_indexes[iel] = el_index;
            iel++;
        }
        c_color++;
        m_first_color_index[c_color] = n_el_per_color + m_first_color_index[c_color-1];
    }
    
    std::cout << "Number of colors = " << color_map.size() << std::endl;
}
