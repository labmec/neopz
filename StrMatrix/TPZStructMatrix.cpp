#include "TPZStructMatrix.h"

#include "pzcmesh.h"
#include "pzcheckmesh.h"
#include "pzerror.h"
#include "pzsubcmesh.h"
#include "TPZLinearAnalysis.h"
#include "pzshtmat.h"
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.TPZStructMatrixOR");
static TPZLogger loggerel("pz.strmatrix.element");
static TPZLogger loggerel2("pz.strmatrix.elementinterface");
static TPZLogger loggerelmat("pz.strmatrix.elementmat");
static TPZLogger loggerCheck("pz.strmatrix.checkconsistency");
static TPZLogger loggerGlobStiff("pz.strmatrix.globalstiffness");
#else
static int logger;
#endif

//this is not in header file to avoid including cmesh.h there
TPZStructMatrix::~TPZStructMatrix(){}

TPZStructMatrix::TPZStructMatrix(TPZCompMesh *mesh)
    : TPZStrMatParInterface(), fEquationFilter(0){
    SetMesh(mesh);
}

TPZStructMatrix::TPZStructMatrix(TPZAutoPointer<TPZCompMesh> mesh)
    : TPZStrMatParInterface(), fEquationFilter(0) {
    SetMesh(mesh);
}

TPZStructMatrix::TPZStructMatrix(const TPZStructMatrix &copy)
    : TPZStrMatParInterface(copy),fMesh(copy.fMesh),
      fCompMesh(copy.fCompMesh),
      fMaterialIds(copy.fMaterialIds), fEquationFilter(copy.fEquationFilter) {
}

void TPZStructMatrix::SetMesh(TPZCompMesh *mesh) {
    fMesh = mesh;
    fEquationFilter.SetNumEq(mesh ? mesh->NEquations() : 0);
#ifdef PZDEBUG
    if (mesh) {
        TPZCheckMesh checkmesh(fMesh, &std::cout);
        if (checkmesh.CheckConnectSeqNumberConsistency() != 0) {
            DebugStop();
        }
    }
#endif
}

void TPZStructMatrix::SetMesh(TPZAutoPointer<TPZCompMesh> mesh) {
    fCompMesh = mesh;
    SetMesh(mesh.operator->());
}

void TPZStructMatrix::SetMaterialIds(const std::set<int> &materialids)
{
    fMaterialIds = materialids;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::set<int>::const_iterator it;
        std::stringstream sout;
        sout << "setting input material ids ";
        for(it=materialids.begin(); it!= materialids.end(); it++)
        {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    if(!fMesh)
    {
        LOGPZ_ERROR(logger,"SetMaterialIds called without mesh")
        return;
    }
    int64_t iel;
    TPZAdmChunkVector<TPZCompEl*> &elvec = fMesh->ElementVec();
    int64_t nel = elvec.NElements();
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = elvec[iel];
        if(!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (cel);
        if(!subcmesh) continue;
        TPZAutoPointer<TPZLinearAnalysis> analysis = subcmesh->Analysis();
        if(!analysis)
        {
            LOGPZ_ERROR(logger,"SetMaterialIds called for substructure without analysis object")
            DebugStop();
        }
        TPZStructMatrix *str = analysis->StructMatrix().operator->();
        if(!str)
        {
            LOGPZ_ERROR(logger,"SetMaterialIds called for substructure without structural matrix")
            continue;
        }
        str->SetMaterialIds(materialids);
    }
}

int TPZStructMatrix::ComputeElementColors(TPZVec<int> &elementcolors)
{
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    if(!fMaterialIds.size())
    {
        fMesh->ComputeElGraph(elgraph, elgraphindex);
    }
    else
    {
        fMesh->ComputeElGraph(elgraph, elgraphindex, fMaterialIds);
    }
    int64_t nconnects = fMesh->NConnects();
    const int color_stored = 10;
    TPZGenMatrix<int64_t> domain(nconnects,color_stored);
    int64_t nel = elgraphindex.size()-1;
    elementcolors.resize(nel);
    elementcolors.Fill(-1);
    bool element_not_colored = true;
    int highest_color = 0;
    int min_color = 0;
    int max_color = min_color+color_stored;
    while(element_not_colored)
    {
        element_not_colored = false;
        for (int64_t el = 0; el<nel; el++) {
            int icolor;
            bool color_found = true;
            for(int icolor = 0; icolor<color_stored; icolor++)
            {
                color_found = true;
                int64_t firstnode = elgraphindex[el];
                int64_t lastnode = elgraphindex[el+1];
                for (int inod = firstnode; inod < lastnode; inod++) {
                    int64_t ic = elgraph[inod];
                    if(domain(ic,icolor) != 0)
                    {
                        color_found = false;
                        break;
                    }
                    domain(ic,icolor) = 1;
                }
                if(color_found == true)
                {
                    elementcolors[el] = min_color+icolor;
                    if(min_color+icolor > highest_color) highest_color = min_color+icolor;
                    break;
                }
            }
            // the element didn't fit any color, the loop will have to be executed again
            if(!color_found) element_not_colored = true;
        }
        min_color += color_stored;
        max_color += color_stored;
        domain.Fill(0);
    }
    return highest_color;
}

int TPZStructMatrix::ClassId() const{
    return Hash("TPZStructMatrix") ^
        TPZStrMatParInterface::ClassId() << 1;
}

void TPZStructMatrix::Read(TPZStream& buf, void* context) {
    //DO NOT READ FROM THE TPZMATPARINTERFACE
    fMesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::GetInstance(&buf));
    fCompMesh = TPZAutoPointerDynamicCast<TPZCompMesh>(TPZPersistenceManager::GetAutoPointer(&buf));
    buf.Read(fMaterialIds);
    fEquationFilter.Read(buf,context);
}

void TPZStructMatrix::Write(TPZStream& buf, int withclassid) const {
    //DO NOT WRITE THE TPZMATPARINTERFACE
    TPZPersistenceManager::WritePointer(fMesh, &buf);
    TPZPersistenceManager::WritePointer(fCompMesh.operator ->(), &buf);
    buf.Write(fMaterialIds);
    fEquationFilter.Write(buf,withclassid);
}



