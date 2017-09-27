
#include "tpzagglomeratemesh.h"

int TPZAgglomerateMesh::ClassId(){
    return TPZFlowCompMesh::ClassId() ^ Hash("TPZAgglomerateMesh");
}