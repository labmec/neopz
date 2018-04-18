
#include "tpzagglomeratemesh.h"

int TPZAgglomerateMesh::ClassId() const{
    return Hash("TPZAgglomerateMesh") ^ TPZFlowCompMesh::ClassId() << 1;
}