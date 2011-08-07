/**
 * @file
 * @brief Contains the implementation of the TPZReadMesh methods. 
 */
/*****************************************************************************
 * O conteudo desse arquivo eh de propriedade do LabMeC-DES-FEC-UNICAMP e do
 * CENPES-Petrobras. 
 * O uso de qualquer parte ou do todo esta condicionado a expressa autorizacao
 * dos proprietarios.
 *****************************************************************************/

#include "pzreadmesh.h"

TPZReadMesh::TPZReadMesh(const char * FileFullPath) :
fInputFile(FileFullPath)
{
}

TPZReadMesh::~TPZReadMesh()
{
}
