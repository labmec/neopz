#include "pzblackoilanalysis.h"


/** Calcula a pressao media nos volumes vizinhos as faces de material matid.
 */
double Pressao(TPZBlackOilAnalysis &an, int matid);

/** Calcula a vazao total que passa pelas faces de material matid
 */
void Vazao(TPZBlackOilAnalysis &an, int matid, double & VazaoAgua, double  & VazaoOleo);
