#include "pzmetis.h"
#include "pzvec.h"

int main() {

	int NElements = 5;
   int NNodes = 6;
	TPZVec<int> elgraph(NElements*2,-1);
   TPZVec<int> elgraphindex(NNodes+1);
   TPZVec<int> perm, iperm;
   int el,elindex = 0;
   for(el=0; el<NElements; el++) {
      elgraphindex[el] = elindex;
      elgraph[elindex++] = el;
      elgraph[elindex++] = el+1;
	}
   elgraphindex[el] = elindex;
   TPZMetis renum(NElements,NNodes);
   renum.SetElementGraph(elgraph,elgraphindex);
   renum.Resequence(perm,iperm);
   int i;
   for(i=0; i<NNodes; i++) {
   	cout << perm[i] << ' ';
   }
   cout << endl;
   for(i=0; i<NNodes; i++) {
   	cout << iperm[i] << ' ';
   }
   cout << endl;
   return 0;
}
