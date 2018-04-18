#include "sloan.h"

/* Subroutine */ 
int iorele_ (int64_t *ior, int64_t *nnn, int64_t *npn, int64_t *xnpn, int64_t *key, int64_t *, int64_t *nel, int64_t *)
 // ior,nnn,npn,xnpn,key,numno,nel,nop
{
 /* System generated locals */
    int64_t     i__1, i__2;
 /* Local variables */
    static int64_t mini, i, j, noeud;

 /* Parameter adjustments */
    --npn;
    --nnn;
    --key;
    --xnpn;
    --ior;

 /* Function Body */
    i__1 = *nel;
    for (i = 1; i <= i__1; ++i)
    {
        mini = 9999999;
        i__2 = xnpn[i + 1] - 1;
        for (j = xnpn[i]; j <= i__2; ++j)
        {
            noeud = nnn[npn[j]];
            mini = (mini < noeud) ? mini : noeud;
/* L2: */
        }
        key[i] = mini;
/* L1: */
    }

    vsrtp1_ (&key[1], &ior[1], nel);

    return 0;
}                               /* iorele_ */
