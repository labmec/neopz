// I DON'T NEED THIS ROUTINE-ARIEL!
#include "sloan.h"

 /* @brief Purpose: Compute the profiles using both original and new numbers */
int profi1_ (int64_t *n, int64_t *nnn, int64_t *, int64_t *adj, int64_t *xadj, int64_t *oldpro, int64_t *newpro)
//n, nnn, e2, adj, xadj, oldpro, newpro
{
 /* System generated locals */
    int64_t     i__1, i__2, i__3, i__4;

 /* Local variables */
    static int64_t i, j, jstop, jstrt, oldmin, newmin;


/*     INPUT: */
/*     ------ */
/*     N      - Number of nodes in graph */
/*     NNN    - List of new node numbers for graph */
/*            - New node number for node I is given by NNN(I) */
/*     E2     - Twice the number of edges in the graph = XADJ(N+1)-1 */
/*     ADJ    - Adjacency list for all nodes in graph */
/*            - List of length 2E where E is the number of edges in */
/*              the graph and 2E = XADJ(N+1)-1 */
/*     XADJ   - Index number for ADJ */
/*            - Nodes adjacent to node I are found in ADJ(I), where */
/*              J=XADJ(I), XADJ(I)+1, ..., XADJ(I+1)-1 */
/*     OLDPRO - Undefined */
/*     NEWPRO - Undefined */

/*     OUTPUT: */
/*     ------- */
/*     N      - Unchanged */
/*     NNN    - Unchanged */
/*     E2     - Unchanged */
/*     ADJ    - Unchanged */
/*     XADJ   - Unchanged */
/*     OLDPRO - Profile with original node numbering */
/*     NEWPRO - Profile with new node numbering */

/*     NOTE:    Profiles include diagonal terms */
/*     ----- */

/*     PROGRAMMER:   Scott Sloan */
/*     ----------- */

/*     LAST MODIFIED:   10 March 1989        Scott Sloan */
/*     -------------- */
/* *********************************************************************** */


/*     Set profiles and loop over each node in graph */

 /* Parameter adjustments */
    --xadj;
    --nnn;
    --adj;

 /* Function Body */
    if (*oldpro)
        *oldpro = 0;
    if (*newpro)
        *newpro = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        jstrt = xadj[i];
        jstop = xadj[i + 1] - 1;
        oldmin = adj[jstrt];
        newmin = nnn[adj[jstrt]];

/*        Find lowest numbered neighbour of node I */
/*        (using both old and new node numbers) */

        i__2 = jstop;
        for (j = jstrt + 1; j <= i__2; ++j)
        {
/* Computing MIN */
            i__3 = oldmin, i__4 = adj[j];
            oldmin = (i__3 < i__4) ? i__3 : i__4;
            //oldmin = min(i__3,i__4);
/* Computing MIN */
            i__3 = newmin, i__4 = nnn[adj[j]];
            newmin = (i__3 < i__4) ? i__3 : i__4;
/* L10: */
        }

/*        Update profiles */

      *oldpro += i_dim(i, oldmin);
      *newpro += i_dim(nnn[i], newmin);

/* L20: */
    }

/*     Add diagonal terms to profiles */

    *oldpro += *n;
    *newpro += *n;
    return 0;
}                               /* profi1_ */

int64_t i_dim(int64_t value_1, int64_t value_2)
{
	if(value_1<=value_2)
	{
		return value_1;
	}
	else
	{
		return value_2;
	}
}
