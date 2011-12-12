#include <iostream>
using namespace std;
#include <stdlib.h> 
#include "sloan.h"


 /* Subroutine */ int 
gegra_ (int *n, int *ne, int *, int *
        npn, int *xnpn, int *iadj, int *adj, int *xadj,
        int *nop)
//n, ne, inpn, npn, xnpn, iadj, adj, xadj, nop
{
 /* Format strings */
    const char   *fmt_10900 = "(\0020%%%E01-GEGRA \002,//,1x,\002CANNOT ASSE\
MBLE NODE ADJACENCY LIST\002,//,1x,\002CHECK NPN AND XNPN ARRAYS\002)";
 /* System generated locals */
    int     i__1, i__2, i__3, i__4;
 /* Local variables */
    int i, j, k, l, m, nodej, nodek, jstop, lstop, mstop, jstrt, jrtst, lstrt, mstrt, nen1;


/*     PURPOSE: */
/*     -------- */

/*     Form adjacency list for a graph corresponding to a finite element
*/
/*     mesh */

/*     INPUT: */
/*     ------ */

/*     N    - Number of nodes in graph (finite element mesh) */
/*     NE   - Number of elements in finte element mesh */
/*     INPN - Length of NPN = XNPN(NE+1)-1 */
/*     NPN  - List of node numbers for each element */
/*     XNPN - Index vector for npn */
/*          - nodes for element I are found in NPN(J), where */
/*            J=XNPN(I), XNPN(I+1), ..., XNPN(I+1)-1 */
/*     IADJ - Length of vector ADJ */
/*          - Set IADJ=NE*NEN*(NEN-1) for a mesh of simple type of */
/*          - element with NEN nodes */
/*          - IADJ=(NEN(1)*(NEN(1)-1)+,.....,+NEN(NE)*(NEN(NE)-1)) */
/*          - for a mesh of elements with varying numbers of nodes */
/*     ADJ  - Undefined */
/*     XADJ - Undefined */

/*     OUTPUT : */
/*     -------- */

/*     N    - Unchanged */
/*     NE   - Unchanged */
/*     INPN - Unchanged */
/*     NPN  - Unchanged */
/*     XNPN - Unchanged */
/*     IADJ - Modified = XADJ(N+1)-1 longueur reelle */
/*     ADJ  - Adjacency list for all nodes in graph */
/*          - List of length 2E where E is the number of edges in */
/*            the graph (note that 2E = XADJ(N+1)-1) */
/*     XADJ - Index vector for ADJ */
/*          - Nodes adjacent to node I are found in ADJ(J), where */
/*          - J = XADJ(I),XADJ(I)+1, ..., XADJ(I+1)-1 */
/*          - Degree of node I give by XADJ(I+1)-XADJ(I) */

/*     NOTES: */
/*     ------ */

/*     This routine typically requires about 25 percent elbow room for */
/*     assembling the ADJ list (i.e. IADJ/2E is typically around 1.25). */
/*     In some cases, the elbow room may be larger (IADJ/2E is slightly */
/*     less than 2 for the 3-nodes triangle) and in order cases it may be
*/
/*     zero (IADJ/2E = 1 for bar elements) */


/*     PROGRAMMER:          Scott Sloan */
/*     ----------- */

/*     LAST MODIFIED:       10 March 1989         Scott Sloan */


/* ***********************************************************************
 */


/*     Initialise the adjacency list and its index vector */

 /* Parameter adjustments */
    
	
	--xadj;
    --xnpn;
    --npn;
    --adj;
	
 /* Function Body */
    i__1 = *iadj;
    for (i = 1; i <= i__1; ++i)
    {
        adj[i] = 0;
/* L5: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        xadj[i] = 0;
/* L10: */
    }

/*     Estimate the degree of each node (always an overestimate) */

    i__1 = *ne;
    for (i = 1; i <= i__1; ++i)
    {
        jstrt = xnpn[i];
        jstop = xnpn[i + 1] - 1;
        nen1 = jstop - jstrt;
        i__2 = jstop;
        for (j = jstrt; j <= i__2; ++j)
        {
            nodej = npn[j];
            xadj[nodej] += nen1;
/* L20: */
        }
/* L30: */
    }

/*     Reconstruct XADJ to point to start of each set of neighbours */


    l=1;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        l += xadj[i];
        xadj[i] = l - xadj[i];
/* L40: */
    }
    xadj[*n + 1] = l;

/*     Form adjency list (which may contain zeros) */

    i__1 = *ne;
    for (i = 1; i <= i__1; ++i)
    {
        jrtst = xnpn[i];
        jstop = xnpn[i + 1] - 1;
        i__2 = jstop - 1;
        for (j = jrtst; j <= i__2; ++j)
        {
            nodej = npn[j];
            lstrt = xadj[nodej];
            lstop = xadj[nodej + 1] - 1;
            i__3 = jstop;
            for (k = j + 1; k <= i__3; ++k)
            {
                nodek = npn[k];
                i__4 = lstop;
                for (l = lstrt; l <= i__4; ++l)
                {
                    if (adj[l] == nodek)
                    {
                        goto L70;
                    }
                    if (adj[l] == 0)
                    {
                        goto L55;
                    }
/* L50: */
                }
                cerr << fmt_10900;
                *nop = 900;
                return 0;
        L55:
                adj[l] = nodek;
                mstrt = xadj[nodek];
                mstop = xadj[nodek + 1] - 1;
                i__4 = mstop;
                for (m = mstrt; m <= i__4; ++m)
                {
                    if (adj[m] == 0)
                    {
                        goto L65;
                    }
/* L60: */
                }
                cerr << fmt_10900;
                *nop = 900;
                return 0;
        L65:
                adj[m] = nodej;
        L70:
                ;
            }
/* L80: */
        }
/* L90: */
    }

/*     Strip any zeros from adjacency list */

    k = 0;
    jstrt = 1;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        jstop = xadj[i + 1] - 1;
        i__2 = jstop;
        for (j = jstrt; j <= i__2; ++j)
        {
            if (adj[j] == 0)
            {
                goto L105;
            }
            ++k;
            adj[k] = adj[j];
/* L100: */
        }
L105:
        xadj[i + 1] = k + 1;
        jstrt = jstop + 1;
/* L110: */
    }
    *iadj = xadj[*n + 1] - 1;
    *iadj = (*iadj > 1) ? *iadj : 1;

    return 0;
}                               /* gegra_ */
