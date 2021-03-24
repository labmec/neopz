#include <stdint.h>

/* Subroutine */ 
int rootls_ (int64_t *, int64_t *root, int64_t *maxwid,
         int64_t *, int64_t *adj, int64_t *xadj, int64_t *mask, int64_t *ls,
         int64_t *xls, int64_t *depth, int64_t *width)
//n, root, maxwid, e2, adj, xadj, mask, ls, xls, depth, width
{
 /* System generated locals */
    int64_t     i__1, i__2;
 /* Local variables */
    static int64_t node, i, j, lwdth, jstop, lstop, jstrt, lstrt, nc, nbr;


/*     PURPOSE: */
/*     -------- */

/*     Generate rooted level structure using a FORTRAN 77  implementation
*/
/*     of the algorithm given by George and Liu */

/*     INPUT: */
/*     ------ */

/*     N      - Numbered of nodes */
/*     ROOT   - Root node for level structure */
/*     MAXWID - Max permissible width ot rooted level structure */
/*            - Abort assembly of level structure if width is ge MAXWID */
/*            - Assembly ensured by setting MAXWID =N+1 */
/*     E2     - Twice the number of edges in the graph = XADJ(N+1)-1 */
/*     ADJ    - Adjacency list for all nodes in graph */
/*            - List of length 2E where E is the number of edges in */
/*              the graph and 2E = XADJ(N+1)-1 */
/*     XADJ   - Index vector for adj */
/*            - Nodes adjacent to node I are found in ADJ(J), where */
/*              J = XADJ(I), XADJ(I)+1, ..., XADJ(I+1)-1 */
/*            - Degree of node I is XADJ(I+1)-XADJ(I) */
/*     MASK   - Masking vector for graph */
/*            - Visible nodes have MASK = 0 */
/*     LS     - Undefined */
/*     XLS    - Undefined */
/*     DEPTH  - Undefined */
/*     WIDTH  - Undefined */

/*     OUTPUT: */
/*     ------- */

/*     N      - Unchanged */
/*     ROOT   - Unchanged */
/*     MAXWID - Unchanged */
/*     E2     - Unchanged */
/*     ADJ    - Unchanged */
/*     XADJ   - Unchanged */
/*     MASK   - Unchanged */
/*     LS     - List containing a rooted level structure */
/*            - List of length NC */
/*     XLS    - Index vector for LS */
/*            - Nodes in level I are found in LS(J), where */
/*              J=XLS(I), XLS(I)+1, ..., XLS(I+1)-1 */
/*            - List of max length NC+1 */
/*     DEPTH  - Number of levels in rooted level structure */
/*     WIDTH  - Width of rooted level structure */

/*     NOTE:  If WIDTH ge MAXWID then assembly has been aborted */
/*     ----- */

/*     PROGRAMMER:     Scott Sloan */
/*     ----------- */

/*     LAST MODIFIED:  10 March 1989      Scott Sloan */
/*     -------------- */

/* ***********************************************************************
 */


/*     Initialisation */

 /* Parameter adjustments */
    --xls;
    --ls;
    --mask;
    --xadj;
    --adj;

 /* Function Body */
    mask[*root] = 1;
    ls[1] = *root;
    nc = 1;
    *width = 1;
    *depth = 0;
    lstop = 0;
    lwdth = 1;
L10:
    if (lwdth > 0)
    {

/*      LWDTH is the width of the current level */
/*      LSTRT points to start of current level */
/*      LSTOP points to end of current level */
/*      NC    counts the nodes in component */

        lstrt = lstop + 1;
        lstop = nc;
        ++(*depth);
        xls[*depth] = lstrt;

/*      Generate next level by finding all visible neighbours */
/*      of nodes in current level */

        i__1 = lstop;
        for (i = lstrt; i <= i__1; ++i)
        {
            node = ls[i];
            jstrt = xadj[node];
            jstop = xadj[node + 1] - 1;
            i__2 = jstop;
            for (j = jstrt; j <= i__2; ++j)
            {
                nbr = adj[j];
                if (mask[nbr] == 0)
                {
                    ++nc;
                    ls[nc] = nbr;
                    mask[nbr] = 1;
                }
/* L20: */
            }
/* L30: */
        }

/*      Compute width of level just assembled and the width of the */
/*      level structure so far */

        lwdth = nc - lstop;
        *width = ((lwdth) > (*width)) ? (lwdth) : (*width);
//      *width = max(lwdth,*width);

/*      Abort assembly if level structure is too wide */

        if (*width >= *maxwid)
        {
            goto L35;
        }
        goto L10;
    }
    xls[*depth + 1] = lstop + 1;

/*     Reset MASK = 0 for nodes in the level structure */

L35:
    i__1 = nc;
    for (i = 1; i <= i__1; ++i)
    {
        mask[ls[i]] = 0;
/* L40: */
    }
    return 0;
}                               /* rootls_ */
