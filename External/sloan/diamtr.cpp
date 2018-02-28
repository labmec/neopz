#include "sloan.h"

/** @brief Purpose: Find nodes which define a pseudo-diameter of a graph and store distances from end node */
int diamtr_ (int64_t *n, int64_t *e2, int64_t *adj, int64_t *
         xadj, int64_t *mask, int64_t *ls, int64_t *xls, int64_t *hlevel, int64_t *snode, int64_t *nc)
{
 /* System generated locals */
    int64_t     i__1, i__2;
 /* Local variables */
    static int64_t node, i, j, enode, depth, width, hsize, istop, jstop, istrt, jstrt, degree, mindeg, ewidth, sdepth;

/*     INPUT: */
/*     ------ */
/*     N       - The total number of nodes in the graph */
/*     E2      - Twice the number of edges in the graph = XADJ(N+1)-1 */
/*     ADJ     - Adjacency list for all nodes in the graph */
/*             - List of length 2E where E is the number of edges in */
/*               the graph and 2E = XADJ(N+1)-1 */
/*     XADJ    - Index vector for ADJ */
/*             - Nodes adjacent to node I are found in ADJ(J), where */
/*               J = XADJ(I),XADJ(I)+1, ..., XADJ(I+1)-1 */
/*             - Degree of node I given by XADJ(I=1)-XADJ(I) */
/*     MASK    - Masking vector for graph */
/*             - Visible nodes have MASK = 0, node invisible otherwise */
/*     LS      - Undefined */
/*     XLS     - Undefined */
/*     HLEVEL  - Undefined */
/*     SNODE   - Undefined */
/*     NC      - Undefined */

/*     OUTPUT: */
/*     ------- */
/*     N       - Unchanged */
/*     E2      - Unchanged */
/*     ADJ     - Unchanged */
/*     XADJ    - Unchanged */
/*     MASK    - List of distances of nodes from the end node */
/*     LS      - List of nodes in the component */
/*     XLS     - Not used */
/*     HLEVEL  - Not used */
/*     SNODE   - Starting node for numbering */
/*     NC      - The number of nodes in this component of graph */

/*     SUBROUTINES CALLED:  ROOTLS, ISORTI */
/*     ------------------- */
/*     NOTE:       SNODE and ENODE define a pseudo-diameter */

/*     PROGRAMMER: Scott Sloan */
/*     ----------- */

/*     LAST MODIFIED:  10 March 1989     Scott Sloan */
/* ***********************************************************************
 */


/*  Choose first guess for starting node by min degree */
/*  Ignore nodes that are invisible (MASK NE 0) */

 /* Parameter adjustments */
    --hlevel;
    --xls;
    --ls;
    --mask;
    --xadj;
    --adj;

 /* Function Body */
    mindeg = *n;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        if (mask[i] == 0)
        {
            degree = xadj[i + 1] - xadj[i];
            if (degree < mindeg)
            {
                *snode = i;
                mindeg = degree;
            }
        }
/* L10: */
    }

/*  Generate level structure for node with min degree */

    i__1 = *n + 1;
    rootls_ (n, snode, &i__1, e2, &adj[1], &xadj[1], &mask[1], &ls[1], &xls[1],
             &sdepth, &width);

/*  Store number of nodes in this component */

    *nc = xls[sdepth + 1] - 1;

/*  Iterate to find start and end nodes */

L15:

/*  Store list of nodes that are at max distance from starting node */
/*  Store their degrees in XLS */
    hsize = 0;
    istrt = xls[sdepth];
    istop = xls[sdepth + 1] - 1;
    i__1 = istop;
    for (i = istrt; i <= i__1; ++i)
    {
        node = ls[i];
        ++hsize;
        hlevel[hsize] = node;
        xls[node] = xadj[node + 1] - xadj[node];
/* L20: */
    }

/*  Sort list of nodes in ascending sequence of their degree */
/*  Use insertion sort algorithm */
    if (hsize > 1)
    {
        isorti_ (&hsize, &hlevel[1], n, &xls[1]);
    }

/*  Remove nodes with duplicate degrees */

    istop = hsize;
    hsize = 1;
    degree = xls[hlevel[1]];
    i__1 = istop;
    for (i = 2; i <= i__1; ++i)
    {
        node = hlevel[i];
        if (xls[node] != degree)
        {
            degree = xls[node];
            ++hsize;
            hlevel[hsize] = node;
        }
/* L25: */
    }

/*  Loop over nodes in skrunken level */

    ewidth = *nc + 1;
    i__1 = hsize;
    for (i = 1; i <= i__1; ++i)
    {
        node = hlevel[i];

/*     Form rooted level structures for each node in skrunken level */

        rootls_ (n, &node, &ewidth, e2, &adj[1], &xadj[1], &mask[1], &ls[1], &
                 xls[1], &depth, &width);
        if (width < ewidth)
        {

/*         Level structure was not aborted during assembly */

            if (depth > sdepth)
            {

/*          Level structure of greater depth found */
/*          Store new starting node, new max depth, and begin */
/*          a new iteration */

                *snode = node;
                sdepth = depth;
                goto L15;
            }
/*         Level struture width for this end node is smallest so far */
/*         store end node and new min width */

            enode = node;
            ewidth = width;
        }
/* L30: */
    }

/* Generate level structure rooted at end node if necessary */

    if (node != enode)
    {
        i__1 = *nc + 1;
        rootls_ (n, &enode, &i__1, e2, &adj[1], &xadj[1], &mask[1], &ls[1], &
                 xls[1], &depth, &width);
    }
/* Store distances of each node from end node */

    i__1 = depth;
    for (i = 1; i <= i__1; ++i)
    {
        jstrt = xls[i];
        jstop = xls[i + 1] - 1;
        i__2 = jstop;
        for (j = jstrt; j <= i__2; ++j)
        {
            mask[ls[j]] = i - 1;
/* L40: */
        }
/* L50: */
    }
    return 0;
}
