#include <stdint.h>

/* @brief Purpose: Number nodes in component of graph for small profile and rms wavefront */
int number_ (int64_t *, int64_t *nc, int64_t *snode, int64_t *lstnum, int64_t *, int64_t *adj, int64_t *xadj,
			 int64_t *s, int64_t *q, int64_t *p)
//n, nc, snode, lstnum,e2, adj, xadj, s, q, p
{
 /* System generated locals */
    int64_t     i__1, i__2;
 /* Local variables */
    static int64_t node, next, prty, i, j, nabor, istop, jstop, istrt, jstrt, nn, addres, maxprt, nbr;

/*     INPUT: */
/*     ------ */
/*     N      - Number of nodes in graph */
/*     NC     - Number of nodes in component of graph */
/*     SNODE  - Node of which numbering starts */
/*     LSTNUM - Count of nodes which have already been numbered */
/*     E2     - Twice the number of edges in the graph = XADJ(N+1)-1 */
/*     ADJ    - Adjacency list for all nodes in graph */
/*            - List of length 2E where E is the number of edges in */
/*              the graph and 2E = XADJ(N+1)-1 */
/*     XADJ   - Index number for ADJ */
/*            - Nodes adjacent to node I are found in ADJ(I), where */
/*              J=XADJ(I), XADJ(I)+1, ..., XADJ(I+1)-1 */
/*     S      - List giving the distance of each node in this */
/*              component */
/*     Q      - List of nodes which are in this component */
/*            - Also used to store queue of active or preactivnodes */
/*     P      - Undefined */

/*     OUTPUT: */
/*     ------- */
/*     N      - Unchanged */
/*     NC     - Unchanged */
/*     SNODE  - Unchanged */
/*     LSTNUM - Count of numbered nodes (input value incremented by NC) */
/*     E2     - Unchanged */
/*     ADJ    - Unchanged */
/*     XADJ   - Unchanged */
/*     S      - List of new node numbers */
/*            - New number for node, I is S(I) */
/*     Q      - Not used */
/*     P      - Not used */

/*     NOTES: */
/*     ------ */
/*     S also serves as a list giving the status of the nodes */
/*     during the numbering process: */
/*     S(I) gt 0 indicates node i is postactive */
/*     S(I) =  0 indicates node i is active */
/*     S(I) = -1 indicates node i is preactive */
/*     S(I) = -2 indicates node i is inactive */
/*     P is used to hold the priorities for each node */

/*     PROGRAMMER:    Scott sloan */
/*     ----------- */

/*     LAST MODIFIED: 10 March 1989    Scott Sloan */
/*     -------------- */
/* *********************************************************************** */


/*     Initialise priorities and status for each node in this component */
/*     Initial priority = W1*DIST - W2*DEGREE    where: */
/*     W1     = a positive weigth */
/*     W2     = a positive weigth */
/*     DEGREE = initial current degree for node */
/*     DIST   = distance of node from end node */

 /* Parameter adjustments */
    --p;
    --s;
    --xadj;
    --q;
    --adj;

 /* Function Body */
    i__1 = *nc;
    for (i = 1; i <= i__1; ++i)
    {
        node = q[i];
        p[node] = s[node] - ((xadj[node + 1] - xadj[node] + 1) << 1);
        s[node] = -2;
/* L10: */
    }

/*     Insert starting node in queue and assign it a preactive status */
/*     NN  is the size of queue */

    nn = 1;
    q[nn] = *snode;
    s[*snode] = -1;

/*     Loop while queue is not empty */

L30:
    if (nn > 0)
    {

/*      Scan queue for node with max prioity */

        addres = 1;
        maxprt = p[q[1]];
        i__1 = nn;
        for (i = 2; i <= i__1; ++i)
        {
            prty = p[q[i]];
            if (prty > maxprt)
            {
                addres = i;
                maxprt = prty;
            }
/* L35: */
        }

/*      NEXT is the node to be numbered next */

        next = q[addres];

/*      Delete node NEXT from queue */

        q[addres] = q[nn];
        --nn;
        istrt = xadj[next];
        istop = xadj[next + 1] - 1;
        if (s[next] == -1)
        {

/*       Node NEXT is preactive, examine its neighbours */

            i__1 = istop;
            for (i = istrt; i <= i__1; ++i)
            {

/*          Decrease current degree of neighbour by -1 */

                nbr = adj[i];
                p[nbr] += 2;

/*          Add neighbour to queue if it is inactive */
/*          assign it a preactive status */

                if (s[nbr] == -2)
                {
                    ++nn;
                    q[nn] = nbr;
                    s[nbr] = -1;
                }
/* L50: */
            }
        }
/*      Store new node number for node NEXT */
/*      Status for node NEXT is now postactive */

        ++(*lstnum);
        s[next] = *lstnum;

/*      Search for preactive neightbours of node NEXT */

        i__1 = istop;
        for (i = istrt; i <= i__1; ++i)
        {
            nbr = adj[i];
            if (s[nbr] == -1)
            {

/*          Decrease current degree of preactive neighbour by -1 */
/*          assign neighbour an active status */

                p[nbr] += 2;
                s[nbr] = 0;

/*          Loop over nodes adjacent to preactive neighbour */

                jstrt = xadj[nbr];
                jstop = xadj[nbr + 1] - 1;
                i__2 = jstop;
                for (j = jstrt; j <= i__2; ++j)
                {
                    nabor = adj[j];

/*             Decrease current degree of adjacent node by -1 */

                    p[nabor] += 2;
                    if (s[nabor] == -2)
                    {

/*              Insert inactive node in queue with a p reactive status */

                        ++nn;
                        q[nn] = nabor;
                        s[nabor] = -1;
                    }
/* L60: */
                }
            }
/* L80: */
        }
        goto L30;
    }
    return 0;
}                               /* number_ */
