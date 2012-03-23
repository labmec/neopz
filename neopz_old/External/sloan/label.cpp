#include "sloan.h" 
 /* Subroutine */ int 
label_ (int *n, int *e2, int *adj, int *
        xadj, int *nnn, int *iw, int *oldpro, int *newpro)
{
 /* System generated locals */
    int     i__1;
 /* Local variables */
    static int i, snode, i1, i2, i3;
    static int nc;
    static int lstnum;


/*     PURPOSE: */
/*     -------- */

/*     Label a graph for small profile and rms wavefront */

/*     INPUT: */
/*     ------ */

/*     N      - Total number of nodes in graph */
/*     E2     - Twice the number of edges in the graph = XADJ(N+1)-1 */
/*     ADJ    - Adjacency list for all nodes in graph */
/*            - List of length 2E where E is the number of edges in */
/*              the graph and 2E = XADJ(N+1)-1 */
/*     XADJ   - Index vector for ADJ */
/*            - Nodes adjacentto node I are found in ADJ(J), where */
/*              J = XADJ(I),XADJ(I)+1, ..., XADJ(I+1)-1 */
/*            - Degreeof node I given by XADJ(I+1)-XADJ(I) */
/*     NNN    - Undefined */
/*     IW     - Undefined */
/*     OLDPRO - Undefined */
/*     NEWPRO - Undefined */

/*     OUTPUT: */
/*     ------- */

/*     N      - Unchanged */
/*     E2     - Unchanged */
/*     ADJ    - Unchanged */
/*     XADJ   - Unchanged */
/*     NNN    - List of new node numbers */
/*            - New number for node I given by NNN(I) */
/*            - If original node numbers give a smaller profile then */
/*              NNN is set so that NNN(I)=I for I=1,N */
/*     IW     - Not used */
/*     OLDPRO - Profile using original node numbering */
/*     NEWPRO - Profile for new node numbering */
/*            - If original profile is smaller than new profile, then */
/*            - original node numbers are used and NEWPRO=OLDPRO */

/*      SUBROUTINES CALLED: DIAMTR, NUMBER, PROFIL */
/*      ------------------- */

/*      PROGRAMMER:     Scott Sloan */
/*      ----------- */

/*      LAST MODIFIED:  10 March 1989     Scott Sloan */
/*      -------------- */

/* ***********************************************************************
 */


/*     Set all new node numbers =0 */
/*     This is used to denote all visible nodes */

 /* Parameter adjustments */
    --iw;
    --nnn;
    --xadj;
    --adj;

 /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        nnn[i] = 0;
/* L10: */
    }

/*     Define offsets */

    i1 = 1;
    i2 = i1 + *n;
    i3 = i2 + *n + 1;

/*     Loop while some nodes remain unnumbered */

    lstnum = 0;
L20:
    if (lstnum < *n)
    {

/*     Find end points of p-diameter for nodes in this component */
/*     Compute distances of nodes from end node */

        diamtr_ (n, e2, &adj[1], &xadj[1], &nnn[1], &iw[i1], &iw[i2], &iw[i3],
                 &snode, &nc);

/*      Number nodes in this component */

        number_ (n, &nc, &snode, &lstnum, e2, &adj[1], &xadj[1], &nnn[1], &iw[
                                                           i1], &iw[i2]);
        goto L20;
    }
/*     Compute profiles for old and new node numbers */

    profi1_ (n, &nnn[1], e2, &adj[1], &xadj[1], oldpro, newpro);

/*     Use original numbering if it gives a smaller profile */

/*    if (*oldpro < *newpro)
    {
        i__1 = *n;
        for (i = 1; i <= i__1; ++i)
        {
            nnn[i] = i;
//  L30: 
        }
        *newpro = *oldpro;
    }
*/
	++adj;
	++xadj;

    return 0;
}                               /* label_ */
