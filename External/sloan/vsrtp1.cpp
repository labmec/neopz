 /* Subroutine */ int 
vsrtp1_ (int *a, int *ir, int *la)
{
 /* System generated locals */
    int     i__1;
 /* Local variables */
    static int i, j, k, l, m, r, t, ij, il[21], it, iu[21], tt, itt;


/*                                  SPECIFICATIONS FOR ARGUMENTS */
/*                                  SPECIFICATIONS FOR LOCAL VARIABLES */
/*                                  FIRST EXECUTABLE STATEMENT */
/*                                  FIND ABSOLUTE VALUES OF ARRAY A */
 /* Parameter adjustments */
    --ir;
    --a;

 /* Function Body */
    if (*la <= 0)
    {
        return 0;
    }
    i__1 = *la;
    for (i = 1; i <= i__1; ++i)
    {
/*         IF (A(I) .LT. 0.0) A(I)=-A(I) */
        ir[i] = i;
/* L5: */
    }
    m = 1;
    i = 1;
    j = *la;
    r = (float) .375;
L10:
    if (i == j)
    {
        goto L55;
    }
/* L15: */
    if ((double) r > (float) .5898437)
    {
        r += (float) -.21875;
    } else
    {
        r += (float) .0390625;
    }
L25:
    k = i;
/*                                  SELECT A CENTRAL ELEMENT OF THE */
/*                                  ARRAY AND SAVE IT IN LOCATION T */
    ij = i + (j - i) * r;
    t = a[ij];
    it = ir[ij];
/*                                  IF FIRST ELEMENT OF ARRAY IS GREATER
*/
/*                                  THAN T, INTERCHANGE WITH T */
    if (a[i] > t)
    {
        a[ij] = a[i];
        a[i] = t;
        t = a[ij];
        ir[ij] = ir[i];
        ir[i] = it;
        it = ir[ij];
    }
    l = j;
/*                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
*/
/*                                  T, INTERCHANGE WITH T */
    if (a[j] < t)
    {
        a[ij] = a[j];
        a[j] = t;
        t = a[ij];
        ir[ij] = ir[j];
        ir[j] = it;
        it = ir[ij];
/*                                  IF FIRST ELEMENT OF ARRAY IS GREAT
ER */
/*                                  THAN T, INTERCHANGE WITH T */
        if (a[i] > t)
        {
            a[ij] = a[i];
            a[i] = t;
            t = a[ij];
            ir[ij] = ir[i];
            ir[i] = it;
            it = ir[ij];
        }
    }
    goto L40;
L35:
    if (a[l] != a[k])
    {
        tt = a[l];
        a[l] = a[k];
        a[k] = tt;
        itt = ir[l];
        ir[l] = ir[k];
        ir[k] = itt;
    }
/*                                  FIND AN ELEMENT IN THE SECOND HALF OF
*/
/*                                  THE ARRAY WHICH IS SMALLER THAN T */
L40:
    --l;
    if (a[l] > t)
    {
        goto L40;
    }
/*                                  FIND AN ELEMENT IN THE FIRST HALF OF
*/
/*                                  THE ARRAY WHICH IS GREATER THAN T */
L45:
    ++k;
    if (a[k] < t)
    {
        goto L45;
    }
/*                                  INTERCHANGE THESE ELEMENTS */
    if (k <= l)
    {
        goto L35;
    }
/*                                  SAVE UPPER AND LOWER SUBSCRIPTS OF */
/*                                  THE ARRAY YET TO BE SORTED */
    if (l - i > j - k)
    {
        il[m - 1] = i;
        iu[m - 1] = l;
        i = k;
        ++m;
    } else
    {
        il[m - 1] = k;
        iu[m - 1] = j;
        j = l;
        ++m;
    }
    goto L60;
/*                                  BEGIN AGAIN ON ANOTHER PORTION OF */
/*                                  THE UNSORTED ARRAY */
L55:
    --m;
    if (m == 0)
    {
        return 0;
    }
    i = il[m - 1];
    j = iu[m - 1];
L60:
    if (j - i >= 11)
    {
        goto L25;
    }
    if (i == 1)
    {
        goto L10;
    }
    --i;
L65:
    ++i;
    if (i == j)
    {
        goto L55;
    }
    t = a[i + 1];
    it = ir[i + 1];
    if (a[i] <= t)
    {
        goto L65;
    }
    k = i;
L70:
    a[k + 1] = a[k];
    ir[k + 1] = ir[k];
    --k;
    if (t < a[k])
    {
        goto L70;
    }
    a[k + 1] = t;
    ir[k + 1] = it;
    goto L65;
}                               /* vsrtp1_ */
