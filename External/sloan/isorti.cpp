#include <stdint.h>

/* Subroutine */ 
int isorti_ (int64_t *nl, int64_t *list, int64_t *, int64_t *key)
//nl, list, nk, key
{
    static int64_t i, j, k, l, m;
    static double r, t;
    static int64_t ij, il[21], iu[21];
    static double tt;
    static int64_t itt;

/*                                  SPECIFICATIONS FOR ARGUMENTS */
/*                                  SPECIFICATIONS FOR LOCAL VARIABLES */
/*                                  FIRST EXECUTABLE STATEMENT */
/*                                  FIND ABSOLUTE VALUES OF ARRAY A */
 
/* Parameter adjustments */
    --list;
    --key;

/* Function Body */
    if (*nl <= 0)
    {
        return 0;
    }
    m = 1;
    i = 1;
    j = *nl;
    r = (float) .375;
L10:
    if (i == j)
    {
        goto L55;
    }
/* L15: */
    if (r > (float) .5898437)
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
    ij = i + (int)((j - i) * r);
    t = (double) key[list[ij]];

/*                                  IF FIRST ELEMENT OF ARRAY IS GREATER */
/*                                  THAN T, INTERCHANGE WITH T */
    if ((double) key[list[i]] > t)
    {
        tt = (double) list[ij];
        list[ij] = list[i];
        list[i] = (int)tt;
        t = (double) key[list[ij]];
/*          IR(IJ)=IR(I) */
/*          IR(I)=IT */
/*          IT=IR(IJ) */
    }
    l = j;
/*                                  IF LAST ELEMENT OF ARRAY IS LESS THAN */
/*                                  T, INTERCHANGE WITH T */
    if ((double) key[list[j]] < t)
    {
        tt = (double) list[ij];
        list[ij] = list[j];
        list[j] = (int)tt;
        t = (double) key[list[ij]];
/*                                  IF FIRST ELEMENT OF ARRAY IS GREATER */
/*                                  THAN T, INTERCHANGE WITH T */
        if ((double) key[list[i]] > t)
        {
            tt = (double) list[ij];
            list[ij] = list[i];
            list[i] = (int)tt;
            t = (double) key[list[ij]];
        }
    }
/*                                  FIND AN ELEMENT IN THE SECOND HALF OF */
/*                                  THE ARRAY WHICH IS SMALLER THAN T */
L40:
    --l;
    if ((double) key[list[l]] > t)
    {
        goto L40;
    }
/*                                  FIND AN ELEMENT IN THE FIRST HALF OF */
/*                                  THE ARRAY WHICH IS GREATER THAN T */
L45:
    ++k;
    if ((double) key[list[k]] < t)
    {
        goto L45;
    }
/*                                  INTERCHANGE THESE ELEMENTS */
    if (k <= l)
    {
        if (key[list[l]] != key[list[k]])
        {
            tt = (double) list[l];
            list[l] = list[k];
            list[k] = (int)tt;
        }
        goto L40;
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
    t = (double) key[list[i + 1]];
    itt = list[i + 1];
    if ((double) key[list[i]] <= t)
    {
        goto L65;
    }
    k = i;
L70:
    list[k + 1] = list[k];
    --k;
    if (t < (double) key[list[k]])
    {
        goto L70;
    }
    list[k + 1] = itt;
    goto L65;
}                               /* isorti_ */
