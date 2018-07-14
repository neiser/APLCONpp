/* condutil.f90.F -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
        on Microsoft Windows system, link with libf2c.lib;
        on Linux or Unix systems, link with .../path/to/libf2c.a -lm
        or, if you install libf2c.a in a standard place, with -lf2c -lm
        -- in that order, at the end of the command line, as in
                cc *.o -lf2c -lm
        Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

                http://www.netlib.org/f2c/libf2c.zip
*/

#include <f2c.h>

integer ijsym_(integer *i__, integer *j) {
  /* System generated locals */
  integer ret_val;

  /* index (I,J)=(J,I) in symmetric matri */
  /* ,NCOUNT */
  if (*i__ <= *j) {
    ret_val = (*j * *j - *j) / 2 + *i__;
  } else {
    ret_val = (*i__ * *i__ - *i__) / 2 + *j;
  }
  return ret_val;
} /* ijsym_ */

/* Subroutine */ int duminv_(doublereal *a, doublereal *w, doublereal *b,
                             integer *nx, integer *nf, integer *mb,
                             integer *nrank, doublereal *aux,
                             doublereal *qnext) {
  /* Initialized data */

  static doublereal eps = 1e-6;

  /* System generated locals */
  integer i__1, i__2, i__3;
  doublereal d__1, d__2, d__3;

  /* Local variables */
  static integer i__, j, k, l, m, n, ia, ij, jk, jj, kk, jl, lk;
  static doublereal vjk, vkk, sum;
  static integer last, nmeas, jlast;
  static logical solve;
  static integer jfirst;

  /*     Obtain solution of a system of linear equations V *  X  =  B  with */
  /*     symmetric matrix V and inverse (for M =  1)  or  matrix  inversion */
  /*     only (for M = 0). */

  /*                   - - - - */
  /*        CALL SMINV(W,B,N,M,NRANK,AUX,QNEXT) */
  /*                   - -     ----- */

  /*           W = symmetric N-by-N matrix in symmetric storage mode */
  /*               W(1) = W11, W(2) = W12, W(3) = W22, W(4) = W13, . . . */
  /*               replaced by inverse matrix */
  /*           B = N-vector   (for M = 0 use a dummy argument) */
  /*               replaced by solution vector */
  /*           M = see above */

  /*     Method of solution is by elimination selecting the  pivot  on  the */
  /*     diagonal each stage. The rank of the matrix is returned in  NRANK. */
  /*     For NRANK ne N, all remaining  rows  and  cols  of  the  resulting */
  /*     matrix V and the corresponding elements of  B  are  set  to  zero. */
  /*     SMINV can be used for a dimension up to 100 (see INVCDR). */

  /* matrix inve */
  /* Parameter adjustments */
  --qnext;
  --aux;
  --b;
  --w;
  --a;

  /* Function Body */
  /*     special entry for partial inversion ****************************** */
  /*     ... */
  *nrank = 0;
  n = *nx + *nf;
  /*     make sure AUX is zero, prevents uninit access */
  /*     in continued execution of DBMINV */
  /* dimension parameter */
  i__1 = n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    aux[i__] = 0.;
  }
  /*     -VX(NX-sym) is already inserted in W(NX+NF-sym) ------------------ */
  ij = (*nx * *nx + *nx) / 2;
  /* number of elements of V(.) */
  i__1 = n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    qnext[i__] = 0.;
    /* reset pointer */
  }
  ia = 0;
  i__1 = *nf;
  for (j = 1; j <= i__1; ++j) {
    i__2 = *nx;
    for (i__ = 1; i__ <= i__2; ++i__) {
      w[ij + i__] = a[ia + i__];
      /* copy A(.) into W_12 */
    }
    i__2 = j;
    for (i__ = 1; i__ <= i__2; ++i__) {
      w[ij + *nx + i__] = 0.f;
      /* reset last submatrix W_22 of W(.) */
    }
    ij = ij + *nx + j;
    ia += *nx;
  }
  /*     distinguish between measured and unmeasured variables ------------ */
  jfirst = 0;
  /* first index of measured variable */
  nmeas = 0;
  /* number of measured variables */
  i__1 = *nx;
  for (i__ = 1; i__ <= i__1; ++i__) {
    if (w[(i__ * i__ + i__) / 2] < 0.f) {
      /* measured variable */
      if (jfirst == 0) {
        jfirst = i__;
        /* first index of measured variable */
      } else {
        qnext[jlast] = (doublereal)i__;
        /* insert index at previous index */
      }
      jlast = i__;
      /* save index */
      ++nmeas;
    }
  }
  if (jlast == 0) {
    goto L10;
  }
  /* nothing to do */
  qnext[jlast] = -1.;
  /*     apply exchange algorithm to sub-matrices ------------------------- */
  /* stop index for last measured variable */
  i__1 = n;
  for (i__ = *nx + 1; i__ <= i__1; ++i__) {
    /* loop I over constraint equations */
    j = jfirst;
    /* first index of unmeasured variable */
    i__2 = nmeas;
    for (m = 1; m <= i__2; ++m) {
      /* already inverted element index J */
      sum = 0.;
      jk = (j * j - j) / 2;
      /* index of diagonal element before */
      i__3 = *nx;
      for (k = 1; k <= i__3; ++k) {
        if (k <= j) {
          ++jk;
        }
        /* index in j column */
        if (qnext[k] != 0.) {
          sum += w[jk] * w[(i__ * i__ - i__) / 2 + k];
        }
        if (k >= j) {
          jk += k;
        }
        /* index in j row */
      }
      aux[j] = sum;
      /* = A-row * VX-row/col */
      j = (integer)qnext[j];
      /* next index of unmeasured variable */
    }
    i__2 = n;
    for (k = i__; k <= i__2; ++k) {
      sum = 0.;
      j = jfirst;
      /* first index of unmeasured variable */
      i__3 = nmeas;
      for (m = 1; m <= i__3; ++m) {
        /* already inverted element index J */
        sum += w[(k * k - k) / 2 + j] * aux[j];
        /* = A-row * H */
        j = (integer)qnext[j];
        /* next index of unmeasured variable */
      }
      w[(k * k - k) / 2 + i__] += sum;
      /* add to diagonal W_22 */
    }
    j = jfirst;
    /* first index of unmeasured variable */
    i__2 = nmeas;
    for (m = 1; m <= i__2; ++m) {
      w[(i__ * i__ - i__) / 2 + j] = -aux[j];
      /* add to off-diagonal W_22 */
      j = (integer)qnext[j];
      /* next index of unmeasured variable */
    }
  }
  /*     set pointer for unmeasured variables ----------------------------- */
  jfirst = 0;
  jlast = 0;
  i__1 = n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    if (qnext[i__] == 0.) {
      /* unmeasured variable */
      if (jfirst == 0) {
        jfirst = i__;
        /* first index of unmeasured variable */
      } else {
        qnext[jlast] = (doublereal)i__;
        /* next index of unmeasured variable */
      }
      jlast = i__;
    } else {
      qnext[i__] = 0.;
      /* reset index for measured variable */
    }
  }
  if (jlast == 0) {
    goto L10;
  }
  /* no unmeasured variable */
  qnext[jlast] = -1.;
/*     common code for inversion and (M=1) solution of matrix equation */
/* end flag */
L10:
  solve = TRUE_;
  if (*mb == 0) {
    solve = FALSE_;
  }
  /*     loop begin (loop on all remaining rows/cols) */
  /* solution flag */
  i__1 = n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* loop on all remaining elements */
    vkk = 0.;
    /* search for pivot element */
    k = 0;
    /* pivot index */
    j = jfirst;
    /* first candidate index */
    last = 0;
  L20:
    if (j > 0) {
      /* test for linearity and zero matrix */
      jj = (j * j + j) / 2;
      /* diagonal index */
      /* Computing MAX */
      d__2 = abs(vkk), d__3 = eps * aux[j];
      if ((d__1 = w[jj], abs(d__1)) > max(d__2, d__3)) {
        vkk = w[jj];
        /* largest pivot candidate so far */
        k = j;
        /* index of largest */
        l = last;
      }
      last = j;
      j = (integer)qnext[j];
      /* index of next candidate */
      goto L20;
    }
    if (k != 0) {
      /* pivot element found - proceed */
      ++(*nrank);
      /* increase rank counter */
      kk = (k * k + k) / 2;
      if (l == 0) {
        jfirst = (integer)qnext[k];
        /* new first index */
      } else {
        qnext[l] = qnext[k];
        /* bridge used index */
      }
      qnext[k] = 0.;
      /* reset used index */
      ++(*nrank);
      /* increase rank */
      vkk = 1.f / vkk;
      /* invert pivot */
      w[kk] = -vkk;
      if (solve) {
        b[k] *= vkk;
      }
      jk = kk - k;
      jl = 0;
      i__2 = n;
      for (j = 1; j <= i__2; ++j) {
        /* elimination */
        if (j == k) {
          jk = kk;
          jl += j;
        } else {
          if (j < k) {
            ++jk;
          } else {
            jk = jk + j - 1;
          }
          vjk = w[jk];
          w[jk] = vkk * vjk;
          if (solve) {
            b[j] -= b[k] * vjk;
          }
          lk = kk - k;
          i__3 = j;
          for (l = 1; l <= i__3; ++l) {
            ++jl;
            if (l == k) {
              lk = kk;
            } else {
              if (l < k) {
                ++lk;
              } else {
                lk = lk + l - 1;
              }
              w[jl] -= w[lk] * vjk;
            }
          }
        }
      }
    } else {
      /* no pivot candadate found - reset */
      i__2 = n;
      for (k = 1; k <= i__2; ++k) {
        if (qnext[k] != 0.) {
          /* undefined variable */
          if (solve) {
            b[k] = 0.;
          }
          /* clear undefined vector element */
          i__3 = k;
          for (j = 1; j <= i__3; ++j) {
            if (qnext[j] != 0.) {
              w[(k * k - k) / 2 + j] = 0.;
            }
            /* clear matri */
          }
        }
      }
      goto L30;
    }
  }
/* end of inversion loop */
L30:
  i__1 = (n * n + n) / 2;
  for (i__ = 1; i__ <= i__1; ++i__) {
    w[i__] = -w[i__];
    /* finally reverse sign */
  }
  return 0;
} /* duminv_ */

doublereal scalxy_(doublereal *x, doublereal *y, integer *n) {
  /* System generated locals */
  integer i__1;
  doublereal ret_val;

  /* Local variables */
  static integer j;
  static doublereal sum;

  /*     Scalar product of two vectors */
  /*                - - -                 T */
  /*        S = VXY(X,Y,N)           S = X  * Y (scalar product) */

  /* scalar vector product */
  /* I,M */
  /*     ... */
  /* Parameter adjustments */
  --y;
  --x;

  /* Function Body */
  sum = 0.;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    sum += x[j] * y[j];
  }
  ret_val = sum;
  return ret_val;
} /* scalxy_ */
