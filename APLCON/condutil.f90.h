#include <cmath>
#include <algorithm>

#include "helper.h"

using namespace std;

inline int ijsym_(int *i, int *j) {

    /* index (I,J)=(J,I) in symmetric matrix */
    if (*i <= *j) {
        return (*j * *j - *j) / 2 + *i;
    } else {
        return (*i * *i - *i) / 2 + *j;
    }
} /* ijsym_ */

/* Subroutine */
inline void duminv_(const vecd& a, vecd& w, vecd& b, int nx, int nf) {
    /* Initialized data */

    static double eps = 1e-6;

    /* Local variables */
    static int i, j, k, l, m, n, ia, ij, jk, jj, kk, jl, lk;
    static double vjk, vkk, sum;
    static int last, nmeas, jlast;
    static int jfirst;

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

    /* Function Body */
    /*     special entry for partial inversion ****************************** */
    /*     ... */
    n = nx + nf;

    vecd aux;
    aux.resize(n);

    vec<int> qnext;
    qnext.resize(n);


    /*     -VX(NX-sym) is already inserted in W(NX+NF-sym) ------------------ */
    ij = (nx * nx + nx) / 2;
    /* number of elements of V(.) */
    ia = 0;
    for (j = 1; j <= nf; ++j) {
        for (i = 1; i <= nx; ++i) {
            w[ij + i] = a[ia + i];
            /* copy A(.) into W_12 */
        }
        for (i = 1; i <= j; ++i) {
            w[ij + nx + i] = 0.f;
            /* reset last submatrix W_22 of W(.) */
        }
        ij = ij + nx + j;
        ia += nx;
    }
    /*     distinguish between measured and unmeasured variables ------------ */
    jfirst = 0;
    /* first index of measured variable */
    nmeas = 0;
    /* number of measured variables */
    for (i = 1; i <= nx; ++i) {
        if (w[(i * i + i) / 2] < 0.f) {
            /* measured variable */
            if (jfirst == 0) {
                jfirst = i;
                /* first index of measured variable */
            } else {
                qnext[jlast] = i;
                /* insert index at previous index */
            }
            jlast = i;
            /* save index */
            ++nmeas;
        }
    }
    if (jlast == 0) {
        goto L10;
    }
    /* nothing to do */
    qnext[jlast] = -1;
    /*     apply exchange algorithm to sub-matrices ------------------------- */
    /* stop index for last measured variable */
    for (i = nx + 1; i <= n; ++i) {
        /* loop I over constraint equations */
        j = jfirst;
        /* first index of unmeasured variable */
        for (m = 1; m <= nmeas; ++m) {
            /* already inverted element index J */
            sum = 0.;
            jk = (j * j - j) / 2;
            /* index of diagonal element before */
            for (k = 1; k <= nx; ++k) {
                if (k <= j) {
                    ++jk;
                }
                /* index in j column */
                if (qnext[k] != 0) {
                    sum += w[jk] * w[(i * i - i) / 2 + k];
                }
                if (k >= j) {
                    jk += k;
                }
                /* index in j row */
            }
            aux[j] = sum;
            /* = A-row * VX-row/col */
            j = qnext[j];
            /* next index of unmeasured variable */
        }
        for (k = i; k <= n; ++k) {
            sum = 0.;
            j = jfirst;
            /* first index of unmeasured variable */
            for (m = 1; m <= nmeas; ++m) {
                /* already inverted element index J */
                sum += w[(k * k - k) / 2 + j] * aux[j];
                /* = A-row * H */
                j = qnext[j];
                /* next index of unmeasured variable */
            }
            w[(k * k - k) / 2 + i] += sum;
            /* add to diagonal W_22 */
        }
        j = jfirst;
        /* first index of unmeasured variable */
        for (m = 1; m <= nmeas; ++m) {
            w[(i * i - i) / 2 + j] = -aux[j];
            /* add to off-diagonal W_22 */
            j = qnext[j];
            /* next index of unmeasured variable */
        }
    }
    /*     set pointer for unmeasured variables ----------------------------- */
    jfirst = 0;
    jlast = 0;
    for (i = 1; i <= n; ++i) {
        if (qnext[i] == 0) {
            /* unmeasured variable */
            if (jfirst == 0) {
                jfirst = i;
                /* first index of unmeasured variable */
            } else {
                qnext[jlast] = i;
                /* next index of unmeasured variable */
            }
            jlast = i;
        } else {
            qnext[i] = 0;
            /* reset index for measured variable */
        }
    }
    if (jlast == 0) {
        goto L10;
    }
    /* no unmeasured variable */
    qnext[jlast] = -1;
    /*     common code for inversion and (M=1) solution of matrix equation */
    /* end flag */
L10:
    /*     loop begin (loop on all remaining rows/cols) */
    /* solution flag */
    for (i = 1; i <= n; ++i) {
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
            if (abs(w[jj]) > max(abs(vkk), eps * aux[j])) {
                vkk = w[jj];
                /* largest pivot candidate so far */
                k = j;
                /* index of largest */
                l = last;
            }
            last = j;
            j = qnext[j];
            /* index of next candidate */
            goto L20;
        }
        if (k != 0) {
            /* pivot element found - proceed */
            kk = (k * k + k) / 2;
            if (l == 0) {
                jfirst = (int)qnext[k];
                /* new first index */
            } else {
                qnext[l] = qnext[k];
                /* bridge used index */
            }
            qnext[k] = 0;
            /* reset used index */
            vkk = 1.f / vkk;
            /* invert pivot */
            w[kk] = -vkk;
            b[k] *= vkk;
            jk = kk - k;
            jl = 0;
            for (j = 1; j <= n; ++j) {
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
                    b[j] -= b[k] * vjk;
                    lk = kk - k;
                    for (l = 1; l <= j; ++l) {
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
            for (k = 1; k <= n; ++k) {
                if (qnext[k] != 0) {
                    /* undefined variable */
                    b[k] = 0;
                    /* clear undefined vector element */
                    for (j = 1; j <= k; ++j) {
                        if (qnext[j] != 0) {
                            w[(k * k - k) / 2 + j] = 0.;
                        }
                        /* clear matrix */
                    }
                }
            }
            goto L30;
        }
    }
    /* end of inversion loop */
L30:
    for (i = 1; i <= (n * n + n) / 2; ++i) {
        w[i] = -w[i];
        /* finally reverse sign */
    }

} /* duminv_ */

inline double scalxy_(double *x, double *y, int *n) {
    /* Local variables */
    static int j;
    static double sum;

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
    for (j = 1; j <= *n; ++j) {
        sum += x[j] * y[j];
    }
    return sum;
} /* scalxy_ */
