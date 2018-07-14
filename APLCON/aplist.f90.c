/* aplist.f90.F -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h> /* For exit() */
#include <f2c.h>

/* Common Block Declarations */

struct {
    doublereal epsf, epschi, chisq, ftest, ftestp, chsqp, frms, frmsp, derfac,
	     decxp, derufc, derlow, weight;
    integer nx, nf, num, iflg, init, ipak, indcf, istat, indst, indlm, ndenda,
	     ndende, ndactl, indas, indvs, indtr, indfc, indhh, indxs, inddx, 
	    indxp, indrh, indwm, india, ndtot, icnt, nxf, mxf, ndf, iunph, 
	    ncst, iter, ncalls, ndpdim, indqn, itermx, nauxc, indpu, nfprim, 
	    ndtotl;
    real tab[10000]	/* was [1000][10] */;
} simcom_;

#define simcom_1 simcom_

struct {
    doublereal aux[125000];
} nauxcm_;

#define nauxcm_1 nauxcm_

/* Subroutine */ int appull_(doublereal *pulls)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*     =========================== MACRO ================================ */
/*     Parameter statement for basic dimensions */
/*     Definition of common for APLCON/ERRPRP/SIM... subroutines */

/*     NX       = number of parameters */
/*     NF       = number of constraint equations */
/*     NUM      = flag for numerical differentiation */
/*                = 0  analytical derivatives */
/*                = 1  numerical derivatives */
/*                = 2  numerical derivatives, compare with analytical */
/*     IFLG     = flag for first case (check derivative) */
/*     INIT     = flag */
/*                = 1  deriving */
/*                = 0  else */

/*     EPSF     = required |F| accuracy */
/*     ISTAT    = status flag */
/*                =0 */
/*                =1  derivative matrix finished */
/*                =2  X(.) corrections applied */

/*     XL(2,.)  = lower and upper values of parameters */
/*     ST(.)    = step sizes for numerical differentiation */
/*     FC(.)    = central values of parameters */
/*     H(.)     = copy */

/*     A        = derivative matrix a/flags during matrix inversion/pulls */
/*     A(NX,NF) */

/*     ****************************************************************** */

/*     ND       = number of degrees of freedom */
/*              = number of constraint equations minus number of */
/*                unmeasured parameters */
/*     CHSQ     = chi square */

/*     DERFAC   = factor for standard deviation in numerical derivatives */

/*     NDENDE   = index of last used word incl. single-precision array */
/*     NDACTL   = index of actual single-precision array */
/*     =========================end=of=macro============================= */
/*     Auxiliary array for all arrays used in APLCON */
/* 1 Mega byte */
/*     ... */
    /* Parameter adjustments */
    --pulls;

    /* Function Body */
    i__1 = simcom_1.nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pulls[i__] = nauxcm_1.aux[simcom_1.indpu + i__ - 1];
    }
    return 0;
} /* appull_ */

