typedef double doublereal;
typedef int integer;
typedef float real;

/* Common Block Declarations */

struct {
  doublereal epsf, epschi, chisq, ftest, ftestp, chsqp, frms, frmsp, derfac,
      decxp, derufc, derlow, weight;
  integer nx, nf, num, iflg, init, ipak, indcf, istat, indst, indlm, ndenda,
      ndende, ndactl, indas, indvs, indtr, indfc, indhh, indxs, inddx, indxp,
      indrh, indwm, india, ndtot, icnt, nxf, mxf, ndf, iunph, ncst, iter,
      ncalls, ndpdim, indqn, itermx, nauxc, indpu, nfprim, ndtotl;
  real tab[10000] /* was [1000][10] */;
} simcom_;

#define simcom_1 simcom_

struct {
  doublereal aux[125000];
} nauxcm_;

#define nauxcm_1 nauxcm_

/* Subroutine */ int appull_(doublereal *pulls) {
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i__;

  /* Parameter adjustments */
  --pulls;

  /* Function Body */
  i__1 = simcom_1.nx;
  for (i__ = 1; i__ <= i__1; ++i__) {
    pulls[i__] = nauxcm_1.aux[simcom_1.indpu + i__ - 1];
  }
  return 0;
} /* appull_ */
