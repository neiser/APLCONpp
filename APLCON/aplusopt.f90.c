typedef double doublereal;
typedef int integer;
typedef float real;
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define abs(x) ((x) >= 0 ? (x) : -(x))

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

/* Subroutine */ int adummy_0_(int n__, integer *lunp, integer *jpr,
                               doublereal *arg, integer *it, integer *i__,
                               doublereal *step, doublereal *xlow,
                               doublereal *xhig, integer *nbinom,
                               doublereal *pow) {
  static integer ntder, ntine, ntrfl, ntvar, ntmes, ntlim, ntprf;
  /*     parameters for fit method */
  /*     __________________________________________________________________ */
  switch (n__) {
  case 2:
    goto L_apdeps;
  case 8:
    goto L_apstep;
  case 12:
    goto L_apoiss;
  }

L_apdeps:
  /* constraint accuracy */
  simcom_1.epsf = *arg;
  /* |F| accuracy */
  return 0;
  /*     __________________________________________________________________ */

L_apstep:
  /* step size for numdif */
  if (*i__ < 1 || *i__ > simcom_1.nx) {
    return 0;
  }
  simcom_1.ipak = *i__;
  /*     unpackfl.inc = code for flag unpacking */
  ntrfl = (integer)nauxcm_1.aux[simcom_1.indtr + simcom_1.ipak - 1];
  /* get packed flags */
  ntvar = ntrfl % 10;
  /* transformation flag */
  ntmes = ntrfl / 10 % 10;
  /* M-estimate flag */
  ntder = ntrfl / 100 % 10;
  /* derivative type flag */
  ntine = ntrfl / 1000 % 10;
  /* inequality flag */
  ntlim = ntrfl / 10000 % 10;
  /* limit flag */
  ntprf = ntrfl / 100000 % 10;
  /* profile flag */
  nauxcm_1.aux[simcom_1.indst + *i__ - 1] = abs(*step);
  /* ST(I)= ... */
  if (*step != 0.) {
    ntine = 0;
    /* variable */
  } else {
    ntine = 1;
    /* fixed by user */
  }
  goto L100;
  /*     __________________________________________________________________ */

L_apoiss:
  /* Poisson distributed variable */
  if (*i__ < 1 || *i__ > simcom_1.nx) {
    return 0;
  }
  simcom_1.ipak = *i__;
  /*     unpackfl.inc = code for flag unpacking */
  ntrfl = (integer)nauxcm_1.aux[simcom_1.indtr + simcom_1.ipak - 1];
  /* get packed flags */
  ntvar = ntrfl % 10;
  /* transformation flag */
  ntmes = ntrfl / 10 % 10;
  /* M-estimate flag */
  ntder = ntrfl / 100 % 10;
  /* derivative type flag */
  ntine = ntrfl / 1000 % 10;
  /* inequality flag */
  ntlim = ntrfl / 10000 % 10;
  /* limit flag */
  ntprf = ntrfl / 100000 % 10;
  /* profile flag */
  ntvar = 2;
  /* Poisson distributed variable */
  ntlim = 1;
  goto L100;
  /*     __________________________________________________________________ */


/*     __________________________________________________________________ */
L100:
  /*     packfl.inc   = code for flag packing */
  /*     explanation see: */
  /*     unpackfl.inc = code for flag unpacking */
  nauxcm_1.aux[simcom_1.indtr + simcom_1.ipak - 1] = (doublereal)(
      ((((ntprf * 10 + ntlim) * 10 + ntine) * 10 + ntder) * 10 + ntmes) * 10 +
      ntvar);
  return 0;
} /* adummy_ */

/* Subroutine */ int apdeps_(doublereal *arg) {
  return adummy_0_(2, (integer *)0, (integer *)0, arg, (integer *)0,
                   (integer *)0, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apstep_(integer *i__, doublereal *step) {
  return adummy_0_(8, (integer *)0, (integer *)0, (doublereal *)0, (integer *)0,
                   i__, step, (doublereal *)0, (doublereal *)0, (integer *)0,
                   (doublereal *)0);
}

/* Subroutine */ int apoiss_(integer *i__) {
  return adummy_0_(12, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}


/* Subroutine */ int apstat_(doublereal *fopt, integer *nfun, integer *niter) {
  /*     __________________________________________________________________ */
  /*     return information after the fit */
  /*     __________________________________________________________________ */
  /* return Fopt and Nfun */
  *fopt = simcom_1.chisq;
  *nfun = simcom_1.ncalls;
  *niter = simcom_1.iter;
  return 0;
} /* apstat_ */
