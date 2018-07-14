typedef double doublereal;
typedef int integer;
typedef float real;
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


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

  /*     __________________________________________________________________ */
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
  /*     flags used in packfl.inc, unpackfl.inc */
  /*     __________________________________________________________________ */

  /*     parameters for fit method */
  /*     __________________________________________________________________ */
  /*     Do nothing...no printing anymore */
  switch (n__) {
  case 1:
    goto L_aprint;
  case 2:
    goto L_apdeps;
  case 3:
    goto L_apepschi;
  case 4:
    goto L_apderf;
  case 5:
    goto L_apderu;
  case 6:
    goto L_apdlow;
  case 7:
    goto L_apiter;
  case 8:
    goto L_apstep;
  case 9:
    goto L_apfix;
  case 10:
    goto L_aplimt;
  case 11:
    goto L_aptrin;
  case 12:
    goto L_apoiss;
  case 13:
    goto L_abinom;
  case 14:
    goto L_aplogn;
  case 15:
    goto L_apsqrt;
  case 16:
    goto L_apower;
  case 17:
    goto L_aposit;
  }

L_aprint:
  return 0;
  /*     __________________________________________________________________ */

L_apdeps:
  /* constraint accuracy */
  simcom_1.epsf = *arg;
  /* |F| accuracy */
  return 0;
  /*     __________________________________________________________________ */

L_apepschi:
  /* chi2 accuracy */
  simcom_1.epschi = *arg;
  return 0;
  /*     __________________________________________________________________ */

L_apderf:
  /* factor for step definition */
  simcom_1.derfac = *arg;
  return 0;
  /*     __________________________________________________________________ */

L_apderu:
  /* factor for step definition */
  simcom_1.derufc = *arg;
  return 0;
  /*     __________________________________________________________________ */

L_apdlow:
  /* factor for step definition */
  simcom_1.derlow = *arg;
  return 0;
  /*     __________________________________________________________________ */

L_apiter:
  /* iteration limit */
  simcom_1.itermx = max(3, *it);
  /* max number of iterations */
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

L_apfix:
  /* fixed parameter */
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
  ntine = 1;
  /* fixed by user */
  goto L100;
  /*     __________________________________________________________________ */

L_aplimt:
  /* range of variable */
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
  nauxcm_1.aux[simcom_1.indlm + (*i__ - 1 << 1)] = min(*xlow, *xhig);
  /* lower limit XL(1,I) */
  nauxcm_1.aux[simcom_1.indlm + (*i__ - 1 << 1) + 1] = max(*xlow, *xhig);
  /* upper limit XL(2,I) */
  ntlim = 4;
  goto L100;
  /*     __________________________________________________________________ */

L_aptrin:
  /* inverse value */
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
  ntvar = 1;
  /* transformation to inverse */
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
  /*      WRITE(*,*) 'APOISS I,IPAK,NTVAR,NTLIM ',I,IPAK,NTVAR,NTLIM */
  goto L100;
  /*     __________________________________________________________________ */

L_abinom:
  /* Binomial distributed variable */
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
  ntvar = 3;
  /* Binomial distributed variable */
  nauxcm_1.aux[simcom_1.indlm + (*i__ << 1) - 1] = (doublereal)(*nbinom);
  ntlim = 1;
  goto L100;
  /*     __________________________________________________________________ */

L_aplogn:
  /* Lognormal distributed variable */
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
  ntvar = 4;
  /* Lognormal distributed variable */
  ntlim = 1;
  goto L100;
  /*     __________________________________________________________________ */

L_apsqrt:
  /* SQRT transformation */
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
  ntvar = 5;
  /* SQRT transformation */
  ntlim = 1;
  goto L100;
  /*     __________________________________________________________________ */

L_apower:
  /* x^power transformation */
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
  ntvar = 6;
  /* x^power transformation */
  nauxcm_1.aux[simcom_1.indlm + (*i__ << 1) - 1] = *pow;
  ntlim = 1;
  goto L100;
  /*     __________________________________________________________________ */

L_aposit:
  /* positive */
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
  ntlim = 1;
  /* positive */
  goto L100;
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

/* Subroutine */ int adummy_(void) {
  return adummy_0_(0, (integer *)0, (integer *)0, (doublereal *)0, (integer *)0,
                   (integer *)0, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int aprint_(integer *lunp, integer *jpr) {
  return adummy_0_(1, lunp, jpr, (doublereal *)0, (integer *)0, (integer *)0,
                   (doublereal *)0, (doublereal *)0, (doublereal *)0,
                   (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apdeps_(doublereal *arg) {
  return adummy_0_(2, (integer *)0, (integer *)0, arg, (integer *)0,
                   (integer *)0, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apepschi_(doublereal *arg) {
  return adummy_0_(3, (integer *)0, (integer *)0, arg, (integer *)0,
                   (integer *)0, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apderf_(doublereal *arg) {
  return adummy_0_(4, (integer *)0, (integer *)0, arg, (integer *)0,
                   (integer *)0, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apderu_(doublereal *arg) {
  return adummy_0_(5, (integer *)0, (integer *)0, arg, (integer *)0,
                   (integer *)0, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apdlow_(doublereal *arg) {
  return adummy_0_(6, (integer *)0, (integer *)0, arg, (integer *)0,
                   (integer *)0, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apiter_(integer *it) {
  return adummy_0_(7, (integer *)0, (integer *)0, (doublereal *)0, it,
                   (integer *)0, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apstep_(integer *i__, doublereal *step) {
  return adummy_0_(8, (integer *)0, (integer *)0, (doublereal *)0, (integer *)0,
                   i__, step, (doublereal *)0, (doublereal *)0, (integer *)0,
                   (doublereal *)0);
}

/* Subroutine */ int apfix_(integer *i__) {
  return adummy_0_(9, (integer *)0, (integer *)0, (doublereal *)0, (integer *)0,
                   i__, (doublereal *)0, (doublereal *)0, (doublereal *)0,
                   (integer *)0, (doublereal *)0);
}

/* Subroutine */ int aplimt_(integer *i__, doublereal *xlow, doublereal *xhig) {
  return adummy_0_(10, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, xlow, xhig, (integer *)0,
                   (doublereal *)0);
}

/* Subroutine */ int aptrin_(integer *i__) {
  return adummy_0_(11, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apoiss_(integer *i__) {
  return adummy_0_(12, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int abinom_(integer *i__, integer *nbinom) {
  return adummy_0_(13, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, nbinom, (doublereal *)0);
}

/* Subroutine */ int aplogn_(integer *i__) {
  return adummy_0_(14, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apsqrt_(integer *i__) {
  return adummy_0_(15, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apower_(integer *i__, doublereal *pow) {
  return adummy_0_(16, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, pow);
}

/* Subroutine */ int aposit_(integer *i__) {
  return adummy_0_(17, (integer *)0, (integer *)0, (doublereal *)0,
                   (integer *)0, i__, (doublereal *)0, (doublereal *)0,
                   (doublereal *)0, (integer *)0, (doublereal *)0);
}

/* Subroutine */ int apstat_(doublereal *fopt, integer *nfun, integer *niter) {
  /*     __________________________________________________________________ */
  /*     return information after the fit */
  /*     __________________________________________________________________ */
  /* return Fopt and Nfun */
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
  /*     ... */
  *fopt = simcom_1.chisq;
  *nfun = simcom_1.ncalls;
  *niter = simcom_1.iter;
  return 0;
} /* apstat_ */
