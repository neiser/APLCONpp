#include <math.h>

typedef double doublereal;
typedef int integer;
typedef integer logical;
typedef float real;
#define TRUE_ (1)
#define FALSE_ (0)
#define dabs(x) (doublereal) abs(x)
#define min(a, b) ((a) <= (b) ? (a) : (b))
#define max(a, b) ((a) >= (b) ? (a) : (b))
#define abs(x) ((x) >= 0 ? (x) : -(x))

/* Common Block Declarations */

struct {
  doublereal epsf, epschi, chisq, ftest, ftestp, chsqp, frms, frmsp, derfac,
      decxp, derufc, derlow, weight;
  integer nx, nf, num, iflg, init, ipak, indcf, istat, indst, indlm, ndenda,
      ndende, ndactl, indas, indvs, indtr, indfc, indhh, indxs, inddx, indxp,
      indrh, indwm, india, ndtot, icnt, nxf, mxf, ndf, iunph, ncst, iter,
      ncalls, ndpdim, indqn, itermx, nauxc, indpu, nfprim, ndtotl;
} simcom_;

struct {
  doublereal aux[125000];
} nauxcm_;

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int aplcon_(integer *nvar, integer *mcst) {
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer ij;

  /*     ================================================================== */
  /*     initialize and define dimension  NVAR/MCST */
  /*     set default parameters */
  /*     define pointer to arrays within aux array */
  /*     clear aux array */
  /*     ================================================================== */
  simcom_.nx = *nvar;
  /* number of variables */
  simcom_.nf = *mcst;
  /* number of constraint equations */
  simcom_.nfprim = simcom_.nf;
  /* primary value of NF */
  simcom_.ndpdim = 125000;
  /* dimension of AUX array */
  simcom_.derfac = .001;
  /* derivative factor */
  simcom_.derufc = 1e-5;
  /* factor or unmeasured variable */
  simcom_.derlow = .01;
  /* factor for lower limit */
  simcom_.epsf = 1e-6;
  /* accuracy limit */
  simcom_.epschi = 1e-5;
  /* chi2 accuracy limit */
  simcom_.itermx = 10;
  /* max number of iterations */
  simcom_.nauxc = 125000;
  /* copy AUX dimension */
  simcom_.init = 0;
  simcom_.istat = 0;
  /* init phase */
  simcom_.nxf = simcom_.nx + simcom_.nf;
  /* total number of fit equations */
  simcom_.mxf = (simcom_.nxf * simcom_.nxf + simcom_.nxf) / 2;
  /* ______________________________________________________________________ */
  /*     Indices for steps, flags and limits are defined here */
  /* reserve NX * (NF + 2) for Jacobian */
  /* number elements symmetric matrix */
  simcom_.indst = simcom_.nx * (simcom_.nf + 2);
  /* steps */
  simcom_.indtr = simcom_.indst + simcom_.nx;
  /* transformation flags */
  simcom_.indlm = simcom_.indtr + simcom_.nx;
  /* 2*limits for variables */
  simcom_.ndtot = simcom_.indlm + (simcom_.nx << 1);
  /*     __________________________________________________________________ */
  /*     storage of initial sub-arrays */
  /*           1 ...    NX * (NF+2)   Jacobian, derivative matrix   A(.) */
  /*     INDST+1 ...    NX            steps                         ST(.) */
  /*     INDTR+1 ...    NX            properties of variables */
  /*     INDLM+1 ...    2 * NX        limits for variables          XL(2,.) */
  /*                    ----------- */
  /*     NDTOT =        NX * (NF+6)   initial memory area */
  /* space used so far */
  i__1 = simcom_.ndtot;
  for (ij = 1; ij <= i__1; ++ij) {
    /* NX*(NF+6) */
    nauxcm_.aux[ij - 1] = 0.;
    /* clear  A(.), ST(.),...,XL(2,.) */
  }
  simcom_.ndf = simcom_.nf;
  /* reset n d f */
  simcom_.ncalls = 0;
  /* reset number of calls */
  return 0;
} /* aplcon_ */

/* Subroutine */ int aploop_(doublereal *x, doublereal *vx, doublereal *f,
                             integer *iret) {
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i__, nff;
  extern /* Subroutine */ int asteps_(doublereal *, doublereal *, doublereal *),
      iploop_(doublereal *, doublereal *, doublereal *, doublereal *,
              doublereal *, doublereal *, doublereal *, doublereal *,
              integer *);

  /* steering routine for loop */

  /* Parameter adjustments */
  --f;
  --vx;
  --x;

  /* Function Body */
  if (simcom_.ncalls != 0) {
    goto L10;
  }
  /*     __________________________________________________________________ */
  /*     indices/pointer etc at first APLOOP entry */
  nff = simcom_.nf;
  /* max. number of constraints */
  simcom_.nxf = simcom_.nx + nff;
  /* total number of fit equations */
  simcom_.mxf = (simcom_.nxf * simcom_.nxf + simcom_.nxf) / 2;
  /* number elements symmetric matrix */
  simcom_.indfc = simcom_.indlm + (simcom_.nx << 1);
  /* pointer to FC(NF) = copy of F(NF) */
  simcom_.indcf = simcom_.indfc + nff;
  /* pointer to FCOPY */
  simcom_.indhh = simcom_.indcf + nff;
  /* pointer to HH(NF) = copy of F(.) */
  simcom_.indxs = simcom_.indhh + nff;
  /* save X(.)        pointer */
  simcom_.inddx = simcom_.indxs + simcom_.nx;
  /* step */
  simcom_.indxp = simcom_.inddx + simcom_.nx;
  /* previos step */
  simcom_.indrh = simcom_.indxp + simcom_.nx;
  /* right-hand side */
  simcom_.indwm = simcom_.indrh + simcom_.nxf;
  /* weight matrix */
  simcom_.india = simcom_.indwm + simcom_.mxf;
  /* matrix diagonal "DIAG" */
  simcom_.indqn = simcom_.india + simcom_.nxf;
  /* next pointer    "QNEXT" */
  simcom_.indas = simcom_.indqn + simcom_.nxf;
  /* X result */
  simcom_.indpu = simcom_.indas + simcom_.nx;
  /* pulls, solution X and Vx */
  simcom_.ndtot = simcom_.indpu + simcom_.nx;
  /* total number of words (so far) */
  simcom_.ndtotl = simcom_.ndtot;
  i__1 = simcom_.ndtot;
  for (i__ = simcom_.indlm + (simcom_.nx << 1) + 1; i__ <= i__1; ++i__) {
    nauxcm_.aux[i__ - 1] = 0.;
    /* reset part of aux, unused so far */
  }
  asteps_(&x[1], &vx[1], &nauxcm_.aux[simcom_.indst]);
/*     __________________________________________________________________ */
/*     internal APLOOP */
/* initial steps ST(.) */
L10:
  *iret = -1;
  /* default status is -1 = continue */
  iploop_(&x[1], &vx[1], &f[1], &nauxcm_.aux[simcom_.indxs],
          &nauxcm_.aux[simcom_.inddx], &nauxcm_.aux[simcom_.indcf],
          &nauxcm_.aux[simcom_.indxp], &nauxcm_.aux[simcom_.indrh], iret);
  return 0;
} /* aploop_ */

/* Subroutine */ int iploop_(doublereal *x, doublereal *vx, doublereal *f,
                             doublereal *xs, doublereal *dx, doublereal *fcopy,
                             doublereal *xp, doublereal *rh, integer *iret) {
  /* System generated locals */
  integer i__1;
  doublereal d__1, d__2;

  /* Local variables */
  static integer j;
  static doublereal fj, fex[100];
  static integer nfit, jret;
  extern /* Subroutine */ int anumde_(doublereal *, doublereal *, doublereal *,
                                      doublereal *, doublereal *, doublereal *,
                                      doublereal *, integer *),
      aniter_(doublereal *, doublereal *, doublereal *, doublereal *,
              doublereal *, doublereal *, doublereal *, doublereal *),
      addtox_(doublereal *, doublereal *, doublereal *, doublereal *),
      antest_(integer *),
      acopxv_(doublereal *, doublereal *, doublereal *, doublereal *,
              doublereal *, doublereal *);
  static integer istatu;

  /*     ================================================================== */
  /*     call IPLDER */
  /*     call IPLCON */
  /*     ================================================================== */
  /* steering routine for loop */
  /* ,FOPT,FAC */
  /*     local variables */
  /* Parameter adjustments */
  --rh;
  --xp;
  --fcopy;
  --dx;
  --xs;
  --f;
  --vx;
  --x;

  /* Function Body */
  ++simcom_.ncalls;
  /*     __________________________________________________________________ */
  /*     initialization */
  /* count calls */
  if (simcom_.ncalls != 1) {
    goto L20;
  }
  istatu = 0;
  /* !! */
  nfit = 0;
  /* reset fit count */
  simcom_.iter = 0;
  i__1 = simcom_.nx;
  for (j = 1; j <= i__1; ++j) {
    xs[j] = x[j];
    /* save initial X values */
    dx[j] = 0.;
    /* reset correction DX */
  }
  /*     __________________________________________________________________ */
  /*     start/restart */
  /* L10: */
  ++nfit;
  /* count fits */
  simcom_.iter = 0;
  simcom_.ncst = 0;
  simcom_.chisq = 0.;
L20:
  i__1 = simcom_.nf;
  for (j = 1; j <= i__1; ++j) {
    fj = f[j];
    fex[j - 1] = fj;
    /* extended F */
  }
  /*     __________________________________________________________________ */
  /*     constraint test summary */
  if (istatu < 0) {
    goto L40;
  }
  simcom_.ftestp = simcom_.ftest;
  /* save previous value */
  simcom_.ftest = 0.;
  /* reset constraint tests */
  simcom_.frmsp = simcom_.frms;
  simcom_.frms = 0.;
  i__1 = simcom_.nf;
  for (j = 1; j <= i__1; ++j) {
    fj = f[j];
    fcopy[j] = fj;
    /* copy constraint vector */
    simcom_.ftest += abs(fj);
    /* sum absolute values */
    /* Computing 2nd power */
    d__1 = fj;
    simcom_.frms += d__1 * d__1;
    /* sum squares */
  }
  /* Computing MAX */
  d__1 = 1e-16, d__2 = simcom_.ftest / (real)simcom_.nf;
  simcom_.ftest = max(d__1, d__2);
  /* average |F| */
  simcom_.frms = sqrt(simcom_.frms / (real)simcom_.nf + 1e-32);
  /* LS mean */
  if (istatu == 1) {
    goto L60;
  }
/*     __________________________________________________________________ */
/*     start numerical derivatives */
L30:
  istatu = -1;
/*     __________________________________________________________________ */
/*     derivative calculation */
/* !! */
L40:
  if (istatu + 1 != 0) {
    goto L50;
  }
  /* ISTATU=-1 */
  anumde_(&x[1], fex, nauxcm_.aux, &nauxcm_.aux[simcom_.indst],
          &nauxcm_.aux[simcom_.indlm], &nauxcm_.aux[simcom_.indfc],
          &nauxcm_.aux[simcom_.indhh], &jret);
  /* derivative matrix A */
  /* steps  ST(.) */
  /* limits XL(2,.) */
  /* copy FC(.) central F(.) */
  /* copy HH(.) shifted F(.) */
  *iret = -1;
  if (jret < 0) {
    return 0;
  }
/*     __________________________________________________________________ */
/*     next iteration */
/* ...for constraint calculation */
L50:
  aniter_(&x[1], &vx[1], &fcopy[1], nauxcm_.aux, &xp[1], &rh[1],
          &nauxcm_.aux[simcom_.indwm], &dx[1]);
  goto L70;
/*     __________________________________________________________________ */
/*     test cutsteps */
L60:
  antest_(iret);
  if (*iret + 1 == 0) {
    goto L30;
  }
  /* numerical derivative:   ISTATU=-1 */
  if (*iret >= 0) {
    goto L80;
  }
/*     __________________________________________________________________ */
/*     apply corrections DX(.) to X(.) with transformations */
/* convergence or failure: ISTATU= 2 */
L70:
  addtox_(&x[1], &xs[1], &dx[1], &xp[1]);
  istatu = 1;
  /* test at next entry */
  return 0;
/*     __________________________________________________________________ */
/*     end-of-primary-fit (NFIT=1) */
L80:
  istatu = 2;
  acopxv_(&x[1], &vx[1], &nauxcm_.aux[simcom_.inddx],
          &nauxcm_.aux[simcom_.indas], &nauxcm_.aux[simcom_.indwm],
          &nauxcm_.aux[simcom_.indpu]);
  *iret = 0;
  return 0;
} /* iploop_ */

/* Subroutine */ int asteps_(doublereal *x, doublereal *vx, doublereal *st) {
  /* System generated locals */
  integer i__1, i__2;
  doublereal d__1, d__2, d__3, d__4, d__5;

  /* Local variables */
  static integer i__, j, ii;
  static doublereal vii;
  static integer ntder, ntine, ntlim, ntrfl, ntmes, ntprf, ntvar;
  extern integer ijsym_(integer *, integer *);

  /*     ================================================================== */
  /*     check transformed variables (e.g. > 0) */
  /*     define initial step size for */
  /*        o  measured variables */
  /*        o  unmeasured variables */
  /*     determine number of degrees of freedom */
  /*     transform steps for transformed variables */
  /*     transform covariance matrix for transformed variables */
  /*     ================================================================== */

  /* define initial steps */

  /* Parameter adjustments */
  --st;
  --vx;
  --x;

  /* Function Body */
  ii = 0;
  i__1 = simcom_.nx;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* loop on all variables */
    simcom_.ipak = i__;
    /*     unpackfl.inc = code for flag unpacking */
    ntrfl = (integer)nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1];
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
    ii += i__;
    vii = (d__1 = vx[ii], abs(d__1));
    /* original diagonal element */
    if (vii == 0.) {
      ntmes = 1;
      /* set unmeasured flag */
    }
    /*      _________________________________________________________________ */
    /*      define step size for derivative calculation */
    if (vii != 0.) {
      /* measured variable */
      if (st[i__] > 0.) {
        /* Computing MIN */
        d__1 = st[i__], d__2 = simcom_.derfac * sqrt(vii);
        st[i__] = min(d__1, d__2);
        /* user step, if smaller */
      } else if (st[i__] <= 0.) {
        /* step is undefined */
        st[i__] = simcom_.derfac * sqrt(vii);
        /* step from cov matrix */
        if (ntvar != 2) {
          /* Computing MIN */
          /* Computing MAX */
          d__4 = 1e-6, d__5 = (d__1 = x[i__], abs(d__1));
          d__2 = st[i__], d__3 = simcom_.derlow * max(d__4, d__5);
          st[i__] = min(d__2, d__3);
        } else {
          /* Poisson */
          /* Computing MIN */
          /* Computing MAX */
          d__4 = 1e-6, d__5 = (d__1 = x[i__] + 1., abs(d__1));
          d__2 = st[i__], d__3 = simcom_.derlow * max(d__4, d__5);
          st[i__] = min(d__2, d__3);
        }
      }
    } else if (vii == 0.) {
      /* unmeasured variable */
      --simcom_.ndf;
      /* reduce degrees of freedom */
      i__2 = simcom_.nx;
      for (j = 1; j <= i__2; ++j) {
        vx[ijsym_(&i__, &j)] = 0.;
        /* clear matrix elements */
      }
      /* Computing MAX */
      d__2 = 1., d__3 = (d__1 = x[i__], abs(d__1));
      st[i__] = simcom_.derufc * max(d__2, d__3);
    }
    if (ntine == 1) {
      st[i__] = 0.;
    }
  }
  return 0;
} /* asteps_ */

/* Subroutine */ int anumde_(doublereal *x, doublereal *f, doublereal *a,
                             doublereal *st, doublereal *xl, doublereal *fc,
                             doublereal *hh, integer *jret) {
  /* Initialized data */

  static logical tinue = FALSE_;

  /* System generated locals */
  integer i__1;
  doublereal d__1, d__2;

  /* Local variables */
  static integer i__, j, ij;
  static doublereal xd[2], xt[2], der;
  static integer ilr;
  static doublereal stm;
  static integer nalz, nzer, nonz, ntder, ntine, ntvar, ntmes, ntlim, ntprf;
  static doublereal xsave;
  static logical limdef;
  static doublereal ratdif, ratmax, derzer;

  /*     ================================================================== */
  /*     calculation of numerical derivatives, one variable at a time */
  /*        check limits for variable */
  /*        calculate displaced value of variable */
  /*        calculate derivative in Jacobian matrix A */
  /*        classify derivative properties of variable */
  /*        return with JRET=-1, unless finished */

  /*     X(.),F(.) = variables and constraints */
  /*     A(.)      = Jacobian */
  /*     ST(.)     = steps */
  /*     XL(2,.)   = limits */
  /*     FC(.)     = Constraints at central variable values */
  /*     HH(.)     = Constraints at first step (variable + step) */
  /*     __________________________________________________________________ */
  /*     Logic: */
  /*     Entry  1: start with I=0 ... */
  /*            2: continue with indices */
  /*     Return 1: continue with constraint evaluation */
  /*            2: Jacobian ready */
  /*     ================================================================== */

  /* numerical derivatives */

  /* Parameter adjustments */
  --hh;
  --fc;
  xl -= 3;
  --st;
  --a;
  --f;
  --x;

  /* Function Body */
  /*     ... */
  /* entry flag */
  *jret = -1;
  /* ...means continue at return */
  if (tinue) {
    goto L30;
  }
  /* continue */
  tinue = TRUE_;
  i__ = 0;
/* initialize derivative loop */
L10:
  if (i__ >= simcom_.nx) {
    /* finished */
    *jret = 0;
    /* ... means differentation finished */
    tinue = FALSE_;
    /* Jacobian ready */
    return 0;
  }
  ++i__;
  /* next variable */
  if (st[i__] == 0.f) {
    goto L10;
  }
  /* skip fixed variable */
  /* skip repeated derivative calculation */
  xsave = x[i__];
  /* save current value of variable */
  ilr = 0;
  /*     __________________________________________________________________ */
  /*     check limits for variable */
  /* define steps */
  limdef = xl[(i__ << 1) + 1] != xl[(i__ << 1) + 2];
  /* true if limits defined */
  if (limdef) {
    if (xsave + st[i__] > xl[(i__ << 1) + 2] ||
        xsave - st[i__] > xl[(i__ << 1) + 1]) {
      /* Computing MIN */
      d__1 = xl[(i__ << 1) + 2] - xsave, d__2 = xsave - xl[(i__ << 1) + 1];
      stm = min(d__1, d__2) * .9999f;
      /* minimal step s */
      if (stm * 3.f > st[i__]) {
        st[i__] = stm;
        /* use smaller symmetric step */
      } else {
        /* Computing MAX */
        d__1 = xl[(i__ << 1) + 2] - xsave, d__2 = xsave - xl[(i__ << 1) + 1];
        stm = max(d__1, d__2) * .4999f;
        /* minimal ste */
        if (stm * 2.f < st[i__]) {
          st[i__] = stm;
        }
        if (st[i__] < xl[(i__ << 1) + 2] - xsave) {
          xd[0] = xsave + st[i__];
          xd[1] = xsave + st[i__] * 2.f;
          /* + one-sided steps */
          ilr = 1;
        } else {
          xd[0] = xsave - st[i__];
          xd[1] = xsave - st[i__] * 2.f;
          /* - one-sided steps */
          ilr = 2;
        }
      }
    }
  }
  /*     __________________________________________________________________ */
  /*     define displaced values for derivative calculation */
  if (ilr == 0) {
    /* L20: */
    if (ntvar == 0 || ntvar == 2 || ntvar == 3) {
      xt[0] = xsave + st[i__];
      /* symmetric (two-sided) steps */
      xt[1] = xsave - st[i__];
      xd[0] = xt[0];
      xd[1] = xt[1];
    }
  }
  /*     __________________________________________________________________ */
  /*     set variable to displaced value and return for calculation */
  x[i__] = xd[0];
  /* first step */
  i__ = -i__;
  return 0;
/*     __________________________________________________________________ */
/*     continue */
L30:
  if (i__ < 0) {
    /* calculation of first step done ... */
    i__1 = simcom_.nf;
    for (j = 1; j <= i__1; ++j) {
      hh[j] = f[j];
      /* save constraint values */
    }
    i__ = -i__;
    /* reverse flag */
    x[i__] = xd[1];
    /* set next step ... */
    return 0;
    /* ... and return for second step */
  }
  /*     __________________________________________________________________ */
  /*     INIT ne 0: second step done - calculate derivative */
  x[i__] = xsave;
  /* restore variable I */
  ij = i__;
  /* derivative calculation */
  nzer = 0;
  /* reset: number of zero derivatives */
  nonz = 0;
  /* reset: number of non-zero derivatives */
  nalz = 0;
  /* flag all derivatives are zero */
  ratmax = 0.;
  /* max of diff-ratio */
  derzer = 0.;
  /* abs of nonzero-derivative */
  i__1 = simcom_.nf;
  for (j = 1; j <= i__1; ++j) {
    /* loop on all constraint functions */
    if (ilr == 0) {
      /* symmetric formula */
      der = (hh[j] - f[j]) / (xt[0] - xt[1]);
      /* !! internal variable */
    } else {
      /* asymmetric formula */
      der = (fc[j] * 3. + f[j] - hh[j] * 4.) * .5 / st[i__];
      if (ilr == 2) {
        der = -der;
      }
      /* sign */
    }
    /*      _________________________________________________________________ */
    /*      classify derivative properties of variable I */
    if (a[ij] != 0. || der != 0.) {
      nalz = 1;
    }
    /* non all zero */
    if (der == 0.) {
      ++nzer;
      /* derivative zero - count */
    } else {
      /* derivative non-zero */
      ratdif = (d__2 = a[ij] - der, abs(d__2)) /
               ((d__1 = a[ij], abs(d__1)) + abs(der));
      ratmax = max(ratmax, ratdif);
      derzer = abs(der);
      /* abs value of this derivative */
      ++nonz;
      /* count non-zero derivative */
    }
    a[ij] = der;
    /* insert into Jacobian matrix A */
    ij += simcom_.nx;
  }
  if (ntder == 0) {
    ntder = 1;
  } else {
    if (nonz == 1 && (d__1 = derzer - 1., abs(d__1)) < 1e-12) {
      /* Computing MIN */
      i__1 = ntder + 1;
      ntder = min(i__1, 7);
    } else if (ratmax < 1e-12) {
      /* Computing MIN */
      i__1 = ntder + 1;
      ntder = min(i__1, 7);
    } else {
      ntder = 2;
      /* reset to 2 */
    }
  }
  ntder = 0;
  /*     packfl.inc   = code for flag packing */
  /*     explanation see: */
  /*     unpackfl.inc = code for flag unpacking */
  nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1] = (doublereal)(
      ((((ntprf * 10 + ntlim) * 10 + ntine) * 10 + ntder) * 10 + ntmes) * 10 +
      ntvar);
  goto L10;
} /* anumde_ */

/* Subroutine */ int aniter_(doublereal *x, doublereal *vx, doublereal *f,
                             doublereal *a, doublereal *xp, doublereal *rh,
                             doublereal *wm, doublereal *dx) {
  /* System generated locals */
  integer i__1;
  doublereal d__1;

  /* Local variables */
  static integer i__, j, ia, ii;
  static doublereal diag[1000];
  static integer nrank, ntrfl, ntvar;
  static doublereal qnext[1000];
  extern /* Subroutine */ int duminv_(doublereal *, doublereal *, doublereal *,
                                      integer *, integer *, integer *,
                                      integer *, doublereal *, doublereal *);
  extern doublereal scalxy_(doublereal *, doublereal *, integer *);

  /* next iteration step */

  /* Parameter adjustments */
  --dx;
  --wm;
  --rh;
  --xp;
  --a;
  --f;
  --vx;
  --x;

  /* Function Body */
  ++simcom_.iter;
  /* start next iteration */
  simcom_.chsqp = simcom_.chisq;
  /*     __________________________________________________________________ */
  /*     right-hand side of equation */
  /* save current chi^2 */
  simcom_.ncst = 0;
  i__1 = simcom_.nx;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* first NX components */
    rh[i__] = 0.;
    /* define right hand side of equation */
  }
  i__1 = simcom_.nf;
  for (j = 1; j <= i__1; ++j) {
    /* next NF components */
    nauxcm_.aux[simcom_.indfc + j - 1] = f[j];
    rh[simcom_.nx + j] = -f[j];
  }
  ia = 0;
  i__1 = simcom_.nf;
  for (j = 1; j <= i__1; ++j) {
    /* "subtract" actual step */
    rh[simcom_.nx + j] += scalxy_(&a[ia + 1], &dx[1], &simcom_.nx);
    ia += simcom_.nx;
    nauxcm_.aux[simcom_.indhh + j - 1] = rh[simcom_.nx + j];
    /* right hand side for chi**2 */
  }
  /*     __________________________________________________________________ */
  /*     form matrix and solve */
  i__1 = (simcom_.nx * simcom_.nx + simcom_.nx) / 2;
  for (i__ = 1; i__ <= i__1; ++i__) {
    wm[i__] = -vx[i__];
    /* copy -VX(.) into W_11 */
  }
  ii = 0;
  /* modify V for Poisson variables */
  i__1 = simcom_.nx;
  for (i__ = 1; i__ <= i__1; ++i__) {
    ii += i__;
    simcom_.ipak = i__;
    /*     unpackfl.inc = code for flag unpacking */
    ntrfl = (integer)nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1];
    /* get packed flags */
    ntvar = ntrfl % 10;
    /* transformation flag */
    /* profile flag */
    if (ntvar == 2) {
      /* Poisson */
      /* Computing 2nd power */
      d__1 = x[i__];
      wm[ii] = -sqrt(d__1 * d__1 + 1.f);
      /* -MAX(ABS(X(I)),1.0D0) */
    }
  }
  duminv_(&a[1], &wm[1], &rh[1], &simcom_.nx, &simcom_.nf, &c__1, &nrank,
          diag, qnext);
  simcom_.chisq = -scalxy_(&nauxcm_.aux[simcom_.indhh], &rh[simcom_.nx + 1],
                            &simcom_.nf);
  /* next chi^2 */
  if (simcom_.chisq < 0.) {
    simcom_.chisq = 0.;
  }
  /*     __________________________________________________________________ */
  /*     handle corrections and cutstep */
  simcom_.weight = 1.;
  /* default weight */
  if (simcom_.iter > 1 && simcom_.chisq >= simcom_.chsqp * 2.) {
    simcom_.weight = .1;
  }
  if (simcom_.iter > 1 && simcom_.chisq >= simcom_.chsqp * 3.) {
    simcom_.weight = .05;
  }
  i__1 = simcom_.nx;
  for (i__ = 1; i__ <= i__1; ++i__) {
    xp[i__] = dx[i__];
    /* save previous corrections */
    dx[i__] = rh[i__];
    /* store new corrections */
  }
  return 0;
} /* aniter_ */

/* Subroutine */ int addtox_(doublereal *x, doublereal *xs, doublereal *dx,
                             doublereal *xp) {
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i__;

  /* Parameter adjustments */
  --xp;
  --dx;
  --xs;
  --x;

  /* Function Body */
  i__1 = simcom_.nx;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dx[i__] = simcom_.weight * dx[i__] + (1. - simcom_.weight) * xp[i__];
    /* reduce step evtl. */
    x[i__] = xs[i__] + dx[i__];
    /* correct x and return to test constraints */
  }
  return 0;
} /* addtox_ */

/* Subroutine */ int antest_(integer *iret) {
  /* System generated locals */
  doublereal d__1, d__2;

  /* Local variables */
  static integer i__;
  static doublereal cm[14];
  extern /* Subroutine */ int lesfcm_(doublereal *);
  static real factor, dchisq;

  *iret = -1;
  /* calculate new Jacobian */
  simcom_.iunph = 0;
  /*     __________________________________________________________________ */
  /*     combined measure? */
  if (simcom_.iter == 1 && simcom_.ncst == 0) {
    factor =
        simcom_.chisq /
        ((d__1 = simcom_.ftestp - simcom_.ftest, abs(d__1)) + simcom_.epsf);
    for (i__ = 1; i__ <= 14; ++i__) {
      cm[i__ - 1] = 0.;
    }
    cm[8] = .5;
    /* damping factor */
    cm[0] = 2.;
    /* Computing 2nd power */
    d__1 = simcom_.frms;
    /* Computing 2nd power */
    d__2 = simcom_.frmsp;
    cm[1] = d__1 * d__1 + d__2 * d__2;
    /* Computing 4th power */
    d__1 = simcom_.frms, d__1 *= d__1;
    /* Computing 4th power */
    d__2 = simcom_.frmsp, d__2 *= d__2;
    cm[2] = d__1 * d__1 + d__2 * d__2;
    cm[3] = simcom_.chisq + simcom_.chsqp;
    /* Computing 2nd power */
    d__1 = simcom_.frms;
    /* Computing 2nd power */
    d__2 = simcom_.frmsp;
    cm[4] = simcom_.chisq * (d__1 * d__1) + simcom_.chsqp * (d__2 * d__2);
    lesfcm_(cm);
    /* fit */
    /* Computing 2nd power */
    d__1 = simcom_.frms;
    cm[12] = cm[5] - cm[6] * (d__1 * d__1);
    /* combined penalty */
  } else if (simcom_.ncst == 0) {
    cm[0] = cm[8] * cm[0] + 1.;
    /* Computing 2nd power */
    d__1 = simcom_.frms;
    cm[1] = cm[8] * cm[1] + d__1 * d__1;
    /* Computing 4th power */
    d__1 = simcom_.frms, d__1 *= d__1;
    cm[2] = cm[8] * cm[2] + d__1 * d__1;
    cm[3] = cm[8] * cm[3] + simcom_.chisq;
    /* Computing 2nd power */
    d__1 = simcom_.frms;
    cm[4] = cm[8] * cm[4] + simcom_.chisq * (d__1 * d__1);
    lesfcm_(cm);
    /* fit */
    cm[13] = cm[12];
    /* Computing 2nd power */
    d__1 = simcom_.frms;
    cm[12] = cm[5] - cm[6] * (d__1 * d__1);
    /* combined penalty */
  }
  /*     __________________________________________________________________ */
  /*     cutstep */
  if (simcom_.ncst < 2 &&
      (simcom_.iunph != 0 ||
       simcom_.iter > 1 &&
           simcom_.ftest > simcom_.ftestp * 2. + simcom_.epsf)) {
    ++simcom_.ncst;
    simcom_.weight = .25;
    simcom_.weight = .5;
    *iret = -2;
    /* cutstep - add corrections */
    return 0;
  }
  /*     __________________________________________________________________ */
  /*     convergent */
  if (simcom_.iter >= 2 && simcom_.ncst == 0) {
    dchisq = simcom_.chisq - simcom_.chsqp;
    if (dabs(dchisq) <= simcom_.epschi && simcom_.ftest < simcom_.epsf) {
      *iret = 0;
      /* convergence */
      return 0;
    }
  }
  /*     __________________________________________________________________ */
  /*     failure */
  if (simcom_.iter > simcom_.itermx) {
    *iret = 2;
  }
  /* non-convergence */
  return 0;
} /* antest_ */

/* Subroutine */ int acopxv_(doublereal *x, doublereal *vx, doublereal *dx,
                             doublereal *as, doublereal *wm, doublereal *pu) {
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i__, ii;
  static doublereal scopy;

  /*     __________________________________________________________________ */
  /*     convergence: pull calculation */
  /* Parameter adjustments */
  --pu;
  --wm;
  --as;
  --dx;
  --vx;
  --x;

  /* Function Body */
  ii = 0;
  i__1 = simcom_.nx;
  for (i__ = 1; i__ <= i__1; ++i__) {
    as[i__] = x[i__];
    ii += i__;
    pu[i__] = 0.f;
    if (vx[ii] > 0.f) {
      if (vx[ii] - wm[ii] > 0.f) {
        pu[i__] = dx[i__] / sqrt(vx[ii] - wm[ii]);
      }
    }
  }
  /*     __________________________________________________________________ */
  /*     copy/exchange result/input covariance matrix */
  i__1 = (simcom_.nx * simcom_.nx + simcom_.nx) / 2;
  for (i__ = 1; i__ <= i__1; ++i__) {
    scopy = vx[i__];
    vx[i__] = wm[i__];
    /* copy fitted covariance matrix */
    wm[i__] = scopy;
    /* ... and save input matrix */
  }
  return 0;
} /* acopxv_ */

/* Subroutine */ int lesfcm_(doublereal *c__) {
  /*     ... */
  /* Parameter adjustments */
  --c__;

  /* Function Body */
  c__[8] = c__[1] * c__[3] - c__[2] * c__[2];
  /* determinant */
  c__[10] = c__[3] / c__[8];
  /* V_11 */
  c__[11] = -c__[2] / c__[8];
  /* V_12 */
  c__[12] = c__[1] / c__[8];
  /* V_22 */
  c__[6] = c__[4] * c__[10] + c__[5] * c__[11];
  /* 1. parameter */
  c__[7] = c__[4] * c__[11] + c__[5] * c__[12];
  /* 2. parameter */
  return 0;
} /* lesfcm_ */

/* Subroutine */ int chndpv_(real *chi2, integer *nd, real *pval) {
  extern doublereal chprob_(doublereal *, integer *);

  *chi2 = simcom_.chisq;
  /* chi^square */
  *nd = simcom_.ndf;
  /* number of degrees of freedom */
  *pval = chprob_(&simcom_.chisq, &simcom_.ndf);
  /* p-value */
  return 0;
} /* chndpv_ */
