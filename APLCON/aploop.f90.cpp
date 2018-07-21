#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <sstream>

typedef double doublereal;
typedef int integer;
typedef bool logical;
typedef float real;

using namespace std;

/* Common Block Declarations */

struct {
    doublereal epsf, epschi, chisq, ftest, ftestp, chsqp, frms, frmsp, derfac,
    derufc, derlow, weight;
    integer nx, nf, ipak, indcf, indst, indlm,
    indas, indtr, indfc, indhh, indxs, inddx, indxp,
    indrh, indwm, ndtot, nxf, mxf, ndf, ncst, iter,
    ncalls, itermx, indpu;
} simcom_;

struct {
    doublereal aux[125000];
} nauxcm_;

class vec {
    vector<double> c;

public:

    void resize(size_t n) {
        c.resize(n);
    }

    double& operator[](size_t i) {
        return c[i-1];
    }
};

struct aplcon {

    /* Subroutine */
    static int aplcon_(integer *nvar, integer *mcst) {
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
        /* primary value of NF */
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

    /* Subroutine */
    static int aploop_(doublereal *x, doublereal *vx, doublereal *f,
                       integer *iret) {

        /* Local variables */
        static integer nff;

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
        simcom_.indas = simcom_.indwm + simcom_.mxf + simcom_.nxf + simcom_.nxf;
        /* X result */
        simcom_.indpu = simcom_.indas + simcom_.nx;
        /* pulls, solution X and Vx */
        simcom_.ndtot = simcom_.indpu + simcom_.nx;
        /* total number of words (so far) */
        asteps_(&x[1], &vx[1], &nauxcm_.aux[simcom_.indst]);
        /*     __________________________________________________________________ */
        /*     internal APLOOP */
        /* initial steps ST(.) */

        static vec x_save;
        x_save.resize(simcom_.nx);
L10:

        *iret = -1;
        /* default status is -1 = continue */
        iploop_(&x[1], &vx[1], &f[1], x_save,
                &nauxcm_.aux[simcom_.inddx], &nauxcm_.aux[simcom_.indcf],
                &nauxcm_.aux[simcom_.indxp], &nauxcm_.aux[simcom_.indrh], iret);
        return 0;
    } /* aploop_ */

    /* Subroutine */
    static int iploop_(doublereal *x, doublereal *vx, doublereal *f,
                       vec& xs, doublereal *dx, doublereal *fcopy,
                       doublereal *xp, doublereal *rh, integer *iret) {
        /* System generated locals */
        integer i__1;
        doublereal d__1, d__2;

        /* Local variables */
        static integer j;
        static doublereal fj, fex[100];
        static integer nfit, jret;
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

    /* Subroutine */
    static int asteps_(doublereal *x, doublereal *vx, doublereal *st) {
        /* System generated locals */
        integer i__1, i__2;
        doublereal d__1, d__2, d__3, d__4, d__5;

        /* Local variables */
        static integer i__, j, ii;
        static doublereal vii;
        static integer ntine, ntrfl, ntvar;
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
            ntine = ntrfl / 1000 % 10;
            /* inequality flag */
            ii += i__;
            vii = abs(vx[ii]);
            /* original diagonal element */
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

    /* Subroutine */
    static int anumde_(doublereal *x, doublereal *f, doublereal *a,
                       doublereal *st, doublereal *xl, doublereal *fc,
                       doublereal *hh, integer *jret) {
        /* Initialized data */

        static logical tinue = false;

        /* System generated locals */
        integer i__1;
        doublereal d__1, d__2;

        /* Local variables */
        static integer i__, j, ij;
        static doublereal xd[2], xt[2], der;
        static integer nzer, nonz, ntine, ntvar, ntmes, ntlim, ntprf;
        static doublereal xsave;
        static doublereal ratdif, ratmax;

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
        tinue = true;
        i__ = 0;
        /* initialize derivative loop */
L10:
        if (i__ >= simcom_.nx) {
            /* finished */
            *jret = 0;
            /* ... means differentation finished */
            tinue = false;
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

        /*     __________________________________________________________________ */
        /*     define displaced values for derivative calculation */
        /* L20: */
        if (ntvar == 0 || ntvar == 2 || ntvar == 3) {
            xt[0] = xsave + st[i__];
            /* symmetric (two-sided) steps */
            xt[1] = xsave - st[i__];
            xd[0] = xt[0];
            xd[1] = xt[1];
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
        ratmax = 0.;
        /* max of diff-ratio */
        i__1 = simcom_.nf;
        for (j = 1; j <= i__1; ++j) {
            /* loop on all constraint functions */

            /* symmetric formula */
            der = (hh[j] - f[j]) / (xt[0] - xt[1]);
            /* !! internal variable */

            /*      _________________________________________________________________ */
            /*      classify derivative properties of variable I */
            if (der == 0.) {
                ++nzer;
                /* derivative zero - count */
            } else {
                /* derivative non-zero */
                ratdif = (d__2 = a[ij] - der, abs(d__2)) /
                         ((d__1 = a[ij], abs(d__1)) + abs(der));
                ratmax = max(ratmax, ratdif);
                ++nonz;
                /* count non-zero derivative */
            }
            a[ij] = der;
            /* insert into Jacobian matrix A */
            ij += simcom_.nx;
        }
        /*     packfl.inc   = code for flag packing */
        /*     explanation see: */
        /*     unpackfl.inc = code for flag unpacking */
        nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1] = (doublereal)(
                                                            ((((ntprf * 10 + ntlim) * 10 + ntine) * 10) * 10 + ntmes) * 10 +
                                                            ntvar);
        goto L10;
    } /* anumde_ */

    /* Subroutine */
    static int aniter_(doublereal *x, doublereal *vx, doublereal *f,
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
                                            integer *, integer *,
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
            if (ntvar == 2) {
                /* Poisson */
                /* Computing 2nd power */
                d__1 = x[i__];
                wm[ii] = -sqrt(d__1 * d__1 + 1.f);
                /* -MAX(ABS(X(I)),1.0D0) */
            }
        }
        duminv_(&a[1], &wm[1], &rh[1], &simcom_.nx, &simcom_.nf, &nrank,
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

    /* Subroutine */
    static int addtox_(doublereal *x, doublereal *xs, doublereal *dx,
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

    /* Subroutine */
    static int antest_(integer *iret) {
        /* System generated locals */
        doublereal d__1, d__2;

        /* Local variables */
        static integer i__;
        static doublereal cm[14];

        *iret = -1;
        /* calculate new Jacobian */
        /*     __________________________________________________________________ */
        /*     combined measure? */
        if (simcom_.iter == 1 && simcom_.ncst == 0) {
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
        if (simcom_.ncst < 2 && simcom_.iter > 1 &&
            simcom_.ftest > simcom_.ftestp * 2. + simcom_.epsf) {
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
            double dchisq = simcom_.chisq - simcom_.chsqp;
            if (abs(dchisq) <= simcom_.epschi && simcom_.ftest < simcom_.epsf) {
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

    /* Subroutine */
    static int acopxv_(doublereal *x, doublereal *vx, doublereal *dx,
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

    /* Subroutine */
    static int lesfcm_(doublereal *c__) {
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

    /* Subroutine */
    static int chndpv_(real *chi2, integer *nd, real *pval) {
        extern doublereal chprob_(doublereal *, integer *);

        *chi2 = simcom_.chisq;
        /* chi^square */
        *nd = simcom_.ndf;
        /* number of degrees of freedom */
        *pval = chprob_(&simcom_.chisq, &simcom_.ndf);
        /* p-value */
        return 0;
    } /* chndpv_ */

    /* Subroutine */
    static int adummy_0_(int n__, doublereal *arg, integer *i__, doublereal *step) {
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
        simcom_.epsf = *arg;
        /* |F| accuracy */
        return 0;
        /*     __________________________________________________________________ */

L_apstep:
        /* step size for numdif */
        if (*i__ < 1 || *i__ > simcom_.nx) {
            return 0;
        }
        simcom_.ipak = *i__;
        /*     unpackfl.inc = code for flag unpacking */
        ntrfl = (integer)nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1];
        /* get packed flags */
        ntvar = ntrfl % 10;
        /* transformation flag */
        ntmes = ntrfl / 10 % 10;
        /* M-estimate flag */
        ntine = ntrfl / 1000 % 10;
        /* inequality flag */
        ntlim = ntrfl / 10000 % 10;
        /* limit flag */
        ntprf = ntrfl / 100000 % 10;
        /* profile flag */
        nauxcm_.aux[simcom_.indst + *i__ - 1] = abs(*step);
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
        if (*i__ < 1 || *i__ > simcom_.nx) {
            return 0;
        }
        simcom_.ipak = *i__;
        /*     unpackfl.inc = code for flag unpacking */
        ntrfl = (integer)nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1];
        /* get packed flags */
        ntvar = ntrfl % 10;
        /* transformation flag */
        ntmes = ntrfl / 10 % 10;
        /* M-estimate flag */
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
        nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1] = (doublereal)(
                                                            ((((ntprf * 10 + ntlim) * 10 + ntine) * 10 + ntder) * 10 + ntmes) * 10 +
                                                            ntvar);
        return 0;
    } /* adummy_ */

    /* Subroutine */
    static int apdeps_(doublereal *arg) {
        return adummy_0_(2, arg, (integer *)0, (doublereal *)0);
    }

    /* Subroutine */
    static int apstep_(integer *i__, doublereal *step) {
        return adummy_0_(8, (doublereal *)0, i__, step);
    }

    /* Subroutine */
    static int apoiss_(integer *i__) {
        return adummy_0_(12, (doublereal *)0, i__, (doublereal*)0);
    }

    /* Subroutine */
    static int apstat_(doublereal *fopt, integer *nfun, integer *niter) {
        /*     __________________________________________________________________ */
        /*     return information after the fit */
        /*     __________________________________________________________________ */
        /* return Fopt and Nfun */
        *fopt = simcom_.chisq;
        *nfun = simcom_.ncalls;
        *niter = simcom_.iter;
        return 0;
    } /* apstat_ */

    /* Subroutine */
    static int appull_(doublereal *pulls) {
        /* System generated locals */
        integer i__1;

        /* Local variables */
        static integer i__;

        /* Parameter adjustments */
        --pulls;

        /* Function Body */
        i__1 = simcom_.nx;
        for (i__ = 1; i__ <= i__1; ++i__) {
            pulls[i__] = nauxcm_.aux[simcom_.indpu + i__ - 1];
        }
        return 0;
    } /* appull_ */

};