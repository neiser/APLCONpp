#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <sstream>

#include "chprob.f90.h"
#include "condutil.f90.h"
#include "helper.h"

using namespace std;


/* Common Block Declarations */

struct {
    double epsf, epschi, chisq, ftest, ftestp, chsqp, derfac,
    derufc, derlow, weight;
    int nx, nf, ipak, indst,
    indas, indtr, indfc, indhh, indxp,
    indrh, indwm, ndtot, nxf, mxf, ndf, ncst, iter,
    ncalls, itermx, indpu;
} simcom_;

struct {
    double aux[125000];
} nauxcm_;




struct aplcon {

    /* Subroutine */
    static int aplcon_(int *nvar, int *mcst) {
        /* Local variables */
        static int ij;

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
        simcom_.indst = 0;
        /* steps */
        simcom_.indtr = simcom_.indst + simcom_.nx;
        /* transformation flags */
        simcom_.ndtot = simcom_.indtr + simcom_.nx;
        /*     __________________________________________________________________ */
        /*     storage of initial sub-arrays */
        /*           1 ...    NX * (NF+2)   Jacobian, derivative matrix   A(.) */
        /*     INDST+1 ...    NX            steps                         ST(.) */
        /*     INDTR+1 ...    NX            properties of variables */
        /*     INDLM+1 ...    2 * NX        limits for variables          XL(2,.) */
        /*                    ----------- */
        /*     NDTOT =        NX * (NF+6)   initial memory area */
        /* space used so far */
        for (ij = 1; ij <= simcom_.ndtot; ++ij) {
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
    static int aploop_(double *x, double *vx, double *f,
                       int *iret) {
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

        simcom_.nxf = simcom_.nx + simcom_.nf;
        /* total number of fit equations */
        simcom_.mxf = (simcom_.nxf * simcom_.nxf + simcom_.nxf) / 2;
        /* number elements symmetric matrix */
        simcom_.indfc = simcom_.ndtot + simcom_.nx;
        /* pointer to FC(NF) = copy of F(NF) */
        simcom_.indhh = simcom_.indfc + simcom_.nf;
        /* step */
        simcom_.indxp = simcom_.indhh + simcom_.nx;
        /* previos step */
        simcom_.indrh = simcom_.indxp + simcom_.nx;
        /* right-hand side */
        simcom_.indwm = simcom_.indrh + simcom_.nxf;
        /* weight matrix */
        simcom_.indas = simcom_.indwm + simcom_.mxf;
        /* X result */
        simcom_.indpu = simcom_.indas + simcom_.nx;
        /* pulls, solution X and Vx */
        asteps_(&x[1], &vx[1], &nauxcm_.aux[simcom_.indst]);
        /*     __________________________________________________________________ */
        /*     internal APLOOP */
        /* initial steps ST(.) */


L10:

        *iret = -1;
        /* default status is -1 = continue */
        iploop_(&x[1], &vx[1], &f[1],
                &nauxcm_.aux[simcom_.indxp], &nauxcm_.aux[simcom_.indrh], iret);
        return 0;
    } /* aploop_ */

    /* Subroutine */
    static int iploop_(double *x, double *vx, double *f,
                       double *xp, double *rh, int *iret) {

        static vecd a;
        a.resize(simcom_.nx * simcom_.nf);

        static vecd xs;
        xs.resize(simcom_.nx);

        static vecd dx;
        dx.resize(simcom_.nx);

        static vecd fcopy;
        fcopy.resize(simcom_.nf);

        /* Local variables */
        static int j;
        static double fj;
        static int nfit, jret;
        static int istatu;

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
        --f;
        --vx;
        --x;


        // asdsa

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
        for (j = 1; j <= simcom_.nx; ++j) {
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
        /*     __________________________________________________________________ */
        /*     constraint test summary */
        if (istatu < 0) {
            goto L30;
        }
        simcom_.ftestp = simcom_.ftest;
        /* save previous value */
        simcom_.ftest = 0.;
        /* reset constraint tests */
        for (j = 1; j <= simcom_.nf; ++j) {
            fj = f[j];
            fcopy[j] = fj;
            /* copy constraint vector */
            simcom_.ftest += abs(fj);
            /* sum absolute values */
        }
        simcom_.ftest = max(1e-16, simcom_.ftest / (double)simcom_.nf);
        /* average |F| */
        if (istatu == 1) {
            goto L60;
        }
        /*     __________________________________________________________________ */
        /*     start numerical derivatives */
L30:
        istatu = -1;
        /*     __________________________________________________________________ */
        /*     derivative calculation */
        anumde_(&x[1], &f[1], a, &nauxcm_.aux[simcom_.indst],
                &nauxcm_.aux[simcom_.indfc],
                &nauxcm_.aux[simcom_.indhh], &jret);
        /* derivative matrix A */
        /* steps  ST(.) */
        /* copy FC(.) central F(.) */
        /* copy HH(.) shifted F(.) */
        *iret = -1;
        if (jret < 0) {
            return 0;
        }
        /*     __________________________________________________________________ */
        /*     next iteration */
        /* ...for constraint calculation */
        aniter_(&x[1], &vx[1], &fcopy[1], a, &xp[1], &rh[1],
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
        acopxv_(&x[1], &vx[1], dx,
                &nauxcm_.aux[simcom_.indas], &nauxcm_.aux[simcom_.indwm],
                &nauxcm_.aux[simcom_.indpu]);
        *iret = 0;
        return 0;
    } /* iploop_ */

    /* Subroutine */
    static int asteps_(double *x, double *vx, double *st) {
        /* Local variables */
        static int i, j, ii;
        static double vii;
        static int ntine, ntrfl, ntvar;

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
        for (i = 1; i <= simcom_.nx; ++i) {
            /* loop on all variables */
            simcom_.ipak = i;
            /*     unpackfl.inc = code for flag unpacking */
            ntrfl = (int)nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1];
            /* get packed flags */
            ntvar = ntrfl % 10;
            /* transformation flag */
            ntine = ntrfl / 1000 % 10;
            /* inequality flag */
            ii += i;
            vii = abs(vx[ii]);
            /* original diagonal element */
            /*      _________________________________________________________________ */
            /*      define step size for derivative calculation */
            if (vii != 0.) {
                /* measured variable */
                if (st[i] > 0.) {
                    /* Computing MIN */
                    st[i] = min(st[i], simcom_.derfac * sqrt(vii));
                    /* user step, if smaller */
                } else if (st[i] <= 0.) {
                    /* step is undefined */
                    st[i] = simcom_.derfac * sqrt(vii);
                    /* step from cov matrix */
                    if (ntvar != 2) {
                        st[i] = min(st[i], simcom_.derlow * max(1e-6, abs(x[i])));
                    } else {
                        /* Poisson */
                        st[i] = min(st[i], simcom_.derlow * max(1e-6, abs(x[i] + 1.)));
                    }
                }
            } else if (vii == 0.) {
                /* unmeasured variable */
                --simcom_.ndf;
                /* reduce degrees of freedom */
                for (j = 1; j <= simcom_.nx; ++j) {
                    vx[ijsym_(&i, &j)] = 0.;
                    /* clear matrix elements */
                }
                /* Computing MAX */
                st[i] = simcom_.derufc * max(1.0, abs(x[i]));
            }
            if (ntine == 1) {
                st[i] = 0.;
            }
        }
        return 0;
    } /* asteps_ */

    /* Subroutine */
    static int anumde_(double *x, double *f, vecd& a,
                       double *st, double *fc,
                       double *hh, int *jret) {
        /* Initialized data */
        static bool tinue = false;

        /* Local variables */
        static int i, j, ij;
        static double xd[2], xt[2], der;
        static int nzer, nonz, ntvar;
        static double xsave;
        static double ratdif, ratmax;

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
        --st;
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
        i = 0;
        /* initialize derivative loop */
L10:
        if (i >= simcom_.nx) {
            /* finished */
            *jret = 0;
            /* ... means differentation finished */
            tinue = false;
            /* Jacobian ready */
            return 0;
        }
        ++i;
        /* next variable */
        if (st[i] == 0.f) {
            goto L10;
        }
        /* skip fixed variable */
        /* skip repeated derivative calculation */
        xsave = x[i];

        /*     __________________________________________________________________ */
        /*     define displaced values for derivative calculation */
        /* L20: */
        if (ntvar == 0 || ntvar == 2 || ntvar == 3) {
            xt[0] = xsave + st[i];
            /* symmetric (two-sided) steps */
            xt[1] = xsave - st[i];
            xd[0] = xt[0];
            xd[1] = xt[1];
        }
        /*     __________________________________________________________________ */
        /*     set variable to displaced value and return for calculation */
        x[i] = xd[0];
        /* first step */
        i = -i;
        return 0;
        /*     __________________________________________________________________ */
        /*     continue */
L30:
        if (i < 0) {
            /* calculation of first step done ... */
            for (j = 1; j <= simcom_.nf; ++j) {
                hh[j] = f[j];
                /* save constraint values */
            }
            i = -i;
            /* reverse flag */
            x[i] = xd[1];
            /* set next step ... */
            return 0;
            /* ... and return for second step */
        }
        /*     __________________________________________________________________ */
        /*     INIT ne 0: second step done - calculate derivative */
        x[i] = xsave;
        /* restore variable I */
        ij = i;
        /* derivative calculation */
        nzer = 0;
        /* reset: number of zero derivatives */
        nonz = 0;
        /* reset: number of non-zero derivatives */
        ratmax = 0.;
        /* max of diff-ratio */
        for (j = 1; j <= simcom_.nf; ++j) {
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
                ratdif = abs(a[ij] - der) / (abs(a[ij]) + abs(der));
                ratmax = max(ratmax, ratdif);
                ++nonz;
                /* count non-zero derivative */
            }
            a[ij] = der;
            /* insert into Jacobian matrix A */
            ij += simcom_.nx;
        }
        goto L10;
    } /* anumde_ */

    /* Subroutine */
    static int aniter_(double *x, double *vx, double *f,
                       vecd& a, double *xp, double *rh,
                       double *wm, double *dx) {
        /* Local variables */
        static int i, j, ia, ii;
        static int ntrfl, ntvar;

        /* next iteration step */

        /* Parameter adjustments */
        --dx;
        --wm;
        --rh;
        --xp;
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
        for (i = 1; i <= simcom_.nx; ++i) {
            /* first NX components */
            rh[i] = 0.;
            /* define right hand side of equation */
        }
        for (j = 1; j <= simcom_.nf; ++j) {
            /* next NF components */
            nauxcm_.aux[simcom_.indfc + j - 1] = f[j];
            rh[simcom_.nx + j] = -f[j];
        }
        ia = 0;
        for (j = 1; j <= simcom_.nf; ++j) {
            /* "subtract" actual step */
            rh[simcom_.nx + j] += scalxy_(&a[ia + 1], &dx[1], &simcom_.nx);
            ia += simcom_.nx;
            nauxcm_.aux[simcom_.indhh + j - 1] = rh[simcom_.nx + j];
            /* right hand side for chi**2 */
        }
        /*     __________________________________________________________________ */
        /*     form matrix and solve */
        for (i = 1; i <= (simcom_.nx * simcom_.nx + simcom_.nx) / 2; ++i) {
            wm[i] = -vx[i];
            /* copy -VX(.) into W_11 */
        }
        ii = 0;
        /* modify V for Poisson variables */
        for (i = 1; i <= simcom_.nx; ++i) {
            ii += i;
            simcom_.ipak = i;
            /*     unpackfl.inc = code for flag unpacking */
            ntrfl = (int)nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1];
            /* get packed flags */
            ntvar = ntrfl % 10;
            /* transformation flag */
            if (ntvar == 2) {
                /* Poisson */
                /* Computing 2nd power */
                wm[ii] = -max(abs(x[i]), 1.0);
            }
        }
        duminv_(&a[1], &wm[1], &rh[1], simcom_.nx, simcom_.nf);
        simcom_.chisq = -scalxy_(&nauxcm_.aux[simcom_.indhh], &rh[simcom_.nx + 1], &simcom_.nf);
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
        for (i = 1; i <= simcom_.nx; ++i) {
            xp[i] = dx[i];
            /* save previous corrections */
            dx[i] = rh[i];
            /* store new corrections */
        }
        return 0;
    } /* aniter_ */

    /* Subroutine */
    static int addtox_(double *x, double *xs, double *dx,
                       double *xp) {
        /* Parameter adjustments */
        --xp;
        --dx;
        --xs;
        --x;

        /* Function Body */
        for (int i = 1; i <= simcom_.nx; ++i) {
            dx[i] = simcom_.weight * dx[i] + (1. - simcom_.weight) * xp[i];
            /* reduce step evtl. */
            x[i] = xs[i] + dx[i];
            /* correct x and return to test constraints */
        }
        return 0;
    } /* addtox_ */

    /* Subroutine */
    static int antest_(int *iret) {
        *iret = -1;
        /* combined penalty */
        //        }
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
    static int acopxv_(double *x, double *vx, vecd& dx,
                       double *as, double *wm, double *pu) {

        /* Local variables */
        static int i, ii;
        static double scopy;

        /*     __________________________________________________________________ */
        /*     convergence: pull calculation */
        /* Parameter adjustments */
        --pu;
        --wm;
        --as;
        --vx;
        --x;

        /* Function Body */
        ii = 0;
        for (i = 1; i <= simcom_.nx; ++i) {
            as[i] = x[i];
            ii += i;
            pu[i] = 0.f;
            if (vx[ii] > 0.f) {
                if (vx[ii] - wm[ii] > 0.f) {
                    pu[i] = dx[i] / sqrt(vx[ii] - wm[ii]);
                }
            }
        }
        /*     __________________________________________________________________ */
        /*     copy/exchange result/input covariance matrix */
        for (i = 1; i <= (simcom_.nx * simcom_.nx + simcom_.nx) / 2; ++i) {
            scopy = vx[i];
            vx[i] = wm[i];
            /* copy fitted covariance matrix */
            wm[i] = scopy;
            /* ... and save input matrix */
        }
        return 0;
    } /* acopxv_ */

    /* Subroutine */
    static int chndpv_(double *chi2, int *nd, double *pval) {

        *chi2 = simcom_.chisq;
        /* chi^square */
        *nd = simcom_.ndf;
        /* number of degrees of freedom */
        *pval = chprob_(simcom_.chisq, simcom_.ndf);
        /* p-value */
        return 0;
    } /* chndpv_ */

    /* Subroutine */
    static int adummy_0_(int n__, double *arg, int *i, double *step) {
        static int ntder, ntine, ntrfl, ntvar, ntmes, ntlim, ntprf;
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
        if (*i < 1 || *i > simcom_.nx) {
            return 0;
        }
        simcom_.ipak = *i;
        /*     unpackfl.inc = code for flag unpacking */
        ntrfl = (int)nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1];
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
        nauxcm_.aux[simcom_.indst + *i - 1] = abs(*step);
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
        if (*i < 1 || *i > simcom_.nx) {
            return 0;
        }
        simcom_.ipak = *i;
        /*     unpackfl.inc = code for flag unpacking */
        ntrfl = (int)nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1];
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
        nauxcm_.aux[simcom_.indtr + simcom_.ipak - 1] = (double)(
                                                            ((((ntprf * 10 + ntlim) * 10 + ntine) * 10 + ntder) * 10 + ntmes) * 10 +
                                                            ntvar);
        return 0;
    } /* adummy_ */

    /* Subroutine */
    static int apdeps_(double *arg) {
        return adummy_0_(2, arg, (int *)0, (double *)0);
    }

    /* Subroutine */
    static int apstep_(int *i, double *step) {
        return adummy_0_(8, (double *)0, i, step);
    }

    /* Subroutine */
    static int apoiss_(int *i) {
        return adummy_0_(12, (double *)0, i, (double*)0);
    }

    /* Subroutine */
    static int apstat_(double *fopt, int *nfun, int *niter) {
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
    static int appull_(double *pulls) {
        /* Local variables */
        static int i;

        /* Parameter adjustments */
        --pulls;

        /* Function Body */
        for (i = 1; i <= simcom_.nx; ++i) {
            pulls[i] = nauxcm_.aux[simcom_.indpu + i - 1];
        }
        return 0;
    } /* appull_ */

};