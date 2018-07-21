#include <cmath>

typedef double doublereal;
typedef int integer;

using namespace std;

inline doublereal dgamml_(doublereal x) {
    /* Initialized data */

    static doublereal cof[6] = {76.18009172947146,   -86.50532032941677,
                                24.01409824083091,   -1.231739572450155,
                                .001208650973866179, -5.395239384953e-7};
    static doublereal stp = 2.5066282746310005;

    /* Local variables */
    static integer j;
    static doublereal yy, ser, tmp;

    /* ln[Gamma(x)] */
    /*     ... */
    yy = x;
    tmp = (x + .5) * log(x + 5.5) - x - 5.5;
    ser = 1.000000000190015;
    for (j = 1; j <= 6; ++j) {
        yy += 1.;
        ser += cof[j - 1] / yy;
    }
    return tmp + log(stp * ser / x);
} /* dgamml_ */

inline doublereal dgamin_(doublereal a, doublereal x) {
    /* Local variables */
    static doublereal b, c__, d__, h__;
    static integer i__, n;
    static doublereal an, ap, del, gln, sum;

    /*     incomplete gamma function P(a,x) */
    /*     returns -1.0 for x < 0 or a =< 0 */
    /*     ... */
    if (x < 0. || a <= 0.) {
        return -1.;
        /* error */
    } else if (x == 0.f) {
        return 0.;
    } else if (x < a + 1.) {
        /* series representation */
        gln = dgamml_(a);
        /* ln[Gamma(a)] */
        ap = a;
        sum = 1. / a;
        del = sum;
        for (n = 1; n <= 300; ++n) {
            ap += 1.f;
            del = del * x / ap;
            sum += del;
            if (abs(del) < abs(sum) * 3e-7) {
                goto L10;
            }
        }
L10:
        return sum * exp(-(x) + a * log(x) - gln);
    } else {
        /* continued fraction representation */
        gln = dgamml_(a);
        /* ln[Gamma(a)] */
        b = x + 1. - a;
        c__ = 9.9999999999999988e29;
        d__ = 1. / b;
        h__ = d__;
        for (i__ = 1; i__ <= 300; ++i__) {
            an = -i__ * (i__ - a);
            b += 2.;
            d__ = an * d__ + b;
            if (abs(d__) < 1e-30) {
                d__ = 1e-30;
            }
            c__ = b + an / c__;
            if (abs(c__) < 1e-30) {
                c__ = 1e-30;
            }
            d__ = 1. / d__;
            del = d__ * c__;
            h__ *= del;
            if (abs(del - 1.) < 3e-7) {
                goto L20;
            }
        }
        return 0;
L20:
        return 1. - exp(-(x) + a * log(x) - gln) * h__;
    }
} /* dgamin_ */

inline doublereal chprob_(doublereal *chisq, integer *n) {
    /*     chi square probability for N degrees of freedom at CHISQ */
    /*     =integral from 0 ... CHISQ */
    /* prob from chisquare */
    /*     ... */
    if (*chisq <= 0.f) {
        return 1.;
    } else {
        return 1. - dgamin_(*n * .5, *chisq * .5);
    }
} /* chprob_ */




